! Bandwidth microbenchmark for OpenMP target map(to:)/map(from:)/map(tofrom:)
! data transfers, to measure the raw host<->device transfer cost on this
! machine's integrated GPU (gfx1103, Radeon 780M in a Ryzen 8845HS APU).
! The iGPU shares system RAM with the host and has no XNACK, so map()
! always does a real copy rather than being a no-op or a page-fault-driven
! migration — this measures that copy.
!
! Five patterns, each swept over array sizes from 4 KiB to 256 MiB:
!
!   map(to)        target region with map(to: a), touching a(1);
!                  device buffer allocated and freed every call (cold)
!   map(from)      as above with map(from: a)/map(delete:) (cold)
!   map(tofrom)    as above with map(tofrom: a) — both directions cross
!                  the bus, so its GB/s is computed against 2*bytes
!   update(to)     device buffer allocated once outside the timed loop;
!                  timed loop is only `target update to(a)` (pure H2D
!                  copy, no alloc/dealloc)
!   update(from)   as above with `target update from(a)` (pure D2H copy)
!
! There's also a kernel-launch-overhead measurement (a single already-
! resident element) to give a floor latency: at small sizes the reported
! bandwidth is dominated by this fixed cost, not by the interconnect.
!
! Each timed kernel dispatch lives in its own subroutine, called from a
! timing loop one level up (never a `!$omp target` construct written
! directly inside a host `do` loop). Writing it the other way — the
! construct lexically inside the timed loop — compiled fine here but
! crashed at run time on every call: amdflang's device linker silently
! dropped the kernel's device-side symbol while leaving the host-side
! dispatch call in place, producing "requested object was not found in
! the binary image". Every other test in this directory (av_rem_test.F90,
! col_norm_test.F90, ...) already puts the target construct inside a
! plain subroutine called from the timing loop, not the other way round;
! this file follows the same shape.
!
! Bandwidths are reported as best-of-n_runs (i.e. from the minimum time),
! in GB/s (10^9 bytes/s).

module bandwidth_mod
  use test_utils_mod, only: dp
  use omp_lib, only: omp_get_wtime
  use iso_c_binding, only: c_int
  implicit none

  ! omp_lib's own omp_is_initial_device() reports wrong (host) results here
  ! when called from inside a teams/parallel-do construct on this compiler,
  ! even though LIBOMPTARGET_INFO=16 confirms the kernel really launched on
  ! the device -- see dc_exit_test.F90's identical workaround in this same
  ! directory for the reproducer (there it's declared `pure` to satisfy
  ! `do concurrent`'s purity rule; that's not needed here, only ordinary
  ! `do` loops are used, so the declaration is left impure to match
  ! omp_lib's own, avoiding an interface-mismatch warning).
  interface
    integer(c_int) function is_initial_device() bind(C, name="omp_is_initial_device")
      import :: c_int
    end function is_initial_device
  end interface

contains

  subroutine dispatch_check(flag)
    integer, intent(out) :: flag
    integer :: i
    flag = 1
    !$omp target teams distribute parallel do map(tofrom: flag)
    do i = 1, 1
      flag = is_initial_device()
    enddo
  end subroutine dispatch_check

  ! Positive check that target regions actually dispatch to the device
  ! rather than silently falling back to the host — see the omp_patterns
  ! README for why a passing run is not otherwise proof of that.
  subroutine check_device_dispatch()
    integer :: flag
    call dispatch_check(flag)
    if (flag == 0) then
      write(*,'(A)') 'Device check: target region executed on the GPU.'
    else
      write(*,'(A)') &
        'WARNING: target region executed on the HOST — the numbers below are not GPU transfers.'
    endif
  end subroutine check_device_dispatch

  function bandwidth_gbs(nbytes, seconds) result(gbs)
    integer(8), intent(in) :: nbytes
    real(dp),   intent(in) :: seconds
    real(dp) :: gbs
    gbs = real(nbytes, dp) / seconds / 1.0e9_dp
  end function bandwidth_gbs

  subroutine dispatch_launch_only(x)
    real(dp), intent(inout) :: x
    integer :: i
    !$omp target teams distribute parallel do map(tofrom: x)
    do i = 1, 1
      x = x + 1.0_dp
    enddo
  end subroutine dispatch_launch_only

  subroutine time_launch_only(times)
    real(dp), intent(out) :: times(:)
    real(dp) :: x, t0, t1
    integer  :: irun

    x = 0.0_dp
    !$omp target enter data map(to: x)

    do irun = 1, size(times)
      t0 = omp_get_wtime()
      call dispatch_launch_only(x)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo

    !$omp target exit data map(delete: x)
  end subroutine time_launch_only

  subroutine dispatch_map_to(n, a)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: a(n)
    integer :: i
    !$omp target teams distribute parallel do map(to: a)
    do i = 1, 1
      a(1) = a(1) + 1.0_dp
    enddo
  end subroutine dispatch_map_to

  subroutine time_map_to(n, times)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: times(:)
    real(dp), allocatable :: a(:)
    real(dp) :: t0, t1
    integer  :: irun

    allocate(a(n))
    a = 1.0_dp

    do irun = 1, size(times)
      t0 = omp_get_wtime()
      call dispatch_map_to(n, a)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo

    deallocate(a)
  end subroutine time_map_to

  subroutine dispatch_map_from(n, a)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: a(n)
    integer :: i
    !$omp target teams distribute parallel do map(from: a)
    do i = 1, 1
      a(1) = 1.0_dp
    enddo
  end subroutine dispatch_map_from

  subroutine time_map_from(n, times)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: times(:)
    real(dp), allocatable :: a(:)
    real(dp) :: t0, t1
    integer  :: irun

    allocate(a(n))

    do irun = 1, size(times)
      t0 = omp_get_wtime()
      call dispatch_map_from(n, a)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo

    deallocate(a)
  end subroutine time_map_from

  subroutine dispatch_map_tofrom(n, a)
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: a(n)
    integer :: i
    !$omp target teams distribute parallel do map(tofrom: a)
    do i = 1, 1
      a(1) = a(1) + 1.0_dp
    enddo
  end subroutine dispatch_map_tofrom

  subroutine time_map_tofrom(n, times)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: times(:)
    real(dp), allocatable :: a(:)
    real(dp) :: t0, t1
    integer  :: irun

    allocate(a(n))
    a = 1.0_dp

    do irun = 1, size(times)
      t0 = omp_get_wtime()
      call dispatch_map_tofrom(n, a)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo

    deallocate(a)
  end subroutine time_map_tofrom

  subroutine time_update_to(n, times)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: times(:)
    real(dp), allocatable :: a(:)
    real(dp) :: t0, t1
    integer  :: irun

    allocate(a(n))
    a = 1.0_dp
    !$omp target enter data map(to: a)

    do irun = 1, size(times)
      t0 = omp_get_wtime()
      !$omp target update to(a)
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo

    !$omp target exit data map(delete: a)
    deallocate(a)
  end subroutine time_update_to

  subroutine time_update_from(n, times)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: times(:)
    real(dp), allocatable :: a(:)
    real(dp) :: t0, t1
    integer  :: irun

    allocate(a(n))
    a = 1.0_dp
    !$omp target enter data map(to: a)

    do irun = 1, size(times)
      t0 = omp_get_wtime()
      !$omp target update from(a)
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo

    !$omp target exit data map(delete: a)
    deallocate(a)
  end subroutine time_update_from

end module bandwidth_mod


program test_bandwidth
  use bandwidth_mod
  use test_utils_mod
  implicit none

  integer,    parameter :: n_sizes = 9
  integer(8), parameter :: byte_sizes(n_sizes) = &
    [4096_8, 16384_8, 65536_8, 262144_8, 1048576_8, 4194304_8, &
     16777216_8, 67108864_8, 268435456_8]

  real(dp) :: times(n_runs)
  real(dp) :: bw_to, bw_from, bw_tofrom, bw_upd_to, bw_upd_from
  integer(8) :: nbytes
  integer :: n, isize

  write(*,'(A)') 'OpenMP target map bandwidth microbenchmark'
  write(*,'(A)') '(integrated GPU, see README.md for the device this was run on)'
  write(*,*)

  call check_device_dispatch()
  write(*,*)

  call time_launch_only(times)
  write(*,'(A,I0,A,F10.3,A)') 'Kernel launch overhead (single resident element, best of ', &
       n_runs, ' runs): ', minval(times) * 1.0e6_dp, ' us'
  write(*,*)

  write(*,'(A,I0,A)') 'All bandwidths are best-of-', n_runs, ' runs, in GB/s (10^9 bytes/s).'
  write(*,'(A)') 'map(tofrom) bandwidth is computed against 2*bytes, since it crosses the bus both ways.'
  write(*,*)
  write(*,'(A)') &
    '       bytes    map(to)   map(from)  map(tofrom)  update(to) update(from)'

  do isize = 1, n_sizes
    nbytes = byte_sizes(isize)
    n = int(nbytes / 8_8)

    call time_map_to(n, times)
    bw_to = bandwidth_gbs(nbytes, minval(times))

    call time_map_from(n, times)
    bw_from = bandwidth_gbs(nbytes, minval(times))

    call time_map_tofrom(n, times)
    bw_tofrom = bandwidth_gbs(2_8 * nbytes, minval(times))

    call time_update_to(n, times)
    bw_upd_to = bandwidth_gbs(nbytes, minval(times))

    call time_update_from(n, times)
    bw_upd_from = bandwidth_gbs(nbytes, minval(times))

    write(*,'(I12,5F12.3)') nbytes, bw_to, bw_from, bw_tofrom, bw_upd_to, bw_upd_from
  enddo

  write(*,*)
  write(*,*) 'Done.'

end program test_bandwidth
