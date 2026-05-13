! Test for GPU/CPU bitwise reproducibility of the av_rem_v accumulation
! pattern from MOM_barotropic.F90 (lines ~1595-1604):
!
!   do concurrent (J=js-1:je)
!     do concurrent (i=is-1:ie+1)
!       av_rem_v(i,J) = 0.0
!     enddo
!     do k=1,nz
!       do concurrent (i=is:ie)
!         av_rem_v(i,J) = av_rem_v(i,J) + frhatv(i,J,k) * visc_rem_v(i,J,k)
!       enddo
!     enddo
!   enddo
!
! Five variants:
!   run_av_rem_omp    – !$omp target teams distribute over J, k serial in
!                       team, parallel do over i inside each k
!   run_av_rem_omp_ji – !$omp target teams distribute parallel do over (J,I)
!                       using nteams, serial k within each thread
!   run_av_rem_dc     – do concurrent, k serial in outer-j body, i concurrent
!                       per k; GPU offload via compiler flag
!   run_av_rem_dc_ji  – fused do concurrent(j,i) for init and accumulation,
!                       serial k within each (j,i) iteration
!   run_av_rem_cpu    – plain do loops, CPU reference
!
! Comparisons: each GPU/DC variant vs CPU.
! Timing: multiple problem sizes (32..256), 5 timed runs each.

#include "omp_macros.inc"

module av_rem_mod
  use test_utils_mod, only: dp
  implicit none

contains

  ! ------------------------------------------------------------------ !
  !  OpenMP target offload.                                             !
  !  Outer J distributed across teams (thread_limit(1) → serial inner  !
  !  loops per team), mirroring the do concurrent (J) outer structure. !
  ! ------------------------------------------------------------------ !
  subroutine run_av_rem_omp(nx, ny, nz, frhatv, visc_rem_v, av_rem_v)
    integer,  intent(in)  :: nx, ny, nz
    real(dp), intent(in)  :: frhatv(nx, ny, nz)
    real(dp), intent(in)  :: visc_rem_v(nx, ny, nz)
    real(dp), intent(out) :: av_rem_v(nx, ny)
    integer :: i, j, k, nthreads

    nthreads = int((nx + WAVEFRONT-1)/WAVEFRONT) * WAVEFRONT

    !$omp target teams TEAMS_OUTER_LOOP num_teams(ny) thread_limit(nthreads) &
    !$omp&   map(to: frhatv, visc_rem_v) map(from: av_rem_v)
    do j = 1, ny
      !$omp PARALLEL_INNER_LOOP
      do i = 1, nx
        av_rem_v(i,j) = 0.0_dp
      enddo
      do k = 1, nz
        !$omp PARALLEL_INNER_LOOP
        do i = 1, nx
          av_rem_v(i,j) = av_rem_v(i,j) + frhatv(i,j,k) * visc_rem_v(i,j,k)
        enddo
      enddo
    enddo
  end subroutine run_av_rem_omp

  ! ------------------------------------------------------------------ !
  !  do concurrent — GPU offload controlled by compiler flag only;     !
  !  source is unchanged.  Mirrors the selection exactly:              !
  !    outer do concurrent(J), inner do concurrent(i) init,           !
  !    then do k { do concurrent(i) accum }.                          !
  ! ------------------------------------------------------------------ !
  subroutine run_av_rem_dc(nx, ny, nz, frhatv, visc_rem_v, av_rem_v)
    integer,  intent(in)  :: nx, ny, nz
    real(dp), intent(in)  :: frhatv(nx, ny, nz)
    real(dp), intent(in)  :: visc_rem_v(nx, ny, nz)
    real(dp), intent(out) :: av_rem_v(nx, ny)
    integer :: i, j, k

    do concurrent (j = 1:ny)
      do concurrent (i = 1:nx)
        av_rem_v(i,j) = 0.0_dp
      enddo
      do k = 1, nz
        do concurrent (i = 1:nx)
          av_rem_v(i,j) = av_rem_v(i,j) + frhatv(i,j,k) * visc_rem_v(i,j,k)
        enddo
      enddo
    enddo
  end subroutine run_av_rem_dc

  ! ------------------------------------------------------------------ !
  !  OpenMP ji variant: j distributed (one team per j-row), i         !
  !  parallelised within the team, k serial within each i-thread.     !
  !  Loop order: distribute(j) → parallel(i) → serial k.             !
  ! ------------------------------------------------------------------ !
  subroutine run_av_rem_omp_ji(nx, ny, nz, nteams, frhatv, visc_rem_v, av_rem_v)
    integer,  intent(in)  :: nx, ny, nz, nteams
    real(dp), intent(in)  :: frhatv(nx, ny, nz)
    real(dp), intent(in)  :: visc_rem_v(nx, ny, nz)
    real(dp), intent(out) :: av_rem_v(nx, ny)
    integer :: i, j, k

    !$omp target teams COMBINED_LOOP collapse(2) num_teams(nteams) &
    !$omp&   map(to: frhatv, visc_rem_v) map(from: av_rem_v)
    do j = 1, ny
      do i = 1, nx
        av_rem_v(i,j) = 0.0_dp
        do k = 1, nz
          av_rem_v(i,j) = av_rem_v(i,j) + frhatv(i,j,k) * visc_rem_v(i,j,k)
        enddo
      enddo
    enddo
  end subroutine run_av_rem_omp_ji

  ! ------------------------------------------------------------------ !
  !  do concurrent ji variant: fused (j,i) concurrent for init and    !
  !  accumulation, serial k within each (j,i) iteration.              !
  !  Loop order: concurrent(j,i) → serial k.                         !
  ! ------------------------------------------------------------------ !
  subroutine run_av_rem_dc_ji(nx, ny, nz, frhatv, visc_rem_v, av_rem_v)
    integer,  intent(in)  :: nx, ny, nz
    real(dp), intent(in)  :: frhatv(nx, ny, nz)
    real(dp), intent(in)  :: visc_rem_v(nx, ny, nz)
    real(dp), intent(out) :: av_rem_v(nx, ny)
    integer :: i, j, k

    do concurrent (j = 1:ny, i = 1:nx)
      av_rem_v(i,j) = 0.0_dp
      do k = 1, nz
        av_rem_v(i,j) = av_rem_v(i,j) + frhatv(i,j,k) * visc_rem_v(i,j,k)
      enddo
    enddo
  end subroutine run_av_rem_dc_ji

  ! ------------------------------------------------------------------ !
  !  CPU reference: plain do loops, identical logic, no directives.   !
  ! ------------------------------------------------------------------ !
  subroutine run_av_rem_cpu(nx, ny, nz, frhatv, visc_rem_v, av_rem_v)
    integer,  intent(in)  :: nx, ny, nz
    real(dp), intent(in)  :: frhatv(nx, ny, nz)
    real(dp), intent(in)  :: visc_rem_v(nx, ny, nz)
    real(dp), intent(out) :: av_rem_v(nx, ny)
    integer :: i, j, k

    do j = 1, ny
      do i = 1, nx
        av_rem_v(i,j) = 0.0_dp
      enddo
      do k = 1, nz
        do i = 1, nx
          av_rem_v(i,j) = av_rem_v(i,j) + frhatv(i,j,k) * visc_rem_v(i,j,k)
        enddo
      enddo
    enddo
  end subroutine run_av_rem_cpu

end module av_rem_mod


program test_av_rem
  use av_rem_mod
  use test_utils_mod
  use omp_lib
  implicit none

  integer,  parameter :: nz           = 100
  integer,  parameter :: n_sizes      = 5
  integer,  parameter :: all_sizes(n_sizes) = [32, 64, 128, 256, 512]

  real(dp), allocatable :: frhatv(:,:,:)
  real(dp), allocatable :: visc_rem_v(:,:,:)
  real(dp), allocatable :: av_rem_omp(:,:)
  real(dp), allocatable :: av_rem_omp_ji(:,:)
  real(dp), allocatable :: av_rem_dc(:,:)
  real(dp), allocatable :: av_rem_dc_ji(:,:)
  real(dp), allocatable :: av_rem_cpu(:,:)

  integer  :: nx, ny, nteams, i, j, k, isize, irun
  real(dp) :: t0, t1
  real(dp) :: times(n_runs)
  logical  :: p1, p2, p3, p4, all_pass

  do isize = 1, n_sizes
    nx = all_sizes(isize) ; ny = all_sizes(isize)
    nteams = (nx * ny + TEAM_SIZE-1) / TEAM_SIZE

    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      'Size: nx=', nx, '  ny=', ny, '  nz=', nz, '  nteams=', nteams
    write(*,'(A)') '========================================================'

    allocate(frhatv(nx, ny, nz), visc_rem_v(nx, ny, nz))
    allocate(av_rem_omp(nx, ny), av_rem_omp_ji(nx, ny))
    allocate(av_rem_dc(nx, ny),  av_rem_dc_ji(nx, ny))
    allocate(av_rem_cpu(nx, ny))

    do k = 1, nz ; do j = 1, ny ; do i = 1, nx
      frhatv(i,j,k)     = 0.5_dp + 0.5_dp * sin(pi * real(i + j + k,     dp) / 20.0_dp)
      visc_rem_v(i,j,k) = 0.8_dp + 0.2_dp * cos(pi * real(2*i + j + k, dp) / 25.0_dp)
    enddo ; enddo ; enddo

    ! -------------------------------------------------------------- !
    !  Correctness (one run per variant)                              !
    ! -------------------------------------------------------------- !
    write(*,*)
    write(*,*) '--- Correctness ---'

    call run_av_rem_cpu(nx, ny, nz, frhatv, visc_rem_v, av_rem_cpu)

    call run_av_rem_omp(nx, ny, nz, frhatv, visc_rem_v, av_rem_omp)
    !$omp taskwait
    call compare_2d('omp    |cpu', av_rem_omp,    av_rem_cpu, p1)

    call run_av_rem_omp_ji(nx, ny, nz, nteams, frhatv, visc_rem_v, av_rem_omp_ji)
    !$omp taskwait
    call compare_2d('omp_ji |cpu', av_rem_omp_ji, av_rem_cpu, p2)

    call run_av_rem_dc(nx, ny, nz, frhatv, visc_rem_v, av_rem_dc)
    !$omp taskwait
    call compare_2d('dc     |cpu', av_rem_dc,     av_rem_cpu, p3)

    call run_av_rem_dc_ji(nx, ny, nz, frhatv, visc_rem_v, av_rem_dc_ji)
    !$omp taskwait
    call compare_2d('dc_ji  |cpu', av_rem_dc_ji,  av_rem_cpu, p4)

    all_pass = p1 .and. p2 .and. p3 .and. p4
    if (all_pass) then
      write(*,*) '  All correct.'
    else
      write(*,*) '  FAILURES — see above.'
    endif

    ! -------------------------------------------------------------- !
    !  Timings (n_runs runs per variant, all times printed)           !
    ! -------------------------------------------------------------- !
    write(*,*)
    write(*,*) '--- Timings ---'

    !$omp target enter data map(to: frhatv, visc_rem_v) &
    !$omp&   map(alloc: av_rem_omp, av_rem_omp_ji, av_rem_dc, av_rem_dc_ji)

    write(*,*) 'OMP (k serial in team, i parallel per k):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_av_rem_omp(nx, ny, nz, frhatv, visc_rem_v, av_rem_omp)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'OMP ji (i parallel in team, k serial per thread):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_av_rem_omp_ji(nx, ny, nz, nteams, frhatv, visc_rem_v, av_rem_omp_ji)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'do concurrent (k serial in team, i concurrent per k):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_av_rem_dc(nx, ny, nz, frhatv, visc_rem_v, av_rem_dc)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'do concurrent ji (i concurrent in outer j, k serial per i):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_av_rem_dc_ji(nx, ny, nz, frhatv, visc_rem_v, av_rem_dc_ji)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'CPU reference:'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_av_rem_cpu(nx, ny, nz, frhatv, visc_rem_v, av_rem_cpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    !$omp target exit data map(release: frhatv, visc_rem_v, &
    !$omp&   av_rem_omp, av_rem_omp_ji, av_rem_dc, av_rem_dc_ji)

    deallocate(frhatv, visc_rem_v)
    deallocate(av_rem_omp, av_rem_omp_ji, av_rem_dc, av_rem_dc_ji, av_rem_cpu)

  enddo ! isize

  write(*,*)
  write(*,*) 'Done.'

end program test_av_rem
