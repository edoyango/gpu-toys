! Column-normalisation / 3-D broadcast pattern from MOM_barotropic.F90.
!
! The pattern (btstep, ~line 1108) is a 4-phase pipeline applied to a 3-D
! weight array (wt):
!
!   Phase 1 – init 2-D column sum from the first layer:
!               col_sum(I,j) = wt(I,j,1)
!   Phase 2 – accumulate remaining layers (serial k, concurrent ij):
!               do k=2,nz; col_sum(I,j) += wt(I,j,k)
!   Phase 3 – conditional invert (mask guard against zero or masked sum):
!               if (abs(col_sum) > 0) col_sum = mask / col_sum
!   Phase 4 – broadcast 2-D result back to 3-D:
!               wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j)
!
! Four GPU variants + one CPU reference:
!   run_colnorm_omp    – single !$omp target teams region; separate
!                        distribute parallel do for each phase; serial k
!                        loop in phase 2
!   run_colnorm_omp_ji – !$omp target teams loop num_teams(nj); one team
!                        per j-row; phases executed in sequence within each
!                        team; parallel do over i per phase, serial k
!   run_colnorm_dc     – separate do concurrent construct per phase;
!                        serial k outer with do concurrent(j,i) inner in
!                        phase 2; do concurrent filter mask in phase 3
!   run_colnorm_dc_ji  – phases 1-3 fused into do concurrent(j,i) with
!                        serial k accumulation and conditional invert inside
!                        each (j,i) iteration; phase 4 separate
!                        do concurrent(k,j,i)
!   run_colnorm_cpu    – plain do loops, CPU reference
!
! Comparisons: each GPU/DC variant vs CPU.
! Timing: multiple problem sizes (32..256), 5 timed runs each.

module col_norm_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)

contains

  ! ------------------------------------------------------------------ !
  !  OMP ij-outer: all four phases run inside a single target teams    !
  !  region using separate distribute parallel do constructs.          !
  !  col_sum persists on device between phases via target enter/exit.  !
  ! ------------------------------------------------------------------ !
  subroutine run_colnorm_omp(ni, nj, nz, nteams, wt, mask, wt_out)
    integer,  intent(in)  :: ni, nj, nz, nteams
    real(dp), intent(in)  :: wt(ni, nj, nz)
    real(dp), intent(in)  :: mask(ni, nj)
    real(dp), intent(out) :: wt_out(ni, nj, nz)
    real(dp) :: col_sum(ni, nj)
    integer  :: I, j, k

    !$omp target enter data map(alloc: col_sum)

    !$omp target teams num_teams(nteams) &
    !$omp&   map(to: wt, mask) map(from: wt_out)

    ! Phase 1: init column sums from k=1
    !$omp distribute parallel do collapse(2) private(I, j)
    do j = 1, nj ; do I = 1, ni
      col_sum(I,j) = wt(I,j,1)
    enddo ; enddo

    ! Phase 2: accumulate layers k=2..nz (serial k, distributed ij per layer)
    do k = 2, nz
      !$omp distribute parallel do collapse(2) private(I, j)
      do j = 1, nj ; do I = 1, ni
        col_sum(I,j) = col_sum(I,j) + wt(I,j,k)
      enddo ; enddo
    enddo

    ! Phase 3: conditional invert — mask guards against zero/masked columns
    !$omp distribute parallel do collapse(2) private(I, j)
    do j = 1, nj ; do I = 1, ni
      if (abs(col_sum(I,j)) > 0.0_dp) col_sum(I,j) = mask(I,j) / col_sum(I,j)
    enddo ; enddo

    ! Phase 4: broadcast 2-D result back to the full 3-D weight array
    !$omp distribute parallel do collapse(3) private(I, j, k)
    do k = 1, nz ; do j = 1, nj ; do I = 1, ni
      wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j)
    enddo ; enddo ; enddo

    !$omp end target teams

    !$omp target exit data map(release: col_sum)
  end subroutine run_colnorm_omp

  ! ------------------------------------------------------------------ !
  !  OMP ji-outer: j distributed (one team per j-row), i parallel     !
  !  within each team, k serial.  All four phases run within the same  !
  !  team in order, separated by implicit parallel-do barriers.        !
  ! ------------------------------------------------------------------ !
  subroutine run_colnorm_omp_ji(ni, nj, nz, wt, mask, wt_out)
    integer,  intent(in)  :: ni, nj, nz
    real(dp), intent(in)  :: wt(ni, nj, nz)
    real(dp), intent(in)  :: mask(ni, nj)
    real(dp), intent(out) :: wt_out(ni, nj, nz)
    real(dp) :: col_sum(ni, nj)
    integer  :: I, j, k

    !$omp target enter data map(alloc: col_sum)

    !$omp target teams loop num_teams(nj) &
    !$omp&   map(to: wt, mask) map(from: wt_out)
    do j = 1, nj

      ! Phase 1: init from k=1
      !$omp loop bind(parallel) private(I)
      do I = 1, ni ; col_sum(I,j) = wt(I,j,1) ; enddo

      ! Phase 2: accumulate k=2..nz
      do k = 2, nz
        !$omp loop bind(parallel) private(I)
        do I = 1, ni ; col_sum(I,j) = col_sum(I,j) + wt(I,j,k) ; enddo
      enddo

      ! Phase 3: conditional invert
      !$omp loop bind(parallel) private(I)
      do I = 1, ni
        if (abs(col_sum(I,j)) > 0.0_dp) col_sum(I,j) = mask(I,j) / col_sum(I,j)
      enddo

      ! Phase 4: apply to 3-D
      do k = 1, nz
        !$omp loop bind(parallel) private(I)
        do I = 1, ni ; wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j) ; enddo
      enddo

    enddo

    !$omp target exit data map(release: col_sum)
  end subroutine run_colnorm_omp_ji

  ! ------------------------------------------------------------------ !
  !  do concurrent, separate phases: mirrors the OMP ij-outer          !
  !  structure.  Phase 2 uses a serial k outer with do concurrent(j,i) !
  !  inner.  Phase 3 uses the do concurrent filter-mask syntax.        !
  ! ------------------------------------------------------------------ !
  subroutine run_colnorm_dc(ni, nj, nz, wt, mask, wt_out)
    integer,  intent(in)  :: ni, nj, nz
    real(dp), intent(in)  :: wt(ni, nj, nz)
    real(dp), intent(in)  :: mask(ni, nj)
    real(dp), intent(out) :: wt_out(ni, nj, nz)
    real(dp) :: col_sum(ni, nj)
    integer  :: I, j, k

    ! Phase 1
    do concurrent (j = 1:nj, I = 1:ni)
      col_sum(I,j) = wt(I,j,1)
    enddo

    ! Phase 2: serial k, concurrent ij per layer
    do k = 2, nz
      do concurrent (j = 1:nj, I = 1:ni)
        col_sum(I,j) = col_sum(I,j) + wt(I,j,k)
      enddo
    enddo

    ! Phase 3: conditional invert via do concurrent filter mask
    do concurrent (j = 1:nj, I = 1:ni, abs(col_sum(I,j)) > 0.0_dp)
      col_sum(I,j) = mask(I,j) / col_sum(I,j)
    enddo

    ! Phase 4
    do concurrent (k = 1:nz, j = 1:nj, I = 1:ni)
      wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j)
    enddo
  end subroutine run_colnorm_dc

  ! ------------------------------------------------------------------ !
  !  do concurrent ji: phases 1-3 fused into a single concurrent(j,i) !
  !  construct.  Each (j,i) thread runs serial k accumulation and      !
  !  conditional invert independently.  Phase 4 is a separate fused   !
  !  concurrent(k,j,i).                                               !
  ! ------------------------------------------------------------------ !
  subroutine run_colnorm_dc_ji(ni, nj, nz, wt, mask, wt_out)
    integer,  intent(in)  :: ni, nj, nz
    real(dp), intent(in)  :: wt(ni, nj, nz)
    real(dp), intent(in)  :: mask(ni, nj)
    real(dp), intent(out) :: wt_out(ni, nj, nz)
    real(dp) :: col_sum(ni, nj)
    integer  :: I, j, k

    ! Phases 1-3 fused: each (j,i) independently accumulates over k then inverts
    do concurrent (j = 1:nj, I = 1:ni)
      col_sum(I,j) = wt(I,j,1)
      do k = 2, nz
        col_sum(I,j) = col_sum(I,j) + wt(I,j,k)
      enddo
      if (abs(col_sum(I,j)) > 0.0_dp) col_sum(I,j) = mask(I,j) / col_sum(I,j)
    enddo

    ! Phase 4
    do concurrent (k = 1:nz, j = 1:nj, I = 1:ni)
      wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j)
    enddo
  end subroutine run_colnorm_dc_ji

  ! ------------------------------------------------------------------ !
  !  CPU reference: plain do loops, identical logic, no directives.   !
  ! ------------------------------------------------------------------ !
  subroutine run_colnorm_cpu(ni, nj, nz, wt, mask, wt_out)
    integer,  intent(in)  :: ni, nj, nz
    real(dp), intent(in)  :: wt(ni, nj, nz)
    real(dp), intent(in)  :: mask(ni, nj)
    real(dp), intent(out) :: wt_out(ni, nj, nz)
    real(dp) :: col_sum(ni, nj)
    integer  :: I, j, k

    do j = 1, nj ; do I = 1, ni
      col_sum(I,j) = wt(I,j,1)
    enddo ; enddo
    do k = 2, nz ; do j = 1, nj ; do I = 1, ni
      col_sum(I,j) = col_sum(I,j) + wt(I,j,k)
    enddo ; enddo ; enddo
    do j = 1, nj ; do I = 1, ni
      if (abs(col_sum(I,j)) > 0.0_dp) col_sum(I,j) = mask(I,j) / col_sum(I,j)
    enddo ; enddo
    do k = 1, nz ; do j = 1, nj ; do I = 1, ni
      wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j)
    enddo ; enddo ; enddo
  end subroutine run_colnorm_cpu

  ! ------------------------------------------------------------------ !
  !  Comparison helper for 3-D arrays (assumed-shape).                !
  ! ------------------------------------------------------------------ !
  subroutine compare_3d(label, a, b, passed)
    character(*), intent(in)  :: label
    real(dp),     intent(in)  :: a(:,:,:), b(:,:,:)
    logical,      intent(out) :: passed
    integer :: I, j, k, ni, nj, nz, ndiff
    ni = size(a,1) ; nj = size(a,2) ; nz = size(a,3)
    ndiff = 0
    do k = 1, nz ; do j = 1, nj ; do I = 1, ni
      if (a(I,j,k) /= b(I,j,k)) ndiff = ndiff + 1
    enddo ; enddo ; enddo
    passed = (ndiff == 0)
    if (passed) then
      write(*,'(A,A)') '  PASS  ', label
    else
      write(*,'(A,A,I0,A,I0,A)') '  FAIL  ', label//' — ', ndiff, ' of ', ni*nj*nz, ' points differ'
      outer: do k = 1, nz ; do j = 1, nj ; do I = 1, ni
        if (a(I,j,k) /= b(I,j,k)) then
          write(*,'(A,3I4)')       '    first diff (I,j,k)=', I, j, k
          write(*,'(A,E28.20)')    '    a = ', a(I,j,k)
          write(*,'(A,E28.20)')    '    b = ', b(I,j,k)
          write(*,'(A,Z16,A,Z16)') '    bits  a=', transfer(a(I,j,k), 0_8), &
                                    '  b=', transfer(b(I,j,k), 0_8)
          exit outer
        endif
      enddo ; enddo ; enddo outer
    endif
  end subroutine compare_3d

end module col_norm_mod


program test_col_norm
  use col_norm_mod
  use omp_lib
  implicit none

  integer,  parameter :: nz           = 20
  integer,  parameter :: n_sizes      = 4
  integer,  parameter :: n_runs       = 5
  integer,  parameter :: all_sizes(n_sizes) = [32, 64, 128, 256]
  real(dp), parameter :: pi           = 3.14159265358979323846_dp

  real(dp), allocatable :: wt(:,:,:)
  real(dp), allocatable :: mask(:,:)
  real(dp), allocatable :: wt_omp(:,:,:), wt_omp_ji(:,:,:)
  real(dp), allocatable :: wt_dc(:,:,:),  wt_dc_ji(:,:,:)
  real(dp), allocatable :: wt_cpu(:,:,:)

  integer  :: ni, nj, nteams, I, j, k, isize, irun
  real(dp) :: t0, t1
  logical  :: p1, p2, p3, p4, all_pass

  do isize = 1, n_sizes
    ni = all_sizes(isize) ; nj = all_sizes(isize)
    nteams = (ni * nj + 255) / 256

    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      'Size: ni=', ni, '  nj=', nj, '  nz=', nz, '  nteams=', nteams
    write(*,'(A)') '========================================================'

    allocate(wt(ni, nj, nz), mask(ni, nj))
    allocate(wt_omp(ni, nj, nz), wt_omp_ji(ni, nj, nz))
    allocate(wt_dc(ni, nj, nz),  wt_dc_ji(ni, nj, nz))
    allocate(wt_cpu(ni, nj, nz))

    ! wt is strictly positive (range [0.5, 1.5]) so column sums are always nonzero.
    ! mask is 0 for ~1-in-7 (I,j) pairs to exercise the conditional branch.
    do k = 1, nz ; do j = 1, nj ; do I = 1, ni
      wt(I,j,k) = 1.0_dp + 0.5_dp * sin(pi * real(I + j + k, dp) / 20.0_dp)
    enddo ; enddo ; enddo
    do j = 1, nj ; do I = 1, ni
      mask(I,j) = merge(0.0_dp, 1.0_dp, mod(I + j, 7) == 0)
    enddo ; enddo

    ! -------------------------------------------------------------- !
    !  Correctness (one run per variant)                              !
    ! -------------------------------------------------------------- !
    write(*,*)
    write(*,*) '--- Correctness ---'

    call run_colnorm_cpu(ni, nj, nz, wt, mask, wt_cpu)

    call run_colnorm_omp(ni, nj, nz, nteams, wt, mask, wt_omp)
    !$omp taskwait
    call compare_3d('omp    |cpu', wt_omp,    wt_cpu, p1)

    call run_colnorm_omp_ji(ni, nj, nz, wt, mask, wt_omp_ji)
    !$omp taskwait
    call compare_3d('omp_ji |cpu', wt_omp_ji, wt_cpu, p2)

    call run_colnorm_dc(ni, nj, nz, wt, mask, wt_dc)
    !$omp taskwait
    call compare_3d('dc     |cpu', wt_dc,     wt_cpu, p3)

    call run_colnorm_dc_ji(ni, nj, nz, wt, mask, wt_dc_ji)
    !$omp taskwait
    call compare_3d('dc_ji  |cpu', wt_dc_ji,  wt_cpu, p4)

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

    write(*,*) 'OMP (ij teams, separate phases):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_omp(ni, nj, nz, nteams, wt, mask, wt_omp)
      !$omp taskwait
      t1 = omp_get_wtime()
      write(*,'(A,I0,A,F10.6,A)') '    run ', irun, ': ', t1-t0, ' s'
    enddo

    write(*,*) 'OMP ji (j teams, parallel i, serial k per team):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_omp_ji(ni, nj, nz, wt, mask, wt_omp_ji)
      !$omp taskwait
      t1 = omp_get_wtime()
      write(*,'(A,I0,A,F10.6,A)') '    run ', irun, ': ', t1-t0, ' s'
    enddo

    write(*,*) 'do concurrent (separate phases, serial k outer in phase 2):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_dc(ni, nj, nz, wt, mask, wt_dc)
      !$omp taskwait
      t1 = omp_get_wtime()
      write(*,'(A,I0,A,F10.6,A)') '    run ', irun, ': ', t1-t0, ' s'
    enddo

    write(*,*) 'do concurrent ji (phases 1-3 fused, serial k within each (j,i)):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_dc_ji(ni, nj, nz, wt, mask, wt_dc_ji)
      !$omp taskwait
      t1 = omp_get_wtime()
      write(*,'(A,I0,A,F10.6,A)') '    run ', irun, ': ', t1-t0, ' s'
    enddo

    write(*,*) 'CPU reference:'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_cpu(ni, nj, nz, wt, mask, wt_cpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      write(*,'(A,I0,A,F10.6,A)') '    run ', irun, ': ', t1-t0, ' s'
    enddo

    deallocate(wt, mask, wt_omp, wt_omp_ji, wt_dc, wt_dc_ji, wt_cpu)

  enddo ! isize

  write(*,*)
  write(*,*) 'Done.'

end program test_col_norm
