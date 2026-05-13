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
! One GPU variant + one CPU reference:
!   run_colnorm_omp    – separate !$omp target loop per phase; serial k
!                        outer with parallel ij inner in phase 2; fused
!                        kji target loop in phase 4
!   run_colnorm_cpu    – plain do loops, CPU reference
!
! Comparisons: GPU variant vs CPU.
! Timing: multiple problem sizes (32..256), 5 timed runs each.

#include "omp_macros.inc"

module col_norm_mod
  use test_utils_mod, only: dp
  implicit none

contains

  ! ------------------------------------------------------------------ !
  !  OMP separate loops: each phase is its own !$omp target loop       !
  !  kernel with its own map clauses.  Phase 2 launches nz-1 kernels  !
  !  (serial k, parallel ij).  col_sum persists on device via         !
  !  target enter/exit data.                                           !
  ! ------------------------------------------------------------------ !
  subroutine run_colnorm_omp(ni, nj, nz, wt, mask, wt_out)
    integer,  intent(in)  :: ni, nj, nz
    real(dp), intent(in)  :: wt(ni, nj, nz)
    real(dp), intent(in)  :: mask(ni, nj)
    real(dp), intent(out) :: wt_out(ni, nj, nz)
    real(dp) :: col_sum(ni, nj)
    integer  :: I, j, k

    !$omp target enter data map(alloc: col_sum)

    ! Phase 1: init column sums from k=1
    !$omp target teams COMBINED_LOOP collapse(2) map(to: wt)
    do j = 1, nj ; do I = 1, ni
      col_sum(I,j) = wt(I,j,1)
    enddo ; enddo

    ! Phase 2: accumulate layers k=2..nz (serial k, parallel ij per layer)
    do k = 2, nz
      !$omp target teams COMBINED_LOOP collapse(2) map(to: wt)
      do j = 1, nj ; do I = 1, ni
        col_sum(I,j) = col_sum(I,j) + wt(I,j,k)
      enddo ; enddo
    enddo

    ! Phase 3: conditional invert
    !$omp target teams COMBINED_LOOP collapse(2) map(to: mask)
    do j = 1, nj ; do I = 1, ni
      if (abs(col_sum(I,j)) > 0.0_dp) col_sum(I,j) = mask(I,j) / col_sum(I,j)
    enddo ; enddo

    ! Phase 4: broadcast 2-D result back to 3-D
    !$omp target teams COMBINED_LOOP collapse(3) map(to: wt) map(from: wt_out)
    do k = 1, nz ; do j = 1, nj ; do I = 1, ni
      wt_out(I,j,k) = wt(I,j,k) * col_sum(I,j)
    enddo ; enddo ; enddo

    !$omp target exit data map(release: col_sum)
  end subroutine run_colnorm_omp

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

end module col_norm_mod


program test_col_norm
  use col_norm_mod
  use test_utils_mod
  use omp_lib
  implicit none

  integer,  parameter :: nz           = 100
  integer,  parameter :: n_sizes      = 5
  integer,  parameter :: all_sizes(n_sizes) = [32, 64, 128, 256, 512]

  real(dp), allocatable :: wt(:,:,:)
  real(dp), allocatable :: mask(:,:)
  real(dp), allocatable :: wt_omp(:,:,:)
  real(dp), allocatable :: wt_cpu(:,:,:)

  integer  :: ni, nj, I, j, k, isize, irun
  real(dp) :: t0, t1
  real(dp) :: times(n_runs)
  logical  :: p1, all_pass

  do isize = 1, n_sizes
    ni = all_sizes(isize) ; nj = all_sizes(isize)

    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I0,A,I0,A,I0)') &
      'Size: ni=', ni, '  nj=', nj, '  nz=', nz
    write(*,'(A)') '========================================================'

    allocate(wt(ni, nj, nz), mask(ni, nj))
    allocate(wt_omp(ni, nj, nz))
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

    call run_colnorm_omp(ni, nj, nz, wt, mask, wt_omp)
    !$omp taskwait
    call compare_3d('omp    |cpu', wt_omp, wt_cpu, p1)

    all_pass = p1
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

    !$omp target enter data map(to: wt, mask) map(alloc: wt_omp)

    write(*,*) 'OMP (separate target loop per phase):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_omp(ni, nj, nz, wt, mask, wt_omp)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'CPU reference:'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_colnorm_cpu(ni, nj, nz, wt, mask, wt_cpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    !$omp target exit data map(release: wt, mask, wt_omp)

    deallocate(wt, mask, wt_omp, wt_cpu)

  enddo ! isize

  write(*,*)
  write(*,*) 'Done.'

end program test_col_norm
