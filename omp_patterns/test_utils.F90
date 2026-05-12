! Shared utilities for omp_patterns test programs.
!
! Provides: dp, pi, n_runs, print_timing_stats, compare_2d, compare_3d.
! Each test module does: use test_utils_mod, only: dp
! Each test program does: use test_utils_mod

module test_utils_mod
  implicit none
  integer,  parameter :: dp     = kind(1.0d0)
  real(dp), parameter :: pi     = 3.14159265358979323846_dp
  integer,  parameter :: n_runs = 5

contains

  subroutine print_timing_stats(times)
    real(dp), intent(in) :: times(:)
    real(dp) :: sorted(size(times)), tmp, tmin, tmax, tavg, tmed
    integer  :: n, i, j
    n = size(times)
    sorted = times
    do i = 2, n
      tmp = sorted(i) ; j = i - 1
      do while (j >= 1 .and. sorted(j) > tmp)
        sorted(j+1) = sorted(j) ; j = j - 1
      enddo
      sorted(j+1) = tmp
    enddo
    tmin = sorted(1)    ; tmax = sorted(n)
    tavg = sum(times) / real(n, dp)
    tmed = sorted((n+1)/2)
    write(*,'(A,F10.6,A,F10.6,A,F10.6,A,F10.6,A)') &
      '    min=', tmin, '  max=', tmax, '  avg=', tavg, '  med=', tmed, ' s'
  end subroutine print_timing_stats

  subroutine compare_2d(label, a, b, passed)
    character(*), intent(in)  :: label
    real(dp),     intent(in)  :: a(:,:), b(:,:)
    logical,      intent(out) :: passed
    integer :: i, j, ni, nj, ndiff
    ni = size(a,1) ; nj = size(a,2)
    ndiff = 0
    do j = 1, nj ; do i = 1, ni
      if (a(i,j) /= b(i,j)) ndiff = ndiff + 1
    enddo ; enddo
    passed = (ndiff == 0)
    if (passed) then
      write(*,'(A,A)') '  PASS  ', label
    else
      write(*,'(A,A,I0,A,I0,A)') '  FAIL  ', label//' — ', ndiff, ' of ', ni*nj, ' points differ'
      outer: do j = 1, nj ; do i = 1, ni
        if (a(i,j) /= b(i,j)) then
          write(*,'(A,2I4)')       '    first diff (i,j)=', i, j
          write(*,'(A,E28.20)')    '    a = ', a(i,j)
          write(*,'(A,E28.20)')    '    b = ', b(i,j)
          write(*,'(A,Z16,A,Z16)') '    bits  a=', transfer(a(i,j), 0_8), &
                                    '  b=', transfer(b(i,j), 0_8)
          exit outer
        endif
      enddo ; enddo outer
    endif
  end subroutine compare_2d

  subroutine compare_3d(label, a, b, passed)
    character(*), intent(in)  :: label
    real(dp),     intent(in)  :: a(:,:,:), b(:,:,:)
    logical,      intent(out) :: passed
    integer :: i, j, k, ni, nj, nz, ndiff
    ni = size(a,1) ; nj = size(a,2) ; nz = size(a,3)
    ndiff = 0
    do k = 1, nz ; do j = 1, nj ; do i = 1, ni
      if (a(i,j,k) /= b(i,j,k)) ndiff = ndiff + 1
    enddo ; enddo ; enddo
    passed = (ndiff == 0)
    if (passed) then
      write(*,'(A,A)') '  PASS  ', label
    else
      write(*,'(A,A,I0,A,I0,A)') '  FAIL  ', label//' — ', ndiff, ' of ', ni*nj*nz, ' points differ'
      outer: do k = 1, nz ; do j = 1, nj ; do i = 1, ni
        if (a(i,j,k) /= b(i,j,k)) then
          write(*,'(A,3I4)')       '    first diff (i,j,k)=', i, j, k
          write(*,'(A,E28.20)')    '    a = ', a(i,j,k)
          write(*,'(A,E28.20)')    '    b = ', b(i,j,k)
          write(*,'(A,Z16,A,Z16)') '    bits  a=', transfer(a(i,j,k), 0_8), &
                                    '  b=', transfer(b(i,j,k), 0_8)
          exit outer
        endif
      enddo ; enddo ; enddo outer
    endif
  end subroutine compare_3d

end module test_utils_mod
