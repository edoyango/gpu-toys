program test

  use nvtx

  implicit none
  double precision, allocatable:: a(:, :, :), b(:, :)
  integer:: n(3), i, j, k, rep
  character(len=4) arg
  integer, parameter:: nreps = 10000
  character(*), parameter:: fmt = "(I4)"

  do i = 1, 3
    call get_command_argument(number=i, value=arg)
    read(arg, *) n(i)
  enddo

  allocate(a(n(1), n(2), n(3)), b(n(1), n(2)), source=0.d0)

  call random_number(a)

  !$omp target enter data map(to: a, b)

  do rep = 1, nreps
    ! zero b so it isn't huge
    do concurrent (i=1:n(1), j=1:n(2))
      b(i, j) = 0.d0
    enddo
    call nvtxStartRange("mainloop")
#ifdef KOUTSIDE
    do k = 1,n(3)
      do concurrent (i=1:n(1), j=1:n(2))
        b(i, j) = b(i, j) + a(i, j, k)
      enddo
    enddo
#elif defined(JKI)
    do concurrent (j=1:n(2))
      do k = 1,n(3)
        do concurrent (i=1:n(1))
          b(i, j) = b(i, j) + a(i, j, k)
        enddo
      enddo
    enddo
#else ! k inside
    do concurrent(i=1:n(1), i=1:n(2))
      do k = 1, nz
        b(i, j) = b(i, j) + a(i, j, k)
      enddo
    enddo
#endif
    call nvtxEndRange
  enddo

  !$omp target exit data map(from: b) map(release: a)

  print *, sum(b)

end program test
