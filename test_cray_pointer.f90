program test_cray_ptr

    double precision :: a(10, 20, 30), a_ptrdata(10, 20, 30), b(10, 20, 30)
    pointer(a_ptr, a_ptrdata)
    integer :: i, j, k, test

    ! ininitialize a
    n = 1
    do k=1,30
        do j=1,20
            do i=1,10
                a(i,j,k) = dble(n)
                n = n + 1
            enddo
        enddo
    enddo

    ! point a_ptrdata to a
    a_ptr = loc(a)

    ! correct value
    print *, sum(a_ptrdata)

    !$omp target teams distribute parallel do map(to: a_ptrdata) map(from: b) collapse(3)
    do k = 1,30
        do j = 1,20
            do i =1,10
                b(i,j,k) = a_ptrdata(i,j,k)
            enddo
        enddo
    enddo

    ! both forms below fail with:
    ! NVFORTRAN-S-0155-Compiler failed to translate accelerator region (see -Minfo messages): Unable to find associated device pointer
    ! do concurrent (k=1:30, j=1:20, i=1:10)
    !     b(i,j,k) = a_ptrdata(i,j,k)
    ! enddo
    ! !$omp target teams loop map(to: a_ptrdata) map(from: b) collapse(3)
    ! do k = 1,30
    !     do j = 1,20
    !         do i =1,10
    !             b(i,j,k) = a_ptrdata(i,j,k)
    !         enddo
    !     enddo
    ! enddo

    print *, sum(b)

end program test_cray_ptr
