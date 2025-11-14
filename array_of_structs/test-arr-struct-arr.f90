program test

    use nvtx

    implicit none
    type mystruct
        double precision, allocatable:: a(:, :, :)
    end type mystruct
    type(mystruct), allocatable, target:: arr_of_structs(:)
    double precision, pointer:: ptr(:, :)
    integer:: n, nx, ny, nz, nargs, i, j, k, m
    character(12):: argc

    nargs = command_argument_count()

    if (nargs < 4) error stop "need 4 args"

    call get_command_argument(1, argc)
    read(argc, "(I12)") n
    call get_command_argument(2, argc)
    read(argc, "(I12)") nx
    call get_command_argument(3, argc)
    read(argc, "(I12)") ny
    call get_command_argument(4, argc)
    read(argc, "(I12)") nz

    allocate(arr_of_structs(n))

    do m = 1, n
        allocate(arr_of_structs(m)%a(nx,ny,nz))
    enddo

    call init

    call nvtxStartRange("kjmi")
    do k = 1, nz
        !$omp target teams loop
        do j = 1, ny
            do concurrent(m=1:n, i=1:nx)
                arr_of_structs(m)%a(i,j,k) = arr_of_structs(m)%a(i,j,k) + 1.d0
            enddo
        enddo
    enddo
    call nvtxEndRange

    print*,mysum(arr_of_structs)

    call init

    call nvtxStartRange("mkji flat")
    do concurrent (m=1:n, k=1:nz, j=1:ny, i=1:nx)
        arr_of_structs(m)%a(i,j,k) = arr_of_structs(m)%a(i,j,k) + 1
    enddo
    call nvtxEndRange

    print*,mysum(arr_of_structs)

    call init

    call nvtxStartRange("mkji with ptr")
    do m=1,n
        do k = 1, nz
            ptr => arr_of_structs(m)%a(:,:,k)
            !$omp target teams loop
            do j = 1, ny
                do concurrent(i=1:nx)
                    ptr(i,j) = ptr(i,j) + 1.d0
                enddo
            enddo
        enddo
    enddo
    call nvtxEndRange

    print*,mysum(arr_of_structs)

contains

    double precision function mysum(arr)
        type(mystruct):: arr(:)

        mysum = 0.d0
        do m=1,n
            do k=1,nz
                do j = 1,ny
                    do i = 1,nx
                        mysum = mysum + arr(m)%a(i,j,k)
                    enddo
                enddo
            enddo
        enddo
    end function mysum

    subroutine init()

        do m = 1,n
            do k = 1,nz
                do j = 1, ny
                    do i = 1, nx
                        arr_of_structs(m)%a(i,j,k) = dble((m-1)*nz*ny*nx + (k-1)*ny*nx + (j-1)*nx + i)/1000.d0
                    enddo
                enddo
            enddo
        enddo

    end subroutine init

end program test