
!/// The i-shape of a dummy argument staggered at h- or v-points.
#  define SZI_(G)     G%isd:G%ied
!/// The j-shape of a dummy argument staggered at h- or u-points.
#  define SZJ_(G)     G%jsd:G%jed
!/// The k-shape of a layer dummy argument.
#  define SZK_(G)     G%ke
!/// The i-shape of a dummy argument staggered at q- or u-points.
#  define SZIB_(G)    G%IsdB:G%IedB
!/// The j-shape of a dummy argument staggered at q- or v-points.
#  define SZJB_(G)    G%JsdB:G%JedB

module prec_m
    use iso_fortran_env, only: real64, int64

    integer, parameter:: fp = real64
    integer, parameter:: chcksum_ip = int64
end module prec_m

module vertical_grid_m
    use prec_m, only: fp
    implicit none
    !> Describes the vertical ocean grid, including unit conversion factors
    type, public :: verticalGrid_type

    ! Commonly used parameters
    integer :: ke     !< The number of layers/levels in the vertical
    real(fp) :: H_subroundoff !< A thickness that is so small that it can be added to a thickness of
                            !! Angstrom or larger without changing it at the bit level [H ~> m or kg m-2].
                            !! If Angstrom is 0 or exceedingly small, this is negligible compared to 1e-17 m.

    end type verticalGrid_type

end module vertical_grid_m

module grid_m

    use prec_m, only: fp

    implicit none
    !> Ocean grid type. See mom_grid for details.
    type, public :: ocean_grid_type

        integer :: isc !< The start i-index of cell centers within the computational domain
        integer :: iec !< The end i-index of cell centers within the computational domain
        integer :: jsc !< The start j-index of cell centers within the computational domain
        integer :: jec !< The end j-index of cell centers within the computational domain

        integer :: isd !< The start i-index of cell centers within the data domain
        integer :: ied !< The end i-index of cell centers within the data domain
        integer :: jsd !< The start j-index of cell centers within the data domain
        integer :: jed !< The end j-index of cell centers within the data domain

        integer :: IscB !< The start i-index of cell vertices within the computational domain
        integer :: IecB !< The end i-index of cell vertices within the computational domain
        integer :: JscB !< The start j-index of cell vertices within the computational domain
        integer :: JecB !< The end j-index of cell vertices within the computational domain

        integer :: IsdB !< The start i-index of cell vertices within the data domain
        integer :: IedB !< The end i-index of cell vertices within the data domain
        integer :: JsdB !< The start j-index of cell vertices within the data domain
        integer :: JedB !< The end j-index of cell vertices within the data domain

        real(fp), allocatable, dimension(:,:) :: &
            mask2dT   !< 0 for land points and 1 for ocean points on the h-grid [nondim].

        real(fp), allocatable, dimension(:,:) :: &
            IdxCu, &     !< 1/dxCu [L-1 ~> m-1].
            IdyCu, &     !< 1/dyCu [L-1 ~> m-1].
            IareaCu   !< The masked inverse areas of u-grid cells [L-2 ~> m-2].

        real(fp), allocatable, dimension(:,:) :: &
            IdxCv, &     !< 1/dxCv [L-1 ~> m-1].
            IdyCv, &     !< 1/dyCv [L-1 ~> m-1].
            IareaCv      !< The masked inverse areas of v-grid cells [L-2 ~> m-2].
        real(fp), allocatable, dimension(:,:) :: &
            mask2dBu  !< 0 for boundary points and 1 for ocean points on the q grid [nondim].

    end type ocean_grid_type

end module grid_m

module hor_visc_m

    use prec_m, only: fp

    use grid_m, only: ocean_grid_type
    use vertical_grid_m, only: verticalGrid_type

    implicit none
    
    !> Control structure for horizontal viscosity
    type, public :: hor_visc_CS
        real(fp)    :: Kh_bg_min       !< The minimum value allowed for Laplacian horizontal
                                    !! viscosity [L2 T-1 ~> m2 s-1]. The default is 0.0.

        real(fp), allocatable, dimension(:,:) :: Kh_bg_xx
                            !< The background Laplacian viscosity at h points [L2 T-1 ~> m2 s-1].
                            !! The actual viscosity may be the larger of this
                            !! viscosity and the Smagorinsky and Leith viscosities.
        real(fp), allocatable, dimension(:,:) :: reduction_xx
                            !< The amount by which stresses through h points are reduced
                            !! due to partial barriers [nondim].
        real(fp), allocatable :: Kh_Max_xx(:,:)     !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].
        real(fp), allocatable, dimension(:,:) :: Kh_bg_xy
                            !< The background Laplacian viscosity at q points [L2 T-1 ~> m2 s-1].
                            !! The actual viscosity may be the larger of this
                            !! viscosity and the Smagorinsky and Leith viscosities.
        real(fp), allocatable, dimension(:,:) :: reduction_xy
                            !< The amount by which stresses through q points are reduced
                            !! due to partial barriers [nondim].
        real(fp), allocatable :: Kh_Max_xy(:,:)  !< The maximum permitted Laplacian viscosity [L2 T-1 ~> m2 s-1].

        real(fp), allocatable, dimension(:,:) :: &
            dx2h,           & !< Pre-calculated dx^2 at h points [L2 ~> m2]
            dy2h,           & !< Pre-calculated dy^2 at h points [L2 ~> m2]
            dx_dyT,         & !< Pre-calculated dx/dy at h points [nondim]
            dy_dxT            !< Pre-calculated dy/dx at h points [nondim]
        real(fp), allocatable, dimension(:,:) :: &
            dx2q,    & !< Pre-calculated dx^2 at q points [L2 ~> m2]
            dy2q,    & !< Pre-calculated dy^2 at q points [L2 ~> m2]
            dx_dyBu, & !< Pre-calculated dx/dy at q points [nondim]
            dy_dxBu    !< Pre-calculated dy/dx at q points [nondim]

        ! The following variables are precalculated time-invariant combinations of
        ! parameters and metric terms.
        real(fp), allocatable :: Laplac2_const_xx(:,:) !< Laplacian metric-dependent constants [L2 ~> m2]

        real(fp), allocatable :: Laplac2_const_xy(:,:) !< Laplacian metric-dependent constants [L2 ~> m2]

    end type hor_visc_CS

contains

    subroutine horizontal_viscosity(u, v, h, diffu, diffv, G, GV, CS)
        type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
        type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
        real(fp), dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                        intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
        real(fp), dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                        intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
        real(fp), dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                        intent(inout) :: h    !< Layer thicknesses [H ~> m or kg m-2].
        real(fp), dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                        intent(out) :: diffu  !< Zonal acceleration due to convergence of
                                                            !! along-coordinate stress tensor [L T-2 ~> m s-2]
        real(fp), dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                        intent(out) :: diffv  !< Meridional acceleration due to convergence
                                                            !! of along-coordinate stress tensor [L T-2 ~> m s-2].
        type(hor_visc_CS),             intent(inout) :: CS   !< Horizontal viscosity control structure
        real(fp), dimension(SZIB_(G),SZJ_(G)) :: &
            h_u        ! Thickness interpolated to u points [H ~> m or kg m-2].
        real(fp), dimension(SZI_(G),SZJB_(G)) :: &
            h_v        ! Thickness interpolated to v points [H ~> m or kg m-2].
        real(fp), dimension(SZI_(G),SZJ_(G)) :: &
            sh_xx, &      ! horizontal tension (du/dx - dv/dy) including metric terms [T-1 ~> s-1]
            str_xx,&      ! str_xx is the diagonal term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
                        ! at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
            dudx, dvdy    ! components in the horizontal tension [T-1 ~> s-1]
        real(fp), dimension(SZIB_(G),SZJB_(G)) :: &
            dvdx, dudy, & ! components in the shearing strain [T-1 ~> s-1]
            sh_xy, &     ! horizontal shearing strain (du/dy + dv/dx) including metric terms [T-1 ~> s-1]
            str_xy, &     ! str_xy is the cross term in the stress tensor [H L2 T-2 ~> m3 s-2 or kg s-2], but
                        ! at some points in the code it is not yet layer integrated, so is in [L2 T-2 ~> m2 s-2].
            hq          ! harmonic mean of the harmonic means of the u- & v point thicknesses [H ~> m or kg m-2]
                        ! This form guarantees that hq/hu < 4.
        real(fp), dimension(SZIB_(G),SZJB_(G)) :: &
            Kh, &           ! Laplacian  viscosity (h or q) [L2 T-1 ~> m2 s-1]
            Shear_mag, &    ! magnitude of the shear (h or q) [T-1 ~> s-1]
            hrat_min     ! h_min divided by the thickness at the stress point (h or q) [nondim]
        real(fp) :: h2uq, h2vq, h_min, h_neglect, sh_xx_sq, sh_xy_sq, h_neglect3
        integer:: i, j, k
        integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
        integer :: js_vort, je_vort, is_vort, ie_vort
        integer :: js_Kh, je_Kh, is_Kh, ie_Kh

        is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
        Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB


        js_Kh = Jsq ; je_Kh = je+1 ; is_Kh = Isq ; ie_Kh = ie+1
        js_vort = js-2 ; je_vort = Jeq+1 ; is_vort = is-2 ; ie_vort = Ieq+1

        h_neglect  = GV%H_subroundoff
        h_neglect3 = h_neglect*h_neglect*h_neglect

        !$omp target enter data map(alloc: &
        !$omp   h_u, &
        !$omp   h_v, &
        !$omp   sh_xx, &
        !$omp   str_xx, &
        !$omp   dudx, &
        !$omp   dvdy, &
        !$omp   dvdx, &
        !$omp   dudy, &
        !$omp   sh_xy, &
        !$omp   str_xy, &
        !$omp   hq, &
        !$omp   Kh, &
        !$omp   shear_mag, &
        !$omp   hrat_min &
        !$Omp )

        do k=1,nz
            ! The following are the forms of the horizontal tension and horizontal
            ! shearing strain advocated by Smagorinsky (1993) and discussed in
            ! Griffies and Hallberg (2000).

            ! Calculate horizontal tension
            do concurrent (j=Jsq-1:Jeq+2, i=Isq-1:Ieq+2)
                dudx(i,j) = CS%DY_dxT(i,j)*((G%IdyCu(I,j) * u(I,j,k)) - &
                                            (G%IdyCu(I-1,j) * u(I-1,j,k)))
            enddo

            do concurrent (j=Jsq-1:Jeq+2, i=Isq-1:Ieq+2)
                dvdy(i,j) = CS%DX_dyT(i,j)*((G%IdxCv(i,J) * v(i,J,k)) - &
                                            (G%IdxCv(i,J-1) * v(i,J-1,k)))
            enddo

            do concurrent (j=Jsq-1:Jeq+2, i=Isq-1:Ieq+2)
                sh_xx(i,j) = dudx(i,j) - dvdy(i,j)
            enddo

            ! Components for the shearing strain
            do concurrent (J=js_vort:je_vort, I=is_vort:ie_vort)
                dvdx(I,J) = CS%DY_dxBu(I,J)*((v(i+1,J,k)*G%IdyCv(i+1,J)) - (v(i,J,k)*G%IdyCv(i,J)))
                dudy(I,J) = CS%DX_dyBu(I,J)*((u(I,j+1,k)*G%IdxCu(I,j+1)) - (u(I,j,k)*G%IdxCu(I,j)))
            enddo
            do concurrent (j=js-2:je+2, I=is-2:Ieq+1)
                h_u(I,j) = 0.5_fp * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k))
            enddo
            do concurrent (J=js-2:Jeq+1, i=is-2:ie+2)
                h_v(i,J) = 0.5_fp * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k))
            enddo
            do concurrent (J=js-2:Jeq+1, I=is-2:Ieq+1)
                sh_xy(I,J) = G%mask2dBu(I,J) * ( dvdx(I,J) + dudy(I,J) )
            enddo
            do concurrent (j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
                sh_xx_sq = sh_xx(i,j)**2
                sh_xy_sq = 0.25_fp * ( ((sh_xy(I-1,J-1)**2) + (sh_xy(I,J)**2)) &
                                + ((sh_xy(I-1,J)**2) + (sh_xy(I,J-1)**2)) )
                Shear_mag(i,j) = sqrt(sh_xx_sq + sh_xy_sq)
            enddo
            do concurrent (j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
                h_min = min(h_u(I,j), h_u(I-1,j), h_v(i,J), h_v(i,J-1))
                hrat_min(i,j) = min(1.0_fp, h_min / (h(i,j,k) + h_neglect))
            enddo

            ! Static (pre-computed) background viscosity
            do concurrent (j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
                Kh(i,j) = CS%Kh_bg_xx(i,j)
            enddo
                do concurrent (j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
                    Kh(i,j) = max(Kh(i,j), CS%Laplac2_const_xx(i,j) * Shear_mag(i,j))
                enddo

            ! Place a floor on the viscosity, if desired.
            do concurrent (j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
                Kh(i,j) = max(Kh(i,j), CS%Kh_bg_min)
            enddo
                do concurrent (j=js_Kh:je_Kh, i=is_Kh:ie_Kh)
                Kh(i,j) = min(Kh(i,j), hrat_min(i,j) * CS%Kh_Max_xx(i,j))
                enddo

            do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
                str_xx(i,j) = -Kh(i,j) * sh_xx(i,j)
            enddo
            do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                sh_xy_sq = sh_xy(I,J)**2
                sh_xx_sq = 0.25_fp * ( ((sh_xx(i,j)**2) + (sh_xx(i+1,j+1)**2)) &
                                + ((sh_xx(i,j+1)**2) + (sh_xx(i+1,j)**2)) )
                Shear_mag(I,J) = sqrt(sh_xy_sq + sh_xx_sq)
            enddo

            do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                h2uq = 4.0_fp * (h_u(I,j) * h_u(I,j+1))
                h2vq = 4.0_fp * (h_v(i,J) * h_v(i+1,J))
                hq(I,J) = (2.0_fp * (h2uq * h2vq)) &
                    / (h_neglect3 + (h2uq + h2vq) * ((h_u(I,j) + h_u(I,j+1)) + (h_v(i,J) + h_v(i+1,J))))
            enddo

            do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                h_min = min(h_u(I,j), h_u(I,j+1), h_v(i,J), h_v(i+1,J))
                hrat_min(I,J) = min(1.0_fp, h_min / (hq(I,J) + h_neglect))
            enddo

            ! Static (pre-computed) background viscosity
            do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                Kh(I,J) = CS%Kh_bg_xy(I,J)
            enddo
                do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                    Kh(I,J) = max(Kh(I,J), CS%Laplac2_const_xy(I,J) * Shear_mag(I,J) )
                enddo
            do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                ! Place a floor on the viscosity, if desired.
                Kh(I,J) = max(Kh(I,J), CS%Kh_bg_min)
            enddo
                do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                    Kh(I,J) = min(Kh(I,J), hrat_min(I,J) * CS%Kh_Max_xy(I,J))
                enddo
                do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                    str_xy(I,J) = -Kh(I,J) * sh_xy(I,J)
                enddo
            ! This changes the units of str_xx from [L2 T-2 ~> m2 s-2] to [H L2 T-2 ~> m3 s-2 or kg s-2].
            do concurrent (j=Jsq:Jeq+1, i=Isq:Ieq+1)
                str_xx(i,j) = str_xx(i,j) * (h(i,j,k) * CS%reduction_xx(i,j))
            enddo

                do concurrent (J=js-1:Jeq, I=is-1:Ieq)
                    str_xy(I,J) = str_xy(I,J) * (hq(I,J) * G%mask2dBu(I,J) * CS%reduction_xy(I,J))
                enddo

            ! Evaluate 1/h x.Div(h Grad u) or the biharmonic equivalent.
            do concurrent (j=js:je, I=Isq:Ieq)
                diffu(I,j,k) = ((G%IdxCu(I,j)*((CS%dx2q(I,J-1)*str_xy(I,J-1)) - (CS%dx2q(I,J)*str_xy(I,J))) + &
                                G%IdyCu(I,j)*((CS%dy2h(i,j)*str_xx(i,j)) - (CS%dy2h(i+1,j)*str_xx(i+1,j)))) * &
                                G%IareaCu(I,j)) / (h_u(I,j) + h_neglect)
            enddo

            ! Evaluate 1/h y.Div(h Grad u) or the biharmonic equivalent.
            do concurrent (J=Jsq:Jeq, i=is:ie)
                diffv(i,J,k) = ((G%IdyCv(i,J)*((CS%dy2q(I-1,J)*str_xy(I-1,J)) - (CS%dy2q(I,J)*str_xy(I,J))) - &
                                G%IdxCv(i,J)*((CS%dx2h(i,j)*str_xx(i,j)) - (CS%dx2h(i,j+1)*str_xx(i,j+1)))) * &
                                G%IareaCv(i,J)) / (h_v(i,J) + h_neglect)
            enddo

        enddo ! end of k loop

        !$omp target exit data map(release: &
        !$omp   h_u, &
        !$omp   h_v, &
        !$omp   sh_xx, &
        !$omp   str_xx, &
        !$omp   dudx, &
        !$omp   dvdy, &
        !$omp   dvdx, &
        !$omp   dudy, &
        !$omp   sh_xy, &
        !$omp   str_xy, &
        !$omp   hq, &
        !$omp   Kh, &
        !$omp   shear_mag, &
        !$omp   hrat_min &
        !$Omp )

    end subroutine horizontal_viscosity

end module hor_visc_m

program test_horizontal_viscosity
    use iso_fortran_env, only: int64
    use prec_m
    use grid_m
    use vertical_grid_m
    use hor_visc_m
    implicit none

    type(ocean_grid_type) :: G
    type(verticalGrid_type) :: GV
    type(hor_visc_CS) :: CS
    real(fp), allocatable :: u(:,:,:), v(:,:,:), h(:,:,:)
    real(fp), allocatable :: diffu(:,:,:), diffv(:,:,:)
    integer :: i, j, k, ni, nj, nk, niter
    integer(int64) :: start_time, end_time, count_rate
    real(fp) :: elapsed_time
    integer(chcksum_ip) :: checksum_u, checksum_v
    real(fp), parameter :: pi = 4.0_fp * atan(1.0_fp)
    character(len=32) :: arg
    integer :: num_args

    ! Get grid dimensions from command line arguments
    num_args = command_argument_count()
    
    if (num_args /= 3) then
        print *, "Usage: test_horvisc <ni> <nj> <nk>"
        print *, "  ni: number of i points"
        print *, "  nj: number of j points"
        print *, "  nk: number of k levels"
        print *, ""
        print *, "Example: test_horvisc 100 100 50"
        stop
    endif
    
    call get_command_argument(1, arg)
    read(arg, *) ni
    call get_command_argument(2, arg)
    read(arg, *) nj
    call get_command_argument(3, arg)
    read(arg, *) nk

    ! Initialize grid structure
    call initialize_grid(G, ni, nj)
    
    ! Initialize vertical grid
    call initialize_vertical_grid(GV, nk)
    
    ! Initialize horizontal viscosity control structure
    call initialize_hor_visc_cs(CS, G, GV)
    
    ! Allocate velocity and thickness arrays
    allocate(u(G%IsdB:G%IedB, G%jsd:G%jed, GV%ke))
    allocate(v(G%isd:G%ied, G%JsdB:G%JedB, GV%ke))
    allocate(h(G%isd:G%ied, G%jsd:G%jed, GV%ke))
    allocate(diffu(G%IsdB:G%IedB, G%jsd:G%jed, GV%ke))
    allocate(diffv(G%isd:G%ied, G%JsdB:G%JedB, GV%ke))
    
    ! Initialize with non-uniform values
    call initialize_velocity_and_thickness(u, v, h, G, GV, pi)
    
    print *, "Grid dimensions: ni=", ni, "nj=", nj, "nk=", nk
    print *, "Running horizontal_viscosity subroutine..."

    !$omp target enter data map(alloc: diffu, diffv)
    
    ! Time the subroutine
    call system_clock(start_time, count_rate)

    do niter = 1, 10
    
        call horizontal_viscosity(u, v, h, diffu, diffv, G, GV, CS)

    enddo
    
    call system_clock(end_time)
    
    elapsed_time = real(end_time - start_time, fp) / real(count_rate, fp)
    
    ! Calculate checksums
    checksum_u = calculate_checksum(diffu)
    checksum_v = calculate_checksum(diffv)
    
    print *, ""
    print *, "=== RESULTS ==="
    print *, "Elapsed time: ", elapsed_time, " seconds"
    print *, "Checksum diffu: ", checksum_u
    print *, "Checksum diffv: ", checksum_v
    print *, ""
    print *, "Sample diffu values:"
    print *, "  diffu(1,1,1) = ", diffu(G%IscB, G%jsc, 1)
    print *, "  diffu(mid)   = ", diffu((G%IscB+G%IecB)/2, (G%jsc+G%jec)/2, nk/2)
    print *, "Sample diffv values:"
    print *, "  diffv(1,1,1) = ", diffv(G%isc, G%JscB, 1)
    print *, "  diffv(mid)   = ", diffv((G%isc+G%iec)/2, (G%JscB+G%JecB)/2, nk/2)

contains

    subroutine initialize_grid(G, ni, nj)
        type(ocean_grid_type), intent(inout) :: G
        integer, intent(in) :: ni, nj
        integer :: i, j
        real(fp) :: x, y

        ! Set up index bounds (with halo of 2)
        G%isd = 1; G%ied = ni
        G%jsd = 1; G%jed = nj
        G%isc = 3; G%iec = ni - 2
        G%jsc = 3; G%jec = nj - 2
        G%IsdB = 0; G%IedB = ni - 1
        G%JsdB = 0; G%JedB = nj - 1
        G%IscB = 2; G%IecB = ni - 3
        G%JscB = 2; G%JecB = nj - 3

        ! Allocate grid arrays
        allocate(G%mask2dT(G%isd:G%ied, G%jsd:G%jed))
        allocate(G%IdxCu(G%IsdB:G%IedB, G%jsd:G%jed))
        allocate(G%IdyCu(G%IsdB:G%IedB, G%jsd:G%jed))
        allocate(G%IareaCu(G%IsdB:G%IedB, G%jsd:G%jed))
        allocate(G%IdxCv(G%isd:G%ied, G%JsdB:G%JedB))
        allocate(G%IdyCv(G%isd:G%ied, G%JsdB:G%JedB))
        allocate(G%IareaCv(G%isd:G%ied, G%JsdB:G%JedB))
        allocate(G%mask2dBu(G%IsdB:G%IedB, G%JsdB:G%JedB))

        ! Initialize masks (all ocean points)
        G%mask2dT = 1.0_fp
        G%mask2dBu = 1.0_fp

        ! Initialize metric terms with spatial variation
        do j = G%jsd, G%jed
            do i = G%IsdB, G%IedB
                x = real(i, fp) / real(ni, fp)
                y = real(j, fp) / real(nj, fp)
                G%IdxCu(i,j) = 0.001_fp * (1.0_fp + 0.1_fp * sin(2.0_fp * pi * x))
                G%IdyCu(i,j) = 0.001_fp * (1.0_fp + 0.1_fp * cos(2.0_fp * pi * y))
                G%IareaCu(i,j) = G%IdxCu(i,j) * G%IdyCu(i,j)
            enddo
        enddo

        do j = G%JsdB, G%JedB
            do i = G%isd, G%ied
                x = real(i, fp) / real(ni, fp)
                y = real(j, fp) / real(nj, fp)
                G%IdxCv(i,j) = 0.001_fp * (1.0_fp + 0.1_fp * cos(2.0_fp * pi * x))
                G%IdyCv(i,j) = 0.001_fp * (1.0_fp + 0.1_fp * sin(2.0_fp * pi * y))
                G%IareaCv(i,j) = G%IdxCv(i,j) * G%IdyCv(i,j)
            enddo
        enddo

        !$omp target enter data map(to: &
        !$omp   G, &
        !$omp   G%mask2dT, &
        !$omp   G%IdxCu, &
        !$omp   G%IdyCu, &
        !$omp   G%IareaCu, &
        !$omp   G%IdxCv, &
        !$omp   G%IdyCv, &
        !$omp   G%IareaCv, &
        !$omp   G%mask2dBu &
        !$omp )
    end subroutine initialize_grid

    subroutine initialize_vertical_grid(GV, nk)
        type(verticalGrid_type), intent(inout) :: GV
        integer, intent(in) :: nk

        GV%ke = nk
        GV%H_subroundoff = 1.0e-30_fp
        !$omp target enter data map(to: GV)
    end subroutine initialize_vertical_grid

    subroutine initialize_hor_visc_cs(CS, G, GV)
        type(hor_visc_CS), intent(inout) :: CS
        type(ocean_grid_type), intent(in) :: G
        type(verticalGrid_type), intent(in) :: GV
        integer :: i, j
        real(fp) :: x, y, pi_local

        pi_local = 4.0_fp * atan(1.0_fp)

        CS%Kh_bg_min = 1.0_fp

        ! Allocate CS arrays
        allocate(CS%Kh_bg_xx(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%reduction_xx(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%Kh_Max_xx(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%Kh_bg_xy(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%reduction_xy(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%Kh_Max_xy(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%dx2h(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%dy2h(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%dx_dyT(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%dy_dxT(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%dx2q(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%dy2q(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%dx_dyBu(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%dy_dxBu(G%IsdB:G%IedB, G%JsdB:G%JedB))
        allocate(CS%Laplac2_const_xx(G%isd:G%ied, G%jsd:G%jed))
        allocate(CS%Laplac2_const_xy(G%IsdB:G%IedB, G%JsdB:G%JedB))

        ! Initialize with spatially varying values
        do j = G%jsd, G%jed
            do i = G%isd, G%ied
                x = real(i, fp) / real(G%ied - G%isd, fp)
                y = real(j, fp) / real(G%jed - G%jsd, fp)
                CS%Kh_bg_xx(i,j) = 100.0_fp * (1.0_fp + 0.5_fp * sin(pi_local * x) * cos(pi_local * y))
                CS%reduction_xx(i,j) = 1.0_fp
                CS%Kh_Max_xx(i,j) = 1000.0_fp
                CS%dx2h(i,j) = 1.0e6_fp
                CS%dy2h(i,j) = 1.0e6_fp
                CS%dx_dyT(i,j) = 1.0_fp + 0.05_fp * sin(2.0_fp * pi_local * x)
                CS%dy_dxT(i,j) = 1.0_fp + 0.05_fp * cos(2.0_fp * pi_local * y)
                CS%Laplac2_const_xx(i,j) = 50.0_fp
            enddo
        enddo

        do j = G%JsdB, G%JedB
            do i = G%IsdB, G%IedB
                x = real(i, fp) / real(G%IedB - G%IsdB, fp)
                y = real(j, fp) / real(G%JedB - G%JsdB, fp)
                CS%Kh_bg_xy(i,j) = 100.0_fp * (1.0_fp + 0.5_fp * cos(pi_local * x) * sin(pi_local * y))
                CS%reduction_xy(i,j) = 1.0_fp
                CS%Kh_Max_xy(i,j) = 1000.0_fp
                CS%dx2q(i,j) = 1.0e6_fp
                CS%dy2q(i,j) = 1.0e6_fp
                CS%dx_dyBu(i,j) = 1.0_fp + 0.05_fp * cos(2.0_fp * pi_local * x)
                CS%dy_dxBu(i,j) = 1.0_fp + 0.05_fp * sin(2.0_fp * pi_local * y)
                CS%Laplac2_const_xy(i,j) = 50.0_fp
            enddo
        enddo
        !$omp target enter data map(to: CS, &
        !$omp                           CS%Kh_bg_xx, &
        !$omp                           CS%reduction_xx, &
        !$omp                           CS%Kh_Max_xx, &
        !$omp                           CS%Kh_bg_xy, &
        !$omp                           CS%reduction_xy, &
        !$omp                           CS%Kh_Max_xy, &
        !$omp                           CS%dx2h, &
        !$omp                           CS%dy2h, &
        !$omp                           CS%dx_dyT, &
        !$omp                           CS%dy_dxT, &
        !$omp                           CS%dx2q, &
        !$omp                           CS%dy2q, &
        !$omp                           CS%dx_dyBu, &
        !$omp                           CS%dy_dxBu, &
        !$omp                           CS%Laplac2_const_xx, &
        !$omp                           CS%Laplac2_const_xy &
        !$omp )
    end subroutine initialize_hor_visc_cs

    subroutine initialize_velocity_and_thickness(u, v, h, G, GV, pi)
        type(ocean_grid_type), intent(in) :: G
        type(verticalGrid_type), intent(in) :: GV
        real(fp), intent(inout) :: u(SZIB_(G),SZJ_(G),SZK_(GV))
        real(fp), intent(inout) :: v(SZI_(G),SZJB_(G),SZK_(GV))
        real(fp), intent(inout) :: h(SZI_(G),SZJ_(G),SZK_(GV))
        real(fp), intent(in) :: pi
        integer :: i, j, k
        real(fp) :: x, y, z, ni, nj

        ni = real(G%ied - G%isd, fp)
        nj = real(G%jed - G%jsd, fp)

        ! Initialize u velocity with vortex-like pattern
        do k = 1, GV%ke
            z = real(k, fp) / real(GV%ke, fp)
            do j = G%jsd, G%jed
                do i = G%IsdB, G%IedB
                    x = real(i, fp) / ni
                    y = real(j, fp) / nj
                    u(i,j,k) = 0.1_fp * sin(2.0_fp * pi * x) * cos(2.0_fp * pi * y) * (1.0_fp - z)
                enddo
            enddo
        enddo

        ! Initialize v velocity with complementary pattern
        do k = 1, GV%ke
            z = real(k, fp) / real(GV%ke, fp)
            do j = G%JsdB, G%JedB
                do i = G%isd, G%ied
                    x = real(i, fp) / ni
                    y = real(j, fp) / nj
                    v(i,j,k) = -0.1_fp * cos(2.0_fp * pi * x) * sin(2.0_fp * pi * y) * (1.0_fp - z)
                enddo
            enddo
        enddo

        ! Initialize thickness with vertical and horizontal variation
        do k = 1, GV%ke
            z = real(k, fp) / real(GV%ke, fp)
            do j = G%jsd, G%jed
                do i = G%isd, G%ied
                    x = real(i, fp) / ni
                    y = real(j, fp) / nj
                    h(i,j,k) = 10.0_fp + 5.0_fp * (1.0_fp - z) * (1.0_fp + 0.3_fp * sin(pi * x) * cos(pi * y))
                enddo
            enddo
        enddo
        !$omp target enter data map(to: u, v, h)
    end subroutine initialize_velocity_and_thickness

    function calculate_checksum(array) result(checksum)
        real(fp), intent(in) :: array(:,:,:)
        integer(chcksum_ip) :: checksum
        integer :: i, j, k
        integer(chcksum_ip) :: temp

        !$omp target update from(array)

        checksum = 0_chcksum_ip
        do k = 1, size(array, 3)
            do j = 1, size(array, 2)
                do i = 1, size(array, 1)
                    temp = transfer(array(i,j,k), checksum)
                    checksum = checksum + temp
                enddo
            enddo
        enddo
    end function calculate_checksum

end program test_horizontal_viscosity