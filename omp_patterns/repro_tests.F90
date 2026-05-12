! Combined GPU/CPU bitwise-reproducibility tests.
!
! Test 1 – continuity flux
!   Pattern: single !$omp target teams region, two sequential
!   !$omp distribute parallel do nests per layer (visc_rem copy then flux).
!
! Test 2 – zonal flux adjustment
!   Pattern: !$omp target enter/exit data for persistent work arrays,
!   single !$omp target teams enclosing a sequential Newton/bisection
!   do-itt loop with multiple !$omp distribute parallel do nests inside.

module repro_mod
  implicit none
  integer, parameter :: dp       = kind(1.0d0)
  integer, parameter :: max_itts = 20           ! flux-adjust iteration cap

contains

  ! ------------------------------------------------------------------ !
  !  Shared flux element: simplified second-order upwind scheme.        !
  !  h / h_p1  : cell-centre thickness left / right of the u-face [m]  !
  !  visc_rem  : viscous-remainder coefficient [nondim]                 !
  !  dy / IdxT / IdxT_p1 : face width and inverse grid spacings        !
  ! ------------------------------------------------------------------ !
  elemental subroutine flux_elem(u, h, h_p1, visc_rem, dy, IdxT, IdxT_p1, dt, uh, duhdu)
    !$omp declare target
    real(dp), intent(in)  :: u, h, h_p1, visc_rem, dy, IdxT, IdxT_p1, dt
    real(dp), intent(out) :: uh, duhdu
    real(dp) :: CFL, h_upwind

    if (u > 0.0_dp) then
      CFL      = u * dt * IdxT
      h_upwind = h + 0.5_dp * (1.0_dp - CFL) * (h_p1 - h)
      uh       = dy * u * h_upwind
    elseif (u < 0.0_dp) then
      CFL      = -u * dt * IdxT_p1
      h_upwind = h_p1 + 0.5_dp * (1.0_dp - CFL) * (h - h_p1)
      uh       = dy * u * h_upwind
    else
      h_upwind = 0.5_dp * (h + h_p1)
      uh       = 0.0_dp
    endif
    duhdu = dy * h_upwind * visc_rem
  end subroutine flux_elem

  ! ================================================================== !
  !  Test 1 – continuity flux                                           !
  ! ================================================================== !

  ! GPU: single target teams, two distribute parallel do per layer —
  ! (1) copy visc_rem_u → block-indexed visc_rem, (2) call flux_elem.
  subroutine run_continuity_gpu(nx, ny, nz, ni, nj,              &
                                 i_start, i_end, j_start, j_end,  &
                                 nteams, u, h_in, visc_rem_u,      &
                                 dy_Cu, IdxT, IdxT_xp1, dt,        &
                                 uh_t, duhdu_out)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end, nteams
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem_u(nx, ny, nz)
    real(dp), intent(in)  :: dy_Cu(nx, ny)
    real(dp), intent(in)  :: IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt
    real(dp), intent(out) :: uh_t(ni, nj, nz)
    real(dp), intent(out) :: duhdu_out(ni, nj, nz)
    real(dp) :: visc_rem(ni, nj, nz)
    integer  :: i, j, k, ii, jj

    !$omp target teams num_teams(nteams)                          &
    !$omp& map(to: u, h_in, visc_rem_u, dy_Cu, IdxT, IdxT_xp1)  &
    !$omp& map(from: uh_t, duhdu_out)                             &
    !$omp& map(alloc: visc_rem)
    do k = 1, nz
      !$omp distribute parallel do collapse(2) private(ii, jj)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        visc_rem(ii,jj,k) = visc_rem_u(i,j,k)
      enddo ; enddo
      !$omp distribute parallel do collapse(2) private(ii, jj)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        call flux_elem(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                       dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,                  &
                       uh_t(ii,jj,k), duhdu_out(ii,jj,k))
      enddo ; enddo
    enddo
    !$omp end target teams
  end subroutine run_continuity_gpu

  ! CPU reference: identical logic, no OpenMP.
  subroutine run_continuity_cpu(nx, ny, nz, ni, nj,              &
                                 i_start, i_end, j_start, j_end,  &
                                 u, h_in, visc_rem_u,              &
                                 dy_Cu, IdxT, IdxT_xp1, dt,        &
                                 uh_t, duhdu_out)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem_u(nx, ny, nz)
    real(dp), intent(in)  :: dy_Cu(nx, ny)
    real(dp), intent(in)  :: IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt
    real(dp), intent(out) :: uh_t(ni, nj, nz)
    real(dp), intent(out) :: duhdu_out(ni, nj, nz)
    real(dp) :: visc_rem(ni, nj, nz)
    integer  :: i, j, k, ii, jj

    do k = 1, nz
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        visc_rem(ii,jj,k) = visc_rem_u(i,j,k)
      enddo ; enddo
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        call flux_elem(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                       dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,                  &
                       uh_t(ii,jj,k), duhdu_out(ii,jj,k))
      enddo ; enddo
    enddo
  end subroutine run_continuity_cpu

  ! ================================================================== !
  !  Test 2 – zonal flux adjustment                                     !
  ! ================================================================== !

  ! GPU: persistent work arrays via target enter/exit data; sequential
  ! Newton/bisection do-itt loop runs entirely inside one target teams
  ! region, with distribute parallel do for all (j,i) work.
  subroutine zonal_flux_adjust_gpu(nx, ny, nz, ni, nj,                    &
                                    i_start, i_end, j_start, j_end, nteams, &
                                    u, h_in, visc_rem,                       &
                                    uhbt, uh_tot_0, duhdu_tot_0,             &
                                    du_max_CFL, du_min_CFL, do_I_in,         &
                                    IareaT, IareaT_xp1,                      &
                                    dy_Cu, IdxT, IdxT_xp1,                   &
                                    dt, tol_eta_base, tol_vel, better_iter,  &
                                    du, uh_3d)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end, nteams
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem(ni, nj, nz)       ! block-indexed input
    real(dp), intent(in)  :: uhbt(ni, nj)               ! target barotropic transport
    real(dp), intent(in)  :: uh_tot_0(ni, nj)           ! initial summed transport
    real(dp), intent(in)  :: duhdu_tot_0(ni, nj)        ! initial summed d(uh)/du
    real(dp), intent(in)  :: du_max_CFL(ni, nj)         ! CFL upper bound on du
    real(dp), intent(in)  :: du_min_CFL(ni, nj)         ! CFL lower bound on du
    logical,  intent(in)  :: do_I_in(ni, nj)            ! active-column mask
    real(dp), intent(in)  :: IareaT(nx, ny)             ! 1/area at (i,j)
    real(dp), intent(in)  :: IareaT_xp1(nx, ny)         ! 1/area at (i+1,j)
    real(dp), intent(in)  :: dy_Cu(nx, ny)
    real(dp), intent(in)  :: IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt, tol_eta_base, tol_vel
    logical,  intent(in)  :: better_iter
    real(dp), intent(out) :: du(ni, nj)                 ! velocity adjustment
    real(dp), intent(out) :: uh_3d(ni, nj, nz)          ! updated layer fluxes

    real(dp) :: uh_err(ni, nj), uh_err_best(ni, nj)
    real(dp) :: duhdu_tot(ni, nj)
    real(dp) :: du_min(ni, nj), du_max(ni, nj)
    logical  :: do_I(ni, nj)
    integer  :: i, j, k, ii, jj, itt
    real(dp) :: tol_eta, u_new, duhdu_loc, ddu, du_prev

    !$omp target enter data &
    !$omp&   map(alloc: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)

    !$omp target teams distribute parallel do collapse(2) private(ii, jj) &
    !$omp&   map(to: do_I_in, du_max_CFL, du_min_CFL, uh_tot_0, uhbt, duhdu_tot_0) &
    !$omp&   map(from: du)
    do j = j_start, j_end ; do i = i_start, i_end
      ii = i - i_start + 1 ; jj = j - j_start + 1
      du(ii,jj)          = 0.0_dp  ;  do_I(ii,jj)      = do_I_in(ii,jj)
      du_max(ii,jj)      = du_max_CFL(ii,jj)
      du_min(ii,jj)      = du_min_CFL(ii,jj)
      uh_err(ii,jj)      = uh_tot_0(ii,jj) - uhbt(ii,jj)
      duhdu_tot(ii,jj)   = duhdu_tot_0(ii,jj)
      uh_err_best(ii,jj) = abs(uh_err(ii,jj))
    enddo ; enddo

    !$omp target teams private(k,itt,tol_eta)                                                             &
    !$omp&   map(to: u, h_in, visc_rem, uhbt, uh_tot_0, duhdu_tot_0,        &
    !$omp&       du_max_CFL, du_min_CFL, do_I_in, IareaT, IareaT_xp1,       &
    !$omp&       dy_Cu, IdxT, IdxT_xp1)                                      &
    !$omp&   map(tofrom: du) map(from: uh_3d)

    do itt = 1, max_itts
      select case (itt)
        case (:1)    ; tol_eta = 1.0e-6_dp * tol_eta_base
        case (2)     ; tol_eta = 1.0e-4_dp * tol_eta_base
        case (3)     ; tol_eta = 1.0e-2_dp * tol_eta_base
        case default ; tol_eta = tol_eta_base
      end select

      !$omp distribute parallel do collapse(2) private(ii, jj)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        if     (uh_err(ii,jj) > 0.0_dp) then ; du_max(ii,jj) = du(ii,jj)
        elseif (uh_err(ii,jj) < 0.0_dp) then ; du_min(ii,jj) = du(ii,jj)
        else                                  ; do_I(ii,jj)   = .false.
        endif
      enddo ; enddo

      !$omp distribute parallel do collapse(2) private(ii, jj, ddu, du_prev)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        if (do_I(ii,jj)) then
          if (dt * min(IareaT(i,j), IareaT_xp1(i,j)) * abs(uh_err(ii,jj)) > tol_eta &
              .or. (better_iter .and.                                                   &
                    (abs(uh_err(ii,jj)) > tol_vel * duhdu_tot(ii,jj) .or.              &
                     abs(uh_err(ii,jj)) > uh_err_best(ii,jj)))) then
            ddu       = -uh_err(ii,jj) / duhdu_tot(ii,jj)
            du_prev   = du(ii,jj)
            du(ii,jj) = du(ii,jj) + ddu
            if (abs(ddu) < 1.0e-15_dp * abs(du(ii,jj))) then
              do_I(ii,jj) = .false.
            elseif (ddu > 0.0_dp) then
              if (du(ii,jj) >= du_max(ii,jj)) then
                du(ii,jj) = 0.5_dp * (du_prev + du_max(ii,jj))
                if (du_max(ii,jj) - du_prev < 1.0e-15_dp * abs(du(ii,jj))) &
                  do_I(ii,jj) = .false.
              endif
            else
              if (du(ii,jj) <= du_min(ii,jj)) then
                du(ii,jj) = 0.5_dp * (du_prev + du_min(ii,jj))
                if (du_prev - du_min(ii,jj) < 1.0e-15_dp * abs(du(ii,jj))) &
                  do_I(ii,jj) = .false.
              endif
            endif
          else
            do_I(ii,jj) = .false.
          endif
        endif
      enddo ; enddo

      !$omp distribute parallel do collapse(2) private(ii, jj)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        uh_err(ii,jj)    = -uhbt(ii,jj)
        duhdu_tot(ii,jj) = 0.0_dp
      enddo ; enddo

      do k = 1, nz
        !$omp distribute parallel do collapse(2) private(ii, jj, u_new, duhdu_loc)
        do j = j_start, j_end ; do i = i_start, i_end
          ii = i - i_start + 1 ; jj = j - j_start + 1
          if (do_I(ii,jj)) then
            u_new = u(i,j,k) + du(ii,jj) * visc_rem(ii,jj,k)
            call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                           dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,               &
                           uh_3d(ii,jj,k), duhdu_loc)
            uh_err(ii,jj)    = uh_err(ii,jj)    + uh_3d(ii,jj,k)
            duhdu_tot(ii,jj) = duhdu_tot(ii,jj) + duhdu_loc
          endif
        enddo ; enddo
      enddo

      !$omp distribute parallel do collapse(2) private(ii, jj)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        uh_err_best(ii,jj) = min(uh_err_best(ii,jj), abs(uh_err(ii,jj)))
      enddo ; enddo

    enddo ! itt

    !$omp end target teams

    !$omp target exit data &
    !$omp&   map(release: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)
  end subroutine zonal_flux_adjust_gpu

  ! CPU reference: identical logic, no OpenMP, domore early-exit enabled.
  subroutine zonal_flux_adjust_cpu(nx, ny, nz, ni, nj,                    &
                                    i_start, i_end, j_start, j_end,         &
                                    u, h_in, visc_rem,                       &
                                    uhbt, uh_tot_0, duhdu_tot_0,             &
                                    du_max_CFL, du_min_CFL, do_I_in,         &
                                    IareaT, IareaT_xp1,                      &
                                    dy_Cu, IdxT, IdxT_xp1,                   &
                                    dt, tol_eta_base, tol_vel, better_iter,  &
                                    du, uh_3d)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem(ni, nj, nz)
    real(dp), intent(in)  :: uhbt(ni, nj), uh_tot_0(ni, nj), duhdu_tot_0(ni, nj)
    real(dp), intent(in)  :: du_max_CFL(ni, nj), du_min_CFL(ni, nj)
    logical,  intent(in)  :: do_I_in(ni, nj)
    real(dp), intent(in)  :: IareaT(nx, ny), IareaT_xp1(nx, ny)
    real(dp), intent(in)  :: dy_Cu(nx, ny), IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt, tol_eta_base, tol_vel
    logical,  intent(in)  :: better_iter
    real(dp), intent(out) :: du(ni, nj), uh_3d(ni, nj, nz)

    real(dp) :: uh_err(ni, nj), uh_err_best(ni, nj)
    real(dp) :: duhdu_tot(ni, nj)
    real(dp) :: du_min(ni, nj), du_max(ni, nj)
    logical  :: do_I(ni, nj), domore
    integer  :: i, j, k, ii, jj, itt
    real(dp) :: tol_eta, u_new, duhdu_loc, ddu, du_prev

    do j = j_start, j_end ; do i = i_start, i_end
      ii = i - i_start + 1 ; jj = j - j_start + 1
      du(ii,jj)          = 0.0_dp  ;  do_I(ii,jj)      = do_I_in(ii,jj)
      du_max(ii,jj)      = du_max_CFL(ii,jj)
      du_min(ii,jj)      = du_min_CFL(ii,jj)
      uh_err(ii,jj)      = uh_tot_0(ii,jj) - uhbt(ii,jj)
      duhdu_tot(ii,jj)   = duhdu_tot_0(ii,jj)
      uh_err_best(ii,jj) = abs(uh_err(ii,jj))
    enddo ; enddo

    do itt = 1, max_itts
      select case (itt)
        case (:1)    ; tol_eta = 1.0e-6_dp * tol_eta_base
        case (2)     ; tol_eta = 1.0e-4_dp * tol_eta_base
        case (3)     ; tol_eta = 1.0e-2_dp * tol_eta_base
        case default ; tol_eta = tol_eta_base
      end select

      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        if     (uh_err(ii,jj) > 0.0_dp) then ; du_max(ii,jj) = du(ii,jj)
        elseif (uh_err(ii,jj) < 0.0_dp) then ; du_min(ii,jj) = du(ii,jj)
        else                                  ; do_I(ii,jj)   = .false.
        endif
      enddo ; enddo

      domore = .false.
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        if (do_I(ii,jj)) then
          if (dt * min(IareaT(i,j), IareaT_xp1(i,j)) * abs(uh_err(ii,jj)) > tol_eta &
              .or. (better_iter .and.                                                   &
                    (abs(uh_err(ii,jj)) > tol_vel * duhdu_tot(ii,jj) .or.              &
                     abs(uh_err(ii,jj)) > uh_err_best(ii,jj)))) then
            ddu       = -uh_err(ii,jj) / duhdu_tot(ii,jj)
            du_prev   = du(ii,jj)
            du(ii,jj) = du(ii,jj) + ddu
            if (abs(ddu) < 1.0e-15_dp * abs(du(ii,jj))) then
              do_I(ii,jj) = .false.
            elseif (ddu > 0.0_dp) then
              if (du(ii,jj) >= du_max(ii,jj)) then
                du(ii,jj) = 0.5_dp * (du_prev + du_max(ii,jj))
                if (du_max(ii,jj) - du_prev < 1.0e-15_dp * abs(du(ii,jj))) &
                  do_I(ii,jj) = .false.
              endif
            else
              if (du(ii,jj) <= du_min(ii,jj)) then
                du(ii,jj) = 0.5_dp * (du_prev + du_min(ii,jj))
                if (du_prev - du_min(ii,jj) < 1.0e-15_dp * abs(du(ii,jj))) &
                  do_I(ii,jj) = .false.
              endif
            endif
            if (do_I(ii,jj)) domore = .true.
          else
            do_I(ii,jj) = .false.
          endif
        endif
      enddo ; enddo
      if (.not. domore) exit

      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        uh_err(ii,jj)    = -uhbt(ii,jj)
        duhdu_tot(ii,jj) = 0.0_dp
      enddo ; enddo

      do k = 1, nz
        do j = j_start, j_end ; do i = i_start, i_end
          ii = i - i_start + 1 ; jj = j - j_start + 1
          if (do_I(ii,jj)) then
            u_new = u(i,j,k) + du(ii,jj) * visc_rem(ii,jj,k)
            call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                           dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,               &
                           uh_3d(ii,jj,k), duhdu_loc)
            uh_err(ii,jj)    = uh_err(ii,jj)    + uh_3d(ii,jj,k)
            duhdu_tot(ii,jj) = duhdu_tot(ii,jj) + duhdu_loc
          endif
        enddo ; enddo
      enddo

      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1
        uh_err_best(ii,jj) = min(uh_err_best(ii,jj), abs(uh_err(ii,jj)))
      enddo ; enddo

    enddo ! itt
  end subroutine zonal_flux_adjust_cpu

  ! ================================================================== !
  !  Test 3 – zonal flux adjustment, ij-outer GPU variant              !
  ! ================================================================== !

  ! GPU: single target teams, outer distribute parallel do over (j,i).
  ! The full do-itt Newton/bisection loop and do-k flux loop run serially
  ! inside each thread.  Per-column accumulators become scalar private
  ! variables — no persistent device allocations are required.
  subroutine zonal_flux_adjust_gpu_ij(nx, ny, nz, ni, nj,                    &
                                       i_start, i_end, j_start, j_end, nteams, &
                                       u, h_in, visc_rem,                       &
                                       uhbt, uh_tot_0, duhdu_tot_0,             &
                                       du_max_CFL, du_min_CFL, do_I_in,         &
                                       IareaT, IareaT_xp1,                      &
                                       dy_Cu, IdxT, IdxT_xp1,                   &
                                       dt, tol_eta_base, tol_vel, better_iter,  &
                                       du, uh_3d)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end, nteams
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem(ni, nj, nz)
    real(dp), intent(in)  :: uhbt(ni, nj), uh_tot_0(ni, nj), duhdu_tot_0(ni, nj)
    real(dp), intent(in)  :: du_max_CFL(ni, nj), du_min_CFL(ni, nj)
    logical,  intent(in)  :: do_I_in(ni, nj)
    real(dp), intent(in)  :: IareaT(nx, ny), IareaT_xp1(nx, ny)
    real(dp), intent(in)  :: dy_Cu(nx, ny), IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt, tol_eta_base, tol_vel
    logical,  intent(in)  :: better_iter
    real(dp), intent(out) :: du(ni, nj), uh_3d(ni, nj, nz)

    integer  :: i, j, k, ii, jj, itt
    ! Per-column scalars — all become private in the distribute loop.
    real(dp) :: uh_err, uh_err_best, duhdu_tot_loc, du_min_loc, du_max_loc, du_loc
    real(dp) :: tol_eta, u_new, duhdu_loc, ddu, du_prev
    logical  :: do_I_loc

    !$omp target teams num_teams(nteams)                                        &
    !$omp&   map(to: u, h_in, visc_rem, uhbt, uh_tot_0, duhdu_tot_0,            &
    !$omp&       du_max_CFL, du_min_CFL, do_I_in, IareaT, IareaT_xp1,           &
    !$omp&       dy_Cu, IdxT, IdxT_xp1)                                          &
    !$omp&   map(from: du, uh_3d)
    !$omp distribute parallel do collapse(2)                                      &
    !$omp&   private(ii, jj, uh_err, uh_err_best, duhdu_tot_loc,                 &
    !$omp&           du_min_loc, du_max_loc, du_loc, do_I_loc,                   &
    !$omp&           tol_eta, u_new, duhdu_loc, ddu, du_prev)
    do j = j_start, j_end ; do i = i_start, i_end
      ii = i - i_start + 1 ; jj = j - j_start + 1

      du_loc        = 0.0_dp
      do_I_loc      = do_I_in(ii,jj)
      du_max_loc    = du_max_CFL(ii,jj)
      du_min_loc    = du_min_CFL(ii,jj)
      uh_err        = uh_tot_0(ii,jj) - uhbt(ii,jj)
      duhdu_tot_loc = duhdu_tot_0(ii,jj)
      uh_err_best   = abs(uh_err)

      do itt = 1, max_itts
        select case (itt)
          case (:1)    ; tol_eta = 1.0e-6_dp * tol_eta_base
          case (2)     ; tol_eta = 1.0e-4_dp * tol_eta_base
          case (3)     ; tol_eta = 1.0e-2_dp * tol_eta_base
          case default ; tol_eta = tol_eta_base
        end select

        if     (uh_err > 0.0_dp) then ; du_max_loc = du_loc
        elseif (uh_err < 0.0_dp) then ; du_min_loc = du_loc
        else                          ; do_I_loc   = .false.
        endif

        if (do_I_loc) then
          if (dt * min(IareaT(i,j), IareaT_xp1(i,j)) * abs(uh_err) > tol_eta &
              .or. (better_iter .and.                                            &
                    (abs(uh_err) > tol_vel * duhdu_tot_loc .or.                 &
                     abs(uh_err) > uh_err_best))) then
            ddu     = -uh_err / duhdu_tot_loc
            du_prev = du_loc
            du_loc  = du_loc + ddu
            if (abs(ddu) < 1.0e-15_dp * abs(du_loc)) then
              do_I_loc = .false.
            elseif (ddu > 0.0_dp) then
              if (du_loc >= du_max_loc) then
                du_loc = 0.5_dp * (du_prev + du_max_loc)
                if (du_max_loc - du_prev < 1.0e-15_dp * abs(du_loc)) do_I_loc = .false.
              endif
            else
              if (du_loc <= du_min_loc) then
                du_loc = 0.5_dp * (du_prev + du_min_loc)
                if (du_prev - du_min_loc < 1.0e-15_dp * abs(du_loc)) do_I_loc = .false.
              endif
            endif
          else
            do_I_loc = .false.
          endif
        endif

        uh_err        = -uhbt(ii,jj)
        duhdu_tot_loc = 0.0_dp

        do k = 1, nz
          if (do_I_loc) then
            u_new = u(i,j,k) + du_loc * visc_rem(ii,jj,k)
            call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                           dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,               &
                           uh_3d(ii,jj,k), duhdu_loc)
            uh_err        = uh_err        + uh_3d(ii,jj,k)
            duhdu_tot_loc = duhdu_tot_loc + duhdu_loc
          endif
        enddo

        uh_err_best = min(uh_err_best, abs(uh_err))

      enddo ! itt

      du(ii,jj) = du_loc

    enddo ; enddo
    !$omp end target teams

  end subroutine zonal_flux_adjust_gpu_ij

  ! ================================================================== !
  !  Test 4 – zonal flux adjustment, fused-ij GPU variant             !
  ! ================================================================== !

  ! GPU: same do-itt outer loop as Test 2, but each iteration uses a
  ! single distribute parallel do collapse(2) that fuses the bounds
  ! update, Newton/bisection, reset, serial do-k flux accumulation, and
  ! best-error update into one loop body.  Shared work arrays still need
  ! target enter/exit data because they persist across do-itt iterations.
  subroutine zonal_flux_adjust_gpu_fused(nx, ny, nz, ni, nj,                    &
                                          i_start, i_end, j_start, j_end, nteams, &
                                          u, h_in, visc_rem,                       &
                                          uhbt, uh_tot_0, duhdu_tot_0,             &
                                          du_max_CFL, du_min_CFL, do_I_in,         &
                                          IareaT, IareaT_xp1,                      &
                                          dy_Cu, IdxT, IdxT_xp1,                   &
                                          dt, tol_eta_base, tol_vel, better_iter,  &
                                          du, uh_3d)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end, nteams
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem(ni, nj, nz)
    real(dp), intent(in)  :: uhbt(ni, nj), uh_tot_0(ni, nj), duhdu_tot_0(ni, nj)
    real(dp), intent(in)  :: du_max_CFL(ni, nj), du_min_CFL(ni, nj)
    logical,  intent(in)  :: do_I_in(ni, nj)
    real(dp), intent(in)  :: IareaT(nx, ny), IareaT_xp1(nx, ny)
    real(dp), intent(in)  :: dy_Cu(nx, ny), IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt, tol_eta_base, tol_vel
    logical,  intent(in)  :: better_iter
    real(dp), intent(out) :: du(ni, nj), uh_3d(ni, nj, nz)

    real(dp) :: uh_err(ni, nj), uh_err_best(ni, nj)
    real(dp) :: duhdu_tot(ni, nj)
    real(dp) :: du_min(ni, nj), du_max(ni, nj)
    logical  :: do_I(ni, nj)
    integer  :: i, j, k, ii, jj, itt
    real(dp) :: tol_eta, u_new, duhdu_loc, ddu, du_prev

    !$omp target enter data &
    !$omp&   map(alloc: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)

    !$omp target teams private(k,itt,tol_eta)                         &
    !$omp&   map(to: u, h_in, visc_rem, uhbt, uh_tot_0, duhdu_tot_0,            &
    !$omp&       du_max_CFL, du_min_CFL, do_I_in, IareaT, IareaT_xp1,           &
    !$omp&       dy_Cu, IdxT, IdxT_xp1)                                          &
    !$omp&   map(from: du, uh_3d)

    !$omp distribute parallel do collapse(2) private(ii, jj)
    do j = j_start, j_end ; do i = i_start, i_end
      ii = i - i_start + 1 ; jj = j - j_start + 1
      du(ii,jj)          = 0.0_dp  ;  do_I(ii,jj)      = do_I_in(ii,jj)
      du_max(ii,jj)      = du_max_CFL(ii,jj)
      du_min(ii,jj)      = du_min_CFL(ii,jj)
      uh_err(ii,jj)      = uh_tot_0(ii,jj) - uhbt(ii,jj)
      duhdu_tot(ii,jj)   = duhdu_tot_0(ii,jj)
      uh_err_best(ii,jj) = abs(uh_err(ii,jj))
    enddo ; enddo

    do itt = 1, max_itts
      select case (itt)
        case (:1)    ; tol_eta = 1.0e-6_dp * tol_eta_base
        case (2)     ; tol_eta = 1.0e-4_dp * tol_eta_base
        case (3)     ; tol_eta = 1.0e-2_dp * tol_eta_base
        case default ; tol_eta = tol_eta_base
      end select

      !$omp distribute parallel do collapse(2) &
      !$omp& private(ii, jj, ddu, du_prev, u_new, duhdu_loc)
      do j = j_start, j_end ; do i = i_start, i_end
        ii = i - i_start + 1 ; jj = j - j_start + 1

        ! Bounds update
        if     (uh_err(ii,jj) > 0.0_dp) then ; du_max(ii,jj) = du(ii,jj)
        elseif (uh_err(ii,jj) < 0.0_dp) then ; du_min(ii,jj) = du(ii,jj)
        else                                  ; do_I(ii,jj)   = .false.
        endif

        ! Newton / bisection step
        if (do_I(ii,jj)) then
          if (dt * min(IareaT(i,j), IareaT_xp1(i,j)) * abs(uh_err(ii,jj)) > tol_eta &
              .or. (better_iter .and.                                                   &
                    (abs(uh_err(ii,jj)) > tol_vel * duhdu_tot(ii,jj) .or.              &
                     abs(uh_err(ii,jj)) > uh_err_best(ii,jj)))) then
            ddu         = -uh_err(ii,jj) / duhdu_tot(ii,jj)
            du_prev     = du(ii,jj)
            du(ii,jj)   = du(ii,jj) + ddu
            if (abs(ddu) < 1.0e-15_dp * abs(du(ii,jj))) then
              do_I(ii,jj) = .false.
            elseif (ddu > 0.0_dp) then
              if (du(ii,jj) >= du_max(ii,jj)) then
                du(ii,jj) = 0.5_dp * (du_prev + du_max(ii,jj))
                if (du_max(ii,jj) - du_prev < 1.0e-15_dp * abs(du(ii,jj))) &
                  do_I(ii,jj) = .false.
              endif
            else
              if (du(ii,jj) <= du_min(ii,jj)) then
                du(ii,jj) = 0.5_dp * (du_prev + du_min(ii,jj))
                if (du_prev - du_min(ii,jj) < 1.0e-15_dp * abs(du(ii,jj))) &
                  do_I(ii,jj) = .false.
              endif
            endif
          else
            do_I(ii,jj) = .false.
          endif
        endif

        ! Reset accumulators
        uh_err(ii,jj)    = -uhbt(ii,jj)
        duhdu_tot(ii,jj) = 0.0_dp

        ! Serial layer flux accumulation
        do k = 1, nz
          if (do_I(ii,jj)) then
            u_new = u(i,j,k) + du(ii,jj) * visc_rem(ii,jj,k)
            call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                           dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,               &
                           uh_3d(ii,jj,k), duhdu_loc)
            uh_err(ii,jj)    = uh_err(ii,jj)    + uh_3d(ii,jj,k)
            duhdu_tot(ii,jj) = duhdu_tot(ii,jj) + duhdu_loc
          endif
        enddo

        ! Best error update
        uh_err_best(ii,jj) = min(uh_err_best(ii,jj), abs(uh_err(ii,jj)))

      enddo ; enddo

    enddo ! itt

    !$omp end target teams

    !$omp target exit data &
    !$omp&   map(release: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)
  end subroutine zonal_flux_adjust_gpu_fused

  ! ================================================================== !
  !  Test 5 – zonal flux adjustment, j-outer / i-parallel GPU variant  !
  ! ================================================================== !

  ! GPU: distribute over j (one team per j), parallel do over i inside.
  ! Loop order: distribute(j) → serial do-itt → serial do-k → parallel do(i).
  ! k is outside the inner i loops: each k level does a full parallel sweep
  ! over i before advancing to the next level.  Work arrays (uh_err etc.)
  ! are 2D (ni,nj) on device via target enter/exit data; each team
  ! accesses its own jj-slice, so there are no cross-team conflicts.
  subroutine zonal_flux_adjust_gpu_ji(nx, ny, nz, ni, nj,                    &
                                       i_start, i_end, j_start, j_end, nteams, &
                                       u, h_in, visc_rem,                       &
                                       uhbt, uh_tot_0, duhdu_tot_0,             &
                                       du_max_CFL, du_min_CFL, do_I_in,         &
                                       IareaT, IareaT_xp1,                      &
                                       dy_Cu, IdxT, IdxT_xp1,                   &
                                       dt, tol_eta_base, tol_vel, better_iter,  &
                                       du, uh_3d)
    integer,  intent(in)  :: nx, ny, nz, ni, nj
    integer,  intent(in)  :: i_start, i_end, j_start, j_end, nteams
    real(dp), intent(in)  :: u(nx, ny, nz)
    real(dp), intent(in)  :: h_in(0:nx+1, ny, nz)
    real(dp), intent(in)  :: visc_rem(ni, nj, nz)
    real(dp), intent(in)  :: uhbt(ni, nj), uh_tot_0(ni, nj), duhdu_tot_0(ni, nj)
    real(dp), intent(in)  :: du_max_CFL(ni, nj), du_min_CFL(ni, nj)
    logical,  intent(in)  :: do_I_in(ni, nj)
    real(dp), intent(in)  :: IareaT(nx, ny), IareaT_xp1(nx, ny)
    real(dp), intent(in)  :: dy_Cu(nx, ny), IdxT(nx, ny), IdxT_xp1(nx, ny)
    real(dp), intent(in)  :: dt, tol_eta_base, tol_vel
    logical,  intent(in)  :: better_iter
    real(dp), intent(out) :: du(ni, nj), uh_3d(ni, nj, nz)

    real(dp) :: uh_err(ni, nj), uh_err_best(ni, nj)
    real(dp) :: duhdu_tot(ni, nj)
    real(dp) :: du_min(ni, nj), du_max(ni, nj)
    logical  :: do_I(ni, nj)
    integer  :: i, j, k, ii, jj, itt
    real(dp) :: tol_eta, u_new, duhdu_loc, ddu, du_prev

    !$omp target enter data &
    !$omp&   map(alloc: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)

    !$omp target teams loop num_teams(nj)                               &
    !$omp&   map(to: u, h_in, visc_rem, uhbt, uh_tot_0, duhdu_tot_0,         &
    !$omp&       du_max_CFL, du_min_CFL, do_I_in, IareaT, IareaT_xp1,        &
    !$omp&       dy_Cu, IdxT, IdxT_xp1)                                       &
    !$omp&   map(from: du, uh_3d)
    do j = j_start, j_end
      jj = j - j_start + 1

      !$omp loop bind(parallel) private(i, ii)
      do i = i_start, i_end
        ii = i - i_start + 1
        du(ii,jj)          = 0.0_dp  ;  do_I(ii,jj)      = do_I_in(ii,jj)
        du_max(ii,jj)      = du_max_CFL(ii,jj)
        du_min(ii,jj)      = du_min_CFL(ii,jj)
        uh_err(ii,jj)      = uh_tot_0(ii,jj) - uhbt(ii,jj)
        duhdu_tot(ii,jj)   = duhdu_tot_0(ii,jj)
        uh_err_best(ii,jj) = abs(uh_err(ii,jj))
      enddo

      do itt = 1, max_itts
        select case (itt)
          case (:1)    ; tol_eta = 1.0e-6_dp * tol_eta_base
          case (2)     ; tol_eta = 1.0e-4_dp * tol_eta_base
          case (3)     ; tol_eta = 1.0e-2_dp * tol_eta_base
          case default ; tol_eta = tol_eta_base
        end select

        !$omp loop bind(parallel) private(i, ii)
        do i = i_start, i_end
          ii = i - i_start + 1
          if     (uh_err(ii,jj) > 0.0_dp) then ; du_max(ii,jj) = du(ii,jj)
          elseif (uh_err(ii,jj) < 0.0_dp) then ; du_min(ii,jj) = du(ii,jj)
          else                                  ; do_I(ii,jj)   = .false.
          endif
        enddo

        !$omp loop bind(parallel) private(i, ii, ddu, du_prev)
        do i = i_start, i_end
          ii = i - i_start + 1
          if (do_I(ii,jj)) then
            if (dt * min(IareaT(i,j), IareaT_xp1(i,j)) * abs(uh_err(ii,jj)) > tol_eta &
                .or. (better_iter .and.                                                   &
                      (abs(uh_err(ii,jj)) > tol_vel * duhdu_tot(ii,jj) .or.              &
                       abs(uh_err(ii,jj)) > uh_err_best(ii,jj)))) then
              ddu         = -uh_err(ii,jj) / duhdu_tot(ii,jj)
              du_prev     = du(ii,jj)
              du(ii,jj)   = du(ii,jj) + ddu
              if (abs(ddu) < 1.0e-15_dp * abs(du(ii,jj))) then
                do_I(ii,jj) = .false.
              elseif (ddu > 0.0_dp) then
                if (du(ii,jj) >= du_max(ii,jj)) then
                  du(ii,jj) = 0.5_dp * (du_prev + du_max(ii,jj))
                  if (du_max(ii,jj) - du_prev < 1.0e-15_dp * abs(du(ii,jj))) &
                    do_I(ii,jj) = .false.
                endif
              else
                if (du(ii,jj) <= du_min(ii,jj)) then
                  du(ii,jj) = 0.5_dp * (du_prev + du_min(ii,jj))
                  if (du_prev - du_min(ii,jj) < 1.0e-15_dp * abs(du(ii,jj))) &
                    do_I(ii,jj) = .false.
                endif
              endif
            else
              do_I(ii,jj) = .false.
            endif
          endif
        enddo

        !$omp loop bind(parallel) private(i, ii)
        do i = i_start, i_end
          ii = i - i_start + 1
          uh_err(ii,jj)    = -uhbt(ii,jj)
          duhdu_tot(ii,jj) = 0.0_dp
        enddo

        do k = 1, nz
          !$omp loop bind(parallel) private(i, ii, u_new, duhdu_loc)
          do i = i_start, i_end
            ii = i - i_start + 1
            if (do_I(ii,jj)) then
              u_new = u(i,j,k) + du(ii,jj) * visc_rem(ii,jj,k)
              call flux_elem(u_new, h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                             dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt,               &
                             uh_3d(ii,jj,k), duhdu_loc)
              uh_err(ii,jj)    = uh_err(ii,jj)    + uh_3d(ii,jj,k)
              duhdu_tot(ii,jj) = duhdu_tot(ii,jj) + duhdu_loc
            endif
          enddo
        enddo

        !$omp loop bind(parallel) private(i, ii)
        do i = i_start, i_end
          ii = i - i_start + 1
          uh_err_best(ii,jj) = min(uh_err_best(ii,jj), abs(uh_err(ii,jj)))
        enddo

      enddo ! itt

    enddo

    !$omp target exit data &
    !$omp&   map(release: uh_err, uh_err_best, duhdu_tot, du_min, du_max, do_I)
  end subroutine zonal_flux_adjust_gpu_ji

  ! ================================================================== !
  !  Shared comparison helpers                                          !
  ! ================================================================== !

  subroutine compare_3d(label, a_gpu, a_cpu, passed)
    character(*), intent(in)  :: label
    real(dp),     intent(in)  :: a_gpu(:,:,:), a_cpu(:,:,:)
    logical,      intent(out) :: passed
    integer :: i, j, k, ni, nj, nz, ndiff
    ni = size(a_gpu,1) ; nj = size(a_gpu,2) ; nz = size(a_gpu,3)
    ndiff = 0
    do k = 1, nz ; do j = 1, nj ; do i = 1, ni
      if (a_gpu(i,j,k) /= a_cpu(i,j,k)) ndiff = ndiff + 1
    enddo ; enddo ; enddo
    passed = (ndiff == 0)
    if (passed) then
      write(*,'(A,A)') '  PASS  ', label
    else
      write(*,'(A,A,I0,A,I0,A)') '  FAIL  ', label//' — ', ndiff, ' of ', ni*nj*nz, ' points differ'
      outer3: do k = 1, nz ; do j = 1, nj ; do i = 1, ni
        if (a_gpu(i,j,k) /= a_cpu(i,j,k)) then
          write(*,'(A,3I4)')       '    first diff (i,j,k)=', i, j, k
          write(*,'(A,E28.20)')    '    gpu = ', a_gpu(i,j,k)
          write(*,'(A,E28.20)')    '    cpu = ', a_cpu(i,j,k)
          write(*,'(A,Z16,A,Z16)') '    bits  gpu=', transfer(a_gpu(i,j,k), 0_8), &
                                    '  cpu=', transfer(a_cpu(i,j,k), 0_8)
          exit outer3
        endif
      enddo ; enddo ; enddo outer3
    endif
  end subroutine compare_3d

  subroutine compare_2d(label, a_gpu, a_cpu, passed)
    character(*), intent(in)  :: label
    real(dp),     intent(in)  :: a_gpu(:,:), a_cpu(:,:)
    logical,      intent(out) :: passed
    integer :: i, j, ni, nj, ndiff
    ni = size(a_gpu,1) ; nj = size(a_gpu,2)
    ndiff = 0
    do j = 1, nj ; do i = 1, ni
      if (a_gpu(i,j) /= a_cpu(i,j)) ndiff = ndiff + 1
    enddo ; enddo
    passed = (ndiff == 0)
    if (passed) then
      write(*,'(A,A)') '  PASS  ', label
    else
      write(*,'(A,A,I0,A,I0,A)') '  FAIL  ', label//' — ', ndiff, ' of ', ni*nj, ' points differ'
      outer2: do j = 1, nj ; do i = 1, ni
        if (a_gpu(i,j) /= a_cpu(i,j)) then
          write(*,'(A,2I4)')       '    first diff (i,j)=', i, j
          write(*,'(A,E28.20)')    '    gpu = ', a_gpu(i,j)
          write(*,'(A,E28.20)')    '    cpu = ', a_cpu(i,j)
          write(*,'(A,Z16,A,Z16)') '    bits  gpu=', transfer(a_gpu(i,j), 0_8), &
                                    '  cpu=', transfer(a_cpu(i,j), 0_8)
          exit outer2
        endif
      enddo ; enddo outer2
    endif
  end subroutine compare_2d

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

end module repro_mod


program test_repro
  use repro_mod
  use omp_lib
  implicit none

  integer,  parameter :: nz           = 100
  integer,  parameter :: n_sizes      = 2
  integer,  parameter :: n_runs       = 5
  integer,  parameter :: all_sizes(n_sizes) = [32, 64]
  real(dp), parameter :: dt           = 900.0_dp
  real(dp), parameter :: tol_eta_base = 1.0e-10_dp
  real(dp), parameter :: tol_vel      = 1.0e-10_dp
  real(dp), parameter :: pi           = 3.14159265358979323846_dp

  ! Allocatable arrays — resized for each problem size in the loop below.
  real(dp), allocatable :: u(:,:,:)
  real(dp), allocatable :: h_in(:,:,:)
  real(dp), allocatable :: visc_rem(:,:,:)
  real(dp), allocatable :: dy_Cu(:,:), IdxT(:,:), IdxT_xp1(:,:)
  real(dp), allocatable :: IareaT(:,:), IareaT_xp1(:,:)
  real(dp), allocatable :: uh_gpu(:,:,:), duhdu_gpu(:,:,:)
  real(dp), allocatable :: uh_cpu(:,:,:), duhdu_cpu(:,:,:)
  real(dp), allocatable :: uh_tot_0(:,:), duhdu_tot_0(:,:), uhbt(:,:)
  real(dp), allocatable :: du_max_CFL(:,:), du_min_CFL(:,:)
  logical,  allocatable :: do_I_in(:,:)
  real(dp), allocatable :: du_gpu(:,:), du_cpu(:,:)
  real(dp), allocatable :: du_gpuij(:,:), du_gpufu(:,:), du_gpuji(:,:)
  real(dp), allocatable :: uh3d_gpu(:,:,:), uh3d_cpu(:,:,:)
  real(dp), allocatable :: uh3d_gpuij(:,:,:), uh3d_gpufu(:,:,:), uh3d_gpuji(:,:,:)

  integer  :: nx, ny, ni, nj, nteams
  integer  :: i_start, i_end, j_start, j_end
  integer  :: i, j, k, ii, jj, isize, irun
  real(dp) :: tmp_uh, tmp_duhdu, t0, t1
  real(dp) :: times(n_runs)
  logical  :: p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, all_pass

  do isize = 1, n_sizes
    nx = all_sizes(isize) ; ny = all_sizes(isize)
    ni = nx ; nj = ny
    i_start = 1 ; i_end = nx
    j_start = 1 ; j_end = ny
    nteams  = (ni * nj + 255) / 256

    write(*,*)
    write(*,'(A)') '========================================================'
    write(*,'(A,I0,A,I0,A,I0,A,I0)') &
      'Size: ni=', ni, '  nj=', nj, '  nz=', nz, '  nteams=', nteams
    write(*,'(A)') '========================================================'

    allocate(u(nx, ny, nz))
    allocate(h_in(0:nx+1, ny, nz))
    allocate(visc_rem(ni, nj, nz))
    allocate(dy_Cu(nx, ny), IdxT(nx, ny), IdxT_xp1(nx, ny))
    allocate(IareaT(nx, ny), IareaT_xp1(nx, ny))
    allocate(uh_gpu(ni, nj, nz), duhdu_gpu(ni, nj, nz))
    allocate(uh_cpu(ni, nj, nz), duhdu_cpu(ni, nj, nz))
    allocate(uh_tot_0(ni, nj), duhdu_tot_0(ni, nj), uhbt(ni, nj))
    allocate(du_max_CFL(ni, nj), du_min_CFL(ni, nj))
    allocate(do_I_in(ni, nj))
    allocate(du_gpu(ni, nj), du_cpu(ni, nj))
    allocate(du_gpuij(ni, nj), du_gpufu(ni, nj), du_gpuji(ni, nj))
    allocate(uh3d_gpu(ni, nj, nz), uh3d_cpu(ni, nj, nz))
    allocate(uh3d_gpuij(ni, nj, nz), uh3d_gpufu(ni, nj, nz), uh3d_gpuji(ni, nj, nz))

    ! Initialise fields
    do k = 1, nz ; do j = 1, ny ; do i = 1, nx
      u(i,j,k)        = 0.5_dp  * sin(pi * real(i + j + k,     dp) / 20.0_dp)
      h_in(i,j,k)     = 10.0_dp + 2.0_dp * cos(pi * real(2*i + j + k, dp) / 30.0_dp)
      dy_Cu(i,j)      = 1000.0_dp + 10.0_dp * sin(pi * real(i + j,     dp) / 40.0_dp)
      IdxT(i,j)       = 1.0_dp / (1000.0_dp + 5.0_dp * sin(pi * real(i   + j, dp) / 50.0_dp))
      IdxT_xp1(i,j)   = 1.0_dp / (1000.0_dp + 5.0_dp * sin(pi * real(i+1 + j, dp) / 50.0_dp))
      IareaT(i,j)     = IdxT(i,j)     / dy_Cu(i,j)
      IareaT_xp1(i,j) = IdxT_xp1(i,j) / dy_Cu(i,j)
    enddo ; enddo ; enddo
    do k = 1, nz
      h_in(0,:,k)    = h_in(1,:,k)
      h_in(nx+1,:,k) = h_in(nx,:,k)
    enddo
    do k = 1, nz ; do j = j_start, j_end ; do i = i_start, i_end
      ii = i - i_start + 1 ; jj = j - j_start + 1
      visc_rem(ii,jj,k) = 0.9_dp + 0.1_dp * cos(pi * real(i + 2*j + k, dp) / 25.0_dp)
    enddo ; enddo ; enddo

    ! Compute initial summed fluxes (needed by tests 2-5)
    uh_tot_0(:,:)    = 0.0_dp
    duhdu_tot_0(:,:) = 0.0_dp
    do k = 1, nz ; do j = j_start, j_end ; do i = i_start, i_end
      ii = i - i_start + 1 ; jj = j - j_start + 1
      call flux_elem(u(i,j,k), h_in(i,j,k), h_in(i+1,j,k), visc_rem(ii,jj,k), &
                     dy_Cu(i,j), IdxT(i,j), IdxT_xp1(i,j), dt, tmp_uh, tmp_duhdu)
      uh_tot_0(ii,jj)    = uh_tot_0(ii,jj)    + tmp_uh
      duhdu_tot_0(ii,jj) = duhdu_tot_0(ii,jj) + tmp_duhdu
    enddo ; enddo ; enddo
    uhbt       = 0.9_dp * uh_tot_0
    du_max_CFL =  0.5_dp
    du_min_CFL = -0.5_dp
    do_I_in    = .true.

    ! -------------------------------------------------------------- !
    !  Correctness (one run per variant)                              !
    ! -------------------------------------------------------------- !
    write(*,*)
    write(*,*) '--- Correctness ---'

    call run_continuity_gpu(nx, ny, nz, ni, nj,                    &
                            i_start, i_end, j_start, j_end, nteams, &
                            u, h_in, visc_rem,                       &
                            dy_Cu, IdxT, IdxT_xp1, dt,               &
                            uh_gpu, duhdu_gpu)
    call run_continuity_cpu(nx, ny, nz, ni, nj,                    &
                            i_start, i_end, j_start, j_end,          &
                            u, h_in, visc_rem,                       &
                            dy_Cu, IdxT, IdxT_xp1, dt,               &
                            uh_cpu, duhdu_cpu)
    !$omp taskwait
    call compare_3d('Test1 uh_t  ', uh_gpu,    uh_cpu,    p1)
    call compare_3d('Test1 duhdu ', duhdu_gpu, duhdu_cpu, p2)

    call zonal_flux_adjust_gpu(nx, ny, nz, ni, nj,                         &
                                i_start, i_end, j_start, j_end, nteams,     &
                                u, h_in, visc_rem,                           &
                                uhbt, uh_tot_0, duhdu_tot_0,                 &
                                du_max_CFL, du_min_CFL, do_I_in,             &
                                IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                dt, tol_eta_base, tol_vel, .false.,          &
                                du_gpu, uh3d_gpu)
    call zonal_flux_adjust_cpu(nx, ny, nz, ni, nj,                         &
                                i_start, i_end, j_start, j_end,             &
                                u, h_in, visc_rem,                           &
                                uhbt, uh_tot_0, duhdu_tot_0,                 &
                                du_max_CFL, du_min_CFL, do_I_in,             &
                                IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                dt, tol_eta_base, tol_vel, .false.,          &
                                du_cpu, uh3d_cpu)
    !$omp taskwait
    call compare_2d('Test2 du    ', du_gpu,   du_cpu,   p3)
    call compare_3d('Test2 uh_3d ', uh3d_gpu, uh3d_cpu, p4)

    call zonal_flux_adjust_gpu_ij(nx, ny, nz, ni, nj,                         &
                                   i_start, i_end, j_start, j_end, nteams,     &
                                   u, h_in, visc_rem,                           &
                                   uhbt, uh_tot_0, duhdu_tot_0,                 &
                                   du_max_CFL, du_min_CFL, do_I_in,             &
                                   IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                   dt, tol_eta_base, tol_vel, .false.,          &
                                   du_gpuij, uh3d_gpuij)
    !$omp taskwait
    call compare_2d('Test3 du ij ', du_gpuij,   du_cpu,   p5)
    call compare_3d('Test3 uh ij ', uh3d_gpuij, uh3d_cpu, p6)

    call zonal_flux_adjust_gpu_fused(nx, ny, nz, ni, nj,                         &
                                      i_start, i_end, j_start, j_end, nteams,     &
                                      u, h_in, visc_rem,                           &
                                      uhbt, uh_tot_0, duhdu_tot_0,                 &
                                      du_max_CFL, du_min_CFL, do_I_in,             &
                                      IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                      dt, tol_eta_base, tol_vel, .false.,          &
                                      du_gpufu, uh3d_gpufu)
    !$omp taskwait
    call compare_2d('Test4 du fu ', du_gpufu,   du_cpu,   p7)
    call compare_3d('Test4 uh fu ', uh3d_gpufu, uh3d_cpu, p8)

    call zonal_flux_adjust_gpu_ji(nx, ny, nz, ni, nj,                         &
                                   i_start, i_end, j_start, j_end, nteams,     &
                                   u, h_in, visc_rem,                           &
                                   uhbt, uh_tot_0, duhdu_tot_0,                 &
                                   du_max_CFL, du_min_CFL, do_I_in,             &
                                   IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                   dt, tol_eta_base, tol_vel, .false.,          &
                                   du_gpuji, uh3d_gpuji)
    !$omp taskwait
    call compare_2d('Test5 du ji ', du_gpuji,   du_cpu,   p9)
    call compare_3d('Test5 uh ji ', uh3d_gpuji, uh3d_cpu, p10)

    all_pass = p1 .and. p2 .and. p3 .and. p4 .and. p5 .and. p6 .and. p7 .and. p8 .and. p9 .and. p10
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

    write(*,*) 'Test 1 GPU (continuity):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_continuity_gpu(nx, ny, nz, ni, nj,                    &
                              i_start, i_end, j_start, j_end, nteams, &
                              u, h_in, visc_rem,                       &
                              dy_Cu, IdxT, IdxT_xp1, dt,               &
                              uh_gpu, duhdu_gpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'Test 1 CPU (continuity):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call run_continuity_cpu(nx, ny, nz, ni, nj,                    &
                              i_start, i_end, j_start, j_end,          &
                              u, h_in, visc_rem,                       &
                              dy_Cu, IdxT, IdxT_xp1, dt,               &
                              uh_cpu, duhdu_cpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'Test 2 GPU (ij teams, separate loops):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call zonal_flux_adjust_gpu(nx, ny, nz, ni, nj,                         &
                                  i_start, i_end, j_start, j_end, nteams,     &
                                  u, h_in, visc_rem,                           &
                                  uhbt, uh_tot_0, duhdu_tot_0,                 &
                                  du_max_CFL, du_min_CFL, do_I_in,             &
                                  IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                  dt, tol_eta_base, tol_vel, .false.,          &
                                  du_gpu, uh3d_gpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'Test 2 CPU:'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call zonal_flux_adjust_cpu(nx, ny, nz, ni, nj,                         &
                                  i_start, i_end, j_start, j_end,             &
                                  u, h_in, visc_rem,                           &
                                  uhbt, uh_tot_0, duhdu_tot_0,                 &
                                  du_max_CFL, du_min_CFL, do_I_in,             &
                                  IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                  dt, tol_eta_base, tol_vel, .false.,          &
                                  du_cpu, uh3d_cpu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'Test 3 GPU (ij-outer, scalar private):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call zonal_flux_adjust_gpu_ij(nx, ny, nz, ni, nj,                         &
                                     i_start, i_end, j_start, j_end, nteams,     &
                                     u, h_in, visc_rem,                           &
                                     uhbt, uh_tot_0, duhdu_tot_0,                 &
                                     du_max_CFL, du_min_CFL, do_I_in,             &
                                     IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                     dt, tol_eta_base, tol_vel, .false.,          &
                                     du_gpuij, uh3d_gpuij)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'Test 4 GPU (fused-ij):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call zonal_flux_adjust_gpu_fused(nx, ny, nz, ni, nj,                         &
                                        i_start, i_end, j_start, j_end, nteams,     &
                                        u, h_in, visc_rem,                           &
                                        uhbt, uh_tot_0, duhdu_tot_0,                 &
                                        du_max_CFL, du_min_CFL, do_I_in,             &
                                        IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                        dt, tol_eta_base, tol_vel, .false.,          &
                                        du_gpufu, uh3d_gpufu)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    write(*,*) 'Test 5 GPU (ji-outer, parallel i):'
    do irun = 1, n_runs
      t0 = omp_get_wtime()
      call zonal_flux_adjust_gpu_ji(nx, ny, nz, ni, nj,                         &
                                     i_start, i_end, j_start, j_end, nteams,     &
                                     u, h_in, visc_rem,                           &
                                     uhbt, uh_tot_0, duhdu_tot_0,                 &
                                     du_max_CFL, du_min_CFL, do_I_in,             &
                                     IareaT, IareaT_xp1, dy_Cu, IdxT, IdxT_xp1,  &
                                     dt, tol_eta_base, tol_vel, .false.,          &
                                     du_gpuji, uh3d_gpuji)
      !$omp taskwait
      t1 = omp_get_wtime()
      times(irun) = t1 - t0
    enddo
    call print_timing_stats(times)

    deallocate(u, h_in, visc_rem, dy_Cu, IdxT, IdxT_xp1, IareaT, IareaT_xp1)
    deallocate(uh_gpu, duhdu_gpu, uh_cpu, duhdu_cpu)
    deallocate(uh_tot_0, duhdu_tot_0, uhbt, du_max_CFL, du_min_CFL, do_I_in)
    deallocate(du_gpu, du_cpu, du_gpuij, du_gpufu, du_gpuji)
    deallocate(uh3d_gpu, uh3d_cpu, uh3d_gpuij, uh3d_gpufu, uh3d_gpuji)

  enddo ! isize

  write(*,*)
  write(*,*) 'Done.'

end program test_repro
