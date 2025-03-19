!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions on nonlinear analysis

module m_fstr_QuasiNewton

  use m_fstr_NonLinearMethod

  implicit none
  ! parameters for line search
  real(kind=kreal), parameter :: C_line_search=2.0, Psi0_line_search=1.0d-2
  real(kind=kreal), parameter :: delta_wolfe=0.2, sigma_wolfe=0.9, eps_wolfe=1.0d-3
  real(kind=kreal), parameter :: omega_wolfe=0.001, Delta_approx_wolfe=0.7
  real(kind=kreal) :: C_wolfe, Q_Wolfe
  integer, parameter :: n_mem_max=10
  integer(kind=kint), parameter :: maxiter_ls = 10
  integer(kind=kint), parameter :: pot_type=1

contains

  !> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
  !> method
  subroutine fstr_Quasi_Newton( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, &
    restrt_step_num, sub_step, ctime, dtime )
    implicit none

    integer, intent(in)                   :: cstep     !< current loading step
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer, intent(in)                   :: sub_step  !< substep number of current loading step
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    type (hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange

    type (hecmwST_local_mesh), pointer :: hecMESHmpc
    type (hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: ndof
    integer(kind=kint) :: i, iter
    integer(kind=kint) :: stepcnt
    integer(kind=kint) :: restrt_step_num
    real(kind=kreal)   :: tt0, tt, res, qnrm, rres, tincr, xnrm, dunrm, rxnrm
    logical :: isLinear = .false.
    integer(kind=kint) :: iterStatus

    real(kind=kreal), allocatable :: z_k(:), s_k(:,:), y_k(:,:), g_prev(:), rho_k(:)
    real(kind=kreal) :: sdoty
    integer :: n_mem
    integer :: len_vector
    integer(kind=kint) :: k

    integer :: u_debug
    integer :: max_iter_bak
    logical :: flag_approx_Wolfe

    max_iter_bak = fstrSOLID%step_ctrl(cstep)%max_iter
    fstrSOLID%step_ctrl(cstep)%max_iter = 100*fstrSOLID%step_ctrl(cstep)%max_iter

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    if(.not. fstrPR%nlgeom)then
      isLinear = .true.
    endif

    hecMAT%NDOF = hecMESH%n_dof
    NDOF = hecMAT%NDOF

    stepcnt = 0

    tincr = dtime
    if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0
    call fstr_init_Newton(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, hecLagMAT, ndof)
    fstrSOLID%GL0(:) = fstrSOLID%GL(:) !store external load at du=0

    !! initialize du for non-zero Dirichlet condition
    call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 1, RHSvector=fstrSOLID%dunode)
    !! update residual vector
    call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)

    len_vector = hecMESH%n_node*ndof
    allocate(z_k(len_vector))
    allocate(s_k(len_vector, n_mem_max))
    allocate(y_k(len_vector, n_mem_max))
    allocate(g_prev(len_vector))
    allocate(rho_k(n_mem_max))
    z_k(:) = 0.0d0
    s_k(:,:) = 0.0d0
    y_k(:,:) = 0.0d0
    g_prev(:) = 0.0d0
    rho_k(:) = 0.0d0
    do i=1,len_vector
      y_k(i,1) = -hecMAT%B(i)
    enddo

    ! parameter to judge Wolfe/approx Wolfe selection
    C_wolfe = 0.0d0
    Q_wolfe = 0.0d0
    flag_approx_Wolfe = .false.
    n_mem = 1
    ! ----- Inner Iteration, lagrange multiplier constant
    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt = stepcnt+1

      ! ----- calculate search direction by limited BFGS method
      do i=1,hecMESH%n_node*ndof
          g_prev(i) = -hecMAT%B(i)
      enddo
      call fstr_calc_direction_LBFGS(hecMesh, g_prev, s_k, y_k, rho_k, z_k, n_mem)

      ! ! ----- Set Boundary condition
      call fstr_AddBC_to_direction_vector(z_k, hecMESH,fstrSOLID, cstep)

      !----- line search of step length
      call fstr_line_search_along_direction(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k)

      ! ----- update the small displacement and the displacement for 1step
      !       \delta u^k => solver's solution
      !       \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i = 1, hecMESH%n_node*ndof
        fstrSOLID%dunode(i) = fstrSOLID%dunode(i) + hecMAT%X(i)
      enddo

      !! set du for non-zero Dirichlet condition
      ! call fstr_AddBC(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 1, RHSvector=fstrSOLID%dunode)
      call fstr_calc_residual_vector(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
 
      ! ----- check convergence
      iterStatus = fstr_check_iteration_converged(hecMESH, hecMAT, fstrSOLID, ndof, iter, sub_step, cstep )
      if (iterStatus == kitrConverged) exit
      if (iterStatus == kitrDiverged .or. iterStatus==kitrFloatingError) return
      ! if (iterStatus == kitrDiverged) exit
      ! if (iterStatus==kitrFloatingError) return

      n_mem = min(n_mem+1, n_mem_max)
      do k=n_mem, 2, -1
        s_k(:,k) = s_k(:,k-1)
        y_k(:,k) = y_k(:,k-1)
        rho_k(k) = rho_k(k-1)
      enddo

      do i = 1, hecMESH%n_node*ndof
        s_k(i,1) = hecMAT%X(i)
        y_k(i,1) = -hecMAT%B(i) - g_prev(i)
      enddo
      call hecmw_innerProduct_R(hecMESH,ndof,s_k(:,1), y_k(:,1), sdoty)
      if (abs(sdoty) < 1.0d-10) then
        rho_k(1) = 0.0d0
      else
        rho_k(1) = 1.0d0/sdoty
      endif
    enddo
    ! ----- end of inner loop

    fstrSOLID%NRstat_i(knstMAXIT) = max(fstrSOLID%NRstat_i(knstMAXIT),iter) ! logging newton iteration(maxtier)
    fstrSOLID%NRstat_i(knstSUMIT) = fstrSOLID%NRstat_i(knstSUMIT) + iter    ! logging newton iteration(sum of iter)

    ! ----- update the total displacement
    ! u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i) = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    fstrSOLID%CutBack_stat = 0
    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    fstrSOLID%step_ctrl(cstep)%max_iter = max_iter_bak
  end subroutine fstr_Quasi_Newton

  subroutine fstr_calc_direction_LBFGS(hecMESH, g_prev, s_k, y_k, rho_k, z_k, n_mem)
    implicit none

    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    real(kind=kreal) :: z_k(:), s_k(:,:), y_k(:,:), g_prev(:), rho_k(:)
    integer :: n_mem

    real(kind=kreal), allocatable :: q(:)
    real(kind=kreal) :: alpha(n_mem), beta
    real(kind=kreal) :: sdotq, ysq, gamma, ydotz, g_max
    
    integer :: len_vector, ndof
    integer(kind=kint) :: k,i

    ndof = hecMesh%n_dof
    len_vector = hecMESH%n_node*hecMesh%n_dof
    allocate(q(len_vector))

    q(1:len_vector) = g_prev(1:len_vector)

    do k=1, n_mem
      call hecmw_innerProduct_R(hecMESH,ndof,s_k(:,k), q, sdotq)
      alpha(k) = rho_k(k) * sdotq
      do i=1, len_vector
        q(i) = q(i) - alpha(k)*y_k(i,k)
      enddo
    enddo
    call hecmw_innerProduct_R(hecMESH,ndof,y_k(:,1), y_k(:,1), ysq)
    if (n_mem==1) then
      call hecmw_absMax_R(hecMESH, ndof, g_prev, g_max)
      ! if (g_max==0.0d0) then
      !   write(6,*) 'gradient of potential is zero-vector'
      !   stop
      ! endif
      ! gamma = 1.0d0/g_max
      gamma = 1.0d0
    else if (abs(rho_k(1)) < 1.0d-10) then
      gamma = 1.0d0
    else
      gamma = 1.0d0/(rho_k(1)*ysq)
    endif

    do i=1, len_vector
      z_k(i) = gamma*q(i)
    enddo

    do k=n_mem, 1, -1
      call hecmw_innerProduct_R(hecMESH,ndof,y_k(:,k), z_k, ydotz)
      beta = rho_k(k)*ydotz
      do i=1, len_vector
        z_k(i)=z_k(i) + s_k(i,k)*(alpha(k)-beta)
      enddo
    enddo
    deallocate(q)
  end subroutine fstr_calc_direction_LBFGS

  subroutine fstr_AddBC_to_direction_vector(z_k, hecMESH,fstrSOLID, cstep)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    integer(kind=kint) :: cstep
    real(kind=kreal) :: z_k(:)

    integer(kind=kint) :: ig0, grpid, ig, ityp, idofS, idofE, iS0, iE0, ik, in, idof, ndof

    ndof = hecMesh%n_dof
    !   ----- Prescibed displacement Boundary Conditions
    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      ityp = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      idofS = ityp/10
      idofE = ityp - idofS*10
      !
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      !
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        !
        do idof = idofS, idofE
            z_k(ndof*(in-1)+idof) = 0.0d0
        enddo
      enddo
    enddo
  end subroutine fstr_AddBC_to_direction_vector

  subroutine fstr_apply_alpha0(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, h_prime, pot)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: z_k(:)
    real(kind=kreal) :: h_prime, pot

    hecMat%X(:) = 0.0d0
    call hecmw_innerProduct_R(hecMESH, hecMAT%NDOF, hecMat%B, z_k, h_prime)
    pot = fstr_get_potential(cstep,hecMESH,hecMAT,fstrSOLID,pot_type)
  end subroutine fstr_apply_alpha0

  subroutine fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha, h_prime, pot)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: z_k(:)
    real(kind=kreal), intent(in) :: alpha
    real(kind=kreal) :: h_prime, pot

    hecMat%X(:) = -alpha*z_k(:)
    call fstr_calc_residual_vector_with_X(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM)
    call hecmw_innerProduct_R(hecMESH, hecMAT%NDOF, hecMat%B, z_k, h_prime)
    pot = fstr_get_potential_with_X(cstep,hecMESH,hecMAT,fstrSOLID,pot_type)
  end subroutine fstr_apply_alpha

  subroutine fstr_init_line_search_range(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, &
     h_prime_0, pot_0, alpha_S, h_prime_S, pot_S,  alpha_E, h_prime_E, pot_E)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: z_k(:)
    real(kind=kreal), intent(in) :: h_prime_0, pot_0
    real(kind=kreal) :: alpha_S, h_prime_S, pot_S
    real(kind=kreal) :: alpha_E, h_prime_E, pot_E

    real(kind=kreal) :: alpha_S_new, h_prime_S_new, pot_S_new
    real(kind=kreal) :: alpha_E_new, h_prime_E_new, pot_E_new
    real(kind=kreal) :: alpha_tmp, h_prime_tmp, pot_tmp
    real(kind=kreal) :: z_max
    real(kind=kreal) :: pot_eps
    pot_eps = eps_wolfe*C_Wolfe

    alpha_S = 0.0d0
    h_prime_S = h_prime_0
    pot_S = pot_0

    if (iter==1) then
      call hecmw_absMax_R(hecMESH, hecMAT%NDOF, z_k, z_max)
      alpha_E = Psi0_line_search/z_max
    else
      alpha_E = 1.0d0
    end if
    call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha_E, h_prime_E, pot_E)

    do while (h_prime_E < 0.0d0)
      ! if (pot_E <= pot_0 + pot_eps) then
        ! alpha_S = alpha_E
        ! h_prime_S = h_prime_E ! so h_prime_S < 0
        ! pot_S = h_prime_S

        alpha_E = alpha_E * C_line_search
        call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, &
                                     &  fstrPARAM, z_k, alpha_E, h_prime_E, pot_E)
      ! else
      !   alpha_tmp = 2.0d0*alpha_E
      !   h_prime_tmp = 0.0d0 ! dummy value
      !   pot_tmp = 0.0d0 ! dummy value
      !   call fstr_get_new_range_with_potential(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, pot_0, &
      !     alpha_S, h_prime_S, pot_S, alpha_tmp, h_prime_tmp, pot_tmp, alpha_E, h_prime_E, pot_E, &
      !     alpha_S_new, h_prime_S_new, pot_S_new, alpha_E_new, h_prime_E_new, pot_E_new)

      !   alpha_S = alpha_S_new
      !   h_prime_S = h_prime_S_new
      !   pot_S = pot_S_new

      !   alpha_E = alpha_E_new
      !   h_prime_E = h_prime_E_new
      !   pot_E = pot_E_new
      !   return
      ! end if
    enddo
  end subroutine fstr_init_line_search_range

  subroutine fstr_update_range(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, &
     alpha_S, h_prime_S, pot_S,  alpha_E, h_prime_E, pot_E,  alpha_c, h_prime_c, pot_c, h_prime_0, pot_0)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal) :: z_k(:)
    real(kind=kreal) :: alpha_S, h_prime_S, pot_S
    real(kind=kreal) :: alpha_E, h_prime_E, pot_E
    real(kind=kreal) :: alpha_c, h_prime_c, pot_c
    real(kind=kreal) :: h_prime_0, pot_0

    real(kind=kreal) :: alpha_A, h_prime_A, pot_A
    real(kind=kreal) :: alpha_B, h_prime_B, pot_B
    real(kind=kreal) :: alpha_a_bar, h_prime_a_bar, pot_a_bar
    real(kind=kreal) :: alpha_b_bar, h_prime_b_bar, pot_b_bar
    real(kind=kreal) :: alpha_c_bar, h_prime_c_bar, pot_c_bar

    call fstr_get_new_range_with_potential(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, pot_0, &
      alpha_S, h_prime_S, pot_S, alpha_E, h_prime_E, pot_E, alpha_c, h_prime_c, pot_c, &
      alpha_A, h_prime_A, pot_A, alpha_B, h_prime_B, pot_B)

    if (alpha_c == alpha_A) then
      call fstr_get_secant(alpha_S, h_prime_S, alpha_A, h_prime_A, alpha_c_bar)
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, &
        alpha_c_bar, h_prime_c_bar, pot_c_bar)
    end if
    if (alpha_c == alpha_B) then
      call fstr_get_secant(alpha_B, h_prime_B, alpha_E, h_prime_E, alpha_c_bar)
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, &
        alpha_c_bar, h_prime_c_bar, pot_c_bar)
    end if

    if ((alpha_c == alpha_A) .or. (alpha_c == alpha_B)) then
      call fstr_get_new_range_with_potential(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, pot_0, &
        alpha_A, h_prime_A, pot_A, alpha_B, h_prime_B, pot_B, alpha_c_bar, h_prime_c_bar, pot_c_bar, &
        alpha_a_bar, h_prime_a_bar, pot_a_bar, alpha_b_bar, h_prime_b_bar, pot_b_bar)
    else
      alpha_a_bar = alpha_A
      h_prime_a_bar = h_prime_A
      pot_a_bar = pot_A

      alpha_b_bar = alpha_B
      h_prime_b_bar = h_prime_B
      pot_b_bar = pot_B
    end if


    if ((alpha_b_bar - alpha_a_bar) > 0.66*(alpha_E - alpha_S)) then
      alpha_c_bar = 0.5*(alpha_a_bar+alpha_b_bar)
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, &
        alpha_c_bar, h_prime_c_bar, pot_c_bar)
      call fstr_get_new_range_with_potential(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, pot_0, &
        alpha_A, h_prime_A, pot_A, alpha_B, h_prime_B, pot_B, alpha_c_bar, h_prime_c_bar, pot_c_bar, &
        alpha_S, h_prime_S, pot_S, alpha_E, h_prime_E, pot_E)
    else
      alpha_S = alpha_a_bar
      h_prime_S = h_prime_a_bar
      pot_S = pot_a_bar

      alpha_E = alpha_b_bar
      h_prime_E = h_prime_b_bar
      pot_E = pot_b_bar
    end if
  end subroutine fstr_update_range

  subroutine fstr_get_new_range_with_potential(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, &
    fstrPARAM, z_k, pot_0, alpha_a, h_prime_a, pot_a, alpha_b, h_prime_b, pot_b, alpha_c, h_prime_c, pot_c, &
     alpha_a_bar, h_prime_a_bar, pot_a_bar, alpha_b_bar, h_prime_b_bar, pot_b_bar)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal), intent(in) :: z_k(:)
    real(kind=kreal), intent(in) :: pot_0
    real(kind=kreal), intent(in) :: alpha_a, h_prime_a, pot_a
    real(kind=kreal), intent(in) :: alpha_b, h_prime_b, pot_b
    real(kind=kreal), intent(in) :: alpha_c, h_prime_c, pot_c
    real(kind=kreal) :: alpha_a_bar, h_prime_a_bar, pot_a_bar
    real(kind=kreal) :: alpha_b_bar, h_prime_b_bar, pot_b_bar

    integer, parameter :: count_max=100
    integer :: count_while
    real(kind=kreal) :: alpha_d, h_prime_d, pot_d
    real(kind=kreal) :: theta_ls = 0.5d0
    real(kind=kreal) :: pot_eps
    pot_eps = eps_wolfe*C_Wolfe

    ! case of NOT (a < c < b)
    if (alpha_c <= alpha_a .or. alpha_b <= alpha_c) then
      alpha_a_bar = alpha_a
      h_prime_a_bar = h_prime_a
      pot_a_bar = pot_a

      alpha_b_bar = alpha_b
      h_prime_b_bar = h_prime_b
      pot_b_bar = pot_b
      return
    end if

    if (h_prime_c >= 0.0d0) then
      alpha_a_bar = alpha_a
      h_prime_a_bar = h_prime_a
      pot_a_bar = pot_a

      alpha_b_bar = alpha_c
      h_prime_b_bar = h_prime_c
      pot_b_bar = pot_c
      return
    end if
    ! below here, it can be assumed that h_prime_c<0

    if(pot_c <= pot_0 + pot_eps) then
      alpha_a_bar = alpha_c
      h_prime_a_bar = h_prime_c
      pot_a_bar = pot_c

      alpha_b_bar = alpha_b
      h_prime_b_bar = h_prime_b
      pot_b_bar = pot_b
      return
    end if

    alpha_a_bar = alpha_a
    h_prime_a_bar = h_prime_a
    pot_a_bar = pot_a

    alpha_b_bar = alpha_b
    h_prime_b_bar = h_prime_b
    pot_b_bar = pot_b

    count_while = 0
    do while(count_while < count_max)
      count_while = count_while + 1

      alpha_d = (1.0d0-theta_ls)*alpha_a_bar + theta_ls*alpha_b_bar
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha_d, h_prime_d, pot_d)

      if (h_prime_d >= 0.0d0) then
        alpha_b_bar = alpha_d
        h_prime_b_bar = h_prime_d
        pot_b_bar = pot_d
        return
      end if

      if(pot_d <= pot_0 + pot_eps) then
        alpha_a_bar = alpha_d
        h_prime_a_bar = h_prime_d
        pot_a_bar = pot_d
      else
        alpha_b_bar = alpha_d
        h_prime_b_bar = h_prime_d
        pot_b_bar = pot_d
      end if
    end do
    write(6,*) 'fstr_get_new_range_with_potential reached loop count max.', hecMESH%my_rank, alpha_a_bar, alpha_b_bar
  end subroutine fstr_get_new_range_with_potential


  subroutine fstr_line_search_along_direction(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k)
    implicit none
    type (hecmwST_local_mesh)             :: hecMESH   !< hecmw mesh
    type (hecmwST_matrix)                 :: hecMAT    !< hecmw matrix
    type (fstr_solid)                     :: fstrSOLID !< fstr_solid
    real(kind=kreal), intent(in)          :: ctime     !< current time
    real(kind=kreal), intent(in) :: tincr
    integer(kind=kint) :: iter
    integer, intent(in)                   :: cstep     !< current loading step
    real(kind=kreal), intent(in)          :: dtime     !< time increment
    type (fstr_param)                     :: fstrPARAM !< type fstr_param
    real(kind=kreal) :: z_k(:)

    real(kind=kreal) :: alpha_S, h_prime_S, pot_S
    real(kind=kreal) :: alpha_E, h_prime_E, pot_E
    real(kind=kreal) :: alpha_c, h_prime_c, pot_c
    real(kind=kreal) :: h_prime_0, pot_0
    logical :: flag_converged
    integer :: ndof, len_vector
    real(kind=kreal) :: res

    integer(kind=kint) :: i, ierr, iter_ls
    real(kind=kreal) :: z_max
    integer :: dummy

    ndof = hecMAT%NDOF
    len_vector = hecMESH%n_node*hecMesh%n_dof

    call fstr_apply_alpha0(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, h_prime_0, pot_0)
    if (h_prime_0 > 0.0d0) then
      write(6,*) 'residual vector is not directed to potential decretion.', h_prime_0
      stop
    endif

    call fstr_init_line_search_range(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, &
    fstrPARAM, z_k, h_prime_0, pot_0, alpha_S, h_prime_S, pot_S,  alpha_E, h_prime_E, pot_E)
    if( hecMESH%my_rank == 0 ) then
      write(6,'(a, 6es27.16e3)') 'range initialized: alpha_S, alpha_E, h_prime_S, h_prime_E, pot_S, pot_E', &
          &  alpha_S, alpha_E, h_prime_S, h_prime_E, pot_S, pot_E
    endif    

		Q_Wolfe = 1 + Delta_approx_wolfe * Q_Wolfe
		C_Wolfe = C_Wolfe + (abs(pot_0)-C_Wolfe) / Q_Wolfe

    flag_converged = .false.
    iter_ls=0
    do while (iter_ls<maxiter_ls)
      call fstr_get_secant(alpha_S, h_prime_S, alpha_E, h_prime_E, alpha_c)
      call fstr_apply_alpha(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, alpha_c, h_prime_c, pot_c)

      if (abs( pot_c - pot_0) <= (omega_Wolfe * C_Wolfe) ) then
        flag_converged = fstr_approx_wolfe_condition(alpha_c, h_prime_c, pot_c, h_prime_0, pot_0)
      else
        flag_converged = fstr_wolfe_condition(alpha_c, h_prime_c, pot_c, h_prime_0, pot_0)
      endif
      if (flag_converged) exit

      call fstr_update_range(hecMESH, hecMAT, fstrSOLID, ctime, tincr, iter, cstep, dtime, fstrPARAM, z_k, &
        alpha_S, h_prime_S, pot_S,  alpha_E, h_prime_E, pot_E,  alpha_c, h_prime_c, pot_c, h_prime_0, pot_0)
      iter_ls = iter_ls +1
    enddo
    if( hecMESH%my_rank == 0 ) then
      write(6,'(a, 1i8, 5es27.16e3)') 'converged: alpha_S, alpha_E, alpha_c, h_prime_c, pot_c', &
          &  iter_ls, &
          &  alpha_S, alpha_E, alpha_c, h_prime_c, pot_c
    endif    
  end subroutine fstr_line_search_along_direction

  subroutine fstr_get_secant(a, Fa, b, Fb, c)
    implicit none
    real(kind=kreal), intent(in) :: a,b, Fa,Fb
    real(kind=kreal), intent(out) :: c
    if(Fb /= Fa) then
      c = (a*Fb  - b*Fa) / ( Fb - Fa )
    else
      c = 0.5*a + 0.5*b
    endif
  end subroutine fstr_get_secant

  function fstr_wolfe_condition(alpha_c, h_prime_c, pot_c, h_prime_0, pot_0) result(flag_converged)
    implicit none
    logical :: flag_converged

    real(kind=kreal), intent(in) :: alpha_c, h_prime_0, h_prime_c, pot_0, pot_c
    real(kind=kreal) :: wolfe1_left, wolfe1_right, wolfe2_left, wolfe2_right

    wolfe1_left = pot_c - pot_0
    wolfe1_right = delta_wolfe*h_prime_0*alpha_c

    wolfe2_left = h_prime_c
    wolfe2_right = sigma_wolfe*h_prime_0

    flag_converged = (wolfe1_left<=wolfe1_right) .and. (wolfe2_left >= wolfe2_right)
    if (flag_converged) write(6,'(a, 4es27.16e3)') 'oWolfe: ', wolfe1_left, wolfe1_right, wolfe2_left, wolfe2_right
  end function fstr_wolfe_condition

  function fstr_approx_wolfe_condition(alpha_c, h_prime_c, pot_c, h_prime_0, pot_0) result(flag_converged)
    implicit none
    logical :: flag_converged

    real(kind=kreal), intent(in) :: alpha_c, h_prime_0, h_prime_c, pot_0, pot_c
    real(kind=kreal) :: wolfe1_left, wolfe1_right, wolfe2_left, wolfe2_right

    real(kind=kreal) :: pot_eps
    pot_eps = eps_wolfe*C_Wolfe

    wolfe1_left =  ( 2.0 * delta_wolfe - 1.0d0 ) * h_prime_0
    wolfe1_right = h_prime_c

    wolfe2_left = h_prime_c
    wolfe2_right = sigma_wolfe*h_prime_0

    flag_converged = &
      (wolfe1_left>=wolfe1_right) &
      .and. (wolfe2_left >= wolfe2_right) &
      .and. (pot_c <= pot_0 + pot_eps)
    if (flag_converged) write(6,'(a, 6es27.16e3)') 'aWolfe: ', wolfe1_left, wolfe1_right, &
    &  wolfe2_left, wolfe2_right, pot_c, pot_0 + (eps_wolfe*C_Wolfe)
  end function fstr_approx_wolfe_condition
end module m_fstr_QuasiNewton
