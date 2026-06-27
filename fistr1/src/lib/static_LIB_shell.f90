!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module ...
module m_static_LIB_shell
  use hecmw, only : kint, kreal
  use elementInfo

  !--------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------

  !--------------------------------------------------------------
  !
  ! (Programmer)
  ! Gaku Hashimoto
  ! Department of Human and Engineered Environmental Studies
  ! Graduate School of Frontier Sciences, The University of Tokyo
  ! 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-8563 JAPAN
  !
  ! (Ref.)
  ! [1] Noguchi, H. and Hisada, T.,
  !     "Sensitivity analysis in post-buckling problems of shell
  !     structures,"
  !     Computers & Structures, Vol.47, No.4, pp.699-710, (1993).
  ! [2] Dvorkin, E.N. and Bathe, K.J.,
  !     "A Continuum Mechanics Based Four-node Shell Element for
  !     General Non-linear Analysis,"
  !     Engineering Computations, Vol.1, pp.77-88, (1984).
  ! [3] Bucalem, M.L. and Bathe, K.J.,
  !     "Higher-order MITC general shell element,"
  !     International Journal for Numerical Methods in
  !     Engineering, Vol.36, pp.3729-3754, (1993).
  ! [4] Lee, P.S. and Bathe, K.J.,
  !     "Development of MITC Isotropic Triangular Shell Finite
  !     Elements,"
  !     Computers & Structures, Vol.82, pp.945-962, (2004).
  !
  ! Xi YUAN
  !   Apr. 13, 2019: Introduce mass matrix calculation
  ! (Ref.)
  ! [5] E.Hinton, T.A.Rock, O.C.Zienkiewicz(1976): A Note on Mass Lumping
  ! and Related Process in FEM. International Journal on Earthquake Eng
  ! and structural dynamics, 4, pp245-249
  !
  !--------------------------------------------------------------

  !--------------------------------------------------------------------

contains

  pure function outer_product3(a, b) result(ab)

    real(kind = kreal), intent(in) :: a(3), b(3)
    real(kind = kreal) :: ab(3, 3)
    integer :: i, j

    do j = 1, 3
      do i = 1, 3
        ab(i, j) = a(i)*b(j)
      end do
    end do

  end function outer_product3


  pure subroutine shell_basis_from_covariant(g1, g2, g3, e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, det)

    real(kind = kreal), intent(in) :: g1(3), g2(3), g3(3)
    real(kind = kreal), intent(out) :: e1_hat(3), e2_hat(3), e3_hat(3)
    real(kind = kreal), intent(out) :: cg1(3), cg2(3), cg3(3)
    real(kind = kreal), intent(out) :: det
    real(kind = kreal) :: det_inv, g3_abs, e1_abs, e2_abs

    det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
      +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
      +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )
    det_inv = 1.0D0/det

    cg1(1) = det_inv*( g2(2)*g3(3)-g2(3)*g3(2) )
    cg1(2) = det_inv*( g2(3)*g3(1)-g2(1)*g3(3) )
    cg1(3) = det_inv*( g2(1)*g3(2)-g2(2)*g3(1) )
    cg2(1) = det_inv*( g3(2)*g1(3)-g3(3)*g1(2) )
    cg2(2) = det_inv*( g3(3)*g1(1)-g3(1)*g1(3) )
    cg2(3) = det_inv*( g3(1)*g1(2)-g3(2)*g1(1) )
    cg3(1) = det_inv*( g1(2)*g2(3)-g1(3)*g2(2) )
    cg3(2) = det_inv*( g1(3)*g2(1)-g1(1)*g2(3) )
    cg3(3) = det_inv*( g1(1)*g2(2)-g1(2)*g2(1) )

    g3_abs = dsqrt( dot_product(g3, g3) )
    e3_hat(:) = g3(:)/g3_abs
    e1_hat(1) = g2(2)*e3_hat(3)-g2(3)*e3_hat(2)
    e1_hat(2) = g2(3)*e3_hat(1)-g2(1)*e3_hat(3)
    e1_hat(3) = g2(1)*e3_hat(2)-g2(2)*e3_hat(1)
    e1_abs = dsqrt( dot_product(e1_hat, e1_hat) )
    e1_hat(:) = e1_hat(:)/e1_abs
    e2_hat(1) = e3_hat(2)*e1_hat(3)-e3_hat(3)*e1_hat(2)
    e2_hat(2) = e3_hat(3)*e1_hat(1)-e3_hat(1)*e1_hat(3)
    e2_hat(3) = e3_hat(1)*e1_hat(2)-e3_hat(2)*e1_hat(1)
    e2_abs = dsqrt( dot_product(e2_hat, e2_hat) )
    e2_hat(:) = e2_hat(:)/e2_abs

  end subroutine shell_basis_from_covariant

  !--------------------------------------------------------------------
  subroutine ShellStressVectorToTensor( stress, tensor )
    implicit none

    real(kind=kreal), intent(in)  :: stress(6)
    real(kind=kreal), intent(out) :: tensor(3, 3)

    tensor(:, :) = 0.0D0
    tensor(1, 1) = stress(1)
    tensor(2, 2) = stress(2)
    tensor(3, 3) = stress(3)
    tensor(1, 2) = stress(4)
    tensor(2, 1) = tensor(1, 2)
    tensor(2, 3) = stress(5)
    tensor(3, 2) = tensor(2, 3)
    tensor(3, 1) = stress(6)
    tensor(1, 3) = tensor(3, 1)

  end subroutine ShellStressVectorToTensor

  subroutine ShellTensorToStressVector( tensor, stress )
    implicit none

    real(kind=kreal), intent(in)  :: tensor(3, 3)
    real(kind=kreal), intent(out) :: stress(6)

    stress(1) = tensor(1, 1)
    stress(2) = tensor(2, 2)
    stress(3) = tensor(3, 3)
    stress(4) = tensor(1, 2)
    stress(5) = tensor(2, 3)
    stress(6) = tensor(3, 1)

  end subroutine ShellTensorToStressVector

  !--------------------------------------------------------------------
  subroutine ShellObjectiveStressIncrement( stress_old, dstrain, dstress_obj, trace_coeff )
    implicit none

    real(kind=kreal), intent(in)  :: stress_old(6), dstrain(6)
    real(kind=kreal), intent(out) :: dstress_obj(6)
    real(kind=kreal), intent(in), optional :: trace_coeff

    call ShellObjectiveTraceStressIncrement( stress_old, dstrain, dstress_obj, trace_coeff )

  end subroutine ShellObjectiveStressIncrement

  !--------------------------------------------------------------------
  subroutine ShellObjectiveTraceStressIncrement( stress_old, dstrain, dstress_trace, trace_coeff )
    implicit none

    real(kind=kreal), intent(in)  :: stress_old(6), dstrain(6)
    real(kind=kreal), intent(out) :: dstress_trace(6)
    real(kind=kreal), intent(in), optional :: trace_coeff

    real(kind=kreal) :: stress_tensor(3, 3), dstress_tensor(3, 3)
    real(kind=kreal) :: trD, coeff

    call ShellStressVectorToTensor( stress_old, stress_tensor )
    coeff = 1.0D0
    if( present( trace_coeff ) ) coeff = trace_coeff
    trD = coeff*(dstrain(1)+dstrain(2))+dstrain(3)
    dstress_tensor(:, :) = -stress_tensor(:, :)*trD
    call ShellTensorToStressVector( dstress_tensor, dstress_trace )

  end subroutine ShellObjectiveTraceStressIncrement

  !--------------------------------------------------------------------
  real(kind=kreal) function ShellPlaneStressTraceCoeff( gauss, n_layer )
    use mMechGauss
    use mMaterial, only: M_POISSON, MC_ISOELASTIC, getElasticType
    implicit none

    type(tGaussStatus), intent(in) :: gauss
    integer(kind=kint), intent(in) :: n_layer

    real(kind=kreal) :: nu, outa(2)
    logical :: ierr

    ! Isotropic plane-stress trace coefficient for the UL stress update.
    ShellPlaneStressTraceCoeff = 1.0D0
    if( .not. associated( gauss%pMaterial ) ) return
    if( getElasticType( gauss%pMaterial%mtype ) == 1 ) then
      stop "MITC4 shell UL orthotropic trace correction is not supported"
    endif
    nu = gauss%pMaterial%variables(M_POISSON)
    call fetch_TableData(MC_ISOELASTIC, gauss%pMaterial%dict, outa, ierr)
    if( associated( gauss%pMaterial%shell_var ) ) then
      if( n_layer >= 1 .and. n_layer <= size( gauss%pMaterial%shell_var ) ) then
        if( gauss%pMaterial%shell_var(n_layer)%ortho == 0 ) then
          if( ierr ) then
            nu = gauss%pMaterial%shell_var(n_layer)%pp
          else
            nu = outa(2)
          endif
        else
          stop "MITC4 shell UL orthotropic trace correction is not supported"
        endif
      else if( .not. ierr ) then
        nu = outa(2)
      endif
    else if( .not. ierr ) then
      nu = outa(2)
    endif

    if( abs(1.0D0-nu) > 1.0D-12 ) then
      ShellPlaneStressTraceCoeff = (1.0D0-2.0D0*nu)/(1.0D0-nu)
    endif

  end function ShellPlaneStressTraceCoeff

  !--------------------------------------------------------------------
  subroutine ShellAddStressVectorToDB( dstress, j, DB )
    implicit none

    integer(kind=kint), intent(in) :: j
    real(kind=kreal), intent(in)    :: dstress(6)
    real(kind=kreal), intent(inout) :: DB(:, :)

    DB(1, j) = DB(1, j)+dstress(1)
    DB(2, j) = DB(2, j)+dstress(2)
    DB(3, j) = DB(3, j)+dstress(4)
    DB(4, j) = DB(4, j)+dstress(5)
    DB(5, j) = DB(5, j)+dstress(6)

  end subroutine ShellAddStressVectorToDB

  !--------------------------------------------------------------------
  subroutine ShellAddULObjectiveTraceTangent( stress_old, ncol, B, DB, trace_coeff )
    implicit none

    integer(kind=kint), intent(in) :: ncol
    real(kind=kreal), intent(in)    :: stress_old(6)
    real(kind=kreal), intent(in)    :: B(5, ncol)
    real(kind=kreal), intent(inout) :: DB(5, ncol)
    real(kind=kreal), intent(in), optional :: trace_coeff

    integer(kind=kint) :: j
    real(kind=kreal) :: dstrain_col(6), dstress_trace(6)

    do j = 1, ncol
      dstrain_col(:) = 0.0D0
      dstrain_col(1) = B(1, j)
      dstrain_col(2) = B(2, j)
      call ShellObjectiveTraceStressIncrement( stress_old, dstrain_col, dstress_trace, trace_coeff )
      call ShellAddStressVectorToDB( dstress_trace, j, DB )
    end do

  end subroutine ShellAddULObjectiveTraceTangent

  !--------------------------------------------------------------------
  subroutine ShellAddULObjectiveTangent( stress_old, ncol, B, DB, trace_coeff )
    implicit none

    integer(kind=kint), intent(in) :: ncol
    real(kind=kreal), intent(in)    :: stress_old(6)
    real(kind=kreal), intent(in)    :: B(5, ncol)
    real(kind=kreal), intent(inout) :: DB(5, ncol)
    real(kind=kreal), intent(in), optional :: trace_coeff

    call ShellAddULObjectiveTraceTangent( stress_old, ncol, B, DB, trace_coeff )

  end subroutine ShellAddULObjectiveTangent


  !####################################################################
  subroutine STF_Shell_MITC                                           &
      (etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, nddisp, element, qf_stress, include_geo_stiff, &
      nddirector, ndrefdirector, nddrill)
    !####################################################################

    use mMechGauss
    use Quadrature
    use m_MatMatrix
    use mMaterial, only: INFINITESIMAL, UPDATELAG, TOTALLAG, isElastic

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn, mixflag
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(out)  :: stiff(:, :)
    real(kind = kreal), intent(in)   :: thick

    real(kind = kreal), intent(in), optional :: nddisp(ndof, nn)
    type(tElement), intent(in), optional :: element
    real(kind = kreal), intent(out), optional :: qf_stress(:)
    logical, intent(in), optional :: include_geo_stiff
    real(kind = kreal), intent(in), optional :: nddirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndrefdirector(3, nn)
    real(kind = kreal), intent(in), optional :: nddrill(nn)

    !--------------------------------------------------------------------

    integer :: flag, flag_dof
    integer :: ndof_shell
    integer(kind=kint) :: ierr_quad
    integer(kind=kint) :: ishell
    integer :: i, j, m, n
    integer :: lx, ly
    integer :: fetype
    integer :: ny
    integer :: ntying
    integer :: npoints_tying(3)
    integer :: it, ip
    integer :: na, nb
    integer :: isize, jsize
    integer :: jsize1, jsize2, jsize3, &
      jsize4, jsize5, jsize6
    integer :: n_layer,n_totlyr, sstable(24)

    real(kind = kreal) :: D(5, 5), B(5, ndof*nn), DB(5, ndof*nn)
    real(kind = kreal) :: stress_old_vec(6)
    real(kind = kreal) :: tmpstiff(ndof*nn, ndof*nn)
    real(kind = kreal) :: qf_tmp(ndof*nn), qf_mix(ndof*nn), qf_disp(ndof*nn)
    real(kind = kreal) :: qf_stress_integrand(ndof*nn), vol_deriv(ndof*nn)
    real(kind = kreal) :: Sv_force(5), geo_term
    real(kind = kreal) :: S_global(3, 3), Smat(9, 9)
    real(kind = kreal) :: BN(9, ndof*nn), SBN(9, ndof*nn)
    real(kind = kreal) :: elem(3, nn)
    real(kind = kreal) :: shell_disp(6, nn)
    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: w_w_lx, w_ly
    real(kind = kreal) :: B_di(5, ndof*nn, 6, 3, 7)
    real(kind = kreal) :: BG1_di(3, ndof*nn, 6, 3, 7)
    real(kind = kreal) :: BG2_di(3, ndof*nn, 6, 3, 7)
    real(kind = kreal) :: BG3_di(3, ndof*nn, 6, 3, 7)
    real(kind = kreal) :: B1(3, ndof*nn), B2(3, ndof*nn), &
      B3(3, ndof*nn)
    real(kind = kreal), allocatable :: B2rot(:, :, :, :)
    real(kind = kreal), allocatable :: B2rot_di(:, :, :, :, :, :, :)
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: tpcoord(6, 2, 3)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: aa1(3), aa2(3), aa3(3)
    real(kind = kreal) :: bb1(3), bb2(3), bb3(3)
    real(kind = kreal) :: cc1(3), cc2(3)
    real(kind = kreal) :: alpha
    real(kind = kreal) :: trace_coeff
    real(kind = kreal) :: xxi_lx, eeta_lx
    real(kind = kreal) :: xxi_di(6, 3), eeta_di(6, 3)
    real(kind = kreal) :: h(nn, 3)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_i(3), v2_i(3), v3_i(3)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: a_over_2_v3_ref(3, nn)
    real(kind = kreal) :: a_over_2_v3_deriv(3, 3, nn)
    real(kind = kreal) :: a_over_2_v3_second(3, 3, 3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
      dudzeta_rot(3, nn)
    real(kind = kreal) :: dudxi_rot_deriv(3, 3), dudeta_rot_deriv(3, 3), &
      dudzeta_rot_deriv(3, 3)
    real(kind = kreal) :: dudxi_rot_second(3, 3, 3), dudeta_rot_second(3, 3, 3), &
      dudzeta_rot_second(3, 3, 3)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: g1_tl(3), g2_tl(3)
    real(kind = kreal) :: g1_weight(3), g2_weight(3), g3_weight(3)
    real(kind = kreal) :: dudxi_trans(3), dudeta_trans(3)
    real(kind = kreal) :: g3_abs
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: cg1(3), cg2(3), cg3(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_weight
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_inv
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: w_w_w_det
    real(kind = kreal) :: e1_hat(3), e2_hat(3), e3_hat(3)
    real(kind = kreal) :: e1_hat_abs, e2_hat_abs
    real(kind = kreal) :: e1_hat_mat(3), e2_hat_mat(3), e3_hat_mat(3)
    real(kind = kreal) :: cg1_mat(3), cg2_mat(3), cg3_mat(3)
    real(kind = kreal) :: det_mat
    real(kind = kreal) :: Cv12(ndof*nn), Cv13(ndof*nn), &
      Cv21(ndof*nn), Cv23(ndof*nn), &
      Cv31(ndof*nn), Cv32(ndof*nn)
    real(kind = kreal) :: Cv_theta(ndof*nn), Cv_w(ndof*nn)
    real(kind = kreal) :: Cv(ndof*nn)
    real(kind = kreal) :: Cv_w_second(3, 3, nn)
    real(kind = kreal) :: Cv_disp, Cv_deriv, Cv_deriv_disp
    real(kind = kreal) :: Cv12_2, Cv13_2, Cv21_2, Cv23_2, Cv31_2, Cv32_2
    real(kind = kreal) :: drill_axis(3), rot_projector(3, 3), hrot(3, 3)
    real(kind = kreal) :: hess_coeff, drill_coeff, axis_norm
    real(kind = kreal) :: director_inc(3), director_ref(3), director_deriv(3, 3)
    real(kind = kreal) :: director_second(3, 3, 3)
    logical :: finite_rotation_director, add_geo_stiff, use_tl_green, use_director_tangent

    sstable = 0
    flag_dof = 0
    ny = 0
    shell_disp(:, :) = 0.0d0
    B_di(:, :, :, :, :) = 0.0D0
    B1(:, :) = 0.0D0
    B2(:, :) = 0.0D0
    B3(:, :) = 0.0D0
    a_over_2_v3_second(:, :, :, :) = 0.0d0
    ndof_shell = min(ndof, 6)

    !--------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      fetype = fe_mitc4_shell

      ny = 2

      ntying = 1
      npoints_tying(1)= 4

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      fetype = fe_mitc9_shell

      ny = 3

      ntying = 3
      npoints_tying(1)= 6
      npoints_tying(2)= 6
      npoints_tying(3)= 4

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      fetype = fe_mitc3_shell

      ny = 2

      ntying = 1
      npoints_tying(1)= 3

    end if

    !--------------------------------------------------------------------

    if( present( nddisp ) ) then

      shell_disp(1:ndof_shell, 1:nn) = nddisp(1:ndof_shell, 1:nn)

    end if

    !--------------------------------------------------------------------

    flag = gausses(1)%pMaterial%nlgeom_flag

    if( .not. present( nddisp ) ) flag = INFINITESIMAL
    finite_rotation_director = ( flag == TOTALLAG .or. flag == UPDATELAG ) .and. ndof_shell >= 6 &
      .and. etype == fe_mitc4_shell .and. nn == 4 &
      .and. isElastic(gausses(1)%pMaterial%mtype)
    use_director_tangent = finite_rotation_director .and. etype == fe_mitc4_shell .and. nn == 4 &
      .and. isElastic(gausses(1)%pMaterial%mtype)
    add_geo_stiff = flag /= INFINITESIMAL
    if( present( include_geo_stiff ) ) add_geo_stiff = add_geo_stiff .and. include_geo_stiff
    ! Green-Lagrange strain for elastic MITC4
    use_tl_green = flag == TOTALLAG .and. etype == fe_mitc4_shell .and. nn == 4 &
      .and. isElastic(gausses(1)%pMaterial%mtype)
    if( use_director_tangent ) then
      allocate( B2rot(5, 3, 3, nn) )
      allocate( B2rot_di(5, 3, 3, nn, 6, 3, ny) )
      B2rot_di(:, :, :, :, :, :, :) = 0.0D0
    endif

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG ) elem(:, :) = elem(:, :) + shell_disp(1:3, :)
    BG1_di(:, :, :, :, :) = 0.0D0
    BG2_di(:, :, :, :, :) = 0.0D0
    BG3_di(:, :, :, :, :) = 0.0D0

    !--------------------------------------------------------------------

    tmpstiff(:, :) = 0.0D0
    qf_tmp(:) = 0.0D0
    qf_disp(:) = 0.0D0
    do nb = 1, nn
      qf_disp(ndof*(nb-1)+1:ndof*(nb-1)+ndof_shell) = shell_disp(1:ndof_shell, nb)
    end do
    if( present( qf_stress ) ) then
      qf_stress(1:ndof*nn) = 0.0D0
    endif

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    !-------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) =  0.0D0
      tpcoord(2, 1, 1) =  1.0D0
      tpcoord(3, 1, 1) =  0.0D0
      tpcoord(4, 1, 1) = -1.0D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0
      tpcoord(2, 2, 1) =  0.0D0
      tpcoord(3, 2, 1) =  1.0D0
      tpcoord(4, 2, 1) =  0.0D0

      !--------------------------------------------------------

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      !--------------------------------------------------------

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 1)  = -1.0D0
      xxi_di(2, 1)  =  1.0D0
      xxi_di(3, 1)  =  1.0D0
      xxi_di(4, 1)  = -1.0D0
      xxi_di(5, 1)  =  1.0D0
      xxi_di(6, 1)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 1) = -1.0D0
      eeta_di(2, 1) = -1.0D0
      eeta_di(3, 1) =  1.0D0
      eeta_di(4, 1) =  1.0D0
      eeta_di(5, 1) =  0.0D0
      eeta_di(6, 1) =  0.0D0

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 2)  = -1.0D0
      xxi_di(2, 2)  =  0.0D0
      xxi_di(3, 2)  =  1.0D0
      xxi_di(4, 2)  =  1.0D0
      xxi_di(5, 2)  =  0.0D0
      xxi_di(6, 2)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 2) = -1.0D0
      eeta_di(2, 2) = -1.0D0
      eeta_di(3, 2) = -1.0D0
      eeta_di(4, 2) =  1.0D0
      eeta_di(5, 2) =  1.0D0
      eeta_di(6, 2) =  1.0D0

      !--------------------------------------------------------

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1)  = 0.5D0
      tpcoord(2, 1, 1)  = 0.0D0
      tpcoord(3, 1, 1)  = 0.5D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = 0.0D0
      tpcoord(2, 2, 1) = 0.5D0
      tpcoord(3, 2, 1) = 0.5D0

      !--------------------------------------------------------

    end if

    !--------------------------------------------------------------------

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    do i = 1, 3

      g1(i) = 0.0D0

      do na = 1, nn

        g1(i) = g1(i)+shapederiv(na, 1) &
          *elem(i, na)

      end do

    end do

    e_0(1) = g1(1)
    e_0(2) = g1(2)
    e_0(3) = g1(3)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1) &
            *elem(i, na)
          g2(i) = g2(i)+shapederiv(na, 2) &
            *elem(i, na)

        end do

      end do

      !------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(1, nb) = det_cg3(1)/det_cg3_abs
      v3(2, nb) = det_cg3(2)/det_cg3_abs
      v3(3, nb) = det_cg3(3)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs .GT. 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else    ! YX: impossible

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !---------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)
      a_over_2_v3_ref(1:3, nb) = a_over_2_v3(1:3, nb)
      if( finite_rotation_director .and. present( ndrefdirector ) ) then
        ! reference nodal director
        a_over_2_v3_ref(1:3, nb) = ndrefdirector(1:3, nb)
      endif

      director_ref(1:3) = a_over_2_v3_ref(1:3, nb)
      if( finite_rotation_director ) then
        if( present( nddirector ) ) then
          a_over_2_v3(1:3, nb) = nddirector(1:3, nb)
        else
          call ShellDirectorIncrement( shell_disp(4:6, nb), director_ref, director_inc )
          a_over_2_v3(1:3, nb) = director_ref(1:3) + director_inc(1:3)
        endif
        if( use_director_tangent ) then
          call ShellDirectorIncrementalDeriv( a_over_2_v3(1:3, nb), director_deriv )
          call ShellDirectorIncrementalSecondDeriv( a_over_2_v3(1:3, nb), director_second )
          a_over_2_v3_deriv(1:3, 1:3, nb) = director_deriv(1:3, 1:3)
          a_over_2_v3_second(1:3, 1:3, 1:3, nb) = director_second(1:3, 1:3, 1:3)
        endif
      endif
      if( .not. use_director_tangent ) then
        a_over_2_v3_deriv(1:3, 1, nb) = (/ 0.0D0, -a_over_2_v3(3, nb), a_over_2_v3(2, nb) /)
        a_over_2_v3_deriv(1:3, 2, nb) = (/ a_over_2_v3(3, nb), 0.0D0, -a_over_2_v3(1, nb) /)
        a_over_2_v3_deriv(1:3, 3, nb) = (/ -a_over_2_v3(2, nb), a_over_2_v3(1, nb), 0.0D0 /)
      endif

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------
    !     MODIFIED to LAMINATED SHELL ANALYSIS
    !--------------------------------------------------------------------
    zeta_ly = 0.0D0

    n_totlyr =  gausses(1)%pMaterial%totallyr
    do n_layer=1,n_totlyr
      do ly = 1, ny

        !--------------------------------------------------------

        ! MITC4
        if( etype .EQ. fe_mitc4_shell ) then

          zeta_ly = 0.0D0

          ! MITC9
        else if( etype .EQ. fe_mitc9_shell ) then

          zeta_ly = gauss1d3(1,ly)

          ! MITC3
        else if( etype .EQ. fe_mitc3_shell )then

          zeta_ly = 0.0D0

        end if

        !--------------------------------------------------------

        do it = 1, ntying

          do ip = 1, npoints_tying(it)

            !-------------------------------------------------

            naturalcoord(1) = tpcoord(ip, 1, it)
            naturalcoord(2) = tpcoord(ip, 2, it)

            call getShapeFunc(fetype, naturalcoord, shapefunc)

            call getShapeDeriv(fetype, naturalcoord, shapederiv)

            !-------------------------------------------------

            do na = 1, nn

              do i = 1, 3

                u_rot(i, na)                      &
                  = shapefunc(na)                   &
                  *( zeta_ly*a_over_2_v3(i, na) )

                dudxi_rot(i, na)                  &
                  = shapederiv(na, 1)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudeta_rot(i, na)                 &
                  = shapederiv(na, 2)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudzeta_rot(i, na)                &
                  = shapefunc(na)                   &
                  *( a_over_2_v3(i, na) )

              end do

            end do

            !-------------------------------------------------

            ! Covariant basis vector
            do i = 1, 3

              g1(i) = 0.0D0
              g2(i) = 0.0D0
              g3(i) = 0.0D0

              do na = 1, nn

                g1(i) = g1(i)+shapederiv(na, 1)  &
                  *elem(i, na)       &
                  +dudxi_rot(i, na)
                g2(i) = g2(i)+shapederiv(na, 2)  &
                  *elem(i, na)       &
                  +dudeta_rot(i, na)
                g3(i) = g3(i)+dudzeta_rot(i, na)

              end do

            end do

            !-------------------------------------------------

            g1_tl(:) = g1(:)
            g2_tl(:) = g2(:)
            if( use_tl_green ) then
              dudxi_trans(:) = 0.0D0
              dudeta_trans(:) = 0.0D0
              do na = 1, nn
                dudxi_trans(:) = dudxi_trans(:)+shapederiv(na, 1)*shell_disp(1:3, na)
                dudeta_trans(:) = dudeta_trans(:)+shapederiv(na, 2)*shell_disp(1:3, na)
              end do
              if( etype == fe_mitc9_shell ) then
                g1_tl(:) = g1_tl(:)+dudxi_trans(:)
                g2_tl(:) = g2_tl(:)+dudeta_trans(:)
              else
                g1(:) = g1(:)+dudxi_trans(:)
                g2(:) = g2(:)+dudeta_trans(:)
                g1_tl(:) = g1(:)
                g2_tl(:) = g2(:)
              endif
            endif

            !-------------------------------------------------

            ! [ B L ] matrix
            do nb = 1, nn

              jsize1 = ndof*(nb-1)+1
              jsize2 = ndof*(nb-1)+2
              jsize3 = ndof*(nb-1)+3
              jsize4 = ndof*(nb-1)+4
              jsize5 = ndof*(nb-1)+5
              jsize6 = ndof*(nb-1)+6

              aa1(1) = dudxi_rot(2, nb)  *g1(3)-dudxi_rot(3, nb)  *g1(2)
              aa1(2) = dudxi_rot(3, nb)  *g1(1)-dudxi_rot(1, nb)  *g1(3)
              aa1(3) = dudxi_rot(1, nb)  *g1(2)-dudxi_rot(2, nb)  *g1(1)

              aa2(1) = dudxi_rot(2, nb)  *g2(3)-dudxi_rot(3, nb)  *g2(2)
              aa2(2) = dudxi_rot(3, nb)  *g2(1)-dudxi_rot(1, nb)  *g2(3)
              aa2(3) = dudxi_rot(1, nb)  *g2(2)-dudxi_rot(2, nb)  *g2(1)

              aa3(1) = dudxi_rot(2, nb)  *g3(3)-dudxi_rot(3, nb)  *g3(2)
              aa3(2) = dudxi_rot(3, nb)  *g3(1)-dudxi_rot(1, nb)  *g3(3)
              aa3(3) = dudxi_rot(1, nb)  *g3(2)-dudxi_rot(2, nb)  *g3(1)

              bb1(1) = dudeta_rot(2, nb) *g1(3)-dudeta_rot(3, nb) *g1(2)
              bb1(2) = dudeta_rot(3, nb) *g1(1)-dudeta_rot(1, nb) *g1(3)
              bb1(3) = dudeta_rot(1, nb) *g1(2)-dudeta_rot(2, nb) *g1(1)

              bb2(1) = dudeta_rot(2, nb) *g2(3)-dudeta_rot(3, nb) *g2(2)
              bb2(2) = dudeta_rot(3, nb) *g2(1)-dudeta_rot(1, nb) *g2(3)
              bb2(3) = dudeta_rot(1, nb) *g2(2)-dudeta_rot(2, nb) *g2(1)

              bb3(1) = dudeta_rot(2, nb) *g3(3)-dudeta_rot(3, nb) *g3(2)
              bb3(2) = dudeta_rot(3, nb) *g3(1)-dudeta_rot(1, nb) *g3(3)
              bb3(3) = dudeta_rot(1, nb) *g3(2)-dudeta_rot(2, nb) *g3(1)

              cc1(1) = dudzeta_rot(2, nb)*g1(3)-dudzeta_rot(3, nb)*g1(2)
              cc1(2) = dudzeta_rot(3, nb)*g1(1)-dudzeta_rot(1, nb)*g1(3)
              cc1(3) = dudzeta_rot(1, nb)*g1(2)-dudzeta_rot(2, nb)*g1(1)

              cc2(1) = dudzeta_rot(2, nb)*g2(3)-dudzeta_rot(3, nb)*g2(2)
              cc2(2) = dudzeta_rot(3, nb)*g2(1)-dudzeta_rot(1, nb)*g2(3)
              cc2(3) = dudzeta_rot(1, nb)*g2(2)-dudzeta_rot(2, nb)*g2(1)

              if( use_director_tangent ) then
                dudxi_rot_deriv(1:3, 1:3) = shapederiv(nb, 1)*zeta_ly &
                  *a_over_2_v3_deriv(1:3, 1:3, nb)
                dudeta_rot_deriv(1:3, 1:3) = shapederiv(nb, 2)*zeta_ly &
                  *a_over_2_v3_deriv(1:3, 1:3, nb)
                dudzeta_rot_deriv(1:3, 1:3) = shapefunc(nb) &
                  *a_over_2_v3_deriv(1:3, 1:3, nb)
                dudxi_rot_second(1:3, 1:3, 1:3) = shapederiv(nb, 1)*zeta_ly &
                  *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
                dudeta_rot_second(1:3, 1:3, 1:3) = shapederiv(nb, 2)*zeta_ly &
                  *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
                dudzeta_rot_second(1:3, 1:3, 1:3) = shapefunc(nb) &
                  *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
              endif

              B_di(1, jsize1, ip, it, ly) = shapederiv(nb, 1)*g1_tl(1)
              B_di(2, jsize1, ip, it, ly) = shapederiv(nb, 2)*g2_tl(1)
              B_di(3, jsize1, ip, it, ly) = shapederiv(nb, 1)*g2_tl(1) &
                +shapederiv(nb, 2)*g1_tl(1)
              B_di(4, jsize1, ip, it, ly) = shapederiv(nb, 2)*g3(1)
              B_di(5, jsize1, ip, it, ly) = shapederiv(nb, 1)*g3(1)

              B_di(1, jsize2, ip, it, ly) = shapederiv(nb, 1)*g1_tl(2)
              B_di(2, jsize2, ip, it, ly) = shapederiv(nb, 2)*g2_tl(2)
              B_di(3, jsize2, ip, it, ly) = shapederiv(nb, 1)*g2_tl(2) &
                +shapederiv(nb, 2)*g1_tl(2)
              B_di(4, jsize2, ip, it, ly) = shapederiv(nb, 2)*g3(2)
              B_di(5, jsize2, ip, it, ly) = shapederiv(nb, 1)*g3(2)

              B_di(1, jsize3, ip, it, ly) = shapederiv(nb, 1)*g1_tl(3)
              B_di(2, jsize3, ip, it, ly) = shapederiv(nb, 2)*g2_tl(3)
              B_di(3, jsize3, ip, it, ly) = shapederiv(nb, 1)*g2_tl(3) &
                +shapederiv(nb, 2)*g1_tl(3)
              B_di(4, jsize3, ip, it, ly) = shapederiv(nb, 2)*g3(3)
              B_di(5, jsize3, ip, it, ly) = shapederiv(nb, 1)*g3(3)

              if( use_director_tangent ) then
                do m = 1, 3
                  jsize = ndof*(nb-1)+3+m
                  B_di(1, jsize, ip, it, ly) = dot_product(dudxi_rot_deriv(1:3, m), g1)
                  B_di(2, jsize, ip, it, ly) = dot_product(dudeta_rot_deriv(1:3, m), g2)
                  B_di(3, jsize, ip, it, ly) = dot_product(dudxi_rot_deriv(1:3, m), g2) &
                    +dot_product(dudeta_rot_deriv(1:3, m), g1)
                  B_di(4, jsize, ip, it, ly) = dot_product(dudeta_rot_deriv(1:3, m), g3) &
                    +dot_product(dudzeta_rot_deriv(1:3, m), g2)
                  B_di(5, jsize, ip, it, ly) = dot_product(dudxi_rot_deriv(1:3, m), g3) &
                    +dot_product(dudzeta_rot_deriv(1:3, m), g1)
                end do
                do n = 1, 3
                  do m = 1, 3
                    B2rot_di(1, m, n, nb, ip, it, ly) = dot_product(dudxi_rot_second(1:3, m, n), g1)
                    B2rot_di(2, m, n, nb, ip, it, ly) = dot_product(dudeta_rot_second(1:3, m, n), g2)
                    B2rot_di(3, m, n, nb, ip, it, ly) = dot_product(dudxi_rot_second(1:3, m, n), g2) &
                      +dot_product(dudeta_rot_second(1:3, m, n), g1)
                    B2rot_di(4, m, n, nb, ip, it, ly) = dot_product(dudeta_rot_second(1:3, m, n), g3) &
                      +dot_product(dudzeta_rot_second(1:3, m, n), g2)
                    B2rot_di(5, m, n, nb, ip, it, ly) = dot_product(dudxi_rot_second(1:3, m, n), g3) &
                      +dot_product(dudzeta_rot_second(1:3, m, n), g1)
                  end do
                end do
              else
                B_di(1, jsize4, ip, it, ly) = aa1(1)
                B_di(2, jsize4, ip, it, ly) = bb2(1)
                B_di(3, jsize4, ip, it, ly) = aa2(1)+bb1(1)
                B_di(4, jsize4, ip, it, ly) = bb3(1)+cc2(1)
                B_di(5, jsize4, ip, it, ly) = aa3(1)+cc1(1)

                B_di(1, jsize5, ip, it, ly) = aa1(2)
                B_di(2, jsize5, ip, it, ly) = bb2(2)
                B_di(3, jsize5, ip, it, ly) = aa2(2)+bb1(2)
                B_di(4, jsize5, ip, it, ly) = bb3(2)+cc2(2)
                B_di(5, jsize5, ip, it, ly) = aa3(2)+cc1(2)

                B_di(1, jsize6, ip, it, ly) = aa1(3)
                B_di(2, jsize6, ip, it, ly) = bb2(3)
                B_di(3, jsize6, ip, it, ly) = aa2(3)+bb1(3)
                B_di(4, jsize6, ip, it, ly) = bb3(3)+cc2(3)
                B_di(5, jsize6, ip, it, ly) = aa3(3)+cc1(3)
              endif

              ! First variations of g1/g2/g3 at tying points, used by the MITC4 shear geometric stiffness.
              BG1_di(1:3, jsize1, ip, it, ly) = (/ shapederiv(nb, 1), 0.0D0, 0.0D0 /)
              BG1_di(1:3, jsize2, ip, it, ly) = (/ 0.0D0, shapederiv(nb, 1), 0.0D0 /)
              BG1_di(1:3, jsize3, ip, it, ly) = (/ 0.0D0, 0.0D0, shapederiv(nb, 1) /)
              BG2_di(1:3, jsize1, ip, it, ly) = (/ shapederiv(nb, 2), 0.0D0, 0.0D0 /)
              BG2_di(1:3, jsize2, ip, it, ly) = (/ 0.0D0, shapederiv(nb, 2), 0.0D0 /)
              BG2_di(1:3, jsize3, ip, it, ly) = (/ 0.0D0, 0.0D0, shapederiv(nb, 2) /)
              if( use_director_tangent ) then
                do m = 1, 3
                  jsize = ndof*(nb-1)+3+m
                  BG1_di(1:3, jsize, ip, it, ly) = dudxi_rot_deriv(1:3, m)
                  BG2_di(1:3, jsize, ip, it, ly) = dudeta_rot_deriv(1:3, m)
                  BG3_di(1:3, jsize, ip, it, ly) = dudzeta_rot_deriv(1:3, m)
                end do
              else
                BG1_di(1:3, jsize4, ip, it, ly) = (/ 0.0D0, -dudxi_rot(3, nb), dudxi_rot(2, nb) /)
                BG1_di(1:3, jsize5, ip, it, ly) = (/ dudxi_rot(3, nb), 0.0D0, -dudxi_rot(1, nb) /)
                BG1_di(1:3, jsize6, ip, it, ly) = (/ -dudxi_rot(2, nb), dudxi_rot(1, nb), 0.0D0 /)
                BG2_di(1:3, jsize4, ip, it, ly) = (/ 0.0D0, -dudeta_rot(3, nb), dudeta_rot(2, nb) /)
                BG2_di(1:3, jsize5, ip, it, ly) = (/ dudeta_rot(3, nb), 0.0D0, -dudeta_rot(1, nb) /)
                BG2_di(1:3, jsize6, ip, it, ly) = (/ -dudeta_rot(2, nb), dudeta_rot(1, nb), 0.0D0 /)
                BG3_di(1:3, jsize4, ip, it, ly) = (/ 0.0D0, -dudzeta_rot(3, nb), dudzeta_rot(2, nb) /)
                BG3_di(1:3, jsize5, ip, it, ly) = (/ dudzeta_rot(3, nb), 0.0D0, -dudzeta_rot(1, nb) /)
                BG3_di(1:3, jsize6, ip, it, ly) = (/ -dudzeta_rot(2, nb), dudzeta_rot(1, nb), 0.0D0 /)
              endif

            end do

            !-------------------------------------------------

          end do

        end do

        !--------------------------------------------------------

        call fstr_shell_layer_quadrature_gauss( etype, gausses(1), n_layer, ly, &
          zeta_ly, w_ly, ierr_quad )
        if( ierr_quad /= 0 ) cycle

        !--------------------------------------------------------

        do lx = 1, NumOfQuadPoints(fetype)

          !--------------------------------------------------

          call getQuadPoint(fetype, lx, naturalcoord)

          xi_lx  = naturalcoord(1)
          eta_lx = naturalcoord(2)

          w_w_lx = getWeight(fetype, lx)

          call getShapeFunc(fetype, naturalcoord, shapefunc)

          call getShapeDeriv(fetype, naturalcoord, shapederiv)

          !--------------------------------------------------

          do i = 1, 3

            v1_i(i) = 0.0D0
            v2_i(i) = 0.0D0
            v3_i(i) = 0.0D0

            do na = 1, nn

              v1_i(i) = v1_i(i)+shapefunc(na)*v1(i, na)
              v2_i(i) = v2_i(i)+shapefunc(na)*v2(i, na)
              v3_i(i) = v3_i(i)+shapefunc(na)*v3(i, na)

            end do

          end do

          !--------------------------------------------------

          do na = 1, nn

            do i = 1, 3

              u_rot(i, na)                      &
                = shapefunc(na)                   &
                *( zeta_ly*a_over_2_v3(i, na) )

              dudxi_rot(i, na)                  &
                = shapederiv(na, 1)               &
                *( zeta_ly*a_over_2_v3(i, na) )
              dudeta_rot(i, na)                 &
                = shapederiv(na, 2)               &
                *( zeta_ly*a_over_2_v3(i, na) )
              dudzeta_rot(i, na)                &
                = shapefunc(na)                   &
                *( a_over_2_v3(i, na) )

            end do

          end do

          !--------------------------------------------------

          ! Covariant basis vector
          do i = 1, 3

            g1(i) = 0.0D0
            g2(i) = 0.0D0
            g3(i) = 0.0D0

            do na = 1, nn

              g1(i) = g1(i)+shapederiv(na, 1)  &
                *elem(i, na)       &
                +dudxi_rot(i, na)
              g2(i) = g2(i)+shapederiv(na, 2)  &
                *elem(i, na)       &
                +dudeta_rot(i, na)
              g3(i) = g3(i)+dudzeta_rot(i, na)

            end do

          end do

          !--------------------------------------------------

          if( use_tl_green .and. etype == fe_mitc4_shell ) then
            ! reference shell geometry for TL quadrature weight
            g1_weight(:) = 0.0D0
            g2_weight(:) = 0.0D0
            g3_weight(:) = 0.0D0
            do na = 1, nn
              g1_weight(:) = g1_weight(:)+shapederiv(na, 1) &
                *(ecoord(1:3, na)+zeta_ly*a_over_2_v3_ref(1:3, na))
              g2_weight(:) = g2_weight(:)+shapederiv(na, 2) &
                *(ecoord(1:3, na)+zeta_ly*a_over_2_v3_ref(1:3, na))
              g3_weight(:) = g3_weight(:)+shapefunc(na)*a_over_2_v3_ref(1:3, na)
            end do
            det_weight = g1_weight(1)*( g2_weight(2)*g3_weight(3)-g2_weight(3)*g3_weight(2) ) &
              +g1_weight(2)*( g2_weight(3)*g3_weight(1)-g2_weight(1)*g3_weight(3) ) &
              +g1_weight(3)*( g2_weight(1)*g3_weight(2)-g2_weight(2)*g3_weight(1) )
          endif

          !--------------------------------------------------

          if( use_tl_green .and. etype == fe_mitc4_shell ) then
            ! Green-Lagrange membrane and bending strain
            dudxi_trans(:) = 0.0D0
            dudeta_trans(:) = 0.0D0
            do na = 1, nn
              dudxi_trans(:) = dudxi_trans(:)+shapederiv(na, 1)*shell_disp(1:3, na)
              dudeta_trans(:) = dudeta_trans(:)+shapederiv(na, 2)*shell_disp(1:3, na)
            end do
            g1(:) = g1(:)+dudxi_trans(:)
            g2(:) = g2(:)+dudeta_trans(:)
          endif

          !--------------------------------------------------

          ! Jacobian
          det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
            +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
            +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

          if( .not. ( use_tl_green .and. etype == fe_mitc4_shell ) ) det_weight = det

          det_inv = 1.0D0/det

          !--------------------------------------------------

          ! Contravariant basis vector
          cg1(1) = det_inv                      &
            *( g2(2)*g3(3)-g2(3)*g3(2) )
          cg1(2) = det_inv                      &
            *( g2(3)*g3(1)-g2(1)*g3(3) )
          cg1(3) = det_inv                      &
            *( g2(1)*g3(2)-g2(2)*g3(1) )
          cg2(1) = det_inv                      &
            *( g3(2)*g1(3)-g3(3)*g1(2) )
          cg2(2) = det_inv                      &
            *( g3(3)*g1(1)-g3(1)*g1(3) )
          cg2(3) = det_inv                      &
            *( g3(1)*g1(2)-g3(2)*g1(1) )
          cg3(1) = det_inv                      &
            *( g1(2)*g2(3)-g1(3)*g2(2) )
          cg3(2) = det_inv                      &
            *( g1(3)*g2(1)-g1(1)*g2(3) )
          cg3(3) = det_inv                      &
            *( g1(1)*g2(2)-g1(2)*g2(1) )

          !--------------------------------------------------

          g3_abs = dsqrt( g3(1)*g3(1)   &
            +g3(2)*g3(2)   &
            +g3(3)*g3(3) )

          !--------------------------------------------------

          ! Orthonormal vectors

          e3_hat(1) = g3(1)/g3_abs
          e3_hat(2) = g3(2)/g3_abs
          e3_hat(3) = g3(3)/g3_abs

          e1_hat(1) = g2(2)*e3_hat(3) &
            -g2(3)*e3_hat(2)
          e1_hat(2) = g2(3)*e3_hat(1) &
            -g2(1)*e3_hat(3)
          e1_hat(3) = g2(1)*e3_hat(2) &
            -g2(2)*e3_hat(1)
          e1_hat_abs = dsqrt( e1_hat(1)*e1_hat(1)   &
            +e1_hat(2)*e1_hat(2)   &
            +e1_hat(3)*e1_hat(3) )
          e1_hat(1) = e1_hat(1)/e1_hat_abs
          e1_hat(2) = e1_hat(2)/e1_hat_abs
          e1_hat(3) = e1_hat(3)/e1_hat_abs

          e2_hat(1) = e3_hat(2)*e1_hat(3) &
            -e3_hat(3)*e1_hat(2)
          e2_hat(2) = e3_hat(3)*e1_hat(1) &
            -e3_hat(1)*e1_hat(3)
          e2_hat(3) = e3_hat(1)*e1_hat(2) &
            -e3_hat(2)*e1_hat(1)
          e2_hat_abs = dsqrt( e2_hat(1)*e2_hat(1)   &
            +e2_hat(2)*e2_hat(2)   &
            +e2_hat(3)*e2_hat(3) )
          e2_hat(1) = e2_hat(1)/e2_hat_abs
          e2_hat(2) = e2_hat(2)/e2_hat_abs
          e2_hat(3) = e2_hat(3)/e2_hat_abs

          !--------------------------------------------------

          e1_hat_mat(:) = e1_hat(:)
          e2_hat_mat(:) = e2_hat(:)
          e3_hat_mat(:) = e3_hat(:)
          cg1_mat(:) = cg1(:)
          cg2_mat(:) = cg2(:)
          cg3_mat(:) = cg3(:)
          if( use_tl_green .and. etype == fe_mitc4_shell ) then
            call shell_basis_from_covariant(g1_weight, g2_weight, g3_weight, &
              e1_hat_mat, e2_hat_mat, e3_hat_mat, cg1_mat, cg2_mat, cg3_mat, det_mat)
          endif

          !--------------------------------------------------

          if( present( element ) ) then
            ishell = fstr_shell_layer_gauss_index( element, lx, n_layer, ly )
          else
            ishell = 0
          endif

          if( ishell > 0 ) then
            call MatlMatrix_Shell                        &
              (element%shell_layer_gausses(ishell), Shell, D, &
              e1_hat_mat, e2_hat_mat, e3_hat_mat, cg1_mat, cg2_mat, cg3_mat, &
              alpha, n_layer)
          else
            call MatlMatrix_Shell                        &
              (gausses(lx), Shell, D,                 &
              e1_hat_mat, e2_hat_mat, e3_hat_mat, cg1_mat, cg2_mat, cg3_mat, &
              alpha, n_layer)
          endif
          !--------------------------------------------------

          g1_tl(:) = g1(:)
          g2_tl(:) = g2(:)
          if( use_tl_green ) then
            dudxi_trans(:) = 0.0D0
            dudeta_trans(:) = 0.0D0
            do na = 1, nn
              dudxi_trans(:) = dudxi_trans(:)+shapederiv(na, 1)*shell_disp(1:3, na)
              dudeta_trans(:) = dudeta_trans(:)+shapederiv(na, 2)*shell_disp(1:3, na)
            end do
            if( etype == fe_mitc9_shell ) then
              g1_tl(:) = g1_tl(:)+dudxi_trans(:)
              g2_tl(:) = g2_tl(:)+dudeta_trans(:)
            else if( etype /= fe_mitc4_shell ) then
              g1(:) = g1(:)+dudxi_trans(:)
              g2(:) = g2(:)+dudeta_trans(:)
              g1_tl(:) = g1(:)
              g2_tl(:) = g2(:)
            endif
          endif

          !--------------------------------------------------

          ! [ B L ] matrix
          B(:, :) = 0.0D0
          if( use_director_tangent ) B2rot(:, :, :, :) = 0.0D0
          do nb = 1, nn

            jsize1 = ndof*(nb-1)+1
            jsize2 = ndof*(nb-1)+2
            jsize3 = ndof*(nb-1)+3
            jsize4 = ndof*(nb-1)+4
            jsize5 = ndof*(nb-1)+5
            jsize6 = ndof*(nb-1)+6

            aa1(1) = dudxi_rot(2, nb)  *g1(3)-dudxi_rot(3, nb)  *g1(2)
            aa1(2) = dudxi_rot(3, nb)  *g1(1)-dudxi_rot(1, nb)  *g1(3)
            aa1(3) = dudxi_rot(1, nb)  *g1(2)-dudxi_rot(2, nb)  *g1(1)

            aa2(1) = dudxi_rot(2, nb)  *g2(3)-dudxi_rot(3, nb)  *g2(2)
            aa2(2) = dudxi_rot(3, nb)  *g2(1)-dudxi_rot(1, nb)  *g2(3)
            aa2(3) = dudxi_rot(1, nb)  *g2(2)-dudxi_rot(2, nb)  *g2(1)

            aa3(1) = dudxi_rot(2, nb)  *g3(3)-dudxi_rot(3, nb)  *g3(2)
            aa3(2) = dudxi_rot(3, nb)  *g3(1)-dudxi_rot(1, nb)  *g3(3)
            aa3(3) = dudxi_rot(1, nb)  *g3(2)-dudxi_rot(2, nb)  *g3(1)

            bb1(1) = dudeta_rot(2, nb) *g1(3)-dudeta_rot(3, nb) *g1(2)
            bb1(2) = dudeta_rot(3, nb) *g1(1)-dudeta_rot(1, nb) *g1(3)
            bb1(3) = dudeta_rot(1, nb) *g1(2)-dudeta_rot(2, nb) *g1(1)

            bb2(1) = dudeta_rot(2, nb) *g2(3)-dudeta_rot(3, nb) *g2(2)
            bb2(2) = dudeta_rot(3, nb) *g2(1)-dudeta_rot(1, nb) *g2(3)
            bb2(3) = dudeta_rot(1, nb) *g2(2)-dudeta_rot(2, nb) *g2(1)

            bb3(1) = dudeta_rot(2, nb) *g3(3)-dudeta_rot(3, nb) *g3(2)
            bb3(2) = dudeta_rot(3, nb) *g3(1)-dudeta_rot(1, nb) *g3(3)
            bb3(3) = dudeta_rot(1, nb) *g3(2)-dudeta_rot(2, nb) *g3(1)

            cc1(1) = dudzeta_rot(2, nb)*g1(3)-dudzeta_rot(3, nb)*g1(2)
            cc1(2) = dudzeta_rot(3, nb)*g1(1)-dudzeta_rot(1, nb)*g1(3)
            cc1(3) = dudzeta_rot(1, nb)*g1(2)-dudzeta_rot(2, nb)*g1(1)

            cc2(1) = dudzeta_rot(2, nb)*g2(3)-dudzeta_rot(3, nb)*g2(2)
            cc2(2) = dudzeta_rot(3, nb)*g2(1)-dudzeta_rot(1, nb)*g2(3)
            cc2(3) = dudzeta_rot(1, nb)*g2(2)-dudzeta_rot(2, nb)*g2(1)

            if( use_director_tangent ) then
              dudxi_rot_deriv(1:3, 1:3) = shapederiv(nb, 1)*zeta_ly &
                *a_over_2_v3_deriv(1:3, 1:3, nb)
              dudeta_rot_deriv(1:3, 1:3) = shapederiv(nb, 2)*zeta_ly &
                *a_over_2_v3_deriv(1:3, 1:3, nb)
              dudzeta_rot_deriv(1:3, 1:3) = shapefunc(nb) &
                *a_over_2_v3_deriv(1:3, 1:3, nb)
              dudxi_rot_second(1:3, 1:3, 1:3) = shapederiv(nb, 1)*zeta_ly &
                *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
              dudeta_rot_second(1:3, 1:3, 1:3) = shapederiv(nb, 2)*zeta_ly &
                *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
              dudzeta_rot_second(1:3, 1:3, 1:3) = shapefunc(nb) &
                *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
            endif

            B(1, jsize1) = shapederiv(nb, 1)*g1_tl(1)
            B(2, jsize1) = shapederiv(nb, 2)*g2_tl(1)
            B(3, jsize1) = shapederiv(nb, 1)*g2_tl(1) &
              +shapederiv(nb, 2)*g1_tl(1)
            B(4, jsize1) = shapederiv(nb, 2)*g3(1)
            B(5, jsize1) = shapederiv(nb, 1)*g3(1)

            B(1, jsize2) = shapederiv(nb, 1)*g1_tl(2)
            B(2, jsize2) = shapederiv(nb, 2)*g2_tl(2)
            B(3, jsize2) = shapederiv(nb, 1)*g2_tl(2) &
              +shapederiv(nb, 2)*g1_tl(2)
            B(4, jsize2) = shapederiv(nb, 2)*g3(2)
            B(5, jsize2) = shapederiv(nb, 1)*g3(2)

            B(1, jsize3) = shapederiv(nb, 1)*g1_tl(3)
            B(2, jsize3) = shapederiv(nb, 2)*g2_tl(3)
            B(3, jsize3) = shapederiv(nb, 1)*g2_tl(3) &
              +shapederiv(nb, 2)*g1_tl(3)
            B(4, jsize3) = shapederiv(nb, 2)*g3(3)
            B(5, jsize3) = shapederiv(nb, 1)*g3(3)

            if( use_director_tangent ) then
              do m = 1, 3
                jsize = ndof*(nb-1)+3+m
                B(1, jsize) = dot_product(dudxi_rot_deriv(1:3, m), g1)
                B(2, jsize) = dot_product(dudeta_rot_deriv(1:3, m), g2)
                B(3, jsize) = dot_product(dudxi_rot_deriv(1:3, m), g2) &
                  +dot_product(dudeta_rot_deriv(1:3, m), g1)
                B(4, jsize) = dot_product(dudeta_rot_deriv(1:3, m), g3) &
                  +dot_product(dudzeta_rot_deriv(1:3, m), g2)
                B(5, jsize) = dot_product(dudxi_rot_deriv(1:3, m), g3) &
                  +dot_product(dudzeta_rot_deriv(1:3, m), g1)
              end do
              do n = 1, 3
                do m = 1, 3
                  B2rot(1, m, n, nb) = dot_product(dudxi_rot_second(1:3, m, n), g1)
                  B2rot(2, m, n, nb) = dot_product(dudeta_rot_second(1:3, m, n), g2)
                  B2rot(3, m, n, nb) = dot_product(dudxi_rot_second(1:3, m, n), g2) &
                    +dot_product(dudeta_rot_second(1:3, m, n), g1)
                  B2rot(4, m, n, nb) = dot_product(dudeta_rot_second(1:3, m, n), g3) &
                    +dot_product(dudzeta_rot_second(1:3, m, n), g2)
                  B2rot(5, m, n, nb) = dot_product(dudxi_rot_second(1:3, m, n), g3) &
                    +dot_product(dudzeta_rot_second(1:3, m, n), g1)
                end do
              end do
            else
              B(1, jsize4) = aa1(1)
              B(2, jsize4) = bb2(1)
              B(3, jsize4) = aa2(1)+bb1(1)
              B(4, jsize4) = bb3(1)+cc2(1)
              B(5, jsize4) = aa3(1)+cc1(1)

              B(1, jsize5) = aa1(2)
              B(2, jsize5) = bb2(2)
              B(3, jsize5) = aa2(2)+bb1(2)
              B(4, jsize5) = bb3(2)+cc2(2)
              B(5, jsize5) = aa3(2)+cc1(2)

              B(1, jsize6) = aa1(3)
              B(2, jsize6) = bb2(3)
              B(3, jsize6) = aa2(3)+bb1(3)
              B(4, jsize6) = bb3(3)+cc2(3)
              B(5, jsize6) = aa3(3)+cc1(3)
            endif

          end do

          !--------------------------------------------------

          ! MITC4
          if( etype .EQ. fe_mitc4_shell ) then

            do jsize = 1, ndof*nn

              B(4, jsize) = 0.0D0
              B(5, jsize) = 0.0D0

              ! B_as(4, jsize)
              B(4, jsize)                                       &
                = 0.5D0*( 1.0D0-xi_lx  )*B_di(4, jsize, 4, 1, ly) &
                +0.5D0*( 1.0D0+xi_lx  )*B_di(4, jsize, 2, 1, ly)
              ! B_as(5, jsize)
              B(5, jsize)                                       &
                = 0.5D0*( 1.0D0-eta_lx )*B_di(5, jsize, 1, 1, ly) &
                +0.5D0*( 1.0D0+eta_lx )*B_di(5, jsize, 3, 1, ly)

            end do
            if( use_director_tangent ) then
              do nb = 1, nn
                do n = 1, 3
                  do m = 1, 3
                    B2rot(4, m, n, nb) &
                      = 0.5D0*( 1.0D0-xi_lx  )*B2rot_di(4, m, n, nb, 4, 1, ly) &
                      +0.5D0*( 1.0D0+xi_lx  )*B2rot_di(4, m, n, nb, 2, 1, ly)
                    B2rot(5, m, n, nb) &
                      = 0.5D0*( 1.0D0-eta_lx )*B2rot_di(5, m, n, nb, 1, 1, ly) &
                      +0.5D0*( 1.0D0+eta_lx )*B2rot_di(5, m, n, nb, 3, 1, ly)
                  end do
                end do
              end do
            endif

            ! MITC9
          else if( etype .EQ. fe_mitc9_shell ) then

            xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
            eeta_lx = eta_lx/dsqrt( 3.0D0/5.0D0 )

            do ip = 1, npoints_tying(1)

              h(ip, 1)                                     &
                = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )   &
                *( ( 0.5D0*eeta_di(ip, 1)*eeta_lx )        &
                *( 1.0D0+eeta_di(ip, 1)*eeta_lx )       &
                +( 1.0D0-eeta_di(ip, 1)*eeta_di(ip, 1) ) &
                *( 1.0D0-eeta_lx*eeta_lx ) )

            end do

            xxi_lx  = xi_lx /dsqrt( 3.0D0/5.0D0 )
            eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

            do ip = 1, npoints_tying(2)

              h(ip, 2)                                      &
                = ( ( 0.5D0*xxi_di(ip, 2) *xxi_lx  )          &
                *( 1.0D0+xxi_di(ip, 2) *xxi_lx  )         &
                +( 1.0D0-xxi_di(ip, 2) *xxi_di(ip, 2)  )   &
                *( 1.0D0-xxi_lx*xxi_lx ) )                &
                *( 0.5D0*( 1.0D0+eeta_di(ip, 2)*eeta_lx ) )

            end do

            xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
            eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

            do ip = 1, npoints_tying(3)

              h(ip, 3)                                      &
                = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )    &
                *( 0.5D0*( 1.0D0+eeta_di(ip, 1)*eeta_lx ) )

            end do

            do jsize = 1, ndof*nn

              B(1, jsize) = 0.0D0
              B(2, jsize) = 0.0D0
              B(3, jsize) = 0.0D0
              B(4, jsize) = 0.0D0
              B(5, jsize) = 0.0D0

              do ip = 1, npoints_tying(1)

                ! B_as(1, jsize)
                B(1, jsize)                                      &
                  = B(1, jsize)+h(ip, 1)*B_di(1, jsize, ip, 1, ly)
                ! B_as(5, jsize)
                B(5, jsize)                                      &
                  = B(5, jsize)+h(ip, 1)*B_di(5, jsize, ip, 1, ly)

              end do

              do ip = 1, npoints_tying(2)

                ! B_as(2, jsize)
                B(2, jsize)                                      &
                  = B(2, jsize)+h(ip, 2)*B_di(2, jsize, ip, 2, ly)
                ! B_as(4, jsize)
                B(4, jsize)                                      &
                  = B(4, jsize)+h(ip, 2)*B_di(4, jsize, ip, 2, ly)

              end do

              do ip = 1, npoints_tying(3)

                ! B_as(3, jsize)
                B(3, jsize)                                      &
                  = B(3, jsize)+h(ip, 3)*B_di(3, jsize, ip, 3, ly)

              end do

            end do

            ! MITC3
          else if( etype .EQ. fe_mitc3_shell ) then

            do jsize = 1, ndof*nn

              B(4, jsize) = 0.0D0
              B(5, jsize) = 0.0D0

              ! B_as(4, jsize)
              B(4, jsize)                                 &
                = ( 1.0D0-xi_lx  )*B_di(4, jsize, 2, 1, ly) &
                +xi_lx *B_di(5, jsize, 1, 1, ly)           &
                +xi_lx *( B_di(4, jsize, 3, 1, ly)         &
                -B_di(5, jsize, 3, 1, ly)  )

              ! B_as(5, jsize)
              B(5, jsize)                                 &
                = eta_lx*B_di(4, jsize, 2, 1, ly)           &
                +( 1.0D0-eta_lx )*B_di(5, jsize, 1, 1, ly) &
                -eta_lx*( B_di(4, jsize, 3, 1, ly)         &
                -B_di(5, jsize, 3, 1, ly)  )

            end do

          end if

          !--------------------------------------------------

          w_w_w_det = w_w_lx*w_ly*det_weight

          !--------------------------------------------------

          if( present( qf_stress ) .or. add_geo_stiff ) then
            if( ishell > 0 ) then
              Sv_force(1) = element%shell_layer_gausses(ishell)%stress(1)
              Sv_force(2) = element%shell_layer_gausses(ishell)%stress(2)
              Sv_force(3) = element%shell_layer_gausses(ishell)%stress(4)
              Sv_force(4) = element%shell_layer_gausses(ishell)%stress(5)
              Sv_force(5) = element%shell_layer_gausses(ishell)%stress(6)
            else
              Sv_force(1) = gausses(lx)%stress(1)
              Sv_force(2) = gausses(lx)%stress(2)
              Sv_force(3) = gausses(lx)%stress(4)
              Sv_force(4) = gausses(lx)%stress(5)
              Sv_force(5) = gausses(lx)%stress(6)
            endif
          endif

          if( present( qf_stress ) ) then
            qf_tmp(1:ndof*nn) = qf_tmp(1:ndof*nn) &
              +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
              *matmul( Sv_force(1:5), B(1:5, 1:ndof*nn) )
          endif

          !--------------------------------------------------

          DB(1:5, 1:ndof*nn) = matmul( D, B(1:5, 1:ndof*nn ) )

          if( flag == UPDATELAG .and. etype == fe_mitc4_shell .and. nn == 4 ) then
            if( ishell > 0 ) then
              stress_old_vec(1:6) = element%shell_layer_gausses(ishell)%stress_bak(1:6)
              trace_coeff = ShellPlaneStressTraceCoeff(element%shell_layer_gausses(ishell), n_layer)
            else
              stress_old_vec(1:6) = gausses(lx)%stress_bak(1:6)
              trace_coeff = ShellPlaneStressTraceCoeff(gausses(lx), n_layer)
            endif
            call ShellAddULObjectiveTangent( stress_old_vec, ndof*nn, B, DB, trace_coeff=trace_coeff )
          endif

          !--------------------------------------------------

          do jsize=1,ndof*nn 
            do isize=1,ndof*nn
              tmpstiff(isize, jsize)                               &
                = tmpstiff(isize, jsize)                             &
                +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight*dot_product( B(:, isize), DB(:, jsize) )
            end do
          end do

          !--------------------------------------------------

          ! [ B_{i} ] matrix
          B1(:, :) = 0.0D0
          B2(:, :) = 0.0D0
          B3(:, :) = 0.0D0
          Cv_w_second(:, :, :) = 0.0D0
          do nb = 1, nn

            jsize1 = ndof*(nb-1)+1
            jsize2 = ndof*(nb-1)+2
            jsize3 = ndof*(nb-1)+3
            jsize4 = ndof*(nb-1)+4
            jsize5 = ndof*(nb-1)+5
            jsize6 = ndof*(nb-1)+6

            B1(1, jsize1) =  shapederiv(nb, 1)
            B1(2, jsize1) =  0.0D0
            B1(3, jsize1) =  0.0D0
            B1(1, jsize2) =  0.0D0
            B1(2, jsize2) =  shapederiv(nb, 1)
            B1(3, jsize2) =  0.0D0
            B1(1, jsize3) =  0.0D0
            B1(2, jsize3) =  0.0D0
            B1(3, jsize3) =  shapederiv(nb, 1)
            B1(1, jsize4) =  0.0D0
            B1(2, jsize4) = -dudxi_rot(3, nb)
            B1(3, jsize4) =  dudxi_rot(2, nb)
            B1(1, jsize5) =  dudxi_rot(3, nb)
            B1(2, jsize5) =  0.0D0
            B1(3, jsize5) = -dudxi_rot(1, nb)
            B1(1, jsize6) = -dudxi_rot(2, nb)
            B1(2, jsize6) =  dudxi_rot(1, nb)
            B1(3, jsize6) =  0.0D0

            B2(1, jsize1) =  shapederiv(nb, 2)
            B2(2, jsize1) =  0.0D0
            B2(3, jsize1) =  0.0D0
            B2(1, jsize2) =  0.0D0
            B2(2, jsize2) =  shapederiv(nb, 2)
            B2(3, jsize2) =  0.0D0
            B2(1, jsize3) =  0.0D0
            B2(2, jsize3) =  0.0D0
            B2(3, jsize3) =  shapederiv(nb, 2)
            B2(1, jsize4) =  0.0D0
            B2(2, jsize4) = -dudeta_rot(3, nb)
            B2(3, jsize4) =  dudeta_rot(2, nb)
            B2(1, jsize5) =  dudeta_rot(3, nb)
            B2(2, jsize5) =  0.0D0
            B2(3, jsize5) = -dudeta_rot(1, nb)
            B2(1, jsize6) = -dudeta_rot(2, nb)
            B2(2, jsize6) =  dudeta_rot(1, nb)
            B2(3, jsize6) =  0.0D0

            B3(1, jsize1) =  0.0D0
            B3(2, jsize1) =  0.0D0
            B3(3, jsize1) =  0.0D0
            B3(1, jsize2) =  0.0D0
            B3(2, jsize2) =  0.0D0
            B3(3, jsize2) =  0.0D0
            B3(1, jsize3) =  0.0D0
            B3(2, jsize3) =  0.0D0
            B3(3, jsize3) =  0.0D0
            B3(1, jsize4) =  0.0D0
            B3(2, jsize4) = -dudzeta_rot(3, nb)
            B3(3, jsize4) =  dudzeta_rot(2, nb)
            B3(1, jsize5) =  dudzeta_rot(3, nb)
            B3(2, jsize5) =  0.0D0
            B3(3, jsize5) = -dudzeta_rot(1, nb)
            B3(1, jsize6) = -dudzeta_rot(2, nb)
            B3(2, jsize6) =  dudzeta_rot(1, nb)
            B3(3, jsize6) =  0.0D0

            if( use_director_tangent ) then
              dudxi_rot_deriv(1:3, 1:3) = shapederiv(nb, 1)*zeta_ly &
                *a_over_2_v3_deriv(1:3, 1:3, nb)
              dudeta_rot_deriv(1:3, 1:3) = shapederiv(nb, 2)*zeta_ly &
                *a_over_2_v3_deriv(1:3, 1:3, nb)
              dudzeta_rot_deriv(1:3, 1:3) = shapefunc(nb) &
                *a_over_2_v3_deriv(1:3, 1:3, nb)
              dudxi_rot_second(1:3, 1:3, 1:3) = shapederiv(nb, 1)*zeta_ly &
                *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
              dudeta_rot_second(1:3, 1:3, 1:3) = shapederiv(nb, 2)*zeta_ly &
                *a_over_2_v3_second(1:3, 1:3, 1:3, nb)
              dudzeta_rot_second(1:3, 1:3, 1:3) = shapefunc(nb) &
                *a_over_2_v3_second(1:3, 1:3, 1:3, nb)

              do m = 1, 3
                jsize = ndof*(nb-1)+3+m
                B1(1:3, jsize) = dudxi_rot_deriv(1:3, m)
                B2(1:3, jsize) = dudeta_rot_deriv(1:3, m)
                B3(1:3, jsize) = dudzeta_rot_deriv(1:3, m)
              end do
              do n = 1, 3
                do m = 1, 3
                  Cv12_2 = ( cg1_mat(1)*dudxi_rot_second(2, m, n) &
                    +cg2_mat(1)*dudeta_rot_second(2, m, n) &
                    +cg3_mat(1)*dudzeta_rot_second(2, m, n) ) &
                    -( cg1_mat(2)*dudxi_rot_second(1, m, n) &
                    +cg2_mat(2)*dudeta_rot_second(1, m, n) &
                    +cg3_mat(2)*dudzeta_rot_second(1, m, n) )
                  Cv13_2 = ( cg1_mat(1)*dudxi_rot_second(3, m, n) &
                    +cg2_mat(1)*dudeta_rot_second(3, m, n) &
                    +cg3_mat(1)*dudzeta_rot_second(3, m, n) ) &
                    -( cg1_mat(3)*dudxi_rot_second(1, m, n) &
                    +cg2_mat(3)*dudeta_rot_second(1, m, n) &
                    +cg3_mat(3)*dudzeta_rot_second(1, m, n) )
                  Cv21_2 = ( cg1_mat(2)*dudxi_rot_second(1, m, n) &
                    +cg2_mat(2)*dudeta_rot_second(1, m, n) &
                    +cg3_mat(2)*dudzeta_rot_second(1, m, n) ) &
                    -( cg1_mat(1)*dudxi_rot_second(2, m, n) &
                    +cg2_mat(1)*dudeta_rot_second(2, m, n) &
                    +cg3_mat(1)*dudzeta_rot_second(2, m, n) )
                  Cv23_2 = ( cg1_mat(2)*dudxi_rot_second(3, m, n) &
                    +cg2_mat(2)*dudeta_rot_second(3, m, n) &
                    +cg3_mat(2)*dudzeta_rot_second(3, m, n) ) &
                    -( cg1_mat(3)*dudxi_rot_second(2, m, n) &
                    +cg2_mat(3)*dudeta_rot_second(2, m, n) &
                    +cg3_mat(3)*dudzeta_rot_second(2, m, n) )
                  Cv31_2 = ( cg1_mat(3)*dudxi_rot_second(1, m, n) &
                    +cg2_mat(3)*dudeta_rot_second(1, m, n) &
                    +cg3_mat(3)*dudzeta_rot_second(1, m, n) ) &
                    -( cg1_mat(1)*dudxi_rot_second(3, m, n) &
                    +cg2_mat(1)*dudeta_rot_second(3, m, n) &
                    +cg3_mat(1)*dudzeta_rot_second(3, m, n) )
                  Cv32_2 = ( cg1_mat(3)*dudxi_rot_second(2, m, n) &
                    +cg2_mat(3)*dudeta_rot_second(2, m, n) &
                    +cg3_mat(3)*dudzeta_rot_second(2, m, n) ) &
                    -( cg1_mat(2)*dudxi_rot_second(3, m, n) &
                    +cg2_mat(2)*dudeta_rot_second(3, m, n) &
                    +cg3_mat(2)*dudzeta_rot_second(3, m, n) )
                  Cv_w_second(m, n, nb) = v1_i(1)*Cv12_2*v2_i(2) &
                    +v1_i(1)*Cv13_2*v2_i(3) &
                    +v1_i(2)*Cv21_2*v2_i(1) &
                    +v1_i(2)*Cv23_2*v2_i(3) &
                    +v1_i(3)*Cv31_2*v2_i(1) &
                    +v1_i(3)*Cv32_2*v2_i(2)
                end do
              end do
            endif

          end do

          !--------------------------------------------------

          if( flag == UPDATELAG ) then
            do jsize = 1, ndof*nn
              vol_deriv(jsize) = dot_product(cg1(1:3), B1(1:3, jsize)) &
                +dot_product(cg2(1:3), B2(1:3, jsize)) &
                +dot_product(cg3(1:3), B3(1:3, jsize))
              qf_stress_integrand(jsize) = dot_product(Sv_force(1:5), B(1:5, jsize))
            end do
          endif

          !--------------------------------------------------

          if( add_geo_stiff ) then
            if( etype == fe_mitc4_shell .and. ( use_tl_green .or. flag == UPDATELAG ) ) then
              ! Match the MITC4 residual integration points.
              do jsize=1,ndof*nn
                do isize=1,ndof*nn
                  geo_term = Sv_force(1)*dot_product(B1(1:3, isize), B1(1:3, jsize)) &
                    +Sv_force(2)*dot_product(B2(1:3, isize), B2(1:3, jsize)) &
                    +Sv_force(3)*(dot_product(B1(1:3, isize), B2(1:3, jsize)) &
                    +dot_product(B2(1:3, isize), B1(1:3, jsize))) &
                    +Sv_force(4)*(0.5D0*(1.0D0-xi_lx) &
                    *(dot_product(BG2_di(1:3, isize, 4, 1, ly), BG3_di(1:3, jsize, 4, 1, ly)) &
                    +dot_product(BG3_di(1:3, isize, 4, 1, ly), BG2_di(1:3, jsize, 4, 1, ly))) &
                    +0.5D0*(1.0D0+xi_lx) &
                    *(dot_product(BG2_di(1:3, isize, 2, 1, ly), BG3_di(1:3, jsize, 2, 1, ly)) &
                    +dot_product(BG3_di(1:3, isize, 2, 1, ly), BG2_di(1:3, jsize, 2, 1, ly)))) &
                    +Sv_force(5)*(0.5D0*(1.0D0-eta_lx) &
                    *(dot_product(BG3_di(1:3, isize, 1, 1, ly), BG1_di(1:3, jsize, 1, 1, ly)) &
                    +dot_product(BG1_di(1:3, isize, 1, 1, ly), BG3_di(1:3, jsize, 1, 1, ly))) &
                    +0.5D0*(1.0D0+eta_lx) &
                    *(dot_product(BG3_di(1:3, isize, 3, 1, ly), BG1_di(1:3, jsize, 3, 1, ly)) &
                    +dot_product(BG1_di(1:3, isize, 3, 1, ly), BG3_di(1:3, jsize, 3, 1, ly))))
                  tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                    +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                    *geo_term
                end do
              end do
            else
              S_global(:, :) = 0.0D0
              S_global(:, :) = S_global(:, :) &
                +Sv_force(1)*outer_product3(e1_hat_mat, e1_hat_mat) &
                +Sv_force(2)*outer_product3(e2_hat_mat, e2_hat_mat) &
                +Sv_force(3)*(outer_product3(e1_hat_mat, e2_hat_mat) &
                +outer_product3(e2_hat_mat, e1_hat_mat)) &
                +Sv_force(4)*(outer_product3(e2_hat_mat, e3_hat_mat) &
                +outer_product3(e3_hat_mat, e2_hat_mat)) &
                +Sv_force(5)*(outer_product3(e1_hat_mat, e3_hat_mat) &
                +outer_product3(e3_hat_mat, e1_hat_mat))

              BN(1:9, 1:ndof*nn) = 0.0D0
              do jsize = 1, ndof*nn
                BN(1, jsize) = B1(1, jsize)*cg1_mat(1) + B2(1, jsize)*cg2_mat(1) &
                  +B3(1, jsize)*cg3_mat(1)
                BN(2, jsize) = B1(2, jsize)*cg1_mat(1) + B2(2, jsize)*cg2_mat(1) &
                  +B3(2, jsize)*cg3_mat(1)
                BN(3, jsize) = B1(3, jsize)*cg1_mat(1) + B2(3, jsize)*cg2_mat(1) &
                  +B3(3, jsize)*cg3_mat(1)
                BN(4, jsize) = B1(1, jsize)*cg1_mat(2) + B2(1, jsize)*cg2_mat(2) &
                  +B3(1, jsize)*cg3_mat(2)
                BN(5, jsize) = B1(2, jsize)*cg1_mat(2) + B2(2, jsize)*cg2_mat(2) &
                  +B3(2, jsize)*cg3_mat(2)
                BN(6, jsize) = B1(3, jsize)*cg1_mat(2) + B2(3, jsize)*cg2_mat(2) &
                  +B3(3, jsize)*cg3_mat(2)
                BN(7, jsize) = B1(1, jsize)*cg1_mat(3) + B2(1, jsize)*cg2_mat(3) &
                  +B3(1, jsize)*cg3_mat(3)
                BN(8, jsize) = B1(2, jsize)*cg1_mat(3) + B2(2, jsize)*cg2_mat(3) &
                  +B3(2, jsize)*cg3_mat(3)
                BN(9, jsize) = B1(3, jsize)*cg1_mat(3) + B2(3, jsize)*cg2_mat(3) &
                  +B3(3, jsize)*cg3_mat(3)
              end do

              Smat(:, :) = 0.0D0
              do j = 1, 3
                Smat(j  , j  ) = S_global(1, 1)
                Smat(j  , j+3) = S_global(1, 2)
                Smat(j  , j+6) = S_global(1, 3)
                Smat(j+3, j  ) = S_global(2, 1)
                Smat(j+3, j+3) = S_global(2, 2)
                Smat(j+3, j+6) = S_global(2, 3)
                Smat(j+6, j  ) = S_global(3, 1)
                Smat(j+6, j+3) = S_global(3, 2)
                Smat(j+6, j+6) = S_global(3, 3)
              end do

              SBN(1:9, 1:ndof*nn) = matmul(Smat(1:9, 1:9), BN(1:9, 1:ndof*nn))
              do jsize=1,ndof*nn
                do isize=1,ndof*nn
                  tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                    +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                    *dot_product(BN(:, isize), SBN(:, jsize))
                end do
              end do
            endif
            if( flag == UPDATELAG ) then
              do jsize=1,ndof*nn
                do isize=1,ndof*nn
                  tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                    +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                    *qf_stress_integrand(isize)*vol_deriv(jsize)
                end do
              end do
            endif
            if( use_director_tangent .and. use_tl_green ) then
              ! Director stress stiffness in the director/drilling split.
              do nb = 1, nn
                drill_axis(1:3) = a_over_2_v3(1:3, nb)
                axis_norm = sqrt(sum(drill_axis(1:3)*drill_axis(1:3)))
                if( axis_norm > 0.0D0 ) then
                  drill_axis(1:3) = drill_axis(1:3)/axis_norm
                else
                  drill_axis(1:3) = (/ 0.0D0, 0.0D0, 1.0D0 /)
                endif
                rot_projector(:, :) = 0.0D0
                do i = 1, 3
                  rot_projector(i, i) = 1.0D0
                end do
                do i = 1, 3
                  do j = 1, 3
                    rot_projector(i, j) = rot_projector(i, j)-drill_axis(i)*drill_axis(j)
                  end do
                end do
                do n = 1, 3
                  do m = 1, 3
                    hrot(m, n) = dot_product(Sv_force(1:5), B2rot(1:5, m, n, nb))
                  end do
                end do
                do n = 1, 3
                  jsize = ndof*(nb-1)+3+n
                  do m = 1, 3
                    isize = ndof*(nb-1)+3+m
                    hess_coeff = 0.0D0
                    do j = 1, 3
                      do i = 1, 3
                        hess_coeff = hess_coeff &
                          +rot_projector(i, m)*hrot(i, j)*rot_projector(j, n)
                      end do
                    end do
                    drill_coeff = 0.0D0
                    do i = 1, 3
                      drill_coeff = drill_coeff+drill_axis(i)*hrot(i, n)
                    end do
                    do j = 1, 3
                      do i = 1, 3
                        drill_coeff = drill_coeff &
                          +rot_projector(i, n)*hrot(i, j)*drill_axis(j)
                      end do
                    end do
                    tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                      +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                      *(hess_coeff+drill_axis(m)*drill_coeff)
                  end do
                end do
              end do
            endif
            if( use_director_tangent .and. flag == UPDATELAG ) then
              ! Current-frame counterpart of the director/drilling split.
              do nb = 1, nn
                drill_axis(1:3) = a_over_2_v3(1:3, nb)
                axis_norm = sqrt(sum(drill_axis(1:3)*drill_axis(1:3)))
                if( axis_norm > 0.0D0 ) then
                  drill_axis(1:3) = drill_axis(1:3)/axis_norm
                else
                  drill_axis(1:3) = (/ 0.0D0, 0.0D0, 1.0D0 /)
                endif
                rot_projector(:, :) = 0.0D0
                do i = 1, 3
                  rot_projector(i, i) = 1.0D0
                end do
                do i = 1, 3
                  do j = 1, 3
                    rot_projector(i, j) = rot_projector(i, j)-drill_axis(i)*drill_axis(j)
                  end do
                end do
                do n = 1, 3
                  do m = 1, 3
                    hrot(m, n) = dot_product(Sv_force(1:5), B2rot(1:5, m, n, nb))
                  end do
                end do
                do n = 1, 3
                  jsize = ndof*(nb-1)+3+n
                  do m = 1, 3
                    isize = ndof*(nb-1)+3+m
                    hess_coeff = 0.0D0
                    do j = 1, 3
                      do i = 1, 3
                        hess_coeff = hess_coeff &
                          +rot_projector(i, m)*hrot(i, j)*rot_projector(j, n)
                      end do
                    end do
                    drill_coeff = 0.0D0
                    do i = 1, 3
                      drill_coeff = drill_coeff+drill_axis(i)*hrot(i, n)
                    end do
                    do j = 1, 3
                      do i = 1, 3
                        drill_coeff = drill_coeff &
                          +rot_projector(i, n)*hrot(i, j)*drill_axis(j)
                      end do
                    end do
                    tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                      +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                      *(hess_coeff+drill_axis(m)*drill_coeff)
                  end do
                end do
              end do
            endif
          endif

          !--------------------------------------------------

          ! { C_{ij} } vector
          do jsize = 1, ndof*nn

            Cv12(jsize) = ( cg1_mat(1)*B1(2, jsize)   &
              +cg2_mat(1)*B2(2, jsize)   &
              +cg3_mat(1)*B3(2, jsize) ) &
              -( cg1_mat(2)*B1(1, jsize)   &
              +cg2_mat(2)*B2(1, jsize)   &
              +cg3_mat(2)*B3(1, jsize) )
            Cv13(jsize) = ( cg1_mat(1)*B1(3, jsize)   &
              +cg2_mat(1)*B2(3, jsize)   &
              +cg3_mat(1)*B3(3, jsize) ) &
              -( cg1_mat(3)*B1(1, jsize)   &
              +cg2_mat(3)*B2(1, jsize)   &
              +cg3_mat(3)*B3(1, jsize) )
            Cv21(jsize) = ( cg1_mat(2)*B1(1, jsize)   &
              +cg2_mat(2)*B2(1, jsize)   &
              +cg3_mat(2)*B3(1, jsize) ) &
              -( cg1_mat(1)*B1(2, jsize)   &
              +cg2_mat(1)*B2(2, jsize)   &
              +cg3_mat(1)*B3(2, jsize) )
            Cv23(jsize) = ( cg1_mat(2)*B1(3, jsize)   &
              +cg2_mat(2)*B2(3, jsize)   &
              +cg3_mat(2)*B3(3, jsize) ) &
              -( cg1_mat(3)*B1(2, jsize)   &
              +cg2_mat(3)*B2(2, jsize)   &
              +cg3_mat(3)*B3(2, jsize) )
            Cv31(jsize) = ( cg1_mat(3)*B1(1, jsize)   &
              +cg2_mat(3)*B2(1, jsize)   &
              +cg3_mat(3)*B3(1, jsize) ) &
              -( cg1_mat(1)*B1(3, jsize)   &
              +cg2_mat(1)*B2(3, jsize)   &
              +cg3_mat(1)*B3(3, jsize) )
            Cv32(jsize) = ( cg1_mat(3)*B1(2, jsize)   &
              +cg2_mat(3)*B2(2, jsize)   &
              +cg3_mat(3)*B3(2, jsize) ) &
              -( cg1_mat(2)*B1(3, jsize)   &
              +cg2_mat(2)*B2(3, jsize)   &
              +cg3_mat(2)*B3(3, jsize) )

          end do

          !--------------------------------------------------

          ! { Cw } vector
          do nb = 1, nn

            do j = 1, ndof

              jsize = ndof*(nb-1)+j

              Cv_w(jsize)                       &
                = v1_i(1)*Cv12(jsize)*v2_i(2) &
                +v1_i(1)*Cv13(jsize)*v2_i(3) &
                +v1_i(2)*Cv21(jsize)*v2_i(1) &
                +v1_i(2)*Cv23(jsize)*v2_i(3) &
                +v1_i(3)*Cv31(jsize)*v2_i(1) &
                +v1_i(3)*Cv32(jsize)*v2_i(2)

            end do

          end do

          ! { Ctheta } vector
          do nb = 1, nn

            jsize1 = ndof*(nb-1)+1
            jsize2 = ndof*(nb-1)+2
            jsize3 = ndof*(nb-1)+3
            jsize4 = ndof*(nb-1)+4
            jsize5 = ndof*(nb-1)+5
            jsize6 = ndof*(nb-1)+6

            Cv_theta(jsize1) = 0.0D0
            Cv_theta(jsize2) = 0.0D0
            Cv_theta(jsize3) = 0.0D0
            Cv_theta(jsize4) = v3_i(1)*shapefunc(nb)
            Cv_theta(jsize5) = v3_i(2)*shapefunc(nb)
            Cv_theta(jsize6) = v3_i(3)*shapefunc(nb)

          end do

          ! { C } vector
          ! drilling stabilization for compatibility rotation
          do jsize = 1, ndof*nn

            Cv(jsize) = Cv_theta(jsize)-0.5D0*Cv_w(jsize)

          end do

          !--------------------------------------------------

          Cv_disp = dot_product( Cv(1:ndof*nn), qf_disp(1:ndof*nn) )
          if( finite_rotation_director .and. present( nddrill ) ) then
            Cv_disp = -0.5D0*dot_product( Cv_w(1:ndof*nn), qf_disp(1:ndof*nn) )
            do nb = 1, nn
              Cv_disp = Cv_disp + shapefunc(nb)*nddrill(nb)
            end do
          endif

          if( present( qf_stress ) ) then
            qf_tmp(1:ndof*nn) = qf_tmp(1:ndof*nn) &
              +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
              *alpha*Cv(1:ndof*nn)*Cv_disp
          endif

          !--------------------------------------------------

          ! [ K L ] matrix
          do jsize = 1, ndof*nn
            do isize = 1, ndof*nn

              tmpstiff(isize, jsize)                &
                = tmpstiff(isize, jsize)              &
                +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight*alpha*Cv(isize)*Cv(jsize)

            end do
          end do

          if( use_director_tangent ) then
            do nb = 1, nn
              do n = 1, 3
                jsize = ndof*(nb-1)+3+n
                Cv_deriv_disp = 0.0D0
                do m = 1, 3
                  isize = ndof*(nb-1)+3+m
                  Cv_deriv = -0.5D0*Cv_w_second(m, n, nb)
                  Cv_deriv_disp = Cv_deriv_disp + Cv_deriv*qf_disp(isize)
                  tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                    +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                    *alpha*Cv_deriv*Cv_disp
                end do
                do isize = 1, ndof*nn
                  tmpstiff(isize, jsize) = tmpstiff(isize, jsize) &
                    +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight &
                    *alpha*Cv(isize)*Cv_deriv_disp
                end do
              end do
            end do
          endif


          !--------------------------------------------------

        end do

        !--------------------------------------------------------

      end do

      !--------------------------------------------------------------

      stiff(1:nn*ndof, 1:nn*ndof) = tmpstiff(1:nn*ndof, 1:nn*ndof)

      !--------------------------------------------------------------------
    end do     !< LAMINATED SHELL ANALYSIS

    !   write(*,"(24E25.16)") stiff

    !******************** Shell-Solid mixed analysis ********************
    !		mixglaf = 0 < natural shell (6 dof)
    !		mixglaf = 1 < mixed 361 (3*2 dof)*8 nod
    !		mixglaf = 2 < mixed 351 (3*2 dof)*6 nod
    if( mixflag == 1 )then

      !       write(*,*) 'convert for shell-solid mixed analysis'
      sstable(1) = 1
      sstable(2) = 2
      sstable(3) = 3
      sstable(4) = 7
      sstable(5) = 8
      sstable(6) = 9
      sstable(7) = 13
      sstable(8) = 14
      sstable(9) = 15
      sstable(10)= 19
      sstable(11)= 20
      sstable(12)= 21
      sstable(13)= 4
      sstable(14)= 5
      sstable(15)= 6
      sstable(16)= 10
      sstable(17)= 11
      sstable(18)= 12
      sstable(19)= 16
      sstable(20)= 17
      sstable(21)= 18
      sstable(22)= 22
      sstable(23)= 23
      sstable(24)= 24

      tmpstiff(1:nn*ndof, 1:nn*ndof) = stiff(1:nn*ndof, 1:nn*ndof)

      do i = 1, nn*ndof
        do j = 1, nn*ndof
          stiff(i,j) = tmpstiff(sstable(i),sstable(j))
        enddo
      enddo

      if( present( qf_stress ) ) then
        qf_mix(1:nn*ndof) = qf_tmp(1:nn*ndof)
        do i = 1, nn*ndof
          qf_stress(i) = qf_mix(sstable(i))
        enddo
      endif

    elseif( mixflag == 2 )then
      !		write(*,*) 'convert for shell-solid mixed analysis 351'
      sstable(1) = 1
      sstable(2) = 2
      sstable(3) = 3
      sstable(4) = 7
      sstable(5) = 8
      sstable(6) = 9
      sstable(7) = 13
      sstable(8) = 14
      sstable(9) = 15
      sstable(10)= 4
      sstable(11)= 5
      sstable(12)= 6
      sstable(13)= 10
      sstable(14)= 11
      sstable(15)= 12
      sstable(16)= 16
      sstable(17)= 17
      sstable(18)= 18

      tmpstiff(1:nn*ndof, 1:nn*ndof) = stiff(1:nn*ndof, 1:nn*ndof)

      do i = 1, nn*ndof
        do j = 1, nn*ndof
          stiff(i,j) = tmpstiff(sstable(i),sstable(j))
        enddo
      enddo

      if( present( qf_stress ) ) then
        qf_mix(1:nn*ndof) = qf_tmp(1:nn*ndof)
        do i = 1, nn*ndof
          qf_stress(i) = qf_mix(sstable(i))
        enddo
      endif

    else
      if( present( qf_stress ) ) qf_stress(1:nn*ndof) = qf_tmp(1:nn*ndof)

    endif

    return

    !####################################################################
  end subroutine STF_Shell_MITC
  !####################################################################


  !####################################################################
  subroutine ShellDirectorIncrement( theta, director_ref, director_inc )
    !####################################################################

    real(kind = kreal), intent(in)  :: theta(3)
    real(kind = kreal), intent(in)  :: director_ref(3)
    real(kind = kreal), intent(out) :: director_inc(3)

    real(kind = kreal) :: theta_norm, theta_norm2
    real(kind = kreal) :: cross1(3), cross2(3)
    real(kind = kreal) :: sin_over_theta, one_minus_cos_over_theta2

    cross1(1) = theta(2)*director_ref(3) - theta(3)*director_ref(2)
    cross1(2) = theta(3)*director_ref(1) - theta(1)*director_ref(3)
    cross1(3) = theta(1)*director_ref(2) - theta(2)*director_ref(1)

    cross2(1) = theta(2)*cross1(3) - theta(3)*cross1(2)
    cross2(2) = theta(3)*cross1(1) - theta(1)*cross1(3)
    cross2(3) = theta(1)*cross1(2) - theta(2)*cross1(1)

    theta_norm2 = dot_product( theta(1:3), theta(1:3) )
    theta_norm = dsqrt( theta_norm2 )

    if( theta_norm < 1.0d-12 ) then
      director_inc(1:3) = cross1(1:3) + 0.5d0*cross2(1:3)
    else
      sin_over_theta = dsin( theta_norm )/theta_norm
      one_minus_cos_over_theta2 = ( 1.0d0-dcos( theta_norm ) )/theta_norm2
      director_inc(1:3) = sin_over_theta*cross1(1:3) &
        +one_minus_cos_over_theta2*cross2(1:3)
    endif

    return

    !####################################################################
  end subroutine ShellDirectorIncrement
  !####################################################################

  !####################################################################
  pure subroutine ShellRotationVectorToMatrix( theta, rotmat )
    !####################################################################

    real(kind = kreal), intent(in)  :: theta(3)
    real(kind = kreal), intent(out) :: rotmat(3, 3)

    real(kind = kreal) :: theta_norm, theta_norm2
    real(kind = kreal) :: sin_over_theta, one_minus_cos_over_theta2
    real(kind = kreal) :: skew(3, 3), skew2(3, 3)
    integer :: i, j

    theta_norm2 = dot_product( theta(1:3), theta(1:3) )
    theta_norm = dsqrt( theta_norm2 )

    skew(:, :) = 0.0D0
    skew(1, 2) = -theta(3)
    skew(1, 3) =  theta(2)
    skew(2, 1) =  theta(3)
    skew(2, 3) = -theta(1)
    skew(3, 1) = -theta(2)
    skew(3, 2) =  theta(1)

    skew2 = matmul( skew, skew )

    if( theta_norm < 1.0D-12 ) then
      sin_over_theta = 1.0D0 - theta_norm2/6.0D0 + theta_norm2*theta_norm2/120.0D0
      one_minus_cos_over_theta2 = 0.5D0 - theta_norm2/24.0D0 + theta_norm2*theta_norm2/720.0D0
    else
      sin_over_theta = dsin( theta_norm )/theta_norm
      one_minus_cos_over_theta2 = ( 1.0D0-dcos( theta_norm ) )/theta_norm2
    endif

    rotmat(:, :) = sin_over_theta*skew(:, :) + one_minus_cos_over_theta2*skew2(:, :)
    do i = 1, 3
      rotmat(i, i) = rotmat(i, i) + 1.0D0
    end do

    return

    !####################################################################
  end subroutine ShellRotationVectorToMatrix
  !####################################################################


  !####################################################################
  pure subroutine ShellRotationMatrixToVector( rotmat, theta )
    !####################################################################

    real(kind = kreal), intent(in)  :: rotmat(3, 3)
    real(kind = kreal), intent(out) :: theta(3)

    real(kind = kreal) :: pi, trace_r, cos_angle, angle, sin_angle
    real(kind = kreal) :: axis(3), axis_abs

    pi = 4.0D0*datan( 1.0D0 )
    trace_r = rotmat(1, 1) + rotmat(2, 2) + rotmat(3, 3)
    cos_angle = 0.5D0*( trace_r - 1.0D0 )
    cos_angle = max( -1.0D0, min( 1.0D0, cos_angle ) )
    angle = dacos( cos_angle )

    theta(1) = rotmat(3, 2) - rotmat(2, 3)
    theta(2) = rotmat(1, 3) - rotmat(3, 1)
    theta(3) = rotmat(2, 1) - rotmat(1, 2)

    if( angle < 1.0D-12 ) then
      theta(1:3) = 0.5D0*theta(1:3)
    else if( pi-angle < 1.0D-8 ) then
      axis(1) = dsqrt( max( 0.0D0, 0.5D0*( rotmat(1, 1)+1.0D0 ) ) )
      axis(2) = dsqrt( max( 0.0D0, 0.5D0*( rotmat(2, 2)+1.0D0 ) ) )
      axis(3) = dsqrt( max( 0.0D0, 0.5D0*( rotmat(3, 3)+1.0D0 ) ) )

      if( rotmat(2, 1)+rotmat(1, 2) < 0.0D0 ) axis(2) = -axis(2)
      if( rotmat(3, 1)+rotmat(1, 3) < 0.0D0 ) axis(3) = -axis(3)
      axis_abs = dsqrt( dot_product( axis(1:3), axis(1:3) ) )
      if( axis_abs > 1.0D-12 ) then
        theta(1:3) = angle*axis(1:3)/axis_abs
      else
        theta(1:3) = 0.0D0
      endif
    else
      sin_angle = dsin( angle )
      theta(1:3) = angle*theta(1:3)/( 2.0D0*sin_angle )
    endif

    return

    !####################################################################
  end subroutine ShellRotationMatrixToVector
  !####################################################################


  !####################################################################
  pure subroutine ShellComposeRotationVector( theta_old, theta_inc, theta_new )
    !####################################################################

    real(kind = kreal), intent(in)  :: theta_old(3)
    real(kind = kreal), intent(in)  :: theta_inc(3)
    real(kind = kreal), intent(out) :: theta_new(3)

    real(kind = kreal) :: rot_old(3, 3), rot_inc(3, 3), rot_new(3, 3)

    call ShellRotationVectorToMatrix( theta_old, rot_old )
    call ShellRotationVectorToMatrix( theta_inc, rot_inc )
    rot_new = matmul( rot_inc, rot_old )
    call ShellRotationMatrixToVector( rot_new, theta_new )

    return

    !####################################################################
  end subroutine ShellComposeRotationVector
  !####################################################################


  !####################################################################
  pure subroutine ShellRelativeRotationVector( theta_old, theta_target, theta_inc )
    !####################################################################

    real(kind = kreal), intent(in)  :: theta_old(3)
    real(kind = kreal), intent(in)  :: theta_target(3)
    real(kind = kreal), intent(out) :: theta_inc(3)

    real(kind = kreal) :: rot_old(3, 3), rot_target(3, 3), rot_inc(3, 3)

    call ShellRotationVectorToMatrix( theta_old, rot_old )
    call ShellRotationVectorToMatrix( theta_target, rot_target )
    rot_inc = matmul( rot_target, transpose( rot_old ) )
    call ShellRotationMatrixToVector( rot_inc, theta_inc )

    return

    !####################################################################
  end subroutine ShellRelativeRotationVector
  !####################################################################


  !####################################################################
  pure subroutine ShellOrthonormalizeTriad( triad_in, triad_out )
    !####################################################################

    real(kind = kreal), intent(in)  :: triad_in(3, 3)
    real(kind = kreal), intent(out) :: triad_out(3, 3)

    real(kind = kreal) :: e1(3), e2(3), e3(3), normv

    e1(1:3) = triad_in(1:3, 1)
    e3(1:3) = triad_in(1:3, 3)

    normv = dsqrt( dot_product( e3(1:3), e3(1:3) ) )
    if( normv < 1.0D-14 ) then
      e3(1:3) = (/ 0.0D0, 0.0D0, 1.0D0 /)
    else
      e3(1:3) = e3(1:3)/normv
    endif

    e1(1:3) = e1(1:3) - dot_product( e1(1:3), e3(1:3) )*e3(1:3)
    normv = dsqrt( dot_product( e1(1:3), e1(1:3) ) )
    if( normv < 1.0D-14 ) then
      if( dabs(e3(1)) < 0.9D0 ) then
        e1(1:3) = (/ 1.0D0, 0.0D0, 0.0D0 /)
      else
        e1(1:3) = (/ 0.0D0, 1.0D0, 0.0D0 /)
      endif
      e1(1:3) = e1(1:3) - dot_product( e1(1:3), e3(1:3) )*e3(1:3)
      normv = dsqrt( dot_product( e1(1:3), e1(1:3) ) )
    endif
    e1(1:3) = e1(1:3)/normv

    e2(1) = e3(2)*e1(3)-e3(3)*e1(2)
    e2(2) = e3(3)*e1(1)-e3(1)*e1(3)
    e2(3) = e3(1)*e1(2)-e3(2)*e1(1)

    triad_out(1:3, 1) = e1(1:3)
    triad_out(1:3, 2) = e2(1:3)
    triad_out(1:3, 3) = e3(1:3)

    return

    !####################################################################
  end subroutine ShellOrthonormalizeTriad
  !####################################################################


  !####################################################################
  pure subroutine ShellUpdateTriadWithIncrement( triad_old, drill_old, theta_inc, triad_new, drill_new )
    !####################################################################

    real(kind = kreal), intent(in)  :: triad_old(3, 3)
    real(kind = kreal), intent(in)  :: drill_old
    real(kind = kreal), intent(in)  :: theta_inc(3)
    real(kind = kreal), intent(out) :: triad_new(3, 3)
    real(kind = kreal), intent(out) :: drill_new

    real(kind = kreal) :: triad_base(3, 3)
    real(kind = kreal) :: director(3), theta_step(3), theta_phys(3)
    real(kind = kreal) :: drill_acc, drill_inc, theta_norm
    real(kind = kreal) :: rot_inc(3, 3), rot_full(3, 3), triad_trial(3, 3)
    integer :: isub, nsub

    call ShellOrthonormalizeTriad( triad_old, triad_base )

    theta_norm = sqrt( sum( theta_inc(1:3)*theta_inc(1:3) ) )
    nsub = max( 1, ceiling( theta_norm/5.0D-2 ) )
    theta_step(1:3) = theta_inc(1:3)/dble(nsub)
    drill_acc = drill_old

    do isub = 1, nsub
      director(1:3) = triad_base(1:3, 3)
      drill_inc = dot_product( theta_step(1:3), director(1:3) )
      theta_phys(1:3) = theta_step(1:3) - drill_inc*director(1:3)

      call ShellRotationVectorToMatrix( theta_phys, rot_inc )
      triad_trial = matmul( rot_inc, triad_base )
      call ShellOrthonormalizeTriad( triad_trial, triad_base )
      drill_acc = drill_acc + drill_inc
    end do

    ! exact director update with separated drilling rotation
    call ShellRotationVectorToMatrix( theta_inc, rot_full )
    triad_trial(1:3, 1:3) = triad_base(1:3, 1:3)
    triad_trial(1:3, 3) = matmul( rot_full, triad_old(1:3, 3) )
    call ShellOrthonormalizeTriad( triad_trial, triad_new )
    drill_new = drill_acc

    return

    !####################################################################
  end subroutine ShellUpdateTriadWithIncrement
  !####################################################################


  !####################################################################
  subroutine ShellComposeNodalDisplacement( ndof, nn, disp_old, disp_inc, disp_new )
    !####################################################################

    integer(kind = kint), intent(in) :: ndof
    integer(kind = kint), intent(in) :: nn
    real(kind = kreal), intent(in)  :: disp_old(:, :)
    real(kind = kreal), intent(in)  :: disp_inc(:, :)
    real(kind = kreal), intent(out) :: disp_new(6, nn)

    integer :: i, ndof_copy

    disp_new(:, :) = 0.0D0
    ndof_copy = min( ndof, 6 )

    do i = 1, nn
      disp_new(1:min(3, ndof_copy), i) = disp_old(1:min(3, ndof_copy), i) &
        + disp_inc(1:min(3, ndof_copy), i)
      if( ndof_copy >= 6 ) then
        call ShellComposeRotationVector( disp_old(4:6, i), disp_inc(4:6, i), disp_new(4:6, i) )
      else if( ndof_copy > 3 ) then
        disp_new(4:ndof_copy, i) = disp_old(4:ndof_copy, i) + disp_inc(4:ndof_copy, i)
      endif
    end do

    return

    !####################################################################
  end subroutine ShellComposeNodalDisplacement
  !####################################################################


  !####################################################################
  pure subroutine ShellDirectorIncrementalDeriv( director_current, director_deriv )
    !####################################################################

    real(kind = kreal), intent(in)  :: director_current(3)
    real(kind = kreal), intent(out) :: director_deriv(3, 3)

    integer :: i
    real(kind = kreal) :: basis(3)

    do i = 1, 3
      basis(1:3) = 0.0D0
      basis(i) = 1.0D0

      director_deriv(1, i) = basis(2)*director_current(3) - basis(3)*director_current(2)
      director_deriv(2, i) = basis(3)*director_current(1) - basis(1)*director_current(3)
      director_deriv(3, i) = basis(1)*director_current(2) - basis(2)*director_current(1)
    end do

    return

    !####################################################################
  end subroutine ShellDirectorIncrementalDeriv
  !####################################################################


  !####################################################################
  pure subroutine ShellDirectorIncrementalSecondDeriv( director_current, director_second )
    !####################################################################

    real(kind = kreal), intent(in)  :: director_current(3)
    real(kind = kreal), intent(out) :: director_second(3, 3, 3)

    integer :: m, n
    real(kind = kreal) :: basis_m(3), basis_n(3)
    real(kind = kreal) :: cross_n(3), cross_mn(3), cross_m(3), cross_nm(3)

    do n = 1, 3
      basis_n(1:3) = 0.0D0
      basis_n(n) = 1.0D0
      cross_n(1) = basis_n(2)*director_current(3) - basis_n(3)*director_current(2)
      cross_n(2) = basis_n(3)*director_current(1) - basis_n(1)*director_current(3)
      cross_n(3) = basis_n(1)*director_current(2) - basis_n(2)*director_current(1)

      do m = 1, 3
        basis_m(1:3) = 0.0D0
        basis_m(m) = 1.0D0
        cross_m(1) = basis_m(2)*director_current(3) - basis_m(3)*director_current(2)
        cross_m(2) = basis_m(3)*director_current(1) - basis_m(1)*director_current(3)
        cross_m(3) = basis_m(1)*director_current(2) - basis_m(2)*director_current(1)

        cross_mn(1) = basis_m(2)*cross_n(3) - basis_m(3)*cross_n(2)
        cross_mn(2) = basis_m(3)*cross_n(1) - basis_m(1)*cross_n(3)
        cross_mn(3) = basis_m(1)*cross_n(2) - basis_m(2)*cross_n(1)
        cross_nm(1) = basis_n(2)*cross_m(3) - basis_n(3)*cross_m(2)
        cross_nm(2) = basis_n(3)*cross_m(1) - basis_n(1)*cross_m(3)
        cross_nm(3) = basis_n(1)*cross_m(2) - basis_n(2)*cross_m(1)

        director_second(1:3, m, n) = 0.5D0*( cross_mn(1:3) + cross_nm(1:3) )
      end do
    end do

    return

    !####################################################################
  end subroutine ShellDirectorIncrementalSecondDeriv
  !####################################################################

  !####################################################################
  subroutine ElementStress_Shell_MITC                  &
      (etype, nn, ndof, ecoord, gausses, edisp, &
      strain, stress, thick, zeta, n_layer, n_totlyr, surface_gauss_points, &
      local_strain, local_stress, local_stress_override, nddirector, ndrefdirector, &
      ndbase_disp)
    !####################################################################

    use m_fstr, only: OPSSTYPE, kOPSS_SOLUTION
    use m_utilities, only: get_principal
    use mMechGauss
    use m_MatMatrix
    use mMaterial, only: INFINITESIMAL, TOTALLAG, UPDATELAG, isElastic

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(in)   :: edisp(6, nn)
    real(kind = kreal), intent(out)  :: strain(:,:)
    real(kind = kreal), intent(out)  :: stress(:,:)
    real(kind = kreal), intent(in)   :: thick
    real(kind = kreal), intent(in)   :: zeta
    logical, intent(in), optional    :: surface_gauss_points
    real(kind = kreal), intent(out), optional :: local_strain(:,:)
    real(kind = kreal), intent(out), optional :: local_stress(:,:)
    real(kind = kreal), intent(in), optional :: local_stress_override(:,:)
    real(kind = kreal), intent(in), optional :: nddirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndrefdirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndbase_disp(6, nn)

    !--------------------------------------------------------------------

    integer :: i
    integer :: flag
    integer(kind=kint) :: ierr_quad
    integer :: lx, npoints
    integer :: fetype
    integer :: ntying
    integer :: npoints_tying(3)
    integer :: it, ip
    integer :: na, nb
    integer :: n_layer, n_totlyr

    real(kind = kreal) :: D(5, 5)
    real(kind = kreal) :: elem(3, nn)
    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: tpcoord(6, 2, 3)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: alpha
    real(kind = kreal) :: xxi_lx, eeta_lx
    real(kind = kreal) :: xxi_di(6, 3), eeta_di(6, 3)
    real(kind = kreal) :: h(nn, 3)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: a_over_2_theta_cross_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: theta(3, nn)
    real(kind = kreal) :: dudxi(3), dudeta(3), dudzeta(3)
    real(kind = kreal) :: dudxi_trans(3), dudeta_trans(3)
    real(kind = kreal) :: g1_cur(3), g2_cur(3), g3_cur(3)
    real(kind = kreal) :: g1_ref(3), g2_ref(3), g3_ref(3)
    real(kind = kreal) :: cg1_ref(3), cg2_ref(3), cg3_ref(3)
    real(kind = kreal) :: det_cur, det_ref, jac
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
      dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: g3_abs
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: cg1(3), cg2(3), cg3(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_inv
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: e1_hat(3), e2_hat(3), e3_hat(3)
    real(kind = kreal) :: e1_hat_abs, e2_hat_abs
    real(kind = kreal) :: e1_hat_mat(3), e2_hat_mat(3), e3_hat_mat(3)
    real(kind = kreal) :: cg1_mat(3), cg2_mat(3), cg3_mat(3)
    real(kind = kreal) :: det_mat
    real(kind = kreal) :: e11, e22, e12_2, e23_2, e31_2
    real(kind = kreal) :: e11_di(6, 3), e22_di(6, 3),     &
      e12_di_2(6, 3), e23_di_2(6, 3), &
      e31_di_2(6, 3)
    real(kind = kreal) :: E(3, 3), Ev(5)
    real(kind = kreal) :: S(3, 3), Sv(5)
    real(kind = kreal) :: stretch_b(3, 3), tensor(6), eigval(3), princ(3, 3)
    real(kind = kreal) :: cauchy(3, 3), logstrain(3, 3), cg_metric(3, 3)
    real(kind = kreal) :: eig_norm
    logical :: use_surface_gauss
    logical :: finite_rotation_director
    logical :: use_gl_strain


    zeta_ly = 0.0d0
    use_surface_gauss = .false.
    if( present( surface_gauss_points ) ) use_surface_gauss = surface_gauss_points
    flag = gausses(1)%pMaterial%nlgeom_flag
    finite_rotation_director = ( flag == TOTALLAG .or. flag == UPDATELAG ) &
      .and. etype == fe_mitc4_shell .and. nn == 4 &
      .and. isElastic(gausses(1)%pMaterial%mtype)
    ! Green-Lagrange strain and spatial output for elastic MITC4
    use_gl_strain = flag == TOTALLAG .and. etype == fe_mitc4_shell .and. nn == 4 &
      .and. isElastic(gausses(1)%pMaterial%mtype)

    !--------------------------------------------------------------------

    ! for lamina stress

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      fetype = fe_mitc4_shell

      ntying = 1
      npoints_tying(1)= 4

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      fetype = fe_mitc9_shell

      ntying = 3
      npoints_tying(1)= 6
      npoints_tying(2)= 6
      npoints_tying(3)= 4

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      fetype = fe_mitc3_shell

      ntying = 1
      npoints_tying(1)= 3

    end if

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)
    if( flag == UPDATELAG .and. present( ndbase_disp ) ) elem(:, :) = elem(:, :) + ndbase_disp(1:3, :)
    if( flag == UPDATELAG ) elem(:, :) = elem(:, :) + 0.5D0*edisp(1:3, :)

    !--------------------------------------------------------------------

    do na = 1, nn

      theta(1, na) = edisp(4, na)
      theta(2, na) = edisp(5, na)
      theta(3, na) = edisp(6, na)

    end do

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    !-------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) =  0.0D0
      tpcoord(2, 1, 1) =  1.0D0
      tpcoord(3, 1, 1) =  0.0D0
      tpcoord(4, 1, 1) = -1.0D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0
      tpcoord(2, 2, 1) =  0.0D0
      tpcoord(3, 2, 1) =  1.0D0
      tpcoord(4, 2, 1) =  0.0D0

      !--------------------------------------------------------

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      !--------------------------------------------------------

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 1)  = -1.0D0
      xxi_di(2, 1)  =  1.0D0
      xxi_di(3, 1)  =  1.0D0
      xxi_di(4, 1)  = -1.0D0
      xxi_di(5, 1)  =  1.0D0
      xxi_di(6, 1)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 1) = -1.0D0
      eeta_di(2, 1) = -1.0D0
      eeta_di(3, 1) =  1.0D0
      eeta_di(4, 1) =  1.0D0
      eeta_di(5, 1) =  0.0D0
      eeta_di(6, 1) =  0.0D0

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 2)  = -1.0D0
      xxi_di(2, 2)  =  0.0D0
      xxi_di(3, 2)  =  1.0D0
      xxi_di(4, 2)  =  1.0D0
      xxi_di(5, 2)  =  0.0D0
      xxi_di(6, 2)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 2) = -1.0D0
      eeta_di(2, 2) = -1.0D0
      eeta_di(3, 2) = -1.0D0
      eeta_di(4, 2) =  1.0D0
      eeta_di(5, 2) =  1.0D0
      eeta_di(6, 2) =  1.0D0

      !--------------------------------------------------------

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1)  = 0.5D0
      tpcoord(2, 1, 1)  = 0.0D0
      tpcoord(3, 1, 1)  = 0.5D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = 0.0D0
      tpcoord(2, 2, 1) = 0.5D0
      tpcoord(3, 2, 1) = 0.5D0

      !--------------------------------------------------------

    end if

    !--------------------------------------------------------------------

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    do i = 1, 3

      g1(i) = 0.0D0

      do na = 1, nn

        g1(i) = g1(i)+shapederiv(na, 1) &
          *elem(i, na)

      end do

    end do

    e_0(1) = g1(1)
    e_0(2) = g1(2)
    e_0(3) = g1(3)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1) &
            *elem(i, na)
          g2(i) = g2(i)+shapederiv(na, 2) &
            *elem(i, na)

        end do

      end do

      !--------------------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(1, nb) = det_cg3(1)/det_cg3_abs
      v3(2, nb) = det_cg3(2)/det_cg3_abs
      v3(3, nb) = det_cg3(3)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs .GT. 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !--------------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)
      if( finite_rotation_director .and. present( ndrefdirector ) ) then
        ! fixed nodal reference director
        a_over_2_v3(1:3, nb) = ndrefdirector(1:3, nb)
      endif

      !--------------------------------------------------------

      if( finite_rotation_director ) then
        if( present( nddirector ) ) then
          a_over_2_theta_cross_v3(1:3, nb) = nddirector(1:3, nb) - a_over_2_v3(1:3, nb)
          if( flag == UPDATELAG ) then
            a_over_2_v3(1:3, nb) = a_over_2_v3(1:3, nb) + 0.5D0*a_over_2_theta_cross_v3(1:3, nb)
          else
            a_over_2_v3(1:3, nb) = nddirector(1:3, nb)
          endif
        else
          call ShellDirectorIncrement( theta(1:3, nb), a_over_2_v3(1:3, nb), &
            a_over_2_theta_cross_v3(1:3, nb) )
          a_over_2_v3(1:3, nb) = a_over_2_v3(1:3, nb) &
            + a_over_2_theta_cross_v3(1:3, nb)
        endif
      else
        a_over_2_theta_cross_v3(1, nb)    &
          = theta(2, nb)*a_over_2_v3(3, nb) &
          -theta(3, nb)*a_over_2_v3(2, nb)
        a_over_2_theta_cross_v3(2, nb)    &
          = theta(3, nb)*a_over_2_v3(1, nb) &
          -theta(1, nb)*a_over_2_v3(3, nb)
        a_over_2_theta_cross_v3(3, nb)    &
          = theta(1, nb)*a_over_2_v3(2, nb) &
          -theta(2, nb)*a_over_2_v3(1, nb)
      endif

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------
    !     Modified stress in laminated shell
    !--------------------------------------------------------------------
    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      zeta_ly = 0.0D0

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      zeta_ly = zeta

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      zeta_ly = 0.0D0

    end if

    !---------------------------------------------------------

    do it = 1, ntying

      do ip = 1, npoints_tying(it)

        !-------------------------------------------------

        naturalcoord(1) = tpcoord(ip, 1, it)
        naturalcoord(2) = tpcoord(ip, 2, it)

        call getShapeFunc(fetype, naturalcoord, shapefunc)

        call getShapeDeriv(fetype, naturalcoord, shapederiv)

        !-------------------------------------------------

        do na = 1, nn

          do i = 1, 3

            u_rot(i, na)                      &
              = shapefunc(na)                   &
              *( zeta_ly*a_over_2_v3(i, na) )

            dudxi_rot(i, na)                  &
              = shapederiv(na, 1)               &
              *( zeta_ly*a_over_2_v3(i, na) )
            dudeta_rot(i, na)                 &
              = shapederiv(na, 2)               &
              *( zeta_ly*a_over_2_v3(i, na) )
            dudzeta_rot(i, na)                &
              = shapefunc(na)                   &
              *( a_over_2_v3(i, na) )

          end do

        end do

        !-------------------------------------------------

        ! Covariant basis vector
        do i = 1, 3

          g1(i) = 0.0D0
          g2(i) = 0.0D0
          g3(i) = 0.0D0

          do na = 1, nn

            g1(i) = g1(i)+shapederiv(na, 1)  &
              *elem(i, na)       &
              +dudxi_rot(i, na)
            g2(i) = g2(i)+shapederiv(na, 2)  &
              *elem(i, na)       &
              +dudeta_rot(i, na)
            g3(i) = g3(i)+dudzeta_rot(i, na)

          end do

        end do

        !---------------------------------------------

        do i = 1, 3

          dudxi(i)   = 0.0D0
          dudeta(i)  = 0.0D0
          dudzeta(i) = 0.0D0
          dudxi_trans(i) = 0.0D0
          dudeta_trans(i) = 0.0D0

          do na = 1, nn

            dudxi(i)                                      &
              = dudxi(i)                                    &
              +shapederiv(na, 1)                           &
              *( edisp(i, na)                             &
              +zeta_ly*a_over_2_theta_cross_v3(i, na) )
            dudeta(i)                                     &
              = dudeta(i)                                   &
              +shapederiv(na, 2)                           &
              *( edisp(i, na)                             &
              +zeta_ly*a_over_2_theta_cross_v3(i, na) )
            dudzeta(i)                                    &
              = dudzeta(i)                                  &
              +shapefunc(na)                               &
              *( a_over_2_theta_cross_v3(i, na) )
            dudxi_trans(i) = dudxi_trans(i)+shapederiv(na, 1)*edisp(i, na)
            dudeta_trans(i) = dudeta_trans(i)+shapederiv(na, 2)*edisp(i, na)

          end do

        end do

        !---------------------------------------------

        if( use_gl_strain ) then
          g1_cur(:) = g1(:)+dudxi_trans(:)
          g2_cur(:) = g2(:)+dudeta_trans(:)
          g3_cur(:) = g3(:)
          g1_ref(:) = g1_cur(:)-dudxi(:)
          g2_ref(:) = g2_cur(:)-dudeta(:)
          g3_ref(:) = g3_cur(:)-dudzeta(:)

          e11_di(ip, it) = 0.5D0*(dot_product(g1_cur, g1_cur)-dot_product(g1_ref, g1_ref))
          e22_di(ip, it) = 0.5D0*(dot_product(g2_cur, g2_cur)-dot_product(g2_ref, g2_ref))
          e12_2 = dot_product(g1_cur, g2_cur)-dot_product(g1_ref, g2_ref)
          e23_di_2(ip, it) = dot_product(g2_cur, g3_cur)-dot_product(g2_ref, g3_ref)
          e31_di_2(ip, it) = dot_product(g3_cur, g1_cur)-dot_product(g3_ref, g1_ref)
        else
          ! Infinitesimal strain tensor
          e11_di(ip, it)                             &
            = 0.5D0                                    &
            *( ( g1(1)*dudxi(1) +dudxi(1) *g1(1) )   &
            +( g1(2)*dudxi(2) +dudxi(2) *g1(2) )   &
            +( g1(3)*dudxi(3) +dudxi(3) *g1(3) ) )
          e22_di(ip, it)                             &
            = 0.5D0                                    &
            *( ( g2(1)*dudeta(1)+dudeta(1)*g2(1) )   &
            +( g2(2)*dudeta(2)+dudeta(2)*g2(2) )   &
            +( g2(3)*dudeta(3)+dudeta(3)*g2(3) ) )
          e12_2                                   &
            = ( g1(1)*dudeta(1) +dudxi(1)  *g2(1) ) &
            +( g1(2)*dudeta(2) +dudxi(2)  *g2(2) ) &
            +( g1(3)*dudeta(3) +dudxi(3)  *g2(3) )
          e23_di_2(ip, it)                        &
            = ( g2(1)*dudzeta(1)+dudeta(1) *g3(1) ) &
            +( g2(2)*dudzeta(2)+dudeta(2) *g3(2) ) &
            +( g2(3)*dudzeta(3)+dudeta(3) *g3(3) )
          e31_di_2(ip, it)                        &
            = ( g3(1)*dudxi(1)  +dudzeta(1)*g1(1) ) &
            +( g3(2)*dudxi(2)  +dudzeta(2)*g1(2) ) &
            +( g3(3)*dudxi(3)  +dudzeta(3)*g1(3) )
        endif

        !-------------------------------------------------

      end do

    end do

    !--------------------------------------------------------

    call fstr_shell_layer_zeta( gausses(1), n_layer, zeta, zeta_ly, ierr_quad )
    if( ierr_quad /= 0 ) stop "Invalid shell layer zeta"

    !--------------------------------------------------------

    npoints = nn
    if( use_surface_gauss ) npoints = NumOfQuadPoints(fetype)

    do lx = 1, npoints

      !--------------------------------------------------

      if( use_surface_gauss ) then
        call getQuadPoint(fetype, lx, naturalcoord)
      else
        naturalcoord(1) = nncoord(lx, 1)
        naturalcoord(2) = nncoord(lx, 2)
      endif

      xi_lx  = naturalcoord(1)
      eta_lx = naturalcoord(2)

      call getShapeFunc(fetype, naturalcoord, shapefunc)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------

      do na = 1, nn

        do i = 1, 3

          u_rot(i, na)                      &
            = shapefunc(na)                   &
            *( zeta_ly*a_over_2_v3(i, na) )

          dudxi_rot(i, na)                  &
            = shapederiv(na, 1)               &
            *( zeta_ly*a_over_2_v3(i, na) )
          dudeta_rot(i, na)                 &
            = shapederiv(na, 2)                &
            *( zeta_ly*a_over_2_v3(i, na) )
          dudzeta_rot(i, na)                &
            = shapefunc(na)                   &
            *( a_over_2_v3(i, na) )

        end do

      end do

      !--------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0
        g3(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1)  &
            *elem(i, na)       &
            +dudxi_rot(i, na)
          g2(i) = g2(i)+shapederiv(na, 2)  &
            *elem(i, na)       &
            +dudeta_rot(i, na)
          g3(i) = g3(i)+dudzeta_rot(i, na)

        end do

      end do

      !--------------------------------------------------

      if( use_gl_strain .and. etype == fe_mitc4_shell ) then
        dudxi_trans(:) = 0.0D0
        dudeta_trans(:) = 0.0D0
        do na = 1, nn
          dudxi_trans(:) = dudxi_trans(:)+shapederiv(na, 1)*edisp(1:3, na)
          dudeta_trans(:) = dudeta_trans(:)+shapederiv(na, 2)*edisp(1:3, na)
        end do
        g1(:) = g1(:)+dudxi_trans(:)
        g2(:) = g2(:)+dudeta_trans(:)
      endif

      !--------------------------------------------------

      ! Jacobian
      det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
        +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
        +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )
      det_cur = det

      if(det == 0.0d0) then
        write(*,*)"ERROR:LIB Shell in l2009 Not Jacobian"
        stop
      endif

      det_inv = 1.0D0/det

      !--------------------------------------------------

      ! Contravariant basis vector
      cg1(1) = det_inv                      &
        *( g2(2)*g3(3)-g2(3)*g3(2) )
      cg1(2) = det_inv                      &
        *( g2(3)*g3(1)-g2(1)*g3(3) )
      cg1(3) = det_inv                      &
        *( g2(1)*g3(2)-g2(2)*g3(1) )
      cg2(1) = det_inv                      &
        *( g3(2)*g1(3)-g3(3)*g1(2) )
      cg2(2) = det_inv                      &
        *( g3(3)*g1(1)-g3(1)*g1(3) )
      cg2(3) = det_inv                      &
        *( g3(1)*g1(2)-g3(2)*g1(1) )
      cg3(1) = det_inv                      &
        *( g1(2)*g2(3)-g1(3)*g2(2) )
      cg3(2) = det_inv                      &
        *( g1(3)*g2(1)-g1(1)*g2(3) )
      cg3(3) = det_inv                      &
        *( g1(1)*g2(2)-g1(2)*g2(1) )

      !--------------------------------------------------

      g3_abs = dsqrt( g3(1)*g3(1)   &
        +g3(2)*g3(2)   &
        +g3(3)*g3(3) )

      !--------------------------------------------------

      ! Orthonormal vectors

      e3_hat(1) = g3(1)/g3_abs
      e3_hat(2) = g3(2)/g3_abs
      e3_hat(3) = g3(3)/g3_abs

      e1_hat(1) = g2(2)*e3_hat(3) &
        -g2(3)*e3_hat(2)
      e1_hat(2) = g2(3)*e3_hat(1) &
        -g2(1)*e3_hat(3)
      e1_hat(3) = g2(1)*e3_hat(2) &
        -g2(2)*e3_hat(1)
      e1_hat_abs = dsqrt( e1_hat(1)*e1_hat(1)   &
        +e1_hat(2)*e1_hat(2)   &
        +e1_hat(3)*e1_hat(3) )
      e1_hat(1) = e1_hat(1)/e1_hat_abs
      e1_hat(2) = e1_hat(2)/e1_hat_abs
      e1_hat(3) = e1_hat(3)/e1_hat_abs

      e2_hat(1) = e3_hat(2)*e1_hat(3) &
        -e3_hat(3)*e1_hat(2)
      e2_hat(2) = e3_hat(3)*e1_hat(1) &
        -e3_hat(1)*e1_hat(3)
      e2_hat(3) = e3_hat(1)*e1_hat(2) &
        -e3_hat(2)*e1_hat(1)
      e2_hat_abs = dsqrt( e2_hat(1)*e2_hat(1)   &
        +e2_hat(2)*e2_hat(2)   &
        +e2_hat(3)*e2_hat(3) )
      e2_hat(1) = e2_hat(1)/e2_hat_abs
      e2_hat(2) = e2_hat(2)/e2_hat_abs
      e2_hat(3) = e2_hat(3)/e2_hat_abs

      !--------------------------------------------------

      do i = 1, 3

        dudxi(i)   = 0.0D0
        dudeta(i)  = 0.0D0
        dudzeta(i) = 0.0D0
        dudxi_trans(i) = 0.0D0
        dudeta_trans(i) = 0.0D0

        do na = 1, nn

          dudxi(i)                                      &
            = dudxi(i)                                    &
            +shapederiv(na, 1)                           &
            *( edisp(i, na)                             &
            +zeta_ly*a_over_2_theta_cross_v3(i, na) )
          dudeta(i)                                     &
            = dudeta(i)                                   &
            +shapederiv(na, 2)                           &
            *( edisp(i, na)                             &
            +zeta_ly*a_over_2_theta_cross_v3(i, na) )
          dudzeta(i)                                    &
            = dudzeta(i)                                  &
            +shapefunc(na)                               &
            *( a_over_2_theta_cross_v3(i, na) )
          dudxi_trans(i) = dudxi_trans(i)+shapederiv(na, 1)*edisp(i, na)
          dudeta_trans(i) = dudeta_trans(i)+shapederiv(na, 2)*edisp(i, na)

        end do

      end do

      !--------------------------------------------------

      !--------------------------------------------------

      if( use_gl_strain ) then
        if( etype == fe_mitc4_shell ) then
          g1_cur(:) = g1(:)
          g2_cur(:) = g2(:)
        else
          g1_cur(:) = g1(:)+dudxi_trans(:)
          g2_cur(:) = g2(:)+dudeta_trans(:)
        endif
        g3_cur(:) = g3(:)
        g1_ref(:) = g1_cur(:)-dudxi(:)
        g2_ref(:) = g2_cur(:)-dudeta(:)
        g3_ref(:) = g3_cur(:)-dudzeta(:)
        det_ref = g1_ref(1)*( g2_ref(2)*g3_ref(3)-g2_ref(3)*g3_ref(2) ) &
          +g1_ref(2)*( g2_ref(3)*g3_ref(1)-g2_ref(1)*g3_ref(3) ) &
          +g1_ref(3)*( g2_ref(1)*g3_ref(2)-g2_ref(2)*g3_ref(1) )

        e11 = 0.5D0*(dot_product(g1_cur, g1_cur)-dot_product(g1_ref, g1_ref))
        e22 = 0.5D0*(dot_product(g2_cur, g2_cur)-dot_product(g2_ref, g2_ref))
        e12_2 = dot_product(g1_cur, g2_cur)-dot_product(g1_ref, g2_ref)
        e23_2 = dot_product(g2_cur, g3_cur)-dot_product(g2_ref, g3_ref)
        e31_2 = dot_product(g3_cur, g1_cur)-dot_product(g3_ref, g1_ref)
      else
        ! Infinitesimal strain tensor
        e11                                        &
          = 0.5D0                                    &
          *( ( g1(1)*dudxi(1) +dudxi(1) *g1(1) )   &
          +( g1(2)*dudxi(2) +dudxi(2) *g1(2) )   &
          +( g1(3)*dudxi(3) +dudxi(3) *g1(3) ) )
        e22                                        &
          = 0.5D0                                    &
          *( ( g2(1)*dudeta(1)+dudeta(1)*g2(1) )   &
          +( g2(2)*dudeta(2)+dudeta(2)*g2(2) )   &
          +( g2(3)*dudeta(3)+dudeta(3)*g2(3) ) )
        e12_2                                   &
          = ( g1(1)*dudeta(1) +dudxi(1)  *g2(1) ) &
          +( g1(2)*dudeta(2) +dudxi(2)  *g2(2) ) &
          +( g1(3)*dudeta(3) +dudxi(3)  *g2(3) )
        e23_2                                   &
          = ( g2(1)*dudzeta(1)+dudeta(1) *g3(1) ) &
          +( g2(2)*dudzeta(2)+dudeta(2) *g3(2) ) &
          +( g2(3)*dudzeta(3)+dudeta(3) *g3(3) )
        e31_2                                   &
          = ( g3(1)*dudxi(1)  +dudzeta(1)*g1(1) ) &
          +( g3(2)*dudxi(2)  +dudzeta(2)*g1(2) ) &
          +( g3(3)*dudxi(3)  +dudzeta(3)*g1(3) )
      endif

      ! MITC4
      if( etype .EQ. fe_mitc4_shell ) then

        e23_2 = 0.0D0
        e31_2 = 0.0D0

        ! e23_as_2
        e23_2                                   &
          = 0.5D0*( 1.0D0-xi_lx  )*e23_di_2(4, 1) &
          +0.5D0*( 1.0D0+xi_lx  )*e23_di_2(2, 1)
        ! e31_as_2
        e31_2                                   &
          = 0.5D0*( 1.0D0-eta_lx )*e31_di_2(1, 1) &
          +0.5D0*( 1.0D0+eta_lx )*e31_di_2(3, 1)

        ! MITC9
      else if( etype .EQ. fe_mitc9_shell ) then

        xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
        eeta_lx = eta_lx/dsqrt( 3.0D0/5.0D0 )

        do ip = 1, npoints_tying(1)

          h(ip, 1)                                     &
            = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )   &
            *( ( 0.5D0*eeta_di(ip, 1)*eeta_lx )        &
            *( 1.0D0+eeta_di(ip, 1)*eeta_lx )       &
            +( 1.0D0-eeta_di(ip, 1)*eeta_di(ip, 1) ) &
            *( 1.0D0-eeta_lx*eeta_lx ) )

        end do

        xxi_lx  = xi_lx /dsqrt( 3.0D0/5.0D0 )
        eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

        do ip = 1, npoints_tying(2)

          h(ip, 2)                                      &
            = ( ( 0.5D0*xxi_di(ip, 2) *xxi_lx  )          &
            *( 1.0D0+xxi_di(ip, 2) *xxi_lx  )         &
            +( 1.0D0-xxi_di(ip, 2) *xxi_di(ip, 2)  )   &
            *( 1.0D0-xxi_lx*xxi_lx ) )                &
            *( 0.5D0*( 1.0D0+eeta_di(ip, 2)*eeta_lx ) )

        end do

        xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
        eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

        do ip = 1, npoints_tying(3)

          h(ip, 3)                                      &
            = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )    &
            *( 0.5D0*( 1.0D0+eeta_di(ip, 1)*eeta_lx ) )

        end do

        e11   = 0.0D0
        e31_2 = 0.0D0

        ! e11_as, e31_as_2
        do ip = 1, npoints_tying(1)

          e11   = e11  +h(ip, 1)*e11_di(ip, 1)
          e31_2 = e31_2+h(ip, 1)*e31_di_2(ip, 1)

        end do

        e22   = 0.0D0
        e23_2 = 0.0D0

        ! e22_as, e23_as_2
        do ip = 1, npoints_tying(2)

          e22   = e22  +h(ip, 2)*e22_di(ip, 2)
          e23_2 = e23_2+h(ip, 2)*e23_di_2(ip, 2)

        end do

        e12_2 = 0.0D0

        ! e12_as_2
        do ip = 1, npoints_tying(3)

          e12_2 = e12_2+h(ip, 3)*e12_di_2(ip, 3)

        end do

        ! MITC3
      else if( etype .EQ. fe_mitc3_shell ) then

        e23_2 = 0.0D0
        e31_2 = 0.0D0

        ! e23_as_2
        e23_2                                       &
          = ( 1.0D0-xi_lx )*e23_di_2(2, 1)           &
          +xi_lx *e31_di_2(1, 1)                    &
          +xi_lx *( e23_di_2(3, 1)-e31_di_2(3, 1) )
        ! e31_as_2
        e31_2                                       &
          = eta_lx*e23_di_2(2, 1)                    &
          +( 1.0D0-eta_lx )*e31_di_2(1, 1)          &
          -eta_lx*( e23_di_2(3, 1)-e31_di_2(3, 1) )

      end if

      !--------------------------------------------------

      ! { E } vector
      Ev(1) = e11
      Ev(2) = e22
      Ev(3) = e12_2
      Ev(4) = e23_2
      Ev(5) = e31_2

      ! Infinitesimal strain tensor
      ! [ E ] matrix
      E(1, 1) = Ev(1)
      E(2, 2) = Ev(2)
      E(3, 3) = 0.0D0
      E(1, 2) = 0.5D0*Ev(3)
      E(2, 1) = 0.5D0*Ev(3)
      E(2, 3) = 0.5D0*Ev(4)
      E(3, 2) = 0.5D0*Ev(4)
      E(3, 1) = 0.5D0*Ev(5)
      E(1, 3) = 0.5D0*Ev(5)

      !--------------------------------------------------
      !   write(*,*) 'Stress_n_layer', n_layer
      e1_hat_mat(:) = e1_hat(:)
      e2_hat_mat(:) = e2_hat(:)
      e3_hat_mat(:) = e3_hat(:)
      cg1_mat(:) = cg1(:)
      cg2_mat(:) = cg2(:)
      cg3_mat(:) = cg3(:)
      if( use_gl_strain .and. etype == fe_mitc4_shell ) then
        call shell_basis_from_covariant(g1_ref, g2_ref, g3_ref, &
          e1_hat_mat, e2_hat_mat, e3_hat_mat, cg1_mat, cg2_mat, cg3_mat, det_mat)
      endif

      call MatlMatrix_Shell                        &
        (gausses(lx), Shell, D,                 &
        e1_hat_mat, e2_hat_mat, e3_hat_mat, cg1_mat, cg2_mat, cg3_mat, &
        alpha, n_layer)

      !--------------------------------------------------

      Sv = matmul( D, Ev )

      if( present( local_stress_override ) ) then
        if( lx <= size(local_stress_override, 1) .and. size(local_stress_override, 2) >= 6 ) then
          Sv(1) = local_stress_override(lx, 1)
          Sv(2) = local_stress_override(lx, 2)
          Sv(3) = local_stress_override(lx, 4)
          Sv(4) = local_stress_override(lx, 5)
          Sv(5) = local_stress_override(lx, 6)
        endif
      endif

      ! Infinitesimal stress tensor
      ! [ S ] matrix
      S(1, 1) = Sv(1)
      S(2, 2) = Sv(2)
      S(3, 3) = 0.0D0
      S(1, 2) = Sv(3)
      S(2, 1) = Sv(3)
      S(2, 3) = Sv(4)
      S(3, 2) = Sv(4)
      S(3, 1) = Sv(5)
      S(1, 3) = Sv(5)

      !------------------------------------------------

      if( present( local_strain ) ) then
        local_strain(lx, 1) = Ev(1)
        local_strain(lx, 2) = Ev(2)
        local_strain(lx, 3) = 0.0d0
        local_strain(lx, 4) = Ev(3)
        local_strain(lx, 5) = Ev(4)
        local_strain(lx, 6) = Ev(5)
      endif

      if( present( local_stress ) ) then
        local_stress(lx, 1) = Sv(1)
        local_stress(lx, 2) = Sv(2)
        local_stress(lx, 3) = 0.0d0
        local_stress(lx, 4) = Sv(3)
        local_stress(lx, 5) = Sv(4)
        local_stress(lx, 6) = Sv(5)
      endif

      if( use_gl_strain .and. etype == fe_mitc4_shell ) then
        g1(:) = g1_ref(:)
        g2(:) = g2_ref(:)
        g3(:) = g3_ref(:)

        det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
          +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
          +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )
        if(det == 0.0d0) then
          write(*,*)"ERROR:LIB Shell in l2009 Not Jacobian"
          stop
        endif

        det_inv = 1.0D0/det
        cg1(1) = det_inv*( g2(2)*g3(3)-g2(3)*g3(2) )
        cg1(2) = det_inv*( g2(3)*g3(1)-g2(1)*g3(3) )
        cg1(3) = det_inv*( g2(1)*g3(2)-g2(2)*g3(1) )
        cg2(1) = det_inv*( g3(2)*g1(3)-g3(3)*g1(2) )
        cg2(2) = det_inv*( g3(3)*g1(1)-g3(1)*g1(3) )
        cg2(3) = det_inv*( g3(1)*g1(2)-g3(2)*g1(1) )
        cg3(1) = det_inv*( g1(2)*g2(3)-g1(3)*g2(2) )
        cg3(2) = det_inv*( g1(3)*g2(1)-g1(1)*g2(3) )
        cg3(3) = det_inv*( g1(1)*g2(2)-g1(2)*g2(1) )
        cg1_ref(:) = cg1(:)
        cg2_ref(:) = cg2(:)
        cg3_ref(:) = cg3(:)
      endif

      stress(lx,1)                  &
        = ( S(1, 1)*g1(1)*g1(1)   &
        +S(1, 2)*g1(1)*g2(1)   &
        +S(1, 3)*g1(1)*g3(1)   &
        +S(2, 1)*g2(1)*g1(1)   &
        +S(2, 2)*g2(1)*g2(1)   &
        +S(2, 3)*g2(1)*g3(1)   &
        +S(3, 1)*g3(1)*g1(1)   &
        +S(3, 2)*g3(1)*g2(1) )
      stress(lx,2)                 &
        = ( S(1, 1)*g1(2)*g1(2)   &
        +S(1, 2)*g1(2)*g2(2)   &
        +S(1, 3)*g1(2)*g3(2)   &
        +S(2, 1)*g2(2)*g1(2)   &
        +S(2, 2)*g2(2)*g2(2)   &
        +S(2, 3)*g2(2)*g3(2)   &
        +S(3, 1)*g3(2)*g1(2)   &
        +S(3, 2)*g3(2)*g2(2) )
      stress(lx,3)                 &
        = ( S(1, 1)*g1(3)*g1(3)   &
        +S(1, 2)*g1(3)*g2(3)   &
        +S(1, 3)*g1(3)*g3(3)   &
        +S(2, 1)*g2(3)*g1(3)   &
        +S(2, 2)*g2(3)*g2(3)   &
        +S(2, 3)*g2(3)*g3(3)   &
        +S(3, 1)*g3(3)*g1(3)   &
        +S(3, 2)*g3(3)*g2(3) )
      stress(lx,4)                 &
        = ( S(1, 1)*g1(1)*g1(2)   &
        +S(1, 2)*g1(1)*g2(2)   &
        +S(1, 3)*g1(1)*g3(2)   &
        +S(2, 1)*g2(1)*g1(2)   &
        +S(2, 2)*g2(1)*g2(2)   &
        +S(2, 3)*g2(1)*g3(2)   &
        +S(3, 1)*g3(1)*g1(2)   &
        +S(3, 2)*g3(1)*g2(2) )
      stress(lx,5)                 &
        = ( S(1, 1)*g1(2)*g1(3)   &
        +S(1, 2)*g1(2)*g2(3)   &
        +S(1, 3)*g1(2)*g3(3)   &
        +S(2, 1)*g2(2)*g1(3)   &
        +S(2, 2)*g2(2)*g2(3)   &
        +S(2, 3)*g2(2)*g3(3)   &
        +S(3, 1)*g3(2)*g1(3)   &
        +S(3, 2)*g3(2)*g2(3) )
      stress(lx,6)                 &
        = ( S(1, 1)*g1(3)*g1(1)   &
        +S(1, 2)*g1(3)*g2(1)   &
        +S(1, 3)*g1(3)*g3(1)   &
        +S(2, 1)*g2(3)*g1(1)   &
        +S(2, 2)*g2(3)*g2(1)   &
        +S(2, 3)*g2(3)*g3(1)   &
        +S(3, 1)*g3(3)*g1(1)   &
        +S(3, 2)*g3(3)*g2(1) )

      strain(lx,1)                   &
        = ( E(1, 1)*cg1(1)*cg1(1)   &
        +E(1, 2)*cg1(1)*cg2(1)   &
        +E(1, 3)*cg1(1)*cg3(1)   &
        +E(2, 1)*cg2(1)*cg1(1)   &
        +E(2, 2)*cg2(1)*cg2(1)   &
        +E(2, 3)*cg2(1)*cg3(1)   &
        +E(3, 1)*cg3(1)*cg1(1)   &
        +E(3, 2)*cg3(1)*cg2(1) )
      strain(lx,2)                   &
        = ( E(1, 1)*cg1(2)*cg1(2)   &
        +E(1, 2)*cg1(2)*cg2(2)   &
        +E(1, 3)*cg1(2)*cg3(2)   &
        +E(2, 1)*cg2(2)*cg1(2)   &
        +E(2, 2)*cg2(2)*cg2(2)   &
        +E(2, 3)*cg2(2)*cg3(2)   &
        +E(3, 1)*cg3(2)*cg1(2)   &
        +E(3, 2)*cg3(2)*cg2(2) )
      strain(lx,3)                   &
        = ( E(1, 1)*cg1(3)*cg1(3)   &
        +E(1, 2)*cg1(3)*cg2(3)   &
        +E(1, 3)*cg1(3)*cg3(3)   &
        +E(2, 1)*cg2(3)*cg1(3)   &
        +E(2, 2)*cg2(3)*cg2(3)   &
        +E(2, 3)*cg2(3)*cg3(3)   &
        +E(3, 1)*cg3(3)*cg1(3)   &
        +E(3, 2)*cg3(3)*cg2(3) )
      strain(lx,4)                   &
        = ( E(1, 1)*cg1(1)*cg1(2)   &
        +E(1, 2)*cg1(1)*cg2(2)   &
        +E(1, 3)*cg1(1)*cg3(2)   &
        +E(2, 1)*cg2(1)*cg1(2)   &
        +E(2, 2)*cg2(1)*cg2(2)   &
        +E(2, 3)*cg2(1)*cg3(2)   &
        +E(3, 1)*cg3(1)*cg1(2)   &
        +E(3, 2)*cg3(1)*cg2(2) )
      strain(lx,5)                   &
        = ( E(1, 1)*cg1(2)*cg1(3)   &
        +E(1, 2)*cg1(2)*cg2(3)   &
        +E(1, 3)*cg1(2)*cg3(3)   &
        +E(2, 1)*cg2(2)*cg1(3)   &
        +E(2, 2)*cg2(2)*cg2(3)   &
        +E(2, 3)*cg2(2)*cg3(3)   &
        +E(3, 1)*cg3(2)*cg1(3)   &
        +E(3, 2)*cg3(2)*cg2(3) )
      strain(lx,6)                   &
        = ( E(1, 1)*cg1(3)*cg1(1)   &
        +E(1, 2)*cg1(3)*cg2(1)   &
        +E(1, 3)*cg1(3)*cg3(1)   &
        +E(2, 1)*cg2(3)*cg1(1)   &
        +E(2, 2)*cg2(3)*cg2(1)   &
        +E(2, 3)*cg2(3)*cg3(1)   &
        +E(3, 1)*cg3(3)*cg1(1)   &
        +E(3, 2)*cg3(3)*cg2(1) )

      if( use_gl_strain .and. etype == fe_mitc4_shell .and. OPSSTYPE == kOPSS_SOLUTION ) then
        jac = det_cur/det_ref
        if( jac == 0.0d0 ) stop "Fail to convert shell stress: detF=0"

        cauchy(:, :) = 0.0d0
        cauchy(:, :) = cauchy(:, :) &
          +S(1, 1)*outer_product3(g1_cur, g1_cur) &
          +S(1, 2)*outer_product3(g1_cur, g2_cur) &
          +S(1, 3)*outer_product3(g1_cur, g3_cur) &
          +S(2, 1)*outer_product3(g2_cur, g1_cur) &
          +S(2, 2)*outer_product3(g2_cur, g2_cur) &
          +S(2, 3)*outer_product3(g2_cur, g3_cur) &
          +S(3, 1)*outer_product3(g3_cur, g1_cur) &
          +S(3, 2)*outer_product3(g3_cur, g2_cur) &
          +S(3, 3)*outer_product3(g3_cur, g3_cur)
        cauchy(:, :) = cauchy(:, :)/jac

        cg_metric(1, 1) = dot_product(cg1_ref, cg1_ref)
        cg_metric(1, 2) = dot_product(cg1_ref, cg2_ref)
        cg_metric(1, 3) = dot_product(cg1_ref, cg3_ref)
        cg_metric(2, 1) = cg_metric(1, 2)
        cg_metric(2, 2) = dot_product(cg2_ref, cg2_ref)
        cg_metric(2, 3) = dot_product(cg2_ref, cg3_ref)
        cg_metric(3, 1) = cg_metric(1, 3)
        cg_metric(3, 2) = cg_metric(2, 3)
        cg_metric(3, 3) = dot_product(cg3_ref, cg3_ref)

        stretch_b(:, :) = 0.0d0
        stretch_b(:, :) = stretch_b(:, :) &
          +cg_metric(1, 1)*outer_product3(g1_cur, g1_cur) &
          +cg_metric(1, 2)*outer_product3(g1_cur, g2_cur) &
          +cg_metric(1, 3)*outer_product3(g1_cur, g3_cur) &
          +cg_metric(2, 1)*outer_product3(g2_cur, g1_cur) &
          +cg_metric(2, 2)*outer_product3(g2_cur, g2_cur) &
          +cg_metric(2, 3)*outer_product3(g2_cur, g3_cur) &
          +cg_metric(3, 1)*outer_product3(g3_cur, g1_cur) &
          +cg_metric(3, 2)*outer_product3(g3_cur, g2_cur) &
          +cg_metric(3, 3)*outer_product3(g3_cur, g3_cur)

        tensor(1) = stretch_b(1, 1)
        tensor(2) = stretch_b(2, 2)
        tensor(3) = stretch_b(3, 3)
        tensor(4) = stretch_b(1, 2)
        tensor(5) = stretch_b(2, 3)
        tensor(6) = stretch_b(3, 1)
        call get_principal(tensor, eigval, princ)

        do i = 1, 3
          if( eigval(i) <= 0.0d0 ) stop "Fail to calc shell log strain: stretch<0"
          eigval(i) = 0.5d0*dlog(eigval(i))
          eig_norm = dsqrt(dot_product(princ(1:3, i), princ(1:3, i)))
          if( eig_norm <= 0.0d0 ) stop "Fail to calc shell log strain: direction vector=0"
          princ(1:3, i) = princ(1:3, i)/eig_norm
        end do

        logstrain(:, :) = 0.0d0
        do i = 1, 3
          logstrain(:, :) = logstrain(:, :) &
            +eigval(i)*outer_product3(princ(1:3, i), princ(1:3, i))
        end do

        stress(lx, 1) = cauchy(1, 1)
        stress(lx, 2) = cauchy(2, 2)
        stress(lx, 3) = cauchy(3, 3)
        stress(lx, 4) = cauchy(1, 2)
        stress(lx, 5) = cauchy(2, 3)
        stress(lx, 6) = cauchy(3, 1)

        strain(lx, 1) = logstrain(1, 1)
        strain(lx, 2) = logstrain(2, 2)
        strain(lx, 3) = logstrain(3, 3)
        strain(lx, 4) = 2.0d0*logstrain(1, 2)
        strain(lx, 5) = 2.0d0*logstrain(2, 3)
        strain(lx, 6) = 2.0d0*logstrain(3, 1)
      endif

      !--------------------------------------------------

    end do
    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine ElementStress_Shell_MITC
  !####################################################################

  !####################################################################
  subroutine UpdateStressShellUL_Elastic( gauss, dstrain, dstress, trace_coeff )
    !####################################################################

    use mMechGauss

    !--------------------------------------------------------------------

    type(tGaussStatus), intent(inout) :: gauss
    real(kind = kreal), intent(in)    :: dstrain(6)
    real(kind = kreal), intent(in)    :: dstress(6)
    real(kind = kreal), intent(in), optional :: trace_coeff

    !--------------------------------------------------------------------

    real(kind = kreal) :: dstress_obj(6)

    !--------------------------------------------------------------------

    call ShellObjectiveStressIncrement( gauss%stress_bak(1:6), dstrain(1:6), dstress_obj, trace_coeff )

    gauss%strain(1:6) = gauss%strain_bak(1:6) + dstrain(1:6)
    gauss%stress(1:6) = gauss%stress_bak(1:6) + dstress_obj(1:6) + dstress(1:6)

    gauss%strain_energy = gauss%strain_energy_bak &
      + dot_product( gauss%stress(1:6), dstrain(1:6) )
    gauss%strain_energy = gauss%strain_energy &
      - 0.5D0*dot_product( dstress(1:6), dstrain(1:6) )
    gauss%strain_out(1:6) = gauss%strain(1:6)
    gauss%stress_out(1:6) = gauss%stress(1:6)

    return

    !####################################################################
  end subroutine UpdateStressShellUL_Elastic
  !####################################################################

  !####################################################################
  subroutine UpdateShellLayerGauss_Shell_MITC                  &
      (etype, nn, ndof, ecoord, element, edisp, thick, nddirector, ndrefdirector, ndbase_disp)
    !####################################################################

    use mMechGauss
    use Quadrature
    use mMaterial, only: UPDATELAG, isElastic

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    type(tElement), intent(inout)    :: element
    real(kind = kreal), intent(in)   :: edisp(6, nn)
    real(kind = kreal), intent(in)   :: thick
    real(kind = kreal), intent(in), optional :: nddirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndrefdirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndbase_disp(6, nn)

    !--------------------------------------------------------------------

    integer(kind = kint) :: ng, ig, ilayer, ithick, ishell, ierr
    integer(kind = kint) :: flag
    real(kind = kreal) :: zeta, weight
    real(kind = kreal) :: trace_coeff
    real(kind = kreal) :: gpstrain(9, 6), gpstress(9, 6)
    real(kind = kreal) :: local_gpstrain(9, 6), local_gpstress(9, 6)

    !--------------------------------------------------------------------

    if( .not. associated( element%gausses ) ) return
    if( .not. associated( element%shell_layer_gausses ) ) return
    if( element%shell_nlayer <= 0 .or. element%shell_nthick <= 0 ) return
    if( .not. ( etype == fe_mitc4_shell .and. nn == 4 &
      .and. isElastic(element%gausses(1)%pMaterial%mtype) ) ) return

    flag = element%gausses(1)%pMaterial%nlgeom_flag
    ng = NumOfQuadPoints( etype )
    if( ng <= 0 .or. ng > 9 ) return

    do ilayer = 1, element%shell_nlayer
      do ithick = 1, element%shell_nthick
        call fstr_shell_thickness_quadrature( etype, ithick, zeta, weight, ierr )
        if( ierr /= 0 ) cycle

        gpstrain(:, :) = 0.0d0
        gpstress(:, :) = 0.0d0
        local_gpstrain(:, :) = 0.0d0
        local_gpstress(:, :) = 0.0d0
        call ElementStress_Shell_MITC( etype, nn, ndof, ecoord, element%gausses, edisp, &
          gpstrain(1:ng, 1:6), gpstress(1:ng, 1:6), thick, zeta, ilayer, element%shell_nlayer, &
          surface_gauss_points=.true., local_strain=local_gpstrain(1:ng, 1:6), &
          local_stress=local_gpstress(1:ng, 1:6), nddirector=nddirector, ndrefdirector=ndrefdirector, &
          ndbase_disp=ndbase_disp )

        do ig = 1, ng
          ishell = fstr_shell_layer_gauss_index( element, ig, ilayer, ithick )
          if( ishell <= 0 ) cycle
          if( flag == UPDATELAG ) then
            trace_coeff = ShellPlaneStressTraceCoeff(element%shell_layer_gausses(ishell), ilayer)
            call UpdateStressShellUL_Elastic( element%shell_layer_gausses(ishell), &
              local_gpstrain(ig, 1:6), local_gpstress(ig, 1:6), trace_coeff=trace_coeff )
          else
            element%shell_layer_gausses(ishell)%strain(1:6) = local_gpstrain(ig, 1:6)
            element%shell_layer_gausses(ishell)%stress(1:6) = local_gpstress(ig, 1:6)
            element%shell_layer_gausses(ishell)%strain_energy = &
              0.5d0*dot_product( local_gpstress(ig, 1:6), local_gpstrain(ig, 1:6) )
            element%shell_layer_gausses(ishell)%strain_out(1:6) = gpstrain(ig, 1:6)
            element%shell_layer_gausses(ishell)%stress_out(1:6) = gpstress(ig, 1:6)
          endif
        enddo
      enddo
    enddo

    return

    !####################################################################
  end subroutine UpdateShellLayerGauss_Shell_MITC
  !####################################################################

  !####################################################################
  subroutine DL_Shell                                  &
      (etype, nn, ndof, xx, yy, zz, rho, thick, &
      ltype, params, vect, nsize, gausses)
    !####################################################################

    use hecmw
    use m_utilities
    use mMechGauss
    use Quadrature

    type(tGaussStatus), intent(in)   :: gausses(:)
    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in) :: xx(*), yy(*), zz(*)
    real(kind = kreal), intent(in) :: rho
    real(kind = kreal), intent(in) :: thick
    real(kind = kreal), intent(in) :: params(*)
    real(kind = kreal), intent(out) :: vect(*)
    integer(kind = kint), intent(out) :: nsize

    !--------------------------------------------------------------------

    integer :: ivol, isurf
    integer :: lx, ly
    integer :: fetype
    integer :: ny
    integer :: i
    integer(kind=kint) :: ierr_quad
    integer :: na, nb
    integer :: isize
    integer :: jsize1, jsize2, jsize3, &
      jsize4, jsize5, jsize6
    integer :: ltype
    integer :: n_totlyr, n_layer

    real(kind = kreal) :: elem(3, nn)
    real(kind = kreal) :: val
    real(kind = kreal) :: ax, ay, az
    real(kind = kreal) :: rx, ry, rz
    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: w_w_lx, w_ly
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
      dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: g1_cross_g2(3)
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: w_w_w_det
    real(kind = kreal) :: N(3, ndof*nn)
    real(kind = kreal) :: hx, hy, hz
    real(kind = kreal) :: phx, phy, phz
    real(kind = kreal) :: coefx, coefy, coefz
    real(kind = kreal) :: x, y, z

    ny = 0

    !--------------------------------------------------------------------

    !   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
    !   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
    !   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
    !   CRAV LTYPE=4  :GRAVITY FORCE
    !   CENT LTYPE=5  :CENTRIFUGAL LOAD
    !   P    LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR SHELL SURFACE

    !--------------------------------------------------------------------

    nsize = ndof*nn

    !--------------------------------------------------------------------

    val = params(1)
    ax  = params(2)
    ay  = params(3)
    az  = params(4)
    rx  = params(5)
    ry  = params(6)
    rz  = params(7)

    !--------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      fetype = fe_mitc4_shell

      ny = 2

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      fetype = fe_mitc9_shell

      ny = 3

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      fetype = fe_mitc3_shell

      ny = 2

    end if

    !--------------------------------------------------------------------

    do na = 1, nn

      elem(1, na) = xx(na)
      elem(2, na) = yy(na)
      elem(3, na) = zz(na)

    end do

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    !--------------------------------------------------------------------

    ! Local load vector
    do isize = 1, ndof*nn

      vect(isize) = 0.0D0

    end do

    !--------------------------------------------------------------------

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    do i = 1, 3

      g1(i) = 0.0D0

      do na = 1, nn

        g1(i) = g1(i)+shapederiv(na, 1) &
          *elem(i, na)

      end do

    end do

    e_0(1) = g1(1)
    e_0(2) = g1(2)
    e_0(3) = g1(3)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1) &
            *elem(i, na)
          g2(i) = g2(i)+shapederiv(na, 2) &
            *elem(i, na)

        end do

      end do

      !--------------------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(1, nb) = det_cg3(1)/det_cg3_abs
      v3(2, nb) = det_cg3(2)/det_cg3_abs
      v3(3, nb) = det_cg3(3)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs .GT. 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !--------------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------

    ! Selection of load type

    ivol = 0
    isurf = 0

    if( ltype .LT. 10 ) then

      ivol = 1

    else if( ltype .GE. 10 ) then

      isurf = 1

    end if

    !--------------------------------------------------------------------

    !** Surface load
    if( isurf .EQ. 1 ) then

      !--------------------------------------------------------

      do lx = 1, NumOfQuadPoints(fetype)

        !--------------------------------------------------

        call getQuadPoint(fetype, lx, naturalcoord)

        xi_lx  = naturalcoord(1)
        eta_lx = naturalcoord(2)

        w_w_lx = getWeight(fetype, lx)

        call getShapeFunc(fetype, naturalcoord, shapefunc)

        call getShapeDeriv(fetype, naturalcoord, shapederiv)

        !--------------------------------------------------

        do na = 1, nn

          do i = 1, 3

            u_rot(i, na)                    &
              = shapefunc(na)                 &
              *( 0.0D0*a_over_2_v3(i, na) )

            dudxi_rot(i, na)                &
              = shapederiv(na, 1)             &
              *( 0.0D0*a_over_2_v3(i, na) )
            dudeta_rot(i, na)               &
              = shapederiv(na, 2)             &
              *( 0.0D0*a_over_2_v3(i, na) )
            dudzeta_rot(i, na)              &
              = shapefunc(na)                 &
              *( a_over_2_v3(i, na) )

          end do

        end do

        !--------------------------------------------------

        ! Covariant basis vector
        do i = 1, 3

          g1(i) = 0.0D0
          g2(i) = 0.0D0
          !g3(i) = 0.0D0

          do na = 1, nn

            g1(i) = g1(i)+shapederiv(na, 1)  &
              *elem(i, na)       &
              +dudxi_rot(i, na)
            g2(i) = g2(i)+shapederiv(na, 2)  &
              *elem(i, na)       &
              +dudeta_rot(i, na)
            !g3(i) = g3(i)+dudzeta_rot(i, na)

          end do

        end do

        !--------------------------------------------------

        !g3_abs = DSQRT( g3(1)*g3(1)   &
          !               +g3(2)*g3(2)   &
          !               +g3(3)*g3(3) )

        !--------------------------------------------------

        !e3_hat(1) = g3(1)/g3_abs
        !e3_hat(2) = g3(2)/g3_abs
        !e3_hat(3) = g3(3)/g3_abs

        !--------------------------------------------------

        ! Jacobian
        !det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
          !     +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
          !     +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

        !--------------------------------------------------

        g1_cross_g2(1) = g1(2)*g2(3)-g1(3)*g2(2)
        g1_cross_g2(2) = g1(3)*g2(1)-g1(1)*g2(3)
        g1_cross_g2(3) = g1(1)*g2(2)-g1(2)*g2(1)

        !--------------------------------------------------

        do nb = 1, nn

          jsize1 = ndof*(nb-1)+1
          jsize2 = ndof*(nb-1)+2
          jsize3 = ndof*(nb-1)+3
          jsize4 = ndof*(nb-1)+4
          jsize5 = ndof*(nb-1)+5
          jsize6 = ndof*(nb-1)+6

          N(1, jsize1) = shapefunc(nb)
          N(1, jsize2) = 0.0D0
          N(1, jsize3) = 0.0D0
          N(1, jsize4) = 0.0D0
          N(1, jsize5) = 0.0D0
          N(1, jsize6) = 0.0D0
          N(2, jsize1) = 0.0D0
          N(2, jsize2) = shapefunc(nb)
          N(2, jsize3) = 0.0D0
          N(2, jsize4) = 0.0D0
          N(2, jsize5) = 0.0D0
          N(2, jsize6) = 0.0D0
          N(3, jsize1) = 0.0D0
          N(3, jsize2) = 0.0D0
          N(3, jsize3) = shapefunc(nb)
          N(3, jsize4) = 0.0D0
          N(3, jsize5) = 0.0D0
          N(3, jsize6) = 0.0D0

        end do

        do isize = 1, ndof*nn

          vect(isize)                                    &
            = vect(isize)                                  &
            +w_w_lx*( N(1, isize)*g1_cross_g2(1)       &
            +N(2, isize)*g1_cross_g2(2)       &
            +N(3, isize)*g1_cross_g2(3) )*val

        end do

        !--------------------------------------------------

      end do

      !--------------------------------------------------------

    end if

    !--------------------------------------------------------------------

    !** Volume load
    if( ivol .EQ. 1 ) then

      !--------------------------------------------------------
      n_totlyr =  gausses(1)%pMaterial%totallyr
      do n_layer=1,n_totlyr
        do ly = 1, ny

          !--------------------------------------------------


          call fstr_shell_layer_quadrature_gauss( etype, gausses(1), n_layer, ly, &
            zeta_ly, w_ly, ierr_quad )
          if( ierr_quad /= 0 ) cycle

          !--------------------------------------------------

          do lx = 1, NumOfQuadPoints(fetype)

            !--------------------------------------------

            call getQuadPoint(fetype, lx, naturalcoord)

            xi_lx  = naturalcoord(1)
            eta_lx = naturalcoord(2)

            w_w_lx = getWeight(fetype, lx)

            call getShapeFunc(fetype, naturalcoord, shapefunc)

            call getShapeDeriv(fetype, naturalcoord, shapederiv)

            !--------------------------------------------

            do na = 1, nn

              do i = 1, 3

                u_rot(i, na)                      &
                  = shapefunc(na)                   &
                  *( zeta_ly*a_over_2_v3(i, na) )

                dudxi_rot(i, na)                  &
                  = shapederiv(na, 1)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudeta_rot(i, na)                 &
                  = shapederiv(na, 2)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudzeta_rot(i, na)                &
                  = shapefunc(na)                   &
                  *( a_over_2_v3(i, na) )

              end do

            end do

            !--------------------------------------------

            ! Covariant basis vector
            do i = 1, 3

              g1(i) = 0.0D0
              g2(i) = 0.0D0
              g3(i) = 0.0D0

              do na = 1, nn

                g1(i) = g1(i)+shapederiv(na, 1)  &
                  *elem(i, na)       &
                  +dudxi_rot(i, na)
                g2(i) = g2(i)+shapederiv(na, 2)  &
                  *elem(i, na)       &
                  +dudeta_rot(i, na)
                g3(i) = g3(i)+dudzeta_rot(i, na)

              end do

            end do

            !--------------------------------------------

            ! Jacobian
            det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
              +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
              +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

            !--------------------------------------------

            ! [ N ] matrix
            do nb = 1, nn

              jsize1 = ndof*(nb-1)+1
              jsize2 = ndof*(nb-1)+2
              jsize3 = ndof*(nb-1)+3
              jsize4 = ndof*(nb-1)+4
              jsize5 = ndof*(nb-1)+5
              jsize6 = ndof*(nb-1)+6

              N(1, jsize1) =  shapefunc(nb)
              N(2, jsize1) =  0.0D0
              N(3, jsize1) =  0.0D0
              N(1, jsize2) =  0.0D0
              N(2, jsize2) =  shapefunc(nb)
              N(3, jsize2) =  0.0D0
              N(1, jsize3) =  0.0D0
              N(2, jsize3) =  0.0D0
              N(3, jsize3) =  shapefunc(nb)
              N(1, jsize4) =  0.0D0
              N(2, jsize4) = -u_rot(3, nb)
              N(3, jsize4) =  u_rot(2, nb)
              N(1, jsize5) =  u_rot(3, nb)
              N(2, jsize5) =  0.0D0
              N(3, jsize5) = -u_rot(1, nb)
              N(1, jsize6) = -u_rot(2, nb)
              N(2, jsize6) =  u_rot(1, nb)
              N(3, jsize6) =  0.0D0

            enddo

            !--------------------------------------------

            w_w_w_det = w_w_lx*w_ly*det

            !--------------------------------------------

            if( ltype .EQ. 1 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*N(1, isize)*val

              end do

            else if( ltype .EQ. 2 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*N(2, isize)*val

              end do

            else if( ltype .EQ. 3 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*N(3, isize)*val

              end do

            else if( ltype .EQ. 4 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*rho*ax*N(1, isize)*val
                vect(isize) = vect(isize)+w_w_w_det*rho*ay*N(2, isize)*val
                vect(isize) = vect(isize)+w_w_w_det*rho*az*N(3, isize)*val

              end do

            else if( ltype .EQ. 5 ) then

              x = 0.0D0
              y = 0.0D0
              z = 0.0D0

              do nb = 1, nn

                x = x+shapefunc(nb)*elem(1, nb)
                y = y+shapefunc(nb)*elem(2, nb)
                z = z+shapefunc(nb)*elem(3, nb)

              end do

              hx = ax+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*rx
              hy = ay+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*ry
              hz = az+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*rz

              phx = x-hx
              phy = y-hy
              phz = z-hz

              coefx = phx*val*rho*val
              coefy = phy*val*rho*val
              coefz = phz*val*rho*val

              do isize = 1, ndof*nn

                vect(isize)                       &
                  = vect(isize)                     &
                  +w_w_w_det*( N(1, isize)*coefx   &
                  +N(2, isize)*coefy   &
                  +N(3, isize)*coefz )

              end do

            end if

            !--------------------------------------------

          end do

          !----------------------------------------------

        end do

        !----------------------------------------------

      end do

      !--------------------------------------------------------

    end if
    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine DL_Shell
  !####################################################################


  !####################################################################
  subroutine DL_Shell_33                               &
      (ic_type, nn, ndof, xx, yy, zz, rho, thick, &
      ltype, params, vect, nsize, gausses)
    !####################################################################

    use hecmw
    use m_utilities
    use mMechGauss

    type(tGaussStatus) :: gausses(:)
    !--------------------------------------------------------------------

    integer(kind = kint) :: ic_type
    integer(kind = kint) :: nn
    integer(kind = kint) :: ndof
    real(kind = kreal) :: xx(*), yy(*), zz(*)
    real(kind = kreal) :: rho
    real(kind = kreal) :: thick
    real(kind = kreal) :: params(*)
    real(kind = kreal) :: vect(*)
    integer(kind = kint) :: nsize
    integer :: ltype, i
    real(kind = kreal) :: tmp(24)
    !--------------------------------------------------------------------

    if(ic_type == 761)then
      !ic_type = 731
      !nn = 3
      !ndof = 6
      call DL_Shell(731, 3, 6, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
      !ic_type = 761
      !nn = 6
      !ndof = 3

      tmp = 0.0
      do i=1,18
        tmp(i) = vect(i)
      enddo

      vect( 1) = tmp(1)
      vect( 2) = tmp(2)
      vect( 3) = tmp(3)
      vect( 4) = tmp(7)
      vect( 5) = tmp(8)
      vect( 6) = tmp(9)
      vect( 7) = tmp(13)
      vect( 8) = tmp(14)
      vect( 9) = tmp(15)
      vect(10) = tmp(4)
      vect(11) = tmp(5)
      vect(12) = tmp(6)
      vect(13) = tmp(10)
      vect(14) = tmp(11)
      vect(15) = tmp(12)
      vect(16) = tmp(16)
      vect(17) = tmp(17)
      vect(18) = tmp(18)

    elseif(ic_type == 781)then
      !ic_type = 741
      !nn = 4
      !ndof = 6
      call DL_Shell(741, 4, 6, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
      !ic_type = 781
      !nn = 8
      !ndof = 3

      tmp = 0.0
      do i=1,24
        tmp(i) = vect(i)
      enddo

      vect( 1) = tmp(1)
      vect( 2) = tmp(2)
      vect( 3) = tmp(3)
      vect( 4) = tmp(7)
      vect( 5) = tmp(8)
      vect( 6) = tmp(9)
      vect( 7) = tmp(13)
      vect( 8) = tmp(14)
      vect( 9) = tmp(15)
      vect(10) = tmp(19)
      vect(11) = tmp(20)
      vect(12) = tmp(21)
      vect(13) = tmp(4)
      vect(14) = tmp(5)
      vect(15) = tmp(6)
      vect(16) = tmp(10)
      vect(17) = tmp(11)
      vect(18) = tmp(12)
      vect(19) = tmp(16)
      vect(20) = tmp(17)
      vect(21) = tmp(18)
      vect(22) = tmp(22)
      vect(23) = tmp(23)
      vect(24) = tmp(24)

    endif

  end subroutine DL_Shell_33

  !####################################################################
  ! Update shell stress and equivalent nodal force.
  subroutine UpdateST_Shell_MITC                                           &
      (etype, nn, ndof, ecoord, u, du, gausses, qf, thick, mixflag, nddisp, element, nddirector, ndrefdirector, &
      ndcurdirector, nddrill)
    !####################################################################

    use mMechGauss
    use m_MatMatrix
    use mMaterial, only: TOTALLAG, UPDATELAG, isElastic

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn, mixflag
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    real(kind = kreal), intent(in)   :: u(:, :)
    real(kind = kreal), intent(in)   :: du(:, :)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(out)  :: qf(:)
    real(kind = kreal), intent(in)   :: thick

    real(kind = kreal), intent(in), optional :: nddisp(ndof, nn)
    type(tElement), intent(inout), optional :: element
    real(kind = kreal), intent(in), optional :: nddirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndrefdirector(3, nn)
    real(kind = kreal), intent(in), optional :: ndcurdirector(3, nn)
    real(kind = kreal), intent(in), optional :: nddrill(nn)
    !--------------------------------------------------------------------

    real(kind = kreal)   :: stiff(nn*ndof, nn*ndof), totaldisp(nn*ndof), edisp(6, nn), qf_direct(nn*ndof)
    real(kind = kreal)   :: incdisp(6, nn), basedisp(6, nn)
    integer(kind = kint) :: i
    integer(kind = kint) :: flag
    logical :: use_stress_force

    flag = gausses(1)%pMaterial%nlgeom_flag
    use_stress_force = present( element ) .and. &
      etype == fe_mitc4_shell .and. nn == 4 .and. &
      ( flag == TOTALLAG .or. flag == UPDATELAG ) .and. &
      isElastic( gausses(1)%pMaterial%mtype )

    totaldisp = 0.d0
    edisp(:, :) = 0.0d0
    incdisp(:, :) = 0.0d0
    basedisp(:, :) = 0.0d0
    do i=1,nn
      totaldisp(ndof*(i-1)+1:ndof*i) = u(1:ndof,i) + du(1:ndof,i)
      basedisp(1:6, i) = u(1:6, i)
      incdisp(1:6, i) = du(1:6, i)
      if( use_stress_force ) then
        edisp(1:3, i) = u(1:3, i) + du(1:3, i)
        call ShellComposeRotationVector( u(4:6, i), du(4:6, i), edisp(4:6, i) )
      else
        edisp(1:6, i) = u(1:6, i) + du(1:6, i)
      endif
    end do

    ! Equivalent nodal force from layer stress.
    qf_direct(:) = 0.0d0

    if( use_stress_force ) then
      if( flag == UPDATELAG ) then
        if( present( ndcurdirector ) ) then
          call UpdateShellLayerGauss_Shell_MITC( &
            etype, nn, ndof, ecoord, element, incdisp, thick, nddirector=nddirector, &
            ndrefdirector=ndcurdirector, ndbase_disp=basedisp )
        else
          call UpdateShellLayerGauss_Shell_MITC( &
            etype, nn, ndof, ecoord, element, incdisp, thick, nddirector=nddirector, &
            ndrefdirector=ndrefdirector, ndbase_disp=basedisp )
        endif
      else
        call UpdateShellLayerGauss_Shell_MITC( &
          etype, nn, ndof, ecoord, element, edisp, thick, nddirector=nddirector, ndrefdirector=ndrefdirector )
      endif
      if( present( nddisp ) ) then
        if( flag == UPDATELAG .and. present( ndcurdirector ) ) then
          call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, &
            mixflag, nddisp=nddisp, element=element, qf_stress=qf_direct, &
            include_geo_stiff=.false., nddirector=nddirector, ndrefdirector=ndcurdirector, nddrill=nddrill)
        else
          call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, &
            mixflag, nddisp=nddisp, element=element, qf_stress=qf_direct, &
            include_geo_stiff=.false., nddirector=nddirector, ndrefdirector=ndrefdirector, nddrill=nddrill)
        endif
      else
        if( flag == UPDATELAG .and. present( ndcurdirector ) ) then
          call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, &
            mixflag, nddisp=edisp, element=element, qf_stress=qf_direct, &
            include_geo_stiff=.false., nddirector=nddirector, ndrefdirector=ndcurdirector, nddrill=nddrill)
        else
          call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, &
            mixflag, nddisp=edisp, element=element, qf_stress=qf_direct, &
            include_geo_stiff=.false., nddirector=nddirector, ndrefdirector=ndrefdirector, nddrill=nddrill)
        endif
      endif
    else
      if( present( nddisp ) ) then
        call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, &
          mixflag, nddisp=nddisp, include_geo_stiff=.false., nddirector=nddirector, ndrefdirector=ndrefdirector, &
          nddrill=nddrill)
      else
        call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag)
      endif
    endif

    if( use_stress_force ) then
      qf = qf_direct
    else
      qf = matmul(stiff,totaldisp)
    endif

  end subroutine UpdateST_Shell_MITC

  !####################################################################
  ! this subroutine can be used only for linear analysis
  subroutine UpdateST_Shell_MITC33                                         &
      (etype, nn, ndof, ecoord, u, du, gausses, qf, thick, mixflag, nddisp)
    !####################################################################

    use mMechGauss
    use m_MatMatrix

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn, mixflag
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    real(kind = kreal), intent(in)   :: u(3, nn*2)
    real(kind = kreal), intent(in)   :: du(3, nn*2)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(out)  :: qf(:)
    real(kind = kreal), intent(in)   :: thick

    real(kind = kreal), intent(in), optional :: nddisp(3, nn)
    !--------------------------------------------------------------------

    real(kind = kreal)   :: stiff(nn*ndof, nn*ndof), totaldisp(nn*ndof)
    integer(kind = kint) :: i

    call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, nddisp)

    totaldisp = 0.d0
    do i=1,nn
      totaldisp(ndof*(i-1)+1:ndof*(i-1)+3) = u(1:3,2*i-1) + du(1:3,2*i-1)
      totaldisp(ndof*(i-1)+4:ndof*(i-1)+6) = u(1:3,2*i) + du(1:3,2*i)
    end do

    qf = matmul(stiff,totaldisp)

  end subroutine UpdateST_Shell_MITC33
  
  !####################################################################
  subroutine mass_Shell(etype, nn, elem, rho, thick, gausses, mass, lumped)
  !####################################################################
    use hecmw
    use m_utilities
    use mMechGauss
    use Quadrature
    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    real(kind = kreal), intent(in)   :: elem(3,nn)
    real(kind = kreal), intent(in)   :: rho
    real(kind = kreal), intent(in)   :: thick
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind=kreal), intent(out)    :: mass(:,:)      !< mass matrix
    real(kind=kreal), intent(out)    :: lumped(:)      !< lumped mass matrix

    !--------------------------------------------------------------------

    integer :: lx, ly, nsize, ndof
    integer :: fetype
    integer :: ny
    integer :: i
    integer(kind=kint) :: ierr_quad
    integer :: na, nb
    integer :: jsize1, jsize2, jsize3, jsize4, jsize5, jsize6
    integer :: n_totlyr, n_layer

    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: w_w_lx, w_ly
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: w_w_w_det
    real(kind = kreal) :: N(3, 6*nn)
    real(kind = kreal) :: totalmass, totdiag

    ny = 0; ndof=6
    nsize = ndof*nn

    !--------------------------------------------------------------------

    ! MITC4
    if( etype == fe_mitc4_shell ) then
      fetype = fe_mitc4_shell
      ny = 2

      ! MITC9
    else if( etype == fe_mitc9_shell ) then
      fetype = fe_mitc9_shell
      ny = 3

      ! MITC3
    else if( etype == fe_mitc3_shell ) then
      fetype = fe_mitc3_shell
      ny = 2

    end if

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    g1(:) = matmul( elem, shapederiv(:,1) )

    e_0(:) = g1(:)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)
      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      g1(:) = matmul( elem, shapederiv(:,1) )
      g2(:) = matmul( elem, shapederiv(:,2) )

      !--------------------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(:, nb) = det_cg3(:)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs > 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !--------------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------

    mass(:,:) = 0.0D0
    totalmass  = 0.d0
    n_totlyr =  gausses(1)%pMaterial%totallyr
    do n_layer=1,n_totlyr
        do ly = 1, ny

          !--------------------------------------------------

          call fstr_shell_layer_quadrature_gauss( etype, gausses(1), n_layer, ly, &
            zeta_ly, w_ly, ierr_quad )
          if( ierr_quad /= 0 ) cycle

          !--------------------------------------------------

          do lx = 1, NumOfQuadPoints(fetype)

            !--------------------------------------------

            call getQuadPoint(fetype, lx, naturalcoord)

            xi_lx  = naturalcoord(1)
            eta_lx = naturalcoord(2)

            w_w_lx = getWeight(fetype, lx)

            call getShapeFunc(fetype, naturalcoord, shapefunc)
            call getShapeDeriv(fetype, naturalcoord, shapederiv)

            !--------------------------------------------

            do na = 1, nn

              do i = 1, 3

                u_rot(i, na) = shapefunc(na)*( zeta_ly*a_over_2_v3(i, na) )

                dudxi_rot(i, na) = shapederiv(na, 1)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudeta_rot(i, na) = shapederiv(na, 2)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudzeta_rot(i, na) = shapefunc(na)                   &
                  *( a_over_2_v3(i, na) )

              end do

            end do

            !--------------------------------------------

            ! Covariant basis vector
            do i = 1, 3
              g1(i) = 0.0D0
              g2(i) = 0.0D0
              g3(i) = 0.0D0
              do na = 1, nn
                g1(i) = g1(i)+shapederiv(na, 1) *elem(i, na)       &
                  +dudxi_rot(i, na)
                g2(i) = g2(i)+shapederiv(na, 2) *elem(i, na)       &
                  +dudeta_rot(i, na)
                g3(i) = g3(i)+dudzeta_rot(i, na)
              end do
            end do

            !--------------------------------------------

            ! Jacobian
            det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
              +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
              +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

            !--------------------------------------------

            ! [ N ] matrix
            do nb = 1, nn

              jsize1 = ndof*(nb-1)+1
              jsize2 = ndof*(nb-1)+2
              jsize3 = ndof*(nb-1)+3
              jsize4 = ndof*(nb-1)+4
              jsize5 = ndof*(nb-1)+5
              jsize6 = ndof*(nb-1)+6

              N(1, jsize1) =  shapefunc(nb)
              N(2, jsize1) =  0.0D0
              N(3, jsize1) =  0.0D0
              N(1, jsize2) =  0.0D0
              N(2, jsize2) =  shapefunc(nb)
              N(3, jsize2) =  0.0D0
              N(1, jsize3) =  0.0D0
              N(2, jsize3) =  0.0D0
              N(3, jsize3) =  shapefunc(nb)
              N(1, jsize4) =  0.0D0
              N(2, jsize4) = -u_rot(3, nb)
              N(3, jsize4) =  u_rot(2, nb)
              N(1, jsize5) =  u_rot(3, nb)
              N(2, jsize5) =  0.0D0
              N(3, jsize5) = -u_rot(1, nb)
              N(1, jsize6) = -u_rot(2, nb)
              N(2, jsize6) =  u_rot(1, nb)
              N(3, jsize6) =  0.0D0

            enddo

            !--------------------------------------------

            w_w_w_det = w_w_lx*w_ly*det*gausses(1)%pMaterial%shell_var(n_layer)%weight
            mass(1:nsize,1:nsize) = mass(1:nsize,1:nsize)+ matmul( transpose(N), N )*w_w_w_det*rho
            totalmass = totalmass + w_w_w_det*rho
            !--------------------------------------------

          end do

          !----------------------------------------------

        end do

        !----------------------------------------------

    end do
    totalmass = totalmass*3.d0

    totdiag=0.d0
    do nb = 1, nn
        DO i = 1, 3
          lx = (nb-1)*ndof+i
          totdiag = totdiag + mass(lx,lx)
        END DO
    ENDDO
    do nb = 1, nn
        DO i = 1, 6
          lx = (nb-1)*ndof+i
          lumped(lx) = mass(lx,lx)/totdiag* totalmass
        END DO
    ENDDO

  !####################################################################
  end subroutine mass_Shell
  !####################################################################

end module m_static_LIB_shell
