!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module ...
module m_static_LIB_shell
  use hecmw, only : kint, kreal, hecmw_abort, hecmw_comm_get_comm
  use elementInfo
  use MITC_Tying, only: NumOfTyingSets, NumOfTyingPoints, getTyingPoint, mitc9_xi_sign, mitc9_eta_sign
  use m_utilities, only: cross_product
  use m_fstr_FiniteRotationKinematics, only: fstr_is_finite_rotation_shell_element, &
    ShellRotationVectorToMatrix, ShellComposeRotationVector, ShellSkewMatrix

  implicit none

  private

  !> Upper bound on the number of terms in the MITC assumed-strain (tying) operator.
  !> MITC9 is the largest case: 6 tying points x 5 target components.
  integer(kind=kint), parameter :: MAX_TYING_TERMS = 64
  !> Natural-coordinate direction indices used in grouped basis arrays.
  integer(kind=kint), parameter :: SHELL_XI = 1_kint
  integer(kind=kint), parameter :: SHELL_ETA = 2_kint
  integer(kind=kint), parameter :: SHELL_ZETA = 3_kint


  public :: STF_Shell_MITC
  public :: ElementStress_Shell_MITC
  public :: DL_Shell
  public :: DL_Shell_33
  public :: UPDATE_Shell_MITC
  public :: UPDATE_Shell_MITC33
  public :: mass_Shell

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
  ! [1] Noguchi, H. and Hisada, T., "Sensitivity analysis in post-buckling problems of shell structures,"
  !     Computers & Structures, Vol.47, No.4, pp.699-710, (1993).
  ! [2] Dvorkin, E.N. and Bathe, K.J., "A Continuum Mechanics Based Four-node Shell Element for General Non-linear Analysis,"
  !     Engineering Computations, Vol.1, pp.77-88, (1984).
  ! [3] Bucalem, M.L. and Bathe, K.J., "Higher-order MITC general shell element,"
  !     International Journal for Numerical Methods in Engineering, Vol.36, pp.3729-3754, (1993).
  ! [4] Lee, P.S. and Bathe, K.J., "Development of MITC Isotropic Triangular Shell Finite Elements,"
  !     Computers & Structures, Vol.82, pp.945-962, (2004).
  !
  ! Xi YUAN
  !   Apr. 13, 2019: Introduce mass matrix calculation
  ! (Ref.)
  ! [5] E. Hinton, T. A. Rock and O. C. Zienkiewicz, "A Note on Mass Lumping and Related Processes in FEM,"
  !     Earthquake Engineering & Structural Dynamics, Vol.4, pp.245-249, (1976).
  !
  !--------------------------------------------------------------

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

  !--------------------------------------------------------------------
  !> Pack the symmetric part of a 3x3 matrix into the five shell strain components.
  pure function ShellMITC_SymmetricComponents(matrix) result(components)
    real(kind=kreal), intent(in) :: matrix(3, 3)
    real(kind=kreal) :: components(5)

    components = (/ matrix(1, 1), matrix(2, 2), matrix(1, 2)+matrix(2, 1), matrix(2, 3)+matrix(3, 2), &
      matrix(3, 1)+matrix(1, 3) /)
  end function ShellMITC_SymmetricComponents

  !--------------------------------------------------------------------
  !> Determinant of the covariant shell basis.
  pure real(kind=kreal) function ShellMITC_CovariantJacobian(covariant_basis)
    real(kind=kreal), intent(in) :: covariant_basis(3, 3)

    ShellMITC_CovariantJacobian = covariant_basis(1, SHELL_XI) &
      *(covariant_basis(2, SHELL_ETA)*covariant_basis(3, SHELL_ZETA) &
      -covariant_basis(3, SHELL_ETA)*covariant_basis(2, SHELL_ZETA)) +covariant_basis(2, SHELL_XI) &
      *(covariant_basis(3, SHELL_ETA)*covariant_basis(1, SHELL_ZETA) &
      -covariant_basis(1, SHELL_ETA)*covariant_basis(3, SHELL_ZETA)) +covariant_basis(3, SHELL_XI) &
      *(covariant_basis(1, SHELL_ETA)*covariant_basis(2, SHELL_ZETA) &
      -covariant_basis(2, SHELL_ETA)*covariant_basis(1, SHELL_ZETA))
  end function ShellMITC_CovariantJacobian

  pure real(kind=kreal) function ShellMITC_TyingZeta(etype, zeta)
    integer(kind=kint), intent(in) :: etype
    real(kind=kreal), intent(in) :: zeta

    ShellMITC_TyingZeta = 0.0D0
    if( etype == fe_mitc9_shell ) ShellMITC_TyingZeta = zeta
  end function ShellMITC_TyingZeta

  !--------------------------------------------------------------------
  !> Resolve the shell formulation once at an entry point.
  !> Results are returned as plain scalars; no persistent formulation object is created.
  subroutine ShellMITC_ResolveFormulation(etype, nn, ndof, gauss, has_nodal_state, &
      has_element_state, require_layer_state, kinematics, ndof_shell, finite_rotation, &
      use_director_tangent, use_green_lagrange, add_geometric_stiffness, update_state)
    use mMechGauss, only: tGaussStatus
    use mMaterial, only: INFINITESIMAL, TOTALLAG, UPDATELAG, isElastic
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof
    type(tGaussStatus), intent(in) :: gauss
    logical, intent(in) :: has_nodal_state, has_element_state, require_layer_state
    integer(kind=kint), intent(out) :: kinematics, ndof_shell
    logical, intent(out) :: finite_rotation, use_director_tangent, use_green_lagrange
    logical, intent(out) :: add_geometric_stiffness, update_state

    kinematics = gauss%pMaterial%nlgeom_flag
    if( .not. has_nodal_state ) kinematics = INFINITESIMAL
    ndof_shell = min(ndof, 6_kint)
    finite_rotation = (kinematics == TOTALLAG .or. kinematics == UPDATELAG) &
      .and. ndof_shell >= 6 .and. fstr_is_finite_rotation_shell_element(etype, nn) .and. isElastic(gauss%pMaterial%mtype)
    use_director_tangent = finite_rotation
    use_green_lagrange = kinematics == TOTALLAG .and. finite_rotation
    add_geometric_stiffness = kinematics /= INFINITESIMAL
    update_state = finite_rotation .and. has_element_state

    if( kinematics /= INFINITESIMAL ) then
      if( .not. finite_rotation ) call ShellMITC_AbortNonlinearUnsupported(etype)
      if( require_layer_state .and. .not. update_state ) then
        call ShellMITC_AbortNonlinearUnsupported(etype)
      endif
    endif
  end subroutine ShellMITC_ResolveFormulation

  !--------------------------------------------------------------------
  !> Prepare the evaluation coordinates and nodal directors shared by STF/UPDATE.
  subroutine ShellMITC_PrepareNodalKinematics(etype, nn, thick, kinematics, finite_rotation, &
      use_director_tangent, need_second_tangent, ecoord, nodal_state, ndtriad, ndreftriad, &
      ndcurtriad, evaluation_coords, v1, v2, v3, director, reference_director, director_tangent, director_second_tangent)
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, kinematics
    real(kind=kreal), intent(in) :: thick, ecoord(3, nn), nodal_state(6, nn)
    logical, intent(in) :: finite_rotation, use_director_tangent, need_second_tangent
    real(kind=kreal), intent(in), optional :: ndtriad(9, nn), ndreftriad(9, nn), ndcurtriad(9, nn)
    real(kind=kreal), intent(out) :: evaluation_coords(3, nn)
    real(kind=kreal), intent(out) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal), intent(out) :: director(3, nn), reference_director(3, nn)
    real(kind=kreal), intent(out) :: director_tangent(3, 3, nn)
    real(kind=kreal), intent(out) :: director_second_tangent(3, 3, 3, nn)

    evaluation_coords = ecoord
    if( kinematics == UPDATELAG ) evaluation_coords = evaluation_coords+nodal_state(1:3, :)
    call ShellMITC_SetupNodalDirectors(etype, nn, thick, kinematics, evaluation_coords, &
      nodal_state, finite_rotation, use_director_tangent, need_second_tangent, ndtriad, &
      ndreftriad, ndcurtriad, v1, v2, v3, director, reference_director, director_tangent, director_second_tangent)
  end subroutine ShellMITC_PrepareNodalKinematics

  !--------------------------------------------------------------------
  !> Evaluate shell covariant basis vectors at one point.
  subroutine ShellMITC_CovariantBasis(nn, coords, director, zeta, shapefunc, shapederiv, covariant_basis)
    implicit none

    integer(kind=kint), intent(in) :: nn
    real(kind=kreal), intent(in) :: coords(3, nn), director(3, nn)
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(out) :: covariant_basis(3, 3)

    covariant_basis(:, SHELL_XI:SHELL_ETA) = matmul(coords+zeta*director, shapederiv)
    covariant_basis(:, SHELL_ZETA) = matmul(director, shapefunc)

  end subroutine ShellMITC_CovariantBasis

  !--------------------------------------------------------------------
  !> Build the nodal director contributions to the point basis variation.
  subroutine ShellMITC_DirectorContributions(nn, director, zeta, shapefunc, shapederiv, director_contribution)
    implicit none

    integer(kind=kint), intent(in) :: nn
    real(kind=kreal), intent(in) :: director(3, nn)
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(out) :: director_contribution(3, nn, 3)

    integer :: node

    do node = 1, nn
      director_contribution(:, node, SHELL_XI) = zeta*shapederiv(node, SHELL_XI)*director(:, node)
      director_contribution(:, node, SHELL_ETA) = zeta*shapederiv(node, SHELL_ETA)*director(:, node)
      director_contribution(:, node, SHELL_ZETA) = shapefunc(node)*director(:, node)
    end do
  end subroutine ShellMITC_DirectorContributions
  pure subroutine ShellMITC_BasisFromCovariant(covariant_basis, local_basis, reciprocal_basis, det)
    real(kind=kreal), intent(in) :: covariant_basis(3, 3)
    real(kind=kreal), intent(out) :: local_basis(3, 3), reciprocal_basis(3, 3)
    real(kind=kreal), intent(out) :: det

    real(kind=kreal) :: det_inv, basis_norm

    det = ShellMITC_CovariantJacobian(covariant_basis)
    det_inv = 1.0D0/det

    reciprocal_basis(1, SHELL_XI) = det_inv *(covariant_basis(2, SHELL_ETA)*covariant_basis(3, SHELL_ZETA) &
      -covariant_basis(3, SHELL_ETA)*covariant_basis(2, SHELL_ZETA))
    reciprocal_basis(2, SHELL_XI) = det_inv *(covariant_basis(3, SHELL_ETA)*covariant_basis(1, SHELL_ZETA) &
      -covariant_basis(1, SHELL_ETA)*covariant_basis(3, SHELL_ZETA))
    reciprocal_basis(3, SHELL_XI) = det_inv *(covariant_basis(1, SHELL_ETA)*covariant_basis(2, SHELL_ZETA) &
      -covariant_basis(2, SHELL_ETA)*covariant_basis(1, SHELL_ZETA))
    reciprocal_basis(1, SHELL_ETA) = det_inv *(covariant_basis(2, SHELL_ZETA)*covariant_basis(3, SHELL_XI) &
      -covariant_basis(3, SHELL_ZETA)*covariant_basis(2, SHELL_XI))
    reciprocal_basis(2, SHELL_ETA) = det_inv *(covariant_basis(3, SHELL_ZETA)*covariant_basis(1, SHELL_XI) &
      -covariant_basis(1, SHELL_ZETA)*covariant_basis(3, SHELL_XI))
    reciprocal_basis(3, SHELL_ETA) = det_inv *(covariant_basis(1, SHELL_ZETA)*covariant_basis(2, SHELL_XI) &
      -covariant_basis(2, SHELL_ZETA)*covariant_basis(1, SHELL_XI))
    reciprocal_basis(1, SHELL_ZETA) = det_inv *(covariant_basis(2, SHELL_XI)*covariant_basis(3, SHELL_ETA) &
      -covariant_basis(3, SHELL_XI)*covariant_basis(2, SHELL_ETA))
    reciprocal_basis(2, SHELL_ZETA) = det_inv *(covariant_basis(3, SHELL_XI)*covariant_basis(1, SHELL_ETA) &
      -covariant_basis(1, SHELL_XI)*covariant_basis(3, SHELL_ETA))
    reciprocal_basis(3, SHELL_ZETA) = det_inv *(covariant_basis(1, SHELL_XI)*covariant_basis(2, SHELL_ETA) &
      -covariant_basis(2, SHELL_XI)*covariant_basis(1, SHELL_ETA))

    basis_norm = dsqrt(dot_product(covariant_basis(:, SHELL_ZETA), covariant_basis(:, SHELL_ZETA)))
    local_basis(:, SHELL_ZETA) = covariant_basis(:, SHELL_ZETA)/basis_norm
    local_basis(1, SHELL_XI) = covariant_basis(2, SHELL_ETA)*local_basis(3, SHELL_ZETA) &
      -covariant_basis(3, SHELL_ETA)*local_basis(2, SHELL_ZETA)
    local_basis(2, SHELL_XI) = covariant_basis(3, SHELL_ETA)*local_basis(1, SHELL_ZETA) &
      -covariant_basis(1, SHELL_ETA)*local_basis(3, SHELL_ZETA)
    local_basis(3, SHELL_XI) = covariant_basis(1, SHELL_ETA)*local_basis(2, SHELL_ZETA) &
      -covariant_basis(2, SHELL_ETA)*local_basis(1, SHELL_ZETA)
    basis_norm = dsqrt(dot_product(local_basis(:, SHELL_XI), local_basis(:, SHELL_XI)))
    local_basis(:, SHELL_XI) = local_basis(:, SHELL_XI)/basis_norm
    local_basis(1, SHELL_ETA) = local_basis(2, SHELL_ZETA)*local_basis(3, SHELL_XI) &
      -local_basis(3, SHELL_ZETA)*local_basis(2, SHELL_XI)
    local_basis(2, SHELL_ETA) = local_basis(3, SHELL_ZETA)*local_basis(1, SHELL_XI) &
      -local_basis(1, SHELL_ZETA)*local_basis(3, SHELL_XI)
    local_basis(3, SHELL_ETA) = local_basis(1, SHELL_ZETA)*local_basis(2, SHELL_XI) &
      -local_basis(2, SHELL_ZETA)*local_basis(1, SHELL_XI)
    basis_norm = dsqrt(dot_product(local_basis(:, SHELL_ETA), local_basis(:, SHELL_ETA)))
    local_basis(:, SHELL_ETA) = local_basis(:, SHELL_ETA)/basis_norm
  end subroutine ShellMITC_BasisFromCovariant

  !--------------------------------------------------------------------
  !> Build point geometry, reference/material frame, and tangent bases once.
  !> All values are returned through ordinary arrays; no point-state object is retained.
  !> Second variations remain outside this routine and are evaluated by STF only.
  subroutine ShellMITC_PreparePointKinematics(etype, nn, use_green_lagrange, reference_coords, &
      evaluation_coords, translation, director, reference_director, zeta, shapefunc, &
      shapederiv, covariant_basis, tangent_basis, reciprocal_basis, local_basis, &
      material_reciprocal_basis, material_local_basis, integration_jacobian, director_contribution)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn
    logical, intent(in) :: use_green_lagrange
    real(kind=kreal), intent(in) :: reference_coords(3, nn), evaluation_coords(3, nn)
    real(kind=kreal), intent(in) :: translation(3, nn), director(3, nn)
    real(kind=kreal), intent(in) :: reference_director(3, nn)
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(out) :: covariant_basis(3, 3), tangent_basis(3, 2)
    real(kind=kreal), intent(out) :: reciprocal_basis(3, 3), local_basis(3, 3)
    real(kind=kreal), intent(out) :: material_reciprocal_basis(3, 3)
    real(kind=kreal), intent(out) :: material_local_basis(3, 3)
    real(kind=kreal), intent(out) :: integration_jacobian
    real(kind=kreal), intent(out) :: director_contribution(3, nn, 3)

    real(kind=kreal) :: reference_basis(3, 3), translation_gradient(3, 2)
    real(kind=kreal) :: current_jacobian, material_jacobian

    call ShellMITC_DirectorContributions(nn, director, zeta, shapefunc, shapederiv, director_contribution)
    call ShellMITC_CovariantBasis(nn, evaluation_coords, director, zeta, shapefunc, shapederiv, covariant_basis)
    if( abs(ShellMITC_CovariantJacobian(covariant_basis)) <= tiny(1.0D0) ) &
      stop "Invalid shell Jacobian"
    if( use_green_lagrange ) then
      translation_gradient = matmul(translation, shapederiv)
      covariant_basis(:, SHELL_XI:SHELL_ETA) = covariant_basis(:, SHELL_XI:SHELL_ETA)+translation_gradient
    endif
    call ShellMITC_BasisFromCovariant(covariant_basis, local_basis, reciprocal_basis, current_jacobian)

    tangent_basis = covariant_basis(:, SHELL_XI:SHELL_ETA)
    material_local_basis = local_basis
    material_reciprocal_basis = reciprocal_basis
    integration_jacobian = current_jacobian

    if( .not. use_green_lagrange ) return

    call ShellMITC_CovariantBasis(nn, reference_coords, reference_director, zeta, shapefunc, shapederiv, reference_basis)
    call ShellMITC_BasisFromCovariant(reference_basis, material_local_basis, material_reciprocal_basis, material_jacobian)
    integration_jacobian = ShellMITC_CovariantJacobian(reference_basis)

    ! Preserve the existing element-specific tangent basis construction.
    if( etype == fe_mitc9_shell ) then
      tangent_basis = tangent_basis + translation_gradient
    else if( etype /= fe_mitc4_shell ) then
      covariant_basis(:, SHELL_XI:SHELL_ETA) = covariant_basis(:, SHELL_XI:SHELL_ETA) + translation_gradient
      tangent_basis = covariant_basis(:, SHELL_XI:SHELL_ETA)
    endif
  end subroutine ShellMITC_PreparePointKinematics

  !--------------------------------------------------------------------
  !> Evaluate the ordinary shell strain at one point before MITC interpolation.
  subroutine ShellMITC_EvaluatePointStrain(nn, zeta, shapefunc, shapederiv, translation, &
      director_increment, covariant_basis, use_green_lagrange, strain, reference_basis, &
      current_basis, reference_jacobian, current_jacobian)
    implicit none

    integer(kind=kint), intent(in) :: nn
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(in) :: translation(3, nn), director_increment(3, nn)
    real(kind=kreal), intent(in) :: covariant_basis(3, 3)
    logical, intent(in) :: use_green_lagrange
    real(kind=kreal), intent(out) :: strain(5)
    real(kind=kreal), intent(out) :: reference_basis(3, 3), current_basis(3, 3)
    real(kind=kreal), intent(out) :: reference_jacobian, current_jacobian

    real(kind=kreal) :: displacement_gradient(3, 3), translation_gradient(3, 2)

    displacement_gradient(:, SHELL_XI:SHELL_ETA) = matmul(translation+zeta*director_increment, shapederiv)
    displacement_gradient(:, SHELL_ZETA) = matmul(director_increment, shapefunc)

    reference_basis = covariant_basis
    current_basis = covariant_basis

    if( use_green_lagrange ) then
      translation_gradient = matmul(translation, shapederiv)
      current_basis(:, SHELL_XI:SHELL_ETA) = covariant_basis(:, SHELL_XI:SHELL_ETA)+translation_gradient
      reference_basis = current_basis-displacement_gradient

      strain = 0.5D0*ShellMITC_SymmetricComponents( matmul(transpose(current_basis), current_basis) &
        -matmul(transpose(reference_basis), reference_basis))
    else
      strain = ShellMITC_SymmetricComponents( matmul(transpose(covariant_basis), displacement_gradient))
    endif

    reference_jacobian = ShellMITC_CovariantJacobian(reference_basis)
    current_jacobian = ShellMITC_CovariantJacobian(current_basis)

  end subroutine ShellMITC_EvaluatePointStrain

  !--------------------------------------------------------------------
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
  !> Construct the reference nodal triads and half-thickness directors.
  subroutine ShellMITC_SetupReferenceDirectors(etype, nn, thick, elem, v1, v2, v3, director)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in) :: thick, elem(3, nn)
    real(kind=kreal), intent(out) :: v1(3, nn), v2(3, nn), v3(3, nn), director(3, nn)

    integer :: nb
    real(kind=kreal) :: naturalcoord(2), nncoord(nn, 2), shapederiv(nn, 2)
    real(kind=kreal) :: g1(3), g2(3), e0(3), normv

    call getNodalNaturalCoord(etype, nncoord)
    naturalcoord = 0.0D0
    call getShapeDeriv(etype, naturalcoord, shapederiv)
    e0 = matmul(elem, shapederiv(:, 1))

    do nb = 1, nn
      naturalcoord = nncoord(nb, :)
      call getShapeDeriv(etype, naturalcoord, shapederiv)
      g1 = matmul(elem, shapederiv(:, 1))
      g2 = matmul(elem, shapederiv(:, 2))

      call cross_product(g1, g2, v3(:, nb))
      normv = dsqrt(dot_product(v3(:, nb), v3(:, nb)))
      v3(:, nb) = v3(:, nb)/normv

      call cross_product(v3(:, nb), e0, v2(:, nb))
      normv = dsqrt(dot_product(v2(:, nb), v2(:, nb)))
      if (normv > 1.0D-15) then
        v2(:, nb) = v2(:, nb)/normv
        call cross_product(v2(:, nb), v3(:, nb), v1(:, nb))
        normv = dsqrt(dot_product(v1(:, nb), v1(:, nb)))
        v1(:, nb) = v1(:, nb)/normv
      else
        v1(:, nb) = (/ 0.0D0, 0.0D0, -1.0D0 /)
        v2(:, nb) = (/ 0.0D0, 1.0D0, 0.0D0 /)
      endif

      call cross_product(v1(:, nb), v2(:, nb), v3(:, nb))
      normv = dsqrt(dot_product(v3(:, nb), v3(:, nb)))
      v3(:, nb) = v3(:, nb)/normv
      director(:, nb) = 0.5D0*thick*v3(:, nb)
    end do
  end subroutine ShellMITC_SetupReferenceDirectors

  !--------------------------------------------------------------------
  !> Build the shell displacement interpolation matrix at one integration point.
  subroutine ShellMITC_BuildInterpolationMatrix(nn, ndof, zeta, shapefunc, director, N)
    implicit none

    integer(kind=kint), intent(in) :: nn, ndof
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), director(3, nn)
    real(kind=kreal), intent(out) :: N(3, ndof*nn)

    integer :: nb, jsize
    real(kind=kreal) :: urot(3)

    N = 0.0D0
    do nb = 1, nn
      jsize = ndof*(nb-1)
      N(1, jsize+1) = shapefunc(nb)
      N(2, jsize+2) = shapefunc(nb)
      N(3, jsize+3) = shapefunc(nb)
      if (ndof >= 6) then
        urot = zeta*shapefunc(nb)*director(:, nb)
        N(:, jsize+4:jsize+6) = ShellSkewMatrix(urot)
      endif
    end do
  end subroutine ShellMITC_BuildInterpolationMatrix


  !--------------------------------------------------------------------
  !> Construct nodal shell frames and the current/reference directors.
  !> Nodal frames (triads) are packed as e1(1:3), e2(4:6), e3=director axis(7:9).
  !> The director itself (0.5*thick*e3) is computed here, once, for every caller.
  subroutine ShellMITC_SetupNodalDirectors(etype, nn, thick, flag, elem, shell_disp, &
      finite_rotation_director, use_director_tangent, need_second_tangent, ndtriad, ndreftriad, ndcurtriad, &
      v1, v2, v3, a_over_2_v3, a_over_2_v3_ref, a_over_2_v3_deriv, a_over_2_v3_second)
    use mMaterial, only: UPDATELAG
    implicit none

    integer, intent(in) :: flag
    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in) :: thick, elem(3, nn), shell_disp(6, nn)
    logical, intent(in) :: finite_rotation_director, use_director_tangent, need_second_tangent
    real(kind=kreal), intent(in), optional :: ndtriad(9, nn), ndreftriad(9, nn), ndcurtriad(9, nn)
    real(kind=kreal), intent(out) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal), intent(out) :: a_over_2_v3(3, nn), a_over_2_v3_ref(3, nn)
    real(kind=kreal), intent(out) :: a_over_2_v3_deriv(3, 3, nn), a_over_2_v3_second(3, 3, 3, nn)

    integer :: nb
    real(kind=kreal) :: director_ref(3), rotmat(3, 3)
    real(kind=kreal) :: nddirector(3, nn), ndrefdirector(3, nn), ndcurdirector(3, nn)
    logical :: has_nddirector, has_ndrefdirector, has_ndcurdirector

    has_nddirector = present(ndtriad)
    has_ndrefdirector = present(ndreftriad)
    has_ndcurdirector = present(ndcurtriad)
    if (has_nddirector)    nddirector(1:3, 1:nn)    = 0.5D0*thick*ndtriad(7:9, 1:nn)
    if (has_ndrefdirector) ndrefdirector(1:3, 1:nn) = 0.5D0*thick*ndreftriad(7:9, 1:nn)
    if (has_ndcurdirector) ndcurdirector(1:3, 1:nn) = 0.5D0*thick*ndcurtriad(7:9, 1:nn)

    a_over_2_v3_deriv = 0.0D0
    a_over_2_v3_second = 0.0D0
    call ShellMITC_SetupReferenceDirectors(etype, nn, thick, elem, v1, v2, v3, a_over_2_v3)

    do nb = 1, nn
      a_over_2_v3_ref(:, nb) = a_over_2_v3(:, nb)
      if (finite_rotation_director .and. flag == UPDATELAG .and. has_ndcurdirector) then
        a_over_2_v3_ref(:, nb) = ndcurdirector(:, nb)
      else if (finite_rotation_director .and. has_ndrefdirector) then
        a_over_2_v3_ref(:, nb) = ndrefdirector(:, nb)
      endif

      director_ref = a_over_2_v3_ref(:, nb)
      if (finite_rotation_director) then
        if (has_nddirector) then
          a_over_2_v3(:, nb) = nddirector(:, nb)
        else
          call ShellRotationVectorToMatrix(shell_disp(4:6, nb), rotmat)
          a_over_2_v3(:, nb) = matmul(rotmat, director_ref)
        endif
        if (use_director_tangent) then
          a_over_2_v3_deriv(:, :, nb) = ShellSkewMatrix(a_over_2_v3(:, nb))
          if( need_second_tangent ) then
            call ShellDirectorIncrementalSecondDeriv(a_over_2_v3(:, nb), a_over_2_v3_second(:, :, :, nb))
          endif
        endif
      endif

      if (.not. use_director_tangent) then
        a_over_2_v3_deriv(:, :, nb) = ShellSkewMatrix(a_over_2_v3(:, nb))
      endif
    end do
  end subroutine ShellMITC_SetupNodalDirectors

  !--------------------------------------------------------------------
  !> Evaluate one MITC tying point. First and second variations are returned separately
  !> so UPDATE can retain only the B matrix while STF also keeps geometric-stiffness data.
  subroutine ShellMITC_EvaluateTyingPointVariation(etype, nn, ndof, tying_set, tying_point, &
      zeta_tying, elem, shell_disp, director, director_tangent, use_green_lagrange, &
      use_director_tangent, point_B, point_basis_variation, director_second_tangent, point_director_second_variation)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof, tying_set, tying_point
    real(kind=kreal), intent(in) :: zeta_tying
    real(kind=kreal), intent(in) :: elem(3, nn), shell_disp(6, nn), director(3, nn)
    real(kind=kreal), intent(in) :: director_tangent(3, 3, nn)
    logical, intent(in) :: use_green_lagrange, use_director_tangent
    real(kind=kreal), intent(out) :: point_B(5, ndof*nn)
    real(kind=kreal), intent(out) :: point_basis_variation(3, ndof*nn, 3)
    real(kind=kreal), intent(in), optional :: director_second_tangent(3, 3, 3, nn)
    real(kind=kreal), intent(out), optional :: point_director_second_variation(5, 3, 3, nn)

    real(kind=kreal) :: naturalcoord(2), shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal) :: director_contribution(3, nn, 3)
    real(kind=kreal) :: kinematic_basis(3, 3), membrane_basis(3, 2)
    real(kind=kreal) :: translation_gradient(3, 2)

    call getTyingPoint(etype, tying_set, tying_point, naturalcoord)
    call getShapeFunc(etype, naturalcoord, shapefunc)
    call getShapeDeriv(etype, naturalcoord, shapederiv)

    call ShellMITC_DirectorContributions(nn, director, zeta_tying, shapefunc, shapederiv, director_contribution)
    call ShellMITC_CovariantBasis(nn, elem, director, zeta_tying, shapefunc, shapederiv, kinematic_basis)
    if( abs(ShellMITC_CovariantJacobian(kinematic_basis)) <= tiny(1.0D0) ) &
      stop "Invalid shell Jacobian"

    membrane_basis = kinematic_basis(:, SHELL_XI:SHELL_ETA)
    if( use_green_lagrange ) then
      translation_gradient = matmul(shell_disp(1:3, :), shapederiv)
      if( etype == fe_mitc9_shell ) then
        membrane_basis = membrane_basis+translation_gradient
      else
        kinematic_basis(:, SHELL_XI:SHELL_ETA) = kinematic_basis(:, SHELL_XI:SHELL_ETA)+translation_gradient
        membrane_basis = kinematic_basis(:, SHELL_XI:SHELL_ETA)
      endif
    endif

    call ShellMITC_BuildFirstStrainVariation(nn, ndof, zeta_tying, shapefunc, shapederiv, &
      kinematic_basis, membrane_basis, director_contribution, director_tangent, &
      use_director_tangent, point_B, point_basis_variation)

    if( use_director_tangent .and. present(director_second_tangent) .and. present(point_director_second_variation) ) then
      call ShellMITC_BuildDirectorSecondVariation(nn, zeta_tying, shapefunc, shapederiv, &
        kinematic_basis, director_second_tangent, point_director_second_variation)
    endif
  end subroutine ShellMITC_EvaluateTyingPointVariation

  !--------------------------------------------------------------------
  !> Evaluate only the assumed-strain B matrices needed by UPDATE at one zeta.
  subroutine ShellMITC_EvaluateTyingBAtZeta(etype, nn, ndof, zeta_tying, elem, shell_disp, &
      director, director_tangent, use_green_lagrange, use_director_tangent, tying_B)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof
    real(kind=kreal), intent(in) :: zeta_tying
    real(kind=kreal), intent(in) :: elem(3, nn), shell_disp(6, nn), director(3, nn)
    real(kind=kreal), intent(in) :: director_tangent(3, 3, nn)
    logical, intent(in) :: use_green_lagrange, use_director_tangent
    real(kind=kreal), intent(out) :: tying_B(5, ndof*nn, 6, 3)

    integer :: tying_set, tying_point
    real(kind=kreal) :: point_B(5, ndof*nn)
    real(kind=kreal) :: point_basis_variation(3, ndof*nn, 3)

    tying_B = 0.0D0
    do tying_set = 1, NumOfTyingSets(etype)
      do tying_point = 1, NumOfTyingPoints(etype, tying_set)
        call ShellMITC_EvaluateTyingPointVariation(etype, nn, ndof, tying_set, tying_point, &
          zeta_tying, elem, shell_disp, director, director_tangent, use_green_lagrange, &
          use_director_tangent, point_B, point_basis_variation)
        tying_B(:, :, tying_point, tying_set) = point_B
      end do
    end do
  end subroutine ShellMITC_EvaluateTyingBAtZeta

  !--------------------------------------------------------------------
  !> Evaluate all first/second variations needed by STF at one zeta.
  subroutine ShellMITC_EvaluateTyingStiffnessDataAtZeta(etype, nn, ndof, zeta_tying, &
      elem, shell_disp, director, director_tangent, director_second_tangent, &
      use_green_lagrange, use_director_tangent, tying_B, tying_basis_variation, tying_director_second_variation)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof
    real(kind=kreal), intent(in) :: zeta_tying
    real(kind=kreal), intent(in) :: elem(3, nn), shell_disp(6, nn), director(3, nn)
    real(kind=kreal), intent(in) :: director_tangent(3, 3, nn)
    real(kind=kreal), intent(in) :: director_second_tangent(3, 3, 3, nn)
    logical, intent(in) :: use_green_lagrange, use_director_tangent
    real(kind=kreal), intent(out) :: tying_B(5, ndof*nn, 6, 3)
    real(kind=kreal), intent(out) :: tying_basis_variation(3, ndof*nn, 3, 6, 3)
    real(kind=kreal), intent(out) :: tying_director_second_variation(5, 3, 3, nn, 6, 3)

    integer :: tying_set, tying_point
    real(kind=kreal) :: point_B(5, ndof*nn)
    real(kind=kreal) :: point_basis_variation(3, ndof*nn, 3)
    real(kind=kreal) :: point_director_second_variation(5, 3, 3, nn)

    tying_B = 0.0D0
    tying_basis_variation = 0.0D0
    tying_director_second_variation = 0.0D0

    do tying_set = 1, NumOfTyingSets(etype)
      do tying_point = 1, NumOfTyingPoints(etype, tying_set)
        if( use_director_tangent ) then
          call ShellMITC_EvaluateTyingPointVariation(etype, nn, ndof, tying_set, tying_point, &
            zeta_tying, elem, shell_disp, director, director_tangent, use_green_lagrange, &
            use_director_tangent, point_B, point_basis_variation, director_second_tangent, point_director_second_variation)
          tying_director_second_variation(:, :, :, :, tying_point, tying_set) = point_director_second_variation
        else
          call ShellMITC_EvaluateTyingPointVariation(etype, nn, ndof, tying_set, tying_point, &
            zeta_tying, elem, shell_disp, director, director_tangent, use_green_lagrange, &
            use_director_tangent, point_B, point_basis_variation)
        endif
        tying_B(:, :, tying_point, tying_set) = point_B
        tying_basis_variation(:, :, :, tying_point, tying_set) = point_basis_variation
      end do
    end do
  end subroutine ShellMITC_EvaluateTyingStiffnessDataAtZeta

  !--------------------------------------------------------------------
  !> Build the first strain variation and the covariant-basis variations.
  subroutine ShellMITC_BuildFirstStrainVariation(nn, ndof, zeta, shapefunc, shapederiv, &
      kinematic_basis, membrane_basis, director_contribution, director_tangent, &
      use_director_tangent, strain_variation, basis_variation)
    implicit none

    integer(kind=kint), intent(in) :: nn, ndof
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(in) :: kinematic_basis(3, 3), membrane_basis(3, 2)
    real(kind=kreal), intent(in) :: director_contribution(3, nn, 3)
    real(kind=kreal), intent(in) :: director_tangent(3, 3, nn)
    logical, intent(in) :: use_director_tangent
    real(kind=kreal), intent(out) :: strain_variation(5, ndof*nn)
    real(kind=kreal), intent(out) :: basis_variation(3, ndof*nn, 3)

    integer :: node, component, dof_index, rotation_start

    basis_variation = 0.0D0
    do node = 1, nn
      do component = 1, min(3, ndof)
        dof_index = ndof*(node-1)+component
        basis_variation(component, dof_index, SHELL_XI:SHELL_ETA) = shapederiv(node, :)
      end do
      if( ndof < 6 ) cycle

      rotation_start = ndof*(node-1)+4
      if( use_director_tangent ) then
        basis_variation(:, rotation_start:rotation_start+2, SHELL_XI) = &
          shapederiv(node, SHELL_XI)*zeta*director_tangent(:, :, node)
        basis_variation(:, rotation_start:rotation_start+2, SHELL_ETA) = &
          shapederiv(node, SHELL_ETA)*zeta*director_tangent(:, :, node)
        basis_variation(:, rotation_start:rotation_start+2, SHELL_ZETA) = shapefunc(node)*director_tangent(:, :, node)
      else
        basis_variation(:, rotation_start:rotation_start+2, SHELL_XI) = &
          ShellSkewMatrix(director_contribution(:, node, SHELL_XI))
        basis_variation(:, rotation_start:rotation_start+2, SHELL_ETA) = &
          ShellSkewMatrix(director_contribution(:, node, SHELL_ETA))
        basis_variation(:, rotation_start:rotation_start+2, SHELL_ZETA) = &
          ShellSkewMatrix(director_contribution(:, node, SHELL_ZETA))
      endif
    end do

    do dof_index = 1, ndof*nn
      strain_variation(:, dof_index) = ShellMITC_SymmetricComponents( &
        matmul(transpose(basis_variation(:, dof_index, :)), kinematic_basis))
    end do

    ! In the TL MITC9 formulation only the translational membrane columns use
    ! the translated metric basis; rotational columns retain the director basis.
    do node = 1, nn
      do component = 1, min(3, ndof)
        dof_index = ndof*(node-1)+component
        strain_variation(1, dof_index) = shapederiv(node, SHELL_XI)*membrane_basis(component, SHELL_XI)
        strain_variation(2, dof_index) = shapederiv(node, SHELL_ETA)*membrane_basis(component, SHELL_ETA)
        strain_variation(3, dof_index) = shapederiv(node, SHELL_XI)*membrane_basis(component, SHELL_ETA) &
          +shapederiv(node, SHELL_ETA)*membrane_basis(component, SHELL_XI)
      end do
    end do
  end subroutine ShellMITC_BuildFirstStrainVariation

  !--------------------------------------------------------------------
  !> Build the director contribution to the second strain variation used by STF.
  subroutine ShellMITC_BuildDirectorSecondVariation(nn, zeta, shapefunc, shapederiv, &
      kinematic_basis, director_second_tangent, director_second_variation)
    implicit none

    integer(kind=kint), intent(in) :: nn
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(in) :: kinematic_basis(3, 3)
    real(kind=kreal), intent(in) :: director_second_tangent(3, 3, 3, nn)
    real(kind=kreal), intent(out) :: director_second_variation(5, 3, 3, nn)

    integer :: node, m, n
    real(kind=kreal) :: second_basis(3, 3)

    director_second_variation = 0.0D0
    do node = 1, nn
      do n = 1, 3
        do m = 1, 3
          second_basis(:, SHELL_XI) = zeta*shapederiv(node, SHELL_XI) *director_second_tangent(:, m, n, node)
          second_basis(:, SHELL_ETA) = zeta*shapederiv(node, SHELL_ETA) *director_second_tangent(:, m, n, node)
          second_basis(:, SHELL_ZETA) = shapefunc(node) *director_second_tangent(:, m, n, node)
          director_second_variation(:, m, n, node) = ShellMITC_SymmetricComponents( &
            matmul(transpose(second_basis), kinematic_basis))
        end do
      end do
    end do
  end subroutine ShellMITC_BuildDirectorSecondVariation

  !--------------------------------------------------------------------
  !> Evaluate the sparse tying operator at one surface point.
  !> E_as(target) = sum_k coefficient(k)*E(source) at the selected tying point.
  subroutine ShellMITC_EvaluateTyingOperator(etype, xi, eta, nterms, target_component, &
      source_component, tying_point, tying_group, coefficient, replaced_component)
    implicit none

    integer(kind=kint), intent(in) :: etype
    real(kind=kreal), intent(in) :: xi, eta
    integer, intent(out) :: nterms
    integer, intent(out) :: target_component(MAX_TYING_TERMS), source_component(MAX_TYING_TERMS)
    integer, intent(out) :: tying_point(MAX_TYING_TERMS), tying_group(MAX_TYING_TERMS)
    real(kind=kreal), intent(out) :: coefficient(MAX_TYING_TERMS)
    logical, intent(out) :: replaced_component(5)

    integer :: ip, k
    real(kind=kreal) :: xxi, eeta, h

    nterms = 0

    if( etype == fe_mitc4_shell ) then
      nterms = 4
      target_component(1:4) = (/ 4, 4, 5, 5 /)
      source_component(1:4) = (/ 4, 4, 5, 5 /)
      tying_point (1:4) = (/ 4, 2, 1, 3 /)
      tying_group (1:4) = (/ 1, 1, 1, 1 /)
      coefficient  (1:4) = (/ 0.5D0*( 1.0D0-xi ), 0.5D0*( 1.0D0+xi ), 0.5D0*( 1.0D0-eta ), 0.5D0*( 1.0D0+eta ) /)

    else if( etype == fe_mitc9_shell ) then
      xxi  = xi /dsqrt( 1.0D0/3.0D0 )
      eeta = eta/dsqrt( 3.0D0/5.0D0 )
      do ip = 1, NumOfTyingPoints(etype, 1)
        h = ( 0.5D0*( 1.0D0+mitc9_xi_sign(ip, 1)*xxi ) )   &
          *( ( 0.5D0*mitc9_eta_sign(ip, 1)*eeta ) *( 1.0D0+mitc9_eta_sign(ip, 1)*eeta )       &
          +( 1.0D0-mitc9_eta_sign(ip, 1)*mitc9_eta_sign(ip, 1) ) *( 1.0D0-eeta*eeta ) )
        target_component(nterms+1:nterms+2) = (/ 1, 5 /)
        source_component(nterms+1:nterms+2) = (/ 1, 5 /)
        tying_point (nterms+1:nterms+2) = ip
        tying_group (nterms+1:nterms+2) = 1
        coefficient  (nterms+1:nterms+2) = h
        nterms = nterms+2
      end do

      xxi  = xi /dsqrt( 3.0D0/5.0D0 )
      eeta = eta/dsqrt( 1.0D0/3.0D0 )
      do ip = 1, NumOfTyingPoints(etype, 2)
        h = ( ( 0.5D0*mitc9_xi_sign(ip, 2)*xxi ) *( 1.0D0+mitc9_xi_sign(ip, 2)*xxi ) &
          +( 1.0D0-mitc9_xi_sign(ip, 2)*mitc9_xi_sign(ip, 2) ) *( 1.0D0-xxi*xxi ) ) &
          *( 0.5D0*( 1.0D0+mitc9_eta_sign(ip, 2)*eeta ) )
        target_component(nterms+1:nterms+2) = (/ 2, 4 /)
        source_component(nterms+1:nterms+2) = (/ 2, 4 /)
        tying_point (nterms+1:nterms+2) = ip
        tying_group (nterms+1:nterms+2) = 2
        coefficient  (nterms+1:nterms+2) = h
        nterms = nterms+2
      end do

      xxi  = xi /dsqrt( 1.0D0/3.0D0 )
      eeta = eta/dsqrt( 1.0D0/3.0D0 )
      do ip = 1, NumOfTyingPoints(etype, 3)
        h = ( 0.5D0*( 1.0D0+mitc9_xi_sign(ip, 1)*xxi ) ) *( 0.5D0*( 1.0D0+mitc9_eta_sign(ip, 1)*eeta ) )
        nterms = nterms+1
        target_component(nterms) = 3
        source_component(nterms) = 3
        tying_point (nterms) = ip
        tying_group (nterms) = 3
        coefficient  (nterms) = h
      end do

    else if( etype == fe_mitc3_shell ) then
      nterms = 8
      target_component(1:8) = (/ 4, 4, 4, 4, 5, 5, 5, 5 /)
      source_component(1:8) = (/ 4, 5, 4, 5, 4, 5, 4, 5 /)
      tying_point (1:8) = (/ 2, 1, 3, 3, 2, 1, 3, 3 /)
      tying_group (1:8) = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)
      coefficient  (1:8) = (/ 1.0D0-xi, xi, xi, -xi, eta, 1.0D0-eta, -eta, eta /)
    endif

    replaced_component(:) = .false.
    do k = 1, nterms
      replaced_component(target_component(k)) = .true.
    end do

  end subroutine ShellMITC_EvaluateTyingOperator

  !--------------------------------------------------------------------
  !> Replace the ordinary shell strain variations by the MITC assumed-strain field.
  subroutine ShellMITC_ApplyAssumedStrain(etype, nn, ndof, xi_lx, eta_lx, tying_B, B, &
      tying_director_second_variation, B2rot)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof
    real(kind=kreal), intent(in) :: xi_lx, eta_lx
    real(kind=kreal), intent(in) :: tying_B(5, ndof*nn, 6, 3)
    real(kind=kreal), intent(inout) :: B(5, ndof*nn)
    real(kind=kreal), intent(in), optional :: tying_director_second_variation(5, 3, 3, nn, 6, 3)
    real(kind=kreal), intent(inout), optional :: B2rot(:, :, :, :)

    integer :: k, c, jsize, m, n, nb, nterms
    integer :: target_component(MAX_TYING_TERMS), source_component(MAX_TYING_TERMS)
    integer :: tying_point(MAX_TYING_TERMS), tying_group(MAX_TYING_TERMS)
    real(kind=kreal) :: coefficient(MAX_TYING_TERMS)
    logical :: replaced_component(5)

    call ShellMITC_EvaluateTyingOperator(etype, xi_lx, eta_lx, nterms, target_component, &
      source_component, tying_point, tying_group, coefficient, replaced_component)
    if( nterms == 0 ) return

    do c = 1, 5
      if( replaced_component(c) ) B(c, 1:ndof*nn) = 0.0D0
    end do
    do k = 1, nterms
      do jsize = 1, ndof*nn
        B(target_component(k), jsize) = B(target_component(k), jsize) &
          +coefficient(k)*tying_B(source_component(k), jsize, tying_point(k), tying_group(k))
      end do
    end do

    if( present(tying_director_second_variation) .and. present(B2rot) ) then
      do c = 1, 5
        if( replaced_component(c) ) B2rot(c, :, :, :) = 0.0D0
      end do
      do k = 1, nterms
        do nb = 1, nn
          do n = 1, 3
            do m = 1, 3
              B2rot(target_component(k), m, n, nb) = B2rot(target_component(k), m, n, nb) &
                +coefficient(k)*tying_director_second_variation(source_component(k), m, n, nb, &
                tying_point(k), tying_group(k))
            end do
          end do
        end do
      end do
    endif

  end subroutine ShellMITC_ApplyAssumedStrain

  !--------------------------------------------------------------------
  !> Second variation of strain component c for the basis variations at dof i and dof j.
  pure real(kind=kreal) function ShellMITC_StrainSecondVariationPair(component, basis_variation_i, basis_variation_j)
    implicit none

    integer, intent(in) :: component
    real(kind=kreal), intent(in) :: basis_variation_i(3, 3), basis_variation_j(3, 3)

    select case( component )
      case( 1 )
        ShellMITC_StrainSecondVariationPair = dot_product( basis_variation_i(:, SHELL_XI), basis_variation_j(:, SHELL_XI))
      case( 2 )
        ShellMITC_StrainSecondVariationPair = dot_product( basis_variation_i(:, SHELL_ETA), basis_variation_j(:, SHELL_ETA))
      case( 3 )
        ShellMITC_StrainSecondVariationPair = dot_product( &
          basis_variation_i(:, SHELL_XI), basis_variation_j(:, SHELL_ETA)) &
          +dot_product(basis_variation_i(:, SHELL_ETA), basis_variation_j(:, SHELL_XI))
      case( 4 )
        ShellMITC_StrainSecondVariationPair = dot_product( &
          basis_variation_i(:, SHELL_ETA), basis_variation_j(:, SHELL_ZETA)) &
          +dot_product(basis_variation_i(:, SHELL_ZETA), basis_variation_j(:, SHELL_ETA))
      case( 5 )
        ShellMITC_StrainSecondVariationPair = dot_product( &
          basis_variation_i(:, SHELL_ZETA), basis_variation_j(:, SHELL_XI)) &
          +dot_product(basis_variation_i(:, SHELL_XI), basis_variation_j(:, SHELL_ZETA))
      case default
        ShellMITC_StrainSecondVariationPair = 0.0D0
    end select

  end function ShellMITC_StrainSecondVariationPair

  !--------------------------------------------------------------------
  !> Assemble initial-stress, volume-change and director geometric stiffness.
  subroutine ShellMITC_AddGeometricStiffness(etype, nn, ndof, kinematics, &
      use_director_tangent, use_green_lagrange, add_geometric_stiffness, &
      xi, eta, integration_weight, layer_weight, stress, B, basis_variation, &
      tying_basis_variation, reciprocal_basis, material_local_basis, material_reciprocal_basis, B2rot, director, stiff)
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof, kinematics
    logical, intent(in) :: use_director_tangent, use_green_lagrange
    logical, intent(in) :: add_geometric_stiffness
    real(kind=kreal), intent(in) :: xi, eta, integration_weight, layer_weight
    real(kind=kreal), intent(in) :: stress(5), B(5, ndof*nn)
    real(kind=kreal), intent(in) :: basis_variation(3, ndof*nn, 3)
    real(kind=kreal), intent(in) :: tying_basis_variation(3, ndof*nn, 3, 6, 3)
    real(kind=kreal), intent(in) :: reciprocal_basis(3, 3)
    real(kind=kreal), intent(in) :: material_local_basis(3, 3)
    real(kind=kreal), intent(in) :: material_reciprocal_basis(3, 3)
    real(kind=kreal), intent(in) :: B2rot(5, 3, 3, nn)
    real(kind=kreal), intent(in) :: director(3, nn)
    real(kind=kreal), intent(inout) :: stiff(ndof*nn, ndof*nn)

    integer :: isize, jsize, k, c, ip, it, nterms
    integer :: target_component(MAX_TYING_TERMS), source_component(MAX_TYING_TERMS)
    integer :: tying_point(MAX_TYING_TERMS), tying_group(MAX_TYING_TERMS)
    real(kind=kreal) :: coefficient(MAX_TYING_TERMS)
    logical :: replaced_component(5)
    real(kind=kreal) :: qf_stress_integrand(ndof*nn), vol_deriv(ndof*nn)
    real(kind=kreal) :: geo_term, local_stress(3, 3), S_global(3, 3)
    real(kind=kreal) :: basis_gradient(3, 3, ndof*nn)
    real(kind=kreal) :: stress_basis_gradient(3, 3, ndof*nn)

    if( kinematics == UPDATELAG ) then
      do jsize = 1, ndof*nn
        vol_deriv(jsize) = sum(reciprocal_basis*basis_variation(:, jsize, :))
        qf_stress_integrand(jsize) = dot_product(stress, B(:, jsize))
      end do
    endif

    if( add_geometric_stiffness ) then
      if( use_green_lagrange .or. kinematics == UPDATELAG ) then
        call ShellMITC_EvaluateTyingOperator(etype, xi, eta, nterms, target_component, &
          source_component, tying_point, tying_group, coefficient, replaced_component)
        do jsize = 1, ndof*nn
          do isize = 1, ndof*nn
            geo_term = 0.0D0
            do c = 1, 5
              if( replaced_component(c) ) cycle
              geo_term = geo_term+stress(c)*ShellMITC_StrainSecondVariationPair(c, &
                basis_variation(:, isize, :), basis_variation(:, jsize, :))
            end do
            do k = 1, nterms
              ip = tying_point(k)
              it = tying_group(k)
              geo_term = geo_term+stress(target_component(k))*coefficient(k) &
                *ShellMITC_StrainSecondVariationPair(source_component(k), tying_basis_variation(:, isize, :, ip, it), &
                tying_basis_variation(:, jsize, :, ip, it))
            end do
            stiff(isize, jsize) = stiff(isize, jsize) +integration_weight*layer_weight*geo_term
          end do
        end do
      else
        local_stress = reshape((/ stress(1), stress(3), stress(5), &
          stress(3), stress(2), stress(4), stress(5), stress(4), 0.0D0 /), (/ 3, 3 /))
        S_global = matmul(material_local_basis, matmul(local_stress, transpose(material_local_basis)))

        do jsize = 1, ndof*nn
          basis_gradient(:, :, jsize) = matmul(basis_variation(:, jsize, :), transpose(material_reciprocal_basis))
          stress_basis_gradient(:, :, jsize) = matmul(basis_gradient(:, :, jsize), S_global)
        end do
        do jsize = 1, ndof*nn
          do isize = 1, ndof*nn
            stiff(isize, jsize) = stiff(isize, jsize) +integration_weight*layer_weight*sum( &
              basis_gradient(:, :, isize)*stress_basis_gradient(:, :, jsize))
          end do
        end do
      endif
      if( kinematics == UPDATELAG ) then
        do jsize = 1, ndof*nn
          do isize = 1, ndof*nn
            stiff(isize, jsize) = stiff(isize, jsize) &
              +integration_weight*layer_weight*qf_stress_integrand(isize)*vol_deriv(jsize)
          end do
        end do
      endif
      if( use_director_tangent .and. (use_green_lagrange .or. kinematics == UPDATELAG) ) then
        call ShellMITC_AddDirectorStressStiffness(nn, ndof, integration_weight, layer_weight, &
          stress, B2rot, director, stiff)
      endif
    endif

  end subroutine ShellMITC_AddGeometricStiffness

  !--------------------------------------------------------------------
  !> Build the director second variation used by the drilling term.
  subroutine ShellMITC_BuildDrillingSecondVariation(nn, zeta, shapefunc, shapederiv, &
      director_second_tangent, reciprocal_basis, local_basis, use_director_tangent, drilling_second_variation)
    implicit none

    integer(kind=kint), intent(in) :: nn
    real(kind=kreal), intent(in) :: zeta, shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal), intent(in) :: director_second_tangent(3, 3, 3, nn)
    real(kind=kreal), intent(in) :: reciprocal_basis(3, 3), local_basis(3, 3)
    logical, intent(in) :: use_director_tangent
    real(kind=kreal), intent(out) :: drilling_second_variation(3, 3, nn)

    integer :: node, m, n
    real(kind=kreal) :: second_basis(3, 3), skew_gradient(3, 3)

    drilling_second_variation = 0.0D0
    if( .not. use_director_tangent ) return
    do node = 1, nn
      do n = 1, 3
        do m = 1, 3
          second_basis(:, SHELL_XI) = shapederiv(node, SHELL_XI)*zeta*director_second_tangent(:, m, n, node)
          second_basis(:, SHELL_ETA) = shapederiv(node, SHELL_ETA)*zeta*director_second_tangent(:, m, n, node)
          second_basis(:, SHELL_ZETA) = shapefunc(node)*director_second_tangent(:, m, n, node)
          skew_gradient = matmul(reciprocal_basis, transpose(second_basis))
          skew_gradient = skew_gradient-transpose(skew_gradient)
          drilling_second_variation(m, n, node) = dot_product( local_basis(:, SHELL_XI), &
            matmul(skew_gradient, local_basis(:, SHELL_ETA)))
        end do
      end do
    end do

  end subroutine ShellMITC_BuildDrillingSecondVariation

  !--------------------------------------------------------------------
  !> Build the drilling compatibility vector and its generalized strain.
  subroutine ShellMITC_BuildDrillingVector(nn, ndof, finite_rotation_director, shapefunc, &
      point_triad, reciprocal_basis, basis_variation, displacement, director, nddrill, Cv, Cv_disp)
    implicit none

    integer(kind=kint), intent(in) :: nn, ndof
    logical, intent(in) :: finite_rotation_director
    real(kind=kreal), intent(in) :: shapefunc(nn), point_triad(3, 3)
    real(kind=kreal), intent(in) :: reciprocal_basis(3, 3)
    real(kind=kreal), intent(in) :: basis_variation(3, ndof*nn, 3)
    real(kind=kreal), intent(in) :: displacement(ndof*nn), director(3, nn)
    real(kind=kreal), intent(in), optional :: nddrill(nn)
    real(kind=kreal), intent(out) :: Cv(ndof*nn), Cv_disp

    integer :: j, nb, jrot
    real(kind=kreal) :: Cv_w(ndof*nn), Cv_theta(ndof*nn)
    real(kind=kreal) :: cmat(3, 3), drill_axis(3), axis_norm

    do j = 1, ndof*nn
      cmat = matmul(reciprocal_basis, transpose(basis_variation(:, j, :)))
      cmat = cmat-transpose(cmat)
      Cv_w(j) = dot_product(point_triad(:, SHELL_XI), matmul(cmat, point_triad(:, SHELL_ETA)))
    end do

    Cv_theta = 0.0D0
    if( ndof >= 6 ) then
      do nb = 1, nn
        jrot = ndof*(nb-1)+4
        drill_axis = point_triad(:, SHELL_ZETA)
        if( finite_rotation_director .and. present(nddrill) ) then
          drill_axis = director(:, nb)
          axis_norm = sqrt(dot_product(drill_axis, drill_axis))
          if( axis_norm > 0.0D0 ) then
            drill_axis = drill_axis/axis_norm
          else
            drill_axis = point_triad(:, SHELL_ZETA)
          endif
        endif
        Cv_theta(jrot:jrot+2) = shapefunc(nb)*drill_axis
      end do
    endif

    Cv = Cv_theta-0.5D0*Cv_w
    Cv_disp = dot_product(Cv, displacement)
    if( finite_rotation_director .and. present(nddrill) ) then
      Cv_disp = -0.5D0*dot_product(Cv_w, displacement)+dot_product(shapefunc, nddrill)
    endif
  end subroutine ShellMITC_BuildDrillingVector

  !--------------------------------------------------------------------
  !> Add the drilling stabilization tangent.
  subroutine ShellMITC_AddDrillingStiffness(nn, ndof, finite_rotation_director, &
      use_director_tangent, shapefunc, Cv, Cv_disp, Cv_w_second, displacement, &
      director, director_deriv, alpha, integration_weight, layer_weight, nddrill, stiff)
    implicit none

    integer(kind=kint), intent(in) :: nn, ndof
    logical, intent(in) :: finite_rotation_director, use_director_tangent
    real(kind=kreal), intent(in) :: shapefunc(nn), Cv(ndof*nn), Cv_disp
    real(kind=kreal), intent(in) :: Cv_w_second(3, 3, nn), displacement(ndof*nn)
    real(kind=kreal), intent(in) :: director(3, nn), director_deriv(3, 3, nn)
    real(kind=kreal), intent(in) :: alpha, integration_weight, layer_weight
    real(kind=kreal), intent(in), optional :: nddrill(nn)
    real(kind=kreal), intent(inout) :: stiff(ndof*nn, ndof*nn)

    integer :: isize, jsize, nb, m, n
    real(kind=kreal) :: scale, Cv_deriv, Cv_deriv_disp
    real(kind=kreal) :: drill_axis(3), drill_coeff, axis_norm

    scale = integration_weight*layer_weight*alpha
    do jsize = 1, ndof*nn
      do isize = 1, ndof*nn
        stiff(isize, jsize) = stiff(isize, jsize)+scale*Cv(isize)*Cv(jsize)
      end do
    end do
    if( .not. use_director_tangent ) return

    do nb = 1, nn
      do n = 1, 3
        jsize = ndof*(nb-1)+3+n
        Cv_deriv_disp = 0.0D0
        do m = 1, 3
          isize = ndof*(nb-1)+3+m
          Cv_deriv = -0.5D0*Cv_w_second(m, n, nb)
          if( finite_rotation_director .and. present(nddrill) ) then
            drill_axis = director(:, nb)
            axis_norm = sqrt(dot_product(drill_axis, drill_axis))
            if( axis_norm > 0.0D0 ) then
              drill_axis = drill_axis/axis_norm
              drill_coeff = dot_product(drill_axis, director_deriv(:, n, nb))
              Cv_deriv = Cv_deriv+shapefunc(nb) *(director_deriv(m, n, nb)-drill_axis(m)*drill_coeff)/axis_norm
            endif
          endif
          Cv_deriv_disp = Cv_deriv_disp+Cv_deriv*displacement(isize)
          stiff(isize, jsize) = stiff(isize, jsize)+scale*Cv_deriv*Cv_disp
        end do
        stiff(:, jsize) = stiff(:, jsize)+scale*Cv*Cv_deriv_disp
      end do
    end do
  end subroutine ShellMITC_AddDrillingStiffness

  !--------------------------------------------------------------------
  !> Add the stress stiffness associated with the nonlinear director map.
  subroutine ShellMITC_AddDirectorStressStiffness(nn, ndof, w_w_w_det, layer_weight, Sv_force, B2rot, a_over_2_v3, stiff)
    implicit none

    integer(kind=kint), intent(in) :: nn, ndof
    real(kind=kreal), intent(in) :: w_w_w_det, layer_weight, Sv_force(5)
    real(kind=kreal), intent(in) :: B2rot(5, 3, 3, nn), a_over_2_v3(3, nn)
    real(kind=kreal), intent(inout) :: stiff(ndof*nn, ndof*nn)

    integer :: i, j, m, n, nb, isize, jsize
    real(kind=kreal) :: drill_axis(3), rot_projector(3, 3), hrot(3, 3)
    real(kind=kreal) :: hess_coeff, drill_coeff, axis_norm

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
              hess_coeff = hess_coeff +rot_projector(i, m)*hrot(i, j)*rot_projector(j, n)
            end do
          end do
          drill_coeff = 0.0D0
          do i = 1, 3
            drill_coeff = drill_coeff+drill_axis(i)*hrot(i, n)
          end do
          do j = 1, 3
            do i = 1, 3
              drill_coeff = drill_coeff +rot_projector(i, n)*hrot(i, j)*drill_axis(j)
            end do
          end do
          stiff(isize, jsize) = stiff(isize, jsize) +w_w_w_det*layer_weight *(hess_coeff+drill_axis(m)*drill_coeff)
        end do
      end do
    end do

  end subroutine ShellMITC_AddDirectorStressStiffness

  !--------------------------------------------------------------------
  !> Return the map from shell-solid mixed ordering to natural shell ordering.
  subroutine ShellMixedDofMap(mixflag, nn, ndof, sstable, is_mixed)
    implicit none

    integer(kind=kint), intent(in) :: mixflag, nn, ndof
    integer, intent(out) :: sstable(ndof*nn)
    logical, intent(out) :: is_mixed

    integer :: i

    is_mixed = .true.
    if( mixflag == 1 .and. ndof*nn == 24 ) then
      sstable = (/ 1, 2, 3, 7, 8, 9, 13, 14, 15, 19, 20, 21, 4, 5, 6, 10, 11, 12, 16, 17, 18, 22, 23, 24 /)
    else if( mixflag == 2 .and. ndof*nn == 18 ) then
      sstable = (/ 1, 2, 3, 7, 8, 9, 13, 14, 15, 4, 5, 6, 10, 11, 12, 16, 17, 18 /)
    else
      is_mixed = .false.
      sstable = [(i, i=1, ndof*nn)]
    endif
  end subroutine ShellMixedDofMap

  !--------------------------------------------------------------------
  !> Convert the natural shell ordering to the shell-solid mixed ordering.
  subroutine ShellApplyMixedDofOrdering(mixflag, nn, ndof, stiff, qf)
    implicit none

    integer(kind=kint), intent(in) :: mixflag, nn, ndof
    real(kind=kreal), intent(inout), optional :: stiff(ndof*nn, ndof*nn), qf(ndof*nn)

    integer :: i, j, sstable(ndof*nn)
    real(kind=kreal) :: stiff_work(ndof*nn, ndof*nn), qf_work(ndof*nn)
    logical :: is_mixed

    call ShellMixedDofMap(mixflag, nn, ndof, sstable, is_mixed)
    if( .not. is_mixed ) return

    if( present(stiff) ) then
      stiff_work = stiff
      do j = 1, ndof*nn
        do i = 1, ndof*nn
          stiff(i, j) = stiff_work(sstable(i), sstable(j))
        end do
      end do
    endif
    if( present(qf) ) then
      qf_work = qf
      do i = 1, ndof*nn
        qf(i) = qf_work(sstable(i))
      end do
    endif
  end subroutine ShellApplyMixedDofOrdering

  !--------------------------------------------------------------------
  !> Evaluate strains at MITC tying points before interpolation.
  subroutine ShellMITC_EvaluateTyingStrains(etype, nn, zeta, elem, edisp, director, &
      director_increment, use_green_lagrange, tying_strain)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in) :: zeta
    real(kind=kreal), intent(in) :: elem(3, nn), edisp(6, nn)
    real(kind=kreal), intent(in) :: director(3, nn), director_increment(3, nn)
    logical, intent(in) :: use_green_lagrange
    real(kind=kreal), intent(out) :: tying_strain(5, 6, 3)

    integer :: tying_set, tying_point
    real(kind=kreal) :: tying_zeta, naturalcoord(2)
    real(kind=kreal) :: shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal) :: covariant_basis(3, 3), reference_basis(3, 3), current_basis(3, 3)
    real(kind=kreal) :: reference_jacobian, current_jacobian, point_strain(5)

    tying_strain = 0.0D0
    tying_zeta = ShellMITC_TyingZeta(etype, zeta)

    do tying_set = 1, NumOfTyingSets(etype)
      do tying_point = 1, NumOfTyingPoints(etype, tying_set)
        call getTyingPoint(etype, tying_set, tying_point, naturalcoord)
        call getShapeFunc(etype, naturalcoord, shapefunc)
        call getShapeDeriv(etype, naturalcoord, shapederiv)
        call ShellMITC_CovariantBasis(nn, elem, director, tying_zeta, shapefunc, shapederiv, covariant_basis)
        if( abs(ShellMITC_CovariantJacobian(covariant_basis)) <= tiny(1.0D0) ) &
          stop "Invalid shell Jacobian"
        call ShellMITC_EvaluatePointStrain(nn, tying_zeta, shapefunc, shapederiv, &
          edisp(1:3, :), director_increment, covariant_basis, use_green_lagrange, &
          point_strain, reference_basis, current_basis, reference_jacobian, current_jacobian)
        tying_strain(:, tying_point, tying_set) = point_strain
      end do
    end do

  end subroutine ShellMITC_EvaluateTyingStrains

  !--------------------------------------------------------------------
  !> Prepare the nodal and director kinematics used by stress evaluation.
  !> Combine director values already computed by the caller's ShellMITC_SetupNodalDirectors
  !> call (a_over_2_v3, a_over_2_v3_ref, a_over_2_v3_deriv) into the stress-evaluation
  !> configuration (elem, director, director_increment). Does not recompute directors itself.
  subroutine ShellMITC_SetupStressKinematics(nn, ecoord, kinematics, finite_rotation, edisp, &
      nddirector, ndrefdirector, nddirector_deriv, ndbase_disp, elem, director, director_increment)
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: nn, kinematics
    logical, intent(in) :: finite_rotation
    real(kind=kreal), intent(in) :: ecoord(3, nn), edisp(6, nn)
    real(kind=kreal), intent(in) :: nddirector(3, nn), ndrefdirector(3, nn)
    real(kind=kreal), intent(in) :: nddirector_deriv(3, 3, nn)
    real(kind=kreal), intent(in), optional :: ndbase_disp(6, nn)
    real(kind=kreal), intent(out) :: elem(3, nn)
    real(kind=kreal), intent(out) :: director(3, nn), director_increment(3, nn)

    integer :: na

    elem = ecoord
    if( kinematics == UPDATELAG .and. present(ndbase_disp) ) then
      elem = elem+ndbase_disp(1:3, :)
    endif
    if( kinematics == UPDATELAG ) elem = elem+0.5D0*edisp(1:3, :)

    do na = 1, nn
      director(:, na) = nddirector(:, na)
      director_increment(:, na) = director(:, na)-ndrefdirector(:, na)
      if( .not. finite_rotation ) then
        director_increment(:, na) = matmul(nddirector_deriv(:, :, na), edisp(4:6, na))
      else if( kinematics == UPDATELAG ) then
        director(:, na) = ndrefdirector(:, na)+0.5D0*director_increment(:, na)
      endif
    end do
  end subroutine ShellMITC_SetupStressKinematics

  !--------------------------------------------------------------------
  !> Evaluate the MITC assumed strain and stress-evaluation bases at one point.
  subroutine ShellMITC_EvaluateAssumedStrain(etype, nn, naturalcoord, zeta, elem, edisp, &
      director, director_increment, use_green_lagrange, tying_strain, dstrain, &
      strain_tensor, covariant_basis, reciprocal_basis, material_local_basis, &
      material_reciprocal_basis, reference_basis, current_basis, reference_jacobian, current_jacobian)
    implicit none

    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in) :: naturalcoord(2), zeta
    real(kind=kreal), intent(in) :: elem(3, nn), edisp(6, nn)
    real(kind=kreal), intent(in) :: director(3, nn), director_increment(3, nn)
    logical, intent(in) :: use_green_lagrange
    real(kind=kreal), intent(in) :: tying_strain(5, 6, 3)
    real(kind=kreal), intent(out) :: dstrain(6), strain_tensor(3, 3)
    real(kind=kreal), intent(out) :: covariant_basis(3, 3), reciprocal_basis(3, 3)
    real(kind=kreal), intent(out) :: material_local_basis(3, 3)
    real(kind=kreal), intent(out) :: material_reciprocal_basis(3, 3)
    real(kind=kreal), intent(out) :: reference_basis(3, 3), current_basis(3, 3)
    real(kind=kreal), intent(out) :: reference_jacobian, current_jacobian

    integer :: k, c, nterms
    integer :: target_component(MAX_TYING_TERMS), source_component(MAX_TYING_TERMS)
    integer :: tying_point(MAX_TYING_TERMS), tying_group(MAX_TYING_TERMS)
    real(kind=kreal) :: coefficient(MAX_TYING_TERMS)
    logical :: replaced_component(5)
    real(kind=kreal) :: shapefunc(nn), shapederiv(nn, 2), strain(5)
    real(kind=kreal) :: local_basis(3, 3), jacobian, material_jacobian

    call getShapeFunc(etype, naturalcoord, shapefunc)
    call getShapeDeriv(etype, naturalcoord, shapederiv)
    call ShellMITC_CovariantBasis(nn, elem, director, zeta, shapefunc, shapederiv, covariant_basis)
    if( abs(ShellMITC_CovariantJacobian(covariant_basis)) <= tiny(1.0D0) ) &
      stop "Invalid shell Jacobian"
    call ShellMITC_BasisFromCovariant(covariant_basis, local_basis, reciprocal_basis, jacobian)
    call ShellMITC_EvaluatePointStrain(nn, zeta, shapefunc, shapederiv, edisp(1:3, :), &
      director_increment, covariant_basis, use_green_lagrange, strain, reference_basis, &
      current_basis, reference_jacobian, current_jacobian)
    call ShellMITC_EvaluateTyingOperator(etype, naturalcoord(1), naturalcoord(2), nterms, &
      target_component, source_component, tying_point, tying_group, coefficient, replaced_component)
    do c = 1, 5
      if( replaced_component(c) ) strain(c) = 0.0D0
    end do
    do k = 1, nterms
      strain(target_component(k)) = strain(target_component(k)) &
        +coefficient(k)*tying_strain(source_component(k), tying_point(k), tying_group(k))
    end do

    strain_tensor = 0.0D0
    strain_tensor(1, 1) = strain(1)
    strain_tensor(2, 2) = strain(2)
    strain_tensor(1, 2) = 0.5D0*strain(3)
    strain_tensor(2, 1) = strain_tensor(1, 2)
    strain_tensor(2, 3) = 0.5D0*strain(4)
    strain_tensor(3, 2) = strain_tensor(2, 3)
    strain_tensor(3, 1) = 0.5D0*strain(5)
    strain_tensor(1, 3) = strain_tensor(3, 1)

    material_local_basis = local_basis
    material_reciprocal_basis = reciprocal_basis
    if( use_green_lagrange ) then
      call ShellMITC_BasisFromCovariant(reference_basis, material_local_basis, material_reciprocal_basis, material_jacobian)
    endif

    dstrain = (/ strain(1), strain(2), 0.0D0, strain(3), strain(4), strain(5) /)
  end subroutine ShellMITC_EvaluateAssumedStrain

  !--------------------------------------------------------------------
  !> Calculate shell stress and optionally update the integration-point state.
  !> The constitutive matrix is local to this routine, as in Update_Stress3D.
  subroutine ShellMITC_UpdateStress(flag, update_state, gauss, n_layer, dstrain, &
      material_local_basis, material_reciprocal_basis, stress, alpha)
    use mMechGauss
    use m_MatMatrix
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: flag
    integer, intent(in) :: n_layer
    logical, intent(in) :: update_state
    type(tGaussStatus), intent(inout) :: gauss
    real(kind=kreal), intent(in) :: dstrain(6)
    real(kind=kreal), intent(in) :: material_local_basis(3, 3)
    real(kind=kreal), intent(in) :: material_reciprocal_basis(3, 3)
    real(kind=kreal), intent(out) :: stress(6), alpha

    real(kind=kreal) :: D(5, 5), strain(5), stress_shell(5)
    real(kind=kreal) :: dstress(6), dstress_trace(6), trace_coeff

    strain = (/ dstrain(1), dstrain(2), dstrain(4), dstrain(5), dstrain(6) /)
    call MatlMatrix_Shell(gauss, Shell, D, material_local_basis(:, SHELL_XI), &
      material_local_basis(:, SHELL_ETA), material_local_basis(:, SHELL_ZETA), material_reciprocal_basis(:, SHELL_XI), &
      material_reciprocal_basis(:, SHELL_ETA), material_reciprocal_basis(:, SHELL_ZETA), alpha, n_layer)

    stress_shell = matmul(D, strain)
    dstress = (/ stress_shell(1), stress_shell(2), 0.0D0, stress_shell(3), stress_shell(4), stress_shell(5) /)
    stress = dstress
    if( .not. update_state ) return

    if( flag == UPDATELAG ) then
      trace_coeff = ShellPlaneStressTraceCoeff(gauss, n_layer)
      call ShellObjectiveTraceStressIncrement(gauss%stress_bak(1:6), dstrain, dstress_trace, trace_coeff)
      gauss%strain(1:6) = gauss%strain_bak(1:6)+dstrain
      gauss%stress(1:6) = gauss%stress_bak(1:6)+dstress_trace+dstress
      gauss%strain_energy = gauss%strain_energy_bak+dot_product(gauss%stress(1:6), dstrain)
      gauss%strain_energy = gauss%strain_energy-0.5D0*dot_product(dstress, dstrain)
    else
      gauss%strain(1:6) = dstrain
      gauss%stress(1:6) = dstress
      gauss%strain_energy = 0.5D0*dot_product(gauss%stress(1:6), gauss%strain(1:6))
    endif
    stress = gauss%stress(1:6)
  end subroutine ShellMITC_UpdateStress

  !--------------------------------------------------------------------
  !> Transform shell-local stress and strain to the requested output measure.
  subroutine ShellMITC_TransformOutput(use_gl_strain, S, E, covariant_basis, &
      reciprocal_basis, reference_basis, current_basis, det_ref, det_cur, strain_out, stress_out)
    use m_fstr, only: OPSSTYPE, kOPSS_SOLUTION
    use m_utilities, only: get_principal
    implicit none

    logical, intent(in) :: use_gl_strain
    real(kind=kreal), intent(in) :: S(3, 3), E(3, 3)
    real(kind=kreal), intent(in) :: covariant_basis(3, 3), reciprocal_basis(3, 3)
    real(kind=kreal), intent(in) :: reference_basis(3, 3), current_basis(3, 3)
    real(kind=kreal), intent(in) :: det_ref, det_cur
    real(kind=kreal), intent(out) :: strain_out(6), stress_out(6)

    integer :: i
    real(kind=kreal) :: output_basis(3, 3), output_reciprocal_basis(3, 3)
    real(kind=kreal) :: reference_reciprocal_basis(3, 3), reference_local_basis(3, 3)
    real(kind=kreal) :: stress_tensor(3, 3), strain_tensor(3, 3)
    real(kind=kreal) :: stretch_b(3, 3), tensor(6), eigval(3), princ(3, 3)
    real(kind=kreal) :: logstrain(3, 3), cg_metric(3, 3)
    real(kind=kreal) :: det, jac, eig_norm

    output_basis = covariant_basis
    output_reciprocal_basis = reciprocal_basis

    if( use_gl_strain ) then
      output_basis = reference_basis
      call ShellMITC_BasisFromCovariant(reference_basis, reference_local_basis, reference_reciprocal_basis, det)
      output_reciprocal_basis = reference_reciprocal_basis
    endif

    stress_tensor = matmul(output_basis, matmul(S, transpose(output_basis)))
    strain_tensor = matmul(output_reciprocal_basis, matmul(E, transpose(output_reciprocal_basis)))
    call ShellTensorToStressVector(stress_tensor, stress_out)
    strain_out(1:3) = (/ strain_tensor(1, 1), strain_tensor(2, 2), strain_tensor(3, 3) /)
    ! The finite-rotation Green-Lagrange output reports the tensor shear (no
    ! engineering factor 2), matching the original ElementStress_Shell_MITC;
    ! only the small-strain/UL output uses the engineering shear.
    if( use_gl_strain ) then
      strain_out(4:6) = (/ strain_tensor(1, 2), strain_tensor(2, 3), strain_tensor(3, 1) /)
    else
      strain_out(4:6) = 2.0D0*(/ strain_tensor(1, 2), strain_tensor(2, 3), strain_tensor(3, 1) /)
    endif

    if( .not. (use_gl_strain .and. OPSSTYPE == kOPSS_SOLUTION) ) return

    jac = det_cur/det_ref
    if( abs(jac) <= tiny(1.0D0) ) stop "Fail to convert shell stress: detF=0"

    stress_tensor = matmul(current_basis, matmul(S, transpose(current_basis)))/jac
    cg_metric = matmul(transpose(reference_reciprocal_basis), reference_reciprocal_basis)
    stretch_b = matmul(current_basis, matmul(cg_metric, transpose(current_basis)))

    call ShellTensorToStressVector(stretch_b, tensor)
    call get_principal(tensor, eigval, princ)
    do i = 1, 3
      if( eigval(i) <= 0.0D0 ) stop "Fail to calc shell log strain: stretch<0"
      eigval(i) = 0.5D0*dlog(eigval(i))
      eig_norm = dsqrt(dot_product(princ(:, i), princ(:, i)))
      if( eig_norm <= 0.0D0 ) stop "Fail to calc shell log strain: direction vector=0"
      princ(:, i) = princ(:, i)/eig_norm
    end do

    logstrain = 0.0D0
    do i = 1, 3
      logstrain = logstrain+eigval(i)*outer_product3(princ(:, i), princ(:, i))
    end do

    call ShellTensorToStressVector(stress_tensor, stress_out)
    strain_out(1:3) = (/ logstrain(1, 1), logstrain(2, 2), logstrain(3, 3) /)
    strain_out(4:6) = 2.0D0*(/ logstrain(1, 2), logstrain(2, 3), logstrain(3, 1) /)

  end subroutine ShellMITC_TransformOutput

  pure subroutine ShellStressVectorToTensor(stress, tensor)
    real(kind=kreal), intent(in) :: stress(6)
    real(kind=kreal), intent(out) :: tensor(3, 3)

    tensor = 0.0D0
    tensor(1, 1) = stress(1)
    tensor(2, 2) = stress(2)
    tensor(3, 3) = stress(3)
    tensor(1, 2) = stress(4)
    tensor(2, 1) = stress(4)
    tensor(2, 3) = stress(5)
    tensor(3, 2) = stress(5)
    tensor(3, 1) = stress(6)
    tensor(1, 3) = stress(6)
  end subroutine ShellStressVectorToTensor

  pure subroutine ShellTensorToStressVector(tensor, stress)
    real(kind=kreal), intent(in) :: tensor(3, 3)
    real(kind=kreal), intent(out) :: stress(6)

    stress(1) = tensor(1, 1)
    stress(2) = tensor(2, 2)
    stress(3) = tensor(3, 3)
    stress(4) = tensor(1, 2)
    stress(5) = tensor(2, 3)
    stress(6) = tensor(3, 1)
  end subroutine ShellTensorToStressVector

  pure subroutine ShellObjectiveTraceStressIncrement(stress_old, dstrain, dstress_trace, trace_coeff)
    real(kind=kreal), intent(in) :: stress_old(6), dstrain(6)
    real(kind=kreal), intent(out) :: dstress_trace(6)
    real(kind=kreal), intent(in), optional :: trace_coeff

    real(kind=kreal) :: stress_tensor(3, 3), dstress_tensor(3, 3)
    real(kind=kreal) :: trD, coeff

    call ShellStressVectorToTensor(stress_old, stress_tensor)
    coeff = 1.0D0
    if (present(trace_coeff)) coeff = trace_coeff
    trD = coeff*(dstrain(1)+dstrain(2))+dstrain(3)
    dstress_tensor = -stress_tensor*trD
    call ShellTensorToStressVector(dstress_tensor, dstress_trace)
  end subroutine ShellObjectiveTraceStressIncrement

  pure subroutine ShellAddULObjectiveTraceTangent(stress_old, ncol, B, DB, trace_coeff)
    integer(kind=kint), intent(in) :: ncol
    real(kind=kreal), intent(in) :: stress_old(6), B(5, ncol)
    real(kind=kreal), intent(inout) :: DB(5, ncol)
    real(kind=kreal), intent(in), optional :: trace_coeff

    integer(kind=kint) :: j
    real(kind=kreal) :: dstrain_col(6), dstress_trace(6)

    do j = 1, ncol
      dstrain_col = 0.0D0
      dstrain_col(1:2) = B(1:2, j)
      call ShellObjectiveTraceStressIncrement(stress_old, dstrain_col, dstress_trace, trace_coeff)
      DB(1, j) = DB(1, j)+dstress_trace(1)
      DB(2, j) = DB(2, j)+dstress_trace(2)
      DB(3, j) = DB(3, j)+dstress_trace(4)
      DB(4, j) = DB(4, j)+dstress_trace(5)
      DB(5, j) = DB(5, j)+dstress_trace(6)
    end do
  end subroutine ShellAddULObjectiveTraceTangent

  pure subroutine ShellDirectorIncrementalSecondDeriv(director_current, director_second)
    real(kind=kreal), intent(in) :: director_current(3)
    real(kind=kreal), intent(out) :: director_second(3, 3, 3)

    integer :: m, n
    real(kind=kreal) :: basis_m(3), basis_n(3)
    real(kind=kreal) :: cross_n(3), cross_mn(3), cross_m(3), cross_nm(3)

    do n = 1, 3
      basis_n = 0.0D0
      basis_n(n) = 1.0D0
      cross_n = matmul(ShellSkewMatrix(basis_n), director_current)
      do m = 1, 3
        basis_m = 0.0D0
        basis_m(m) = 1.0D0
        cross_m = matmul(ShellSkewMatrix(basis_m), director_current)
        cross_mn = matmul(ShellSkewMatrix(basis_m), cross_n)
        cross_nm = matmul(ShellSkewMatrix(basis_n), cross_m)
        director_second(:, m, n) = 0.5D0*(cross_mn+cross_nm)
      end do
    end do
  end subroutine ShellDirectorIncrementalSecondDeriv


  subroutine ShellMITC_AbortNonlinearUnsupported(etype)
    integer(kind=kint), intent(in) :: etype

    write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
    write(*,*) ' ic_type = ', etype
    call hecmw_abort(hecmw_comm_get_comm())
  end subroutine ShellMITC_AbortNonlinearUnsupported

  !--------------------------------------------------------------------
  !> Integrate one physical shell layer for the material and geometric stiffness.
  !> Tying data are local to this call, so workspace size is independent of nlayer.
  subroutine ShellMITC_IntegrateStiffnessLayer(etype, nn, ndof, ilayer, kinematics, &
      finite_rotation, use_director_tangent, use_green_lagrange, &
      add_geometric_stiffness, ecoord, elem, shell_disp, gausses, element, &
      v1, v2, v3, director, reference_director, director_tangent, &
      director_second_tangent, nodal_kinematic_dofs, nddrill, stiff_work)
    use mMechGauss
    use m_MatMatrix
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof, ilayer, kinematics
    logical, intent(in) :: finite_rotation, use_director_tangent
    logical, intent(in) :: use_green_lagrange, add_geometric_stiffness
    real(kind=kreal), intent(in) :: ecoord(3, nn), elem(3, nn), shell_disp(6, nn)
    type(tGaussStatus), intent(in) :: gausses(:)
    type(tElement), intent(in), optional :: element
    real(kind=kreal), intent(in) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal), intent(in) :: director(3, nn), reference_director(3, nn)
    real(kind=kreal), intent(in) :: director_tangent(3, 3, nn)
    real(kind=kreal), intent(in) :: director_second_tangent(3, 3, 3, nn)
    real(kind=kreal), intent(in) :: nodal_kinematic_dofs(ndof*nn)
    real(kind=kreal), intent(in), optional :: nddrill(nn)
    real(kind=kreal), intent(inout) :: stiff_work(ndof*nn, ndof*nn)

    integer :: lx, ly, ny, isize, jsize, ishell
    integer(kind=kint) :: ierr
    real(kind=kreal) :: zeta_ly, tying_zeta, xi_lx, eta_lx
    real(kind=kreal) :: surface_weight, thickness_weight, integration_jacobian
    real(kind=kreal) :: integration_weight, layer_weight
    real(kind=kreal) :: naturalcoord(2), shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal) :: D(5, 5), B(5, ndof*nn), DB(5, ndof*nn)
    real(kind=kreal) :: basis_variation(3, ndof*nn, 3)
    real(kind=kreal) :: stress_vec(5), stress_old_vec(6), alpha, trace_coeff
    real(kind=kreal) :: point_triad(3, 3), covariant_basis(3, 3)
    real(kind=kreal) :: tangent_basis(3, 2), reciprocal_basis(3, 3), local_basis(3, 3)
    real(kind=kreal) :: material_reciprocal_basis(3, 3), material_local_basis(3, 3)
    real(kind=kreal) :: director_contribution(3, nn, 3)
    real(kind=kreal) :: director_second_variation(5, 3, 3, nn)
    real(kind=kreal) :: tying_B(5, ndof*nn, 6, 3)
    real(kind=kreal) :: tying_basis_variation(3, ndof*nn, 3, 6, 3)
    real(kind=kreal) :: tying_director_second_variation(5, 3, 3, nn, 6, 3)
    real(kind=kreal) :: Cv(ndof*nn), Cv_disp, Cv_w_second(3, 3, nn)

    ny = NumOfShellThicknessQuadPoints(etype)
    layer_weight = gausses(1)%pMaterial%shell_var(ilayer)%weight
    stress_vec = 0.0D0
    director_second_variation = 0.0D0

    ! MITC3/MITC4 use the same midsurface tying data for every thickness point.
    ! Re-evaluating it once per layer keeps the workspace layer-local by design.
    if( etype /= fe_mitc9_shell ) then
      tying_zeta = ShellMITC_TyingZeta(etype, 0.0D0)
      call ShellMITC_EvaluateTyingStiffnessDataAtZeta(etype, nn, ndof, tying_zeta, &
        elem, shell_disp, director, director_tangent, director_second_tangent, &
        use_green_lagrange, use_director_tangent, tying_B, tying_basis_variation, tying_director_second_variation)
    endif

    do ly = 1, ny
      call fstr_shell_layer_quadrature_gauss(etype, gausses(1), ilayer, ly, zeta_ly, thickness_weight, ierr)
      if( ierr /= 0 ) cycle

      ! MITC9 tying data depend on the physical thickness coordinate.
      if( etype == fe_mitc9_shell ) then
        tying_zeta = ShellMITC_TyingZeta(etype, zeta_ly)
        call ShellMITC_EvaluateTyingStiffnessDataAtZeta(etype, nn, ndof, tying_zeta, &
          elem, shell_disp, director, director_tangent, director_second_tangent, &
          use_green_lagrange, use_director_tangent, tying_B, tying_basis_variation, tying_director_second_variation)
      endif

      do lx = 1, NumOfQuadPoints(etype)
        call getQuadPoint(etype, lx, naturalcoord)
        xi_lx = naturalcoord(SHELL_XI)
        eta_lx = naturalcoord(SHELL_ETA)
        surface_weight = getWeight(etype, lx)
        call getShapeFunc(etype, naturalcoord, shapefunc)
        call getShapeDeriv(etype, naturalcoord, shapederiv)

        point_triad(:, SHELL_XI) = matmul(v1, shapefunc)
        point_triad(:, SHELL_ETA) = matmul(v2, shapefunc)
        point_triad(:, SHELL_ZETA) = matmul(v3, shapefunc)

        call ShellMITC_PreparePointKinematics(etype, nn, use_green_lagrange, ecoord, elem, &
          shell_disp(1:3, :), director, reference_director, zeta_ly, shapefunc, shapederiv, &
          covariant_basis, tangent_basis, reciprocal_basis, local_basis, &
          material_reciprocal_basis, material_local_basis, integration_jacobian, director_contribution)

        if( add_geometric_stiffness .and. present(element) ) then
          ishell = fstr_shell_layer_gauss_index(element, lx, ilayer, ly)
        else
          ishell = 0
        endif

        if( ishell > 0 ) then
          call MatlMatrix_Shell(element%shell_layer_gausses(ishell), Shell, D, &
            material_local_basis(:, SHELL_XI), material_local_basis(:, SHELL_ETA), material_local_basis(:, SHELL_ZETA), &
            material_reciprocal_basis(:, SHELL_XI), material_reciprocal_basis(:, SHELL_ETA), &
            material_reciprocal_basis(:, SHELL_ZETA), alpha, ilayer)
        else
          call MatlMatrix_Shell(gausses(lx), Shell, D, &
            material_local_basis(:, SHELL_XI), material_local_basis(:, SHELL_ETA), material_local_basis(:, SHELL_ZETA), &
            material_reciprocal_basis(:, SHELL_XI), material_reciprocal_basis(:, SHELL_ETA), &
            material_reciprocal_basis(:, SHELL_ZETA), alpha, ilayer)
        endif

        call ShellMITC_BuildFirstStrainVariation(nn, ndof, zeta_ly, shapefunc, shapederiv, &
          covariant_basis, tangent_basis, director_contribution, director_tangent, use_director_tangent, B, basis_variation)
        if( use_director_tangent ) then
          call ShellMITC_BuildDirectorSecondVariation(nn, zeta_ly, shapefunc, shapederiv, &
            covariant_basis, director_second_tangent, director_second_variation)
        endif

        if( use_director_tangent ) then
          call ShellMITC_ApplyAssumedStrain(etype, nn, ndof, xi_lx, eta_lx, tying_B, B, &
            tying_director_second_variation, director_second_variation)
        else
          call ShellMITC_ApplyAssumedStrain(etype, nn, ndof, xi_lx, eta_lx, tying_B, B)
        endif

        integration_weight = surface_weight*thickness_weight*integration_jacobian

        if( add_geometric_stiffness ) then
          if( ishell > 0 ) then
            stress_vec = (/ element%shell_layer_gausses(ishell)%stress(1), element%shell_layer_gausses(ishell)%stress(2), &
              element%shell_layer_gausses(ishell)%stress(4), element%shell_layer_gausses(ishell)%stress(5), &
              element%shell_layer_gausses(ishell)%stress(6) /)
          else
            stress_vec = (/ gausses(lx)%stress(1), gausses(lx)%stress(2), &
              gausses(lx)%stress(4), gausses(lx)%stress(5), gausses(lx)%stress(6) /)
          endif
        endif

        DB = matmul(D, B)
        if( kinematics == UPDATELAG .and. etype == fe_mitc4_shell .and. nn == 4 ) then
          if( ishell > 0 ) then
            stress_old_vec = element%shell_layer_gausses(ishell)%stress_bak(1:6)
            trace_coeff = ShellPlaneStressTraceCoeff(element%shell_layer_gausses(ishell), ilayer)
          else
            stress_old_vec = gausses(lx)%stress_bak(1:6)
            trace_coeff = ShellPlaneStressTraceCoeff(gausses(lx), ilayer)
          endif
          call ShellAddULObjectiveTraceTangent(stress_old_vec, ndof*nn, B, DB, trace_coeff=trace_coeff)
        endif

        do jsize = 1, ndof*nn
          do isize = 1, ndof*nn
            stiff_work(isize, jsize) = stiff_work(isize, jsize) &
              +integration_weight*layer_weight*dot_product(B(:, isize), DB(:, jsize))
          end do
        end do

        call ShellMITC_AddGeometricStiffness(etype, nn, ndof, kinematics, &
          use_director_tangent, use_green_lagrange, add_geometric_stiffness, &
          xi_lx, eta_lx, integration_weight, layer_weight, stress_vec, B, &
          basis_variation, tying_basis_variation, reciprocal_basis, material_local_basis, &
          material_reciprocal_basis, director_second_variation, director, stiff_work)

        call ShellMITC_BuildDrillingSecondVariation(nn, zeta_ly, shapefunc, shapederiv, &
          director_second_tangent, material_reciprocal_basis, point_triad, use_director_tangent, Cv_w_second)
        call ShellMITC_BuildDrillingVector(nn, ndof, finite_rotation, shapefunc, &
          point_triad, material_reciprocal_basis, basis_variation, nodal_kinematic_dofs, director, nddrill, Cv, Cv_disp)
        call ShellMITC_AddDrillingStiffness(nn, ndof, finite_rotation, &
          use_director_tangent, shapefunc, Cv, Cv_disp, Cv_w_second, nodal_kinematic_dofs, &
          director, director_tangent, alpha, integration_weight, layer_weight, nddrill, stiff_work)
      end do
    end do
  end subroutine ShellMITC_IntegrateStiffnessLayer

  !--------------------------------------------------------------------
  !> Integrate one physical shell layer for stress update and internal force.
  !> Only one zeta's tying B matrices are retained at a time.
  subroutine ShellMITC_IntegrateInternalForceLayer(etype, nn, ndof, ilayer, kinematics, &
      finite_rotation, use_director_tangent, use_green_lagrange, update_state, &
      ecoord, evaluation_coords, total_nodal_state, strain_nodal_state, gausses, element, &
      v1, v2, v3, director, reference_director, director_tangent, stress_elem, &
      stress_director, stress_director_increment, nodal_kinematic_dofs, nddrill, qf_work)
    use mMechGauss
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof, ilayer, kinematics
    logical, intent(in) :: finite_rotation, use_director_tangent
    logical, intent(in) :: use_green_lagrange, update_state
    real(kind=kreal), intent(in) :: ecoord(3, nn), evaluation_coords(3, nn)
    real(kind=kreal), intent(in) :: total_nodal_state(6, nn), strain_nodal_state(6, nn)
    type(tGaussStatus), intent(in) :: gausses(:)
    type(tElement), intent(inout), optional :: element
    real(kind=kreal), intent(in) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal), intent(in) :: director(3, nn), reference_director(3, nn)
    real(kind=kreal), intent(in) :: director_tangent(3, 3, nn)
    real(kind=kreal), intent(in) :: stress_elem(3, nn), stress_director(3, nn)
    real(kind=kreal), intent(in) :: stress_director_increment(3, nn)
    real(kind=kreal), intent(in) :: nodal_kinematic_dofs(ndof*nn)
    real(kind=kreal), intent(in), optional :: nddrill(nn)
    real(kind=kreal), intent(inout) :: qf_work(ndof*nn)

    integer :: lx, ly, ny, ishell
    integer(kind=kint) :: ierr
    logical :: store_state
    type(tGaussStatus), target :: gauss_work
    type(tGaussStatus), pointer :: gauss
    real(kind=kreal) :: zeta_ly, tying_zeta, xi, eta
    real(kind=kreal) :: surface_weight, thickness_weight, integration_jacobian
    real(kind=kreal) :: integration_weight, layer_weight
    real(kind=kreal) :: naturalcoord(2), shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal) :: B(5, ndof*nn), strain_vec(5), stress_vec(5), alpha
    real(kind=kreal) :: basis_variation(3, ndof*nn, 3), point_triad(3, 3)
    real(kind=kreal) :: covariant_basis(3, 3), tangent_basis(3, 2)
    real(kind=kreal) :: reciprocal_basis(3, 3), local_basis(3, 3)
    real(kind=kreal) :: material_reciprocal_basis(3, 3), material_local_basis(3, 3)
    real(kind=kreal) :: director_contribution(3, nn, 3)
    real(kind=kreal) :: tying_B(5, ndof*nn, 6, 3), tying_strain(5, 6, 3)
    real(kind=kreal) :: dstrain(6), stress(6), strain_out(6), stress_out(6)
    real(kind=kreal) :: strain_tensor(3, 3), stress_tensor(3, 3)
    real(kind=kreal) :: stress_covariant_basis(3, 3), stress_reciprocal_basis(3, 3)
    real(kind=kreal) :: stress_material_local_basis(3, 3)
    real(kind=kreal) :: stress_material_reciprocal_basis(3, 3)
    real(kind=kreal) :: reference_basis(3, 3), current_basis(3, 3)
    real(kind=kreal) :: reference_jacobian, current_jacobian
    real(kind=kreal) :: Cv(ndof*nn), Cv_disp

    ny = NumOfShellThicknessQuadPoints(etype)
    layer_weight = gausses(1)%pMaterial%shell_var(ilayer)%weight
    ! MITC3/MITC4 tying is evaluated at zeta=0 once for this layer.
    if( etype /= fe_mitc9_shell ) then
      tying_zeta = ShellMITC_TyingZeta(etype, 0.0D0)
      call ShellMITC_EvaluateTyingBAtZeta(etype, nn, ndof, tying_zeta, evaluation_coords, &
        total_nodal_state, director, director_tangent, use_green_lagrange, use_director_tangent, tying_B)
    endif

    do ly = 1, ny
      call fstr_shell_layer_quadrature_gauss(etype, gausses(1), ilayer, ly, zeta_ly, thickness_weight, ierr)
      if( ierr /= 0 ) cycle

      if( etype == fe_mitc9_shell ) then
        tying_zeta = ShellMITC_TyingZeta(etype, zeta_ly)
        call ShellMITC_EvaluateTyingBAtZeta(etype, nn, ndof, tying_zeta, evaluation_coords, &
          total_nodal_state, director, director_tangent, use_green_lagrange, use_director_tangent, tying_B)
      endif

      if( update_state ) then
        call ShellMITC_EvaluateTyingStrains(etype, nn, zeta_ly, stress_elem, &
          strain_nodal_state, stress_director, stress_director_increment, use_green_lagrange, tying_strain)
      endif

      do lx = 1, NumOfQuadPoints(etype)
        call getQuadPoint(etype, lx, naturalcoord)
        xi = naturalcoord(SHELL_XI)
        eta = naturalcoord(SHELL_ETA)
        surface_weight = getWeight(etype, lx)
        call getShapeFunc(etype, naturalcoord, shapefunc)
        call getShapeDeriv(etype, naturalcoord, shapederiv)

        ishell = 0
        store_state = .false.
        if( update_state .and. present(element) ) then
          ishell = fstr_shell_layer_gauss_index(element, lx, ilayer, ly)
          store_state = ishell > 0
        endif
        if( update_state .and. .not. store_state ) then
          stop "Missing shell layer Gauss state"
        endif
        if( store_state ) then
          gauss => element%shell_layer_gausses(ishell)
        else
          gauss_work = gausses(lx)
          gauss => gauss_work
        endif

        if( update_state ) then
          call ShellMITC_EvaluateAssumedStrain(etype, nn, naturalcoord, zeta_ly, stress_elem, &
            strain_nodal_state, stress_director, stress_director_increment, &
            use_green_lagrange, tying_strain, dstrain, strain_tensor, stress_covariant_basis, stress_reciprocal_basis, &
            stress_material_local_basis, stress_material_reciprocal_basis, &
            reference_basis, current_basis, reference_jacobian, current_jacobian)
        endif

        point_triad(:, SHELL_XI) = matmul(v1, shapefunc)
        point_triad(:, SHELL_ETA) = matmul(v2, shapefunc)
        point_triad(:, SHELL_ZETA) = matmul(v3, shapefunc)
        call ShellMITC_PreparePointKinematics(etype, nn, use_green_lagrange, ecoord, &
          evaluation_coords, total_nodal_state(1:3, :), director, reference_director, &
          zeta_ly, shapefunc, shapederiv, covariant_basis, tangent_basis, reciprocal_basis, &
          local_basis, material_reciprocal_basis, material_local_basis, integration_jacobian, director_contribution)

        call ShellMITC_BuildFirstStrainVariation(nn, ndof, zeta_ly, shapefunc, shapederiv, &
          covariant_basis, tangent_basis, director_contribution, director_tangent, use_director_tangent, B, basis_variation)
        call ShellMITC_ApplyAssumedStrain(etype, nn, ndof, xi, eta, tying_B, B)

        if( .not. update_state ) then
          strain_vec = matmul(B, nodal_kinematic_dofs)
          dstrain = (/ strain_vec(1), strain_vec(2), 0.0D0, strain_vec(3), strain_vec(4), strain_vec(5) /)
          stress_material_local_basis = material_local_basis
          stress_material_reciprocal_basis = material_reciprocal_basis
        endif

        call ShellMITC_UpdateStress(kinematics, store_state, gauss, ilayer, dstrain, &
          stress_material_local_basis, stress_material_reciprocal_basis, stress, alpha)

        if( store_state ) then
          if( kinematics == UPDATELAG ) then
            gauss%strain_out(1:6) = gauss%strain(1:6)
            gauss%stress_out(1:6) = gauss%stress(1:6)
          else
            call ShellStressVectorToTensor(stress, stress_tensor)
            call ShellMITC_TransformOutput(use_green_lagrange, stress_tensor, strain_tensor, &
              stress_covariant_basis, stress_reciprocal_basis, reference_basis, &
              current_basis, reference_jacobian, current_jacobian, strain_out, stress_out)
            gauss%strain_out(1:6) = strain_out
            gauss%stress_out(1:6) = stress_out
          endif
        endif

        stress_vec = (/ stress(1), stress(2), stress(4), stress(5), stress(6) /)
        integration_weight = surface_weight*thickness_weight*integration_jacobian
        qf_work = qf_work+integration_weight*layer_weight*matmul(stress_vec, B)

        call ShellMITC_BuildDrillingVector(nn, ndof, finite_rotation, shapefunc, &
          point_triad, material_reciprocal_basis, basis_variation, nodal_kinematic_dofs, director, nddrill, Cv, Cv_disp)
        qf_work = qf_work+integration_weight*layer_weight*alpha*Cv*Cv_disp
      end do
    end do
  end subroutine ShellMITC_IntegrateInternalForceLayer

  !=====================================================================
  ! Public shell element entry points
  !=====================================================================

  !--------------------------------------------------------------------
  !> Calculate the tangent stiffness matrix of a MITC shell element.
  subroutine STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, &
      nddisp, element, ndtriad, ndreftriad, ndcurtriad, nddrill)

    use mMechGauss
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, mixflag, ndof
    real(kind=kreal), intent(in) :: ecoord(3, nn), thick
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: stiff(:, :)
    real(kind=kreal), intent(in), optional :: nddisp(ndof, nn)
    type(tElement), intent(in), optional :: element
    !> Nodal frames (triads), packed as e1(1:3), e2(4:6), e3=director axis(7:9).
    !> The director itself (0.5*thick*e3) is computed inside this routine.
    real(kind=kreal), intent(in), optional :: ndtriad(9, nn)
    real(kind=kreal), intent(in), optional :: ndreftriad(9, nn), ndcurtriad(9, nn)
    real(kind=kreal), intent(in), optional :: nddrill(nn)

    integer :: ndof_shell, nb, ilayer, nlayer
    integer(kind=kint) :: kinematics
    logical :: finite_rotation, use_director_tangent, use_green_lagrange
    logical :: add_geometric_stiffness, update_state
    real(kind=kreal) :: elem(3, nn), shell_disp(6, nn)
    real(kind=kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal) :: director(3, nn), reference_director(3, nn)
    real(kind=kreal) :: director_tangent(3, 3, nn)
    real(kind=kreal) :: director_second_tangent(3, 3, 3, nn)
    real(kind=kreal) :: nodal_kinematic_dofs(ndof*nn)
    real(kind=kreal) :: stiff_work(ndof*nn, ndof*nn)

    shell_disp = 0.0D0
    director_second_tangent = 0.0D0
    ndof_shell = min(ndof, 6)
    if( present(nddisp) ) shell_disp(1:ndof_shell, :) = nddisp(1:ndof_shell, :)

    call ShellMITC_ResolveFormulation(etype, nn, ndof, gausses(1), present(nddisp), &
      present(element), .true., kinematics, ndof_shell, finite_rotation, &
      use_director_tangent, use_green_lagrange, add_geometric_stiffness, update_state)

    call ShellMITC_PrepareNodalKinematics(etype, nn, thick, kinematics, finite_rotation, &
      use_director_tangent, use_director_tangent, ecoord, shell_disp, ndtriad, ndreftriad, &
      ndcurtriad, elem, v1, v2, v3, director, reference_director, director_tangent, director_second_tangent)

    stiff = 0.0D0
    stiff_work = 0.0D0
    nodal_kinematic_dofs = 0.0D0
    do nb = 1, nn
      nodal_kinematic_dofs(ndof*(nb-1)+1:ndof*(nb-1)+ndof_shell) = shell_disp(1:ndof_shell, nb)
    end do

    nlayer = gausses(1)%pMaterial%totallyr
    do ilayer = 1, nlayer
      call ShellMITC_IntegrateStiffnessLayer(etype, nn, ndof, ilayer, kinematics, &
        finite_rotation, use_director_tangent, use_green_lagrange, &
        add_geometric_stiffness, ecoord, elem, shell_disp, gausses, element, &
        v1, v2, v3, director, reference_director, director_tangent, &
        director_second_tangent, nodal_kinematic_dofs, nddrill, stiff_work)
    end do

    stiff(1:nn*ndof, 1:nn*ndof) = stiff_work
    call ShellApplyMixedDofOrdering(mixflag, nn, ndof, stiff=stiff)
  end subroutine STF_Shell_MITC

  !--------------------------------------------------------------------
  !> Evaluate MITC shell stress and strain for result output.
  subroutine ElementStress_Shell_MITC(etype, nn, ndof, ecoord, gausses, edisp, &
      strain, stress, thick, zeta, n_layer, surface_gauss_points, &
      local_strain, local_stress, local_stress_override, ndtriad, ndreftriad, ndbase_disp)
    use mMechGauss
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof
    integer, intent(in) :: n_layer
    real(kind=kreal), intent(in) :: ecoord(3, nn), edisp(6, nn), thick, zeta
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: strain(:, :), stress(:, :)
    logical, intent(in), optional :: surface_gauss_points
    real(kind=kreal), intent(out), optional :: local_strain(:, :), local_stress(:, :)
    real(kind=kreal), intent(in), optional :: local_stress_override(:, :)
    !> Nodal frames (triads), packed as e1(1:3), e2(4:6), e3=director axis(7:9).
    !> The director itself (0.5*thick*e3) is computed inside ShellMITC_SetupNodalDirectors.
    real(kind=kreal), intent(in), optional :: ndtriad(9, nn), ndreftriad(9, nn)
    real(kind=kreal), intent(in), optional :: ndbase_disp(6, nn)

    integer :: lx, npoints
    integer(kind=kint) :: ierr_quad, kinematics, ndof_shell
    logical :: finite_rotation, use_director_tangent, use_green_lagrange
    logical :: add_geometric_stiffness, update_state, use_surface_gauss
    type(tGaussStatus), target :: gauss_work
    real(kind=kreal) :: elem(3, nn), naturalcoord(2), nncoord(nn, 2), zeta_ly
    real(kind=kreal) :: director(3, nn), director_increment(3, nn)
    real(kind=kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal) :: a_over_2_v3(3, nn), a_over_2_v3_ref(3, nn)
    real(kind=kreal) :: a_over_2_v3_deriv(3, 3, nn), a_over_2_v3_second(3, 3, 3, nn)
    real(kind=kreal) :: tying_strain(5, 6, 3)
    real(kind=kreal) :: point_strain(6), point_stress(6), alpha
    real(kind=kreal) :: strain_tensor(3, 3), stress_tensor(3, 3)
    real(kind=kreal) :: covariant_basis(3, 3), reciprocal_basis(3, 3)
    real(kind=kreal) :: material_local_basis(3, 3), material_reciprocal_basis(3, 3)
    real(kind=kreal) :: reference_basis(3, 3), current_basis(3, 3)
    real(kind=kreal) :: reference_jacobian, current_jacobian

    use_surface_gauss = .false.
    if( present(surface_gauss_points) ) use_surface_gauss = surface_gauss_points

    ! Reproduce the same stress-evaluation configuration used by UPDATE.
    call ShellMITC_ResolveFormulation(etype, nn, ndof, gausses(1), .true., .false., &
      .false., kinematics, ndof_shell, finite_rotation, use_director_tangent, &
      use_green_lagrange, add_geometric_stiffness, update_state)
    call getNodalNaturalCoord(etype, nncoord)
    elem = ecoord
    if( kinematics == UPDATELAG .and. present(ndbase_disp) ) elem = elem+ndbase_disp(1:3, :)
    if( kinematics == UPDATELAG ) elem = elem+0.5D0*edisp(1:3, :)
    call ShellMITC_SetupNodalDirectors(etype, nn, thick, kinematics, elem, edisp, &
      finite_rotation, use_director_tangent, .false., ndtriad, ndreftriad, &
      v1=v1, v2=v2, v3=v3, a_over_2_v3=a_over_2_v3, a_over_2_v3_ref=a_over_2_v3_ref, &
      a_over_2_v3_deriv=a_over_2_v3_deriv, a_over_2_v3_second=a_over_2_v3_second)

    call ShellMITC_SetupStressKinematics(nn, ecoord, kinematics, finite_rotation, edisp, &
      a_over_2_v3, a_over_2_v3_ref, a_over_2_v3_deriv, ndbase_disp, elem, director, director_increment)

    call ShellMITC_EvaluateTyingStrains(etype, nn, zeta, elem, edisp, director, &
      director_increment, use_green_lagrange, tying_strain)
    call fstr_shell_layer_zeta(gausses(1), n_layer, zeta, zeta_ly, ierr_quad)
    if( ierr_quad /= 0 ) stop "Invalid shell layer zeta"

    npoints = nn
    if( use_surface_gauss ) npoints = NumOfQuadPoints(etype)
    do lx = 1, npoints
      if( use_surface_gauss ) then
        call getQuadPoint(etype, lx, naturalcoord)
      else
        naturalcoord = nncoord(lx, :)
      endif

      call ShellMITC_EvaluateAssumedStrain(etype, nn, naturalcoord, zeta_ly, elem, edisp, &
        director, director_increment, use_green_lagrange, tying_strain, point_strain, &
        strain_tensor, covariant_basis, reciprocal_basis, material_local_basis, &
        material_reciprocal_basis, reference_basis, current_basis, reference_jacobian, current_jacobian)

      if( present(local_stress_override) .and. lx <= size(local_stress_override, 1) &
          .and. size(local_stress_override, 2) >= 6 ) then
        point_stress = (/ local_stress_override(lx, 1), local_stress_override(lx, 2), 0.0D0, local_stress_override(lx, 4), &
          local_stress_override(lx, 5), local_stress_override(lx, 6) /)
      else
        gauss_work = gausses(lx)
        call ShellMITC_UpdateStress(kinematics, update_state, gauss_work, n_layer, &
          point_strain, material_local_basis, material_reciprocal_basis, point_stress, alpha)
      endif

      call ShellStressVectorToTensor(point_stress, stress_tensor)
      call ShellMITC_TransformOutput(use_green_lagrange, stress_tensor, strain_tensor, &
        covariant_basis, reciprocal_basis, reference_basis, current_basis, &
        reference_jacobian, current_jacobian, strain(lx, 1:6), stress(lx, 1:6))
      if( present(local_strain) ) local_strain(lx, 1:6) = point_strain
      if( present(local_stress) ) local_stress(lx, 1:6) = point_stress
    end do

  end subroutine ElementStress_Shell_MITC

  !--------------------------------------------------------------------
  !> Calculate the distributed load vector of a MITC shell element.
  subroutine DL_Shell(etype, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
    use mMechGauss
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof
    real(kind=kreal), intent(in) :: xx(*), yy(*), zz(*), rho, thick, params(*)
    integer, intent(in) :: ltype
    real(kind=kreal), intent(out) :: vect(*)
    integer(kind=kint), intent(out) :: nsize
    type(tGaussStatus), intent(in) :: gausses(:)

    integer :: ny, lx, ly, nb, ilayer, nlayer, jsize
    integer(kind=kint) :: ierr_quad
    real(kind=kreal) :: elem(3, nn)
    real(kind=kreal) :: v1(3, nn), v2(3, nn), v3(3, nn), director(3, nn)
    real(kind=kreal) :: naturalcoord(2), shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal) :: covariant_basis(3, 3), normal_jac(3), N(3, ndof*nn)
    real(kind=kreal) :: body_force(3), position(3), origin(3), axis(3), radial(3)
    real(kind=kreal) :: val, zeta, surface_weight, thickness_weight, integration_weight, projection

    nsize = ndof*nn
    vect(1:nsize) = 0.0D0
    elem(1, :) = xx(1:nn)
    elem(2, :) = yy(1:nn)
    elem(3, :) = zz(1:nn)
    val = params(1)
    origin = params(2:4)
    axis = params(5:7)
    ny = NumOfShellThicknessQuadPoints(etype)

    if (ltype >= 10) then
      do lx = 1, NumOfQuadPoints(etype)
        call getQuadPoint(etype, lx, naturalcoord)
        surface_weight = getWeight(etype, lx)
        call getShapeFunc(etype, naturalcoord, shapefunc)
        call getShapeDeriv(etype, naturalcoord, shapederiv)
        covariant_basis(:, SHELL_XI:SHELL_ETA) = matmul(elem, shapederiv)
        call cross_product(covariant_basis(:, SHELL_XI), covariant_basis(:, SHELL_ETA), normal_jac)
        do nb = 1, nn
          jsize = ndof*(nb-1)
          vect(jsize+1:jsize+3) = vect(jsize+1:jsize+3) + surface_weight*shapefunc(nb)*normal_jac*val
        end do
      end do
    else
      call ShellMITC_SetupReferenceDirectors(etype, nn, thick, elem, v1, v2, v3, director)
      nlayer = gausses(1)%pMaterial%totallyr
      do ilayer = 1, nlayer
        do ly = 1, ny
          call fstr_shell_layer_quadrature_gauss(etype, gausses(1), ilayer, ly, zeta, thickness_weight, ierr_quad)
          if (ierr_quad /= 0) cycle
          do lx = 1, NumOfQuadPoints(etype)
            call getQuadPoint(etype, lx, naturalcoord)
            surface_weight = getWeight(etype, lx)
            call getShapeFunc(etype, naturalcoord, shapefunc)
            call getShapeDeriv(etype, naturalcoord, shapederiv)
            call ShellMITC_CovariantBasis(nn, elem, director, zeta, shapefunc, shapederiv, covariant_basis)
            integration_weight = surface_weight*thickness_weight *ShellMITC_CovariantJacobian(covariant_basis)
            call ShellMITC_BuildInterpolationMatrix(nn, ndof, zeta, shapefunc, director, N)

            body_force = 0.0D0
            select case (ltype)
            case (1)
              body_force(1) = val
            case (2)
              body_force(2) = val
            case (3)
              body_force(3) = val
            case (4)
              body_force = rho*val*origin
            case (5)
              position = matmul(elem, shapefunc)
              projection = dot_product(position-origin, axis)/dot_product(axis, axis)
              radial = position-(origin+projection*axis)
              body_force = rho*val*val*radial
            end select
            vect(1:nsize) = vect(1:nsize)+integration_weight*matmul(body_force, N)
          end do
        end do
      end do
    endif
  end subroutine DL_Shell

  !--------------------------------------------------------------------
  !> Adapt the natural shell load ordering to shell-solid mixed elements.
  subroutine DL_Shell_33(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
    use mMechGauss
    implicit none

    integer(kind=kint), intent(in) :: ic_type, nn, ndof
    real(kind=kreal), intent(in) :: xx(*), yy(*), zz(*), rho, thick, params(*)
    integer, intent(in) :: ltype
    real(kind=kreal), intent(out) :: vect(*)
    integer(kind=kint), intent(out) :: nsize
    type(tGaussStatus), intent(in) :: gausses(:)

    select case (ic_type)
    case (761)
      call DL_Shell(fe_mitc3_shell, 3_kint, 6_kint, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
      call ShellApplyMixedDofOrdering(2_kint, 3_kint, 6_kint, qf=vect(1:18))
    case (781)
      call DL_Shell(fe_mitc4_shell, 4_kint, 6_kint, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
      call ShellApplyMixedDofOrdering(1_kint, 4_kint, 6_kint, qf=vect(1:24))
    case default
      nsize = nn*ndof
      vect(1:nsize) = 0.0D0
    end select
  end subroutine DL_Shell_33

  !> Update shell stress and assemble the equivalent nodal force.
  subroutine UPDATE_Shell_MITC(etype, nn, ndof, ecoord, u, du, gausses, qf, thick, mixflag, &
      nddisp, element, ndtriad, ndreftriad, ndcurtriad, nddrill)

    use mMechGauss
    use mMaterial, only: UPDATELAG
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof, mixflag
    real(kind=kreal), intent(in) :: ecoord(3, nn), u(:, :), du(:, :), thick
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: qf(:)
    real(kind=kreal), intent(in), optional :: nddisp(ndof, nn)
    type(tElement), intent(inout), optional :: element
    !> Nodal frames (triads), packed as e1(1:3), e2(4:6), e3=director axis(7:9).
    !> The director itself (0.5*thick*e3) is computed inside this routine.
    real(kind=kreal), intent(in), optional :: ndtriad(9, nn), ndreftriad(9, nn)
    real(kind=kreal), intent(in), optional :: ndcurtriad(9, nn), nddrill(nn)

    integer :: ndof_shell, nb, ilayer, nlayer
    integer(kind=kint) :: kinematics
    logical :: finite_rotation, use_director_tangent, use_green_lagrange
    logical :: add_geometric_stiffness, update_state
    real(kind=kreal) :: total_nodal_state(6, nn), strain_nodal_state(6, nn)
    real(kind=kreal) :: step_base_nodal_state(6, nn)
    real(kind=kreal) :: nodal_kinematic_dofs(ndof*nn), qf_work(ndof*nn)
    real(kind=kreal) :: evaluation_coords(3, nn), stress_elem(3, nn)
    real(kind=kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind=kreal) :: director(3, nn), reference_director(3, nn)
    real(kind=kreal) :: director_tangent(3, 3, nn)
    real(kind=kreal) :: director_second_tangent(3, 3, 3, nn)
    real(kind=kreal) :: stress_director(3, nn), stress_director_increment(3, nn)

    call ShellMITC_ResolveFormulation(etype, nn, ndof, gausses(1), .true., present(element), &
      .true., kinematics, ndof_shell, finite_rotation, use_director_tangent, &
      use_green_lagrange, add_geometric_stiffness, update_state)

    total_nodal_state = 0.0D0
    strain_nodal_state = 0.0D0
    step_base_nodal_state = 0.0D0
    do nb = 1, nn
      total_nodal_state(1:ndof_shell, nb) = u(1:ndof_shell, nb)+du(1:ndof_shell, nb)
      strain_nodal_state(1:ndof_shell, nb) = total_nodal_state(1:ndof_shell, nb)
      step_base_nodal_state(1:ndof_shell, nb) = u(1:ndof_shell, nb)
      if( finite_rotation ) then
        call ShellComposeRotationVector(u(4:6, nb), du(4:6, nb), total_nodal_state(4:6, nb))
        strain_nodal_state(4:6, nb) = total_nodal_state(4:6, nb)
      endif
    end do
    if( present(nddisp) ) total_nodal_state(1:ndof_shell, :) = nddisp(1:ndof_shell, :)
    if( kinematics == UPDATELAG .and. update_state ) then
      strain_nodal_state = du(1:6, 1:nn)
    endif

    nodal_kinematic_dofs = 0.0D0
    do nb = 1, nn
      nodal_kinematic_dofs(ndof*(nb-1)+1:ndof*(nb-1)+ndof_shell) = total_nodal_state(1:ndof_shell, nb)
    end do
    qf_work = 0.0D0
    director_second_tangent = 0.0D0
    stress_elem = 0.0D0
    stress_director = 0.0D0
    stress_director_increment = 0.0D0

    call ShellMITC_PrepareNodalKinematics(etype, nn, thick, kinematics, finite_rotation, &
      use_director_tangent, .false., ecoord, total_nodal_state, ndtriad, ndreftriad, &
      ndcurtriad, evaluation_coords, v1, v2, v3, director, reference_director, director_tangent, director_second_tangent)

    if( update_state ) then
      call ShellMITC_SetupStressKinematics(nn, ecoord, kinematics, finite_rotation, &
        strain_nodal_state, director, reference_director, director_tangent, &
        step_base_nodal_state, stress_elem, stress_director, stress_director_increment)
    endif

    nlayer = gausses(1)%pMaterial%totallyr
    do ilayer = 1, nlayer
      call ShellMITC_IntegrateInternalForceLayer(etype, nn, ndof, ilayer, kinematics, &
        finite_rotation, use_director_tangent, use_green_lagrange, update_state, &
        ecoord, evaluation_coords, total_nodal_state, strain_nodal_state, gausses, element, &
        v1, v2, v3, director, reference_director, director_tangent, stress_elem, &
        stress_director, stress_director_increment, nodal_kinematic_dofs, nddrill, qf_work)
    end do

    qf(1:ndof*nn) = qf_work
    call ShellApplyMixedDofOrdering(mixflag, nn, ndof, qf=qf(1:ndof*nn))
  end subroutine UPDATE_Shell_MITC

  !> Linear shell adapter for the split translational/rotational node layout.
  subroutine UPDATE_Shell_MITC33(etype, nn, ndof, ecoord, u, du, gausses, qf, thick, mixflag, nddisp)

    use mMechGauss
    implicit none

    integer(kind=kint), intent(in) :: etype, nn, ndof, mixflag
    real(kind=kreal), intent(in) :: ecoord(3, nn), u(3, nn*2), du(3, nn*2), thick
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: qf(:)
    real(kind=kreal), intent(in), optional :: nddisp(3, nn)

    integer :: i, sstable(ndof*nn)
    real(kind=kreal) :: mixed_disp(ndof*nn), natural_disp(ndof*nn)
    real(kind=kreal) :: shell_u(6, nn), shell_du(6, nn), shell_nddisp(6, nn)
    logical :: is_mixed

    mixed_disp = 0.0D0
    do i = 1, nn
      mixed_disp(ndof*(i-1)+1:ndof*(i-1)+3) = u(1:3, 2*i-1)+du(1:3, 2*i-1)
      mixed_disp(ndof*(i-1)+4:ndof*(i-1)+6) = u(1:3, 2*i)+du(1:3, 2*i)
    end do

    call ShellMixedDofMap(mixflag, nn, ndof, sstable, is_mixed)
    natural_disp = mixed_disp
    if( is_mixed ) then
      do i = 1, ndof*nn
        natural_disp(sstable(i)) = mixed_disp(i)
      end do
    endif

    shell_u = 0.0D0
    do i = 1, nn
      shell_du(:, i) = natural_disp(ndof*(i-1)+1:ndof*i)
    end do
    if( present(nddisp) ) then
      shell_nddisp = shell_du
      shell_nddisp(1:3, :) = nddisp
      call UPDATE_Shell_MITC(etype, nn, ndof, ecoord, shell_u, shell_du, gausses, qf, thick, mixflag, nddisp=shell_nddisp)
    else
      call UPDATE_Shell_MITC(etype, nn, ndof, ecoord, shell_u, shell_du, gausses, qf, thick, mixflag)
    endif
  end subroutine UPDATE_Shell_MITC33
  

  !--------------------------------------------------------------------
  !> Calculate the consistent and lumped mass matrices of a MITC shell.
  subroutine mass_Shell(etype, nn, elem, rho, thick, gausses, mass, lumped)
    use mMechGauss
    implicit none

    integer(kind=kint), intent(in) :: etype, nn
    real(kind=kreal), intent(in) :: elem(3, nn), rho, thick
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: mass(:, :), lumped(:)

    integer(kind=kint), parameter :: ndof = 6_kint
    integer :: ny, nsize, lx, ly, nb, i, ilayer, nlayer, idof
    integer(kind=kint) :: ierr_quad
    real(kind=kreal) :: v1(3, nn), v2(3, nn), v3(3, nn), director(3, nn)
    real(kind=kreal) :: naturalcoord(2), shapefunc(nn), shapederiv(nn, 2)
    real(kind=kreal) :: covariant_basis(3, 3), N(3, ndof*nn)
    real(kind=kreal) :: zeta, surface_weight, thickness_weight, integration_weight, layer_weight
    real(kind=kreal) :: totalmass, totdiag

    nsize = ndof*nn
    mass(1:nsize, 1:nsize) = 0.0D0
    lumped(1:nsize) = 0.0D0
    totalmass = 0.0D0
    ny = NumOfShellThicknessQuadPoints(etype)
    call ShellMITC_SetupReferenceDirectors(etype, nn, thick, elem, v1, v2, v3, director)

    nlayer = gausses(1)%pMaterial%totallyr
    do ilayer = 1, nlayer
      layer_weight = gausses(1)%pMaterial%shell_var(ilayer)%weight
      do ly = 1, ny
        call fstr_shell_layer_quadrature_gauss(etype, gausses(1), ilayer, ly, zeta, thickness_weight, ierr_quad)
        if (ierr_quad /= 0) cycle
        do lx = 1, NumOfQuadPoints(etype)
          call getQuadPoint(etype, lx, naturalcoord)
          surface_weight = getWeight(etype, lx)
          call getShapeFunc(etype, naturalcoord, shapefunc)
          call getShapeDeriv(etype, naturalcoord, shapederiv)
          call ShellMITC_CovariantBasis(nn, elem, director, zeta, shapefunc, shapederiv, covariant_basis)
          integration_weight = surface_weight*thickness_weight *ShellMITC_CovariantJacobian(covariant_basis)*layer_weight
          call ShellMITC_BuildInterpolationMatrix(nn, ndof, zeta, shapefunc, director, N)
          mass(1:nsize, 1:nsize) = mass(1:nsize, 1:nsize) + rho*integration_weight*matmul(transpose(N), N)
          totalmass = totalmass+rho*integration_weight
        end do
      end do
    end do

    totalmass = 3.0D0*totalmass
    totdiag = 0.0D0
    do nb = 1, nn
      do i = 1, 3
        idof = ndof*(nb-1)+i
        totdiag = totdiag+mass(idof, idof)
      end do
    end do
    do nb = 1, nn
      do i = 1, ndof
        idof = ndof*(nb-1)+i
        lumped(idof) = mass(idof, idof)*totalmass/totdiag
      end do
    end do
  end subroutine mass_Shell

end module m_static_LIB_shell
