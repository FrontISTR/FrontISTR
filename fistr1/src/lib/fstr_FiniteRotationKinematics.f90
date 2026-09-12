!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Shared finite-rotation nodal kinematics and rotation algebra.

module m_fstr_FiniteRotationKinematics
  use hecmw
  use elementInfo, only: fe_mitc4_shell
  use m_fstr, only: fstr_solid
  use mMaterial, only: tMaterial, TOTALLAG, UPDATELAG, isElastic
  implicit none

  private

  public :: fstr_uses_finite_rotation_kinematics
  public :: fstr_has_finite_rotation_kinematics
  public :: fstr_mark_finite_rotation_nodes
  public :: fstr_is_finite_rotation_shell_element
  public :: ShellRotationVectorToMatrix
  public :: ShellRotationMatrixToVector
  public :: ShellSkewMatrix
  public :: ShellComposeRotationVector
  public :: ShellOrthonormalizeTriad
  public :: ShellUpdateTriadWithIncrement
  public :: ShellComposeNodalDisplacement

contains

  logical function fstr_is_finite_rotation_shell_element( etype, nn )
    integer(kind=kint), intent(in) :: etype
    integer(kind=kint), intent(in) :: nn

    fstr_is_finite_rotation_shell_element = etype == fe_mitc4_shell .and. nn == 4
  end function fstr_is_finite_rotation_shell_element

  logical function fstr_uses_finite_rotation_kinematics( etype, nn, material )
    integer(kind=kint), intent(in) :: etype
    integer(kind=kint), intent(in) :: nn
    type(tMaterial), intent(in)    :: material

    fstr_uses_finite_rotation_kinematics = ( fstr_is_finite_rotation_shell_element( etype, nn ) &
      .and. ( material%nlgeom_flag == TOTALLAG .or. material%nlgeom_flag == UPDATELAG ) &
      .and. isElastic( material%mtype ) )
  end function fstr_uses_finite_rotation_kinematics

  logical function fstr_has_finite_rotation_kinematics( hecMESH, fstrSOLID )
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(in)         :: fstrSOLID

    integer(kind=kint) :: itype, is, iE, ic_type, icel, iiS, nn

    fstr_has_finite_rotation_kinematics = .false.
    if( hecMESH%n_dof < 6 ) return
    if( .not. associated( fstrSOLID%elements ) ) return

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype)
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type /= fe_mitc4_shell ) cycle

      do icel = is, iE
        iiS = hecMESH%elem_node_index(icel-1)
        nn = hecMESH%elem_node_index(icel) - iiS
        if( .not. associated( fstrSOLID%elements(icel)%gausses ) ) cycle
        if( fstr_uses_finite_rotation_kinematics( ic_type, nn, &
            fstrSOLID%elements(icel)%gausses(1)%pMaterial ) ) then
          fstr_has_finite_rotation_kinematics = .true.
          return
        endif
      end do
    end do
  end function fstr_has_finite_rotation_kinematics

  subroutine fstr_mark_finite_rotation_nodes( hecMESH, fstrSOLID, ndof, shell_node_mode )
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(in)         :: fstrSOLID
    integer(kind=kint), intent(in)        :: ndof
    integer(kind=kint), intent(out)       :: shell_node_mode(:)

    integer(kind=kint) :: itype, is, iE, ic_type, icel, iiS, nn, j, node_id

    shell_node_mode(:) = 0
    if( ndof < 6 ) return
    if( .not. associated( fstrSOLID%elements ) ) return

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype)
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type /= fe_mitc4_shell ) cycle

      do icel = is, iE
        iiS = hecMESH%elem_node_index(icel-1)
        nn = hecMESH%elem_node_index(icel) - iiS
        if( .not. associated( fstrSOLID%elements(icel)%gausses ) ) cycle
        if( .not. fstr_uses_finite_rotation_kinematics( ic_type, nn, &
            fstrSOLID%elements(icel)%gausses(1)%pMaterial ) ) cycle

        do j = 1, nn
          node_id = hecMESH%elem_node_item(iiS+j)
          if( node_id > 0 .and. node_id <= size(shell_node_mode) ) shell_node_mode(node_id) = 1
        end do
      end do
    end do
  end subroutine fstr_mark_finite_rotation_nodes

  pure subroutine ShellRotationVectorToMatrix(theta, rotmat)
    real(kind=kreal), intent(in) :: theta(3)
    real(kind=kreal), intent(out) :: rotmat(3, 3)

    integer :: i
    real(kind=kreal) :: theta_norm, theta_norm2
    real(kind=kreal) :: sin_over_theta, one_minus_cos_over_theta2
    real(kind=kreal) :: skew(3, 3), skew2(3, 3)

    theta_norm2 = dot_product(theta, theta)
    theta_norm = dsqrt(theta_norm2)
    skew = ShellSkewMatrix(theta)
    skew2 = matmul(skew, skew)

    if (theta_norm < 1.0D-12) then
      sin_over_theta = 1.0D0-theta_norm2/6.0D0+theta_norm2*theta_norm2/120.0D0
      one_minus_cos_over_theta2 = 0.5D0-theta_norm2/24.0D0+theta_norm2*theta_norm2/720.0D0
    else
      sin_over_theta = dsin(theta_norm)/theta_norm
      one_minus_cos_over_theta2 = (1.0D0-dcos(theta_norm))/theta_norm2
    endif

    rotmat = sin_over_theta*skew+one_minus_cos_over_theta2*skew2
    do i = 1, 3
      rotmat(i, i) = rotmat(i, i)+1.0D0
    end do
  end subroutine ShellRotationVectorToMatrix

  pure subroutine ShellRotationMatrixToVector(rotmat, theta)
    real(kind=kreal), intent(in) :: rotmat(3, 3)
    real(kind=kreal), intent(out) :: theta(3)

    real(kind=kreal) :: pi, trace_r, cos_angle, angle, sin_angle
    real(kind=kreal) :: axis(3), axis_abs

    pi = 4.0D0*datan(1.0D0)
    trace_r = rotmat(1, 1)+rotmat(2, 2)+rotmat(3, 3)
    cos_angle = max(-1.0D0, min(1.0D0, 0.5D0*(trace_r-1.0D0)))
    angle = dacos(cos_angle)

    theta(1) = rotmat(3, 2)-rotmat(2, 3)
    theta(2) = rotmat(1, 3)-rotmat(3, 1)
    theta(3) = rotmat(2, 1)-rotmat(1, 2)

    if (angle < 1.0D-12) then
      theta = 0.5D0*theta
    else if (pi-angle < 1.0D-8) then
      axis(1) = dsqrt(max(0.0D0, 0.5D0*(rotmat(1, 1)+1.0D0)))
      axis(2) = dsqrt(max(0.0D0, 0.5D0*(rotmat(2, 2)+1.0D0)))
      axis(3) = dsqrt(max(0.0D0, 0.5D0*(rotmat(3, 3)+1.0D0)))
      if (rotmat(2, 1)+rotmat(1, 2) < 0.0D0) axis(2) = -axis(2)
      if (rotmat(3, 1)+rotmat(1, 3) < 0.0D0) axis(3) = -axis(3)
      axis_abs = dsqrt(dot_product(axis, axis))
      if (axis_abs > 1.0D-12) then
        theta = angle*axis/axis_abs
      else
        theta = 0.0D0
      endif
    else
      sin_angle = dsin(angle)
      theta = angle*theta/(2.0D0*sin_angle)
    endif
  end subroutine ShellRotationMatrixToVector

  pure function ShellSkewMatrix(vector) result(matrix)
    real(kind=kreal), intent(in) :: vector(3)
    real(kind=kreal) :: matrix(3, 3)

    matrix = 0.0D0
    matrix(1, 2) = -vector(3)
    matrix(1, 3) = vector(2)
    matrix(2, 1) = vector(3)
    matrix(2, 3) = -vector(1)
    matrix(3, 1) = -vector(2)
    matrix(3, 2) = vector(1)
  end function ShellSkewMatrix

  pure subroutine ShellComposeRotationVector(theta_old, theta_inc, theta_new)
    real(kind=kreal), intent(in) :: theta_old(3), theta_inc(3)
    real(kind=kreal), intent(out) :: theta_new(3)

    real(kind=kreal) :: rot_old(3, 3), rot_inc(3, 3), rot_new(3, 3)

    call ShellRotationVectorToMatrix(theta_old, rot_old)
    call ShellRotationVectorToMatrix(theta_inc, rot_inc)
    rot_new = matmul(rot_inc, rot_old)
    call ShellRotationMatrixToVector(rot_new, theta_new)
  end subroutine ShellComposeRotationVector


  pure subroutine ShellOrthonormalizeTriad(triad_in, triad_out)
    real(kind=kreal), intent(in) :: triad_in(3, 3)
    real(kind=kreal), intent(out) :: triad_out(3, 3)

    real(kind=kreal) :: e1(3), e2(3), e3(3), normv

    e1 = triad_in(:, 1)
    e3 = triad_in(:, 3)
    normv = dsqrt(dot_product(e3, e3))
    if (normv < 1.0D-14) then
      e3 = (/ 0.0D0, 0.0D0, 1.0D0 /)
    else
      e3 = e3/normv
    endif

    e1 = e1-dot_product(e1, e3)*e3
    normv = dsqrt(dot_product(e1, e1))
    if (normv < 1.0D-14) then
      if (dabs(e3(1)) < 0.9D0) then
        e1 = (/ 1.0D0, 0.0D0, 0.0D0 /)
      else
        e1 = (/ 0.0D0, 1.0D0, 0.0D0 /)
      endif
      e1 = e1-dot_product(e1, e3)*e3
      normv = dsqrt(dot_product(e1, e1))
    endif
    e1 = e1/normv
    e2(1) = e3(2)*e1(3)-e3(3)*e1(2)
    e2(2) = e3(3)*e1(1)-e3(1)*e1(3)
    e2(3) = e3(1)*e1(2)-e3(2)*e1(1)

    triad_out(:, 1) = e1
    triad_out(:, 2) = e2
    triad_out(:, 3) = e3
  end subroutine ShellOrthonormalizeTriad

  pure subroutine ShellUpdateTriadWithIncrement(triad_old, drill_old, theta_inc, triad_new, drill_new)
    real(kind=kreal), intent(in) :: triad_old(3, 3), drill_old, theta_inc(3)
    real(kind=kreal), intent(out) :: triad_new(3, 3), drill_new

    integer :: isub, nsub
    real(kind=kreal) :: triad_base(3, 3)
    real(kind=kreal) :: director(3), theta_step(3), theta_phys(3)
    real(kind=kreal) :: drill_acc, drill_inc, theta_norm
    real(kind=kreal) :: rot_inc(3, 3), rot_full(3, 3), triad_trial(3, 3)

    call ShellOrthonormalizeTriad(triad_old, triad_base)
    theta_norm = dsqrt(dot_product(theta_inc, theta_inc))
    nsub = max(1, ceiling(theta_norm/5.0D-2))
    theta_step = theta_inc/dble(nsub)
    drill_acc = drill_old

    do isub = 1, nsub
      director = triad_base(:, 3)
      drill_inc = dot_product(theta_step, director)
      theta_phys = theta_step-drill_inc*director
      call ShellRotationVectorToMatrix(theta_phys, rot_inc)
      triad_trial = matmul(rot_inc, triad_base)
      call ShellOrthonormalizeTriad(triad_trial, triad_base)
      drill_acc = drill_acc+drill_inc
    end do

    call ShellRotationVectorToMatrix(theta_inc, rot_full)
    triad_trial = triad_base
    triad_trial(:, 3) = matmul(rot_full, triad_old(:, 3))
    call ShellOrthonormalizeTriad(triad_trial, triad_new)
    drill_new = drill_acc
  end subroutine ShellUpdateTriadWithIncrement

  pure subroutine ShellComposeNodalDisplacement(ndof, nn, disp_old, disp_inc, disp_new)
    integer(kind=kint), intent(in) :: ndof, nn
    real(kind=kreal), intent(in) :: disp_old(:, :), disp_inc(:, :)
    real(kind=kreal), intent(out) :: disp_new(6, nn)

    integer :: i, ndof_copy

    disp_new = 0.0D0
    ndof_copy = min(ndof, 6)
    do i = 1, nn
      disp_new(1:min(3, ndof_copy), i) = disp_old(1:min(3, ndof_copy), i)+disp_inc(1:min(3, ndof_copy), i)
      if (ndof_copy >= 6) then
        call ShellComposeRotationVector(disp_old(4:6, i), disp_inc(4:6, i), disp_new(4:6, i))
      else if (ndof_copy > 3) then
        disp_new(4:ndof_copy, i) = disp_old(4:ndof_copy, i)+disp_inc(4:ndof_copy, i)
      endif
    end do
  end subroutine ShellComposeNodalDisplacement

end module m_fstr_FiniteRotationKinematics
