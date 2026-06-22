!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Shared predicates for finite-rotation nodal kinematics.

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

contains

  logical function fstr_uses_finite_rotation_kinematics( etype, nn, material )
    integer(kind=kint), intent(in) :: etype
    integer(kind=kint), intent(in) :: nn
    type(tMaterial), intent(in)    :: material

    fstr_uses_finite_rotation_kinematics = ( etype == fe_mitc4_shell .and. nn == 4 &
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

end module m_fstr_FiniteRotationKinematics
