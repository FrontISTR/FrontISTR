!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact damping module for CONTACTNEAR state
!!
!! Provides smoothstep weight function w(g), and routines to compute
!! damping nodal force and tangent stiffness for NEAR-state slave nodes.
!!
!! Damping residual:   r_damp = alpha * w(g) * (Bn^T * ddu) * Bn
!! Damping stiffness:  K_damp = alpha * w(g) * Bn (x) Bn
!!
!! where g = cstate%distance (frozen during Newton iteration),
!! alpha = DAMP_ALPHA * refStiff, and w(g) is a smoothstep
!! that transitions from 1 at g=0 to 0 at g=gact.
module m_fstr_contact_damping
  use hecmw
  use elementInfo
  use mContactDef
  use m_fstr_contact_elem_common
  implicit none

  private
  public :: is_damping_enabled
  public :: contact_damping_weight
  public :: getDampingStiffness
  public :: getDampingNodalForce

contains

  !> Check if damping is enabled for a contact pair
  pure logical function is_damping_enabled(contact)
    type(tContact), intent(in) :: contact
    is_damping_enabled = (contact%damp_alpha > 0.0d0 .and. contact%damp_gact > 0.0d0)
  end function is_damping_enabled

  !> Compute smoothstep weight w(g) for contact damping
  !!
  !! w(g) = 1                              if g <= 0       (penetrating)
  !! w(g) = 1 - 3t^2 + 2t^3               if 0 < g < gact (near zone)
  !! w(g) = 0                              if g >= gact    (far away)
  !!
  !! where t = g / gact
  pure real(kind=kreal) function contact_damping_weight(g, gact)
    real(kind=kreal), intent(in) :: g     !< signed distance (+ = gap, - = penetration)
    real(kind=kreal), intent(in) :: gact  !< activation distance (NEAR_DIST)
    real(kind=kreal) :: t

    if (g <= 0.0d0) then
      contact_damping_weight = 1.0d0
    else if (g >= gact) then
      contact_damping_weight = 0.0d0
    else
      t = g / gact
      contact_damping_weight = 1.0d0 - 3.0d0*t*t + 2.0d0*t*t*t
    end if
  end function contact_damping_weight

  !> Compute damping stiffness matrix K_damp = alpha * w(g) * Bn (x) Bn
  !!
  !! Same interface pattern as getContactStiffness_Alag:
  !!   compute stiffness matrix, caller assembles via hecmw_mat_ass_elem.
  subroutine getDampingStiffness(cstate, surf, alpha, gact, stiff, smoothing_type)
    type(tContactState), intent(in)    :: cstate        !< contact state (NEAR)
    type(tSurfElement), intent(in)     :: surf          !< master surface element
    real(kind=kreal), intent(in)       :: alpha         !< damping coefficient (DAMP_ALPHA * refStiff)
    real(kind=kreal), intent(in)       :: gact          !< activation distance (NEAR_DIST)
    real(kind=kreal), intent(out)      :: stiff(:,:)    !< damping stiffness matrix
    integer(kind=kint), optional, intent(in) :: smoothing_type

    integer(kind=kint) :: nnode, ndof, i, j
    real(kind=kreal) :: Bn(l_max_surface_node*3+3)
    real(kind=kreal) :: Tm(3, 3*(l_max_surface_node+1))
    real(kind=kreal) :: Tt(3, 3*(l_max_surface_node+1))
    real(kind=kreal) :: w, coeff

    nnode = size(surf%nodes)
    ndof = nnode*3 + 3

    stiff = 0.0d0

    ! Compute Bn directly via computeTm_Tt (no need for metric/Ht/Gt)
    call computeTm_Tt(cstate, surf, 0.0d0, Tm, Tt, smoothing_type, Bn=Bn)

    ! Compute weight from frozen distance
    w = contact_damping_weight(cstate%distance, gact)
    coeff = alpha * w

    if (coeff <= 0.0d0) return

    ! Build stiffness: K_damp = coeff * Bn (x) Bn
    do j = 1, ndof
      do i = 1, ndof
        stiff(i,j) = coeff * Bn(i) * Bn(j)
      enddo
    enddo

  end subroutine getDampingStiffness

  !> Compute damping nodal force vector f_damp = alpha * w(g) * (Bn^T * ddu) * Bn
  !!
  !! Same interface pattern as getContactNodalForce_Alag:
  !!   compute force vector, caller assembles via assemble_contact_force_residual.
  !!   ctNForce receives the damping force; ctTForce is zeroed.
  subroutine getDampingNodalForce(cstate, surf, ndDu, alpha, gact, ctNForce, ctTForce, smoothing_type)
    type(tContactState), intent(in)    :: cstate        !< contact state (NEAR)
    type(tSurfElement), intent(in)     :: surf          !< master surface element
    real(kind=kreal), intent(in)       :: ndDu(:)       !< nodal displacement increment (ddisp for this contact element)
    real(kind=kreal), intent(in)       :: alpha         !< damping coefficient (DAMP_ALPHA * refStiff)
    real(kind=kreal), intent(in)       :: gact          !< activation distance (NEAR_DIST)
    real(kind=kreal), intent(out)      :: ctNForce(:)   !< damping force vector (normal-like)
    real(kind=kreal), intent(out)      :: ctTForce(:)   !< tangential force vector (zeroed)
    integer(kind=kint), optional, intent(in) :: smoothing_type

    integer(kind=kint) :: nnode, ndof
    real(kind=kreal) :: Bn(l_max_surface_node*3+3)
    real(kind=kreal) :: Tm(3, 3*(l_max_surface_node+1))
    real(kind=kreal) :: Tt(3, 3*(l_max_surface_node+1))
    real(kind=kreal) :: w, coeff, delta_gn

    nnode = size(surf%nodes)
    ndof = nnode*3 + 3

    ctNForce = 0.0d0
    ctTForce = 0.0d0

    ! Compute Bn directly via computeTm_Tt (no need for metric/Ht/Gt/elemcrd)
    call computeTm_Tt(cstate, surf, 0.0d0, Tm, Tt, smoothing_type, Bn=Bn)

    ! Compute weight from frozen distance
    w = contact_damping_weight(cstate%distance, gact)
    coeff = alpha * w

    if (coeff <= 0.0d0) return

    ! Compute normal gap increment: delta_gn = Bn^T * ddu
    delta_gn = dot_product(Bn(1:ndof), ndDu(1:ndof))

    ! Damping force: f_damp = coeff * delta_gn * Bn
    ctNForce(1:ndof) = coeff * delta_gn * Bn(1:ndof)

    ! Lagrange row (not used, set to 0 for interface compatibility)
    ctNForce((nnode+1)*3+1) = 0.0d0
    ctTForce((nnode+1)*3+1) = 0.0d0

  end subroutine getDampingNodalForce

end module m_fstr_contact_damping
