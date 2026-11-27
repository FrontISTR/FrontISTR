!-------------------------------------------------------------------------------
! Copyright (c) 2021 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module manage the parameters for contact calculation
!!
module mContactParam

  use hecmw

  implicit none

  private
  public :: tContactParam
  public :: init_ContactParam
  public :: tContactInterference
  public :: init_Contact_IF
  public :: tInteractionProperty
  public :: init_InteractionProperty

  type tContactParam
    character(HECMW_NAME_LEN) :: name    !< name of contact parameter set
    real(kind=kreal) :: CLEARANCE        !< ordinary clearance
    real(kind=kreal) :: CLR_SAME_ELEM    !< clearance for already-in-contct elems (loosen to avoid moving too easily)
    real(kind=kreal) :: CLR_DIFFLPOS     !< clearance to be recognized as different position (loosen to avoid oscillation)
    real(kind=kreal) :: CLR_CAL_NORM     !< clearance used when calculating surface normal
    real(kind=kreal) :: DISTCLR_INIT     !< dist clearance for initial scan
    real(kind=kreal) :: DISTCLR_FREE     !< dist clearance for free nodes (wait until little penetration to be judged as contact)
    real(kind=kreal) :: DISTCLR_NOCHECK  !< dist clearance for skipping distance check for nodes already in contact
                                         !< (big value to keep contact because contact-to-free is judged by tensile force)
    real(kind=kreal) :: DISTCLR_C2F      !< dist clearance for the contact-to-free judgment of MORTAR=YES
                                         !< (scaled by master reflen; forms hysteresis band with DISTCLR_INIT/FREE on OFF->ON side)
    real(kind=kreal) :: PENCLR_FREE      !< penetration clearance for free nodes (deeper nodes are not judged as contact)
    real(kind=kreal) :: PENCLR_NOCHECK   !< penetration clearance for nodes already in contact
                                         !< (big value to keep contact because contact-to-free is not judged by depth)
    real(kind=kreal) :: TENSILE_FORCE    !< tensile force to be judged as free node
    real(kind=kreal) :: NEAR_DIST        !< distance threshold for CONTACTNEAR detection [length] (<=0: disabled)
  end type tContactParam

  type tContactInterference
    integer(kind=kint) :: if_type
    real(kind=kreal)   :: etime
    real(kind=kreal)   :: initial_pos
    real(kind=kreal)   :: end_pos
    character(HECMW_NAME_LEN) :: cp_name
  end type tContactInterference

  !> Interaction property: the friction and the penalties that contact pairs refer to by name
  type tInteractionProperty
    character(HECMW_NAME_LEN) :: name    !< name of interaction property
    real(kind=kreal) :: fcoeff           !< friction coefficient
    real(kind=kreal) :: nPenalty         !< normal penalty
    real(kind=kreal) :: tPenalty         !< tangential penalty
  end type tInteractionProperty

contains

  subroutine init_ContactParam( cparam )
    type(tContactParam), intent(out) :: cparam
    cparam%name            = ''
    cparam%CLEARANCE       = 1.d-2
    cparam%CLR_SAME_ELEM   = 5.d-3
    cparam%CLR_DIFFLPOS    = 1.d-2
    cparam%CLR_CAL_NORM    = 1.d-4
    cparam%DISTCLR_INIT    = 1.d-2
    cparam%DISTCLR_FREE    =-1.d-6
    cparam%DISTCLR_NOCHECK = 1.d0
    cparam%DISTCLR_C2F     = 1.d-3
    cparam%PENCLR_FREE     = 1.d0
    cparam%PENCLR_NOCHECK  = 1.d10
    cparam%TENSILE_FORCE   =-1.d-2
    cparam%NEAR_DIST       = 0.0d0
  end subroutine init_ContactParam

  subroutine init_Contact_IF( contact_if )
    type(tContactInterference), intent(out) :: contact_if
    contact_if%if_type     = 1
    contact_if%etime       = 1.0d0
    contact_if%initial_pos = 0.d0
    contact_if%end_pos     = 0.d0
    contact_if%cp_name = ''
  end subroutine init_Contact_IF

  !> The penalties default to those of !CONTACT without NPENALTY, TPENALTY and factor
  subroutine init_InteractionProperty( iprop )
    type(tInteractionProperty), intent(out) :: iprop
    iprop%name         = ''
    iprop%fcoeff       = 0.d0
    iprop%nPenalty     = 5.d0
    iprop%tPenalty     = 0.5d0
  end subroutine init_InteractionProperty

end module mContactParam
