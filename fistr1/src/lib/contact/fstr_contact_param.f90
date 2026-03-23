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
    real(kind=kreal) :: TENSILE_FORCE    !< tensile force to be judged as free node
    real(kind=kreal) :: BOX_EXP_RATE     !< recommended: (1.0..2.0] (the smaller the faster, the bigger the safer)
  end type tContactParam

  type tContactInterference
    integer(kind=kint) :: if_type
    real(kind=kreal)   :: etime
    real(kind=kreal)   :: initial_pos
    real(kind=kreal)   :: end_pos
    character(HECMW_NAME_LEN) :: cp_name
  end type tContactInterference

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
    cparam%TENSILE_FORCE   =-1.d-2
    cparam%BOX_EXP_RATE    = 1.05d0
  end subroutine init_ContactParam

  subroutine init_Contact_IF( contact_if )
    type(tContactInterference), intent(out) :: contact_if
    contact_if%if_type     = 1
    contact_if%etime       = 1.0d0
    contact_if%initial_pos = 0.d0
    contact_if%end_pos     = 0.d0
    contact_if%cp_name = ''
  end subroutine init_Contact_IF

end module mContactParam
