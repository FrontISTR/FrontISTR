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

contains

  subroutine init_ContactParam( cparam )
    type(tContactParam), intent(out) :: cparam
    cparam%name            = ''
    cparam%CLEARANCE       = 1.d-4
    cparam%CLR_SAME_ELEM   = 5.d-3
    cparam%CLR_DIFFLPOS    = 1.d-2
    cparam%CLR_CAL_NORM    = 1.d-4
    cparam%DISTCLR_INIT    = 1.d-6
    cparam%DISTCLR_FREE    =-1.d-6
    cparam%DISTCLR_NOCHECK = 1.d0
    cparam%TENSILE_FORCE   =-1.d-8
    cparam%BOX_EXP_RATE    = 1.05d0
  end subroutine init_ContactParam

end module mContactParam
