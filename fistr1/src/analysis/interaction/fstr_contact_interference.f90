!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provides interference fit (shrink) functions for contact
module m_fstr_contact_interference
  use hecmw
  use mContactDef
  implicit none

  public :: set_shrink_factor

contains

  subroutine set_shrink_factor(ctime, cstate, etime, if_type)
    real(kind=kreal),intent(in)         :: ctime, etime
    type(tContactState), intent(inout)  :: cstate        !< Recorde of contact information
    integer, intent(in)                 :: if_type

    if (if_type == C_IF_SLAVE .and. cstate%init_pos == 0.d0) then
      cstate%shrink_factor = 0.0d0; return
    end if
    cstate%shrink_factor = cstate%time_factor*ctime + cstate%init_pos
    if(ctime >= etime) cstate%shrink_factor = cstate%end_pos

  end subroutine set_shrink_factor

  subroutine get_shrink_elemact_surf(cstate, coords, nnode)
    use m_fstr_TimeInc
    integer(kind=kint)              :: nnode, i
    type(tContactState)             :: cstate
    real(kind=kreal)                :: coords(:)
    real(kind=kreal)   :: d(3)

    d = - cstate%shrink_factor * cstate%direction(1:3)

    do i = 1, nnode
      coords(1+i*3:(i+1)*3) = coords(1+i*3:(i+1)*3) + d(1:3)
    enddo

  end subroutine get_shrink_elemact_surf

end module m_fstr_contact_interference
