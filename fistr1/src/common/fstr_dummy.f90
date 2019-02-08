!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide a function to dummy elements
module m_fstr_dummy
  use hecmw
  use m_fstr
  use m_dummy

contains

  subroutine fstr_update_dummy_solid( hecMESH, fstrSOLID, cstep, ctime )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(fstr_solid), intent(inout)        :: fstrSOLID    !< fstr_solid
    integer(kind=kint), intent(in)         :: cstep        !< current step number
    real(kind=kreal), intent(in)           :: ctime        !< current analysis time

    integer(kind=kint) :: idum, amp_id, gid
    real(kind=kreal)   :: amp_val

    call clear_dummy_flag_all( hecMESH, fstrSOLID%dummy, fstrSOLID%elements )

    do idum = 1, fstrSOLID%dummy%DUMMY_egrp_tot
      gid = fstrSOLID%dummy%DUMMY_egrp_GRPID(idum)
      if( .not. fstr_isDummyActive( fstrSOLID, gid, cstep ) ) cycle

      amp_id = fstrSOLID%dummy%DUMMY_egrp_amp(idum)
      amp_val = 1.d0
      if( amp_id > 0 ) then
        call hecmw_get_amplitude_value(hecMESH%amp, amp_id, ctime, amp_val)
        if( amp_val < 1.d0 ) cycle
      end if

      call activate_dummy_flag( hecMESH, fstrSOLID%dummy, idum, fstrSOLID%elements )
    end do

  end subroutine

  subroutine fstr_update_dummy_heat( hecMESH, dummy, ctime, elements )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tDummy), intent(in)               :: dummy        !< dummy info
    real(kind=kreal), intent(in)           :: ctime        !< current analysis time
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(dummy flags will be updated)

    integer(kind=kint) :: idum, amp_id
    real(kind=kreal)   :: amp_val

    call clear_dummy_flag_all( hecMESH, dummy, elements )

    do idum = 1, dummy%DUMMY_egrp_tot
      amp_id = dummy%DUMMY_egrp_amp(idum)
      amp_val = 1.d0
      if( amp_id > 0 ) then
        call hecmw_get_amplitude_value(hecMESH%amp, amp_id, ctime, amp_val)
        if( amp_val < 1.d0 ) cycle
      end if

      call activate_dummy_flag( hecMESH, dummy, idum, elements )
    end do

  end subroutine

  subroutine activate_dummy_flag( hecMESH, dummy, dumid, elements )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tDummy), intent(in)               :: dummy        !< dummy info
    integer(kind=kint), intent(in)         :: dumid        !< dummy id
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(dummy flags will be updated)

    integer(kind=kint) :: ig, iS0, iE0, ik, icel

    if( dumid < 0 .or. dumid > dummy%DUMMY_egrp_tot ) return

    ig = dummy%DUMMY_egrp_ID(dumid)
    iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
    iE0 = hecMESH%elem_group%grp_index(ig  )

    do ik=iS0,iE0
      icel = hecMESH%elem_group%grp_item(ik)
      elements(icel)%dummy_flag = kDUM_ACTIVE
      elements(icel)%dummy_coeff = dummy%DUMMY_egrp_eps(dumid)
    end do

  end subroutine

  subroutine clear_dummy_flag_all( hecMESH, dummy, elements )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tDummy), intent(in)               :: dummy        !< dummy info
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(dummy flags will be updated)

    integer(kind=kint) :: idum, ig, iS0, iE0, ik, icel

    do idum = 1, dummy%DUMMY_egrp_tot
      ig = dummy%DUMMY_egrp_GRPID(idum)
      iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
      iE0 = hecMESH%elem_group%grp_index(ig  )

      do ik=iS0,iE0
        icel = hecMESH%elem_group%grp_item(ik)
        elements(ik)%dummy_flag = kDUM_INACTIVE
      end do
    end do

  end subroutine

end module m_fstr_dummy
