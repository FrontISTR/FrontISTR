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

<<<<<<< HEAD
  subroutine fstr_updatedof_dummy( ndof, hecMESH, dummy, elements, vec_old, vec_new )
    integer(kind=kint), intent(in)            :: ndof
    type(hecmwST_local_mesh), intent(in)      :: hecMESH      !< mesh information
    type(tDummy), intent(in)                  :: dummy        !< dummy info
    type(tElement), pointer, intent(inout)    :: elements(:)  !< elements info(dummy flags will be updated)
    real(kind=kreal), pointer, intent(in)     :: vec_old(:)
    real(kind=kreal), pointer, intent(inout)  :: vec_new(:)

    integer(kind=kint) :: icel, in0
    integer(kind=kint) :: iS, iE, ic_type, nodlocal, i
    real(kind=kreal), pointer   :: active(:)

    allocate(active(hecMESH%n_node))
    active = -1.d0

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)

      if (hecmw_is_etype_link(ic_type)) cycle
      if(ic_type == 3414) cycle

      do icel = iS, iE
        if( elements(icel)%dummy_flag > 0 ) cycle
        in0 = hecMESH%elem_node_index(icel-1)
        nn = hecmw_get_max_node(ic_type)
        do i = 1, nn
          nodlocal = hecMESH%elem_node_item(in0+i)
          active(nodlocal) = 1.d0
        enddo
      enddo
    enddo

    call hecmw_update_1_R(hecMESH,active,hecMESH%n_node)

    do i = 1, hecMESH%n_node
      if( active(i) > 0.d0 ) cycle
      vec_new(ndof*(i-1)+1:ndof*i) = vec_old(ndof*(i-1)+1:ndof*i)
    end do

    deallocate(active)

  end subroutine

  subroutine output_dummy_flag( hecMESH, elements, outval )
    type(hecmwST_local_mesh), intent(in)     :: hecMESH      !< mesh information
    type(tElement), pointer, intent(in)      :: elements(:)  !< elements info(dummy flags will be updated)
    real(kind=kreal), pointer, intent(inout) :: outval(:)    !< active dummy flag

    integer(kind=kint) :: icel

    outval = 0.d0
    do icel = 1, hecMESH%n_elem
      if( elements(icel)%dummy_flag > 0 ) outval(icel) = 1.d0
    end do

  end subroutine

=======
>>>>>>> 6b55241c ([dummy element] Added dummy data and reading !DUMMY keyword)
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
<<<<<<< HEAD
        elements(icel)%dummy_flag = kDUM_INACTIVE
=======
        elements(ik)%dummy_flag = kDUM_INACTIVE
>>>>>>> 6b55241c ([dummy element] Added dummy data and reading !DUMMY keyword)
      end do
    end do

  end subroutine

end module m_fstr_dummy
