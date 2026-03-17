!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide a function to elemact elements
module m_fstr_elemact
  use hecmw
  use m_fstr
  use m_elemact

contains

  !> Apply amplitude-based element activation control
  !> amp_val > 0.5: ACTIVE, amp_val <= 0.5: INACTIVE
  subroutine apply_amplitude_control( hecMESH, elemact, dumid, elements, amp_id, ctime )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH
    type(tElemact), intent(in)             :: elemact
    integer(kind=kint), intent(in)         :: dumid, amp_id
    type(tElement), pointer, intent(inout) :: elements(:)
    real(kind=kreal), intent(in)           :: ctime
    
    real(kind=kreal)   :: amp_val
    integer(kind=kint) :: amp_state
    
    call hecmw_get_amplitude_value(hecMESH%amp, amp_id, ctime, amp_val)
    
    ! Simple rule: amp_val > 0.5 -> ACTIVE, otherwise INACTIVE
    if( amp_val > 0.5d0 ) then
      amp_state = kELACT_ACTIVE
    else
      amp_state = kELACT_INACTIVE
    endif
    
    call set_elemact_flag( hecMESH, elemact, dumid, elements, amp_state, .false. )
  end subroutine apply_amplitude_control

  subroutine fstr_update_elemact_solid( hecMESH, fstrSOLID, cstep, ctime )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(fstr_solid), intent(inout)        :: fstrSOLID    !< fstr_solid
    integer(kind=kint), intent(in)         :: cstep        !< current step number
    real(kind=kreal), intent(in)           :: ctime        !< current analysis time

    integer(kind=kint) :: idum, amp_id, gid
    real(kind=kreal)   :: amp_val
    integer(kind=kint) :: target_state  ! Target state from control file (STATE=ON/OFF)

    do idum = 1, fstrSOLID%elemact%ELEMACT_egrp_tot
      gid = fstrSOLID%elemact%ELEMACT_egrp_GRPID(idum)
      if( .not. fstr_isElemActivationActive( fstrSOLID, gid, cstep ) ) then
        call set_elemact_flag( hecMESH, fstrSOLID%elemact, idum, fstrSOLID%elements, kELACT_ACTIVE, .false. )
        cycle
      endif

      amp_id = fstrSOLID%elemact%ELEMACT_egrp_amp(idum)
      if( amp_id > 0 ) then
        ! Amplitude control overrides STATE
        call apply_amplitude_control( hecMESH, fstrSOLID%elemact, idum, fstrSOLID%elements, amp_id, ctime )
        cycle
      end if

      ! No amplitude: use STATE from control file
      target_state = fstrSOLID%elemact%ELEMACT_egrp_state(idum)

      if( fstrSOLID%elemact%ELEMACT_egrp_depends(idum) == kELACTD_NONE ) then
        call set_elemact_flag( hecMESH, fstrSOLID%elemact, idum, fstrSOLID%elements, target_state, .false. )
      else
        call set_elemact_flag( hecMESH, fstrSOLID%elemact, idum, fstrSOLID%elements, target_state, .true. )
      endif
    end do

  end subroutine

  subroutine fstr_update_elemact_solid_by_value( hecMESH, fstrSOLID, cstep, ctime )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(fstr_solid), intent(inout)        :: fstrSOLID    !< fstr_solid
    integer(kind=kint), intent(in)         :: cstep        !< current step number
    real(kind=kreal), intent(in)           :: ctime        !< current analysis time

    integer(kind=kint) :: idum, amp_id, gid
    integer(kind=kint) :: n_changed_local, n_changed_total
    real(kind=kreal)   :: amp_val

    n_changed_total = 0

    do idum = 1, fstrSOLID%elemact%ELEMACT_egrp_tot
      gid = fstrSOLID%elemact%ELEMACT_egrp_GRPID(idum)
      if( .not. fstr_isElemActivationActive( fstrSOLID, gid, cstep ) ) cycle

      amp_id = fstrSOLID%elemact%ELEMACT_egrp_amp(idum)
      if( amp_id > 0 ) then
        ! Amplitude control overrides stress-based control
        call apply_amplitude_control( hecMESH, fstrSOLID%elemact, idum, fstrSOLID%elements, amp_id, ctime )
        cycle
      end if

      ! No amplitude: use stress-based control
      n_changed_local = 0
      call activate_elemact_flag_by_value( hecMESH, fstrSOLID%elemact, idum, fstrSOLID%elements, n_changed_local )
      n_changed_total = n_changed_total + n_changed_local
    end do

    fstrSOLID%elemact%ELEMACT_n_changed = n_changed_total

    if( hecMESH%my_rank == 0 .and. n_changed_total > 0 ) then
      write(*,'(a,i6,a)') ' *** ELEMACT: ', n_changed_total, ' element(s) changed state (ACTIVE -> INACTIVE)'
    endif

  end subroutine

  subroutine fstr_update_elemact_heat( hecMESH, elemact, ctime, elements )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tElemact), intent(in)               :: elemact        !< elemact info
    real(kind=kreal), intent(in)           :: ctime        !< current analysis time
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(elemact flags will be updated)

    integer(kind=kint) :: idum, amp_id
    real(kind=kreal)   :: amp_val

    call clear_elemact_flag_all( hecMESH, elemact, elements )

    do idum = 1, elemact%ELEMACT_egrp_tot
      amp_id = elemact%ELEMACT_egrp_amp(idum)
      amp_val = 1.d0
      if( amp_id > 0 ) then
        call hecmw_get_amplitude_value(hecMESH%amp, amp_id, ctime, amp_val)
        if( amp_val < 1.d0 ) cycle
      end if

      call activate_elemact_flag( hecMESH, elemact, idum, elements )
    end do

  end subroutine

  subroutine fstr_updatedof_elemact( ndof, hecMESH, elemact, elements, vec_old, vec_new )
    integer(kind=kint), intent(in)            :: ndof
    type(hecmwST_local_mesh), intent(in)      :: hecMESH      !< mesh information
    type(tElemact), intent(in)                  :: elemact        !< elemact info
    type(tElement), pointer, intent(inout)    :: elements(:)  !< elements info(elemact flags will be updated)
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
        if( elements(icel)%elemact_flag == kELACT_INACTIVE ) cycle
        in0 = hecMESH%elem_node_index(icel-1)
        nn = hecmw_get_max_node(ic_type)
        do i = 1, nn
          nodlocal = hecMESH%elem_node_item(in0+i)
          active(nodlocal) = 1.d0
        enddo
      enddo
    enddo

    call hecmw_update_R(hecMESH,active,hecMESH%n_node,1)

    do i = 1, hecMESH%n_node
      if( active(i) > 0.d0 ) cycle
      vec_new(ndof*(i-1)+1:ndof*i) = vec_old(ndof*(i-1)+1:ndof*i)
    end do

    deallocate(active)

  end subroutine

  subroutine output_elemact_flag( hecMESH, elements, outval )
    type(hecmwST_local_mesh), intent(in)     :: hecMESH      !< mesh information
    type(tElement), pointer, intent(in)      :: elements(:)  !< elements info(elemact flags will be updated)
    real(kind=kreal), pointer, intent(inout) :: outval(:)    !< active elemact flag

    integer(kind=kint) :: icel

    outval = 0.d0
    do icel = 1, hecMESH%n_elem
      if( elements(icel)%elemact_flag == kELACT_INACTIVE ) outval(icel) = 1.d0
    end do

  end subroutine

  subroutine activate_elemact_flag( hecMESH, elemact, dumid, elements )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tElemact), intent(in)               :: elemact        !< elemact info
    integer(kind=kint), intent(in)         :: dumid        !< elemact id
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(elemact flags will be updated)

    integer(kind=kint) :: ig, iS0, iE0, ik, icel

    if( dumid < 0 .or. dumid > elemact%ELEMACT_egrp_tot ) return

    ig = elemact%ELEMACT_egrp_ID(dumid)
    iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
    iE0 = hecMESH%elem_group%grp_index(ig  )

    do ik=iS0,iE0
      icel = hecMESH%elem_group%grp_item(ik)
      elements(icel)%elemact_flag = kELACT_INACTIVE
      elements(icel)%elemact_coeff = elemact%ELEMACT_egrp_eps(dumid)
    end do

  end subroutine

  subroutine activate_elemact_flag_by_value( hecMESH, elemact, dumid, elements, n_changed )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tElemact), intent(in)               :: elemact        !< elemact info
    integer(kind=kint), intent(in)         :: dumid        !< elemact id
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(elemact flags will be updated)
    integer(kind=kint), intent(out)        :: n_changed    !< number of elements that changed state

    integer(kind=kint) :: ig, iS0, iE0, ik, icel, dtype, ig0
    integer(kind=kint) :: old_flag
    real(kind=kreal)   :: thlow, thup, stress(6), mises, ps
    integer(kind=kint) :: target_state  ! Target state from control file (STATE=ON/OFF)

    n_changed = 0
    if( dumid < 0 .or. dumid > elemact%ELEMACT_egrp_tot ) return
    if( elemact%ELEMACT_egrp_depends(dumid) == kELACTD_NONE ) return 

    ig = elemact%ELEMACT_egrp_ID(dumid)
    iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
    iE0 = hecMESH%elem_group%grp_index(ig  )

    thlow = elemact%ELEMACT_egrp_ts_lower(dumid)
    thup = elemact%ELEMACT_egrp_ts_upper(dumid)

    ! Get the target state from control file
    ! STATE=ON  -> kELACT_ACTIVE   (0): Elements start active
    ! STATE=OFF -> kELACT_INACTIVE (1): Elements start inactive
    target_state = kELACT_INACTIVE  ! Default is INACTIVE
    if (associated(elemact%ELEMACT_egrp_state)) then
      target_state = elemact%ELEMACT_egrp_state(dumid)
    endif

    do ik=iS0,iE0
      icel = hecMESH%elem_group%grp_item(ik)
      
      ! Case 1: target_state is INACTIVE (STATE=OFF in control file)
      ! Set to INACTIVE regardless of stress (no stress dependency check needed)
      if( target_state == kELACT_INACTIVE ) then
        elements(icel)%elemact_flag = kELACT_INACTIVE
        elements(icel)%elemact_coeff = elemact%ELEMACT_egrp_eps(dumid)
        cycle
      endif
      
      ! Case 2: Element is already INACTIVE
      ! Once deactivated, element stays INACTIVE (no reactivation: one-way transition)
      if( elements(icel)%elemact_flag == kELACT_INACTIVE ) cycle
      
      ! Case 3: target_state is ACTIVE (STATE=ON) and element is currently ACTIVE
      ! Check stress to decide whether to deactivate: ACTIVE -> INACTIVE transition
      do ig0=1,size(elements(icel)%gausses)
        ! get mises
        if( elemact%ELEMACT_egrp_depends(dumid) == kELACTD_STRESS ) then
          stress(1:6) = elements(icel)%gausses(ig0)%stress(1:6)
        elseif( elemact%ELEMACT_egrp_depends(dumid) == kELACTD_STRAIN ) then
          stress(1:6) = elements(icel)%gausses(ig0)%strain(1:6)
        else
          return
        endif
        ps = ( stress(1) + stress(2) + stress(3) ) / 3.0d0
        mises = 0.5d0 * ( (stress(1)-ps)**2 + (stress(2)-ps)**2 + (stress(3)-ps)**2 )
        mises = mises + stress(4)**2 + stress(5)**2 + stress(6)**2
        mises = dsqrt( 3.0d0 * mises )
        
        ! If stress is OUT OF RANGE [thlow, thup], deactivate the element
        old_flag = elements(icel)%elemact_flag
        if( .not. (thlow <= mises .and. mises <= thup) ) then
          elements(icel)%elemact_flag = kELACT_INACTIVE
          elements(icel)%elemact_coeff = elemact%ELEMACT_egrp_eps(dumid)
        endif
        ! Count state change
        if( elements(icel)%elemact_flag /= old_flag ) n_changed = n_changed + 1
        
        exit
      enddo
    end do

  end subroutine

  subroutine set_elemact_flag( hecMESH, elemact, dumid, elements, flag, init_only )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tElemact), intent(in)               :: elemact        !< elemact info
    integer(kind=kint), intent(in)         :: dumid        !< elemact id
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(elemact flags will be updated)
    integer(kind=kint), intent(in)         :: flag        !< elemact id
    logical, intent(in)                    :: init_only   !< set flag if only initial elemact flag is set

    integer(kind=kint) :: ig, iS0, iE0, ik, icel

    if( dumid < 0 .or. dumid > elemact%ELEMACT_egrp_tot ) return

    ig = elemact%ELEMACT_egrp_ID(dumid)
    iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
    iE0 = hecMESH%elem_group%grp_index(ig  )

    do ik=iS0,iE0
      icel = hecMESH%elem_group%grp_item(ik)
      if( init_only ) then
        if( elements(icel)%elemact_flag == kELACT_UNDEFINED ) elements(icel)%elemact_flag = flag
      else
        elements(icel)%elemact_flag = flag
      endif
      elements(icel)%elemact_coeff = elemact%ELEMACT_egrp_eps(dumid)
    end do

  end subroutine

  subroutine clear_elemact_flag_all( hecMESH, elemact, elements )
    type(hecmwST_local_mesh), intent(in)   :: hecMESH      !< mesh information
    type(tElemact), intent(in)               :: elemact        !< elemact info
    type(tElement), pointer, intent(inout) :: elements(:)  !< elements info(elemact flags will be updated)

    integer(kind=kint) :: idum, ig, iS0, iE0, ik, icel

    do idum = 1, elemact%ELEMACT_egrp_tot
      ig = elemact%ELEMACT_egrp_GRPID(idum)
      iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
      iE0 = hecMESH%elem_group%grp_index(ig  )

      do ik=iS0,iE0
        icel = hecMESH%elem_group%grp_item(ik)
        elements(icel)%elemact_flag = kELACT_UNDEFINED
      end do
    end do

  end subroutine

end module m_fstr_elemact
