!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief  This module manages the data structure for contact calculation
!!
!!  Contact calculation takes into act after calling the following three
!!  subrotuines provided in this module
!!-#      Reading contact definition with subroutine: fstr_ctrl_get_CONTACT
!!-#      Check its consistency with mesh definition: fstr_contact_check
!!-#      Initializing the contact calculation      : fstr_contact_init
module mContactDef

  use hecmw
  use elementInfo
  use mSurfElement
  use m_hecmw_contact_comm
  use bucket_search
  use mContactParam

  implicit none

  real(kind=kreal), save :: cgn=1.d-5 !< convergent condition of penetration
  real(kind=kreal), save :: cgt=1.d-3 !< convergent condition of relative tangent disp

  real(kind=kreal), save :: gnt(2)    !< 1:current average penetration;
  !< 2:current relative tangent displacement
  real(kind=kreal), save :: bakgnt(2) !< 1:current average penetration;
  !< 2:current relative tangent displacement!

  integer, parameter :: CONTACTUNKNOWN = -1
  !> contact state definition
  integer, parameter :: CONTACTFREE = -1
  integer, parameter :: CONTACTSTICK = 1
  integer, parameter :: CONTACTSLIP = 2

  !> contact type or algorithm definition
  integer, parameter :: CONTACTTIED = 1
  integer, parameter :: CONTACTGLUED = 2
  integer, parameter :: CONTACTSSLID = 3
  integer, parameter :: CONTACTFSLID = 4

  !> contact interference type
  integer, parameter :: C_IF_SLAVE = 1
  integer, parameter :: C_IF_MASTER = 2

  !> This structure records contact status
  type tContactState
    integer          :: state !< -1:free, 1:in contact, or other needed
    integer          :: surface !< contacting surface number
    real(kind=kreal) :: distance !< penetration value
    real(kind=kreal) :: wkdist !< copy of penetration value
    real(kind=kreal) :: lpos(3) !< contact position(local coordinate)
    real(kind=kreal) :: gpos(3) !< contact position(global coordinate)
    real(kind=kreal) :: direction(3) !< contact direction
    real(kind=kreal) :: multiplier(3) !< Lagrangian multiplier or contact force
    !< 1: normal 2:tangent component
    real(kind=kreal) :: tangentForce(3) !< friction force
    real(kind=kreal) :: tangentForce1(3) !< friction force rotated by element(for trial friction force)
    real(kind=kreal) :: tangentForce_trial(3) !< trial friction force
    real(kind=kreal) :: tangentForce_final(3) !< final friction force
    real(kind=kreal)    :: reldisp(3)
    !
    real(kind=kreal)    :: shrink_factor
    real(kind=kreal)    :: time_factor
    real(kind=kreal)    :: init_pos
    real(kind=kreal)    :: end_pos
    integer             :: interference_flag
  end type

  !> Structure to includes all info needed by contact calculation
  type tContact
    ! following contact definition
    character(len=HECMW_NAME_LEN) :: name                    !< name
    integer                       :: ctype                   !< 1:node-surface 2: surface-surface
    integer                       :: group                   !< group number
    character(len=HECMW_NAME_LEN) :: pair_name               !< name of contact pair
    integer                       :: surf_id1, surf_id2      !< slave surface, master surface
    integer                       :: surf_id1_sgrp           !< surface group id of slave surface
    type(tSurfElement), pointer   :: master(:)=>null()       !< master surface (element )
    integer, pointer              :: slave(:)=>null()        !< slave surface (node)
    real(kind=kreal)              :: fcoeff                  !< coeeficient of friction
    real(kind=kreal)              :: nPenalty                !< normal penalty coefficient
    real(kind=kreal)              :: tPenalty                !< tangential penalty coefficient
    real(kind=kreal)              :: refStiff                !< reference stiffness for penalty calculation
    
    real(kind=kreal)    :: ctime
    integer(kind=kint)  :: if_type
    real(kind=kreal)    :: if_etime
    real(kind=kreal)    :: initial_pos
    real(kind=kreal)    :: end_pos
    ! following algorithm
    ! -1: not initialized
    ! 1: TIED-Just rigidly fixed the two surfaces
    ! 2: GLUED-Distance between the two surfaces to zero and glue them
    ! 3: SSLID-Small sliding contact( no position but with contact state change)
    ! 4: FSLID-Finite sliding contact (both changes in contact state and position possible)
    integer                       :: algtype                 !< algorithm flag

    logical                       :: mpced                   !< if turns into mpc condition
    logical                       :: symmetric               !< if symmetrizalized in friction calculation

    ! following contact state
    type(tContactState), pointer  :: states(:)=>null()       !< contact states of each slave nodes

    type(hecmwST_contact_comm)    :: comm                    !< contact communication table
    type(bucketDB)                :: master_bktDB            !< bucket DB for master surface

    type(tContactParam), pointer  :: cparam=>null()          !< contact parameter
  end type tContact

  type fstr_info_contactChange
    logical            :: active
    integer(kind=kint) :: contact2free           !< counter: contact to free state change
    integer(kind=kint) :: contact2neighbor       !< counter: contact to neighbor state change
    integer(kind=kint) :: contact2diffLpos       !< counter: contact to different local position state change
    integer(kind=kint) :: free2contact           !< counter: free to contact state change
    integer(kind=kint) :: contactNode_previous   !< previous number of nodes in contact
    integer(kind=kint) :: contactNode_current    !< current number of nodes in contact
  end type fstr_info_contactChange

  private :: is_MPC_available
  private :: is_active_contact

contains

  !> Initializer
  subroutine contact_state_init(cstate)
    type(tContactState), intent(inout) :: cstate !< contact state
    cstate%state = -1
    cstate%surface = -1
    cstate%distance = 0.0d0
    cstate%wkdist = 0.0d0
    cstate%lpos(:) = 0.0d0
    cstate%gpos(:) = 0.0d0
    cstate%direction(:) = 0.0d0
    cstate%multiplier(:) = 0.0d0
    cstate%tangentForce(:) = 0.0d0
    cstate%tangentForce1(:) = 0.0d0
    cstate%tangentForce_trial(:) = 0.0d0
    cstate%tangentForce_final(:) = 0.0d0
    cstate%reldisp(:) = 0.0d0
    cstate%shrink_factor = 0.0d0
    cstate%time_factor = 0.0d0
    cstate%init_pos = 0.0d0
    cstate%end_pos = 0.0d0
    cstate%interference_flag = 0
  end subroutine

  !> Copy
  subroutine contact_state_copy(cstate1, cstate2)
    type(tContactState), intent(in)    :: cstate1 !< contact state
    type(tContactState), intent(inout) :: cstate2 !< contact state
    cstate2 = cstate1
  end subroutine

  !> Print out contact state
  subroutine print_contact_state(fnum, cstate)
    integer, intent(in)             :: fnum !< file number
    type(tContactState), intent(in) :: cstate !< contact state
    write(fnum, *) "--Contact state=",cstate%state
    write(fnum, *) cstate%surface, cstate%distance
    write(fnum, *) cstate%lpos
    write(fnum, *) cstate%direction
    write(fnum, *) cstate%multiplier
  end subroutine

  !> If contact to mpc condiitions
  logical function is_MPC_available( contact )
    type(tContact), intent(in)        :: contact   !< contact definition
    is_MPC_available = .true.
    if( contact%fcoeff/=0.d0 ) is_MPC_available = .false.
  end function

  !> Write out contact definition
  subroutine fstr_write_contact( file, contact )
    integer(kind=kint), intent(in)    :: file      !< file number
    type(tContact), intent(in)        :: contact   !< contact definition
    integer :: i
    write(file,*) "CONTACT:", contact%ctype,contact%group,trim(contact%pair_name),contact%fcoeff
    write(file,*) "---Slave----"
    write(file,*) "num.slave",size(contact%slave)
    if( associated(contact%slave) ) then
      do i=1,size(contact%slave)
        write(file, *) contact%slave(i)
      enddo
    endif
    write(file,*) "----master---"
    write(file,*) "num.master",size(contact%master)
    if( associated(contact%master) ) then
      do i=1,size(contact%master)
        call write_surf( file, contact%master(i) )
      enddo
    endif
  end subroutine

  !> Finalizer
  subroutine fstr_contact_finalize( contact )
    type(tContact), intent(inout)     :: contact !< contact definition
    integer  :: i
    if( associated( contact%slave ) ) deallocate(contact%slave)
    if( associated( contact%master ) ) then
      do i=1,size( contact%master )
        call  finalize_surf( contact%master(i) )
      enddo
      deallocate(contact%master)
    endif
    if( associated(contact%states) ) deallocate(contact%states)
    call hecmw_contact_comm_finalize(contact%comm)
    call bucketDB_finalize( contact%master_bktDB )
  end subroutine

  !>  Check the consistency with given mesh of contact definition
  logical function fstr_contact_check( contact, hecMESH )
    type(tContact), intent(inout)     :: contact  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition

    integer  :: i
    logical  :: isfind

    fstr_contact_check = .false.

    ! if contact pair exist?
    isfind = .false.
    do i=1,hecMESH%contact_pair%n_pair
      if( hecMESH%contact_pair%name(i) == contact%pair_name ) then
        contact%ctype = hecMESH%contact_pair%type(i)
        contact%surf_id1 = hecMESH%contact_pair%slave_grp_id(i)
        contact%surf_id2 = hecMESH%contact_pair%master_grp_id(i)
        contact%surf_id1_sgrp = hecMESH%contact_pair%slave_orisgrp_id(i)
        isfind = .true.
      endif
    enddo
    if( .not. isfind ) return;
    if( contact%fcoeff<=0.d0 ) contact%fcoeff=0.d0
    if( contact%ctype < 1 .and. contact%ctype > 3 ) return
    if( contact%group<=0 ) return

    fstr_contact_check = .true.
  end function

  !>  Initializer of tContactState
  logical function fstr_contact_init( contact, hecMESH, cparam, myrank )
    type(tContact), intent(inout)     :: contact  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    type(tContactParam), target       :: cparam   !< contact parameter
    integer(kint),intent(in),optional :: myrank

    integer  :: i, j, is, ie, cgrp, nsurf, nslave, ic, ic_type, iss, nn, ii
    integer  :: count,nodeID

    fstr_contact_init = .false.

    contact%cparam => cparam

    !  master surface
    cgrp = contact%surf_id2
    if( cgrp<=0 ) return
    is= hecMESH%surf_group%grp_index(cgrp-1) + 1
    ie= hecMESH%surf_group%grp_index(cgrp  )

    if(present(myrank)) then
      ! PARA_CONTACT
      count = 0
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        if(hecMESH%elem_ID(ic*2) /= myrank) cycle
        count = count + 1
      enddo
      allocate( contact%master(count) )
      count = 0
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        if(hecMESH%elem_ID(ic*2) /= myrank) cycle
        count = count + 1
        nsurf = hecMESH%surf_group%grp_item(2*i)
        ic_type = hecMESH%elem_type(ic)
        call initialize_surf( ic, ic_type, nsurf, contact%master(count) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( contact%master(count)%nodes )
          nn = contact%master(count)%nodes(j)
          contact%master(count)%nodes(j) = hecMESH%elem_node_item( iss+nn )
        enddo
      enddo

    else
      ! not PARA_CONTACT
      allocate( contact%master(ie-is+1) )
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        nsurf = hecMESH%surf_group%grp_item(2*i)
        ic_type = hecMESH%elem_type(ic)
        call initialize_surf( ic, ic_type, nsurf, contact%master(i-is+1) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( contact%master(i-is+1)%nodes )
          nn = contact%master(i-is+1)%nodes(j)
          contact%master(i-is+1)%nodes(j) = hecMESH%elem_node_item( iss+nn )
        enddo
      enddo

    endif

    call update_surface_reflen( contact%master, hecMESH%node )

    cgrp = contact%surf_id1
    if( cgrp<=0 ) return
    is= hecMESH%node_group%grp_index(cgrp-1) + 1
    ie= hecMESH%node_group%grp_index(cgrp  )
    nslave = 0
    do i=is,ie
      nodeID = hecMESH%global_node_ID(hecMESH%node_group%grp_item(i))
      if(present(myrank)) then
        ! PARA_CONTACT
        nslave = nslave + 1
      else
        ! not PARA_CONTACT
        if( hecMESH%node_group%grp_item(i) <= hecMESH%nn_internal) then
          nslave = nslave + 1
        endif
      endif
    enddo
    allocate( contact%slave(nslave) )
    ii = 0
    do i=is,ie
      if(.not.present(myrank)) then
        ! not PARA_CONTACT
        if( hecMESH%node_group%grp_item(i) > hecMESH%nn_internal) cycle
      endif
      ii = ii + 1
      contact%slave(ii) = hecMESH%node_group%grp_item(i)
    enddo

    ! contact state
    allocate( contact%states(nslave) )
    do i=1,nslave
      call contact_state_init( contact%states(i) )
    enddo

    ! neighborhood of surface group
    call update_surface_box_info( contact%master, hecMESH%node )
    call bucketDB_init( contact%master_bktDB )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )
    call find_surface_neighbor( contact%master, contact%master_bktDB )

    ! initialize contact communication table
    call hecmw_contact_comm_init( contact%comm, hecMESH, 1, nslave, contact%slave )

    contact%symmetric = .true.
    fstr_contact_init = .true.
  end function

  !>  Initializer of tContactState for embed case
  logical function fstr_embed_init( embed, hecMESH, cparam, myrank )
    type(tContact), intent(inout)     :: embed  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    type(tContactParam), target       :: cparam   !< contact parameter
    integer(kint),intent(in),optional :: myrank

    integer  :: i, j, is, ie, cgrp, nsurf, nslave, ic, ic_type, iss, nn, ii
    integer  :: count,nodeID

    fstr_embed_init = .false.

    embed%cparam => cparam

    !  master surface
    cgrp = embed%surf_id2
    if( cgrp<=0 ) return
    is= hecMESH%elem_group%grp_index(cgrp-1) + 1
    ie= hecMESH%elem_group%grp_index(cgrp  )

    if(present(myrank)) then
      ! PARA_CONTACT
      count = 0
      do i=is,ie
        ic   = hecMESH%elem_group%grp_item(i)
        if(hecMESH%elem_ID(ic*2) /= myrank) cycle
        count = count + 1
      enddo
      allocate( embed%master(count) )
      count = 0
      do i=is,ie
        ic   = hecMESH%elem_group%grp_item(i)
        if(hecMESH%elem_ID(ic*2) /= myrank) cycle
        count = count + 1
        ic_type = hecMESH%elem_type(ic)
        call initialize_surf( ic, ic_type, 0, embed%master(count) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( embed%master(count)%nodes )
          nn = embed%master(count)%nodes(j)
          embed%master(count)%nodes(j) = hecMESH%elem_node_item( iss+nn )
        enddo
      enddo

    else
      ! not PARA_CONTACT
      allocate( embed%master(ie-is+1) )
      do i=is,ie
        ic   = hecMESH%elem_group%grp_item(i)
        ic_type = hecMESH%elem_type(ic)
        call initialize_surf( ic, ic_type, 0, embed%master(i-is+1) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( embed%master(i-is+1)%nodes )
          nn = embed%master(i-is+1)%nodes(j)
          embed%master(i-is+1)%nodes(j) = hecMESH%elem_node_item( iss+nn )
        enddo
      enddo

    endif

    ! slave surface
    cgrp = embed%surf_id1
    if( cgrp<=0 ) return
    is= hecMESH%node_group%grp_index(cgrp-1) + 1
    ie= hecMESH%node_group%grp_index(cgrp  )
    nslave = 0
    do i=is,ie
      nodeID = hecMESH%global_node_ID(hecMESH%node_group%grp_item(i))
      if(present(myrank)) then
        ! PARA_CONTACT
        nslave = nslave + 1
      else
        ! not PARA_CONTACT
        if( hecMESH%node_group%grp_item(i) <= hecMESH%nn_internal) then
          nslave = nslave + 1
        endif
      endif
    enddo
    allocate( embed%slave(nslave) )
    ii = 0
    do i=is,ie
      if(.not.present(myrank)) then
        ! not PARA_CONTACT
        if( hecMESH%node_group%grp_item(i) > hecMESH%nn_internal) cycle
      endif
      ii = ii + 1
      embed%slave(ii) = hecMESH%node_group%grp_item(i)
    enddo

    ! embed state
    allocate( embed%states(nslave) )
    do i=1,nslave
      call contact_state_init( embed%states(i) )
    enddo

    ! neighborhood of surface group
    call update_surface_box_info( embed%master, hecMESH%node )
    call bucketDB_init( embed%master_bktDB )
    call update_surface_bucket_info( embed%master, embed%master_bktDB )
    call find_surface_neighbor( embed%master, embed%master_bktDB )

    ! initialize contact communication table
    call hecmw_contact_comm_init( embed%comm, hecMESH, 1, nslave, embed%slave )

    ! initialize penalty coefficients
    embed%nPenalty = 1.0d0     ! default normal penalty coefficient
    embed%tPenalty = 0.1d0     ! default tangential penalty coefficient
    embed%refStiff = 0.0d0     ! will be calculated after first stiffness assembly

    embed%symmetric = .true.
    fstr_embed_init = .true.
  end function

  function check_apply_Contact_IF( contact_if, contacts )
    type(tContactInterference), intent(inout)     :: contact_if  !< contact definition
    type(tContact)     :: contacts(:) !< type fstr_solid
    
    integer  :: i, j
    logical  :: isfind
    integer(kind=kint)            :: check_apply_Contact_IF

    check_apply_Contact_IF = -1
    ! if contact pair exist?
    isfind = .false.
    do i = 1, size(contacts)
      if( contacts(i)%pair_name == contact_if%cp_name ) then
        contacts(i)%if_type     = contact_if%if_type
        contacts(i)%if_etime    = contact_if%etime
        contacts(i)%initial_pos = contact_if%initial_pos
        contacts(i)%end_pos     = contact_if%end_pos
        do j = 1, size(contacts(i)%states)
          contacts(i)%states(j)%interference_flag = contact_if%if_type
          contacts(i)%states(j)%init_pos = contact_if%initial_pos
          contacts(i)%states(j)%end_pos  = contact_if%end_pos
          if( contact_if%if_type /= C_IF_SLAVE )then
            contacts(i)%states(j)%time_factor = (contact_if%end_pos - contact_if%initial_pos) / contact_if%etime
          else
            contacts(i)%states(j)%time_factor = contact_if%etime
          end if
        end do
        isfind = .true.
        check_apply_Contact_IF = 0; return
      endif
    enddo
    if( .not. isfind ) return;
    check_apply_Contact_IF = 0

  end function

  !> Reset contact state all to free
  subroutine clear_contact_state( contact )
    type(tContact), intent(inout) :: contact    !< contact definition
    integer :: i
    if( .not. associated(contact%states) ) return
    do i=1,size( contact%states )
      contact%states(i)%state = -1
    enddo
  end subroutine

  !> if contact is active is curr step
  logical function is_active_contact( acgrp, contact )
    integer, intent(in)        :: acgrp(:)      !< active contact group numbers
    type(tContact), intent(in) :: contact       !< contact definition
    if( any( acgrp==contact%group ) ) then
      is_active_contact = .true.
    else
      is_active_contact = .false.
    endif
  end function

  !> Write out the contact definition read from mesh file
  subroutine print_contatct_pair( file, pair )
    integer(kind=kint), intent(in)           :: file
    type( hecmwST_contact_pair ), intent(in) :: pair

    integer(kind=kint) :: i
    write(file,*) "Number of contact pair", pair%n_pair
    do i=1,pair%n_pair
      write(file,*) trim(pair%name(i)), pair%type(i), pair%slave_grp_id(i)  &
        ,pair%master_grp_id(i), pair%slave_orisgrp_id(i)
    enddo
  end subroutine


end module mContactDef
