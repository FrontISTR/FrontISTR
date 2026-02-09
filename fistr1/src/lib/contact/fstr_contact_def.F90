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
  integer, parameter :: CONTACTACTIVE = 10
  integer, parameter :: CONTACTINACTIVE = -1 !!!! dbg-later
  
  !> contact type or algorithm definition
  integer, parameter :: CONTACTTIED = 1
  integer, parameter :: CONTACTGLUED = 2
  integer, parameter :: CONTACTSSLID = 3
  integer, parameter :: CONTACTFSLID = 4

  !> contact method
  integer, parameter :: CONTACTN2S = 1
  integer, parameter :: CONTACTS2S = 2

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

  !> Structure to define contact surface
  type tContactSurf
    integer(kind=kint)              :: eid                  !< elemental index(global)
    integer(kind=kint)              :: etype                !< type of surface element
    integer(kind=kint), pointer     :: nodes(:)=>null()     !< nodes index(global)
    integer(kind=kint)              :: state
    integer(kind=kint)              :: n_intp               !< num of surface integral point
    type(tContactState), pointer    :: states(:)=>null()    !< contact states of slave surf
    integer(kind=kint), pointer     :: nslave_index(:)=>null()     !< nodes index(global)
    real(kind=kreal), pointer       :: phi(:,:)=>null()  
    integer(kind=kint), pointer     :: msurf_list(:)=>null()     !< nodes index(global)
  end type tContactSurf

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
    
    type(tContactSurf), pointer   :: slave_surf(:)=>null()       

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
    integer                       :: method                  !< ndoe-surf or surf-surf

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
    integer, allocatable  :: slave_index(:)
    allocate( slave_index(hecMESH%n_node) )
    slave_index(:) = 0
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
      slave_index(contact%slave(ii)) = ii
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

    if(contact%method == CONTACTS2S) then
      !  slave surface
      cgrp = contact%surf_id1_sgrp
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
        allocate( contact%slave_surf(count) )
        count = 0
        do i=is,ie
          ic   = hecMESH%surf_group%grp_item(2*i-1)
          if(hecMESH%elem_ID(ic*2) /= myrank) cycle
          count = count + 1
          nsurf = hecMESH%surf_group%grp_item(2*i)
          ic_type = hecMESH%elem_type(ic)
          call initialize_csurf( ic, ic_type, nsurf, contact%slave_surf(count) )
          iss = hecMESH%elem_node_index(ic-1)
          do j=1, size( contact%slave_surf(count)%nodes )
            nn = contact%slave_surf(count)%nodes(j)
            contact%slave_surf(count)%nodes(j) = hecMESH%elem_node_item( iss+nn )
            contact%slave_surf(count)%nslave_index(j) = slave_index(hecMESH%elem_node_item( iss+nn ))
          enddo
        enddo

      else
        ! not PARA_CONTACT
        allocate( contact%slave_surf(ie-is+1) )
        do i=is,ie
          ic   = hecMESH%surf_group%grp_item(2*i-1)
          nsurf = hecMESH%surf_group%grp_item(2*i)
          ic_type = hecMESH%elem_type(ic)
          call initialize_csurf( ic, ic_type, nsurf, contact%slave_surf(i-is+1) )
          iss = hecMESH%elem_node_index(ic-1)
          do j=1, size( contact%slave_surf(i-is+1)%nodes )
            nn = contact%slave_surf(i-is+1)%nodes(j)
            contact%slave_surf(i-is+1)%nodes(j) = hecMESH%elem_node_item( iss+nn )
            contact%slave_surf(i-is+1)%nslave_index(j) = slave_index(hecMESH%elem_node_item( iss+nn ))
          enddo
        enddo

      endif

      ! state for each integration points
      do i=1, size( contact%slave_surf )
        nn = contact%slave_surf(i)%n_intp
        allocate( contact%slave_surf(i)%states(nn) )
        do j = 1, contact%slave_surf(i)%n_intp
          contact%slave_surf(i)%states(j)%state = -1
          contact%slave_surf(i)%states(j)%multiplier(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce1(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce_trial(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce_final(:) = 0.d0
          contact%slave_surf(i)%states(j)%reldisp(:) = 0.d0
          contact%slave_surf(i)%states(j)%time_factor = 0.d0
          contact%slave_surf(i)%states(j)%interference_flag = 0
        enddo
      enddo
    endif

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

  !> Initializer
  subroutine initialize_csurf( eid, etype, nsurf, surf )
    use elementInfo
    integer(kind=kint), intent(in)    :: eid    !< element ID
    integer(kind=kint), intent(in)    :: etype  !< element type
    integer(kind=kint), intent(in)    :: nsurf  !< surface ID
    type(tContactSurf), intent(inout) :: surf   !< surface element
    integer(kind=kint) :: i, n, outtype, nodes(100)
    surf%eid = eid

    call getSubFace( etype, nsurf, outtype, nodes )
    surf%etype = outtype
    n=getNumberOfNodes( outtype )
    if(surf%etype == fe_quad4n )then
      ! surf%n_intp = 1
      ! surf%n_intp = 4
      surf%n_intp = 16
      ! surf%n_intp = 100
    else if(surf%etype == fe_tri3n )then
      surf%n_intp = 7
    end if
    allocate( surf%nodes(n) )
    allocate( surf%nslave_index(n) )
    surf%nodes(1:n)=nodes(1:n)
    surf%nslave_index(:)= 0

    allocate( surf%phi(surf%n_intp,n) )
    allocate( surf%msurf_list(surf%n_intp) )
    surf%msurf_list(:) = -1
  end subroutine
  
  subroutine getDualShapfunc(slave_surf, nn, n_intp, elecoord, weight, N, success)
    type(tContactSurf), intent(inout):: slave_surf
    integer, intent(in)           :: nn, n_intp
    real(kind=kreal), intent(in)  :: elecoord(:,:), weight(:), N(:,:)
    integer :: intp, i, j, k, info,counter
    real(kind=kreal) :: Me(nn,nn), Ae(nn,nn), ncoord(2)
    real(kind=kreal) :: De(nn), det(n_intp)
    real(kind=kreal) :: biorth_check(nn,nn), sum_check
    logical, intent(out) :: success
    logical, parameter :: debug_check = .true.
    
    ! 積分点での形状関数とヤコビアンを計算
    Ae = 0.d0

    ! De = ∫ N_i dA を計算（対角ベクトル）
    De = 0.0d0
    success = .true.
    do j = 1, nn
      do intp = 1, n_intp
        if( slave_surf%states(intp)%state == CONTACTFREE ) cycle
        De(j) = De(j) + N(intp,j) * weight(intp)
      end do
    end do
    ! write(*,*) 'De = ∫ N_i dA:'
    ! do i = 1, nn
    !   write(*,'(A,I2,A,4E15.6)') '  De(', i, ') = ', De(i)
    ! end do
    ! Me = ∫ N_j N_k dA を計算
    Me = 0.0d0
    do j = 1, nn
      do k = 1, nn
        do intp = 1, n_intp
          if( slave_surf%states(intp)%state == CONTACTFREE ) cycle
          Me(j,k) = Me(j,k) + N(intp,j)*N(intp,k)* weight(intp)
        end do
      end do
    end do
    ! write(*,*) 'Me_jk = ∫ N_j N_k dA (m_els):'
    ! do i = 1, nn
    !   write(*,'(A,I2,A,4E15.6)') '  Me(', i, ') = ', Me(i,:)
    ! end do

    ! 3. Me の逆行列を計算（対称行列用 Cholesky 分解）
    !    Me は対称正定値なので dpotrf/dpotri を使用
    call dpotrf('U', nn, Me, nn, info)
    if (info /= 0) then
      write(*,*) 'Error in dpotrf (Cholesky factorization), info=', info
      write(*,*) 'Matrix may not be positive definite'
      success = .false.
      return
      stop
    endif
    
    call dpotri('U', nn, Me, nn, info)
    if (info /= 0) then
      write(*,*) 'Error in dpotri (inverse calculation), info=', info
      success = .false.
      return
      stop
    endif
    
    ! dpotri は上三角のみ計算するので、下三角を埋める
    do j = 1, nn
      do k = j+1, nn
        Me(k,j) = Me(j,k)
      end do
    end do

    ! Ae = W_j × (Me^-1)_jk を計算
    !    これが Dual basis 関数の展開係数
    do j = 1, nn
      do k = 1, nn
        Ae(j,k) = De(j) * Me(j,k)
      end do
    end do

    ! Ae(1,:) = Ae(1,:)+Ae(2,:)
    ! Ae(4,:) = Ae(4,:)+Ae(3,:)

    ! Ae(2,:) = 0.0d0
    ! Ae(3,:) = 0.0d0

    ! write(*,*) 'Ae (pslavdual):'
    ! do i = 1, nn
    !   write(*,'(A,I2,A,4E15.6)') '  Ae(', i, ') = ', Ae(i,:)
    ! end do

    ! φ_j(積分点) = Σ_k Ae_jk N_k(積分点) を計算
    !    各積分点で Dual basis 関数を評価
    slave_surf%phi = 0.0d0
    do intp = 1, n_intp
      if( slave_surf%states(intp)%state == CONTACTFREE ) cycle
      do j = 1, nn
        do k = 1, nn
          slave_surf%phi(intp,j) = slave_surf%phi(intp,j) + Ae(j,k) * N(intp,k)
        end do
      end do
    end do

    ! 6. Biorthogonality 条件の検証（デバッグ用）
    !    ∫ φ_i N_j dA ≈ W_i δ_ij を確認
    if (debug_check) then
      ! write(*,*) '=== Dual Shape Function Verification ==='
      ! write(*,*) 'Number of nodes:', nn
      ! write(*,*) 'Number of integration points:', n_intp
      
      ! Biorthogonality check: ∫ φ_i N_j dA
      ! W と M の計算と同じ積分点のみを使用
      biorth_check = 0.0d0
      do intp = 1, n_intp
        do i = 1, nn
          do j = 1, nn
            biorth_check(i,j) = biorth_check(i,j) + &
                                slave_surf%phi(intp,i) * N(intp,j) * weight(intp)
          end do
        end do
      end do
      
      ! write(*,*) 'Biorthogonality check (∫ φ_i N_j dA):'
      ! write(*,*) '  (Should be W_i for i=j, zero for i≠j)'
      ! do i = 1, nn
      !   write(*,'(A,I2,A)', advance='no') '  Row ', i, ':'
      !   do j = 1, nn
      !     write(*,'(E12.4)', advance='no') biorth_check(i,j)
      !   end do
      !   write(*,*)
      ! end do
      
      ! 最大誤差を計算
      ! write(*,*) 'Maximum errors:'
      do i = 1, nn
        ! 対角成分の誤差（W_i との差）
        ! write(*,'(A,I2,A,E15.6)') '  Diagonal(', i, ') error: ', &
        !                           abs(biorth_check(i,i) - De(i))
        ! 非対角成分の誤差（ゼロとの差）
        do j = 1, nn
          if (i /= j) then
            if (abs(biorth_check(i,j)) > 1.0d-8) then
              write(*,'(A,I2,A,I2,A,E15.6)') '  Off-diagonal(', i, ',', j, &
                                             ') error: ', abs(biorth_check(i,j))
              success = .false.
              return
            endif
          endif
        end do
      end do
      
      ! Partition of unity check: Σ_i φ_i = 1
      ! W と M の計算と同じ積分点のみを使用
      ! write(*,*) 'Partition of unity check (Σ φ_i at each integration point):'
      ! do intp = 1, n_intp
      !   if( slave_surf%states(intp)%state == CONTACTFREE ) cycle
      !   sum_check = 0.0d0
      !   do i = 1, nn
      !     sum_check = sum_check + slave_surf%phi(intp,i)
      !   end do
      !   write(*,'(A,I3,A,E15.6,A,E15.6)') '  Intp ', intp, ': sum = ', &
      !                                     sum_check, ', error = ', abs(sum_check - 1.0d0)
      ! end do
      
      ! write(*,*) '========================================='
    endif

  end subroutine

  subroutine switch_intp_contact(slave_surf, num, nstate)
    type(tContactSurf), intent(inout):: slave_surf
    integer, intent(in)           :: num, nstate
    integer :: n_intp, i, intp_list(25), state
    logical, parameter :: debug_check = .true.

    n_intp = 4
    if (nstate == CONTACTFREE) then
      state = CONTACTINACTIVE
    else
      state = CONTACTACTIVE
      slave_surf%state = CONTACTACTIVE
    end if

    select case (slave_surf%etype)
      case ( fe_tri3n )
        select case (num)
          case (1)
            intp_list(1:4) = [1,4,5,6]
          case (2)
            intp_list(1:4) = [1,2,5,7]
          case (3)
            intp_list(1:4) = [1,3,6,7]
        end select
      case ( fe_quad4n )
        select case (num)
          case (1)
            intp_list(1:4) = [1,2,5,6]
          case (2)
            intp_list(1:4) = [3,4,7,8]
          case (3)
            intp_list(1:4) = [11,12,15,16]
          case (4)
            intp_list(1:4) = [9,10,13,14]
        end select
      case default
        ! error message
        stop "element type not defined-qp"
    end select

    do i = 1, n_intp
      if (slave_surf%states(intp_list(i))%state > 0 .and. state > 0) cycle !!!! dbg-later
      slave_surf%states(intp_list(i))%state = state
    end do

  end subroutine switch_intp_contact

end module mContactDef
