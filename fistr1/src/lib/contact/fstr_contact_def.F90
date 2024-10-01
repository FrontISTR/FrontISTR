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
  use mSurfElement
  use m_contact_lib
  use m_hecmw_contact_comm
  use bucket_search
  use mContactParam

  implicit none

  include 'fstr_ctrl_util_f.inc'

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
    real(kind=kreal)              :: tPenalty                !< tangential penalty

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
  private :: cal_node_normal
  private :: track_contact_position
  private :: reset_contact_force

contains

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
    allocate( contact%states(nslave) )
    ii = 0
    do i=is,ie
      if(.not.present(myrank)) then
        ! not PARA_CONTACT
        if( hecMESH%node_group%grp_item(i) > hecMESH%nn_internal) cycle
      endif
      ii = ii + 1
      contact%slave(ii) = hecMESH%node_group%grp_item(i)
      contact%states(ii)%state = -1
      contact%states(ii)%multiplier(:) = 0.d0
      contact%states(ii)%tangentForce(:) = 0.d0
      contact%states(ii)%tangentForce1(:) = 0.d0
      contact%states(ii)%tangentForce_trial(:) = 0.d0
      contact%states(ii)%tangentForce_final(:) = 0.d0
      contact%states(ii)%reldisp(:) = 0.d0
    enddo
    !    endif

    ! contact state
    !      allocate( contact%states(nslave) )
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

    !call update_surface_reflen( embed%master, hecMESH%node )

    ! slave surface
    !    if( contact%ctype==1 ) then
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
    allocate( embed%states(nslave) )
    ii = 0
    do i=is,ie
      if(.not.present(myrank)) then
        ! not PARA_CONTACT
        if( hecMESH%node_group%grp_item(i) > hecMESH%nn_internal) cycle
      endif
      ii = ii + 1
      embed%slave(ii) = hecMESH%node_group%grp_item(i)
      embed%states(ii)%state = -1
      embed%states(ii)%multiplier(:) = 0.d0
      embed%states(ii)%tangentForce(:) = 0.d0
      embed%states(ii)%tangentForce1(:) = 0.d0
      embed%states(ii)%tangentForce_trial(:) = 0.d0
      embed%states(ii)%tangentForce_final(:) = 0.d0
      embed%states(ii)%reldisp(:) = 0.d0
    enddo
    !    endif

    ! contact state
    !      allocate( contact%states(nslave) )
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

    embed%symmetric = .true.
    fstr_embed_init = .true.
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

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_contact_state( flag_ctAlgo, contact, currpos, currdisp, ndforce, infoCTChange, &
      nodeID, elemID, is_init, active, mu, B )
    character(len=9), intent(in)                     :: flag_ctAlgo  !< contact analysis algorithm flag
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)  !< current displacement of each nodes
    real(kind=kreal), intent(in)                     :: ndforce(:)   !< nodal force
    integer(kind=kint), intent(in)                   :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                   :: elemID(:)    !< global elemental ID, just for print out
    logical, intent(in)                              :: is_init      !< whether initial scan or not
    logical, intent(out)                             :: active       !< if any in contact
    real(kind=kreal), intent(in)                     :: mu           !< penalty
    real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node )
    real(kind=kreal)    :: nlforce, slforce(3)
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_cand, is_present_B
    real(kind=kreal), pointer :: Bp(:)

    if( is_init ) then
      distclr = contact%cparam%DISTCLR_INIT
    else
      distclr = contact%cparam%DISTCLR_FREE
      if( contact%algtype == CONTACTTIED ) then
        active = .false.
        do i= 1, size(contact%slave)
          if( contact%states(i)%state==CONTACTSTICK ) then
            active = .true.
            exit
          endif
        enddo
      endif
    endif

    allocate(contact_surf(size(nodeID)))
    allocate(states_prev(size(contact%slave)))
    contact_surf(:) = huge(0)
    do i = 1, size(contact%slave)
      states_prev(i) = contact%states(i)%state
    enddo

    call update_surface_box_info( contact%master, currpos )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )

    ! for gfortran-10: optional parameter seems not allowed within omp parallel
    is_present_B = present(B)
    if( is_present_B ) Bp => B

    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,slforce,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,elem,is_cand,idm,etype,isin, &
      !$omp&         bktID,nCand,indexCand) &
      !$omp& firstprivate(nMasterMax,is_present_B) &
      !$omp& shared(contact,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,mu,nodeID,elemID,Bp,distclr,contact_surf,is_init) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      slave = contact%slave(i)
      if( contact%states(i)%state==CONTACTSTICK .or. contact%states(i)%state==CONTACTSLIP ) then
        slforce(1:3)=ndforce(3*slave-2:3*slave)
        id = contact%states(i)%surface
        nlforce = contact%states(i)%multiplier(1)

        ! update direction of TIED contact
        if( contact%algtype == CONTACTTIED ) then
          call update_direction( i, contact, currpos )
          if (.not.is_init) cycle
        endif

        if( nlforce < contact%cparam%TENSILE_FORCE ) then
          contact%states(i)%state = CONTACTFREE
          contact%states(i)%multiplier(:) = 0.d0
          write(*,'(A,i10,A,i10,A,e12.3)') "Node",nodeID(slave)," free from contact with element", &
            elemID(contact%master(id)%eid), " with tensile force ", nlforce
          cycle
        endif
        if( contact%algtype /= CONTACTFSLID .or. (.not. is_present_B) ) then   ! small slide problem
          contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
        else
          call track_contact_position( flag_ctAlgo, i, contact, currpos, currdisp, mu, infoCTChange, nodeID, elemID, Bp )
          if( contact%states(i)%state /= CONTACTFREE ) then
            id = contact%states(i)%surface
            contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
          endif
        endif

      else if( contact%states(i)%state==CONTACTFREE ) then
        if( contact%algtype == CONTACTTIED .and. .not. is_init ) cycle 
        coord(:) = currpos(3*slave-2:3*slave)

        ! get master candidates from bucketDB
        bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
        nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
        if (nCand == 0) cycle
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)

        nMasterMax = nCand
        allocate(indexMaster(nMasterMax))
        nMaster = 0

        ! narrow down candidates
        do idm= 1, nCand
          id = indexCand(idm)
          if (.not. is_in_surface_box( contact%master(id), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          nMaster = nMaster + 1
          indexMaster(nMaster) = id
        enddo
        deallocate(indexCand)

        if(nMaster == 0) then
          deallocate(indexMaster)
          cycle
        endif

        do idm = 1,nMaster
          id = indexMaster(idm)
          etype = contact%master(id)%etype
          nn = size( contact%master(id)%nodes )
          do j=1,nn
            iSS = contact%master(id)%nodes(j)
            elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          enddo
          call project_Point2Element( coord,etype,nn,elem,contact%master(id)%reflen,contact%states(i), &
            isin,distclr,localclr=contact%cparam%CLEARANCE )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          iSS = isInsideElement( etype, contact%states(i)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
          if( iSS>0 ) &
            call cal_node_normal( id, iSS, contact%master, currpos, contact%states(i)%lpos(1:2), &
              contact%states(i)%direction(:) )
          contact_surf(contact%slave(i)) = elemID(contact%master(id)%eid)
          write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3,A,i6)') "Node",nodeID(slave)," contact with element", &
            elemID(contact%master(id)%eid),       &
            " with distance ", contact%states(i)%distance," at ",contact%states(i)%lpos(1:2), &
            " along direction ", contact%states(i)%direction," rank=",hecmw_comm_get_rank()
          exit
        enddo
        deallocate(indexMaster)
      endif
    enddo
    !$omp end parallel do

    if( contact%algtype == CONTACTTIED .and. .not. is_init ) then
      deallocate(contact_surf)
      deallocate(states_prev)
      return
    endif

    call hecmw_contact_comm_allreduce_i(contact%comm, contact_surf, HECMW_MIN)
    nactive = 0
    do i = 1, size(contact%slave)
      if (contact%states(i)%state /= CONTACTFREE) then                    ! any slave in contact
        id = contact%states(i)%surface
        if (abs(contact_surf(contact%slave(i))) /= elemID(contact%master(id)%eid)) then ! that is in contact with other surface
          contact%states(i)%state = CONTACTFREE                           ! should be freed
          write(*,'(A,i10,A,i10,A,i6,A,i6,A)') "Node",nodeID(contact%slave(i))," contact with element", &
            &  elemID(contact%master(id)%eid), " in rank",hecmw_comm_get_rank()," freed due to duplication"
        else
          nactive = nactive + 1
        endif
      endif
      if (states_prev(i) == CONTACTFREE .and. contact%states(i)%state /= CONTACTFREE) then
        infoCTChange%free2contact = infoCTChange%free2contact + 1
      elseif (states_prev(i) /= CONTACTFREE .and. contact%states(i)%state == CONTACTFREE) then
        infoCTChange%contact2free = infoCTChange%contact2free + 1
      endif
    enddo
    active = (nactive > 0)
    deallocate(contact_surf)
    deallocate(states_prev)
  end subroutine scan_contact_state

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_embed_state( flag_ctAlgo, embed, currpos, currdisp, ndforce, infoCTChange, &
    nodeID, elemID, is_init, active, mu, B )
  character(len=9), intent(in)                     :: flag_ctAlgo  !< contact analysis algorithm flag
  type( tContact ), intent(inout)                  :: embed      !< contact info
  type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
  real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
  real(kind=kreal), intent(in)                     :: currdisp(:)  !< current displacement of each nodes
  real(kind=kreal), intent(in)                     :: ndforce(:)   !< nodal force
  integer(kind=kint), intent(in)                   :: nodeID(:)    !< global nodal ID, just for print out
  integer(kind=kint), intent(in)                   :: elemID(:)    !< global elemental ID, just for print out
  logical, intent(in)                              :: is_init      !< whether initial scan or not
  logical, intent(out)                             :: active       !< if any in contact
  real(kind=kreal), intent(in)                     :: mu           !< penalty
  real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

  real(kind=kreal)    :: distclr
  integer(kind=kint)  :: slave, id, etype
  integer(kind=kint)  :: nn, i, j, iSS, nactive
  real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node )
  real(kind=kreal)    :: nlforce, slforce(3)
  logical             :: isin
  integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
  !
  integer, pointer :: indexMaster(:),indexCand(:)
  integer   ::  nMaster,idm,nMasterMax,bktID,nCand
  logical :: is_cand, is_present_B
  real(kind=kreal), pointer :: Bp(:)

  if( is_init ) then
    distclr = embed%cparam%DISTCLR_INIT
  else
    distclr = embed%cparam%DISTCLR_FREE
    active = .false.
    do i= 1, size(embed%slave)
      if( embed%states(i)%state==CONTACTSTICK ) then
        active = .true.
        exit
      endif
    enddo
    return
  endif

  allocate(contact_surf(size(nodeID)))
  allocate(states_prev(size(embed%slave)))
  contact_surf(:) = huge(0)
  do i = 1, size(embed%slave)
    states_prev(i) = embed%states(i)%state
  enddo

  call update_surface_box_info( embed%master, currpos )
  call update_surface_bucket_info( embed%master, embed%master_bktDB )

  ! for gfortran-10: optional parameter seems not allowed within omp parallel
  is_present_B = present(B)
  if( is_present_B ) Bp => B

  !$omp parallel do &
    !$omp& default(none) &
    !$omp& private(i,slave,slforce,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,elem,is_cand,idm,etype,isin, &
    !$omp&         bktID,nCand,indexCand) &
    !$omp& firstprivate(nMasterMax,is_present_B) &
    !$omp& shared(embed,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,mu,nodeID,elemID,Bp,distclr,contact_surf) &
    !$omp& reduction(.or.:active) &
    !$omp& schedule(dynamic,1)
  do i= 1, size(embed%slave)
    slave = embed%slave(i)
    if( embed%states(i)%state==CONTACTFREE ) then
      coord(:) = currpos(3*slave-2:3*slave)

      ! get master candidates from bucketDB
      bktID = bucketDB_getBucketID(embed%master_bktDB, coord)
      nCand = bucketDB_getNumCand(embed%master_bktDB, bktID)
      if (nCand == 0) cycle
      allocate(indexCand(nCand))
      call bucketDB_getCand(embed%master_bktDB, bktID, nCand, indexCand)

      nMasterMax = nCand
      allocate(indexMaster(nMasterMax))
      nMaster = 0

      ! narrow down candidates
      do idm= 1, nCand
        id = indexCand(idm)
        if (.not. is_in_surface_box( embed%master(id), coord(1:3), embed%cparam%BOX_EXP_RATE )) cycle
        nMaster = nMaster + 1
        indexMaster(nMaster) = id
      enddo
      deallocate(indexCand)

      if(nMaster == 0) then
        deallocate(indexMaster)
        cycle
      endif

      do idm = 1,nMaster
        id = indexMaster(idm)
        etype = embed%master(id)%etype
        if( mod(etype,10) == 2 ) etype = etype - 1 !search by 1st-order shape function
        nn = getNumberOfNodes(etype)
        do j=1,nn
          iSS = embed%master(id)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call project_Point2SolidElement( coord,etype,nn,elem,embed%master(id)%reflen,embed%states(i), &
          isin,distclr,localclr=embed%cparam%CLEARANCE )
        if( .not. isin ) cycle
        embed%states(i)%surface = id
        embed%states(i)%multiplier(:) = 0.d0
        contact_surf(embed%slave(i)) = elemID(embed%master(id)%eid)
        write(*,'(A,i10,A,i10,A,3f7.3,A,i6)') "Node",nodeID(slave)," embeded to element", &
          elemID(embed%master(id)%eid), " at ",embed%states(i)%lpos(:)," rank=",hecmw_comm_get_rank()
        exit
      enddo
      deallocate(indexMaster)
    endif
  enddo
  !$omp end parallel do

  call hecmw_contact_comm_allreduce_i(embed%comm, contact_surf, HECMW_MIN)
  nactive = 0
  do i = 1, size(embed%slave)
    if (embed%states(i)%state /= CONTACTFREE) then                    ! any slave in contact
      id = embed%states(i)%surface
      if (abs(contact_surf(embed%slave(i))) /= elemID(embed%master(id)%eid)) then ! that is in contact with other surface
        embed%states(i)%state = CONTACTFREE                           ! should be freed
        write(*,'(A,i10,A,i10,A,i6,A,i6,A)') "Node",nodeID(embed%slave(i))," contact with element", &
          &  elemID(embed%master(id)%eid), " in rank",hecmw_comm_get_rank()," freed due to duplication"
      else
        nactive = nactive + 1
      endif
    endif
    if (states_prev(i) == CONTACTFREE .and. embed%states(i)%state /= CONTACTFREE) then
      infoCTChange%free2contact = infoCTChange%free2contact + 1
    elseif (states_prev(i) /= CONTACTFREE .and. embed%states(i)%state == CONTACTFREE) then
      infoCTChange%contact2free = infoCTChange%contact2free + 1
    endif
  enddo
  active = (nactive > 0)
  deallocate(contact_surf)
  deallocate(states_prev)
end subroutine scan_embed_state


  !> Calculate averaged nodal normal
  subroutine cal_node_normal( csurf, isin, surf, currpos, lpos, normal )
    use elementInfo, only:getVertexCoord, SurfaceNormal
    integer, intent(in)            :: csurf       !< current surface element
    integer, intent(in)            :: isin        !< return value from isInsideElement()
    type(tSurfElement), intent(in) :: surf(:)     !< surface elements
    real(kind=kreal), intent(in)   :: currpos(:)  !< current coordinate of each nodes
    real(kind=kreal), intent(in)   :: lpos(:)     !< local coordinate of contact position
    real(kind=kreal), intent(out)  :: normal(3)   !< averaged node nomral
    integer(kind=kint) :: cnode, i, j, cnt, nd1, gn, etype, iSS, nn,cgn
    real(kind=kreal) :: cnpos(2), elem(3, l_max_elem_node )
    integer(kind=kint) :: cnode1, cnode2, gn1, gn2, nsurf, cgn1, cgn2, isin_n
    real(kind=kreal) :: x=0, normal_n(3), lpos_n(2)

    if( 1 <= isin .and. isin <= 4 ) then  ! corner
      cnode = isin
      gn = surf(csurf)%nodes(cnode)
      etype = surf(csurf)%etype
      call getVertexCoord( etype, cnode, cnpos )
      nn = size( surf(csurf)%nodes )
      do j=1,nn
        iSS = surf(csurf)%nodes(j)
        elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      enddo
      normal = SurfaceNormal( etype, nn, cnpos, elem )
      cnt = 1
      do i=1,surf(csurf)%n_neighbor
        nd1 = surf(csurf)%neighbor(i)
        nn = size( surf(nd1)%nodes )
        etype = surf(nd1)%etype
        cgn = 0
        do j=1,nn
          iSS = surf(nd1)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          if( iSS==gn ) cgn=j
        enddo
        if( cgn>0 ) then
          call getVertexCoord( etype, cgn, cnpos )
          !normal = normal+SurfaceNormal( etype, nn, cnpos, elem )
          normal_n = SurfaceNormal( etype, nn, cnpos, elem )
          normal = normal+normal_n
          cnt = cnt+1
        endif
      enddo
      !normal = normal/cnt                                        !!-???
    elseif( 12 <= isin .and. isin <= 41 ) then  ! edge
      cnode1 = isin / 10
      cnode2 = mod(isin, 10)
      gn1 = surf(csurf)%nodes(cnode1)
      gn2 = surf(csurf)%nodes(cnode2)
      etype = surf(csurf)%etype
      nn = size( surf(csurf)%nodes )
      do j=1,nn
        iSS = surf(csurf)%nodes(j)
        elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      enddo
      normal = SurfaceNormal( etype, nn, lpos, elem )
      select case (etype)
      case (fe_tri3n, fe_tri6n, fe_tri6nc)
        if    ( isin==12 ) then; x=lpos(2)-lpos(1)
        elseif( isin==23 ) then; x=1.d0-2.d0*lpos(2)
        elseif( isin==31 ) then; x=2.d0*lpos(1)-1.d0
        else; stop "Error: cal_node_normal: invalid isin"
        endif
      case (fe_quad4n, fe_quad8n)
        if    ( isin==12 ) then; x=lpos(1)
        elseif( isin==23 ) then; x=lpos(2)
        elseif( isin==34 ) then; x=-lpos(1)
        elseif( isin==41 ) then; x=-lpos(2)
        else; stop "Error: cal_node_normal: invalid isin"
        endif
      end select
      ! find neighbor surf that includes cnode1 and cnode2
      nsurf = 0
      NEIB_LOOP: do i=1, surf(csurf)%n_neighbor
        nd1 = surf(csurf)%neighbor(i)
        nn = size( surf(nd1)%nodes )
        etype = surf(nd1)%etype
        cgn1 = 0
        cgn2 = 0
        do j=1,nn
          iSS = surf(nd1)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          if( iSS==gn1 ) cgn1=j
          if( iSS==gn2 ) cgn2=j
        enddo
        if( cgn1>0 .and. cgn2>0 ) then
          nsurf = nd1
          isin_n = 10*cgn2 + cgn1
          x = -x
          select case (etype)
          case (fe_tri3n, fe_tri6n, fe_tri6nc)
            if    ( isin_n==12 ) then; lpos_n(1)=0.5d0*(1.d0-x); lpos_n(2)=0.5d0*(1.d0+x)
            elseif( isin_n==23 ) then; lpos_n(1)=0.d0;           lpos_n(2)=0.5d0*(1.d0-x)
            elseif( isin_n==31 ) then; lpos_n(1)=0.5d0*(1.d0+x); lpos_n(2)=0.d0
            else; stop "Error: cal_node_normal: invalid isin_n"
            endif
          case (fe_quad4n, fe_quad8n)
            if    ( isin_n==12 ) then; lpos_n(1)= x;    lpos_n(2)=-1.d0
            elseif( isin_n==23 ) then; lpos_n(1)= 1.d0; lpos_n(2)= x
            elseif( isin_n==34 ) then; lpos_n(1)=-x;    lpos_n(2)= 1.d0
            elseif( isin_n==41 ) then; lpos_n(1)=-1.d0; lpos_n(2)=-x
            else; stop "Error: cal_node_normal: invalid isin_n"
            endif
          end select
          !normal = normal + SurfaceNormal( etype, nn, lpos_n, elem )
          normal_n = SurfaceNormal( etype, nn, lpos_n, elem )
          normal = normal+normal_n
          exit NEIB_LOOP
        endif
      enddo NEIB_LOOP
      !if( nsurf==0 ) write(0,*) "Warning: cal_node_normal: neighbor surf not found"
      !normal = normal/2
    endif
    normal = normal/ dsqrt( dot_product( normal, normal ) )
  end subroutine cal_node_normal

  !> This subroutine tracks down next contact position after a finite slide
  subroutine track_contact_position( flag_ctAlgo, nslave, contact, currpos, currdisp, mu, infoCTChange, nodeID, elemID, B )
    character(len=9), intent(in)                    :: flag_ctAlgo  !< contact analysis algorithm flag
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)    !< current displacement of each nodes
    real(kind=kreal), intent(in)                     :: mu           !< penalty
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out
    real(kind=kreal), intent(inout)                  :: B(:)         !< nodal force residual

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node ), elem0(3, l_max_elem_node )
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)

    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    odirec = contact%states(nslave)%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2Element( coord,etype,nn,elem,contact%master(sid0)%reflen,contact%states(nslave), &
      isin,contact%cparam%DISTCLR_NOCHECK,contact%states(nslave)%lpos(1:2),contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        etype = contact%master(sid)%etype
        nn = getNumberOfNodes( etype )
        do j=1,nn
          iSS = contact%master(sid)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call project_Point2Element( coord,etype,nn,elem,contact%master(sid)%reflen,contact%states(nslave), &
          isin,contact%cparam%DISTCLR_NOCHECK,localclr=contact%cparam%CLEARANCE )
        if( isin ) then
          contact%states(nslave)%surface = sid
          exit
        endif
      enddo
    endif

    if( .not. isin ) then   ! such case is considered to rarely or never occur
      write(*,*) 'Warning: contact moved beyond neighbor elements'
      ! get master candidates from bucketDB
      bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
      nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
      if (nCand > 0) then
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)
        do idm= 1, nCand
          sid = indexCand(idm)
          if( sid==sid0 ) cycle
          if( associated(contact%master(sid0)%neighbor) ) then
            if( any(sid==contact%master(sid0)%neighbor(:)) ) cycle
          endif
          if (.not. is_in_surface_box( contact%master(sid), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          etype = contact%master(sid)%etype
          nn = size( contact%master(sid)%nodes )
          do j=1,nn
            iSS = contact%master(sid)%nodes(j)
            elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          enddo
          call project_Point2Element( coord,etype,nn,elem,contact%master(sid)%reflen,contact%states(nslave), &
               isin,contact%cparam%DISTCLR_NOCHECK,localclr=contact%cparam%CLEARANCE )
          if( isin ) then
            contact%states(nslave)%surface = sid
            exit
          endif
        enddo
        deallocate(indexCand)
      endif
    endif

    if( isin ) then
      if( contact%states(nslave)%surface==sid0 ) then
        if(any(dabs(contact%states(nslave)%lpos(1:2)-opos(:)) >= contact%cparam%CLR_DIFFLPOS))  then
          !$omp atomic
          infoCTChange%contact2difflpos = infoCTChange%contact2difflpos + 1
        endif
      else
        write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3)') "Node",nodeID(slave)," move to contact with", &
          elemID(contact%master(sid)%eid), " with distance ",      &
          contact%states(nslave)%distance," at ",contact%states(nslave)%lpos(1:2)
        !$omp atomic
        infoCTChange%contact2neighbor = infoCTChange%contact2neighbor + 1
        if( flag_ctAlgo=='ALagrange' )  &
          call reset_contact_force( contact, currpos, nslave, sid0, opos, odirec, B )
      endif
      if( flag_ctAlgo=='SLagrange' ) call update_TangentForce(etype,nn,elem0,elem,contact%states(nslave))
      iSS = isInsideElement( etype, contact%states(nslave)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) &
        call cal_node_normal( contact%states(nslave)%surface, iSS, contact%master, currpos, &
          contact%states(nslave)%lpos(1:2), contact%states(nslave)%direction(:) )
    else if( .not. isin ) then
      write(*,'(A,i10,A)') "Node",nodeID(slave)," move out of contact"
      contact%states(nslave)%state = CONTACTFREE
      contact%states(nslave)%multiplier(:) = 0.d0
    endif

  end subroutine track_contact_position

  !> This subroutine tracks down next contact position after a finite slide
  subroutine update_direction( nslave, contact, currpos )
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                  :: contact      !< contact info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node ), elem0(3, l_max_elem_node )
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)
    type(tContactState) :: cstate_tmp  !< Recorde of contact information

    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    odirec = contact%states(nslave)%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem(1:3,j)=currpos(3*iSS-2:3*iSS)
    enddo

    cstate_tmp%state = contact%states(nslave)%state
	cstate_tmp%surface = sid0
    call project_Point2Element( coord,etype,nn,elem,contact%master(sid0)%reflen,cstate_tmp, &
      isin,contact%cparam%DISTCLR_NOCHECK,contact%states(nslave)%lpos,contact%cparam%CLR_SAME_ELEM )

    if( isin ) then
      iSS = isInsideElement( etype, cstate_tmp%lpos, contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) &
        call cal_node_normal( cstate_tmp%surface, iSS, contact%master, currpos, &
        cstate_tmp%lpos, cstate_tmp%direction(:) )
    endif

    contact%states(nslave)%direction = cstate_tmp%direction

  end subroutine 


  !>\brief This subroutine update contact force in case that contacting element
  !> is changed
  subroutine reset_contact_force( contact, currpos, lslave, omaster, opos, odirec, B )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: currpos(:)     !< current coordinate of each nodes
    integer, intent(in)               :: lslave         !< slave node
    integer, intent(in)               :: omaster        !< former master element
    real(kind=kreal), intent(in)      :: opos(2)        !< former contact pos
    real(kind=kreal), intent(in)      :: odirec(3)      !< former contact direction
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, j, iSS
    real(kind=kreal)    :: nrlforce, fcoeff, tangent(3,2)
    real(kind=kreal)    :: elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(l_max_elem_node*3+3)
    integer(kind=kint)  :: i, idx0

    slave = contact%slave(lslave)
    fcoeff = contact%fcoeff
    ! clear contact force in former contacted element
    nrlforce = contact%states(lslave)%multiplier(1)
    B(3*slave-2:3*slave) = B(3*slave-2:3*slave)+nrlforce*odirec
    nn = size( contact%master(omaster)%nodes )
    etype = contact%master(omaster)%etype
    call getShapeFunc( etype, opos(:), shapefunc )
    do j=1,nn
      iSS = contact%master(omaster)%nodes(j)
      ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-nrlforce*shapefunc(j)*odirec
      idx0 = 3*(iSS-1)
      do i=1,3
        !$omp atomic
        B(idx0+i) = B(idx0+i)-nrlforce*shapefunc(j)*odirec(i)
      enddo
    enddo
    if( fcoeff/=0.d0 ) then
      do j=1,nn
        iSS = contact%master(omaster)%nodes(j)
        elemcrd(:,j) = currpos(3*iSS-2:3*iSS)
      enddo
      call DispIncreMatrix( opos(:), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      fric(1:2) = contact%states(lslave)%multiplier(2:3)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)+f3(1:3)
      do j=1,nn
        iSS = contact%master(omaster)%nodes(j)
        ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+f3(3*j+1:3*j+3)
        idx0 = 3*(iSS-1)
        do i=1,3
          !$omp atomic
          B(idx0+i) = B(idx0+i)+f3(3*j+i)
        enddo
      enddo
    endif

    ! reset contact force in new contacting element
    master = contact%states(lslave)%surface
    nn = size( contact%master(master)%nodes )
    etype = contact%master(master)%etype
    call getShapeFunc( etype, contact%states(lslave)%lpos(1:2), shapefunc )
    B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce*contact%states(lslave)%direction
    do j=1,nn
      iSS = contact%master(master)%nodes(j)
      ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+nrlforce*        &
        !            shapefunc(j)*contact%states(lslave)%direction
      idx0 = 3*(iSS-1)
      do i=1,3
        !$omp atomic
        B(idx0+i) = B(idx0+i)+nrlforce*        &
          shapefunc(j)*contact%states(lslave)%direction(i)
      enddo
    enddo
    if( fcoeff/=0.d0 ) then
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemcrd(:,j) = currpos(3*iSS-2:3*iSS)
      enddo
      call DispIncreMatrix( contact%states(lslave)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      fric(1:2) = contact%states(lslave)%multiplier(2:3)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-f3(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        ! B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-f3(3*j+1:3*j+3)
        idx0 = 3*(iSS-1)
        do i=1,3
          !$omp atomic
          B(idx0+i) = B(idx0+i)-f3(3*j+i)
        enddo
      enddo
    endif

  end subroutine reset_contact_force

  !>\brief This subroutine updates contact condition as follows:
  !!-# Contact force from multiplier and disp increment
  !!-# Update nodal force residual
  subroutine calcu_contact_force0( contact, coord, disp, ddisp, fcoeff, mu,     &
      mut, B )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: fcoeff         !< frictional coeff
    real(kind=kreal), intent(in)      :: mu, mut        !< penalty
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: nrlforce, elemdisp(3,l_max_elem_node), tangent(3,2)
    real(kind=kreal)    :: dg(3), elemg(3), elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(3*l_max_elem_node+3), edisp(3*l_max_elem_node+3)

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      edisp(1:3) = ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemdisp(:,j) = ddisp(3*iSS-2:3*iSS)
        edisp(3*j+1:3*j+3) = ddisp(3*iSS-2:3*iSS)
        elemcrd(:,j) = coord(3*iSS-2:3*iSS)+disp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      elemg = 0.d0
      do j=1,nn
        elemg(:) = elemg(:)+shapefunc(j)*elemdisp(:,j)
      enddo
      dg(:) = ddisp(3*slave-2:3*slave) -  elemg(:)
      dgn = dot_product( contact%states(i)%direction, dg )
      nrlforce = contact%states(i)%multiplier(1)- mu*(contact%states(i)%wkdist-dgn)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce*contact%states(i)%direction

      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+nrlforce*        &
          shapefunc(j)*contact%states(i)%direction
      enddo

      if( fcoeff==0.d0 ) cycle
      ! tangent component
      call DispIncreMatrix( contact%states(i)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      dxi(1) = dot_product( dispmat(1,1:nn*3+3), edisp(1:nn*3+3) )
      dxi(2) = dot_product( dispmat(2,1:nn*3+3), edisp(1:nn*3+3) )
      dxy(:) = matmul( metric, dxi )
      fric(1:2) = contact%states(i)%multiplier(2:3) + mut*dxy(1:2)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)

      if(  contact%states(i)%state==CONTACTSLIP ) then
        dgn = dsqrt( f3(1)*f3(1)+f3(2)*f3(2)+f3(3)*f3(3) )
        f3(:) = f3(:)*fcoeff*contact%states(i)%multiplier(1)/dgn
      endif
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-f3(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-f3(3*j+1:3*j+3)
      enddo
    enddo
  end subroutine calcu_contact_force0

  !>\brief This subroutine updates contact condition as follows:
  !!-# Contact force from multiplier and disp increment
  !!-# Update nodal force residual
  subroutine calcu_tied_force0( contact, disp, ddisp, mu, B )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: mu             !< penalty
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: nrlforce(3)
    real(kind=kreal)    :: dg(3)
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: edisp(3*l_max_elem_node+3)

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      edisp(1:3) = disp(3*slave-2:3*slave)+ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        edisp(3*j+1:3*j+3) = disp(3*iSS-2:3*iSS)+ddisp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      dg(1:3) = edisp(1:3)
      do j=1,nn
        dg(1:3) = dg(1:3)-shapefunc(j)*edisp(3*j+1:3*j+3)
      enddo

      nrlforce(1:3) = contact%states(i)%multiplier(1:3)+mu*dg(1:3)

      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS) + shapefunc(j)*nrlforce(1:3)
      enddo
    enddo

  end subroutine 

  !> This subroutine update lagrangian multiplier and the
  !> distance between contacting nodes
  subroutine update_contact_multiplier( contact, coord, disp, ddisp, fcoeff, mu, mut,   &
      gnt, ctchanged )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: fcoeff         !< frictional coeff
    real(kind=kreal), intent(in)      :: mu, mut        !< penalty
    real(kind=kreal), intent(out)     :: gnt(2)         !< convergency information
    logical, intent(inout)            :: ctchanged      !< if contact state changes

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS, cnt
    real(kind=kreal)    :: elemdisp(3,l_max_elem_node), tangent(3,2)
    real(kind=kreal)    :: dg(3), elemg(3), elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: lgnt(2), fric(2), f3(3*l_max_elem_node+3), edisp(3*l_max_elem_node+3)

    cnt =0; lgnt(:)=0.d0
    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      edisp(1:3) = ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemdisp(:,j) = ddisp(3*iSS-2:3*iSS)
        edisp(3*j+1:3*j+3) = ddisp(3*iSS-2:3*iSS)
        elemcrd(:,j) = coord(3*iSS-2:3*iSS)+disp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      elemg = 0.d0
      do j=1,nn
        elemg(:) = elemg(:)+shapefunc(j)*elemdisp(:,j)
      enddo
      dg(:) = ddisp(3*slave-2:3*slave) -  elemg(:)
      dgn = dot_product( contact%states(i)%direction, dg )
      contact%states(i)%wkdist = contact%states(i)%wkdist-dgn
      contact%states(i)%multiplier(1) = contact%states(i)%multiplier(1)- mu*contact%states(i)%wkdist
      contact%states(i)%distance = contact%states(i)%distance - dgn
      cnt = cnt+1
      lgnt(1) = lgnt(1)- contact%states(i)%wkdist

      if( fcoeff==0.d0 ) cycle
      ! tangent component
      call DispIncreMatrix( contact%states(i)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )
      dxi(1) = dot_product( dispmat(1,1:nn*3+3), edisp(1:nn*3+3) )
      dxi(2) = dot_product( dispmat(2,1:nn*3+3), edisp(1:nn*3+3) )
      dxy(:) = matmul( metric, dxi )
      fric(1:2) = contact%states(i)%multiplier(2:3) + mut*dxy(1:2)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      dgn = dsqrt( f3(1)*f3(1)+f3(2)*f3(2)+f3(3)*f3(3) )
      if( contact%states(i)%multiplier(1)>0.d0 ) then
        if(  dgn > fcoeff*contact%states(i)%multiplier(1) ) then
          if( contact%states(i)%state==CONTACTSTICK ) then
            ctchanged= .true.
            print *, "Node", slave, "to slip state",dgn, fcoeff*contact%states(i)%multiplier(1)
          endif
          contact%states(i)%state = CONTACTSLIP
          fric(:) = fric(:)*fcoeff*contact%states(i)%multiplier(1)/dgn
        else
          if( contact%states(i)%state==CONTACTSLIP ) then
            ctchanged= .true.
            print *, "Node", slave, "to stick state",dgn, fcoeff*contact%states(i)%multiplier(1)
          endif
          contact%states(i)%state = CONTACTSTICK
        endif
      endif
      contact%states(i)%multiplier(2:3) = fric(:)
      dxy(:) = matmul( dg, tangent )
      lgnt(2) = lgnt(2)+dsqrt( dxy(1)*dxy(1)+dxy(2)*dxy(2) )
    enddo
    if(cnt>0) lgnt(:) = lgnt(:)/cnt
    gnt = gnt + lgnt
  end subroutine update_contact_multiplier

  !> This subroutine update lagrangian multiplier and the
  !> distance between contacting nodes
  subroutine update_tied_multiplier( contact, disp, ddisp, mu, ctchanged )
    type( tContact ), intent(inout)   :: contact        !< contact info
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current step
    real(kind=kreal), intent(in)      :: ddisp(:)       !< disp till current substep
    real(kind=kreal), intent(in)      :: mu             !< penalty
    logical, intent(inout)            :: ctchanged      !< if contact state changes

    integer(kind=kint)  :: slave, etype, master
    integer(kind=kint)  :: nn, i, j, iSS, cnt
    real(kind=kreal)    :: dg(3), dgmax
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: edisp(3*l_max_elem_node+3)

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      edisp(1:3) = disp(3*slave-2:3*slave)+ddisp(3*slave-2:3*slave)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        edisp(3*j+1:3*j+3) = disp(3*iSS-2:3*iSS)+ddisp(3*iSS-2:3*iSS)
      enddo
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      dg(1:3) = edisp(1:3)
      do j=1,nn
        dg(1:3) = dg(1:3)-shapefunc(j)*edisp(3*j+1:3*j+3)
      enddo

      contact%states(i)%multiplier(1:3) = contact%states(i)%multiplier(1:3) + mu*dg(1:3)

      ! check if tied constraint converged
      dgmax = 0.d0
      do j=1,(nn+1)*3
        dgmax = dgmax + dabs(edisp(j))
      enddo
      dgmax = dgmax/dble((nn+1)*3)
      do j=1,3
        if( dabs(dg(j))/dmax1(1.d0,dgmax) > 1.d-3 ) ctchanged = .true.
      enddo

    enddo
  end subroutine 

  !>\brief This subroutine assemble contact force into contacing nodes
  subroutine ass_contact_force( contact, coord, disp, B )
    type( tContact ), intent(in)      :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current now
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, iSS
    real(kind=kreal)    :: fcoeff, nrlforce, tangent(3,2)
    real(kind=kreal)    :: elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(l_max_elem_node*3+3)
    fcoeff = contact%fcoeff
    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      call getShapeFunc( etype, contact%states(i)%lpos(1:2), shapefunc )

      ! normal component
      nrlforce = contact%states(i)%multiplier(1)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-nrlforce*contact%states(i)%direction

      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)+nrlforce*        &
          shapefunc(j)*contact%states(i)%direction
      enddo

      if( fcoeff==0.d0 ) cycle
      ! tangent component
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        elemcrd(:,j) = coord(3*iSS-2:3*iSS)+disp(3*iSS-2:3*iSS)
      enddo
      call DispIncreMatrix( contact%states(i)%lpos(1:2), etype, nn, elemcrd, tangent,   &
        metric, dispmat )

      fric(1:2) = contact%states(i)%multiplier(2:3)
      f3(:) = fric(1)*dispmat(1,:)+fric(2)*dispmat(2,:)
      B(3*slave-2:3*slave) = B(3*slave-2:3*slave)-f3(1:3)
      do j=1,nn
        iSS = contact%master(master)%nodes(j)
        B(3*iSS-2:3*iSS) = B(3*iSS-2:3*iSS)-f3(3*j+1:3*j+3)
      enddo
    enddo

  end subroutine ass_contact_force

  !>\brief This subroutine setup contact output nodal vectors
  subroutine set_contact_state_vector( contact, dt, relvel_vec, state_vec )
      type( tContact ), intent(in)      :: contact        !< contact info
      real(kind=kreal), intent(in)      :: dt
      real(kind=kreal), intent(inout)   :: relvel_vec(:)       !< mesh coordinate
      real(kind=kreal), intent(inout)   :: state_vec(:)        !< disp till current now

      integer(kind=kint)  :: i, slave

      do i= 1, size(contact%slave)
        slave = contact%slave(i)
        if( state_vec(slave) < 0.1d0 .or. contact%states(i)%state > 0 ) &
          &  state_vec(slave) = dble(contact%states(i)%state)

        if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
        if( dt < 1.d-16 ) cycle ! too small delta t
        relvel_vec(3*slave-2:3*slave) = contact%states(i)%reldisp(1:3)/dt
      enddo

  end subroutine set_contact_state_vector

  subroutine update_contact_TangentForce( contact )
    type( tContact ), intent(inout)   :: contact        !< contact info

    integer(kind=kint)  :: i

    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) then
        contact%states(i)%tangentForce(1:3) = 0.d0
        contact%states(i)%tangentForce_trial(1:3) = 0.d0
        contact%states(i)%tangentForce_final(1:3) = 0.d0
      else
        contact%states(i)%tangentForce(1:3) = contact%states(i)%tangentForce_final(1:3)
      end if
      contact%states(i)%tangentForce1(1:3) = contact%states(i)%tangentForce(1:3)
    enddo
  end subroutine update_contact_TangentForce

  !> This subroutine tracks down next contact position after a finite slide
  subroutine track_contact_position_exp( nslave, contact, currpos, currdisp, infoCTChange, nodeID, elemID )
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                 :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)  :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                    :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                    :: currdisp(:)  !< current displacement of each nodes
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node ), elem0(3, l_max_elem_node )
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)

    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    odirec = contact%states(nslave)%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2Element( coord,etype,nn,elem,contact%master(sid0)%reflen,contact%states(nslave), &
      isin,contact%cparam%DISTCLR_NOCHECK,contact%states(nslave)%lpos(1:2),contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        etype = contact%master(sid)%etype
        nn = getNumberOfNodes( etype )
        do j=1,nn
          iSS = contact%master(sid)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call project_Point2Element( coord,etype,nn,elem,contact%master(sid)%reflen,contact%states(nslave), &
          isin,contact%cparam%DISTCLR_NOCHECK,localclr=contact%cparam%CLEARANCE )
        if( isin ) then
          contact%states(nslave)%surface = sid
          exit
        endif
      enddo
    endif

    if( .not. isin ) then   ! such case is considered to rarely or never occur
      write(*,*) 'Warning: contact moved beyond neighbor elements'
      ! get master candidates from bucketDB
      bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
      nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
      if (nCand > 0) then
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)
        do idm= 1, nCand
          sid = indexCand(idm)
          if( sid==sid0 ) cycle
          if( any(sid==contact%master(sid0)%neighbor(:)) ) cycle
          if (.not. is_in_surface_box( contact%master(sid), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          etype = contact%master(sid)%etype
          nn = size( contact%master(sid)%nodes )
          do j=1,nn
            iSS = contact%master(sid)%nodes(j)
            elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          enddo
          call project_Point2Element( coord,etype,nn,elem,contact%master(sid)%reflen,contact%states(nslave), &
               isin,contact%cparam%DISTCLR_NOCHECK,localclr=contact%cparam%CLEARANCE )
          if( isin ) then
            contact%states(nslave)%surface = sid
            exit
          endif
        enddo
        deallocate(indexCand)
      endif
    endif

    if( isin ) then
      if( contact%states(nslave)%surface==sid0 ) then
        if(any(dabs(contact%states(nslave)%lpos(1:2)-opos(:)) >= contact%cparam%CLR_DIFFLPOS))  then
          !$omp atomic
          infoCTChange%contact2difflpos = infoCTChange%contact2difflpos + 1
        endif
      else
        write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3)') "Node",nodeID(slave)," move to contact with", &
          elemID(contact%master(sid)%eid), " with distance ",      &
          contact%states(nslave)%distance," at ",contact%states(nslave)%lpos(1:2)
        !$omp atomic
        infoCTChange%contact2neighbor = infoCTChange%contact2neighbor + 1
      endif
      iSS = isInsideElement( etype, contact%states(nslave)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) then
        call cal_node_normal( contact%states(nslave)%surface, iSS, contact%master, currpos, &
          contact%states(nslave)%lpos(1:2), contact%states(nslave)%direction(:) )
      endif
    else if( .not. isin ) then
      write(*,'(A,i10,A)') "Node",nodeID(slave)," move out of contact"
      contact%states(nslave)%state = CONTACTFREE
      contact%states(nslave)%multiplier(:) = 0.d0
    endif

  end subroutine track_contact_position_exp

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_contact_state_exp( contact, currpos, currdisp, infoCTChange, &
      nodeID, elemID, is_init, active )
    type( tContact ), intent(inout)                 :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)  :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                    :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                    :: currdisp(:)  !< current displacement of each nodes
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out
    logical, intent(in)                             :: is_init      !< whether initial scan or not
    logical, intent(out)                            :: active       !< if any in contact

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node )
    real(kind=kreal)    :: nlforce
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_cand

    if( is_init ) then
      distclr = contact%cparam%DISTCLR_INIT
    else
      distclr = contact%cparam%DISTCLR_FREE
      if( contact%algtype == CONTACTTIED ) then
        active = .false.
        do i= 1, size(contact%slave)
          if( contact%states(i)%state==CONTACTSTICK ) then
            active = .true.
            exit
          endif
        enddo
        return
      endif
    endif

    allocate(contact_surf(size(nodeID)))
    allocate(states_prev(size(contact%slave)))
    contact_surf(:) = size(elemID)+1
    do i = 1, size(contact%slave)
      states_prev(i) = contact%states(i)%state
    enddo

    call update_surface_box_info( contact%master, currpos )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )

    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,elem,is_cand,idm,etype,isin, &
      !$omp&         bktID,nCand,indexCand) &
      !$omp& firstprivate(nMasterMax) &
      !$omp& shared(contact,infoCTChange,currpos,currdisp,nodeID,elemID,distclr,contact_surf) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      slave = contact%slave(i)
      if( contact%states(i)%state==CONTACTSTICK .or. contact%states(i)%state==CONTACTSLIP ) then
        call track_contact_position_exp( i, contact, currpos, currdisp, infoCTChange, nodeID, elemID )
        if( contact%states(i)%state /= CONTACTFREE ) then
          contact_surf(contact%slave(i)) = -contact%states(i)%surface
        endif
      else if( contact%states(i)%state==CONTACTFREE ) then
        coord(:) = currpos(3*slave-2:3*slave)

        ! get master candidates from bucketDB
        bktID = bucketDB_getBucketID(contact%master_bktDB, coord)
        nCand = bucketDB_getNumCand(contact%master_bktDB, bktID)
        if (nCand == 0) cycle
        allocate(indexCand(nCand))
        call bucketDB_getCand(contact%master_bktDB, bktID, nCand, indexCand)

        nMasterMax = nCand
        allocate(indexMaster(nMasterMax))
        nMaster = 0

        ! narrow down candidates
        do idm= 1, nCand
          id = indexCand(idm)
          if (.not. is_in_surface_box( contact%master(id), coord(1:3), contact%cparam%BOX_EXP_RATE )) cycle
          nMaster = nMaster + 1
          indexMaster(nMaster) = id
        enddo
        deallocate(indexCand)

        if(nMaster == 0) then
          deallocate(indexMaster)
          cycle
        endif

        do idm = 1,nMaster
          id = indexMaster(idm)
          etype = contact%master(id)%etype
          nn = size( contact%master(id)%nodes )
          do j=1,nn
            iSS = contact%master(id)%nodes(j)
            elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          enddo
          call project_Point2Element( coord,etype,nn,elem,contact%master(id)%reflen,contact%states(i), &
            isin,distclr,localclr=contact%cparam%CLEARANCE )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          iSS = isInsideElement( etype, contact%states(i)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
          if( iSS>0 ) &
            call cal_node_normal( id, iSS, contact%master, currpos, contact%states(i)%lpos(1:2), &
              contact%states(i)%direction(:) )
          contact_surf(contact%slave(i)) = id
          write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3)') "Node",nodeID(slave)," contact with element", &
            elemID(contact%master(id)%eid),       &
            " with distance ", contact%states(i)%distance," at ",contact%states(i)%lpos(1:2), &
            " along direction ", contact%states(i)%direction
          exit
        enddo
        deallocate(indexMaster)
      endif
    enddo
    !$omp end parallel do

    call hecmw_contact_comm_allreduce_i(contact%comm, contact_surf, HECMW_MIN)
    nactive = 0
    do i = 1, size(contact%slave)
      if (contact%states(i)%state /= CONTACTFREE) then                    ! any slave in contact
        if (abs(contact_surf(contact%slave(i))) /= contact%states(i)%surface) then ! that is in contact with other surface
          contact%states(i)%state = CONTACTFREE                           ! should be freed
          write(*,'(A,i10,A,i6,A,i6,A)') "Node",nodeID(contact%slave(i)), &
            " in rank",hecmw_comm_get_rank()," freed due to duplication"
        else
          nactive = nactive + 1
        endif
      endif
      if (states_prev(i) == CONTACTFREE .and. contact%states(i)%state /= CONTACTFREE) then
        infoCTChange%free2contact = infoCTChange%free2contact + 1
      elseif (states_prev(i) /= CONTACTFREE .and. contact%states(i)%state == CONTACTFREE) then
        infoCTChange%contact2free = infoCTChange%contact2free + 1
      endif
    enddo
    active = (nactive > 0)
    deallocate(contact_surf)
    deallocate(states_prev)
  end subroutine scan_contact_state_exp

end module mContactDef
