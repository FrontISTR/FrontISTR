!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module manage the data structure for contact calculation
!!
!!  Contact calculation takes into act after calling the following three
!!  subrotuines provided in this module
!!-#      Reading contact definition with subroutine: fstr_ctrl_get_CONTACT
!!-#      Check its consistency with mesh definition: fstr_contact_check
!!-#      Initilizing the contact calculation       : fstr_contact_init
module mContactDef

  use hecmw
  use mSurfElement
  use m_contact_lib
  use m_fstr_contact_comm

  implicit none

  include 'fstr_ctrl_util_f.inc'


  !> Structure to includes all info needed by contact calculation
  type tContact
    ! following contact definition
    character(len=HECMW_NAME_LEN) :: name                    !< name
    integer                       :: ctype                   !< 1:node-surface 2: surface-surface
    integer                       :: group                   !< group number
    character(len=HECMW_NAME_LEN) :: pair_name               !< name of contact pair
    integer                       :: surf_id1, surf_id2      !< slave surafce, master surface
    type(tSurfElement), pointer   :: master(:)=>null()       !< master surface (element )
    integer, pointer              :: slave(:)=>null()        !< slave surface (node)
    real(kind=kreal)              :: fcoeff                  !< coeeficient of friction
    real(kind=kreal)              :: tPenalty                !< tangential penalty

    ! following algorithm
    ! -1: not initialized
    ! 1: TIED-Just rigidly fixed the two surfaces
    ! 2: GLUED-Distance between the two surfaces to zero and glue them
    ! 3: SSLID-Small sliding contact( no position but with contact state change)
    ! 4: FSLID-Finite sliding contact (both changes in contact state and poisition possible)
    integer                       :: algtype                 !< algorithm flag

    logical                       :: mpced                   !< if turns into mpc condition
    logical                       :: symmetric               !< if symmetrizalized in friction calculation

    ! following contact state
    type(tContactState), pointer  :: states(:)=>null()       !< contact states of each slave nodes

    type(fstrST_contact_comm)     :: comm                    !< contact communication table
  end type tContact

  type fstr_info_contactChange
    integer(kind=kint) :: contact2free           !< counter: contact to free state change
    integer(kind=kint) :: contact2neighbor       !< counter: contact to neighbor state change
    integer(kind=kint) :: contact2diffLpos       !< counter: contact to different local position state change
    integer(kind=kint) :: free2contact           !< counter: free to contact state change
    integer(kind=kint) :: contactNode_previous   !< previous number of nodes in contact
    integer(kind=kint) :: contactNode_current    !< current number of nodes in contact
  end type fstr_info_contactChange

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
    if( associated(contact%slave) ) then
      do i=1,size(contact%slave)
        write(file, *) contact%slave(i)
      enddo
    endif
    write(file,*) "----master---"
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
    call fstr_contact_comm_finalize(contact%comm)
  end subroutine

  !>  Check the consistency with given mesh of contact defintiion
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
        isfind = .true.
      endif
    enddo
    if( .not. isfind ) return;
    if( contact%fcoeff<=0.d0 ) contact%fcoeff=0.d0
    if( contact%ctype/=1 .and. contact%ctype/=2 ) return
    if( contact%group<=0 ) return

    fstr_contact_check = .true.
  end function

  !>  Initializer of tContactState
  logical function fstr_contact_init( contact, hecMESH,myrank )
    type(tContact), intent(inout)     :: contact  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    integer(kint),intent(in),optional :: myrank

    integer  :: i, j, is, ie, cgrp, nsurf, nslave, ic, ic_type, iss, nn, ii
    integer  :: count,nodeID

    fstr_contact_init = .false.
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

    ! slave surface
    !    if( contact%ctype==1 ) then
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
    call find_surface_neighbor( contact%master )

    ! initialize contact communication table
    call fstr_contact_comm_init( contact%comm, hecMESH, 1, nslave, contact%slave )

    contact%symmetric = .true.
    fstr_contact_init = .true.
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
  !!-# Free to contact ot contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_contact_state( flag_ctAlgo, contact, currpos, currdisp, ndforce, infoCTChange, &
      nodeID, elemID, active, mu, B )
    character(len=9), intent(in)                    :: flag_ctAlgo  !< contact analysis algorithm flag
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)    !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)    !< current displacement of each nodes
    real(kind=kreal), intent(in)                     :: ndforce(:)    !< nodal force
    integer(kind=kint), intent(in)                  :: nodeID(:)     !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)     !< global elemental ID, just for print out
    logical, intent(out)                            :: active        !< if any in contact
    real(kind=kreal), intent(in)                     :: mu            !< penalty
    real(kind=kreal), optional                       :: B(:)          !< nodal force residual

    real(kind=kreal)    :: clearance
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node )
    real(kind=kreal)    :: nlforce, slforce(3)
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    !#ifdef CONTACT_FILTER
    !integer   ::  indexMaster(100),nMaster=0,minID(3),maxID(3),idm
    integer, pointer :: indexMaster(:),itmp(:)
    integer   ::  nMaster=0,minID(3),maxID(3),idm,nMasterMax=100
    real(kreal) :: width = 4.0D0,x0(3)
    !#endif

    clearance = 1.d-6
    if( contact%algtype<=2 ) return

    allocate(contact_surf(size(nodeID)))
    allocate(states_prev(size(contact%slave)))
    contact_surf(:) = size(elemID)+1
    do i = 1, size(contact%slave)
      states_prev(i) = contact%states(i)%state
    enddo

    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,slforce,id,nlforce,coord,indexMaster,nMaster,nn,j,iSS,elem,x0,minID,maxID,itmp,idm,etype,isin) &
      !$omp& firstprivate(nMasterMax) &
      !$omp& shared(contact,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,mu,nodeID,elemID,B,width,clearance,contact_surf) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      slave = contact%slave(i)
      if( contact%states(i)%state==CONTACTSTICK .or. contact%states(i)%state==CONTACTSLIP ) then
        slforce(1:3)=ndforce(3*slave-2:3*slave)
        id = contact%states(i)%surface
        nlforce = contact%states(i)%multiplier(1)

        if( nlforce < -1.0d-8 ) then
          contact%states(i)%state = CONTACTFREE
          contact%states(i)%multiplier(:) = 0.d0
          write(*,'(A,i10,A,i10,A,e12.3)') "Node",nodeID(slave)," free from contact with element", &
            elemID(contact%master(id)%eid), " with tensile force ", nlforce
          cycle
        endif
        if( contact%algtype /= CONTACTFSLID .or. (.not. present(B)) ) then   ! small slide problem
          contact_surf(contact%slave(i)) = -id
        else
          call track_contact_position( flag_ctAlgo, i, contact, currpos, currdisp, mu, infoCTChange, nodeID, elemID, B )
          if( contact%states(i)%state /= CONTACTFREE ) then
            contact_surf(contact%slave(i)) = -contact%states(i)%surface
          endif
        endif

      else if( contact%states(i)%state==CONTACTFREE ) then
        coord(:) = currpos(3*slave-2:3*slave)
        allocate(indexMaster(nMasterMax))
        !#ifdef CONTACT_FILTER
        nMaster = 0
        !#endif
        do id= 1, size(contact%master)

          !#ifdef CONTACT_FILTER
          nn = size( contact%master(id)%nodes )
          do j=1,nn
            iSS = contact%master(id)%nodes(j)
            elem(1:3,j) = currpos(3*iSS-2:3*iSS)
          enddo
          x0(:) = coord(:)
          call getMinMaxBoxIDPassedByMultiPoint(elem(1:3,1:nn),x0,width,minID,maxID)
          !if(minID(1) <= 1.and.maxID(1) >= 1.and. &
            !   minID(2) <= 1.and.maxID(2) >= 1.and. &
            !   minID(3) <= 1.and.maxID(3) >= 1) then
          if(any(minID > 1).or.any(maxID < 0)) cycle
          nMaster = nMaster + 1
          if(nMaster > size(indexMaster)) then
            !stop 'Error: Too many master faces are possibly in contact!'
            itmp => indexMaster
            allocate(indexMaster(nMasterMax*2))
            indexMaster(1:nMasterMax) = itmp(1:nMasterMax)
            deallocate(itmp)
            nMasterMax = nMasterMax*2
            write(*,*) 'Info: increased nMasterMax to ', nMasterMax
          endif
          indexMaster(nMaster) = id
          !endif
        enddo

        if(nMaster == 0) then
          deallocate(indexMaster)
          cycle
        endif
        !          do id= 1, size(contact%master)
        do idm = 1,nMaster
          id = indexMaster(idm)
          !#endif
          etype = contact%master(id)%etype
          nn = size( contact%master(id)%nodes )
          do j=1,nn
            iSS = contact%master(id)%nodes(j)
            elem(1:3,j)=currpos(3*iSS-2:3*iSS)
          enddo
          call project_Point2Element(coord,etype,nn,elem,contact%states(i), isin, clearance )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          iSS = isInsideElement( etype, contact%states(i)%lpos, clearance )
          if( iSS>0 ) call cal_node_normal( id, iSS, contact%master, currpos, contact%states(i)%direction(:) )
          contact_surf(contact%slave(i)) = id
          write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3)') "Node",nodeID(slave)," contact with element", &
            elemID(contact%master(id)%eid),       &
            " with distance ", contact%states(i)%distance," at ",contact%states(i)%lpos(:), &
            " along direction ", contact%states(i)%direction
          exit
        enddo
        deallocate(indexMaster)
      endif
    enddo
    !$omp end parallel do

    call fstr_contact_comm_allreduce_i(contact%comm, contact_surf, HECMW_MIN)
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
  end subroutine scan_contact_state

  !> Calculate averaged nodal normal
  subroutine cal_node_normal( csurf,cnode, surf, currpos, normal )
    use elementInfo, only:getVertexCoord, SurfaceNormal
    integer, intent(in)            :: csurf       !< current surface element
    integer, intent(in)            :: cnode       !< current node position
    type(tSurfElement), intent(in) :: surf(:)     !< surface elements
    real(kind=kreal), intent(in)   :: currpos(:)  !< current coordinate of each nodes
    real(kind=kreal), intent(out)  :: normal(3)   !< averaged node nomral
    integer :: i, j, cnt, nd1, gn, etype, iSS, nn,cgn
    real(kind=kreal) :: cnpos(2), elem(3, l_max_elem_node )

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
        if( iSS==gn ) cgn=iSS
      enddo
      if( cgn>0 ) then
        call getVertexCoord( etype, cgn, cnpos )
        normal = normal+SurfaceNormal( etype, nn, cnpos, elem )
        cnt = cnt+1
      endif
    enddo
    normal = normal/cnt                                        !!-???
    normal = normal/ dsqrt( dot_product( normal, normal ) )
  end subroutine

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

    real(kind=kreal)    :: distclr, clearance
    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node ), elem0(3, l_max_elem_node )
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)

    distclr = 1.0d0          !1.d-1
    clearance = 1.d-6
    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos
    odirec = contact%states(nslave)%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem(1:3,j)=currpos(3*iSS-2:3*iSS)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2Element(coord,etype,nn,elem,contact%states(nslave), isin, distclr, &
      contact%states(nslave)%lpos, clearance )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        etype = contact%master(sid)%etype
        nn = getNumberOfNodes( etype )
        do j=1,nn
          iSS = contact%master(sid)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call project_Point2Element(coord,etype,nn,elem,contact%states(nslave), isin, distclr )
        if( isin ) then
          contact%states(nslave)%surface = sid
          exit
        endif
      enddo
    endif

    if( .not. isin ) then   ! such case is considered to rarely or never occur
      do sid= 1, size(contact%master)
        if( sid==sid0 ) cycle
        if( any(sid==contact%master(sid0)%neighbor(:)) ) cycle
        etype = contact%master(sid)%etype
        nn = size( contact%master(sid)%nodes )
        do j=1,nn
          iSS = contact%master(sid)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call project_Point2Element(coord,etype,nn,elem,contact%states(nslave), isin, distclr )
        if( isin ) then
          contact%states(nslave)%surface = sid
          exit
        endif
      enddo
    endif

    if( isin ) then
      if( contact%states(nslave)%surface==sid0 ) then
        if(any(dabs(contact%states(nslave)%lpos(:)-opos(:)) >= 1.0d-3))  then
          !$omp atomic
          infoCTChange%contact2difflpos = infoCTChange%contact2difflpos + 1
        endif
      else
        write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3)') "Node",nodeID(slave)," move to contact with", &
          elemID(contact%master(sid)%eid), " with distance ",      &
          contact%states(nslave)%distance," at ",contact%states(nslave)%lpos(:)
        !$omp atomic
        infoCTChange%contact2neighbor = infoCTChange%contact2neighbor + 1
        if( flag_ctAlgo=='ALagrange' )  &
          call reset_contact_force( contact, currpos, nslave, sid0, opos, odirec, B )
      endif
      if( flag_ctAlgo=='SLagrange' ) call update_TangentForce(etype,nn,elem0,elem,contact%states(nslave))
      iSS = isInsideElement( etype, contact%states(nslave)%lpos, clearance )
      if( iSS>0 ) &
        call cal_node_normal( contact%states(nslave)%surface, iSS, contact%master, currpos, contact%states(nslave)%direction(:) )
    else if( .not. isin ) then
      write(*,'(A,i10,A)') "Node",nodeID(slave)," move out of contact"
      contact%states(nslave)%state = CONTACTFREE
      contact%states(nslave)%multiplier(:) = 0.d0
    endif

  end subroutine track_contact_position

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
    real(kind=kreal)    :: dum, dxi(2), shapefunc(l_max_surface_node)
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
    call getShapeFunc( etype, contact%states(lslave)%lpos(:), shapefunc )
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
      call DispIncreMatrix( contact%states(lslave)%lpos, etype, nn, elemcrd, tangent,   &
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

  !>\brief This subroutine update contact condition as follows:
  !!-# Contact force from multipler and disp increment
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
    real(kind=kreal)    :: dum, dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
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
      call getShapeFunc( etype, contact%states(i)%lpos(:), shapefunc )

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
      call DispIncreMatrix( contact%states(i)%lpos, etype, nn, elemcrd, tangent,   &
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
    logical, intent(inout)            :: ctchanged       !< if contact state changes

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
      call getShapeFunc( etype, contact%states(i)%lpos(:), shapefunc )

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
      call DispIncreMatrix( contact%states(i)%lpos, etype, nn, elemcrd, tangent,   &
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

  !>\brief This subroutine assemble contact force into contacing nodes
  subroutine ass_contact_force( contact, coord, disp, B )
    type( tContact ), intent(in)      :: contact        !< contact info
    real(kind=kreal), intent(in)      :: coord(:)       !< mesh coordinate
    real(kind=kreal), intent(in)      :: disp(:)        !< disp till current now
    real(kind=kreal), intent(inout)   :: B(:)           !< nodal force residual

    integer(kind=kint)  :: slave,  etype, master
    integer(kind=kint)  :: nn, i, j, k, iSS
    real(kind=kreal)    :: fcoeff, nrlforce, tangent(3,2)
    real(kind=kreal)    :: dg(3), elemg(3), elemcrd(3, l_max_elem_node )
    real(kind=kreal)    :: dum, dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
    real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
    real(kind=kreal)    :: fric(2), f3(l_max_elem_node*3+3)
    fcoeff = contact%fcoeff
    do i= 1, size(contact%slave)
      if( contact%states(i)%state==CONTACTFREE ) cycle   ! not in contact
      slave = contact%slave(i)
      master = contact%states(i)%surface

      nn = size( contact%master(master)%nodes )
      etype = contact%master(master)%etype
      call getShapeFunc( etype, contact%states(i)%lpos(:), shapefunc )

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
      call DispIncreMatrix( contact%states(i)%lpos, etype, nn, elemcrd, tangent,   &
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

  subroutine getMinMaxBoxIDPassedByMultiPoint(x,x0,width,minID,maxID)
    real(kreal),intent(in)      ::  x(:,:)  ! Multi-Points' coordinate in space
    real(kreal),intent(in)      ::  x0(3)
    real(kreal),intent(in)      ::  width
    integer(kint),intent(out)   ::  minID(3)
    integer(kint),intent(out)   ::  maxID(3)
    !
    integer(kint) ::  i,j,boxID(3)
    !
    do i=1,size(x,2)
      !boxID(:) = (x(:,i)-x0(:))/width + 1
      boxID(:) = ceiling((x(:,i)-x0(:))/width)
      if(i == 1) then
        minID(:) = boxID(:)
        maxID(:) = boxID(:)
      endif
      do j=1,3
        if(boxID(j) < minID(j)) minID(j) = boxID(j)
        if(boxID(j) > maxID(j)) maxID(j) = boxID(j)
      enddo
    enddo
  end subroutine getMinMaxBoxIDPassedByMultiPoint

  !>\brief This subroutine setup contact output nodal vectors
  subroutine set_contact_state_vector( contact, dt, relvel_vec, state_vec )
      type( tContact ), intent(in)      :: contact        !< contact info
      real(kind=kreal), intent(in)      :: dt
      real(kind=kreal), intent(inout)   :: relvel_vec(:)       !< mesh coordinate
      real(kind=kreal), intent(inout)   :: state_vec(:)        !< disp till current now

      integer(kind=kint)  :: slave,  etype, master
      integer(kind=kint)  :: nn, i, j, k, iSS
      real(kind=kreal)    :: fcoeff, nrlforce, tangent(3,2)
      real(kind=kreal)    :: dg(3), elemg(3), elemcrd(3, l_max_elem_node )
      real(kind=kreal)    :: dum, dgn, dxi(2), dxy(2), shapefunc(l_max_surface_node)
      real(kind=kreal)    :: metric(2,2), dispmat(2,l_max_elem_node*3+3)
      real(kind=kreal)    :: fric(2), f3(l_max_elem_node*3+3)

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

end module mContactDef
