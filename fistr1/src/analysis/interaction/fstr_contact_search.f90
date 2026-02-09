!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact search and state scanning at search level (all pairs in one tContact)
!>
module m_fstr_contact_search
  use m_fstr
  use hecmw
  use elementInfo
  use mContactDef
  use m_fstr_contact_geom
  use m_fstr_contact_element
  use m_fstr_contact_interference
  implicit none

contains

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

  !> This subroutine tracks down next contact position after a finite slide
  subroutine track_contact_position( flag_ctAlgo, nslave, contact, currpos, currdisp, infoCTChange, nodeID, elemID, B )
    character(len=9), intent(in)                    :: flag_ctAlgo  !< contact analysis algorithm flag
    integer, intent(in)                             :: nslave       !< slave node
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)    !< current displacement of each nodes
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out
    real(kind=kreal), intent(inout)                  :: B(:)         !< nodal force residual

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node), elem0(3, l_max_elem_node)
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
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
      contact%states(nslave)%lpos(1:2), contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        call project_Point2SurfElement( coord, contact%master(sid), currpos, &
          contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
          localclr=contact%cparam%CLEARANCE )
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
          call project_Point2SurfElement( coord, contact%master(sid), currpos, &
            contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
            localclr=contact%cparam%CLEARANCE )
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
      if( flag_ctAlgo=='SLagrange' ) then
        ! Setup elem array for update_TangentForce
        etype = contact%master(contact%states(nslave)%surface)%etype
        nn = size(contact%master(contact%states(nslave)%surface)%nodes)
        do j=1,nn
          iSS = contact%master(contact%states(nslave)%surface)%nodes(j)
          elem(1:3,j)=currpos(3*iSS-2:3*iSS)
        enddo
        call update_TangentForce(etype,nn,elem0,elem,contact%states(nslave))
      endif
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
    real(kind=kreal)    :: coord(3), elem0(3, l_max_elem_node )
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
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
      contact%states(nslave)%lpos(1:2), contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        call project_Point2SurfElement( coord, contact%master(sid), currpos, &
          contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
          localclr=contact%cparam%CLEARANCE )
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
          call project_Point2SurfElement( coord, contact%master(sid), currpos, &
            contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
            localclr=contact%cparam%CLEARANCE )
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
  subroutine scan_contact_state( flag_ctAlgo, contact, currpos, currdisp, ndforce, infoCTChange, &
    nodeID, elemID, is_init, active, B )
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
    real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: i, iSS, nactive
    real(kind=kreal)    :: coord(3)
    real(kind=kreal)    :: nlforce, slforce(3)
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_present_B
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
    !$omp& private(i,slave,slforce,id,nlforce,coord,indexMaster,nMaster,iSS,idm,etype,isin, &
    !$omp&         bktID,nCand,indexCand) &
    !$omp& firstprivate(nMasterMax,is_present_B) &
    !$omp& shared(contact,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,nodeID,elemID,Bp,distclr,contact_surf,is_init) &
    !$omp& reduction(.or.:active) &
    !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      if(contact%if_type /= 0) call set_shrink_factor(contact%ctime, contact%states(i), contact%if_etime, contact%if_type)
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
          call track_contact_position( flag_ctAlgo, i, contact, currpos, currdisp, infoCTChange, nodeID, elemID, Bp )
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
          call project_Point2SurfElement( coord, contact%master(id), currpos, &
            contact%states(i), isin, distclr, localclr=contact%cparam%CLEARANCE )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          etype = contact%master(id)%etype
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
    integer(kind=kint)  :: i, iSS, nactive
    real(kind=kreal)    :: coord(3)
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand

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
    !$omp& private(i,slave,id,coord,indexMaster,nMaster,iSS,idm,etype,isin, &
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
          call project_Point2SurfElement( coord, contact%master(id), currpos, &
            contact%states(i), isin, distclr, localclr=contact%cparam%CLEARANCE )
          if( .not. isin ) cycle
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          etype = contact%master(id)%etype
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

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  subroutine scan_embed_state( flag_ctAlgo, embed, currpos, currdisp, ndforce, infoCTChange, &
    nodeID, elemID, is_init, active, B )
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
    real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nn, i, j, iSS, nactive
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node )
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_present_B
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
    !$omp& private(i,slave,id,coord,indexMaster,nMaster,nn,j,iSS,elem,idm,etype,isin, &
    !$omp&         bktID,nCand,indexCand) &
    !$omp& firstprivate(nMasterMax,is_present_B) &
    !$omp& shared(embed,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,nodeID,elemID,Bp,distclr,contact_surf) &
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

  !> Scanning contact state
  subroutine remove_duplication_tiedcontact( cstep, hecMESH, fstrSOLID, infoCTChange )
    integer(kind=kint), intent(in)         :: cstep      !< current step number
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
    type(fstr_info_contactChange), intent(inout):: infoCTChange   !<

    integer(kind=kint) :: i, j, grpid, slave
    integer(kind=kint) :: k, id, iSS
    integer(kind=kint) :: ig0, ig, iS0, iE0
    integer(kind=kint), allocatable :: states(:)

    allocate(states(hecMESH%n_node))
    states(:) = CONTACTFREE

    ! if a boundary condition is given, the slave
    do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      iS0= hecMESH%node_group%grp_index(ig-1) + 1
      iE0= hecMESH%node_group%grp_index(ig  )
      do k= iS0, iE0
        iSS = hecMESH%node_group%grp_item(k)
        !states(iSS) = CONTACTSTICK
      enddo
    enddo

    do i=1,fstrSOLID%n_contacts
      if( fstrSOLID%contacts(i)%algtype /= CONTACTTIED ) cycle
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      do j=1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state==CONTACTFREE ) cycle   ! free
        slave = fstrSOLID%contacts(i)%slave(j)
        if( states(slave) == CONTACTFREE ) then
          states(slave) = fstrSOLID%contacts(i)%states(j)%state
          id = fstrSOLID%contacts(i)%states(j)%surface
          do k=1,size( fstrSOLID%contacts(i)%master(id)%nodes )
            iSS = fstrSOLID%contacts(i)%master(id)%nodes(k)
            states(iSS) = fstrSOLID%contacts(i)%states(j)%state
          enddo
        else !found duplicate tied contact slave node
          fstrSOLID%contacts(i)%states(j)%state = CONTACTFREE
          infoCTChange%free2contact = infoCTChange%free2contact - 1
          write(*,'(A,i10,A,i6,A,i6,A)') "Node",hecMESH%global_node_ID(slave), &
            " in rank",hecmw_comm_get_rank()," freed due to duplication"
        endif
      enddo
    enddo

  end subroutine

  subroutine scan_contact_state_ss( flag_ctAlgo, contact, currpos, currdisp, ndforce, infoCTChange, &
      nodeID, elemID, is_init, active, B )
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
    real(kind=kreal), optional, target               :: B(:)         !< nodal force residual

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: nnode_s, nnode_m, i, j, k, iSS, nactive, state, ii, n_intp
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node ), ncoord(2), normal(3)
    real(kind=kreal)    :: nlforce, nlforces(4)
    real(kind=kreal)    :: surf_node_pos(3,4), sfunc(4), N(128,4), weight(128)
    logical             :: isin, success
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
    endif
    !$omp parallel do &
      !$omp& default(none) &
      !$omp& private(i,slave,id,nlforce,nlforces,coord,ncoord,surf_node_pos,sfunc,indexMaster,nMaster,nnode_s,nnode_m, &
      !$omp&         normal,j,iSS,elem,is_cand,idm,etype,isin,bktID,nCand,indexCand) &
      !$omp& firstprivate(nMasterMax,is_present_B) &
      !$omp& shared(contact,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,nodeID,elemID,Bp,distclr,contact_surf,is_init) &
      !$omp& reduction(.or.:active) &
      !$omp& schedule(dynamic,1)

    ! loop element; if any node(nslave_index) in contact, call intp point projection
    do i = 1, size(contact%slave_surf)
      surf_node_pos(:,:) = 0.d0
      nnode_s = size(contact%slave_surf(i)%nodes)
      contact%slave_surf(i)%state = CONTACTINACTIVE !!!! dbg-later

      do j = 1, nnode_s
        slave = contact%slave_surf(i)%nslave_index(j)
        call switch_intp_contact(contact%slave_surf(i),j,contact%states(slave)%state)
      enddo
      if (contact%slave_surf(i)%state == CONTACTINACTIVE) cycle

      do j = 1, nnode_s
        slave = contact%slave_surf(i)%nodes(j)
        surf_node_pos(:,j) = currpos(3*slave-2:3*slave)
      enddo
      do j = 1, contact%slave_surf(i)%n_intp
        if(contact%slave_surf(i)%states(j)%state==CONTACTINACTIVE) cycle
        if(contact%if_type /= 0) &
         call set_shrink_factor(contact%ctime, contact%slave_surf(i)%states(j), contact%if_etime, contact%if_type)
        
        if( contact%slave_surf(i)%states(j)%state==CONTACTSTICK .or. contact%slave_surf(i)%states(j)%state==CONTACTSLIP ) then
          id = contact%slave_surf(i)%states(j)%surface
          call getIntPoint4ss( contact%slave_surf(i)%etype, j, ncoord, contact%slave_surf(i)%n_intp, sfunc )
          coord(:) = matmul( surf_node_pos(:,1:nnode_s), sfunc(1:nnode_s) )

          ! update direction of TIED contact
          if( contact%algtype == CONTACTTIED ) then
            call update_direction( i, contact, currpos ) !要修正
            if (.not.is_init) cycle
          endif

          ! if( contact%algtype /= CONTACTFSLID .or. (.not. is_present_B) ) then   ! small slide problem
          !   continue
          ! !   contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid) ! ??? -> contact_surf
          ! else
          call track_contact_position_ss( flag_ctAlgo, contact%slave_surf(i)%eid, coord, contact%slave_surf(i)%states(j),&
            contact, currpos, currdisp, infoCTChange, nodeID, elemID, Bp )

          if( contact%slave_surf(i)%states(j)%state /= CONTACTFREE ) then
            id = contact%slave_surf(i)%states(j)%surface
            ! contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
            ! Calculate and store normal direction in states structure
            call getIntPoint4ss( contact%slave_surf(i)%etype, j, ncoord, contact%slave_surf(i)%n_intp, sfunc )
            normal(:) = - SurfaceNormal( contact%slave_surf(i)%etype, nnode_s, ncoord, surf_node_pos )
            contact%slave_surf(i)%states(j)%direction = normal(:)/dsqrt( dot_product(normal, normal) )
          endif
          ! endif

        else if( contact%slave_surf(i)%states(j)%state==CONTACTACTIVE ) then
          if( contact%algtype == CONTACTTIED .and. .not. is_init ) cycle 
          call getIntPoint4ss( contact%slave_surf(i)%etype, j, ncoord, contact%slave_surf(i)%n_intp, sfunc )
          coord(:) = matmul( surf_node_pos(:,1:nnode_s), sfunc(1:nnode_s) )

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
            call project_Point2SurfElement( coord, contact%master(id), currpos, &
              contact%slave_surf(i)%states(j), isin, distclr, localclr=contact%cparam%CLEARANCE )
            if( .not. isin ) cycle
            contact%slave_surf(i)%states(j)%surface = id
            contact%slave_surf(i)%states(j)%multiplier(:) = 0.d0
            etype = contact%master(id)%etype

            ! Calculate and store normal direction for new contact
            normal(:) = - SurfaceNormal( contact%slave_surf(i)%etype, nnode_s, ncoord, surf_node_pos )
            contact%slave_surf(i)%states(j)%direction = normal(:)/dsqrt( dot_product(normal, normal) )
            iSS = isInsideElement( etype, contact%slave_surf(i)%states(j)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
            if( iSS>0 ) &
              call cal_node_normal( id, iSS, contact%master, currpos, contact%slave_surf(i)%states(j)%lpos(1:2), &
                contact%slave_surf(i)%states(j)%direction(:) )
            ! contact_surf(contact%slave(i)) = elemID(contact%master(id)%eid)
            contact%slave_surf(i)%msurf_list(j) = id
            write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3,A,i6)') &
              "Element",elemID(contact%slave_surf(i)%eid)," contact with element", elemID(contact%master(id)%eid),       &
              " with distance ", contact%slave_surf(i)%states(j)%distance," at ",contact%slave_surf(i)%states(j)%lpos(1:2), &
              " along direction ", contact%slave_surf(i)%states(j)%direction," rank=",hecmw_comm_get_rank()
            exit
          enddo
          deallocate(indexMaster)
        endif
      enddo
    enddo
    !$omp end parallel do

    if( contact%algtype == CONTACTTIED .and. .not. is_init ) then
      deallocate(contact_surf)
      return
    endif

    ! calc dual basis
    do i = 1, size(contact%slave_surf)
      if (contact%slave_surf(i)%state == CONTACTFREE) cycle
      nnode_s = size(contact%slave_surf(i)%nodes)
      n_intp = contact%slave_surf(i)%n_intp
      do j = 1, nnode_s
        slave = contact%slave_surf(i)%nodes(j)
        surf_node_pos(:,j) = currpos(3*slave-2:3*slave)
      enddo
      call getWeightAndN4ss(contact%slave_surf(i)%etype, nnode_s, n_intp, surf_node_pos, &
      & weight(1:n_intp), N(1:n_intp,1:nnode_s))
      call getDualShapfunc(contact%slave_surf(i), nnode_s, n_intp, surf_node_pos, weight, N, success)
      if(.not. success) contact%slave_surf(i)%state = CONTACTFREE
    enddo

  end subroutine scan_contact_state_ss

  subroutine track_contact_position_ss( flag_ctAlgo, eid, coord, state, contact, currpos, currdisp, infoCTChange,&
     nodeID, elemID, B )
    character(len=9), intent(in)                    :: flag_ctAlgo  !< contact analysis algorithm flag
    integer, intent(in)                             :: eid       !< slave elem
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( tContactState ), intent(inout)             :: state      !< slave intp state
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)    !< current displacement of each nodes
    integer(kind=kint), intent(in)                  :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                  :: elemID(:)    !< global elemental ID, just for print out
    real(kind=kreal), intent(inout)                  :: B(:)         !< nodal force residual
    real(kind=kreal), intent(inout)                  :: coord(:)         !< position of int point

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: elem(3, l_max_elem_node ), elem0(3, l_max_elem_node )
    logical            :: isin
    real(kind=kreal)    :: opos(2), odirec(3)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)

    sid = 0

    !> checking the contact element of last step
    sid0 = state%surface
    opos = state%lpos(1:2)
    odirec = state%direction
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo

    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      state, isin, contact%cparam%DISTCLR_NOCHECK, state%lpos(1:2), contact%cparam%CLR_SAME_ELEM )
    if( .not. isin ) then ! search neighbor master surf
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        etype = contact%master(sid)%etype
        call project_Point2SurfElement( coord, contact%master(sid), currpos, &
          state, isin, contact%cparam%DISTCLR_NOCHECK, &
          localclr=contact%cparam%CLEARANCE )
        if( isin ) then
          state%surface = sid
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
          call project_Point2SurfElement( coord, contact%master(sid), currpos, &
            state, isin, contact%cparam%DISTCLR_NOCHECK, &
            localclr=contact%cparam%CLEARANCE )
          if( isin ) then
            state%surface = sid
            exit
          endif
        enddo
        deallocate(indexCand)
      endif
    endif

    if( isin ) then
      if( state%surface==sid0 ) then
        if(any(dabs(state%lpos(1:2)-opos(:)) >= contact%cparam%CLR_DIFFLPOS))  then
          !$omp atomic
          infoCTChange%contact2difflpos = infoCTChange%contact2difflpos + 1
        endif
      else
        write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3)') "Element",elemID(eid)," move to contact with", &
          elemID(contact%master(sid)%eid), " with distance ",      &
          state%distance," at ",state%lpos(1:2)
        !$omp atomic
        infoCTChange%contact2neighbor = infoCTChange%contact2neighbor + 1
        ! if( flag_ctAlgo=='ALagrange' )  &
        !   call reset_contact_force( contact, currpos, nslave, sid0, opos, odirec, B )
      endif
      !ここはあとで修正
      write(*,*)'need to fix'
      if( flag_ctAlgo=='SLagrange' ) call update_TangentForce(etype,nn,elem0,elem,state)
      iSS = isInsideElement( etype, state%lpos(1:2), contact%cparam%CLR_CAL_NORM )
      if( iSS>0 ) &
        call cal_node_normal( state%surface, iSS, contact%master, currpos, &
          state%lpos(1:2), state%direction(:) )
    else if( .not. isin ) then
      write(*,'(A,i10,A)') "Node",nodeID(slave)," move out of contact"
      state%state = CONTACTFREE
      state%multiplier(:) = 0.d0
    endif

  end subroutine track_contact_position_ss

end module m_fstr_contact_search
