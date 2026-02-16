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
  use m_fstr_contact_smoothing
  implicit none

contains

  !> This subroutine tracks down next contact position after a finite slide
  !! When flag_ctAlgo is present, implicit-specific processing
  !! (update_TangentForce for SLagrange) is performed.
  subroutine track_contact_position( nslave, contact, currpos, currdisp, infoCTChange, nodeID, elemID, flag_ctAlgo )
    integer, intent(in)                              :: nslave       !< slave node
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)  !< current displacement of each nodes
    integer(kind=kint), intent(in)                   :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                   :: elemID(:)    !< global elemental ID, just for print out
    character(len=9), intent(in), optional           :: flag_ctAlgo  !< contact analysis algorithm flag (implicit only)

    integer(kind=kint) :: slave, sid0, sid, etype
    integer(kind=kint) :: nn, i, j, iSS
    real(kind=kreal)    :: coord(3), elem(3, l_max_elem_node), elem0(3, l_max_elem_node)
    logical            :: isin
    real(kind=kreal)    :: opos(2)
    integer(kind=kint) :: bktID, nCand, idm
    integer(kind=kint), allocatable :: indexCand(:)
    logical            :: is_implicit

    is_implicit = present(flag_ctAlgo)

    sid = 0

    slave = contact%slave(nslave)
    coord(:) = currpos(3*slave-2:3*slave)
    !> checking the contact element of last step
    sid0 = contact%states(nslave)%surface
    opos = contact%states(nslave)%lpos(1:2)
    etype = contact%master(sid0)%etype
    nn = getNumberOfNodes( etype )
    do j=1,nn
      iSS = contact%master(sid0)%nodes(j)
      elem0(1:3,j)=currpos(3*iSS-2:3*iSS)-currdisp(3*iSS-2:3*iSS)
    enddo
    call project_Point2SurfElement( coord, contact%master(sid0), currpos, &
      contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
      contact%states(nslave)%lpos(1:2), contact%cparam%CLR_SAME_ELEM, smoothing=contact%smoothing )
    if( .not. isin ) then
      do i=1, contact%master(sid0)%n_neighbor
        sid = contact%master(sid0)%neighbor(i)
        call project_Point2SurfElement( coord, contact%master(sid), currpos, &
          contact%states(nslave), isin, contact%cparam%DISTCLR_NOCHECK, &
          localclr=contact%cparam%CLEARANCE, smoothing=contact%smoothing )
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
            localclr=contact%cparam%CLEARANCE, smoothing=contact%smoothing )
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
      if( is_implicit .and. flag_ctAlgo=='SLagrange' ) then
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

  !> This subroutine update contact states, which include
  !!-# Free to contact or contact to free state changes
  !!-# Clear lagrangian multipliers when free to contact
  !! When flag_ctAlgo is present, implicit-specific processing is performed:
  !!   - TIED direction update with cycle (instead of early return)
  !!   - Tensile force check using ndforce
  !!   - SLagrange update_TangentForce in track_contact_position
  subroutine scan_contact_state( contact, currpos, currdisp, infoCTChange, &
    nodeID, elemID, is_init, active, hecMESH, flag_ctAlgo, ndforce )
    type( tContact ), intent(inout)                  :: contact      !< contact info
    type( fstr_info_contactChange ), intent(inout)   :: infoCTChange !< contact change info
    real(kind=kreal), intent(in)                     :: currpos(:)   !< current coordinate of each nodes
    real(kind=kreal), intent(in)                     :: currdisp(:)  !< current displacement of each nodes
    integer(kind=kint), intent(in)                   :: nodeID(:)    !< global nodal ID, just for print out
    integer(kind=kint), intent(in)                   :: elemID(:)    !< global elemental ID, just for print out
    logical, intent(in)                              :: is_init      !< whether initial scan or not
    logical, intent(out)                             :: active       !< if any in contact
    type( hecmwST_local_mesh ), intent(in), optional :: hecMESH      !< mesh for MPI communication
    character(len=9), intent(in), optional           :: flag_ctAlgo  !< contact analysis algorithm flag (implicit only)
    real(kind=kreal), intent(in), optional           :: ndforce(:)   !< nodal force (implicit only)

    real(kind=kreal)    :: distclr
    integer(kind=kint)  :: slave, id, etype
    integer(kind=kint)  :: i, iSS, nactive
    real(kind=kreal)    :: coord(3)
    real(kind=kreal)    :: nlforce
    logical             :: isin
    integer(kind=kint), allocatable :: contact_surf(:), states_prev(:)
    real(kind=kreal)    :: effective_near_dist, distclr_use
    !
    integer, pointer :: indexMaster(:),indexCand(:)
    integer   ::  nMaster,idm,nMasterMax,bktID,nCand
    logical :: is_implicit

    is_implicit = present(flag_ctAlgo)

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

    if( contact%smoothing == kcsNAGATA ) then
      call update_surface_normal( contact%master, currpos, hecMESH )
      call create_intermediate_points( contact%master, currpos, contact%pair_name )
    endif

    call update_surface_box_info( contact%master, currpos )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )

    ! Compute effective near distance for CONTACTNEAR detection
    effective_near_dist = contact%cparam%NEAR_DIST

    !$omp parallel do &
    !$omp& default(none) &
    !$omp& private(i,slave,id,nlforce,coord,indexMaster,nMaster,iSS,idm,etype,isin, &
    !$omp&         bktID,nCand,indexCand,distclr_use) &
    !$omp& firstprivate(nMasterMax,is_implicit,effective_near_dist) &
    !$omp& shared(contact,ndforce,flag_ctAlgo,infoCTChange,currpos,currdisp,nodeID,elemID,distclr,contact_surf,is_init) &
    !$omp& reduction(.or.:active) &
    !$omp& schedule(dynamic,1)
    do i= 1, size(contact%slave)
      if(contact%if_type /= 0) call set_shrink_factor(contact%ctime, contact%states(i), contact%if_etime, contact%if_type)
      slave = contact%slave(i)
      if( contact%states(i)%state==CONTACTSTICK .or. contact%states(i)%state==CONTACTSLIP ) then
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
        if( contact%algtype /= CONTACTFSLID ) then   ! small slide problem
          contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
        else
          if( is_implicit ) then
            call track_contact_position( i, contact, currpos, currdisp, infoCTChange, nodeID, elemID, &
              flag_ctAlgo=flag_ctAlgo )
          else
            call track_contact_position( i, contact, currpos, currdisp, infoCTChange, nodeID, elemID )
          endif
          if( contact%states(i)%state /= CONTACTFREE ) then
            id = contact%states(i)%surface
            contact_surf(contact%slave(i)) = -elemID(contact%master(id)%eid)
          endif
        endif

      else if( contact%states(i)%state==CONTACTNEAR ) then
        ! NEAR tracking: re-project on current master surface
        coord(:) = currpos(3*slave-2:3*slave)
        id = contact%states(i)%surface
        if (effective_near_dist > 0.0d0) then
          distclr_use = max(distclr, effective_near_dist / contact%master(id)%reflen)
        else
          distclr_use = distclr
        end if
        call project_Point2SurfElement( coord, contact%master(id), currpos, &
          contact%states(i), isin, distclr_use, &
          contact%states(i)%lpos(1:2), contact%cparam%CLR_SAME_ELEM, smoothing=contact%smoothing )
        if (isin) then
          if (contact%states(i)%distance <= distclr * contact%master(id)%reflen) then
            ! Within contact threshold -> upgrade to STICK
            contact%states(i)%state = CONTACTSTICK
            contact%states(i)%multiplier(:) = 0.d0
            write(*,'(A,i10,A,i10,A,i6)') "Node",nodeID(slave)," upgraded NEAR->STICK on element", &
              elemID(contact%master(id)%eid)," rank=",hecmw_comm_get_rank()
          else if (contact%states(i)%distance > effective_near_dist) then
            ! Beyond NEAR range -> free
            contact%states(i)%state = CONTACTFREE
            contact%states(i)%multiplier(:) = 0.d0
          end if
          ! else: stay NEAR
          if (.not. is_contact_free(contact%states(i)%state)) then
            etype = contact%master(id)%etype
            iSS = isInsideElement( etype, contact%states(i)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
            if( iSS>0 ) &
              call cal_node_normal( id, iSS, contact%master, currpos, &
              contact%states(i)%lpos(1:2), contact%states(i)%direction(:) )
            contact_surf(contact%slave(i)) = elemID(contact%master(id)%eid)
          end if
        else
          ! Lost projection -> free
          contact%states(i)%state = CONTACTFREE
          contact%states(i)%multiplier(:) = 0.d0
        end if

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
          ! Expand distclr for NEAR detection
          if (effective_near_dist > 0.0d0) then
            distclr_use = max(distclr, effective_near_dist / contact%master(id)%reflen)
          else
            distclr_use = distclr
          end if
          call project_Point2SurfElement( coord, contact%master(id), currpos, &
            contact%states(i), isin, distclr_use, localclr=contact%cparam%CLEARANCE, smoothing=contact%smoothing )
          if( .not. isin ) cycle
          ! Classify: STICK or NEAR
          if (effective_near_dist > 0.0d0 .and. &
              contact%states(i)%distance > distclr * contact%master(id)%reflen) then
            contact%states(i)%state = CONTACTNEAR
          end if
          contact%states(i)%surface = id
          contact%states(i)%multiplier(:) = 0.d0
          etype = contact%master(id)%etype
          iSS = isInsideElement( etype, contact%states(i)%lpos(1:2), contact%cparam%CLR_CAL_NORM )
          if( iSS>0 ) &
            call cal_node_normal( id, iSS, contact%master, currpos, contact%states(i)%lpos(1:2), &
            contact%states(i)%direction(:) )
          contact_surf(contact%slave(i)) = elemID(contact%master(id)%eid)
          if (contact%states(i)%state == CONTACTNEAR) then
            write(*,'(A,i10,A,i10,A,f7.3,A,i6)') "Node",nodeID(slave)," near element", &
              elemID(contact%master(id)%eid), &
              " with distance ", contact%states(i)%distance," rank=",hecmw_comm_get_rank()
          else
            write(*,'(A,i10,A,i10,A,f7.3,A,2f7.3,A,3f7.3,A,i6)') "Node",nodeID(slave)," contact with element", &
              elemID(contact%master(id)%eid),       &
              " with distance ", contact%states(i)%distance," at ",contact%states(i)%lpos(1:2), &
              " along direction ", contact%states(i)%direction," rank=",hecmw_comm_get_rank()
          end if
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
      if (.not. is_contact_free(contact%states(i)%state)) then             ! any slave in contact or near
        id = contact%states(i)%surface
        if (abs(contact_surf(contact%slave(i))) /= elemID(contact%master(id)%eid)) then ! that is in contact with other surface
          contact%states(i)%state = CONTACTFREE                           ! should be freed
          write(*,'(A,i10,A,i10,A,i6,A,i6,A)') "Node",nodeID(contact%slave(i))," contact with element", &
          &  elemID(contact%master(id)%eid), " in rank",hecmw_comm_get_rank()," freed due to duplication"
        else if (is_contact_active(contact%states(i)%state)) then
          nactive = nactive + 1
        endif
      endif
      if (is_contact_free(states_prev(i)) .and. .not. is_contact_free(contact%states(i)%state)) then
        infoCTChange%free2contact = infoCTChange%free2contact + 1
      elseif (.not. is_contact_free(states_prev(i)) .and. is_contact_free(contact%states(i)%state)) then
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
        if( is_contact_free(fstrSOLID%contacts(i)%states(j)%state) ) cycle   ! free
        slave = fstrSOLID%contacts(i)%slave(j)
        if( is_contact_free(states(slave)) ) then
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

end module m_fstr_contact_search

