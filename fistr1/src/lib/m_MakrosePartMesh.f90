!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Jun Yin (ASTOM)                                !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines to create distributed mesh information
!   for parallel contact simulation

module m_MakrosePartMesh
  use m_Makrose
  implicit none
  integer(kint),private,save,allocatable   ::  helpNode2(:,:)
  integer(kint),private,save,allocatable   ::  helpElmt2(:,:)

!  interface Mak_GetLocalMesh_NodeBase1
!    module procedure Mak_GetLocalMesh_NodeBase1_1
!    module procedure Mak_GetLocalMesh_NodeBase1_all
!  end interface
  
contains

!< Node-based partition
subroutine Mak_GetNeighborIndex_NodeBase(mak,nparts,part,partID,nbIndex)
! Get neighboring partition ID of the current one(partID),
  type(MAKROSE_STRUCT),intent(in)   ::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),intent(in)          ::  partID
  integer(kint),pointer             ::  nbIndex(:)  
!
  integer(kint)   ::  i,j,k,n,istat,nodeID
  integer(kint),allocatable   ::  help(:),idx(:)
!
    if(.not.allocated(helpNode2)) then
      allocate(helpNode2(nparts,nparts),stat=istat)
      helpNode2(:,:) = 0
      allocate(help(nparts),stat=istat)
      allocate(idx(nparts),stat=istat)
      do i=1,mak%ne
        idx(:) = 0; n = 0; help(:) = 0
        do j=mak%eptr(i),mak%eptr(i+1)-1
          nodeID = mak%eind(j)
          if(help(part(nodeID)) /= 0) cycle
          n = n + 1
          idx(n) = part(nodeID)
          help(part(nodeID)) = 1
        enddo
        if(n == 1) cycle
        do j=1,n
          do k=1,n
            if(k == j) cycle
            if(helpNode2(idx(k),idx(j)) /= 0) cycle
            helpNode2(idx(k),idx(j)) = 1
          enddo
        enddo
      enddo
      deallocate(help,idx,stat=istat)
    endif
    n = 0
    do i=1,nparts
      if(i == partID) cycle
      if(helpNode2(i,partID) /= 0) n = n + 1
    enddo
    if(associated(nbIndex)) deallocate(nbIndex,stat=istat)
    allocate(nbIndex(n),stat=istat)
    n = 0
    do i=1,nparts
      if(i == partID) cycle
      if(helpNode2(i,partID) /= 0) then
        n = n + 1
        nbIndex(n) = i
      endif
    enddo    
end subroutine Mak_GetNeighborIndex_NodeBase

subroutine Mak_GetLocalMeshMapping_NodeBase(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
  type(MAKROSE_STRUCT),intent(in)   ::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),pointer             ::  ni_loc(:)
  integer(kint),pointer             ::  nx_loc(:)
  integer(kint),pointer             ::  ei_loc(:)
  integer(kint),pointer             ::  ex_loc(:)
  integer(kint),pointer             ::  indexNodeG2L(:,:)
  integer(kint),pointer             ::  indexElmtG2L(:,:)
!
  integer(kint)   ::  i,j,k,n,istat,nodeID,nodeID1,nodeID2,partID
  integer(kint),allocatable   ::  helpNode(:,:),helpElmt(:,:),help(:)
!
!   Initiate pointers
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)
!
    allocate(ni_loc(nparts),stat=istat);  ni_loc(:) = 0
    allocate(nx_loc(nparts),stat=istat);  nx_loc(:) = 0
    allocate(ei_loc(nparts),stat=istat);  ei_loc(:) = 0
    allocate(ex_loc(nparts),stat=istat);  ex_loc(:) = 0
    allocate(indexNodeG2L(mak%nn,nparts),stat=istat);  indexNodeG2L(:,:) = 0
    allocate(indexElmtG2L(mak%ne,nparts),stat=istat);  indexElmtG2L(:,:) = 0
!   Internal nodes
    do i=1,mak%nn
      ni_loc(part(i)) = ni_loc(part(i)) + 1
      indexNodeG2L(i,part(i)) = ni_loc(part(i))
    enddo
!   External nodes
    if(allocated(helpNode)) deallocate(helpNode,stat=istat)
    allocate(helpNode(mak%nn,nparts),stat=istat); helpNode(:,:) = 0
    if(allocated(helpElmt)) deallocate(helpElmt,stat=istat)
    allocate(helpElmt(mak%ne,nparts),stat=istat); helpElmt(:,:) = 0
    allocate(help(nparts),stat=istat)
!   Internal elements first
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) then
        ei_loc(partID) = ei_loc(partID) + 1
        indexElmtG2L(i,partID) = ei_loc(partID)
        cycle
      endif
    enddo
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) cycle
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID1 = mak%eind(j)
        do k=mak%eptr(i),mak%eptr(i+1)-1
          nodeID2 = mak%eind(k)
          if(nodeID1 == nodeID2.or.part(nodeID1) == part(nodeID2)) cycle
          if(helpNode(nodeID2,part(nodeID1)) == 0) then
            helpNode(nodeID2,part(nodeID1)) = part(nodeID2)
            nx_loc(part(nodeID1)) = nx_loc(part(nodeID1)) + 1
            indexNodeG2L(nodeID2,part(nodeID1)) = ni_loc(part(nodeID1)) + nx_loc(part(nodeID1))
          endif
          if(helpNode(nodeID1,part(nodeID2)) == 0) then
            helpNode(nodeID1,part(nodeID2)) = part(nodeID1)
            nx_loc(part(nodeID2)) = nx_loc(part(nodeID2)) + 1
            indexNodeG2L(nodeID1,part(nodeID2)) = ni_loc(part(nodeID2)) + nx_loc(part(nodeID2))
          endif
          if(helpElmt(i,part(nodeID1)) == 0) then
            helpElmt(i,part(nodeID1)) = part(nodeID2)
            ex_loc(part(nodeID1)) = ex_loc(part(nodeID1)) + 1
            indexElmtG2L(i,part(nodeID1)) = ei_loc(part(nodeID1)) + ex_loc(part(nodeID1))
          endif
          if(helpElmt(i,part(nodeID2)) == 0) then
            helpElmt(i,part(nodeID2)) = part(nodeID1)
            ex_loc(part(nodeID2)) = ex_loc(part(nodeID2)) + 1
            indexElmtG2L(i,part(nodeID2)) = ei_loc(part(nodeID2)) + ex_loc(part(nodeID2))
          endif
        enddo
      enddo      
    enddo
    deallocate(helpNode,helpElmt,help,stat=istat) 
end subroutine Mak_GetLocalMeshMapping_NodeBase

subroutine Mak_GetLocalMeshMapping_NodeBase_InOrder(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
  type(MAKROSE_STRUCT),intent(inout)::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),pointer             ::  ni_loc(:)
  integer(kint),pointer             ::  nx_loc(:)
  integer(kint),pointer             ::  ei_loc(:)
  integer(kint),pointer             ::  ex_loc(:)
  integer(kint),pointer             ::  indexNodeG2L(:,:)
  integer(kint),pointer             ::  indexElmtG2L(:,:)
!
  integer(kint)   ::  i,j,k,n,istat,nodeID,nodeID1,nodeID2,partID,elmtID,pid
  integer(kint),allocatable   ::  helpNode(:,:),helpElmt(:,:),help(:)
!
!   Initiate pointers
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)
!
    allocate(ni_loc(nparts),stat=istat);  ni_loc(:) = 0
    allocate(nx_loc(nparts),stat=istat);  nx_loc(:) = 0
    allocate(ei_loc(nparts),stat=istat);  ei_loc(:) = 0
    allocate(ex_loc(nparts),stat=istat);  ex_loc(:) = 0
    allocate(indexNodeG2L(mak%nn,nparts),stat=istat);  indexNodeG2L(:,:) = 0
    allocate(indexElmtG2L(mak%ne,nparts),stat=istat);  indexElmtG2L(:,:) = 0
!   Internal nodes
    do i=1,mak%nn
      ni_loc(part(i)) = ni_loc(part(i)) + 1
      indexNodeG2L(i,part(i)) = ni_loc(part(i))
    enddo
!   External nodes
    if(allocated(helpNode)) deallocate(helpNode,stat=istat)
    allocate(helpNode(mak%nn,nparts),stat=istat); helpNode(:,:) = 0
    if(allocated(helpElmt)) deallocate(helpElmt,stat=istat)
    allocate(helpElmt(mak%ne,nparts),stat=istat); helpElmt(:,:) = 0
    allocate(help(nparts),stat=istat)
!   Internal elements first
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) then
        ei_loc(partID) = ei_loc(partID) + 1
        indexElmtG2L(i,partID) = ei_loc(partID)
        cycle
      endif
    enddo
!   External nodes
    if(.not.associated(mak%nptr)) then
      call Mak_GetElementOnNode(mak)
    endif
    do i=1,mak%nn
      partID = part(i)
      help(:) = 0
      do j=mak%nptr(i),mak%nptr(i+1)-1
        elmtID = mak%nind(j)
        if(any(indexElmtG2L(elmtID,1:nparts) /= 0)) cycle
        do k=mak%eptr(elmtID),mak%eptr(elmtID+1)-1
          nodeID = mak%eind(k)
          if(nodeID == i) cycle
          if(part(nodeID) == partID) cycle
          if(help(part(nodeID)) /= 0) cycle
          help(part(nodeID)) = 1
          nx_loc(part(nodeID)) = nx_loc(part(nodeID)) + 1
          indexNodeG2L(i,part(nodeID)) = ni_loc(part(nodeID)) + nx_loc(part(nodeID))
        enddo
      enddo
    enddo

!   External elements     
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
      enddo
      if(n == 1) cycle
!     Find min nodal partition #
      pid = nparts
      help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) == 1) cycle
        help(part(nodeID)) = 1
        pid = min(pid,part(nodeID))
      enddo
      ei_loc(pid) = ei_loc(pid) + 1
      indexElmtG2L(i,pid) = ei_loc(pid) + ex_loc(pid)
      do j=1,nparts
        if(j /= pid.and.help(j) == 1) then
          ex_loc(j) = ex_loc(j) + 1
          if(indexElmtG2L(i,j) /= 0) then
            print *,'again',pid,j,ei_loc(j),ex_loc(j)
          endif
          indexElmtG2L(i,j) = ei_loc(j) + ex_loc(j)
        endif
      enddo
      
!      do j=mak%eptr(i),mak%eptr(i+1)-1
!        nodeID1 = mak%eind(j)
!        do k=mak%eptr(i),mak%eptr(i+1)-1
!          nodeID2 = mak%eind(k)
!          if(nodeID1 == nodeID2.or.part(nodeID1) == part(nodeID2)) cycle
!!          if(helpNode(nodeID2,part(nodeID1)) == 0) then
!!            helpNode(nodeID2,part(nodeID1)) = part(nodeID2)
!!            nx_loc(part(nodeID1)) = nx_loc(part(nodeID1)) + 1
!!            indexNodeG2L(nodeID2,part(nodeID1)) = ni_loc(part(nodeID1)) + nx_loc(part(nodeID1))
!!          endif
!!          if(helpNode(nodeID1,part(nodeID2)) == 0) then
!!            helpNode(nodeID1,part(nodeID2)) = part(nodeID1)
!!            nx_loc(part(nodeID2)) = nx_loc(part(nodeID2)) + 1
!!            indexNodeG2L(nodeID1,part(nodeID2)) = ni_loc(part(nodeID2)) + nx_loc(part(nodeID2))
!!          endif
!          if(helpElmt(i,part(nodeID1)) == 0) then
!            helpElmt(i,part(nodeID1)) = part(nodeID2)
!            ex_loc(part(nodeID1)) = ex_loc(part(nodeID1)) + 1
!            indexElmtG2L(i,part(nodeID1)) = ei_loc(part(nodeID1)) + ex_loc(part(nodeID1))
!          endif
!          if(helpElmt(i,part(nodeID2)) == 0) then
!            helpElmt(i,part(nodeID2)) = part(nodeID1)
!            ex_loc(part(nodeID2)) = ex_loc(part(nodeID2)) + 1
!            indexElmtG2L(i,part(nodeID2)) = ei_loc(part(nodeID2)) + ex_loc(part(nodeID2))
!          endif
!        enddo
!      enddo      
    enddo
    deallocate(helpNode,helpElmt,help,stat=istat) 
end subroutine Mak_GetLocalMeshMapping_NodeBase_InOrder

subroutine Mak_GetLocalMeshMapping_NodeBase_InOrder_FrontISTR(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
  type(MAKROSE_STRUCT),intent(inout)::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),pointer             ::  ni_loc(:)
  integer(kint),pointer             ::  nx_loc(:)
  integer(kint),pointer             ::  ei_loc(:)
  integer(kint),pointer             ::  ex_loc(:)
  integer(kint),pointer             ::  indexNodeG2L(:,:)
  integer(kint),pointer             ::  indexElmtG2L(:,:)
!
  integer(kint)   ::  i,j,k,n,istat,nodeID,nodeID1,nodeID2,partID,elmtID,pid
  integer(kint),allocatable   ::  helpNode(:,:),helpElmt(:,:),help(:)
!
!   Initiate pointers
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
!    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
!    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)
    indexNodeG2L => null()
    indexElmtG2L => null()
!
    allocate(ni_loc(nparts),stat=istat);  ni_loc(:) = 0
    allocate(nx_loc(nparts),stat=istat);  nx_loc(:) = 0
    allocate(ei_loc(nparts),stat=istat);  ei_loc(:) = 0
    allocate(ex_loc(nparts),stat=istat);  ex_loc(:) = 0
    allocate(indexNodeG2L(mak%nn,nparts),stat=istat);  indexNodeG2L(:,:) = 0
    allocate(indexElmtG2L(mak%ne,nparts),stat=istat);  indexElmtG2L(:,:) = 0
!   Internal nodes
    do i=1,mak%nn
      ni_loc(part(i)) = ni_loc(part(i)) + 1
      indexNodeG2L(i,part(i)) = ni_loc(part(i))
    enddo
!   External nodes
    if(allocated(helpNode)) deallocate(helpNode,stat=istat)
    allocate(helpNode(mak%nn,nparts),stat=istat); helpNode(:,:) = 0
    if(allocated(helpElmt)) deallocate(helpElmt,stat=istat)
    allocate(helpElmt(mak%ne,nparts),stat=istat); helpElmt(:,:) = 0
    allocate(help(nparts),stat=istat)
!   Internal elements first
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) then
        ei_loc(partID) = ei_loc(partID) + 1
        indexElmtG2L(i,partID) = ei_loc(partID)
        cycle
      endif
    enddo
!   External nodes
    if(.not.associated(mak%nptr)) then
      call Mak_GetElementOnNode(mak)
    endif
    do i=1,mak%nn
      partID = part(i)
      help(:) = 0
      do j=mak%nptr(i),mak%nptr(i+1)-1
        elmtID = mak%nind(j)
        if(any(indexElmtG2L(elmtID,1:nparts) /= 0)) cycle
        do k=mak%eptr(elmtID),mak%eptr(elmtID+1)-1
          nodeID = mak%eind(k)
          if(nodeID == i) cycle
          if(part(nodeID) == partID) cycle
          if(help(part(nodeID)) /= 0) cycle
          help(part(nodeID)) = 1
          nx_loc(part(nodeID)) = nx_loc(part(nodeID)) + 1
          indexNodeG2L(i,part(nodeID)) = ni_loc(part(nodeID)) + nx_loc(part(nodeID))
        enddo
      enddo
    enddo

!   Internal again, then external elements
    ei_loc(:) = 0; ex_loc(:) = 0
    indexElmtG2L(:,:) = 0
    do i=1,mak%ne
      n = 0; help(:) = 0
      pid = nparts
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        pid = min(pid,part(nodeID))
      enddo
      if(n == 1) then
        ei_loc(pid) = ei_loc(pid) + 1
        indexElmtG2L(i,pid) = ei_loc(pid) + ex_loc(pid)
        mak%egrp(i) = pid + 1
      else
        ! internal
        ei_loc(pid) = ei_loc(pid) + 1
        indexElmtG2L(i,pid) = ei_loc(pid) + ex_loc(pid)
        mak%egrp(i) = pid + 1
        ! external
        do j=1,nparts
          if(help(j) == 1.and.j /= pid) then
            ex_loc(j) = ex_loc(j) + 1
            indexElmtG2L(i,j) = ei_loc(j) + ex_loc(j)
          endif
        enddo
      endif
    enddo
    deallocate(helpNode,helpElmt,help,stat=istat) 
end subroutine Mak_GetLocalMeshMapping_NodeBase_InOrder_FrontISTR

subroutine Mak_GetLocalMesh_NodeBase1(mak,nparts,part,partID,mak_loc,indexNodeG2L_me,indexElmtG2L_me)
  type(MAKROSE_STRUCT),intent(in)   ::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),intent(in)          ::  partID  ! PartID of local mesh
  type(MAKROSE_STRUCT),intent(out)  ::  mak_loc ! Local mesh
  integer(kint),pointer,optional    ::  indexNodeG2L_me(:)
  integer(kint),pointer,optional    ::  indexElmtG2L_me(:)
!
  integer(kint)   ::  i,istat
  integer(kint),pointer     ::  ni_loc(:)=>null(),nx_loc(:)=>null(),ei_loc(:)=>null(),ex_loc(:)=>null()
  integer(kint),pointer     ::  indexNodeG2L(:,:)=>null(),indexElmtG2L(:,:)=>null()
!
!   Get local mesh mapping info
    call Mak_GetLocalMeshMapping_NodeBase(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
!   Initiate local mesh data
!#ifdef VIA_ELMT_TYPE
    call Mak_ResetIndexElmtG2L_Type(mak,indexElmtG2L(:,partID))
!#endif    
    call Mak_SetLocalMesh(mak,ni_loc(partID),nx_loc(partID),ei_loc(partID),ex_loc(partID),  &
                          indexNodeG2L(:,partID),indexElmtG2L(:,partID),partID,mak_loc)
!   Get index global to local of current partID
    if(present(indexNodeG2L_me)) then
      allocate(indexNodeG2L_me(size(indexNodeG2L,1)),stat=istat)
      indexNodeG2L_me(:) = indexNodeG2L(:,partID)
    endif
    if(present(indexElmtG2L_me)) then
      allocate(indexElmtG2L_me(size(indexElmtG2L,1)),stat=istat)
      indexElmtG2L_me(:) = indexElmtG2L(:,partID)
    endif
!   Release temp memory
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)    
end subroutine Mak_GetLocalMesh_NodeBase1

subroutine Mak_GetLocalMesh_NodeBase1_all(mak,nparts,part,partID,mak_loc,indexNodeG2L_me,indexElmtG2L_me)
  type(MAKROSE_STRUCT),intent(inout)::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),intent(in)          ::  partID  ! PartID of local mesh
  type(MAKROSE_STRUCT),intent(out)  ::  mak_loc ! Local mesh
  integer(kint),pointer,optional    ::  indexNodeG2L_me(:,:)
  integer(kint),pointer,optional    ::  indexElmtG2L_me(:,:)
!
  integer(kint)   ::  i,istat
  integer(kint),pointer     ::  ni_loc(:)=>null(),nx_loc(:)=>null(),ei_loc(:)=>null(),ex_loc(:)=>null()
  integer(kint),pointer     ::  indexNodeG2L(:,:)=>null(),indexElmtG2L(:,:)=>null()
!
!   Get local mesh mapping info
    call Mak_GetLocalMeshMapping_NodeBase_InOrder(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
!   Initiate local mesh data
!#ifdef VIA_ELMT_TYPE
    call Mak_ResetIndexElmtG2L_Type(mak,indexElmtG2L(:,partID))
!#endif
    print *,partID,'local',ei_loc(partID),ex_loc(partID)
    call Mak_SetLocalMesh(mak,ni_loc(partID),nx_loc(partID),ei_loc(partID),ex_loc(partID),  &
                          indexNodeG2L(:,partID),indexElmtG2L(:,partID),partID,mak_loc)
!   Get index global to local of current partID
    if(present(indexNodeG2L_me)) then
      if(associated(indexNodeG2L_me)) deallocate(indexNodeG2L_me,stat=istat)
      allocate(indexNodeG2L_me(size(indexNodeG2L,1),size(indexNodeG2L,2)),stat=istat)
      indexNodeG2L_me(:,:) = indexNodeG2L(:,:)
    endif
    if(present(indexElmtG2L_me)) then
      if(associated(indexElmtG2L_me)) deallocate(indexElmtG2L_me,stat=istat)
      allocate(indexElmtG2L_me(size(indexElmtG2L,1),size(indexElmtG2L,2)),stat=istat)
      indexElmtG2L_me(:,:) = indexElmtG2L(:,:)
    endif
!   Release temp memory
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)    
end subroutine Mak_GetLocalMesh_NodeBase1_all

subroutine Mak_GetLocalMesh_NodeBase1_FrontISTR(mak,nparts,part,partID,mak_loc,indexNodeG2L_me,indexElmtG2L_me)
  type(MAKROSE_STRUCT),intent(inout)::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),intent(in)          ::  partID  ! PartID of local mesh
  type(MAKROSE_STRUCT),intent(out)  ::  mak_loc ! Local mesh
  integer(kint),pointer,optional    ::  indexNodeG2L_me(:,:)
  integer(kint),pointer,optional    ::  indexElmtG2L_me(:,:)
!
  integer(kint)   ::  i,istat
  integer(kint),pointer     ::  ni_loc(:)=>null(),nx_loc(:)=>null(),ei_loc(:)=>null(),ex_loc(:)=>null()
  integer(kint),pointer     ::  indexNodeG2L(:,:)=>null(),indexElmtG2L(:,:)=>null()
!
!   Get local mesh mapping info
    call Mak_GetLocalMeshMapping_NodeBase_InOrder_FrontISTR(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
!   Initiate local mesh data
!#ifdef VIA_ELMT_TYPE
    call Mak_ResetIndexElmtG2L_Type(mak,indexElmtG2L(:,partID))
!#endif
!    print *,partID,'local',ei_loc(partID),ex_loc(partID)
    call Mak_SetLocalMesh(mak,ni_loc(partID),nx_loc(partID),ei_loc(partID),ex_loc(partID),  &
                          indexNodeG2L(:,partID),indexElmtG2L(:,partID),partID,mak_loc)
!   Get index global to local of current partID
    if(present(indexNodeG2L_me)) then
      if(associated(indexNodeG2L_me)) deallocate(indexNodeG2L_me,stat=istat)
!      allocate(indexNodeG2L_me(size(indexNodeG2L,1),size(indexNodeG2L,2)),stat=istat)
!      indexNodeG2L_me(:,:) = indexNodeG2L(:,:)
      indexNodeG2L_me => indexNodeG2L
    endif
    if(present(indexElmtG2L_me)) then
      if(associated(indexElmtG2L_me)) deallocate(indexElmtG2L_me,stat=istat)
!      allocate(indexElmtG2L_me(size(indexElmtG2L,1),size(indexElmtG2L,2)),stat=istat)
!      indexElmtG2L_me(:,:) = indexElmtG2L(:,:)
      indexElmtG2L_me => indexElmtG2L
    endif
!   Release temp memory
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
!    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
!    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)    
end subroutine Mak_GetLocalMesh_NodeBase1_FrontISTR

subroutine Mak_GetLocalMesh_NodeBaseAll(mak,nparts,part,mak_loc)
  type(MAKROSE_STRUCT),intent(in)   ::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  type(MAKROSE_STRUCT),intent(out)  ::  mak_loc(:) ! Local mesh
!
  integer(kint)   ::  i,istat
  integer(kint),pointer     ::  ni_loc(:)=>null(),nx_loc(:)=>null(),ei_loc(:)=>null(),ex_loc(:)=>null()
  integer(kint),pointer     ::  indexNodeG2L(:,:)=>null(),indexElmtG2L(:,:)=>null()
!
!   Get local mesh mapping info
    call Mak_GetLocalMeshMapping_NodeBase(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
!   Initiate local mesh data
    do i=1,nparts
!#ifdef VIA_ELMT_TYPE
      call Mak_ResetIndexElmtG2L_Type(mak,indexElmtG2L(:,i))
!#endif
      call Mak_SetLocalMesh(mak,ni_loc(i),nx_loc(i),ei_loc(i),ex_loc(i),    &
                            indexNodeG2L(:,i),indexElmtG2L(:,i),i,mak_loc(i))
    enddo
!   Release temp memory
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)
end subroutine Mak_GetLocalMesh_NodeBaseAll

subroutine Mak_ResetIndexElmtG2L_Type(mak,indexElmtG2L)
  type(MAKROSE_STRUCT),intent(in)   ::  mak     ! Total global mesh
  integer(kint),intent(inout)       ::  indexElmtG2L(:) ! Reorder the local mesh elementID by element types
!
  integer(kint)   ::  i,j,n,istat
  integer(kint),allocatable   ::  indexType(:),counter(:),indexElmtG2L_tmp(:),ptr(:)
!
    allocate(indexType(20),stat=istat)
    indexType(:) = 0
    allocate(counter(20),stat=istat)
    counter(:) = 0
    allocate(indexElmtG2L_tmp(mak%ne),stat=istat)
    indexElmtG2L_tmp(:) = 0
!
    n = 0
    next_elmt : do i=1,mak%ne
      if(indexElmtG2L(i) == 0) cycle
      do j=1,n
        if(mak%etyp(i) == indexType(j)) then
          counter(j) = counter(j) + 1
          cycle next_elmt
        endif
      enddo
      n = n + 1
      indexType(n) = mak%etyp(i)
      counter(n) = counter(n) + 1
    enddo next_elmt
    if(n == 1.or.n == 0) then
      continue
      return
    endif
    if(n == 0) then
      print *,'number of element types',n
    endif
    allocate(ptr(n),stat=istat)
    ptr(1) = 1
    do i=1,n-1
      ptr(i+1) = counter(i) + ptr(i)
    enddo
    counter(:) = 0
    do i=1,mak%ne
      if(indexElmtG2L(i) == 0) cycle
      do j=1,n
        indexType(j) = mak%etyp(i)
        exit
      enddo
      counter(j) = counter(j) + 1
      indexElmtG2L_tmp(i) = ptr(j)-1+counter(j)
    enddo
    indexElmtG2L(:) = indexElmtG2L_tmp(:)
!
    deallocate(indexElmtG2L_tmp,stat=istat)
    deallocate(indexType,stat=istat)
    deallocate(counter,stat=istat)
    deallocate(ptr,stat=istat)   
end subroutine Mak_ResetIndexElmtG2L_Type

subroutine Mak_SetLocalMesh(mak,ni,nx,ei,ex,indexNodeG2L,indexElmtG2L,partID,mak_loc)
  type(MAKROSE_STRUCT),intent(in)   ::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  ni,nx,ei,ex
  integer(kint),intent(in)          ::  indexNodeG2L(:),indexElmtG2L(:)
  integer(kint),intent(in)          ::  partID  ! PartID of local mesh
  type(MAKROSE_STRUCT),intent(out)  ::  mak_loc ! Local mesh
!
  integer(kint)   ::  i,j,istat,n,egid
  integer(kint),allocatable         ::  indexNodeL2G(:),indexElmtL2G(:)
!
    call Mak_Init(mak_loc,mak%ndim,ni+nx,ei+ex)
    mak_loc%nn_i = ni
    mak_loc%ne_i = ei
    do i=1,mak%nn
      if(indexNodeG2L(i) == 0) cycle
      mak_loc%x(1:3,indexNodeG2L(i)) = mak%x(1:3,i)
!      if(associated(mak%ngid)) then
!        mak_loc%ngid(indexNodeG2L(i)) = mak%ngid(i)
!      else
        mak_loc%ngid(indexNodeG2L(i)) = i
!      endif
    enddo
    allocate(mak_loc%eptr(mak_loc%ne+1),stat=istat)
    do i=1,mak%ne
      if(indexElmtG2L(i) == 0) cycle
!      if(associated(mak%egid)) then
!        mak_loc%egid(indexElmtG2L(i)) = mak%egid(i)
!      else
        mak_loc%egid(indexElmtG2L(i)) = i
!      endif
      mak_loc%etyp(indexElmtG2L(i)) = mak%etyp(i)
      mak_loc%egrp(indexElmtG2L(i)) = mak%egrp(i)
      mak_loc%emat(indexElmtG2L(i)) = mak%emat(i)
    enddo
    mak_loc%eptr(1) = 1
    do i=1,mak_loc%ne
      egid = mak_loc%egid(i)
      if(egid < 1.or.egid > mak%ne) then
        print *,partID,i,egid
      endif
      mak_loc%eptr(i+1) = mak%eptr(egid+1) - mak%eptr(egid) + mak_loc%eptr(i)
    enddo
    allocate(mak_loc%eind(mak_loc%eptr(mak_loc%ne+1)-1),stat=istat)
    n = 0
    do i=1,mak_loc%ne
      egid = mak_loc%egid(i)
      do j=mak%eptr(egid),mak%eptr(egid+1)-1
        n = n + 1
        if(i == 280) then
          continue
        endif
        if(indexNodeG2L(mak%eind(j)) > mak_loc%nn.or.indexNodeG2L(mak%eind(j)) == 0) then
          continue
        endif
        mak_loc%eind(n) = indexNodeG2L(mak%eind(j))
      enddo
    enddo    
end subroutine Mak_SetLocalMesh

!< Element-based partition
subroutine Mak_GetNeighborIndex_ElementBase(mak,nparts,part,partID,nbIndex)
  type(MAKROSE_STRUCT),intent(inout)::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%nn))
  integer(kint),intent(in)          ::  partID
  integer(kint),pointer             ::  nbIndex(:)  
!
  integer(kint)   ::  i,j,k,n,istat,elmtID
  integer(kint),allocatable   ::  help(:),idx(:)
!
    if(.not.(associated(mak%nptr).and.associated(mak%nind))) then
      call Mak_GetElementOnNode(mak)
    endif
    if(.not.allocated(helpElmt2)) then
      allocate(helpElmt2(nparts,nparts),stat=istat)
      helpElmt2(:,:) = 0
      allocate(help(nparts),stat=istat)
      allocate(idx(nparts),stat=istat)
      do i=1,mak%nn
        idx(:) = 0; n = 0; help(:) = 0
        do j=mak%nptr(i),mak%nptr(i+1)-1
          elmtID = mak%nind(j)
          if(help(part(elmtID)) /= 0) cycle
          n = n + 1
          idx(n) = part(elmtID)
          help(part(elmtID)) = 1
        enddo
        if(n == 1) cycle
        do j=1,n
          do k=1,n
            if(k == j) cycle
            if(helpElmt2(idx(k),idx(j)) /= 0) cycle
            helpElmt2(idx(k),idx(j)) = 1
          enddo
        enddo
      enddo
      deallocate(help,idx,stat=istat)
    endif
    n = 0
    do i=1,nparts
      if(i == partID) cycle
      if(helpElmt2(i,partID) /= 0) n = n + 1
    enddo
    if(associated(nbIndex)) deallocate(nbIndex,stat=istat)
    allocate(nbIndex(n),stat=istat)
    n = 0
    do i=1,nparts
      if(i == partID) cycle
      if(helpElmt2(i,partID) /= 0) then
        n = n + 1
        nbIndex(n) = i
      endif
    enddo    
end subroutine Mak_GetNeighborIndex_ElementBase

subroutine Mak_GetLocalMeshMapping_ElementBase(mak,nparts,part,  &
            ni_loc,nx_loc,ei_loc,ex_loc,indexNodeG2L,indexElmtG2L)
  type(MAKROSE_STRUCT),intent(inout)::  mak     ! Total global mesh
  integer(kint),intent(in)          ::  nparts  ! Number of partitions
  integer(kint),intent(in)          ::  part(:) ! Nodal partition IDs (part(mak%ne))
  integer(kint),pointer             ::  ni_loc(:) ! Number of internal nodes in each part
  integer(kint),pointer             ::  nx_loc(:) ! Number of nodes on partition interface in each part
  integer(kint),pointer             ::  ei_loc(:) ! Number of elements whose nodes are internal 
  integer(kint),pointer             ::  ex_loc(:) ! Number of elements on partition interface
  integer(kint),pointer             ::  indexNodeG2L(:,:)
  integer(kint),pointer             ::  indexElmtG2L(:,:)
!
  integer(kint)   ::  i,j,k,n,istat,nodeID,nodeID1,nodeID2,partID,elmtID
  integer(kint),allocatable   ::  helpNode(:,:),helpElmt(:,:),help(:),marker(:)
!
!   Initiate pointers
    if(associated(ni_loc)) deallocate(ni_loc,stat=istat)
    if(associated(nx_loc)) deallocate(nx_loc,stat=istat)
    if(associated(ei_loc)) deallocate(ei_loc,stat=istat)
    if(associated(ex_loc)) deallocate(ex_loc,stat=istat)
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)
!
    allocate(ni_loc(nparts),stat=istat);  ni_loc(:) = 0
    allocate(nx_loc(nparts),stat=istat);  nx_loc(:) = 0
    allocate(ei_loc(nparts),stat=istat);  ei_loc(:) = 0
    allocate(ex_loc(nparts),stat=istat);  ex_loc(:) = 0
    allocate(indexNodeG2L(mak%nn,nparts),stat=istat);  indexNodeG2L(:,:) = 0
    allocate(indexElmtG2L(mak%ne,nparts),stat=istat);  indexElmtG2L(:,:) = 0
!
    if(.not.(associated(mak%nptr).and.associated(mak%nind))) then
      call Mak_GetElementOnNode(mak)
    endif
!   Internal nodes
    allocate(help(nparts),stat=istat)
    allocate(marker(mak%nn),stat=istat)
    marker(:) = 0
    do i=1,mak%nn
      n = 0; help(:) = 0
      do j=mak%nptr(i),mak%nptr(i+1)-1
        elmtID = mak%nind(j)
        if(help(part(elmtID)) /= 0) cycle
        help(part(elmtID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) then
        marker(i) = 1
        ni_loc(partID) = ni_loc(partID) + 1
        indexNodeG2L(i,partID) = ni_loc(partID)
        cycle
      endif
    enddo
!   External nodes
    do i=1,mak%nn
      if(marker(i) /= 0) cycle
      n = 0; help(:) = 0
      do j=mak%nptr(i),mak%nptr(i+1)-1
        elmtID = mak%nind(j)
        if(help(part(elmtID)) /= 0) cycle
        help(part(elmtID)) = 1
        partID = part(nodeID)
        nx_loc(partID) = nx_loc(partID) + 1
        indexNodeG2L(i,partID) = nx_loc(partID) + ni_loc(partID)
      enddo
    enddo
!   External nodes
    allocate(helpNode(mak%nn,nparts),stat=istat); helpNode(:,:) = 0
    allocate(helpElmt(mak%ne,nparts),stat=istat); helpElmt(:,:) = 0
    allocate(help(nparts),stat=istat)
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) then
        ei_loc(partID) = ei_loc(partID) + 1
        indexElmtG2L(i,partID) = ei_loc(partID)
        cycle
      endif
    enddo
    do i=1,mak%ne
      n = 0; help(:) = 0
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID = mak%eind(j)
        if(help(part(nodeID)) /= 0) cycle
        help(part(nodeID)) = 1
        n = n + 1
        partID = part(nodeID)
      enddo
      if(n == 1) cycle
      do j=mak%eptr(i),mak%eptr(i+1)-1
        nodeID1 = mak%eind(j)
        do k=mak%eptr(i),mak%eptr(i+1)-1
          nodeID2 = mak%eind(k)
          if(nodeID1 == nodeID2.or.part(nodeID1) == part(nodeID2)) cycle
          if(helpNode(nodeID2,part(nodeID1)) == 0) then
            helpNode(nodeID2,part(nodeID1)) = part(nodeID2)
            nx_loc(part(nodeID1)) = nx_loc(part(nodeID1)) + 1
            indexNodeG2L(nodeID2,part(nodeID1)) = ni_loc(part(nodeID1)) + nx_loc(part(nodeID1))
          endif
          if(helpNode(nodeID1,part(nodeID2)) == 0) then
            helpNode(nodeID1,part(nodeID2)) = part(nodeID1)
            nx_loc(part(nodeID2)) = nx_loc(part(nodeID2)) + 1
            indexNodeG2L(nodeID1,part(nodeID2)) = ni_loc(part(nodeID2)) + nx_loc(part(nodeID2))
          endif
          if(helpElmt(i,part(nodeID1)) == 0) then
            helpElmt(i,part(nodeID1)) = part(nodeID2)
            ex_loc(part(nodeID1)) = ex_loc(part(nodeID1)) + 1
            indexElmtG2L(i,part(nodeID1)) = ei_loc(part(nodeID1)) + ex_loc(part(nodeID1))
          endif
          if(helpElmt(i,part(nodeID2)) == 0) then
            helpElmt(i,part(nodeID2)) = part(nodeID1)
            ex_loc(part(nodeID2)) = ex_loc(part(nodeID2)) + 1
            indexElmtG2L(i,part(nodeID2)) = ei_loc(part(nodeID2)) + ex_loc(part(nodeID2))
          endif
        enddo
      enddo      
    enddo
    deallocate(helpNode,helpElmt,help,stat=istat) 
end subroutine Mak_GetLocalMeshMapping_ElementBase

end module m_MakrosePartMesh
