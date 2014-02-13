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
!> \brief This module contains subroutines for mesh partition and
!   frontISTR distributed mesh structure generation

module m_fstr_para_contact
  use m_fstr
  use m_Metis403API
  use m_MakrosePartMesh
  implicit none
  
  integer,parameter   ::  PARTITION_DEFAULT   = 0
  integer,parameter   ::  PARTITION_PMETIS    = 1
  integer,parameter   ::  PARTITION_KMETIS    = 2
  integer,parameter   ::  PARTITION_MCPMETIS  = 11
  integer,parameter   ::  PARTITION_MCKMETIS  = 21
  
  integer(kint),save              ::  npe_import = 0
  integer(kint),allocatable,save  ::  import_pe(:)
  integer(kint),allocatable,save  ::  import_index(:)
  integer(kint),allocatable,save  ::  import_item(:)
  integer(kint),save              ::  npe_export = 0
  integer(kint),allocatable,save  ::  export_pe(:)
  integer(kint),allocatable,save  ::  export_index(:)
  integer(kint),allocatable,save  ::  export_item(:)
  
contains

subroutine paraContact_DomainPartition(hecMESH_G,hecMESH_L)
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_G
  type(hecmwST_local_mesh),intent(out)  ::  hecMESH_L
!
  type(MAKROSE_STRUCT),save    ::  mak  ,mak_loc
  integer(kint),pointer   ::  part(:)=>null()
  integer(idx_t),pointer  ::  xadj(:)=>null(),adjncy(:)=>null()
  integer(idx_t),pointer  ::  vwgt(:)=>null(),adjwgt(:)=>null(),vsize(:)=>null()
  integer(kint),pointer   ::  options(:)=>null()
  real(real_t),pointer    ::  tpwgts(:)=>null(),ubvec(:)=>null()
  integer(idx_t)   ::  nvtxs,ncon,wgtflag=0,numflag,objval
!
  integer(kint)           ::  i,j,n,nparts,istat,partID
  integer(kint),allocatable ::  ind(:),help(:),metisPartResult(:)
  integer(kint),pointer   ::  indexNodeG2L(:,:)=>null(),indexElmtG2L(:,:)=>null()
  integer(kint)   ::  ierr,maxSize,ppid,count

  character(len=256)      ::  filename,lmsh_file,name
  integer(kint)   ::  partMethod = PARTITION_DEFAULT
  character(len=6)        ::  actualPartMethod
!
    call paraContact_GetFSTR2Makrose(hecMESH_G,mak)

    if(myrank == 0) then
!    call Mak_WriteFile(mak,'mesh')
      call Mak_MeshToNodal(mak,xadj,adjncy)
      allocate(part(hecMESH_G%n_node),stat=istat)
      nparts = nprocs
      nvtxs = mak%nn
      allocate(options(5),stat=istat)
      options(:) = 0
      numflag = 1
      wgtflag = 0
!      ncon = 1
!      goto 133
      ncon = 2
      allocate(vwgt(hecMESH_G%n_node*ncon),stat=istat)
      allocate(help(hecMESH_G%n_node),stat=istat)
      help(:) = 0
      call paraContact_MarkMasterNode(hecMESH_G,help)
      do i=1,hecMESH_G%n_node
        vwgt((i-1)*ncon+1) = 1
        vwgt((i-1)*ncon+2) = 0
        if(help(i) /= 0) then
          vwgt((i-1)*ncon+2) = 1
        endif
      enddo
      deallocate(help,stat=istat)
133   continue      
!
!      if(myrank == 0) print *,'Part Method',partMethod
      select case(partMethod)
      case(PARTITION_DEFAULT)
        if(nparts < 8) then
          if(ncon < 2) then
            ierr = METIS_PartGraphRecursive                   &
                        (nvtxs,   xadj,   adjncy,             &
                         vwgt,    adjwgt, wgtflag,numflag,    &
                         nparts,  options,objval, part)
            actualPartMethod = 'PMETIS'
          else
            ierr = METIS_mCPartGraphRecursive                 &
                        (nvtxs,   ncon,   xadj,   adjncy,     &
                         vwgt,    adjwgt, wgtflag,numflag,    &
                         nparts,  options,objval, part)
            actualPartMethod = 'MCPMETIS'        
          endif
        else
          if(ncon < 2) then
            ierr = METIS_PartGraphKway                        &
                        (nvtxs,   xadj,   adjncy,             &
                         vwgt,    adjwgt, wgtflag,numflag,    &
                         nparts,  options,objval, part)
            actualPartMethod = 'KMETIS'
          else
            allocate(ubvec(ncon),stat=istat)
            ubvec(:) = 1.03D0
            ierr = METIS_mCPartGraphKway                      &
                        (nvtxs,   ncon,   xadj,   adjncy,     &
                         vwgt,    adjwgt, wgtflag,numflag,    &
                         nparts,  ubvec,  options,objval, part) 
            actualPartMethod = 'MCKMETIS'     
          endif  
        endif
      case(PARTITION_PMETIS)
          ierr = METIS_PartGraphRecursive                   &
                      (nvtxs,   xadj,   adjncy,             &
                       vwgt,    adjwgt, wgtflag,numflag,    &
                       nparts,  options,objval, part)
          actualPartMethod = 'PMETIS'    
      case(PARTITION_KMETIS)
          ierr = METIS_PartGraphKway                        &
                      (nvtxs,   xadj,   adjncy,             &
                       vwgt,    adjwgt, wgtflag,numflag,    &
                       nparts,  options,objval, part) 
          actualPartMethod = 'KMETIS'
      case(PARTITION_MCPMETIS)
          ierr = METIS_mCPartGraphRecursive                 &
                      (nvtxs,   ncon,   xadj,   adjncy,     &
                       vwgt,    adjwgt, wgtflag,numflag,    &
                       nparts,  options,objval, part)
          actualPartMethod = 'MCPMETIS'    
      case(PARTITION_MCKMETIS)
          if(ncon <= 1) stop
          allocate(ubvec(ncon),stat=istat)
          ubvec(:) = 1.03D0
          ierr = METIS_mCPartGraphKway                      &
                      (nvtxs,   ncon,   xadj,   adjncy,     &
                       vwgt,    adjwgt, wgtflag,numflag,    &
                       nparts,  ubvec,  options,objval, part) 
          actualPartMethod = 'MCKMETIS'  
      case default
        stop 'Error: Undefined partition method!'
      end select
!
!      call MPI_BARRIER(hecMESH_G%MPI_COMM,ierr)
      print *,myrank,'Metis Over'
    else
!     Restart and read existing partition for local mesh construction
      allocate(part(hecMESH_G%n_node),stat=istat)
!     Read existing partitions
      nparts = nprocs
    endif
    if(nprocs > 1) then
      call hecmw_bcast_I(hecMESH_G,part,hecMESH_G%n_node,0)
!      call MPI_BCAST(part,hecMESH_G%n_node,MPI_INTEGER,0,hecMESH_G%MPI_COMM,istat)
    endif
    call hecmw_BARRIER(hecMESH_G)
    
!#ifdef OUTPUT_PARTITION
!    if(myrank == 0) then
!      do i=1,mak%ne
!        ppid = 0
!        count = 0
!        do j=mak%eptr(i),mak%eptr(i+1)-1
!          if(part(mak%eind(j)) /= ppid) then
!            ppid = part(mak%eind(j))
!            count = count + 1
!          endif
!        enddo
!        if(count > 1) then
!          mak%egrp(i) = 1
!        elseif(count == 1) then
!          mak%egrp(i) = ppid + 1
!        else
!          stop 'error:'
!        endif
!      enddo
!      write(name,'(A,I2.2)')'full_mesh_init_np',nparts
!      call Mak_WriteFile(mak,name)
!    endif
!#endif

!   Get local hecMESH from global one and partition info
    partID = myrank + 1
!    call Mak_GetLocalMesh_NodeBase1_all(mak,nparts,part,partID,mak_loc,indexNodeG2L,indexElmtG2L)
!    do partID=1,nparts
!    print *,myrank,'mak',mak%nn,mak%ne,nparts,partID
    call Mak_GetLocalMesh_NodeBase1_FrontISTR(mak,nparts,part,partID,mak_loc,indexNodeG2L,indexElmtG2L)

    
!    call rtri_GetLocalMesh_all(hecMESH_G,mak_loc,part,partID,hecMESH_L,indexNodeG2L,indexElmtG2L)
    call paraContact_GetLocalMesh_all_new(hecMESH_G,mak_loc,part,partID,hecMESH_L,indexNodeG2L,indexElmtG2L)
    call hecmw_BARRIER(hecMESH_G)
    call paraContact_CreateExportImport(hecMESH_L)
    if(myrank == 0) then
!      write(*,'(A,I15,A,I15)')'number of edgecut :',objval,'/',size(adjncy)/2
      write(*,'(A,A)')'Partition Method: ',actualPartMethod(1:6)
      write(*,'(A8,4A15)')' rankID','          nodes',' internal nodes','          elems',' internal elems'
    endif
!
    call hecmw_BARRIER(hecMESH_G)
    write(*,'(I8,5I15)')myrank,hecMESH_L%n_node,hecMesh_L%nn_internal,hecMESH_L%n_elem,hecMESH_L%ne_internal,hecMesh_L%nn_middle
!    enddo
!
    call hecmw_BARRIER(hecMESH_G)
!    stop
     
    if(associated(indexNodeG2L)) deallocate(indexNodeG2L,stat=istat)
    if(associated(indexElmtG2L)) deallocate(indexElmtG2L,stat=istat)
        
!   deallocation pointers and pointer-related structures
    if(associated(xadj)) deallocate(xadj,stat=istat)
    if(associated(adjncy)) deallocate(adjncy,stat=istat)

    if(associated(options)) deallocate(options,stat=istat)
    if(allocated(ind)) deallocate(ind,stat=istat)
    if(allocated(help)) deallocate(help,stat=istat)
    if(allocated(metisPartResult)) deallocate(metisPartResult,stat=istat)
    if(associated(part)) deallocate(part,stat=istat)
!
    call Mak_Init(mak)
    call Mak_Init(mak_loc)
end subroutine paraContact_DomainPartition

subroutine paraContact_GetFSTR2Makrose(hecMESH,mak)
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH
  type(MAKROSE_STRUCT),intent(out)      ::  mak
!
  integer   ::  i,istat
!   Check if hecMESH is global mesh
    if(hecMESH%n_node /= hecMESH%nn_internal.or.  &
       hecMESH%n_elem /= hecMESH%ne_internal) then
      print *,'myrank',myrank,' is not a global mesh!'
      stop
    endif
    do i=1,hecMESH%n_node
      if(i /= hecMESH%global_node_ID(i)) then
        print *,'myrank',myrank,' Global nodeIDs are not consecutive!'
        stop
      endif
    enddo
    call Mak_init(mak,ndim=hecMESH%n_dof,nn=hecMESH%n_node,ne=hecMESH%n_elem)
    do i=1,hecMESH%n_node
      mak%x(1:hecMESH%n_dof,i) = hecMESH%node((i-1)*hecMESH%n_dof+1:i*hecMESH%n_dof)
      mak%ngid(i) = hecMESH%global_node_ID(i)
    enddo
    allocate(mak%eptr(hecMESH%n_elem+1),stat=istat)
    allocate(mak%eind(hecMESH%elem_node_index(hecMESH%n_elem)),stat=istat)
    mak%eptr(1) = 1
    do i=1,hecMESH%n_elem
      mak%eptr(i+1) = hecMESH%elem_node_index(i)+1
      mak%eind(mak%eptr(i):mak%eptr(i+1)-1) = hecMESH%elem_node_item(hecMESH%elem_node_index(i-1)+1:hecMESH%elem_node_index(i))
      mak%egid(i) = hecMESH%global_elem_ID(i)
      mak%egrp(i) = hecMESH%section_ID(i)
      mak%emat(i) = hecMESH%elem_mat_ID_item(i)
      select case(hecMESH%elem_type(i))
      case(361)
        mak%etyp(i) = MAKROSE_HEXA8
      case(351)
        mak%etyp(i) = MAKROSE_PRI6
      case default
        print *,'Element type ID ',hecMESH%elem_type(i),' is not yet implemented for RTRI applications!'
        stop
      end select
    enddo
end subroutine paraContact_GetFSTR2Makrose

subroutine paraContact_GetLocalElementType(mak,ntype,ptrtype,itemtype)
! Elements of same types must be continuously ordered in mak%etype
  type(MAKROSE_STRUCT),intent(in)   ::  mak
  integer(kint),intent(out)         ::  ntype
  integer(kint),pointer             ::  ptrtype(:)
  integer(kint),pointer             ::  itemtype(:)

!
  integer(kint)   ::  i,istat,etype
  integer(kint),allocatable   ::  counter(:)
!
    allocate(counter(20),stat=istat)
    counter(:) = 0
    ntype = 1
    etype = mak%etyp(1)
    counter(1) = 1
    do i=2,mak%ne
      if(etype /= mak%etyp(i)) then
          ntype = ntype + 1
          etype = mak%etyp(i)
      else
        counter(ntype) = counter(ntype) + 1
      endif
    enddo
    if(associated(ptrtype)) deallocate(ptrtype,stat=istat)
    allocate(ptrtype(0:ntype),stat=istat)
    ptrtype(0) = 0
    do i=1,ntype
      ptrtype(i) = ptrtype(i-1) + counter(i)
    enddo
    if(associated(itemtype)) deallocate(itemtype,stat=istat)
    allocate(itemtype(ntype),stat=istat)
    do i=1,ntype
      select case(mak%etyp(ptrtype(i)))
      case(MAKROSE_PRI6)
        itemtype(i) = 351
      case(MAKROSE_HEXA8)
        itemtype(i) = 361
      case default
        print *,'Element type ID ',mak%etyp(ptrtype(i)),' is not yet implemented for RTRI applications!'
        stop
      end select
    enddo
end subroutine paraContact_GetLocalElementType

subroutine paraContact_MarkMasterNode(hecMESH,help)
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH
  integer(kint),intent(inout)           ::  help(:)
!
  integer(kint)   ::  i,j,ic,ic_type,iss,mm,surf_grp_ID,node_grp_ID
!
    help(:) = 0
!
    do i=1,hecMESH%contact_pair%n_pair
      surf_grp_ID = hecMESH%contact_pair%master_grp_id(i)
      if(hecMESH%surf_group%grp_index(surf_grp_ID) == hecMESH%surf_group%grp_index(surf_grp_ID-1)+1) cycle
!     This partation has contact pair master group
      do j=hecMESH%surf_group%grp_index(surf_grp_ID-1)+1,hecMESH%surf_group%grp_index(surf_grp_ID)
        ic   = hecMESH%surf_group%grp_item(2*j-1)
        ic_type = hecMESH%elem_type(ic)
        iss = hecMESH%elem_node_index(ic-1)
        do mm=1,getNumberOfNodes( ic_type )
          if(help(hecMESH%elem_node_item(iss+mm)) == 0) then
            help(hecMESH%elem_node_item(iss+mm)) = 1
          endif
        enddo
      enddo
    enddo
    
end subroutine paraContact_MarkMasterNode

subroutine paraContact_GetLocalMesh_all_new(hecMESH_G,mak_loc,part,partID,hecMESH_L,  &
                                        indexNodeG2L,indexElmtG2L)
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_G
  type(MAKROSE_STRUCT),intent(inout)    ::  mak_loc
  integer(kint),intent(in)              ::  part(:)
  integer(kint),intent(in)              ::  partID    ! partition ID
  type(hecmwST_local_mesh),intent(inout)::  hecMESH_L
  integer(kint),intent(in)              ::  indexNodeG2L(:,:),indexElmtG2L(:,:)
!
  integer(kint) ::  i,j,n,istat,pid,nodeID,elemID,nodeID1,ii,jj,maxNum,k,l
  integer(kint),allocatable   ::  neighbor(:),help(:),helpnode(:,:),helpelem(:,:),  &
                                  helpCountNode(:),helpCountElem(:)
  integer(kint),allocatable   ::  markelem(:),marknode(:),temp(:),markPE(:),indexNodeG2L_Hanging(:),exportSlavePE(:)
  integer(kint),pointer   ::  indexHangingNode(:)=>null()
  integer(kint)           ::  nHangingNode = 0,nOtherSlaveNode = 0,count
  logical ::  sameflag
  integer ::  buf(2),rev(2),ierr
!
    call hecmw_nullify_mesh(hecMESH_L)
!-- General data
    hecMESH_L%gridfile              = hecMESH_G%gridfile
    hecMESH_L%hecmw_n_file          = hecMESH_G%hecmw_n_file
    !allocate(hecMESH_L%files(size(hecMESH_G%files)),stat=istat)
    !hecMESH_L%files                 = hecMESH_G%files
    hecMESH_L%header                = hecMESH_G%header
    hecMESH_L%hecmw_flag_adapt      = hecMESH_G%hecmw_flag_adapt
    hecMESH_L%hecmw_flag_initcon    = hecMESH_G%hecmw_flag_initcon
    
    hecMESH_L%hecmw_flag_parttype   = hecMESH_G%hecmw_flag_parttype
    hecMESH_L%hecmw_flag_partdepth  = hecMESH_G%hecmw_flag_partdepth
    hecMESH_L%hecmw_flag_version    = hecMESH_G%hecmw_flag_version
    hecMESH_L%zero_temp             = hecMESH_G%zero_temp
!
!-- Surface group data
    maxNum = 0
    do i=1,hecMESH_G%surf_group%n_grp
      maxNum = max(maxNum,hecMESH_G%surf_group%grp_index(i)- hecMESH_G%surf_group%grp_index(i-1))
    enddo
    call paraContact_getHecmwLocalSurfaceGroup(hecMESH_G%surf_group,hecMESH_L%surf_group, &
         indexElmtG2L(:,partID),maxNum,hecMESH_G,part,partID)

!-- Contact pair data
    call paraContact_copyHecmwContactPair(hecMESH_G%contact_pair,hecMESH_L%contact_pair)

!-- Get temp hanging slave nodes which are not specified by partition info.
    call paraContact_GetHangingSlaveNode(hecMESH_G,hecMesh_L,indexNodeG2L(:,partID),indexHangingNode,nHangingNode,nOtherSlaveNode)
!    print *,myrank,'nHangingNode',nHangingNode,nOtherSlaveNode
!    print *,myrank,'nHangingNode',nHangingNode,'index',indexHangingNode(1:nHangingNode)
!    print *,myrank,'nHangingNode',nHangingNode,'part ',part(indexHangingNode(1:nHangingNode))
!    nHangingNode = 0    ! debug
!

!-- Node data
    hecMESH_L%n_node                = mak_loc%nn
    hecMESH_L%nn_middle             = mak_loc%nn
    hecMESH_L%n_node                = mak_loc%nn + nHangingNode
    hecMESH_L%n_node_gross          = hecMESH_L%n_node !hecMESH_G%n_node_gross    !?
    hecMESH_L%nn_internal           = mak_loc%nn_i
    hecMESH_L%n_dof                 = hecMESH_G%n_dof
    hecMESH_L%n_dof_grp             = hecMESH_G%n_dof_grp
    hecMESH_L%n_dof_tot             = hecMESH_G%n_dof_tot       !?
    allocate(hecMESH_L%node(hecMESH_L%n_node*hecMESH_L%n_dof),stat=istat)
    do i=1,hecMESH_L%n_node
      if(i <= mak_loc%nn) then
        hecMESH_L%node((i-1)*hecMESH_L%n_dof+1:i*hecMESH_L%n_dof) = mak_loc%x(1:hecMESH_L%n_dof,i)
      else
        j = i - mak_loc%nn
        hecMESH_L%node((i-1)*hecMESH_L%n_dof+1:i*hecMESH_L%n_dof)   &
          = hecMESH_G%node((indexHangingNode(j)-1)*hecMESH_L%n_dof+1:indexHangingNode(j)*hecMESH_L%n_dof)
      endif
    enddo
    allocate(hecMESH_L%node_ID(hecMESH_L%n_node*2),stat=istat)
    do i=1,hecMESH_L%n_node
      if(i <= mak_loc%nn) then
        if(part(mak_loc%ngid(i)) == partID) then
          hecMESH_L%node_ID((i-1)*2+1)  = i
        else
          hecMESH_L%node_ID((i-1)*2+1)  = indexNodeG2L(mak_loc%ngid(i),part(mak_loc%ngid(i)))
        endif
        hecMESH_L%node_ID(i*2)        = part(mak_loc%ngid(i)) - 1
      else
        j = i - mak_loc%nn
        if(part(indexHangingNode(j)) == partID) then
          stop
        else
          hecMESH_L%node_ID((i-1)*2+1)  = indexNodeG2L(indexHangingNode(j),part(indexHangingNode(j)))
        endif
        hecMESH_L%node_ID(i*2)        = part(indexHangingNode(j)) - 1
      endif
    enddo
    allocate(hecMESH_L%global_node_ID(hecMESH_L%n_node),stat=istat)
    hecMESH_L%global_node_ID(1:hecMESH_L%nn_middle) = mak_loc%ngid(1:mak_loc%nn)
!   debug
    if(nHangingNode > 0.and.size(indexHangingNode) == nHangingNode) then
      hecMESH_L%global_node_ID(hecMESH_L%nn_middle+1:hecMESH_L%n_node) = indexHangingNode(1:nHangingNode)
    endif
    allocate(hecMESH_L%node_dof_index(0:hecMESH_L%n_dof_grp),stat=istat)
    hecMESH_L%node_dof_index(0) = 0
    hecMESH_L%node_dof_index(1) = hecMESH_L%n_node
    allocate(hecMESH_L%node_dof_item(hecMESH_L%n_dof_grp),stat=istat)
    hecMESH_L%node_dof_item(:) = hecMESH_L%n_dof
    allocate(hecMESH_L%node_internal_list(hecMESH_L%nn_internal),stat=istat)
    n = 0
    do i=1,hecMESH_L%n_node
      if(i > mak_loc%nn) exit
      if(part(mak_loc%ngid(i)) == partID) then
        n = n + 1
        hecMESH_L%node_internal_list(n) = i
      endif
    enddo
    
!    write(myrank+30,*)hecMESH_L%nn_internal,hecMESH_L%n_node
!    do i=1,hecMESH_L%n_node
!      write(myrank+30,*)i,hecMESH_L%node((i-1)*3+1:i*3)
!    enddo
!    do i=1,hecMESH_L%n_node 
!      if(i > hecMESH_L%nn_internal) then
!        write(myrank+30,*)hecMESH_L%node_ID((i-1)*2+1:i*2), i
!      else
!        write(myrank+30,*)hecMESH_L%node_ID((i-1)*2+1:i*2)
!      endif
!    enddo
!    write(myrank+30,'(10I6)')hecMESH_L%global_node_ID(:)
    
!-- Node group data
    maxNum = 0
    do i=1,hecMESH_G%node_group%n_grp
      maxNum = max(maxNum,hecMESH_G%node_group%grp_index(i)- hecMESH_G%node_group%grp_index(i-1))
    enddo
!   debug
    allocate(indexNodeG2L_Hanging(hecMESH_G%n_node),stat=istat)
    indexNodeG2L_Hanging(:) = indexNodeG2L(:,partID)
    do i=1,nHangingNode
      if(indexNodeG2L_Hanging(indexHangingNode(i)) /= 0) then
        stop
      else
        indexNodeG2L_Hanging(indexHangingNode(i)) = hecMESH_L%nn_middle + i
      endif
    enddo
    call paraContact_getHecmwLocalNodeGroup(hecMESH_G%node_group,hecMESH_L%node_group,indexNodeG2L_Hanging,maxNum)
    deallocate(indexNodeG2L_Hanging,stat=istat)
    
!-- Element data 
    hecMESH_L%n_elem                = mak_loc%ne
    hecMESH_L%n_elem_gross          = hecMESH_L%n_elem  !hecMESH_G%n_elem_gross      ! ?
    hecMESH_L%ne_internal           = mak_loc%ne_i
    call paraContact_GetLocalElementType(mak_loc,hecMESH_L%n_elem_type,hecMESH_L%elem_type_index,hecMESH_L%elem_type_item)
    allocate(hecMESH_L%elem_type(hecMESH_L%n_elem),stat=istat)
    do i=1,hecMESH_L%n_elem_type
      hecMESH_L%elem_type(hecMESH_L%elem_type_index(i-1)+1:hecMESH_L%elem_type_index(i)) = hecMESH_L%elem_type_item(i)
    enddo
    hecMESH_L%n_elem_mat_ID = hecMESH_L%n_elem
    allocate(hecMESH_L%elem_mat_ID_index(0:hecMESH_L%n_elem),stat=istat)
    do i=0,hecMESH_L%n_elem
      hecMESH_L%elem_mat_ID_index(i) = i
    enddo
    allocate(hecMESH_L%elem_mat_ID_item(hecMESH_L%n_elem),stat=istat)
    hecMESH_L%elem_mat_ID_item(:) = mak_loc%emat(:)
    allocate(hecMESH_L%elem_node_index(0:hecMESH_L%n_elem),stat=istat)
    do i=0,hecMESH_L%n_elem
      hecMESH_L%elem_node_index(i) = mak_loc%eptr(i+1)-1
    enddo
    allocate(hecMESH_L%elem_node_item(size(mak_loc%eind)),stat=istat)
    hecMESH_L%elem_node_item(:) = mak_loc%eind(:)
    allocate(hecMESH_L%global_elem_ID(hecMESH_L%n_elem),stat=istat)
    do i=1,hecMESH_L%n_elem
      hecMESH_L%global_elem_ID(i) = hecMESH_G%global_elem_ID(mak_loc%egid(i))
    enddo
    allocate(hecMESH_L%elem_ID(hecMESH_L%n_elem*2),stat=istat)
!    allocate(hecMESH_L%elem_internal_list(hecMESH_L%ne_internal),stat=istat)
    allocate(temp(hecMESH_L%n_elem),stat=istat)
    n = 0
    next_elem : do i=1,hecMESH_L%n_elem
       sameflag = .true.
       pid = part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(hecMESH_L%elem_node_index(i-1)+1)))
       do j=hecMESH_L%elem_node_index(i-1)+1,hecMESH_L%elem_node_index(i)
         if(pid /= part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(j)))) then
           if(sameflag) sameflag = .false.
           pid = min(pid,part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(j))))
         endif
       enddo
!       if(sameflag) then
!         n = n + 1
!         hecMESH_L%elem_internal_list(n) = i
!       endif
       hecMESH_L%elem_ID(i*2) = pid - 1
       if(hecMESH_L%elem_ID(i*2) + 1 == partID) then
         n = n + 1
         temp(n) = i
         hecMESH_L%elem_ID((i-1)*2+1) = n
       else
         hecMESH_L%elem_ID((i-1)*2+1) = indexElmtG2L(mak_loc%egid(i),pid)
!         hecMESH_L%elem_internal_list(n) = i
       endif
       
    enddo next_elem
    hecMESH_L%ne_internal = n
    allocate(hecMESH_L%elem_internal_list(hecMESH_L%ne_internal),stat=istat)
    hecMESH_L%elem_internal_list(1:n) = temp(1:n)
    deallocate(temp,stat=istat)
!    print *,'my_rank',partID-1, hecMESH_L%ne_internal,hecMESH_L%n_elem,hecMESH_L%nn_internal,hecMESH_L%n_node
        
    allocate(hecMESH_L%section_ID(hecMESH_L%n_elem),stat=istat)
    do i=1,hecMESH_L%n_elem
      hecMESH_L%section_ID(i) = hecMESH_G%section_ID(mak_loc%egid(i))
    enddo

!    write(myrank+30,*)'elements'
!    write(myrank+30,*)hecMESH_L%ne_internal,hecMESH_L%n_elem
!    do i=1,hecMESH_L%n_elem
!      write(myrank+30,'(20I10)')i,hecMESH_L%elem_node_item(hecMESH_L%elem_node_index(i-1)+1:hecMESH_L%elem_node_index(i))
!    enddo
!    do i=1,hecMESH_L%n_elem
!      write(myrank+30,*)hecMESH_L%elem_ID((i-1)*2+1:i*2)
!    enddo
!    write(myrank+30,'(10I6)')hecMESH_L%global_elem_ID(:)

!-- PE data
    hecMESH_L%PETOT = hecMESH_G%PETOT 
    if(partID == 1) then
      hecMESH_L%zero = 1
    else
      hecMESH_L%zero = 0
    endif
    hecMESH_L%MPI_COMM = hecMESH_G%MPI_COMM
    hecMESH_L%my_rank = myrank
    hecMESH_L%n_subdomain = hecMESH_L%PETOT
    hecMESH_L%n_neighbor_pe = 0
    allocate(neighbor(hecMESH_L%PETOT-1),stat=istat)
    i = hecMESH_L%PETOT
    allocate(help(i),stat=istat)
    allocate(helpnode(hecMESH_L%n_node - hecMESH_L%nn_internal,hecMESH_L%PETOT),stat=istat)
    helpnode(:,:) = 0
    help(:) = 0
!   Count Neighboring Processors
!   import_index, import_item
    do i=1,hecMESH_L%n_node
      pid = part(hecMESH_L%global_node_ID(i))
      if(pid /= partID) then
        help(pid) = help(pid) + 1
        helpnode(help(pid),pid) = i
      endif
    enddo
    do i=1,size(help)
      if(help(i) /= 0) then
!        print *,myrank,'neighborID',i-1,help(i)
        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
      endif
    enddo
!    do i=hecMESH_L%nn_middle+1,hecMESH_L%n_node
!      pid = part(hecMESH_L%global_node_ID(i))
!      if(pid /= partID) then
!        help(pid) = help(pid) + 1
!        helpnode(help(pid),pid) = i
!      endif
!    enddo
!    do i=1,size(help)
!      if(help(i) /= 0) then
!        print *,myrank,'neighborID',i-1,help(i)
!        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
!      endif
!    enddo
!    print *,myrank,'n_neighbor_pe',hecMESH_L%n_neighbor_pe
    hecMESH_L%n_neighbor_pe = hecMESH_L%PETOT - 1
    allocate(hecMESH_L%neighbor_pe(hecMESH_L%n_neighbor_pe),stat=istat)
!    hecMESH_L%n_neighbor_pe = 0
    count = 0
    do i=1,size(help)
!      if(help(i) /= 0) then
!        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
!        hecMESH_L%neighbor_pe(hecMESH_L%n_neighbor_pe) = i - 1  !start from 0
!      endif
      if(i == partID) then
        count = 1
        cycle
      endif
      hecMESH_L%neighbor_pe(i-count) = i - 1
    enddo
    allocate(hecMESH_L%import_index(0:hecMESH_L%n_neighbor_pe),stat=istat)
!    hecMESH_L%n_neighbor_pe = 0
    hecMESH_L%import_index(0) = 0
    count = 0
    do i=1,size(help)
!      if(help(i) /= 0) then
!        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1        
!        hecMESH_L%import_index(hecMESH_L%n_neighbor_pe) = hecMESH_L%import_index(hecMESH_L%n_neighbor_pe-1) + help(i)
!      endif
      if(i == partID) then
        count = 1
        cycle
      endif
      hecMESH_L%import_index(i-count) = hecMESH_L%import_index(i-count-1) + help(i)
    enddo
!    print *,myrank,'item',hecMESH_L%import_index(hecMESH_L%n_neighbor_pe)
    allocate(hecMESH_L%import_item(hecMESH_L%import_index(hecMESH_L%n_neighbor_pe)),stat=istat)   
!    hecMESH_L%n_neighbor_pe = 0
    count = 0
    do i=1,size(help)
!      if(help(i) /= 0) then
!        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
!        hecMESH_L%import_item(hecMESH_L%import_index(hecMESH_L%n_neighbor_pe-1)+1:hecMESH_L%import_index(hecMESH_L%n_neighbor_pe)) = helpnode(1:help(i),i)
!      endif
      if(i == partID) then
        count = 1
        cycle
      endif
      if(hecMESH_L%import_index(i-count-1) == hecMESH_L%import_index(i-count)) cycle
      hecMESH_L%import_item(hecMESH_L%import_index(i-count-1)+1:hecMESH_L%import_index(i-count)) = helpnode(1:help(i),i)
    enddo
!!  export_index, export_item
    if(.not.(associated(mak_loc%nptr).and.associated(mak_loc%nind))) then
      call Mak_GetElementOnNode(mak_loc)
    endif
!
    deallocate(helpnode)
    allocate(helpnode(hecMESH_L%nn_internal,hecMESH_L%n_neighbor_pe),stat=istat)
!    print *,'myrank',myrank,size(helpnode,1),hecMESH_L%nn_internal
    helpnode(:,:) = 0
    allocate(helpCountNode(hecMESH_L%n_neighbor_pe),stat=istat)
    helpCountNode(:) = 0
    allocate(markPE(hecMESH_L%PETOT),stat=istat)
    allocate(exportSlavePE(hecMESH_L%PETOT),stat=istat)
    
    do i=1,hecMESH_L%n_node
      if(hecMESH_L%node_ID(i*2) /= partID - 1) cycle
      markPE(:) = 0
      do j=mak_loc%nptr(i),mak_loc%nptr(i+1)-1
        elemID = mak_loc%nind(j)
        help(:) = 0
        n = 0
        do k=mak_loc%eptr(elemID),mak_loc%eptr(elemID+1)-1
          if(help(part(mak_loc%ngid(mak_loc%eind(k)))) == 0) then
            help(part(mak_loc%ngid(mak_loc%eind(k)))) = 1
            n = n + 1
          endif
        enddo
        if(n == 1) cycle
        do k=1,hecMESH_L%PETOT
          if(help(k) == 0) cycle
          if(k == partID) cycle
          if(markPE(k) /= 0) cycle
          markPE(k) = 1
          do l=1,hecMESH_L%n_neighbor_pe
            if(hecMESH_L%neighbor_pe(l) == k-1) goto 110
          enddo
          print *,'not found cycle'
          cycle
110       helpCountNode(l) = helpCountNode(l) + 1
          if( helpCountNode(l) > size(helpnode,1)) then
            print *,'PPP ',helpCountNode(l),size(helpnode,1),hecMESH_L%nn_internal
          endif
          helpnode(helpCountNode(l),l) = i
        enddo
      enddo
    enddo
    
    do i=1,hecMESH_L%n_node
      if(hecMESH_L%node_ID(i*2) /= partID - 1) cycle
      markPE(:) = 0
      do j=mak_loc%nptr(i),mak_loc%nptr(i+1)-1
        elemID = mak_loc%nind(j)
        help(:) = 0
        n = 0
        do k=mak_loc%eptr(elemID),mak_loc%eptr(elemID+1)-1
          if(help(part(mak_loc%ngid(mak_loc%eind(k)))) == 0) then
            help(part(mak_loc%ngid(mak_loc%eind(k)))) = 1
            n = n + 1
          endif
        enddo
        if(n == 1) cycle
        do k=1,hecMESH_L%PETOT
          if(help(k) == 0) cycle
          if(k == partID) cycle
          if(markPE(k) /= 0) cycle
          markPE(k) = 1
        enddo
      enddo
      if(paraContact_isContactSlaveNode(hecMESH_G,hecMesh_L,part,i,exportSlavePE)) then
!        write(myrank+30,'(A,   32I5)')'nodeLID ',i,markPE(:)
!        write(myrank+30,'(A,A5,32I5)')'        ','     ',exportSlavePE(:)
!        print *,myrank,i,'expPE',exportSlavePE(:)
!        cycle ! debug
        do k=1,hecMESH_L%PETOT
          if(k == partID) cycle
          if(markPE(k) == 0.and.exportSlavePE(k) /= 0) then
            do l=1,hecMESH_L%n_neighbor_pe
              if(hecMESH_L%neighbor_pe(l) == k-1) goto 111
            enddo
!            write(myrank+60,'(I10,40I5)')i,k
!            write(myrank+60,'(10X,40I5)')hecMESH_L%neighbor_pe(:)
!            write(myrank+60,'(10X,40I5)')markPE(:)
!            write(myrank+60,'(10X,40I5)')exportSlavePE(:)
            print *,'cycle'
            cycle
111         helpCountNode(l) = helpCountNode(l) + 1
            if( helpCountNode(l) > size(helpnode,1)) then
              print *,'PPP+',helpCountNode(l),size(helpnode,1),hecMESH_L%nn_internal
            endif
            helpnode(helpCountNode(l),l) = i
          endif
        enddo
      endif
    enddo

    allocate(hecMESH_L%export_index(0:hecMESH_L%n_neighbor_pe),stat=istat)
    hecMESH_L%export_index(0) = 0
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%export_index(i) = hecMESH_L%export_index(i-1) + helpCountNode(i)
    enddo
    allocate(hecMESH_L%export_item(hecMESH_L%export_index(hecMESH_L%n_neighbor_pe)),stat=istat)
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%export_item(hecMESH_L%export_index(i-1)+1:hecMESH_L%export_index(i)) = helpnode(1:helpCountNode(i),i)
!      print *,myrank,'exp_index',hecMESH_L%export_index(i)
!      print *,myrank,'exp_item ',hecMESH_L%export_item(hecMESH_L%export_index(i-1)+1:hecMESH_L%export_index(i))
    enddo
!   shared_index, shared_item
!    allocate(helpelem(hecMESH_L%ne_internal,hecMESH_L%n_neighbor_pe),stat=istat)
    allocate(helpelem(hecMESH_L%n_elem,hecMESH_L%n_neighbor_pe),stat=istat)
    helpelem(:,:) = 0
    allocate(helpCountElem(hecMESH_L%n_neighbor_pe),stat=istat)
    helpCountElem(:) = 0
!
    do i=1,hecMESH_L%n_elem
      n = 0; help(:) = 0
      do j=mak_loc%eptr(i),mak_loc%eptr(i+1)-1
        if(help(part(mak_loc%ngid(mak_loc%eind(j)))) /= 0) cycle
        n = n + 1
        help(part(mak_loc%ngid(mak_loc%eind(j)))) = 1
      enddo
      if(n == 1) cycle
      do j=1,hecMESH_L%PETOT
        if(help(j) == 0) cycle
        if(j == partID) cycle
        do k=1,hecMESH_L%n_neighbor_pe
          if(hecMESH_L%neighbor_pe(k) == j-1) goto 112
        enddo
        cycle
112     helpCountElem(k) = helpCountElem(k) + 1
        helpelem(helpCountElem(k),k) = i
      enddo
    enddo
    allocate(hecMESH_L%shared_index(0:hecMESH_L%n_neighbor_pe),stat=istat)
    hecMESH_L%shared_index(0) = 0
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%shared_index(i) = hecMESH_L%shared_index(i-1) + helpCountElem(i)
    enddo
    allocate(hecMESH_L%shared_item(hecMESH_L%shared_index(hecMESH_L%n_neighbor_pe)),stat=istat)
    hecMESH_L%shared_item(:) = 0
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%shared_item(hecMESH_L%shared_index(i-1)+1:hecMESH_L%shared_index(i)) = helpelem(1:helpCountElem(i),i)
!      print *,myrank,'shared_index',hecMESH_L%shared_index(i)
!      print *,myrank,'shared_item',hecMESH_L%shared_item(hecMESH_L%shared_index(i-1)+1:hecMESH_L%shared_index(i))
    enddo
    
!    write(myrank+30,*)'PE'
!    write(myrank+30,*)hecMESH_L%n_neighbor_pe
!    write(myrank+30,*)hecMESH_L%neighbor_pe(:)
!    write(myrank+30,*)'import'
!    write(myrank+30,'(10I10)')hecMESH_L%import_index(0:)
!    write(myrank+30,'(10I10)')hecMESH_L%import_item(:)
!    write(myrank+30,*)'export'
!    write(myrank+30,'(10I10)')hecMESH_L%export_index(0:)
!    write(myrank+30,'(10I10)')hecMESH_L%export_item(:)
!    write(myrank+30,*)'shared'
!    write(myrank+30,'(10I10)')hecMESH_L%shared_index(0:)
!    write(myrank+30,'(10I10)')hecMESH_L%shared_item(:)

!-- Material data
    call paraContact_copyHecmwMaterial(hecMESH_G%material,hecMESH_L%material)
    
!-- Section data
    call paraContact_copyHecmwSection(hecMESH_G%section,hecMESH_L%section)

!-- MPC data (is not supported)
    hecMESH_L%mpc%n_mpc = 0

!-- Amplitude data
    call paraContact_copyHecmwAmplitude(hecMESH_G%amp,hecMESH_L%amp)

!-- Node group data
!    maxNum = 0
!    do i=1,hecMESH_G%node_group%n_grp
!      maxNum = max(maxNum,hecMESH_G%node_group%grp_index(i)- hecMESH_G%node_group%grp_index(i-1))
!    enddo
!    allocate(indexNodeG2L_Hanging(hecMESH_G%n_node),stat=istat)
!    indexNodeG2L_Hanging(:) = indexNodeG2L(:,partID)
!    do i=1,nHangingNode
!      if(indexNodeG2L_Hanging(indexHangingNode(i)) /= 0) then
!        stop
!      else
!        indexNodeG2L_Hanging(indexHangingNode(i)) = hecMESH_L%nn_middle + i
!      endif
!    enddo
!    call paraContact_getHecmwLocalNodeGroup(hecMESH_G%node_group,hecMESH_L%node_group,indexNodeG2L_Hanging,maxNum)
!    deallocate(indexNodeG2L_Hanging,stat=istat)

!-- Element group data
    maxNum = 0
    do i=1,hecMESH_G%elem_group%n_grp
      maxNum = max(maxNum,hecMESH_G%elem_group%grp_index(i)- hecMESH_G%elem_group%grp_index(i-1))
    enddo
    call paraContact_getHecmwLocalElementGroup(hecMESH_G%elem_group,hecMESH_L%elem_group,indexElmtG2L(:,partID),maxNum)
    if(associated(indexHangingNode)) deallocate(indexHangingNode,stat=istat)
    
end subroutine paraContact_GetLocalMesh_all_new

subroutine rtri_GetLocalMesh_all(hecMESH_G,mak_loc,part,partID,hecMESH_L,indexNodeG2L,indexElmtG2L)
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_G
  type(MAKROSE_STRUCT),intent(inout)    ::  mak_loc
  integer(kint),intent(in)              ::  part(:)
  integer(kint),intent(in)              ::  partID    ! partition ID
  type(hecmwST_local_mesh),intent(inout)::  hecMESH_L
  integer(kint),intent(in)              ::  indexNodeG2L(:,:),indexElmtG2L(:,:)
!
  integer(kint) ::  i,j,n,istat,pid,nodeID,elemID,nodeID1,ii,jj,maxNum,k,l
  integer(kint),allocatable   ::  neighbor(:),help(:),helpnode(:,:),helpelem(:,:),helpCountNode(:),helpCountElem(:)
  integer(kint),allocatable   ::  markelem(:),marknode(:),temp(:),markPE(:)
  logical ::  sameflag
  integer ::  buf(2),rev(2),ierr
!
    call hecmw_nullify_mesh(hecMESH_L)
!-- General data
    hecMESH_L%gridfile              = hecMESH_G%gridfile
    hecMESH_L%hecmw_n_file          = hecMESH_G%hecmw_n_file
    allocate(hecMESH_L%files(size(hecMESH_G%files)),stat=istat)
    hecMESH_L%files                 = hecMESH_G%files
    hecMESH_L%header                = hecMESH_G%header
    hecMESH_L%hecmw_flag_adapt      = hecMESH_G%hecmw_flag_adapt
    hecMESH_L%hecmw_flag_initcon    = hecMESH_G%hecmw_flag_initcon
    
    hecMESH_L%hecmw_flag_parttype   = hecMESH_G%hecmw_flag_parttype
    hecMESH_L%hecmw_flag_partdepth  = hecMESH_G%hecmw_flag_partdepth
    hecMESH_L%hecmw_flag_version    = hecMESH_G%hecmw_flag_version
    hecMESH_L%zero_temp             = hecMESH_G%zero_temp
!-- Node data
!    hecMESH_L%n_node                = mak_loc%nn
!    hecMESH_L%n_node_gross          = hecMESH_L%n_node !hecMESH_G%n_node_gross    !?
!    hecMESH_L%nn_internal           = mak_loc%nn_i
!    hecMESH_L%n_dof                 = hecMESH_G%n_dof
!    hecMESH_L%n_dof_grp             = hecMESH_G%n_dof_grp
!    hecMESH_L%n_dof_tot             = hecMESH_G%n_dof_tot       !?
    
!-- Surface group data
    maxNum = 0
    do i=1,hecMESH_G%surf_group%n_grp
      maxNum = max(maxNum,hecMESH_G%surf_group%grp_index(i)- hecMESH_G%surf_group%grp_index(i-1))
    enddo
    call paraContact_getHecmwLocalSurfaceGroup(hecMESH_G%surf_group,hecMESH_L%surf_group, &
         indexElmtG2L(:,partID),maxNum,hecMESH_G,part,partID)

!-- Contact pair data
      call paraContact_copyHecmwContactPair(hecMESH_G%contact_pair,hecMESH_L%contact_pair)
      
!-- Node data
    hecMESH_L%n_node                = mak_loc%nn
    hecMESH_L%n_node_gross          = hecMESH_L%n_node !hecMESH_G%n_node_gross    !?
    hecMESH_L%nn_internal           = mak_loc%nn_i
    hecMESH_L%n_dof                 = hecMESH_G%n_dof
    hecMESH_L%n_dof_grp             = hecMESH_G%n_dof_grp
    hecMESH_L%n_dof_tot             = hecMESH_G%n_dof_tot       !?

    allocate(hecMESH_L%node(hecMESH_L%n_node*hecMESH_L%n_dof),stat=istat)
    do i=1,hecMESH_L%n_node
      hecMESH_L%node((i-1)*hecMESH_L%n_dof+1:i*hecMESH_L%n_dof) = mak_loc%x(1:hecMESH_L%n_dof,i)
    enddo
    allocate(hecMESH_L%node_ID(hecMESH_L%n_node*2),stat=istat)
    do i=1,hecMESH_L%n_node
      if(part(mak_loc%ngid(i)) == partID) then
        hecMESH_L%node_ID((i-1)*2+1)  = i
      else
        hecMESH_L%node_ID((i-1)*2+1)  = indexNodeG2L(mak_loc%ngid(i),part(mak_loc%ngid(i)))
      endif
      hecMESH_L%node_ID(i*2)        = part(mak_loc%ngid(i)) - 1
    enddo
    allocate(hecMESH_L%global_node_ID(hecMESH_L%n_node),stat=istat)
    hecMESH_L%global_node_ID(:) = mak_loc%ngid(:)
    allocate(hecMESH_L%node_dof_index(0:hecMESH_L%n_dof_grp),stat=istat)
    hecMESH_L%node_dof_index(0) = 0
    hecMESH_L%node_dof_index(1) = hecMESH_L%n_node
    allocate(hecMESH_L%node_dof_item(hecMESH_L%n_dof_grp),stat=istat)
    hecMESH_L%node_dof_item(:) = hecMESH_L%n_dof
    allocate(hecMESH_L%node_internal_list(hecMESH_L%nn_internal),stat=istat)
    n = 0
    do i=1,hecMESH_L%n_node
      if(part(mak_loc%ngid(i)) == partID) then
        n = n + 1
        hecMESH_L%node_internal_list(n) = i
      endif
    enddo
    
!    write(myrank+40,*)hecMESH_L%nn_internal,hecMESH_L%n_node
!    do i=1,hecMESH_L%n_node
!      write(myrank+40,*)i,hecMESH_L%node((i-1)*3+1:i*3)
!    enddo
!    do i=1,hecMESH_L%n_node
!      write(myrank+40,*)hecMESH_L%node_ID((i-1)*2+1:i*2)
!    enddo
!    write(myrank+40,'(10I6)')hecMESH_L%global_node_ID(:)
!
!-- Node group data
    maxNum = 0
    do i=1,hecMESH_G%node_group%n_grp
      maxNum = max(maxNum,hecMESH_G%node_group%grp_index(i)- hecMESH_G%node_group%grp_index(i-1))
    enddo
    call paraContact_getHecmwLocalNodeGroup(hecMESH_G%node_group,hecMESH_L%node_group,indexNodeG2L(:,partID),maxNum)
    
!-- Element data 
    hecMESH_L%n_elem                = mak_loc%ne
    hecMESH_L%n_elem_gross          = hecMESH_L%n_elem  !hecMESH_G%n_elem_gross      ! ?
    hecMESH_L%ne_internal           = mak_loc%ne_i
    call paraContact_GetLocalElementType(mak_loc,hecMESH_L%n_elem_type,hecMESH_L%elem_type_index,hecMESH_L%elem_type_item)
    allocate(hecMESH_L%elem_type(hecMESH_L%n_elem),stat=istat)
    do i=1,hecMESH_L%n_elem_type
      hecMESH_L%elem_type(hecMESH_L%elem_type_index(i-1)+1:hecMESH_L%elem_type_index(i)) = hecMESH_L%elem_type_item(i)
    enddo
    hecMESH_L%n_elem_mat_ID = hecMESH_L%n_elem
    allocate(hecMESH_L%elem_mat_ID_index(0:hecMESH_L%n_elem),stat=istat)
    do i=0,hecMESH_L%n_elem
      hecMESH_L%elem_mat_ID_index(i) = i
    enddo
    allocate(hecMESH_L%elem_mat_ID_item(hecMESH_L%n_elem),stat=istat)
    hecMESH_L%elem_mat_ID_item(:) = mak_loc%emat(:)
    allocate(hecMESH_L%elem_node_index(0:hecMESH_L%n_elem),stat=istat)
    do i=0,hecMESH_L%n_elem
      hecMESH_L%elem_node_index(i) = mak_loc%eptr(i+1)-1
    enddo
    allocate(hecMESH_L%elem_node_item(size(mak_loc%eind)),stat=istat)
    hecMESH_L%elem_node_item(:) = mak_loc%eind(:)
    allocate(hecMESH_L%global_elem_ID(hecMESH_L%n_elem),stat=istat)
    do i=1,hecMESH_L%n_elem
      hecMESH_L%global_elem_ID(i) = hecMESH_G%global_elem_ID(mak_loc%egid(i))
    enddo
    allocate(hecMESH_L%elem_ID(hecMESH_L%n_elem*2),stat=istat)
!    allocate(hecMESH_L%elem_internal_list(hecMESH_L%ne_internal),stat=istat)
    allocate(temp(hecMESH_L%n_elem),stat=istat)
    n = 0
    next_elem : do i=1,hecMESH_L%n_elem
!      hecMESH_L%elem_ID((i-1)*2+1) = i
!      pid = part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(hecMESH_L%elem_node_index(i-1)+1)))
!      if( pid /= partID) then
!        hecMESH_L%elem_ID(i*2) = 0
!      else
!        do j=hecMESH_L%elem_node_index(i-1)+1,hecMESH_L%elem_node_index(i)
!          if(pid /= part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(j)))) then
!            hecMESH_L%elem_ID(i*2) = 0
!            cycle next_elem
!          endif
!        enddo
!        hecMESH_L%elem_ID(i*2) = pid - 1
!        n = n + 1
!        hecMESH_L%elem_internal_list(n) = i
!      endif
       sameflag = .true.
       pid = part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(hecMESH_L%elem_node_index(i-1)+1)))
       do j=hecMESH_L%elem_node_index(i-1)+1,hecMESH_L%elem_node_index(i)
         if(pid /= part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(j)))) then
           if(sameflag) sameflag = .false.
           pid = min(pid,part(hecMESH_L%global_node_ID(hecMESH_L%elem_node_item(j))))
         endif
       enddo
!       if(sameflag) then
!         n = n + 1
!         hecMESH_L%elem_internal_list(n) = i
!       endif
       hecMESH_L%elem_ID(i*2) = pid - 1
       if(hecMESH_L%elem_ID(i*2) + 1 == partID) then
         n = n + 1
         temp(n) = i
         hecMESH_L%elem_ID((i-1)*2+1) = n
       else
         hecMESH_L%elem_ID((i-1)*2+1) = indexElmtG2L(mak_loc%egid(i),pid)
!         hecMESH_L%elem_internal_list(n) = i
       endif
       
    enddo next_elem
    hecMESH_L%ne_internal = n
    allocate(hecMESH_L%elem_internal_list(hecMESH_L%ne_internal),stat=istat)
    hecMESH_L%elem_internal_list(1:n) = temp(1:n)
    deallocate(temp,stat=istat)
!    print *,'my_rank',partID-1, hecMESH_L%ne_internal,hecMESH_L%n_elem,hecMESH_L%nn_internal,hecMESH_L%n_node
        
    allocate(hecMESH_L%section_ID(hecMESH_L%n_elem),stat=istat)
    do i=1,hecMESH_L%n_elem
      hecMESH_L%section_ID(i) = hecMESH_G%section_ID(mak_loc%egid(i))
    enddo
    
!    write(myrank+40,*)'elements'
!    write(myrank+40,*)hecMESH_L%ne_internal,hecMESH_L%n_elem
!    do i=1,hecMESH_L%n_elem
!      write(myrank+40,'(20I10)')i,hecMESH_L%elem_node_item(hecMESH_L%elem_node_index(i-1)+1:hecMESH_L%elem_node_index(i))
!    enddo
!    do i=1,hecMESH_L%n_elem
!      write(myrank+40,*)hecMESH_L%elem_ID((i-1)*2+1:i*2)
!    enddo
!    write(myrank+40,'(10I6)')hecMESH_L%global_elem_ID(:)

!-- PE data
    hecMESH_L%PETOT = hecMESH_G%PETOT 
    if(partID == 1) then
      hecMESH_L%zero = 1
    else
      hecMESH_L%zero = 0
    endif
    hecMESH_L%MPI_COMM = hecMESH_G%MPI_COMM
    hecMESH_L%my_rank = myrank
    hecMESH_L%n_subdomain = hecMESH_L%PETOT
    hecMESH_L%n_neighbor_pe = 0
    allocate(neighbor(hecMESH_L%PETOT-1),stat=istat)
    allocate(help(hecMESH_L%PETOT),stat=istat)
    allocate(helpnode(hecMESH_L%n_node - hecMESH_L%nn_internal,hecMESH_L%PETOT),stat=istat)
    helpnode(:,:) = 0
    help(:) = 0
!   Count Neighboring Processors
!   import_index, import_item
    do i=1,hecMESH_L%n_node
      pid = part(hecMESH_L%global_node_ID(i))
      if(pid /= partID) then
        help(pid) = help(pid) + 1
        helpnode(help(pid),pid) = i
      endif
    enddo
    do i=1,size(help)
      if(help(i) /= 0) hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
    enddo
    allocate(hecMESH_L%neighbor_pe(hecMESH_L%n_neighbor_pe),stat=istat)
    hecMESH_L%n_neighbor_pe = 0
    do i=1,size(help)
      if(help(i) /= 0) then
        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
        hecMESH_L%neighbor_pe(hecMESH_L%n_neighbor_pe) = i - 1  !start from 0
      endif
    enddo
    allocate(hecMESH_L%import_index(0:hecMESH_L%n_neighbor_pe),stat=istat)
    hecMESH_L%n_neighbor_pe = 0
    hecMESH_L%import_index(0) = 0
    do i=1,size(help)
      if(help(i) /= 0) then
        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
        hecMESH_L%import_index(hecMESH_L%n_neighbor_pe) = hecMESH_L%import_index(hecMESH_L%n_neighbor_pe-1) + help(i)
      endif
    enddo
    allocate(hecMESH_L%import_item(hecMESH_L%import_index(hecMESH_L%n_neighbor_pe)),stat=istat)   
    hecMESH_L%n_neighbor_pe = 0
    do i=1,size(help)
      if(help(i) /= 0) then
        hecMESH_L%n_neighbor_pe = hecMESH_L%n_neighbor_pe + 1
        hecMESH_L%import_item(hecMESH_L%import_index(hecMESH_L%n_neighbor_pe-1)+  &
                              1:hecMESH_L%import_index(hecMESH_L%n_neighbor_pe)) = helpnode(1:help(i),i)
      endif
    enddo
!!  export_index, export_item
    if(.not.(associated(mak_loc%nptr).and.associated(mak_loc%nind))) then
      call Mak_GetElementOnNode(mak_loc)
    endif
!
    deallocate(helpnode)
    allocate(helpnode(hecMESH_L%nn_internal,hecMESH_L%n_neighbor_pe),stat=istat)
!    print *,'myrank',myrank,size(helpnode,1),hecMESH_L%nn_internal
    helpnode(:,:) = 0
    allocate(helpCountNode(hecMESH_L%n_neighbor_pe),stat=istat)
    helpCountNode(:) = 0
    allocate(markPE(hecMESH_L%PETOT),stat=istat)
    
    do i=1,hecMESH_L%n_node
      if(hecMESH_L%node_ID(i*2) /= partID - 1) cycle
      markPE(:) = 0
      do j=mak_loc%nptr(i),mak_loc%nptr(i+1)-1
        elemID = mak_loc%nind(j)
        help(:) = 0
        n = 0
        do k=mak_loc%eptr(elemID),mak_loc%eptr(elemID+1)-1
          if(help(part(mak_loc%ngid(mak_loc%eind(k)))) == 0) then
            help(part(mak_loc%ngid(mak_loc%eind(k)))) = 1
            n = n + 1
          endif
        enddo
        if(n == 1) cycle
        do k=1,hecMESH_L%PETOT
          if(help(k) == 0) cycle
          if(k == partID) cycle
          if(markPE(k) /= 0) cycle
          markPE(k) = 1
          do l=1,hecMESH_L%n_neighbor_pe
            if(hecMESH_L%neighbor_pe(l) == k-1) exit
          enddo
          helpCountNode(l) = helpCountNode(l) + 1
          if( helpCountNode(l) > size(helpnode,1)) then
!            print *,helpCountNode(l),size(helpnode,1),hecMESH_L%nn_internal
          endif
          helpnode(helpCountNode(l),l) = i
        enddo
      enddo
    enddo
!
    allocate(hecMESH_L%export_index(0:hecMESH_L%n_neighbor_pe),stat=istat)
    hecMESH_L%export_index(0) = 0
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%export_index(i) = hecMESH_L%export_index(i-1) + helpCountNode(i)
    enddo
    allocate(hecMESH_L%export_item(hecMESH_L%export_index(hecMESH_L%n_neighbor_pe)),stat=istat)
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%export_item(hecMESH_L%export_index(i-1)+1:hecMESH_L%export_index(i)) = helpnode(1:helpCountNode(i),i)
    enddo
!   shared_index, shared_item
!    allocate(helpelem(hecMESH_L%ne_internal,hecMESH_L%n_neighbor_pe),stat=istat)
    allocate(helpelem(hecMESH_L%n_elem,hecMESH_L%n_neighbor_pe),stat=istat)
    helpelem(:,:) = 0
    allocate(helpCountElem(hecMESH_L%n_neighbor_pe),stat=istat)
    helpCountElem(:) = 0
!
    do i=1,hecMESH_L%n_elem
      n = 0; help(:) = 0
      do j=mak_loc%eptr(i),mak_loc%eptr(i+1)-1
        if(help(part(mak_loc%ngid(mak_loc%eind(j)))) /= 0) cycle
        n = n + 1
        help(part(mak_loc%ngid(mak_loc%eind(j)))) = 1
      enddo
      if(n == 1) cycle
      do j=1,hecMESH_L%PETOT
        if(help(j) == 0) cycle
        if(j == partID) cycle
        do k=1,hecMESH_L%n_neighbor_pe
          if(hecMESH_L%neighbor_pe(k) == j-1) exit
        enddo
        helpCountElem(k) = helpCountElem(k) + 1
        helpelem(helpCountElem(k),k) = i
      enddo
    enddo
    allocate(hecMESH_L%shared_index(0:hecMESH_L%n_neighbor_pe),stat=istat)
    hecMESH_L%shared_index(0) = 0
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%shared_index(i) = hecMESH_L%shared_index(i-1) + helpCountElem(i)
    enddo
    allocate(hecMESH_L%shared_item(hecMESH_L%shared_index(hecMESH_L%n_neighbor_pe)),stat=istat)
    do i=1,hecMESH_L%n_neighbor_pe
      hecMESH_L%shared_item(hecMESH_L%shared_index(i-1)+1:hecMESH_L%shared_index(i)) = helpelem(1:helpCountElem(i),i)
    enddo
    
!    write(myrank+40,*)'PE'
!    write(myrank+40,*)hecMESH_L%n_neighbor_pe
!    write(myrank+40,*)hecMESH_L%neighbor_pe(:)
!    write(myrank+40,*)'import'
!    write(myrank+40,'(10I10)')hecMESH_L%import_index(0:)
!    write(myrank+40,'(10I10)')hecMESH_L%import_item(:)
!    write(myrank+40,*)'import'
!    write(myrank+40,'(10I10)')hecMESH_L%export_index(0:)
!    write(myrank+40,'(10I10)')hecMESH_L%export_item(:)
!    write(myrank+40,*)'shared'
!    write(myrank+40,'(10I10)')hecMESH_L%shared_index(0:)
!    write(myrank+40,'(10I10)')hecMESH_L%shared_item(:)

!-- Material data
    call paraContact_copyHecmwMaterial(hecMESH_G%material,hecMESH_L%material)
    
!-- Section data
    call paraContact_copyHecmwSection(hecMESH_G%section,hecMESH_L%section)

!-- MPC data (is not supported)

!-- Amplitude data
    call paraContact_copyHecmwAmplitude(hecMESH_G%amp,hecMESH_L%amp)
    

!-- Node group data
!    maxNum = 0
!    do i=1,hecMESH_G%node_group%n_grp
!      maxNum = max(maxNum,hecMESH_G%node_group%grp_index(i)- hecMESH_G%node_group%grp_index(i-1))
!    enddo
!    call paraContact_getHecmwLocalNodeGroup(hecMESH_G%node_group,hecMESH_L%node_group,indexNodeG2L(:,partID),maxNum)

!-- Element group data
    maxNum = 0
    do i=1,hecMESH_G%elem_group%n_grp
      maxNum = max(maxNum,hecMESH_G%elem_group%grp_index(i)- hecMESH_G%elem_group%grp_index(i-1))
    enddo
    call paraContact_getHecmwLocalElementGroup(hecMESH_G%elem_group,hecMESH_L%elem_group,indexElmtG2L(:,partID),maxNum)
    
!-- Surface group data
!    maxNum = 0
!    do i=1,hecMESH_G%surf_group%n_grp
!      maxNum = max(maxNum,hecMESH_G%surf_group%grp_index(i)- hecMESH_G%surf_group%grp_index(i-1))
!    enddo
!    call paraContact_getHecmwLocalSurfaceGroup(hecMESH_G%surf_group,hecMESH_L%surf_group,indexElmtG2L(:,partID),maxNum)

!-- Contact pair data
!      call paraContact_copyHecmwContactPair(hecMESH_G%contact_pair,hecMESH_L%contact_pair)


end subroutine rtri_GetLocalMesh_all

subroutine paraContact_GetHangingSlaveNode(hecMESH_G,hecMesh_L,indexNodeG2L,indexHangingNode,nHangingNode,nOtherSlaveNode)
! Get hanging slave nodes with existing contact pair and surface group data
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_G
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_L
  integer(kint),intent(in)              ::  indexNodeG2L(:)
  integer(kint),pointer                 ::  indexHangingNode(:)
  integer(kint),intent(out)             ::  nHangingNode
  integer(kint),intent(out),optional    ::  nOtherSlaveNode
!
  integer   ::  istat
  integer(kint)   ::  i,j,surf_grp_ID,node_grp_ID
  integer(kint),allocatable ::  help(:),helpOther(:)
!
    nHangingNode = 0
    nOtherSlaveNode = 0
    if(associated(indexHangingNode)) deallocate(indexHangingNode,stat=istat)
    allocate(help(hecMESH_G%n_node),stat=istat)
    help(:) = 0
    allocate(helpOther(hecMESH_G%n_node),stat=istat)
    helpOther(:) = 0
!
    do i=1,hecMESH_L%contact_pair%n_pair
      surf_grp_ID = hecMESH_L%contact_pair%master_grp_id(i)
      if(hecMESH_L%surf_group%grp_index(surf_grp_ID) == hecMESH_L%surf_group%grp_index(surf_grp_ID-1)+1) cycle
!     This partation has contact pair master group
      node_grp_ID = hecMESH_L%contact_pair%slave_grp_id(i)
      do j=hecMESH_G%node_group%grp_index(node_grp_ID-1)+1,hecMESH_G%node_group%grp_index(node_grp_ID)
        if(indexNodeG2L(hecMESH_G%node_group%grp_item(j)) == 0.and.help(hecMESH_G%node_group%grp_item(j)) == 0) then
          help(hecMESH_G%node_group%grp_item(j)) = 1
          nHangingNode = nHangingNode + 1
        elseif(indexNodeG2L(hecMESH_G%node_group%grp_item(j)) /= 0.and.helpOther(hecMESH_G%node_group%grp_item(j)) == 0) then
          helpOther(hecMESH_G%node_group%grp_item(j)) = 1
          if(present(nOtherSlaveNode)) then
            nOtherSlaveNode = nOtherSlaveNode + 1
          endif
        endif
      enddo
    enddo
    deallocate(helpOther,stat=istat)
    allocate(indexHangingNode(nHangingNode),stat=istat)
    nHangingNode = 0
    help(:) = 0
    do i=1,hecMESH_L%contact_pair%n_pair
      surf_grp_ID = hecMESH_L%contact_pair%master_grp_id(i)
      if(hecMESH_L%surf_group%grp_index(surf_grp_ID) == hecMESH_L%surf_group%grp_index(surf_grp_ID-1)+1) cycle
!     This partation has contact pair master group
      node_grp_ID = hecMESH_L%contact_pair%slave_grp_id(i)
      do j=hecMESH_G%node_group%grp_index(node_grp_ID-1)+1,hecMESH_G%node_group%grp_index(node_grp_ID)
        if(indexNodeG2L(hecMESH_G%node_group%grp_item(j)) == 0.and.help(hecMESH_G%node_group%grp_item(j)) == 0) then
          help(hecMESH_G%node_group%grp_item(j)) = 1
          nHangingNode = nHangingNode + 1
          indexHangingNode(nHangingNode) = hecMESH_G%node_group%grp_item(j)
        endif
      enddo
    enddo
    deallocate(help,stat=istat)
end subroutine paraContact_GetHangingSlaveNode

function paraContact_isContactSlaveNode(hecMESH_G,hecMesh_L,part,nodeLID,exportSlavePE) result(yes)
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_G
  type(hecmwST_local_mesh),intent(in)   ::  hecMESH_L
  integer(kint),intent(in)              ::  part(:)
  integer(kint),intent(in)              ::  nodeLID   ! Local internal nodeID
  integer,intent(out)                   ::  exportSlavePE(:)
  logical   ::  yes
!
  integer(kint)   ::  i,j,surf_grp_ID,node_grp_ID,nodeGID
  integer(kint)   ::  ic,nsurf,ic_type,outtype,node_index(20),nodeID,nn,iss,m
!
    yes = .false.
    exportSlavePE(:) = 0
!
    nodeGID = hecMESH_L%global_node_ID(nodeLID)
    if(hecMESH_L%contact_pair%n_pair /= hecMESH_G%contact_pair%n_pair)  &
       stop 'hecMESH_L%contact_pair%n_pair /= hecMESH_G%contact_pair%n_pair'
    do i=1,hecMESH_L%contact_pair%n_pair
      surf_grp_ID = hecMESH_L%contact_pair%master_grp_id(i)
!     This partation has contact pair master group
      node_grp_ID = hecMESH_L%contact_pair%slave_grp_id(i)
      if(hecMESH_L%node_group%grp_index(node_grp_ID) == hecMESH_L%node_group%grp_index(node_grp_ID-1)+1) cycle
      do j=hecMESH_L%node_group%grp_index(node_grp_ID-1)+1,hecMESH_L%node_group%grp_index(node_grp_ID)
        if(hecMESH_L%node_group%grp_item(j) == nodeLID) then
          yes = .true.
          exit
        endif
      enddo
      if(.not.yes) cycle
      do j=hecMESH_G%surf_group%grp_index(surf_grp_ID-1)+1,hecMESH_G%surf_group%grp_index(surf_grp_ID)
      
        ic   = hecMESH_G%surf_group%grp_item(2*j-1)
        nsurf = hecMESH_G%surf_group%grp_item(2*j)
        ic_type = hecMESH_G%elem_type(ic)
        call getSubFace( ic_type, nsurf, outtype, node_index )  
        nn=getNumberOfNodes( outtype )
        iss = hecMESH_G%elem_node_index(ic-1)
!       Master surface nodes only
        do m=1, nn
          nodeID = hecMESH_G%elem_node_item(iss + node_index(m))
          if(part(nodeID) /= part(nodeGID).and.exportSlavePE(part(nodeID)) == 0) then
            exportSlavePE(part(nodeID)) = 1
          endif
        enddo
      enddo     
    enddo
    
end function paraContact_isContactSlaveNode

subroutine paraContact_copyHecmwMaterial(matIn,matOut)
  type(hecmwST_material),intent(in)   ::  matIn
  type(hecmwST_material),intent(inout)::  matOut
!
  integer   ::  istat
!
    call hecmw_nullify_material(matOut)
    matOut%n_mat          = matIn%n_mat
    matOut%n_mat_item     = matIn%n_mat_item
    matOut%n_mat_subitem  = matIn%n_mat_subitem
    matOut%n_mat_table    = matIn%n_mat_table
    allocate(matOut%mat_name(matOut%n_mat),stat=istat)
    matOut%mat_name(:) = matIn%mat_name(:)
    allocate(matOut%mat_item_index(0:matOut%n_mat),stat=istat)
    matOut%mat_item_index(:) = matIn%mat_item_index(:)
    allocate(matOut%mat_subitem_index(0:matOut%n_mat_item),stat=istat)
    matOut%mat_subitem_index(:) = matIn%mat_subitem_index(:)
    allocate(matOut%mat_table_index(0:matOut%n_mat_subitem),stat=istat)
    matOut%mat_table_index(:) = matIn%mat_table_index(:)
    allocate(matOut%mat_val(matOut%mat_table_index(matOut%n_mat_subitem)),stat=istat)
    matOut%mat_val(:) = matIn%mat_val(:)
    allocate(matOut%mat_temp(matOut%mat_table_index(matOut%n_mat_subitem)),stat=istat)
    matOut%mat_temp(:) = matIn%mat_temp(:)
end subroutine paraContact_copyHecmwMaterial

subroutine paraContact_copyHecmwSection(sectIn,sectOut)
  type(hecmwST_section),intent(in)    ::  sectIn
  type(hecmwST_section),intent(inout) ::  sectOut
!
  integer   ::  istat
!
    call hecmw_nullify_section(sectOut)
    sectOut%n_sect          = sectIn%n_sect
    if(sectOut%n_sect == 0) return
    allocate(sectOut%sect_type(sectOut%n_sect),stat=istat)
    sectOut%sect_type(:) = sectIn%sect_type(:)
    allocate(sectOut%sect_opt(sectOut%n_sect),stat=istat)
    sectOut%sect_opt(:) = sectIn%sect_opt(:)
    
    allocate(sectOut%sect_mat_ID_index(0:sectOut%n_sect),stat=istat)
    sectOut%sect_mat_ID_index(:) = sectIn%sect_mat_ID_index(:)
    allocate(sectOut%sect_mat_ID_item(sectOut%sect_mat_ID_index(sectOut%n_sect)),stat=istat)
    sectOut%sect_mat_ID_item(:) = sectIn%sect_mat_ID_item(:)

    if(associated(sectIn%sect_I_index)) then
      allocate(sectOut%sect_I_index(0:size(sectIn%sect_I_index)-1),stat=istat)
      sectOut%sect_I_index(:) = sectIn%sect_I_index(:)
    endif
    if(associated(sectIn%sect_I_item)) then
      allocate(sectOut%sect_I_item(size(sectIn%sect_I_item)),stat=istat)
      sectOut%sect_I_index(:) = sectIn%sect_I_index(:)
    endif
    if(associated(sectIn%sect_R_index)) then
      allocate(sectOut%sect_R_index(0:size(sectIn%sect_R_index)-1),stat=istat)
      sectOut%sect_R_index(:) = sectIn%sect_R_index(:)
    endif
    if(associated(sectIn%sect_R_item)) then
      allocate(sectOut%sect_R_item(size(sectIn%sect_R_item)),stat=istat)
      sectOut%sect_R_item(:) = sectIn%sect_R_item(:)
    endif
    if(associated(sectIn%sect_orien_ID)) then
      allocate(sectOut%sect_orien_ID(size(sectIn%sect_orien_ID)),stat=istat)
      sectOut%sect_orien_ID(:) = sectIn%sect_orien_ID(:)
    endif

end subroutine paraContact_copyHecmwSection

subroutine paraContact_copyHecmwAmplitude(ampIn,ampOut)
  type(hecmwST_amplitude),intent(in)    ::  ampIn
  type(hecmwST_amplitude),intent(inout) ::  ampOut
!
  integer   ::  istat
!
    call hecmw_nullify_amplitude(ampOut)
    ampOut%n_amp          = ampIn%n_amp
    if(ampOut%n_amp == 0) return
    allocate(ampOut%amp_name(ampOut%n_amp),stat=istat)
    ampOut%amp_name(:) = ampIn%amp_name(:)
    allocate(ampOut%amp_type_definition(ampOut%n_amp),stat=istat)
    ampOut%amp_type_definition(:) = ampIn%amp_type_definition(:)
    allocate(ampOut%amp_type_time(ampOut%n_amp),stat=istat)
    ampOut%amp_type_time(:) = ampIn%amp_type_time(:)
    allocate(ampOut%amp_type_value(ampOut%n_amp),stat=istat)
    ampOut%amp_type_value(:) = ampIn%amp_type_value(:)
    allocate(ampOut%amp_index(0:size(ampIn%amp_index)-1),stat=istat)
    ampOut%amp_index(:) = ampIn%amp_index(:)
    allocate(ampOut%amp_val(size(ampIn%amp_val)),stat=istat)
    ampOut%amp_val(:) = ampIn%amp_val(:)
    allocate(ampOut%amp_table(size(ampIn%amp_table)),stat=istat)
    ampOut%amp_table(:) = ampIn%amp_table(:)
end subroutine paraContact_copyHecmwAmplitude

subroutine paraContact_getHecmwLocalNodeGroup(ngrpIn,ngrpOut,indexNodeG2L,maxNodeInGroup)
  type(hecmwST_node_grp),intent(in)   ::  ngrpIn
  type(hecmwST_node_grp),intent(inout)::  ngrpOut
  integer(kint),intent(in)            ::  indexNodeG2L(:)
  integer(kint),intent(in)            ::  maxNodeInGroup
!
  integer   ::  istat,count,i,j,ngid
  integer(kint),allocatable ::  nodeCounter(:),nodeIndex(:,:)
!
    call hecmw_nullify_node_grp(ngrpOut)
    ngrpOut%n_grp = ngrpIn%n_grp
    if(ngrpOut%n_grp > 0) then
      allocate(ngrpOut%grp_name(ngrpOut%n_grp),stat=istat)
      ngrpOut%grp_name(:) = ngrpIn%grp_name(:)
      allocate(nodeCounter(ngrpOut%n_grp),stat=istat)
      nodeCounter(:) = 0
      allocate(nodeIndex(maxNodeInGroup,ngrpOut%n_grp),stat=istat)
      nodeIndex(:,:) = 0
      do i=1,ngrpIn%n_grp
        do j=ngrpIn%grp_index(i-1)+1,ngrpIn%grp_index(i)
          ngid = ngrpIn%grp_item(j)
          if(indexNodeG2L(ngid) /= 0) then
            nodeCounter(i) = nodeCounter(i) + 1
            nodeIndex(nodeCounter(i),i) = indexNodeG2L(ngid)
          endif
        enddo
      enddo
      allocate(ngrpOut%grp_index(0:ngrpOut%n_grp),stat=istat)
      ngrpOut%grp_index(0) = 0
      do i=1,ngrpOut%n_grp
        ngrpOut%grp_index(i) = ngrpOut%grp_index(i-1) + nodeCounter(i)
      enddo
      allocate(ngrpOut%grp_item(ngrpOut%grp_index(ngrpOut%n_grp)),stat=istat)
      count = 0
      do i=1,ngrpOut%n_grp
        do j=1,nodeCounter(i)
          count = count + 1
          ngrpOut%grp_item(count) = nodeIndex(j,i)
        enddo
      enddo      
    endif
    deallocate(nodeCounter,stat=istat)
    deallocate(nodeIndex,stat=istat)
    
    ngrpOut%n_bc = ngrpIn%n_bc
    if(ngrpOut%n_bc == 0) return
    print *,'node_group%bc_... are not implemented!'
!    pause
end subroutine paraContact_getHecmwLocalNodeGroup

subroutine paraContact_getHecmwLocalElementGroup(egrpIn,egrpOut,indexElmtG2L,maxElemInGroup)
  type(hecmwST_elem_grp),intent(in)   ::  egrpIn
  type(hecmwST_elem_grp),intent(inout)::  egrpOut
  integer(kint),intent(in)            ::  indexElmtG2L(:)
  integer(kint),intent(in)            ::  maxElemInGroup
!
  integer   ::  istat,count,i,j,egid
  integer(kint),allocatable ::  nodeCounter(:),nodeIndex(:,:)
!
    call hecmw_nullify_elem_grp(egrpOut)
    egrpOut%n_grp = egrpIn%n_grp
    if(egrpOut%n_grp > 0) then
      allocate(egrpOut%grp_name(egrpOut%n_grp),stat=istat)
      egrpOut%grp_name(:) = egrpIn%grp_name(:)
      allocate(nodeCounter(egrpOut%n_grp),stat=istat)
      nodeCounter(:) = 0
      allocate(nodeIndex(maxElemInGroup,egrpOut%n_grp),stat=istat)
      nodeIndex(:,:) = 0
      do i=1,egrpIn%n_grp
        do j=egrpIn%grp_index(i-1)+1,egrpIn%grp_index(i)
          egid = egrpIn%grp_item(j)
          if(indexElmtG2L(egid) /= 0) then
            nodeCounter(i) = nodeCounter(i) + 1
            nodeIndex(nodeCounter(i),i) = indexElmtG2L(egid)
          endif
        enddo
      enddo
      allocate(egrpOut%grp_index(0:egrpOut%n_grp),stat=istat)
      egrpOut%grp_index(0) = 0
      do i=1,egrpOut%n_grp
        egrpOut%grp_index(i) = egrpOut%grp_index(i-1) + nodeCounter(i)
      enddo
      allocate(egrpOut%grp_item(egrpOut%grp_index(egrpOut%n_grp)),stat=istat)
      count = 0
      do i=1,egrpOut%n_grp
        do j=1,nodeCounter(i)
          count = count + 1
          egrpOut%grp_item(count) = nodeIndex(j,i)
        enddo
      enddo      
    endif
    deallocate(nodeCounter,stat=istat)
    deallocate(nodeIndex,stat=istat)
    
    egrpOut%n_bc = egrpIn%n_bc
    if(egrpOut%n_bc == 0) return
    print *,'elem_group%bc_... are not implemented!'
!    pause
end subroutine paraContact_getHecmwLocalElementGroup

subroutine paraContact_getHecmwLocalSurfaceGroup(egrpIn,egrpOut,indexElmtG2L,maxSurfInGroup,hecMESH_G,part,partID)
  type(hecmwST_surf_grp),intent(in)   ::  egrpIn
  type(hecmwST_surf_grp),intent(inout)::  egrpOut
  integer(kint),intent(in)            ::  indexElmtG2L(:)
  integer(kint),intent(in)            ::  maxSurfInGroup
  type(hecmwST_local_mesh),intent(in) ::  hecMESH_G
  integer(kint),intent(in)            ::  part(:)
  integer(kint),intent(in)            ::  partID
!
  integer   ::  istat,count,i,j,egid,pid,k
  integer(kint),allocatable ::  nodeCounter(:),nodeIndex(:,:),surfIndex(:,:)
!
    call hecmw_nullify_surf_grp(egrpOut)
    egrpOut%n_grp = egrpIn%n_grp
    if(egrpOut%n_grp > 0) then
      allocate(egrpOut%grp_name(egrpOut%n_grp),stat=istat)
      egrpOut%grp_name(:) = egrpIn%grp_name(:)
      allocate(nodeCounter(egrpOut%n_grp),stat=istat)
      nodeCounter(:) = 0
      allocate(nodeIndex(maxSurfInGroup,egrpOut%n_grp),stat=istat)
      allocate(surfIndex(maxSurfInGroup,egrpOut%n_grp),stat=istat)
      nodeIndex(:,:) = 0
      do i=1,egrpIn%n_grp
        do j=egrpIn%grp_index(i-1)+1,egrpIn%grp_index(i)
          egid = egrpIn%grp_item(j*2-1)
          if(indexElmtG2L(egid) /= 0) then
!!           Check if element is internal elements
!            pid = part(hecMESH_G%elem_node_item(hecMESH_G%elem_node_index(egid-1)+1))
!            do k=hecMESH_G%elem_node_index(egid-1)+1,hecMESH_G%elem_node_index(egid)
!              if(pid /= part(hecMESH_G%elem_node_item(k))) then
!                pid = min(pid,part(hecMESH_G%elem_node_item(k)))
!              endif
!            enddo
!            if(pid == partID) then
              nodeCounter(i) = nodeCounter(i) + 1
              nodeIndex(nodeCounter(i),i) = indexElmtG2L(egid)
              surfIndex(nodeCounter(i),i) = egrpIn%grp_item(j*2)
!            endif
          endif
        enddo
      enddo
      allocate(egrpOut%grp_index(0:egrpOut%n_grp),stat=istat)
      egrpOut%grp_index(0) = 0
      do i=1,egrpOut%n_grp
        egrpOut%grp_index(i) = egrpOut%grp_index(i-1) + nodeCounter(i)
      enddo
      allocate(egrpOut%grp_item(2*egrpOut%grp_index(egrpOut%n_grp)),stat=istat)
      count = 0
      do i=1,egrpOut%n_grp
        do j=1,nodeCounter(i)
          count = count + 1
          egrpOut%grp_item(2*count-1) = nodeIndex(j,i)
          egrpOut%grp_item(2*count  ) = surfIndex(j,i)
        enddo
      enddo      
    endif
    deallocate(nodeCounter,stat=istat)
    deallocate(nodeIndex,stat=istat)
    deallocate(surfIndex,stat=istat)
    
    egrpOut%n_bc = egrpIn%n_bc
    if(egrpOut%n_bc == 0) return
    print *,'elem_group%bc_... are not implemented!'
!    pause
end subroutine paraContact_getHecmwLocalSurfaceGroup

subroutine paraContact_copyHecmwContactPair(contactIn,contactOut)
  type(hecmwST_contact_pair),intent(in)   ::  contactIn
  type(hecmwST_contact_pair),intent(inout)::  contactOut
!
  integer   ::  istat
!
    call hecmw_nullify_contact_pair(contactOut)
    contactOut%n_pair          = contactIn%n_pair
    allocate(contactOut%name(contactOut%n_pair),stat=istat)
    contactOut%name(:) = contactIn%name(:)
    allocate(contactOut%type(contactOut%n_pair),stat=istat)
    contactOut%type(:) = contactIn%type(:)
    allocate(contactOut%slave_grp_id(contactOut%n_pair),stat=istat)
    contactOut%slave_grp_id(:) = contactIn%slave_grp_id(:)
    allocate(contactOut%master_grp_id(contactOut%n_pair),stat=istat)
    contactOut%master_grp_id(:) = contactIn%master_grp_id(:)
end subroutine paraContact_copyHecmwContactPair

subroutine copyClearMatrix(hecMAT,conMAT)
  type(hecmwST_matrix),intent(in)     :: hecMAT
  type(hecmwST_matrix),intent(inout)  :: conMAT
!
  integer ::  ierr
    call hecmw_nullify_matrix( conMAT )

    conMAT%N = hecMAT%N
    conMAT%NP = hecMAT%NP
    conMAT%ndof = hecMAT%ndof
    allocate(conMAT%AL(size(hecMAT%AL)), stat=ierr) 
        if ( ierr /= 0 ) stop " Allocation error, conMAT%AL "
        conMAT%AL = 0.0D0  

    allocate(conMAT%AU(size(hecMAT%AU)), stat=ierr) 
        if ( ierr /= 0 ) stop " Allocation error, conMAT%AU "
        conMAT%AU = 0.0D0   
    
    allocate(conMAT%B(size(hecMAT%B)), stat=ierr)      
        conMAT%B = 0.0D0                                                              
    
    allocate(conMAT%X(size(hecMAT%X)), stat=ierr)    
        conMAT%X = 0.0D0

    allocate(conMAT%D(size(hecMAT%D)), stat=ierr)    
        conMAT%D = 0.0D0
end subroutine copyClearMatrix

! ----------------------------------------------------------------------------
        subroutine hecMAT_clear( hecMAT )
! ----------------------------------------------------------------------------
!  Purpose: clear hecMAT matrix 
!           Dec.2, 2009  YUAN Xi
!  Notes: 1. Should be a function of middleware !
!  This subroutine should completely rewritten next year
! ----------------------------------------------------------------------------
        type( hecmwST_matrix ) :: hecMAT
        hecMAT%AL= 0.d0
        hecMAT%AU= 0.d0
        hecMAT%X = 0.d0
        hecMAT%D = 0.d0
        call hecmw_cmat_clear( hecMAT%cmat )
        end subroutine hecMAT_clear

!Update external nodal data
subroutine paraContact_CreateExportImport(hecMESH)
  use hecmw_util
  type (hecmwST_local_mesh),intent(in)  ::  hecMESH
!
  integer(kint)   ::  i,j,pid,count,istat,ierr,inum
  integer(kint),allocatable ::  help(:)
  integer(kint),allocatable ::  expMap(:)
  integer(kint),allocatable ::  req_s(:),  req_r(:)
  integer(kint),allocatable ::  sta_s(:,:),sta_r(:,:)
!
    allocate(expMap(nprocs),stat=istat)
    expMap(:) = 0
    allocate(help(0:nprocs-1),stat=istat)
    
    help(:) = 0
    npe_import = 0
    do i=hecMESH%nn_internal+1,hecMESH%n_node
      pid = hecMESH%node_ID(i*2)
      if(help(pid) == 0) then
        npe_import = npe_import + 1
      endif
      help(pid) = help(pid) + 1
    enddo
    allocate(import_pe(npe_import),stat=istat)
    allocate(import_index(0:npe_import),stat=istat)
    allocate(import_item(hecMESH%n_node - hecMESH%nn_internal),stat=istat)
    
    import_index(0) = 0
    count = 0
    do i=0,nprocs-1
      if(help(i) /= 0) then
        count = count + 1
        import_pe(count) = i
        import_index(count) = import_index(count-1) + help(i)
      endif
    enddo
    help(:) = 0
    do i=hecMESH%nn_internal+1,hecMESH%n_node
      pid = hecMESH%node_ID(i*2)
      help(pid) = help(pid) + 1
      do j=1,npe_import
        if(import_pe(j) == pid) exit
      enddo
      import_item(import_index(j-1) + help(pid)) = hecMESH%node_ID(i*2-1)
    enddo
    call hecmw_BARRIER(hecMESH)
!
!   Send to PEs which contain the import nodes as internal nodes
    do i=0,nprocs-1
      call hecmw_gather_int_1(help(i), expMap, i, hecMESH%MPI_COMM)
!      call MPI_GATHER(help(i),1,MPI_INTEGER,expMap,1,MPI_INTEGER,i,hecMESH%MPI_COMM,ierr)
    enddo
!    print *,myrank,'expMap',expMap(:)
    call hecmw_BARRIER(hecMESH)
!
    npe_export = 0
    count = 0
    do i=1,nprocs
      if(expMap(i) > 0) then
        npe_export = npe_export + 1
        count = count + expMap(i)
      endif
    enddo
    if(npe_export > 0) then
      allocate(export_pe(npe_export),stat=istat)
      allocate(export_index(0:npe_export),stat=istat)
      allocate(export_item(count),stat=istat)
      export_index(0) = 0
      count = 0
      do i=1,nprocs
        if(expMap(i) > 0) then
          count = count + 1
          export_pe(count) = i - 1
          export_index(count) = export_index(count-1) + expMap(i)
        endif
      enddo
!      print *,myrank,'ext_idx',export_index(:)
    endif
    call hecmw_BARRIER(hecMESH)
!    
    allocate(req_s(npe_import),stat=istat)
    if(npe_export > 0) allocate(req_r(npe_export),stat=istat)
    allocate(sta_s(hecmw_status_size,npe_import),stat=istat)
    if(npe_export > 0) allocate(sta_r(hecmw_status_size,npe_export),stat=istat)
    
    do i=1,npe_import
      inum = import_index(i) - import_index(i-1)
      call hecmw_isend_int(import_item(import_index(i-1)+1),inum,import_pe(i),0,hecMESH%MPI_COMM,req_s(i))
!      call MPI_ISEND(import_item(import_index(i-1)+1),inum,MPI_INTEGER,import_pe(i),0,hecMESH%MPI_COMM,req_s(i),ierr)
    enddo
    
    if(npe_export > 0) then
      do i=1,npe_export
        inum = export_index(i) - export_index(i-1)
        call hecmw_irecv_int(export_item(export_index(i-1)+1),inum,export_pe(i),0,hecMESH%MPI_COMM,req_r(i))
!        call MPI_IRECV(export_item(export_index(i-1)+1),inum,MPI_INTEGER,export_pe(i),0,hecMESH%MPI_COMM,req_r(i),ierr)
      enddo
      call hecmw_waitall(npe_export,req_r,sta_r)
!      call MPI_WAITALL(npe_export,req_r,sta_r,ierr)
    endif
    call hecmw_waitall(npe_import,req_s,sta_s)
!    call MPI_WAITALL(npe_import,req_s,sta_s,ierr)
    
    call hecmw_BARRIER(hecMESH) 
!
    help(:) = 0
    do i=hecMESH%nn_internal+1,hecMESH%n_node
      pid = hecMESH%node_ID(i*2)
      help(pid) = help(pid) + 1
      do j=1,npe_import
        if(import_pe(j) == pid) exit
      enddo
      import_item(import_index(j-1) + help(pid)) = i
    enddo
!    write(myrank+40,*)hecMESH%nn_internal,hecMESH%n_node
!    write(myrank+40,*)'external node'
!    write(myrank+40,'(2I10)')hecMESH%node_ID(hecMESH%nn_internal*2+1:)
!    write(myrank+40,*)'import'
!    write(myrank+40,*)'import_pe'
!    write(myrank+40,'(10I10)')import_pe(:)
!    write(myrank+40,*)'import_index'
!    write(myrank+40,'(10I10)')import_index(:)
!    write(myrank+40,*)'import_item'
!    do i=1,npe_import
!      write(myrank+40,'(10I10)')import_item(import_index(i-1)+1:import_index(i))
!    enddo
!
!    write(myrank+40,*)'export'
!    if(npe_export > 0) then
!      write(myrank+40,*)'export_pe'
!      write(myrank+40,'(10I10)')export_pe(:)
!      write(myrank+40,*)'export_index'
!      write(myrank+40,'(10I10)')export_index(:)
!      write(myrank+40,*)'export_item'
!      do i=1,npe_export
!        write(myrank+40,'(10I10)')export_item(export_index(i-1)+1:export_index(i))
!      enddo
!    endif
!    
    call hecmw_BARRIER(hecMESH)
    if(allocated(help)) deallocate(help,stat=istat)
    if(allocated(expMap)) deallocate(expMap,stat=istat)
    if(allocated(req_s)) deallocate(req_s,stat=istat)
    if(allocated(req_r)) deallocate(req_r,stat=istat)
    if(allocated(sta_s)) deallocate(sta_s,stat=istat)
    if(allocated(sta_r)) deallocate(sta_r,stat=istat)
end subroutine paraContact_CreateExportImport
! 
subroutine paraContact_Update_3_R(hecMESH,val)
  type (hecmwST_local_mesh),intent(in)  ::  hecMESH
  real(kreal),intent(inout)             ::  val(3*hecMESH%n_node)
!
  integer(kint)   ::   n,istat
  real(kreal),allocatable ::  WS(:), WR(:) 
!
    n = hecMESH%n_node
    allocate(WS(3*n),stat=istat)
    allocate(WR(3*n),stat=istat)
    WS(:) = 0.0D0
    WR(:) = 0.0D0
!
    call paraContact_send_recv_33(n,WS,WR,val,hecMESH%MPI_COMM,myrank)
    deallocate(WS,stat=istat)
    deallocate(WR,stat=istat)
end subroutine paraContact_Update_3_R

subroutine paraContact_send_recv_33(n,WS,WR,X,SOLVER_COMM,my_rank)
  integer(kint),intent(in)  ::  n
  real(kreal),intent(inout) ::  WS(3*n)
  real(kreal),intent(inout) ::  WR(3*n)
  real(kreal),intent(inout) ::  X(3*n)
  integer(kind=kint ),intent(in)  ::  SOLVER_COMM
  integer(kind=kint ),intent(in)  ::  my_rank
!
  integer(kint)   ::  neib,istart,inum,k,ii,ierr
  integer(kint), dimension(:,:), allocatable :: sta1
  integer(kint), dimension(:,:), allocatable :: sta2
  integer(kint), dimension(:  ), allocatable :: req1
  integer(kint), dimension(:  ), allocatable :: req2  
!
!C
!C-- INIT.
    if(npe_export > 0) allocate (sta1(hecmw_status_size,npe_export),stat=ierr)
    allocate (sta2(hecmw_status_size,npe_import),stat=ierr)
    if(npe_export > 0) allocate (req1(npe_export),stat=ierr)
    allocate (req2(npe_import),stat=ierr)
       
!C
!C-- SEND
    do neib= 1, npe_export
      istart= export_index(neib-1)
      inum  = export_index(neib  ) - istart
      do k= istart+1, istart+inum
        if(export_item(k) > n) print *,myrank,neib,k,export_item(k)
        ii = 3*export_item(k)
        WS(3*(k-istart)-2)= X(ii-2)
        WS(3*(k-istart)-1)= X(ii-1)
        WS(3*(k-istart)  )= X(ii  )
      enddo
      call hecmw_isend_r(WS(1),3*inum,export_pe(neib),0,SOLVER_COMM,req1(neib))
!      call MPI_ISEND (WS(1), 3*inum,MPI_DOUBLE_PRECISION,    &
!                      export_pe(neib), 0, SOLVER_COMM, req1(neib), ierr)
    enddo

!C
!C-- RECEIVE
    do neib= 1, npe_import
      istart= import_index(neib-1)
      inum  = import_index(neib  ) - istart
      call hecmw_irecv_r(WR(3*istart+1),3*inum,import_pe(neib),0,SOLVER_COMM,req2(neib))
!      call MPI_IRECV (WR(3*istart+1), 3*inum, MPI_DOUBLE_PRECISION,   &
!                      import_pe(neib), 0, SOLVER_COMM, req2(neib), ierr)
    enddo

    call hecmw_waitall(npe_import, req2, sta2)
!    call MPI_WAITALL (npe_import, req2, sta2, ierr)
   
    do neib= 1, npe_import
      istart= import_index(neib-1)
      inum  = import_index(neib  ) - istart
      do k= istart+1, istart+inum
        ii = 3*import_item(k)
        X(ii-2)= WR(3*k-2)
        X(ii-1)= WR(3*k-1)
        X(ii  )= WR(3*k  )
      enddo
    enddo

    if(npe_export > 0) then
      call hecmw_waitall(npe_export, req1, sta1)
!      call MPI_WAITALL (npe_export, req1, sta1, ierr)
    endif
    if(allocated(sta1)) deallocate (sta1)
    if(allocated(sta2)) deallocate (sta2)
    if(allocated(req1)) deallocate (req1)
    if(allocated(req2)) deallocate (req2)
end subroutine paraContact_send_recv_33

function useEntireMesh() result(yes)
  logical ::  yes
!
  logical ::  log
  character(len=256)  ::  tx
  integer ::  istat
    yes = .false.
    inquire(file='hecmw_ctrl.dat',exist=log)
    if(.not.log) then
      stop 'Error: hecmw_ctrl.dat file does not exist!'
    else
      open(14,file='hecmw_ctrl.dat',status='old')
      do
        read(14,'(a)',iostat=istat)tx
        if(istat < 0) exit
        if(istat > 0) print *, 'Read Error!'
        if(tx == '') cycle
        if(tx(1:5) == '!MESH') then
          if(foundToken(tx,'!MESH').and.    &
             foundToken(tx,'fstrMSH').and.  &
             foundToken(tx,'HECMW-ENTIRE')) then
            yes = .true.
            return
          endif
        endif
      enddo
    endif
end function useEntireMesh

function foundToken(string,token) result(yes)
  character(*),intent(in) ::  string,token
  logical :: yes
!
  integer ::  i,j,ns,nt
    yes = .false.
    ns = len_trim(string)
    nt = len_trim(token)
    next_char : do i=1,ns
      if(string(i:i) /= token(1:1)) cycle
      do j=1,nt-1
        if(i+j > ns) cycle next_char
        if(string(i+j:i+j) /= token(j+1:j+1)) cycle next_char
      enddo
      yes = .true.
      return
    enddo next_char
end function foundToken
      
end module m_fstr_para_contact
