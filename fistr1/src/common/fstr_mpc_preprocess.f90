!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provides functions to calculate mpc coeff
module mMPCPreprocess

  use mContactDef
  use hecmw
  use m_fstr
  use m_utilities
  implicit none

  private
  public :: fstr_create_coeff_tiedcontact
  public :: resize_structures
  public :: resize_structures_static

  integer(kind=kint), parameter :: LOGLVL=1
  integer(kind=kint), parameter :: DEBUG=0

contains

!========================================================================
! Create MPC structure from CONTACT, INTERACTION=TIED.
!======================================================================== 
  !> create mpc setup
  subroutine fstr_create_coeff_tiedcontact( cstep, hecMESH, hecMAT, fstrSOLID, &
      &  infoCTChange, tied_method, dump_equation )
    integer(kind=kint), intent(in)         :: cstep      !< current step number
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
    type(fstr_info_contactChange), intent(inout):: infoCTChange   !<
    integer(kind=kint), intent(in)         :: tied_method  !< tiecontact processing method
    integer(kind=kint), intent(in)         :: dump_equation

    integer(kind=kint) :: i, j, grpid, n_tied_slave, n_tied_slave_total
    type(tMPCCond), allocatable :: mpcs_old(:), mpcs_new(:), mpcs_all(:)
    real(kind=kreal) :: T(6)

    T(1) = hecmw_Wtime()
    ! create original mpc coeff
    n_tied_slave_total = 0
    do i = 1, fstrSOLID%n_contacts
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
      if( fstrSOLID%contacts(i)%algtype /= CONTACTTIED ) cycle

      call create_coeff_tiedcontact( fstrSOLID%contacts(i), n_tied_slave )
      n_tied_slave_total = n_tied_slave_total + n_tied_slave
    enddo

    do i = 1, fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      call create_coeff_tiedcontact( fstrSOLID%embeds(i), n_tied_slave )
      n_tied_slave_total = n_tied_slave_total + n_tied_slave
    enddo

    ! extract mpc coeff of dof 1
    allocate(mpcs_old(n_tied_slave_total))
    n_tied_slave_total = 0
    do i = 1, fstrSOLID%n_contacts
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
      if( fstrSOLID%contacts(i)%algtype /= CONTACTTIED ) cycle

      call extract_coeff_tiedcontact( fstrSOLID%contacts(i), mpcs_old, n_tied_slave_total )

      ! disable contact for tied
      if( tied_method == ktMETHOD_MPC ) fstrSOLID%contacts(i)%group = -1
      do j=1,size(fstrSOLID%contacts(i)%states)
        fstrSOLID%contacts(i)%states(j)%state = CONTACTFREE
      enddo
    enddo

    do i = 1, fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      call extract_coeff_tiedcontact( fstrSOLID%embeds(i), mpcs_old, n_tied_slave_total )

      ! disable embed for tied
      if( tied_method == ktMETHOD_MPC ) fstrSOLID%embeds(i)%group = -1
      do j=1,size(fstrSOLID%embeds(i)%states)
        fstrSOLID%embeds(i)%states(j)%state = CONTACTFREE
      enddo
    enddo

    T(2) = hecmw_Wtime()
    if( myrank == 0 .and. LOGLVL > 0 ) then
      write(*,'(A)') "create_coeff_tiedcontact timelog"
      write(*,'(A20,f8.2)') "create mpc coeff", T(2)-T(1)
    endif

    ! update mpcs
    call get_newmpc( hecMESH, hecMAT, mpcs_old, mpcs_new, dump_equation )

    T(3) = hecmw_Wtime()
    ! create mpc
    if( tied_method == ktMETHOD_MPC ) then
      call set_hecmwST_mpc( hecMESH, mpcs_new )
    else if ( tied_method == ktMETHOD_CONTACT ) then
      call set_contact_structure( cstep, fstrSOLID, mpcs_new )
    endif

    T(4) = hecmw_Wtime()
    if( myrank == 0 .and. LOGLVL > 0 ) write(*,'(A20,f8.2)') "set data structure", T(4)-T(3)
  end subroutine fstr_create_coeff_tiedcontact

  !> create mpc coeff 
  subroutine create_coeff_tiedcontact( contact, n_tied_slave )
    type( tContact ), intent(inout)         :: contact      !< contact info
    integer(kind=kint), intent(out)         :: n_tied_slave

    integer(kind=kint) :: ctsurf, etype, nnode, ndLocal(21) !< contents of type tContact
    integer(kind=kint) :: i, j, idx, dof(l_max_elem_node+1)
    real(kind=kreal)   :: shapefunc(l_max_elem_node), coeff(l_max_elem_node+1)

    n_tied_slave = 0
    do i=1,size(contact%slave)
      if( contact%states(i)%state == CONTACTFREE ) cycle

      ctsurf = contact%states(i)%surface
      etype = contact%master(ctsurf)%etype
      nnode = size(contact%master(ctsurf)%nodes)
      ndLocal(1) = contact%slave(i)
      ndLocal(2:nnode+1) = contact%master(ctsurf)%nodes(1:nnode)

      call getShapeFunc( etype, contact%states(i)%lpos, shapefunc )
      coeff(1) = 1
      coeff(2:nnode+1) = -shapefunc(1:nnode)

      !cut-off small coeff
      idx = 1
      do j=2,nnode+1
        if( dabs(coeff(j)) < 1.d-6 ) cycle
        idx = idx + 1
        coeff(idx) = coeff(j)
        ndLocal(idx) = ndLocal(j)
      enddo

      if( ndLocal(1) == ndLocal(2) ) then
        contact%states(i)%state = CONTACTFREE
        cycle ! skip meaningless constraint
      endif

      do j=1,3
        dof(1:idx+1) = j
        call init_mpc_cond(contact%states(i)%mpc_cond(j),idx)
        call set_mpc_cond(contact%states(i)%mpc_cond(j), idx, &
          &  ndLocal(1:idx), dof(1:idx), coeff(1:idx))
        !call print_mpc_cond(103, contact%states(i)%mpc_cond(j))
      enddo

      n_tied_slave = n_tied_slave + 1
    enddo
  end subroutine 

  !> extract mpc coeff of dof 1 for singular value decomposition
  subroutine extract_coeff_tiedcontact( contact, mpcs_old, n_tied_slave_total )
    type( tContact ), intent(inout)            :: contact      !< contact info
    type(tMPCCond), allocatable, intent(inout) :: mpcs_old(:)
    integer(kind=kint), intent(inout)          :: n_tied_slave_total

    integer(kind=kint) :: ctsurf, etype, nnode, ndLocal(21) !< contents of type tContact
    integer(kind=kint) :: i, j, dof(l_max_elem_node+1)
    real(kind=kreal)   :: shapefunc(l_max_elem_node), coeff(l_max_elem_node+1)

    do i=1,size(contact%slave)
      if( contact%states(i)%state == CONTACTFREE ) cycle
      n_tied_slave_total = n_tied_slave_total + 1
      mpcs_old(n_tied_slave_total) = contact%states(i)%mpc_cond(1)
    enddo
  end subroutine

!========================================================================
! Gather global node id
!======================================================================== 
  subroutine gather_all_gnid( hecMESH, global_node_ID_all )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    integer(kind=kint), allocatable, intent(inout) :: global_node_ID_all(:)

    integer(kind=kint), allocatable :: displs_nid(:), recvcount_nid(:)
    integer(kind=kint) :: N_loc, N, offset_node, i
    integer(kind=kint), allocatable :: item_sendbuf(:), item_recvbuf(:)

    ! create offset of node id
    N_loc = hecMESH%nn_internal
    call create_offset( hecMESH, N_loc, recvcount_nid, displs_nid )
    offset_node = displs_nid(myrank)
    N = displs_nid(nprocs)

    ! initialize gather array
    if( myrank == 0 ) then !root process
      allocate(item_recvbuf(N))
    else
      allocate(item_recvbuf(1))
    endif
    item_recvbuf(:) = 0
    allocate(item_sendbuf(N_loc))

    do i=1,N_loc
      item_sendbuf(i) = hecMESH%global_node_ID(i)
    enddo

    call hecmw_gatherv_int(item_sendbuf, N_loc, item_recvbuf, recvcount_nid, &
         & displs_nid(0:nprocs-1), 0, hecMESH%MPI_COMM)

    ! only rank 0 processes the gathered data
    if( myrank == 0 ) then
      if(allocated(global_node_ID_all)) deallocate(global_node_ID_all)
      allocate(global_node_ID_all(N))
      global_node_ID_all = item_recvbuf
    endif

    deallocate(item_sendbuf, recvcount_nid, displs_nid)
  end subroutine

!========================================================================
! Remove the excessive constraints of MPC.
!======================================================================== 
  subroutine get_newmpc( hecMESH, hecMAT, mpcs_old, mpcs_new, dump_equation )
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    type(hecmwST_matrix), intent(in)         :: hecMAT
    type(tMPCCond), allocatable, intent(inout) :: mpcs_old(:)
    type(tMPCCond), allocatable, intent(inout) :: mpcs_new(:)
    integer(kind=kint), intent(in)         :: dump_equation

    type(tMPCCond), allocatable :: mpcs_all(:), mpcs_all_new(:), mpcs_new_dofremoved(:)
    integer(kind=kint), allocatable :: displs_nid(:)
    integer(kind=kint), allocatable :: global_node_ID_all(:)
    real(kind=kreal) :: T(6)

    T(1) = hecmw_Wtime()

    ! gather all mpcs_old to rank 0 process
    call gather_all_mpcs( hecMESH, mpcs_old, mpcs_all, displs_nid )

    ! dump original mpc as equation
    if( dump_equation == ktDUMP_ASIS ) then
      call gather_all_gnid( hecMESH, global_node_ID_all )
      if( myrank == 0 ) call print_full_mpc_conditions_3d( mpcs_all, global_node_ID_all )
    endif

    T(2) = hecmw_Wtime()

    ! get regular mpc conditions by singular value decomposition
    if( myrank == 0 ) call get_newmpc_by_svd( mpcs_all, mpcs_all_new )

    T(3) = hecmw_Wtime()

    ! dump updated mpc as equation
    if( dump_equation == ktDUMP_REGULARIZED ) then
      call gather_all_gnid( hecMESH, global_node_ID_all )
      if( myrank == 0 ) call print_full_mpc_conditions_3d( mpcs_all_new, global_node_ID_all )
    endif

    ! distribute mpc conditions(mpcs_new is still specified by global node id )
    call distribute_mpcs( hecMESH, mpcs_all_new, mpcs_new, displs_nid )

    T(4) = hecmw_Wtime()

    ! create new connectivity by dof remove
    call create_new_connectivity_by_slave_dofremove( hecMESH, hecMAT, displs_nid, mpcs_new, mpcs_new_dofremoved )

    T(5) = hecmw_Wtime()

    ! re-construct communication table
    call reconst_commtable( hecMESH, mpcs_new, displs_nid )

    T(6) = hecmw_Wtime()
    if( myrank == 0 .and. LOGLVL > 0 ) then
      write(*,'(A20,f8.2)') "gather_all_mpcs", T(2)-T(1)
      write(*,'(A20,f8.2)') "get_newmpc_by_svd", T(3)-T(2)
      write(*,'(A20,f8.2)') "distribute_mpcs", T(4)-T(3)
      write(*,'(A20,f8.2)') "create_new_connectivity", T(5)-T(4)
      write(*,'(A20,f8.2)') "reconst_commtable", T(6)-T(5)
    endif
  end subroutine

!========================================================================
! Gather and distribute MPC functions.
!======================================================================== 
  subroutine gather_all_mpcs( hecMESH, mpcs_old, mpcs_all, displs_nid )
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    type(tMPCCond), allocatable, intent(in)    :: mpcs_old(:)
    type(tMPCCond), allocatable, intent(inout) :: mpcs_all(:)
    integer(kind=kint), allocatable :: displs_nid(:)

    integer(kind=kint) ::  i,j,N,i0,N_loc,nndof,N_mpc_loc
    integer(kind=kint) :: offset_node, offset_mpc, offset_mpcitem, pid, lid, dof(1000)
    integer(kind=kint), allocatable :: recvcount_nid(:)
    integer(kind=kint), allocatable :: recvcount_mpc(:), recvcount_mpc_item(:)
    integer(kind=kint), allocatable :: displs_mpc(:), displs_mpc_item(:)
    integer(kind=kint), allocatable :: nitemall(:), itemall(:), nitem_sendbuf(:), item_sendbuf(:)
    real(kind=kreal), allocatable :: coeffall(:), coeff_sendbuf(:)
    integer(kind=kint) :: n_mpcall, n_mpcitemall, idx, local_nid, global_nid, tmprank, nitem

    ! create offset of node id
    call create_offset( hecMESH, hecMESH%nn_internal, recvcount_nid, displs_nid )
    offset_node = displs_nid(myrank)
    N = displs_nid(nprocs)

    N_mpc_loc = size(mpcs_old)

    ! total number of mpcs
    ! create offset of mpc
    call create_offset( hecMESH, N_mpc_loc, recvcount_mpc, displs_mpc )
    offset_mpc = displs_mpc(myrank)
    n_mpcall = displs_mpc(nprocs)

    ! create offset of mpc items
    N_loc = 0
    do i=1,N_mpc_loc
      N_loc = N_loc + mpcs_old(i)%nitem
    enddo
    call create_offset( hecMESH, N_loc, recvcount_mpc_item, displs_mpc_item )
    offset_mpcitem = displs_mpc_item(myrank)
    n_mpcitemall = displs_mpc_item(nprocs)

    ! initialize gather array
    if( myrank == 0 ) then !root process
      allocate(nitemall(n_mpcall))
      allocate(itemall(n_mpcitemall),coeffall(n_mpcitemall))
    else
      allocate(nitemall(1))
      allocate(itemall(1),coeffall(1))
    endif
    nitemall(:) = 0
    itemall(:) = 0
    coeffall(:) = 0.d0
    allocate(nitem_sendbuf(N_mpc_loc))
    allocate(item_sendbuf(N_loc),coeff_sendbuf(N_loc))

    idx = 0
    do i=1,N_mpc_loc
      nitem_sendbuf(i) = mpcs_old(i)%nitem
      do j=1,mpcs_old(i)%nitem
        local_nid = hecMESH%node_ID(2*mpcs_old(i)%pid(j)-1)
        tmprank = hecMESH%node_ID(2*mpcs_old(i)%pid(j))
        global_nid = displs_nid(tmprank) + local_nid

        item_sendbuf(idx+j) = global_nid
        coeff_sendbuf(idx+j) = mpcs_old(i)%coeff(j)
      enddo
      idx = idx + nitem_sendbuf(i)
    enddo

    call hecmw_gatherv_int(nitem_sendbuf, N_mpc_loc, nitemall, recvcount_mpc, &
         & displs_mpc(0:nprocs-1), 0, hecMESH%MPI_COMM)
    call hecmw_gatherv_int(item_sendbuf, N_loc, itemall, recvcount_mpc_item, &
         & displs_mpc_item(0:nprocs-1), 0, hecMESH%MPI_COMM)
    call hecmw_gatherv_real(coeff_sendbuf, N_loc, coeffall, recvcount_mpc_item, &
         & displs_mpc_item(0:nprocs-1), 0, hecMESH%MPI_COMM)

    if( myrank == 0 ) then
      allocate(mpcs_all(n_mpcall))

      idx = 0
      dof(:) = 1
      do i=1,n_mpcall
        nitem = nitemall(i)
        call init_mpc_cond(mpcs_all(i), nitem)
        call set_mpc_cond(mpcs_all(i), nitem, itemall(idx+1:idx+nitem), dof(1:nitem), coeffall(idx+1:idx+nitem))
        idx = idx + nitem
      enddo
    endif

  end subroutine

  subroutine create_offset( hecMESH, N_loc, recvcounts, displs )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    integer(kind=kint), intent(in) :: N_loc
    integer(kind=kint), allocatable, intent(inout) :: recvcounts(:)
    integer(kind=kint), allocatable, intent(inout) :: displs(:)

    integer(kind=kint) :: i

    allocate(recvcounts(nprocs))
    allocate(displs(0:nprocs))
    recvcounts(:) = 0
    displs(:) = 0
    displs(myrank+1) = N_loc
    call hecmw_allreduce_I(hecMESH, displs, nprocs+1, hecmw_sum)
    do i=1,nprocs
      recvcounts(i) = displs(i)
      displs(i) = displs(i-1) + displs(i)
    end do
  end subroutine

  subroutine distribute_mpcs( hecMESH, mpcs_all, mpcs_new, displs_nid )
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    type(tMPCCond), allocatable, intent(in) :: mpcs_all(:)
    type(tMPCCond), allocatable, intent(inout) :: mpcs_new(:)
    integer(kind=kint), allocatable, intent(in) :: displs_nid(:)

    integer(kind=kint) ::  i,j,pid,nmpc_loc,nmpc_item_loc,count,nitem_tmp,ierr
    integer(kind=kint) :: iS, iE, offset_mpc, offset_mpcitem, nid, dof(1000), mpcid
    type(hecmwST_varray_int), allocatable :: mpclist(:)
    integer(kind=kint), allocatable :: pidlist_nid(:), item_offset(:)
    real(kind=kreal), allocatable :: sendbuf_size_nitem(:), sendbuf_size_item(:)
    integer(kind=kint), allocatable :: nitem_recvbuf(:), item_recvbuf(:)
    
    integer(kind=kint), allocatable :: sendcount_mpc(:), sendcount_mpc_item(:)
    integer(kind=kint), allocatable :: displs_mpc(:), displs_mpc_item(:)
    integer(kind=kint), allocatable :: nitemall(:), itemall(:), nitem_sendbuf(:), item_sendbuf(:)
    real(kind=kreal), allocatable :: coeffall(:), coeff_sendbuf(:), coeff_recvbuf(:)
    integer(kind=kint) :: n_mpcall, n_mpcitemall, idx, local_nid, global_nid, tmprank, nitem

    if( myrank == 0 ) then
      call HECMW_varray_int_initialize_all(mpclist,nprocs,1000)

      ! create node id - process id list
      allocate(pidlist_nid(displs_nid(nprocs)))
      do i=1,nprocs
        pidlist_nid(displs_nid(i-1)+1:displs_nid(i)) = i
      enddo
      
      ! get mpcs in each rank
      count = 0
      do i=1,size(mpcs_all)
        do j=1,mpcs_all(i)%nitem
          nid = mpcs_all(i)%pid(j)
          pid = pidlist_nid(nid)
          call HECMW_varray_int_add( mpclist(pid), i )
        enddo
      enddo
      do i=1,nprocs
        call HECMW_varray_int_unique_items(mpclist(i))
      enddo
      !! get offset of mpc items
      allocate(item_offset(0:size(mpcs_all)))
      item_offset(:) = 0
      do i=1,size(mpcs_all)
        item_offset(i) = item_offset(i-1) + mpcs_all(i)%nitem
      enddo
  
      deallocate(pidlist_nid)
    endif

    ! 1. scatter nitems
    allocate(displs_mpc(0:nprocs),sendcount_mpc(nprocs))
    displs_mpc(:) = 0
    sendcount_mpc(:) = 0
    if( myrank == 0 ) then
      !! get sendbuff size
      n_mpcall = 0
      do i=1,nprocs
        n_mpcall = n_mpcall + HECMW_varray_int_get_nitem(mpclist(i))
      enddo
      !! set nitem sendbuff
      allocate(nitem_sendbuf(n_mpcall))
      count = 0
      do i=1,nprocs
        do j=1,HECMW_varray_int_get_nitem(mpclist(i))
          mpcid = HECMW_varray_int_get_item(mpclist(i),j)
          count = count + 1
          nitem_sendbuf(count) = mpcs_all(mpcid)%nitem
        enddo
        displs_mpc(i) = count
        sendcount_mpc(i) = displs_mpc(i) - displs_mpc(i-1)
      enddo
    else
      allocate(nitem_sendbuf(1))
    endif
    call hecmw_allreduce_I(hecMESH, sendcount_mpc, nprocs, hecmw_sum)
    nmpc_loc = sendcount_mpc(myrank+1)
    allocate(nitem_recvbuf(nmpc_loc))

    call hecmw_scatterv_int(nitem_sendbuf, sendcount_mpc, displs_mpc(0:nprocs-1), &
          &     nitem_recvbuf, nmpc_loc, 0, hecMESH%MPI_COMM)
    deallocate(nitem_sendbuf)

    ! 2. scatter pids
    allocate(displs_mpc_item(0:nprocs),sendcount_mpc_item(nprocs))
    displs_mpc_item(:) = 0
    sendcount_mpc_item(:) = 0
    if( myrank == 0 ) then
      !! get sendbuff size
      count = 0
      do i=1,nprocs
        do j=1,HECMW_varray_int_get_nitem(mpclist(i))
          mpcid = HECMW_varray_int_get_item(mpclist(i),j)
          count = count + mpcs_all(mpcid)%nitem
        enddo
        displs_mpc_item(i) = count
        sendcount_mpc_item(i) = displs_mpc_item(i) - displs_mpc_item(i-1)
      enddo
      n_mpcitemall = displs_mpc_item(nprocs)

      !! set item sendbuff
      allocate(item_sendbuf(n_mpcitemall))
      count = 0
      do i=1,nprocs
        do j=1,HECMW_varray_int_get_nitem(mpclist(i))
          mpcid = HECMW_varray_int_get_item(mpclist(i),j)
          nitem_tmp = mpcs_all(mpcid)%nitem
          item_sendbuf(count+1:count+nitem_tmp) = mpcs_all(mpcid)%pid(1:nitem_tmp)
          count = count + nitem_tmp
        enddo
      enddo
    else
      allocate(item_sendbuf(1))
    endif
    call hecmw_allreduce_I(hecMESH, sendcount_mpc_item, nprocs, hecmw_sum)
    nmpc_item_loc = sendcount_mpc_item(myrank+1)
    allocate(item_recvbuf(nmpc_item_loc))
    call hecmw_scatterv_int(item_sendbuf, sendcount_mpc_item, displs_mpc_item(0:nprocs-1), &
          &     item_recvbuf, nmpc_item_loc, 0, hecMESH%MPI_COMM)
    deallocate(item_sendbuf)

    ! 3. scatter coeffs
    if( myrank == 0 ) then
      !! set item sendbuff
      allocate(coeff_sendbuf(n_mpcitemall))
      count = 0
      do i=1,nprocs
        do j=1,HECMW_varray_int_get_nitem(mpclist(i))
          mpcid = HECMW_varray_int_get_item(mpclist(i),j)
          nitem_tmp = mpcs_all(mpcid)%nitem
          coeff_sendbuf(count+1:count+nitem_tmp) = mpcs_all(mpcid)%coeff(1:nitem_tmp)
          count = count + nitem_tmp
        enddo
      enddo
    else
      allocate(coeff_sendbuf(1))
    endif
    allocate(coeff_recvbuf(nmpc_item_loc))
    call hecmw_scatterv_real(coeff_sendbuf, sendcount_mpc_item, displs_mpc_item(0:nprocs-1), &
          &     coeff_recvbuf, nmpc_item_loc, 0, hecMESH%MPI_COMM)
    deallocate(coeff_sendbuf)

    !set up local mpc
    allocate(mpcs_new(nmpc_loc))
    dof(:) = 1
    count = 0
    do i=1,nmpc_loc
      nitem = nitem_recvbuf(i)
      call init_mpc_cond(mpcs_new(i), nitem)
      call set_mpc_cond(mpcs_new(i), nitem, item_recvbuf(count+1:count+nitem), &
         &  dof(1:nitem), coeff_recvbuf(count+1:count+nitem))
      count = count + nitem
    enddo

  end subroutine

!========================================================================
! Update connectivity.
!======================================================================== 
  subroutine create_new_connectivity_by_slave_dofremove( hecMESH, hecMAT, displs_nid, mpcs_new, mpcs_new_dofremoved )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(hecmwST_matrix), intent(in)         :: hecMAT
    integer(kind=kint), allocatable :: displs_nid(:)
    type(tMPCCond), allocatable, intent(in) :: mpcs_new(:)
    type(tMPCCond), allocatable, intent(inout) :: mpcs_new_dofremoved(:)

    integer(kind=kint) :: i, j, iS, iE, pid, lid, iiS, iiE, offset, global_nid
    integer(kind=kint) :: tmplid, tmppid, nslave, nslave_total
    integer(kind=kint), allocatable :: ncon, slavecon_index(:), slavecon_item1(:), slavecon_item2(:) 

    !create local slave connectivity
    iS = displs_nid(myrank)
    iE = displs_nid(myrank+1)

    !! count number of slave nodes and connectivity
    nslave = 0
    do i=1,size(mpcs_new)
      pid = mpcs_new(i)%pid(1) !slave node
      if( pid <= iS ) cycle
      if( pid > iE ) cycle
      lid = pid-iS
      nslave = nslave + 1
    enddo

    !! create slavecon_index
    allocate(slavecon_index(0:nprocs))
    slavecon_index(0:nprocs) = 0
    slavecon_index(myrank+1) = nslave
    call hecmw_allreduce_I(hecMESH, slavecon_index, nprocs+1, hecmw_sum)
    do i=1,nprocs
      slavecon_index(i) = slavecon_index(i) + slavecon_index(i-1)
    enddo
    nslave_total = slavecon_index(nprocs)

    !! create slavecon_item1: index of slave con node id
    allocate(slavecon_item1(0:nslave_total))
    slavecon_item1(0:nslave_total) = 0
    offset = slavecon_index(myrank)
    nslave = 0
    do i=1,size(mpcs_new)
      pid = mpcs_new(i)%pid(1) !slave node
      if( pid <= iS ) cycle
      if( pid > iE ) cycle
      lid = pid-iS
      nslave = nslave + 1
      ncon = 1 !diag
      ncon = ncon + hecMAT%indexU(lid  )-hecMAT%indexU(lid-1) !upper
      ncon = ncon + hecMAT%indexL(lid  )-hecMAT%indexL(lid-1) !lower
      slavecon_item1(offset+nslave) = ncon
    enddo
    nslave = slavecon_index(nprocs)
    call hecmw_allreduce_I(hecMESH, slavecon_item1, nslave_total+1, hecmw_sum)
    do i=1,nslave_total
      slavecon_item1(i) = slavecon_item1(i) + slavecon_item1(i-1)
    enddo
    ncon = slavecon_item1(nslave_total)

    !! create slavecon_item2: index of slave con node id
    allocate(slavecon_item2(ncon))
    slavecon_item2(1:ncon) = 0
    nslave = 0
    ncon = slavecon_item1(slavecon_index(myrank))
    do i=1,size(mpcs_new)
      pid = mpcs_new(i)%pid(1) !slave node
      if( pid <= iS ) cycle
      if( pid > iE ) cycle
      lid = pid-iS
      !diag
      ncon = ncon + 1 
      slavecon_item2(ncon) = pid
      !lower
      do j=hecMAT%indexL(lid-1)+1,hecMAT%indexL(lid  ) !upper
        tmplid = hecMAT%itemL(j) !local id in myrank
        tmppid = hecMESH%node_ID(tmplid*2)  !process id of internal rank
        tmplid = hecMESH%node_ID(tmplid*2-1) !local id in internal rank
        global_nid = displs_nid(tmppid) + tmplid
        ncon = ncon + 1
        slavecon_item2(ncon) = global_nid
      enddo
      !upper
      do j=hecMAT%indexU(lid-1)+1,hecMAT%indexU(lid  ) !upper
        tmplid = hecMAT%itemU(j) !local id in myrank
        tmppid = hecMESH%node_ID(tmplid*2)  !process id of internal rank
        tmplid = hecMESH%node_ID(tmplid*2-1) !local id in internal rank
        global_nid = displs_nid(tmppid) + tmplid
        ncon = ncon + 1
        slavecon_item2(ncon) = global_nid
      enddo
    enddo
    ncon = slavecon_item1(nslave_total)
    call hecmw_allreduce_I(hecMESH, slavecon_item2, ncon, hecmw_sum)
  end subroutine

!========================================================================
! Update Communication table.
!======================================================================== 
  subroutine reconst_commtable( hecMESH, mpcs_new, displs_nid )
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    type(tMPCCond), allocatable, intent(inout) :: mpcs_new(:)
    integer(kind=kint), allocatable, intent(in) :: displs_nid(:)

    integer(kind=kint) :: i, j, nid, pid, lid, flag, count, nitem, iS, iE
    type(hecmwST_varray_int) :: external_new
    type(hecmwST_varray_int), allocatable :: import_gids(:)
    integer(kind=kint), allocatable :: requests(:), statuses(:,:)
    integer(kind=kint) :: n_newnodes, n_neighbor_pe_import, n_neighbor_pe_export, n_neighbor_pe_new

    integer(kind=kint), allocatable :: pids(:) !< 0:internal, 1:external(already defined), 2:external(newly defined by mpc)
    integer(kind=kint), allocatable :: node_ID_new(:), neighbor_pe_list(:), import_index_new(:), import_item_new(:)
    integer(kind=kint), allocatable :: export_index_new(:), export_item_new(:), send_index(:), recv_index(:)

    if( nprocs == 1 ) return

    allocate(pids(displs_nid(nprocs)))
    pids(:) = -1
    do i=1,hecMESH%n_node
      ! set internal and already defined external node flag
      pid = hecMESH%node_ID(i*2)
      lid = hecMESH%node_ID(i*2-1)
      nid = displs_nid(pid) + lid
      pids(nid) = pid
    enddo

    ! newly defined external node flag ad nprocs
    do i=1,size(mpcs_new)
      do j=1,mpcs_new(i)%nitem
        nid = mpcs_new(i)%pid(j)
        if( pids(nid) > -1 ) cycle
        pids(nid) = nprocs
      enddo
    enddo

    ! init additional nodes
    call HECMW_varray_int_initialize( external_new, 100 )
    do nid=1,displs_nid(nprocs)
      if( pids(nid) == nprocs ) call HECMW_varray_int_add( external_new, nid )
    enddo
    n_newnodes = HECMW_varray_int_get_nitem( external_new )

    ! append new import nodes to node_ID_new
    allocate( node_ID_new(2*(hecMESH%n_node+n_newnodes)) )
    node_ID_new(1:2*hecMESH%n_node) = hecMESH%node_ID(1:2*hecMESH%n_node)
    do i=1,n_newnodes
      nid = HECMW_varray_int_get_item( external_new, i )
      call get_lid_and_pid_from_gid( nid, displs_nid, pid, lid )
      node_ID_new(2*hecMESH%n_node+i*2  ) = pid
      node_ID_new(2*hecMESH%n_node+i*2-1) = lid
    enddo

    !!! create import_index_new
    allocate(import_index_new(0:nprocs))
    import_index_new(:) = 0
    do i=hecMESH%nn_internal+1,hecMESH%n_node+n_newnodes
      pid = node_ID_new(i*2)
      import_index_new(pid+1) = import_index_new(pid+1) + 1
    enddo

    !!! create export_index_new
    allocate(export_index_new(0:nprocs))
    export_index_new(:) = 0
    allocate(requests(nprocs-1), statuses(HECMW_STATUS_SIZE, nprocs-1))
    requests = 0
    statuses = 0
    count = 0
    do pid=0,nprocs-1
      if( pid == myrank ) cycle
      count = count + 1
      call hecmw_isend_int(import_index_new(pid+1:pid+1), 1, pid, 0, hecmw_comm_get_comm(), requests(count))
      !write(*,*) "send",pid, requests(pid+1)
    enddo
    call HECMW_Waitall(count, requests, statuses)
    count = 0
    do pid=0,nprocs-1
      if( pid == myrank ) cycle
      count = count + 1
      call hecmw_recv_int(export_index_new(pid+1:pid+1), 1, pid, 0, hecmw_comm_get_comm(), statuses(:, count))
    enddo

    !!! create n_neighbor_pe_new
    n_neighbor_pe_new = 0
    do i=1,nprocs
      pid = i-1
      if( export_index_new(i)+import_index_new(i) == 0 ) cycle
      n_neighbor_pe_new = n_neighbor_pe_new + 1
    enddo
    allocate(neighbor_pe_list(n_neighbor_pe_new))
    n_neighbor_pe_new = 0
    do i=1,nprocs
      pid = i-1
      if( export_index_new(i)+import_index_new(i) == 0 ) cycle
      n_neighbor_pe_new = n_neighbor_pe_new + 1
      neighbor_pe_list(n_neighbor_pe_new) = pid
      import_index_new(n_neighbor_pe_new) = import_index_new(n_neighbor_pe_new-1) + import_index_new(i)
      export_index_new(n_neighbor_pe_new) = export_index_new(n_neighbor_pe_new-1) + export_index_new(i)
    enddo

    !! create import_item_new
    call HECMW_varray_int_initialize_all(import_gids,nprocs,100)
    !!! already exist import nodes
    do i=hecMESH%nn_internal+1,hecMESH%n_node
      pid = node_ID_new(i*2)
      lid = hecMESH%node_ID(i*2-1)
      nid = displs_nid(pid) + lid
      call HECMW_varray_int_add( import_gids(pid+1), nid )
    enddo
    !!! newly added import nodes
    do i=1,n_newnodes
      nid = HECMW_varray_int_get_item( external_new, i )
      call get_lid_and_pid_from_gid( nid, displs_nid, pid, lid )
      call HECMW_varray_int_add( import_gids(pid+1), nid )
    enddo
    !!! allocate import_item_new
    allocate(import_item_new(import_index_new(n_neighbor_pe_new)))
    count = 0
    do i=1,nprocs
      nitem = HECMW_varray_int_get_nitem(import_gids(i))
      if( nitem == 0 ) cycle
      call HECMW_varray_int_get_item_all( import_gids(i), import_item_new(count+1:count+nitem))
      count = count + nitem
    enddo
    if( count /= import_index_new(n_neighbor_pe_new) ) &
      & write(*,*) "count,import_index_new(n_neighbor_pe_import)",count,import_index_new(n_neighbor_pe_import)

    !! create export_item_new
    allocate(export_item_new(export_index_new(n_neighbor_pe_new)))
    requests = 0
    statuses = 0
    count = 0
    do i=1,n_neighbor_pe_new
      pid = neighbor_pe_list(i)
      iS = import_index_new(i-1)+1
      iE = import_index_new(i)
      if( iS == iE ) cycle
      count = count + 1
      call hecmw_isend_int(import_item_new(iS:iE), iE-iS+1, pid, 0, hecmw_comm_get_comm(), requests(count))
    enddo
    call HECMW_Waitall(count, requests, statuses)
    count = 0
    do i=1,n_neighbor_pe_new
      pid = neighbor_pe_list(i)
      iS = export_index_new(i-1)+1
      iE = export_index_new(i)
      if( iS == iE ) cycle
      count = count + 1
      call hecmw_recv_int(export_item_new(iS:iE), iE-iS+1, pid, 0, hecmw_comm_get_comm(), statuses(:,count))
    enddo

    !!! change item id from global to local
    pids(:) = -1
    do i=1,hecMESH%n_node
      pid = node_ID_new(i*2)
      lid = node_ID_new(i*2-1)
      nid = displs_nid(pid) + lid
      if( nid > size(pids) ) write(*,*) "oversize:",nid,size(pids)
      pids(nid) = i
    enddo
    do i=1,n_newnodes
      pid = node_ID_new(2*hecMESH%n_node+i*2)
      nid = HECMW_varray_int_get_item( external_new, i )
      if( nid > size(pids) ) write(*,*) "oversize:",nid,size(pids)
      pids(nid) = hecMESH%n_node+i
    enddo
    do i=1,import_index_new(n_neighbor_pe_new)
      if( import_item_new(i) < 0 ) write(*,*) "Error: import_item_new(i) < 0 at",i
      import_item_new(i) = pids(import_item_new(i))
    enddo
    do i=1,export_index_new(n_neighbor_pe_new)
      if( export_item_new(i) < 0 ) write(*,*) "Error: export_item_new(i) < 0 at",i
      export_item_new(i) = pids(export_item_new(i))
    enddo
    do i=1,size(mpcs_new)
      do j=1,mpcs_new(i)%nitem
        nid = mpcs_new(i)%pid(j)
        if( pids(nid) == -1 ) write(*,*) "Error: pids(mpcs_new(i)%pid(j)) < 0 at",i,j 
        mpcs_new(i)%pid(j) = pids(nid)
      enddo
    enddo

    call update_hecMESH_commtable( hecMESH, n_newnodes, n_neighbor_pe_new, neighbor_pe_list, &
      & import_index_new, import_item_new, export_index_new, export_item_new, node_ID_new )

  end subroutine

  subroutine get_lid_and_pid_from_gid( gid, displs, pid, lid )
    integer(kind=kint), intent(in) :: gid
    integer(kind=kint), allocatable, intent(in) :: displs(:)
    integer(kind=kint), intent(out) :: pid
    integer(kind=kint), intent(out) :: lid

    integer(kind=kint) :: i
    do i=1,nprocs
      if( gid < displs(i-1)+1 ) cycle
      if( gid > displs(i)     ) cycle
      pid = i-1
      lid = gid - displs(i-1)
      exit
    enddo
  end subroutine

  subroutine update_hecMESH_commtable( hecMESH, n_newnodes, n_neighbor_pe_new, neighbor_pe_list, &
    & import_index_new, import_item_new, export_index_new, export_item_new, node_ID_new )
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    integer(kind=kint) :: n_newnodes, n_neighbor_pe_new
    integer(kind=kint), allocatable :: neighbor_pe_list(:), import_index_new(:), import_item_new(:)
    integer(kind=kint), allocatable :: export_index_new(:), export_item_new(:), node_ID_new(:)
    integer(kind=kint) :: nitem

    if( nprocs == 1 ) return

    hecMESH%n_node = hecMESH%n_node + n_newnodes

    hecMESH%n_neighbor_pe = n_neighbor_pe_new

    if( associated(hecMESH%neighbor_pe)) deallocate(hecMESH%neighbor_pe)
    allocate(hecMESH%neighbor_pe(n_neighbor_pe_new))
    hecMESH%neighbor_pe(1:n_neighbor_pe_new) = neighbor_pe_list(1:n_neighbor_pe_new)

    if( associated(hecMESH%import_index)) deallocate(hecMESH%import_index)
    allocate(hecMESH%import_index(0:n_neighbor_pe_new))
    hecMESH%import_index(0:n_neighbor_pe_new) = import_index_new(0:n_neighbor_pe_new)

    nitem = import_index_new(n_neighbor_pe_new)
    if( associated(hecMESH%import_item)) deallocate(hecMESH%import_item)
    allocate(hecMESH%import_item(nitem))
    hecMESH%import_item(1:nitem) = import_item_new(1:nitem)

    if( associated(hecMESH%export_index)) deallocate(hecMESH%export_index)
    allocate(hecMESH%export_index(0:n_neighbor_pe_new))
    hecMESH%export_index(0:n_neighbor_pe_new) = export_index_new(0:n_neighbor_pe_new)

    nitem = export_index_new(n_neighbor_pe_new)
    if( associated(hecMESH%export_item)) deallocate(hecMESH%export_item)
    allocate(hecMESH%export_item(nitem))
    hecMESH%export_item(1:nitem) = export_item_new(1:nitem)

    if( associated(hecMESH%node_ID)) deallocate(hecMESH%node_ID)
    allocate(hecMESH%node_ID(size(node_ID_new)))
    hecMESH%node_ID(:) = node_ID_new(:)

    if( associated(hecMESH%shared_index)) deallocate(hecMESH%shared_index)
    allocate(hecMESH%shared_index(0:n_neighbor_pe_new))
    hecMESH%shared_index(0:n_neighbor_pe_new) = export_index_new(0:n_neighbor_pe_new) !set dummy to shared_index
    if( associated(hecMESH%shared_item)) deallocate(hecMESH%shared_item)
    allocate(hecMESH%shared_item(hecMESH%shared_index(n_neighbor_pe_new)))
    hecMESH%shared_item(1:hecMESH%shared_index(n_neighbor_pe_new)) = export_item_new(1:hecMESH%shared_index(n_neighbor_pe_new))

    call resize_int_array( hecMESH%n_node, hecMESH%global_node_ID )
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, hecMESH%node )

  end subroutine

  subroutine set_hecmwST_mpc( hecMESH, mpcs )
    type( hecmwST_local_mesh ), intent(inout) :: hecMESH     !< type mesh
    type(tMPCCond), allocatable :: mpcs(:)

    integer(kind=kint) :: i, j, idof, count, idx

    if( hecMESH%mpc%n_mpc > 0 ) stop "currently cannot use mpc and tied contact together"

    hecMESH%mpc%n_mpc = 3*size(mpcs)
    !if( associated(hecMESH%mpc%mpc_index) ) deallocate(hecMESH%mpc%mpc_index)
    allocate(hecMESH%mpc%mpc_index(0:hecMESH%mpc%n_mpc))
    hecMESH%mpc%mpc_index(0) = 0
    idx = 0
    count = 0
    do idof=1,3
      do i=1,size(mpcs)
        idx = idx + 1
        count = count + mpcs(i)%nitem
        hecMESH%mpc%mpc_index(idx) = count
      enddo
    enddo
    allocate(hecMESH%mpc%mpc_item(count))
    allocate(hecMESH%mpc%mpc_dof(count))
    allocate(hecMESH%mpc%mpc_val(count))
    allocate(hecMESH%mpc%mpc_const(count))
    !set value from mpcs
    idx = 0
    count = 0
    do idof=1,3
      do i=1,size(mpcs)
        do j=1,mpcs(i)%nitem
          idx = idx + 1
          hecMESH%mpc%mpc_item(idx) = mpcs(i)%pid(j)
          hecMESH%mpc%mpc_dof(idx) = idof
          hecMESH%mpc%mpc_val(idx) = mpcs(i)%coeff(j)
          hecMESH%mpc%mpc_const(idx) = 0.d0
        enddo
      enddo
    enddo
  end subroutine

  subroutine set_contact_structure( cstep, fstrSOLID, mpcs )
    integer(kind=kint), intent(in)         :: cstep      !< current step number
    type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
    type(tMPCCond), allocatable :: mpcs(:)

    integer(kind=kint) :: i, j, grpid, n_mpc, nitem, idof, dof(1:1000)

    do i = 1, size(fstrSOLID%contacts)
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
      if( fstrSOLID%contacts(i)%algtype /= CONTACTTIED ) cycle

      deallocate(fstrSOLID%contacts(i)%slave,fstrSOLID%contacts(i)%states)
      deallocate(fstrSOLID%contacts(i)%master)
      n_mpc = size(mpcs)
      allocate(fstrSOLID%contacts(i)%slave(n_mpc))
      allocate(fstrSOLID%contacts(i)%states(n_mpc))
      allocate(fstrSOLID%contacts(i)%master(n_mpc))

      do j=1,n_mpc
        nitem = mpcs(j)%nitem
        fstrSOLID%contacts(i)%slave(j) = mpcs(j)%pid(1)
        fstrSOLID%contacts(i)%states(j)%surface = j
        fstrSOLID%contacts(i)%master(j)%etype = fe_tri3n
        allocate(fstrSOLID%contacts(i)%master(j)%nodes(nitem))
        fstrSOLID%contacts(i)%master(j)%nodes(1:nitem) = mpcs(j)%pid(1:nitem)
        do idof=1,3
          fstrSOLID%contacts(i)%states(j)%state = CONTACTSTICK
          dof(1:nitem) = idof
          call init_mpc_cond(fstrSOLID%contacts(i)%states(j)%mpc_cond(idof), mpcs(j)%nitem)
          call set_mpc_cond(fstrSOLID%contacts(i)%states(j)%mpc_cond(idof), &
            &  nitem, mpcs(j)%pid(1:nitem), dof(1:nitem), mpcs(j)%coeff(1:nitem))
        enddo
      enddo
      exit

    enddo
  end subroutine


!========================================================================
! Get new mpcs coeff by singular value decomposition
!======================================================================== 
  subroutine get_newmpc_by_svd( mpcs_old, mpcs_new )
    type(tMPCCond), allocatable, intent(in)  :: mpcs_old(:)
    type(tMPCCond), allocatable, intent(out) :: mpcs_new(:)


    ! variables for work
    type(tMPCCond), allocatable :: mpcs_work(:)
    integer(kind=kint) :: i, nmpc_remian
    integer(kind=kint), parameter :: IMAX = 100

    allocate(mpcs_work(size(mpcs_old)))
    do i=1,size(mpcs_old)
      call copy_mpc_cond( mpcs_old(i), mpcs_work(i) )
    enddo

    do i=1,IMAX
      call get_newmpc_by_svd_main( mpcs_work, nmpc_remian )
      if( nmpc_remian == 0 ) exit
      if( i == IMAX ) then
        write(*,*) "cannot remove duplication for all mpcs, iter, nmpc_remain=",i,nmpc_remian
        stop
      endif
    enddo

    allocate(mpcs_new(size(mpcs_work)))
    do i=1,size(mpcs_work)
      !! cutoff small coeff
      call cutoff_small_coeff( mpcs_work(i), 1.d-4 )
      call copy_mpc_cond( mpcs_work(i), mpcs_new(i) )
    enddo

    do i=1,size(mpcs_work)
      call finalize_mpc_cond( mpcs_work(i) )
    enddo

  end subroutine

  !> get new mpcs coeff by singular value decomposition
  subroutine get_newmpc_by_svd_main( mpcs_work, nmpc_remian )
    type(tMPCCond), allocatable, intent(inout)  :: mpcs_work(:)
    integer(kind=kint), intent(out)  :: nmpc_remian

    ! variables for local to global index conversion
    integer(kind=kint) :: n_index, n_expanded
    integer(kind=kint), allocatable :: index(:)
    ! variables for grouping
    integer(kind=kint) :: n_mpc_group
    type(tMPCGroup), allocatable :: mpc_group(:) !< mpc group divieded by slave pids
    type(tMPCGroup), allocatable :: mpc_group_intermaster(:) !< mpc group between master nodes
    type(tMPCGroup) :: mpc_group_intermaster_all !< mpc group between master nodes

    integer(kind=kint) :: i, j, pid, idx
    integer(kind=kint) :: n_mpcs_new
    logical, parameter :: make_group = .true.

    ! compress pids
    call compress_pids( mpcs_work, n_index, index, .false. )

    ! expand obvious slave chain
    do i=1,10
      call expand_slave_chain( mpcs_work, n_index, n_expanded )
      if( n_expanded == 0 ) exit
    enddo

    ! create mpc groups
    if( make_group ) then !! divide mpc group by slave connectivity
      call grouping_mps( mpcs_work, n_index, n_mpc_group, mpc_group )
      allocate(mpc_group_intermaster(n_mpc_group))
    else !! do not divide mpcs
      n_mpc_group = 1
      allocate(mpc_group(n_mpc_group))
      mpc_group(1)%nitem = size(mpcs_work(:))
      mpc_group(1)%mpcs = mpcs_work(:)
      allocate(mpc_group_intermaster(1))
    endif

    ! singular value decomposition by group
    do i=1,n_mpc_group
      mpc_group_intermaster(i)%nitem = 0
      call regularize_mpc_bysvd( mpc_group(i), mpc_group_intermaster(i) )
    enddo
    call merge_mpc_group_intermaster( n_mpc_group, mpc_group_intermaster, mpc_group_intermaster_all )
    deallocate(mpc_group_intermaster)

    call choose_slave_intermaster_mpcs( n_index, n_mpc_group, mpc_group, mpc_group_intermaster_all )

    do i=1,n_mpc_group
      call expand_pids( mpc_group(i)%mpcs(1:mpc_group(i)%nitem), index )
    enddo
    call expand_pids( mpc_group_intermaster_all%mpcs(:), index )

    do i=1,size(mpcs_work)
      call finalize_mpc_cond( mpcs_work(i) )
    enddo
    deallocate(mpcs_work)

    idx = 0
    do i=1,n_mpc_group
      idx = idx + mpc_group(i)%nslaves
    enddo
    idx = idx + mpc_group_intermaster_all%nitem
    allocate(mpcs_work(idx))
    idx = 0
    do i=1,n_mpc_group
      do j=1,mpc_group(i)%nslaves
        idx = idx + 1
        mpcs_work(idx) = mpc_group(i)%mpcs(j)
      enddo
    enddo
    do i=1,mpc_group_intermaster_all%nitem
      idx = idx + 1
      mpcs_work(idx) = mpc_group_intermaster_all%mpcs(i)
    enddo
    nmpc_remian = mpc_group_intermaster_all%nitem

    call remove_tinycoeff( mpcs_work, nmpc_remian )
  end subroutine

  subroutine remove_tinycoeff(mpcs_work, nmpc_remian)
    implicit none
    type(tMPCCond), allocatable, intent(inout) :: mpcs_work(:)
    integer(kind=kint), intent(inout)          :: nmpc_remian

    integer(kind=kint) :: i, j, k, n
    real(kind=kreal)   :: max_val, factor, tmp_real
    integer(kind=kint) :: tmp_int

    do i = 1, size(mpcs_work)
      n = mpcs_work(i)%nitem
      if (n <= 1) cycle

      max_val = 0.d0
      do j = 2, n
        if ( dabs(mpcs_work(i)%coeff(j)) > max_val ) then
          max_val = dabs(mpcs_work(i)%coeff(j))
        end if
      end do

      if ( dabs(mpcs_work(i)%coeff(1)) > 1.d-6 * max_val ) cycle

      do j = 1, n - 1
        mpcs_work(i)%coeff(j) = mpcs_work(i)%coeff(j+1)
        mpcs_work(i)%pid(j)   = mpcs_work(i)%pid(j+1)
        mpcs_work(i)%dof(j)   = mpcs_work(i)%dof(j+1)
      end do
      mpcs_work(i)%nitem = n - 1
      n = mpcs_work(i)%nitem

      k = 1
      max_val = abs(mpcs_work(i)%coeff(1))
      do j = 2, n
        if ( abs(mpcs_work(i)%coeff(j)) > max_val ) then
          max_val = abs(mpcs_work(i)%coeff(j))
          k = j
        end if
      end do

      if (k /= 1) then
        tmp_real = mpcs_work(i)%coeff(1)
        mpcs_work(i)%coeff(1) = mpcs_work(i)%coeff(k)
        mpcs_work(i)%coeff(k) = tmp_real
        tmp_int = mpcs_work(i)%pid(1)
        mpcs_work(i)%pid(1) = mpcs_work(i)%pid(k)
        mpcs_work(i)%pid(k) = tmp_int
        tmp_int = mpcs_work(i)%dof(1)
        mpcs_work(i)%dof(1) = mpcs_work(i)%dof(k)
        mpcs_work(i)%dof(k) = tmp_int
      end if

      factor = mpcs_work(i)%coeff(1)
      if (factor /= 0.d0) then
        do j = 1, n
          mpcs_work(i)%coeff(j) = mpcs_work(i)%coeff(j) / factor
        end do
      end if

      nmpc_remian = nmpc_remian + 1
    end do
  end subroutine remove_tinycoeff

  subroutine compress_pids( mpcs, n_index, index, slave_first, n_slave )
    type(tMPCCond), allocatable, intent(inout)     :: mpcs(:)  !< mpc condition
    integer(kind=kint), intent(inout)              :: n_index  !< number of pids
    integer(kind=kint), allocatable, intent(inout) :: index(:) !< pid convert table
    logical, intent(in)                            :: slave_first
    integer(kind=kint), intent(inout), optional    :: n_slave  !< number of pids

    integer(kind=kint) :: i, j, pid, idx
    integer(kind=kint) :: n_mpcs, nitem
    integer(kind=kint), allocatable :: slave_first_ids(:), slave_first_ids_inv(:)
    integer(kind=kint), allocatable :: iwork(:)

    n_mpcs = size(mpcs)

    ! create index table
    n_index = 0
    do i=1,n_mpcs
      nitem = mpcs(i)%nitem
      n_index = n_index + nitem
    enddo
    allocate(index(n_index))
    n_index = 0
    do i=1,n_mpcs
      nitem = mpcs(i)%nitem
      index(n_index+1:n_index+nitem) = mpcs(i)%pid(1:nitem)
      n_index = n_index + nitem
    enddo
    !! sort index
    call qsort_int_array(index, 1, n_index)
    call uniq_int_array(index, n_index, n_index)

    ! compress pids
    do i=1,n_mpcs
      do j=1,mpcs(i)%nitem
        pid = mpcs(i)%pid(j)
        call bsearch_int_array(index, 1, n_index, pid, idx)
        mpcs(i)%pid(j) = idx
      enddo
    enddo

    ! slave_first
    if( slave_first ) then
      allocate(slave_first_ids(n_index),iwork(n_index))
      iwork(:) = 0
      do i=1,n_mpcs
        iwork(mpcs(i)%pid(1)) = 1
      enddo
      idx = 0
      !get slave
      do i=1,n_index
        if( iwork(i) == 0 ) cycle
        idx = idx + 1
        slave_first_ids(idx) = i 
      enddo
      if( present(n_slave) ) n_slave = idx
      !get master
      do i=1,n_index
        if( iwork(i) == 1 ) cycle
        idx = idx + 1
        slave_first_ids(idx) = i 
      enddo
      ! crate local pid to original pid
      allocate(slave_first_ids_inv(n_index))
      do i=1,n_index
        idx = slave_first_ids(i)
        slave_first_ids_inv(idx) = i 
      enddo
      ! update pids by new pid
      do i=1,n_mpcs
        do j=1,mpcs(i)%nitem
          pid = mpcs(i)%pid(j)
          mpcs(i)%pid(j) = slave_first_ids_inv(pid)
        enddo
      enddo
      iwork(1:n_index) = index(1:n_index)
      do i=1,n_index
        index(i) = iwork(slave_first_ids(i))
      enddo
    endif
  end subroutine

  subroutine expand_pids( mpcs, index )
    type(tMPCCond), intent(inout)  :: mpcs(:)  !< mpc condition
    integer(kind=kint), allocatable, intent(in) :: index(:) !< pid convert table

    integer(kind=kint) :: i, j, n_mpcs

    n_mpcs = size(mpcs)
    ! expand pids
    do i=1,n_mpcs
      do j=1,mpcs(i)%nitem
        mpcs(i)%pid(j) = index(mpcs(i)%pid(j))
      enddo
    enddo
  end subroutine

  subroutine cutoff_small_coeff( mpc, threshold )
    type(tMPCCond), intent(inout) :: mpc
    real(kind=kreal), intent(in)  :: threshold

    integer(kind=kint) :: i, nitem
    real(kind=kreal)   :: coeff_sum

    ! suppose coeff(1) == 1. if not, do nothing
    if( dabs(mpc%coeff(1)-1.d0) > 1.d-6 ) return

    nitem = 1
    do i=2,mpc%nitem
      if( dabs(mpc%coeff(i)) < threshold ) cycle
      nitem = nitem + 1
      mpc%pid(nitem) = mpc%pid(i) 
      mpc%dof(nitem) = mpc%dof(i) 
      mpc%coeff(nitem) = mpc%coeff(i) 
    enddo
    mpc%nitem = nitem

  end subroutine

  subroutine expand_slave_chain( mpcs, n_index, n_expanded )
    type(tMPCCond), allocatable, intent(inout)     :: mpcs(:)  !< mpc condition
    integer(kind=kint), intent(in)                 :: n_index  !< number of pids
    integer(kind=kint), intent(out)                :: n_expanded

    logical, allocatable :: is_slave(:)
    integer(kind=kint), allocatable :: expanding_mpcs(:)
    logical :: expanding_flag
    integer(kind=kint) :: i, j, slave, pid, n_expsl, expsl(15)
    integer(kind=kint) :: n_mpcs


    integer(kind=kint) :: iter, jter, mpcid, mpcid_tmp
    integer(kind=kint), allocatable :: group_id(:), group_id_mpc(:)
    type(hecmwST_varray_int), allocatable :: table_pid2mpc(:), group_mpcs(:)
    type(hecmwST_varray_int) :: slave_ids
    integer(kind=kint) :: pid_finished
    logical :: target_found
    logical, parameter :: DEBUG=.false.

    n_mpcs = size(mpcs)

    !! create slave pid flag
    allocate(is_slave(n_index))
    is_slave = .false.
    do i=1,n_mpcs
      is_slave(mpcs(i)%pid(1)) = .true.
    enddo

    !! create expand mpclist
    allocate(expanding_mpcs(n_index))
    expanding_mpcs(:) = -1
    do i=1,n_mpcs
      !! if no master pid are not slave for another mpc, use this mpc to expand
      expanding_flag = .true.
      do j=2,mpcs(i)%nitem !check masters
        if( is_slave(mpcs(i)%pid(j)) ) then ! if slave chain found, skip expanding
          expanding_flag = .false.
          exit
        endif
      enddo
      if( expanding_flag ) then
        slave = mpcs(i)%pid(1)
        if( expanding_mpcs(slave) < 0 ) then
          expanding_mpcs(slave) = i
        else !if expanding_mpcs is already found
          mpcid = expanding_mpcs(slave)
          if( mpcs(i)%nitem < mpcs(mpcid)%nitem ) then !if shorter mpc found, use it to expand
            expanding_mpcs(slave) = i
          endif
        endif
      endif
    enddo

    !! expand mpcs
    n_expanded = 0
    do i=1,n_mpcs
      !!! get expand slave lists
      n_expsl = 0
      do j=2,mpcs(i)%nitem !check masters
        pid = mpcs(i)%pid(j)
        if( expanding_mpcs(pid) > 0 ) then
          n_expsl = n_expsl + 1
          if( n_expsl > size(expsl) ) stop ""
          expsl(n_expsl) = expanding_mpcs(pid)
        endif
      enddo
      do j=1,n_expsl
        call expand_mpc_cond( mpcs(i), mpcs(expsl(j)) )
      enddo
      if( n_expsl > 0 ) n_expanded = n_expanded + 1
    enddo

  end subroutine

  subroutine grouping_mps( mpcs, n_index, n_mpc_group, mpc_group )
    type(tMPCCond), allocatable, intent(in)     :: mpcs(:)  !< mpc condition
    integer(kind=kint), intent(in)              :: n_index  !< number of pids
    integer(kind=kint), intent(inout)           :: n_mpc_group
    type(tMPCGroup), allocatable, intent(inout) :: mpc_group(:) !< mpc group divieded by slave pids

    integer(kind=kint) :: i, j, iter, jter, pid, mpcid, mpcid_tmp
    integer(kind=kint) :: n_mpcs
    integer(kind=kint), allocatable :: group_id(:), group_id_mpc(:)
    type(hecmwST_varray_int), allocatable :: table_pid2mpc(:), group_mpcs(:)
    type(hecmwST_varray_int) :: slave_ids
    integer(kind=kint) :: pid_finished
    logical :: target_found
    logical, parameter :: DEBUG=.false.

    n_mpcs = size(mpcs)

    !! initialize group_id: set 0 to slave pid and 1 to master pid
    allocate(group_id(n_index))
    group_id = -1
    do i=1,n_mpcs
      group_id(mpcs(i)%pid(1)) = 0
    enddo

    ! create group id list
    !! create slave pid to mpc cond list
    call HECMW_varray_int_initialize_all(table_pid2mpc,n_index,4)
    do i=1,n_mpcs
      do j=1,mpcs(i)%nitem
        pid = mpcs(i)%pid(j)
        if( group_id(pid) /= 0 ) cycle
        call HECMW_varray_int_add_if_not_exits(table_pid2mpc(pid),i)
      enddo
    enddo

    !! get connectivity
    n_mpc_group = 0
    pid_finished = 0  ! cache starting mpc id for speed
    call HECMW_varray_int_initialize_all(group_mpcs,n_mpcs,10)
    do iter=1,n_mpcs
      !!! search ungrouped mpc
      !!! scan pid to initialize mpc group list
      target_found = .false.
      do pid=pid_finished+1,n_index
        if( group_id(pid) /= 0 ) cycle
        n_mpc_group = n_mpc_group + 1  ! increment group id
        !! add mpc ids to group_mpcs from table_pid2mpc using slave pid as key  
        do i=1,HECMW_varray_int_get_nitem(table_pid2mpc(pid))
          mpcid = HECMW_varray_int_get_item(table_pid2mpc(pid),i)
          call HECMW_varray_int_add( group_mpcs(n_mpc_group), mpcid )
        enddo
        pid_finished = pid
        target_found = .true.
        exit
      enddo
      !!! exit if ungrouped mpc is not found
      if( .not. target_found ) exit

      !! depth-first search by slave-slave connectivity
      do jter=1,n_mpcs
        !!! exit if all mpcs are scanned
        if( jter > HECMW_varray_int_get_nitem(group_mpcs(n_mpc_group)) ) exit
        !!! pop mpc
        mpcid = HECMW_varray_int_get_item(group_mpcs(n_mpc_group),jter)
        do i=1,mpcs(mpcid)%nitem !scan pid
          pid = mpcs(mpcid)%pid(i)
          if( group_id(pid) < 0 ) then
            cycle !!! do nothing to master
          else if( group_id(pid) == 0 ) then
            !!! if slave pid is found, add mpcs which share the slave pid
            do j=1,HECMW_varray_int_get_nitem(table_pid2mpc(pid))
              mpcid_tmp = HECMW_varray_int_get_item(table_pid2mpc(pid),j)
              call HECMW_varray_int_add_if_not_exits( group_mpcs(n_mpc_group), mpcid_tmp )
            enddo
            !!! set current group id to list
            group_id(pid) = n_mpc_group
          else if( group_id(pid) == n_mpc_group ) then
            cycle  !!! do nothing if current group_id is re-founded
          else
            write(*,*) "error in grouping_mps: group_id mismatch found"
            write(*,*) "group_id(pid),n_mpc_group",group_id(pid),n_mpc_group
            stop
          endif
        enddo
      enddo

      !! check group_mpcs(n_mpc_group) if all slave found
      if( DEBUG ) then
        write(*,*) "## check group mpc: n_mpc_group=",n_mpc_group
        do i=1,HECMW_varray_int_get_nitem(group_mpcs(n_mpc_group))
          mpcid = HECMW_varray_int_get_item(group_mpcs(n_mpc_group),i)
          write(*,*) "mpcid",mpcid
          do j=1,mpcs(mpcid)%nitem
            pid = mpcs(mpcid)%pid(j)
            if( group_id(pid) < 0 ) then
              cycle !!! do nothing to master
            else if( group_id(pid) == 0 ) then
              write(*,*) "error in grouping_mps: ungrouped slave found",i,pid
            else if( group_id(pid) == n_mpc_group ) then
              cycle  !!! do nothing if current group_id is re-founded
            else
              write(*,*) "error in grouping_mps: group_id mismatch found"
              write(*,*) "group_id(pid),n_mpc_group",group_id(pid),n_mpc_group
            endif
          enddo
        enddo
      endif
    enddo

    !! check grouping
    if( .true. ) then
      do pid=1,n_index
        if( group_id(pid) == 0 ) write(*,*) "founded ungrouped pid:",pid
      enddo
      allocate(group_id_mpc(n_mpcs))
      group_id_mpc(:) = 0
      do i=1,n_mpc_group
        do j=1,HECMW_varray_int_get_nitem(group_mpcs(i))
          mpcid = HECMW_varray_int_get_item(group_mpcs(i),j)
          if( group_id_mpc(mpcid) /= 0 ) &
            & write(*,*) "founded double-grouped mpc,groupid:",mpcid,i,group_id_mpc(mpcid)
            group_id_mpc(mpcid) = i
        enddo
      enddo
      do mpcid=1,n_mpcs
        if( group_id_mpc(mpcid) == 0 ) write(*,*) "founded ungrouped mpc:",mpcid
      enddo
    endif

    allocate(mpc_group(n_mpc_group))
    call HECMW_varray_int_initialize(slave_ids,10)
    do i=1,n_mpc_group
      mpc_group(i)%nitem = HECMW_varray_int_get_nitem(group_mpcs(i))
      call HECMW_varray_int_clear(slave_ids)
      allocate(mpc_group(i)%mpcs(mpc_group(i)%nitem))
      do j=1,mpc_group(i)%nitem
        mpcid = HECMW_varray_int_get_item(group_mpcs(i),j)
        mpc_group(i)%mpcs(j) = mpcs(mpcid)
        call HECMW_varray_int_add_if_not_exits(slave_ids,mpcs(mpcid)%pid(1))
      enddo
      mpc_group(i)%nslaves = HECMW_varray_int_get_nitem(slave_ids)
    enddo
  end subroutine

  subroutine solve_dof_with_row_pivoting(Amat,nrow,mcol,index,nrow_to_be_solved)
    real(kind=kreal), allocatable, intent(inout)  :: Amat(:,:)
    integer(kind=kint), intent(in)                :: nrow
    integer(kind=kint), intent(in)                :: mcol
    integer(kind=kint), allocatable, intent(inout) :: index(:)
    integer(kind=kint), intent(in)                :: nrow_to_be_solved

    integer(kind=kint) :: i, j, k, i_maxpiv, j_maxpiv, tmpidx
    real(kind=kreal) :: pivot, maxpivot, coeff, diff
    real(kind=kreal) :: tmprow(mcol), tmpcol(nrow) ! for swap

    ! forward elimination
    do i=1,nrow_to_be_solved
      ! row pivoting
      maxpivot = 0.d0
      do k=i,nrow
        pivot = Amat(k,i)
        if( dabs(pivot) < 1.d-10 ) cycle
        if( dabs(pivot) - maxpivot > 1.d-3*maxpivot ) then
          maxpivot = dabs(pivot)
          i_maxpiv = k
        endif
      enddo

      if( i_maxpiv == 0 .and. j_maxpiv == 0 ) continue

      !swap row i and row i_maxpiv
      if( i_maxpiv > 0 .and. i /= i_maxpiv ) then
        tmprow(i:mcol) = Amat(i,i:mcol)
        Amat(i,i:mcol) = Amat(i_maxpiv,i:mcol)
        Amat(i_maxpiv,i:mcol) = tmprow(i:mcol)
      endif

      if( dabs( Amat(i,i) ) < 1.d-10 ) cycle

      ! gauss elim
      coeff = 1.d0/Amat(i,i)
      !! normalize first line
      do j=i,mcol
        Amat(i,j) = coeff*Amat(i,j)
      enddo
      !! delete second line and after
      do j=i+1,nrow
        pivot = Amat(j,i)
        if( dabs(pivot) < 1.d-10 ) cycle
        do k=i,mcol
          Amat(j,k) = Amat(j,k)-pivot*Amat(i,k)
        enddo

        !! normalize line
        pivot = Amat(j,i)
        do k=i+1,mcol
          if( dabs(Amat(j,k)) - dabs(pivot) > 1.d-4 ) pivot = Amat(j,k)
        enddo
        pivot = 1.d0/pivot
        do k=i,mcol
          Amat(j,k) = pivot*Amat(j,k)
        enddo
      enddo
    enddo

    ! backward substitution
    do i=nrow_to_be_solved,2,-1
      do j=i-1,1,-1
        pivot = Amat(j,i)
        do k=j,mcol
          Amat(j,k) = Amat(j,k)-pivot*Amat(i,k)
        enddo
      enddo
    enddo
  end subroutine

  subroutine create_Amat_from_mpcs( mpcs, Amat, n_mpcs, n_index )
    type(tMPCCond), allocatable, intent(inout) :: mpcs(:)
    real(kind=kreal), allocatable, intent(inout) :: Amat(:,:)
    integer(kind=kint), intent(inout) :: n_mpcs
    integer(kind=kint), intent(inout) :: n_index

    integer(kind=kint) :: i, j, pid

    allocate(Amat(n_mpcs,n_index))
    Amat(:,:) = 0.d0

    do i=1,n_mpcs
      do j=1,mpcs(i)%nitem
        pid = mpcs(i)%pid(j)
        Amat(i,pid) = mpcs(i)%coeff(j)
      enddo
    enddo
  end subroutine

  subroutine set_mpc_from_row_Amat( mpc, Arow, n_index )
    type(tMPCCond), intent(inout) :: mpc
    real(kind=kreal), intent(in)   :: Arow(:)
    integer(kind=kint), intent(in)              :: n_index

    integer(kind=kint) :: i, j, nitem
    real(kind=kreal), parameter :: threshold = 1.d-4
    integer(kind=kint) :: pid(n_index),dof(n_index)
    real(kind=kreal)   :: coeff(n_index)

    nitem = 0
    do i=1,n_index
      if( dabs(Arow(i)) < threshold ) cycle
      nitem = nitem + 1
    enddo
    call init_mpc_cond(mpc,nitem)

    dof(1:nitem) = 1

    nitem = 0
    do i=1,n_index
      if( dabs(Arow(i)) < threshold ) cycle
      nitem = nitem + 1
      pid(nitem) = i
      coeff(nitem) = Arow(i)
    enddo
    call set_mpc_cond(mpc, nitem, pid(1:nitem), dof(1:nitem), coeff(1:nitem))
  end subroutine

  subroutine regularize_mpc_bysvd( mpc_group, mpc_group_intermaster )
    type(tMPCGroup), intent(inout) :: mpc_group !< mpc group divieded by slave pids
    type(tMPCGroup), intent(inout) :: mpc_group_intermaster !< mpc group between master nodes

    integer(kind=kint) :: i, j, n_mpcs, pid
    ! variables for local to global index conversion
    integer(kind=kint) :: n_index_local
    integer(kind=kint), allocatable :: index_local(:)
    ! variables for svd
    real(kind=kreal), allocatable :: Amat(:,:), Amat_ori(:,:), sval(:), Vs_inv(:,:)
    integer(kind=kint) :: n_slave, n_active_mpcs

    n_mpcs = mpc_group%nitem
    if( n_mpcs == 1 ) return !do nothing for single mpc

    call compress_pids( mpc_group%mpcs, n_index_local, index_local, .true., n_slave )

    !if( n_mpcs > 1 ) then !debug
    !  write(109,*) "mpc_compressed"
    !  call print_mpc_conditions(109, mpc_group%mpcs)
    !endif

    call create_Amat_from_mpcs( mpc_group%mpcs, Amat, n_mpcs, n_index_local )
    ! save original
    allocate(Amat_ori(n_mpcs,n_index_local))
    Amat_ori(:,:) = Amat(:,:)

    if( DEBUG > 0 .and. n_mpcs > 1 ) then !debug
      write(109,*) "Amat_ori,n_mpcs,n_index_local,n_slave",n_mpcs,n_index_local,n_slave
      write(109,*) index_local(1:n_index_local)
      call print_densemat(Amat_ori(1:n_mpcs,1:n_index_local),n_mpcs,n_index_local,109)
    endif

    ! singular value decomposition
    allocate(sval(n_mpcs))
    call svd_coeffmat(Amat,sval,n_mpcs,n_index_local)

    if( DEBUG > 0 .and. n_mpcs > 1 ) then !debug
      write(109,*) "Amat_svd,n_mpcs,n_index_local,n_slave",n_mpcs,n_index_local,n_slave
      call print_densemat(Amat(1:n_mpcs,1:n_index_local),n_mpcs,n_index_local,109)
    endif

    !! count active slave
    n_active_mpcs = n_mpcs
    do i=1,n_mpcs
      if( sval(i) > 0.01d0 ) cycle
      n_active_mpcs = i-1
      exit
    enddo
    n_slave = min(n_slave,n_active_mpcs)
    if( n_active_mpcs == n_mpcs ) Amat(:,:) = Amat_ori(:,:) ! do not use
    call solve_dof_with_row_pivoting(Amat,n_active_mpcs,n_index_local,index_local,n_slave)

    if( DEBUG > 0 .and. n_mpcs > 1 ) then !debug
      write(109,*) "Amat,n_active_mpcs,n_index_local,n_slave",n_active_mpcs,n_index_local,n_slave
      write(109,*) index_local(1:n_index_local)
      write(109,*) "sval",sval
      call print_densemat(Amat(1:n_active_mpcs,1:n_index_local),n_active_mpcs,n_index_local,109)
    endif

    ! set new coeff to mpc_group
    mpc_group%nitem = n_slave
    do i=1,n_slave
      call finalize_mpc_cond( mpc_group%mpcs(i) )
      call set_mpc_from_row_Amat( mpc_group%mpcs(i), Amat(i,:), n_index_local )
    enddo
    call expand_pids( mpc_group%mpcs(1:n_slave), index_local )

    if( DEBUG > 0 .and. n_mpcs > 1 ) then !debug
      write(109,*) "mpc_compressed"
      call print_mpc_conditions(109, mpc_group%mpcs(1:n_slave))
    endif

    mpc_group_intermaster%nitem = n_active_mpcs - n_slave
    if( mpc_group_intermaster%nitem > 0 ) then
      allocate(mpc_group_intermaster%mpcs(mpc_group_intermaster%nitem))
      do i=n_slave+1,n_active_mpcs
        call set_mpc_from_row_Amat( mpc_group_intermaster%mpcs(i-n_slave), Amat(i,:), n_index_local )
      enddo
      call expand_pids( mpc_group_intermaster%mpcs(1:mpc_group_intermaster%nitem), index_local )

      if( DEBUG > 0 .and. n_mpcs > 1 ) then !debug
        write(109,*) "mpc_compressed_intermaster"
        call print_mpc_conditions(109, mpc_group_intermaster%mpcs(1:mpc_group_intermaster%nitem))
      endif
  
    endif

  end subroutine

  subroutine merge_mpc_group_intermaster( n_mpc_group, mpc_group_intermaster, mpc_group_intermaster_all )
    integer(kind=kint), intent(in) :: n_mpc_group  !< number of pids
    type(tMPCGroup), intent(inout) :: mpc_group_intermaster(:)  !< mpc group between master nodes
    type(tMPCGroup), intent(inout) :: mpc_group_intermaster_all !< mpc group between master nodes(merged)

    integer(kind=kint) :: i, j, nitem_all

    nitem_all = 0
    do i=1,n_mpc_group
      if( mpc_group_intermaster(i)%nitem == 0 ) cycle
      nitem_all = nitem_all + mpc_group_intermaster(i)%nitem
    enddo

    mpc_group_intermaster_all%nitem = nitem_all
    mpc_group_intermaster_all%nslaves = nitem_all
    allocate(mpc_group_intermaster_all%mpcs(nitem_all))
    nitem_all = 0
    do i=1,n_mpc_group
      if( mpc_group_intermaster(i)%nitem == 0 ) cycle
      do j=1,mpc_group_intermaster(i)%nitem
        nitem_all = nitem_all + 1
        call copy_mpc_cond( mpc_group_intermaster(i)%mpcs(j), mpc_group_intermaster_all%mpcs(nitem_all) )
        call finalize_mpc_cond( mpc_group_intermaster(i)%mpcs(j) )
      enddo
    enddo
  end subroutine

  subroutine choose_slave_intermaster_mpcs( n_index, n_mpc_group, mpc_group, mpc_group_intermaster_all )
    integer(kind=kint), intent(in) :: n_index  !< number of pids
    integer(kind=kint), intent(in) :: n_mpc_group  !< number of mpc_group
    type(tMPCGroup), intent(inout) :: mpc_group(:)  !< mpc group
    type(tMPCGroup), intent(inout) :: mpc_group_intermaster_all !< mpc group between master nodes(merged)

    integer(kind=kint) :: i, j, nitem, pid, minidx, tmppid
    real(kind=kreal) :: tmpcoeff, tmpcoeff0

    do i=1,mpc_group_intermaster_all%nitem
      nitem = mpc_group_intermaster_all%mpcs(i)%nitem
      ! choose maxcoeff
      minidx = 1
      tmpcoeff = dabs(mpc_group_intermaster_all%mpcs(i)%coeff(1))
      do j=2,nitem
        tmpcoeff0 = dabs(mpc_group_intermaster_all%mpcs(i)%coeff(j))
        if( tmpcoeff0 > tmpcoeff ) then
          tmpcoeff = tmpcoeff0
          minidx = j
        endif
      enddo
      do j=1,nitem
        tmpcoeff0 = dabs(mpc_group_intermaster_all%mpcs(i)%coeff(j))
        if( dabs(tmpcoeff0-tmpcoeff) < 1.d-4*tmpcoeff ) then
          minidx = j
          exit
        endif
      enddo

      !swap j and 1
      tmppid = mpc_group_intermaster_all%mpcs(i)%pid(minidx)
      tmpcoeff = mpc_group_intermaster_all%mpcs(i)%coeff(minidx)
      mpc_group_intermaster_all%mpcs(i)%pid(minidx) = mpc_group_intermaster_all%mpcs(i)%pid(1)
      mpc_group_intermaster_all%mpcs(i)%coeff(minidx) = mpc_group_intermaster_all%mpcs(i)%coeff(1)
      mpc_group_intermaster_all%mpcs(i)%pid(1) = tmppid
      mpc_group_intermaster_all%mpcs(i)%coeff(1) = tmpcoeff
      do j=1,nitem
        mpc_group_intermaster_all%mpcs(i)%coeff(j) = mpc_group_intermaster_all%mpcs(i)%coeff(j)/tmpcoeff
      enddo
    enddo

  end subroutine

  subroutine svd_coeffmat(Amat,sval,m,n)
    real(kind=kreal), intent(inout) :: Amat(m,n)
    real(kind=kreal), intent(inout) :: sval(m)
    integer(kind=kint), intent(in)  :: m !< nrow
    integer(kind=kint), intent(in)  :: n !< ncol

    integer, parameter :: nb = 64
    integer :: i, info, lda, ldu, ldvt, lwork, msize
    real(kind=kreal), allocatable :: b(:), s(:), u(:,:), vt(:,:), work(:)
    real(kind=kreal) :: dummy(1,1)

    lda = m
    ldu = m
    ldvt = n
    msize = min(m,n)
    allocate (s(msize), vt(ldvt,n), u(ldu,m), b(msize))

    lwork = -1
    call dgesvd('A', 'A', m, n, Amat, lda, s, u, ldu, vt, ldvt, dummy, lwork, info)

    lwork = max(m+4*n+nb*(m+n), int(dummy(1,1)))
    allocate (work(lwork))

    call dgesvd('A', 'A', m, n, Amat, lda, s, u, ldu, vt, ldvt, work, lwork, info)

    sval(:) = 0.d0
    sval(1:msize) = s(1:msize)

    if (info/=0) Then
      Write (*,*) 'Failure in DGESVD. INFO =', info
    endif

    Amat(:,:) = 0.d0
    Amat(1:msize,1:n) = vt(1:msize,1:n)
  end subroutine

  subroutine lufact_coeffmat_trans(Amat,m,n)
    real(kind=kreal), intent(inout) :: Amat(m,n)
    integer(kind=kint), intent(in)  :: m !< nrow
    integer(kind=kint), intent(in)  :: n !< ncol

    integer(kind=kint) :: i, j, ifail, info, lda, mt, nt
    
    integer(kind=kint), allocatable :: ipiv(:), ipiv2(:)
    real(kind=kreal) :: AmatT(n,m)

    AmatT = transpose(Amat)

    mt = n
    nt = m

    !debug
    call print_densemat(AmatT,mt,nt,201)

    lda = mt
    Allocate (ipiv(nt),ipiv2(mt))
    call dgetrf(mt, nt, AmatT, lda, ipiv, info)

    if (info/=0) Then
      Write (*,*) 'Failure in DGETRF. INFO =', info
    endif

    !get permutation from ipiv
    do i=1,mt
      ipiv2(i) = i
    enddo
    do i=1,nt
      j = ipiv(i)
      !swap ipiv2(i) and ipiv2(j)
      lda = ipiv2(j)
      ipiv2(j) = ipiv2(i)
      ipiv2(i) = lda
    enddo

    !clear upper A
    do j=2,m
      do i=1,j-1
        AmatT(i,j) = 0.d0
      enddo
      AmatT(j,j) = 1.d0
    enddo

    !reorder
    do i=1,n
      lda = ipiv2(i)
      do j=1,m
        Amat(j,lda) = AmatT(i,j)
      enddo
    enddo

    return

    !solve upper A
    do i=m-1,1
      do j=1,i
        Amat(i,j) = 0.d0
      enddo
    enddo

    write(*,*) "m,n",m,n
    stop
  end subroutine

  subroutine create_mpcs_by_coeffmat(Amat,n_mpcs_new,n_index,index,mpcs_new)
    real(kind=kreal), intent(inout) :: Amat(:,:)
    integer(kind=kint), intent(in)  :: n_mpcs_new
    integer(kind=kint), intent(in)  :: n_index !< nrow
    integer(kind=kint), intent(in)  :: index(:)
    type(tMPCCond), allocatable, intent(inout) :: mpcs_new(:)

    integer(kind=kint) :: i, j, idx
    integer(kind=kint) :: pid(n_index), dof(n_index)
    real(kind=kreal)   :: coeff(n_index)

    allocate(mpcs_new(n_mpcs_new))
    dof(:) = 1

    do i=1,n_mpcs_new
      idx = 0
      do j=1,n_index
        if( dabs(Amat(i,j)) < 1.d-4 ) cycle
        idx = idx + 1
        pid(idx) = index(j)
        coeff(idx) = Amat(i,j)
      enddo

      call init_mpc_cond(mpcs_new(i),idx)
      call set_mpc_cond(mpcs_new(i), idx, pid(1:idx), dof(1:idx), coeff(1:idx))
    enddo
  end subroutine

  !> for debug
  subroutine print_mpc_use( mpcs )
    type(tMPCCond), allocatable, intent(in) :: mpcs(:)

    integer(kind=kint) :: i, j, maxi, c1
    integer(kind=kint), allocatable :: counts(:)

    maxi = 0
    do i=1,size(mpcs)
      do j=1,mpcs(i)%nitem
        if( mpcs(i)%pid(j) > maxi ) maxi = mpcs(i)%pid(j)
      enddo
    enddo

    allocate(counts(maxi))
    counts(:) = 0
    do i=1,size(mpcs)
      do j=1,mpcs(i)%nitem
        counts(mpcs(i)%pid(j)) = counts(mpcs(i)%pid(j)) + 1
      enddo
    enddo

    write(*,*) "size(mpcs),",size(mpcs)
    c1 = 0
    do i=1,maxi
      if( counts(i) == 1 ) c1 = c1 + 1
    enddo
    write(*,*) "counts 1,",c1
  end

  !> for debug
  subroutine print_index(idx,n,fnum)
    integer(kind=kint), intent(in)  :: idx(n)
    integer(kind=kint), intent(in)  :: n
    integer(kind=kint), intent(in)  :: fnum

    integer(kind=kint) :: i
    character(len=128) :: fmt

    write(fmt,'(A,I0,A)') '(',n,'I7)'
    write(fnum,trim(fmt)) idx(1:n)
  end subroutine

  !> for debug
  subroutine print_densemat(A,m,n,fnum)
    real(kind=kreal), intent(inout) :: A(m,n)
    integer(kind=kint), intent(in)  :: m !< nrow
    integer(kind=kint), intent(in)  :: n !< ncol
    integer(kind=kint), intent(in)  :: fnum

    integer(kind=kint) :: i
    character(len=128) :: fmt

    write(fmt,'(A,I0,A)') '(',n,'f7.3)'
    do i=1,m
      write(fnum,trim(fmt)) A(i,1:n)
    enddo
  end subroutine

  subroutine debug_svd(m,n,u,s,vt)
    integer(kind=kint), intent(in)  :: m !< nrow
    integer(kind=kint), intent(in)  :: n !< ncol
    real(kind=kreal), intent(inout) :: u(m,m)
    real(kind=kreal), intent(inout) :: s(n)
    real(kind=kreal), intent(inout) :: vt(n,n)

    integer(kind=kint) :: i,j
    integer(kind=kint), parameter :: fnum = 121
    real(kind=kreal) :: amat(m,n),tmp(m,n)

    write(fnum,*) ""
    write(fnum,*) "umat"
    call print_densemat(u,m,m,fnum)    

    write(fnum,*) ""
    write(fnum,*) "s"
    write(fnum,'(f12.6)') s(1:n)

    write(fnum,*) ""
    write(fnum,*) "vtmat"
    call print_densemat(vt,n,n,fnum)

    write(123,*) ""
    write(123,*) "u*s*vt"
    tmp(1:m,1:n) = vt(1:m,1:n)
    do i=1,m
      tmp(i,1:n) = s(i)*vt(i,1:n)
    enddo
    amat(1:m,1:n) = matmul(u(1:m,1:m),tmp(1:m,1:n))
    call print_densemat(amat,m,n,123)

  end subroutine

!========================================================================
! With the update of hecMESH, we will expand the data storage structure.
!======================================================================== 
  subroutine resize_structures( hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrEIG )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid), intent(inout)      :: fstrSOLID   !< type fstr_solid
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_eigen)                     :: fstrEIG

    ! resize fstrSOLID
    call resize_structures_static( hecMESH, hecMAT, fstrSOLID )
    call resize_real_array( 3*hecMESH%n_node, fstrEIG%mass )
    call resize_real_array2( hecMESH%n_node*hecMESH%n_dof, size(fstrDYNAMIC%DISP(1,:)), fstrDYNAMIC%DISP )
    call resize_real_array2( hecMESH%n_node*hecMESH%n_dof, size(fstrDYNAMIC%VEL(1,:)), fstrDYNAMIC%VEL )
    call resize_real_array2( hecMESH%n_node*hecMESH%n_dof, size(fstrDYNAMIC%ACC(1,:)), fstrDYNAMIC%ACC )
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, fstrDYNAMIC%VEC1 )
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, fstrDYNAMIC%VEC2 )
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, fstrDYNAMIC%VEC3 )
  end subroutine

  subroutine resize_structures_static( hecMESH, hecMAT, fstrSOLID )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid), intent(inout)      :: fstrSOLID   !< type fstr_solid

    type(fstr_solid_physic_val), pointer :: solid_physic_val

    if (associated(fstrSOLID%SOLID) ) then
      solid_physic_val => fstrSOLID%SOLID
    elseif (associated(fstrSOLID%SHELL) ) then
      solid_physic_val => fstrSOLID%SHELL
    endif

    ! resize fstrSOLID
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%QFORCE )
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%GL )
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%EFORCE )
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%REACTION )
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%CONT_NFORCE )
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%CONT_FRIC )
    call resize_real_array( 3*hecMESH%n_node, fstrSOLID%CONT_RELVEL )
    call resize_real_array( hecMESH%n_node, fstrSOLID%CONT_STATE )
    call resize_real_array( 6*hecMESH%n_node, solid_physic_val%STRESS )
    call resize_real_array( 6*hecMESH%n_node, solid_physic_val%STRAIN )
    call resize_real_array( hecMESH%n_node, solid_physic_val%MISES )
    call resize_real_array( 6*hecMESH%n_elem, solid_physic_val%ESTRESS )
    call resize_real_array( 6*hecMESH%n_elem, solid_physic_val%ESTRAIN )
    call resize_real_array( hecMESH%n_elem, solid_physic_val%EMISES )
    call resize_real_array( 12*hecMESH%n_elem, solid_physic_val%ENQM )
    fstrSOLID%STRAIN => solid_physic_val%STRAIN
    fstrSOLID%STRESS => solid_physic_val%STRESS
    fstrSOLID%MISES  => solid_physic_val%MISES
    fstrSOLID%ESTRAIN => solid_physic_val%ESTRAIN
    fstrSOLID%ESTRESS => solid_physic_val%ESTRESS
    fstrSOLID%EMISES  => solid_physic_val%EMISES
    fstrSOLID%ENQM    => solid_physic_val%ENQM
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, fstrSOLID%unode )
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, fstrSOLID%dunode )
    call resize_real_array( hecMESH%n_dof*hecMESH%n_node, fstrSOLID%ddunode )
    hecMAT%NP = hecMESH%n_node
    call resize_real_array( 3*hecMESH%n_node, hecMAT%X )
    call resize_real_array( 3*hecMESH%n_node, hecMAT%B )
  end subroutine

end module mMPCPreprocess
