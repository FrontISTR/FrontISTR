!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains control file data obtaining functions for dynamic analysis

module hecmw_es_mesh_connectivity
  use hecmw_util
  use hecmw_varray_int
  use hecmw_etype
  implicit none

  private
  public :: hecmw_create_smoothing_element_connectivity

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  common functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine concat_int_list(inlist, n_in, outlist)
    integer, allocatable, intent(in) :: inlist(:)
    integer, intent(in)              :: n_in
    integer, pointer, intent(inout)  :: outlist(:)

    integer :: i, n_ori, n_new, lbd
    integer, allocatable :: tmp(:)

    n_ori = size(outlist)
    lbd = lbound(outlist,1)
    allocate(tmp(n_ori))

    do i=1,n_ori
      tmp(i) = outlist(lbd+i-1)
    end do

    n_new = n_ori + n_in

    deallocate(outlist)
    allocate(outlist(lbd:lbd+n_new-1))
    do i=1,n_ori
      outlist(lbd+i-1) = tmp(i)
    end do
    do i=1,n_in
      outlist(lbd+i-1+n_ori) = inlist(i)
    end do

    deallocate(tmp)
  end subroutine

  subroutine add_newelements_to_hecmesh( etype, elems, sectionID, hecMESH )
    integer, intent(in)                       :: etype
    type( hecmwST_varray_int ), allocatable, intent(in) :: elems(:)
    integer, allocatable, intent(in)          :: sectionID(:)
    type(hecmwST_local_mesh), intent(inout)   :: hecMESH

    integer(kind=kint) :: i
    integer(kind=kint) :: n_newelem, n_newelemnodes
    integer(kind=kint), allocatable :: iwork(:)

    n_newelem = size(elems)
    n_newelemnodes = 0
    do i=1,n_newelem
      n_newelemnodes = n_newelemnodes + elems(i)%nitem
    end do
    allocate(iwork(0:n_newelemnodes))

    ! elem_type_index(:)
    iwork(1) = hecMESH%elem_type_index(hecMESH%n_elem_type)+n_newelem
    call concat_int_list(iwork, 1, hecMESH%elem_type_index)

    ! elem_type_item(:)
    iwork(1) = etype
    call concat_int_list(iwork, 1, hecMESH%elem_type_item)

    ! elem_type(:)
    iwork(1:n_newelem) = etype
    call concat_int_list(iwork, n_newelem, hecMESH%elem_type)

    ! section_ID(:)
    call concat_int_list(sectionID, n_newelem, hecMESH%section_ID)

    ! elem_mat_ID_index(:)
    iwork(0) = hecMESH%elem_mat_ID_index(ubound(hecMESH%elem_mat_ID_index,1))
    do i=1,n_newelem
      iwork(i) = iwork(i-1) + 1
    end do
    call concat_int_list(iwork, n_newelem, hecMESH%elem_mat_ID_index)

    ! elem_mat_ID_item(:)

    ! elem_node_index(:)
    iwork(0) = hecMESH%elem_node_index(ubound(hecMESH%elem_node_index,1))
    do i=1,n_newelem
      iwork(i) = iwork(i-1) + elems(i)%nitem
    end do
    call concat_int_list(iwork, n_newelem, hecMESH%elem_node_index)

    ! elem_node_item(:)
    n_newelemnodes = 0
    do i=1,n_newelem
      iwork(n_newelemnodes+1:n_newelemnodes+elems(i)%nitem) = elems(i)%items(1:elems(i)%nitem)
      n_newelemnodes = n_newelemnodes + elems(i)%nitem
    end do
    call concat_int_list(iwork, n_newelemnodes, hecMESH%elem_node_item)

    ! elem_ID(:)
    do i=1,n_newelem
      iwork(2*i-1) = hecMESH%n_elem + 1
      iwork(2*i  ) = 0
    end do
    call concat_int_list(iwork, 2*n_newelem, hecMESH%elem_ID)

    ! global_elem_ID(:)
    do i=1,n_newelem
      iwork(i) = 0
    end do
    call concat_int_list(iwork, n_newelem, hecMESH%global_elem_ID)

    ! elem_internal_list(:)
    ! elem_mat_int_index(:)
    ! elem_mat_int_val(:)
    ! elem_val_index(:)
    ! elem_val_item(:)

    ! n_elem
    !hecMESH%n_elem = hecMESH%n_elem + n_newelem

    ! n_elem_gross
    ! ne_internal
    !hecMESH%ne_internal = hecMESH%ne_internal + n_newelem

    ! n_elem_type
    hecMESH%n_elem_type = hecMESH%n_elem_type + 1

    ! n_elem_mat_ID

  end subroutine

  subroutine reorder_tet1( nid, nodlocal )
    integer, intent(in)    :: nid
    integer, intent(inout) :: nodlocal(4)

    integer :: i, idx
    integer :: nodlocal_ori(4)

    nodlocal_ori(1:4) = nodlocal(1:4)
    idx = 0
    do i=1,4
      if( nodlocal(i) == nid ) then
        idx = i
        exit
      end if
    end do

    if( idx == 1 ) then
      return
    else if( idx == 2 ) then
      nodlocal(1) = nodlocal_ori(2)
      nodlocal(2) = nodlocal_ori(3)
      nodlocal(3) = nodlocal_ori(1)
    else if( idx == 3 ) then
      nodlocal(1) = nodlocal_ori(3)
      nodlocal(2) = nodlocal_ori(1)
      nodlocal(3) = nodlocal_ori(2)
    else if( idx == 4 ) then
      nodlocal(1) = nodlocal_ori(4)
      nodlocal(2) = nodlocal_ori(2)
      nodlocal(3) = nodlocal_ori(1)
      nodlocal(4) = nodlocal_ori(3)
    end if

  end subroutine

  subroutine reorder_tet1_edge( nid1, nid2, nodlocal )
    integer, intent(in)    :: nid1
    integer, intent(in)    :: nid2
    integer, intent(inout) :: nodlocal(4)

    integer :: i, idx
    integer :: nodlocal_ori(4)

    call reorder_tet1( nid1, nodlocal )

    nodlocal_ori(1:4) = nodlocal(1:4)
    idx = 1
    do i=2,4
      if( nodlocal(i) == nid2 ) then
        idx = i
        exit
      end if
    end do

    if( idx == 2 ) then
      return
    else if( idx == 3 ) then
      nodlocal(2) = nodlocal_ori(3)
      nodlocal(3) = nodlocal_ori(4)
      nodlocal(4) = nodlocal_ori(2)
    else if( idx == 4 ) then
      nodlocal(2) = nodlocal_ori(4)
      nodlocal(3) = nodlocal_ori(2)
      nodlocal(4) = nodlocal_ori(3)
    else
      stop "error in reorder_tet1_edge"
    end if

  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  create node-smoothed elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_node_elem_list( hecMESH, n_elem_ori, is_selem_list, node_elem_list )
    type(hecmwST_local_mesh),target              :: hecMESH
    integer(kind=kint), intent(in)               :: n_elem_ori
    logical, allocatable, intent(in)             :: is_selem_list(:)
    type( hecmwST_varray_int ), allocatable, intent(inout) :: node_elem_list(:)

    integer(kind=kint) :: i, icel, iS, iE, nid

    call HECMW_varray_int_initialize_all( hecMESH%n_node, 8, node_elem_list )

    do icel=1,n_elem_ori
      if( .not. is_selem_list(icel) ) cycle

      iS = hecMESH%elem_node_index(icel-1)+1
      iE = hecMESH%elem_node_index(icel)

      do i=iS,iE
        nid = hecMESH%elem_node_item(i)
        call HECMW_varray_int_add( icel, node_elem_list(nid) )
      end do
    end do

  end subroutine

  subroutine create_node_smoothed_elements( hecMESH, node_elem_list, nselems, sectionID )
    type(hecmwST_local_mesh), intent(in)         :: hecMESH
    type( hecmwST_varray_int ), allocatable, intent(in)    :: node_elem_list(:)
    type( hecmwST_varray_int ), allocatable, intent(inout) :: nselems(:)
    integer, allocatable, intent(inout)          :: sectionID(:)

    integer(kind=kint) :: n_nselems
    type( hecmwST_varray_int ), allocatable :: n_sections(:)
    integer(kind=kint) :: i, j, k, icel, isect, iS, iE
    integer(kind=kint) :: nodlocal(4)

    !get num of node-smoothed elements
    call HECMW_varray_int_initialize_all( hecMESH%n_node, 2, n_sections )

    n_nselems = 0
    do i=1,hecMESH%n_node
      do j=1,node_elem_list(i)%nitem
        icel = node_elem_list(i)%items(j)
        isect = hecMESH%section_ID(icel)
        call HECMW_varray_int_add_if_not_exits( isect, n_sections(i) )
      enddo
      n_nselems = n_nselems + n_sections(i)%nitem
    end do

    allocate(sectionID(n_nselems))

    !create node-smoothed elements
    !elements which belongs to different sections must not be connected
    call HECMW_varray_int_initialize_all( n_nselems, 16, nselems )

    n_nselems = 0
    do i=1,hecMESH%n_node
      do j=1,n_sections(i)%nitem
        isect = n_sections(i)%items(j)
        n_nselems = n_nselems + 1
        sectionID(n_nselems) = isect
        do k=1,node_elem_list(i)%nitem
          icel = node_elem_list(i)%items(k)
          if ( hecMESH%section_ID(icel) == isect ) then
            iS = hecMESH%elem_node_index(icel-1)+1
            iE = hecMESH%elem_node_index(icel)
            if( iE-iS+1 /= 4 ) stop "error in add_elemments_smoothed_by_node(1)"
            nodlocal(1:4) = hecMESH%elem_node_item(is:iE)
            call reorder_tet1( i, nodlocal )
            if( nselems(n_nselems)%nitem == 0 ) call HECMW_varray_int_add( nodlocal(1), nselems(n_nselems) )
            call HECMW_varray_int_expand( 3, nodlocal(2:4), nselems(n_nselems) )
          end if
        end do
      enddo
    end do

    !finalize n_sections
    call HECMW_varray_int_finalize_all( n_sections )

  end subroutine

  subroutine add_elemments_smoothed_by_node( hecMESH, n_elem_ori, is_selem_list )
    type(hecmwST_local_mesh),target  :: hecMESH
    integer(kind=kint), intent(in)   :: n_elem_ori
    logical, allocatable, intent(in) :: is_selem_list(:)

    type( hecmwST_varray_int ), allocatable :: node_elem_list(:)
    type( hecmwST_varray_int ), allocatable :: nselems(:)
    integer, allocatable :: sectionID(:)


    !initialize node_elem_list
    call initialize_node_elem_list( hecMESH, n_elem_ori, is_selem_list, node_elem_list )

    !create node-smoothed elements and store them in list data
    call create_node_smoothed_elements( hecMESH, node_elem_list, nselems, sectionID )

    !finalize node_elem_list
    call HECMW_varray_int_finalize_all( node_elem_list )

    !add node-smoothed elements to hecMESH
    call add_newelements_to_hecmesh( 881, nselems, sectionID, hecMESH )

    !finalize nselems
    call HECMW_varray_int_finalize_all( nselems )

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  create edge-smoothed elements
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_edge_elem_list( hecMESH, n_elem_ori, is_selem_list, edges, edge_elem_list )
    type(hecmwST_local_mesh),target              :: hecMESH
    integer(kind=kint), intent(in)               :: n_elem_ori
    logical, allocatable, intent(in)             :: is_selem_list(:)
    integer, allocatable, intent(inout)          :: edges(:)
    type( hecmwST_varray_int ), allocatable, intent(inout) :: edge_elem_list(:)

    type( hecmwST_varray_int ), allocatable :: edge_con(:)
    type( hecmwST_varray_int ), allocatable :: edge_id(:)
    integer(kind=kint) :: ndlocal(4)
    integer(kind=kint), parameter :: leid(2,6) = reshape((/1,2,2,3,3,1,1,4,2,4,3,4/),shape(leid))

    integer(kind=kint) :: i, j, icel, iS, iE, nid1, nid2
    integer(kind=kint) :: n_edges, edid

    !create edge connectivity list in node num order
    !edge_con(i)%items : list of nodes connecting to i by edge
    call HECMW_varray_int_initialize_all( hecMESH%n_node, 8, edge_con )

    do icel=1,n_elem_ori
      if( .not. is_selem_list(icel) ) cycle

      iS = hecMESH%elem_node_index(icel-1)+1
      iE = hecMESH%elem_node_index(icel)

      ndlocal(1:4) = hecMESH%elem_node_item(iS:iE)
      do i=1,6
        nid1 = ndlocal(leid(1,i))
        nid2 = ndlocal(leid(2,i))
        if( nid1 < nid2 ) then
          call HECMW_varray_int_add_if_not_exits( nid2, edge_con(nid1) )
        else
          call HECMW_varray_int_add_if_not_exits( nid1, edge_con(nid2) )
        end if
      end do
    end do

    !create edges and edge_id
    n_edges = 0
    do i=1,hecMESH%n_node
      n_edges = n_edges + edge_con(i)%nitem
    end do

    allocate(edges(2*n_edges))
    call HECMW_varray_int_initialize_all( hecMESH%n_node, 8, edge_id )
    n_edges = 0
    do i=1,hecMESH%n_node
      do j=1,edge_con(i)%nitem
        n_edges = n_edges + 1
        edges(2*n_edges-1) = i
        edges(2*n_edges  ) = edge_con(i)%items(j)
        call HECMW_varray_int_add( n_edges, edge_id(i) )
      end do
    end do

    !create edge-elem list
    call HECMW_varray_int_initialize_all( n_edges, 8, edge_elem_list )

    do icel=1,n_elem_ori
      if( .not. is_selem_list(icel) ) cycle

      iS = hecMESH%elem_node_index(icel-1)+1
      iE = hecMESH%elem_node_index(icel)

      ndlocal(1:4) = hecMESH%elem_node_item(iS:iE)
      do i=1,6
        nid1 = ndlocal(leid(1,i))
        nid2 = ndlocal(leid(2,i))
        if( nid1 < nid2 ) then
          edid = edge_id(nid1)%items( HECMW_varray_int_find( nid2, edge_con(nid1) ) )
        else
          edid = edge_id(nid2)%items( HECMW_varray_int_find( nid1, edge_con(nid2) ) )
        end if
        call HECMW_varray_int_add( icel, edge_elem_list(edid) )
      end do
    end do

    call HECMW_varray_int_finalize_all( edge_con )
    call HECMW_varray_int_finalize_all( edge_id )

  end subroutine

  subroutine create_edge_smoothed_elements( hecMESH, edges, edge_elem_list, eselems, sectionID )
    type(hecmwST_local_mesh), intent(in)         :: hecMESH
    integer(kind=kint), allocatable, intent(in)  :: edges(:)
    type( hecmwST_varray_int ), allocatable, intent(in)    :: edge_elem_list(:)
    type( hecmwST_varray_int ), allocatable, intent(inout) :: eselems(:)
    integer, allocatable, intent(inout)          :: sectionID(:)

    integer(kind=kint) :: n_edges, n_eselems
    type( hecmwST_varray_int ), allocatable :: n_sections(:)
    integer(kind=kint) :: i, j, k, icel, isect, iS, iE
    integer(kind=kint) :: nodlocal(4)

    n_edges = size(edge_elem_list)

    !get num of edge-smoothed elements
    call HECMW_varray_int_initialize_all( n_edges, 2, n_sections )

    n_eselems = 0
    do i=1,n_edges
      do j=1,edge_elem_list(i)%nitem
        icel = edge_elem_list(i)%items(j)
        isect = hecMESH%section_ID(icel)
        call HECMW_varray_int_add_if_not_exits( isect, n_sections(i) )
      enddo
      n_eselems = n_eselems + n_sections(i)%nitem
    end do

    allocate(sectionID(n_eselems))

    !create edge-smoothed elements
    !elements which belongs to different sections must not be connected
    call HECMW_varray_int_initialize_all( n_eselems, 16, eselems )

    n_eselems = 0
    do i=1,n_edges
      do j=1,n_sections(i)%nitem
        isect = n_sections(i)%items(j)
        n_eselems = n_eselems + 1
        sectionID(n_eselems) = isect
        do k=1,edge_elem_list(i)%nitem
          icel = edge_elem_list(i)%items(k)
          if ( hecMESH%section_ID(icel) == isect ) then
            iS = hecMESH%elem_node_index(icel-1)+1
            iE = hecMESH%elem_node_index(icel)
            if( iE-iS+1 /= 4 ) stop "error in add_elemments_smoothed_by_node(1)"
            nodlocal(1:4) = hecMESH%elem_node_item(is:iE)
            call reorder_tet1_edge( edges(2*i-1), edges(2*i), nodlocal )
            if( eselems(n_eselems)%nitem == 0 ) call HECMW_varray_int_expand( 2, nodlocal(1:2), eselems(n_eselems) )
            call HECMW_varray_int_expand( 2, nodlocal(3:4), eselems(n_eselems) )
          end if
        end do
      enddo
    end do

    !finalize n_sections
    call HECMW_varray_int_finalize_all( n_sections )

  end subroutine

  subroutine add_elemments_smoothed_by_edge( hecMESH, n_elem_ori, is_selem_list )
    type(hecmwST_local_mesh),target :: hecMESH
    integer(kind=kint), intent(in)  :: n_elem_ori
    logical, allocatable, intent(in) :: is_selem_list(:)

    integer, allocatable          :: edges(:)
    type( hecmwST_varray_int ), allocatable :: edge_elem_list(:)
    type( hecmwST_varray_int ), allocatable :: eselems(:)
    integer, allocatable :: sectionID(:)


    !initialize edge_elem_list
    call initialize_edge_elem_list( hecMESH, n_elem_ori, is_selem_list, edges, edge_elem_list )

    !create edge-smoothed elements and store them in list data
    call create_edge_smoothed_elements( hecMESH, edges, edge_elem_list, eselems, sectionID )

    !finalize edges and edge_elem_list
    deallocate(edges)
    call HECMW_varray_int_finalize_all( edge_elem_list )

    !add edge-smoothed elements to hecMESH
    call add_newelements_to_hecmesh( 891, eselems, sectionID, hecMESH )

    !finalize eselems
    call HECMW_varray_int_finalize_all( eselems )

  end subroutine

  !> setup selective es/ns smoothing FEM connectivity
  subroutine hecmw_create_smoothing_element_connectivity( hecMESH, is_selem_list )
    type(hecmwST_local_mesh),target :: hecMESH
    logical, allocatable, intent(in) :: is_selem_list(:)

    integer(kind=kint) :: n_elem_ori

    n_elem_ori = hecMESH%n_elem

    call add_elemments_smoothed_by_node( hecMESH, n_elem_ori, is_selem_list )

    call add_elemments_smoothed_by_edge( hecMESH, n_elem_ori, is_selem_list )

  end subroutine

end module
