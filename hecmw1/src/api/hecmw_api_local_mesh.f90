!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_api_local_mesh
  use iso_c_binding
  use hecmw_util, only : hecmwST_local_mesh
  implicit none

contains

  !> @brief Generating mesh handler
  !! @return Mesh handler of type(c_ptr)
  function hecmw_api_mesh_new() bind(C,name='hecmw_api_mesh_new')
    use hecmw_util, only : hecmw_nullify_mesh
    implicit none
    type(c_ptr) :: hecmw_api_mesh_new
    type(hecmwST_local_mesh), pointer :: hecMESH => null()
    allocate(hecMESH)
    call hecmw_nullify_mesh(hecMESH)
    allocate( hecMESH%node_group%grp_index(0:0) )
    hecMESH%node_group%grp_index(0) = 0
    allocate( hecMESH%surf_group%grp_index(0:0) )
    hecMESH%surf_group%grp_index(0) = 0
    allocate( hecMESH%elem_group%grp_index(0:0) )
    hecMESH%elem_group%grp_index(0) = 0
    hecmw_api_mesh_new = c_loc(hecMESH)
  end function

  !> @brief Destroying the mesh handler
  !! @param[in] mesh : Mesh handler
  subroutine hecmw_api_mesh_delete(mesh) bind(C,name='hecmw_api_mesh_delete')
    use hecmw_dist_free_f, only : hecmw_dist_free
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call hecmw_dist_free(hecMESH)
    deallocate(hecMESH)
  end subroutine

  !> @brief Setting up the node array
  !! @param[in] mesh : Mesh handler
  !! @param[in] nnode : Number of nodes
  !! @param[in] node : Node array
  subroutine hecmw_api_mesh_set_node(mesh,nnode,node) bind(C,name='hecmw_api_mesh_set_node')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: nnode
    real(c_double), intent(in) :: node(3*nnode)

    integer :: i
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=mesh, fptr=hecMESH)

    hecMESH%n_node = nnode
    hecMESH%n_node_gross = nnode
    hecMESH%nn_middle = nnode
    hecMESH%nn_internal = nnode
    allocate(hecMESH%node(3*nnode))
    allocate(hecMESH%node_ID(2*nnode))
    allocate(hecMESH%node_internal_list(nnode))
    allocate(hecMESH%global_node_ID(nnode))

    do i=1, nnode
      hecMESH%node_internal_list(i) = i
      hecMESH%node_ID(2*i-1) = i
      hecMESH%node_ID(2*i) = 0
      hecMESH%global_node_ID(i) = i
      hecMESH%node(3*i-2) = node(3*i-2)
      hecMESH%node(3*i-1) = node(3*i-1)
      hecMESH%node(3*i-0) = node(3*i-0)
    enddo

  end subroutine

  !> @brief Retrieving the number of nodes.
  !! @param[in] mesh : Mesh handler
  !! @return Number of nodes
  function hecmw_api_mesh_n_node(mesh) bind(C,name='hecmw_api_mesh_n_node')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_node
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_node = hecMESH%n_node
  end function

  !> @brief Retrieving the array of nodal coordinates.
  !! @param[in] mesh : Mesh handler
  !! @return Pointer to an array of nodes
  !! @remark The recipient must not deallocate the pointer obtained from this function.
  function hecmw_api_mesh_get_node(mesh) bind(C,name='hecmw_api_mesh_get_node')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_node
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_node = c_loc(hecMESH%node)
  end function

  !> @brief Setting of Nodal Degrees of Freedom
  !! @param[in] mesh : Mesh handler
  !! @param[in] ndof : Nodal degrees of freedom
  subroutine hecmw_api_mesh_set_n_dof(mesh,ndof) bind(C,name='hecmw_api_mesh_set_n_dof')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: ndof
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecMESH%n_dof = ndof
  end subroutine

  !> @brief Retrieving nodal degrees of freedom
  !! @param[in] mesh : Mesh handler
  !! @return Nodal degrees of freedom
  function hecmw_api_mesh_get_n_dof(mesh) bind(C,name='hecmw_api_mesh_get_n_dof')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) :: hecmw_api_mesh_get_n_dof
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_n_dof = hecMESH%n_dof
  end function

  !> @brief Element Settings
  !! @param[in] mesh : Mesh handler
  !! @param[in] nelem : Number of elements
  !! @param[in] elemtype : Array of element types
  !! @param[in] element : Array of node tables (1-indexed)
  !! @param[in] sectionID : Array of section numbers (1-indexed)
  subroutine hecmw_api_mesh_set_element(mesh,nelem,elemtype,element,sectionID) bind(C,name='hecmw_api_mesh_set_element')
    use hecmw_etype, only : hecmw_get_max_node
    implicit none
    type(c_ptr), value                :: mesh
    integer(c_int), value, intent(in) :: nelem
    integer(c_int), intent(in)        :: elemtype(nelem)
    integer(c_int), intent(in)        :: element(*)
    integer(c_int), intent(in)        :: sectionID(nelem)

    type(hecmwST_local_mesh), pointer :: hecMESH
    integer :: i, j, ncon, n, ii, off
    integer, parameter :: etypes(35) = [ &
      111, 112, &
      231, 232, 2322, 241, 242, &
      301, 341, 3414, 342, 3422, 351, 352, 361, 362, 363, &
      511, &
      611, 612, 641, &
      732, 733, 731, 741, 742, 743, 761, 781, &
      881, 891, &
      1031, 1032, 1041, 1042 &
    ]
    integer :: etype_counter(35)

    call c_f_pointer(cptr=mesh, fptr=hecMESH)

    hecMESH%n_elem = nelem
    hecMESH%n_elem_gross = nelem
    hecMESH%ne_internal = nelem
    allocate(hecMESH%elem_internal_list(nelem))
    allocate(hecMESH%elem_ID(2*nelem))
    allocate(hecMESH%global_elem_ID(nelem))
    allocate(hecMESH%elem_node_index(0:nelem))
    allocate(hecMESH%elem_type(nelem))
    allocate(hecMESH%section_ID(nelem))

    ! Determine the count for each element type.
    etype_counter(:) = 0
    do i=1, nelem
      do j=1, 35
        if (elemtype(i) == etypes(j)) then
          etype_counter(j) = etype_counter(j) + 1
        end if
      end do
    end do
    ! Determine how many types exist.
    n = 0
    do j=1, 35
      if (etype_counter(j) > 0) then
        n = n + 1
      end if
    end do

    hecMESH%n_elem_type = n
    allocate(hecMESH%elem_type_index(0:n))
    allocate(hecMESH%elem_type_item(n))

    ! Record the number of elements for each type.
    i = 1
    n = 0
    hecMESH%elem_type_index(0) = n
    do j=1, 35
      if (etype_counter(j) > 0) then
        n = n + etype_counter(j)
        hecMESH%elem_type_item(i) = etypes(j)
        hecMESH%elem_type_index(i) = n
        i = i + 1
      end if
    end do

    ! Determine the number of elem_node_item entries.
    ncon = 0
    do i=1, nelem
      n = HECMW_get_max_node(elemtype(i))
      ncon = ncon + n
    end do
    allocate(hecMESH%elem_node_item(ncon))

    ! Sort by element type and assign
    ncon = 0 ! Offset of elem_node_item
    hecMESH%elem_node_index(0) = ncon
    ii = 0
    do j=1, 35
      if (etype_counter(j) == 0) continue
      n = HECMW_get_max_node(etypes(j))
      off = 0 ! Offset of the input array
      do i=1, nelem
        if (elemtype(i)==etypes(j)) then
          hecMESH%elem_node_item(ncon+1:ncon+n) = element(off+1:off+n)
          ncon = ncon + n
          ii = ii + 1
          hecMESH%elem_node_index(ii)    = ncon
          hecMESH%elem_ID(2*ii-1)        = i
          hecMESH%elem_ID(2*ii)          = 0
          hecMESH%global_elem_ID(ii)     = i
          hecMESH%elem_internal_list(ii) = i
          hecMESH%elem_type(ii)          = elemtype(i)
          hecMESH%section_ID(ii)         = sectionID(i)
        end if
        off = off + HECMW_get_max_node(elemtype(i))
      end do
    end do

    ! Examine n_sect
    hecMESH%section%n_sect = 0
    do i=1, nelem
      if (hecMESH%section%n_sect < hecMESH%section_ID(i)) then
        hecMESH%section%n_sect = hecMESH%section_ID(i)
      end if
    end do

  end subroutine

  !> @brief Retrieving element types
  !! @param[in] mesh : Mesh handler
  !! @return Pointer to an array of elements
  function hecmw_api_mesh_get_elem_type(mesh) bind(C,name='hecmw_api_mesh_get_elem_type')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_elem_type
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_elem_type = c_loc(hecMESH%elem_type(1))
  end function

  !> @brief Retrieving the node table
  !! @param[in] mesh : Mesh handler
  !! @return Pointer to the node table array
  function hecmw_api_mesh_get_elem_node_item(mesh) bind(C,name='hecmw_api_mesh_get_elem_node_item')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_elem_node_item
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_elem_node_item = c_loc(hecMESH%elem_node_item(1))
  end function

  !> @brief Retrieving the section number array
  !! @param[in] mesh : Mesh handler
  !! @return Pointer to an array of section numbers
  function hecmw_api_mesh_get_section_ID(mesh) bind(C,name='hecmw_api_mesh_get_section_id')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_section_ID
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_section_ID = c_loc(hecMESH%section_ID(1))
  end function

  !> @brief Getting the number of elements
  !! @param[in] mesh : Mesh handler
  !! @return Number of elements
  function hecmw_api_mesh_n_elem(mesh) bind(C,name='hecmw_api_mesh_n_elem')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_elem
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_elem = hecMESH%n_elem
  end function

  !> @brief Retrieving the size of the element-node array
  !! @param[in] mesh : Mesh handler
  !! @return Size of the element-node array
  function hecmw_api_mesh_n_elem_node_item(mesh) bind(C,name='hecmw_api_mesh_n_elem_node_item')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_elem_node_item
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_elem_node_item = hecMESH%elem_node_index(hecMESH%n_elem)
  end function

  !> @brief Number of node groups
  !! @param[in] mesh : Mesh handler
  !! @return Number of groups
  !! @remark NALL for all nodes is not included unless explicitly added.
  function hecmw_api_mesh_n_ngrp(mesh) bind(C,name='hecmw_api_mesh_n_ngrp')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_ngrp
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_ngrp = hecMESH%node_group%n_grp
  end function

  !> @brief Add node group
  !! @param[in] mesh : Mesh handler
  !! @param[in] grp_name : Node group name
  !! @param[in] count : Number of nodes
  !! @param[in] list : Array of node numbers
  subroutine hecmw_api_mesh_append_ngrp(mesh,grp_name,count,list) bind(C,name='hecmw_api_mesh_append_ngrp')
    use hecmw_api_common, only : c_f_str_copy
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw, only: kint
    use hecmw_setup_util, only : append_new_group 
    implicit none
    type(c_ptr), value, intent(in) :: mesh
    type(c_ptr), value, intent(in) :: grp_name
    integer(c_int), value, intent(in) :: count
    integer(c_int), intent(in) :: list(count)
    type(hecmwST_local_mesh), pointer :: hecMESH
    character(kind=c_char), pointer :: char_array(:) => null()
    character(len=HECMW_NAME_LEN) :: grp_name_f
    integer(kind=kint) :: grp_id

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call c_f_pointer(grp_name, char_array, [HECMW_NAME_LEN])
    call c_f_str_copy(char_array, grp_name_f, HECMW_NAME_LEN)
    call append_new_group(hecMESH, 'node_grp', grp_name_f, count, list, grp_id)
  end subroutine

  !> @brief Node group name
  !! @param[in] mesh : Mesh handler
  !! @param[in] i : Index of the node group (in order of addition)
  !! @param[inout] buf : Node group name
  !! @param[in] buflen : Size of the allocated buffer
  subroutine hecmw_api_mesh_get_ngrp_name(mesh,i,buf,buflen) bind(C,name='hecmw_api_mesh_get_ngrp_name')
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int), value, intent(in) :: i
    character(kind=c_char), intent(inout) :: buf(*)
    integer(c_int), value, intent(in) :: buflen
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call f_c_str_copy(hecMESH%node_group%grp_name(i),buf,buflen)
  end subroutine

  !> @brief Array of node numbers for a node group
  !! @param[in] mesh : Mesh handler
  !! @param[in] i : Index of the node group (in order of addition)
  !! @param[out] array : Pointer to the node number array
  !! @param[out] count : Array size
  subroutine hecmw_api_mesh_get_ngrp(mesh,i,array,count) bind(C,name='hecmw_api_mesh_get_ngrp')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: i
    type(c_ptr), intent(out) :: array
    integer(c_int), intent(out) :: count
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer :: is, ie

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    is = hecMESH%node_group%grp_index(i-1)+1
    ie = hecMESH%node_group%grp_index(i)
    array = c_loc(hecMESH%node_group%grp_item(is))
    count = ie - is + 1
  end subroutine

  !> @brief Number of face groups
  !! @param[in] mesh : Mesh hander
  !! @return Number of groups
  function hecmw_api_mesh_n_sgrp(mesh) bind(C,name='hecmw_api_mesh_n_sgrp')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_sgrp
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_sgrp = hecMESH%surf_group%n_grp
  end function

  !> @brief Add surface group
  !! @param[in] mesh : Mesh handler
  !! @param[in] grp_name : Name of the surface group
  !! @param[in] count : Number of surfaces
  !! @param[in] list : Array of surface indices
  subroutine hecmw_api_mesh_append_sgrp(mesh,grp_name,count,list) bind(C,name='hecmw_api_mesh_append_sgrp')
    use hecmw_api_common, only : c_f_str_copy
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw, only: kint
    use hecmw_setup_util, only : append_new_group 
    implicit none
    type(c_ptr), value, intent(in) :: mesh
    type(c_ptr), value, intent(in) :: grp_name
    integer(c_int), value, intent(in) :: count
    integer(c_int), intent(in) :: list(count)
    type(hecmwST_local_mesh), pointer :: hecMESH
    character(kind=c_char), pointer :: char_array(:) => null()
    character(len=HECMW_NAME_LEN) :: grp_name_f
    integer(kind=kint) :: grp_id

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call c_f_pointer(grp_name, char_array, [HECMW_NAME_LEN])
    call c_f_str_copy(char_array, grp_name_f, HECMW_NAME_LEN)

    call append_new_group(hecMESH, 'surf_grp', grp_name_f, count, list, grp_id)
  end subroutine

  !> @brief Name of the surface group
  !! @param[in] mesh : Mesh handler
  !! @param[in] i : Surface group index (in order of addition)
  !! @param[inout] buf : Name of the surface group
  !! @param[in] buflen : Size of the allocated buffer
  subroutine hecmw_api_mesh_get_sgrp_name(mesh,i,buf,buflen) bind(C,name='hecmw_api_mesh_get_sgrp_name')
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int), value, intent(in) :: i
    character(kind=c_char), intent(inout) :: buf(*)
    integer(c_int), value, intent(in) :: buflen
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call f_c_str_copy(hecMESH%surf_group%grp_name(i),buf,buflen)
  end subroutine

  !> @brief Array of element indices and face indices for a face group
  !! @param[in] mesh : Mesh handler
  !! @param[in] i : Surface group index (in order of addition)
  !! @param[out] array : Pointer to an array of element numbers and face numbers
  !! @param[out] count : Array size
  subroutine hecmw_api_mesh_get_sgrp(mesh,i,array,count) bind(C,name='hecmw_api_mesh_get_sgrp')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: i
    type(c_ptr), intent(out) :: array
    integer(c_int), intent(out) :: count
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer :: is, ie

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    is = hecMESH%surf_group%grp_index(i-1)+1
    ie = hecMESH%surf_group%grp_index(i)
    array = c_loc(hecMESH%surf_group%grp_item(is))
    count = ie - is + 1
  end subroutine

  !> @brief Number of element groups
  !! @param[in] mesh : Mesh handler
  !! @return Number of groups
  !! @remark EALL for all elements is not included unless explicitly added.
  function hecmw_api_mesh_n_egrp(mesh) bind(C,name='hecmw_api_mesh_n_egrp')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_egrp
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_egrp = hecMESH%elem_group%n_grp
  end function

  !> @brief Add element group
  !! @param[in] mesh : Mesh handler
  !! @param[in] grp_name : Name of the element group
  !! @param[in] count : Number of elements
  !! @param[in] list : Array of element numbers
  subroutine hecmw_api_mesh_append_egrp(mesh,grp_name,count,list) bind(C,name='hecmw_api_mesh_append_egrp')
    use hecmw_api_common, only : c_f_str_copy
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw, only: kint
    use hecmw_setup_util, only : append_new_group 
    implicit none
    type(c_ptr), value, intent(in) :: mesh
    type(c_ptr), value, intent(in) :: grp_name
    integer(c_int), value, intent(in) :: count
    integer(c_int), intent(in) :: list(count)
    type(hecmwST_local_mesh), pointer :: hecMESH
    character(kind=c_char), pointer :: char_array(:) => null()
    character(len=HECMW_NAME_LEN) :: grp_name_f
    integer(kind=kint) :: grp_id

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call c_f_pointer(grp_name, char_array, [HECMW_NAME_LEN])
    call c_f_str_copy(char_array, grp_name_f, HECMW_NAME_LEN)

    call append_new_group(hecMESH, 'elem_grp', grp_name_f, count, list, grp_id)
  end subroutine

  !> @brief Name of the element group
  !! @param[in] mesh : Mesh handler
  !! @param[in] i : Index of the element group (in order of addition)
  !! @param[inout] buf : Element group name
  !! @param[in] buflen : Size of the allocated buffer
  subroutine hecmw_api_mesh_get_egrp_name(mesh,i,buf,buflen) bind(C,name='hecmw_api_mesh_get_egrp_name')
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int), value, intent(in) :: i
    character(kind=c_char), intent(inout) :: buf(*)
    integer(c_int), value, intent(in) :: buflen
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call f_c_str_copy(hecMESH%elem_group%grp_name(i),buf,buflen)
  end subroutine
  
  !> @brief An array of element numbers for a group of elements
  !! @param[in] mesh : Mesh handler
  !! @param[in] i : Index of the element group (in order of addition)
  !! @param[out] array : Element number array
  !! @param[out] count : Array size
  subroutine hecmw_api_mesh_get_egrp(mesh,i,array,count) bind(C,name='hecmw_api_mesh_get_egrp')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: i
    type(c_ptr), intent(out) :: array
    integer(c_int), intent(out) :: count
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer :: is, ie

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    is = hecMESH%elem_group%grp_index(i-1)+1
    ie = hecMESH%elem_group%grp_index(i)
    array = c_loc(hecMESH%elem_group%grp_item(is))
    count = ie - is + 1
  end subroutine

end module
