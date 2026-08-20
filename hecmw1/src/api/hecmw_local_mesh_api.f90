module hecmw_local_mesh_api
  use iso_c_binding
  use hecmw_util, only : hecmwST_local_mesh
  implicit none

contains

  !> @brief メッシュハンドラの生成
  !> @return メッシュ構造体のハンドラ type(c_ptr)
  function hecmw_api_mesh_new() bind(C,name='hecmw_api_mesh_new')
    use hecmw_util, only : hecmw_nullify_mesh
    implicit none
    type(c_ptr) :: hecmw_api_mesh_new
    type(hecmwST_local_mesh), target, save :: hecMESH
    call hecmw_nullify_mesh(hecMESH)
    hecmw_api_mesh_new = c_loc(hecMESH)
  end function

  !> @brief メッシュハンドラの破棄
  !> @param[in] mesh メッシュ構造体のハンドラ
  subroutine hecmw_api_mesh_delete(mesh) bind(C,name='hecmw_api_mesh_delete')
    use hecmw_dist_free_f, only : hecmw_dist_free
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call hecmw_dist_free(hecMESH)
  end subroutine

  !> @brief 節点配列の設定
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] nnode 節点数
  !> @param[in] node 節点配列
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

  !> @brief 節点数の取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 節点数
  function hecmw_api_mesh_n_node(mesh) bind(C,name='hecmw_api_mesh_n_node')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_node
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_node = hecMESH%n_node
  end function

  !> @brief 節点座標配列の取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 節点配列のポインタ
  !> @details この関数で得られたポインタを受取側は解放してはならない
  function hecmw_api_mesh_get_node(mesh) bind(C,name='hecmw_api_mesh_get_node')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_node
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_node = c_loc(hecMESH%node)
  end function

  subroutine hecmw_api_mesh_set_n_dof(mesh,ndof) bind(C,name='hecmw_api_mesh_set_n_dof')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: ndof
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecMESH%n_dof = ndof
  end subroutine

  function hecmw_api_mesh_get_n_dof(mesh) bind(C,name='hecmw_api_mesh_get_n_dof')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) :: hecmw_api_mesh_get_n_dof
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_n_dof = hecMESH%n_dof
  end function

  subroutine hecmw_api_mesh_set_element(mesh,nelem,elemtype,element,sectionID) bind(C,name='hecmw_api_mesh_set_element')
    use hecmw_etype, only : hecmw_get_max_node
    implicit none
    type(c_ptr), value                :: mesh
    integer(c_int), value, intent(in) :: nelem
    integer(c_int), intent(in)        :: element(*)
    integer(c_int), intent(in)        :: elemtype(nelem)
    integer(c_int), intent(in)        :: sectionID(nelem)

    type(hecmwST_local_mesh), pointer :: hecMESH
    integer :: i, j, ncon, n, start

    call c_f_pointer(cptr=mesh, fptr=hecMESH)

    hecMESH%n_elem = nelem
    hecMESH%n_elem_gross = nelem
    hecMESH%ne_internal = nelem
    hecMESH%n_elem_type = nelem
    hecMESH%section%n_sect = 0
    allocate(hecMESH%elem_internal_list(nelem))
    allocate(hecMESH%elem_ID(2*nelem))
    allocate(hecMESH%global_elem_ID(nelem))
    allocate(hecMESH%elem_type(nelem))
    allocate(hecMESH%elem_node_index(0:nelem))
    allocate(hecMESH%elem_type_index(0:nelem))
    allocate(hecMESH%elem_type_item(nelem))
    allocate(hecMESH%section_ID(nelem))

    ncon = 0
    hecMESH%elem_node_index(0) = ncon
    hecMESH%elem_type_index(0) = 0
    do i=1, nelem
      n = HECMW_get_max_node(elemtype(i))
      ncon = ncon + n
      hecMESH%elem_node_index(i) = ncon
      hecMESH%elem_type_index(i) = i
    end do

    allocate(hecMESH%elem_node_item(ncon))
    do i=1, nelem
      n = hecMESH%elem_node_index(i) - hecMESH%elem_node_index(i-1)
      start = hecMESH%elem_node_index(i-1)
      do j=1, n
        hecMESH%elem_node_item(start+j) = element(start+j) + 1
      end do
      hecMESH%elem_ID(2*i)          = i
      hecMESH%elem_ID(2*i+1)        = 0
      hecMESH%global_elem_ID(i)     = i
      hecMESH%elem_internal_list(i) = i
      hecMESH%elem_type(i)          = elemtype(i)
      hecMESH%elem_type_item(i)     = elemtype(i)
      hecMESH%section_ID(i)         = sectionID(i)+1
      if (hecMESH%section%n_sect < hecMESH%section_ID(i)) then
        hecMESH%section%n_sect = hecMESH%section_ID(i)
      end if
    end do

  end subroutine

  function hecmw_api_mesh_get_elem_type(mesh) bind(C,name='hecmw_api_mesh_get_elem_type')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_elem_type
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_elem_type = c_loc(hecMESH%elem_type(1))
  end function

  function hecmw_api_mesh_get_elem_node_item(mesh) bind(C,name='hecmw_api_mesh_get_elem_node_item')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_elem_node_item
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_elem_node_item = c_loc(hecMESH%elem_node_item(1))
  end function

  function hecmw_api_mesh_get_section_ID(mesh) bind(C,name='hecmw_api_mesh_get_section_id')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_section_ID
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_section_ID = c_loc(hecMESH%section_ID(1))
  end function

  function hecmw_api_mesh_n_elem(mesh) bind(C,name='hecmw_api_mesh_n_elem')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_elem
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_elem = hecMESH%n_elem
  end function

  function hecmw_api_mesh_n_elem_node_item(mesh) bind(C,name='hecmw_api_mesh_n_elem_node_item')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_elem_node_item
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_elem_node_item = hecMESH%elem_node_index(hecMESH%n_elem)
  end function

end module
