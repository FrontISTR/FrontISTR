module hecmw_local_mesh_api
  use iso_c_binding
  use hecmw_util, only : hecmwST_local_mesh, hecmw_nullify_mesh
  use hecmw_dist_free_f, only : hecmw_dist_free
  implicit none

contains

  !> @brief メッシュハンドラの生成
  !> @return メッシュ構造体のハンドラ type(c_ptr)
  function hecmw_api_mesh_new() bind(C,name='hecmw_api_mesh_new')
    implicit none
    type(c_ptr) :: hecmw_api_mesh_new
    type(hecmwST_local_mesh), target, save :: hecMESH
    call hecmw_nullify_mesh(hecMESH)
    hecmw_api_mesh_new = c_loc(hecMESH)
  end function

  !> @brief メッシュハンドラの破棄
  !> @param[in] mesh メッシュ構造体のハンドラ
  subroutine hecmw_api_mesh_delete(mesh) bind(C,name='hecmw_api_mesh_delete')
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

end module
