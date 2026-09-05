!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_api_local_mesh
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
    allocate( hecMESH%node_group%grp_index(0:0) )
    hecMESH%node_group%grp_index(0) = 0
    allocate( hecMESH%surf_group%grp_index(0:0) )
    hecMESH%surf_group%grp_index(0) = 0
    allocate( hecMESH%elem_group%grp_index(0:0) )
    hecMESH%elem_group%grp_index(0) = 0
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

  !> @brief 節点自由度の設定
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] ndof 節点自由度
  subroutine hecmw_api_mesh_set_n_dof(mesh,ndof) bind(C,name='hecmw_api_mesh_set_n_dof')
    implicit none
    type(c_ptr), value :: mesh
    integer(c_int), value, intent(in) :: ndof
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecMESH%n_dof = ndof
  end subroutine

  !> @brief 節点自由度の取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 節点自由度
  function hecmw_api_mesh_get_n_dof(mesh) bind(C,name='hecmw_api_mesh_get_n_dof')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) :: hecmw_api_mesh_get_n_dof
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_n_dof = hecMESH%n_dof
  end function

  !> @brief 要素の設定
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] nelem 要素数
  !> @param[in] elemtype 要素型の配列
  !> @param[in] element 節点テーブルの配列 (1始まり)
  !> @param[in] sectionID セクション番号の配列 (1始まり)
  subroutine hecmw_api_mesh_set_element(mesh,nelem,elemtype,element,sectionID) bind(C,name='hecmw_api_mesh_set_element')
    use hecmw_etype, only : hecmw_get_max_node
    implicit none
    type(c_ptr), value                :: mesh
    integer(c_int), value, intent(in) :: nelem
    integer(c_int), intent(in)        :: element(*)
    integer(c_int), intent(in)        :: elemtype(nelem)
    integer(c_int), intent(in)        :: sectionID(nelem)

    type(hecmwST_local_mesh), pointer :: hecMESH
    integer :: i, j, ncon, n, ii, jj, off
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

    ! 要素の型ごとの個数を調べる
    etype_counter(:) = 0
    do i=1, nelem
      do j=1, 35
        if (elemtype(i) == etypes(j)) then
          etype_counter(j) = etype_counter(j) + 1
        end if
      end do
    end do
    ! 存在する型は何種類あるか
    n = 0
    do j=1, 35
      if (etype_counter(j) > 0) then
        n = n + 1
      end if
    end do

    hecMESH%n_elem_type = n
    allocate(hecMESH%elem_type_index(0:n))
    allocate(hecMESH%elem_type_item(n))

    ! 型ごとに要素がいくつかるかを記録
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

    ! elem_node_item の個数を調べる
    ncon = 0
    do i=1, nelem
      n = HECMW_get_max_node(elemtype(i))
      ncon = ncon + n
    end do
    allocate(hecMESH%elem_node_item(ncon))

    ! 要素の型ごとに並び替えて代入
    ncon = 0 ! elem_node_item　のオフセット
    hecMESH%elem_node_index(0) = ncon
    ii = 0
    do j=1, 35
      if (etype_counter(j) == 0) continue
      n = HECMW_get_max_node(etypes(j))
      off = 0 ! 入力配列のオフセット
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

    ! n_sect を調べる
    hecMESH%section%n_sect = 0
    do i=1, nelem
      if (hecMESH%section%n_sect < hecMESH%section_ID(i)) then
        hecMESH%section%n_sect = hecMESH%section_ID(i)
      end if
    end do

  end subroutine

  !> @brief 要素型の取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 要素型の配列のポインタ
  function hecmw_api_mesh_get_elem_type(mesh) bind(C,name='hecmw_api_mesh_get_elem_type')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_elem_type
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_elem_type = c_loc(hecMESH%elem_type(1))
  end function

  !> @brief 節点テーブルの取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 節点テーブル配列のポインタ
  function hecmw_api_mesh_get_elem_node_item(mesh) bind(C,name='hecmw_api_mesh_get_elem_node_item')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_elem_node_item
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_elem_node_item = c_loc(hecMESH%elem_node_item(1))
  end function

  !> @brief セクション番号配列の取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return セクション番号配列のポインタ
  function hecmw_api_mesh_get_section_ID(mesh) bind(C,name='hecmw_api_mesh_get_section_id')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(c_ptr) :: hecmw_api_mesh_get_section_ID
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_get_section_ID = c_loc(hecMESH%section_ID(1))
  end function

  !> @brief 要素数の取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 要素数
  function hecmw_api_mesh_n_elem(mesh) bind(C,name='hecmw_api_mesh_n_elem')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_elem
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_elem = hecMESH%n_elem
  end function

  !> @brief 要素節点配列の大きさの取得
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return 要素節点配列の大きさ
  function hecmw_api_mesh_n_elem_node_item(mesh) bind(C,name='hecmw_api_mesh_n_elem_node_item')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_elem_node_item
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_elem_node_item = hecMESH%elem_node_index(hecMESH%n_elem)
  end function

  !> @brief 節点グループのグループ数
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return グループの個数
  function hecmw_api_mesh_n_ngrp(mesh) bind(C,name='hecmw_api_mesh_n_ngrp')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_ngrp
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_ngrp = hecMESH%node_group%n_grp
  end function

  !> @brief 節点グループを追加
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] grp_name 節点グループの名前
  !> @param[in] count 節点の個数
  !> @param[in] list 節点番号の配列
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

  !> @brief 節点グループの名前
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] i 節点グループのインデックス（追加順）
  !> @param[out] buf 節点グループの名前
  !> @param[in] buflen 確保済みの buf の大きさ
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

  !> @brief 節点グループの節点番号の配列
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] i 節点グループのインデックス（追加順）
  !> @param[out] array 節点番号配列
  !> @param[out] count 配列の大きさ
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

  !> @brief 面グループのグループ数
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return グループの個数
  function hecmw_api_mesh_n_sgrp(mesh) bind(C,name='hecmw_api_mesh_n_sgrp')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_sgrp
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_sgrp = hecMESH%surf_group%n_grp
  end function

  !> @brief 面グループを追加
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] grp_name 面グループの名前
  !> @param[in] count 面の個数
  !> @param[in] list 面番号の配列
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

  !> @brief 面グループの名前
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] i 面グループのインデックス（追加順）
  !> @param[out] buf 節点グループの名前
  !> @param[in] buflen 確保済みの buf の大きさ
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

  !> @brief 面グループの要素番号、面番号の配列
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] i 面グループのインデックス（追加順）
  !> @param[out] array 要素番号、面番号配列
  !> @param[out] count 配列の大きさ
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

  !> @brief 要素グループのグループ数
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @return グループの個数
  function hecmw_api_mesh_n_egrp(mesh) bind(C,name='hecmw_api_mesh_n_egrp')
    implicit none
    type(c_ptr), value :: mesh
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(c_int) hecmw_api_mesh_n_egrp
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    hecmw_api_mesh_n_egrp = hecMESH%elem_group%n_grp
  end function

  !> @brief 要素グループを追加
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] grp_name 要素グループの名前
  !> @param[in] count 要素の個数
  !> @param[in] list 要素番号の配列
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

  !> @brief 要素グループの名前
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] i 要素グループのインデックス（追加順）
  !> @param[out] buf 節点グループの名前
  !> @param[in] buflen 確保済みの buf の大きさ
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
  
  !> @brief 要素グループの要素番号の配列
  !> @param[in] mesh メッシュ構造体のハンドラ
  !> @param[in] i 要素グループのインデックス（追加順）
  !> @param[out] array 要素番号配列
  !> @param[out] count 配列の大きさ
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
