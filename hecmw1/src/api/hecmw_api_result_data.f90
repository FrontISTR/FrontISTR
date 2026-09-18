!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_api_result_data
  use iso_c_binding
  use hecmw_result, only : hecmwST_result_data
  implicit none

contains

  !> @brief 計算結果ハンドラの生成
  !! @return 計算結果構造体のハンドラ type(c_ptr)
  function hecmw_api_result_new() bind(C,name='hecmw_api_result_new')
    implicit none
    type(c_ptr) :: hecmw_api_result_new
    type(hecmwST_result_data), target, save :: hecRESULT
    hecmw_api_result_new = c_loc(hecRESULT)
  end function

  !> @brief 計算結果ハンドラの破棄
  !! @param[in] result 計算結果構造体のハンドラ
  subroutine hecmw_api_result_delete(result) bind(C,name='hecmw_api_result_delete')
    implicit none
    type(c_ptr), value :: result
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    call hecmw_result_free(hecRESULT)
  end subroutine

  !> @brief 計算結果に含まれる大域データの個数の取得
  !! @param[in] result 計算結果構造体のハンドラ
  !! @return 大域データの個数
  function hecmw_api_result_ng_component(result) bind(C,name='hecmw_api_result_ng_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_ng_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_ng_component = hecRESULT%ng_component
  end function

  !> @brief 計算結果に含まれる大域データのラベルと値を取得
  !! @param[in] result 計算結果構造体のハンドラ
  !! @param[in] i 格納されている大域データのインデックス
  !! @param[out] dof 大域データの自由度
  !! @param[inout] label 大域データの名前（呼び出し側で確保済み）
  !! @param[in] label_len 確保済みの label の大きさ
  !! @param[out] value 大域データの先頭ポインタ
  !! @remark value(1:dof) に格納されている
  subroutine hecmw_api_result_global_val(result,i,dof,label,label_len,value) bind(C,name='hecmw_api_result_global_val')
    use hecmw_api_common, only : f_c_str_copy
    use hecmw_result
    implicit none
    type(c_ptr), value :: result
    integer(c_int), value :: i
    integer(c_int), intent(out) :: dof
    character(kind=c_char), intent(inout) :: label(*)
    integer(c_int), value, intent(in) :: label_len
    type(c_ptr), intent(out) :: value

    type(hecmwST_result_data), pointer :: hecRESULT
    integer :: index, j

    call c_f_pointer(cptr=result, fptr=hecRESULT)
    call f_c_str_copy(hecRESULT%global_label(i), label, label_len)
    index = 1
    do j=1, i-1
      index = index + hecRESULT%ng_dof(j)
    end do
    
    dof = hecRESULT%ng_dof(i)
    value = c_loc(hecRESULT%global_val_item(index))
  end subroutine

  !> @brief 計算結果に含まれる節点データの個数の取得
  !! @param[in] result 計算結果構造体のハンドラ
  !! @return 節点データの個数
  function hecmw_api_result_nn_component(result) bind(C,name='hecmw_api_result_nn_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_nn_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_nn_component = hecRESULT%nn_component
  end function

  !> @brief 計算結果に含まれる節点データのラベルと値を取得
  !! @param[in] result 計算結果構造体のハンドラ
  !! @param[in] i 格納されている節点データのインデックス
  !! @param[out] dim 節点ごとに格納されている全データの自由度
  !! @param[out] index 節点ごとにデータが格納されているインデックス
  !! @param[out] dof 節点データの自由度
  !! @param[inout] label 節点データの名前（呼び出し側で確保済み）
  !! @param[in] label_len 確保済みの label の大きさ
  !! @param[out] value 節点データの先頭ポインタ
  !! @remark value を (1:dim,1:n_node) に reshape した (index:index+dof,:) に格納されている
  subroutine hecmw_api_result_node_val(result,i,dim,index,dof,label,label_len,value) bind(C,name='hecmw_api_result_node_val')
    use hecmw_api_common, only : f_c_str_copy
    use hecmw_result
    implicit none
    type(c_ptr), value :: result
    integer(c_int), value :: i
    integer(c_int), intent(out) :: dim
    integer(c_int), intent(out) :: index
    integer(c_int), intent(out) :: dof
    character(kind=c_char), intent(inout) :: label(*)
    integer(c_int), value, intent(in) :: label_len
    type(c_ptr), intent(out) :: value

    type(hecmwST_result_data), pointer :: hecRESULT
    integer :: j

    call c_f_pointer(cptr=result, fptr=hecRESULT)
    call f_c_str_copy(hecRESULT%node_label(i), label, label_len)

    dof = hecRESULT%nn_dof(i)
    index = 0
    dim = 0
    do j=1, i-1
      index = index + hecRESULT%nn_dof(j)
      dim = dim + hecRESULT%nn_dof(j)
    end do
    do j=i, hecRESULT%nn_component
      dim = dim + hecRESULT%nn_dof(j)
    end do
    value = c_loc(hecRESULT%node_val_item)
  end subroutine

  !> @brief 計算結果に含まれる要素データの個数の取得
  !! @param[in] result 計算結果構造体のハンドラ
  !! @return 要素データの個数
  function hecmw_api_result_ne_component(result) bind(C,name='hecmw_api_result_ne_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_ne_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_ne_component = hecRESULT%ne_component
  end function

  !> @brief 計算結果に含まれる要素データのラベルと値を取得
  !! @param[in] result 計算結果構造体のハンドラ
  !! @param[in] i 格納されている要素データのインデックス
  !! @param[out] dim 要素ごとに格納されている全データの自由度
  !! @param[out] index 要素ごとにデータが格納されているインデックス
  !! @param[out] dof 要素データの自由度
  !! @param[inout] label 要素データの名前（呼び出し側で確保済み）
  !! @param[in] label_len 確保済みの label の大きさ
  !! @param[out] value 要素データの先頭ポインタ
  !! @remark value を (1:dim,1:n_elem) に reshape した (index:index+dof,:) に格納されている
  subroutine hecmw_api_result_elem_val(result,i,dim,index,dof,label,label_len,value) bind(C,name='hecmw_api_result_elem_val')
    use hecmw_api_common, only : f_c_str_copy
    use hecmw_result
    implicit none
    type(c_ptr), value :: result
    integer(c_int), value :: i
    integer(c_int), intent(out) :: dim
    integer(c_int), intent(out) :: index
    integer(c_int), intent(out) :: dof
    character(kind=c_char), intent(inout) :: label(*)
    integer(c_int), value, intent(in) :: label_len
    type(c_ptr), intent(out) :: value

    integer :: j
    type(hecmwST_result_data), pointer :: hecRESULT

    call c_f_pointer(cptr=result, fptr=hecRESULT)
    call f_c_str_copy(hecRESULT%elem_label(i), label, label_len)

    dof = hecRESULT%ne_dof(i)
    index = 0
    dim = 0
    do j=1, i-1
      index = index + hecRESULT%ne_dof(j)
      dim = dim + hecRESULT%ne_dof(j)
    end do
    do j=i, hecRESULT%ne_component
      dim = dim + hecRESULT%ne_dof(j)
    end do
    value = c_loc(hecRESULT%elem_val_item)
  end subroutine

end module