!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_api_result_data
  use iso_c_binding
  use hecmw_result, only : hecmwST_result_data
  implicit none

contains

  function hecmw_api_result_new() bind(C,name='hecmw_api_result_new')
    implicit none
    type(c_ptr) :: hecmw_api_result_new
    type(hecmwST_result_data), target, save :: hecRESULT
    hecmw_api_result_new = c_loc(hecRESULT)
  end function

  subroutine hecmw_api_result_delete(result) bind(C,name='hecmw_api_result_delete')
    implicit none
    type(c_ptr), value :: result
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    call hecmw_result_free(hecRESULT)
  end subroutine

  ! result に格納されている大域データの個数
  function hecmw_api_result_ng_component(result) bind(C,name='hecmw_api_result_ng_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_ng_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_ng_component = hecRESULT%ng_component
  end function

  ! i 番目の大域データのラベルと値を取得
  subroutine hecmw_api_result_global_val(result,i,dof,label,label_len,value) bind(C,name='hecmw_api_result_global_val')
    use hecmw_api_common, only : f_c_str_copy
    use hecmw_result
    implicit none
    type(c_ptr), value :: result
    integer(c_int), value :: i
    integer(c_int), intent(out) :: dof
    character(kind=c_char), intent(out) :: label(*)
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

  ! result に格納されている節点データの個数
  function hecmw_api_result_nn_component(result) bind(C,name='hecmw_api_result_nn_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_nn_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_nn_component = hecRESULT%nn_component
  end function

  !
  ! i 番目の節点データを取得 
  ! 節点ごとに dim 個の値が並んでいる index から始まって dof 個
  ! nv(1:dim,1:n_node) に reshape して nv(index:index+dof,:) で取り出す
  ! 
  subroutine hecmw_api_result_node_val(result,i,dim,index,dof,label,label_len,value) bind(C,name='hecmw_api_result_node_val')
    use hecmw_api_common, only : f_c_str_copy
    use hecmw_result
    implicit none
    type(c_ptr), value :: result
    integer(c_int), value :: i
    integer(c_int), intent(out) :: dim
    integer(c_int), intent(out) :: index
    integer(c_int), intent(out) :: dof
    character(kind=c_char), intent(out) :: label(*)
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

  ! result に格納されている要素データの個数
  function hecmw_api_result_ne_component(result) bind(C,name='hecmw_api_result_ne_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_ne_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_ne_component = hecRESULT%ne_component
  end function

  !
  ! i 番目の要素データを取得 
  ! 要素ごとに dim 個の値が並んでいる index から始まって dof 個
  ! ev(1:dim,1:n_elem) に reshape して ev(index:index+dof,:) で取り出す
  ! 
  subroutine hecmw_api_result_elem_val(result,i,dim,index,dof,label,label_len,value) bind(C,name='hecmw_api_result_elem_val')
    use hecmw_api_common, only : f_c_str_copy
    use hecmw_result
    implicit none
    type(c_ptr), value :: result
    integer(c_int), value :: i
    integer(c_int), intent(out) :: dim
    integer(c_int), intent(out) :: index
    integer(c_int), intent(out) :: dof
    character(kind=c_char), intent(out) :: label(*)
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