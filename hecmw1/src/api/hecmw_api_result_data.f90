!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_api_result_data
  use iso_c_binding
  use hecmw_result, only : hecmwST_result_data
  implicit none

contains

  !> @brief Generating result handler
  !! @return Result handler of type(c_ptr)
  function hecmw_api_result_new() bind(C,name='hecmw_api_result_new')
    use hecmw_result, only : hecmw_nullify_result_data
    implicit none
    type(c_ptr) :: hecmw_api_result_new
    type(hecmwST_result_data), pointer :: hecRESULT => null()
    allocate(hecRESULT)
    call hecmw_nullify_result_data(hecRESULT)
    hecRESULT%ng_component = 0
    hecRESULT%nn_component = 0
    hecRESULT%ne_component = 0
    hecmw_api_result_new = c_loc(hecRESULT)
  end function

  !> @brief Destroying the result handler
  !! @param[in] result : Result handler
  subroutine hecmw_api_result_delete(result) bind(C,name='hecmw_api_result_delete')
    use hecmw_result, only : hecmw_result_free
    implicit none
    type(c_ptr), value :: result
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    call hecmw_result_free(hecRESULT)
    deallocate(hecRESULT)
  end subroutine

  !> @brief Retrieving the number of global data items included in the results
  !! @param[in] result : Result handler
  !! @return Number of global data items
  function hecmw_api_result_ng_component(result) bind(C,name='hecmw_api_result_ng_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_ng_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_ng_component = hecRESULT%ng_component
  end function

  !> @brief Retrieve the labels and values ​​of the global data contained in the results.
  !! @param[in] result : Result handler
  !! @param[in] i : Index of the stored global data
  !! @param[out] dof : Degrees of freedom of global data
  !! @param[inout] label : Name of the global data (already allocated by the caller)
  !! @param[in] label_len : Size of the allocated label
  !! @param[out] value : Pointer to the beginning of global data
  !! @remark stored in value(1:dof)
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

  !> @brief Retrieving the number of node data entries in the results
  !! @param[in] result : Result handler
  !! @return Number of nodal data points
  function hecmw_api_result_nn_component(result) bind(C,name='hecmw_api_result_nn_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_nn_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_nn_component = hecRESULT%nn_component
  end function

  !> @brief Retrieve the labels and values ​​of the node data contained in the results.
  !! @param[in] result : Result handler
  !! @param[in] i : Index of the stored node data
  !! @param[out] dim : Degrees of freedom for all data stored at each node
  !! @param[out] index : An index in which data is stored for each node
  !! @param[out] dof : Degrees of freedom of nodal data
  !! @param[inout] label : Name of the node data (allocated by the caller)
  !! @param[in] label_len : Reserved label size
  !! @param[out] value : Pointer to the beginning of the node data
  !! @remark Stored in (index:index+dof, :) after reshaping value into (1:dim, 1:n_node)
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

  !> @brief Retrieving the number of element data items contained in the result
  !! @param[in] result : Result handler
  !! @return Number of element data items
  function hecmw_api_result_ne_component(result) bind(C,name='hecmw_api_result_ne_component')
    implicit none
    type(c_ptr), value :: result
    integer(c_int) :: hecmw_api_result_ne_component
    type(hecmwST_result_data), pointer :: hecRESULT
    call c_f_pointer(cptr=result, fptr=hecRESULT)
    hecmw_api_result_ne_component = hecRESULT%ne_component
  end function

  !> @brief Retrieve the labels and values ​​of the element data contained in the result.
  !! @param[in] result : Result handler
  !! @param[in] i : Index of the stored element data
  !! @param[out] dim : Degrees of freedom for all data stored for each element
  !! @param[out] index : An index in which data is stored for each element
  !! @param[out] dof : Degrees of freedom of element data
  !! @param[inout] label : Name of the element data (already allocated by the caller)
  !! @param[in] label_len : Reserved label size
  !! @param[out] value : Pointer to the beginning of the element data
  !! @remark Stored in (index:index+dof, :) after reshaping value to (1:dim, 1:n_elem)
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
