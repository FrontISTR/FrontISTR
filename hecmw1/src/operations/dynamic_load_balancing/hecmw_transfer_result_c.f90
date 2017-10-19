!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Dynamic Load Balancing

subroutine hecmw_dist_result_copy_f2c(adapRES)
  use   hecmw_util
  use   hecmw_result
  use   hecmw_io
  use   hecmw_dist_copy_f2c_f
  type (hecmwST_result_data):: adapRES



  if(adapRES%nn_component .gt.0) then
    call  hecmw_set_result_node(adapRES%nn_component,adapRES%nn_dof,adapRES%node_label, &
      adapRES%node_val_item)
    deallocate (adapRES%node_val_item)
  endif
  if(adapRES%ne_component .gt.0) then
    call  hecmw_set_result_elem(adapRES%ne_component,adapRES%ne_dof, adapRES%elem_label, &
      adapRES%elem_val_item)
    deallocate (adapRES%elem_val_item)
  endif

end subroutine  hecmw_dist_result_copy_f2c

