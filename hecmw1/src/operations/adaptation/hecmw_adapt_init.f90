!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Adaptive Mesh Refinement

subroutine hecmw_adapt_init (hecMESH)

  use  hecmw_util
  type      (hecmwST_local_mesh) :: hecMESH

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#allocate (n_node, nn_int, n_elem, ne_int)'
  call hecmw_adapt_allocate   (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#ACTIVE'
  call hecmw_adapt_ACTIVE (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#E-COMM.TAB.'
  call hecmw_adapt_EDGE_COMM_TABLE (hecMESH)
  if (hecMESH%my_rank.eq.0) write (*,'(  a)') '#C-COMM.TAB.'
  call hecmw_adapt_CELL_COMM_TABLE (hecMESH)

end subroutine hecmw_adapt_init
