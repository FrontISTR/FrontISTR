!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Adaptive Mesh Refinement

subroutine hecmw_adapt_proc (hecMESH)

  use  hecmw_util
  type      (hecmwST_local_mesh) :: hecMESH

  hecMESH%n_adapt= hecMESH%n_adapt + 1

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#EXTEND EMB.'
  call hecmw_adapt_EXTEMB      (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(a)') '#GRID smoothing'
  call hecmw_adapt_GRID_SMOOTH (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#create NEW nodes'
  call hecmw_adapt_NEW_NODE (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(  a)') '#create NEW CELL'
  call hecmw_adapt_NEW_CELL (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(  a)') '#create CELL INFO.'
  call hecmw_adapt_GET_NEW_CELL_INFO (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#create NEW BC pointer'
  call hecmw_adapt_BC_POINTER (hecMESH)

  if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#create NEW comm. table'
  call hecmw_adapt_REPRO_COMM_TABLE (hecMESH)

end subroutine hecmw_adapt_proc
