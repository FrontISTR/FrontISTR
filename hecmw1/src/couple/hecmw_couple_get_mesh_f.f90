!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

module hecmw_couple_get_mesh_f

  use hecmw_util
  use hecmw_io
  use hecmw_couple_struct_f
  use hecmw_couple_info_f

  implicit none
  private
  public :: hecmw_couple_get_mesh

contains

subroutine hecmw_couple_get_mesh(name_ID, unit_ID, mesh)

  character(len=HECMW_NAME_LEN), intent(in)  :: name_ID
  character(len=HECMW_NAME_LEN), intent(in)  :: unit_ID
  type(hecmwST_local_mesh),      intent(out) :: mesh

  call hecmw_get_mesh(name_ID, mesh)

  mesh%PETOT    = hecmw_intracomm_get_size_u(unit_ID)
  mesh%my_rank  = hecmw_intracomm_get_rank_u(unit_ID)
  mesh%MPI_COMM = hecmw_intracomm_get_comm_u(unit_ID)

end subroutine hecmw_couple_get_mesh

end module hecmw_couple_get_mesh_f
