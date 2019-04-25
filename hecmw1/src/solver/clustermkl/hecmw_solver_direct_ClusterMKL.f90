!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface for PARADISO
module hecmw_solver_direct_ClusterMKL
  use hecmw_util
  use m_hecmw_ClusterMKL_wrapper

  private
  public :: hecmw_solve_direct_ClusterMKL

contains

  subroutine hecmw_solve_direct_ClusterMKL(hecMESH,hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT

    call hecmw_clustermkl_wrapper(hecMESH, hecMAT)

  end subroutine hecmw_solve_direct_ClusterMKL

end module hecmw_solver_direct_ClusterMKL
