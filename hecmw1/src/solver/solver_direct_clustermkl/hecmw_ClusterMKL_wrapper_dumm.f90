!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides wrapper for parallel sparse direct solver PARADISO
module m_hecmw_ClusterMKL_wrapper
  use hecmw_util
  use m_sparse_matrix

  private
  public :: hecmw_clustermkl_wrapper

contains

  subroutine hecmw_clustermkl_wrapper(hecMESH, hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT

    stop "PARADISO not available"
  end subroutine hecmw_clustermkl_wrapper

end module m_hecmw_ClusterMKL_wrapper
