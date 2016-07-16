!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.8                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Naoki Morita (Univ. of Tokyo)                  !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
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
