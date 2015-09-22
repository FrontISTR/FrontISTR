!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.7                                                !
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
module hecmw_solver_direct_PARADISO
  use hecmw_util
  use m_hecmw_PARADISO_wrapper

  private
  public :: hecmw_solve_direct_PARADISO

contains

  subroutine hecmw_solve_direct_PARADISO(hecMESH,hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT

    call hecmw_paradiso_wrapper(hecMESH, hecMAT)

  end subroutine hecmw_solve_direct_PARADISO

end module hecmw_solver_direct_PARADISO
