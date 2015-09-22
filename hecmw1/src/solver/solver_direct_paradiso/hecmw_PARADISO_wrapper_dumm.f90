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
!> This module provides wrapper for parallel sparse direct solver PARADISO
module m_hecmw_PARADISO_wrapper
  use hecmw_util
  use m_sparse_matrix

  private
  public :: hecmw_paradiso_wrapper

contains

  subroutine hecmw_mumps_wrapper(hecMESH, hecMAT, ii)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    integer(kind=kint) :: ii

    stop "PARADISO not available"
  end subroutine hecmw_paradiso_wrapper

end module m_hecmw_PARADISO_wrapper
