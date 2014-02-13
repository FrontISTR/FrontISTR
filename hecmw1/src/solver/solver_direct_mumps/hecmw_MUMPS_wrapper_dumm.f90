!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
!> This module provides wrapper for parallel sparse direct solver MUMPS
module m_hecmw_MUMPS_wrapper
  use hecmw_util
  use m_sparse_matrix

  private
  public :: hecmw_mumps_wrapper

contains

  subroutine hecmw_mumps_wrapper(spMAT, job, istat)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    integer(kind=kint), intent(out) :: istat
    stop "MUMPS not available"
  end subroutine hecmw_mumps_wrapper

end module m_hecmw_MUMPS_wrapper
