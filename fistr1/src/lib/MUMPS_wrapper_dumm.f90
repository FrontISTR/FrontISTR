!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.4                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by K. Goto (VINAS)                                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module provides wrapper for parallel sparse direct solver MUMPS
module m_MUMPS_wrapper
  use m_fstr
  use m_sparse_matrix

  private
  public :: mumps_wrapper

contains

  subroutine mumps_wrapper(spMAT, job, istat)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    integer(kind=kint), intent(out) :: istat
    stop "MUMPS not available"
  end subroutine mumps_wrapper

end module m_MUMPS_wrapper
