!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Coupling Interface                                 !
!                                                                      !
!            Written by Shin'ichi Ezure (RIST)                         !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!


module hecmw_couple_finalize_f

  use hecmw_util

  implicit none
  private
  public :: hecmw_couple_finalize

contains

subroutine hecmw_couple_finalize(boundary_id)

  character(len=HECMW_NAME_LEN), intent(in) :: boundary_id
  integer(kind=kint)                        :: ierr

  call hecmw_couple_finalize_if(boundary_id, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

end subroutine hecmw_couple_finalize

end module hecmw_couple_finalize_f
