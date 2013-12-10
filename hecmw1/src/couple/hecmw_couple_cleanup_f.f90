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


module hecmw_couple_cleanup_f

  use hecmw_util
  use hecmw_couple_define_f
  use hecmw_couple_struct_f

  implicit none
  private
  public :: hecmw_couple_cleanup

contains

subroutine hecmw_couple_cleanup(couple_value)

  type(hecmw_couple_value), intent(inout) :: couple_value
  integer(kind=kint)                      :: ista

  if(associated(couple_value%item)) then
    deallocate(couple_value%item, stat=ista)
    if(ista > 0) call hecmw_abort(hecmw_comm_get_comm())
  endif

  if(associated(couple_value%value)) then
    deallocate(couple_value%value, stat=ista)
    if(ista > 0) call hecmw_abort(hecmw_comm_get_comm())
  endif

  couple_value%n         = 0
  couple_value%item_type = HECMW_COUPLE_GROUP_UNDEF
  couple_value%n_dof     = 0
  nullify(couple_value%item)
  nullify(couple_value%value)

end subroutine hecmw_couple_cleanup

end module hecmw_couple_cleanup_f
