!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

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
