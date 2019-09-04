!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

module hecmw_couple_copy_f2c_f

  use hecmw_util
  use hecmw_couple_struct_f

  implicit none
  private
  public :: hecmw_couple_copy_f2c

contains

subroutine hecmw_couple_copy_f2c(couple_value, ierr)

  type(hecmw_couple_value), intent(inout) :: couple_value
  integer(kind=kint),       intent(inout) :: ierr
  character(len=100)                      :: sname, vname

  sname = "hecmw_couple_value"

  vname = "n"
  call hecmw_cpl_copy_f2c_set_if(sname, vname, couple_value%n, ierr)
  if(ierr /= 0) return

  vname = "item_type"
  call hecmw_cpl_copy_f2c_set_if(sname, vname, couple_value%item_type, ierr)
  if(ierr /= 0) return

  vname = "n_dof"
  call hecmw_cpl_copy_f2c_set_if(sname, vname, couple_value%n_dof, ierr)
  if(ierr /= 0) return

  if(couple_value%n > 0) then
    if(associated(couple_value%item)) then
      vname = "item"
      call hecmw_cpl_copy_f2c_set_if(sname, vname, couple_value%item, ierr)
      if(ierr /= 0) return
    endif
  endif

  if(couple_value%n > 0 .AND. couple_value%n_dof > 0) then
    if(associated(couple_value%value)) then
      vname = "value"
      call hecmw_cpl_copy_f2c_set_if(sname, vname, couple_value%value, ierr)
      if(ierr /= 0) return
    endif
  endif

end subroutine hecmw_couple_copy_f2c

end module hecmw_couple_copy_f2c_f
