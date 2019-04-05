!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Coupling Interface

module hecmw_couple_startup_f

  use hecmw_util
  use hecmw_couple_copy_c2f_f
  use hecmw_couple_define_f
  use hecmw_couple_struct_f

  implicit none
  private
  public :: hecmw_couple_startup
  public :: hecmw_couple_print_couple_value

contains

subroutine hecmw_couple_startup(boundary_id, couple_value)

  character(len=HECMW_NAME_LEN), intent(in)  :: boundary_id
  type(hecmw_couple_value),      intent(out) :: couple_value
  integer(kind=kint)                         :: ierr

  call hecmw_couple_startup_init_if(boundary_id, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_copy_c2f(couple_value, ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

  call hecmw_couple_startup_final_if(ierr)
  if(ierr /= 0) call hecmw_abort(hecmw_comm_get_comm())

end subroutine hecmw_couple_startup


subroutine hecmw_couple_print_couple_value(couple_value)

  type(hecmw_couple_value), intent(in) :: couple_value
  integer(kind=kint)                   :: ierr, i, j

  write(unit=*,fmt="(a)") "*** Value of coupling area"

  write(unit=*,fmt="(a,i)") "number of item: ", couple_value%n

  if(couple_value%item_type == HECMW_COUPLE_NODE_GROUP) then
    write(unit=*,fmt="(a)") "item type: NODE GROUP"
  else if(couple_value%item_type == HECMW_COUPLE_ELEMENT_GROUP) then
    write(unit=*,fmt="(a)") "item type: ELEMENT GROUP"
  else if(couple_value%item_type == HECMW_COUPLE_SURFACE_GROUP) then
    write(unit=*,fmt="(a)") "item type: SURFACE GROUP"
  else
    write(unit=*,fmt="(a)") "item type: UNKNOWN"
  endif

  write(unit=*,fmt="(a,i)") "number of DOF: ", couple_value%n_dof

  write(unit=*,fmt="(a)") "ID & value:"
  do i= 1, couple_value%n
    if(couple_value%item_type == HECMW_COUPLE_SURFACE_GROUP) then
      write(unit=*,fmt="(2i)",advance="NO") couple_value%item(2*i-1), couple_value%item(2*i)
    else
      write(unit=*,fmt="(i)",advance="NO") couple_value%item(i)
    endif

    do j= 1, couple_value%n_dof
      write(unit=*,fmt="(e15.7)",advance="NO") couple_value%value(couple_value%n_dof*(i-1)+j)
    enddo
    write(unit=*,fmt="(/)",advance="NO")
  enddo

end subroutine hecmw_couple_print_couple_value

end module hecmw_couple_startup_f
