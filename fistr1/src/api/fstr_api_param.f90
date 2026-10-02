!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module fstr_api_param
  use iso_c_binding
  use m_fstr, only : fstr_param
  implicit none

contains

  function fstr_api_param_new() bind(C,name='fstr_api_param_new')
    implicit none
    type(c_ptr) :: fstr_api_param_new
    type(fstr_param), target, save :: fstrPARAM
    call fstr_nullify_fstr_param(fstrPARAM)
    fstr_api_param_new = c_loc(fstrPARAM)
  end function

  subroutine fstr_api_param_delete(param) bind(C,name='fstr_api_param_delete')
    use m_timepoint, only : time_points
    implicit none
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    integer :: i
    call c_f_pointer(cptr=param, fptr=fstrPARAM)

    if( associated(fstrPARAM%dtime) ) deallocate(fstrPARAM%dtime)
    if( associated(fstrPARAM%etime) ) deallocate(fstrPARAM%etime)
    if( associated(fstrPARAM%dtmin) ) deallocate(fstrPARAM%dtmin)
    if( associated(fstrPARAM%delmax) ) deallocate(fstrPARAM%delmax)
    if( associated(fstrPARAM%itmax) ) deallocate(fstrPARAM%itmax)
    if( associated(fstrPARAM%eps) ) deallocate(fstrPARAM%eps)
    if( associated(fstrPARAM%global_local_ID) ) deallocate(fstrPARAM%global_local_ID)
    if( associated(fstrPARAM%contactparam) ) deallocate(fstrPARAM%contactparam)
    if( associated(fstrPARAM%contact_if) ) deallocate(fstrPARAM%contact_if)
    if( associated(fstrPARAM%ainc) ) deallocate(fstrPARAM%ainc)
    if( associated(fstrPARAM%timepoints) ) then
      do i=1, size(fstrPARAM%timepoints)
        if( associated(fstrPARAM%timepoints(i)%points) ) deallocate(fstrPARAM%timepoints(i)%points)
      end do
      deallocate(fstrPARAM%timepoints)
    end if
    if( associated(fstrPARAM%cnvparam) ) deallocate(fstrPARAM%cnvparam)
    deallocate(fstrPARAM)
  end subroutine

  subroutine fstr_api_param_init(param,mesh) bind(C,name='fstr_api_param_init')
    use m_fstr, only : fstr_param_init
    use hecmw, only : hecmwST_local_mesh
    implicit none
    type(c_ptr), value :: param
    type(c_ptr), value :: mesh
    type(fstr_param), pointer :: fstrPARAM
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call fstr_param_init(fstrPARAM,hecMESH)
  end subroutine

  function fstr_api_param_solutuin_type(param) bind(C,name='fstr_api_param_solution_type')
    implicit none
    integer(c_int) :: fstr_api_param_solutuin_type
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_solutuin_type = fstrPARAM%solution_type
  end function

end module