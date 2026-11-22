!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module fstr_api_solid
  use iso_c_binding
  use m_fstr
  implicit none

contains

  type(c_ptr) function fstr_api_solid_new() &
    bind(C,name='fstr_api_solid_new')
    implicit none
    type(fstr_solid), pointer :: fstrSOLID => null()
    allocate(fstrSOLID)
    call fstr_nullify_fstr_solid(fstrSOLID)
    fstr_api_solid_new = c_loc(fstrSOLID)
  end function

  subroutine fstr_api_solid_delete(solid) &
    bind(C,name='fstr_api_solid_delete')
    use m_fstr_setup
    implicit none
    type(c_ptr), value :: solid
    type(fstr_solid), pointer :: fstrSOLID
    call c_f_pointer(cptr=solid, fptr=fstrSOLID)
    call fstr_solid_finalize( fstrSOLID )
  end subroutine

  subroutine fstr_api_solid_init(solid,mesh,param) &
    bind(C,name='fstr_api_solid_init')
    use m_fstr_setup
    implicit none

    type(c_ptr), value :: solid
    type(c_ptr), value :: mesh
    type(c_ptr), value :: param
    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_solid), pointer :: fstrSOLID
    type(fstr_param), pointer :: fstrPARAM

    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call c_f_pointer(cptr=solid, fptr=fstrSOLID)
    call c_f_pointer(cptr=param, fptr=fstrPARAM)

    call fstr_solid_init(fstrSOLID)
    call fstr_solid_init_section(hecMESH,fstrSOLID,fstrPARAM)
    call fstr_solid_init_output(fstrSOLID)
  end subroutine

  integer(c_int) function fstr_api_solid_nmaterial_tot(solid) &
    bind(C,name='fstr_api_solid_nmaterial_tot')
    implicit none
    type(c_ptr), value :: solid
    type(fstr_solid), pointer :: fstrSOLID
    call c_f_pointer(cptr=solid, fptr=fstrSOLID)
    fstr_api_solid_nmaterial_tot = 0
    if( associated(fstrSOLID%materials) ) then
      fstr_api_solid_nmaterial_tot = size(fstrSOLID%materials)
    endif
  end function

  subroutine fstr_api_solid_init_material(solid,n) &
    bind(C,name='fstr_api_solid_init_material')
    implicit none
    type(c_ptr), value :: solid
    integer(c_int), value, intent(in) :: n
    type(fstr_solid), pointer :: fstrSOLID
    integer :: i

    call c_f_pointer(cptr=solid, fptr=fstrSOLID)
    allocate( fstrSOLID%materials( n ) )
    do i = 1, n
      call initMaterial(fstrSOLID%materials(i))
    enddo
  end subroutine

end module
