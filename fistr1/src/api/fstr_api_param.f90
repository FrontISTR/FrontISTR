!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module fstr_api_param
  use iso_c_binding
  use m_fstr, only : fstr_param
  implicit none

contains

  !> @brief Generating parameter handler
  !! @return Parameter handler of type(c_ptr)
  function fstr_api_param_new() bind(C,name='fstr_api_param_new')
    implicit none
    type(c_ptr) :: fstr_api_param_new
    type(fstr_param), target, save :: fstrPARAM
    call fstr_nullify_fstr_param(fstrPARAM)
    fstr_api_param_new = c_loc(fstrPARAM)
  end function

  !> @brief Destroying the parameter handler
  !! @param[in] param : Parameter handler
  subroutine fstr_api_param_delete(param) bind(C,name='fstr_api_param_delete')
    use m_fstr, only : fstr_param_finalize
    implicit none
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    integer :: i
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    call fstr_param_finalize(fstrPARAM)
    deallocate(fstrPARAM)
  end subroutine

  !> @brief Initialize of the parameter sturcture
  !! @param[in] param : Parameter handler
  !! @param[in] mesh : Mesh handler
  subroutine fstr_api_param_init(param,mesh) bind(C,name='fstr_api_param_init')
    use m_fstr, only : fstr_param_init
    use hecmw, only : hecmwST_local_mesh
    use fstr_setup_util, only : reallocate_real, reallocate_integer
    implicit none
    type(c_ptr), value :: param
    type(c_ptr), value :: mesh
    type(fstr_param), pointer :: fstrPARAM
    type(hecmwST_local_mesh), pointer :: hecMESH
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    call c_f_pointer(cptr=mesh, fptr=hecMESH)
    call fstr_param_init(fstrPARAM,hecMESH)

    ! to run it without 'setup'
    call reallocate_real( fstrPARAM%dtime, 1 )
    call reallocate_real( fstrPARAM%etime, 1 )
    call reallocate_real( fstrPARAM%dtmin, 1 )
    call reallocate_real( fstrPARAM%delmax, 1 )
    call reallocate_integer( fstrPARAM%itmax, 1 )
    call reallocate_real( fstrPARAM%eps, 1 )
    fstrPARAM%dtime = 0.0d0
    fstrPARAM%etime = 0.0d0
    fstrPARAM%dtmin = 0.0d0
    fstrPARAM%delmax = 0.0d0
    fstrPARAM%itmax = 20
    fstrPARAM%eps = 1.0e-6

  end subroutine

  !> @brief Retrieve the solution type
  !! @param[in] param : Parameter handler
  !! @return solution type
  function fstr_api_param_solution_type(param) bind(C,name='fstr_api_param_solution_type')
    implicit none
    integer(c_int) :: fstr_api_param_solution_type
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_solution_type = fstrPARAM%solution_type
  end function

  !> @brief Set the solution type
  !! @param[in] param : Parameter handler
  !! @param[in] stype : solution type
  subroutine fstr_api_param_set_solution_type(param,stype) bind(C,name='fstr_api_param_set_solution_type')
    implicit none
    type(c_ptr), value :: param
    integer(c_int), value :: stype
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%solution_type = stype
  end subroutine

  !> @brief Retrieve the solver method
  !! @param[in] param : Parameter handler
  !! @return Solver method
  function fstr_api_param_solver_method(param) bind(C,name='fstr_api_param_solver_method')
    implicit none
    integer(c_int) :: fstr_api_param_solver_method
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_solver_method = fstrPARAM%solver_method
  end function

  !> @brief Set the solver method
  !! @param[in] param : Parameter handler
  !! @param[in] smethod : Solver method
  subroutine fstr_api_param_set_solver_method(param,smethod) bind(C,name='fstr_api_param_set_solver_method')
    implicit none
    type(c_ptr), value :: param
    integer(c_int) :: smethod
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%solver_method = smethod
  end subroutine

  !> @brief Retrieve the flag of nonliniear analysis
  !! @param[in] param : Parameter handler
  !! @return flag of nonlinear analysis
  function fstr_api_param_nlgeom(param) bind(C,name='fstr_api_param_nlgeom')
    implicit none
    logical(c_bool) :: fstr_api_param_nlgeom
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_nlgeom = fstrPARAM%nlgeom
  end function

  !> @brief Set the flag of nonlinear analysis
  !! @param[in] param : Parameter handler
  !! @param[in] nlgeom flag of nonlinear analysis
  subroutine fstr_api_param_set_nlgeom(param,nlgeom) bind(C,name='fstr_api_param_set_nlgeom')
    implicit none
    type(c_ptr), value :: param
    logical(c_bool), value :: nlgeom
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%nlgeom = nlgeom
  end subroutine

  !> @brief Retrieve the method of nonlinear solver
  !! @param[in] param : Parameter handler
  !! @return method of nonlinear solver
  function fstr_api_param_nlsolver_method(param) bind(C,name='fstr_api_param_nlsolver_method')
    implicit none
    integer(c_int) :: fstr_api_param_nlsolver_method
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_nlsolver_method = fstrPARAM%nlsolver_method
  end function

  !> @brief Set the method of nonlinear solver
  !! @param[in] param : Parameter handler
  !! @param[in] nlmethod : method of nonlinear solver
  subroutine fstr_api_param_set_nlsolver_method(param,nlmethod) bind(C,name='fstr_api_param_set_nlsolver_method')
    implicit none
    type(c_ptr), value :: param
    integer(c_int), value :: nlmethod
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%nlsolver_method = nlmethod
  end subroutine

  !> @brief Retrieve the flag of result output
  !! @param[in] param : Parameter handler
  !! @return flag of result output
  function fstr_api_param_fg_result(param) bind(C,name='fstr_api_param_fg_result')
    implicit none
    integer(c_int) :: fstr_api_param_fg_result
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_fg_result = fstrPARAM%fg_result
  end function

  !> @brief Set the flag of result output
  !! @param[in] param : Paramter handler
  !! @param[in] fg_result : flag of result output
  subroutine fstr_api_param_set_fg_result(param,fg_result) bind(C,name='fstr_api_param_set_fg_result')
    implicit none
    type(c_ptr), value :: param
    integer(c_int), value :: fg_result
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%fg_result = fg_result
  end subroutine

  !> @brief Retrieve the flag of visual output
  !! @param[in] param : Parameter handler
  !! @return flag of visual output
  function fstr_api_param_fg_visual(param) bind(C,name='fstr_api_param_fg_visual')
    implicit none
    integer(c_int) :: fstr_api_param_fg_visual
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_fg_visual = fstrPARAM%fg_visual
  end function

  !> @brief Set the flag of visual output
  !! @param[in] param : Parameter handler
  !! @return[in] fg_visual : flag of visual output
  subroutine fstr_api_param_set_fg_visual(param,fg_visual) bind(C,name='fstr_api_param_set_fg_visual')
    implicit none
    type(c_ptr), value :: param
    integer(c_int), value :: fg_visual
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%fg_visual = fg_visual
  end subroutine

  !> @brief Retrieve the algorithm of contact analysis
  !! @param[in] param : Parameter handler
  !! @return algorithm of contact analysis
  function fstr_api_param_contact_algo(param) bind(C,name='fstr_api_param_contact_algo')
    implicit none
    integer(c_int) :: fstr_api_param_contact_algo
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_contact_algo = fstrPARAM%contact_algo
  end function

  !> @brief Set the algorithm of contact analysis
  !! @param[in] param : Parameter handler
  !! @param[in] calgo : algorithm of contact analysis
  subroutine fstr_api_param_set_contact_algo(param,calgo) bind(C,name='fstr_api_param_set_contact_algo')
    implicit none
    type(c_ptr), value :: param
    integer(c_int), value :: calgo
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstrPARAM%contact_algo = calgo
  end subroutine

end module