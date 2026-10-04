!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module fstr_api_param
  use iso_c_binding
  use m_fstr, only : fstr_param
  implicit none

contains

  !> @brief パラメータハンドラの生成
  !! @return パラメータ構造体のハンドラ type(c_ptr)
  function fstr_api_param_new() bind(C,name='fstr_api_param_new')
    implicit none
    type(c_ptr) :: fstr_api_param_new
    type(fstr_param), target, save :: fstrPARAM
    call fstr_nullify_fstr_param(fstrPARAM)
    fstr_api_param_new = c_loc(fstrPARAM)
  end function

  !> @brief パラメータハンドラの破棄
  !! @param[in] param パラメータ構造体のハンドラ
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

  !> @brief パラメータ構造体の初期化
  !! @param[in] param パラメータ構造体のハンドラ
  !! @param[in] mesh メッシュ構造体のハンドラ
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

    ! setup なしでも動作させるため
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

  !> @brief 解析の種別の取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return 解析の種別
  function fstr_api_param_solution_type(param) bind(C,name='fstr_api_param_solution_type')
    implicit none
    integer(c_int) :: fstr_api_param_solution_type
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_solution_type = fstrPARAM%solution_type
  end function

  !> @brief ソルバーの解法の取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return ソルバーの解法
  function fstr_api_param_solver_method(param) bind(C,name='fstr_api_param_solver_method')
    implicit none
    integer(c_int) :: fstr_api_param_solver_method
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_solver_method = fstrPARAM%solver_method
  end function

  !> @brief 非線形を考慮するかのフラグの取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return 非線形を考慮するか
  function fstr_api_param_nlgeom(param) bind(C,name='fstr_api_param_nlgeom')
    implicit none
    logical(c_bool) :: fstr_api_param_nlgeom
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_nlgeom = fstrPARAM%nlgeom
  end function

  !> @brief 非線形ソルバーの解法の取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return 非線形ソルバーの解法
  function fstr_api_param_nlsolver_method(param) bind(C,name='fstr_api_param_nlsolver_method')
    implicit none
    integer(c_int) :: fstr_api_param_nlsolver_method
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_nlsolver_method = fstrPARAM%nlsolver_method
  end function

  !> @brief 結果を出力するかのフラグの取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return 結果を出力するかのフラグ
  function fstr_api_param_fg_result(param) bind(C,name='fstr_api_param_fg_result')
    implicit none
    integer(c_int) :: fstr_api_param_fg_result
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_fg_result = fstrPARAM%fg_result
  end function

  !> @brief 可視化出力するかのフラグの取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return 可視化出力するかのフラグ
  function fstr_api_param_fg_visual(param) bind(C,name='fstr_api_param_fg_visual')
    implicit none
    integer(c_int) :: fstr_api_param_fg_visual
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_fg_visual = fstrPARAM%fg_visual
  end function

  !> @brief 接触解析アルゴリズムの取得
  !! @param[in] param パラメータ構造体のハンドラ
  !! @return 接触解析アルゴリズム
  function fstr_api_param_contact_algo(param) bind(C,name='fstr_api_param_contact_algo')
    implicit none
    integer(c_int) :: fstr_api_param_contact_algo
    type(c_ptr), value :: param
    type(fstr_param), pointer :: fstrPARAM
    call c_f_pointer(cptr=param, fptr=fstrPARAM)
    fstr_api_param_contact_algo = fstrPARAM%contact_algo
  end function

end module