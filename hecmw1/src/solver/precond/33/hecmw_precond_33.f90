!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_33
  use hecmw_util
  use hecmw_matrix_misc
  use hecmw_precond_BILU_33
  use hecmw_precond_DIAG_33
  use hecmw_precond_SSOR_33
  use hecmw_precond_ML_33
  use hecmw_precond_SAINV_33
  use hecmw_precond_RIF_33
  use hecmw_precond_nn
  use hecmw_solver_las_33
  implicit none

  private

  public :: hecmw_precond_33_setup
  public :: hecmw_precond_33_clear
  public :: hecmw_precond_33_apply

contains

  subroutine hecmw_precond_33_setup(hecMAT, hecMESH, sym)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: sym

    select case(hecmw_mat_get_precond( hecMAT ))
      case(1,2)
        call hecmw_precond_SSOR_33_setup(hecMAT)
      case(3)
        call hecmw_precond_DIAG_33_setup(hecMAT)
      case(5)
        call hecmw_precond_ML_33_setup(hecMAT, hecMESH, sym)
      case(10,11,12)
        call hecmw_precond_BILU_33_setup(hecMAT)
      case(20)
        call hecmw_precond_33_SAINV_setup(hecMAT)
      case(21)
        call hecmw_precond_RIF_33_setup(hecMAT)
      case default
        call hecmw_precond_nn_setup(hecMAT, hecMESH, sym)
    end select

  end subroutine hecmw_precond_33_setup

  subroutine hecmw_precond_33_clear(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecmw_mat_get_precond( hecMAT ))
      case(1,2)
        call hecmw_precond_SSOR_33_clear(hecMAT)
      case(3)
        call hecmw_precond_DIAG_33_clear()
      case(5)
        call hecmw_precond_ML_33_clear()
      case(10:12)
        call hecmw_precond_BILU_33_clear()
      case(20)
        call hecmw_precond_33_SAINV_clear()
      case(21)
        call hecmw_precond_RIF_33_clear()
      case default
        call hecmw_precond_nn_clear(hecMAT)
    end select

  end subroutine hecmw_precond_33_clear

  subroutine hecmw_precond_33_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: R(:)
    real(kind=kreal), intent(out) :: Z(:), ZP(:)
    real(kind=kreal), intent(inout) :: time_precond
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint ) :: i, iterPRE, iterPREmax
    real(kind=kreal) :: START_TIME, END_TIME

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    do iterPRE= 1, iterPREmax
      START_TIME = hecmw_Wtime()
      select case(hecmw_mat_get_precond( hecMAT ))
        case(1,2)
          call hecmw_precond_SSOR_33_apply(ZP)
        case(3)
          call hecmw_precond_DIAG_33_apply(ZP)
        case(5)
          call hecmw_precond_ML_33_apply(ZP)
        case(10:12)
          call hecmw_precond_BILU_33_apply(ZP)
        case(20)
          call hecmw_precond_33_SAINV_apply(R,ZP)
        case(21)
          call hecmw_precond_RIF_33_apply(ZP)
        case default
          call hecmw_precond_nn_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
          return
      end select
      END_TIME = hecmw_Wtime()
      time_precond = time_precond + END_TIME - START_TIME

      !C-- additive Schwartz
      do i= 1, hecMAT%N * hecMAT%NDOF
        Z(i)= Z(i) + ZP(i)
      enddo
      if (iterPRE.eq.iterPREmax) exit

      !C--    {ZP} = {R} - [A] {Z}
      call hecmw_matresid_33 (hecMESH, hecMAT, Z, R, ZP, COMMtime)
    enddo
  end subroutine hecmw_precond_33_apply
end module hecmw_precond_33
