!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_ML
  use hecmw_util
  use hecmw_precond_ML_33
  use hecmw_precond_ML_nn
  implicit none

  private
  public :: hecmw_precond_ML_setup
  public :: hecmw_precond_ML_clear
  public :: hecmw_precond_ML_apply

contains

  subroutine hecmw_precond_ML_setup(hecMAT, hecMESH, sym)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: sym

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_ML_33_setup(hecMAT, hecMESH, sym)
      case default
        call hecmw_precond_ML_nn_setup(hecMAT, hecMESH, sym)
    end select
  end subroutine hecmw_precond_ML_setup

  subroutine hecmw_precond_ML_clear(NDOF)
    implicit none
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_ML_33_clear()
      case default
        call hecmw_precond_ML_nn_clear()
    end select
  end subroutine hecmw_precond_ML_clear

  subroutine hecmw_precond_ML_apply(ZP, NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_ML_33_apply(ZP)
      case default
        call hecmw_precond_ML_nn_apply(ZP)
    end select
  end subroutine hecmw_precond_ML_apply

end module hecmw_precond_ML
