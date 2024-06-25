!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_SAINV
  use hecmw_util
  use hecmw_precond_SAINV_33
  use hecmw_precond_SAINV_nn
  implicit none

  private
  public :: hecmw_precond_SAINV_setup
  public :: hecmw_precond_SAINV_clear
  public :: hecmw_precond_SAINV_apply

contains

  subroutine hecmw_precond_SAINV_setup(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_SAINV_33_setup(hecMAT)
      case default
        call hecmw_precond_SAINV_nn_setup(hecMAT)
    end select
  end subroutine hecmw_precond_SAINV_setup

  subroutine hecmw_precond_SAINV_clear(NDOF)
    implicit none
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_SAINV_33_clear()
      case default
        call hecmw_precond_SAINV_nn_clear()
    end select
  end subroutine hecmw_precond_SAINV_clear

  subroutine hecmw_precond_SAINV_apply(R, ZP, NDOF)
    implicit none
    real(kind=kreal), intent(in) :: R(:)
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_SAINV_33_apply(R,ZP)
      case default
        call hecmw_precond_SAINV_nn_apply(R,ZP)
    end select
  end subroutine hecmw_precond_SAINV_apply

end module hecmw_precond_SAINV
