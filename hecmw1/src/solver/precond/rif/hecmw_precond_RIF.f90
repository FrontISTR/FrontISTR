!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_RIF
  use hecmw_util
  use hecmw_precond_RIF_33
  use hecmw_precond_RIF_nn
  implicit none

  private
  public :: hecmw_precond_RIF_setup
  public :: hecmw_precond_RIF_clear
  public :: hecmw_precond_RIF_apply

contains

  subroutine hecmw_precond_RIF_setup(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_RIF_33_setup(hecMAT)
      case default
        call hecmw_precond_RIF_nn_setup(hecMAT)
    end select
  end subroutine hecmw_precond_RIF_setup

  subroutine hecmw_precond_RIF_clear(NDOF)
    implicit none
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_RIF_33_clear()
      case default
        call hecmw_precond_RIF_nn_clear()
    end select
  end subroutine hecmw_precond_RIF_clear

  subroutine hecmw_precond_RIF_apply(ZP, NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_RIF_33_apply(ZP)
      case default
        call hecmw_precond_RIF_nn_apply(ZP,NDOF)
    end select
  end subroutine hecmw_precond_RIF_apply

end module hecmw_precond_RIF
