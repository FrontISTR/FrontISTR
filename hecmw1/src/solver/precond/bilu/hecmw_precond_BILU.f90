!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_BILU
  use hecmw_util
  use hecmw_precond_BILU_33
  use hecmw_precond_BILU_44
  use hecmw_precond_BILU_66
  use hecmw_precond_BILU_nn
  implicit none

  private
  public :: hecmw_precond_BILU_setup
  public :: hecmw_precond_BILU_clear
  public :: hecmw_precond_BILU_apply

contains

  subroutine hecmw_precond_BILU_setup(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_BILU_33_setup(hecMAT)
      case(4)
        call hecmw_precond_BILU_44_setup(hecMAT)
      case(6)
        call hecmw_precond_BILU_66_setup(hecMAT)
      case default
        call hecmw_precond_BILU_nn_setup(hecMAT)
    end select
  end subroutine hecmw_precond_BILU_setup

  subroutine hecmw_precond_BILU_clear(NDOF)
    implicit none
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_BILU_33_clear()
      case(4)
        call hecmw_precond_BILU_44_clear()
      case(6)
        call hecmw_precond_BILU_66_clear()
      case default
        call hecmw_precond_BILU_nn_clear()
    end select
  end subroutine hecmw_precond_BILU_clear

  subroutine hecmw_precond_BILU_apply(ZP, NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_BILU_33_apply(ZP)
      case(4)
        call hecmw_precond_BILU_44_apply(ZP)
      case(6)
        call hecmw_precond_BILU_66_apply(ZP)
      case default
        call hecmw_precond_BILU_nn_apply(ZP,NDOF)
    end select
  end subroutine hecmw_precond_BILU_apply

end module hecmw_precond_BILU
