!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_SSOR
  use hecmw_util
  use hecmw_precond_SSOR_11
  use hecmw_precond_SSOR_22
  use hecmw_precond_SSOR_33
  use hecmw_precond_SSOR_44
  use hecmw_precond_SSOR_66
  use hecmw_precond_SSOR_nn
  implicit none

  private
  public :: hecmw_precond_SSOR_setup
  public :: hecmw_precond_SSOR_clear
  public :: hecmw_precond_SSOR_apply

contains

  subroutine hecmw_precond_SSOR_setup(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_SSOR_33_setup(hecMAT)
      case(4)
        call hecmw_precond_SSOR_44_setup(hecMAT)
      case(6)
        call hecmw_precond_SSOR_66_setup(hecMAT)
      case(1)
        call hecmw_precond_SSOR_11_setup(hecMAT)
      case(2)
        call hecmw_precond_SSOR_22_setup(hecMAT)
      case default
        call hecmw_precond_SSOR_nn_setup(hecMAT)
    end select
  end subroutine hecmw_precond_SSOR_setup

  subroutine hecmw_precond_SSOR_clear(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_SSOR_33_clear(hecMAT)
      case(4)
        call hecmw_precond_SSOR_44_clear(hecMAT)
      case(6)
        call hecmw_precond_SSOR_66_clear(hecMAT)
      case(1)
        call hecmw_precond_SSOR_11_clear(hecMAT)
      case(2)
        call hecmw_precond_SSOR_22_clear(hecMAT)
      case default
        call hecmw_precond_SSOR_nn_clear(hecMAT)
    end select
  end subroutine hecmw_precond_SSOR_clear

  subroutine hecmw_precond_SSOR_apply(ZP, NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_SSOR_33_apply(ZP)
      case(4)
        call hecmw_precond_SSOR_44_apply(ZP)
      case(6)
        call hecmw_precond_SSOR_66_apply(ZP)
      case(1)
        call hecmw_precond_SSOR_11_apply(ZP)
      case(2)
        call hecmw_precond_SSOR_22_apply(ZP)
      case default
        call hecmw_precond_SSOR_nn_apply(ZP,NDOF)
    end select
  end subroutine hecmw_precond_SSOR_apply

end module hecmw_precond_SSOR
