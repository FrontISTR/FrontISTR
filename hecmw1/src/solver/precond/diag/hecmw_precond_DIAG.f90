!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_DIAG
  use hecmw_util
  use hecmw_precond_DIAG_11
  use hecmw_precond_DIAG_22
  use hecmw_precond_DIAG_33
  use hecmw_precond_DIAG_44
  use hecmw_precond_DIAG_66
  use hecmw_precond_DIAG_nn
  implicit none

  private
  public :: hecmw_precond_DIAG_setup
  public :: hecmw_precond_DIAG_clear
  public :: hecmw_precond_DIAG_apply

contains

  subroutine hecmw_precond_DIAG_setup(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_DIAG_33_setup(hecMAT)
      case(4)
        call hecmw_precond_DIAG_44_setup(hecMAT)
      case(6)
        call hecmw_precond_DIAG_66_setup(hecMAT)
      case(1)
        call hecmw_precond_DIAG_11_setup(hecMAT)
      case(2)
        call hecmw_precond_DIAG_22_setup(hecMAT)
      case default
        call hecmw_precond_DIAG_nn_setup(hecMAT)
    end select
  end subroutine hecmw_precond_DIAG_setup

  subroutine hecmw_precond_DIAG_clear(NDOF)
    implicit none
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_DIAG_33_clear()
      case(4)
        call hecmw_precond_DIAG_44_clear()
      case(6)
        call hecmw_precond_DIAG_66_clear()
      case(1)
        call hecmw_precond_DIAG_11_clear()
      case(2)
        call hecmw_precond_DIAG_22_clear()
      case default
        call hecmw_precond_DIAG_nn_clear()
    end select
  end subroutine hecmw_precond_DIAG_clear

  subroutine hecmw_precond_DIAG_apply(ZP, NDOF)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF

    select case(NDOF)
      case(3)
        call hecmw_precond_DIAG_33_apply(ZP)
      case(4)
        call hecmw_precond_DIAG_44_apply(ZP)
      case(6)
        call hecmw_precond_DIAG_66_apply(ZP)
      case(1)
        call hecmw_precond_DIAG_11_apply(ZP)
      case(2)
        call hecmw_precond_DIAG_22_apply(ZP)
      case default
        call hecmw_precond_DIAG_nn_apply(ZP,NDOF)
    end select
  end subroutine hecmw_precond_DIAG_apply

end module hecmw_precond_DIAG
