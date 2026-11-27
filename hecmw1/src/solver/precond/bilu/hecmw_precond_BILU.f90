!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_BILU
  use hecmw_util
  use hecmw_matrix_misc
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
        ! BILU_44 forms only ILU(0); ILU(1) and ILU(2) of 4x4 blocks are formed by BILU_nn.
        ! clear and apply must choose the same module, since each module keeps its own factors.
        if (hecmw_mat_get_precond(hecMAT) == 10) then
          call hecmw_precond_BILU_44_setup(hecMAT)
        else
          call hecmw_precond_BILU_nn_setup(hecMAT)
        endif
      case(6)
        call hecmw_precond_BILU_66_setup(hecMAT)
      case default
        call hecmw_precond_BILU_nn_setup(hecMAT)
    end select
  end subroutine hecmw_precond_BILU_setup

  subroutine hecmw_precond_BILU_clear(NDOF, PRECOND)
    implicit none
    integer(kind=kint), intent(in) :: NDOF, PRECOND

    select case(NDOF)
      case(3)
        call hecmw_precond_BILU_33_clear()
      case(4)
        if (PRECOND == 10) then
          call hecmw_precond_BILU_44_clear()
        else
          call hecmw_precond_BILU_nn_clear()
        endif
      case(6)
        call hecmw_precond_BILU_66_clear()
      case default
        call hecmw_precond_BILU_nn_clear()
    end select
  end subroutine hecmw_precond_BILU_clear

  subroutine hecmw_precond_BILU_apply(ZP, NDOF, PRECOND)
    implicit none
    real(kind=kreal), intent(inout) :: ZP(:)
    integer(kind=kint), intent(in) :: NDOF, PRECOND

    select case(NDOF)
      case(3)
        call hecmw_precond_BILU_33_apply(ZP)
      case(4)
        if (PRECOND == 10) then
          call hecmw_precond_BILU_44_apply(ZP)
        else
          call hecmw_precond_BILU_nn_apply(ZP,NDOF)
        endif
      case(6)
        call hecmw_precond_BILU_66_apply(ZP)
      case default
        call hecmw_precond_BILU_nn_apply(ZP,NDOF)
    end select
  end subroutine hecmw_precond_BILU_apply

end module hecmw_precond_BILU
