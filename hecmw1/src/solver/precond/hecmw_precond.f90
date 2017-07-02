!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond
  use hecmw_util
  implicit none

  private

  public :: hecmw_precond_setup
  public :: hecmw_precond_clear
  public :: hecmw_precond_apply
  public :: hecmw_precond_clear_timer
  public :: hecmw_precond_get_timer

  real(kind=kreal) :: time_precond = 0.d0

contains

  !C
  !C***
  !C*** hecmw_precond_setup
  !C***
  !C
  subroutine hecmw_precond_setup(hecMAT, hecMESH, sym)
    use hecmw_util
    use hecmw_precond_33
    use hecmw_precond_44
    use hecmw_precond_66
    use hecmw_precond_nn
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint) :: sym

    SELECT CASE(hecMAT%NDOF)
      CASE(3)
        call hecmw_precond_33_setup(hecMAT, hecMESH, sym)
      CASE(4)
        call hecmw_precond_44_setup(hecMAT, hecMESH, sym)
!      CASE(6)
!        call hecmw_precond_66_setup(hecMAT)
      case default
        call hecmw_precond_nn_setup(hecMAT, hecMESH, sym)
    END SELECT
  end subroutine hecmw_precond_setup

  !C
  !C***
  !C*** hecmw_precond_clear
  !C***
  !C
  subroutine hecmw_precond_clear(hecMAT)
    use hecmw_util
    use hecmw_precond_33
    use hecmw_precond_44
    use hecmw_precond_66
    use hecmw_precond_nn
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
 
    SELECT CASE(hecMAT%NDOF)
      CASE(3)
        call hecmw_precond_33_clear(hecMAT)
      CASE(4)
        call hecmw_precond_44_clear(hecMAT)
!      CASE(6)
!        call hecmw_precond_66_clear(hecMAT)
      case default
        call hecmw_precond_nn_clear(hecMAT)
    END SELECT

  end subroutine hecmw_precond_clear

  !C
  !C***
  !C*** hecmw_precond_apply
  !C***
  !C
  subroutine hecmw_precond_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
    use hecmw_util
    use hecmw_precond_33
    use hecmw_precond_44
    use hecmw_precond_66
    use hecmw_precond_nn
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(inout)     :: hecMAT
    real(kind=kreal), intent(inout) :: R(:)
    real(kind=kreal), intent(inout) :: Z(:), ZP(:)
    real(kind=kreal), intent(inout) :: COMMtime


    call hecmw_precond_nn_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
    SELECT CASE(hecMAT%NDOF)
      CASE(3)
        call hecmw_precond_33_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
      CASE(4)
        call hecmw_precond_44_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
!      CASE(6)
!        call hecmw_precond_66_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
      case default
        call hecmw_precond_nn_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
    END SELECT

  end subroutine hecmw_precond_apply

  !C
  !C***
  !C*** hecmw_precond_clear_timer
  !C***
  !C
  subroutine hecmw_precond_clear_timer
    implicit none
    time_precond = 0.d0
  end subroutine hecmw_precond_clear_timer

  !C
  !C***
  !C*** hecmw_precond_get_timer
  !C***
  !C
  function hecmw_precond_get_timer()
    implicit none
    real(kind=kreal) :: hecmw_precond_get_timer
    hecmw_precond_get_timer = time_precond
  end function hecmw_precond_get_timer

end module hecmw_precond
