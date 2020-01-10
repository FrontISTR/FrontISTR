!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond
  use hecmw_util
  use hecmw_precond_11
  use hecmw_precond_22
  use hecmw_precond_33
  use hecmw_precond_44
  use hecmw_precond_66
  use hecmw_precond_nn
  use hecmw_matrix_misc
  implicit none

  private
  public :: hecmw_precond_setup
  public :: hecmw_precond_clear
  public :: hecmw_precond_apply
  public :: hecmw_precond_clear_timer
  public :: hecmw_precond_get_timer

  real(kind=kreal) :: time_precond = 0.d0

contains

  subroutine hecmw_precond_setup(hecMAT, hecMESH, sym)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint) :: sym

    if (hecmw_mat_get_iterpremax( hecMAT ).le.0) return

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_33_setup(hecMAT, hecMESH, sym)
      case(4)
        call hecmw_precond_44_setup(hecMAT, hecMESH, sym)
      case(6)
        call hecmw_precond_66_setup(hecMAT, hecMESH, sym)
      case(1)
        call hecmw_precond_11_setup(hecMAT, hecMESH, sym)
      case(2)
        call hecmw_precond_22_setup(hecMAT, hecMESH, sym)
      case default
        call hecmw_precond_nn_setup(hecMAT, hecMESH, sym)
    end select
  end subroutine hecmw_precond_setup

  subroutine hecmw_precond_clear(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    if (hecmw_mat_get_iterpremax( hecMAT ).le.0) return

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_33_clear(hecMAT)
      case(4)
        call hecmw_precond_44_clear(hecMAT)
      case(6)
        call hecmw_precond_66_clear(hecMAT)
      case(1)
        call hecmw_precond_11_clear(hecMAT)
      case(2)
        call hecmw_precond_22_clear(hecMAT)
      case default
        call hecmw_precond_nn_clear(hecMAT)
    end select

  end subroutine hecmw_precond_clear

  subroutine hecmw_precond_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(inout)     :: hecMAT
    real(kind=kreal), intent(in)    :: R(hecMAT%N*hecMAT%NDOF)
    real(kind=kreal), intent(inout) :: Z(hecMAT%NP*hecMAT%NDOF)
    real(kind=kreal), intent(inout) :: ZP(hecMAT%NP*hecMAT%NDOF)
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint ) :: N, NP, NNDOF, NPNDOF
    integer(kind=kint ) :: i, iterPRE
    real(kind=kreal) :: START_TIME, END_TIME
    N = hecMAT%N
    NP = hecMAT%NP
    NNDOF = N * hecMAT%NDOF
    NPNDOF = NP * hecMAT%NDOF

    if (hecmw_mat_get_iterpremax( hecMAT ).le.0) then
      do i= 1, NNDOF
        Z(i)= R(i)
      enddo
      return
    endif

    !C {z}= [Minv]{r}
    do i= 1, NNDOF
      ZP(i)= R(i)
    enddo
    do i= NNDOF+1, NPNDOF
      ZP(i) = 0.d0
    enddo
    do i= 1, NPNDOF
      Z(i)= 0.d0
    enddo

    select case(hecMAT%NDOF)
      case(3)
        call hecmw_precond_33_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
      case(4)
        call hecmw_precond_44_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
      case(6)
        call hecmw_precond_66_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
      case(1)
        call hecmw_precond_11_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
      case(2)
        call hecmw_precond_22_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
      case default
        call hecmw_precond_nn_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
    end select

  end subroutine hecmw_precond_apply

  subroutine hecmw_precond_clear_timer
    implicit none
    time_precond = 0.d0
  end subroutine hecmw_precond_clear_timer
  function hecmw_precond_get_timer()
    implicit none
    real(kind=kreal) :: hecmw_precond_get_timer
    hecmw_precond_get_timer = time_precond
  end function hecmw_precond_get_timer
end module hecmw_precond
