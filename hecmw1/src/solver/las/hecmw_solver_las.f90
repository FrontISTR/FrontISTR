!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_las
  use hecmw_util
  use hecmw_solver_las_11
  use hecmw_solver_las_22
  use hecmw_solver_las_33
  use hecmw_solver_las_44
  use hecmw_solver_las_66
  use hecmw_solver_las_nn

  implicit none

  private

  public :: hecmw_matvec
  public :: hecmw_matvec_set_async
  public :: hecmw_matvec_unset_async
  public :: hecmw_matresid
  public :: hecmw_rel_resid_L2
  public :: hecmw_Tvec
  public :: hecmw_Ttvec
  public :: hecmw_TtmatTvec
  public :: hecmw_matvec_clear_timer
  public :: hecmw_matvec_get_timer
  public :: hecmw_mat_diag_sr
  public :: hecmw_mat_add
  public :: hecmw_mat_multiple

  real(kind=kreal), save :: time_Ax = 0.d0

contains

  !C
  !C***
  !C*** hecmw_matvec
  !C***
  !C
  subroutine hecmw_matvec (hecMESH, hecMAT, X, Y, COMMtime)
    use hecmw_util

    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout), optional :: COMMtime
    select case(hecMAT%NDOF)
      case (3)
        call hecmw_matvec_33(hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
      case (4)
        call hecmw_matvec_44(hecMESH, hecMAT, X, Y, time_Ax,COMMtime)
      case (6)
        call hecmw_matvec_66(hecMESH, hecMAT, X, Y, time_Ax,COMMtime)
      case (1)
        call hecmw_matvec_11(hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
      case (2)
        call hecmw_matvec_22(hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
      case default
        call hecmw_matvec_nn(hecMESH, hecMAT, X, Y, time_Ax, COMMtime)
    end select

  end subroutine hecmw_matvec

  !C
  !C***
  !C*** hecmw_matvec_set_async
  !C***
  !C
  subroutine hecmw_matvec_set_async (hecMAT)
    use hecmw_util
    implicit none
    type (hecmwST_matrix), intent(in) :: hecMAT

  end subroutine hecmw_matvec_set_async

  !C
  !C***
  !C*** hecmw_matvec_unset_async
  !C***
  !C
  subroutine hecmw_matvec_unset_async
    implicit none
  end subroutine hecmw_matvec_unset_async

  !C
  !C***
  !C*** hecmw_matresid
  !C***
  !C
  subroutine hecmw_matresid (hecMESH, hecMAT, X, B, R, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:), B(:)
    real(kind=kreal), intent(out) :: R(:)
    real(kind=kreal), intent(inout), optional :: COMMtime

    select case(hecMAT%NDOF)
      case (3)
        call hecmw_matresid_33(hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
      case (4)
        call hecmw_matresid_44(hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
      case (6)
        call hecmw_matresid_66(hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
      case (1)
        call hecmw_matresid_11(hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
      case (2)
        call hecmw_matresid_22(hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
      case default
        call hecmw_matresid_nn(hecMESH, hecMAT, X, B, R, time_Ax, COMMtime)
    end select
  end subroutine hecmw_matresid

  !C
  !C***
  !C*** hecmw_rel_resid_L2
  !C***
  !C
  function hecmw_rel_resid_L2 (hecMESH, hecMAT, COMMtime)
    use hecmw_util
    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH
    type ( hecmwST_matrix     ), intent(in) :: hecMAT
    real(kind=kreal), intent(inout), optional :: COMMtime

    select case(hecMAT%NDOF)
      case (3)
        hecmw_rel_resid_L2 = hecmw_rel_resid_L2_33(hecMESH, hecMAT, time_Ax, COMMtime)
      case (4)
        hecmw_rel_resid_L2 = hecmw_rel_resid_L2_44(hecMESH, hecMAT, time_Ax, COMMtime)
      case (6)
        hecmw_rel_resid_L2 = hecmw_rel_resid_L2_66(hecMESH, hecMAT, time_Ax, COMMtime)
      case (1)
        hecmw_rel_resid_L2 = hecmw_rel_resid_L2_11(hecMESH, hecMAT, time_Ax, COMMtime)
      case (2)
        hecmw_rel_resid_L2 = hecmw_rel_resid_L2_22(hecMESH, hecMAT, time_Ax, COMMtime)
      case default
        hecmw_rel_resid_L2 = hecmw_rel_resid_L2_nn(hecMESH, hecMAT, time_Ax, COMMtime)
    end select
  end function hecmw_rel_resid_L2

  !C
  !C***
  !C*** hecmw_Tvec
  !C***
  !C
  subroutine hecmw_Tvec (hecMESH, ndof, X, Y, COMMtime)
    use hecmw_util
    use hecmw_solver_las_33
    use hecmw_solver_las_nn

    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    select case(ndof)
      case (3)
        call hecmw_Tvec_33(hecMESH, X, Y, COMMtime)
      case default
        call hecmw_Tvec_nn(hecMESH, ndof, X, Y, COMMtime)
    end select

  end subroutine hecmw_Tvec

  !C
  !C***
  !C*** hecmw_Ttvec
  !C***
  !C
  subroutine hecmw_Ttvec (hecMESH, ndof, X, Y, COMMtime)
    use hecmw_util
    use hecmw_solver_las_33
    use hecmw_solver_las_nn
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in) :: ndof
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime

    select case(ndof)
      case (3)
        call hecmw_Ttvec_33(hecMESH, X, Y, COMMtime)
      case default
        call hecmw_Ttvec_nn(hecMESH, ndof, X, Y, COMMtime)
    end select

  end subroutine hecmw_Ttvec

  !C
  !C***
  !C*** hecmw_TtmatTvec
  !C***
  !C
  subroutine hecmw_TtmatTvec (hecMESH, hecMAT, X, Y, W, COMMtime)
    use hecmw_util
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:), W(:)
    real(kind=kreal), intent(inout) :: COMMtime

    !    call hecmw_Tvec(hecMESH, X, Y, COMMtime)
    !    call hecmw_matvec(hecMESH, hecMAT, Y, W, COMMtime)
    !    call hecmw_Ttvec(hecMESH, W, Y, COMMtime)
    select case(hecMESH%n_dof)
      case (3)
        call hecmw_TtmatTvec_33 (hecMESH, hecMAT, X, Y, W, time_Ax, COMMtime)
      case default
        call hecmw_TtmatTvec_nn (hecMESH, hecMAT, X, Y, W, time_Ax, COMMtime)
    end select

  end subroutine hecmw_TtmatTvec

  !C
  !C***
  !C*** hecmw_matvec_clear_timer
  !C***
  !C
  subroutine hecmw_matvec_clear_timer
    implicit none
    time_Ax = 0.d0
  end subroutine hecmw_matvec_clear_timer

  !C
  !C***
  !C*** hecmw_matvec_get_timer
  !C***
  !C
  function hecmw_matvec_get_timer()
    implicit none
    real(kind=kreal) :: hecmw_matvec_get_timer
    hecmw_matvec_get_timer = time_Ax
  end function hecmw_matvec_get_timer

  !C
  !C***
  !C*** hecmw_mat_diag_sr
  !C***
  !C
  subroutine hecmw_mat_diag_sr(hecMESH, hecMAT, COMMtime)
    use hecmw_util
    use hecmw_solver_las_33
    use hecmw_solver_las_nn
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(inout), target :: hecMAT
    real(kind=kreal), intent(inout), optional :: COMMtime

    select case(hecMESH%n_dof)
      case (3)
        call hecmw_mat_diag_sr_33(hecMESH, hecMAT, COMMtime)
      case default
        call hecmw_mat_diag_sr_nn(hecMESH, hecMAT, COMMtime)
    end select

  end subroutine hecmw_mat_diag_sr

  subroutine hecmw_mat_add(hecMAT1, hecMAT2, hecMAT3)
    use hecmw_util
    use hecmw_solver_las_33
    use hecmw_solver_las_nn
    implicit none
    type (hecmwST_matrix), intent(inout), target :: hecMAT1, hecMAT2, hecMAT3

    call hecmw_mat_add_nn(hecMAT1, hecMAT2, hecMAT3)
  end subroutine hecmw_mat_add

  subroutine hecmw_mat_multiple(hecMAT, alpha)
    use hecmw_util
    use hecmw_solver_las_33
    use hecmw_solver_las_nn
    implicit none
    type (hecmwST_matrix), intent(inout), target :: hecMAT
    real(kind=kreal), intent(in) :: alpha

    call hecmw_mat_multiple_nn(hecMAT, alpha)
  end subroutine hecmw_mat_multiple
end module hecmw_solver_las
