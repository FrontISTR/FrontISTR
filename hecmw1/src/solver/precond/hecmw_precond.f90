!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond
  use hecmw_util
  use hecmw_precond_SSOR
  use hecmw_precond_DIAG
  use hecmw_precond_ML
  use hecmw_precond_BILU
  use hecmw_precond_SAINV
  use hecmw_precond_RIF
  use hecmw_matrix_misc
  use hecmw_solver_las
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

    select case(hecmw_mat_get_precond( hecMAT ))
      case(1,2)
        call hecmw_precond_SSOR_setup(hecMAT)
      case(3)
        call hecmw_precond_DIAG_setup(hecMAT)
      case(5)
        call hecmw_precond_ML_setup(hecMAT, hecMESH, sym)
      case(10,11,12)
        call hecmw_precond_BILU_setup(hecMAT)
      case(20)
        call hecmw_precond_SAINV_setup(hecMAT)
      case(21)
        call hecmw_precond_RIF_setup(hecMAT)
      case default
        write (*,'(/a )')'#### HEC-MW-SOLVER-E-1001'
        write (*,'( a/)')'    inconsistent solver/preconditioning'
        call hecmw_abort( hecmw_comm_get_comm())
    end select
  end subroutine hecmw_precond_setup

  subroutine hecmw_precond_clear(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    if (hecmw_mat_get_iterpremax( hecMAT ).le.0) return

    select case(hecmw_mat_get_precond( hecMAT ))
      case(1,2)
        call hecmw_precond_SSOR_clear(hecMAT)
      case(3)
        call hecmw_precond_DIAG_clear(hecMAT%NDOF)
      case(5)
        call hecmw_precond_ML_clear(hecMAT%NDOF)
      case(10:12)
        call hecmw_precond_BILU_clear(hecMAT%NDOF)
      case(20)
        call hecmw_precond_SAINV_clear(hecMAT%NDOF)
      case(21)
        call hecmw_precond_RIF_clear(hecMAT%NDOF)
      case default
    end select

  end subroutine hecmw_precond_clear

  subroutine hecmw_precond_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix), intent(inout)     :: hecMAT
    real(kind=kreal), intent(inout) :: R(:)
    real(kind=kreal), intent(inout) :: Z(:), ZP(:)
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint ) :: i, N, NP, NNDOF, NPNDOF
    integer(kind=kint) :: iterPREmax, iterPRE
    real(kind=kreal) :: START_TIME, END_TIME

    START_TIME = hecmw_Wtime()

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

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    do iterPRE= 1, iterPREmax

      select case(hecmw_mat_get_precond( hecMAT ))
        case(1,2)
          call hecmw_precond_SSOR_apply(ZP,hecMAT%NDOF)
        case(3)
          call hecmw_precond_DIAG_apply(ZP,hecMAT%NDOF)
        case(5)
          call hecmw_precond_ML_apply(ZP,hecMAT%NDOF)
        case(10:12)
          call hecmw_precond_BILU_apply(ZP,hecMAT%NDOF)
        case(20)
          call hecmw_precond_SAINV_apply(R,ZP,hecMAT%NDOF)
        case(21)
          call hecmw_precond_RIF_apply(ZP,hecMAT%NDOF)
        case default
      end select

      !C-- additive Schwartz
      do i= 1, hecMAT%N * hecMAT%NDOF
        Z(i)= Z(i) + ZP(i)
      enddo
      if (iterPRE.eq.iterPREmax) exit

      !C--    {ZP} = {R} - [A] {Z}
      call hecmw_matresid (hecMESH, hecMAT, Z, R, ZP, COMMtime)
    enddo

    END_TIME = hecmw_Wtime()
    time_precond = time_precond + END_TIME - START_TIME
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
