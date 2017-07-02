!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_11
  use hecmw_util
  use hecmw_matrix_misc
  !use hecmw_precond_BILU_11
  use hecmw_precond_DIAG_11
  use hecmw_precond_SSOR_11
  use hecmw_precond_nn
  use hecmw_solver_las_11
  implicit none

  private

  public :: hecmw_precond_11_setup
  public :: hecmw_precond_11_clear
  public :: hecmw_precond_11_apply

contains

  subroutine hecmw_precond_11_setup(hecMAT, hecMESH, sym)
    implicit none
    type (hecmwST_matrix),     intent(inout) :: hecMAT
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint), intent(in) :: sym

    SELECT CASE(hecmw_mat_get_precond( hecMAT ))
      CASE(1,2)
        call hecmw_precond_SSOR_11_setup(hecMAT)
      CASE(3)
        call hecmw_precond_DIAG_11_setup(hecMAT)
      !CASE(10,11,12)
      !  call hecmw_precond_BILU_11_setup(hecMAT)
      CASE DEFAULT
        call hecmw_precond_nn_setup(hecMAT, hecMESH, sym)
    END SELECT

  end subroutine hecmw_precond_11_setup

  subroutine hecmw_precond_11_clear(hecMAT)
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    SELECT CASE(hecmw_mat_get_precond( hecMAT ))
      CASE(1,2)
        call hecmw_precond_SSOR_11_clear(hecMAT)
      CASE(3)
        call hecmw_precond_DIAG_11_clear()
      !CASE(10,11,12)
      !  call hecmw_precond_BILU_11_clear()
      CASE DEFAULT
        call hecmw_precond_nn_clear(hecMAT)
    END SELECT
    
  end subroutine hecmw_precond_11_clear

  subroutine hecmw_precond_11_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: R(:)
    real(kind=kreal), intent(out) :: Z(:), ZP(:)
    real(kind=kreal), intent(inout) :: time_precond
    real(kind=kreal), intent(inout) :: COMMtime
    integer(kind=kint ) :: i, iterPRE, iterPREmax
    real(kind=kreal) :: START_TIME, END_TIME

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    do iterPRE= 1, iterPREmax
      START_TIME = hecmw_Wtime()
      SELECT CASE(hecmw_mat_get_precond( hecMAT ))
        CASE(1,2)
          call hecmw_precond_SSOR_11_apply(ZP)
        CASE(3)
          call hecmw_precond_DIAG_11_apply(ZP)
        !CASE(10,11,12)
        !  call hecmw_precond_BILU_11_apply(ZP)
        CASE DEFAULT
          call hecmw_precond_nn_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
          return 
      END SELECT
      END_TIME = hecmw_Wtime()
      time_precond = time_precond + END_TIME - START_TIME

      !C-- additive Schwartz
      do i= 1, hecMAT%N
        Z(i)= Z(i) + ZP(i)
      enddo
      if (iterPRE.eq.iterPREmax) exit

      !C--    {ZP} = {R} - [A] {Z}
      call hecmw_matresid_11 (hecMESH, hecMAT, Z, R, ZP, COMMtime)
    enddo
  end subroutine hecmw_precond_11_apply
end module hecmw_precond_11
