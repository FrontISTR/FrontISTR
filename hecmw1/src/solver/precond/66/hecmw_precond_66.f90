!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_precond_66
  use hecmw_util
  implicit none

  private

  public :: hecmw_precond_66_setup
  public :: hecmw_precond_66_clear
  public :: hecmw_precond_66_apply

contains

  !C
  !C***
  !C*** hecmw_precond_66_setup
  !C***
  !C
  subroutine hecmw_precond_66_setup(hecMAT)
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_precond_BILU_66
    use hecmw_precond_DIAG_66
    use hecmw_precond_SSOR_66
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    integer(kind=kint ) :: PRECOND, iterPREmax

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    if (iterPREmax.le.0) return

    PRECOND = hecmw_mat_get_precond( hecMAT )

    if (PRECOND.le.2) then
      call hecmw_precond_SSOR_66_setup(hecMAT)
    else if (PRECOND.eq.3) then
      call hecmw_precond_DIAG_66_setup(hecMAT)
    else if (PRECOND.eq.10 .or. PRECOND.eq.11 .or. PRECOND.eq.12) then
      call hecmw_precond_BILU_66_setup(hecMAT)
    endif

  end subroutine hecmw_precond_66_setup

  !C
  !C***
  !C*** hecmw_precond_66_clear
  !C***
  !C
  subroutine hecmw_precond_66_clear(hecMAT)
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_precond_BILU_66
    use hecmw_precond_DIAG_66
    use hecmw_precond_SSOR_66
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    integer(kind=kint ) :: PRECOND, iterPREmax

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    if (iterPREmax.le.0) return

    PRECOND = hecmw_mat_get_precond( hecMAT )

    if (PRECOND.le.2) then
      call hecmw_precond_SSOR_66_clear(hecMAT)
    else if (PRECOND.eq.3) then
      call hecmw_precond_DIAG_66_clear()
    else if (PRECOND.eq.10 .or. PRECOND.eq.11 .or. PRECOND.eq.12) then
      call hecmw_precond_BILU_66_clear()
    endif

  end subroutine hecmw_precond_66_clear

  !C
  !C***
  !C*** hecmw_precond_66_apply
  !C***
  !C
  subroutine hecmw_precond_66_apply(hecMESH, hecMAT, R, Z, ZP, time_precond, COMMtime)
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_precond_BILU_66
    use hecmw_precond_DIAG_66
    use hecmw_precond_SSOR_66
    use hecmw_solver_las_66
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix) :: hecMAT
    real(kind=kreal), intent(inout) :: R(:)
    real(kind=kreal), intent(out) :: Z(:), ZP(:)
    real(kind=kreal), intent(inout) :: time_precond
    real(kind=kreal), intent(inout) :: COMMtime

    integer(kind=kint ) :: N, NP, NNDOF, NPNDOF
    integer(kind=kint ) :: PRECOND, iterPREmax
    integer(kind=kint ) :: i, iterPRE
    real(kind=kreal) :: START_TIME, END_TIME

    N = hecMAT%N
    NP = hecMAT%NP
    NNDOF = N * 6
    NPNDOF = NP * 6
    PRECOND = hecmw_mat_get_precond( hecMAT )
    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )

    if (iterPREmax.le.0) then
      do i= 1, NNDOF
        Z(i)= R(i)
      enddo
      return
    endif

    do i= 1, NNDOF
      ZP(i)= R(i)
    enddo
    !C===
    !C +----------------+
    !C | {z}= [Minv]{r} |
    !C +----------------+
    !C===

    do i= NNDOF+1, NPNDOF
      ZP(i) = 0.d0
    enddo

    do i= 1, NPNDOF
      Z(i)= 0.d0
    enddo

    do iterPRE= 1, iterPREmax

      START_TIME = hecmw_Wtime()

      if (PRECOND.le.2) then
        call hecmw_precond_SSOR_66_apply(ZP)
      else if (PRECOND.eq.3) then
        call hecmw_precond_DIAG_66_apply(ZP)
      else if (PRECOND.eq.10 .or. PRECOND.eq.11 .or. PRECOND.eq.12) then
        call hecmw_precond_BILU_66_apply(ZP)
      endif

      END_TIME = hecmw_Wtime()
      time_precond = time_precond + END_TIME - START_TIME

      !C
      !C-- additive Schwartz
      !C
      do i= 1, NNDOF
        Z(i)= Z(i) + ZP(i)
      enddo

      if (iterPRE.eq.iterPREmax) exit

      !C--    {ZP} = {R} - [A] {Z}

      call hecmw_matresid_66 (hecMESH, hecMAT, Z, R, ZP, COMMtime)

    enddo

  end subroutine hecmw_precond_66_apply

end module hecmw_precond_66
