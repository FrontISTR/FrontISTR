!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/11/02                                         
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                       Tatushiro Shono (Univ. of Tokyo)               !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
module hecmw_precond_44
  use hecmw_util
  implicit none

  private

  public :: hecmw_precond_44_setup
  public :: hecmw_precond_44_clear
  public :: hecmw_precond_44_apply
  public :: hecmw_precond_44_clear_timer
  public :: hecmw_precond_44_get_timer

  real(kind=kreal) :: time_precond = 0.d0

contains

  !C
  !C***
  !C*** hecmw_precond_44_setup
  !C***
  !C
  subroutine hecmw_precond_44_setup(hecMAT)
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_precond_BILU_44
!    use hecmw_precond_DIAG_44
    use hecmw_precond_SSOR_44
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    integer(kind=kint ) :: PRECOND, iterPREmax

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    if (iterPREmax.le.0) return

    PRECOND = hecmw_mat_get_precond( hecMAT )

    if (PRECOND.le.2) then
      call hecmw_precond_SSOR_44_setup(hecMAT)
!    else if (PRECOND.eq.3) then
!      call hecmw_precond_DIAG_44_setup(hecMAT)
    else if (PRECOND.eq.10.or.PRECOND.eq.11.or.PRECOND.eq.12) then
      call hecmw_precond_BILU_44_setup(hecMAT)
    endif

  end subroutine hecmw_precond_44_setup

  !C
  !C***
  !C*** hecmw_precond_44_clear
  !C***
  !C
  subroutine hecmw_precond_44_clear(hecMAT)
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_precond_BILU_44
!    use hecmw_precond_DIAG_44
    use hecmw_precond_SSOR_44
    implicit none
    type (hecmwST_matrix), intent(inout) :: hecMAT

    integer(kind=kint ) :: PRECOND, iterPREmax

    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )
    if (iterPREmax.le.0) return

    PRECOND = hecmw_mat_get_precond( hecMAT )

    if (PRECOND.le.2) then
      call hecmw_precond_SSOR_44_clear(hecMAT)
!    else if (PRECOND.eq.3) then
!      call hecmw_precond_DIAG_44_clear()
    else if (PRECOND.eq.10.or.PRECOND.eq.11.or.PRECOND.eq.12) then
      call hecmw_precond_BILU_44_clear()
    endif

  end subroutine hecmw_precond_44_clear

  !C
  !C***
  !C*** hecmw_precond_44_apply
  !C***
  !C
  subroutine hecmw_precond_44_apply(hecMESH, hecMAT, R, Z, ZP, COMMtime)
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_precond_BILU_44
!    use hecmw_precond_DIAG_44
    use hecmw_precond_SSOR_44
    use hecmw_solver_las_44
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    real(kind=kreal), intent(in) :: R(:)
    real(kind=kreal), intent(out) :: Z(:), ZP(:)
    real(kind=kreal), intent(inout) :: COMMtime

    integer(kind=kint ) :: N, NP, NNDOF, NPNDOF
    integer(kind=kint ) :: PRECOND, iterPREmax
    integer(kind=kint ) :: i, iterPRE
    real(kind=kreal) :: START_TIME, END_TIME

    N = hecMAT%N
    NP = hecMAT%NP
    NNDOF = N * 4
    NPNDOF = NP * 4
    PRECOND = hecmw_mat_get_precond( hecMAT )
    iterPREmax = hecmw_mat_get_iterpremax( hecMAT )

    if (iterPREmax.le.0) then
      do i= 1, NNDOF
        Z(i)= R(i)
      enddo
      return
    endif

    !C===
    !C +----------------+
    !C | {z}= [Minv]{r} |
    !C +----------------+
    !C===

    do i= 1, NNDOF
      ZP(i)= R(i)
    enddo
    do i= NNDOF+1, NPNDOF
      ZP(i) = 0.d0
    enddo

    do i= 1, NPNDOF
      Z(i)= 0.d0
    enddo

    do iterPRE= 1, iterPREmax

      START_TIME = hecmw_Wtime()

      if (PRECOND.le.2) then
        call hecmw_precond_SSOR_44_apply(ZP)
!      else if (PRECOND.eq.3) then
!        call hecmw_precond_DIAG_44_apply(ZP)
      else if (PRECOND.eq.10.or.PRECOND.eq.11.or.PRECOND.eq.12) then
        call hecmw_precond_BILU_44_apply(ZP)
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

      call hecmw_matresid_44 (hecMESH, hecMAT, Z, R, ZP, COMMtime)

    enddo

  end subroutine hecmw_precond_44_apply

  !C
  !C***
  !C*** hecmw_precond_44_clear_timer
  !C***
  !C
  subroutine hecmw_precond_44_clear_timer
    implicit none
    time_precond = 0.d0
  end subroutine hecmw_precond_44_clear_timer

  !C
  !C***
  !C*** hecmw_precond_44_get_timer
  !C***
  !C
  function hecmw_precond_44_get_timer()
    implicit none
    real(kind=kreal) :: hecmw_precond_44_get_timer
    hecmw_precond_44_get_timer = time_precond
  end function hecmw_precond_44_get_timer

end module hecmw_precond_44
