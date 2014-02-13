!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                       Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_solver_las_11

  private

  public :: hecmw_matvec_11
  public :: hecmw_matresid_11
  public :: hecmw_rel_resid_L2_11

contains

  !C
  !C***
  !C*** hecmw_matvec_11
  !C***
  !C
  subroutine hecmw_matvec_11 (hecMESH, hecMAT, X, Y, n, COMMtime)
    use hecmw_util
    use m_hecmw_comm_f

    implicit none
    integer(kind=kint):: n
    real(kind=kreal), dimension(n) :: X, Y
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    real(kind=kreal), optional :: COMMtime

    real(kind=kreal) :: START_TIME, END_TIME
    integer(kind=kint) :: i, j, jS, jE, in
    real(kind=kreal) :: YV

    START_TIME= HECMW_WTIME()
    call hecmw_update_1_R (hecMESH, X, hecMAT%NP)
    END_TIME= HECMW_WTIME()
    if (present(COMMtime)) COMMtime = COMMtime + END_TIME - START_TIME

    do i= 1, hecMAT%N
      YV= hecMAT%D(i) * X(i)
      jS= hecMAT%indexL(i-1) + 1
      jE= hecMAT%indexL(i  )
      do j= jS, jE
        in= hecMAT%itemL(j)
        YV= YV + hecMAT%AL(j) * X(in)
      enddo
      jS= hecMAT%indexU(i-1) + 1
      jE= hecMAT%indexU(i  )
      do j= jS, jE
        in= hecMAT%itemU(j)
        YV= YV + hecMAT%AU(j) * X(in)
      enddo
      Y(i)= YV
    enddo

  end subroutine hecmw_matvec_11

  !C
  !C***
  !C*** hecmw_matresid_11
  !C***
  !C
  subroutine hecmw_matresid_11 (hecMESH, hecMAT, X, B, R, COMMtime)
    use hecmw_util

    implicit none
    real(kind=kreal) :: X(:), B(:), R(:)
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    real(kind=kreal), optional :: COMMtime

    integer(kind=kint) :: i
    real(kind=kreal) :: Tcomm

    Tcomm = 0.d0
    call hecmw_matvec_11 (hecMESH, hecMAT, X, R, hecMAT%NP, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm

    do i = 1, hecMAT%N
      R(i) = B(i) - R(i)
    enddo

  end subroutine hecmw_matresid_11

  !C
  !C***
  !C*** hecmw_rel_resid_L2_11
  !C***
  !C
  function hecmw_rel_resid_L2_11 (hecMESH, hecMAT, COMMtime)
    use hecmw_util
    use hecmw_solver_misc

    implicit none
    real(kind=kreal) :: hecmw_rel_resid_L2_11
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH
    type ( hecmwST_matrix     ), intent(in) :: hecMAT
    real(kind=kreal), optional :: COMMtime

    real(kind=kreal) :: r(hecMAT%NDOF*hecMAT%NP)
    real(kind=kreal) :: bnorm2, rnorm2
    real(kind=kreal) :: Tcomm

    Tcomm = 0.d0
    call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnorm2, Tcomm)
      call hecmw_matresid_11(hecMESH, hecMAT, hecMAT%X, hecMAT%B, r, Tcomm)
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, r, r, rnorm2, Tcomm)
    if (present(COMMtime)) COMMtime = COMMtime + Tcomm
    hecmw_rel_resid_L2_11 = sqrt(rnorm2 / bnorm2)

  end function hecmw_rel_resid_L2_11

end module hecmw_solver_las_11
