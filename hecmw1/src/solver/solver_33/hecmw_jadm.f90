!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
!> Jagged Diagonal Matrix storage for vector processors.
!> Original code was provided by JAMSTEC.

module hecmw_JAD_TYPE
  use hecmw_util
  use m_hecmw_comm_f
  implicit none

  private

  public :: hecmw_JAD_INIT
  public :: hecmw_JAD_FINALIZE
  public :: hecmw_JAD_MATVEC

  !C---------------------- AU&AL
  real(kind=kreal), allocatable      :: AJAD(:)
  integer(kind=kint), allocatable    :: JAJAD(:)
  integer(kind=kint), allocatable    :: JADORD(:)
  integer(kind=kint), allocatable    :: IAJAD(:)
  integer(kind=kint) :: MJAD
  real(kind=kreal), allocatable  :: WP1(:), WP2(:), WP3(:)

contains

  subroutine hecmw_JAD_INIT(hecMAT)
    type(hecmwST_matrix) :: hecMAT
    ALLOCATE(WP1(hecMAT%NP), WP2(hecMAT%NP), WP3(hecMAT%NP))
    ALLOCATE(AJAD((hecMAT%NPL+hecMAT%NPU)*9))
    ALLOCATE(JAJAD(hecMAT%NPL+hecMAT%NPU))
    ALLOCATE(JADORD(hecMAT%NP))
    ALLOCATE(IAJAD(hecMAT%NP+1))
    CALL REPACK(hecMAT%N, hecMAT, MJAD, AJAD, JAJAD, IAJAD, JADORD)
  end subroutine hecmw_JAD_INIT

  subroutine hecmw_JAD_FINALIZE()
    DEALLOCATE(AJAD)
    DEALLOCATE(JAJAD)
    DEALLOCATE(JADORD)
    DEALLOCATE(IAJAD)
    DEALLOCATE(WP1,WP2,WP3)
  end subroutine hecmw_JAD_FINALIZE

  subroutine hecmw_JAD_MATVEC(hecMESH, hecMAT, X, Y, COMMtime)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime
    real(kind=kreal) :: START_TIME, END_TIME
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: i
    real(kind=kreal) :: X1, X2, X3

    START_TIME= HECMW_WTIME()
    call hecmw_update_3_R (hecMESH, X, hecMAT%NP)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    D => hecMAT%D

    do i= 1, hecMAT%N
      X1= X(3*i-2)
      X2= X(3*i-1)
      X3= X(3*i  )
      Y(3*i -2)= D(9*i-8)*X1 + D(9*i-7)*X2 + D(9*i-6)*X3
      Y(3*i -1)= D(9*i-5)*X1 + D(9*i-4)*X2 + D(9*i-3)*X3
      Y(3*i   )= D(9*i-2)*X1 + D(9*i-1)*X2 + D(9*i  )*X3
    enddo
    CALL MATJAD(hecMAT%N, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, WP1, WP2, WP3)
  end subroutine hecmw_JAD_MATVEC

  subroutine REPACK(N, hecMAT, MJAD,  AJAD, JAJAD, IAJAD, JADORD)
    use hecmw_util
    !C---------------------------------
    type (hecmwST_matrix)     :: hecMAT
    !C----------------------
    integer(kind = kint) :: N, MJAD
    real(kind = kreal), dimension(*) :: AJAD
    integer(kind = kint), dimension(*) :: JAJAD
    integer(kind = kint), dimension(*) :: IAJAD
    integer(kind = kint), dimension(*) :: JADORD

    integer(kind = kint) :: IJAD, MAXNZ, MINNZ
    integer(kind = kint) :: I, J, JS, JE, IN, JC
    integer(kind = kint), allocatable :: LEN(:), LENZ(:), JADREORD(:)

    ALLOCATE(LEN(N))
    ALLOCATE(JADREORD(N))
    DO I=1,N
      LEN(I)= hecMAT%indexL(I) - hecMAT%indexL(I-1) &
           &        + hecMAT%indexU(I) - hecMAT%indexU(I-1)
    END DO
    MAXNZ=MAXVAL(LEN(1:N))
    MINNZ=MINVAL(LEN(1:N))
    MJAD =MAXNZ
    ALLOCATE(LENZ(0:MJAD))
    LENZ = 0
    DO I=1,N
      LENZ(LEN(I))=LENZ(LEN(I))+1
    ENDDO
    DO I=MAXNZ-1,MINNZ,-1
      LENZ(I)=LENZ(I)+LENZ(I+1)
    ENDDO
    DO I=1,N
      JADORD(I)=LENZ(LEN(I))
      LENZ(LEN(I))=LENZ(LEN(I))-1
    ENDDO
    DO I=1,N
      JADREORD(JADORD(I))=I
    ENDDO
    DO I=1,N
      LENZ(LEN(JADREORD(I)))=I
    ENDDO
    DO I=MAXNZ-1,1,-1
      LENZ(I)=MAX(LENZ(I+1),LENZ(I))
    ENDDO
    IAJAD(1)=1
    DO I=1,MAXNZ
      IAJAD(I+1)=IAJAD(I)+LENZ(I)
    ENDDO
    LEN=0
    DO I= 1, N
      IJAD=JADORD(I)
      JS= hecMAT%indexL(I-1) + 1
      JE= hecMAT%indexL(I  )
      DO J=JS,JE
        IN  = hecMAT%itemL(J)
        LEN(IJAD)=LEN(IJAD)+1
        JC=IAJAD(LEN(IJAD))+IJAD-1
        AJAD(JC*9-8:JC*9) = hecMAT%AL(9*J-8:9*J)
        JAJAD(JC) = IN
      END DO
    END DO
    DO I= 1, N
      IJAD=JADORD(I)
      JS= hecMAT%indexU(I-1) + 1
      JE= hecMAT%indexU(I  )
      DO J=JS,JE
        IN  = hecMAT%itemU(J)
        LEN(IJAD)=LEN(IJAD)+1
        JC=IAJAD(LEN(IJAD))+IJAD-1
        AJAD(JC*9-8:JC*9) = hecMAT%AU(9*J-8:9*J)
        JAJAD(JC) = IN
      END DO
    END DO
    DEALLOCATE(LEN)
    DEALLOCATE(JADREORD)
    DEALLOCATE(LENZ)
  END SUBROUTINE REPACK

  SUBROUTINE MATJAD(N, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, W1, W2, W3)
    use hecmw_util
    integer(kind=kint) :: N, MJAD
    integer(kind=kint) :: IAJAD(*), JAJAD(*), JADORD(*)
    real(kind=kreal)   :: AJAD(*), X(*), Y(*), W1(*), W2(*), W3(*)

    integer(kind=kint) :: I, K, NZ, IXX
    real(kind=kreal)   :: X1, X2, X3

    DO I=1,N
      W1(I)=0.D0
      W2(I)=0.D0
      W3(I)=0.D0
    ENDDO

    DO NZ=1,MJAD
      DO K=IAJAD(NZ),IAJAD(NZ+1)-1
        X1=X(JAJAD(K)*3-2)
        X2=X(JAJAD(K)*3-1)
        X3=X(JAJAD(K)*3  )
        IXX = K-IAJAD(NZ)+1
        W1(IXX)=W1(IXX) + AJAD(K*9-8)*X1 + AJAD(K*9-7)*X2 + AJAD(K*9-6)*X3
        W2(IXX)=W2(IXX) + AJAD(K*9-5)*X1 + AJAD(K*9-4)*X2 + AJAD(K*9-3)*X3
        W3(IXX)=W3(IXX) + AJAD(K*9-2)*X1 + AJAD(K*9-1)*X2 + AJAD(K*9-0)*X3
      ENDDO
    ENDDO

    DO I=1,N
      Y(3*I-2)=Y(3*I-2)+W1(JADORD(I))
      Y(3*I-1)=Y(3*I-1)+W2(JADORD(I))
      Y(3*I  )=Y(3*I  )+W3(JADORD(I))
    ENDDO
  END SUBROUTINE MATJAD

end module hecmw_JAD_TYPE
