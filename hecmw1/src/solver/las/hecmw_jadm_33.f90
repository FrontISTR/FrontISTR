!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Jagged Diagonal Matrix storage for vector processors.
!> Original code was provided by JAMSTEC.

module hecmw_JAD_TYPE_33
  use hecmw_util
  use m_hecmw_comm_f
  implicit none

  private

  public :: hecmw_JAD_INIT_33
  public :: hecmw_JAD_FINALIZE_33
  public :: hecmw_JAD_IS_INITIALIZED_33
  public :: hecmw_JAD_MATVEC_33

  !C---------------------- AU&AL
  real(kind=kreal), allocatable      :: AJAD(:)
  integer(kind=kint), allocatable    :: JAJAD(:)
  integer(kind=kint), allocatable    :: JADORD(:)
  integer(kind=kint), allocatable    :: IAJAD(:)
  integer(kind=kint) :: MJAD
  real(kind=kreal), allocatable  :: WP1(:), WP2(:), WP3(:)
  integer(kind=kint) :: INITIALIZED = 0

contains

  subroutine hecmw_JAD_INIT_33(hecMAT)
    type(hecmwST_matrix) :: hecMAT
    allocate(WP1(hecMAT%NP), WP2(hecMAT%NP), WP3(hecMAT%NP))
    allocate(AJAD((hecMAT%NPL+hecMAT%NPU)*9))
    allocate(JAJAD(hecMAT%NPL+hecMAT%NPU))
    allocate(JADORD(hecMAT%NP))
    allocate(IAJAD(hecMAT%NP+1))
    call REPACK(hecMAT%N, hecMAT, MJAD, AJAD, JAJAD, IAJAD, JADORD)
    INITIALIZED = 1
  end subroutine hecmw_JAD_INIT_33

  subroutine hecmw_JAD_FINALIZE_33()
    deallocate(AJAD)
    deallocate(JAJAD)
    deallocate(JADORD)
    deallocate(IAJAD)
    deallocate(WP1,WP2,WP3)
    INITIALIZED = 0
  end subroutine hecmw_JAD_FINALIZE_33

  function hecmw_JAD_IS_INITIALIZED_33()
    integer(kind=kint) :: hecmw_JAD_IS_INITIALIZED_33
    hecmw_JAD_IS_INITIALIZED_33 = INITIALIZED
  end function hecmw_JAD_IS_INITIALIZED_33

  subroutine hecmw_JAD_MATVEC_33(hecMESH, hecMAT, X, Y, COMMtime)
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
    call hecmw_update_R (hecMESH, X, hecMAT%NP, 3)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    D => hecMAT%D

    !$OMP PARALLEL PRIVATE(i)
    !$OMP DO
    do i= 1, hecMAT%N
      X1= X(3*i-2)
      X2= X(3*i-1)
      X3= X(3*i  )
      Y(3*i -2)= D(9*i-8)*X1 + D(9*i-7)*X2 + D(9*i-6)*X3
      Y(3*i -1)= D(9*i-5)*X1 + D(9*i-4)*X2 + D(9*i-3)*X3
      Y(3*i   )= D(9*i-2)*X1 + D(9*i-1)*X2 + D(9*i  )*X3
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call MATJAD(hecMAT%N, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, WP1, WP2, WP3)
  end subroutine hecmw_JAD_MATVEC_33

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
    integer(kind = kint) :: I, J, JS, JE, in, JC
    integer(kind = kint), allocatable :: len(:), LENZ(:), JADREORD(:)

    allocate(len(N))
    allocate(JADREORD(N))
    do I=1,N
      len(I)= hecMAT%indexL(I) - hecMAT%indexL(I-1) &
        &        + hecMAT%indexU(I) - hecMAT%indexU(I-1)
    end do
    MAXNZ=maxval(len(1:N))
    MINNZ=minval(len(1:N))
    MJAD =MAXNZ
    allocate(LENZ(0:MJAD))
    LENZ = 0
    do I=1,N
      LENZ(len(I))=LENZ(len(I))+1
    enddo
    do I=MAXNZ-1,MINNZ,-1
      LENZ(I)=LENZ(I)+LENZ(I+1)
    enddo
    do I=1,N
      JADORD(I)=LENZ(len(I))
      LENZ(len(I))=LENZ(len(I))-1
    enddo
    do I=1,N
      JADREORD(JADORD(I))=I
    enddo
    do I=1,N
      LENZ(len(JADREORD(I)))=I
    enddo
    do I=MAXNZ-1,1,-1
      LENZ(I)=max(LENZ(I+1),LENZ(I))
    enddo
    IAJAD(1)=1
    do I=1,MAXNZ
      IAJAD(I+1)=IAJAD(I)+LENZ(I)
    enddo
    len=0
    do I= 1, N
      IJAD=JADORD(I)
      JS= hecMAT%indexL(I-1) + 1
      JE= hecMAT%indexL(I  )
      do J=JS,JE
        in  = hecMAT%itemL(J)
        len(IJAD)=len(IJAD)+1
        JC=IAJAD(len(IJAD))+IJAD-1
        AJAD(JC*9-8:JC*9) = hecMAT%AL(9*J-8:9*J)
        JAJAD(JC) = in
      end do
    end do
    do I= 1, N
      IJAD=JADORD(I)
      JS= hecMAT%indexU(I-1) + 1
      JE= hecMAT%indexU(I  )
      do J=JS,JE
        in  = hecMAT%itemU(J)
        len(IJAD)=len(IJAD)+1
        JC=IAJAD(len(IJAD))+IJAD-1
        AJAD(JC*9-8:JC*9) = hecMAT%AU(9*J-8:9*J)
        JAJAD(JC) = in
      end do
    end do
    deallocate(len)
    deallocate(JADREORD)
    deallocate(LENZ)
  end subroutine REPACK

  subroutine MATJAD(N, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, W1, W2, W3)
    use hecmw_util
    integer(kind=kint) :: N, MJAD
    integer(kind=kint) :: IAJAD(*), JAJAD(*), JADORD(*)
    real(kind=kreal)   :: AJAD(*), X(*), Y(*), W1(*), W2(*), W3(*)

    integer(kind=kint) :: I, K, NZ, IXX
    real(kind=kreal)   :: X1, X2, X3

    !$OMP PARALLEL PRIVATE(I,K,X1,X2,X3,IXX)
    !$OMP DO
    do I=1,N
      W1(I)=0.D0
      W2(I)=0.D0
      W3(I)=0.D0
    enddo
    !$OMP END DO

    do NZ=1,MJAD
      !$OMP DO
      do K=IAJAD(NZ),IAJAD(NZ+1)-1
        X1=X(JAJAD(K)*3-2)
        X2=X(JAJAD(K)*3-1)
        X3=X(JAJAD(K)*3  )
        IXX = K-IAJAD(NZ)+1
        W1(IXX)=W1(IXX) + AJAD(K*9-8)*X1 + AJAD(K*9-7)*X2 + AJAD(K*9-6)*X3
        W2(IXX)=W2(IXX) + AJAD(K*9-5)*X1 + AJAD(K*9-4)*X2 + AJAD(K*9-3)*X3
        W3(IXX)=W3(IXX) + AJAD(K*9-2)*X1 + AJAD(K*9-1)*X2 + AJAD(K*9-0)*X3
      enddo
      !$OMP END DO
    enddo

    !$OMP DO
    do I=1,N
      Y(3*I-2)=Y(3*I-2)+W1(JADORD(I))
      Y(3*I-1)=Y(3*I-1)+W2(JADORD(I))
      Y(3*I  )=Y(3*I  )+W3(JADORD(I))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine MATJAD

end module hecmw_JAD_TYPE_33
