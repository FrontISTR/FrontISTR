!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Jagged Diagonal Matrix storage for vector processors.
!> Original code was provided by JAMSTEC.

module hecmw_JAD_TYPE_44
  use hecmw_util
  use m_hecmw_comm_f
  implicit none

  private

  public :: hecmw_JAD_INIT_44
  public :: hecmw_JAD_FINALIZE_44
  public :: hecmw_JAD_IS_INITIALIZED_44
  public :: hecmw_JAD_MATVEC_44

  !C---------------------- AU&AL
  real(kind=kreal), allocatable      :: AJAD(:)
  integer(kind=kint), allocatable    :: JAJAD(:)
  integer(kind=kint), allocatable    :: JADORD(:)
  integer(kind=kint), allocatable    :: IAJAD(:)
  integer(kind=kint) :: MJAD
  real(kind=kreal), allocatable  :: WP1(:), WP2(:), WP3(:), WP4(:)
  integer(kind=kint) :: INITIALIZED = 0

contains

  subroutine hecmw_JAD_INIT_44(hecMAT)
    type(hecmwST_matrix) :: hecMAT
    allocate(WP1(hecMAT%NP), WP2(hecMAT%NP), WP3(hecMAT%NP), WP4(hecMAT%NP))
    allocate(AJAD((hecMAT%NPL+hecMAT%NPU)*16))
    allocate(JAJAD(hecMAT%NPL+hecMAT%NPU))
    allocate(JADORD(hecMAT%NP))
    allocate(IAJAD(hecMAT%NP+1))
    call REPACK(hecMAT%N, hecMAT, MJAD, AJAD, JAJAD, IAJAD, JADORD)
    INITIALIZED = 1
  end subroutine hecmw_JAD_INIT_44

  subroutine hecmw_JAD_FINALIZE_44()
    deallocate(AJAD)
    deallocate(JAJAD)
    deallocate(JADORD)
    deallocate(IAJAD)
    deallocate(WP1,WP2,WP3,WP4)
    INITIALIZED = 0
  end subroutine hecmw_JAD_FINALIZE_44

  function hecmw_JAD_IS_INITIALIZED_44()
    integer(kind=kint) :: hecmw_JAD_IS_INITIALIZED_44
    hecmw_JAD_IS_INITIALIZED_44 = INITIALIZED
  end function hecmw_JAD_IS_INITIALIZED_44

  subroutine hecmw_JAD_MATVEC_44(hecMESH, hecMAT, X, Y, COMMtime)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime
    real(kind=kreal) :: START_TIME, END_TIME
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: i
    real(kind=kreal) :: X1, X2, X3, X4

    START_TIME= HECMW_WTIME()
    call hecmw_update_4_R (hecMESH, X, hecMAT%NP)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    D => hecMAT%D

    !$OMP PARALLEL PRIVATE(I,K,X1,X2,X3,IXX)
    !$OMP DO
    do i= 1, hecMAT%N
      X1= X(4*i-3)
      X2= X(4*i-2)
      X3= X(4*i-1)
      X4= X(4*i  )
      Y(4*i-3)= D(16*i-15)*X1 + D(16*i-14)*X2 + D(16*i-13)*X3 + D(16*i-12)*X4
      Y(4*i-2)= D(16*i-11)*X1 + D(16*i-10)*X2 + D(16*i- 9)*X3 + D(16*i- 8)*X4
      Y(4*i-1)= D(16*i- 7)*X1 + D(16*i- 6)*X2 + D(16*i- 5)*X3 + D(16*i- 4)*X4
      Y(4*i  )= D(16*i- 3)*X1 + D(16*i- 2)*X2 + D(16*i- 1)*X3 + D(16*i- 0)*X4
    enddo
    !$OMP END DO
    call MATJAD(hecMAT%N, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, WP1, WP2, WP3, WP4)
  end subroutine hecmw_JAD_MATVEC_44

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
        AJAD(JC*16-15:JC*16) = hecMAT%AL(16*J-15:16*J)
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
        AJAD(JC*16-15:JC*16) = hecMAT%AU(16*J-15:16*J)
        JAJAD(JC) = in
      end do
    end do
    deallocate(len)
    deallocate(JADREORD)
    deallocate(LENZ)
  end subroutine REPACK

  subroutine MATJAD(N, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, W1, W2, W3, W4)
    use hecmw_util
    integer(kind=kint) :: N, MJAD
    integer(kind=kint) :: IAJAD(*), JAJAD(*), JADORD(*)
    real(kind=kreal)   :: AJAD(*), X(*), Y(*), W1(*), W2(*), W3(*), W4(*)

    integer(kind=kint) :: I, K, NZ, IXX
    real(kind=kreal)   :: X1, X2, X3, X4

    !$OMP PARALLEL PRIVATE(I,K,X1,X2,X3,IXX)
    !$OMP DO
    do I=1,N
      W1(I)=0.D0
      W2(I)=0.D0
      W3(I)=0.D0
      W4(I)=0.D0
    enddo
    !$OMP END DO

    do NZ=1,MJAD
      !$OMP DO
      do K=IAJAD(NZ),IAJAD(NZ+1)-1
        X1=X(JAJAD(K)*4-3)
        X2=X(JAJAD(K)*4-2)
        X3=X(JAJAD(K)*4-1)
        X4=X(JAJAD(K)*4  )
        IXX = K-IAJAD(NZ)+1
        W1(IXX)=W1(IXX) + AJAD(K*16-15)*X1 + AJAD(K*16-14)*X2 + AJAD(K*16-13)*X3 + AJAD(K*16-12)*X4
        W2(IXX)=W2(IXX) + AJAD(K*16-11)*X1 + AJAD(K*16-10)*X2 + AJAD(K*16- 9)*X3 + AJAD(K*16- 8)*X4
        W3(IXX)=W3(IXX) + AJAD(K*16- 7)*X1 + AJAD(K*16- 6)*X2 + AJAD(K*16- 5)*X3 + AJAD(K*16- 4)*X4
        W4(IXX)=W4(IXX) + AJAD(K*16- 3)*X1 + AJAD(K*16- 2)*X2 + AJAD(K*16- 1)*X3 + AJAD(K*16- 0)*X4
      enddo
      !$OMP END DO
    enddo

    !$OMP DO
    do I=1,N
      Y(4*I-3)=Y(4*I-3)+W1(JADORD(I))
      Y(4*I-2)=Y(4*I-2)+W2(JADORD(I))
      Y(4*I-1)=Y(4*I-1)+W3(JADORD(I))
      Y(4*I  )=Y(4*I  )+W4(JADORD(I))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine MATJAD

end module hecmw_JAD_TYPE_44
