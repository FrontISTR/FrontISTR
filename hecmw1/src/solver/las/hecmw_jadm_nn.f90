!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Jagged Diagonal Matrix storage for vector processors.
!> Original code was provided by JAMSTEC.

module hecmw_JAD_TYPE_nn
  use hecmw_util
  use m_hecmw_comm_f
  implicit none

  private

  public :: hecmw_JAD_INIT_nn
  public :: hecmw_JAD_FINALIZE_nn
  public :: hecmw_JAD_IS_INITIALIZED_nn
  public :: hecmw_JAD_MATVEC_nn

  !C---------------------- AU&AL
  real(kind=kreal), allocatable      :: AJAD(:)
  integer(kind=kint), allocatable    :: JAJAD(:)
  integer(kind=kint), allocatable    :: JADORD(:)
  integer(kind=kint), allocatable    :: IAJAD(:)
  integer(kind=kint) :: MJAD
  real(kind=kreal), allocatable  :: WP(:,:)
  integer(kind=kint) :: INITIALIZED = 0

contains

  subroutine hecmw_JAD_INIT_nn(hecMAT)
    type(hecmwST_matrix) :: hecMAT
    allocate(WP(hecMAT%NDOF,hecMAT%NP))
    allocate(AJAD((hecMAT%NPL+hecMAT%NPU)*hecMAT%NDOF*hecMAT%NDOF))
    allocate(JAJAD(hecMAT%NPL+hecMAT%NPU))
    allocate(JADORD(hecMAT%NP))
    allocate(IAJAD(hecMAT%NP+1))
    call REPACK(hecMAT%N, hecMAT, MJAD, AJAD, JAJAD, IAJAD, JADORD)
    INITIALIZED = 1
  end subroutine hecmw_JAD_INIT_nn

  subroutine hecmw_JAD_FINALIZE_nn()
    deallocate(AJAD)
    deallocate(JAJAD)
    deallocate(JADORD)
    deallocate(IAJAD)
    deallocate(WP)
    INITIALIZED = 0
  end subroutine hecmw_JAD_FINALIZE_nn

  function hecmw_JAD_IS_INITIALIZED_nn()
    integer(kind=kint) :: hecmw_JAD_IS_INITIALIZED_nn
    hecmw_JAD_IS_INITIALIZED_nn = INITIALIZED
  end function hecmw_JAD_IS_INITIALIZED_nn

  subroutine hecmw_JAD_MATVEC_nn(hecMESH, hecMAT, X, Y, COMMtime)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in), target :: hecMAT
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Y(:)
    real(kind=kreal), intent(inout) :: COMMtime
    real(kind=kreal) :: START_TIME, END_TIME
    real(kind=kreal), pointer :: D(:)
    integer(kind=kint) :: i,idof,jdof,NDOF,NDOF2

    START_TIME= HECMW_WTIME()
    call hecmw_update_R (hecMESH, X, hecMAT%NP,hecMAT%NDOF)
    END_TIME= HECMW_WTIME()
    COMMtime = COMMtime + END_TIME - START_TIME

    D => hecMAT%D
    NDOF = hecMAT%NDOF
    NDOF2 = NDOF*NDOF
    Y=0.0d0

    !$OMP PARALLEL PRIVATE(i, idof, jdof)
    !$OMP DO
    do i= 1, hecMAT%N
      do idof=1,hecMAT%NDOF
        do jdof=1,hecMAT%NDOF
          Y(NDOF*(i-1)+idof) = Y(NDOF*(i-1)+idof) + D(NDOF2*(i-1)+NDOF*(idof-1)+jdof)*X(NDOF*(i-1)+jdof)
        end do
      end do
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call MATJAD(hecMAT%N,hecMAT%NDOF, MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, WP)
  end subroutine hecmw_JAD_MATVEC_nn

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

    integer(kind = kint) :: IJAD, MAXNZ, MINNZ,NDOF,NDOF2
    integer(kind = kint) :: I, J, JS, JE, in, JC
    integer(kind = kint), allocatable :: len(:), LENZ(:), JADREORD(:)
    NDOF = hecMAT%NDOF;NDOF2=NDOF*NDOF
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
        AJAD(NDOF2*(JC-1)+1:NDOF2*(JC)) = hecMAT%AL(NDOF2*(J-1)+1:NDOF2*(J))
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
        AJAD(NDOF2*(JC-1)+1:NDOF2*(JC)) = hecMAT%AU(NDOF2*(J-1)+1:NDOF2*(J))
        JAJAD(JC) = in
      end do
    end do
    deallocate(len)
    deallocate(JADREORD)
    deallocate(LENZ)
  end subroutine REPACK

  subroutine MATJAD(N, NDOF,  MJAD, IAJAD, JAJAD, AJAD, JADORD, X, Y, W)
    use hecmw_util
    integer(kind=kint) :: N,NDOF, MJAD,NDOF2
    integer(kind=kint) :: IAJAD(*), JAJAD(*), JADORD(*)
    real(kind=kreal)   :: AJAD(*), X(*), Y(*),  W(NDOF,N)
    integer(kind=kint) :: I, K, NZ, IXX, idof,jdof
    NDOF2=NDOF*NDOF

    W=0.0d0

    do NZ=1,MJAD
      !$OMP PARALLEL PRIVATE(K,IXX,idof,jdof)
      !$OMP DO
      do K=IAJAD(NZ),IAJAD(NZ+1)-1
        IXX = K-IAJAD(NZ)+1
        do idof = 1, NDOF
          do jdof = 1, NDOF
            W(idof,IXX)=W(idof,IXX)+AJAD(NDOF2*(K-1)+NDOF*(idof-1)+jdof)*X(NDOF*(JAJAD(K)-1)+jdof)
          end do
        end do
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    enddo

    !$OMP PARALLEL PRIVATE(I,idof)
    !$OMP DO
    do I=1,N
      do idof = 1, NDOF
        Y(NDOF*(I-1)+idof)=Y(NDOF*(I-1)+idof)+W(idof,JADORD(I))
      end do
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine MATJAD

end module hecmw_JAD_TYPE_nn
