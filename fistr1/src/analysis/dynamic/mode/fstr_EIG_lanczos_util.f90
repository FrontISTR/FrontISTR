!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Lanczos iteration calculation
module m_fstr_EIG_lanczos_util
  contains

!> Initialize Lanczos iterations
  subroutine lanczos_set_initial_value(hecMESH, hecMAT, fstrEIG, EVEC, WK, LVECP, q0, q1, beta, maxiter)
    USE m_fstr
    USE hecmw_util
    implicit none
    type(fstr_eigen) :: fstrEIG
    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF



    REAL(kind=kreal) :: EVEC(:), WK(:,:)

    REAL(kind=kreal) :: LVECP(:), beta(maxiter), chk
    REAL(kind=kreal), POINTER :: q0(:), q1(:), mass(:), filter(:)
    INTEGER(kind=kint) :: GTOT, IRANK, numnp, iov, IRTN, NNN, NN, NTOT, maxiter, i, ierror, j
    INTEGER(kind=kint) :: IXVEC(0:NPROCS-1), IDISP(0:NPROCS-1)
    TYPE (hecmwST_local_mesh) :: hecMESH
    TYPE (hecmwST_matrix    ) :: hecMAT



    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    mass   => fstrEIG%mass
    filter => fstrEIG%filter

    IRANK = myrank
    NTOT = NNDOF
    NN=NTOT


    CALL URAND1(NN,EVEC,IRTN)

     do i = 1,NNDOF
       EVEC(i) = EVEC(i)*filter(i)
     end do

    CALL VECPRO1(chk, EVEC, EVEC, NN)
      CALL hecmw_allreduce_R1(hecMESH,chk,hecmw_sum)

    EVEC(:) = EVEC(:)/sqrt(chk)

    do J=1,NN
      q0(j) = 0.0D0
    enddo

    CALL MATPRO(WK,mass,EVEC,NN,1)

    CALL VECPRO(beta(1),EVEC,WK,NN,1)
      CALL hecmw_allreduce_R(hecMESH,beta(1),1,hecmw_sum)

    beta(1) = SQRT(beta(1))

    IF( beta(1).EQ.0.0D0 ) THEN
      CALL hecmw_finalize()
      STOP "EL1 Self-orthogonal r0!: file Lanczos.f"
    ENDIF

    do J=1,NN
      q1(j) = EVEC(J)/beta(1)
    enddo

    CALL MATPRO(LVECP,mass,q1,NN,1)

  end subroutine lanczos_set_initial_value



!> Sort eigenvalues
      SUBROUTINE EVSORT(EIG,NEW,NEIG)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION EIG(NEIG),NEW(NEIG)

      DO 10 I=1,NEIG
        NEW(I)=I
   10 CONTINUE

      NM=NEIG-1
      DO 20 I=1,NM
        MINLOC=I
        EMIN=ABS(EIG(NEW(I)))
        IP=I+1
        DO 30 J=IP,NEIG
          IF(ABS(EIG(NEW(J))).LT.EMIN) THEN
            MINLOC=J
            EMIN=ABS(EIG(NEW(J)))
          END IF
   30   CONTINUE
        IBAF=NEW(I)
        NEW(I)=NEW(MINLOC)
        NEW(MINLOC)=IBAF
   20 CONTINUE

      RETURN
      END SUBROUTINE EVSORT

!> Scalar product of two vectors
      SUBROUTINE VECPRO(Z,X,Y,NN,NEIG)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION Z(NEIG),X(NN,NEIG),Y(NN,NEIG)
      DO 10 I=1,NEIG
        S=0.0D0
        DO 20 J=1,NN
          S=S+X(J,I)*Y(J,I)
   20   CONTINUE
          Z(I)=S
   10 CONTINUE
      RETURN
      END SUBROUTINE VECPRO

!> Scalar product of two vectors
      SUBROUTINE VECPRO1(Z,X,Y,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)

      DIMENSION X(NN),Y(NN)
        Z=0.0D0
        DO 20 J=1,NN
          Z=Z+X(J)*Y(J)
   20   CONTINUE

      RETURN
      END SUBROUTINE VECPRO1

!> Copy vector
      SUBROUTINE DUPL(X,Y,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION X(NN),Y(NN)
        DO 20 J=1,NN
          X(J) = Y(J)
   20   CONTINUE
      RETURN
      END SUBROUTINE DUPL

!> Scalar shift vector
      SUBROUTINE SCSHFT(X,Y,A,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION X(NN),Y(NN)
        DO 20 J=1,NN
          X(J) = X(J) - A*Y(J)
   20   CONTINUE
      RETURN
      END SUBROUTINE SCSHFT

!> Product of diagonal matrix and vector
      SUBROUTINE MATPRO(Y,A,X,MM,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION Y(MM,NN),A(MM),X(MM,NN)
      DO 10 I=1,NN
        DO 20 J=1,MM
          Y(J,I)=A(J)*X(J,I)
   20   CONTINUE
   10 CONTINUE
      RETURN
      END SUBROUTINE MATPRO

!> Update Lanczos vectors
      SUBROUTINE UPLCZ(X,Y,U,V,A,NN)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION X(NN), Y(NN), U(NN), V(NN)
        DO 20 J=1,NN
          X(J) = U(J)/A
          Y(J) = V(J)/A
   20   CONTINUE
      RETURN
      END SUBROUTINE UPLCZ

!> Random number generator
      SUBROUTINE URAND1(N,X,IR)
      use hecmw
      REAL(kind=kreal) X(N), INVM
      PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
      PARAMETER (INVM = 1.0D0 / M)
      INTEGER(kind=kint) IR

      DO 10 I = 1, N
        IR = MOD( LAMBDA * IR + MU, M)
         X(I) = IR * INVM
   10 CONTINUE
      RETURN
      END SUBROUTINE URAND1

!> Random number generator
      SUBROUTINE URAND0(N,IR)
      use hecmw
      real(kind=kreal) INVM
      PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
      PARAMETER (INVM = 1.0D0 / M)

      DO 10 I = 1, N
        IR = MOD( LAMBDA * IR + MU, M)
   10 CONTINUE
      RETURN
      END SUBROUTINE URAND0


end module m_fstr_EIG_lanczos_util

