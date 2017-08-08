!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Lanczos iteration calculation
module m_fstr_EIG_lanczos_util
  contains

!> Initialize Lanczos iterations
  SUBROUTINE SETIVL( GMASS, EVEC, EFILT, WK, LVECP, q0, q1, BTA, &
                   & NTOT, NEIG, hecMESH, hecMAT, NDOF, GTOT )
    USE m_fstr
    USE hecmw_util
    implicit none
    REAL(kind=kreal) :: GMASS(NTOT), EVEC(NTOT), EFILT(NTOT), WK(NTOT,1)

    REAL(kind=kreal) :: LVECP(NTOT), BTA(NEIG), chk
    REAL(kind=kreal), POINTER :: xvec(:), q0(:), q1(:)
    INTEGER(kind=kint) :: GTOT, IRANK, NDOF, numnp, iov, IRTN, NNN, NN, NTOT, NEIG, i, ierror, j
    INTEGER(kind=kint) :: IXVEC(0:NPROCS-1), IDISP(0:NPROCS-1)
    TYPE (hecmwST_local_mesh) :: hecMESH
    TYPE (hecmwST_matrix    ) :: hecMAT

    IRANK = myrank
    NN=NTOT

    CALL URAND1(NN,EVEC,IRTN)

     do i = 1,NTOT
       EVEC(i) = EVEC(i)*EFILT(i)
     end do

    CALL VECPRO1(chk,EVEC(1),EVEC(1), NN)
      CALL hecmw_allreduce_R1(hecMESH,chk,hecmw_sum)

    EVEC(:) = EVEC(:)/sqrt(chk)

    do J=1,NN
      q0(j) = 0.0D0
    enddo

    CALL MATPRO(WK,GMASS,EVEC,NN,1)

    CALL VECPRO(BTA(1),EVEC(1),WK(1,1), NN,1)
      CALL hecmw_allreduce_R(hecMESH,BTA,1,hecmw_sum)

    BTA(1) = SQRT(BTA(1))

    IF( BTA(1).EQ.0.0D0 ) THEN
      CALL hecmw_finalize()
      STOP "EL1 Self-orthogonal r0!: file Lanczos.f"
    ENDIF

    do J=1,NN
      q1(j) = EVEC(J)/BTA(1)
    enddo

    CALL MATPRO(LVECP,GMASS,q1,NN,1)
    RETURN
  END SUBROUTINE SETIVL

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

!> Eigenvector regularization
      SUBROUTINE REGVEC(EVEC,GMASS,XMODE,NTOT,NEIG,NORMAL)
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION EVEC(NTOT,NEIG),GMASS(NTOT),XMODE(2,NEIG)

      IF( NORMAL.EQ.0 ) THEN
        DO 100 I = 1, NEIG
          SCALE = 0.0
          DO 200 J = 1, NTOT
            SCALE = SCALE + GMASS(J)*EVEC(J,I)**2
  200     CONTINUE

          DO 300 J = 1, NTOT
            EVEC(J,I) = EVEC(J,I)/SQRT(SCALE)
  300     CONTINUE
          XMODE(1,I) = XMODE(1,I)/SCALE
          XMODE(2,I) = XMODE(2,I)/SCALE
  100   CONTINUE

      ELSE IF( NORMAL.EQ.1 ) THEN

        DO 1000 I = 1, NEIG
          EMAX = 0.0
          DO 1100 J = 1, NTOT
            IF(ABS( EVEC(J,I)).GT.EMAX ) EMAX = ABS(EVEC(J,I))
 1100     CONTINUE
          DO 1200 J = 1, NTOT
            EVEC(J,I) = EVEC(J,I)/EMAX
 1200     CONTINUE

          XMODE(1,I)=XMODE(1,I)/EMAX**2
          XMODE(2,I)=XMODE(2,I)/EMAX**2
 1000   CONTINUE
      END IF

      RETURN
      END SUBROUTINE REGVEC

end module m_fstr_EIG_lanczos_util

