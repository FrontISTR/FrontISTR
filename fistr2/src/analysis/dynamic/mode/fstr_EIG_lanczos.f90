!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Eigen Analysis                                    !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> Lanczos iteration calculation
module m_fstr_EIG_lanczos

USE hecmw

contains

!C=====================================================================!
!                       Description                                    !
!C=====================================================================!
!> Initialize Lanczos iterations
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE SETIVL( GMASS,EVEC,EFILT,WK,LVECP,lvecq,BTA, &
     &                         NTOT,NEIG,ISHF,hecMESH,hecMAT,NDOF,GTOT )
!C---------------------------------------------------------------------*
      USE m_fstr
      use lczparm
!C
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION GMASS(NTOT),EVEC(NTOT),EFILT(NTOT),WK(NTOT,1)

      TYPE(lczvec) :: lvecq(0:lvecq_size)
!C
!C --- Lanczos vector & coefficient ---
      REAL(kind=kreal) :: LVECP(NTOT), BTA(NEIG)
      REAL(kind=kreal), POINTER, DIMENSION(:) :: xvec
!C
      INTEGER(kind=kint) GTOT, IRANK, NDOF, numnp, iov, IRTN, NNN, NN
      INTEGER(kind=kint) IXVEC(0:NPROCS-1), ISHF(0:NPROCS-1), IDISP(0:NPROCS-1)
!C
      TYPE (hecmwST_local_mesh) :: hecMESH
      TYPE (hecmwST_matrix    ) :: hecMAT
!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if

!C
      IRANK = myrank
!C
!C** DEFINE NN
      NN=NTOT
!C
!C** SET 1 TO ALL DOF OF FISRT VECTOR and 0 to first q vector (q0)
      CALL URAND1(NN,EVEC,IRTN)
      GO TO 10
!C
!*Shift to appropriate seed if in parallel mode
      IRTN = 1
      NNN  = 0
      DO I = 0,IRANK-1
        NNN = NNN + ISHF(I)
      END DO
!C
      IF(IRANK.GT.0) THEN
        call urand0(NNN,IRTN)
      ENDIF
!C
!C*Then call URAND again from shift
      ALLOCATE(xvec(gtot),STAT=ierror)
      IF(ierror.NE.0)  STOP "Allocation error, SETIVL"
!C
      IF(IRANK.EQ.0) THEN
        call urand1(GTOT,XVEC,IRTN)
      ENDIF
!C
      IDISP(0) = 0
      DO I = 1,NPROCS-1 
        IDISP(I) = IDISP(I-1)+ISHF(I-1)
      ENDDO
!C
      if (.not. ds) then ! in case of Direct Solver prevent MPI.
   !   CALL HECMW_scatterv_DP( XVEC,ISHF(0),IDISP(0), &
!     &                   EVEC,ISHF(IRANK), &
 !    &                   0,hecMESH%MPI_COMM )
      end if
!C      CALL MPI_scatterv( XVEC,ISHF(0),IDISP(0),MPI_DOUBLE_PRECISION,
!C     &                   EVEC(1),ISHF(IRANK),MPI_DOUBLE_PRECISION,
!C     &                   0,hecMESH%MPI_COMM,IERROR )
!C
      if( associated(xvec) ) DEALLOCATE(xvec)
!C
!C*Update EVEC on boundaries
      numnp = NN/NDOF
      if (.not. ds) then ! in case of Direct Solver prevent MPI.
      IF    ( NDOF.eq.3 ) THEN
    !     CALL hecmw_update_m_R( hecMESH,EVEC(1),numnp,NDOF )
      ELSEIF( NDOF.eq.2 ) THEN
    !     CALL hecmw_update_2_R( hecMESH,EVEC(1),numnp )
      ELSEIF( NDOF.EQ.6 ) THEN
     !    CALL hecmw_update_m_R( hecMESH,EVEC(1),numnp,NDOF )
      ENDIF
      end if

 10   CONTINUE
!SPC Node
       do i = 1,NTOT
         EVEC(i) = EVEC(i)*EFILT(i)
       end do
!C
!*Dot product only over nonoverlapping vectors
!C
      CALL VECPRO1(chk,EVEC(1),EVEC(1),ISHF(IRANK))
      if (.not. ds) then ! in case of Direct Solver prevent MPI.
    !    CALL hecmw_allreduce_R1(hecMESH,chk,hecmw_sum)
      end if
      EVEC(:) = EVEC(:)/sqrt(chk)
!C
      do J=1,NN
        lvecq(0)%q(j) = 0.0D0
      enddo
!C
!C** {WK}={GMASS}T{EVEC}
      CALL MATPRO(WK,GMASS,EVEC,NN,1)
!C*
!C* BTA(1)={EV}T[GM]{X}={EV}T{WK}
      CALL VECPRO(BTA(1),EVEC(1),WK(1,1),ISHF(IRANK),1)
      if (.not. ds) then ! in case of Direct Solver prevent MPI.
    !    CALL hecmw_allreduce_R(hecMESH,BTA,1,hecmw_sum)
      end if
!C
!Calculate the first beta value
      BTA(1) = SQRT(BTA(1))
!C
!Calculate q1
      IF( BTA(1).EQ.0.0D0 ) THEN
     !   CALL hecmw_finalize()
        STOP "EL1 Self-orthogonal r0!: file Lanczos.f"
      ENDIF
!C
      do J=1,NN
        lvecq(1)%q(j) = EVEC(J)/BTA(1)
      enddo
!C
!Calculate p1
      CALL MATPRO(LVECP,GMASS,lvecq(1)%q,NN,1)
!C
      RETURN
      END SUBROUTINE SETIVL
!C=====================================================================!
!                       Description                                    !
!C=====================================================================!
!> Sort eigenvalues
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE EVSORT(EIG,NEW,NEIG) 
!C---------------------------------------------------------------------*
!C*
!C* REORDER EIGEN VALUE
!C*
      use hecmw
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION EIG(NEIG),NEW(NEIG)

      DO 10 I=1,NEIG
        NEW(I)=I
   10 CONTINUE
!C
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
!C
      RETURN
      END SUBROUTINE EVSORT
!C=====================================================================!
!                       Description                                    !
!C=====================================================================!
!> Output eigenvalues and vectors
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE EGLIST( hecMESH,hecMAT,IOUT )
!C---------------------------------------------------------------------*
!C*
!C* DISPLAY RESULTS OF EIGEN VALUE ANALYSIS
!C*                         
      use lczparm
      use lczeigen
!C
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
!C
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
!C
!*For parallel part reduction
      INTEGER(kind=kint) nglobal,istt,ied, GID,gmyrank,groot
      INTEGER(kind=kint) groupcount, GROUP, XDIFF, hecGROUP, IOUT
      INTEGER(kind=kint), POINTER :: istarray(:,:), grouping(:), gmem(:) 
      INTEGER(kind=kint), POINTER :: counts(:),disps(:)
      REAL(kind=kreal), POINTER :: xevec(:),xsend(:)
!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if
!C
!C*--- Create local communication groups for contiguous eigenvector ---*
      ALLOCATE( istarray(2,nprocs) )
      ALLOCATE( grouping(0:nprocs-1) )
      ALLOCATE( gmem(nprocs) )
!C
      DO i = 1,nprocs
        istarray(:,i) = 0
        grouping(i-1) = 0
      ENDDO
!C
      istt = hecMESH%global_node_ID(1)
      ied  = hecMESH%global_NODE_ID(numn)
      istarray(1,myrank+1) = istt
      istarray(2,myrank+1) = ied
      if (.not. ds) then
     !   CALL hecmw_allreduce_I( hecMESH,istarray,2*nprocs,hecmw_sum )
      end if

!C
      nglobal = numn
      if (.not. ds) then
    !    CALL hecmw_allreduce_I1( hecMESH,nglobal,hecmw_sum )
      end if
!C
      IF( myrank.EQ.0 ) ALLOCATE( xevec(nglobal*NDOF) )
      ALLOCATE( xsend(Gntotal) )
!C
      GID = 0
      DO i = 1,nprocs
        IF(istt.EQ.istarray(1,i)) GO TO 10
        GID = GID + 1
      ENDDO
   10 CONTINUE
!C
      GROUP = GID
      grouping(myrank) = GID
      if (.not. ds) then
    !    CALL hecmw_allreduce_I( hecMESH,grouping(0:),nprocs,hecmw_sum )
      end if
      groupcount = 1
      j = grouping(0)
      gmem(1) = j
!C
      DO i = 1, nprocs-1
        IF( j.NE.grouping(i) ) THEN
          j = grouping(i)
          groupcount = groupcount+1
          gmem(groupcount) = j
        ENDIF
      ENDDO
!C
      ALLOCATE( counts(groupcount) )
      ALLOCATE( disps(groupcount) )
!C
      DO i = 1, groupcount
        counts(i) = 0
        disps(i)  = 0
      ENDDO
!C
      disps(1) = 0
      DO i = 1,groupcount
        counts(i) = my_ntotal(gmem(i))
        if( i.lt.groupcount ) then
          disps(i+1) = disps(i) + counts(i)
        endif
      ENDDO
!C
!C comment out by imai 2005/08/29
!C      CALL MPI_comm_group ( hecMESH%MPI_COMM,hecGROUP,ierror )
!C      CALL MPI_group_incl ( hecGROUP,groupcount,gmem,GROUP,ierror )
!C      CALL MPI_comm_create( hecMESH%MPI_COMM,GROUP,XDIFF,ierror )
!C 
      PI = 4.0*ATAN(1.0) 
!C
!C*EIGEN VALUE SORTING 
      CALL EVSORT(EVAL,NEW,LTRIAL)
      IF(myrank==0) THEN
        WRITE(IOUT,*) ''
        WRITE(IOUT,'(''********************************'')')
        WRITE(IOUT,'(''*RESULT OF EIGEN VALUE ANALYSIS*'')')
        WRITE(IOUT,'(''********************************'')')
        WRITE(IOUT,'('' '')')
        WRITE(IOUT,*) 'NUMBER OF ITERATIONS =',LTRIAL
        WRITE(IOUT,'('' '')')
        WRITE(IOUT,'(''   NO.   EIGENVALUE    ANGL.FREQUENCY'',&
     &                 '' FREQUENCY(HZ)'')')
        WRITE(IOUT,'(''   ---   ----------    --------------'',&
     &                 '' -------------'')')
        kcount = 0
        DO 40 i=1,LTRIAL
          II=NEW(I)
          if( modal(ii).eq.1 ) then
            kcount = kcount + 1
            EEE=EVAL(II)
            IF(EEE.LT.0.0) EEE=0.0
            WWW=DSQRT(EEE)
            FFF=WWW*0.5/PI
            WRITE(IOUT,'(I5,5E15.6)') kcount,EEE,WWW,FFF
            if( kcount.EQ.NGET ) go to 41
          endif
   40   CONTINUE
   41   continue
        WRITE(IOUT,*)
      ENDIF
      RETURN
      END SUBROUTINE EGLIST
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Scalar product of two vectors
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE VECPRO(Z,X,Y,NN,NEIG)
!C---------------------------------------------------------------------*
!C*
!C* PRODUCT OF VECTOR {X} ,{Y}: {Z}={X}T{Y}
!C*                         
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
!C
      DIMENSION Z(NEIG),X(NN,NEIG),Y(NN,NEIG)
      DO 10 I=1,NEIG
        S=0.0D0
        DO 20 J=1,NN
          S=S+X(J,I)*Y(J,I)
   20   CONTINUE
          Z(I)=S
   10 CONTINUE
!C
      RETURN
      END SUBROUTINE VECPRO
!> Scalar product of two vectors
!C---------------------------------------------------------------------*
      SUBROUTINE VECPRO1(Z,X,Y,NN)
!C---------------------------------------------------------------------*
!C*
!C* PRODUCT OF VECTOR {X} ,{Y}: {Z}={X}T{Y}
!C*                         
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
!C
      DIMENSION X(NN),Y(NN)
        Z=0.0D0
        DO 20 J=1,NN
          Z=Z+X(J)*Y(J)
   20   CONTINUE
!C
      RETURN
      END SUBROUTINE VECPRO1
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Copy vector 
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE DUPL(X,Y,NN)
!C---------------------------------------------------------------------*
!C*
!C* DUPLICATE VECTORS {Y} INTO {X}: {X} = {Y}
!C*                         
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION X(NN),Y(NN)
        DO 20 J=1,NN
          X(J) = Y(J)
   20   CONTINUE
      RETURN
      END SUBROUTINE DUPL
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Scalar shift vector
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE SCSHFT(X,Y,A,NN)
!C---------------------------------------------------------------------*
!C*
!C* SCALAR SHIFT {X} BY A*{Y}: {X} := {X} - A*{Y}
!C*                           
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION X(NN),Y(NN)
        DO 20 J=1,NN
          X(J) = X(J) - A*Y(J)
   20   CONTINUE
      RETURN
      END SUBROUTINE SCSHFT
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Product of diagonal matrix and vector
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE MATPRO(Y,A,X,MM,NN)
!C---------------------------------------------------------------------*
!C*
!C*  MATRIX PRODUCT [Y]=[A][X]
!C*                         
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION Y(MM,NN),A(MM),X(MM,NN)
      DO 10 I=1,NN
        DO 20 J=1,MM
          Y(J,I)=A(J)*X(J,I)
   20   CONTINUE
   10 CONTINUE
      RETURN
      END SUBROUTINE MATPRO
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Update Lanczos vectors
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE UPLCZ(X,Y,U,V,A,NN)
!C---------------------------------------------------------------------*
!C*
!C* DUPLICATE VECTORS {Y} INTO {X}: {X} = {Y}
!C*                         
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION X(NN), Y(NN), U(NN), V(NN) 
        DO 20 J=1,NN
          X(J) = U(J)/A
          Y(J) = V(J)/A
   20   CONTINUE
      RETURN
      END SUBROUTINE UPLCZ
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Random number generator
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE URAND1(N,X,IR)
!C---------------------------------------------------------------------*
      REAL(kind=kreal) X(N), INVM
      PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
      PARAMETER (INVM = 1.0D0 / M)
      INTEGER(kind=kint) IR
!C
      DO 10 I = 1, N
        IR = MOD( LAMBDA * IR + MU, M)
         X(I) = IR * INVM
   10 CONTINUE
      RETURN
      END SUBROUTINE URAND1
!> Random number generator
!C*--------------------------------------------------------------------*
      SUBROUTINE URAND0(N,IR)
!C*--------------------------------------------------------------------*
      real(kind=kreal) INVM
      PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
      PARAMETER (INVM = 1.0D0 / M)
!C*
      DO 10 I = 1, N
        IR = MOD( LAMBDA * IR + MU, M)
   10 CONTINUE
      RETURN
      END SUBROUTINE URAND0
!C=====================================================================!
!C                      Description                                    !
!C=====================================================================!
!> Eigenvector regularization
!C======================================================================
!C---------------------------------------------------------------------*
      SUBROUTINE REGVEC(EVEC,GMASS,XMODE,NTOT,NEIG,NORMAL)
!C---------------------------------------------------------------------*
!C*
!C* EIGEN VECTOR REGULARIZATION (MASS NORMAL)
!C*
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      DIMENSION EVEC(NTOT,NEIG),GMASS(NTOT),XMODE(2,NEIG)
!C
      IF( NORMAL.EQ.0 ) THEN
        DO 100 I = 1, NEIG
          SCALE = 0.0
          DO 200 J = 1, NTOT
            SCALE = SCALE + GMASS(J)*EVEC(J,I)**2
  200     CONTINUE
!C
          DO 300 J = 1, NTOT
            EVEC(J,I) = EVEC(J,I)/SQRT(SCALE)
  300     CONTINUE
          XMODE(1,I) = XMODE(1,I)/SCALE
          XMODE(2,I) = XMODE(2,I)/SCALE
  100   CONTINUE
!C
      ELSE IF( NORMAL.EQ.1 ) THEN
!C**
!C** MASS NORMALIZATION
!C***
        DO 1000 I = 1, NEIG
          EMAX = 0.0
          DO 1100 J = 1, NTOT
            IF(ABS( EVEC(J,I)).GT.EMAX ) EMAX = ABS(EVEC(J,I))
 1100     CONTINUE
          DO 1200 J = 1, NTOT
            EVEC(J,I) = EVEC(J,I)/EMAX
 1200     CONTINUE
!C***
!C*** modal mass & stiffness
!C***
          XMODE(1,I)=XMODE(1,I)/EMAX**2
          XMODE(2,I)=XMODE(2,I)/EMAX**2
 1000   CONTINUE
      END IF
!C
      RETURN
      END SUBROUTINE REGVEC

end module m_fstr_EIG_lanczos

