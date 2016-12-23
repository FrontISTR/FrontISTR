!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

MODULE HECMW_SOLVER_DIRECT
  IMPLICIT NONE

  INTEGER(KIND=4), PRIVATE :: len_colno
  INTEGER(KIND=4), PRIVATE :: nstop
  INTEGER(KIND=4), PRIVATE :: stage
  INTEGER(KIND=4), PRIVATE :: neqns
  INTEGER(KIND=4), PRIVATE :: nttbr
  INTEGER(KIND=4), PRIVATE :: isym
  INTEGER(KIND=4), PRIVATE :: ndeg
  INTEGER(KIND=4), PRIVATE :: irr
  INTEGER(KIND=4), PRIVATE :: len_dsln
  INTEGER(KIND=4), PRIVATE :: len_iv
  INTEGER(KIND=4), PRIVATE :: total

  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: jcol
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: irow
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: zpiv
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: iperm
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: invp
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: parent
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: nch
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: xlnzr
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: colno
  !*Work arrays
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: jcpt
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: jcolno
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: ia
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: ja
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: deg
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: marker
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: rchset
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: nbrhd
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: qsize
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: qlink
  INTEGER(KIND=4), PRIVATE :: nofsub
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: adjncy
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: btree
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: pordr
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: adjncp
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: xleaf
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: leaf
  INTEGER(KIND=4), PRIVATE, DIMENSION(:), POINTER :: indx

  REAL(KIND=8), PRIVATE, DIMENSION(:), POINTER :: val
  REAL(KIND=8), PRIVATE, DIMENSION(:), POINTER :: temp
  REAL(KIND=8), PRIVATE, DIMENSION(:), POINTER :: diag
  REAL(KIND=8), PRIVATE, DIMENSION(:), POINTER :: zln
  REAL(KIND=8), PRIVATE, DIMENSION(:), POINTER :: dsln
  !*Timing
  REAL(KIND=8), PRIVATE, DIMENSION(10) :: tom
  !*Allocation variables
  INTEGER(KIND=4), PRIVATE :: ialoc
  INTEGER(KIND=4), PRIVATE :: raloc
  INTEGER(KIND=4), PRIVATE :: ierror

CONTAINS
  !----------------------------------------------------------------------
  !> @brief HECMW_SOLVE_DIRECT is a program for the matrix solver
  !     and is also available for performance tests
  !
  !     to solve Ax=b for x, nozero pattern and values must be give.
  !
  !   arrays
  !     jcol     column entrys
  !     irow     row entrys
  !     val      its value
  !     v        interface array used through matini, staijx, nufctx,
  !              nusolx
  !
  !     #coded by t.arakawa of RIST on 040316
  !     Modified G. Prabhakar, RIST, June 7 2004
  !----------------------------------------------------------------------
  SUBROUTINE HECMW_SOLVE_DIRECT(Hecmesh,Hecmat,Ifmsg)
    USE HECMW_UTIL
    USE HECMW_MATRIX_ASS
    USE HECMW_MATRIX_DUMP
    IMPLICIT NONE

    TYPE (HECMWST_LOCAL_MESH), INTENT(IN)::Hecmesh
    INTEGER, INTENT(IN):: Ifmsg
    TYPE (HECMWST_MATRIX), INTENT(OUT)::Hecmat

    INTEGER:: ISEed
    INTEGER:: IXXx
    INTEGER:: LRAtio
    INTEGER:: i98
    INTEGER:: i97
    INTEGER:: ir
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio
    COMMON /QAZ   / ISEed, IXXx

    RMAx = 8.988D+307
    RMIn = 4.941D-300
    EPSm = 2.220D-16
    LRAtio = 2
    ISEed = 1
    ir = 0

    CALL HECMW_MAT_ASS_EQUATION(Hecmesh,Hecmat)
    CALL HECMW_MAT_DUMP(Hecmat,Hecmesh)

    CALL PTIME(t1)

    !*EHM HECMW June 7 2004
    i98 = Hecmat%IARRAY(98)
    IF ( Hecmat%IARRAY(98)==1 ) THEN
       !* Interface to symbolic factorization
       CALL SETIJ(Hecmesh,Hecmat)

       !* Symbolic factorization
       CALL MATINI(ir)
       Hecmat%IARRAY(98) = 0
       WRITE (6,*) "symbolic fct done"
    ENDIF
    CALL PTIME(t2)
    t3 = t2


    i97 = Hecmat%IARRAY(97)
    IF ( Hecmat%IARRAY(97)==1 ) THEN
       !* Interface to numeric factorization
       CALL NUFORM(Hecmesh,Hecmat,ir)
       CALL PTIME(t3)

       !* Numeric factorization
       CALL NUFCT0(ir)
       Hecmat%IARRAY(97) = 0

       !*Memory Details
       WRITE (*,*) '*-----------------------------------*'
       WRITE (*,*) '|   Direct  Solver  Memory  Usage   |'
       WRITE (*,*) '*-----------------------------------*'
       WRITE (*,*) 'INTEGER memory: ', REAL(IALoc*4)/REAL(1048576), 'MB'
       WRITE (*,*) 'REAL*8  memory: ', REAL(RALoc*8)/REAL(1048576), 'MB'
       WRITE (*,*) 'TOTAL   memory: ', REAL((RALoc*2+IALoc)*4)/REAL(1048576), 'MB'
       WRITE (*,*) '*-----------------------------------*'
    ENDIF
    CALL PTIME(t4)

    !* Finalize
    !*  Errors 1
    IF ( i98/=0 .AND. i98/=1 ) THEN
       WRITE (Ifmsg,*) 'ERROR in symb. fact. flag: Should be 1 or 0'
       STOP 'ERROR in symb. fact. flag: Should be 1 or 0'
    ENDIF
    IF ( i97/=0 .AND. i97/=1 ) THEN
       WRITE (Ifmsg,*) 'ERROR in numer. fact. flag: Should be 1 or 0'
       STOP 'ERROR in numer. fact. flag: Should be 1 or 0'
    ENDIF
    IF ( i98==1 .AND. i97==0 ) THEN
       WRITE (Ifmsg,*) 'WARNING: Numeric factorization not performed!'
       STOP 'WARNING: Numeric factorization not performed! Solve will not be performed'
    ENDIF
    !*  Errors 2
    IF ( ir/=0 ) THEN
       !WINDEBUG
       WRITE (Ifmsg,*) 'ERROR in nufct0. ir = ', ir
       STOP
    ENDIF

    TOM(1) = t2-t1
    TOM(2) = t3-t2
    TOM(3) = t4-t3

    !* Solve
    !* Backsubstitute
    CALL NUSOL0(Hecmat%B,ir)
    CALL PTIME(t5)
    !* Errors 4
    IF ( ir/=0 ) THEN
       !WINDEBUG
       WRITE (Ifmsg,*) 'error in nusol0. irr = ', ir
       STOP
    ENDIF
    CALL HECMW_MAT_DUMP_SOLUTION(Hecmat)
  END SUBROUTINE HECMW_SOLVE_DIRECT

  !======================================================================!
  !> @brief ADDR0
  !======================================================================!
  SUBROUTINE ADDR0(Isw,I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ndeg2,Ndeg2l,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: I
    INTEGER, INTENT(IN):: Isw
    INTEGER, INTENT(IN):: J
    INTEGER, INTENT(IN):: Ndeg2
    INTEGER, INTENT(IN):: Ndeg2l
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Invp(*)
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Aij(Ndeg2)
    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg2,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(Ndeg2l,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(Ndeg2,*)

    INTEGER:: i0
    INTEGER:: idbg
    INTEGER:: ii
    INTEGER:: itrans
    INTEGER:: j0
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l

    DATA idbg/0/

    Ir = 0
    ii = Invp(I)
    jj = Invp(J)
    IF ( idbg/=0 ) WRITE (6,*) ii, jj, Aij
    IF ( ii==jj ) THEN
       IF ( Ndeg2==1 ) THEN
          IF ( Isw==0 ) THEN
             Diag(1,ii) = Aij(1)
          ELSE
             Diag(1,ii) = Diag(1,ii) + Aij(1)
          ENDIF
       ELSEIF ( Ndeg2==4 ) THEN
          IF ( Isw==0 ) THEN
             Diag(1,ii) = Aij(1)
             Diag(2,ii) = Aij(2)
             Diag(3,ii) = Aij(4)
          ELSE
             Diag(1,ii) = Diag(1,ii) + Aij(1)
             Diag(2,ii) = Diag(2,ii) + Aij(2)
             Diag(3,ii) = Diag(3,ii) + Aij(4)
          ENDIF
       ENDIF
       RETURN
    ENDIF
    itrans = 0
    IF ( jj>ii ) THEN
       k = jj
       jj = ii
       ii = k
       itrans = 1
    ENDIF
    IF ( jj>=Nstop ) THEN
       i0 = ii-Nstop
       j0 = jj-Nstop + 1
       k = i0*(i0-1)/2 + j0
       IF ( Ndeg2==1 ) THEN
          Dsln(1,k) = Aij(1)
          RETURN
       ELSEIF ( Ndeg2==4 ) THEN
          IF ( itrans==0 ) THEN
             DO l = 1, Ndeg2
                Dsln(l,k) = Aij(l)
             ENDDO
          ELSE
             Dsln(1,k) = Aij(1)
             Dsln(2,k) = Aij(3)
             Dsln(3,k) = Aij(2)
             Dsln(4,k) = Aij(4)
          ENDIF
          RETURN
       ENDIF
    ENDIF
    ks = Xlnzr(ii)
    ke = Xlnzr(ii+1) - 1
    DO k = ks, ke
       IF ( Colno(k)==jj ) THEN
          IF ( Isw==0 ) THEN
             IF ( Ndeg2==1 ) THEN
                Zln(1,k) = Aij(1)
             ELSEIF ( Ndeg2==4 ) THEN
                IF ( itrans==0 ) THEN
                   DO l = 1, Ndeg2
                      Zln(l,k) = Aij(l)
                   ENDDO
                ELSE
                   Zln(1,k) = Aij(1)
                   Zln(2,k) = Aij(3)
                   Zln(3,k) = Aij(2)
                   Zln(4,k) = Aij(4)
                ENDIF
             ENDIF
          ELSEIF ( Ndeg2==1 ) THEN
             Zln(1,k) = Zln(1,k) + Aij(1)
          ELSEIF ( Ndeg2==4 ) THEN
             IF ( itrans==0 ) THEN
                DO l = 1, Ndeg2
                   Zln(l,k) = Zln(l,k) + Aij(l)
                ENDDO
             ELSE
                Zln(1,k) = Zln(1,k) + Aij(1)
                Zln(2,k) = Zln(2,k) + Aij(3)
                Zln(3,k) = Zln(3,k) + Aij(2)
                Zln(4,k) = Zln(4,k) + Aij(4)
             ENDIF
          ENDIF
          RETURN
       ENDIF
    ENDDO
    Ir = 20
  END SUBROUTINE ADDR0

  !======================================================================!
  !> @brief ADDR3
  !======================================================================!
  SUBROUTINE ADDR3(I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: I
    INTEGER, INTENT(IN):: J
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Invp(*)
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Aij(9)
    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Zln(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(9,*)

    INTEGER:: i0
    INTEGER:: idbg
    INTEGER:: ii
    INTEGER:: itrans
    INTEGER:: j0
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: ndeg2
    INTEGER:: ndeg2l
    DATA idbg, ndeg2, ndeg2l/0, 9, 6/

    Ir = 0
    ii = Invp(I)
    jj = Invp(J)
    IF ( idbg/=0 ) WRITE (6,*) ii, jj, Aij
    IF ( ii==jj ) THEN
       Diag(1,ii) = Aij(1)
       Diag(2,ii) = Aij(2)
       Diag(3,ii) = Aij(5)
       Diag(4,ii) = Aij(3)
       Diag(5,ii) = Aij(6)
       Diag(6,ii) = Aij(9)
       RETURN
    ENDIF
    itrans = 0
    IF ( jj>ii ) THEN
       k = jj
       jj = ii
       ii = k
       itrans = 1
    ENDIF
    IF ( jj>=Nstop ) THEN
       i0 = ii - Nstop
       j0 = jj - Nstop + 1
       k = i0*(i0-1)/2 + j0
       IF ( itrans==0 ) THEN
          DO l = 1, ndeg2
             Dsln(l,k) = Aij(l)
          ENDDO
       ELSE
          Dsln(1,k) = Aij(1)
          Dsln(2,k) = Aij(4)
          Dsln(3,k) = Aij(7)
          Dsln(4,k) = Aij(2)
          Dsln(5,k) = Aij(5)
          Dsln(6,k) = Aij(8)
          Dsln(7,k) = Aij(3)
          Dsln(8,k) = Aij(6)
          Dsln(9,k) = Aij(9)
       ENDIF
       RETURN
    ENDIF
    ks = Xlnzr(ii)
    ke = Xlnzr(ii+1) - 1
    DO k = ks, ke
       IF ( Colno(k)==jj ) THEN
          IF ( itrans==0 ) THEN
             DO l = 1, ndeg2
                Zln(l,k) = Aij(l)
             ENDDO
          ELSE
             Zln(1,k) = Aij(1)
             Zln(2,k) = Aij(4)
             Zln(3,k) = Aij(7)
             Zln(4,k) = Aij(2)
             Zln(5,k) = Aij(5)
             Zln(6,k) = Aij(8)
             Zln(7,k) = Aij(3)
             Zln(8,k) = Aij(6)
             Zln(9,k) = Aij(9)
          ENDIF
          RETURN
       ENDIF
    ENDDO
    Ir = 20
  END SUBROUTINE ADDR3

  !======================================================================!
  !> @brief ADDRX
  !======================================================================!
  SUBROUTINE ADDRX(I,J,Aij,Invp,Xlnzr,Colno,Diag,Zln,Dsln,Nstop,Ndeg,Ndeg2l,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: I
    INTEGER, INTENT(IN):: J
    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndeg2l
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Invp(*)
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Aij(Ndeg,Ndeg)
    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(Ndeg2l,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(Ndeg,Ndeg,*)

    INTEGER:: i0
    INTEGER:: idbg
    INTEGER:: ii
    INTEGER:: itrans
    INTEGER:: j0
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: m
    INTEGER:: n
    DATA idbg/0/

    Ir = 0
    ii = Invp(I)
    jj = Invp(J)
    IF ( idbg/=0 ) WRITE (6,*) ii, jj, Aij
    IF ( ii==jj ) THEN
       l = 0
       DO n = 1, Ndeg
          DO m = 1, n
             l = l + 1
             Diag(l,ii) = Aij(n,m)
          ENDDO
       ENDDO
       RETURN
    ENDIF
    itrans = 0
    IF ( jj>ii ) THEN
       k = jj
       jj = ii
       ii = k
       itrans = 1
    ENDIF
    IF ( jj>=Nstop ) THEN
       i0 = ii - Nstop
       j0 = jj - Nstop + 1
       k = i0*(i0-1)/2 + j0
       IF ( itrans==0 ) THEN
          DO m = 1, Ndeg
             DO n = 1, Ndeg
                Dsln(n,m,k) = Aij(n,m)
             ENDDO
          ENDDO
       ELSE
          DO m = 1, Ndeg
             DO n = 1, Ndeg
                Dsln(n,m,k) = Aij(m,n)
             ENDDO
          ENDDO
       ENDIF
       RETURN
    ENDIF
    ks = Xlnzr(ii)
    ke = Xlnzr(ii+1) - 1
    DO k = ks, ke
       IF ( Colno(k)==jj ) THEN
          IF ( itrans==0 ) THEN
             DO m = 1, Ndeg
                DO n = 1, Ndeg
                   Zln(n,m,k) = Aij(n,m)
                ENDDO
             ENDDO
          ELSE
             DO m = 1, Ndeg
                DO n = 1, Ndeg
                   Zln(n,m,k) = Aij(m,n)
                ENDDO
             ENDDO
          ENDIF
          RETURN
       ENDIF
    ENDDO
    Ir = 20
  END SUBROUTINE ADDRX

  !======================================================================!
  !> @brief BRINGU  brings up zero pivots from bottom of the elimination tree to higher nodes
  !
  !      irr = 0     complete
  !          = 1     impossible
  !======================================================================!
  SUBROUTINE BRINGU(Zpiv,Iperm,Invp,Parent,Izz,Neqns,Irr)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Izz
    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Zpiv(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Irr
    INTEGER, INTENT(OUT):: Iperm(*)
    INTEGER, INTENT(OUT):: Invp(*)

    INTEGER:: i
    INTEGER:: ib
    INTEGER:: ib0
    INTEGER:: ibp
    INTEGER:: idbg
    INTEGER:: izzp

    idbg = 0
    Irr = 0
    ib0 = Invp(Izz)
    ib = ib0
    DO WHILE ( ib>0 )
       ibp = Parent(ib)
       izzp = Iperm(ibp)
       IF ( Zpiv(izzp)==0 ) THEN
          Invp(Izz) = ibp
          Invp(izzp) = ib0
          Iperm(ibp) = Izz
          Iperm(ib0) = izzp
          IF ( idbg/=0 ) THEN
             DO i = 1, Neqns
                IF ( Invp(Iperm(i))/=i .OR. Iperm(Invp(i))/=i) THEN
                   WRITE (6,*) 'permutation error'
                   STOP
                ENDIF
             ENDDO
             RETURN
          ENDIF
          RETURN
       ELSE
          ib = ibp
       ENDIF
    ENDDO
    Irr = 1
  END SUBROUTINE BRINGU

  !======================================================================!
  !> @brief D2DOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE D2DOT(T,A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(4,*)
    DOUBLE PRECISION, INTENT(IN):: B(4,*)
    DOUBLE PRECISION, INTENT(OUT):: T(4)

    INTEGER:: jj

    T(1) = 0.0D0
    T(2) = 0.0D0
    T(3) = 0.0D0
    T(4) = 0.0D0
    DO jj = 1, N
       T(1) = T(1) + A(1,jj)*B(1,jj) + A(3,jj)*B(3,jj)
       T(2) = T(2) + A(2,jj)*B(1,jj) + A(4,jj)*B(3,jj)
       T(3) = T(3) + A(1,jj)*B(2,jj) + A(3,jj)*B(4,jj)
       T(4) = T(4) + A(2,jj)*B(2,jj) + A(4,jj)*B(4,jj)
    ENDDO
  END SUBROUTINE D2DOT

  !======================================================================!
  !> @brief D2SDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE D2SDOT(Wi,A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(2,*)
    DOUBLE PRECISION, INTENT(IN):: B(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Wi(2)

    INTEGER:: jj

    DO jj = 1, N
       Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(3,jj)
       Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(4,jj)
    ENDDO
  END SUBROUTINE D2SDOT

  !======================================================================!
  !> @brief D3DOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE D3DOT(T,A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(9,*)
    DOUBLE PRECISION, INTENT(IN):: B(9,*)
    DOUBLE PRECISION, INTENT(OUT):: T(9)

    INTEGER:: jj
    INTEGER:: l

    !$dir max_trips(9)
    DO l = 1, 9
       T(l) = 0.0D0
    ENDDO
    DO jj = 1, N
       T(1) = T(1) + A(1,jj)*B(1,jj) + A(4,jj)*B(4,jj) + A(7,jj)*B(7,jj)
       T(2) = T(2) + A(2,jj)*B(1,jj) + A(5,jj)*B(4,jj) + A(8,jj)*B(7,jj)
       T(3) = T(3) + A(3,jj)*B(1,jj) + A(6,jj)*B(4,jj) + A(9,jj)*B(7,jj)
       T(4) = T(4) + A(1,jj)*B(2,jj) + A(4,jj)*B(5,jj) + A(7,jj)*B(8,jj)
       T(5) = T(5) + A(2,jj)*B(2,jj) + A(5,jj)*B(5,jj) + A(8,jj)*B(8,jj)
       T(6) = T(6) + A(3,jj)*B(2,jj) + A(6,jj)*B(5,jj) + A(9,jj)*B(8,jj)
       T(7) = T(7) + A(1,jj)*B(3,jj) + A(4,jj)*B(6,jj) + A(7,jj)*B(9,jj)
       T(8) = T(8) + A(2,jj)*B(3,jj) + A(5,jj)*B(6,jj) + A(8,jj)*B(9,jj)
       T(9) = T(9) + A(3,jj)*B(3,jj) + A(6,jj)*B(6,jj) + A(9,jj)*B(9,jj)
    ENDDO
  END SUBROUTINE D3DOT

  !======================================================================!
  !> @brief D3DOTL performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE D3DOTL(T,A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(9,*)
    DOUBLE PRECISION, INTENT(IN):: B(9,*)
    DOUBLE PRECISION, INTENT(OUT):: T(6)

    INTEGER:: jj
    INTEGER:: l

    DO l = 1, 6
       T(l) = 0.0D0
    ENDDO
    DO jj = 1, N
       T(1) = T(1) + A(1,jj)*B(1,jj) + A(4,jj)*B(4,jj) + A(7,jj)*B(7,jj)
       T(2) = T(2) + A(2,jj)*B(1,jj) + A(5,jj)*B(4,jj) + A(8,jj)*B(7,jj)
       T(3) = T(3) + A(2,jj)*B(2,jj) + A(5,jj)*B(5,jj) + A(8,jj)*B(8,jj)
       T(4) = T(4) + A(3,jj)*B(1,jj) + A(6,jj)*B(4,jj) + A(9,jj)*B(7,jj)
       T(5) = T(5) + A(3,jj)*B(2,jj) + A(6,jj)*B(5,jj) + A(9,jj)*B(8,jj)
       T(6) = T(6) + A(3,jj)*B(3,jj) + A(6,jj)*B(6,jj) + A(9,jj)*B(9,jj)
    ENDDO
  END SUBROUTINE D3DOTL

  !======================================================================!
  !> @brief D3SDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE D3SDOT(Wi,A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(3,*)
    DOUBLE PRECISION, INTENT(IN):: B(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Wi(3)

    INTEGER:: jj

    DO jj = 1, N
       Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(4,jj) - A(3,jj)*B(7,jj)
       Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(5,jj) - A(3,jj)*B(8,jj)
       Wi(3) = Wi(3) - A(1,jj)*B(3,jj) - A(2,jj)*B(6,jj) - A(3,jj)*B(9,jj)
    ENDDO
  END SUBROUTINE D3SDOT

  !======================================================================!
  !> @brief DDOT performs inner product of sparse vectors
  !======================================================================!
  DOUBLE PRECISION FUNCTION DDOT(A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(N)
    DOUBLE PRECISION, INTENT(IN):: B(N)

    INTEGER:: i
    DOUBLE PRECISION:: s

    s = 0.0D0
    DO i = 1, N
       s = s + A(i)*B(i)
    ENDDO
    DDOT = s
  END FUNCTION DDOT

  !======================================================================!
  !> @brief DXDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE DXDOT(Ndeg,T,A,B,L)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: L
    INTEGER, INTENT(IN):: Ndeg
    DOUBLE PRECISION, INTENT(IN):: A(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(IN):: B(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: T(Ndeg,Ndeg)

    INTEGER:: jj
    INTEGER:: k
    INTEGER:: m
    INTEGER:: n

    DO n = 1, Ndeg
       DO m = 1, Ndeg
          T(n,m) = 0.0D0
          DO k = 1, Ndeg
             DO jj = 1, L
                T(n,m) = T(n,m) + A(n,k,jj)*B(m,k,jj)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE DXDOT

  !======================================================================!
  !> @brief DXDOTL performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE DXDOTL(Ndeg,T,A,B,L)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: L
    INTEGER, INTENT(IN):: Ndeg
    DOUBLE PRECISION, INTENT(IN):: A(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(IN):: B(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: T(Ndeg,Ndeg)

    INTEGER:: jj
    INTEGER:: k
    INTEGER:: m
    INTEGER:: n

    DO n = 1, Ndeg
       DO m = 1, n
          T(n,m) = 0.0D0
          DO k = 1, Ndeg
             DO jj = 1, L
                T(n,m) = T(n,m) + A(n,k,jj)*B(m,k,jj)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE DXDOTL

  !======================================================================!
  !> @brief DXSDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE DXSDOT(Ndeg,Wi,A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    DOUBLE PRECISION, INTENT(IN):: A(Ndeg,*)
    DOUBLE PRECISION, INTENT(IN):: B(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Wi(Ndeg)
    INTEGER, INTENT(INOUT):: N
 
    INTEGER:: jj
    INTEGER:: m

    DO jj = 1, N
       DO m = 1, Ndeg
          DO N = 1, Ndeg
             Wi(N) = Wi(N) - B(N,m,jj)*A(m,jj)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE DXSDOT

  !======================================================================!
  !> @brief S3PDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE S3PDOT(Bi,B,Zln,Colno,Ks,Ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ke
    INTEGER, INTENT(IN):: Ks
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(9,*)
    DOUBLE PRECISION, INTENT(IN):: B(3,*)
    DOUBLE PRECISION, INTENT(OUT):: Bi(3)

    INTEGER:: j
    INTEGER:: jj

    DO jj = Ks, Ke
       j = Colno(jj)
       Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(4,jj)*B(2,j) - Zln(7,jj)*B(3,j)
       Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(5,jj)*B(2,j) - Zln(8,jj)*B(3,j)
       Bi(3) = Bi(3) - Zln(3,jj)*B(1,j) - Zln(6,jj)*B(2,j) - Zln(9,jj)*B(3,j)
    ENDDO
  END SUBROUTINE S3PDOT

  !======================================================================!
  !> @brief S2PDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE S2PDOT(Bi,B,Zln,Colno,Ks,Ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ke
    INTEGER, INTENT(IN):: Ks
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(4,*)
    DOUBLE PRECISION, INTENT(IN):: B(2,*)
    DOUBLE PRECISION, INTENT(OUT):: Bi(2)

    INTEGER:: j
    INTEGER:: jj

    DO jj = Ks, Ke
       j = Colno(jj)
       Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(3,jj)*B(2,j)
       Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(4,jj)*B(2,j)
    ENDDO
  END SUBROUTINE S2PDOT

  !======================================================================!
  !> @brief S6PDOT  performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE S6PDOT(Bi,B,Zln,Colno,Ks,Ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ke
    INTEGER, INTENT(IN):: Ks
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(36,*)
    DOUBLE PRECISION, INTENT(IN):: B(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Bi(6)

    INTEGER:: j
    INTEGER:: jj

    DO jj = Ks, Ke
       j = Colno(jj)
       Bi(1) = Bi(1) - Zln(1,jj)*B(1,j) - Zln(7,jj)*B(2,j)&
            - Zln(13,jj)*B(3,j) - Zln(19,jj)*B(4,j) - Zln(25,jj)&
            *B(5,j) - Zln(31,jj)*B(6,j)
       Bi(2) = Bi(2) - Zln(2,jj)*B(1,j) - Zln(8,jj)*B(2,j)&
            - Zln(14,jj)*B(3,j) - Zln(20,jj)*B(4,j) - Zln(26,jj)&
            *B(5,j) - Zln(32,jj)*B(6,j)
       Bi(3) = Bi(3) - Zln(3,jj)*B(1,j) - Zln(9,jj)*B(2,j)&
            - Zln(15,jj)*B(3,j) - Zln(21,jj)*B(4,j) - Zln(27,jj)&
            *B(5,j) - Zln(33,jj)*B(6,j)
       Bi(4) = Bi(4) - Zln(4,jj)*B(1,j) - Zln(10,jj)*B(2,j)&
            - Zln(16,jj)*B(3,j) - Zln(22,jj)*B(4,j) - Zln(28,jj)&
            *B(5,j) - Zln(34,jj)*B(6,j)
       Bi(5) = Bi(5) - Zln(5,jj)*B(1,j) - Zln(11,jj)*B(2,j)&
            - Zln(17,jj)*B(3,j) - Zln(23,jj)*B(4,j) - Zln(29,jj)&
            *B(5,j) - Zln(35,jj)*B(6,j)
       Bi(6) = Bi(6) - Zln(6,jj)*B(1,j) - Zln(12,jj)*B(2,j)&
            - Zln(18,jj)*B(3,j) - Zln(25,jj)*B(4,j) - Zln(30,jj)&
            *B(5,j) - Zln(36,jj)*B(6,j)
    ENDDO
  END SUBROUTINE S6PDOT

  !======================================================================!
  !> @brief SPDOT2 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  DOUBLE PRECISION FUNCTION SPDOT2(B,Zln,Colno,Ks,Ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ke
    INTEGER, INTENT(IN):: Ks
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(*)
    DOUBLE PRECISION, INTENT(IN):: B(*)

    INTEGER:: j
    INTEGER:: jj
    DOUBLE PRECISION:: s

    s = 0.0D0
    DO jj = Ks, Ke
       j = Colno(jj)
       s = s + Zln(jj)*B(j)
    ENDDO
    SPDOT2 = s
  END FUNCTION SPDOT2

  !======================================================================!
  !> @brief SXPDOT performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE SXPDOT(Ndeg,Bi,B,Zln,Colno,Ks,Ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ke
    INTEGER, INTENT(IN):: Ks
    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Colno(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(IN):: B(Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Bi(Ndeg)

    INTEGER:: j
    INTEGER:: jj
    INTEGER:: m
    INTEGER:: n

    DO jj = Ks, Ke
       j = Colno(jj)
       DO m = 1, Ndeg
          DO n = 1, Ndeg
             Bi(n) = Bi(n) - Zln(n,m,jj)*B(m,j)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE SXPDOT

  !======================================================================!
  !> @brief FORPAR
  !======================================================================!
  SUBROUTINE FORPAR(Neqns,Parent,Nch,Nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Nstop
    INTEGER, INTENT(OUT):: Nch(*)

    INTEGER:: i
    INTEGER:: idens
    INTEGER:: ii

    DO i = 1, Neqns
       Nch(i) = 0
    ENDDO
    Nch(Neqns+1) = 0
    DO i = 1, Neqns
       ii = Parent(i)
       Nch(ii) = Nch(ii) + 1
    ENDDO
    DO i = Neqns, 1, -1
       IF ( Nch(i)/=1 ) EXIT
    ENDDO

    idens = 0
    IF ( idens==1 ) THEN
       Nstop = i
    ELSE
       Nstop = Neqns + 1
    ENDIF
  END SUBROUTINE FORPAR

  !======================================================================!
  !> @brief SETIJ
  !======================================================================!
  SUBROUTINE SETIJ(Hecmesh,Hecmat)
    USE HECMW_UTIL
    IMPLICIT NONE

    TYPE (HECMWST_LOCAL_MESH), INTENT(IN)::Hecmesh
    TYPE (HECMWST_MATRIX), INTENT(IN)::Hecmat

    INTEGER:: i
    INTEGER:: ierr
    INTEGER:: j
    INTEGER:: k
    INTEGER:: kk
    INTEGER:: ntotal
    INTEGER(KIND=kint):: numnp
    INTEGER(KIND=kint):: ndof
    INTEGER(KIND=kint):: ndof2

    numnp = Hecmat%NP
    ndof = Hecmesh%N_DOF
    ntotal = numnp*ndof

    !*NUFACT variables
    NEQns = numnp
    NDEg = ndof
    NTTbr = Hecmat%NP + Hecmat%NPL
    !+hecMAT%NPU if unsymmetric
    ISYm = 0

    !*Allocations
    ALLOCATE (IROw(NTTbr),STAT=ierr)
    ALLOCATE (JCOl(NTTbr),STAT=ierr)
    IF ( ierr/=0 ) STOP "Allocation error: irow/jcol"

    kk = 0
    ndof2 = ndof*ndof
    DO j = 1, numnp
       !*Diagonal
       kk = kk + 1
       IROw(kk) = j
       JCOl(kk) = j
       !*Lower
       DO k = Hecmat%INDEXL(j-1) + 1, Hecmat%INDEXL(j)
          i = Hecmat%ITEML(k)
          kk = kk + 1
          IROw(kk) = j
          JCOl(kk) = i
       ENDDO
    ENDDO
  END SUBROUTINE SETIJ

  !======================================================================!
  !> @brief NUFORM
  !======================================================================!
  SUBROUTINE NUFORM(Hecmesh,Hecmat,Ir)
    USE HECMW_UTIL
    IMPLICIT NONE

    TYPE (HECMWST_LOCAL_MESH), INTENT(IN)::Hecmesh
    TYPE (HECMWST_MATRIX), INTENT(IN)::Hecmat
    INTEGER(KIND=kint), INTENT(OUT):: Ir

    INTEGER:: i
    INTEGER:: idbg
    INTEGER:: ierr
    INTEGER:: j
    INTEGER:: k
    INTEGER:: kk
    INTEGER:: ntotal
    INTEGER(KIND=kint):: numnp
    INTEGER(KIND=kint):: ndof
    INTEGER(KIND=kint):: ndof2

    idbg = 0
    numnp = Hecmat%NP
    ndof = Hecmesh%N_DOF
    ntotal = numnp*ndof

    !*NUFACT variables
    NEQns = numnp
    NDEg = ndof
    NTTbr = Hecmat%NP + Hecmat%NPL
    !+hecMAT%NPU if unsymmetric
    ISYm = 0

    !*Allocations
    ALLOCATE (VAL(NDEg*NDEg),STAT=ierr)
    IF ( ierr/=0 ) STOP "Allocation error:val"
    WRITE (6,*) "nuform:stage = ", STAge
    kk = 0
    ndof2 = ndof*ndof
    DO j = 1, numnp
       !*Diagonal
       kk = kk + 1
       CALL VLCPY(VAL,Hecmat%D(ndof2*(j-1)+1:ndof2*j),ndof)
       CALL STAIJ1(0,j,j,VAL,Ir)

       DO i = 1, ndof
          !           PAUSE 'Error?'
          IF ( VAL((i-1)*ndof+i)<=0 ) WRITE (idbg,*) 'j,j,val:', j, i, VAL((i-1)*ndof+i)
       ENDDO

       !*Lower
       DO k = Hecmat%INDEXL(j-1) + 1, Hecmat%INDEXL(j)
          i = Hecmat%ITEML(k)
          kk = kk + 1
          CALL VLCPY(VAL,Hecmat%AL(ndof2*(k-1)+1:ndof2*k),ndof)
          CALL STAIJ1(0,j,i,VAL,Ir)
       ENDDO
    ENDDO

    DEALLOCATE (VAL)
  END SUBROUTINE NUFORM

  !======================================================================!
  !> @brief VLCPY
  !======================================================================!
  SUBROUTINE VLCPY(A,B,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: B(N,N)
    DOUBLE PRECISION, INTENT(OUT):: A(N,N)

    INTEGER:: i
    INTEGER:: j

    DO i = 1, N
       DO j = 1, N
          A(j,i) = B(i,j)
       ENDDO
    ENDDO
  END SUBROUTINE VLCPY

  !======================================================================!
  !> @brief GENBTQ
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE GENBTQ(Invp,Parent,Btree,Zpiv,Izz,Neqns)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(IN):: Invp(*)
    INTEGER, INTENT(IN):: Zpiv(*)
    INTEGER, INTENT(OUT):: Btree(2,*)

    INTEGER:: i
    INTEGER:: ib
    INTEGER:: IDBg1
    INTEGER:: inext
    INTEGER:: ip
    INTEGER:: Izz
    COMMON /DEBUG / IDBg1

    DO i = 1, Neqns + 1
       Btree(1,i) = 0
       Btree(2,i) = 0
    ENDDO
    DO i = 1, Neqns + 1
       ip = Parent(i)
       IF ( ip>0 ) THEN
          ib = Btree(1,ip)
          IF ( ib==0 ) THEN
             Btree(1,ip) = i
          ELSE
             DO
                inext = Btree(2,ib)
                IF ( inext==0 ) THEN
                   Btree(2,ib) = i
                ELSE
                   ib = inext
                   CYCLE
                ENDIF
                EXIT
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    !
    ! find zeropivot
    !
    DO i = 1, Neqns
       IF ( Zpiv(i)/=0 ) THEN
          IF ( Btree(1,Invp(i))==0 ) THEN
             Izz = i
             GOTO 100
          ENDIF
       ENDIF
    ENDDO
    Izz = 0
100 IF ( IDBg1/=0 ) WRITE (6,"(' binary tree')")
    IF ( IDBg1/=0 ) WRITE (6,"(i6,'(',2I6,')')") (i,Btree(1,i),Btree(2,i),i=1,Neqns)
    IF ( IDBg1/=0 ) WRITE (6,"(' the first zero pivot is ',i4)") Izz
  END SUBROUTINE GENBTQ

  !======================================================================!
  !> @brief GENPAQ
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE GENPAQ(Xadj,Adjncy,Invp,Iperm,Parent,Neqns,Ancstr)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(IN):: Adjncy(*)
    INTEGER, INTENT(IN):: Invp(*)
    INTEGER, INTENT(IN):: Iperm(*)
    INTEGER, INTENT(OUT):: Parent(*)
    INTEGER, INTENT(OUT):: Ancstr(*)

    INTEGER:: i
    INTEGER:: IDBg1
    INTEGER:: ip
    INTEGER:: it
    INTEGER:: k
    INTEGER:: l
    COMMON /DEBUG / IDBg1

    DO i = 1, Neqns
       Parent(i) = 0
       Ancstr(i) = 0
       ip = Iperm(i)
       DO k = Xadj(ip), Xadj(ip+1) - 1
          l = Invp(Adjncy(k))
          IF ( l<i ) THEN
             DO WHILE ( Ancstr(l)/=0 )
                IF ( Ancstr(l)==i ) GOTO 50
                it = Ancstr(l)
                Ancstr(l) = i
                l = it
             ENDDO
             Ancstr(l) = i
             Parent(l) = i
          ENDIF
50     ENDDO
    ENDDO
    DO i = 1, Neqns
       IF ( Parent(i)==0 ) Parent(i) = Neqns + 1
    ENDDO
    Parent(Neqns+1) = 0
    IF ( IDBg1/=0 ) WRITE (6,"(' parent')")
    IF ( IDBg1/=0 ) WRITE (6,"(2I6)") (i,Parent(i),i=1,Neqns)
  END SUBROUTINE GENPAQ

  !======================================================================!
  !> @brief GENQMD
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE GENQMD(Neqns,Xadj,Adj0,Perm,Invp,Deg,Marker,Rchset,Nbrhd,Qsize,Qlink,Nofsub,Adjncy)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Adj0(*)
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(OUT):: Nofsub
    INTEGER, INTENT(OUT):: Adjncy(*)
    INTEGER, INTENT(OUT):: Perm(*)
    INTEGER, INTENT(OUT):: Invp(*)
    INTEGER, INTENT(OUT):: Deg(*)
    INTEGER, INTENT(OUT):: Marker(*)
    INTEGER, INTENT(OUT):: Rchset(*)
    INTEGER, INTENT(OUT):: Nbrhd(*)
    INTEGER, INTENT(OUT):: Qsize(*)
    INTEGER, INTENT(OUT):: Qlink(*)

    INTEGER:: i
    INTEGER:: inode
    INTEGER:: ip
    INTEGER:: irch
    INTEGER:: j
    INTEGER:: mindeg
    INTEGER:: ndeg
    INTEGER:: nhdsze
    INTEGER:: node
    INTEGER:: np
    INTEGER:: num
    INTEGER:: nump1
    INTEGER:: nxnode
    INTEGER:: rchsze
    INTEGER:: search
    INTEGER:: thresh

    mindeg = Neqns
    Nofsub = 0
    DO i = 1, Xadj(Neqns+1) - 1
       Adjncy(i) = Adj0(i)
    ENDDO
    DO node = 1, Neqns
       Perm(node) = node
       Invp(node) = node
       Marker(node) = 0
       Qsize(node) = 1
       Qlink(node) = 0
       ndeg = Xadj(node+1) - Xadj(node)
       Deg(node) = ndeg
       IF ( ndeg<mindeg ) mindeg = ndeg
    ENDDO

    num = 0
100 search = 1
    thresh = mindeg
    mindeg = Neqns
200 nump1 = num + 1
    IF ( nump1>search ) search = nump1
    DO j = search, Neqns
       node = Perm(j)
       IF ( Marker(node)>=0 ) THEN
          ndeg = Deg(node)
          IF ( ndeg<=thresh ) GOTO 300
          IF ( ndeg<mindeg ) mindeg = ndeg
       ENDIF
    ENDDO
    GOTO 100

300 search = j
    Nofsub = Nofsub + Deg(node)
    Marker(node) = 1
    CALL QMDRCH(node,Xadj,Adjncy,Deg,Marker,rchsze,Rchset,nhdsze,Nbrhd)
    nxnode = node
    DO
       num = num + 1
       np = Invp(nxnode)
       ip = Perm(num)
       Perm(np) = ip
       Invp(ip) = np
       Perm(num) = nxnode
       Invp(nxnode) = num
       Deg(nxnode) = -1
       nxnode = Qlink(nxnode)
       IF ( nxnode<=0 ) THEN
          IF ( rchsze>0 ) THEN
             CALL QMDUPD(Xadj,Adjncy,rchsze,Rchset,Deg,Qsize,Qlink,Marker,Rchset(rchsze+1),Nbrhd(nhdsze+1))
             Marker(node) = 0
             DO irch = 1, rchsze
                inode = Rchset(irch)
                IF ( Marker(inode)>=0 ) THEN
                   Marker(inode) = 0
                   ndeg = Deg(inode)
                   IF ( ndeg<mindeg ) mindeg = ndeg
                   IF ( ndeg<=thresh ) THEN
                      mindeg = thresh
                      thresh = ndeg
                      search = Invp(inode)
                   ENDIF
                ENDIF
             ENDDO
             IF ( nhdsze>0 ) CALL QMDOT(node,Xadj,Adjncy,Marker,rchsze,Rchset,Nbrhd)
          ENDIF
          IF ( num>=Neqns ) EXIT
          GOTO 200
       ENDIF
    ENDDO
  END SUBROUTINE GENQMD

  !======================================================================!
  !> @brief GNCLNO
  !======================================================================!
  SUBROUTINE GNCLNO(Parent,Xleaf,Leaf,Xlnzr,Colno,Neqns,Nstop,Lncol,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(IN):: Xleaf(*)
    INTEGER, INTENT(IN):: Leaf(*)
    INTEGER, INTENT(OUT):: Ir
    INTEGER, INTENT(OUT):: Lncol
    INTEGER, INTENT(OUT):: Xlnzr(*)
    INTEGER, INTENT(OUT):: Colno(*)

    INTEGER:: i
    INTEGER:: IDBg1
    INTEGER:: j
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: nc
    INTEGER:: nxleaf
    COMMON /DEBUG / IDBg1

    nc = 0
    Ir = 0
    l = 1
    DO i = 1, Neqns
       Xlnzr(i) = l
       ks = Xleaf(i)
       ke = Xleaf(i+1) - 1
       IF ( ke>=ks ) THEN
          nxleaf = Leaf(ks)
          DO k = ks, ke - 1
             j = nxleaf
             nxleaf = Leaf(k+1)
             DO WHILE ( j<nxleaf )
                IF ( j>=Nstop ) GOTO 100
                Colno(l) = j
                l = l + 1
                j = Parent(j)
             ENDDO
          ENDDO
          j = Leaf(ke)
          DO WHILE ( j<Nstop )
             IF ( j>=i .OR. j==0 ) EXIT
             Colno(l) = j
             l = l + 1
             j = Parent(j)
          ENDDO
       ENDIF
100 ENDDO
    Xlnzr(Neqns+1) = l
    Lncol = l - 1
    IF ( IDBg1/=0 ) WRITE (6,"(' xlnzr')")
    IF ( IDBg1/=0 ) WRITE (6,"(' colno (lncol =',i10,')')") Lncol
    IF ( IDBg1/=0 ) THEN
       DO k = 1, Neqns
          WRITE (6,"(/' row = ',i6)") k
          WRITE (6,"(10I4)") (Colno(i),i=Xlnzr(k),Xlnzr(k+1)-1)
       ENDDO
    ENDIF
  END SUBROUTINE GNCLNO

  !======================================================================!
  !> @brief GNLEAF
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE GNLEAF(Xadj,Adjncy,Invp,Iperm,Nch,Adjncp,Xleaf,Leaf,Neqns,Lnleaf)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(IN):: Adjncy(*)
    INTEGER, INTENT(IN):: Nch(*)
    INTEGER, INTENT(IN):: Invp(*)
    INTEGER, INTENT(IN):: Iperm(*)
    INTEGER, INTENT(OUT):: Lnleaf
    INTEGER, INTENT(OUT):: Adjncp(*)
    INTEGER, INTENT(OUT):: Xleaf(*)
    INTEGER, INTENT(OUT):: Leaf(*)

    INTEGER:: i
    INTEGER:: IDBg1
    INTEGER:: ik
    INTEGER:: ip
    INTEGER:: iq
    INTEGER:: istart
    INTEGER:: k
    INTEGER:: l
    INTEGER:: lc
    INTEGER:: lc1
    INTEGER:: m
    COMMON /DEBUG / IDBg1

    l = 1
    ik = 0
    istart = 0
    DO i = 1, Neqns
       Xleaf(i) = l
       ip = Iperm(i)
       DO k = Xadj(ip), Xadj(ip+1) - 1
          iq = Invp(Adjncy(k))
          IF ( iq<i ) THEN
             ik = ik + 1
             Adjncp(ik) = iq
          ENDIF
       ENDDO
       m = ik - istart
       IF ( m/=0 ) THEN
          CALL QQSORT(Adjncp(istart+1),m)
          lc1 = Adjncp(istart+1)
          IF ( lc1<i ) THEN
             Leaf(l) = lc1
             l = l + 1
             DO k = istart + 2, ik
                lc = Adjncp(k)

                IF ( lc1<lc-Nch(lc) ) THEN
                   Leaf(l) = lc
                   l = l + 1
                ENDIF

                lc1 = lc
             ENDDO
             ik = 1
             istart = ik
          ENDIF
       ENDIF
    ENDDO
    Xleaf(Neqns+1) = l
    Lnleaf = l - 1
    IF ( IDBg1/=0 ) WRITE (6,"(' xleaf')")
    IF ( IDBg1/=0 ) WRITE (6,"(10I6)") (Xleaf(i),i=1,Neqns+1)
    IF ( IDBg1/=0 ) WRITE (6,"(' leaf (len = ',i6,')')") Lnleaf
    IF ( IDBg1/=0 ) WRITE (6,"(10I6)") (Leaf(i),i=1,Lnleaf)
    RETURN
  END SUBROUTINE GNLEAF

  !======================================================================!
  !> @brief INV2
  !======================================================================!
  SUBROUTINE INV2(Dsln,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Dsln(3)

    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    Ir = 0
    IF ( DABS(Dsln(1))<RMIn ) THEN
       Ir = 10
       RETURN
    ENDIF
    Dsln(1) = 1.0D0/Dsln(1)
    t = Dsln(2)*Dsln(1)
    Dsln(3) = Dsln(3) - t*Dsln(2)
    Dsln(2) = t
    IF ( DABS(Dsln(3))<RMIn ) THEN
       Ir = 10
       RETURN
    ENDIF
    Dsln(3) = 1.0D0/Dsln(3)
  END SUBROUTINE INV2

  !======================================================================!
  !> @brief INV22
  !======================================================================!
  SUBROUTINE INV22(Zln,Zz,Diag)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: Diag(3)
    DOUBLE PRECISION, INTENT(IN):: Zz(4)
    DOUBLE PRECISION, INTENT(OUT):: Zln(4)

    Zln(3) = Zz(3) - Zz(1)*Diag(2)
    Zln(1) = Zz(1)*Diag(1)
    Zln(3) = Zln(3)*Diag(3)
    Zln(1) = Zln(1) - Zln(3)*Diag(2)

    Zln(4) = Zz(4) - Zz(2)*Diag(2)
    Zln(2) = Zz(2)*Diag(1)
    Zln(4) = Zln(4)*Diag(3)
    Zln(2) = Zln(2) - Zln(4)*Diag(2)
  END SUBROUTINE INV22

  !======================================================================!
  !> @brief INV3
  !======================================================================!
  SUBROUTINE INV3(Dsln,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Dsln(6)

    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t(2)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    Ir = 0
    IF ( DABS(Dsln(1))<RMIn ) GOTO 100
    Dsln(1) = 1.0D0/Dsln(1)
    t(1) = Dsln(2)*Dsln(1)
    Dsln(3) = Dsln(3) - t(1)*Dsln(2)
    Dsln(2) = t(1)
    IF ( DABS(Dsln(3))<RMIn ) GOTO 100
    Dsln(3) = 1.0D0/Dsln(3)
    t(1) = Dsln(4)*Dsln(1)
    Dsln(5) = Dsln(5) - Dsln(2)*Dsln(4)
    t(2) = Dsln(5)*Dsln(3)
    Dsln(6) = Dsln(6) - t(1)*Dsln(4) - t(2)*Dsln(5)
    Dsln(4) = t(1)
    Dsln(5) = t(2)
    IF ( DABS(Dsln(6))<RMIn ) GOTO 100
    Dsln(6) = 1.0D0/Dsln(6)
    RETURN
100 Dsln(1) = 1.0D0
    Dsln(2) = 0.0D0
    Dsln(3) = 1.0D0
    Dsln(4) = 0.0D0
    Dsln(5) = 0.0D0
    Dsln(6) = 1.0D0
  END SUBROUTINE INV3

  !======================================================================!
  !> @brief INV33
  !======================================================================!
  SUBROUTINE INV33(Zln,Zz,Diag)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: Diag(6)
    DOUBLE PRECISION, INTENT(IN):: Zz(9)
    DOUBLE PRECISION, INTENT(OUT):: Zln(9)

    Zln(4) = Zz(4) - Zz(1)*Diag(2)
    Zln(7) = Zz(7) - Zz(1)*Diag(4) - Zln(4)*Diag(5)
    Zln(1) = Zz(1)*Diag(1)
    Zln(4) = Zln(4)*Diag(3)
    Zln(7) = Zln(7)*Diag(6)
    Zln(4) = Zln(4) - Zln(7)*Diag(5)
    Zln(1) = Zln(1) - Zln(4)*Diag(2) - Zln(7)*Diag(4)

    Zln(5) = Zz(5) - Zz(2)*Diag(2)
    Zln(8) = Zz(8) - Zz(2)*Diag(4) - Zln(5)*Diag(5)
    Zln(2) = Zz(2)*Diag(1)
    Zln(5) = Zln(5)*Diag(3)
    Zln(8) = Zln(8)*Diag(6)
    Zln(5) = Zln(5) - Zln(8)*Diag(5)
    Zln(2) = Zln(2) - Zln(5)*Diag(2) - Zln(8)*Diag(4)

    Zln(6) = Zz(6) - Zz(3)*Diag(2)
    Zln(9) = Zz(9) - Zz(3)*Diag(4) - Zln(6)*Diag(5)
    Zln(3) = Zz(3)*Diag(1)
    Zln(6) = Zln(6)*Diag(3)
    Zln(9) = Zln(9)*Diag(6)
    Zln(6) = Zln(6) - Zln(9)*Diag(5)
    Zln(3) = Zln(3) - Zln(6)*Diag(2) - Zln(9)*Diag(4)
  END SUBROUTINE INV33

  !======================================================================!
  !> @brief INV6
  !======================================================================!
  SUBROUTINE INV6(Dsln,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Dsln(21)

    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t(5)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    Ir = 0
    Dsln(1) = 1.0D0/Dsln(1)
    t(1) = Dsln(2)*Dsln(1)
    Dsln(3) = 1.0D0/(Dsln(3)-t(1)*Dsln(2))
    Dsln(2) = t(1)
    Dsln(5) = Dsln(5) - Dsln(4)*Dsln(2)
    t(1) = Dsln(4)*Dsln(1)
    t(2) = Dsln(5)*Dsln(3)
    Dsln(6) = 1.0D0/(Dsln(6)-t(1)*Dsln(4)-t(2)*Dsln(5))
    Dsln(4) = t(1)
    Dsln(5) = t(2)
    Dsln(8) = Dsln(8) - Dsln(7)*Dsln(2)
    Dsln(9) = Dsln(9) - Dsln(7)*Dsln(4) - Dsln(8)*Dsln(5)
    t(1) = Dsln(7)*Dsln(1)
    t(2) = Dsln(8)*Dsln(3)
    t(3) = Dsln(9)*Dsln(6)
    Dsln(10) = 1.0D0/(Dsln(10)-t(1)*Dsln(7)-t(2)*Dsln(8)-t(3)*Dsln(9))
    Dsln(7) = t(1)
    Dsln(8) = t(2)
    Dsln(9) = t(3)
    Dsln(12) = Dsln(12) - Dsln(11)*Dsln(2)
    Dsln(13) = Dsln(13) - Dsln(11)*Dsln(4) - Dsln(12)*Dsln(5)
    Dsln(14) = Dsln(14) - Dsln(11)*Dsln(7) - Dsln(12)*Dsln(8) - Dsln(13)*Dsln(9)
    t(1) = Dsln(11)*Dsln(1)
    t(2) = Dsln(12)*Dsln(3)
    t(3) = Dsln(13)*Dsln(6)
    t(4) = Dsln(14)*Dsln(10)
    Dsln(15) = 1.0D0/(Dsln(15)-t(1)*Dsln(11)-t(2)*Dsln(12)-t(3)*Dsln(13)-t(4)*Dsln(14))
    Dsln(11) = t(1)
    Dsln(12) = t(2)
    Dsln(13) = t(3)
    Dsln(14) = t(4)
    Dsln(17) = Dsln(17) - Dsln(16)*Dsln(2)
    Dsln(18) = Dsln(18) - Dsln(16)*Dsln(4) - Dsln(17)*Dsln(5)
    Dsln(19) = Dsln(19) - Dsln(16)*Dsln(7) - Dsln(17)*Dsln(8) - Dsln(18)*Dsln(9)
    Dsln(20) = Dsln(20) - Dsln(16)*Dsln(11) - Dsln(17)*Dsln(12) - Dsln(18)*Dsln(13) - Dsln(19)*Dsln(14)
    t(1) = Dsln(16)*Dsln(1)
    t(2) = Dsln(17)*Dsln(3)
    t(3) = Dsln(18)*Dsln(6)
    t(4) = Dsln(19)*Dsln(10)
    t(5) = Dsln(20)*Dsln(15)
    Dsln(21) = 1.0D0/(Dsln(21)-t(1)*Dsln(16)-t(2)*Dsln(17)-t(3) *Dsln(18)-t(4)*Dsln(19)-t(5)*Dsln(20))
    Dsln(16) = t(1)
    Dsln(17) = t(2)
    Dsln(18) = t(3)
    Dsln(19) = t(4)
    Dsln(20) = t(5)
  END SUBROUTINE INV6

  !======================================================================!
  !> @brief INV66
  !======================================================================!
  SUBROUTINE INV66(Zln,Zz,Diag)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: Diag(21)
    DOUBLE PRECISION, INTENT(IN):: Zz(36)
    DOUBLE PRECISION, INTENT(OUT):: Zln(36)

    INTEGER:: i

    DO i = 0, 5
       Zln(i+7) = Zz(i+7) - Zz(i+1)*Diag(2)
       Zln(i+13) = Zz(i+13) - Zz(i+1)*Diag(4) - Zln(i+7)*Diag(5)
       Zln(i+19) = Zz(i+19) - Zz(i+1)*Diag(7) - Zln(i+7)*Diag(8) - Zln(i+13)*Diag(9)
       Zln(i+25) = Zz(i+25) - Zz(i+1)*Diag(11) - Zln(i+7)*Diag(12) - Zln(i+13)*Diag(13) - Zln(i+19)*Diag(14)
       Zln(i+31) = Zz(i+31) - Zz(i+1)*Diag(16) - Zln(i+7)*Diag(17) - Zln(i+13)*Diag(18) - Zln(i+19)*Diag(19) - Zln(i+25)*Diag(20)
       Zln(i+1) = Zz(i+1)*Diag(1)
       Zln(i+7) = Zln(i+7)*Diag(3)
       Zln(i+13) = Zln(i+13)*Diag(6)
       Zln(i+19) = Zln(i+19)*Diag(10)
       Zln(i+25) = Zln(i+25)*Diag(15)
       Zln(i+31) = Zln(i+31)*Diag(21)
       Zln(i+25) = Zln(i+25) - Zln(i+31)*Diag(20)
       Zln(i+19) = Zln(i+19) - Zln(i+31)*Diag(19) - Zln(i+25)*Diag(14)
       Zln(i+13) = Zln(i+13) - Zln(i+31)*Diag(18) - Zln(i+25)*Diag(13)- Zln(i+19)*Diag(9)
       Zln(i+7) = Zln(i+7) - Zln(i+31)*Diag(17) - Zln(i+25)*Diag(12)- Zln(i+19)*Diag(8) - Zln(i+13)*Diag(5)
       Zln(i+1) = Zln(i+1) - Zln(i+31)*Diag(16) - Zln(i+25)*Diag(11)- Zln(i+19)*Diag(7) - Zln(i+13)*Diag(4) - Zln(i+7)*Diag(2)
    ENDDO
  END SUBROUTINE INV66

  !======================================================================!
  !> @brief INVX
  !======================================================================!
  SUBROUTINE INVX(Dsln,Ndeg,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: Dsln(*)

    INTEGER:: i
    INTEGER:: j
    INTEGER:: k
    INTEGER:: k0
    INTEGER:: l
    INTEGER:: l0
    INTEGER:: ld
    INTEGER:: ll
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t
    DOUBLE PRECISION:: tem
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    Ir = 0
    l = 1
    Dsln(1) = 1.0D0/Dsln(1)
    DO i = 2, Ndeg
       ld = 0
       l0 = l
       DO j = 1, i - 1
          l = l + 1
          DO k = 1, j - 1
             ld = ld + 1
             Dsln(l) = Dsln(l) - Dsln(l0+k)*Dsln(ld)
          ENDDO
          ld = ld + 1
       ENDDO
       t = 0.0D0
       k0 = 0
       ll = 0
       DO k = l - i + 2, l
          ll = ll + 1
          k0 = k0 + ll
          tem = Dsln(k)*Dsln(k0)
          t = t + tem*Dsln(k)
          Dsln(k) = tem
       ENDDO
       l = l + 1
       Dsln(l) = Dsln(l) - t
       Dsln(l) = 1.0D0/Dsln(l)
    ENDDO
  END SUBROUTINE INVX

  !======================================================================!
  !> @brief INVXX
  !======================================================================!
  SUBROUTINE INVXX(Zln,Zz,Diag,Ndeg)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    DOUBLE PRECISION, INTENT(IN):: Diag(*)
    DOUBLE PRECISION, INTENT(IN):: Zz(Ndeg,Ndeg)
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg,Ndeg)

    INTEGER:: joc
    INTEGER:: l
    INTEGER:: loc1
    INTEGER:: m
    INTEGER:: n

    Zln = Zz
    DO l = 1, Ndeg, 2
       joc = 0
       DO m = 1, Ndeg - 1
          joc = joc + m
          loc1 = joc + m
          DO n = m + 1, Ndeg
             Zln(l,n) = Zln(l,n) - Zln(l,m)*Diag(loc1)
             Zln(l+1,n) = Zln(l+1,n) - Zln(l+1,m)*Diag(loc1)
             loc1 = loc1 + n
          ENDDO
       ENDDO
       joc = 0
       DO m = 1, Ndeg
          joc = joc + m
          Zln(l,m) = Zln(l,m)*Diag(joc)
          Zln(l+1,m) = Zln(l+1,m)*Diag(joc)
       ENDDO
       DO n = Ndeg, 2, -1
          joc = joc - 1
          DO m = n - 1, 1, -1
             Zln(l,m) = Zln(l,m) - Zln(l,n)*Diag(joc)
             Zln(l+1,m) = Zln(l+1,m) - Zln(l+1,n)*Diag(joc)
             joc = joc - 1
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE INVXX

  !======================================================================!
  !> @brief MATINI initializes storage for sparse matrix solver.
  !     this routine is used for both symmetric and asymmetric matrices
  !     and must be called once at the beginning
  !
  !    (i)
  !        neqns     number of unknowns
  !        nttbr     number of non0s, pattern of non-zero elements are
  !                  given like following.
  !                  nonz(A)={(i,j);i=irow(l),j=jcol(l); 1<= l <= nttbr}
  !        irow
  !        jcol      to define non-zero pattern
  !        lenv      length of the array v (iv)
  !
  !    (o)
  !        iv        comunication array. v is the original name
  !        ir        return code
  !                              =0    normal
  !                              =-1   non positive index
  !                              =1    too big index
  !                              =10   insufficient storage
  !
  !        contents of iv
  !               pointers 1 &zpiv(1)  2 &iperm(1)  3 &invp(1)
  !                        4 &parent(1)5 &nch(1)    6 &xlnzr(1)
  !                        7 &colno(1) 8 &diag(1)   9 &zln(1)
  !                       10 &dsln(1)
  !
  !               scalars 21 len(colno)  22 nstop     23 stage
  !                       24 neqns       25 len(iv) 26 len(dsln)
  !                       27 total
  !
  !        stage   10  after initialization
  !                20  building up matrix
  !                30  after LU decomposition
  !                40  after solving
  !
  !         #coded by t.arakawa of RIST on 040329
  !======================================================================!
  SUBROUTINE MATINI(Ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: Ir

    INTEGER:: IDBg
    INTEGER:: ir1
    INTEGER:: iv1
    INTEGER:: izz
    INTEGER:: izz0
    INTEGER:: ladp
    INTEGER:: last
    INTEGER:: lbtree
    INTEGER:: lcolno
    INTEGER:: lcpt
    INTEGER:: left
    INTEGER:: lenv
    INTEGER:: lenv2
    INTEGER:: lia
    INTEGER:: lja
    INTEGER:: lleaf
    INTEGER:: lncol
    INTEGER:: lnleaf
    INTEGER:: lpordr
    INTEGER:: LRAtio
    INTEGER:: lwk1
    INTEGER:: lwk2
    INTEGER:: lwk3
    INTEGER:: lwk4
    INTEGER:: lwk5
    INTEGER:: lwk6
    INTEGER:: lwk7
    INTEGER:: lwk8
    INTEGER:: lxleaf
    INTEGER:: maxl
    INTEGER:: neqns1
    INTEGER:: neqnsz
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    COMMON /DEBUG / IDBg
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    IDBg = 0
    izz0 = 0

    Ir = 0
    lenv = 1000000000
    lenv2 = LRAtio*lenv
    neqns1 = NEQns + 2
    LEN_dsln = lenv2
    LEN_iv = lenv2
    iv1 = 51

    !*Initialize allocation measure variables
    IALoc = 0
    RALoc = 0
    !
    !  set z pivot
    !
    ALLOCATE (ZPIv(NEQns),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, zpiv: SUB. matini"
    !ialoc = ialoc + neqns
    CALL ZPIVOT(NEQns,neqnsz,NTTbr,JCOl,IROw,ZPIv,ir1)
    IF ( ir1/=0 ) THEN
       Ir = ir1
       RETURN
    ENDIF

    !
    !  build jcpt,jcolno
    !
    !rmiv
    lcpt = iv1 + neqns1
    lcolno = lcpt + 2*NTTbr
    left = lcolno + 2*NTTbr
    last = lenv2
    !rmem
    ALLOCATE (JCPt(2*NTTbr),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, jcpt: SUB. matini"
    ALLOCATE (JCOlno(2*NTTbr),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, jcolno: SUB. matini"
    CALL STSMAT(NEQns,NTTbr,IROw,JCOl,JCPt,JCOlno)
    !
    !  build ia,ja
    !
    lia = last - neqns1
    lja = lia - NTTbr*2
    last = lja
    !rmem
    ALLOCATE (IA(NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, ia: SUB. matini"
    !WINDEBUG
    ALLOCATE (JA(2*NTTbr),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, ja: SUB. matini"
    CALL STIAJA(NEQns,IA,JA,JCPt,JCOlno)

    !*Deallocation of work array
    DEALLOCATE (JCPt)
    DEALLOCATE (JCOlno)
    !
    !  get permutation vector iperm,invp
    !
    lwk1 = lja - neqns1
    lwk2 = lwk1 - neqns1
    lwk3 = lwk2 - neqns1
    lwk4 = lwk3 - neqns1
    lwk5 = lwk4 - neqns1
    lwk6 = lwk5 - neqns1
    lwk7 = lwk6 - neqns1
    lwk8 = lwk7 - 2*NTTbr
    last = lwk8

    left = iv1 + 5*neqns1

    ALLOCATE (IPErm(NEQns),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, iperm: SUB. matini"
    ALLOCATE (INVp(NEQns),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, invp: SUB. matini"
    ALLOCATE (DEG(NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, deg: SUB. matini"
    ALLOCATE (MARker(NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, marker: SUB. matini"
    ALLOCATE (RCHset(0:NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, rchset: SUB. matini"
    ALLOCATE (NBRhd(NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, nbrhd: SUB. matini"
    ALLOCATE (QSIze(NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, qsize: SUB. matini"
    ALLOCATE (QLInk(NEQns+1),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, qlink: SUB. matini"
    ALLOCATE (ADJncy(2*NTTbr),STAT=IERror)
    IF ( IERror/=0 ) STOP "ALLOCATION ERROR, adjncy: SUB. matini"
    CALL GENQMD(neqnsz,IA,JA,IPErm,INVp,DEG,MARker,RCHset,NBRhd,QSIze,QLInk,NOFsub,ADJncy)
    DO
       !   build up the parent vector parent vector will be saved in
       !   work2 for a while
       CALL GENPAQ(IA,JA,INVp,IPErm,MARker,NEQns,RCHset)
       !
       !   build up the binary tree
       !
       lbtree = lwk3 - 2*NEQns
       last = lbtree
       !rmem
       ALLOCATE (BTRee(2*(NEQns+1)),STAT=IERror)
       IF ( IERror/=0 ) STOP "ALLOCATION ERROR, btree: SUB. matini"
       CALL GENBTQ(INVp,MARker,BTRee,ZPIv,izz,NEQns)
       !
       !   rotate the binary tree to avoid a zero pivot
       !
       IF ( izz==0 ) THEN
          !
          !   post ordering
          !
          lpordr = last - neqns1
          last = lpordr
          !rmem
          ALLOCATE (PARent(NEQns),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, parent: SUB. matini.f"
          ALLOCATE (NCH(NEQns+1),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, nch: SUB. matini.f"
          ALLOCATE (PORdr(NEQns+1),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, pordr: SUB. matini.f"
          CALL POSORD(PARent,BTRee,INVp,IPErm,PORdr,NCH,NEQns,DEG,MARker,RCHset)
          !
          !   generate skelton graph
          !
          lleaf = last - NTTbr
          lxleaf = lleaf - neqns1
          ladp = lxleaf - neqns1
          last = ladp
          !rmem
          ALLOCATE (ADJncp(NEQns+1),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, adjncp: SUB. matini.f"
          ALLOCATE (XLEaf(NEQns+1),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, xleaf: SUB. matini.f"
          ALLOCATE (LEAf(NTTbr),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, leaf: SUB. matini.f"
          CALL GNLEAF(IA,JA,INVp,IPErm,NCH,ADJncp,XLEaf,LEAf,NEQns,lnleaf)
          CALL FORPAR(NEQns,PARent,NCH,NSTop)
          !*Deallocation of work arrays
          DEALLOCATE (IA)
          DEALLOCATE (JA)
          DEALLOCATE (DEG)
          DEALLOCATE (MARker)
          DEALLOCATE (RCHset)
          DEALLOCATE (NBRhd)
          DEALLOCATE (QSIze)
          DEALLOCATE (QLInk)
          DEALLOCATE (ADJncy)
          DEALLOCATE (ZPIv)
          !*Nullify pointers
          NULLIFY (IA)
          NULLIFY (JA)
          NULLIFY (DEG)
          NULLIFY (MARker)
          NULLIFY (RCHset)
          NULLIFY (NBRhd)
          NULLIFY (QSIze)
          NULLIFY (QLInk)
          NULLIFY (ADJncy)
          NULLIFY (ZPIv)
          !
          !   build up xlnzr,colno  (this is the symbolic fct.)
          !
          maxl = lxleaf - (left+neqns1)
          ALLOCATE (XLNzr(NEQns+1),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, xlnzr: SUB. matini.f"
          CALL PRE_GNCLNO(PARent,XLEaf,LEAf,XLNzr,NEQns,NSTop,lncol,ir1)
          ALLOCATE (COLno(lncol),STAT=IERror)
          IF ( IERror/=0 ) STOP "ALLOCATION ERROR, colno: SUB. matini.f"
          CALL GNCLNO(PARent,XLEaf,LEAf,XLNzr,COLno,NEQns,NSTop,lncol,ir1)

          !*Deallocate work arrays
          DEALLOCATE (PORdr)
          DEALLOCATE (ADJncp)
          DEALLOCATE (XLEaf)
          DEALLOCATE (LEAf)
          DEALLOCATE (BTRee)
          !*Nullify pointers
          NULLIFY (PORdr)
          NULLIFY (ADJncp)
          NULLIFY (XLEaf)
          NULLIFY (LEAf)
          NULLIFY (BTRee)
          !rmem
          left = (left+neqns1) + lncol
          !rmiv
          LEN_dsln = (NEQns-NSTop+1)*(NEQns-NSTop)/2

          !Scalar assignments
          LEN_colno = lncol
          !
          !   area for double precision values
          !
          IF ( MOD(left,2)==0 ) left = left + 1

          !rmiv
          TOTal = left
          !rmiv
          STAge = 10
          IALoc = 5*NEQns + lncol + 1
          EXIT
       ELSE
          IF ( izz0==0 ) izz0 = izz
          IF ( izz0/=izz ) THEN

             CALL BRINGU(ZPIv,IPErm,INVp,MARker,izz,NEQns,IRR)
          ELSE
             lwk4 = last - neqns1
             lwk5 = lwk4 - neqns1
             last = lwk5
             CALL ROTATE(IA,JA,INVp,IPErm,MARker,BTRee,izz,NEQns,NBRhd,QSIze,IRR)
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE MATINI

  !======================================================================!
  !> @brief NUFCT performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !         #coded by t.arakawa of RIST on 040329
  !======================================================================!
  SUBROUTINE NUFCT(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Ir
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(*)

    INTEGER:: ic
    INTEGER:: ISEm
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio
    COMMON ISEm

    ISEm = 1
    !
    ! phase I
    !
    CALL PTIME(t1)
    Diag(1) = 1.0D0/Diag(1)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    DO ic = 2, Nstop - 1
       CALL SUM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    ENDDO
    !
    ! phase II
    !
    CALL PTIME(t2)
    DO ic = Nstop, Neqns
       CALL SUM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    ENDDO
    !
    ! phase III
    !
    CALL PTIME(t3)
    CALL SUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    CALL PTIME(t4)
    CALL SUM3(Neqns-Nstop+1,Dsln,Diag(Nstop),Indx,Temp)
    CALL PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
    RETURN
    Ir = 30
  END SUBROUTINE NUFCT

  !======================================================================!
  !> @brief NUFCT0 performs Cholesky factorization
  !          if(iv(22).eq.0)    normal type
  !          if(iv(22).gt.0)    code generation type
  !======================================================================!
  SUBROUTINE NUFCT0(Ir)
    USE HECMW_UTIL
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: Ir

    INTEGER:: ISEed
    INTEGER:: IXXx
    INTEGER:: LRAtio
    INTEGER:: ndeg2
    INTEGER:: ndegl
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio
    COMMON /QAZ   / ISEed, IXXx

    IF ( STAge/=20 ) THEN
       PRINT *, '*********Setting Stage 40!*********'
       Ir = 40
       RETURN
    ELSE
       Ir = 0
    ENDIF
    !
    ALLOCATE (TEMp(NDEg*NDEg*NEQns),STAT=IRR)
    IF ( IRR/=0 ) THEN
       WRITE (*,*) '##Error : Not enough memory'
       CALL HECMW_ABORT(HECMW_COMM_GET_COMM())
       !stop
    ENDIF
    ALLOCATE (INDx(NEQns),STAT=IRR)
    IF ( IRR/=0 ) THEN
       WRITE (*,*) '##Error : Not enough memory'
       CALL HECMW_ABORT(HECMW_COMM_GET_COMM())
       !stop
    ENDIF
    !
    !rmiv
    ndegl = NDEg*(NDEg+1)
    ndegl = ndegl/2
    ndeg2 = NDEg*NDEg
    !rmiv
    IF ( NDEg==1 ) THEN
       CALL NUFCT(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
    ELSEIF ( NDEg==2 ) THEN
       CALL NUFCT2(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
    ELSEIF ( NDEg==3 ) THEN
       CALL NUFCT3(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
    ELSEIF ( NDEg==6 ) THEN
       CALL NUFCT6(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,Ir)
    ELSE
       CALL NUFCTX(XLNzr,COLno,DSLn,ZLN,DIAg,INDx,TEMp,NEQns,PARent,NCH,NSTop,NDEg,ndegl,Ir)
    ENDIF
    STAge = 30
    DEALLOCATE (TEMp)
    DEALLOCATE (INDx)
  END SUBROUTINE NUFCT0

  !======================================================================!
  !> @brief NUFCT2  performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !         #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE NUFCT2(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Ir
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(3,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(4,*)

    INTEGER:: ic
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    !
    ! phase I
    !
    Ir = 0
    CALL PTIME(t1)
    IF ( Nstop>1 ) CALL INV2(Diag(1,1),Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    DO ic = 2, Nstop - 1
       CALL S2UM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    ENDDO
    !
    ! phase II
    !
    CALL PTIME(t2)
    DO ic = Nstop, Neqns
       CALL S2UM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    ENDDO
    !
    ! phase III
    !
    CALL PTIME(t3)
    CALL S2UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    CALL PTIME(t4)
    CALL S2UM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp)
    CALL PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  END SUBROUTINE NUFCT2

  !======================================================================!
  !> @brief NUFCT3 performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !         #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE NUFCT3(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Ir
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(9,*)

    INTEGER:: ic
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    !
    ! phase I
    !
    CALL PTIME(t1)
    IF ( Nstop>1 ) CALL INV3(Diag(1,1),Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    DO ic = 2, Nstop - 1
       CALL S3UM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    ENDDO
    !
    ! phase II
    !
    CALL PTIME(t2)
    DO ic = Nstop, Neqns
       CALL S3UM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    ENDDO
    !
    ! phase III
    !
    CALL PTIME(t3)
    CALL S3UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    CALL PTIME(t4)
    CALL S3UM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp)
    CALL PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  END SUBROUTINE NUFCT3

  !======================================================================!
  !> @brief NUFCT6 performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !         #coded by t.arakawa of RIST 040331
  !======================================================================!
  SUBROUTINE NUFCT6(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(36,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(21,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(36,*)
    INTEGER:: ic
    INTEGER:: Ir
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    !
    ! phase I
    !
    CALL PTIME(t1)
    IF ( Nstop>1 ) CALL INV6(Diag(1,1),Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    DO ic = 2, Nstop - 1
       CALL S6UM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx)
    ENDDO
    !
    ! phase II
    !
    CALL PTIME(t2)
    DO ic = Nstop, Neqns
       CALL S6UM1(ic,Xlnzr,Colno,Zln,Temp,Indx)
    ENDDO
    !
    ! phase III
    !
    CALL PTIME(t3)
    CALL S6UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    !
    ! phase IV
    !
    CALL PTIME(t4)
    CALL S6UM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp)
    CALL PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  END SUBROUTINE NUFCT6

  !======================================================================!
  !> @brief NUFCTX performs cholesky factorization in row order
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !     (o) zln,diag,dsln
  !         #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE NUFCTX(Xlnzr,Colno,Dsln,Zln,Diag,Indx,Temp,Neqns,Parent,Nch,Nstop,Ndeg,Ndegl,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndegl
    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg*Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(Ndegl,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(Ndeg*Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(Ndeg*Ndeg,*)

    INTEGER:: ic
    INTEGER:: Ir
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    DOUBLE PRECISION:: zz(100)
    DOUBLE PRECISION:: t(100)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    !
    ! phase I
    !
    CALL PTIME(t1)
    IF ( Nstop>1 ) CALL INVX(Diag(1,1),Ndeg,Ir)
    l = Parent(1)
    Nch(l) = Nch(l) - 1
    Nch(1) = -1
    DO ic = 2, Nstop - 1
       CALL SXUM(ic,Xlnzr,Colno,Zln,Diag,Nch,Parent,Temp,Indx,Ndeg,Ndegl,zz,t)
    ENDDO
    !
    ! phase II
    !
    CALL PTIME(t2)
    DO ic = Nstop, Neqns
       CALL SXUM1(ic,Xlnzr,Colno,Zln,Temp,Indx,Ndeg,t)
    ENDDO
    !
    ! phase III
    !
    CALL PTIME(t3)
    CALL SXUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx,Ndeg,Ndegl)
    !
    ! phase IV
    !
    CALL PTIME(t4)
    CALL SXUM3(Neqns-Nstop+1,Dsln,Diag(1,Nstop),Indx,Temp,Ndeg,Ndegl,t)
    CALL PTIME(t5)
    tt = t5 - t1
    t1 = t2 - t1
    t2 = t3 - t2
    t3 = t4 - t3
    t4 = t5 - t4
  END SUBROUTINE NUFCTX

  !======================================================================!
  !> @brief NUSOL0 performs forward elimination and backward substitution
  !     (i/o)
  !           r_h_s    on entry     right hand side vector
  !                    on exit      solution vector
  !           iv       communication array
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE NUSOL0(R_h_s,Ir)
    USE HECMW_UTIL
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: Ir
    DOUBLE PRECISION, INTENT(OUT):: R_h_s(*)

    INTEGER:: ISEed
    INTEGER:: IXXx
    INTEGER:: LRAtio
    INTEGER:: lwk
    INTEGER:: ndegl
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION, POINTER :: wk(:)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio
    COMMON /QAZ   / ISEed, IXXx

    IF ( STAge/=30 .AND. STAge/=40 ) THEN
       Ir = 50
       RETURN
    ELSE
       Ir = 0
    ENDIF
    lwk = TOTal

    ALLOCATE (wk(NDEg*NEQns),STAT=IERror)
    IF ( IERror/=0 ) THEN
       WRITE (*,*) "##Error: not enough memory"
       CALL HECMW_ABORT(HECMW_COMM_GET_COMM())
    ENDIF
    !rmiv
    ndegl = NDEg*(NDEg+1)
    ndegl = ndegl/2
    IF ( NDEg==1 ) THEN
       CALL NUSOL1(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop)
    ELSEIF ( NDEg==2 ) THEN
       CALL NUSOL2(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop)
    ELSEIF ( NDEg==3 ) THEN
       CALL NUSOL3(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop)
    ELSEIF ( NDEg==6 ) THEN
       CALL NUSOLX(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop,NDEg,ndegl)
    ELSE
       CALL NUSOLX(XLNzr,COLno,DSLn,ZLN,DIAg,IPErm,R_h_s,wk,NEQns,NSTop,NDEg,ndegl)
    ENDIF
    STAge = 40
    DEALLOCATE (wk)
  END SUBROUTINE NUSOL0

  !======================================================================!
  !> @brief NUSOL1 performs forward elimination and backward substitution
  !======================================================================!
  SUBROUTINE NUSOL1(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Iperm(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(*)
    DOUBLE PRECISION, INTENT(IN):: Diag(*)
    DOUBLE PRECISION, INTENT(IN):: Dsln(*)
    DOUBLE PRECISION, INTENT(OUT):: B(*)
    DOUBLE PRECISION, INTENT(OUT):: Wk(*)

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    ! forward
    DO i = 1, Neqns
       Wk(i) = B(Iperm(i))
    ENDDO
    joc = 1
    DO i = 1, Neqns
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) Wk(i) = Wk(i) - SPDOT2(Wk,Zln,Colno,ks,ke)
       IF ( i>Nstop ) THEN
          Wk(i) = Wk(i) - DDOT(Wk(Nstop),Dsln(joc),i-Nstop)
          joc = joc + i - Nstop
       ENDIF
    ENDDO
    DO i = 1, Neqns
       Wk(i) = Wk(i)*Diag(i)
    ENDDO
    ! back ward
    DO i = Neqns, 1, -1
       IF ( i>=Nstop ) THEN
          DO j = i - 1, Nstop, -1
             joc = joc - 1
             Wk(j) = Wk(j) - Wk(i)*Dsln(joc)
          ENDDO
       ENDIF
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) THEN
          DO k = ks, ke
             j = Colno(k)
             Wk(j) = Wk(j) - Wk(i)*Zln(k)
          ENDDO
       ENDIF
    ENDDO
    ! permutaion
    DO i = 1, Neqns
       B(Iperm(i)) = Wk(i)
    ENDDO
  END SUBROUTINE NUSOL1

  !======================================================================!
  !> @brief NUSOL2 performs forward elimination and backward substitution
  !======================================================================!
  SUBROUTINE NUSOL2(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Iperm(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(4,*)
    DOUBLE PRECISION, INTENT(IN):: Diag(3,*)
    DOUBLE PRECISION, INTENT(IN):: Dsln(4,*)
    DOUBLE PRECISION, INTENT(OUT):: B(2,*)
    DOUBLE PRECISION, INTENT(OUT):: Wk(2,*)

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    ! forward
    DO i = 1, Neqns
       Wk(1,i) = B(1,Iperm(i))
       Wk(2,i) = B(2,Iperm(i))
    ENDDO
    joc = 1
    DO i = 1, Neqns
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) CALL S2PDOT(Wk(1,i),Wk,Zln,Colno,ks,ke)
       IF ( i>Nstop ) THEN
          CALL D2SDOT(Wk(1,i),Wk(1,Nstop),Dsln(1,joc),i-Nstop)
          joc = joc + i - Nstop
       ENDIF
    ENDDO
    DO i = 1, Neqns
       Wk(2,i) = Wk(2,i) - Wk(1,i)*Diag(2,i)
       Wk(1,i) = Wk(1,i)*Diag(1,i)
       Wk(2,i) = Wk(2,i)*Diag(3,i)
       Wk(1,i) = Wk(1,i) - Wk(2,i)*Diag(2,i)
    ENDDO
    ! back ward
    DO i = Neqns, 1, -1
       IF ( i>=Nstop ) THEN
          DO j = i - 1, Nstop, -1
             joc = joc - 1
             Wk(1,j) = Wk(1,j) - Wk(1,i)*Dsln(1,joc) - Wk(2,i)*Dsln(2,joc)
             Wk(2,j) = Wk(2,j) - Wk(1,i)*Dsln(3,joc) - Wk(2,i)*Dsln(4,joc)
          ENDDO
       ENDIF
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) THEN
          DO k = ks, ke
             j = Colno(k)
             Wk(1,j) = Wk(1,j) - Wk(1,i)*Zln(1,k) - Wk(2,i)*Zln(2,k)
             Wk(2,j) = Wk(2,j) - Wk(1,i)*Zln(3,k) - Wk(2,i)*Zln(4,k)
          ENDDO
       ENDIF
    ENDDO
    ! permutaion
    DO i = 1, Neqns
       B(1,Iperm(i)) = Wk(1,i)
       B(2,Iperm(i)) = Wk(2,i)
    ENDDO
  END SUBROUTINE NUSOL2

  !======================================================================!
  !> @brief NUSOL3 performs forward elimination and backward substitution
  !======================================================================!
  SUBROUTINE NUSOL3(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Iperm(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(9,*)
    DOUBLE PRECISION, INTENT(IN):: Diag(6,*)
    DOUBLE PRECISION, INTENT(IN):: Dsln(9,*)
    DOUBLE PRECISION, INTENT(OUT):: B(3,*)
    DOUBLE PRECISION, INTENT(OUT):: Wk(3,*)

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    ! forward
    DO i = 1, Neqns
       Wk(1,i) = B(1,Iperm(i))
       Wk(2,i) = B(2,Iperm(i))
       Wk(3,i) = B(3,Iperm(i))
    ENDDO
    joc = 1
    DO i = 1, Neqns
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) CALL S3PDOT(Wk(1,i),Wk,Zln,Colno,ks,ke)
       IF ( i>Nstop ) THEN
          CALL D3SDOT(Wk(1,i),Wk(1,Nstop),Dsln(1,joc),i-Nstop)
          joc = joc + i - Nstop
       ENDIF
    ENDDO
    DO i = 1, Neqns
       Wk(2,i) = Wk(2,i) - Wk(1,i)*Diag(2,i)
       Wk(3,i) = Wk(3,i) - Wk(1,i)*Diag(4,i) - Wk(2,i)*Diag(5,i)
       Wk(1,i) = Wk(1,i)*Diag(1,i)
       Wk(2,i) = Wk(2,i)*Diag(3,i)
       Wk(3,i) = Wk(3,i)*Diag(6,i)
       Wk(2,i) = Wk(2,i) - Wk(3,i)*Diag(5,i)
       Wk(1,i) = Wk(1,i) - Wk(2,i)*Diag(2,i) - Wk(3,i)*Diag(4,i)
    ENDDO
    ! back ward
    DO i = Neqns, 1, -1
       IF ( i>=Nstop ) THEN
          DO j = i - 1, Nstop, -1
             joc = joc - 1
             Wk(1,j) = Wk(1,j) - Wk(1,i)*Dsln(1,joc) - Wk(2,i)*Dsln(2,joc) - Wk(3,i)*Dsln(3,joc)
             Wk(2,j) = Wk(2,j) - Wk(1,i)*Dsln(4,joc) - Wk(2,i)*Dsln(5,joc) - Wk(3,i)*Dsln(6,joc)
             Wk(3,j) = Wk(3,j) - Wk(1,i)*Dsln(7,joc) - Wk(2,i)*Dsln(8,joc) - Wk(3,i)*Dsln(9,joc)
          ENDDO
       ENDIF
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) THEN
          DO k = ks, ke
             j = Colno(k)
             Wk(1,j) = Wk(1,j) - Wk(1,i)*Zln(1,k) - Wk(2,i)*Zln(2,k) - Wk(3,i)*Zln(3,k)
             Wk(2,j) = Wk(2,j) - Wk(1,i)*Zln(4,k) - Wk(2,i)*Zln(5,k) - Wk(3,i)*Zln(6,k)
             Wk(3,j) = Wk(3,j) - Wk(1,i)*Zln(7,k) - Wk(2,i)*Zln(8,k) - Wk(3,i)*Zln(9,k)
          ENDDO
       ENDIF
    ENDDO
    ! permutaion
    DO i = 1, Neqns
       B(1,Iperm(i)) = Wk(1,i)
       B(2,Iperm(i)) = Wk(2,i)
       B(3,Iperm(i)) = Wk(3,i)
    ENDDO
  END SUBROUTINE NUSOL3

  !======================================================================!
  !> @brief NUSOLX performs forward elimination and backward substitution
  !======================================================================!
  SUBROUTINE NUSOLX(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop,Ndeg,Ndegl)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndegl
    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Iperm(*)
    DOUBLE PRECISION, INTENT(IN):: Zln(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(IN):: Diag(Ndegl,*)
    DOUBLE PRECISION, INTENT(OUT):: B(Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Wk(Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(Ndeg,Ndeg,*)

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: joc1
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: loc1
    INTEGER:: locd
    INTEGER:: m
    INTEGER:: n

    ! forward
    DO l = 1, Ndeg
       DO i = 1, Neqns
          Wk(l,i) = B(l,Iperm(i))
       ENDDO
    ENDDO
    joc = 1
    DO i = 1, Neqns
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) CALL SXPDOT(Ndeg,Wk(1,i),Wk,Zln,Colno,ks,ke)
       IF ( i>Nstop ) THEN
          joc1 = i - Nstop
          CALL DXSDOT(Ndeg,Wk(1,i),Wk(1,Nstop),Dsln(1,1,joc),joc1)
          joc = joc + joc1
       ENDIF
    ENDDO
    DO i = 1, Neqns
       locd = 0
       DO m = 1, Ndeg - 1
          locd = locd + m
          loc1 = locd + m
          DO n = m + 1, Ndeg
             Wk(n,i) = Wk(n,i) - Wk(m,i)*Diag(loc1,i)
             loc1 = loc1 + n
          ENDDO
       ENDDO
       locd = 0
       DO m = 1, Ndeg
          locd = locd + m
          Wk(m,i) = Wk(m,i)*Diag(locd,i)
       ENDDO
       DO n = Ndeg, 2, -1
          locd = locd - 1
          DO m = n - 1, 1, -1
             Wk(m,i) = Wk(m,i) - Wk(n,i)*Diag(locd,i)
             locd = locd - 1
          ENDDO
       ENDDO
    ENDDO
    ! back ward
    DO i = Neqns, 1, -1
       IF ( i>=Nstop ) THEN
          DO j = i - 1, Nstop, -1
             joc = joc - 1
             DO m = 1, Ndeg
                DO n = 1, Ndeg
                   Wk(m,j) = Wk(m,j) - Wk(n,i)*Dsln(n,m,joc)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ks = Xlnzr(i)
       ke = Xlnzr(i+1) - 1
       IF ( ke>=ks ) THEN
          DO k = ks, ke
             j = Colno(k)
             DO m = 1, Ndeg
                DO n = 1, Ndeg
                   Wk(m,j) = Wk(m,j) - Wk(n,i)*Zln(n,m,k)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! permutaion
    DO l = 1, Ndeg
       DO i = 1, Neqns
          B(l,Iperm(i)) = Wk(l,i)
       ENDDO
    ENDDO
  END SUBROUTINE NUSOLX

  !======================================================================!
  !> @brief POSORD
  !======================================================================!
  SUBROUTINE POSORD(Parent,Btree,Invp,Iperm,Pordr,Nch,Neqns,Iw,Qarent,Mch)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Btree(2,*)
    INTEGER, INTENT(IN):: Qarent(*)
    INTEGER, INTENT(OUT):: Parent(*)
    INTEGER, INTENT(OUT):: Pordr(*)
    INTEGER, INTENT(OUT):: Nch(*)
    INTEGER, INTENT(OUT):: Invp(*)
    INTEGER, INTENT(OUT):: Iperm(*)
    INTEGER, INTENT(OUT):: Iw(*)
    INTEGER, INTENT(OUT):: Mch(0:Neqns+1)

    INTEGER:: i
    INTEGER:: IDBg1
    INTEGER:: ii
    INTEGER:: invpos
    INTEGER:: ipinv
    INTEGER:: joc
    INTEGER:: l
    INTEGER:: locc
    INTEGER:: locp
    COMMON /DEBUG / IDBg1

    DO i = 1, Neqns
       Mch(i) = 0
       Pordr(i) = 0
    ENDDO
    l = 1
    locc = Neqns + 1
    DO
       joc = locc
       locc = Btree(1,joc)
       IF ( locc==0 ) THEN
          locp = Qarent(joc)
          Mch(locp) = Mch(locp) + 1
          DO
             Pordr(joc) = l
             IF ( l>=Neqns ) THEN
                DO i = 1, Neqns
                   ipinv = Pordr(Invp(i))
                   Invp(i) = ipinv
                   Iperm(ipinv) = i
                   Iw(Pordr(i)) = i
                ENDDO
                DO i = 1, Neqns
                   invpos = Iw(i)
                   Nch(i) = Mch(invpos)
                   ii = Qarent(invpos)
                   IF ( ii>0 .AND. ii<=Neqns ) THEN
                      Parent(i) = Pordr(ii)
                   ELSE
                      Parent(i) = Qarent(invpos)
                   ENDIF
                ENDDO
                IF ( IDBg1/=0 ) THEN
                   WRITE (6,"(' post order')")
                   WRITE (6,"(10I6)") (Pordr(i),i=1,Neqns)
                   WRITE (6,"(/' invp ')")
                   WRITE (6,"(/' parent')")
                   WRITE (6,"(10I6)") (Parent(i),i=1,Neqns)
                   WRITE (6,"(10I6)") (Invp(i),i=1,Neqns)
                   WRITE (6,"(/' iperm ')")
                   WRITE (6,"(10I6)") (Iperm(i),i=1,Neqns)
                   WRITE (6,"(' nch')")
                   WRITE (6,"(10I6)") (Nch(i),i=1,Neqns)
                ENDIF
                RETURN
             ELSE
                l = l + 1
                locc = Btree(2,joc)
                IF ( locc/=0 ) EXIT
                joc = Qarent(joc)
                locp = Qarent(joc)
                Mch(locp) = Mch(locp) + Mch(joc) + 1
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE POSORD

  !======================================================================!
  !> @brief PRE_GNCLNO
  !======================================================================!
  SUBROUTINE PRE_GNCLNO(Parent,Xleaf,Leaf,Xlnzr,Neqns,Nstop,Lncol,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(IN):: Xleaf(*)
    INTEGER, INTENT(IN):: Leaf(*)
    INTEGER, INTENT(OUT):: Lncol
    INTEGER, INTENT(OUT):: Xlnzr(*)

    INTEGER:: i
    INTEGER:: IDBg1
    INTEGER:: Ir
    INTEGER:: j
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: nc
    INTEGER:: nxleaf
    COMMON /DEBUG / IDBg1

    nc = 0
    Ir = 0
    l = 1
    DO i = 1, Neqns
       Xlnzr(i) = l
       ks = Xleaf(i)
       ke = Xleaf(i+1) - 1
       IF ( ke>=ks ) THEN
          nxleaf = Leaf(ks)
          DO k = ks, ke - 1
             j = nxleaf
             nxleaf = Leaf(k+1)
             DO WHILE ( j<nxleaf )
                IF ( j>=Nstop ) GOTO 100
                l = l + 1
                j = Parent(j)
             ENDDO
          ENDDO
          j = Leaf(ke)
          DO WHILE ( j<Nstop )
             IF ( j>=i .OR. j==0 ) EXIT
             l = l + 1
             j = Parent(j)
          ENDDO
       ENDIF
100 ENDDO
    Xlnzr(Neqns+1) = l
    Lncol = l - 1
  END SUBROUTINE PRE_GNCLNO

  !======================================================================!
  !> @brief PTIME
  !======================================================================!
  SUBROUTINE PTIME(Cputim)
    USE HECMW_UTIL
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(OUT):: Cputim

    ! cpu time by hour
    Cputim = HECMW_WTIME()
  END SUBROUTINE PTIME

  !======================================================================!
  !> @brief QMDMRG
  !======================================================================!
  SUBROUTINE QMDMRG(Xadj,Adjncy,Deg,Qsize,Qlink,Marker,Deg0,Nhdsze,Nbrhd,Rchset,Ovrlp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Deg0
    INTEGER, INTENT(IN):: Nhdsze
    INTEGER, INTENT(IN):: Adjncy(*)
    INTEGER, INTENT(IN):: Nbrhd(*)
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(OUT):: Deg(*)
    INTEGER, INTENT(OUT):: Qsize(*)
    INTEGER, INTENT(OUT):: Qlink(*)
    INTEGER, INTENT(OUT):: Marker(*)
    INTEGER, INTENT(OUT):: Rchset(*)
    INTEGER, INTENT(OUT):: Ovrlp(*)

    INTEGER:: deg1
    INTEGER:: head
    INTEGER:: inhd
    INTEGER:: iov
    INTEGER:: irch
    INTEGER:: j
    INTEGER:: jstrt
    INTEGER:: jstop
    INTEGER:: link
    INTEGER:: lnode
    INTEGER:: mark
    INTEGER:: mrgsze
    INTEGER:: nabor
    INTEGER:: node
    INTEGER:: novrlp
    INTEGER:: rchsze
    INTEGER:: root

    IF ( Nhdsze<=0 ) RETURN
    DO inhd = 1, Nhdsze
       root = Nbrhd(inhd)
       Marker(root) = 0
    ENDDO
    DO inhd = 1, Nhdsze
       root = Nbrhd(inhd)
       Marker(root) = -1
       rchsze = 0
       novrlp = 0
       deg1 = 0
       DO
          jstrt = Xadj(root)
          jstop = Xadj(root+1) - 1
          DO j = jstrt, jstop
             nabor = Adjncy(j)
             root = -nabor
             IF ( nabor<0 ) GOTO 50
             IF ( nabor==0 ) EXIT
             mark = Marker(nabor)

             IF ( mark>=0 ) THEN
                IF ( mark<=0 ) THEN
                   rchsze = rchsze + 1
                   Rchset(rchsze) = nabor
                   deg1 = deg1 + Qsize(nabor)
                   Marker(nabor) = 1
                ELSEIF ( mark<=1 ) THEN
                   novrlp = novrlp + 1
                   Ovrlp(novrlp) = nabor
                   Marker(nabor) = 2
                ENDIF
             ENDIF
          ENDDO
          EXIT
50     ENDDO
       head = 0
       mrgsze = 0
       DO iov = 1, novrlp
          node = Ovrlp(iov)
          jstrt = Xadj(node)
          jstop = Xadj(node+1) - 1
          DO j = jstrt, jstop
             nabor = Adjncy(j)
             IF ( Marker(nabor)==0 ) THEN
                Marker(node) = 1
                GOTO 100
             ENDIF
          ENDDO
          mrgsze = mrgsze + Qsize(node)
          Marker(node) = -1
          lnode = node
          DO
             link = Qlink(lnode)
             IF ( link<=0 ) THEN
                Qlink(lnode) = head
                head = node
                EXIT
             ELSE
                lnode = link
             ENDIF
          ENDDO
100    ENDDO
       IF ( head>0 ) THEN
          Qsize(head) = mrgsze
          Deg(head) = Deg0 + deg1 - 1
          Marker(head) = 2
       ENDIF
       root = Nbrhd(inhd)
       Marker(root) = 0
       IF ( rchsze>0 ) THEN
          DO irch = 1, rchsze
             node = Rchset(irch)
             Marker(node) = 0
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE QMDMRG

  !======================================================================!
  !> @brief QMDOT
  !======================================================================!
  SUBROUTINE QMDOT(Root,Xadj,Adjncy,Marker,Rchsze,Rchset,Nbrhd)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Rchsze
    INTEGER, INTENT(IN):: Root
    INTEGER, INTENT(IN):: Marker(*)
    INTEGER, INTENT(IN):: Rchset(*)
    INTEGER, INTENT(IN):: Nbrhd(*)
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(OUT):: Adjncy(*)

    INTEGER:: inhd
    INTEGER:: irch
    INTEGER:: j
    INTEGER:: jstrt
    INTEGER:: jstop
    INTEGER:: link
    INTEGER:: nabor
    INTEGER:: node

    irch = 0
    inhd = 0
    node = Root
100 jstrt = Xadj(node)
    jstop = Xadj(node+1) - 2
    IF ( jstop>=jstrt ) THEN
       DO j = jstrt, jstop
          irch = irch + 1
          Adjncy(j) = Rchset(irch)
          IF ( irch>=Rchsze ) GOTO 200
       ENDDO
    ENDIF
    link = Adjncy(jstop+1)
    node = -link
    IF ( link>=0 ) THEN
       inhd = inhd + 1
       node = Nbrhd(inhd)
       Adjncy(jstop+1) = -node
    ENDIF
    GOTO 100
200 Adjncy(j+1) = 0
    DO irch = 1, Rchsze
       node = Rchset(irch)
       IF ( Marker(node)>=0 ) THEN
          jstrt = Xadj(node)
          jstop = Xadj(node+1) - 1
          DO j = jstrt, jstop
             nabor = Adjncy(j)
             IF ( Marker(nabor)<0 ) THEN
                Adjncy(j) = Root
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE QMDOT

  !======================================================================!
  !> @brief QMDRCH
  !======================================================================!
  SUBROUTINE QMDRCH(Root,Xadj,Adjncy,Deg,Marker,Rchsze,Rchset,Nhdsze,Nbrhd)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Root
    INTEGER, INTENT(IN):: Adjncy(*)
    INTEGER, INTENT(IN):: Deg(*)
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(OUT):: Nhdsze
    INTEGER, INTENT(OUT):: Rchsze
    INTEGER, INTENT(OUT):: Marker(*)
    INTEGER, INTENT(OUT):: Rchset(*)
    INTEGER, INTENT(OUT):: Nbrhd(*)

    INTEGER:: i
    INTEGER:: istrt
    INTEGER:: istop
    INTEGER:: j
    INTEGER:: jstrt
    INTEGER:: jstop
    INTEGER:: nabor
    INTEGER:: node

    Nhdsze = 0
    Rchsze = 0
    istrt = Xadj(Root)
    istop = Xadj(Root+1) - 1
    IF ( istop<istrt ) RETURN
    DO i = istrt, istop
       nabor = Adjncy(i)
       IF ( nabor==0 ) RETURN
       IF ( Marker(nabor)==0 ) THEN
          IF ( Deg(nabor)<0 ) THEN
             Marker(nabor) = -1
             Nhdsze = Nhdsze + 1
             Nbrhd(Nhdsze) = nabor
             DO
                jstrt = Xadj(nabor)
                jstop = Xadj(nabor+1) - 1
                DO j = jstrt, jstop
                   node = Adjncy(j)
                   nabor = -node
                   IF ( node<0 ) GOTO 10
                   IF ( node==0 ) EXIT
                   IF ( Marker(node)==0 ) THEN
                      Rchsze = Rchsze + 1
                      Rchset(Rchsze) = node
                      Marker(node) = 1
                   ENDIF
                ENDDO
                EXIT
10           ENDDO
          ELSE
             Rchsze = Rchsze + 1
             Rchset(Rchsze) = nabor
             Marker(nabor) = 1
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE QMDRCH

  !======================================================================!
  !> @brief QMDUPD
  !======================================================================!
  SUBROUTINE QMDUPD(Xadj,Adjncy,Nlist,List,Deg,Qsize,Qlink,Marker,Rchset,Nbrhd)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Nlist
    INTEGER, INTENT(IN):: Adjncy(*)
    INTEGER, INTENT(IN):: List(*)
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(OUT):: Deg(*)
    INTEGER, INTENT(OUT):: Marker(*)
    INTEGER, INTENT(OUT):: Rchset(*)
    INTEGER, INTENT(OUT):: Nbrhd(*)
    INTEGER, INTENT(OUT):: Qsize(*)
    INTEGER, INTENT(OUT):: Qlink(*)

    INTEGER:: deg0
    INTEGER:: deg1
    INTEGER:: il
    INTEGER:: inhd
    INTEGER:: inode
    INTEGER:: irch
    INTEGER:: j
    INTEGER:: jstrt
    INTEGER:: jstop
    INTEGER:: mark
    INTEGER:: nabor
    INTEGER:: nhdsze
    INTEGER:: node
    INTEGER:: rchsze

    IF ( Nlist<=0 ) RETURN
    deg0 = 0
    nhdsze = 0
    DO il = 1, Nlist
       node = List(il)
       deg0 = deg0 + Qsize(node)
       jstrt = Xadj(node)
       jstop = Xadj(node+1) - 1
       DO j = jstrt, jstop
          nabor = Adjncy(j)
          IF ( Marker(nabor)==0 .AND. Deg(nabor)<0 ) THEN
             Marker(nabor) = -1
             nhdsze = nhdsze + 1
             Nbrhd(nhdsze) = nabor
          ENDIF
       ENDDO
    ENDDO

    IF ( nhdsze>0 ) CALL QMDMRG(Xadj,Adjncy,Deg,Qsize,Qlink,Marker,deg0,nhdsze,Nbrhd,Rchset,Nbrhd(nhdsze+1))
    DO il = 1, Nlist
       node = List(il)
       mark = Marker(node)
       IF ( mark<=1 .AND. mark>=0 ) THEN
          CALL QMDRCH(node,Xadj,Adjncy,Deg,Marker,rchsze,Rchset,nhdsze,Nbrhd)
          deg1 = deg0
          IF ( rchsze>0 ) THEN
             DO irch = 1, rchsze
                inode = Rchset(irch)
                deg1 = deg1 + Qsize(inode)
                Marker(inode) = 0
             ENDDO
          ENDIF
          Deg(node) = deg1 - 1
          IF ( nhdsze>0 ) THEN
             DO inhd = 1, nhdsze
                inode = Nbrhd(inhd)
                Marker(inode) = 0
             ENDDO
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE QMDUPD

  !======================================================================!
  !> @brief QQSORT
  ! sort in increasing order up to i
  !     iw   array
  !     ik   number of input/output
  !     i    deal with numbers less than this numberi
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE QQSORT(Iw,Ik)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ik
    INTEGER, INTENT(OUT):: Iw(*)

    INTEGER:: IDBg1
    INTEGER:: itemp
    INTEGER:: l
    INTEGER:: m
    COMMON /DEBUG / IDBg1

    IF ( Ik<=1 ) RETURN
    DO l = 1, Ik - 1
       DO m = l + 1, Ik
          IF ( Iw(l)>=Iw(m) ) THEN
             itemp = Iw(l)
             Iw(l) = Iw(m)
             Iw(m) = itemp
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE QQSORT

  !======================================================================!
  !> @brief ROTATE
  ! irr return code irr=0 node izz is not a bottom node
  !                     irr=1          is a bottom node then rotation is
  !                                    performed
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE ROTATE(Xadj,Adjncy,Invp,Iperm,Parent,Btree,Izz,Neqns,Anc,Adjt,Irr)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Izz
    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Xadj(*)
    INTEGER, INTENT(IN):: Adjncy(*)
    INTEGER, INTENT(IN):: Parent(*)
    INTEGER, INTENT(IN):: Btree(2,*)
    INTEGER, INTENT(OUT):: Irr
    INTEGER, INTENT(OUT):: Invp(*)
    INTEGER, INTENT(OUT):: Iperm(*)
    INTEGER, INTENT(OUT):: Anc(*)
    INTEGER, INTENT(OUT):: Adjt(*)

    INTEGER:: i
    INTEGER:: IDBg1
    INTEGER:: iy
    INTEGER:: izzz
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: kk
    INTEGER:: l
    INTEGER:: ll
    INTEGER:: locc
    INTEGER:: nanc
    COMMON /DEBUG / IDBg1

    IF ( Izz==0 ) THEN
       Irr = 0
       RETURN
    ENDIF
    izzz = Invp(Izz)

    IF ( Btree(1,izzz)/=0 ) Irr = 0
    Irr = 1
    !
    !  ancestors of izzz
    !
    nanc = 0
    joc = izzz
    DO
       nanc = nanc + 1
       Anc(nanc) = joc
       joc = Parent(joc)
       IF ( joc==0 ) THEN
          !
          !  to find the eligible node from ancestors of izz
          !
          l = 1
          EXIT
       ENDIF
    ENDDO
100 DO i = 1, Neqns
       Adjt(i) = 0
    ENDDO
    locc = Anc(l)
    DO
       joc = locc
       locc = Btree(1,joc)
       IF ( locc==0 ) THEN
          DO
             DO k = Xadj(Iperm(joc)), Xadj(Iperm(joc)+1) - 1
                Adjt(Invp(Adjncy(k))) = 1
             ENDDO
             IF ( joc>=Anc(l) ) THEN
                DO ll = l + 1, nanc
                   IF ( Adjt(Anc(ll))==0 ) THEN
                      l = l + 1
                      GOTO 100
                   ENDIF
                ENDDO
                IF ( l==1 ) THEN
                   !
                   ! izz can be numbered last
                   !
                   k = 0
                   DO i = 1, Neqns
                      IF ( i/=izzz ) THEN
                         k = k + 1
                         Invp(Iperm(i)) = k
                      ENDIF
                   ENDDO
                   Invp(Iperm(izzz)) = Neqns
                ELSE
                   !
                   !  anc(l-1) is the eligible node
                   !
                   ! (1) number the node not in Ancestor(iy)
                   iy = Anc(l-1)
                   DO i = 1, Neqns
                      Adjt(i) = 0
                   ENDDO
                   DO ll = l, nanc
                      Adjt(Anc(ll)) = 1
                   ENDDO
                   k = 0
                   DO ll = 1, Neqns
                      IF ( Adjt(ll)==0 ) THEN
                         k = k + 1
                         Invp(Iperm(ll)) = k
                      ENDIF
                   ENDDO
                   ! (2) followed by nodes in Ancestor(iy)-Adj(T(iy))
                   DO i = 1, Neqns
                      Adjt(i) = 0
                   ENDDO
                   locc = iy
                   DO
                      joc = locc
                      locc = Btree(1,joc)
                      IF ( locc==0 ) THEN
                         DO
                            DO kk = Xadj(Iperm(joc)), Xadj(Iperm(joc)+1) - 1
                               Adjt(Invp(Adjncy(kk))) = 1
                            ENDDO
                            IF ( joc>=iy ) THEN
                               DO ll = l, nanc
                                  IF ( Adjt(Anc(ll))==0 ) THEN
                                     k = k + 1
                                     Invp(Iperm(Anc(ll))) = k
                                  ENDIF
                               ENDDO
                               ! (3) and finally number the node in Adj(t(iy))
                               DO ll = l, nanc
                                  IF ( Adjt(Anc(ll))/=0 ) THEN
                                     k = k + 1
                                     Invp(Iperm(Anc(ll))) = k
                                  ENDIF
                               ENDDO
                               GOTO 105
                            ELSE
                               locc = Btree(2,joc)
                               IF ( locc/=0 ) EXIT
                               joc = Parent(joc)
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
                !
                ! set iperm
                !
105             DO i = 1, Neqns
                   Iperm(Invp(i)) = i
                ENDDO
                IF ( IDBg1/=0 ) WRITE (6,"(10I6)") (Invp(i),i=1,Neqns)
                RETURN
             ELSE
                locc = Btree(2,joc)
                IF ( locc/=0 ) EXIT
                joc = Parent(joc)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE ROTATE

  !======================================================================!
  !> @brief S2UM
  !======================================================================!
  SUBROUTINE S2UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Par(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(3,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(4,*)

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: s(4)
    DOUBLE PRECISION:: t(3)
    DOUBLE PRECISION:: zz(4)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t(1) = 0.0D0
    t(2) = 0.0D0
    t(3) = 0.0D0
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       DO l = 1, 4
          s(l) = 0.0D0
          zz(l) = Zln(l,k)
       ENDDO
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             zz(1) = zz(1) - Temp(1,j)*Zln(1,jj) - Temp(3,j)*Zln(3,jj)
             zz(2) = zz(2) - Temp(2,j)*Zln(1,jj) - Temp(4,j)*Zln(3,jj)
             zz(3) = zz(3) - Temp(1,j)*Zln(2,jj) - Temp(3,j)*Zln(4,jj)
             zz(4) = zz(4) - Temp(2,j)*Zln(2,jj) - Temp(4,j)*Zln(4,jj)
          ENDIF
       ENDDO
       CALL INV22(Zln(1,k),zz,Diag(1,jc))
       DO l = 1, 4
          Temp(l,jc) = zz(l)
       ENDDO
       t(1) = t(1) + zz(1)*Zln(1,k) + zz(3)*Zln(3,k)
       t(2) = t(2) + zz(1)*Zln(2,k) + zz(3)*Zln(4,k)
       t(3) = t(3) + zz(2)*Zln(2,k) + zz(4)*Zln(4,k)
    ENDDO
    Diag(1,Ic) = Diag(1,Ic) - t(1)
    Diag(2,Ic) = Diag(2,Ic) - t(2)
    Diag(3,Ic) = Diag(3,Ic) - t(3)
    CALL INV2(Diag(1,Ic),ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  END SUBROUTINE S2UM

  !======================================================================!
  !> @brief S2UM1
  !======================================================================!
  SUBROUTINE S2UM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(4,*)

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: s(4)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    DO l = 1, 4
       s(l) = 0.0D0
    ENDDO
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             s(1) = s(1) + Temp(1,j)*Zln(1,jj) + Temp(3,j)*Zln(3,jj)
             s(2) = s(2) + Temp(2,j)*Zln(1,jj) + Temp(4,j)*Zln(3,jj)
             s(3) = s(3) + Temp(1,j)*Zln(2,jj) + Temp(3,j)*Zln(4,jj)
             s(4) = s(4) + Temp(2,j)*Zln(2,jj) + Temp(4,j)*Zln(4,jj)
          ENDIF
       ENDDO
       DO l = 1, 4
          Temp(l,jc) = Zln(l,k) - s(l)
          Zln(l,k) = Temp(l,jc)
          s(l) = 0.0D0
       ENDDO
    ENDDO
  END SUBROUTINE S2UM1

  !======================================================================!
  !> @brief S2UM2
  !======================================================================!
  SUBROUTINE S2UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(3,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(4,*)

    INTEGER:: ic
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    joc = 0
    DO ic = Nstop, Neqns
       ks = Xlnzr(ic)
       ke = Xlnzr(ic+1) - 1
       DO k = ks, ke
          jj = Colno(k)
          Temp(1,jj) = Zln(1,k)
          Temp(2,jj) = Zln(2,k)
          Temp(3,jj) = Zln(3,k)
          Temp(4,jj) = Zln(4,k)

          Zln(3,k) = Temp(3,jj) - Temp(1,jj)*Diag(2,jj)
          Zln(1,k) = Temp(1,jj)*Diag(1,jj)
          Zln(3,k) = Zln(3,k)*Diag(3,jj)
          Zln(1,k) = Zln(1,k) - Zln(3,k)*Diag(2,jj)

          Zln(4,k) = Temp(4,jj) - Temp(2,jj)*Diag(2,jj)
          Zln(2,k) = Temp(2,jj)*Diag(1,jj)
          Zln(4,k) = Zln(4,k)*Diag(3,jj)
          Zln(2,k) = Zln(2,k) - Zln(4,k)*Diag(2,jj)

          Diag(1,ic) = Diag(1,ic) - (Temp(1,jj)*Zln(1,k)+Temp(3,jj)*Zln(3,k))
          Diag(2,ic) = Diag(2,ic) - (Temp(1,jj)*Zln(2,k)+Temp(3,jj)*Zln(4,k))
          Diag(3,ic) = Diag(3,ic) - (Temp(2,jj)*Zln(2,k)+Temp(4,jj)*Zln(4,k))
          Indx(jj) = ic
       ENDDO
       DO jc = Nstop, ic - 1
          joc = joc + 1
          DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
             j = Colno(jj)
             IF ( Indx(j)==ic ) THEN
                Dsln(1,joc) = Dsln(1,joc) - (Temp(1,j)*Zln(1,jj)+Temp(3,j)*Zln(3,jj))
                Dsln(2,joc) = Dsln(2,joc) - (Temp(2,j)*Zln(1,jj)+Temp(4,j)*Zln(3,jj))
                Dsln(3,joc) = Dsln(3,joc) - (Temp(1,j)*Zln(2,jj)+Temp(3,j)*Zln(4,jj))
                Dsln(4,joc) = Dsln(4,joc) - (Temp(2,j)*Zln(2,jj)+Temp(4,j)*Zln(4,jj))
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE S2UM2

  !======================================================================!
  !> @brief S2UM3
  !======================================================================!
  SUBROUTINE S2UM3(N,Dsln,Diag,Indx,Temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(3,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(4,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(4,*)

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    DOUBLE PRECISION:: t(4)

    IF ( N>0 ) THEN
       Indx(1) = 0
       joc = 1
       CALL INV2(Diag(1,1),ir)
       DO i = 2, N
          Indx(i) = joc
          DO j = 1, i - 1
             CALL D2DOT(t,Dsln(1,Indx(i)),Dsln(1,Indx(j)),j-1)
             Dsln(1,joc) = Dsln(1,joc) - t(1)
             Dsln(2,joc) = Dsln(2,joc) - t(2)
             Dsln(3,joc) = Dsln(3,joc) - t(3)
             Dsln(4,joc) = Dsln(4,joc) - t(4)
             joc = joc + 1
          ENDDO
          CALL V2PROD(Dsln(1,Indx(i)),Diag,Temp,i-1)
          CALL D2DOT(t,Temp,Dsln(1,Indx(i)),i-1)
          Diag(1,i) = Diag(1,i) - t(1)
          Diag(2,i) = Diag(2,i) - t(2)
          Diag(3,i) = Diag(3,i) - t(4)
          CALL VCOPY(Temp,Dsln(1,Indx(i)),4*(i-1))
          CALL INV2(Diag(1,i),ir)
       ENDDO
    ENDIF
  END SUBROUTINE S2UM3

  !======================================================================!
  !> @brief S3UM
  !======================================================================!
  SUBROUTINE S3UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Par(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(9,*)

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t(6)
    DOUBLE PRECISION:: zz(9)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    !$dir max_trips(6)
    DO l = 1, 6
       t(l) = 0.0D0
    ENDDO
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       !$dir max_trips(9)
       DO l = 1, 9
          zz(l) = Zln(l,k)
       ENDDO
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             zz(1) = zz(1) - Temp(1,j)*Zln(1,jj) - Temp(4,j)*Zln(4,jj) - Temp(7,j)*Zln(7,jj)
             zz(2) = zz(2) - Temp(2,j)*Zln(1,jj) - Temp(5,j)*Zln(4,jj) - Temp(8,j)*Zln(7,jj)
             zz(3) = zz(3) - Temp(3,j)*Zln(1,jj) - Temp(6,j)*Zln(4,jj) - Temp(9,j)*Zln(7,jj)
             zz(4) = zz(4) - Temp(1,j)*Zln(2,jj) - Temp(4,j)*Zln(5,jj) - Temp(7,j)*Zln(8,jj)
             zz(5) = zz(5) - Temp(2,j)*Zln(2,jj) - Temp(5,j)*Zln(5,jj) - Temp(8,j)*Zln(8,jj)
             zz(6) = zz(6) - Temp(3,j)*Zln(2,jj) - Temp(6,j)*Zln(5,jj) - Temp(9,j)*Zln(8,jj)
             zz(7) = zz(7) - Temp(1,j)*Zln(3,jj) - Temp(4,j)*Zln(6,jj) - Temp(7,j)*Zln(9,jj)
             zz(8) = zz(8) - Temp(2,j)*Zln(3,jj) - Temp(5,j)*Zln(6,jj) - Temp(8,j)*Zln(9,jj)
             zz(9) = zz(9) - Temp(3,j)*Zln(3,jj) - Temp(6,j)*Zln(6,jj) - Temp(9,j)*Zln(9,jj)
          ENDIF
       ENDDO
       CALL INV33(Zln(1,k),zz,Diag(1,jc))
       !$dir max_trips(9)
       DO l = 1, 9
          Temp(l,jc) = zz(l)
       ENDDO
       t(1) = t(1) + zz(1)*Zln(1,k) + zz(4)*Zln(4,k) + zz(7)*Zln(7,k)
       t(2) = t(2) + zz(1)*Zln(2,k) + zz(4)*Zln(5,k) + zz(7)*Zln(8,k)
       t(3) = t(3) + zz(2)*Zln(2,k) + zz(5)*Zln(5,k) + zz(8)*Zln(8,k)
       t(4) = t(4) + zz(1)*Zln(3,k) + zz(4)*Zln(6,k) + zz(7)*Zln(9,k)
       t(5) = t(5) + zz(2)*Zln(3,k) + zz(5)*Zln(6,k) + zz(8)*Zln(9,k)
       t(6) = t(6) + zz(3)*Zln(3,k) + zz(6)*Zln(6,k) + zz(9)*Zln(9,k)
    ENDDO
    !$dir max_trips(6)
    DO l = 1, 6
       Diag(l,Ic) = Diag(l,Ic) - t(l)
    ENDDO
    CALL INV3(Diag(1,Ic),ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  END SUBROUTINE S3UM

  !======================================================================!
  !> @brief S3UM1
  !======================================================================!
  SUBROUTINE S3UM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(9,*)

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: s(9)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    !$dir max_trip(9)
    DO l = 1, 9
       s(l) = 0.0D0
    ENDDO
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             s(1) = s(1) + Temp(1,j)*Zln(1,jj) + Temp(4,j)*Zln(4,jj) + Temp(7,j)*Zln(7,jj)
             s(2) = s(2) + Temp(2,j)*Zln(1,jj) + Temp(5,j)*Zln(4,jj) + Temp(8,j)*Zln(7,jj)
             s(3) = s(3) + Temp(3,j)*Zln(1,jj) + Temp(6,j)*Zln(4,jj) + Temp(9,j)*Zln(7,jj)
             s(4) = s(4) + Temp(1,j)*Zln(2,jj) + Temp(4,j)*Zln(5,jj) + Temp(7,j)*Zln(8,jj)
             s(5) = s(5) + Temp(2,j)*Zln(2,jj) + Temp(5,j)*Zln(5,jj) + Temp(8,j)*Zln(8,jj)
             s(6) = s(6) + Temp(3,j)*Zln(2,jj) + Temp(6,j)*Zln(5,jj) + Temp(9,j)*Zln(8,jj)
             s(7) = s(7) + Temp(1,j)*Zln(3,jj) + Temp(4,j)*Zln(6,jj) + Temp(7,j)*Zln(9,jj)
             s(8) = s(8) + Temp(2,j)*Zln(3,jj) + Temp(5,j)*Zln(6,jj) + Temp(8,j)*Zln(9,jj)
             s(9) = s(9) + Temp(3,j)*Zln(3,jj) + Temp(6,j)*Zln(6,jj) + Temp(9,j)*Zln(9,jj)
          ENDIF
       ENDDO
       !$dir max_trip(9)
       DO l = 1, 9
          Temp(l,jc) = Zln(l,k) - s(l)
          Zln(l,k) = Temp(l,jc)
          s(l) = 0.0D0
       ENDDO
    ENDDO
  END SUBROUTINE S3UM1

  !======================================================================!
  !> @brief S3UM1
  !======================================================================!
  SUBROUTINE S3UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(Neqns,9)
    DOUBLE PRECISION, INTENT(OUT):: Zln(9,*)

    INTEGER:: ic
    INTEGER:: j
    INTEGER:: j1
    INTEGER:: j2
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    joc = 0
    DO ic = Nstop, Neqns
       ks = Xlnzr(ic)
       ke = Xlnzr(ic+1) - 1
       DO k = ks, ke
          jj = Colno(k)
          Temp(jj,1) = Zln(1,k)
          Temp(jj,2) = Zln(2,k)
          Temp(jj,3) = Zln(3,k)
          Temp(jj,4) = Zln(4,k)
          Temp(jj,5) = Zln(5,k)
          Temp(jj,6) = Zln(6,k)
          Temp(jj,7) = Zln(7,k)
          Temp(jj,8) = Zln(8,k)
          Temp(jj,9) = Zln(9,k)
          Indx(jj) = ic
       ENDDO
       DO k = ks, ke
          jj = Colno(k)
          Zln(4,k) = Temp(jj,4) - Temp(jj,1)*Diag(2,jj)
          Zln(7,k) = Temp(jj,7) - Temp(jj,1)*Diag(4,jj) - Zln(4,k)*Diag(5,jj)
          Zln(1,k) = Temp(jj,1)*Diag(1,jj)
          Zln(4,k) = Zln(4,k)*Diag(3,jj)
          Zln(7,k) = Zln(7,k)*Diag(6,jj)
          Zln(4,k) = Zln(4,k) - Zln(7,k)*Diag(5,jj)
          Zln(1,k) = Zln(1,k) - Zln(4,k)*Diag(2,jj) - Zln(7,k)*Diag(4,jj)

          Zln(5,k) = Temp(jj,5) - Temp(jj,2)*Diag(2,jj)
          Zln(8,k) = Temp(jj,8) - Temp(jj,2)*Diag(4,jj) - Zln(5,k)*Diag(5,jj)
          Zln(2,k) = Temp(jj,2)*Diag(1,jj)
          Zln(5,k) = Zln(5,k)*Diag(3,jj)
          Zln(8,k) = Zln(8,k)*Diag(6,jj)
          Zln(5,k) = Zln(5,k) - Zln(8,k)*Diag(5,jj)
          Zln(2,k) = Zln(2,k) - Zln(5,k)*Diag(2,jj) - Zln(8,k)*Diag(4,jj)

          Zln(6,k) = Temp(jj,6) - Temp(jj,3)*Diag(2,jj)
          Zln(9,k) = Temp(jj,9) - Temp(jj,3)*Diag(4,jj) - Zln(6,k)*Diag(5,jj)
          Zln(3,k) = Temp(jj,3)*Diag(1,jj)
          Zln(6,k) = Zln(6,k)*Diag(3,jj)
          Zln(9,k) = Zln(9,k)*Diag(6,jj)
          Zln(6,k) = Zln(6,k) - Zln(9,k)*Diag(5,jj)
          Zln(3,k) = Zln(3,k) - Zln(6,k)*Diag(2,jj) - Zln(9,k)*Diag(4,jj)
       ENDDO

       DO k = ks, ke
          jj = Colno(k)
          Diag(1,ic) = Diag(1,ic) - Temp(jj,1)*Zln(1,k) - Temp(jj,4)*Zln(4,k) - Temp(jj,7)*Zln(7,k)
          Diag(2,ic) = Diag(2,ic) - Temp(jj,1)*Zln(2,k) - Temp(jj,4)*Zln(5,k) - Temp(jj,7)*Zln(8,k)
          Diag(3,ic) = Diag(3,ic) - Temp(jj,2)*Zln(2,k) - Temp(jj,5)*Zln(5,k) - Temp(jj,8)*Zln(8,k)
          Diag(4,ic) = Diag(4,ic) - Temp(jj,1)*Zln(3,k) - Temp(jj,4)*Zln(6,k) - Temp(jj,7)*Zln(9,k)
          Diag(5,ic) = Diag(5,ic) - Temp(jj,2)*Zln(3,k) - Temp(jj,5)*Zln(6,k) - Temp(jj,8)*Zln(9,k)
          Diag(6,ic) = Diag(6,ic) - Temp(jj,3)*Zln(3,k) - Temp(jj,6)*Zln(6,k) - Temp(jj,9)*Zln(9,k)
       ENDDO
       DO jc = Nstop, ic - 1
          joc = joc + 1
          j1 = Xlnzr(jc)
          j2 = Xlnzr(jc+1)
          DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
             j = Colno(jj)
             IF ( Indx(j)==ic ) THEN
                Dsln(1,joc) = Dsln(1,joc) - Temp(j,1)*Zln(1,jj) - Temp(j,4)*Zln(4,jj) - Temp(j,7)*Zln(7,jj)
                Dsln(2,joc) = Dsln(2,joc) - Temp(j,2)*Zln(1,jj) - Temp(j,5)*Zln(4,jj) - Temp(j,8)*Zln(7,jj)
                Dsln(3,joc) = Dsln(3,joc) - Temp(j,3)*Zln(1,jj) - Temp(j,6)*Zln(4,jj) - Temp(j,9)*Zln(7,jj)
                Dsln(4,joc) = Dsln(4,joc) - Temp(j,1)*Zln(2,jj) - Temp(j,4)*Zln(5,jj) - Temp(j,7)*Zln(8,jj)
                Dsln(5,joc) = Dsln(5,joc) - Temp(j,2)*Zln(2,jj) - Temp(j,5)*Zln(5,jj) - Temp(j,8)*Zln(8,jj)
                Dsln(6,joc) = Dsln(6,joc) - Temp(j,3)*Zln(2,jj) - Temp(j,6)*Zln(5,jj) - Temp(j,9)*Zln(8,jj)
                Dsln(7,joc) = Dsln(7,joc) - Temp(j,1)*Zln(3,jj) - Temp(j,4)*Zln(6,jj) - Temp(j,7)*Zln(9,jj)
                Dsln(8,joc) = Dsln(8,joc) - Temp(j,2)*Zln(3,jj) - Temp(j,5)*Zln(6,jj) - Temp(j,8)*Zln(9,jj)
                Dsln(9,joc) = Dsln(9,joc) - Temp(j,3)*Zln(3,jj) - Temp(j,6)*Zln(6,jj) - Temp(j,9)*Zln(9,jj)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE S3UM2

  !======================================================================!
  !> @brief S3UM3
  !======================================================================!
  SUBROUTINE S3UM3(N,Dsln,Diag,Indx,Temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(9,*)

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    DOUBLE PRECISION:: t(9)

    IF ( N>0 ) THEN
       Indx(1) = 0
       joc = 1
       CALL INV3(Diag(1,1),ir)
       DO i = 2, N
          Indx(i) = joc
          DO j = 1, i - 1
             CALL D3DOT(t,Dsln(1,Indx(i)),Dsln(1,Indx(j)),j-1)
             Dsln(:,joc) = Dsln(:,joc) - t(:)
             joc = joc + 1
          ENDDO
          CALL V3PROD(Dsln(1,Indx(i)),Diag,Temp,i-1)
          CALL D3DOTL(t,Temp,Dsln(1,Indx(i)),i-1)
          Diag(:,i) = Diag(:,i) - t(1:6)
          CALL VCOPY(Temp,Dsln(1,Indx(i)),9*(i-1))
          CALL INV3(Diag(1,i),ir)
       ENDDO
    ENDIF
  END SUBROUTINE S3UM3

  !======================================================================!
  !> @brief S6UM
  !======================================================================!
  SUBROUTINE S6UM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Par(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(21,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(36,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(36,*)

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: t(21)
    DOUBLE PRECISION:: zz(36)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    DO l = 1, 21
       t(l) = 0.0D0
    ENDDO
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       DO l = 1, 36
          zz(l) = Zln(l,k)
       ENDDO
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             zz(1) = zz(1) - Temp(1,j)*Zln(1,jj) - Temp(7,j)*Zln(7,jj)&
                  - Temp(13,j)*Zln(13,jj) - Temp(19,j)*Zln(19,jj)&
                  - Temp(25,j)*Zln(25,jj) - Temp(31,j)*Zln(31,jj)
             zz(2) = zz(2) - Temp(2,j)*Zln(1,jj) - Temp(8,j)*Zln(7,jj)&
                  - Temp(14,j)*Zln(13,jj) - Temp(20,j)*Zln(19,jj)&
                  - Temp(26,j)*Zln(25,jj) - Temp(32,j)*Zln(31,jj)
             zz(3) = zz(3) - Temp(3,j)*Zln(1,jj) - Temp(9,j)*Zln(7,jj)&
                  - Temp(15,j)*Zln(13,jj) - Temp(21,j)*Zln(19,jj)&
                  - Temp(27,j)*Zln(25,jj) - Temp(33,j)*Zln(31,jj)
             zz(4) = zz(4) - Temp(4,j)*Zln(1,jj) - Temp(10,j)&
                  *Zln(7,jj) - Temp(16,j)*Zln(13,jj) - Temp(22,j)&
                  *Zln(19,jj) - Temp(28,j)*Zln(25,jj) - Temp(34,j)&
                  *Zln(31,jj)
             zz(5) = zz(5) - Temp(5,j)*Zln(1,jj) - Temp(11,j)&
                  *Zln(7,jj) - Temp(17,j)*Zln(13,jj) - Temp(23,j)&
                  *Zln(19,jj) - Temp(29,j)*Zln(25,jj) - Temp(35,j)&
                  *Zln(31,jj)
             zz(6) = zz(6) - Temp(6,j)*Zln(1,jj) - Temp(12,j)&
                  *Zln(7,jj) - Temp(18,j)*Zln(13,jj) - Temp(24,j)&
                  *Zln(19,jj) - Temp(30,j)*Zln(25,jj) - Temp(36,j)&
                  *Zln(31,jj)
             zz(7) = zz(7) - Temp(1,j)*Zln(2,jj) - Temp(7,j)*Zln(8,jj)&
                  - Temp(13,j)*Zln(14,jj) - Temp(19,j)*Zln(20,jj)&
                  - Temp(25,j)*Zln(26,jj) - Temp(31,j)*Zln(32,jj)
             zz(8) = zz(8) - Temp(2,j)*Zln(2,jj) - Temp(8,j)*Zln(8,jj)&
                  - Temp(14,j)*Zln(14,jj) - Temp(20,j)*Zln(20,jj)&
                  - Temp(26,j)*Zln(26,jj) - Temp(32,j)*Zln(32,jj)
             zz(9) = zz(9) - Temp(3,j)*Zln(2,jj) - Temp(9,j)*Zln(8,jj)&
                  - Temp(15,j)*Zln(14,jj) - Temp(21,j)*Zln(20,jj)&
                  - Temp(27,j)*Zln(26,jj) - Temp(33,j)*Zln(32,jj)
             zz(10) = zz(10) - Temp(4,j)*Zln(2,jj) - Temp(10,j)&
                  *Zln(8,jj) - Temp(16,j)*Zln(14,jj) - Temp(22,j)&
                  *Zln(20,jj) - Temp(28,j)*Zln(26,jj) - Temp(34,j)&
                  *Zln(32,jj)
             zz(11) = zz(11) - Temp(5,j)*Zln(2,jj) - Temp(11,j)&
                  *Zln(8,jj) - Temp(17,j)*Zln(14,jj) - Temp(23,j)&
                  *Zln(20,jj) - Temp(29,j)*Zln(26,jj) - Temp(35,j)&
                  *Zln(32,jj)
             zz(12) = zz(12) - Temp(6,j)*Zln(2,jj) - Temp(12,j)&
                  *Zln(8,jj) - Temp(18,j)*Zln(14,jj) - Temp(24,j)&
                  *Zln(20,jj) - Temp(30,j)*Zln(26,jj) - Temp(36,j)&
                  *Zln(32,jj)
             zz(13) = zz(13) - Temp(1,j)*Zln(3,jj) - Temp(7,j)&
                  *Zln(9,jj) - Temp(13,j)*Zln(15,jj) - Temp(19,j)&
                  *Zln(21,jj) - Temp(25,j)*Zln(27,jj) - Temp(31,j)&
                  *Zln(33,jj)
             zz(14) = zz(14) - Temp(2,j)*Zln(3,jj) - Temp(8,j)&
                  *Zln(9,jj) - Temp(14,j)*Zln(15,jj) - Temp(20,j)&
                  *Zln(21,jj) - Temp(26,j)*Zln(27,jj) - Temp(32,j)&
                  *Zln(33,jj)
             zz(15) = zz(15) - Temp(3,j)*Zln(3,jj) - Temp(9,j)&
                  *Zln(9,jj) - Temp(15,j)*Zln(15,jj) - Temp(21,j)&
                  *Zln(21,jj) - Temp(27,j)*Zln(27,jj) - Temp(33,j)&
                  *Zln(33,jj)
             zz(16) = zz(16) - Temp(4,j)*Zln(3,jj) - Temp(10,j)&
                  *Zln(9,jj) - Temp(16,j)*Zln(15,jj) - Temp(22,j)&
                  *Zln(21,jj) - Temp(28,j)*Zln(27,jj) - Temp(34,j)&
                  *Zln(33,jj)
             zz(17) = zz(17) - Temp(5,j)*Zln(3,jj) - Temp(11,j)&
                  *Zln(9,jj) - Temp(17,j)*Zln(15,jj) - Temp(23,j)&
                  *Zln(21,jj) - Temp(29,j)*Zln(27,jj) - Temp(35,j)&
                  *Zln(33,jj)
             zz(18) = zz(18) - Temp(6,j)*Zln(3,jj) - Temp(12,j)&
                  *Zln(9,jj) - Temp(18,j)*Zln(15,jj) - Temp(24,j)&
                  *Zln(21,jj) - Temp(30,j)*Zln(27,jj) - Temp(36,j)&
                  *Zln(33,jj)
             zz(19) = zz(19) - Temp(1,j)*Zln(4,jj) - Temp(7,j)&
                  *Zln(10,jj) - Temp(13,j)*Zln(16,jj) - Temp(19,j)&
                  *Zln(22,jj) - Temp(25,j)*Zln(28,jj) - Temp(31,j)&
                  *Zln(34,jj)
             zz(20) = zz(20) - Temp(2,j)*Zln(4,jj) - Temp(8,j)&
                  *Zln(10,jj) - Temp(14,j)*Zln(16,jj) - Temp(20,j)&
                  *Zln(22,jj) - Temp(26,j)*Zln(28,jj) - Temp(32,j)&
                  *Zln(34,jj)
             zz(21) = zz(21) - Temp(3,j)*Zln(4,jj) - Temp(9,j)&
                  *Zln(10,jj) - Temp(15,j)*Zln(16,jj) - Temp(21,j)&
                  *Zln(22,jj) - Temp(27,j)*Zln(28,jj) - Temp(33,j)&
                  *Zln(34,jj)
             zz(22) = zz(22) - Temp(4,j)*Zln(4,jj) - Temp(10,j)&
                  *Zln(10,jj) - Temp(16,j)*Zln(16,jj) - Temp(22,j)&
                  *Zln(22,jj) - Temp(28,j)*Zln(28,jj) - Temp(34,j)&
                  *Zln(34,jj)
             zz(23) = zz(23) - Temp(5,j)*Zln(4,jj) - Temp(11,j)&
                  *Zln(10,jj) - Temp(17,j)*Zln(16,jj) - Temp(23,j)&
                  *Zln(22,jj) - Temp(29,j)*Zln(28,jj) - Temp(35,j)&
                  *Zln(34,jj)
             zz(24) = zz(24) - Temp(6,j)*Zln(4,jj) - Temp(12,j)&
                  *Zln(10,jj) - Temp(18,j)*Zln(16,jj) - Temp(24,j)&
                  *Zln(22,jj) - Temp(30,j)*Zln(28,jj) - Temp(36,j)&
                  *Zln(34,jj)
             zz(25) = zz(25) - Temp(1,j)*Zln(5,jj) - Temp(7,j)&
                  *Zln(11,jj) - Temp(13,j)*Zln(17,jj) - Temp(19,j)&
                  *Zln(23,jj) - Temp(25,j)*Zln(29,jj) - Temp(31,j)&
                  *Zln(35,jj)
             zz(26) = zz(26) - Temp(2,j)*Zln(5,jj) - Temp(8,j)&
                  *Zln(11,jj) - Temp(14,j)*Zln(17,jj) - Temp(20,j)&
                  *Zln(23,jj) - Temp(26,j)*Zln(29,jj) - Temp(32,j)&
                  *Zln(35,jj)
             zz(27) = zz(27) - Temp(3,j)*Zln(5,jj) - Temp(9,j)&
                  *Zln(11,jj) - Temp(15,j)*Zln(17,jj) - Temp(21,j)&
                  *Zln(23,jj) - Temp(27,j)*Zln(29,jj) - Temp(33,j)&
                  *Zln(35,jj)
             zz(28) = zz(28) - Temp(4,j)*Zln(5,jj) - Temp(10,j)&
                  *Zln(11,jj) - Temp(16,j)*Zln(17,jj) - Temp(22,j)&
                  *Zln(23,jj) - Temp(28,j)*Zln(29,jj) - Temp(34,j)&
                  *Zln(35,jj)
             zz(29) = zz(29) - Temp(5,j)*Zln(5,jj) - Temp(11,j)&
                  *Zln(11,jj) - Temp(17,j)*Zln(17,jj) - Temp(23,j)&
                  *Zln(23,jj) - Temp(29,j)*Zln(29,jj) - Temp(35,j)&
                  *Zln(35,jj)
             zz(30) = zz(30) - Temp(6,j)*Zln(5,jj) - Temp(12,j)&
                  *Zln(11,jj) - Temp(18,j)*Zln(17,jj) - Temp(24,j)&
                  *Zln(23,jj) - Temp(30,j)*Zln(29,jj) - Temp(36,j)&
                  *Zln(35,jj)
             zz(31) = zz(31) - Temp(1,j)*Zln(6,jj) - Temp(7,j)&
                  *Zln(12,jj) - Temp(13,j)*Zln(18,jj) - Temp(19,j)&
                  *Zln(24,jj) - Temp(25,j)*Zln(30,jj) - Temp(31,j)&
                  *Zln(36,jj)
             zz(32) = zz(32) - Temp(2,j)*Zln(6,jj) - Temp(8,j)&
                  *Zln(12,jj) - Temp(14,j)*Zln(18,jj) - Temp(20,j)&
                  *Zln(24,jj) - Temp(26,j)*Zln(30,jj) - Temp(32,j)&
                  *Zln(36,jj)
             zz(33) = zz(33) - Temp(3,j)*Zln(6,jj) - Temp(9,j)&
                  *Zln(12,jj) - Temp(15,j)*Zln(18,jj) - Temp(21,j)&
                  *Zln(24,jj) - Temp(27,j)*Zln(30,jj) - Temp(33,j)&
                  *Zln(36,jj)
             zz(34) = zz(34) - Temp(4,j)*Zln(6,jj) - Temp(10,j)&
                  *Zln(12,jj) - Temp(16,j)*Zln(18,jj) - Temp(22,j)&
                  *Zln(24,jj) - Temp(28,j)*Zln(30,jj) - Temp(34,j)&
                  *Zln(36,jj)
             zz(35) = zz(35) - Temp(5,j)*Zln(6,jj) - Temp(11,j)&
                  *Zln(12,jj) - Temp(17,j)*Zln(18,jj) - Temp(23,j)&
                  *Zln(24,jj) - Temp(29,j)*Zln(30,jj) - Temp(35,j)&
                  *Zln(36,jj)
             zz(36) = zz(36) - Temp(6,j)*Zln(6,jj) - Temp(12,j)&
                  *Zln(12,jj) - Temp(18,j)*Zln(18,jj) - Temp(24,j)&
                  *Zln(24,jj) - Temp(30,j)*Zln(30,jj) - Temp(36,j)&
                  *Zln(36,jj)
          ENDIF
       ENDDO
       CALL INV66(Zln(1,k),zz,Diag(1,jc))
       !$dir max_trips(9)
       DO l = 1, 36
          Temp(l,jc) = zz(l)
       ENDDO
       t(1) = t(1) + zz(1)*Zln(1,k) + zz(7)*Zln(7,k) + zz(13)&
            *Zln(13,k) + zz(19)*Zln(19,k) + zz(25)*Zln(25,k)&
            + zz(31)*Zln(31,k)
       t(2) = t(2) + zz(1)*Zln(2,k) + zz(7)*Zln(8,k) + zz(13)&
            *Zln(14,k) + zz(19)*Zln(20,k) + zz(25)*Zln(26,k)&
            + zz(31)*Zln(32,k)
       t(3) = t(3) + zz(2)*Zln(2,k) + zz(8)*Zln(8,k) + zz(14)&
            *Zln(14,k) + zz(20)*Zln(20,k) + zz(26)*Zln(26,k)&
            + zz(32)*Zln(32,k)
       t(4) = t(4) + zz(1)*Zln(3,k) + zz(7)*Zln(9,k) + zz(13)&
            *Zln(15,k) + zz(19)*Zln(21,k) + zz(25)*Zln(27,k)&
            + zz(31)*Zln(33,k)
       t(5) = t(5) + zz(2)*Zln(3,k) + zz(8)*Zln(9,k) + zz(14)&
            *Zln(15,k) + zz(20)*Zln(21,k) + zz(26)*Zln(27,k)&
            + zz(32)*Zln(33,k)
       t(6) = t(6) + zz(3)*Zln(3,k) + zz(9)*Zln(9,k) + zz(15)&
            *Zln(15,k) + zz(21)*Zln(21,k) + zz(27)*Zln(27,k)&
            + zz(33)*Zln(33,k)
       t(7) = t(7) + zz(1)*Zln(4,k) + zz(7)*Zln(10,k) + zz(13)&
            *Zln(16,k) + zz(19)*Zln(22,k) + zz(25)*Zln(28,k)&
            + zz(31)*Zln(34,k)
       t(8) = t(8) + zz(2)*Zln(4,k) + zz(8)*Zln(10,k) + zz(14)&
            *Zln(16,k) + zz(20)*Zln(22,k) + zz(26)*Zln(28,k)&
            + zz(32)*Zln(34,k)
       t(9) = t(9) + zz(3)*Zln(4,k) + zz(9)*Zln(10,k) + zz(15)&
            *Zln(16,k) + zz(21)*Zln(22,k) + zz(27)*Zln(28,k)&
            + zz(33)*Zln(34,k)
       t(10) = t(10) + zz(4)*Zln(4,k) + zz(10)*Zln(10,k) + zz(16)&
            *Zln(16,k) + zz(22)*Zln(22,k) + zz(28)*Zln(28,k)&
            + zz(34)*Zln(34,k)
       t(11) = t(11) + zz(1)*Zln(5,k) + zz(7)*Zln(11,k) + zz(13)&
            *Zln(17,k) + zz(19)*Zln(23,k) + zz(25)*Zln(29,k)&
            + zz(31)*Zln(35,k)
       t(12) = t(12) + zz(2)*Zln(5,k) + zz(8)*Zln(11,k) + zz(14)&
            *Zln(17,k) + zz(20)*Zln(23,k) + zz(26)*Zln(29,k)&
            + zz(32)*Zln(35,k)
       t(13) = t(13) + zz(3)*Zln(5,k) + zz(9)*Zln(11,k) + zz(15)&
            *Zln(17,k) + zz(21)*Zln(23,k) + zz(27)*Zln(29,k)&
            + zz(33)*Zln(35,k)
       t(14) = t(14) + zz(4)*Zln(5,k) + zz(10)*Zln(11,k) + zz(16)&
            *Zln(17,k) + zz(22)*Zln(23,k) + zz(28)*Zln(29,k)&
            + zz(34)*Zln(35,k)
       t(15) = t(15) + zz(5)*Zln(5,k) + zz(11)*Zln(11,k) + zz(17)&
            *Zln(17,k) + zz(23)*Zln(23,k) + zz(29)*Zln(29,k)&
            + zz(35)*Zln(35,k)
       t(16) = t(16) + zz(1)*Zln(6,k) + zz(7)*Zln(12,k) + zz(13)&
            *Zln(18,k) + zz(19)*Zln(24,k) + zz(25)*Zln(30,k)&
            + zz(31)*Zln(36,k)
       t(17) = t(17) + zz(2)*Zln(6,k) + zz(8)*Zln(12,k) + zz(14)&
            *Zln(18,k) + zz(20)*Zln(24,k) + zz(26)*Zln(30,k)&
            + zz(32)*Zln(36,k)
       t(18) = t(18) + zz(3)*Zln(6,k) + zz(9)*Zln(12,k) + zz(15)&
            *Zln(18,k) + zz(21)*Zln(24,k) + zz(27)*Zln(30,k)&
            + zz(33)*Zln(36,k)
       t(19) = t(19) + zz(4)*Zln(6,k) + zz(10)*Zln(12,k) + zz(16)&
            *Zln(18,k) + zz(22)*Zln(24,k) + zz(28)*Zln(30,k)&
            + zz(34)*Zln(36,k)
       t(20) = t(20) + zz(5)*Zln(6,k) + zz(11)*Zln(12,k) + zz(17)&
            *Zln(18,k) + zz(23)*Zln(24,k) + zz(29)*Zln(30,k)&
            + zz(35)*Zln(36,k)
       t(21) = t(21) + zz(6)*Zln(6,k) + zz(12)*Zln(12,k) + zz(18)&
            *Zln(18,k) + zz(24)*Zln(24,k) + zz(30)*Zln(30,k)&
            + zz(36)*Zln(36,k)
    ENDDO
    !$dir max_trips(6)
    DO l = 1, 21
       Diag(l,Ic) = Diag(l,Ic) - t(l)
    ENDDO
    CALL INV6(Diag(1,Ic),ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  END SUBROUTINE S6UM

  !======================================================================!
  !> @brief S6UM1
  !======================================================================!
  SUBROUTINE S6UM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(9,*)

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: s(9)
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    !$dir max_trip(9)
    DO l = 1, 9
       s(l) = 0.0D0
    ENDDO
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             s(1) = s(1) + Temp(1,j)*Zln(1,jj) + Temp(4,j)*Zln(4,jj) + Temp(7,j)*Zln(7,jj)
             s(2) = s(2) + Temp(2,j)*Zln(1,jj) + Temp(5,j)*Zln(4,jj) + Temp(8,j)*Zln(7,jj)
             s(3) = s(3) + Temp(3,j)*Zln(1,jj) + Temp(6,j)*Zln(4,jj) + Temp(9,j)*Zln(7,jj)
             s(4) = s(4) + Temp(1,j)*Zln(2,jj) + Temp(4,j)*Zln(5,jj) + Temp(7,j)*Zln(8,jj)
             s(5) = s(5) + Temp(2,j)*Zln(2,jj) + Temp(5,j)*Zln(5,jj) + Temp(8,j)*Zln(8,jj)
             s(6) = s(6) + Temp(3,j)*Zln(2,jj) + Temp(6,j)*Zln(5,jj) + Temp(9,j)*Zln(8,jj)
             s(7) = s(7) + Temp(1,j)*Zln(3,jj) + Temp(4,j)*Zln(6,jj) + Temp(7,j)*Zln(9,jj)
             s(8) = s(8) + Temp(2,j)*Zln(3,jj) + Temp(5,j)*Zln(6,jj) + Temp(8,j)*Zln(9,jj)
             s(9) = s(9) + Temp(3,j)*Zln(3,jj) + Temp(6,j)*Zln(6,jj) + Temp(9,j)*Zln(9,jj)
          ENDIF
       ENDDO
       !$dir max_trip(9)
       DO l = 1, 9
          Temp(l,jc) = Zln(l,k) - s(l)
          Zln(l,k) = Temp(l,jc)
          s(l) = 0.0D0
       ENDDO
    ENDDO
  END SUBROUTINE S6UM1

  !======================================================================!
  !> @brief S6UM2
  !======================================================================!
  SUBROUTINE S6UM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(21,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(36,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(36,Neqns)
    DOUBLE PRECISION, INTENT(OUT):: Zln(36,*)

    INTEGER:: ic
    INTEGER:: j
    INTEGER:: j1
    INTEGER:: j2
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: m

    joc = 0
    DO ic = Nstop, Neqns
       DO m = 1, 36
          DO jj = 1, Nstop
             Temp(jj,m) = 0.0D0
          ENDDO
       ENDDO
       ks = Xlnzr(ic)
       ke = Xlnzr(ic+1) - 1
       DO k = ks, ke
          jj = Colno(k)
          Temp(:,jj) = Zln(:,k)
          Indx(jj) = ic
       ENDDO
       DO k = ks, ke
          jj = Colno(k)
          CALL INV66(Zln(1,k),Temp,Diag(1,jj))
       ENDDO

       DO k = ks, ke
          jj = Colno(k)
          Diag(1,ic) = Diag(1,ic) - Temp(jj,1)*Zln(1,k) - Temp(jj,4)*Zln(4,k) - Temp(jj,7)*Zln(7,k)
          Diag(2,ic) = Diag(2,ic) - Temp(jj,1)*Zln(2,k) - Temp(jj,4)*Zln(5,k) - Temp(jj,7)*Zln(8,k)
          Diag(3,ic) = Diag(3,ic) - Temp(jj,2)*Zln(2,k) - Temp(jj,5)*Zln(5,k) - Temp(jj,8)*Zln(8,k)
          Diag(4,ic) = Diag(4,ic) - Temp(jj,1)*Zln(3,k) - Temp(jj,4)*Zln(6,k) - Temp(jj,7)*Zln(9,k)
          Diag(5,ic) = Diag(5,ic) - Temp(jj,2)*Zln(3,k) - Temp(jj,5)*Zln(6,k) - Temp(jj,8)*Zln(9,k)
          Diag(6,ic) = Diag(6,ic) - Temp(jj,3)*Zln(3,k) - Temp(jj,6)*Zln(6,k) - Temp(jj,9)*Zln(9,k)
       ENDDO
       DO jc = Nstop, ic - 1
          joc = joc + 1
          j1 = Xlnzr(jc)
          j2 = Xlnzr(jc+1)
          DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
             j = Colno(jj)
             IF ( Indx(j)==ic ) THEN
                Dsln(1,joc) = Dsln(1,joc) - Temp(j,1)*Zln(1,jj)&
                     - Temp(j,4)*Zln(4,jj) - Temp(j,7)&
                     *Zln(7,jj)
                Dsln(2,joc) = Dsln(2,joc) - Temp(j,2)*Zln(1,jj)&
                     - Temp(j,5)*Zln(4,jj) - Temp(j,8)&
                     *Zln(7,jj)
                Dsln(3,joc) = Dsln(3,joc) - Temp(j,3)*Zln(1,jj)&
                     - Temp(j,6)*Zln(4,jj) - Temp(j,9)&
                     *Zln(7,jj)
                Dsln(4,joc) = Dsln(4,joc) - Temp(j,1)*Zln(2,jj)&
                     - Temp(j,4)*Zln(5,jj) - Temp(j,7)&
                     *Zln(8,jj)
                Dsln(5,joc) = Dsln(5,joc) - Temp(j,2)*Zln(2,jj)&
                     - Temp(j,5)*Zln(5,jj) - Temp(j,8)&
                     *Zln(8,jj)
                Dsln(6,joc) = Dsln(6,joc) - Temp(j,3)*Zln(2,jj)&
                     - Temp(j,6)*Zln(5,jj) - Temp(j,9)&
                     *Zln(8,jj)
                Dsln(7,joc) = Dsln(7,joc) - Temp(j,1)*Zln(3,jj)&
                     - Temp(j,4)*Zln(6,jj) - Temp(j,7)&
                     *Zln(9,jj)
                Dsln(8,joc) = Dsln(8,joc) - Temp(j,2)*Zln(3,jj)&
                     - Temp(j,5)*Zln(6,jj) - Temp(j,8)&
                     *Zln(9,jj)
                Dsln(9,joc) = Dsln(9,joc) - Temp(j,3)*Zln(3,jj)&
                     - Temp(j,6)*Zln(6,jj) - Temp(j,9)&
                     *Zln(9,jj)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE S6UM2

  !======================================================================!
  !> @brief S6UM3
  !======================================================================!
  SUBROUTINE S6UM3(N,Dsln,Diag,Indx,Temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(6,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(9,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(9,*)

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: l
    DOUBLE PRECISION:: t(9)

    IF ( N>0 ) THEN
       Indx(1) = 0
       joc = 1
       CALL INV3(Diag(1,1),ir)
       DO i = 2, N
          Indx(i) = joc
          DO j = 1, i - 1
             CALL D3DOT(t,Dsln(1,Indx(i)),Dsln(1,Indx(j)),j-1)
             !$dir max_trips(9)
             DO l = 1, 9
                Dsln(l,joc) = Dsln(l,joc) - t(l)
             ENDDO
             joc = joc + 1
          ENDDO
          CALL V3PROD(Dsln(1,Indx(i)),Diag,Temp,i-1)
          CALL D3DOTL(t,Temp,Dsln(1,Indx(i)),i-1)
          !$dir max_trips(6)
          DO l = 1, 6
             Diag(l,i) = Diag(l,i) - t(l)
          ENDDO
          CALL VCOPY(Temp,Dsln(1,Indx(i)),9*(i-1))
          CALL INV3(Diag(1,i),ir)
       ENDDO
    ENDIF
  END SUBROUTINE S6UM3

  !======================================================================!
  !> @brief STAIJ1 routine sets an non-zero entry  of the matrix.
  !      (symmetric version)
  !      (i)
  !          isw      =0    set the value
  !                   =1    add the value
  !          i        row entry
  !          j        column entry
  !          aij      value
  !      (o)
  !          iv       communication array
  !        #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE STAIJ1(Isw,I,J,Aij,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: I
    INTEGER, INTENT(IN):: Isw
    INTEGER, INTENT(IN):: J
    DOUBLE PRECISION, INTENT(IN):: Aij(NDEg*NDEg)
    INTEGER, INTENT(OUT):: Ir

    INTEGER:: IDBg
    INTEGER:: LRAtio
    INTEGER:: ndeg2
    INTEGER:: ndeg2l
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    COMMON /DEBUG / IDBg
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    Ir = 0
    ndeg2 = NDEg*NDEg
    ndeg2l = NDEg*(NDEg+1)/2
    IF ( STAge==30 ) WRITE (6,*) 'warning a matrix was build up '//'but never solved.'
    IF ( STAge==10 ) THEN
       ALLOCATE (DIAg(NEQns*ndeg2l),STAT=IERror)
       RALoc = RALoc + NEQns*ndeg2l
       ALLOCATE (ZLN(LEN_colno*ndeg2),STAT=IERror)

       RALoc = RALoc + LEN_colno*ndeg2
       ALLOCATE (DSLn(LEN_dsln*ndeg2),STAT=IERror)

       IF ( IERror/=0 ) STOP "Allocation error dsln"

       RALoc = RALoc + LEN_dsln*ndeg2
    ENDIF
    IF ( STAge/=20 ) THEN
       !
       ! for diagonal
       !
       DIAg = 0.
       !
       ! for lower triangle
       !
       ZLN = 0.
       !
       ! for dense window
       !
       DSLn = 0.

       STAge = 20
    ENDIF
    !         Print *,'********Set Stage 20 *********'
    !
    IF ( NDEg<=2 ) THEN
       CALL ADDR0(Isw,I,J,Aij,INVp,XLNzr,COLno,DIAg,ZLN,DSLn,NSTop,ndeg2,ndeg2l,Ir)
    ELSEIF ( NDEg==3 ) THEN
       CALL ADDR3(I,J,Aij,INVp,XLNzr,COLno,DIAg,ZLN,DSLn,NSTop,Ir)
    ELSE
       CALL ADDRX(I,J,Aij,INVp,XLNzr,COLno,DIAg,ZLN,DSLn,NSTop,NDEg,ndeg2l,Ir)
    ENDIF
  END SUBROUTINE STAIJ1

  !======================================================================!
  !> @brief STIAJA routine sets an non-zero entry  of the matrix.
  !      (asymmetric version)
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE STIAJA(Neqns,Ia,Ja,Jcpt,Jcolno)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Jcpt(*)
    INTEGER, INTENT(IN):: Jcolno(*)
    INTEGER, INTENT(OUT):: Ia(*)
    INTEGER, INTENT(OUT):: Ja(*)

    INTEGER:: i
    INTEGER:: IDBg
    INTEGER:: ii
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: l
    COMMON /DEBUG / IDBg
    !
    Ia(1) = 1
    l = 0
    DO k = 1, Neqns
       joc = Jcpt(k)
       DO WHILE ( joc/=0 )
          ii = Jcolno(joc)
          IF ( ii/=k ) THEN
             l = l + 1
             Ja(l) = ii
          ENDIF
          joc = Jcpt(joc)
       ENDDO
       Ia(k+1) = l + 1
    ENDDO
    IF ( IDBg/=0 ) THEN
       WRITE (6,*) 'ia '
       WRITE (6,"(10I7)") (Ia(i),i=1,Neqns)
       WRITE (6,*) 'ja '
       WRITE (6,"(10I7)") (Ja(i),i=1,Ia(Neqns+1))
    ENDIF
  END SUBROUTINE STIAJA

  !======================================================================!
  !> @brief STSMAT
  !     coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE STSMAT(Neqns,Nttbr,Irow,Jcol,Jcpt,Jcolno)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nttbr
    INTEGER, INTENT(IN):: Irow(*)
    INTEGER, INTENT(IN):: Jcol(*)
    INTEGER, INTENT(OUT):: Jcpt(*)
    INTEGER, INTENT(OUT):: Jcolno(*)

    INTEGER:: i
    INTEGER:: IDBg
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: l
    INTEGER:: locr
    COMMON /DEBUG / IDBg

    DO i = 1, 2*Nttbr
       Jcpt(i) = 0
       Jcolno(i) = 0
    ENDDO
    DO i = 1, Neqns
       Jcpt(i) = i + Neqns
       Jcolno(i+Neqns) = i
    ENDDO

    k = 2*Neqns

    DO l = 1, Nttbr
       i = Irow(l)
       j = Jcol(l)
       IF ( i/=j ) THEN
          joc = Jcpt(i)
          locr = i
          DO WHILE ( joc/=0 )
             IF ( Jcolno(joc)==j ) GOTO 100
             IF ( Jcolno(joc)>j ) THEN
                k = k + 1
                Jcpt(locr) = k
                Jcpt(k) = joc
                Jcolno(k) = j
                GOTO 20
             ENDIF
             locr = joc
             joc = Jcpt(joc)
          ENDDO
          k = k + 1
          Jcpt(locr) = k
          Jcolno(k) = j
20        joc = Jcpt(j)
          locr = j
          DO WHILE ( joc/=0 )
             IF ( Jcolno(joc)==i ) GOTO 100
             IF ( Jcolno(joc)>i ) THEN
                k = k + 1
                Jcpt(locr) = k
                Jcpt(k) = joc
                Jcolno(k) = i
                GOTO 100
             ENDIF
             locr = joc
             joc = Jcpt(joc)
          ENDDO
          k = k + 1
          Jcpt(locr) = k
          Jcolno(k) = i
       ENDIF
100 ENDDO
    IF ( IDBg/=0 ) THEN
       WRITE (6,*) 'jcolno'
       WRITE (6,"(10I7)") (Jcolno(i),i=1,k)
       WRITE (6,*) 'jcpt'
       WRITE (6,"(10I7)") (Jcpt(i),i=1,k)
    ENDIF
  END SUBROUTINE STSMAT

  !======================================================================!
  !> @brief SUM
  !======================================================================!
  SUBROUTINE SUM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Par(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(*)

    INTEGER:: ISEm
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: piv
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: s
    DOUBLE PRECISION:: t
    DOUBLE PRECISION:: zz
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio
    COMMON ISEm

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t = 0.0D0

    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       s = 0.0D0
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) s = s + Temp(j)*Zln(jj)
       ENDDO

       zz = Zln(k) - s
       Zln(k) = zz*Diag(jc)
       Temp(jc) = zz
       t = t + zz*Zln(k)
    ENDDO
    piv = Diag(Ic) - t
    IF ( DABS(piv)>RMIn ) Diag(Ic) = 1.0D0/piv
    DO WHILE ( ISEm/=1 )
    ENDDO
    ISEm = 0
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
    ISEm = 1
  END SUBROUTINE SUM

  !======================================================================!
  !> @brief SUM1
  !======================================================================!
  SUBROUTINE SUM1(Ic,Xlnzr,Colno,Zln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(*)

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: LRAtio
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    DOUBLE PRECISION:: s
    DOUBLE PRECISION:: t
    DOUBLE PRECISION:: zz
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    t = 0.0D0

    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       s = 0.0D0
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) s = s + Temp(j)*Zln(jj)
       ENDDO
       zz = Zln(k) - s

       Zln(k) = zz
       Temp(jc) = zz
    ENDDO
  END SUBROUTINE SUM1

  !======================================================================!
  !> @brief SUM2
  !======================================================================!
  SUBROUTINE SUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(*)

    INTEGER:: i
    INTEGER:: ic
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    DOUBLE PRECISION:: s

    joc = 0
    DO ic = Nstop, Neqns
       DO i = 1, Nstop
          Temp(i) = 0.0D0
       ENDDO
       ks = Xlnzr(ic)
       ke = Xlnzr(ic+1) - 1
       DO k = ks, ke
          jj = Colno(k)
          Temp(jj) = Zln(k)
          Zln(k) = Temp(jj)*Diag(jj)
          Indx(jj) = ic
          Diag(ic) = Diag(ic) - Temp(jj)*Zln(k)
       ENDDO
       DO jc = Nstop, ic - 1
          s = 0.0D0
          joc = joc + 1
          DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
             j = Colno(jj)
             IF ( Indx(j)==ic ) s = s + Temp(j)*Zln(jj)
          ENDDO
          IF ( s==0.0D0 ) WRITE (16,*) ic, jc
          Dsln(joc) = Dsln(joc) - s
       ENDDO
    ENDDO
  END SUBROUTINE SUM2

  !======================================================================!
  !> @brief SUM3
  !======================================================================!
  SUBROUTINE SUM3(N,Dsln,Diag,Indx,Temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(*)

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc

    IF ( N>0 ) THEN
       Indx(1) = 0
       joc = 1
       Diag(1) = 1.0D0/Diag(1)
       DO i = 2, N
          Indx(i) = joc
          DO j = 1, i - 1
             Dsln(joc) = Dsln(joc) - DDOT(Dsln(Indx(i)),Dsln(Indx(j)),j-1)
             joc = joc + 1
          ENDDO
          CALL VPROD(Dsln(Indx(i)),Diag,Temp,i-1)
          Diag(i) = Diag(i) - DDOT(Temp,Dsln(Indx(i)),i-1)
          CALL VCOPY(Temp,Dsln(Indx(i)),i-1)
          Diag(i) = 1.0D0/Diag(i)
       ENDDO
    ENDIF
  END SUBROUTINE SUM3

  !======================================================================!
  !> @brief SXUM
  !======================================================================!
  SUBROUTINE SXUM(Ic,Xlnzr,Colno,Zln,Diag,Nch,Par,Temp,Indx,Ndeg,Ndegl,Zz,T)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ic
    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndegl
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(IN):: Par(*)
    INTEGER, INTENT(OUT):: Indx(*)
    INTEGER, INTENT(OUT):: Nch(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(Ndegl,*)
    DOUBLE PRECISION, INTENT(OUT):: T(Ndegl)
    DOUBLE PRECISION, INTENT(OUT):: Temp(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Zz(Ndeg,Ndeg)

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: LRAtio
    INTEGER:: m
    INTEGER:: n
    INTEGER:: ndeg22
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ndeg22 = Ndeg*Ndeg
    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    T = 0.0
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       Zz = Zln(:,:,k)
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             DO m = 1, Ndeg, 2
                DO n = 1, Ndeg, 2
                   DO kk = 1, Ndeg, 2
                      Zz(n,m) = Zz(n,m) - Temp(n,kk,j)*Zln(m,kk,jj) - Temp(n,kk+1,j)*Zln(m,kk+1,jj)
                      Zz(n,m+1) = Zz(n,m+1) - Temp(n,kk,j)*Zln(m+1,kk,jj) - Temp(n,kk+1,j)*Zln(m+1,kk+1,jj)
                      Zz(n+1,m) = Zz(n+1,m) - Temp(n+1,kk,j)*Zln(m,kk,jj) - Temp(n+1,kk+1,j)*Zln(m,kk+1,jj)
                      Zz(n+1,m+1) = Zz(n+1,m+1) - Temp(n+1,kk,j)*Zln(m+1,kk,jj) - Temp(n+1,kk+1,j)*Zln(m+1,kk+1,jj)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       CALL INVXX(Zln(1,1,k),Zz,Diag(1,jc),Ndeg)

       Temp(:,:,jc) = Zz
       joc = 0
       DO n = 1, Ndeg
          DO m = 1, n
             joc = joc + 1
             DO kk = 1, Ndeg, 2
                T(joc) = T(joc) + Zz(n,kk)*Zln(m,kk,k) + Zz(n,kk+1)*Zln(m,kk+1,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    Diag(:,Ic) = Diag(:,Ic) - T
    CALL INVX(Diag(1,Ic),Ndeg,ir)
    Nch(Ic) = -1
    kk = Par(Ic)
    Nch(kk) = Nch(kk) - 1
  END SUBROUTINE SXUM

  !======================================================================!
  !> @brief SXUM1
  !======================================================================!
  SUBROUTINE SXUM1(Ic,Xlnzr,Colno,Zln,Temp,Indx,Ndeg,S)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: S(Ndeg,Ndeg)
    DOUBLE PRECISION, INTENT(OUT):: Temp(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg,Ndeg,*)

    INTEGER:: Ic
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: LRAtio
    INTEGER:: m
    INTEGER:: n
    DOUBLE PRECISION:: EPSm
    DOUBLE PRECISION:: RMAx
    DOUBLE PRECISION:: RMIn
    COMMON /MCHDPN/ RMAx, RMIn, EPSm, LRAtio

    ks = Xlnzr(Ic)
    ke = Xlnzr(Ic+1)
    DO m = 1, Ndeg
       DO n = 1, Ndeg
          S(n,m) = 0.0D0
       ENDDO
    ENDDO
    DO k = ks, ke - 1
       jc = Colno(k)
       Indx(jc) = Ic
       DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
          j = Colno(jj)
          IF ( Indx(j)==Ic ) THEN
             DO m = 1, Ndeg
                DO n = 1, Ndeg
                   DO kk = 1, Ndeg
                      S(n,m) = S(n,m) + Temp(n,kk,j)*Zln(m,kk,jj)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       DO m = 1, Ndeg
          DO n = 1, Ndeg
             Temp(n,m,jc) = Zln(n,m,k) - S(n,m)
             Zln(n,m,k) = Temp(n,m,jc)
             S(n,m) = 0.0D0
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE SXUM1

  !======================================================================!
  !> @brief SXUM2
  !======================================================================!
  SUBROUTINE SXUM2(Neqns,Nstop,Xlnzr,Colno,Zln,Diag,Dsln,Temp,Indx,Ndeg,Ndegl)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndegl
    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nstop
    INTEGER, INTENT(IN):: Xlnzr(*)
    INTEGER, INTENT(IN):: Colno(*)
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(Ndegl,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Temp(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg,Ndeg,*)

    INTEGER:: ic
    INTEGER:: j
    INTEGER:: j1
    INTEGER:: j2
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: locd
    INTEGER:: m
    INTEGER:: n

    joc = 0
    DO ic = Nstop, Neqns
       ks = Xlnzr(ic)
       ke = Xlnzr(ic+1) - 1
       DO k = ks, ke
          jj = Colno(k)
          DO m = 1, Ndeg
             DO n = 1, Ndeg
                Temp(n,m,jj) = Zln(n,m,k)
                Indx(jj) = ic
             ENDDO
          ENDDO
       ENDDO
       DO k = ks, ke
          jj = Colno(k)
          CALL INVXX(Zln(1,1,k),Temp(1,1,jj),Diag(1,jj),Ndeg)
       ENDDO

       locd = 0
       DO n = 1, Ndeg
          DO m = 1, n
             locd = locd + 1
             DO k = ks, ke
                jj = Colno(k)
                DO kk = 1, Ndeg
                   Diag(locd,ic) = Diag(locd,ic) - Temp(n,kk,jj)*Zln(m,kk,k)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO jc = Nstop, ic - 1
          joc = joc + 1
          j1 = Xlnzr(jc)
          j2 = Xlnzr(jc+1)
          DO jj = Xlnzr(jc), Xlnzr(jc+1) - 1
             j = Colno(jj)
             IF ( Indx(j)==ic ) THEN
                DO m = 1, Ndeg
                   DO n = 1, Ndeg
                      DO k = 1, Ndeg
                         Dsln(n,m,joc) = Dsln(n,m,joc) - Temp(n,k,j)*Zln(m,k,jj)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE SXUM2

  !======================================================================!
  !> @brief SXUM3
  !======================================================================!
  SUBROUTINE SXUM3(Nn,Dsln,Diag,Indx,Temp,Ndeg,Ndegl,T)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndegl
    INTEGER, INTENT(IN):: Nn
    INTEGER, INTENT(OUT):: Indx(*)
    DOUBLE PRECISION, INTENT(OUT):: Diag(Ndegl,*)
    DOUBLE PRECISION, INTENT(OUT):: Dsln(Ndeg,Ndeg,*)
    DOUBLE PRECISION, INTENT(OUT):: T(Ndeg,Ndeg)
    DOUBLE PRECISION, INTENT(OUT):: Temp(Ndeg,Ndeg,*)

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: locd
    INTEGER:: m
    INTEGER:: n

    IF ( Nn>0 ) THEN
       Indx(1) = 0
       joc = 1
       CALL INVX(Diag(1,1),Ndeg,ir)
       DO i = 2, Nn
          Indx(i) = joc
          DO j = 1, i - 1
             CALL DXDOT(Ndeg,T,Dsln(1,1,Indx(i)),Dsln(1,1,Indx(j)),j-1)
             DO m = 1, Ndeg
                DO n = 1, Ndeg
                   Dsln(n,m,joc) = Dsln(n,m,joc) - T(n,m)
                ENDDO
             ENDDO
             joc = joc + 1
          ENDDO
          CALL VXPROD(Ndeg,Ndegl,Dsln(1,1,Indx(i)),Diag,Temp,i-1)
          CALL DXDOTL(Ndeg,T,Temp,Dsln(1,1,Indx(i)),i-1)
          locd = 0
          DO n = 1, Ndeg
             DO m = 1, n
                locd = locd + 1
                Diag(locd,i) = Diag(locd,i) - T(n,m)
             ENDDO
          ENDDO
          CALL VCOPY(Temp,Dsln(1,1,Indx(i)),Ndeg*Ndeg*(i-1))
          CALL INVX(Diag(1,i),Ndeg,ir)
       ENDDO
    ENDIF
  END SUBROUTINE SXUM3

  !======================================================================!
  !> @brief V2PROD
  !======================================================================!
  SUBROUTINE V2PROD(A,B,C,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(4,N)
    DOUBLE PRECISION, INTENT(IN):: B(3,N)
    DOUBLE PRECISION, INTENT(OUT):: C(4,N)

    INTEGER:: i

    DO i = 1, N
       C(3,i) = A(3,i) - A(1,i)*B(2,i)
       C(1,i) = A(1,i)*B(1,i)
       C(3,i) = C(3,i)*B(3,i)
       C(1,i) = C(1,i) - C(3,i)*B(2,i)

       C(4,i) = A(4,i) - A(2,i)*B(2,i)
       C(2,i) = A(2,i)*B(1,i)
       C(4,i) = C(4,i)*B(3,i)
       C(2,i) = C(2,i) - C(4,i)*B(2,i)
    ENDDO
  END SUBROUTINE V2PROD

  !======================================================================!
  !> @brief V3PROD
  !======================================================================!
  SUBROUTINE V3PROD(Zln,Diag,Zz,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: Diag(6,N)
    DOUBLE PRECISION, INTENT(IN):: Zln(9,N)
    DOUBLE PRECISION, INTENT(OUT):: Zz(9,N)

    INTEGER:: i

    DO i = 1, N
       Zz(4,i) = Zln(4,i) - Zln(1,i)*Diag(2,i)
       Zz(7,i) = Zln(7,i) - Zln(1,i)*Diag(4,i) - Zz(4,i)*Diag(5,i)
       Zz(1,i) = Zln(1,i)*Diag(1,i)
       Zz(4,i) = Zz(4,i)*Diag(3,i)
       Zz(7,i) = Zz(7,i)*Diag(6,i)
       Zz(4,i) = Zz(4,i) - Zz(7,i)*Diag(5,i)
       Zz(1,i) = Zz(1,i) - Zz(4,i)*Diag(2,i) - Zz(7,i)*Diag(4,i)

       Zz(5,i) = Zln(5,i) - Zln(2,i)*Diag(2,i)
       Zz(8,i) = Zln(8,i) - Zln(2,i)*Diag(4,i) - Zz(5,i)*Diag(5,i)
       Zz(2,i) = Zln(2,i)*Diag(1,i)
       Zz(5,i) = Zz(5,i)*Diag(3,i)
       Zz(8,i) = Zz(8,i)*Diag(6,i)
       Zz(5,i) = Zz(5,i) - Zz(8,i)*Diag(5,i)
       Zz(2,i) = Zz(2,i) - Zz(5,i)*Diag(2,i) - Zz(8,i)*Diag(4,i)

       Zz(6,i) = Zln(6,i) - Zln(3,i)*Diag(2,i)
       Zz(9,i) = Zln(9,i) - Zln(3,i)*Diag(4,i) - Zz(6,i)*Diag(5,i)
       Zz(3,i) = Zln(3,i)*Diag(1,i)
       Zz(6,i) = Zz(6,i)*Diag(3,i)
       Zz(9,i) = Zz(9,i)*Diag(6,i)
       Zz(6,i) = Zz(6,i) - Zz(9,i)*Diag(5,i)
       Zz(3,i) = Zz(3,i) - Zz(6,i)*Diag(2,i) - Zz(9,i)*Diag(4,i)
    ENDDO
  END SUBROUTINE V3PROD

  !======================================================================!
  !> @brief VCOPY
  !======================================================================!
  SUBROUTINE VCOPY(A,C,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(N)
    DOUBLE PRECISION, INTENT(OUT):: C(N)

    C = A
  END SUBROUTINE VCOPY

  !======================================================================!
  !> @brief VPROD
  !======================================================================!
  SUBROUTINE VPROD(A,B,C,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: N
    DOUBLE PRECISION, INTENT(IN):: A(N)
    DOUBLE PRECISION, INTENT(IN):: B(N)
    DOUBLE PRECISION, INTENT(OUT):: C(N)

    INTEGER:: i

    DO i = 1, N
       C(i) = A(i)*B(i)
    ENDDO
  END SUBROUTINE VPROD

  !======================================================================!
  !> @brief VXPROD
  !======================================================================!
  SUBROUTINE VXPROD(Ndeg,Ndegl,Zln,Diag,Zz,N)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Ndeg
    INTEGER, INTENT(IN):: Ndegl
    DOUBLE PRECISION, INTENT(IN):: Diag(Ndegl,N)
    DOUBLE PRECISION, INTENT(OUT):: Zln(Ndeg*Ndeg,N)
    DOUBLE PRECISION, INTENT(OUT):: Zz(Ndeg*Ndeg,N)

    INTEGER:: i
    INTEGER:: N

    DO i = 1, N
       CALL INVXX(Zz(1,i),Zln(1,i),Diag(1,i),Ndeg)
    ENDDO
  END SUBROUTINE VXPROD

  !======================================================================!
  !> @brief ZPIVOT
  !======================================================================!
  SUBROUTINE ZPIVOT(Neqns,Neqnsz,Nttbr,Jcol,Irow,Zpiv,Ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Neqns
    INTEGER, INTENT(IN):: Nttbr
    INTEGER, INTENT(IN):: Jcol(*)
    INTEGER, INTENT(IN):: Irow(*)
    INTEGER, INTENT(OUT):: Ir
    INTEGER, INTENT(OUT):: Neqnsz
    INTEGER, INTENT(OUT):: Zpiv(*)

    INTEGER:: i
    INTEGER:: IDBg
    INTEGER:: j
    INTEGER:: l
    COMMON /DEBUG / IDBg
    !
    Ir = 0
    DO l = 1, Neqns
       Zpiv(l) = 1
    ENDDO

    DO l = 1, Nttbr
       i = Irow(l)
       j = Jcol(l)
       IF ( i<=0 .OR. j<=0 ) THEN
          Ir = -1
          GOTO 100
       ELSEIF ( i>Neqns .OR. j>Neqns ) THEN
          Ir = 1
          GOTO 100
       ENDIF
       IF ( i==j ) Zpiv(i) = 0
    ENDDO
    DO i = Neqns, 1, -1
       IF ( Zpiv(i)==0 ) THEN
          Neqnsz = i
          EXIT
       ENDIF
    ENDDO
100 IF ( IDBg/=0 ) WRITE (6,"(20I3)") (Zpiv(i),i=1,Neqns)
  END SUBROUTINE ZPIVOT
END MODULE HECMW_SOLVER_DIRECT
