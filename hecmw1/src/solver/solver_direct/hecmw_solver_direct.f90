!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

MODULE hecmw_solver_direct
  IMPLICIT NONE
  INTEGER (KIND=4),PRIVATE :: len_colno
  INTEGER (KIND=4),PRIVATE :: nstop
  INTEGER (KIND=4),PRIVATE :: stage
  INTEGER (KIND=4),PRIVATE :: neqns
  INTEGER (KIND=4),PRIVATE :: nttbr
  INTEGER (KIND=4),PRIVATE :: isym
  INTEGER (KIND=4),PRIVATE :: ndeg
  INTEGER (KIND=4),PRIVATE :: irr
  INTEGER (KIND=4),PRIVATE :: len_dsln
  INTEGER (KIND=4),PRIVATE :: len_iv
  INTEGER (KIND=4),PRIVATE :: total

  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: jcol
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: irow
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: zpiv
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: iperm
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: invp
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: parent
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: nch
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: xlnzr
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: colno
  !*Work arrays
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: jcpt
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: jcolno
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: ia
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: ja
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: deg
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: marker
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: rchset
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: nbrhd
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: qsize
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: qlink
  !WINDEBUG
  !        INTEGER (KIND=4),private, POINTER  :: nofsub(:)
  INTEGER (KIND=4),PRIVATE :: nofsub
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: adjncy
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: btree
  !DEBUGWIN
  !INTEGER (KIND=4),private :: btree(100000)
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: pordr
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: adjncp
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: xleaf
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: leaf
  INTEGER (KIND=4),DIMENSION(:), PRIVATE, POINTER  :: indx

  REAL(KIND=8),DIMENSION(:), PRIVATE, POINTER  :: val
  REAL(KIND=8),DIMENSION(:), PRIVATE, POINTER  :: temp
  REAL(KIND=8),DIMENSION(:), PRIVATE, POINTER  :: diag
  REAL(KIND=8),DIMENSION(:), PRIVATE, POINTER  :: zln
  REAL(KIND=8),DIMENSION(:), PRIVATE, POINTER  :: dsln
  !*Timing
  REAL(KIND=8),DIMENSION(10), PRIVATE  :: tom
  !*Allocation variables
  INTEGER(KIND=4),PRIVATE :: ialoc
  INTEGER(KIND=4),PRIVATE :: raloc
  INTEGER(KIND=4),PRIVATE :: ierror

CONTAINS
  !======================================================================!
  !> hecmw_solve_direct
  !     this is a sample program for the matrix solver
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
  !======================================================================!
  SUBROUTINE hecmw_solve_direct(hecMESH,hecMAT,ifmsg)
    USE hecmw_util
    USE hecmw_matrix_ass
    USE hecmw_matrix_dump
    IMPLICIT NONE

    TYPE (hecmwST_local_mesh), INTENT(IN):: hecMESH
    INTEGER, INTENT(IN):: ifmsg

    TYPE (hecmwST_matrix    ), INTENT(OUT):: hecMAT

    INTEGER:: iseed
    INTEGER:: ixxx
    INTEGER:: lratio
    INTEGER:: i98
    INTEGER:: i97
    INTEGER:: ir
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    COMMON /mchdpn/ rmax,rmin,epsm,lratio
    COMMON /qaz/ iseed,ixxx

    rmax=8.988d+307
    rmin=4.941d-300
    epsm=2.220d-16
    lratio=2
    iseed=1
    ir=0

    CALL hecmw_mat_ass_equation( hecMESH, hecMAT )
    CALL hecmw_mat_dump(hecMAT, hecMESH)

    !WRITE(IFMSG,*) "Interface to ADS from HECMW..."
    CALL ptime(t1)

    !*EHM HECMW June 7 2004

    i98 = hecMAT%Iarray(98)
    IF(hecMAT%Iarray(98)==1) THEN
      !* Interface to symbolic factorization
      CALL setij(hecMESH,hecMAT)

      !* Symbolic factorization
      CALL matini(ir)
      hecMAT%Iarray(98) = 0
      WRITE(6,*)"symbolic fct done"
    ENDIF
    CALL ptime(t2)
    t3=t2

    i97 = hecMAT%Iarray(97)
    IF(hecMAT%Iarray(97)==1) THEN
      !* Interface to numeric factorization
      CALL nuform(hecMESH,hecMAT,ir)
      CALL ptime(t3)

      !* Numeric factorization
      CALL nufct0(ir)
      hecMAT%Iarray(97) = 0

      !*Memory Details
      WRITE(*,*) '*-----------------------------------*'
      WRITE(*,*) '|   Direct  Solver  Memory  Usage   |'
      WRITE(*,*) '*-----------------------------------*'
      WRITE(*,*) 'INTEGER memory: ',REAL(ialoc*4)/REAL(1048576),'MB'
      WRITE(*,*) 'REAL*8  memory: ',REAL(raloc*8)/REAL(1048576),'MB'
      WRITE(*,*) 'TOTAL   memory: ',REAL((raloc*2+ialoc)*4)/REAL(1048576),'MB'
      WRITE(*,*) '*-----------------------------------*'
    ENDIF
    CALL ptime(t4)

    !* Finalize
    !      tom(5)=tt2-tt1
    !*  Errors 1
    IF(i98/=0.AND.i98/=1) THEN
      WRITE(IFMSG,*) 'ERROR in symb. fact. flag: Should be 1 or 0'
      STOP 'ERROR in symb. fact. flag: Should be 1 or 0'
    ENDIF
    IF(i97/=0.AND.i97/=1) THEN
      WRITE(IFMSG,*) 'ERROR in numer. fact. flag: Should be 1 or 0'
      STOP 'ERROR in numer. fact. flag: Should be 1 or 0'
    ENDIF
    IF(i98==1.AND.i97==0) THEN
      WRITE(IFMSG,*) 'WARNING: Numeric factorization not performed!'
      STOP 'WARNING: Numeric factorization not performed! Solve will not be performed'
      RETURN
    ENDIF
    !*  Errors 2
    IF(ir/=0) THEN
      !WINDEBUG
      WRITE(IFMSG,*) 'ERROR in nufct0. ir = ',ir
      STOP
    ENDIF

    tom(1)=t2-t1
    tom(2)=t3-t2
    tom(3)=t4-t3

    !* Solve
    !* Backsubstitute
    CALL nusol0(hecMAT%B,ir)
    CALL ptime(t5)
    !* Errors 4
    IF(ir/=0) THEN
      !WINDEBUG
      WRITE(IFMSG,*) 'error in nusol0. irr = ',ir
      STOP
    ENDIF
    CALL hecmw_mat_dump_solution(hecMAT)
    RETURN
  END SUBROUTINE hecmw_solve_direct

  !======================================================================!
  !> ddot
  !======================================================================!
  DOUBLE PRECISION FUNCTION ddot(a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n):: b

    INTEGER:: i
    DOUBLE PRECISION:: s

    s=0.0d0
    DO i=1,n
      s=s+a(i)*b(i)
    ENDDO
    ddot=s
    RETURN
  END FUNCTION ddot

  !======================================================================!
  !> spdot2
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  DOUBLE PRECISION FUNCTION spdot2(b,ks,ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ke
    INTEGER, INTENT(IN):: ks
    DOUBLE PRECISION, INTENT(IN), DIMENSION(*):: b

    INTEGER:: j
    INTEGER:: jj
    DOUBLE PRECISION:: s

    s=0.0d0
    DO jj=ks,ke
      j=colno(jj)
      s=s+zln(jj)*b(j)
    ENDDO
    spdot2=s
  END FUNCTION spdot2

  !======================================================================!
  !> addr0
  !======================================================================!
  SUBROUTINE addr0(isw,i,j,aij,diag,zln,dsln,nstop,ndeg2,ndeg2l,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: i
    INTEGER, INTENT(IN):: isw
    INTEGER, INTENT(IN):: j
    INTEGER, INTENT(IN):: ndeg2
    INTEGER, INTENT(IN):: ndeg2l
    INTEGER, INTENT(IN):: nstop
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg2):: aij

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg2,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg2l,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg2,*):: dsln

    INTEGER:: i0
    INTEGER:: ii
    INTEGER:: itrans
    INTEGER:: j0
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    !INTEGER:: l
    INTEGER, PARAMETER:: idbg = 0

    ir=0
    ii=invp(i)
    jj=invp(j)
    IF(idbg/=0) WRITE(6,*) ii,jj,aij
    IF(ii==jj) THEN
      IF(ndeg2==1) THEN
        IF(isw==0) THEN
          diag(1,ii)=aij(1)
        ELSE
          diag(1,ii)=diag(1,ii)+aij(1)
        ENDIF
      ELSEIF(ndeg2==4) THEN
        IF(isw==0) THEN
          diag(1,ii)=aij(1)
          diag(2,ii)=aij(2)
          diag(3,ii)=aij(4)
        ELSE
          diag(1,ii)=diag(1,ii)+aij(1)
          diag(2,ii)=diag(2,ii)+aij(2)
          diag(3,ii)=diag(3,ii)+aij(4)
        ENDIF
      ENDIF
      RETURN
    ENDIF
    itrans=0
    IF(jj>ii) THEN
      k=jj
      jj=ii
      ii=k
      itrans=1
    ENDIF
    IF(jj>=nstop) THEN
      i0=ii-nstop
      j0=jj-nstop+1
      k=i0*(i0-1)/2+j0
      IF(ndeg2==1) THEN
        dsln(1,k)=aij(1)
        RETURN
      ELSEIF(ndeg2==4) THEN
        IF(itrans==0) THEN
          dsln(1:ndeg2,k)=aij(1:ndeg2)
          RETURN
        ELSE
          dsln(1,k)=aij(1)
          dsln(2,k)=aij(3)
          dsln(3,k)=aij(2)
          dsln(4,k)=aij(4)
          RETURN
        ENDIF
      ENDIF
    ENDIF
    ks=xlnzr(ii)
    ke=xlnzr(ii+1)-1
    DO k=ks,ke
      IF(colno(k)==jj) THEN
        IF(isw==0) THEN
          IF(ndeg2==1) THEN
            zln(1,k)=aij(1)
          ELSEIF(ndeg2==4) THEN
            IF(itrans==0) THEN
              zln(1:ndeg2,k)=aij(1:ndeg2)
            ELSE
              zln(1:4,k)=aij(1:4)
            ENDIF
          ENDIF
        ELSE
          IF(ndeg2==1) THEN
            zln(1,k)=zln(1,k)+aij(1)
          ELSEIF(ndeg2==4) THEN
            IF(itrans==0) THEN
              zln(1:ndeg2,k)=zln(1:ndeg2,k)+aij(1:ndeg2)
            ELSE
              zln(1:4,k)=zln(1:4,k)+aij(1:4)
            ENDIF
          ENDIF
        ENDIF
        RETURN
      ENDIF
    ENDDO
    ir=20
    RETURN
  END SUBROUTINE addr0

  !======================================================================!
  !> addr3
  !======================================================================!
  SUBROUTINE addr3(i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: i
    INTEGER, INTENT(IN):: j
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: invp
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9):: aij

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: dsln

    INTEGER:: i0
    INTEGER:: ii
    INTEGER:: itrans
    INTEGER:: j0
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER, PARAMETER:: idbg = 0
    INTEGER, PARAMETER:: ndeg2 = 9

    ir=0
    ii=invp(i)
    jj=invp(j)
    IF(idbg/=0) WRITE(6,*) ii,jj,aij
    IF(ii==jj) THEN
      diag(1,ii)=aij(1)
      diag(2,ii)=aij(2)
      diag(3,ii)=aij(5)
      diag(4,ii)=aij(3)
      diag(5,ii)=aij(6)
      diag(6,ii)=aij(9)
      RETURN
    ENDIF
    itrans=0
    IF(jj>ii) THEN
      k=jj
      jj=ii
      ii=k
      itrans=1
    ENDIF
    IF(jj>=nstop) THEN
      i0=ii-nstop
      j0=jj-nstop+1
      k=i0*(i0-1)/2+j0
      IF(itrans==0) THEN
        dsln(1:ndeg2,k)=aij(1:ndeg2)
        !DO l=1,ndeg2
        !  dsln(l,k)=aij(l)
        !ENDDO
        RETURN
      ELSE
        dsln(1,k)=aij(1)
        dsln(2,k)=aij(4)
        dsln(3,k)=aij(7)
        dsln(4,k)=aij(2)
        dsln(5,k)=aij(5)
        dsln(6,k)=aij(8)
        dsln(7,k)=aij(3)
        dsln(8,k)=aij(6)
        dsln(9,k)=aij(9)
        RETURN
      ENDIF
    ENDIF
    ks=xlnzr(ii)
    ke=xlnzr(ii+1)-1
    DO k=ks,ke
      IF(colno(k)==jj) THEN
        IF(itrans==0) THEN
          zln(1:ndeg2,k)=aij(1:ndeg2)
        ELSE
          zln(1,k)=aij(1)
          zln(2,k)=aij(4)
          zln(3,k)=aij(7)
          zln(4,k)=aij(2)
          zln(5,k)=aij(5)
          zln(6,k)=aij(8)
          zln(7,k)=aij(3)
          zln(8,k)=aij(6)
          zln(9,k)=aij(9)
        ENDIF
        RETURN
      ENDIF
    ENDDO
    ir=20

    RETURN
  END SUBROUTINE addr3

  !======================================================================!
  !> addrx
  !======================================================================!
  SUBROUTINE addrx(i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ndeg2l,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: i
    INTEGER, INTENT(IN):: j
    INTEGER, INTENT(IN):: ndeg2l
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: invp
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg):: aij

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg2l,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: zln

    INTEGER:: i0
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
    INTEGER, PARAMETER:: idbg = 0

    ir=0
    ii=invp(i)
    jj=invp(j)
    IF(idbg/=0) WRITE(6,*) ii,jj,aij
    IF(ii==jj) THEN
      l=0
      DO n=1,ndeg
        DO m=1,n
          l=l+1
          diag(l,ii)=aij(n,m)
        ENDDO
      ENDDO
      RETURN
    ENDIF
    itrans=0
    IF(jj>ii) THEN
      k=jj
      jj=ii
      ii=k
      itrans=1
    ENDIF
    IF(jj>=nstop) THEN
      i0=ii-nstop
      j0=jj-nstop+1
      k=i0*(i0-1)/2+j0
      IF(itrans==0) THEN
        dsln(1:ndeg,1:ndeg,k)=aij(1:ndeg,1:ndeg)
        RETURN
      ELSE
        dsln(1:ndeg,1:ndeg,k)=aij(1:ndeg,1:ndeg)
        RETURN
      ENDIF
    ENDIF
    ks=xlnzr(ii)
    ke=xlnzr(ii+1)-1
    DO k=ks,ke
      IF(colno(k)==jj) THEN
        IF(itrans==0) THEN
          zln(1:ndeg,1:ndeg,k)=aij(1:ndeg,1:ndeg)
        ELSE
          zln(1:ndeg,1:ndeg,k)=aij(1:ndeg,1:ndeg)
        ENDIF
        RETURN
      ENDIF
    ENDDO
    ir=20
    RETURN
  END SUBROUTINE addrx

  !======================================================================!
  !> @brief bringu brings up zero pivots from bottom of the elimination tree to higher nodes
  !      irr = 0     complete
  !          = 1     impossible
  !======================================================================!
  SUBROUTINE bringu(zpiv,iperm,invp,parent,izz,neqns,irr)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: izz
    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: zpiv
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: irr
    INTEGER, INTENT(OUT), DIMENSION(*):: iperm
    INTEGER, INTENT(OUT), DIMENSION(*):: invp

    INTEGER:: i
    INTEGER:: ib
    INTEGER:: ib0
    INTEGER:: ibp
    INTEGER:: idbg
    INTEGER:: izzp

    idbg=0
    irr=0
    ib0=invp(izz)
    ib=ib0
    DO
      IF(ib<=0) THEN
        irr=1
        RETURN
      ENDIF
      ibp=parent(ib)
      izzp=iperm(ibp)
      IF(zpiv(izzp)==0) EXIT
      ib=ibp
    ENDDO

    invp(izz)=ibp
    invp(izzp)=ib0
    iperm(ibp)=izz
    iperm(ib0)=izzp
    IF(idbg/=0) THEN
      DO i=1,neqns
        IF((invp(iperm(i))/=i).OR.(iperm(invp(i))/=i)) THEN
          WRITE(6,*) 'permutation error'
          STOP
        ENDIF
      ENDDO
      RETURN
    ENDIF
    RETURN
  END SUBROUTINE bringu

  !======================================================================!
  !> d2dot
  !      spdot performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE d2dot(t,a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4):: t

    INTEGER:: jj

    t(1:4)=0.0d0
    DO jj=1,n
      t(1)=t(1)+a(1,jj)*b(1,jj)+a(3,jj)*b(3,jj)
      t(2)=t(2)+a(2,jj)*b(1,jj)+a(4,jj)*b(3,jj)
      t(3)=t(3)+a(1,jj)*b(2,jj)+a(3,jj)*b(4,jj)
      t(4)=t(4)+a(2,jj)*b(2,jj)+a(4,jj)*b(4,jj)
    ENDDO
    RETURN
  END SUBROUTINE d2dot

  !======================================================================!
  !> d2sdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE d2sdot(wi,a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(2,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2):: wi

    INTEGER:: jj

    DO jj=1,n
      wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(3,jj)
      wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(4,jj)
    ENDDO
    RETURN
  END SUBROUTINE d2sdot

  !======================================================================!
  !> d3dot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE d3dot(t,a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9):: t

    INTEGER:: jj
    !INTEGER:: l

    !$dir max_trips(9)
    t(1:9)=0.0d0
    DO jj=1,n
      t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
      t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
      t(3)=t(3)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
      t(4)=t(4)+a(1,jj)*b(2,jj)+a(4,jj)*b(5,jj)+a(7,jj)*b(8,jj)
      t(5)=t(5)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
      t(6)=t(6)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
      t(7)=t(7)+a(1,jj)*b(3,jj)+a(4,jj)*b(6,jj)+a(7,jj)*b(9,jj)
      t(8)=t(8)+a(2,jj)*b(3,jj)+a(5,jj)*b(6,jj)+a(8,jj)*b(9,jj)
      t(9)=t(9)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
    ENDDO
    RETURN
  END SUBROUTINE d3dot

  !======================================================================!
  !> d3dotl
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE d3dotl(t,a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6):: t

    INTEGER:: jj

    t(1:6)=0.0d0

    DO jj=1,n
      t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
      t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
      t(3)=t(3)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
      t(4)=t(4)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
      t(5)=t(5)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
      t(6)=t(6)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
    ENDDO
    RETURN
  END SUBROUTINE d3dotl

  !======================================================================!
  !> d3sdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE d3sdot(wi,a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3):: wi

    INTEGER:: jj

    DO jj=1,n
      wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(4,jj)-a(3,jj)*b(7,jj)
      wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(5,jj)-a(3,jj)*b(8,jj)
      wi(3)=wi(3)-a(1,jj)*b(3,jj)-a(2,jj)*b(6,jj)-a(3,jj)*b(9,jj)
    ENDDO
    RETURN
  END SUBROUTINE d3sdot

  !======================================================================!
  !> dxdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE dxdot(t,a,b,l)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg):: t

    INTEGER:: jj
    INTEGER:: k
    INTEGER:: l
    INTEGER:: m
    INTEGER:: n

    DO n=1,ndeg
      DO  m=1,ndeg
        t(n,m)=0.0d0
        DO k=1,ndeg
          DO jj=1,l
            t(n,m)=t(n,m)+a(n,k,jj)*b(m,k,jj)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE dxdot

  !======================================================================!
  !> dxdotl
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE dxdotl(t,a,b,l)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg):: t

    INTEGER:: jj
    INTEGER:: k
    INTEGER:: l
    INTEGER:: m
    INTEGER:: n

    DO n=1,ndeg
      DO m=1,n
        t(n,m)=0.0d0
        DO k=1,ndeg
          DO jj=1,l
            t(n,m)=t(n,m)+a(n,k,jj)*b(m,k,jj)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE dxdotl

  !======================================================================!
  !> dxsdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE dxsdot(wi,a,b,n)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,*):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg):: wi

    INTEGER:: jj
    INTEGER:: m
    INTEGER:: n

    DO jj=1,n
      DO m=1,ndeg
        wi(1:ndeg)=wi(1:ndeg)-b(1:ndeg,m,jj)*a(m,jj)
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE dxsdot

  !======================================================================!
  !> forpar
  !======================================================================!
  SUBROUTINE forpar(neqns,parent,nch,nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: nstop
    INTEGER, INTENT(OUT), DIMENSION(*):: nch

    INTEGER:: i
    INTEGER:: idens
    INTEGER:: ii

    nch(1:neqns)=0
    nch(neqns+1)=0
    DO i=1,neqns
      ii=parent(i)
      nch(ii)=nch(ii)+1
    ENDDO
    DO i=neqns,1,-1
      IF(nch(i)/=1) EXIT
    ENDDO

    idens = 0
    IF(idens==1) THEN
      nstop=i
    ELSE
      nstop=neqns+1
    ENDIF
    RETURN
  END SUBROUTINE forpar

  !======================================================================!
  !> genbtq
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE genbtq(invp,parent,btree,zpiv,izz,neqns)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: parent
    INTEGER, INTENT(IN), DIMENSION(*):: invp
    INTEGER, INTENT(IN), DIMENSION(*):: zpiv

    INTEGER, INTENT(OUT):: izz
    INTEGER, INTENT(OUT), DIMENSION(2,*):: btree

    INTEGER:: i
    INTEGER:: ib
    INTEGER:: idbg1
    INTEGER:: inext
    INTEGER:: ip
    COMMON /debug/ idbg1

    btree(1,1:neqns+1)=0
    btree(2,1:neqns+1)=0
    DO i=1,neqns+1
      ip=parent(i)
      IF(ip<=0) CYCLE
      ib=btree(1,ip)
      IF(ib==0) THEN
        btree(1,ip)=i
      ELSE
        DO
          inext=btree(2,ib)
          IF(inext/=0) THEN
            ib=inext
          ELSE
            btree(2,ib)=i
            EXIT
          ENDIF
        ENDDO

      ENDIF
    ENDDO
    !
    ! find zeropivot
    !
    DO i=1,neqns
      IF(zpiv(i)/=0) THEN
        IF(btree(1,invp(i))==0) THEN
          izz=i
          GOTO 210
        ENDIF
      ENDIF
    ENDDO
    izz=0
210 CONTINUE
    IF(idbg1/=0) WRITE(6,"(' binary tree')")
    IF(idbg1/=0) WRITE(6,"(i6,'(',2i6,')')") (i,btree(1,i),btree(2,i),i=1,neqns)
    IF(idbg1/=0) WRITE(6,"(' the first zero pivot is ',i4)") izz
    RETURN
  END SUBROUTINE genbtq

  !======================================================================!
  !> genpaq
  !     #oded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE genpaq(xadj,adjncy,invp,iperm,parent,neqns,ancstr)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: xadj
    INTEGER, INTENT(IN), DIMENSION(*):: adjncy
    INTEGER, INTENT(IN), DIMENSION(*):: invp
    INTEGER, INTENT(IN), DIMENSION(*):: iperm

    INTEGER, INTENT(OUT), DIMENSION(*):: parent
    INTEGER, INTENT(OUT), DIMENSION(*):: ancstr

    INTEGER:: i
    INTEGER:: idbg1
    INTEGER:: ip
    INTEGER:: it
    INTEGER:: k
    INTEGER:: l
    COMMON /debug/ idbg1

    DO i=1,neqns
      parent(i)=0
      ancstr(i)=0
      ip=iperm(i)
      DO k=xadj(ip),xadj(ip+1)-1
        l=invp(adjncy(k))
        IF(l>=i) GOTO 110
        DO
          IF(ancstr(l)==0) EXIT
          IF(ancstr(l)==i) GOTO 110
          it=ancstr(l)
          ancstr(l)=i
          l=it
        ENDDO
        ancstr(l)=i
        parent(l)=i
110     CONTINUE
      ENDDO
    ENDDO
    DO i=1,neqns
      IF(parent(i)==0) parent(i)=neqns+1
    ENDDO
    parent(neqns+1)=0
    IF(idbg1/=0) WRITE(6,"(' parent')")
    IF(idbg1/=0) WRITE(6,"(2i6)") (i,parent(i),i=1,neqns)
    RETURN
  END SUBROUTINE genpaq

  !======================================================================!
  !> genqmd
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE genqmd(neqns,xadj,adj0,perm,invp,deg,marker,rchset,nbrhd,qsize,qlink,nofsub,adjncy)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: adj0
    INTEGER, INTENT(IN), DIMENSION(*):: xadj

    INTEGER, INTENT(OUT):: nofsub
    INTEGER, INTENT(OUT), DIMENSION(*):: rchset
    INTEGER, INTENT(OUT), DIMENSION(*):: nbrhd
    INTEGER, INTENT(OUT), DIMENSION(*):: adjncy
    INTEGER, INTENT(OUT), DIMENSION(*):: perm
    INTEGER, INTENT(OUT), DIMENSION(*):: invp
    INTEGER, INTENT(OUT), DIMENSION(*):: deg
    INTEGER, INTENT(OUT), DIMENSION(*):: marker
    INTEGER, INTENT(OUT), DIMENSION(*):: qsize
    INTEGER, INTENT(OUT), DIMENSION(*):: qlink

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

    mindeg=neqns
    nofsub=0
    adjncy(1:xadj(neqns+1)-1)=adj0(1:xadj(neqns+1)-1)

    DO node=1,neqns
      perm(node)=node
      invp(node)=node
      marker(node)=0
      qsize(node)=1
      qlink(node)=0
      ndeg=xadj(node+1)-xadj(node)
      deg(node)=ndeg
      IF(ndeg<mindeg) mindeg=ndeg
    ENDDO

    num=0
    DO
      search=1
      thresh=mindeg
      mindeg=neqns

300   CONTINUE
      nump1=num+1
      IF(nump1>search) search=nump1
      DO j=search,neqns
        node=perm(j)
        IF(marker(node)<0) CYCLE
        ndeg=deg(node)
        IF(ndeg<=thresh) GOTO 500
        IF(ndeg<mindeg) mindeg=ndeg
      ENDDO
    ENDDO

500 CONTINUE
    search=j
    nofsub=nofsub+deg(node)
    marker(node)=1
    CALL qmdrch(node,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze,nbrhd)
    nxnode=node

    DO
      num=num+1
      np=invp(nxnode)
      ip=perm(num)
      perm(np)=ip
      invp(ip)=np
      perm(num)=nxnode
      invp(nxnode)=num
      deg(nxnode)=-1
      nxnode=qlink(nxnode)
      IF(nxnode<=0) EXIT
    ENDDO

    IF(rchsze<=0) GOTO 800
    !
    CALL qmdupd(xadj,adjncy,rchsze,rchset,deg,qsize,qlink,marker,rchset(rchsze+1),nbrhd(nhdsze+1))
    marker(node)=0
    DO irch=1,rchsze
      inode=rchset(irch)
      IF(marker(inode)<0) CYCLE
      marker(inode)=0
      ndeg=deg(inode)
      IF(ndeg<mindeg) mindeg=ndeg
      IF(ndeg>thresh) CYCLE
      mindeg=thresh
      thresh=ndeg
      search=invp(inode)
    ENDDO
    IF(nhdsze>0) CALL qmdot(node,xadj,adjncy,marker,rchsze,rchset,nbrhd)

800 CONTINUE
    IF(num<neqns) GOTO 300

    RETURN
  END SUBROUTINE genqmd

  !======================================================================!
  !> gnclno
  !======================================================================!
  SUBROUTINE gnclno(parent,xleaf,leaf,xlnzr,colno,neqns,nstop,lncol,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: parent
    INTEGER, INTENT(IN), DIMENSION(*):: xleaf
    INTEGER, INTENT(IN), DIMENSION(*):: leaf

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT):: lncol
    INTEGER, INTENT(OUT), DIMENSION(*):: xlnzr
    INTEGER, INTENT(OUT), DIMENSION(*):: colno

    INTEGER:: i
    INTEGER:: idbg1
    INTEGER:: j
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: nc
    INTEGER:: nxleaf
    COMMON /debug/ idbg1

    nc=0
    ir=0
    l=1
    DO i=1,neqns
      xlnzr(i)=l
      ks=xleaf(i)
      ke=xleaf(i+1)-1
      IF(ke<ks) CYCLE
      nxleaf=leaf(ks)
      DO k=ks,ke-1
        j=nxleaf
        nxleaf=leaf(k+1)
        DO
          IF(j>=nxleaf) EXIT
          IF(j>=nstop) THEN
            GOTO 100
          ENDIF
          colno(l)=j
          l=l+1
          j=parent(j)
        ENDDO
      ENDDO
      j=leaf(ke)

      DO
        IF(j>=nstop) GOTO 100
        IF(j>=i.OR.j==0) GOTO 100
        colno(l)=j
        l=l+1
        j=parent(j)
      ENDDO
100   CONTINUE
    ENDDO

    xlnzr(neqns+1)=l
    lncol=l-1
    IF(idbg1/=0) WRITE(6,"(' xlnzr')")
    IF(idbg1/=0) WRITE(6,"(' colno (lncol =',i10,')')") lncol
    IF(idbg1/=0) THEN
      DO k=1,neqns
        WRITE(6,"(/' row = ',i6)") k
        WRITE(6,"(10i4)") (colno(i),i=xlnzr(k),xlnzr(k+1)-1)
      ENDDO
    ENDIF
    RETURN
  END SUBROUTINE gnclno

  !======================================================================!
  !> gnleaf
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE gnleaf(xadj,adjncy,invp,iperm,nch,adjncp,xleaf,leaf,neqns,lnleaf)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: xadj
    INTEGER, INTENT(IN), DIMENSION(*):: adjncy
    INTEGER, INTENT(IN), DIMENSION(*):: nch
    INTEGER, INTENT(IN), DIMENSION(*):: invp
    INTEGER, INTENT(IN), DIMENSION(*):: iperm
    INTEGER, INTENT(OUT):: lnleaf
    INTEGER, INTENT(OUT), DIMENSION(*):: adjncp
    INTEGER, INTENT(OUT), DIMENSION(*):: xleaf
    INTEGER, INTENT(OUT), DIMENSION(*):: leaf

    INTEGER:: i
    INTEGER:: idbg1
    INTEGER:: ik
    INTEGER:: ip
    INTEGER:: iq
    INTEGER:: istart
    INTEGER:: k
    INTEGER:: l
    INTEGER:: lc
    INTEGER:: lc1
    INTEGER:: m
    COMMON /debug/ idbg1

    l=1
    ik=0
    istart=0
    DO i=1,neqns
      xleaf(i)=l
      ip=iperm(i)
      DO k=xadj(ip),xadj(ip+1)-1
        iq=invp(adjncy(k))
        IF(iq<i) THEN
          ik=ik+1
          adjncp(ik)=iq
        ENDIF
      ENDDO
      m=ik-istart

      IF(m==0) CYCLE
      CALL qqsort(adjncp(istart+1),m)
      lc1=adjncp(istart+1)
      IF(lc1>=i) CYCLE
      leaf(l)=lc1
      l=l+1
      DO k=istart+2,ik
        lc=adjncp(k)
        IF(lc1<lc-nch(lc)) THEN
          leaf(l)=lc
          l=l+1
        ENDIF

        lc1=lc
      ENDDO
      ik=1
      istart=ik

    ENDDO
    xleaf(neqns+1)=l
    lnleaf=l-1
    IF(idbg1/=0) WRITE(6,"(' xleaf')")
    IF(idbg1/=0) WRITE(6,"(10i6)") (xleaf(i),i=1,neqns+1)
    IF(idbg1/=0) WRITE(6,"(' leaf (len = ',i6,')')") lnleaf
    IF(idbg1/=0) WRITE(6,"(10i6)") (leaf(i),i=1,lnleaf)
    RETURN
  END SUBROUTINE gnleaf

  !======================================================================!
  !> inv2
  !======================================================================!
  SUBROUTINE inv2(dsln,ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3):: dsln

    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ir=0
    IF(dabs(dsln(1))<rmin) THEN
      ir=10
      RETURN
    ENDIF
    dsln(1)=1.0d0/dsln(1)
    t=dsln(2)*dsln(1)
    dsln(3)=dsln(3)-t*dsln(2)
    dsln(2)=t
    IF(dabs(dsln(3))<rmin) THEN
      ir=10
      RETURN
    ENDIF
    dsln(3)=1.0d0/dsln(3)
    RETURN
  END SUBROUTINE inv2

  !======================================================================!
  !> inv22
  !======================================================================!
  SUBROUTINE inv22(zln,zz,diag)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(3):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4):: zz

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4):: zln

    zln(3)=zz(3)-zz(1)*diag(2)
    zln(1)=zz(1)*diag(1)
    zln(3)=zln(3)*diag(3)
    zln(1)=zln(1)-zln(3)*diag(2)

    zln(4)=zz(4)-zz(2)*diag(2)
    zln(2)=zz(2)*diag(1)
    zln(4)=zln(4)*diag(3)
    zln(2)=zln(2)-zln(4)*diag(2)

    RETURN
  END SUBROUTINE inv22

  !======================================================================!
  !> inv3
  !======================================================================!
  SUBROUTINE inv3(dsln,ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6):: dsln

    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(2):: t
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ir=0

    IF(dabs(dsln(1))<rmin) GOTO 1000

    dsln(1)=1.0d0/dsln(1)
    t(1)=dsln(2)*dsln(1)
    dsln(3)=dsln(3)-t(1)*dsln(2)
    dsln(2)=t(1)

    IF(dabs(dsln(3))<rmin) GOTO 1000

    dsln(3)=1.0d0/dsln(3)
    t(1)=dsln(4)*dsln(1)
    dsln(5)=dsln(5)-dsln(2)*dsln(4)
    t(2)=dsln(5)*dsln(3)
    dsln(6)=dsln(6)-t(1)*dsln(4)-t(2)*dsln(5)
    dsln(4)=t(1)
    dsln(5)=t(2)

    IF(dabs(dsln(6))<rmin) GOTO 1000

    dsln(6)=1.0d0/dsln(6)
    RETURN

1000 CONTINUE
    dsln(1)=1.0d0
    dsln(2)=0.0d0
    dsln(3)=1.0d0
    dsln(4)=0.0d0
    dsln(5)=0.0d0
    dsln(6)=1.0d0
    RETURN
  END SUBROUTINE inv3

  !======================================================================!
  !> inv33
  !======================================================================!
  SUBROUTINE inv33(zln,zz,diag)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(6):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9):: zz

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9):: zln

    zln(4)=zz(4)-zz(1)*diag(2)
    zln(7)=zz(7)-zz(1)*diag(4)-zln(4)*diag(5)
    zln(1)=zz(1)*diag(1)
    zln(4)=zln(4)*diag(3)
    zln(7)=zln(7)*diag(6)
    zln(4)=zln(4)-zln(7)*diag(5)
    zln(1)=zln(1)-zln(4)*diag(2)-zln(7)*diag(4)

    zln(5)=zz(5)-zz(2)*diag(2)
    zln(8)=zz(8)-zz(2)*diag(4)-zln(5)*diag(5)
    zln(2)=zz(2)*diag(1)
    zln(5)=zln(5)*diag(3)
    zln(8)=zln(8)*diag(6)
    zln(5)=zln(5)-zln(8)*diag(5)
    zln(2)=zln(2)-zln(5)*diag(2)-zln(8)*diag(4)

    zln(6)=zz(6)-zz(3)*diag(2)
    zln(9)=zz(9)-zz(3)*diag(4)-zln(6)*diag(5)
    zln(3)=zz(3)*diag(1)
    zln(6)=zln(6)*diag(3)
    zln(9)=zln(9)*diag(6)
    zln(6)=zln(6)-zln(9)*diag(5)
    zln(3)=zln(3)-zln(6)*diag(2)-zln(9)*diag(4)
    RETURN
  END SUBROUTINE inv33

  !======================================================================!
  !> inv6
  !======================================================================!
  SUBROUTINE inv6(dsln,ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(21):: dsln

    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(5):: t
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ir=0
    dsln(1)=1.0d0/dsln(1)
    t(1)=dsln(2)*dsln(1)
    dsln(3)=1.0d0/(dsln(3)-t(1)*dsln(2))
    dsln(2)=t(1)
    dsln(5)=dsln(5)-dsln(4)*dsln(2)
    t(1)=dsln(4)*dsln(1)
    t(2)=dsln(5)*dsln(3)
    dsln(6)=1.0d0/(dsln(6)-t(1)*dsln(4)-t(2)*dsln(5))
    dsln(4)=t(1)
    dsln(5)=t(2)
    dsln(8)=dsln(8)-dsln(7)*dsln(2)
    dsln(9)=dsln(9)-dsln(7)*dsln(4)-dsln(8)*dsln(5)
    t(1)=dsln(7)*dsln(1)
    t(2)=dsln(8)*dsln(3)
    t(3)=dsln(9)*dsln(6)
    dsln(10)=1.0d0  /(dsln(10)-t(1)*dsln(7)-t(2)*dsln(8)-t(3)*dsln(9))
    dsln(7)=t(1)
    dsln(8)=t(2)
    dsln(9)=t(3)
    dsln(12)=dsln(12)-dsln(11)*dsln(2)
    dsln(13)=dsln(13)-dsln(11)*dsln(4)-dsln(12)*dsln(5)
    dsln(14)=dsln(14)-dsln(11)*dsln(7)-dsln(12)*dsln(8)-dsln(13)*dsln(9)
    t(1)=dsln(11)*dsln(1)
    t(2)=dsln(12)*dsln(3)
    t(3)=dsln(13)*dsln(6)
    t(4)=dsln(14)*dsln(10)
    dsln(15)=1.0d0/(dsln(15)-t(1)*dsln(11)-t(2)*dsln(12)-t(3)*dsln(13)-t(4)*dsln(14))
    dsln(11)=t(1)
    dsln(12)=t(2)
    dsln(13)=t(3)
    dsln(14)=t(4)
    dsln(17)=dsln(17)-dsln(16)*dsln(2)
    dsln(18)=dsln(18)-dsln(16)*dsln(4)-dsln(17)*dsln(5)
    dsln(19)=dsln(19)-dsln(16)*dsln(7)-dsln(17)*dsln(8)-dsln(18)*dsln(9)
    dsln(20)=dsln(20)-dsln(16)*dsln(11)-dsln(17)*dsln(12)-dsln(18)*dsln(13)-dsln(19)*dsln(14)
    t(1)=dsln(16)*dsln(1)
    t(2)=dsln(17)*dsln(3)
    t(3)=dsln(18)*dsln(6)
    t(4)=dsln(19)*dsln(10)
    t(5)=dsln(20)*dsln(15)
    dsln(21)=1.0d0/(dsln(21)-t(1)*dsln(16)-t(2)*dsln(17)-t(3)*dsln(18)-t(4)*dsln(19)-t(5)*dsln(20))
    dsln(16)=t(1)
    dsln(17)=t(2)
    dsln(18)=t(3)
    dsln(19)=t(4)
    dsln(20)=t(5)
    RETURN
  END SUBROUTINE inv6

  !======================================================================!
  !> inv66
  !======================================================================!
  SUBROUTINE inv66(zln,zz,diag)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(21):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(36):: zz

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36):: zln

    INTEGER:: i

    DO i=0,5
      zln(i+7)=zz(i+7)-zz(i+1)*diag(2)
      zln(i+13)=zz(i+13)-zz(i+1)*diag(4)-zln(i+7)*diag(5)
      zln(i+19)=zz(i+19)-zz(i+1)*diag(7)-zln(i+7)*diag(8)-zln(i+13)*diag(9)
      zln(i+25)=zz(i+25)-zz(i+1)*diag(11)-zln(i+7)*diag(12)-zln(i+13)*diag(13)-zln(i+19)*diag(14)
      zln(i+31)=zz(i+31)-zz(i+1)*diag(16)-zln(i+7)*diag(17)-zln(i+13)*diag(18)-zln(i+19)*diag(19)-zln(i+25)*diag(20)
      zln(i+1)=zz(i+1)*diag(1)
      zln(i+7)=zln(i+7)*diag(3)
      zln(i+13)=zln(i+13)*diag(6)
      zln(i+19)=zln(i+19)*diag(10)
      zln(i+25)=zln(i+25)*diag(15)
      zln(i+31)=zln(i+31)*diag(21)
      zln(i+25)=zln(i+25)-zln(i+31)*diag(20)
      zln(i+19)=zln(i+19)-zln(i+31)*diag(19)-zln(i+25)*diag(14)
      zln(i+13)=zln(i+13)-zln(i+31)*diag(18)-zln(i+25)*diag(13)-zln(i+19)*diag(9)
      zln(i+7)=zln(i+7)-zln(i+31)*diag(17)-zln(i+25)*diag(12)-zln(i+19)*diag(8)-zln(i+13)*diag(5)
      zln(i+1)=zln(i+1)-zln(i+31)*diag(16)-zln(i+25)*diag(11)-zln(i+19)*diag(7)-zln(i+13)*diag(4)-zln(i+7)*diag(2)
    ENDDO
    RETURN
  END SUBROUTINE inv66

  !======================================================================!
  !> invx
  !======================================================================!
  SUBROUTINE invx(dsln,ndeg,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ndeg

    INTEGER, INTENT(OUT):: ir
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: dsln

    INTEGER:: i
    INTEGER:: j
    INTEGER:: k
    INTEGER:: k0
    INTEGER:: l
    INTEGER:: l0
    INTEGER:: ld
    INTEGER:: ll
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t
    DOUBLE PRECISION:: tem
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ir=0
    l=1
    dsln(1)=1.0d0/dsln(1)
    DO i=2,ndeg
      ld=0
      l0=l
      DO j=1,i-1
        l=l+1
        DO k=1,j-1
          ld=ld+1
          dsln(l)=dsln(l)-dsln(l0+k)*dsln(ld)
        ENDDO
        ld=ld+1
      ENDDO
      t=0.0d0
      k0=0
      ll=0
      DO k=l-i+2,l
        ll=ll+1
        k0=k0+ll
        tem=dsln(k)*dsln(k0)
        t=t+tem*dsln(k)
        dsln(k)=tem
      ENDDO
      l=l+1
      dsln(l)=dsln(l)-t
      dsln(l)=1.0d0/dsln(l)
    ENDDO
    RETURN
  END SUBROUTINE invx

  !======================================================================!
  !> invxx
  !======================================================================!
  SUBROUTINE invxx(zln,zz,diag,ndeg)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ndeg
    DOUBLE PRECISION, INTENT(IN), DIMENSION(*):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg):: zz

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg):: zln

    INTEGER:: joc
    INTEGER:: l
    INTEGER:: loc1
    INTEGER:: m
    INTEGER:: n

    zln=zz
    DO l=1,ndeg,2
      joc=0
      DO m=1,ndeg-1
        joc=joc+m
        loc1=joc+m
        DO n=m+1,ndeg
          zln(l,n)=zln(l,n)-zln(l,m)*diag(loc1)
          zln(l+1,n)=zln(l+1,n)-zln(l+1,m)*diag(loc1)
          loc1=loc1+n
        ENDDO
      ENDDO
      joc=0
      DO m=1,ndeg
        joc=joc+m
        zln(l,m)=zln(l,m)*diag(joc)
        zln(l+1,m)=zln(l+1,m)*diag(joc)
      ENDDO
      DO n=ndeg,2,-1
        joc=joc-1
        DO m=n-1,1,-1
          zln(l,m)=zln(l,m)-zln(l,n)*diag(joc)
          zln(l+1,m)=zln(l+1,m)-zln(l+1,n)*diag(joc)
          joc=joc-1
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE invxx

  !======================================================================!
  !> matini
  !     matini initializes storage for sparse matrix solver.
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
  SUBROUTINE matini(ir)
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: ir

    INTEGER:: idbg
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
    INTEGER:: lratio
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
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    COMMON /debug/ idbg
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    idbg=0

    ir=0
    izz=0
    izz0=0
    lenv = 1000000000
    lenv2=lratio*lenv
    neqns1=neqns+2
    len_dsln = lenv2

    !
    !rmiv
    !      iv(1)=51
    !      iv(24)=neqns
    !      iv(25)=lenv2

    len_iv = lenv2

    iv1 = 51

    !*Initialize allocation measure variables
    ialoc = 0
    raloc = 0
    !
    !  set z pivot
    !
    ALLOCATE(zpiv(neqns),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, zpiv: SUB. matini"
    CALL zpivot(neqns,neqnsz,nttbr,jcol,irow,zpiv,ir1)
    IF(ir1/=0) THEN
      ir=ir1
      RETURN
    ENDIF
    !
    !  build jcpt,jcolno
    !
    !rmiv
    lcpt=iv1+neqns1
    lcolno=lcpt+2*nttbr
    left=lcolno+2*nttbr
    last=lenv2

    !rmem
    ALLOCATE(jcpt(2*nttbr),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, jcpt: SUB. matini"
    ALLOCATE(jcolno(2*nttbr),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, jcolno: SUB. matini"
    CALL stsmat(neqns,nttbr,irow,jcol,jcpt,jcolno)

    !
    !  build ia,ja
    !
    lia=last-neqns1
    lja=lia-nttbr*2
    last=lja
    !rmem
    ALLOCATE(ia(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, ia: SUB. matini"
    !ialoc = ialoc + neqns+1
    !WINDEBUG
    ALLOCATE(ja(2*nttbr),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, ja: SUB. matini"
    CALL stiaja(neqns,ia,ja,jcpt,jcolno)

    !*Deallocation of work array
    DEALLOCATE(jcpt)
    DEALLOCATE(jcolno)
    !
    !  get permutation vector iperm,invp
    !
    lwk1=lja-neqns1
    lwk2=lwk1-neqns1
    lwk3=lwk2-neqns1
    lwk4=lwk3-neqns1
    lwk5=lwk4-neqns1
    lwk6=lwk5-neqns1
    lwk7=lwk6-neqns1
    lwk8=lwk7-2*nttbr
    last=lwk8

    left = iv1+5*neqns1

    ALLOCATE(iperm(neqns),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, iperm: SUB. matini"
    ALLOCATE(invp(neqns),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, invp: SUB. matini"
    ALLOCATE(deg(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, deg: SUB. matini"
    ALLOCATE(marker(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, marker: SUB. matini"
    !WINDEBUG
    !      ALLOCATE(rchset(neqns+1),STAT=ierror)
    ALLOCATE(rchset(0:neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, rchset: SUB. matini"
    ALLOCATE(nbrhd(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, nbrhd: SUB. matini"
    ALLOCATE(qsize(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, qsize: SUB. matini"
    ALLOCATE(qlink(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, qlink: SUB. matini"
    !WINDEBUG
    ALLOCATE(adjncy(2*nttbr),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, adjncy: SUB. matini"

    CALL genqmd(neqnsz,ia,ja,iperm,invp,deg,marker,rchset,nbrhd,qsize,qlink,nofsub,adjncy)
    !rmiv
    !
    !   build up the parent vector parent vector will be saved in
    !   work2 for a while
    !
    DO
      CALL genpaq(ia,ja,invp,iperm,marker, neqns,rchset)
      !
      !   build up the binary tree
      !
      lbtree=lwk3-2*neqns
      last=lbtree
      !rmem
      !WINDEBUG
      ALLOCATE(btree(2*(neqns+1)),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, btree: SUB. matini"

      CALL genbtq(invp,marker,btree,zpiv,izz,neqns)
      !
      !   rotate the binary tree to avoid a zero pivot
      !
      IF(izz==0) EXIT
      IF(izz0==0) izz0=izz
      IF(izz0/=izz) THEN
        CALL bringu(zpiv,iperm,invp,marker,izz,neqns,irr)
        CYCLE
      ENDIF
      lwk4=last-neqns1
      lwk5=lwk4-neqns1
      last=lwk5
      CALL rotate(ia,ja,invp,iperm,marker,btree,izz,neqns,nbrhd,qsize,irr)
    ENDDO

    !
    !   post ordering
    !
    lpordr=last-neqns1
    last=lpordr
    !rmem
    ALLOCATE(parent(neqns),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, parent: SUB. matini.f"
    ALLOCATE(nch(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, nch: SUB. matini.f"
    ALLOCATE(pordr(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, pordr: SUB. matini.f"
    CALL posord(parent,btree,invp,iperm,pordr,nch,neqns,deg,marker,rchset)
    !
    !   generate skelton graph
    !
    lleaf=last-nttbr
    lxleaf=lleaf-neqns1
    ladp=lxleaf-neqns1
    last=ladp
    !rmem
    ALLOCATE(adjncp(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, adjncp: SUB. matini.f"
    ALLOCATE(xleaf(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, xleaf: SUB. matini.f"
    ALLOCATE(leaf(nttbr),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, leaf: SUB. matini.f"
    CALL gnleaf(ia,ja,invp,iperm,nch,adjncp,xleaf,leaf,neqns,lnleaf)
    CALL forpar(neqns,parent,nch,nstop)

    !*Deallocation of work arrays
    DEALLOCATE(ia)
    DEALLOCATE(ja)
    DEALLOCATE(deg)
    DEALLOCATE(marker)
    DEALLOCATE(rchset)
    DEALLOCATE(nbrhd)
    DEALLOCATE(qsize)
    DEALLOCATE(qlink)
    DEALLOCATE(adjncy)
    DEALLOCATE(zpiv)
    !*Nullify pointers
    NULLIFY(ia)
    NULLIFY(ja)
    NULLIFY(deg)
    NULLIFY(marker)
    NULLIFY(rchset)
    NULLIFY(nbrhd)
    NULLIFY(qsize)
    NULLIFY(qlink)
    NULLIFY(adjncy)
    NULLIFY(zpiv)

    !
    !   build up xlnzr,colno  (this is the symbolic fct.)
    !

    !rmiv
    maxl=lxleaf-(left+neqns1)
    ALLOCATE(xlnzr(neqns+1),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, xlnzr: SUB. matini.f"
    CALL pre_gnclno(parent,xleaf,leaf,xlnzr,neqns,nstop,lncol,ir1)
    ALLOCATE(colno(lncol),STAT=ierror)
    IF(ierror/=0) STOP "ALLOCATION ERROR, colno: SUB. matini.f"
    CALL gnclno(parent,xleaf,leaf,xlnzr,colno,neqns,nstop,lncol,ir1)

    !*Deallocate work arrays
    DEALLOCATE(pordr)
    DEALLOCATE(adjncp)
    DEALLOCATE(xleaf)
    DEALLOCATE(leaf)
    DEALLOCATE(btree)
    !*Nullify pointers
    NULLIFY(pordr)
    NULLIFY(adjncp)
    NULLIFY(xleaf)
    NULLIFY(leaf)
    NULLIFY(btree)
    !rmem
    !rmiv
    left=(left+neqns1)+lncol
    !rmiv
    len_dsln =(neqns-nstop+1)*(neqns-nstop)/2

    !Scalar assignments
    len_colno = lncol
    !
    !   area for double precision values
    !
    IF(MOD(left,2)==0) left=left+1

    !rmiv
    total = left

    !rmiv
    stage = 10
    ialoc = 5*neqns+lncol+1

    RETURN
  END SUBROUTINE matini

  !======================================================================!
  !> nufct
  !     nufct performs cholesky factorization in row order
  !
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !
  !     (o) zln,diag,dsln
  !
  !         #coded by t.arakawa of RIST on 040329
  !======================================================================!
  SUBROUTINE nufct(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: zln

    INTEGER:: ic
    INTEGER:: isem
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: diag
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /mchdpn/ rmax,rmin,epsm,lratio
    COMMON isem

    isem=1
    !
    ! phase I
    !
    CALL ptime(t1)
    diag(1)=1.0d0/diag(1)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    DO ic=2,nstop-1
      CALL SUM(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
    ENDDO
    !
    ! phase II
    !
    CALL ptime(t2)
    DO ic=nstop,neqns
      CALL sum1(ic,xlnzr,colno,zln,temp,indx)
    ENDDO
    !
    ! phase III
    !
    CALL ptime(t3)
    CALL sum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    !
    ! phase IV
    !
    CALL ptime(t4)
    CALL sum3(neqns-nstop+1,dsln,diag(nstop),indx,temp)
    CALL ptime(t5)
    tt=t5-t1
    t1=t2-t1
    t2=t3-t2
    t3=t4-t3
    t4=t5-t4
    RETURN
    ir=30
    RETURN
  END SUBROUTINE nufct

  !======================================================================!
  !> nufct0
  !     this performs Cholesky factorization
  !
  !          if(iv(22)==0)    normal type
  !          if(iv(22)>0)    code generation type
  !======================================================================!
  SUBROUTINE nufct0(ir)
    USE hecmw_util
    IMPLICIT NONE

    INTEGER, INTENT(OUT):: ir

    INTEGER:: iseed
    INTEGER:: ixxx
    INTEGER:: lratio
    INTEGER:: ndeg2
    INTEGER:: ndegl
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    COMMON /mchdpn/ rmax,rmin,epsm,lratio
    COMMON /qaz/iseed,ixxx
    !dimension iv(*)
    ! caution) under def. may cause stack overflow.
    !dimension temp(ndeg,ndeg,neqns),indx(neqns)

    IF(stage/=20) THEN
      PRINT *,'*********Setting Stage 40!*********'
      ir=40
      RETURN
    ELSE
      ir=0
    ENDIF

    ALLOCATE(temp(ndeg*ndeg*neqns), STAT=irr)
    IF(irr /= 0) THEN
      WRITE(*,*) '##Error : Not enough memory'
      CALL hecmw_abort( hecmw_comm_get_comm())
      !stop
    ENDIF
    ALLOCATE(indx(neqns), STAT=irr)
    IF(irr /= 0) THEN
      WRITE(*,*) '##Error : Not enough memory'
      CALL hecmw_abort( hecmw_comm_get_comm())
      !stop
    ENDIF
    !rmiv
    ndegl=ndeg*(ndeg+1)
    ndegl=ndegl/2
    ndeg2=ndeg*ndeg
    !rmiv
    IF(ndeg==1) THEN
      CALL nufct(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    ELSEIF(ndeg==2) THEN
      CALL nufct2(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    ELSEIF(ndeg==3) THEN
      CALL nufct3(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    ELSEIF(ndeg==6) THEN
      CALL nufct6(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    ELSE
      CALL nufctx(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ndeg,ndegl,ir)
    ENDIF
    stage=30
    DEALLOCATE( temp )
    DEALLOCATE( indx )
    RETURN
  END SUBROUTINE nufct0

  !======================================================================!
  !> nufct2
  !     nufct performs cholesky factorization in row order
  !
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !
  !     (o) zln,diag,dsln
  !
  !         #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE nufct2(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: temp

    INTEGER:: ic
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    !
    ! phase I
    !
    ir=0
    CALL ptime(t1)
    IF(nstop>1) CALL inv2(diag(1,1),ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    DO ic=2,nstop-1
      CALL s2um(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
    ENDDO
    !
    ! phase II
    !
    CALL ptime(t2)
    DO ic=nstop,neqns
      CALL s2um1(ic,xlnzr,colno,zln,temp,indx)
    ENDDO
    !
    ! phase III
    !
    CALL ptime(t3)
    CALL s2um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    !
    ! phase IV
    !
    CALL ptime(t4)
    CALL s2um3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp)
    CALL ptime(t5)
    tt=t5-t1
    t1=t2-t1
    t2=t3-t2
    t3=t4-t3
    t4=t5-t4
    RETURN
  END SUBROUTINE nufct2

  !======================================================================!
  !> nufct3
  !     nufct performs cholesky factorization in row order
  !
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !
  !     (o) zln,diag,dsln
  !
  !         #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE nufct3(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: zln

    INTEGER:: ic
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    !
    ! phase I
    !
    CALL ptime(t1)
    IF(nstop>1) CALL inv3(diag(1,1),ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    DO ic=2,nstop-1
      CALL s3um(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
    ENDDO
    !
    ! phase II
    !
    CALL ptime(t2)
    DO ic=nstop,neqns
      CALL s3um1(ic,xlnzr,colno,zln,temp,indx)
    ENDDO
    !
    ! phase III
    !
    CALL ptime(t3)
    CALL s3um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    !
    ! phase IV
    !
    CALL ptime(t4)
    CALL s3um3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp)
    CALL ptime(t5)
    tt=t5-t1
    t1=t2-t1
    t2=t3-t2
    t3=t4-t3
    t4=t5-t4
    RETURN
  END SUBROUTINE nufct3

  !======================================================================!
  !> nufct6
  !     nufct performs cholesky factorization in row order
  !
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !
  !     (o) zln,diag,dsln
  !
  !         #coded by t.arakawa of RIST 040331
  !======================================================================!
  SUBROUTINE nufct6(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(21,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,*):: dsln

    INTEGER:: ic
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    !
    ! phase I
    !
    CALL ptime(t1)
    IF(nstop>1) CALL inv6(diag(1,1),ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    DO ic=2,nstop-1
      CALL s6um(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
    ENDDO
    !
    ! phase II
    !
    CALL ptime(t2)
    DO ic=nstop,neqns
      CALL s6um1(ic,xlnzr,colno,zln,temp,indx)
    ENDDO
    !
    ! phase III
    !
    CALL ptime(t3)
    CALL s6um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    !
    ! phase IV
    !
    CALL ptime(t4)
    CALL s6um3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp)
    CALL ptime(t5)
    tt=t5-t1
    t1=t2-t1
    t2=t3-t2
    t3=t4-t3
    t4=t5-t4
    RETURN
  END SUBROUTINE nufct6

  !======================================================================!
  !> nufctx
  !     nufct performs cholesky factorization in row order
  !
  !     (i) xlnzr,colno,zln,diag
  !         symbolicaly factorized
  !
  !     (o) zln,diag,dsln
  !
  !         #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE nufctx(xlnzr,colno,dsln,zln,diag,indx,temp,neqns,parent,nch,nstop,ndeg,ndegl,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ndeg
    INTEGER, INTENT(IN):: ndegl
    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: parent

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg*ndeg,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndegl,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg*ndeg,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg*ndeg,*):: zln

    INTEGER:: ic
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: t1
    DOUBLE PRECISION:: t2
    DOUBLE PRECISION:: t3
    DOUBLE PRECISION:: t4
    DOUBLE PRECISION:: t5
    DOUBLE PRECISION:: tt
    DOUBLE PRECISION, DIMENSION(100):: zz
    DOUBLE PRECISION, DIMENSION(100):: t
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    !
    ! phase I
    !
    CALL ptime(t1)
    IF(nstop>1) CALL invx(diag(1,1),ndeg,ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    DO ic=2,nstop-1
      CALL sxum(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx,ndegl,zz,t)
    ENDDO
    !
    ! phase II
    !
    CALL ptime(t2)
    DO ic=nstop,neqns
      CALL sxum1(ic,xlnzr,colno,zln,temp,indx,t)
    ENDDO
    !
    ! phase III
    !
    CALL ptime(t3)
    CALL sxum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx,ndegl)
    !
    ! phase IV
    !
    CALL ptime(t4)
    CALL sxum3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp,ndegl,t)
    CALL ptime(t5)
    tt=t5-t1
    t1=t2-t1
    t2=t3-t2
    t3=t4-t3
    t4=t5-t4
    RETURN
  END SUBROUTINE nufctx

  !======================================================================!
  !> nuform
  !======================================================================!
  SUBROUTINE nuform(hecMESH,hecMAT,ir)
    USE hecmw_util
    IMPLICIT NONE

    TYPE (hecmwST_local_mesh), INTENT(IN):: hecMESH
    TYPE (hecmwST_matrix    ), INTENT(IN):: hecMAT

    INTEGER(KIND=kint), INTENT(OUT):: ir

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

    NUMNP = hecMAT%NP
    NDOF  = hecMESH%n_dof
    ntotal = NUMNP*NDOF
    idbg=0


    !*NUFACT variables
    neqns = NUMNP
    ndeg = NDOF
    nttbr = hecMAT%NP+hecMAT%NPL !+hecMAT%NPU if unsymmetric
    isym = 0

    !*Allocations
    ALLOCATE(val(ndeg*ndeg),STAT=ierr)
    IF(ierr/=0) STOP "Allocation error:val"
    WRITE(6,*) "nuform:stage = ",stage
    kk = 0
    ndof2 = NDOF*NDOF
    DO j= 1, NUMNP
      !*Diagonal
      kk = kk + 1
      CALL vlcpy(val,hecMAT%D(ndof2*(j-1)+1:ndof2*j),ndof)
      CALL staij1(0,j,j,val,ir)

      DO i = 1,NDOF
        IF(val((i-1)*ndof+i)<=0) THEN
          WRITE(IDBG,*)'j,j,val:',j,i,val((i-1)*ndof+i)
          !          PAUSE 'Error?'
        ENDIF
      END DO

      !*Lower
      DO k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i= hecMAT%itemL(k)
        kk = kk + 1
        CALL vlcpy(val,hecMAT%AL(ndof2*(k-1)+1:ndof2*k),ndof)

        CALL staij1(0,j,i,val,ir)
      ENDDO
    ENDDO

    DEALLOCATE(val)
    RETURN
  END SUBROUTINE nuform

  !======================================================================!
  !> nusol0
  !     this performs forward elimination and backward substitution
  !
  !     (i/o)
  !           r_h_s        on entry     right hand side vector
  !                    on exit      solution vector
  !           iv       communication array
  !
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE nusol0(r_h_s,ir)
    USE hecmw_util
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: r_h_s
    INTEGER, INTENT(OUT):: ir

    INTEGER:: iseed
    INTEGER:: ixxx
    INTEGER:: lratio
    INTEGER:: lwk
    INTEGER:: ndegl
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(:), POINTER :: wk
    COMMON /mchdpn/ rmax,rmin,epsm,lratio
    COMMON /qaz/iseed,ixxx

    IF(stage/=30.AND.stage/=40) THEN
      ir=50
      RETURN
    ELSE
      ir=0
    ENDIF
    lwk=total

    ALLOCATE( wk(ndeg*neqns), STAT=ierror )
    IF( ierror /= 0 ) THEN
      WRITE(*,*) "##Error: not enough memory"
      CALL hecmw_abort( hecmw_comm_get_comm() )
    ENDIF
    !rmiv
    ndegl=ndeg*(ndeg+1)
    ndegl=ndegl/2
    IF(ndeg==1) THEN
      CALL nusol1(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop)
    ELSEIF(ndeg==2) THEN
      CALL nusol2(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop)
    ELSEIF(ndeg==3) THEN
      CALL nusol3(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop)
    ELSEIF(ndeg==6) THEN
      CALL nusolx(xlnzr,colno,dsln,zln,diag, iperm,r_h_s,wk,neqns,nstop,ndegl)
    ELSE
      CALL nusolx(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop,ndegl)
    ENDIF
    stage=40
    DEALLOCATE( wk )
    RETURN
  END SUBROUTINE nusol0

  !======================================================================!
  !> nusol1
  !======================================================================!
  SUBROUTINE nusol1(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns,nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: iperm
    DOUBLE PRECISION, INTENT(IN), DIMENSION(*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(*):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(*):: dsln

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: wk

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    ! forward
    DO i=1,neqns
      wk(i)=b(iperm(i))
    ENDDO
    joc=1
    DO i=1,neqns
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) GOTO 110
      wk(i)=wk(i)-spdot2(wk,ks,ke)

110   CONTINUE
      IF(i<=nstop) CYCLE
      wk(i)=wk(i)-ddot(wk(nstop),dsln(joc),i-nstop)
      joc=joc+i-nstop
    ENDDO
    DO i=1,neqns
      wk(i)=wk(i)*diag(i)
    ENDDO
    ! back ward
    DO i=neqns,1,-1
      IF(i<nstop) GOTO 206
      DO j=i-1,nstop,-1
        joc=joc-1
        wk(j)=wk(j)-wk(i)*dsln(joc)
      ENDDO

206   CONTINUE
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) CYCLE
      DO k=ks,ke
        j=colno(k)
        wk(j)=wk(j)-wk(i)*zln(k)
      ENDDO
    ENDDO
    ! permutaion
    DO i=1,neqns
      b(iperm(i))=wk(i)
    ENDDO
    RETURN
  END SUBROUTINE nusol1

  !======================================================================!
  !> nusol2
  !======================================================================!
  SUBROUTINE nusol2(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns,nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: iperm
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,*):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,*):: dsln

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2,*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2,*):: wk

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    ! forward
    DO i=1,neqns
      wk(1,i)=b(1,iperm(i))
      wk(2,i)=b(2,iperm(i))
    ENDDO
    joc=1
    DO i=1,neqns
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) GOTO 110
      CALL s2pdot(wk(1,i),wk,zln,colno,ks,ke)
110   CONTINUE
      IF(i<=nstop) CYCLE
      CALL d2sdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
      joc=joc+i-nstop
    ENDDO
    DO i=1,neqns
      wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
      wk(1,i)=wk(1,i)*diag(1,i)
      wk(2,i)=wk(2,i)*diag(3,i)
      wk(1,i)=wk(1,i)-wk(2,i)*diag(2,i)
    ENDDO
    ! back ward
    DO i=neqns,1,-1
      IF(i<nstop) GOTO 206
      DO j=i-1,nstop,-1
        joc=joc-1
        wk(1,j)=wk(1,j)-wk(1,i)*dsln(1,joc)-wk(2,i)*dsln(2,joc)
        wk(2,j)=wk(2,j)-wk(1,i)*dsln(3,joc)-wk(2,i)*dsln(4,joc)
      ENDDO

206   CONTINUE
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) CYCLE
      DO k=ks,ke
        j=colno(k)
        wk(1,j)=wk(1,j)-wk(1,i)*zln(1,k)-wk(2,i)*zln(2,k)
        wk(2,j)=wk(2,j)-wk(1,i)*zln(3,k)-wk(2,i)*zln(4,k)
      ENDDO
    ENDDO
    ! permutaion
    DO i=1,neqns
      b(1,iperm(i))=wk(1,i)
      b(2,iperm(i))=wk(2,i)
    ENDDO
    RETURN
  END SUBROUTINE nusol2

  !======================================================================!
  !> nusol3
  !======================================================================!
  SUBROUTINE nusol3(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns,nstop)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: iperm
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: dsln

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,*):: wk

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks

    ! forward
    DO i=1,neqns
      wk(1,i)=b(1,iperm(i))
      wk(2,i)=b(2,iperm(i))
      wk(3,i)=b(3,iperm(i))
    ENDDO
    joc=1
    DO i=1,neqns
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) GOTO 110
      CALL s3pdot(wk(1,i),wk,zln,colno,ks,ke)

110   CONTINUE
      IF(i<=nstop) CYCLE
      CALL d3sdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
      joc=joc+i-nstop
    ENDDO
    DO i=1,neqns
      wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
      wk(3,i)=wk(3,i)-wk(1,i)*diag(4,i)-wk(2,i)*diag(5,i)
      wk(1,i)=wk(1,i)*diag(1,i)
      wk(2,i)=wk(2,i)*diag(3,i)
      wk(3,i)=wk(3,i)*diag(6,i)
      wk(2,i)=wk(2,i)-wk(3,i)*diag(5,i)
      wk(1,i)=wk(1,i)-wk(2,i)*diag(2,i)-wk(3,i)*diag(4,i)
    ENDDO
    ! back ward
    DO i=neqns,1,-1
      IF(i<nstop) GOTO 206
      DO j=i-1,nstop,-1
        joc=joc-1
        wk(1,j)=wk(1,j)-wk(1,i)*dsln(1,joc)-wk(2,i)*dsln(2,joc)-wk(3,i)*dsln(3,joc)
        wk(2,j)=wk(2,j)-wk(1,i)*dsln(4,joc)-wk(2,i)*dsln(5,joc)-wk(3,i)*dsln(6,joc)
        wk(3,j)=wk(3,j)-wk(1,i)*dsln(7,joc)-wk(2,i)*dsln(8,joc)-wk(3,i)*dsln(9,joc)
      ENDDO

206   CONTINUE
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) CYCLE
      DO k=ks,ke
        j=colno(k)
        wk(1,j)=wk(1,j)-wk(1,i)*zln(1,k)-wk(2,i)*zln(2,k)-wk(3,i)*zln(3,k)
        wk(2,j)=wk(2,j)-wk(1,i)*zln(4,k)-wk(2,i)*zln(5,k)-wk(3,i)*zln(6,k)
        wk(3,j)=wk(3,j)-wk(1,i)*zln(7,k)-wk(2,i)*zln(8,k)-wk(3,i)*zln(9,k)
      ENDDO
    ENDDO
    ! permutaion
    DO i=1,neqns
      b(1,iperm(i))=wk(1,i)
      b(2,iperm(i))=wk(2,i)
      b(3,iperm(i))=wk(3,i)
    ENDDO
    RETURN
  END SUBROUTINE nusol3

  !======================================================================!
  !> nusolx
  !======================================================================!
  SUBROUTINE nusolx(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns,nstop,ndegl)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: ndegl
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: iperm
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndegl,*):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: dsln

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,*):: wk

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: loc1
    INTEGER:: locd
    INTEGER:: m
    INTEGER:: n

    ! forward
    DO l=1,ndeg
      DO i=1,neqns
        wk(l,i)=b(l,iperm(i))
      ENDDO
    ENDDO
    joc=1
    DO i=1,neqns
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) GOTO 110
      CALL sxpdot(wk(1,i),wk,zln,colno,ks,ke)
110   CONTINUE
      IF(i<=nstop) EXIT
      CALL dxsdot(wk(1,i),wk(1,nstop),dsln(1,1,joc),i-nstop)
      joc=joc+i-nstop
    ENDDO
    DO i=1,neqns
      locd=0
      DO m=1,ndeg-1
        locd=locd+m
        loc1=locd+m
        DO n=m+1,ndeg
          wk(n,i)=wk(n,i)-wk(m,i)*diag(loc1,i)
          loc1=loc1+n
        ENDDO
      ENDDO
      locd=0
      DO m=1,ndeg
        locd=locd+m
        wk(m,i)=wk(m,i)*diag(locd,i)
      ENDDO
      DO n=ndeg,2,-1
        locd=locd-1
        DO m=n-1,1,-1
          wk(m,i)=wk(m,i)-wk(n,i)*diag(locd,i)
          locd=locd-1
        ENDDO
      ENDDO
    ENDDO
    ! back ward
    DO i=neqns,1,-1
      IF(i<nstop) GOTO 208
      DO j=i-1,nstop,-1
        joc=joc-1
        DO m=1,ndeg
          DO n=1,ndeg
            wk(m,j)=wk(m,j)-wk(n,i)*dsln(n,m,joc)
          ENDDO
        ENDDO
      ENDDO

208   CONTINUE
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      IF(ke<ks) CYCLE
      DO k=ks,ke
        j=colno(k)
        DO m=1,ndeg
          DO n=1,ndeg
            wk(m,j)=wk(m,j)-wk(n,i)*zln(n,m,k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! permutaion
    DO l=1,ndeg
      DO i=1,neqns
        b(l,iperm(i))=wk(l,i)
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE nusolx

  !======================================================================!
  !> posord
  !======================================================================!
  SUBROUTINE posord(parent,btree,invp,iperm,pordr,nch,neqns,iw,qarent,mch)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(2,*):: btree
    INTEGER, INTENT(IN), DIMENSION(*):: qarent

    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    INTEGER, INTENT(OUT), DIMENSION(*):: parent
    INTEGER, INTENT(OUT), DIMENSION(*):: pordr
    INTEGER, INTENT(OUT), DIMENSION(*):: invp
    INTEGER, INTENT(OUT), DIMENSION(*):: iperm
    INTEGER, INTENT(OUT), DIMENSION(*):: iw
    INTEGER, INTENT(OUT), DIMENSION(0:neqns+1):: mch

    INTEGER:: i
    INTEGER:: idbg1
    INTEGER:: ii
    INTEGER:: invpos
    INTEGER:: ipinv
    INTEGER:: joc
    INTEGER:: l
    INTEGER:: locc
    INTEGER:: locp
    COMMON /debug/ idbg1
    !
    DO i=1,neqns
      mch(i)=0
      pordr(i)=0
    ENDDO
    l=1
    locc=neqns+1

10  CONTINUE
    DO
      joc=locc
      locc=btree(1,joc)
      IF(locc==0) EXIT
    ENDDO

    locp=qarent(joc)
    mch(locp)=mch(locp)+1

    DO
      pordr(joc)=l
      IF(l>=neqns) GOTO 1000
      l=l+1
      locc=btree(2,joc)
      IF(locc/=0) GOTO 10
      joc=qarent(joc)
      locp=qarent(joc)
      mch(locp)=mch(locp)+mch(joc)+1
    ENDDO

1000 CONTINUE
    DO i=1,neqns
      ipinv=pordr(invp(i))
      invp(i)=ipinv
      iperm(ipinv)=i
      iw(pordr(i))=i
    ENDDO
    DO i=1,neqns
      invpos=iw(i)
      nch(i)=mch(invpos)
      ii=qarent(invpos)
      IF(ii>0.AND.ii<=neqns) THEN
        parent(i)=pordr(ii)
      ELSE
        parent(i)=qarent(invpos)
      ENDIF
    ENDDO
    IF(idbg1/=0) WRITE(6,"(' post order')")
    IF(idbg1/=0) WRITE(6,"(10i6)") (pordr(i),i=1,neqns)
    IF(idbg1/=0) WRITE(6,"(/' invp ')")
    IF(idbg1/=0) WRITE(6,"(/' parent')")
    IF(idbg1/=0) WRITE(6,"(10i6)") (parent(i),i=1,neqns)
    IF(idbg1/=0) WRITE(6,"(10i6)") (invp(i),i=1,neqns)
    IF(idbg1/=0) WRITE(6,"(/' iperm ')")
    IF(idbg1/=0) WRITE(6,"(10i6)") (iperm(i),i=1,neqns)
    IF(idbg1/=0) WRITE(6,"(' nch')")
    IF(idbg1/=0) WRITE(6,"(10i6)") (nch(i),i=1,neqns)
    RETURN
  END SUBROUTINE posord

  !======================================================================!
  !> pre_gnclno
  !======================================================================!
  SUBROUTINE pre_gnclno(parent,xleaf,leaf,xlnzr,neqns,nstop,lncol,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: parent
    INTEGER, INTENT(IN), DIMENSION(*):: xleaf
    INTEGER, INTENT(IN), DIMENSION(*):: leaf
    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT):: lncol
    INTEGER, INTENT(OUT), DIMENSION(*):: xlnzr

    INTEGER:: i
    INTEGER:: j
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: nc
    INTEGER:: nxleaf
    INTEGER:: idbg1
    COMMON /debug/ idbg1

    nc=0
    ir=0
    l=1
    DO i=1,neqns
      xlnzr(i)=l
      ks=xleaf(i)
      ke=xleaf(i+1)-1
      IF(ke<ks) CYCLE
      nxleaf=leaf(ks)
      DO k=ks,ke-1
        j=nxleaf
        nxleaf=leaf(k+1)
        DO
          IF(j>=nxleaf) EXIT
          IF(j>=nstop) GOTO 100
          l=l+1
          j=parent(j)
        ENDDO
      ENDDO

      j=leaf(ke)
      DO
        IF(j>=i.OR.j==0.OR.j>=nstop) EXIT
        l=l+1
        j=parent(j)
      ENDDO
100   CONTINUE
    ENDDO
    xlnzr(neqns+1)=l
    lncol=l-1
    RETURN
  END SUBROUTINE pre_gnclno

  !======================================================================!
  !> ptime
  !======================================================================!
  SUBROUTINE ptime(cputim)
    USE hecmw_util
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(OUT):: cputim

    cputim=hecmw_Wtime()
    RETURN
  END SUBROUTINE ptime
  !======================================================================!
  !> qmdmrg
  !======================================================================!
  SUBROUTINE qmdmrg(xadj,adjncy,deg,qsize,qlink,marker,deg0,nhdsze,nbrhd,rchset,ovrlp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: deg0
    INTEGER, INTENT(IN):: nhdsze
    INTEGER, INTENT(IN), DIMENSION(*):: adjncy
    INTEGER, INTENT(IN), DIMENSION(*):: nbrhd
    INTEGER, INTENT(IN), DIMENSION(*):: xadj

    INTEGER, INTENT(OUT), DIMENSION(*):: deg
    INTEGER, INTENT(OUT), DIMENSION(*):: qsize
    INTEGER, INTENT(OUT), DIMENSION(*):: qlink
    INTEGER, INTENT(OUT), DIMENSION(*):: marker
    INTEGER, INTENT(OUT), DIMENSION(*):: rchset
    INTEGER, INTENT(OUT), DIMENSION(*):: ovrlp

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

    IF(nhdsze<=0) RETURN
    DO inhd=1,nhdsze
      root=nbrhd(inhd)
      marker(root)=0
    ENDDO

    DO inhd=1,nhdsze
      root=nbrhd(inhd)
      marker(root)=-1
      rchsze=0
      novrlp=0
      deg1=0

200   jstrt=xadj(root)
      jstop=xadj(root+1)-1

      DO j=jstrt,jstop
        nabor=adjncy(j)
        root=-nabor
        IF(nabor<0) GOTO 200
        IF(nabor==0) EXIT
        mark=marker(nabor)

        IF(mark<0) CYCLE
        IF(mark>0) THEN
          !GOTO 500
          IF(mark>1) CYCLE
          novrlp=novrlp+1
          ovrlp(novrlp)=nabor
          marker(nabor)=2
          CYCLE
        ENDIF
        rchsze=rchsze+1
        rchset(rchsze)=nabor
        deg1=deg1+qsize(nabor)
        marker(nabor)=1
      ENDDO

      head=0
      mrgsze=0

      DO iov=1,novrlp
        node=ovrlp(iov)
        jstrt=xadj(node)
        jstop=xadj(node+1)-1
        DO j=jstrt,jstop
          nabor=adjncy(j)
          IF(marker(nabor)/=0) CYCLE
          marker(node)=1
          GOTO 1100
        ENDDO

        mrgsze=mrgsze+qsize(node)
        marker(node)=-1
        lnode=node

        DO
          link=qlink(lnode)
          link=qlink(lnode)
          IF(link<=0) EXIT
          lnode=link
        ENDDO
        qlink(lnode)=head
        head=node
1100    CONTINUE
      ENDDO

      IF(head<=0) GOTO 1200
      qsize(head)=mrgsze
      deg(head)=deg0+deg1-1
      marker(head)=2

1200  CONTINUE
      root=nbrhd(inhd)
      marker(root)=0
      IF(rchsze<=0) CYCLE
      DO irch=1,rchsze
        node=rchset(irch)
        marker(node)=0
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE qmdmrg
  !======================================================================!
  !> qmdupd
  !======================================================================!
  SUBROUTINE qmdupd(xadj,adjncy,nlist,list,deg,qsize,qlink,marker,rchset,nbrhd)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: nlist
    INTEGER, INTENT(IN), DIMENSION(*):: adjncy
    INTEGER, INTENT(IN), DIMENSION(*):: list
    INTEGER, INTENT(IN), DIMENSION(*):: xadj

    INTEGER, INTENT(OUT), DIMENSION(*):: deg
    INTEGER, INTENT(OUT), DIMENSION(*):: marker
    INTEGER, INTENT(OUT), DIMENSION(*):: rchset
    INTEGER, INTENT(OUT), DIMENSION(*):: nbrhd
    INTEGER, INTENT(OUT), DIMENSION(*):: qsize
    INTEGER, INTENT(OUT), DIMENSION(*):: qlink

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

    IF(nlist<=0) RETURN
    deg0=0
    nhdsze=0
    DO il=1,nlist
      node=list(il)
      deg0=deg0+qsize(node)
      jstrt=xadj(node)
      jstop=xadj(node+1)-1
      DO j=jstrt,jstop
        nabor=adjncy(j)
        IF(marker(nabor)/=0.OR.deg(nabor)>=0) CYCLE
        marker(nabor)=-1
        nhdsze=nhdsze+1
        nbrhd(nhdsze)=nabor
      ENDDO
    ENDDO

    IF(nhdsze>0) CALL qmdmrg(xadj,adjncy,deg,qsize,qlink,marker,deg0,nhdsze,nbrhd,rchset,nbrhd(nhdsze+1))
    DO il=1,nlist
      node=list(il)
      mark=marker(node)
      IF(mark>1.OR.mark<0) CYCLE
      CALL qmdrch(node,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze,nbrhd)
      deg1=deg0
      IF(rchsze<=0) GOTO 400

      DO irch=1,rchsze
        inode=rchset(irch)
        deg1=deg1+qsize(inode)
        marker(inode)=0
      ENDDO

400   CONTINUE
      deg(node)=deg1-1
      IF(nhdsze<=0) CYCLE
      DO inhd=1,nhdsze
        inode=nbrhd(inhd)
        marker(inode)=0
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE qmdupd

  !======================================================================!
  !> qmdot
  !======================================================================!
  SUBROUTINE qmdot(root,xadj,adjncy,marker,rchsze,rchset,nbrhd)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: rchsze
    INTEGER, INTENT(IN):: root
    INTEGER, INTENT(IN), DIMENSION(*):: marker
    INTEGER, INTENT(IN), DIMENSION(*):: rchset
    INTEGER, INTENT(IN), DIMENSION(*):: nbrhd
    INTEGER, INTENT(IN), DIMENSION(*):: xadj

    INTEGER, INTENT(OUT), DIMENSION(*):: adjncy

    INTEGER:: inhd
    INTEGER:: irch
    INTEGER:: j
    INTEGER:: jstrt
    INTEGER:: jstop
    INTEGER:: link
    INTEGER:: nabor
    INTEGER:: node

    irch=0
    inhd=0
    node=root

    DO
      jstrt=xadj(node)
      jstop=xadj(node+1)-2

      IF(jstop<jstrt) GOTO 300

      DO j=jstrt,jstop
        irch=irch+1
        adjncy(j)=rchset(irch)
        IF(irch>=rchsze) GOTO 400
      ENDDO

300   CONTINUE
      link=adjncy(jstop+1)
      node=-link
      IF(link<0) CYCLE

      inhd=inhd+1
      node=nbrhd(inhd)
      adjncy(jstop+1)=-node
    ENDDO

400 CONTINUE
    adjncy(j+1)=0
    DO irch=1,rchsze
      node=rchset(irch)
      IF(marker(node)<0) CYCLE
      jstrt=xadj(node)
      jstop=xadj(node+1)-1
      DO j=jstrt,jstop
        nabor=adjncy(j)
        IF(marker(nabor)>=0) CYCLE
        adjncy(j)=root
        EXIT
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE qmdot

  !======================================================================!
  !> qmdrch
  !======================================================================!
  SUBROUTINE qmdrch(root,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze,nbrhd)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: root
    INTEGER, INTENT(IN), DIMENSION(*):: adjncy
    INTEGER, INTENT(IN), DIMENSION(*):: deg
    INTEGER, INTENT(IN), DIMENSION(*):: xadj

    INTEGER, INTENT(OUT), DIMENSION(*):: rchset
    INTEGER, INTENT(OUT):: nhdsze
    INTEGER, INTENT(OUT):: rchsze
    INTEGER, INTENT(OUT), DIMENSION(*):: marker
    INTEGER, INTENT(OUT), DIMENSION(*):: nbrhd

    INTEGER:: i
    INTEGER:: istrt
    INTEGER:: istop
    INTEGER:: j
    INTEGER:: jstrt
    INTEGER:: jstop
    INTEGER:: nabor
    INTEGER:: node

    nhdsze=0
    rchsze=0
    istrt=xadj(root)
    istop=xadj(root+1)-1
    IF(istop<istrt) RETURN
    DO i=istrt,istop
      nabor=adjncy(i)
      IF(nabor==0) RETURN
      IF(marker(nabor)/=0) CYCLE
      IF(deg(nabor)<0) GOTO 200
      rchsze=rchsze+1
      rchset(rchsze)=nabor
      marker(nabor)=1
      CYCLE

200   marker(nabor)=-1
      nhdsze=nhdsze+1
      nbrhd(nhdsze)=nabor

300   jstrt=xadj(nabor)
      jstop=xadj(nabor+1)-1
      DO j=jstrt,jstop
        node=adjncy(j)
        nabor=-node
        IF(node<0) GOTO 300
        IF(node==0) EXIT
        IF(marker(node)/=0) CYCLE
        rchsze=rchsze+1
        rchset(rchsze)=node
        marker(node)=1
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE qmdrch


  !======================================================================!
  !> qqsort
  !     sort in increasing order up to i
  !
  !     iw   array
  !     ik   number of input/output
  !     i    deal with numbers less than this numberi
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE qqsort(iw,ik)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ik

    INTEGER, INTENT(OUT), DIMENSION(*):: iw

    INTEGER:: idbg1
    INTEGER:: itemp
    INTEGER:: l
    INTEGER:: m
    COMMON /debug/ idbg1

    IF(ik<=1) RETURN
    DO l=1,ik-1
      DO m=l+1,ik
        IF(iw(l)<iw(m)) CYCLE
        itemp=iw(l)
        iw(l)=iw(m)
        iw(m)=itemp
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE qqsort

  !======================================================================!
  !> rotate
  !     irr return code irr=0 node izz is not a bottom node
  !                     irr=1          is a bottom node then rotation is
  !                                    performed
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE rotate(xadj,adjncy,invp,iperm,parent,btree,izz,neqns,anc,adjt,irr)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: izz
    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: xadj
    INTEGER, INTENT(IN), DIMENSION(*):: adjncy
    INTEGER, INTENT(IN), DIMENSION(*):: parent
    INTEGER, INTENT(IN), DIMENSION(2,*):: btree

    INTEGER, INTENT(OUT):: irr
    INTEGER, INTENT(OUT), DIMENSION(*):: invp
    INTEGER, INTENT(OUT), DIMENSION(*):: iperm
    INTEGER, INTENT(OUT), DIMENSION(*):: anc
    INTEGER, INTENT(OUT), DIMENSION(*):: adjt

    INTEGER:: i
    INTEGER:: idbg1
    INTEGER:: iy
    INTEGER:: izzz
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: kk
    INTEGER:: l
    INTEGER:: ll
    INTEGER:: locc
    INTEGER:: nanc
    COMMON /debug/ idbg1

    IF(izz==0) THEN
      irr=0
      RETURN
    ENDIF
    izzz=invp(izz)
    IF(btree(1,izzz)/=0) THEN
      irr=0
    ENDIF
    irr=1
    !
    !  ancestors of izzz
    !
    nanc=0
    joc=izzz

    DO
      nanc=nanc+1
      anc(nanc)=joc
      joc=parent(joc)
      IF(joc==0) EXIT
    ENDDO
    !
    !  to find the eligible node from ancestors of izz
    !
    !     adjt = Adj(Tree(y))
    l=1
200 CONTINUE
    DO i=1,neqns
      adjt(i)=0
    ENDDO
    locc=anc(l)

220 CONTINUE
    DO
      joc=locc
      locc=btree(1,joc)
      IF(locc==0) EXIT
    ENDDO

    DO
      DO k=xadj(iperm(joc)),xadj(iperm(joc)+1)-1
        adjt(invp(adjncy(k)))=1
      ENDDO

      IF(joc>=anc(l)) GOTO 250
      locc=btree(2,joc)
      IF(locc/=0) GOTO 220
      joc=parent(joc)
    ENDDO

250 CONTINUE
    DO ll=l+1,nanc
      IF(adjt(anc(ll))==0) THEN
        l=l+1
        GOTO 200
      ENDIF
    ENDDO
    IF(l==1) GOTO 500

    !
    !  anc(l-1) is the eligible node

    ! (1) number the node not in Ancestor(iy)
    iy=anc(l-1)
    DO i=1,neqns
      adjt(i)=0
    ENDDO
    DO ll=l,nanc
      adjt(anc(ll))=1
    ENDDO
    k=0
    DO ll=1,neqns
      IF(adjt(ll)==0) THEN
        k=k+1
        invp(iperm(ll))=k
      ENDIF
    ENDDO

    ! (2) followed by nodes in Ancestor(iy)-Adj(T(iy))
    adjt(1:neqns)=0
    locc=iy

350 CONTINUE
    DO
      joc=locc
      locc=btree(1,joc)
      IF(locc==0) EXIT
    ENDDO

    DO
      DO kk=xadj(iperm(joc)),xadj(iperm(joc)+1)-1
        adjt(invp(adjncy(kk)))=1
      ENDDO
      IF(joc>=iy) GOTO 380
      locc=btree(2,joc)
      IF(locc/=0) GOTO 350
      joc=parent(joc)
    ENDDO

380 CONTINUE
    DO ll=l,nanc
      IF(adjt(anc(ll))==0) THEN
        k=k+1
        invp(iperm(anc(ll)))=k
      ENDIF
    ENDDO

    ! (3) and finally number the node in Adj(t(iy))
    DO ll=l,nanc
      IF(adjt(anc(ll))/=0) THEN
        k=k+1
        invp(iperm(anc(ll)))=k
      ENDIF
    ENDDO
    GOTO 600

    ! izz can be numbered last
500 CONTINUE
    k=0
    DO i=1,neqns
      IF(i==izzz) GOTO 510
      k=k+1
      invp(iperm(i))=k
510   CONTINUE
    ENDDO
    invp(iperm(izzz))=neqns

    ! set iperm
600 CONTINUE
    DO i=1,neqns
      iperm(invp(i))=i
    ENDDO
    IF(idbg1/=0) WRITE(6,"(10i6)") (invp(i),i=1,neqns)
    RETURN
  END SUBROUTINE rotate

  !======================================================================!
  !> s2pdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE s2pdot(bi,b,zln,colno,ks,ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ke
    INTEGER, INTENT(IN):: ks
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(2,*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2):: bi

    INTEGER:: j
    INTEGER:: jj

    DO jj=ks,ke
      j=colno(jj)
      bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(3,jj)*b(2,j)
      bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(4,jj)*b(2,j)
    ENDDO
    RETURN
  END SUBROUTINE s2pdot

  !======================================================================!
  !> s2um
  !======================================================================!
  SUBROUTINE s2um(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: par

    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,*):: diag

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: lratio
    INTEGER, DIMENSION(*):: indx
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(4):: s
    DOUBLE PRECISION, DIMENSION(3):: t
    DOUBLE PRECISION, DIMENSION(4,*):: temp
    DOUBLE PRECISION, DIMENSION(4):: zz
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t(1)=0.0d0
    t(2)=0.0d0
    t(3)=0.0d0
    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      DO l=1,4
        s(l)=0.0d0
        zz(l)=zln(l,k)
      ENDDO
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(3,j)*zln(3,jj)
          zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(4,j)*zln(3,jj)
          zz(3)=zz(3)-temp(1,j)*zln(2,jj)-temp(3,j)*zln(4,jj)
          zz(4)=zz(4)-temp(2,j)*zln(2,jj)-temp(4,j)*zln(4,jj)
        ENDIF
      ENDDO
      CALL inv22(zln(1,k),zz,diag(1,jc))
      DO l=1,4
        temp(l,jc)=zz(l)
      ENDDO
      t(1)=t(1)+zz(1)*zln(1,k)+zz(3)*zln(3,k)
      t(2)=t(2)+zz(1)*zln(2,k)+zz(3)*zln(4,k)
      t(3)=t(3)+zz(2)*zln(2,k)+zz(4)*zln(4,k)
    ENDDO
    diag(1,ic)=diag(1,ic)-t(1)
    diag(2,ic)=diag(2,ic)-t(2)
    diag(3,ic)=diag(3,ic)-t(3)
    CALL inv2(diag(1,ic),ir)
    nch(ic)=-1
    kk=par(ic)
    nch(kk)=nch(kk)-1
    RETURN
  END SUBROUTINE s2um

  !======================================================================!
  !> s2um1
  !======================================================================!
  SUBROUTINE s2um1(ic,xlnzr,colno,zln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: zln

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(4):: s
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    DO l=1,4
      s(l)=0.0d0
    ENDDO

    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(3,j)*zln(3,jj)
          s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(4,j)*zln(3,jj)
          s(3)=s(3)+temp(1,j)*zln(2,jj)+temp(3,j)*zln(4,jj)
          s(4)=s(4)+temp(2,j)*zln(2,jj)+temp(4,j)*zln(4,jj)
        ENDIF
      ENDDO
      DO l=1,4
        temp(l,jc)=zln(l,k)-s(l)
        zln(l,k)=temp(l,jc)
        s(l)=0.0d0
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE s2um1

  !======================================================================!
  !> s2um2
  !======================================================================!
  SUBROUTINE s2um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: zln

    INTEGER:: ic
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno

    joc=0
    DO ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      DO k=ks,ke
        jj=colno(k)
        temp(1,jj)=zln(1,k)
        temp(2,jj)=zln(2,k)
        temp(3,jj)=zln(3,k)
        temp(4,jj)=zln(4,k)

        zln(3,k)=temp(3,jj)-temp(1,jj)*diag(2,jj)
        zln(1,k)=temp(1,jj)*diag(1,jj)
        zln(3,k)=zln(3,k)*diag(3,jj)
        zln(1,k)=zln(1,k)-zln(3,k)*diag(2,jj)

        zln(4,k)=temp(4,jj)-temp(2,jj)*diag(2,jj)
        zln(2,k)=temp(2,jj)*diag(1,jj)
        zln(4,k)=zln(4,k)*diag(3,jj)
        zln(2,k)=zln(2,k)-zln(4,k)*diag(2,jj)

        diag(1,ic)=diag(1,ic)-(temp(1,jj)*zln(1,k)+temp(3,jj)*zln(3,k))
        diag(2,ic)=diag(2,ic)-(temp(1,jj)*zln(2,k)+temp(3,jj)*zln(4,k))
        diag(3,ic)=diag(3,ic)-(temp(2,jj)*zln(2,k)+temp(4,jj)*zln(4,k))
        indx(jj)=ic
      ENDDO
      DO jc=nstop,ic-1
        joc=joc+1
        DO jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          IF(indx(j)==ic) THEN
            dsln(1,joc)=dsln(1,joc)-(temp(1,j)*zln(1,jj)+temp(3,j)*zln(3,jj))
            dsln(2,joc)=dsln(2,joc)-(temp(2,j)*zln(1,jj)+temp(4,j)*zln(3,jj))
            dsln(3,joc)=dsln(3,joc)-(temp(1,j)*zln(2,jj)+temp(3,j)*zln(4,jj))
            dsln(4,joc)=dsln(4,joc)-(temp(2,j)*zln(2,jj)+temp(4,j)*zln(4,jj))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE s2um2

  !======================================================================!
  !> s2um3
  !======================================================================!
  SUBROUTINE s2um3(n,dsln,diag,indx,temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,*):: dsln

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    DOUBLE PRECISION, DIMENSION(4):: t

    IF(n<=0) RETURN

    indx(1)=0
    joc=1
    CALL inv2(diag(1,1),ir)
    DO i=2,n
      indx(i)=joc
      DO j=1,i-1
        CALL d2dot(t,dsln(1,indx(i)),dsln(1,indx(j)),j-1)
        dsln(1,joc)=dsln(1,joc)-t(1)
        dsln(2,joc)=dsln(2,joc)-t(2)
        dsln(3,joc)=dsln(3,joc)-t(3)
        dsln(4,joc)=dsln(4,joc)-t(4)
        joc=joc+1
      ENDDO
      CALL v2prod(dsln(1,indx(i)),diag,temp,i-1)
      CALL d2dot(t,temp,dsln(1,indx(i)),i-1)
      diag(1,i)=diag(1,i)-t(1)
      diag(2,i)=diag(2,i)-t(2)
      diag(3,i)=diag(3,i)-t(4)
      CALL vcopy(temp,dsln(1,indx(i)),4*(i-1))
      CALL inv2(diag(1,i),ir)
    ENDDO
    RETURN
  END SUBROUTINE s2um3

  !======================================================================!
  !> s3pdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE s3pdot(bi,b,zln,colno,ks,ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ke
    INTEGER, INTENT(IN):: ks
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,*):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3):: bi

    INTEGER:: j
    INTEGER:: jj

    DO jj=ks,ke
      j=colno(jj)
      bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(4,jj)*b(2,j)-zln(7,jj)*b(3,j)
      bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(5,jj)*b(2,j)-zln(8,jj)*b(3,j)
      bi(3)=bi(3)-zln(3,jj)*b(1,j)-zln(6,jj)*b(2,j)-zln(9,jj)*b(3,j)
    ENDDO
    RETURN
  END SUBROUTINE s3pdot

  !======================================================================!
  !> s3um
  !======================================================================!
  SUBROUTINE s3um(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: par

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: zln

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(6):: t
    DOUBLE PRECISION, DIMENSION(9):: zz
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    !$dir max_trips(6)
    DO l=1,6
      t(l)=0.0d0
    ENDDO

    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      !$dir max_trips(9)
      DO l=1,9
        zz(l)=zln(l,k)
      ENDDO
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(4,j)*zln(4,jj)-temp(7,j)*zln(7,jj)
          zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(5,j)*zln(4,jj)-temp(8,j)*zln(7,jj)
          zz(3)=zz(3)-temp(3,j)*zln(1,jj)-temp(6,j)*zln(4,jj)-temp(9,j)*zln(7,jj)
          zz(4)=zz(4)-temp(1,j)*zln(2,jj)-temp(4,j)*zln(5,jj)-temp(7,j)*zln(8,jj)
          zz(5)=zz(5)-temp(2,j)*zln(2,jj)-temp(5,j)*zln(5,jj)-temp(8,j)*zln(8,jj)
          zz(6)=zz(6)-temp(3,j)*zln(2,jj)-temp(6,j)*zln(5,jj)-temp(9,j)*zln(8,jj)
          zz(7)=zz(7)-temp(1,j)*zln(3,jj)-temp(4,j)*zln(6,jj)-temp(7,j)*zln(9,jj)
          zz(8)=zz(8)-temp(2,j)*zln(3,jj)-temp(5,j)*zln(6,jj)-temp(8,j)*zln(9,jj)
          zz(9)=zz(9)-temp(3,j)*zln(3,jj)-temp(6,j)*zln(6,jj)-temp(9,j)*zln(9,jj)
        ENDIF
      ENDDO
      CALL inv33(zln(1,k),zz,diag(1,jc))
      !$dir max_trips(9)
      DO l=1,9
        temp(l,jc)=zz(l)
      ENDDO
      t(1)=t(1)+zz(1)*zln(1,k)+zz(4)*zln(4,k)+zz(7)*zln(7,k)
      t(2)=t(2)+zz(1)*zln(2,k)+zz(4)*zln(5,k)+zz(7)*zln(8,k)
      t(3)=t(3)+zz(2)*zln(2,k)+zz(5)*zln(5,k)+zz(8)*zln(8,k)
      t(4)=t(4)+zz(1)*zln(3,k)+zz(4)*zln(6,k)+zz(7)*zln(9,k)
      t(5)=t(5)+zz(2)*zln(3,k)+zz(5)*zln(6,k)+zz(8)*zln(9,k)
      t(6)=t(6)+zz(3)*zln(3,k)+zz(6)*zln(6,k)+zz(9)*zln(9,k)
    ENDDO
    !$dir max_trips(6)
    DO l=1,6
      diag(l,ic)=diag(l,ic)-t(l)
    ENDDO

    CALL inv3(diag(1,ic),ir)
    nch(ic)=-1
    kk=par(ic)
    nch(kk)=nch(kk)-1
    RETURN
  END SUBROUTINE s3um

  !======================================================================!
  !> s3um1
  !======================================================================!
  SUBROUTINE s3um1(ic,xlnzr,colno,zln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: zln

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(9):: s
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    !$dir max_trip(9)
    DO l=1,9
      s(l)=0.0d0
    ENDDO

    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(4,j)*zln(4,jj)+temp(7,j)*zln(7,jj)
          s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(5,j)*zln(4,jj)+temp(8,j)*zln(7,jj)
          s(3)=s(3)+temp(3,j)*zln(1,jj)+temp(6,j)*zln(4,jj)+temp(9,j)*zln(7,jj)
          s(4)=s(4)+temp(1,j)*zln(2,jj)+temp(4,j)*zln(5,jj)+temp(7,j)*zln(8,jj)
          s(5)=s(5)+temp(2,j)*zln(2,jj)+temp(5,j)*zln(5,jj)+temp(8,j)*zln(8,jj)
          s(6)=s(6)+temp(3,j)*zln(2,jj)+temp(6,j)*zln(5,jj)+temp(9,j)*zln(8,jj)
          s(7)=s(7)+temp(1,j)*zln(3,jj)+temp(4,j)*zln(6,jj)+temp(7,j)*zln(9,jj)
          s(8)=s(8)+temp(2,j)*zln(3,jj)+temp(5,j)*zln(6,jj)+temp(8,j)*zln(9,jj)
          s(9)=s(9)+temp(3,j)*zln(3,jj)+temp(6,j)*zln(6,jj)+temp(9,j)*zln(9,jj)
        ENDIF
      ENDDO
      !$dir max_trip(9)
      DO l=1,9
        temp(l,jc)=zln(l,k)-s(l)
        zln(l,k)=temp(l,jc)
        s(l)=0.0d0
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE s3um1

  !======================================================================!
  !> s3um2
  !======================================================================!
  SUBROUTINE s3um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: zln

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

    joc=0
    DO ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      DO k=ks,ke
        jj=colno(k)
        temp(jj,1)=zln(1,k)
        temp(jj,2)=zln(2,k)
        temp(jj,3)=zln(3,k)
        temp(jj,4)=zln(4,k)
        temp(jj,5)=zln(5,k)
        temp(jj,6)=zln(6,k)
        temp(jj,7)=zln(7,k)
        temp(jj,8)=zln(8,k)
        temp(jj,9)=zln(9,k)
        indx(jj)=ic
      ENDDO
      DO k=ks,ke
        jj=colno(k)
        zln(4,k)=temp(jj,4)-temp(jj,1)*diag(2,jj)
        zln(7,k)=temp(jj,7)-temp(jj,1)*diag(4,jj)-zln(4,k)*diag(5,jj)
        zln(1,k)=temp(jj,1)*diag(1,jj)
        zln(4,k)=zln(4,k)*diag(3,jj)
        zln(7,k)=zln(7,k)*diag(6,jj)
        zln(4,k)=zln(4,k)-zln(7,k)*diag(5,jj)
        zln(1,k)=zln(1,k)-zln(4,k)*diag(2,jj)-zln(7,k)*diag(4,jj)
        !
        zln(5,k)=temp(jj,5)-temp(jj,2)*diag(2,jj)
        zln(8,k)=temp(jj,8)-temp(jj,2)*diag(4,jj)-zln(5,k)*diag(5,jj)
        zln(2,k)=temp(jj,2)*diag(1,jj)
        zln(5,k)=zln(5,k)*diag(3,jj)
        zln(8,k)=zln(8,k)*diag(6,jj)
        zln(5,k)=zln(5,k)-zln(8,k)*diag(5,jj)
        zln(2,k)=zln(2,k)-zln(5,k)*diag(2,jj)-zln(8,k)*diag(4,jj)
        !
        zln(6,k)=temp(jj,6)-temp(jj,3)*diag(2,jj)
        zln(9,k)=temp(jj,9)-temp(jj,3)*diag(4,jj)-zln(6,k)*diag(5,jj)
        zln(3,k)=temp(jj,3)*diag(1,jj)
        zln(6,k)=zln(6,k)*diag(3,jj)
        zln(9,k)=zln(9,k)*diag(6,jj)
        zln(6,k)=zln(6,k)-zln(9,k)*diag(5,jj)
        zln(3,k)=zln(3,k)-zln(6,k)*diag(2,jj)-zln(9,k)*diag(4,jj)
      ENDDO

      DO k=ks,ke
        jj=colno(k)
        diag(1,ic)=diag(1,ic)-temp(jj,1)*zln(1,k)-temp(jj,4)*zln(4,k)-temp(jj,7)*zln(7,k)
        diag(2,ic)=diag(2,ic)-temp(jj,1)*zln(2,k)-temp(jj,4)*zln(5,k)-temp(jj,7)*zln(8,k)
        diag(3,ic)=diag(3,ic)-temp(jj,2)*zln(2,k)-temp(jj,5)*zln(5,k)-temp(jj,8)*zln(8,k)
        diag(4,ic)=diag(4,ic)-temp(jj,1)*zln(3,k)-temp(jj,4)*zln(6,k)-temp(jj,7)*zln(9,k)
        diag(5,ic)=diag(5,ic)-temp(jj,2)*zln(3,k)-temp(jj,5)*zln(6,k)-temp(jj,8)*zln(9,k)
        diag(6,ic)=diag(6,ic)-temp(jj,3)*zln(3,k)-temp(jj,6)*zln(6,k)-temp(jj,9)*zln(9,k)
      ENDDO
      DO jc=nstop,ic-1
        joc=joc+1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        DO jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          IF(indx(j)==ic) THEN
            dsln(1,joc)=dsln(1,joc)-temp(j,1)*zln(1,jj)-temp(j,4)*zln(4,jj)-temp(j,7)*zln(7,jj)
            dsln(2,joc)=dsln(2,joc)-temp(j,2)*zln(1,jj)-temp(j,5)*zln(4,jj)-temp(j,8)*zln(7,jj)
            dsln(3,joc)=dsln(3,joc)-temp(j,3)*zln(1,jj)-temp(j,6)*zln(4,jj)-temp(j,9)*zln(7,jj)
            dsln(4,joc)=dsln(4,joc)-temp(j,1)*zln(2,jj)-temp(j,4)*zln(5,jj)-temp(j,7)*zln(8,jj)
            dsln(5,joc)=dsln(5,joc)-temp(j,2)*zln(2,jj)-temp(j,5)*zln(5,jj)-temp(j,8)*zln(8,jj)
            dsln(6,joc)=dsln(6,joc)-temp(j,3)*zln(2,jj)-temp(j,6)*zln(5,jj)-temp(j,9)*zln(8,jj)
            dsln(7,joc)=dsln(7,joc)-temp(j,1)*zln(3,jj)-temp(j,4)*zln(6,jj)-temp(j,7)*zln(9,jj)
            dsln(8,joc)=dsln(8,joc)-temp(j,2)*zln(3,jj)-temp(j,5)*zln(6,jj)-temp(j,8)*zln(9,jj)
            dsln(9,joc)=dsln(9,joc)-temp(j,3)*zln(3,jj)-temp(j,6)*zln(6,jj)-temp(j,9)*zln(9,jj)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE s3um2

  !======================================================================!
  !> s3um3
  !======================================================================!
  SUBROUTINE s3um3(n,dsln,diag,indx,temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: dsln

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    DOUBLE PRECISION, DIMENSION(9):: t

    IF(n<=0) RETURN

    indx(1)=0
    joc=1
    CALL inv3(diag(1,1),ir)
    DO i=2,n
      indx(i)=joc
      DO j=1,i-1
        CALL d3dot(t,dsln(1,indx(i)),dsln(1,indx(j)),j-1)
        !$dir max_trips(9)
        dsln(:,joc)=dsln(:,joc)-t(:)
        joc=joc+1
      ENDDO
      CALL v3prod(dsln(1,indx(i)),diag,temp,i-1)
      CALL d3dotl(t,temp,dsln(1,indx(i)),i-1)
      !$dir max_trips(6)
      diag(:,i)=diag(:,i)-t(1:6)
      CALL vcopy(temp,dsln(1,indx(i)),9*(i-1))
      CALL inv3(diag(1,i),ir)
    ENDDO
    RETURN
  END SUBROUTINE s3um3

  !======================================================================!
  !> s6pdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE s6pdot(bi,b,zln,colno,ks,ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ke
    INTEGER, INTENT(IN):: ks
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    DOUBLE PRECISION, INTENT(IN), DIMENSION(36,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(6,*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6):: bi

    INTEGER:: j
    INTEGER:: jj

    DO jj=ks,ke
      j=colno(jj)
      bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(7,jj)*b(2,j)-zln(13,jj)*b(3,j)-zln(19,jj)*b(4,j)-zln(25,jj)*b(5,j)-zln(31,jj)*b(6,j)
      bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(8,jj)*b(2,j)-zln(14,jj)*b(3,j)-zln(20,jj)*b(4,j)-zln(26,jj)*b(5,j)-zln(32,jj)*b(6,j)
      bi(3)=bi(3)-zln(3,jj)*b(1,j)-zln(9,jj)*b(2,j)-zln(15,jj)*b(3,j)-zln(21,jj)*b(4,j)-zln(27,jj)*b(5,j)-zln(33,jj)*b(6,j)
      bi(4)=bi(4)-zln(4,jj)*b(1,j)-zln(10,jj)*b(2,j)-zln(16,jj)*b(3,j)-zln(22,jj)*b(4,j)-zln(28,jj)*b(5,j)-zln(34,jj)*b(6,j)
      bi(5)=bi(5)-zln(5,jj)*b(1,j)-zln(11,jj)*b(2,j)-zln(17,jj)*b(3,j)-zln(23,jj)*b(4,j)-zln(29,jj)*b(5,j)-zln(35,jj)*b(6,j)
      bi(6)=bi(6)-zln(6,jj)*b(1,j)-zln(12,jj)*b(2,j)-zln(18,jj)*b(3,j)-zln(25,jj)*b(4,j)-zln(30,jj)*b(5,j)-zln(36,jj)*b(6,j)
    ENDDO
    RETURN
  END SUBROUTINE s6pdot

  !======================================================================!
  !> s6um
  !======================================================================!
  SUBROUTINE s6um(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: par

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(21,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,*):: zln

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: lratio

    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(21):: t
    DOUBLE PRECISION, DIMENSION(36):: zz
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    DO l=1,21
      t(l)=0.0d0
    ENDDO
    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      DO l=1,36
        zz(l)=zln(l,k)
      ENDDO
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(7,j)*zln(7,jj)-temp(13,j)*zln(13,jj)-temp(19,j)*zln(19,jj) &
               -temp(25,j)*zln(25,jj)-temp(31,j)*zln(31,jj)
          zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(8,j)*zln(7,jj)-temp(14,j)*zln(13,jj)-temp(20,j)*zln(19,jj) &
               -temp(26,j)*zln(25,jj)-temp(32,j)*zln(31,jj)
          zz(3)=zz(3)-temp(3,j)*zln(1,jj)-temp(9,j)*zln(7,jj)-temp(15,j)*zln(13,jj)-temp(21,j)*zln(19,jj) &
               -temp(27,j)*zln(25,jj)-temp(33,j)*zln(31,jj)
          zz(4)=zz(4)-temp(4,j)*zln(1,jj)-temp(10,j)*zln(7,jj)-temp(16,j)*zln(13,jj)-temp(22,j)*zln(19,jj) &
               -temp(28,j)*zln(25,jj)-temp(34,j)*zln(31,jj)
          zz(5)=zz(5)-temp(5,j)*zln(1,jj)-temp(11,j)*zln(7,jj)-temp(17,j)*zln(13,jj)-temp(23,j)*zln(19,jj) &
               -temp(29,j)*zln(25,jj)-temp(35,j)*zln(31,jj)
          zz(6)=zz(6)-temp(6,j)*zln(1,jj)-temp(12,j)*zln(7,jj)-temp(18,j)*zln(13,jj)-temp(24,j)*zln(19,jj) &
               -temp(30,j)*zln(25,jj)-temp(36,j)*zln(31,jj)
          zz(7)=zz(7)-temp(1,j)*zln(2,jj)-temp(7,j)*zln(8,jj)-temp(13,j)*zln(14,jj)-temp(19,j)*zln(20,jj) &
               -temp(25,j)*zln(26,jj)-temp(31,j)*zln(32,jj)
          zz(8)=zz(8)-temp(2,j)*zln(2,jj)-temp(8,j)*zln(8,jj)-temp(14,j)*zln(14,jj)-temp(20,j)*zln(20,jj) &
               -temp(26,j)*zln(26,jj)-temp(32,j)*zln(32,jj)
          zz(9)=zz(9)-temp(3,j)*zln(2,jj)-temp(9,j)*zln(8,jj)-temp(15,j)*zln(14,jj)-temp(21,j)*zln(20,jj) &
               -temp(27,j)*zln(26,jj)-temp(33,j)*zln(32,jj)
          zz(10)=zz(10)-temp(4,j)*zln(2,jj)-temp(10,j)*zln(8,jj)-temp(16,j)*zln(14,jj)-temp(22,j)*zln(20,jj) &
               -temp(28,j)*zln(26,jj)-temp(34,j)*zln(32,jj)
          zz(11)=zz(11)-temp(5,j)*zln(2,jj)-temp(11,j)*zln(8,jj)-temp(17,j)*zln(14,jj)-temp(23,j)*zln(20,jj) &
               -temp(29,j)*zln(26,jj)-temp(35,j)*zln(32,jj)
          zz(12)=zz(12)-temp(6,j)*zln(2,jj)-temp(12,j)*zln(8,jj)-temp(18,j)*zln(14,jj)-temp(24,j)*zln(20,jj) &
               -temp(30,j)*zln(26,jj)-temp(36,j)*zln(32,jj)
          zz(13)=zz(13)-temp(1,j)*zln(3,jj)-temp(7,j)*zln(9,jj)-temp(13,j)*zln(15,jj)-temp(19,j)*zln(21,jj) &
               -temp(25,j)*zln(27,jj)-temp(31,j)*zln(33,jj)
          zz(14)=zz(14)-temp(2,j)*zln(3,jj)-temp(8,j)*zln(9,jj)-temp(14,j)*zln(15,jj)-temp(20,j)*zln(21,jj) &
               -temp(26,j)*zln(27,jj)-temp(32,j)*zln(33,jj)
          zz(15)=zz(15)-temp(3,j)*zln(3,jj)-temp(9,j)*zln(9,jj)-temp(15,j)*zln(15,jj)-temp(21,j)*zln(21,jj) &
               -temp(27,j)*zln(27,jj)-temp(33,j)*zln(33,jj)
          zz(16)=zz(16)-temp(4,j)*zln(3,jj)-temp(10,j)*zln(9,jj)-temp(16,j)*zln(15,jj)-temp(22,j)*zln(21,jj) &
               -temp(28,j)*zln(27,jj)-temp(34,j)*zln(33,jj)
          zz(17)=zz(17)-temp(5,j)*zln(3,jj)-temp(11,j)*zln(9,jj)-temp(17,j)*zln(15,jj)-temp(23,j)*zln(21,jj) &
               -temp(29,j)*zln(27,jj)-temp(35,j)*zln(33,jj)
          zz(18)=zz(18)-temp(6,j)*zln(3,jj)-temp(12,j)*zln(9,jj)-temp(18,j)*zln(15,jj)-temp(24,j)*zln(21,jj) &
               -temp(30,j)*zln(27,jj)-temp(36,j)*zln(33,jj)
          zz(19)=zz(19)-temp(1,j)*zln(4,jj)-temp(7,j)*zln(10,jj)-temp(13,j)*zln(16,jj)-temp(19,j)*zln(22,jj) &
               -temp(25,j)*zln(28,jj)-temp(31,j)*zln(34,jj)
          zz(20)=zz(20)-temp(2,j)*zln(4,jj)-temp(8,j)*zln(10,jj)-temp(14,j)*zln(16,jj)-temp(20,j)*zln(22,jj) &
               -temp(26,j)*zln(28,jj)-temp(32,j)*zln(34,jj)
          zz(21)=zz(21)-temp(3,j)*zln(4,jj)-temp(9,j)*zln(10,jj)-temp(15,j)*zln(16,jj)-temp(21,j)*zln(22,jj) &
               -temp(27,j)*zln(28,jj)-temp(33,j)*zln(34,jj)
          zz(22)=zz(22)-temp(4,j)*zln(4,jj)-temp(10,j)*zln(10,jj)-temp(16,j)*zln(16,jj)-temp(22,j)*zln(22,jj) &
               -temp(28,j)*zln(28,jj)-temp(34,j)*zln(34,jj)
          zz(23)=zz(23)-temp(5,j)*zln(4,jj)-temp(11,j)*zln(10,jj)-temp(17,j)*zln(16,jj)-temp(23,j)*zln(22,jj) &
               -temp(29,j)*zln(28,jj)-temp(35,j)*zln(34,jj)
          zz(24)=zz(24)-temp(6,j)*zln(4,jj)-temp(12,j)*zln(10,jj)-temp(18,j)*zln(16,jj)-temp(24,j)*zln(22,jj) &
               -temp(30,j)*zln(28,jj)-temp(36,j)*zln(34,jj)
          zz(25)=zz(25)-temp(1,j)*zln(5,jj)-temp(7,j)*zln(11,jj)-temp(13,j)*zln(17,jj)-temp(19,j)*zln(23,jj) &
               -temp(25,j)*zln(29,jj)-temp(31,j)*zln(35,jj)
          zz(26)=zz(26)-temp(2,j)*zln(5,jj)-temp(8,j)*zln(11,jj)-temp(14,j)*zln(17,jj)-temp(20,j)*zln(23,jj) &
               -temp(26,j)*zln(29,jj)-temp(32,j)*zln(35,jj)
          zz(27)=zz(27)-temp(3,j)*zln(5,jj)-temp(9,j)*zln(11,jj)-temp(15,j)*zln(17,jj)-temp(21,j)*zln(23,jj) &
               -temp(27,j)*zln(29,jj)-temp(33,j)*zln(35,jj)
          zz(28)=zz(28)-temp(4,j)*zln(5,jj)-temp(10,j)*zln(11,jj)-temp(16,j)*zln(17,jj)-temp(22,j)*zln(23,jj) &
               -temp(28,j)*zln(29,jj)-temp(34,j)*zln(35,jj)
          zz(29)=zz(29)-temp(5,j)*zln(5,jj)-temp(11,j)*zln(11,jj)-temp(17,j)*zln(17,jj)-temp(23,j)*zln(23,jj) &
               -temp(29,j)*zln(29,jj)-temp(35,j)*zln(35,jj)
          zz(30)=zz(30)-temp(6,j)*zln(5,jj)-temp(12,j)*zln(11,jj)-temp(18,j)*zln(17,jj)-temp(24,j)*zln(23,jj) &
               -temp(30,j)*zln(29,jj)-temp(36,j)*zln(35,jj)
          zz(31)=zz(31)-temp(1,j)*zln(6,jj)-temp(7,j)*zln(12,jj)-temp(13,j)*zln(18,jj)-temp(19,j)*zln(24,jj) &
               -temp(25,j)*zln(30,jj)-temp(31,j)*zln(36,jj)
          zz(32)=zz(32)-temp(2,j)*zln(6,jj)-temp(8,j)*zln(12,jj)-temp(14,j)*zln(18,jj)-temp(20,j)*zln(24,jj) &
               -temp(26,j)*zln(30,jj)-temp(32,j)*zln(36,jj)
          zz(33)=zz(33)-temp(3,j)*zln(6,jj)-temp(9,j)*zln(12,jj)-temp(15,j)*zln(18,jj)-temp(21,j)*zln(24,jj) &
               -temp(27,j)*zln(30,jj)-temp(33,j)*zln(36,jj)
          zz(34)=zz(34)-temp(4,j)*zln(6,jj)-temp(10,j)*zln(12,jj)-temp(16,j)*zln(18,jj)-temp(22,j)*zln(24,jj) &
               -temp(28,j)*zln(30,jj)-temp(34,j)*zln(36,jj)
          zz(35)=zz(35)-temp(5,j)*zln(6,jj)-temp(11,j)*zln(12,jj)-temp(17,j)*zln(18,jj)-temp(23,j)*zln(24,jj) &
               -temp(29,j)*zln(30,jj)-temp(35,j)*zln(36,jj)
          zz(36)=zz(36)-temp(6,j)*zln(6,jj)-temp(12,j)*zln(12,jj)-temp(18,j)*zln(18,jj)-temp(24,j)*zln(24,jj) &
               -temp(30,j)*zln(30,jj)-temp(36,j)*zln(36,jj)
        ENDIF
      ENDDO
      CALL inv66(zln(1,k),zz,diag(1,jc))
      !$dir max_trips(9)
      DO l=1,36
        temp(l,jc)=zz(l)
      ENDDO
      t(1)=t(1)+zz(1)*zln(1,k)+zz(7)*zln(7,k)+zz(13)*zln(13,k)+zz(19)*zln(19,k)+zz(25)*zln(25,k)+zz(31)*zln(31,k)
      t(2)=t(2)+zz(1)*zln(2,k)+zz(7)*zln(8,k)+zz(13)*zln(14,k)+zz(19)*zln(20,k)+zz(25)*zln(26,k)+zz(31)*zln(32,k)
      t(3)=t(3)+zz(2)*zln(2,k)+zz(8)*zln(8,k)+zz(14)*zln(14,k)+zz(20)*zln(20,k)+zz(26)*zln(26,k)+zz(32)*zln(32,k)
      t(4)=t(4)+zz(1)*zln(3,k)+zz(7)*zln(9,k)+zz(13)*zln(15,k)+zz(19)*zln(21,k)+zz(25)*zln(27,k)+zz(31)*zln(33,k)
      t(5)=t(5)+zz(2)*zln(3,k)+zz(8)*zln(9,k)+zz(14)*zln(15,k)+zz(20)*zln(21,k)+zz(26)*zln(27,k)+zz(32)*zln(33,k)
      t(6)=t(6)+zz(3)*zln(3,k)+zz(9)*zln(9,k)+zz(15)*zln(15,k)+zz(21)*zln(21,k)+zz(27)*zln(27,k)+zz(33)*zln(33,k)
      t(7)=t(7)+zz(1)*zln(4,k)+zz(7)*zln(10,k)+zz(13)*zln(16,k)+zz(19)*zln(22,k)+zz(25)*zln(28,k)+zz(31)*zln(34,k)
      t(8)=t(8)+zz(2)*zln(4,k)+zz(8)*zln(10,k)+zz(14)*zln(16,k)+zz(20)*zln(22,k)+zz(26)*zln(28,k)+zz(32)*zln(34,k)
      t(9)=t(9)+zz(3)*zln(4,k)+zz(9)*zln(10,k)+zz(15)*zln(16,k)+zz(21)*zln(22,k)+zz(27)*zln(28,k)+zz(33)*zln(34,k)
      t(10)=t(10)+zz(4)*zln(4,k)+zz(10)*zln(10,k)+zz(16)*zln(16,k)+zz(22)*zln(22,k)+zz(28)*zln(28,k)+zz(34)*zln(34,k)
      t(11)=t(11)+zz(1)*zln(5,k)+zz(7)*zln(11,k)+zz(13)*zln(17,k)+zz(19)*zln(23,k)+zz(25)*zln(29,k)+zz(31)*zln(35,k)
      t(12)=t(12)+zz(2)*zln(5,k)+zz(8)*zln(11,k)+zz(14)*zln(17,k)+zz(20)*zln(23,k)+zz(26)*zln(29,k)+zz(32)*zln(35,k)
      t(13)=t(13)+zz(3)*zln(5,k)+zz(9)*zln(11,k)+zz(15)*zln(17,k)+zz(21)*zln(23,k)+zz(27)*zln(29,k)+zz(33)*zln(35,k)
      t(14)=t(14)+zz(4)*zln(5,k)+zz(10)*zln(11,k)+zz(16)*zln(17,k)+zz(22)*zln(23,k)+zz(28)*zln(29,k)+zz(34)*zln(35,k)
      t(15)=t(15)+zz(5)*zln(5,k)+zz(11)*zln(11,k)+zz(17)*zln(17,k)+zz(23)*zln(23,k)+zz(29)*zln(29,k)+zz(35)*zln(35,k)
      t(16)=t(16)+zz(1)*zln(6,k)+zz(7)*zln(12,k)+zz(13)*zln(18,k)+zz(19)*zln(24,k)+zz(25)*zln(30,k)+zz(31)*zln(36,k)
      t(17)=t(17)+zz(2)*zln(6,k)+zz(8)*zln(12,k)+zz(14)*zln(18,k)+zz(20)*zln(24,k)+zz(26)*zln(30,k)+zz(32)*zln(36,k)
      t(18)=t(18)+zz(3)*zln(6,k)+zz(9)*zln(12,k)+zz(15)*zln(18,k)+zz(21)*zln(24,k)+zz(27)*zln(30,k)+zz(33)*zln(36,k)
      t(19)=t(19)+zz(4)*zln(6,k)+zz(10)*zln(12,k)+zz(16)*zln(18,k)+zz(22)*zln(24,k)+zz(28)*zln(30,k)+zz(34)*zln(36,k)
      t(20)=t(20)+zz(5)*zln(6,k)+zz(11)*zln(12,k)+zz(17)*zln(18,k)+zz(23)*zln(24,k)+zz(29)*zln(30,k)+zz(35)*zln(36,k)
      t(21)=t(21)+zz(6)*zln(6,k)+zz(12)*zln(12,k)+zz(18)*zln(18,k)+zz(24)*zln(24,k)+zz(30)*zln(30,k)+zz(36)*zln(36,k)
    ENDDO
    !$dir max_trips(6)
    DO l=1,21
      diag(l,ic)=diag(l,ic)-t(l)
    ENDDO
    CALL inv6(diag(1,ic),ir)
    nch(ic)=-1
    kk=par(ic)
    nch(kk)=nch(kk)-1
    RETURN
  END SUBROUTINE s6um

  !======================================================================!
  !> s6um1
  !======================================================================!
  SUBROUTINE s6um1(ic,xlnzr,colno,zln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: zln

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: l
    INTEGER:: lratio
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION, DIMENSION(9):: s
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    !$dir max_trip(9)
    DO l=1,9
      s(l)=0.0d0
    ENDDO
    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(4,j)*zln(4,jj)+temp(7,j)*zln(7,jj)
          s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(5,j)*zln(4,jj)+temp(8,j)*zln(7,jj)
          s(3)=s(3)+temp(3,j)*zln(1,jj)+temp(6,j)*zln(4,jj)+temp(9,j)*zln(7,jj)
          s(4)=s(4)+temp(1,j)*zln(2,jj)+temp(4,j)*zln(5,jj)+temp(7,j)*zln(8,jj)
          s(5)=s(5)+temp(2,j)*zln(2,jj)+temp(5,j)*zln(5,jj)+temp(8,j)*zln(8,jj)
          s(6)=s(6)+temp(3,j)*zln(2,jj)+temp(6,j)*zln(5,jj)+temp(9,j)*zln(8,jj)
          s(7)=s(7)+temp(1,j)*zln(3,jj)+temp(4,j)*zln(6,jj)+temp(7,j)*zln(9,jj)
          s(8)=s(8)+temp(2,j)*zln(3,jj)+temp(5,j)*zln(6,jj)+temp(8,j)*zln(9,jj)
          s(9)=s(9)+temp(3,j)*zln(3,jj)+temp(6,j)*zln(6,jj)+temp(9,j)*zln(9,jj)
        ENDIF
      ENDDO
      !$dir max_trip(9)
      DO l=1,9
        temp(l,jc)=zln(l,k)-s(l)
        zln(l,k)=temp(l,jc)
        s(l)=0.0d0
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE s6um1

  !======================================================================!
  !> s6um2
  !======================================================================!
  SUBROUTINE s6um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(21,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,neqns):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(36,*):: zln

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

    joc=0
    DO ic=nstop,neqns
      DO m=1,36
        DO jj=1,nstop
          temp(jj,m)=0.0d0
        ENDDO
      ENDDO
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      DO k=ks,ke
        jj=colno(k)
        temp(:,jj)=zln(:,k)
        indx(jj)=ic
      ENDDO
      DO k=ks,ke
        jj=colno(k)
        CALL inv66(zln(1,k),temp,diag(1,jj))
      ENDDO

      DO k=ks,ke
        jj=colno(k)
        diag(1,ic)=diag(1,ic)-temp(jj,1)*zln(1,k)-temp(jj,4)*zln(4,k)-temp(jj,7)*zln(7,k)
        diag(2,ic)=diag(2,ic)-temp(jj,1)*zln(2,k)-temp(jj,4)*zln(5,k)-temp(jj,7)*zln(8,k)
        diag(3,ic)=diag(3,ic)-temp(jj,2)*zln(2,k)-temp(jj,5)*zln(5,k)-temp(jj,8)*zln(8,k)
        diag(4,ic)=diag(4,ic)-temp(jj,1)*zln(3,k)-temp(jj,4)*zln(6,k)-temp(jj,7)*zln(9,k)
        diag(5,ic)=diag(5,ic)-temp(jj,2)*zln(3,k)-temp(jj,5)*zln(6,k)-temp(jj,8)*zln(9,k)
        diag(6,ic)=diag(6,ic)-temp(jj,3)*zln(3,k)-temp(jj,6)*zln(6,k)-temp(jj,9)*zln(9,k)
      ENDDO
      DO jc=nstop,ic-1
        joc=joc+1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        DO jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          IF(indx(j)==ic) THEN
            dsln(1,joc)=dsln(1,joc)-temp(j,1)*zln(1,jj)-temp(j,4)*zln(4,jj)-temp(j,7)*zln(7,jj)
            dsln(2,joc)=dsln(2,joc)-temp(j,2)*zln(1,jj)-temp(j,5)*zln(4,jj)-temp(j,8)*zln(7,jj)
            dsln(3,joc)=dsln(3,joc)-temp(j,3)*zln(1,jj)-temp(j,6)*zln(4,jj)-temp(j,9)*zln(7,jj)
            dsln(4,joc)=dsln(4,joc)-temp(j,1)*zln(2,jj)-temp(j,4)*zln(5,jj)-temp(j,7)*zln(8,jj)
            dsln(5,joc)=dsln(5,joc)-temp(j,2)*zln(2,jj)-temp(j,5)*zln(5,jj)-temp(j,8)*zln(8,jj)
            dsln(6,joc)=dsln(6,joc)-temp(j,3)*zln(2,jj)-temp(j,6)*zln(5,jj)-temp(j,9)*zln(8,jj)
            dsln(7,joc)=dsln(7,joc)-temp(j,1)*zln(3,jj)-temp(j,4)*zln(6,jj)-temp(j,7)*zln(9,jj)
            dsln(8,joc)=dsln(8,joc)-temp(j,2)*zln(3,jj)-temp(j,5)*zln(6,jj)-temp(j,8)*zln(9,jj)
            dsln(9,joc)=dsln(9,joc)-temp(j,3)*zln(3,jj)-temp(j,6)*zln(6,jj)-temp(j,9)*zln(9,jj)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE s6um2

  !======================================================================!
  !> s6um3
  !======================================================================!
  SUBROUTINE s6um3(n,dsln,diag,indx,temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n

    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,*):: dsln

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: l
    DOUBLE PRECISION, DIMENSION(9):: t

    IF(n<=0) RETURN

    indx(1)=0
    joc=1
    CALL inv3(diag(1,1),ir)
    DO i=2,n
      indx(i)=joc
      DO j=1,i-1
        CALL d3dot(t,dsln(1,indx(i)),dsln(1,indx(j)),j-1)
        !$dir max_trips(9)
        DO l=1,9
          dsln(l,joc)=dsln(l,joc)-t(l)
        ENDDO
        joc=joc+1
      ENDDO
      CALL v3prod(dsln(1,indx(i)),diag,temp,i-1)
      CALL d3dotl(t,temp,dsln(1,indx(i)),i-1)
      !$dir max_trips(6)
      DO l=1,6
        diag(l,i)=diag(l,i)-t(l)
      ENDDO
      CALL vcopy(temp,dsln(1,indx(i)),9*(i-1))
      CALL inv3(diag(1,i),ir)
    ENDDO
    RETURN
  END SUBROUTINE s6um3

  !======================================================================!
  !> setij
  !======================================================================!
  SUBROUTINE setij(hecMESH,hecMAT)
    USE hecmw_util
    IMPLICIT NONE

    TYPE (hecmwST_local_mesh), INTENT(IN):: hecMESH
    TYPE (hecmwST_matrix    ), INTENT(IN):: hecMAT

    INTEGER:: i
    INTEGER:: ierr
    INTEGER:: j
    INTEGER:: k
    INTEGER:: kk
    INTEGER:: ntotal
    INTEGER(KIND=kint):: numnp
    INTEGER(KIND=kint):: ndof
    INTEGER(KIND=kint):: ndof2

    NUMNP = hecMAT%NP
    NDOF  = hecMESH%n_dof
    ntotal = NUMNP*NDOF

    !*NUFACT variables
    neqns = NUMNP
    ndeg = NDOF
    nttbr = hecMAT%NP+hecMAT%NPL !+hecMAT%NPU if unsymmetric
    isym = 0

    !*Allocations
    ALLOCATE(irow(nttbr),STAT=ierr)
    ALLOCATE(jcol(nttbr),STAT=ierr)
    IF(ierr/=0) STOP "Allocation error: irow/jcol"

    kk = 0
    ndof2 = NDOF*NDOF
    DO j= 1, NUMNP
      !*Diagonal
      kk = kk + 1
      irow(kk) = j
      jcol(kk) = j
      !*Lower
      DO k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i= hecMAT%itemL(k)
        kk = kk + 1
        irow(kk) = j
        jcol(kk) = i
      ENDDO

      !!*Upper if unsymmetric
    ENDDO

    RETURN
  END SUBROUTINE setij

  !======================================================================!
  !> staij1
  !      this routine sets an non-zero entry  of the matrix.
  !      (symmetric version)
  !
  !      (i)
  !          isw      =0    set the value
  !                   =1    add the value
  !          i        row entry
  !          j        column entry
  !          aij      value
  !
  !      (o)
  !          iv       communication array
  !
  !        #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE staij1(isw,i,j,aij,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: i
    INTEGER, INTENT(IN):: isw
    INTEGER, INTENT(IN):: j
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg*ndeg):: aij

    INTEGER, INTENT(OUT):: ir

    INTEGER:: idbg
    INTEGER:: lratio
    INTEGER:: ndeg2
    INTEGER:: ndeg2l
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    COMMON /debug/ idbg
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ir=0
    ndeg2=ndeg*ndeg
    ndeg2l=ndeg*(ndeg+1)/2
    IF(stage==30) WRITE(6,*) 'warning a matrix was build up '//'but never solved.'
    IF(stage==10) THEN
      ALLOCATE(diag(neqns*ndeg2l),STAT = ierror)
      raloc = raloc + neqns*ndeg2l
      ALLOCATE(zln(len_colno*ndeg2),STAT = ierror)
      raloc = raloc + len_colno*ndeg2
      ALLOCATE(dsln(len_dsln*ndeg2),STAT = ierror)
      IF(ierror/=0) STOP "Allocation error dsln"
      raloc = raloc + len_dsln*ndeg2
    ENDIF
    IF(stage/=20) THEN

      ! for diagonal
      diag = 0.

      ! for lower triangle
      zln = 0.

      ! for dense window
      dsln = 0.
      stage=20
    ENDIF
    !         Print *,'********Set Stage 20 *********'
    !
    IF(ndeg<=2) THEN
      CALL addr0(isw,i,j,aij,diag,zln,dsln,nstop,ndeg2,ndeg2l,ir)
    ELSEIF(ndeg==3) THEN
      CALL addr3(i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ir)
    ELSE
      CALL addrx(i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ndeg2l,ir)
    ENDIF
    RETURN
  END SUBROUTINE staij1

  !======================================================================!
  !> stiaja
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE stiaja(neqns,ia,ja,jcpt,jcolno)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN), DIMENSION(*):: jcpt
    INTEGER, INTENT(IN), DIMENSION(*):: jcolno
    INTEGER, INTENT(OUT), DIMENSION(*):: ia
    INTEGER, INTENT(OUT), DIMENSION(*):: ja

    INTEGER:: i
    INTEGER:: idbg
    INTEGER:: ii
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: l
    COMMON /debug/idbg

    ia(1)=1
    l=0
    DO k=1,neqns
      joc=jcpt(k)
      DO
        IF(joc==0) EXIT
        ii=jcolno(joc)
        IF(ii==k) GOTO 130
        l=l+1
        ja(l)=ii

130     CONTINUE
        joc=jcpt(joc)
      ENDDO

      ia(k+1)=l+1
    ENDDO
    IF(idbg/=0) THEN
      WRITE(6,*) 'ia '
      WRITE(6,"(10i7)") (ia(i),i=1,neqns)
      WRITE(6,*) 'ja '
      WRITE(6,"(10i7)") (ja(i),i=1,ia(neqns+1))
    ENDIF
    RETURN
  END SUBROUTINE stiaja

  !======================================================================!
  !> stsmat
  !     #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE stsmat(neqns,nttbr,irow,jcol,jcpt,jcolno)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nttbr
    INTEGER, INTENT(IN), DIMENSION(*):: irow
    INTEGER, INTENT(IN), DIMENSION(*):: jcol
    INTEGER, INTENT(OUT), DIMENSION(*):: jcpt
    INTEGER, INTENT(OUT), DIMENSION(*):: jcolno

    INTEGER:: i
    INTEGER:: idbg
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: l
    INTEGER:: locr
    COMMON /debug/idbg

    jcpt(1:2*nttbr)=0
    jcolno(1:2*nttbr)=0
    DO i=1,neqns
      jcpt(i)=i+neqns
      jcolno(i+neqns)=i
    ENDDO

    k=2*neqns

    DO l=1,nttbr
      i=irow(l)
      j=jcol(l)
      IF(i==j) CYCLE
      joc=jcpt(i)
      locr=i

      DO
        IF(joc==0) THEN
          k=k+1
          jcpt(locr)=k
          jcolno(k)=j
          GOTO 150
        ENDIF

        IF(jcolno(joc)==j) THEN
          GOTO 100
        ELSEIF(jcolno(joc)>j) THEN
          k=k+1
          jcpt(locr)=k
          jcpt(k)=joc
          jcolno(k)=j
          GOTO 150
        ENDIF
        locr=joc
        joc=jcpt(joc)
      ENDDO

150   CONTINUE
      joc=jcpt(j)
      locr=j

      DO
        IF(joc==0) THEN
          k=k+1
          jcpt(locr)=k
          jcolno(k)=i
          GOTO 100
        ENDIF

        IF(jcolno(joc)==i) THEN
          GOTO 100
        ELSEIF(jcolno(joc)>i) THEN
          k=k+1
          jcpt(locr)=k
          jcpt(k)=joc
          jcolno(k)=i
          GOTO 100
        ENDIF
        locr=joc
        joc=jcpt(joc)
      ENDDO
100   CONTINUE
    ENDDO

    IF(idbg/=0) THEN
      WRITE(6,*) 'jcolno'
      WRITE(6,60) (jcolno(i),i=1,k)
      WRITE(6,*) 'jcpt'
      WRITE(6,60) (jcpt(i),i=1,k)
    ENDIF
    RETURN

60  FORMAT(10i7)
  END SUBROUTINE stsmat

  !======================================================================!
  !> sum
  !======================================================================!
  SUBROUTINE SUM(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
    IMPLICIT NONE


    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: par
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: zln

    INTEGER:: isem
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: lratio

    DOUBLE PRECISION:: piv
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: s
    DOUBLE PRECISION:: t
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: zz
    COMMON /mchdpn/ rmax,rmin,epsm,lratio
    COMMON isem

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t=0.0d0
    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      s=0.0d0
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          s=s+temp(j)*zln(jj)
        ENDIF
      ENDDO
      zz=zln(k)-s
      zln(k)=zz*diag(jc)
      temp(jc)=zz
      t=t+zz*zln(k)
    ENDDO
    piv=diag(ic)-t
    IF(dabs(piv)>rmin) THEN
      diag(ic)=1.0d0/piv
    ENDIF

    DO
      IF(isem==1) THEN
        isem=0
        nch(ic)=-1
        kk=par(ic)
        nch(kk)=nch(kk)-1
        isem=1
        EXIT
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE sum

  !======================================================================!
  !> sum1
  !======================================================================!
  SUBROUTINE sum1(ic,xlnzr,colno,zln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: zln

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: ks
    INTEGER:: lratio
    INTEGER, DIMENSION(*):: colno
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    DOUBLE PRECISION:: s
    DOUBLE PRECISION:: t
    DOUBLE PRECISION:: zz
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t=0.0d0

    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      s=0.0d0
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          s=s+temp(j)*zln(jj)
        ENDIF
      ENDDO
      zz=zln(k)-s

      zln(k)=zz
      temp(jc)=zz
    ENDDO
    RETURN
  END SUBROUTINE sum1

  !======================================================================!
  !> sum2
  !======================================================================!
  SUBROUTINE sum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: zln

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

    joc=0
    DO ic=nstop,neqns
      DO i=1,nstop
        temp(i)=0.0d0
      ENDDO
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      DO k=ks,ke
        jj=colno(k)
        temp(jj)=zln(k)
        zln(k)=temp(jj)*diag(jj)
        indx(jj)=ic
        diag(ic)=diag(ic)-temp(jj)*zln(k)
      ENDDO
      DO jc=nstop,ic-1
        s=0.0d0
        joc=joc+1
        DO jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          IF(indx(j)==ic) THEN
            s=s+temp(j)*zln(jj)
          ENDIF
        ENDDO
        IF(s==0.0d0) WRITE(16,*) ic,jc
        dsln(joc)=dsln(joc)-s
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE sum2

  !======================================================================!
  !> sum3
  !======================================================================!
  SUBROUTINE sum3(n,dsln,diag,indx,temp)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(*):: dsln

    INTEGER:: i
    INTEGER:: j
    INTEGER:: joc

    IF(n<=0) RETURN
    indx(1)=0
    joc=1
    diag(1)=1.0d0/diag(1)
    DO i=2,n
      indx(i)=joc
      DO j=1,i-1
        dsln(joc)=dsln(joc)-ddot(dsln(indx(i)),dsln(indx(j)),j-1)
        joc=joc+1
      ENDDO
      CALL vprod(dsln(indx(i)),diag,temp,i-1)
      diag(i)=diag(i)-ddot(temp,dsln(indx(i)),i-1)
      CALL vcopy(temp,dsln(indx(i)),i-1)
      diag(i)=1.0d0/diag(i)
    ENDDO
    RETURN
  END SUBROUTINE sum3

  !======================================================================!
  !> sxpdot
  !      spdot1 performs inner product of sparse vectors
  !      #coded by t.arakawa of RIST on 040510
  !======================================================================!
  SUBROUTINE sxpdot(bi,b,zln,colno,ks,ke)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ke
    INTEGER, INTENT(IN):: ks
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: zln
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,*):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg):: bi

    INTEGER:: j
    INTEGER:: jj
    INTEGER:: m
    INTEGER:: n

    DO jj=ks,ke
      j=colno(jj)
      DO m=1,ndeg
        DO n=1,ndeg
          bi(n)=bi(n)-zln(n,m,jj)*b(m,j)
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE sxpdot

  !======================================================================!
  !> sxum
  !======================================================================!
  SUBROUTINE sxum(ic,xlnzr,colno,zln,diag,nch,par,temp,indx,ndegl,zz,t)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN):: ndegl
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(IN), DIMENSION(*):: par
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    INTEGER, INTENT(OUT), DIMENSION(*):: nch
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndegl,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndegl):: t
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg):: zz

    INTEGER:: ir
    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: joc
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: lratio
    INTEGER:: m
    INTEGER:: n
    INTEGER:: ndeg22
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ndeg22=ndeg*ndeg
    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t=0.0
    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      zz=zln(:,:,k)
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          DO m=1,ndeg,2
            DO n=1,ndeg,2
              DO kk=1,ndeg,2
                zz(n,m)=zz(n,m)-temp(n,kk,j)*zln(m,kk,jj)-temp(n,kk+1,j)*zln(m,kk+1,jj)
                zz(n,m+1)=zz(n,m+1)-temp(n,kk,j)*zln(m+1,kk,jj)-temp(n,kk+1,j)*zln(m+1,kk+1,jj)
                zz(n+1,m)=zz(n+1,m)-temp(n+1,kk,j)*zln(m,kk,jj)-temp(n+1,kk+1,j)*zln(m,kk+1,jj)
                zz(n+1,m+1)=zz(n+1,m+1)-temp(n+1,kk,j)*zln(m+1,kk,jj)-temp(n+1,kk+1,j)*zln(m+1,kk+1,jj)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO

      CALL invxx(zln(1,1,k),zz,diag(1,jc),ndeg)

      temp(:,:,jc)=zz
      joc=0
      DO n=1,ndeg
        DO m=1,n
          joc=joc+1
          DO kk=1,ndeg,2
            t(joc)=t(joc)+zz(n,kk)*zln(m,kk,k)+zz(n,kk+1)*zln(m,kk+1,k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    diag(:,ic)=diag(:,ic)-t
    CALL invx(diag(1,ic),ndeg,ir)
    nch(ic)=-1
    kk=par(ic)
    nch(kk)=nch(kk)-1
    RETURN
  END SUBROUTINE sxum

  !======================================================================!
  !> sxum1
  !======================================================================!
  SUBROUTINE sxum1(ic,xlnzr,colno,zln,temp,indx,s)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ic
    INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
    INTEGER, INTENT(IN), DIMENSION(*):: colno
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg):: s
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: temp

    INTEGER:: j
    INTEGER:: jc
    INTEGER:: jj
    INTEGER:: k
    INTEGER:: ke
    INTEGER:: kk
    INTEGER:: ks
    INTEGER:: lratio
    INTEGER:: m
    INTEGER:: n
    DOUBLE PRECISION:: epsm
    DOUBLE PRECISION:: rmax
    DOUBLE PRECISION:: rmin
    COMMON /mchdpn/ rmax,rmin,epsm,lratio

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    DO m=1,ndeg
      DO n=1,ndeg
        s(n,m)=0.0d0
      ENDDO
    ENDDO

    DO k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      DO jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        IF(indx(j)==ic) THEN
          DO m=1,ndeg
            DO n=1,ndeg
              DO kk=1,ndeg
                s(n,m)=s(n,m)+temp(n,kk,j)*zln(m,kk,jj)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      DO m=1,ndeg
        DO n=1,ndeg
          temp(n,m,jc)=zln(n,m,k)-s(n,m)
          zln(n,m,k)=temp(n,m,jc)
          s(n,m)=0.0d0
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE sxum1

  !======================================================================!
  !> sxum2
  !======================================================================!
  SUBROUTINE sxum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx,ndegl)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nstop
    INTEGER, INTENT(IN):: ndegl
    INTEGER, DIMENSION(*), INTENT(IN):: xlnzr
    INTEGER, DIMENSION(*), INTENT(IN):: colno
    INTEGER, DIMENSION(*), INTENT(OUT):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndegl,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: zln

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

    joc=0
    DO ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      DO k=ks,ke
        jj=colno(k)
        DO m=1,ndeg
          DO n=1,ndeg
            temp(n,m,jj)=zln(n,m,k)
            indx(jj)=ic
          ENDDO
        ENDDO
      ENDDO
      DO k=ks,ke
        jj=colno(k)
        CALL invxx(zln(1,1,k),temp(1,1,jj),diag(1,jj),ndeg)
      ENDDO

      locd=0
      DO n=1,ndeg
        DO m=1,n
          locd=locd+1
          DO k=ks,ke
            jj=colno(k)
            DO kk=1,ndeg
              diag(locd,ic)=diag(locd,ic)-temp(n,kk,jj)*zln(m,kk,k)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO jc=nstop,ic-1
        joc=joc+1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        DO jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          IF(indx(j)==ic) THEN
            DO m=1,ndeg
              DO n=1,ndeg
                DO k=1,ndeg
                  dsln(n,m,joc)=dsln(n,m,joc)-temp(n,k,j)*zln(m,k,jj)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE sxum2

  !======================================================================!
  !> sxum3
  !======================================================================!
  SUBROUTINE sxum3(nn,dsln,diag,indx,temp,ndegl,t)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: ndegl
    INTEGER, INTENT(IN):: nn
    INTEGER, INTENT(OUT), DIMENSION(*):: indx
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: temp
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndegl,*):: diag
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg,*):: dsln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,ndeg):: t

    INTEGER:: i
    INTEGER:: ir
    INTEGER:: j
    INTEGER:: joc
    INTEGER:: locd
    INTEGER:: m
    INTEGER:: n

    IF(nn<=0) RETURN

    indx(1)=0
    joc=1
    CALL invx(diag(1,1),ndeg,ir)
    DO i=2,nn
      indx(i)=joc
      DO j=1,i-1
        CALL dxdot(t,dsln(1,1,indx(i)),dsln(1,1,indx(j)),j-1)
        DO m=1,ndeg
          DO n=1,ndeg
            dsln(n,m,joc)=dsln(n,m,joc)-t(n,m)
          ENDDO
        ENDDO
        joc=joc+1
      ENDDO
      CALL vxprod(ndegl,dsln(1,1,indx(i)),diag,temp,i-1)
      CALL dxdotl(t,temp,dsln(1,1,indx(i)),i-1)
      locd=0
      DO n=1,ndeg
        DO m=1,n
          locd=locd+1
          diag(locd,i)=diag(locd,i)-t(n,m)
        ENDDO
      ENDDO
      CALL vcopy(temp,dsln(1,1,indx(i)),ndeg*ndeg*(i-1))
      CALL invx(diag(1,i),ndeg,ir)
    ENDDO
    RETURN
  END SUBROUTINE sxum3

  !======================================================================!
  !> v2prod
  !======================================================================!
  SUBROUTINE v2prod(a,b,c,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(4,n):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,n):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4,n):: c

    INTEGER:: i

    DO i=1,n
      c(3,i)=a(3,i)-a(1,i)*b(2,i)
      c(1,i)=a(1,i)*b(1,i)
      c(3,i)=c(3,i)*b(3,i)
      c(1,i)=c(1,i)-c(3,i)*b(2,i)

      c(4,i)=a(4,i)-a(2,i)*b(2,i)
      c(2,i)=a(2,i)*b(1,i)
      c(4,i)=c(4,i)*b(3,i)
      c(2,i)=c(2,i)-c(4,i)*b(2,i)
    ENDDO
    RETURN
  END SUBROUTINE v2prod

  !======================================================================!
  !> v3prod
  !======================================================================!
  SUBROUTINE v3prod(zln,diag,zz,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(6,n):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(9,n):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,n):: zz

    INTEGER:: i

    DO i=1,n
      zz(4,i)=zln(4,i)-zln(1,i)*diag(2,i)
      zz(7,i)=zln(7,i)-zln(1,i)*diag(4,i)-zz(4,i)*diag(5,i)
      zz(1,i)=zln(1,i)*diag(1,i)
      zz(4,i)=zz(4,i)*diag(3,i)
      zz(7,i)=zz(7,i)*diag(6,i)
      zz(4,i)=zz(4,i)-zz(7,i)*diag(5,i)
      zz(1,i)=zz(1,i)-zz(4,i)*diag(2,i)-zz(7,i)*diag(4,i)

      zz(5,i)=zln(5,i)-zln(2,i)*diag(2,i)
      zz(8,i)=zln(8,i)-zln(2,i)*diag(4,i)-zz(5,i)*diag(5,i)
      zz(2,i)=zln(2,i)*diag(1,i)
      zz(5,i)=zz(5,i)*diag(3,i)
      zz(8,i)=zz(8,i)*diag(6,i)
      zz(5,i)=zz(5,i)-zz(8,i)*diag(5,i)
      zz(2,i)=zz(2,i)-zz(5,i)*diag(2,i)-zz(8,i)*diag(4,i)

      zz(6,i)=zln(6,i)-zln(3,i)*diag(2,i)
      zz(9,i)=zln(9,i)-zln(3,i)*diag(4,i)-zz(6,i)*diag(5,i)
      zz(3,i)=zln(3,i)*diag(1,i)
      zz(6,i)=zz(6,i)*diag(3,i)
      zz(9,i)=zz(9,i)*diag(6,i)
      zz(6,i)=zz(6,i)-zz(9,i)*diag(5,i)
      zz(3,i)=zz(3,i)-zz(6,i)*diag(2,i)-zz(9,i)*diag(4,i)
    ENDDO
    RETURN
  END SUBROUTINE v3prod

  !======================================================================!
  !> vcopy
  !======================================================================!
  SUBROUTINE vcopy(a,c,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n):: a
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(n):: c

    c=a
    RETURN
  END SUBROUTINE vcopy

  !======================================================================!
  !> vlcpy
  !======================================================================!
  SUBROUTINE vlcpy(a,b,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n):: b
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,n):: a

    INTEGER:: i
    INTEGER:: j

    DO i = 1,n
      DO j = 1,n
        a(j,i) = b(i,j)
      END DO
    END DO
    RETURN
  END SUBROUTINE vlcpy

  !======================================================================!
  !> vprod
  !======================================================================!
  SUBROUTINE vprod(a,b,c,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n):: a
    DOUBLE PRECISION, INTENT(IN), DIMENSION(n):: b

    DOUBLE PRECISION, INTENT(OUT), DIMENSION(n):: c

    c(1:n)=a(1:n)*b(1:n)
    RETURN
  END SUBROUTINE vprod

  !======================================================================!
  !> vxprod
  !======================================================================!
  SUBROUTINE vxprod(ndegl,zln,diag,zz,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    INTEGER, INTENT(IN):: ndegl

    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndegl,n):: diag
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg*ndeg,n):: zln
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg*ndeg,n):: zz

    INTEGER:: i

    DO i=1,n
      CALL invxx(zz(1,i),zln(1,i),diag(1,i),ndeg)
    ENDDO
    RETURN
  END SUBROUTINE vxprod

  !======================================================================!
  !> zpivot
  !======================================================================!
  SUBROUTINE zpivot(neqns,neqnsz,nttbr,jcol,irow,zpiv,ir)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nttbr
    INTEGER, INTENT(IN), DIMENSION(*):: jcol
    INTEGER, INTENT(IN), DIMENSION(*):: irow

    INTEGER, INTENT(OUT):: ir
    INTEGER, INTENT(OUT):: neqnsz
    INTEGER, INTENT(OUT), DIMENSION(*):: zpiv

    INTEGER:: i
    INTEGER:: j
    INTEGER:: l
    INTEGER:: idbg
    COMMON /debug/idbg

    ir=0
    zpiv(1:neqns)=1

    DO l=1,nttbr
      i=irow(l)
      j=jcol(l)
      IF(i<=0.OR.j<=0) THEN
        ir=-1
        GOTO 1000
      ELSEIF(i>neqns.OR.j>neqns) THEN
        ir=1
        GOTO 1000
      ENDIF
      IF(i==j) zpiv(i)=0
    ENDDO
    DO i=neqns,1,-1
      IF(zpiv(i)==0) THEN
        neqnsz=i
        EXIT
      ENDIF
    ENDDO
1000 CONTINUE
    IF(idbg/=0) WRITE(6,"(20i3)") (zpiv(i),i=1,neqns)
    RETURN

  END SUBROUTINE zpivot

END MODULE hecmw_solver_direct
