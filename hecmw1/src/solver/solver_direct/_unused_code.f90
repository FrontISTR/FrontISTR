!======================================================================!
!> d6dot
!      spdot1 performs inner product of sparse vectors
!      #coded by t.arakawa of RIST on 040510
!======================================================================!
SUBROUTINE d6dot(t,a,b,n)
  IMPLICIT NONE

  INTEGER, INTENT(IN):: n
  DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: a
  DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: b

  DOUBLE PRECISION, INTENT(OUT), DIMENSION(9):: t

  INTEGER:: jj

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
END SUBROUTINE d6dot

!======================================================================!
!> d6dotl
!      spdot1 performs inner product of sparse vectors
!      #coded by t.arakawa of RIST on 040510
!======================================================================!
SUBROUTINE d6dotl(t,a,b,n)
  IMPLICIT NONE

  INTEGER, INTENT(IN):: n
  DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: a
  DOUBLE PRECISION, INTENT(IN), DIMENSION(9,*):: b

  DOUBLE PRECISION, INTENT(OUT), DIMENSION(6):: t

  INTEGER:: jj

  !$dir max_trips(6)
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
END SUBROUTINE d6dotl

!======================================================================!
!> d6sdot
!      spdot1 performs inner product of sparse vectors
!      #coded by t.arakawa of RIST on 040510
!======================================================================!
SUBROUTINE d6sdot(wi,a,b,n)
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
END SUBROUTINE d6sdot

!======================================================================!
!> idntty
!======================================================================!
SUBROUTINE idntty(neqns,invp,iperm)
  IMPLICIT NONE

  INTEGER, INTENT(IN):: neqns

  INTEGER, INTENT(OUT), DIMENSION(*):: invp
  INTEGER, INTENT(OUT), DIMENSION(*):: iperm

  INTEGER:: i
  INTEGER:: idbg1
  COMMON /debug/ idbg1

  i=1
  DO
    IF(i>neqns) THEN
      DO i=1,neqns
        iperm(invp(i))=i
      ENDDO
      RETURN
    ENDIF

    WRITE(6,*) 'invp(',i,')'
    READ(5,*) invp(i)

    IF(invp(i)==0) THEN
      DO i=1,neqns
        invp(i)=i
        iperm(i)=i
      ENDDO
      RETURN
    ENDIF

    IF(invp(i)<0) THEN
      READ(11,*) (invp(i),i=1,neqns)
      DO i=1,neqns
        iperm(invp(i))=i
      ENDDO
      RETURN
    ENDIF
    i=i+1
  ENDDO
  RETURN
END SUBROUTINE idntty

!======================================================================!
!> nusol6
!======================================================================!
SUBROUTINE nusol6(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns,nstop)
  IMPLICIT NONE

  INTEGER, INTENT(IN):: neqns
  INTEGER, INTENT(IN):: nstop
  INTEGER, INTENT(IN), DIMENSION(*):: xlnzr
  INTEGER, INTENT(IN), DIMENSION(*):: colno
  INTEGER, INTENT(IN), DIMENSION(*):: iperm
  !GP: DEBUG 13May04  wk(3 ---> wk(6, b(3 ---> b(6, diag(6 ---> diag(21,
  !GP: DEBUG 13May04  zln(9 ---> zln(36, dsln(9 ---> dsln(36
  !      double precision zln(9,*),diag(6,*),b(3,*),wk(3,*),dsln(9,*)
  DOUBLE PRECISION, INTENT(IN), DIMENSION(36,*):: zln
  DOUBLE PRECISION, INTENT(IN), DIMENSION(21,*):: diag
  DOUBLE PRECISION, INTENT(IN), DIMENSION(36,*):: dsln
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: b
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(6,*):: wk

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
    wk(4,i)=b(4,iperm(i))
    wk(5,i)=b(5,iperm(i))
    wk(6,i)=b(6,iperm(i))
  ENDDO
  joc=1
  DO i=1,neqns
    ks=xlnzr(i)
    ke=xlnzr(i+1)-1
    IF(ke<ks) GOTO 110
    CALL s6pdot(wk(1,i),wk,zln,colno,ks,ke)

110 CONTINUE
    IF(i<=nstop) EXIT
    CALL dxsdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
    joc=joc+i-nstop
  ENDDO
  DO i=1,neqns
    wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
    wk(3,i)=wk(3,i)-wk(1,i)*diag(4,i)-wk(2,i)*diag(5,i)
    wk(4,i)=wk(4,i)-wk(1,i)*diag(7,i)-wk(2,i)*diag(8,i)-wk(3,i)*diag(9,i)
    wk(5,i)=wk(5,i)-wk(1,i)*diag(11,i)-wk(2,i)*diag(12,i)-wk(3,i)*diag(13,i)-wk(4,i)*diag(14,i)
    wk(6,i)=wk(6,i)-wk(1,i)*diag(16,i)-wk(2,i)*diag(17,i)-wk(3,i)*diag(18,i)-wk(4,i)*diag(19,i)-wk(6,i)*diag(20,i)
    wk(1,i)=wk(1,i)*diag(1,i)
    wk(2,i)=wk(2,i)*diag(3,i)
    wk(3,i)=wk(3,i)*diag(6,i)
    wk(4,i)=wk(4,i)*diag(10,i)
    wk(5,i)=wk(5,i)*diag(15,i)
    wk(6,i)=wk(6,i)*diag(21,i)
    wk(5,i)=wk(5,i)-wk(6,i)*diag(20,i)
    wk(4,i)=wk(4,i)-wk(6,i)*diag(19,i)-wk(5,i)*diag(14,i)
    wk(3,i)=wk(3,i)-wk(6,i)*diag(18,i)-wk(5,i)*diag(13,i)-wk(4,i)*diag(9,i)
    wk(2,i)=wk(2,i)-wk(6,i)*diag(17,i)-wk(5,i)*diag(12,i)-wk(4,i)*diag(8,i)-wk(3,i)*diag(5,i)
    wk(1,i)=wk(1,i)-wk(6,i)*diag(16,i)-wk(5,i)*diag(11,i)-wk(4,i)*diag(7,i)-wk(3,i)*diag(4,i)-wk(2,i)*diag(2,i)
  ENDDO
  ! back ward
  DO i=neqns,1,-1
    IF(i<nstop) EXIT
    DO j=i-1,nstop,-1
      joc=joc-1
      wk(1,j)=wk(1,j)-wk(1,i)*dsln(1,joc)-wk(2,i)*dsln(2,joc) &
           -wk(3,i)*dsln(3,joc)-wk(4,i)*dsln(4,joc)-wk(5,i)*dsln(5,joc)-wk(6,i)*dsln(6,joc)
      wk(2,j)=wk(2,j)-wk(1,i)*dsln(7,joc)-wk(2,i)*dsln(8,joc) &
           -wk(3,i)*dsln(9,joc)-wk(4,i)*dsln(10,joc)-wk(5,i)*dsln(11,joc)-wk(6,i)*dsln(12,joc)
      wk(3,j)=wk(3,j)-wk(1,i)*dsln(13,joc)-wk(2,i)*dsln(14,joc) &
           -wk(3,i)*dsln(15,joc)-wk(4,i)*dsln(16,joc)-wk(5,i)*dsln(17,joc)-wk(6,i)*dsln(18,joc)
      wk(4,j)=wk(4,j)-wk(1,i)*dsln(19,joc)-wk(2,i)*dsln(20,joc) &
           -wk(3,i)*dsln(21,joc)-wk(4,i)*dsln(22,joc)-wk(5,i)*dsln(23,joc)-wk(6,i)*dsln(24,joc)
      wk(5,j)=wk(5,j)-wk(1,i)*dsln(25,joc)-wk(2,i)*dsln(26,joc) &
           -wk(3,i)*dsln(27,joc)-wk(4,i)*dsln(28,joc)-wk(5,i)*dsln(29,joc)-wk(6,i)*dsln(30,joc)
      wk(6,j)=wk(6,j)-wk(1,i)*dsln(31,joc)-wk(2,i)*dsln(32,joc) &
           -wk(3,i)*dsln(33,joc)-wk(4,i)*dsln(34,joc)-wk(5,i)*dsln(35,joc)-wk(6,i)*dsln(36,joc)
    ENDDO

    ks=xlnzr(i)
    ke=xlnzr(i+1)-1
    IF(ke<ks) EXIT
    DO k=ks,ke
      j=colno(k)
      wk(1,j)=wk(1,j)-wk(1,i)*zln(1,joc)-wk(2,i)*zln(2,joc) &
           -wk(3,i)*zln(3,joc)-wk(4,i)*zln(4,joc)-wk(5,i)*zln(5,joc)-wk(6,i)*zln(6,joc)
      wk(2,j)=wk(2,j)-wk(1,i)*zln(7,joc)-wk(2,i)*zln(8,joc) &
           -wk(3,i)*zln(9,joc)-wk(4,i)*zln(10,joc)-wk(5,i)*zln(11,joc)-wk(6,i)*zln(12,joc)
      wk(3,j)=wk(3,j)-wk(1,i)*zln(13,joc)-wk(2,i)*zln(14,joc) &
           -wk(3,i)*zln(15,joc)-wk(4,i)*zln(16,joc)-wk(5,i)*zln(17,joc)-wk(6,i)*zln(18,joc)
      wk(4,j)=wk(4,j)-wk(1,i)*zln(19,joc)-wk(2,i)*zln(20,joc) &
           -wk(3,i)*zln(21,joc)-wk(4,i)*zln(22,joc)-wk(5,i)*zln(23,joc)-wk(6,i)*zln(24,joc)
      wk(5,j)=wk(5,j)-wk(1,i)*zln(25,joc)-wk(2,i)*zln(26,joc) &
           -wk(3,i)*zln(27,joc)-wk(4,i)*zln(28,joc)-wk(5,i)*zln(29,joc)-wk(6,i)*zln(30,joc)
      wk(6,j)=wk(6,j)-wk(1,i)*zln(31,joc)-wk(2,i)*zln(32,joc) &
           -wk(3,i)*zln(33,joc)-wk(4,i)*zln(34,joc)-wk(5,i)*zln(35,joc)-wk(6,i)*zln(36,joc)
    ENDDO
  ENDDO
  ! permutaion
  DO i=1,neqns
    b(1,iperm(i))=wk(1,i)
    b(2,iperm(i))=wk(2,i)
    b(3,iperm(i))=wk(3,i)
    b(4,iperm(i))=wk(4,i)
    b(5,iperm(i))=wk(5,i)
    b(6,iperm(i))=wk(6,i)
  ENDDO
  RETURN
END SUBROUTINE nusol6

!======================================================================!
!> prt
!======================================================================!
SUBROUTINE prt(ip,n)
  IMPLICIT NONE

  INTEGER, INTENT(IN):: n
  INTEGER, INTENT(IN), DIMENSION(n):: ip

  INTEGER:: i

  WRITE(6,"(10(2x,i4))") (ip(i),i=1,n)
  RETURN
END SUBROUTINE prt

  !======================================================================!
  !> vlcpy1
  !======================================================================!
  SUBROUTINE vlcpy1(a,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: n
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(n):: a

    a(n)=0
    RETURN
  END SUBROUTINE vlcpy1

    !======================================================================!
  !> verif0
  !======================================================================!
  SUBROUTINE verif0(neqns,nttbr,irow,jcol,val,rhs,x)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: neqns
    INTEGER, INTENT(IN):: nttbr
    INTEGER, INTENT(IN), DIMENSION(*):: irow
    INTEGER, INTENT(IN), DIMENSION(*):: jcol
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,ndeg,*):: val
    DOUBLE PRECISION, INTENT(IN), DIMENSION(ndeg,*):: x
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(ndeg,*):: rhs

    INTEGER:: i
    INTEGER:: j
    INTEGER:: k
    INTEGER:: l
    INTEGER:: m
    DOUBLE PRECISION:: err
    DOUBLE PRECISION:: rel

    !----------------------------------------------------------------------
    !
    !     verify the solution(symmetric matrix)
    !
    !----------------------------------------------------------------------
    rel=0.0d0
    DO i=1,neqns
      DO l=1,ndeg
        rel=rel+dabs(rhs(l,i))
      ENDDO
    ENDDO
    DO k=1,nttbr
      i=irow(k)
      j=jcol(k)
      DO l=1,ndeg
        DO m=1,ndeg
          rhs(l,i)=rhs(l,i)-val(l,m,k)*x(m,j)
          IF(i/=j) rhs(l,j)=rhs(l,j)-val(m,l,k)*x(m,i)
        ENDDO
      ENDDO
    ENDDO
    err=0.0d0
    DO i=1,neqns
      DO l=1,ndeg
        err=err+dabs(rhs(l,i))
      ENDDO
    ENDDO
    WRITE(6,6000) err,rel,err/rel
    !WINDEBUG
    !      write(16,6000) err,rel,err/rel
    RETURN

6000 FORMAT(' ***verification***(symmetric)'/ &
         'norm(Ax-b)            =  ',1pd20.10/ &
         'norm(b)               =  ',1pd20.10/ &
         'norm(Ax-b)/norm(b)    =  ',1pd20.10)
  END SUBROUTINE verif0

  !======================================================================!
  !> v6prod
  !======================================================================!
  SUBROUTINE v6prod(zln,diag,zz,n)
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
  END SUBROUTINE v6prod

