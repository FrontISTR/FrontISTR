!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE D6DOT(T,A,B,N)
  IMPLICIT NONE

  INTEGER :: jj
  INTEGER :: l
  INTEGER :: N
  DOUBLE PRECISION :: T(9)
  DOUBLE PRECISION :: A(9,*)
  DOUBLE PRECISION :: B(9,*)
  !----------------------------------------------------------------------
  !
  !      spdot1 performs inner product of sparse vectors
  !
  !
  !      #coded by t.arakawa of RIST on 040510
  !
  !----------------------------------------------------------------------
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
END SUBROUTINE D6DOT
!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE D6DOTL(T,A,B,N)
  IMPLICIT NONE

  INTEGER :: jj
  INTEGER :: l
  INTEGER :: N
  DOUBLE PRECISION :: T(6)
  DOUBLE PRECISION :: A(9,*)
  DOUBLE PRECISION :: B(9,*)
  !----------------------------------------------------------------------
  !
  !      spdot1 performs inner product of sparse vectors
  !
  !
  !      #coded by t.arakawa of RIST on 040510
  !
  !----------------------------------------------------------------------
  !$dir max_trips(6)
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
END SUBROUTINE D6DOTL
!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE D6SDOT(Wi,A,B,N)
  IMPLICIT NONE

  INTEGER :: jj
  INTEGER :: N
  DOUBLE PRECISION :: Wi(3)
  DOUBLE PRECISION :: A(3,*)
  DOUBLE PRECISION :: B(9,*)
  !----------------------------------------------------------------------
  !
  !      spdot1 performs inner product of sparse vectors
  !
  !
  !      #coded by t.arakawa of RIST on 040510
  !
  !----------------------------------------------------------------------
  DO jj = 1, N
     Wi(1) = Wi(1) - A(1,jj)*B(1,jj) - A(2,jj)*B(4,jj) - A(3,jj)*B(7,jj)
     Wi(2) = Wi(2) - A(1,jj)*B(2,jj) - A(2,jj)*B(5,jj) - A(3,jj)*B(8,jj)
     Wi(3) = Wi(3) - A(1,jj)*B(3,jj) - A(2,jj)*B(6,jj) - A(3,jj)*B(9,jj)
  ENDDO
END SUBROUTINE D6SDOT
!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE IDNTTY(Neqns,Invp,Iperm)
  IMPLICIT NONE

  INTEGER :: i
  INTEGER :: IDBg1
  INTEGER :: Neqns
  INTEGER :: Invp(*)
  INTEGER :: Iperm(*)
  COMMON /DEBUG / IDBg1

  i = 1
  DO WHILE ( i<=Neqns )
     WRITE (6,*) 'invp(', i, ')'
     READ (5,*) Invp(i)
     IF ( Invp(i)==0 ) THEN
        DO i = 1, Neqns
           Invp(i) = i
           Iperm(i) = i
        ENDDO
        RETURN
     ELSEIF ( Invp(i)<0 ) THEN
        READ (11,*) (Invp(i),i=1,Neqns)
        DO i = 1, Neqns
           Iperm(Invp(i)) = i
        ENDDO
        GOTO 99999
     ELSE
        i = i + 1
     ENDIF
  ENDDO
  DO i = 1, Neqns
     Iperm(Invp(i)) = i
  ENDDO
  RETURN
99999 END SUBROUTINE IDNTTY
  !======================================================================!
  !                                                                      !
  !======================================================================!
SUBROUTINE NUSOL6(Xlnzr,Colno,Dsln,Zln,Diag,Iperm,B,Wk,Neqns,Nstop)
  IMPLICIT NONE

  INTEGER :: i
  INTEGER :: j
  INTEGER :: joc
  INTEGER :: k
  INTEGER :: ke
  INTEGER :: ks
  INTEGER :: Neqns
  INTEGER :: Nstop
  INTEGER :: Xlnzr(*)
  INTEGER :: Colno(*)
  INTEGER :: Iperm(*)
  !GP: DEBUG 13May04  wk(3 ---> wk(6, b(3 ---> b(6, diag(6 ---> diag(21,
  !GP: DEBUG 13May04  zln(9 ---> zln(36, dsln(9 ---> dsln(36
  !      double precision zln(9,*),diag(6,*),b(3,*),wk(3,*),dsln(9,*)
  DOUBLE PRECISION :: Zln(36,*)
  DOUBLE PRECISION :: Diag(21,*)
  DOUBLE PRECISION :: B(6,*)
  DOUBLE PRECISION :: Wk(6,*)
  DOUBLE PRECISION :: Dsln(36,*)
  ! forward
  DO i = 1, Neqns
     Wk(1,i) = B(1,Iperm(i))
     Wk(2,i) = B(2,Iperm(i))
     Wk(3,i) = B(3,Iperm(i))
     Wk(4,i) = B(4,Iperm(i))
     Wk(5,i) = B(5,Iperm(i))
     Wk(6,i) = B(6,Iperm(i))
  ENDDO
  joc = 1
  DO i = 1, Neqns
     ks = Xlnzr(i)
     ke = Xlnzr(i+1) - 1
     IF ( ke>=ks ) CALL S6PDOT(Wk(1,i),Wk,Zln,Colno,ks,ke)
     IF ( i>Nstop ) THEN
        !        call d6sdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
        CALL DXSDOT(6,Wk(1,i),Wk(1,Nstop),Dsln(1,joc),i-Nstop)
        joc = joc + i - Nstop
     ENDIF
  ENDDO
  DO i = 1, Neqns
     Wk(2,i) = Wk(2,i) - Wk(1,i)*Diag(2,i)
     Wk(3,i) = Wk(3,i) - Wk(1,i)*Diag(4,i) - Wk(2,i)*Diag(5,i)
     Wk(4,i) = Wk(4,i) - Wk(1,i)*Diag(7,i) - Wk(2,i)*Diag(8,i) - Wk(3,i)*Diag(9,i)
     Wk(5,i) = Wk(5,i) - Wk(1,i)*Diag(11,i) - Wk(2,i)*Diag(12,i) - Wk(3,i)*Diag(13,i) - Wk(4,i)*Diag(14,i)
     Wk(6,i) = Wk(6,i) - Wk(1,i)*Diag(16,i) - Wk(2,i)*Diag(17,i)&
          - Wk(3,i)*Diag(18,i) - Wk(4,i)*Diag(19,i) - Wk(6,i)&
          *Diag(20,i)
     Wk(1,i) = Wk(1,i)*Diag(1,i)
     Wk(2,i) = Wk(2,i)*Diag(3,i)
     Wk(3,i) = Wk(3,i)*Diag(6,i)
     Wk(4,i) = Wk(4,i)*Diag(10,i)
     Wk(5,i) = Wk(5,i)*Diag(15,i)
     Wk(6,i) = Wk(6,i)*Diag(21,i)
     Wk(5,i) = Wk(5,i) - Wk(6,i)*Diag(20,i)
     Wk(4,i) = Wk(4,i) - Wk(6,i)*Diag(19,i) - Wk(5,i)*Diag(14,i)
     Wk(3,i) = Wk(3,i) - Wk(6,i)*Diag(18,i) - Wk(5,i)*Diag(13,i) - Wk(4,i)*Diag(9,i)
     Wk(2,i) = Wk(2,i) - Wk(6,i)*Diag(17,i) - Wk(5,i)*Diag(12,i) - Wk(4,i)*Diag(8,i) - Wk(3,i)*Diag(5,i)
     Wk(1,i) = Wk(1,i) - Wk(6,i)*Diag(16,i) - Wk(5,i)*Diag(11,i)&
          - Wk(4,i)*Diag(7,i) - Wk(3,i)*Diag(4,i) - Wk(2,i)&
          *Diag(2,i)
  ENDDO
  ! back ward
  DO i = Neqns, 1, -1
     IF ( i>=Nstop ) THEN
        DO j = i - 1, Nstop, -1
           joc = joc - 1
           Wk(1,j) = Wk(1,j) - Wk(1,i)*Dsln(1,joc) - Wk(2,i)&
                *Dsln(2,joc) - Wk(3,i)*Dsln(3,joc) - Wk(4,i)&
                *Dsln(4,joc) - Wk(5,i)*Dsln(5,joc) - Wk(6,i)&
                *Dsln(6,joc)
           Wk(2,j) = Wk(2,j) - Wk(1,i)*Dsln(7,joc) - Wk(2,i)&
                *Dsln(8,joc) - Wk(3,i)*Dsln(9,joc) - Wk(4,i)&
                *Dsln(10,joc) - Wk(5,i)*Dsln(11,joc) - Wk(6,i)&
                *Dsln(12,joc)
           Wk(3,j) = Wk(3,j) - Wk(1,i)*Dsln(13,joc) - Wk(2,i)&
                *Dsln(14,joc) - Wk(3,i)*Dsln(15,joc) - Wk(4,i)&
                *Dsln(16,joc) - Wk(5,i)*Dsln(17,joc) - Wk(6,i)&
                *Dsln(18,joc)
           Wk(4,j) = Wk(4,j) - Wk(1,i)*Dsln(19,joc) - Wk(2,i)&
                *Dsln(20,joc) - Wk(3,i)*Dsln(21,joc) - Wk(4,i)&
                *Dsln(22,joc) - Wk(5,i)*Dsln(23,joc) - Wk(6,i)&
                *Dsln(24,joc)
           Wk(5,j) = Wk(5,j) - Wk(1,i)*Dsln(25,joc) - Wk(2,i)&
                *Dsln(26,joc) - Wk(3,i)*Dsln(27,joc) - Wk(4,i)&
                *Dsln(28,joc) - Wk(5,i)*Dsln(29,joc) - Wk(6,i)&
                *Dsln(30,joc)
           Wk(6,j) = Wk(6,j) - Wk(1,i)*Dsln(31,joc) - Wk(2,i)&
                *Dsln(32,joc) - Wk(3,i)*Dsln(33,joc) - Wk(4,i)&
                *Dsln(34,joc) - Wk(5,i)*Dsln(35,joc) - Wk(6,i)&
                *Dsln(36,joc)
        ENDDO
     ENDIF
     ks = Xlnzr(i)
     ke = Xlnzr(i+1) - 1
     IF ( ke>=ks ) THEN
        DO k = ks, ke
           j = Colno(k)
           Wk(1,j) = Wk(1,j) - Wk(1,i)*Zln(1,joc) - Wk(2,i)&
                *Zln(2,joc) - Wk(3,i)*Zln(3,joc) - Wk(4,i)&
                *Zln(4,joc) - Wk(5,i)*Zln(5,joc) - Wk(6,i)&
                *Zln(6,joc)
           Wk(2,j) = Wk(2,j) - Wk(1,i)*Zln(7,joc) - Wk(2,i)&
                *Zln(8,joc) - Wk(3,i)*Zln(9,joc) - Wk(4,i)&
                *Zln(10,joc) - Wk(5,i)*Zln(11,joc) - Wk(6,i)&
                *Zln(12,joc)
           Wk(3,j) = Wk(3,j) - Wk(1,i)*Zln(13,joc) - Wk(2,i)&
                *Zln(14,joc) - Wk(3,i)*Zln(15,joc) - Wk(4,i)&
                *Zln(16,joc) - Wk(5,i)*Zln(17,joc) - Wk(6,i)&
                *Zln(18,joc)
           Wk(4,j) = Wk(4,j) - Wk(1,i)*Zln(19,joc) - Wk(2,i)&
                *Zln(20,joc) - Wk(3,i)*Zln(21,joc) - Wk(4,i)&
                *Zln(22,joc) - Wk(5,i)*Zln(23,joc) - Wk(6,i)&
                *Zln(24,joc)
           Wk(5,j) = Wk(5,j) - Wk(1,i)*Zln(25,joc) - Wk(2,i)&
                *Zln(26,joc) - Wk(3,i)*Zln(27,joc) - Wk(4,i)&
                *Zln(28,joc) - Wk(5,i)*Zln(29,joc) - Wk(6,i)&
                *Zln(30,joc)
           Wk(6,j) = Wk(6,j) - Wk(1,i)*Zln(31,joc) - Wk(2,i)&
                *Zln(32,joc) - Wk(3,i)*Zln(33,joc) - Wk(4,i)&
                *Zln(34,joc) - Wk(5,i)*Zln(35,joc) - Wk(6,i)&
                *Zln(36,joc)
        ENDDO
     ENDIF
  ENDDO
  ! permutation
  DO i = 1, Neqns
     B(1,Iperm(i)) = Wk(1,i)
     B(2,Iperm(i)) = Wk(2,i)
     B(3,Iperm(i)) = Wk(3,i)
     B(4,Iperm(i)) = Wk(4,i)
     B(5,Iperm(i)) = Wk(5,i)
     B(6,Iperm(i)) = Wk(6,i)
  ENDDO
END SUBROUTINE NUSOL6
!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE PRT(Ip,N)
  IMPLICIT NONE

  INTEGER :: i
  INTEGER :: Ip
  INTEGER :: N
  DIMENSION Ip(N)
  WRITE (6,99001) (Ip(i),i=1,N)
99001 FORMAT (10(2x,i4))
END SUBROUTINE PRT

SUBROUTINE VLCPY1(A,N)
  IMPLICIT NONE

  DOUBLE PRECISION :: A
  INTEGER :: N
  DIMENSION A(N)
  INTEGER :: i
  INTEGER :: j
  A(N) = 0
END SUBROUTINE VLCPY1

!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE VERIF0(Neqns,Ndeg,Nttbr,Irow,Jcol,Val,Rhs,X)
  IMPLICIT NONE

  DOUBLE PRECISION :: err
  DOUBLE PRECISION :: rel
  DOUBLE PRECISION :: Rhs
  DOUBLE PRECISION :: Val
  DOUBLE PRECISION :: X
  INTEGER :: i
  INTEGER :: Irow
  INTEGER :: j
  INTEGER :: Jcol
  INTEGER :: k
  INTEGER :: l
  INTEGER :: m
  INTEGER :: Ndeg
  INTEGER :: Neqns
  INTEGER :: Nttbr
  DIMENSION Irow(*), Jcol(*), Val(Ndeg,Ndeg,*), Rhs(Ndeg,*), X(Ndeg,*)
  !----------------------------------------------------------------------
  !
  !     verify the solution(symmetric matrix)
  !
  !----------------------------------------------------------------------
  rel = 0.0D0
  DO i = 1, Neqns
     DO l = 1, Ndeg
        rel = rel + DABS(Rhs(l,i))
     ENDDO
  ENDDO
  DO k = 1, Nttbr
     i = Irow(k)
     j = Jcol(k)
     DO l = 1, Ndeg
        DO m = 1, Ndeg
           Rhs(l,i) = Rhs(l,i) - Val(l,m,k)*X(m,j)
           IF ( i/=j ) Rhs(l,j) = Rhs(l,j) - Val(m,l,k)*X(m,i)
        ENDDO
     ENDDO
  ENDDO
  err = 0.0D0
  DO i = 1, Neqns
     DO l = 1, Ndeg
        err = err + DABS(Rhs(l,i))
     ENDDO
  ENDDO
  WRITE (6,99001) err, rel, err/rel
  !WINDEBUG
  !      write(16,6000) err,rel,err/rel
99001 FORMAT (' ***verification***(symmetric)'/&
       'norm(Ax-b)            =  ',&
       1pd20.10/'norm(b)               =  ',&
       1pd20.10/'norm(Ax-b)/norm(b)    =  ',1pd20.10)
END SUBROUTINE VERIF0

!======================================================================!
!                                                                      !
!======================================================================!
SUBROUTINE V6PROD(Zln,Diag,Zz,N)
  IMPLICIT NONE

  DOUBLE PRECISION :: Diag
  DOUBLE PRECISION :: Zln
  DOUBLE PRECISION :: Zz
  INTEGER :: i
  INTEGER :: N
  DIMENSION Zln(9,N), Diag(6,N), Zz(9,N)
  DO i = 1, N
     Zz(4,i) = Zln(4,i) - Zln(1,i)*Diag(2,i)
     Zz(7,i) = Zln(7,i) - Zln(1,i)*Diag(4,i) - Zz(4,i)*Diag(5,i)
     Zz(1,i) = Zln(1,i)*Diag(1,i)
     Zz(4,i) = Zz(4,i)*Diag(3,i)
     Zz(7,i) = Zz(7,i)*Diag(6,i)
     Zz(4,i) = Zz(4,i) - Zz(7,i)*Diag(5,i)
     Zz(1,i) = Zz(1,i) - Zz(4,i)*Diag(2,i) - Zz(7,i)*Diag(4,i)
     !
     Zz(5,i) = Zln(5,i) - Zln(2,i)*Diag(2,i)
     Zz(8,i) = Zln(8,i) - Zln(2,i)*Diag(4,i) - Zz(5,i)*Diag(5,i)
     Zz(2,i) = Zln(2,i)*Diag(1,i)
     Zz(5,i) = Zz(5,i)*Diag(3,i)
     Zz(8,i) = Zz(8,i)*Diag(6,i)
     Zz(5,i) = Zz(5,i) - Zz(8,i)*Diag(5,i)
     Zz(2,i) = Zz(2,i) - Zz(5,i)*Diag(2,i) - Zz(8,i)*Diag(4,i)
     !
     Zz(6,i) = Zln(6,i) - Zln(3,i)*Diag(2,i)
     Zz(9,i) = Zln(9,i) - Zln(3,i)*Diag(4,i) - Zz(6,i)*Diag(5,i)
     Zz(3,i) = Zln(3,i)*Diag(1,i)
     Zz(6,i) = Zz(6,i)*Diag(3,i)
     Zz(9,i) = Zz(9,i)*Diag(6,i)
     Zz(6,i) = Zz(6,i) - Zz(9,i)*Diag(5,i)
     Zz(3,i) = Zz(3,i) - Zz(6,i)*Diag(2,i) - Zz(9,i)*Diag(4,i)
  ENDDO
END SUBROUTINE V6PROD
