!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takamichi Arakawa (Univ. of Tokyo)             !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
      module hecmw_solver_direct
        INTEGER (KIND=4),private :: len_colno
        INTEGER (KIND=4),private :: nstop
        INTEGER (KIND=4),private :: stage
        INTEGER (KIND=4),private :: neqns
        INTEGER (KIND=4),private :: nttbr
        INTEGER (KIND=4),private :: isym
        INTEGER (KIND=4),private :: ndeg
        INTEGER (KIND=4),private :: irr
        INTEGER (KIND=4),private :: len_dsln
        INTEGER (KIND=4),private :: len_iv
        INTEGER (KIND=4),private :: tot_dof
        INTEGER (KIND=4),private :: total
        INTEGER (KIND=4),private :: tmpmax
!        INTEGER (KIND=4),private :: maxdeg

        INTEGER (KIND=4),private, POINTER  :: jcol(:)
        INTEGER (KIND=4),private, POINTER  :: irow(:)
        INTEGER (KIND=4),private, POINTER  :: zpiv(:)
        INTEGER (KIND=4),private, POINTER  :: iperm(:)
        INTEGER (KIND=4),private, POINTER  :: invp(:)
        INTEGER (KIND=4),private, POINTER  :: parent(:)
        INTEGER (KIND=4),private, POINTER  :: nch(:)
        INTEGER (KIND=4),private, POINTER  :: xlnzr(:)
        INTEGER (KIND=4),private, POINTER  :: colno(:)

!*Work arrays
        INTEGER (KIND=4),private, POINTER  :: jcpt(:)
        INTEGER (KIND=4),private, POINTER  :: jcolno(:)
        INTEGER (KIND=4),private, POINTER  :: ia(:)
        INTEGER (KIND=4),private, POINTER  :: ja(:)
        INTEGER (KIND=4),private, POINTER  :: deg(:)
        INTEGER (KIND=4),private, POINTER  :: marker(:)
        INTEGER (KIND=4),private, POINTER  :: rchset(:)
        INTEGER (KIND=4),private, POINTER  :: nbrhd(:)
        INTEGER (KIND=4),private, POINTER  :: qsize(:)
        INTEGER (KIND=4),private, POINTER  :: qlink(:)
!WINDEBUG
!        INTEGER (KIND=4),private, POINTER  :: nofsub(:)
        INTEGER (KIND=4),private :: nofsub
        INTEGER (KIND=4),private, POINTER  :: adjncy(:)
        INTEGER (KIND=4),private, POINTER  :: btree(:)
!DEBUGWIN
        !INTEGER (KIND=4),private :: btree(100000)
        INTEGER (KIND=4),private, POINTER  :: pordr(:)
        INTEGER (KIND=4),private, POINTER  :: adjncp(:)
        INTEGER (KIND=4),private, POINTER  :: xleaf(:)
        INTEGER (KIND=4),private, POINTER  :: leaf(:)
        INTEGER (KIND=4),private, POINTER  :: indx(:)
!        INTEGER (KIND=4),private :: nndegs(:)
!        INTEGER (KIND=4),private :: istdig(:)
!        INTEGER (KIND=4),private :: istind(:)
!        INTEGER (KIND=4),private :: nd_dof(:)

        REAL(KIND=8),private, POINTER  :: v(:)
        REAL(KIND=8),private, POINTER  :: val(:)
        REAL(KIND=8),private, POINTER  :: rhs(:)
        REAL(KIND=8),private, POINTER  :: b(:)
        REAL(KIND=8),private, POINTER  :: temp(:)
!       REAL(KIND=8),private, POINTER  :: wk(:)
        REAL(KIND=8),private, POINTER  :: diag(:)
        REAL(KIND=8),private, POINTER  :: zln(:)
        REAL(KIND=8),private, POINTER  :: dsln(:)


!*Timing
        REAL(KIND=8),private  :: tom(10)
        REAL(KIND=8),private  :: iv(1)
!*Allocation variables
        INTEGER(KIND=4),private :: ialoc
        INTEGER(KIND=4),private :: raloc
        INTEGER(KIND=4),private :: ierror
      contains

      SUBROUTINE hecmw_solve_direct(hecMESH,hecMAT,ifmsg)
!----------------------------------------------------------------------
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
!----------------------------------------------------------------------
      USE hecmw_util
!     USE sp_direct_solver
      implicit double precision(a-h, o-z)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      INTEGER i98,i97,ir,ifmsg
      common /mchdpn/ rmax,rmin,epsm,lratio
      common /qaz/ iseed,ixxx
      INTEGER, POINTER :: memchk(:)
      rmax=8.988d+307
      rmin=4.941d-300
      epsm=2.220d-16
      lratio=2
      iseed=1
      ir = 0


      !WRITE(IFMSG,*) "Interface to ADS from HECMW..."
      call ptime(t1)

!*EHM HECMW June 7 2004

      i98 = hecMAT%Iarray(98) 
      IF(hecMAT%Iarray(98).EQ.1) THEN
!* Interface to symbolic factorization
       CALL setij(hecMESH,hecMAT)

!* Symbolic factorization
       call matini(ir)
       hecMAT%Iarray(98) = 0
      write(6,*)"symbolic fct done"
      ENDIF
      call ptime(t2)
      t3 = t2


      i97 = hecMAT%Iarray(97) 
      IF(hecMAT%Iarray(97).EQ.1) THEN
!* Interface to numeric factorization
       CALL nuform(hecMESH,hecMAT,ir)
       call ptime(t3)

!* Numeric factorization
       call nufct0(ir)
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
      call ptime(t4)



!* Finalize

!      tom(5)=tt2-tt1

!*  Errors 1
      IF(i98.NE.0.AND.i98.NE.1) THEN
       WRITE(IFMSG,*) 'ERROR in symb. fact. flag: Should be 1 or 0'
       STOP 'ERROR in symb. fact. flag: Should be 1 or 0'
      ENDIF
      IF(i97.NE.0.AND.i97.NE.1) THEN
       WRITE(IFMSG,*) 'ERROR in numer. fact. flag: Should be 1 or 0'
       STOP 'ERROR in numer. fact. flag: Should be 1 or 0'
      ENDIF
      IF(i98.EQ.1.AND.i97.EQ.0) THEN
       WRITE(IFMSG,*) 'WARNING: Numeric factorization not performed!'
       STOP 'WARNING: Numeric factorization not performed! Solve will not be performed'
       GO TO 100
      ENDIF
!*  Errors 2
      if(ir.ne.0) then
!WINDEBUG
         WRITE(IFMSG,*) 'ERROR in nufct0. ir = ',ir
         STOP
      endif

      tom(1)=t2-t1
      tom(2)=t3-t2
      tom(3)=t4-t3

!* Solve
!* Backsubstitute
         call nusol0(hecMAT%B,ir)
         call ptime(t5)
!* Errors 4
      if(ir.ne.0) then
!WINDEBUG
         write(IFMSG,*) 'error in nusol0. irr = ',ir
         STOP
      endif


 100  RETURN

      END subroutine hecmw_solve_direct

!======================================================================!
!                                                                      !
!======================================================================!
      subroutine addr0(isw,i,j,aij,invp,xlnzr,colno,diag,zln, &
                      dsln,nstop,ndeg2,ndeg2l,ir)
 
 
      implicit double precision(a-h,o-z)
      integer invp(*),xlnzr(*),colno(*)
      double precision zln(ndeg2,*),diag(ndeg2l,*),dsln(ndeg2,*), aij(ndeg2)
 
      data idbg/0/
      ir=0
      ii=invp(i)
      jj=invp(j)
      if(idbg.ne.0) write(6,*) ii,jj,aij
      if(ii.eq.jj) then
         if(ndeg2.eq.1) then
            if(isw.eq.0) then
               diag(1,ii)=aij(1)
            else
               diag(1,ii)=diag(1,ii)+aij(1)
            endif
         elseif(ndeg2.eq.4) then
            if(isw.eq.0) then
               diag(1,ii)=aij(1)
               diag(2,ii)=aij(2)
               diag(3,ii)=aij(4)
            else
               diag(1,ii)=diag(1,ii)+aij(1)
               diag(2,ii)=diag(2,ii)+aij(2)
               diag(3,ii)=diag(3,ii)+aij(4)
            endif
         endif
         goto 1000
      endif
      itrans=0
      if(jj.gt.ii) then
         k=jj
         jj=ii
         ii=k
         itrans=1
      endif
      if(jj.ge.nstop) then
         i0=ii-nstop
         j0=jj-nstop+1
         k=i0*(i0-1)/2+j0
         if(ndeg2.eq.1) then
            dsln(1,k)=aij(1)
            goto 1000
         elseif(ndeg2.eq.4) then
            if(itrans.eq.0) then
               do 3 l=1,ndeg2
                  dsln(l,k)=aij(l)
   3           continue
               goto 1000
            else
               dsln(1,k)=aij(1)
               dsln(2,k)=aij(3)
               dsln(3,k)=aij(2)
               dsln(4,k)=aij(4)
               goto 1000
            endif
         endif
      endif
      ks=xlnzr(ii)
      ke=xlnzr(ii+1)-1
      do 100 k=ks,ke
         if(colno(k).eq.jj) then
            if(isw.eq.0) then
               if(ndeg2.eq.1) then
                  zln(1,k)=aij(1)
               elseif(ndeg2.eq.4) then
                  if(itrans.eq.0) then
                     do 4 l=1,ndeg2
                        zln(l,k)=aij(l)
    4                continue
                  else
                     zln(1,k)=aij(1)
                     zln(2,k)=aij(3)
                     zln(3,k)=aij(2)
                     zln(4,k)=aij(4)
                  endif
               endif
            else
               if(ndeg2.eq.1) then
                  zln(1,k)=zln(1,k)+aij(1)
               elseif(ndeg2.eq.4) then
                  if(itrans.eq.0) then
                     do 5 l=1,ndeg2
                        zln(l,k)=zln(l,k)+aij(l)
    5                continue
                  else
                     zln(1,k)=zln(1,k)+aij(1)
                     zln(2,k)=zln(2,k)+aij(3)
                     zln(3,k)=zln(3,k)+aij(2)
                     zln(4,k)=zln(4,k)+aij(4)
                  endif
               endif
            endif
            goto 1000
         endif
  100 continue
      ir=20
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine addr3(isw,i,j,aij,invp,xlnzr,colno,diag,zln, &
                      dsln,nstop,ir)
!
!
      implicit double precision(a-h,o-z)
      integer invp(*),xlnzr(*),colno(*)
      double precision zln(9,*),diag(6,*),dsln(9,*),aij(9)
!
      data idbg,ndeg2,ndeg2l/0,9,6/
      ir=0
      ii=invp(i)
      jj=invp(j)
      if(idbg.ne.0) write(6,*) ii,jj,aij
      if(ii.eq.jj) then
         diag(1,ii)=aij(1)
         diag(2,ii)=aij(2)
         diag(3,ii)=aij(5)
         diag(4,ii)=aij(3)
         diag(5,ii)=aij(6)
         diag(6,ii)=aij(9)
         goto 1000
      endif
      itrans=0
      if(jj.gt.ii) then
         k=jj
         jj=ii
         ii=k
         itrans=1
      endif
      if(jj.ge.nstop) then
         i0=ii-nstop
         j0=jj-nstop+1
         k=i0*(i0-1)/2+j0
         if(itrans.eq.0) then
            do 110 l=1,ndeg2
               dsln(l,k)=aij(l)
  110       continue
            goto 1000
         else
            dsln(1,k)=aij(1)
            dsln(2,k)=aij(4)
            dsln(3,k)=aij(7)
            dsln(4,k)=aij(2)
            dsln(5,k)=aij(5)
            dsln(6,k)=aij(8)
            dsln(7,k)=aij(3)
            dsln(8,k)=aij(6)
            dsln(9,k)=aij(9)
            goto 1000
         endif
      endif
      ks=xlnzr(ii)
      ke=xlnzr(ii+1)-1
      do 100 k=ks,ke
         if(colno(k).eq.jj) then
            if(itrans.eq.0) then
               do 120 l=1,ndeg2
                  zln(l,k)=aij(l)
  120          continue
            else
               zln(1,k)=aij(1)
               zln(2,k)=aij(4)
               zln(3,k)=aij(7)
               zln(4,k)=aij(2)
               zln(5,k)=aij(5)
               zln(6,k)=aij(8)
               zln(7,k)=aij(3)
               zln(8,k)=aij(6)
               zln(9,k)=aij(9)
            endif
            goto 1000
         endif
  100 continue
      ir=20
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine addrx(isw,i,j,aij,invp,xlnzr,colno,diag,zln, &
                      dsln,nstop,ndeg,ndeg2,ndeg2l,ir)
!
!
      implicit double precision(a-h,o-z)
      integer invp(*),xlnzr(*),colno(*)
      double precision zln(ndeg,ndeg,*),diag(ndeg2l,*), &
                       dsln(ndeg,ndeg,*),aij(ndeg,ndeg)
!
      data idbg/0/
      ir=0
      ii=invp(i)
      jj=invp(j)
      if(idbg.ne.0) write(6,*) ii,jj,aij
      if(ii.eq.jj) then
         l=0
         do 100 n=1,ndeg
            do 110 m=1,n
               l=l+1
               diag(l,ii)=aij(n,m)
  110       continue
  100    continue
         goto 1000
      endif
      itrans=0
      if(jj.gt.ii) then
         k=jj
         jj=ii
         ii=k
         itrans=1
      endif
      if(jj.ge.nstop) then
         i0=ii-nstop
         j0=jj-nstop+1
         k=i0*(i0-1)/2+j0
         if(itrans.eq.0) then
            do 120 m=1,ndeg
               do 130 n=1,ndeg
                  dsln(n,m,k)=aij(n,m)
  130          continue
  120       continue
            goto 1000
         else
            do 140 m=1,ndeg
               do 150 n=1,ndeg
                  dsln(n,m,k)=aij(m,n)
  150          continue
  140       continue
            goto 1000
         endif
      endif
      ks=xlnzr(ii)
      ke=xlnzr(ii+1)-1
      do 200 k=ks,ke
         if(colno(k).eq.jj) then
            if(itrans.eq.0) then
               do 160 m=1,ndeg
                  do 170 n=1,ndeg
                     zln(n,m,k)=aij(n,m)
  170             continue
  160          continue
            else
               do 180 m=1,ndeg
                  do 190 n=1,ndeg
                     zln(n,m,k)=aij(m,n)
  190             continue
  180          continue
            endif
            goto 1000
         endif
  200 continue
      ir=20
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine bringu(zpiv,iperm,invp,parent,izz, neqns,irr)
      integer zpiv(*),iperm(*),invp(*),parent(*)
!
!----------------------------------------------------------------------
!
!      bringu brings up zero pivots from bottom of the elimination tree
!      to higher nodes
!
!      irr = 0     complete
!          = 1     impossible
!
!
!----------------------------------------------------------------------
!
      idbg=0
      irr=0
      ib0=invp(izz)
      ib=ib0
  100 continue
      if(ib.le.0) goto 1000
      ibp=parent(ib)
      izzp=iperm(ibp)
      if(zpiv(izzp).eq.0) goto 110
      ib=ibp
      goto 100
  110 continue
      invp(izz)=ibp
      invp(izzp)=ib0
      iperm(ibp)=izz
      iperm(ib0)=izzp
      if(idbg.ne.0) then
         do 200 i=1,neqns
            if(invp(iperm(i)).ne.i) goto 210
            if(iperm(invp(i)).ne.i) goto 210
  200    continue
         goto 220
  210    continue
         write(6,*) 'permutation error'
         stop
      endif
  220 continue
      return
 1000 continue
      irr=1
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d2dot(t,a,b,n)
      implicit double precision(a-h,o-z)
      double precision t(4),a(4,*),b(4,*)
!
!----------------------------------------------------------------------
!
!      spdot1 performs inner product of sparse vectors
!
!
!      #coded by t.arakawa of RIST on 040510
!
!----------------------------------------------------------------------
!
      t(1)=0.0d0
      t(2)=0.0d0
      t(3)=0.0d0
      t(4)=0.0d0
      do 100 jj=1,n
         t(1)=t(1)+a(1,jj)*b(1,jj)+a(3,jj)*b(3,jj)
         t(2)=t(2)+a(2,jj)*b(1,jj)+a(4,jj)*b(3,jj)
         t(3)=t(3)+a(1,jj)*b(2,jj)+a(3,jj)*b(4,jj)
         t(4)=t(4)+a(2,jj)*b(2,jj)+a(4,jj)*b(4,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d2sdot(wi,a,b,n)
      implicit double precision(a-h,o-z)
      double precision wi(2),a(2,*),b(4,*)
!
!----------------------------------------------------------------------
!
!      spdot1 performs inner product of sparse vectors
!
!
!      #coded by t.arakawa of RIST on 040510
!
!----------------------------------------------------------------------
!
      do 100 jj=1,n
         wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(3,jj)
         wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(4,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d3dot(t,a,b,n)
      implicit double precision(a-h,o-z)
      double precision t(9),a(9,*),b(9,*)
!
!----------------------------------------------------------------------
!
!      spdot1 performs inner product of sparse vectors
!
!
!      #coded by t.arakawa of RIST on 040510
!
!----------------------------------------------------------------------
!
!$dir max_trips(9)
      do 10 l=1,9
         t(l)=0.0d0
   10 continue
      do 100 jj=1,n
         t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
         t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
         t(3)=t(3)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
         t(4)=t(4)+a(1,jj)*b(2,jj)+a(4,jj)*b(5,jj)+a(7,jj)*b(8,jj)
         t(5)=t(5)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
         t(6)=t(6)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
         t(7)=t(7)+a(1,jj)*b(3,jj)+a(4,jj)*b(6,jj)+a(7,jj)*b(9,jj)
         t(8)=t(8)+a(2,jj)*b(3,jj)+a(5,jj)*b(6,jj)+a(8,jj)*b(9,jj)
         t(9)=t(9)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d3dotl(t,a,b,n)
      implicit double precision(a-h,o-z)
      double precision t(6),a(9,*),b(9,*)
!
!----------------------------------------------------------------------
!
!      spdot1 performs inner product of sparse vectors
!
!
!      #coded by t.arakawa of RIST on 040510
!
!----------------------------------------------------------------------
!
      do 10 l=1,6
         t(l)=0.0d0
   10 continue
      do 100 jj=1,n
         t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
         t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
         t(3)=t(3)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
         t(4)=t(4)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
         t(5)=t(5)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
         t(6)=t(6)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d3sdot(wi,a,b,n)
      implicit double precision(a-h,o-z)
      double precision wi(3),a(3,*),b(9,*)
!
!----------------------------------------------------------------------
!
!      spdot1 performs inner product of sparse vectors
!
!
!      #coded by t.arakawa of RIST on 040510
!
!----------------------------------------------------------------------
!
      do 100 jj=1,n
         wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(4,jj)-a(3,jj)*b(7,jj)
         wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(5,jj)-a(3,jj)*b(8,jj)
         wi(3)=wi(3)-a(1,jj)*b(3,jj)-a(2,jj)*b(6,jj)-a(3,jj)*b(9,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d6dot(t,a,b,n)
      implicit double precision(a-h,o-z)
      double precision t(9),a(9,*),b(9,*)
!
!----------------------------------------------------------------------
!
!      spdot1 performs inner product of sparse vectors
!
!
!      #coded by t.arakawa of RIST on 040510
!
!----------------------------------------------------------------------
!
      do 10 l=1,9
         t(l)=0.0d0
   10 continue
      do 100 jj=1,n
         t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
         t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
         t(3)=t(3)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
         t(4)=t(4)+a(1,jj)*b(2,jj)+a(4,jj)*b(5,jj)+a(7,jj)*b(8,jj)
         t(5)=t(5)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
         t(6)=t(6)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
         t(7)=t(7)+a(1,jj)*b(3,jj)+a(4,jj)*b(6,jj)+a(7,jj)*b(9,jj)
         t(8)=t(8)+a(2,jj)*b(3,jj)+a(5,jj)*b(6,jj)+a(8,jj)*b(9,jj)
         t(9)=t(9)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d6dotl(t,a,b,n)
      implicit double precision(a-h,o-z)
      double precision t(6),a(9,*),b(9,*)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
!c$dir max_trips(6)
      do 10 l=1,6
         t(l)=0.0d0
   10 continue
      do 100 jj=1,n
         t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
         t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
         t(3)=t(3)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
         t(4)=t(4)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
         t(5)=t(5)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
         t(6)=t(6)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine d6sdot(wi,a,b,n)
      implicit double precision(a-h,o-z)
      double precision wi(3),a(3,*),b(9,*)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      do 100 jj=1,n
         wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(4,jj)-a(3,jj)*b(7,jj)
         wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(5,jj)-a(3,jj)*b(8,jj)
         wi(3)=wi(3)-a(1,jj)*b(3,jj)-a(2,jj)*b(6,jj)-a(3,jj)*b(9,jj)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      double precision function ddot(a,b,n)
      double precision a(n),b(n),s
      s=0.0d0
      do 100 i=1,n
         s=s+a(i)*b(i)
  100 continue
      ddot=s
      return
      end function
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine dxdot(ndeg,t,a,b,l)
      implicit double precision(a-h,o-z)
      double precision t(ndeg,ndeg),a(ndeg,ndeg,*),b(ndeg,ndeg,*)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
         do 223 n=1,ndeg
         do 222 m=1,ndeg
         t(n,m)=0.0d0
         do 221 k=1,ndeg
         do 100 jj=1,l
            t(n,m)=t(n,m)+a(n,k,jj)*b(m,k,jj)
  100    continue
  221    continue
  222    continue
  223    continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine dxdotl(ndeg,t,a,b,l)
      implicit double precision(a-h,o-z)
      double precision t(ndeg,ndeg),a(ndeg,ndeg,*),b(ndeg,ndeg,*)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
         do 223 n=1,ndeg
         do 222 m=1,n
         t(n,m)=0.0d0
         do 221 k=1,ndeg
         do 100 jj=1,l
            t(n,m)=t(n,m)+a(n,k,jj)*b(m,k,jj)
  100    continue
  221    continue
  222    continue
  223    continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine dxsdot(ndeg,wi,a,b,n)
      implicit double precision(a-h,o-z)
      double precision wi(ndeg),a(ndeg,*),b(ndeg,ndeg,*)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      do 102 jj=1,n
         do 101 m=1,ndeg
         do 100 n=1,ndeg
            wi(n)=wi(n)-b(n,m,jj)*a(m,jj)
  100 continue
  101 continue
  102 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine forpar(neqns,parent,nch,nstop)
!WINDEBUG
      !integer parent(neqns),nch(neqns)
      integer parent(*),nch(*)
      do 100 i=1,neqns
         nch(i)=0
  100 continue
      nch(neqns+1)=0
      do 200 i=1,neqns
         ii=parent(i)
         nch(ii)=nch(ii)+1
  200 continue
      do 300 i=neqns,1,-1
         if(nch(i).ne.1) goto 310
  300 continue
  310 continue
      !write(6,*) 'denseform 1 no 0'
      !read(5,*) idens
      idens = 0
      if(idens.eq.1) then
         nstop=i
      else
         nstop=neqns+1
      endif
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!


      SUBROUTINE setij(hecMESH,hecMAT)

      USE hecmw_util 
!     USE sp_direct_solver
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      integer(kind=kint)  l,ir,numnp,ndof
      integer(kind=kint) iiS,iiE,kki,kkj,ndof2

      NUMNP = hecMAT%NP
      NDOF  = hecMESH%n_dof
      ntotal = NUMNP*NDOF

!*NUFACT variables
      neqns = NUMNP
      ndeg = NDOF
      nttbr = hecMAT%NP+hecMAT%NPL !+hecMAT%NPU if unsymmetric
      isym = 0

!*Allocations
!!!@      ALLOCATE(irow(nttbr),STAT=ierror)
!!!@      ALLOCATE(jcol(nttbr),STAT=ierror)
!!!@      IF(ierror/=0) STOP "Allocation error: irow/jcol"
      ALLOCATE(irow(nttbr),STAT=ierr)
      ALLOCATE(jcol(nttbr),STAT=ierr)
      IF(ierr/=0) STOP "Allocation error: irow/jcol"

      kk = 0
      ndof2 = NDOF*NDOF
      do j= 1, NUMNP
!*Diagonal
       kk = kk + 1
       irow(kk) = j
       jcol(kk) = j
!*Lower
       do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i= hecMAT%itemL(k)
        kk = kk + 1
        irow(kk) = j
        jcol(kk) = i
       enddo

!!*Upper if unsymmetric
!       do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
!        i= hecMAT%itemU(k)
!        kk = kk + 1
!        irow(kk) = j
!        jcol(kk) = i
!       enddo

      enddo

      RETURN
      END subroutine

      SUBROUTINE nuform(hecMESH,hecMAT,ir)

      USE hecmw_util 
!     USE sp_direct_solver
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      integer(kind=kint)  l,ir,numnp,ndof
      integer(kind=kint) iiS,iiE,kki,kkj,ndof2

      NUMNP = hecMAT%NP
      NDOF  = hecMESH%n_dof
      ntotal = NUMNP*NDOF

     
!*NUFACT variables
      neqns = NUMNP
      ndeg = NDOF
      nttbr = hecMAT%NP+hecMAT%NPL !+hecMAT%NPU if unsymmetric
      isym = 0

!*Allocations
!!!@      ALLOCATE(val(ndeg*ndeg),STAT=ierror)
!!!@      IF(ierror/=0) STOP "Allocation error:val"
      ALLOCATE(val(ndeg*ndeg),STAT=ierr)
      IF(ierr/=0) STOP "Allocation error:val"
      write(6,*) "nuform:stage = ",stage
!     if(stage.ge.30)then
!         DEALLOCATE(diag,stat = ierror1)
!         DEALLOCATE(zln,STAT = ierror2)
!         DEALLOCATE(dsln,STAT = ierror3)
!         raloc = raloc - neqns*ndeg*(ndeg+1)/2
!         raloc = raloc - len_colno*ndeg*ndeg
!         raloc = raloc - len_dsln*ndeg*ndeg
!     write(6,*) "nuform:ierror = ",ierror1,ierror2,ierror3
!     endif
      kk = 0
      ndof2 = NDOF*NDOF
      do j= 1, NUMNP
!*Diagonal
       kk = kk + 1
       CALL vlcpy(val,hecMAT%D(ndof2*(j-1)+1:ndof2*j),ndof)
       !call staij1(0,j,j,val,v,ndof,ir)
       call staij1(0,j,j,val,ir)

        do i = 1,NDOF
          IF(val((i-1)*ndof+i).LE.0) THEN
           WRITE(IDBG,*)'j,j,val:',j,i,val((i-1)*ndof+i)
!!!           PAUSE 'Error?'
          ENDIF
        end do

!*Lower
       do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i= hecMAT%itemL(k)
        kk = kk + 1
       CALL vlcpy(val,hecMAT%AL(ndof2*(k-1)+1:ndof2*k),ndof)

       !call staij1(0,j,i,val,v,ndof,ir)
       call staij1(0,j,i,val,ir)
       enddo
      enddo

      DEALLOCATE(val)
      RETURN
      END subroutine

      SUBROUTINE vlcpy1(a,n)
      IMPLICIT double precision (a-h,o-z)
      dimension a(n)
      INTEGER i,j
      a(n)=0
      return
      END subroutine
      SUBROUTINE vlcpy(a,b,n)
      IMPLICIT double precision (a-h,o-z)
      INTEGER i,j
      DIMENSION a(n,n), b(n,n)

      DO i = 1,n
       DO j = 1,n
        a(j,i) = b(i,j)
       END DO
      END DO
      RETURN
      END subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine genbtq(xadj,adjncy,invp,parent,btree, zpiv,izz,neqns)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer  xadj(*),adjncy(*),parent(*),btree(2,*),invp(*), zpiv(*)
      common /debug/ idbg1
      do 10 i=1,neqns+1
         btree(1,i)=0
         btree(2,i)=0
   10 continue
      do 100 i=1,neqns+1
         ip=parent(i)
         if(ip.le.0) goto 100
         ib=btree(1,ip)
         if(ib.eq.0) then
            btree(1,ip)=i
         else
  101       continue
            inext=btree(2,ib)
            if(inext.eq.0) then
               btree(2,ib)=i
            else
               ib=inext
               goto 101
            endif
          endif
  100 continue
!c
!c find zeropivot
!c
      do 200 i=1,neqns
         if(zpiv(i).ne.0) then
            if(btree(1,invp(i)).eq.0) then
               izz=i
               goto 210
            endif
         endif
  200 continue
      izz=0
  210 continue
      if(idbg1.ne.0) write(6,6010)
      if(idbg1.ne.0) write(6,6000) (i,btree(1,i),btree(2,i),i=1,neqns)
      if(idbg1.ne.0) write(6,6020) izz
 6000 format(i6,'(',2i6,')')
 6010 format(' binary tree')
 6020 format(' the first zero pivot is ',i4)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine genpaq(xadj,adjncy,invp,iperm,parent,neqns,ancstr)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer xadj(*),adjncy(*),parent(*),invp(*),iperm(*),ancstr(*)
      common /debug/ idbg1
!c
      do 100 i=1,neqns
         parent(i)=0
         ancstr(i)=0
         ip=iperm(i)
         do 110 k=xadj(ip),xadj(ip+1)-1
            l=invp(adjncy(k))
            if(l.ge.i) goto 110
  112       continue
            if(ancstr(l).eq.0) goto 111
            if(ancstr(l).eq.i) goto 110
               it=ancstr(l)
               ancstr(l)=i
               l=it
               goto 112
  111       continue
            ancstr(l)=i
            parent(l)=i
  110    continue
  100 continue
      do 200 i=1,neqns
         if(parent(i).eq.0) parent(i)=neqns+1
  200 continue
      parent(neqns+1)=0
      if(idbg1.ne.0) write(6,6010)
      if(idbg1.ne.0) write(6,6000) (i,parent(i),i=1,neqns)
 6000 format(2i6)
 6010 format(' parent')
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine genqmd(neqns,xadj,adj0,perm,invp,deg,&
                       marker,rchset,nbrhd,qsize,qlink,&
                       nofsub,adjncy)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer adjncy(*),perm(*),invp(*),deg(*),marker(*),&
             rchset(*),nbrhd(*),qsize(*),qlink(*),adj0(*)
      integer xadj(*),inode,ip,irch,j,mindeg,ndeg, &
             neqns,nhdsze,node,nofsub,np,num,nump1,&
             nxnode,rchsze,search,thresh
!c
      mindeg=neqns
      nofsub=0
      do 10 i=1,xadj(neqns+1)-1
         adjncy(i)=adj0(i)
   10 continue
      do 100 node=1,neqns
         perm(node)=node
         invp(node)=node
         marker(node)=0
         qsize(node)=1
         qlink(node)=0
         ndeg=xadj(node+1)-xadj(node)
         deg(node)=ndeg
         if(ndeg.lt.mindeg) mindeg=ndeg
  100 continue
!c
      num=0
  200 search=1
      thresh=mindeg
      mindeg=neqns
  300 nump1=num+1
      if(nump1.gt.search) search=nump1
      do 400 j=search,neqns
         node=perm(j)
         if(marker(node).lt.0) goto 400
         ndeg=deg(node)
         if(ndeg.le.thresh) goto 500
         if(ndeg.lt.mindeg) mindeg=ndeg
  400 continue
      goto 200
!c
  500 search=j
      nofsub=nofsub+deg(node)
      marker(node)=1
      call qmdrch(node,xadj,adjncy,deg,marker,&
                  rchsze,rchset,nhdsze,nbrhd)
      nxnode=node
  600 num=num+1
      np=invp(nxnode)
      ip=perm(num)
      perm(np)=ip
      invp(ip)=np
      perm(num)=nxnode
      invp(nxnode)=num
      deg(nxnode)=-1
      nxnode=qlink(nxnode)
      if(nxnode.gt.0) goto 600
      if(rchsze.le.0) goto 800
!c
      call qmdupd(xadj,adjncy,rchsze,rchset,deg,&
                 qsize,qlink,marker,&
                 rchset(rchsze+1),nbrhd(nhdsze+1))
      marker(node)=0
      do 700 irch=1,rchsze
         inode=rchset(irch)
         if(marker(inode).lt.0) goto 700
         marker(inode)=0
         ndeg=deg(inode)
         if(ndeg.lt.mindeg) mindeg=ndeg
         if(ndeg.gt.thresh) goto 700
         mindeg=thresh
         thresh=ndeg
         search=invp(inode)
  700 continue
      if(nhdsze.gt.0) call qmdot(node,xadj, &
                 adjncy,marker,rchsze,rchset,nbrhd)
  800 if(num.lt.neqns) goto 300
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine gnclno(parent,pordr,xleaf,leaf,xlnzr,colno, &
                       neqns,nstop,lncol,ir)
!c
!c
      integer parent(*),pordr(*),xleaf(*),leaf(*),xlnzr(*),colno(*)
      common /debug/ idbg1
!c
      nc=0
      ir=0
      l=1
      do 100 i=1,neqns
         xlnzr(i)=l
         ks=xleaf(i)
         ke=xleaf(i+1)-1
         if(ke.lt.ks) goto 100
         nxleaf=leaf(ks)
         do 110 k=ks,ke-1
            j=nxleaf
            nxleaf=leaf(k+1)
  105       continue
            if(j.ge.nxleaf) goto 110
            if(j.ge.nstop) then
               goto 100
            endif
            colno(l)=j
            l=l+1
            j=parent(j)
            goto 105
  110    continue
         j=leaf(ke)
  115    continue
         if(j.ge.nstop) goto 100
         if(j.ge.i.or.j.eq.0) goto 100
         colno(l)=j
         l=l+1
         j=parent(j)
         goto 115
  100 continue
      xlnzr(neqns+1)=l
      lncol=l-1
      if(idbg1.ne.0) write(6,6010)
!c     if(idbg1.ne.0) write(6,6000) (xlnzr(i),i=1,neqns+1)
      if(idbg1.ne.0) write(6,6020) lncol
      if(idbg1.ne.0) then
         do 200 k=1,neqns
         write(6,6100) k
         write(6,6000) (colno(i),i=xlnzr(k),xlnzr(k+1)-1)
  200    continue
      endif
 6000 format(10i4)
 6010 format(' xlnzr')
 6020 format(' colno (lncol =',i10,')')
 6100 format(/' row = ',i6)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine gnleaf(xadj,adjncy,invp,iperm,pordr,nch,&
                       adjncp,xleaf,leaf,neqns,lnleaf)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer xadj(*),adjncy(*),pordr(*),nch(*),&
             adjncp(*),xleaf(*),leaf(*),invp(*),iperm(*)
      common /debug/ idbg1
!c
      l=1
      ik=0
      istart=0
      do 100 i=1,neqns
         xleaf(i)=l
         ip=iperm(i)
         do 105 k=xadj(ip),xadj(ip+1)-1
            iq=invp(adjncy(k))
            if(iq.lt.i) then
               ik=ik+1
               adjncp(ik)=iq
            endif
  105    continue
         m=ik-istart
         if(m.eq.0) goto 131
         call qqsort(adjncp(istart+1),m)
         lc1=adjncp(istart+1)
         if(lc1.ge.i) goto 100
         leaf(l)=lc1
         l=l+1
         do 130 k=istart+2,ik
            lc=adjncp(k)
!           if(lc.ge.i) goto 125
            if(lc1.lt.lc-nch(lc)) then
               leaf(l)=lc
               l=l+1
            endif
! 125       continue
            lc1=lc
  130    continue
         ik=1
         istart=ik
  131    continue
  100 continue
      xleaf(neqns+1)=l
      lnleaf=l-1
      if(idbg1.ne.0) write(6,6020)
      if(idbg1.ne.0) write(6,6000) (xleaf(i),i=1,neqns+1)
      if(idbg1.ne.0) write(6,6010) lnleaf
      if(idbg1.ne.0) write(6,6000) (leaf(i),i=1,lnleaf)
      return
 6000 format(10i6)
 6010 format(' leaf (len = ',i6,')')
 6020 format(' xleaf')
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine idntty(neqns,invp,iperm)
      integer invp(*),iperm(*)
      common /debug/ idbg1
!c
      i=1
   10 continue
      if(i.gt.neqns) goto 200
      write(6,*) 'invp(',i,')'
      read(5,*) invp(i)
      if(invp(i).eq.0) goto 120
      if(invp(i).lt.0) goto 220
      i=i+1
      goto 10
  200 continue
      do 210 i=1,neqns
         iperm(invp(i))=i
  210 continue
      return
  120 continue
      do 100 i=1,neqns
         invp(i)=i
         iperm(i)=i
  100 continue
      return
  220 continue
      read(11,*) (invp(i),i=1,neqns)
      do 230 i=1,neqns
         iperm(invp(i))=i
  230 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine inv2(dsln,ir)
      implicit double precision (a-h,o-z)
      dimension dsln(3)
      common /mchdpn/ rmax,rmin,epsm,lratio
      ir=0
      if(dabs(dsln(1)).lt.rmin) then
         ir=10
         return
      endif
      dsln(1)=1.0d0/dsln(1)
      t=dsln(2)*dsln(1)
      dsln(3)=dsln(3)-t*dsln(2)
      dsln(2)=t
      if(dabs(dsln(3)).lt.rmin) then
         ir=10
         return
      endif
      dsln(3)=1.0d0/dsln(3)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine inv22(zln,zz,diag)
      implicit double precision(a-h,o-z)
      dimension zln(4),zz(4),diag(3)
      zln(3)=zz(3)-zz(1)*diag(2)
      zln(1)=zz(1)*diag(1)
      zln(3)=zln(3)*diag(3)
      zln(1)=zln(1)-zln(3)*diag(2)
!c
      zln(4)=zz(4)-zz(2)*diag(2)
      zln(2)=zz(2)*diag(1)
      zln(4)=zln(4)*diag(3)
      zln(2)=zln(2)-zln(4)*diag(2)
!c
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine inv3(dsln,ir)
      implicit double precision (a-h,o-z)
      dimension dsln(6),t(2)
      common /mchdpn/ rmax,rmin,epsm,lratio
      ir=0
      if(dabs(dsln(1)).lt.rmin) then
         goto 999
      endif
      dsln(1)=1.0d0/dsln(1)
      t(1)=dsln(2)*dsln(1)
      dsln(3)=dsln(3)-t(1)*dsln(2)
      dsln(2)=t(1)
      if(dabs(dsln(3)).lt.rmin) then
         goto 999
      endif
      dsln(3)=1.0d0/dsln(3)
      t(1)=dsln(4)*dsln(1)
      dsln(5)=dsln(5)-dsln(2)*dsln(4)
      t(2)=dsln(5)*dsln(3)
      dsln(6)=dsln(6)-t(1)*dsln(4)-t(2)*dsln(5)
      dsln(4)=t(1)
      dsln(5)=t(2)
      if(dabs(dsln(6)).lt.rmin) then
         goto 999
      endif
      dsln(6)=1.0d0/dsln(6)
      return
  999 continue
      dsln(1)=1.0d0
      dsln(2)=0.0d0
      dsln(3)=1.0d0
      dsln(4)=0.0d0
      dsln(5)=0.0d0
      dsln(6)=1.0d0
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine inv33(zln,zz,diag)
      implicit double precision(a-h,o-z)
      dimension zln(9),zz(9),diag(6)
      zln(4)=zz(4)-zz(1)*diag(2)
      zln(7)=zz(7)-zz(1)*diag(4)-zln(4)*diag(5)
      zln(1)=zz(1)*diag(1)
      zln(4)=zln(4)*diag(3)
      zln(7)=zln(7)*diag(6)
      zln(4)=zln(4)-zln(7)*diag(5)
      zln(1)=zln(1)-zln(4)*diag(2)-zln(7)*diag(4)
!c
      zln(5)=zz(5)-zz(2)*diag(2)
      zln(8)=zz(8)-zz(2)*diag(4)-zln(5)*diag(5)
      zln(2)=zz(2)*diag(1)
      zln(5)=zln(5)*diag(3)
      zln(8)=zln(8)*diag(6)
      zln(5)=zln(5)-zln(8)*diag(5)
      zln(2)=zln(2)-zln(5)*diag(2)-zln(8)*diag(4)
!c
      zln(6)=zz(6)-zz(3)*diag(2)
      zln(9)=zz(9)-zz(3)*diag(4)-zln(6)*diag(5)
      zln(3)=zz(3)*diag(1)
      zln(6)=zln(6)*diag(3)
      zln(9)=zln(9)*diag(6)
      zln(6)=zln(6)-zln(9)*diag(5)
      zln(3)=zln(3)-zln(6)*diag(2)-zln(9)*diag(4)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine inv6(dsln,ir)
      implicit double precision (a-h,o-z)
      dimension dsln(21),t(5)
      common /mchdpn/ rmax,rmin,epsm,lratio
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
      dsln(14)=dsln(14)-dsln(11)*dsln(7)-dsln(12)*dsln(8) &
                       -dsln(13)*dsln(9)
      t(1)=dsln(11)*dsln(1)
      t(2)=dsln(12)*dsln(3)
      t(3)=dsln(13)*dsln(6)
      t(4)=dsln(14)*dsln(10)
      dsln(15)=1.0d0 &
            /(dsln(15)-t(1)*dsln(11)-t(2)*dsln(12) &
             -t(3)*dsln(13)-t(4)*dsln(14))
      dsln(11)=t(1)
      dsln(12)=t(2)
      dsln(13)=t(3)
      dsln(14)=t(4)
      dsln(17)=dsln(17)-dsln(16)*dsln(2)
      dsln(18)=dsln(18)-dsln(16)*dsln(4)-dsln(17)*dsln(5)
      dsln(19)=dsln(19)-dsln(16)*dsln(7)-dsln(17)*dsln(8) &
                       -dsln(18)*dsln(9)
      dsln(20)=dsln(20)-dsln(16)*dsln(11)-dsln(17)*dsln(12) &
                       -dsln(18)*dsln(13)-dsln(19)*dsln(14)
      t(1)=dsln(16)*dsln(1)
      t(2)=dsln(17)*dsln(3)
      t(3)=dsln(18)*dsln(6)
      t(4)=dsln(19)*dsln(10)
      t(5)=dsln(20)*dsln(15)
      dsln(21)=1.0d0 &
              /(dsln(21)-t(1)*dsln(16)-t(2)*dsln(17) &
                       -t(3)*dsln(18)-t(4)*dsln(19)-t(5)*dsln(20))
      dsln(16)=t(1)
      dsln(17)=t(2)
      dsln(18)=t(3)
      dsln(19)=t(4)
      dsln(20)=t(5)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine inv66(zln,zz,diag)
      implicit double precision(a-h,o-z)
      dimension zln(36),zz(36),diag(21)
      do 100 i=0,5
      zln(i+7)=zz(i+7)-zz(i+1)*diag(2)
      zln(i+13)=zz(i+13)-zz(i+1)*diag(4)-zln(i+7)*diag(5)
      zln(i+19)=zz(i+19)-zz(i+1)*diag(7)-zln(i+7)*diag(8) &
                       -zln(i+13)*diag(9)
      zln(i+25)=zz(i+25)-zz(i+1)*diag(11)-zln(i+7)*diag(12) &
                       -zln(i+13)*diag(13)-zln(i+19)*diag(14)
      zln(i+31)=zz(i+31)-zz(i+1)*diag(16)-zln(i+7)*diag(17) &
                       -zln(i+13)*diag(18)-zln(i+19)*diag(19) &
                       -zln(i+25)*diag(20)
      zln(i+1)=zz(i+1)*diag(1)
      zln(i+7)=zln(i+7)*diag(3)
      zln(i+13)=zln(i+13)*diag(6)
      zln(i+19)=zln(i+19)*diag(10)
      zln(i+25)=zln(i+25)*diag(15)
      zln(i+31)=zln(i+31)*diag(21)
      zln(i+25)=zln(i+25)-zln(i+31)*diag(20)
      zln(i+19)=zln(i+19)-zln(i+31)*diag(19)-zln(i+25)*diag(14)
      zln(i+13)=zln(i+13)-zln(i+31)*diag(18)-zln(i+25)*diag(13)&
                       -zln(i+19)*diag(9)
      zln(i+7)=zln(i+7)-zln(i+31)*diag(17)-zln(i+25)*diag(12) &
                       -zln(i+19)*diag(8)-zln(i+13)*diag(5)
      zln(i+1)=zln(i+1)-zln(i+31)*diag(16)-zln(i+25)*diag(11) &
                       -zln(i+19)*diag(7)-zln(i+13)*diag(4) &
                       -zln(i+7)*diag(2)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine invx(dsln,ndeg,ir)
      implicit double precision (a-h,o-z)
      dimension dsln(*)
      common /mchdpn/ rmax,rmin,epsm,lratio
      ir=0
      l=1
      dsln(1)=1.0d0/dsln(1)
      do 100 i=2,ndeg
         ld=0
         l0=l
         do 110 j=1,i-1
            l=l+1
            do 120 k=1,j-1
               ld=ld+1
               dsln(l)=dsln(l)-dsln(l0+k)*dsln(ld)
  120       continue
            ld=ld+1
  110    continue
         t=0.0d0
         k0=0
         ll=0
         do 130 k=l-i+2,l
            ll=ll+1
            k0=k0+ll
            tem=dsln(k)*dsln(k0)
            t=t+tem*dsln(k)
            dsln(k)=tem
  130    continue
         l=l+1
         dsln(l)=dsln(l)-t
         dsln(l)=1.0d0/dsln(l)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine invxx(zln,zz,diag,ndeg)
      implicit double precision(a-h,o-z)
      dimension zln(ndeg,ndeg),zz(ndeg,ndeg),diag(*)
      zln=zz
      do 100 l=1,ndeg,2
         joc=0
         do 120 m=1,ndeg-1
            joc=joc+m
            loc1=joc+m
            do 121 n=m+1,ndeg
               zln(l,n)=zln(l,n)-zln(l,m)*diag(loc1)
               zln(l+1,n)=zln(l+1,n)-zln(l+1,m)*diag(loc1)
               loc1=loc1+n
  121       continue
  120    continue
         joc=0
         do 130 m=1,ndeg
            joc=joc+m
            zln(l,m)=zln(l,m)*diag(joc)
            zln(l+1,m)=zln(l+1,m)*diag(joc)
  130    continue
         do 140 n=ndeg,2,-1
         joc=joc-1
         do 141 m=n-1,1,-1
            zln(l,m)=zln(l,m)-zln(l,n)*diag(joc)
            zln(l+1,m)=zln(l+1,m)-zln(l+1,n)*diag(joc)
            joc=joc-1
  141    continue
  140    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      !subroutine matini(neqs,nttbr,irow,jcol,lenv,iv,ir)
      subroutine matini(ir)
!c
!c
!     USE sp_direct_solver
      implicit double precision(a-h,o-z)
!c
!c----------------------------------------------------------------------
!c
!c     matini initializes storage for sparse matrix solver.
!c     this routine is used for both symmetric and asymmetric matrices 
!c     and must be called once at the beginning
!c
!c    (i)
!c        neqns     number of unknowns
!c        nttbr     number of non0s, pattern of non-zero elements are 
!c                  given like following.
!c                  nonz(A)={(i,j);i=irow(l),j=jcol(l); 1<= l <= nttbr}
!c        irow
!c        jcol      to define non-zero pattern
!c        lenv      length of the array v (iv)
!c
!c    (o)
!c        iv        comunication array. v is the original name
!c        ir        return code
!c                              =0    normal
!c                              =-1   non positive index
!c                              =1    too big index
!c                              =10   insufficient storage
!c
!c        contents of iv
!c               pointers 1 &zpiv(1)  2 &iperm(1)  3 &invp(1)
!c                        4 &parent(1)5 &nch(1)    6 &xlnzr(1) 
!c                        7 &colno(1) 8 &diag(1)   9 &zln(1)
!c                       10 &dsln(1)
!c
!c               scalars 21 len(colno)  22 nstop     23 stage
!c                       24 neqns       25 len(iv) 26 len(dsln)
!c                       27 total 
!c
!c        stage   10  after initialization
!c                20  building up matrix
!c                30  after LU decomposition
!c                40  after solving
!c
!c         #coded by t.arakawa of RIST on 040329
!c
!c----------------------------------------------------------------------
!c
      !integer irow(*),jcol(*),iv(*),neqs
      common /debug/ idbg
      common /mchdpn/ rmax,rmin,epsm,lratio
!c     rmax=8.988d+307
!c     rmin=4.941d-324
!c     epsm=2.220d-16
!c     lratio=2

      idbg=0
!c
      ir=0
      lenv = 1000000000
      lenv2=lratio*lenv
      neqns1=neqns+2
      len_dsln = lenv2

!c
!rmiv
!      iv(1)=51
!      iv(24)=neqns
!      iv(25)=lenv2

      len_iv = lenv2

      iv1 = 51

!*Initialize allocation measure variables
      ialoc = 0
      raloc = 0
!c
!c  set z pivot
!c
      ALLOCATE(zpiv(neqns),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, zpiv: SUB. matini"
      !ialoc = ialoc + neqns
      call zpivot(neqns,neqnsz,nttbr,jcol,irow,zpiv,ir1)
      if(ir1.ne.0) then
         ir=ir1
         goto 1000
      endif

!c
!c  build jcpt,jcolno
!c
!rmiv
!      lcpt=iv(1)+neqns1
      lcpt=iv1+neqns1
      lcolno=lcpt+2*nttbr
      left=lcolno+2*nttbr
      last=lenv2

!rmem
!      if(left.gt.last) then
!         ir=10
!         goto 1000
!      endif
      ALLOCATE(jcpt(2*nttbr),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, jcpt: SUB. matini"
      ALLOCATE(jcolno(2*nttbr),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, jcolno: SUB. matini"
      !ialoc = ialoc + 4*nttbr
      call stsmat(neqns,nttbr,irow,jcol,jcpt,jcolno)

!c
!c  build ia,ja
!c
      lia=last-neqns1
      lja=lia-nttbr*2
      last=lja
!rmem
!      if(left.gt.last) then
!        ir=10
!         goto 1000
!      endif
      ALLOCATE(ia(neqns+1),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, ia: SUB. matini"
      !ialoc = ialoc + neqns+1
!WINDEBUG
      ALLOCATE(ja(2*nttbr),STAT=ierror)
      !ALLOCATE(ja(ndeg*ndeg*nttbr),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, ja: SUB. matini"
      !ialoc = ialoc + 2*nttbr
      call stiaja(neqns,ia,ja,jcpt,jcolno)

!*Deallocation of work array
      DEALLOCATE(jcpt)
      DEALLOCATE(jcolno)


!c
!c  get permutation vector iperm,invp
!c
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
!      ALLOCATE(nofsub(neqns+1),STAT=ierror)
!      IF(ierror/=0) STOP "ALLOCATION ERROR, nofsub: SUB. matini"
      ALLOCATE(adjncy(2*nttbr),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, adjncy: SUB. matini"
      !ialoc = ialoc + 8*neqns + 2*nttbr + 7
      call genqmd(neqnsz,ia,ja,iperm,invp, &
                  deg,marker,rchset,nbrhd,qsize, &
                  qlink,nofsub,adjncy)


!c     else
!c     call idntty(neqns,invp,iv(iv(2)))
!c     endif

!rmiv
!      iv(iv(2)+neqns)=0
!      iv(iv(3)+neqns)=0
!c
!c   build up the parent vector parent vector will be saved in
!c   work2 for a while
!c
   10 continue
      call genpaq(ia,ja,invp,iperm,marker, neqns,rchset)
!c
!c   build up the binary tree
!c
      lbtree=lwk3-2*neqns
      last=lbtree
!rmem
!      if(left.gt.last) then
!         ir=10
!         goto 1000
!      endif
!WINDEBUG
      !ALLOCATE(btree(2*neqns),STAT=ierror)
      ALLOCATE(btree(2*(neqns+1)),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, btree: SUB. matini"
      !ialoc = ialoc + 2*neqns
      call genbtq(ia,ja,invp,marker, btree,zpiv,izz,neqns)
!c
!c   rotate the binary tree to avoid a zero pivot
!c
      if(izz.eq.0) goto 20
!c     write(6,*) 'zero pivot at ',izz
      if(izz0.eq.0) izz0=izz
      if(izz0.ne.izz) goto 30
      lwk4=last-neqns1
      lwk5=lwk4-neqns1
      last=lwk5
      call rotate(ia,ja,invp,iperm,marker, &
                 btree,izz,neqns,nbrhd,qsize,irr)
      goto 10
   30 continue
      call bringu(zpiv,iperm,invp,marker,izz, neqns,irr)

      goto 10
!c
!c   post ordering
!c
   20 continue
      lpordr=last-neqns1
      last=lpordr
!rmem
!      if(left.gt.last) then
!         ir=10
!         goto 1000
!      endif
      ALLOCATE(parent(neqns),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, parent: SUB. matini.f"
!WINDEBUG
!      ALLOCATE(nch(neqns),STAT=ierror)
      ALLOCATE(nch(neqns+1),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, nch: SUB. matini.f"
      ALLOCATE(pordr(neqns+1),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, pordr: SUB. matini.f"
      !ialoc = ialoc + 3*neqns + 1
      call posord(parent,btree,invp,iperm,pordr, &
                  nch,neqns,deg,marker,rchset)
!c
!c   generate skelton graph
!c
      lleaf=last-nttbr
      lxleaf=lleaf-neqns1
      ladp=lxleaf-neqns1
      last=ladp
!rmem
!      if(left.gt.last) then
!         ir=10
!         goto 1000
!      endif
      ALLOCATE(adjncp(neqns+1),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, adjncp: SUB. matini.f"
      ALLOCATE(xleaf(neqns+1),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, xleaf: SUB. matini.f"
      ALLOCATE(leaf(nttbr),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, leaf: SUB. matini.f"
      !ialoc = ialoc + 2*neqns+nttbr+2
      call gnleaf(ia,ja,invp,iperm,pordr,&
                  nch,adjncp,xleaf,leaf,neqns,lnleaf)
      call forpar(neqns,parent,nch,nstop)


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

!c
!c   build up xlnzr,colno  (this is the symbolic fct.)
!c

!rmiv
!      iv(6)=left
!      iv(7)=iv(6)+neqns1
!      maxl=lxleaf-iv(7)
 
      maxl=lxleaf-(left+neqns1)
      ALLOCATE(xlnzr(neqns+1),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, xlnzr: SUB. matini.f"
      !ialoc = ialoc + neqns + 1
      call pre_gnclno(parent,pordr,xleaf,leaf, &
                 xlnzr,neqns,nstop,lncol,ir1)
      ALLOCATE(colno(lncol),STAT=ierror)
      IF(ierror/=0) STOP "ALLOCATION ERROR, colno: SUB. matini.f"
      !ialoc = ialoc + lncol
      call gnclno(parent,pordr,xleaf,leaf, &
                 xlnzr,colno,neqns,nstop,lncol,ir1)

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
!      if(ir1.ne.0) then
!         ir=10
!         goto 1000
!      endif
!rmiv
!      left=iv(7)+lncol
      left=(left+neqns1)+lncol
!rmiv
!      iv(21)=lncol
!      iv(22)=nstop
!      iv(26)=(neqns-nstop+1)*(neqns-nstop)/2
      len_dsln =(neqns-nstop+1)*(neqns-nstop)/2

!Scalar assignments
      len_colno = lncol
!c
!c   area for double precision values
!c
      if(mod(left,2).eq.0) left=left+1

!rmiv
!      iv(8)=left
      total = left
!c     iv(9)=iv(8)+lratio*neqns
!c     iv(10)=iv(9)+lratio*lncol
!c     iv(27)=iv(10)+lratio*iv(26)

!rmiv
!      iv(23)=10
      stage = 10      
      ialoc = 5*neqns+lncol+1
 1000 continue

      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nufct(xlnzr,colno,dsln,zln,diag,indx,temp,neqns, &
                       parent,nch,nstop,ir)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),indx(*),parent(*),nch(*)
      double precision zln(*),diag(*),temp(*),dsln(*)
      common /mchdpn/ rmax,rmin,epsm,lratio
      common isem
!c
!c----------------------------------------------------------------------
!c
!c     nufct performs cholesky factorization in row order
!c
!c     (i) xlnzr,colno,zln,diag
!c         symbolicaly factorized
!c
!c     (o) zln,diag,dsln
!c
!c         #coded by t.arakawa of RIST on 040329
!c
!c----------------------------------------------------------------------
!c
!c
      isem=1
!c
!c phase I
!c
      call ptime(t1)
      diag(1)=1.0d0/diag(1)
      l=parent(1)
      nch(l)=nch(l)-1
      nch(1)=-1
      do 100 ic=2,nstop-1
         call sum(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
  100 continue
!c
!c phase II
!c
      call ptime(t2)
      do 200 ic=nstop,neqns
         call sum1(ic,xlnzr,colno,zln,diag,parent,temp,indx)
  200 continue
!c
!c phase III
!c
      call ptime(t3)
      call sum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
!c
!c phase IV
!c
      call ptime(t4)
      call sum3(neqns-nstop+1,dsln,diag(nstop),indx,temp)
      call ptime(t5)
      tt=t5-t1
      t1=t2-t1
      t2=t3-t2
      t3=t4-t3
      t4=t5-t4
      return
      ir=30
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nufct0(ir)
!c
!c
!     USE sp_direct_solver
      use hecmw_util
      implicit double precision(a-h,o-z)
      !dimension iv(*)
      ! caution) under def. may cause stack overflow.
      !dimension temp(ndeg,ndeg,neqns),indx(neqns)
!c
!c----------------------------------------------------------------------
!c
!c     this performs Cholesky factorization
!c
!c          if(iv(22).eq.0)    normal type
!c          if(iv(22).gt.0)    code generation type
!c
!c----------------------------------------------------------------------
!c
      common /mchdpn/ rmax,rmin,epsm,lratio
      common /qaz/iseed,ixxx
!c
      if(stage.ne.20) then
         Print *,'*********Setting Stage 40!*********'
         ir=40
         goto 1000
      else
         ir=0
      endif
!c
      allocate( temp(ndeg*ndeg*neqns), stat=irr )
      if( irr /= 0 ) then
            write(*,*) '##Error : Not enough memory'
            call hecmw_abort( hecmw_comm_get_comm())
            !stop
      endif
      allocate( indx(neqns), stat=irr )
      if( irr /= 0 ) then
            write(*,*) '##Error : Not enough memory'
            call hecmw_abort( hecmw_comm_get_comm())
            !stop
      endif
!c
!rmiv
!      ndeg=iv(28)
      ndegl=ndeg*(ndeg+1)
      ndegl=ndegl/2
      ndeg2=ndeg*ndeg
!rmiv
!      ndeg2=iv(28)*iv(28)
!     lndx=total
!      ltemp=lndx+neqns
!      if(mod(ltemp,2).eq.0) ltemp=ltemp+1
!      last=ltemp+lratio*neqns*(ndeg2+1)
!      if(last.gt.len_iv) then
!         write(6,*) "last,len_iv",last,len_iv
!         ir=10
!         goto 1000
!      endif
!     ALLOCATE(temp(ndeg2*neqns),STAT=ierror)
!     ALLOCATE(indx(neqns),STAT=ierror)
!     IF(ierror/=0) STOP "ALLOCATION ERROR, temp: SUB.nufct0.f"
      if(ndeg.eq.1) then
      call nufct(xlnzr,colno,dsln,zln,diag,  &
                 indx,temp,neqns,parent,nch, &
                 nstop,ir)
      elseif(ndeg.eq.2) then
      call nufct2(xlnzr,colno,dsln,zln,diag, &
                 indx,temp,neqns,parent,nch, &
                 nstop,ir)
      elseif(ndeg.eq.3) then
      call nufct3(xlnzr,colno,dsln,zln,diag, &
                 indx,temp,neqns,parent,nch, &
                 nstop,ir)
      elseif(ndeg.eq.6) then
      call nufct6(xlnzr,colno,dsln,zln,diag, &
                 indx,temp,neqns,parent,nch, &
                 nstop,ir)
      else
      call nufctx(xlnzr,colno,dsln,zln,diag, &
                 indx,temp,neqns,parent,nch, &
                 nstop,ndeg,ndegl,ir)
      endif
      stage=30
      deallocate( temp )
      deallocate( indx )
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nufct2(xlnzr,colno,dsln,zln,diag,indx,temp,neqns, &
                       parent,nch,nstop,ir)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),indx(*),parent(*),nch(*)
      double precision zln(4,*),diag(3,*),temp(4,*),dsln(4,*)
      common /mchdpn/ rmax,rmin,epsm,lratio
!c
!c----------------------------------------------------------------------
!c
!c     nufct performs cholesky factorization in row order
!c
!c     (i) xlnzr,colno,zln,diag
!c         symbolicaly factorized
!c
!c     (o) zln,diag,dsln
!c
!c         #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
!c
!c
!c phase I
!c
      ir=0
      call ptime(t1)
      if(nstop.gt.1) call inv2(diag(1,1),ir)
      l=parent(1)
      nch(l)=nch(l)-1
      nch(1)=-1
      do 100 ic=2,nstop-1
         call s2um(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
  100 continue
!c
!c phase II
!c
      call ptime(t2)
      do 200 ic=nstop,neqns
         call s2um1(ic,xlnzr,colno,zln,diag,parent,temp,indx)
  200 continue
!c
!c phase III
!c
      call ptime(t3)
      call s2um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
!c
!c phase IV
!c
      call ptime(t4)
      call s2um3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp)
      call ptime(t5)
      tt=t5-t1
      t1=t2-t1
      t2=t3-t2
      t3=t4-t3
      t4=t5-t4
!      ttt=tt5-tt1
!      tt1=tt2-tt1
!      tt2=tt3-tt2
!      tt3=tt4-tt3
!      tt4=tt5-tt4
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nufct3(xlnzr,colno,dsln,zln,diag,indx,temp,neqns, &
                       parent,nch,nstop,ir)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),indx(*),parent(*),nch(*)
      double precision zln(9,*),diag(6,*),temp(*),dsln(9,*)
      common /mchdpn/ rmax,rmin,epsm,lratio
!c
!c----------------------------------------------------------------------
!c
!c     nufct performs cholesky factorization in row order
!c
!c     (i) xlnzr,colno,zln,diag
!c         symbolicaly factorized
!c
!c     (o) zln,diag,dsln
!c
!c         #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
!c
!c
!c phase I
!c
      call ptime(t1)
      if(nstop.gt.1) call inv3(diag(1,1),ir)
      l=parent(1)
      nch(l)=nch(l)-1
      nch(1)=-1
      do 100 ic=2,nstop-1
         call s3um(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
  100 continue
!c
!c phase II
!c
      call ptime(t2)
      do 200 ic=nstop,neqns
         call s3um1(ic,xlnzr,colno,zln,diag,parent,temp,indx)
  200 continue
!c
!c phase III
!c
      call ptime(t3)
      call s3um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
!c
!c phase IV
!c
      call ptime(t4)
      call s3um3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp)
      call ptime(t5)
      tt=t5-t1
      t1=t2-t1
      t2=t3-t2
      t3=t4-t3
      t4=t5-t4
!      ttt=tt5-tt1
!      tt1=tt2-tt1
!      tt2=tt3-tt2
!      tt3=tt4-tt3
!      tt4=tt5-tt4
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nufct6(xlnzr,colno,dsln,zln,diag,indx,temp,neqns, &
                       parent,nch,nstop,ir)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),indx(*),parent(*),nch(*)
      double precision zln(36,*),diag(21,*),temp(*),dsln(36,*)
      common /mchdpn/ rmax,rmin,epsm,lratio
!c
!c----------------------------------------------------------------------
!c
!c     nufct performs cholesky factorization in row order
!c
!c     (i) xlnzr,colno,zln,diag
!c         symbolicaly factorized
!c
!c     (o) zln,diag,dsln
!c
!c         #coded by t.arakawa of RIST 040331
!c
!c----------------------------------------------------------------------
!c
!c
!c
!c phase I
!c
      call ptime(t1)
      if(nstop.gt.1) call inv6(diag(1,1),ir)
      l=parent(1)
      nch(l)=nch(l)-1
      nch(1)=-1
      do 100 ic=2,nstop-1
         call s6um(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx)
  100 continue
!c
!c phase II
!c
      call ptime(t2)
      do 200 ic=nstop,neqns
         call s6um1(ic,xlnzr,colno,zln,diag,parent,temp,indx)
  200 continue
!c
!c phase III
!c
      call ptime(t3)
      call s6um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
!c
!c phase IV
!c
      call ptime(t4)
      call s6um3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp)
      call ptime(t5)
      tt=t5-t1
      t1=t2-t1
      t2=t3-t2
      t3=t4-t3
      t4=t5-t4
!      ttt=tt5-tt1
!      tt1=tt2-tt1
!      tt2=tt3-tt2
!      tt3=tt4-tt3
!      tt4=tt5-tt4
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nufctx(xlnzr,colno,dsln,zln,diag,indx,temp,neqns, &
                       parent,nch,nstop,ndeg,ndegl,ir)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),indx(*),parent(*),nch(*)
      double precision zln(ndeg*ndeg,*),diag(ndegl,*), &
                       temp(ndeg*ndeg,*),dsln(ndeg*ndeg,*), &
                       zz(100),t(100)
      common /mchdpn/ rmax,rmin,epsm,lratio
!c
!c----------------------------------------------------------------------
!c
!c     nufct performs cholesky factorization in row order
!c
!c     (i) xlnzr,colno,zln,diag
!c         symbolicaly factorized
!c
!c     (o) zln,diag,dsln
!c
!c         #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
!c
!c
!c phase I
!c
      call ptime(t1)
      if(nstop.gt.1) call invx(diag(1,1),ndeg,ir)
      l=parent(1)
      nch(l)=nch(l)-1
      nch(1)=-1
      do 100 ic=2,nstop-1
         call sxum(ic,xlnzr,colno,zln,diag,nch,parent,temp,indx,ndeg, ndegl,zz,t) 
  100 continue
!c
!c phase II
!c
      call ptime(t2)
      do 200 ic=nstop,neqns
         call sxum1(ic,xlnzr,colno,zln,diag,parent,temp,indx, ndeg,ndegl,t)
  200 continue
!c
!c phase III
!c
      call ptime(t3)
      call sxum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx,ndeg,ndegl)
!c
!c phase IV
!c
      call ptime(t4)
      call sxum3(neqns-nstop+1,dsln,diag(1,nstop),indx,temp,ndeg,ndegl,t)
      call ptime(t5)
      tt=t5-t1
      t1=t2-t1
      t2=t3-t2
      t3=t4-t3
      t4=t5-t4
!      ttt=tt5-tt1
!      tt1=tt2-tt1
!      tt2=tt3-tt2
!      tt3=tt4-tt3
!      tt4=tt5-tt4
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      !subroutine nusol0(b,iv,ir) 
      subroutine nusol0(r_h_s,ir) 
!c
!c
!     USE sp_direct_solver
      use hecmw_util
      implicit double precision(a-h,o-z) 
      dimension r_h_s(*)
      ! dimention wk(ndeg*neqns) -> leading stack overflow
      double precision,pointer :: wk(:)
      INTEGER ir
      !dimension b(*),iv(*) 
!c 
!c---------------------------------------------------------------------- 
!c 
!c     this performs forward elimination and backward substitution 
!c 
!c     (i/o)
!c           r_h_s        on entry     right hand side vector
!c                    on exit      solution vector
!c           iv       communication array
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      common /mchdpn/ rmax,rmin,epsm,lratio
      common /qaz/iseed,ixxx
      if(stage.ne.30.and.stage.ne.40) then
         ir=50
         goto 1000
      else
         ir=0
      endif
      lwk=total

      allocate( wk(ndeg*neqns), stat=ierror )
      if( ierror /= 0 ) then
           write(*,*) "##Error: not enough memory"
           call hecmw_abort( hecmw_comm_get_comm() )
      endif
!rmiv
!      ndeg=iv(28)
      ndegl=ndeg*(ndeg+1)
      ndegl=ndegl/2
      if(ndeg.eq.1) then
      call nusol1(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop)
      elseif(ndeg.eq.2) then
      call nusol2(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop)
      elseif(ndeg.eq.3) then
      call nusol3(xlnzr,colno,dsln,zln,diag,iperm,r_h_s,wk,neqns,nstop)
      !elseif(ndeg.eq.6.and.ixxx.eq.1) then
      elseif(ndeg.eq.6) then
      call nusolx(xlnzr,colno,dsln,zln,diag, iperm,r_h_s,wk,neqns,nstop,ndeg,ndegl)
!      call nusol6(xlnzr,colno,dsln,zln,diag,
!     *                  iv(iv(2)),r_h_s,wk,neqns,nstop)
!      call nusol6(xlnzr,colno,dsln,zln,diag,
!     *                  iperm,r_h_s,wk,neqns,nstop)
      else
      call nusolx(xlnzr,colno,dsln,zln,diag, &
                       iperm,r_h_s,wk,neqns,nstop,ndeg,ndegl)
      endif
      stage=40
      deallocate( wk )
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nusol1(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns,nstop)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),iperm(*)
      double precision zln(*),diag(*),b(*),wk(*),dsln(*)
!c forward
      do 10 i=1,neqns
         wk(i)=b(iperm(i))
   10 continue
      joc=1
      do 100 i=1,neqns
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 110
         wk(i)=wk(i)-spdot2(wk,zln,colno,ks,ke)
  110    continue
         if(i.le.nstop) goto 100
         wk(i)=wk(i)-ddot(wk(nstop),dsln(joc),i-nstop)
         joc=joc+i-nstop
  100 continue
      do 120 i=1,neqns
         wk(i)=wk(i)*diag(i)
  120 continue
!c back ward
      do 200 i=neqns,1,-1
         if(i.lt.nstop) goto 206
         do 205 j=i-1,nstop,-1
            joc=joc-1
            wk(j)=wk(j)-wk(i)*dsln(joc)
  205    continue
  206    continue
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 200
         do 210 k=ks,ke
            j=colno(k)
            wk(j)=wk(j)-wk(i)*zln(k)
  210    continue
  200 continue
!c permutaion
      do 300 i=1,neqns
         b(iperm(i))=wk(i)
  300 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nusol2(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns, &
                     nstop)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),iperm(*)
      double precision zln(4,*),diag(3,*),b(2,*),wk(2,*),dsln(4,*)
!c forward
      do 10 i=1,neqns
         wk(1,i)=b(1,iperm(i))
         wk(2,i)=b(2,iperm(i))
   10 continue
      joc=1
      do 100 i=1,neqns
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 110
         call s2pdot(wk(1,i),wk,zln,colno,ks,ke)
  110    continue
         if(i.le.nstop) goto 100
         call d2sdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
         joc=joc+i-nstop
  100 continue
      do 120 i=1,neqns
         wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
         wk(1,i)=wk(1,i)*diag(1,i)
         wk(2,i)=wk(2,i)*diag(3,i)
         wk(1,i)=wk(1,i)-wk(2,i)*diag(2,i)
  120 continue
!c back ward
      do 200 i=neqns,1,-1
         if(i.lt.nstop) goto 206
         do 205 j=i-1,nstop,-1
            joc=joc-1
            wk(1,j)=wk(1,j)-wk(1,i)*dsln(1,joc)-wk(2,i)*dsln(2,joc)
            wk(2,j)=wk(2,j)-wk(1,i)*dsln(3,joc)-wk(2,i)*dsln(4,joc)
  205    continue
  206    continue
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 200
         do 210 k=ks,ke
            j=colno(k)
            wk(1,j)=wk(1,j)-wk(1,i)*zln(1,k)-wk(2,i)*zln(2,k)
            wk(2,j)=wk(2,j)-wk(1,i)*zln(3,k)-wk(2,i)*zln(4,k)
  210    continue
  200 continue
!c permutaion
      do 300 i=1,neqns
         b(1,iperm(i))=wk(1,i)
         b(2,iperm(i))=wk(2,i)
  300 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nusol3(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns, &
                       nstop)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),iperm(*)
      double precision zln(9,*),diag(6,*),b(3,*),wk(3,*),dsln(9,*)
!c forward
      do 10 i=1,neqns
         wk(1,i)=b(1,iperm(i))
         wk(2,i)=b(2,iperm(i))
         wk(3,i)=b(3,iperm(i))
   10 continue
      joc=1
      do 100 i=1,neqns
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 110
         call s3pdot(wk(1,i),wk,zln,colno,ks,ke)
  110    continue
         if(i.le.nstop) goto 100
         call d3sdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
         joc=joc+i-nstop
  100 continue
      do 120 i=1,neqns
         wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
         wk(3,i)=wk(3,i)-wk(1,i)*diag(4,i)-wk(2,i)*diag(5,i)
         wk(1,i)=wk(1,i)*diag(1,i)
         wk(2,i)=wk(2,i)*diag(3,i)
         wk(3,i)=wk(3,i)*diag(6,i)
         wk(2,i)=wk(2,i)-wk(3,i)*diag(5,i)
         wk(1,i)=wk(1,i)-wk(2,i)*diag(2,i)-wk(3,i)*diag(4,i)
  120 continue
!c back ward
      do 200 i=neqns,1,-1
         if(i.lt.nstop) goto 206
         do 205 j=i-1,nstop,-1
            joc=joc-1
            wk(1,j)=wk(1,j)-wk(1,i)*dsln(1,joc)-wk(2,i)*dsln(2,joc) &
                   -wk(3,i)*dsln(3,joc)
            wk(2,j)=wk(2,j)-wk(1,i)*dsln(4,joc)-wk(2,i)*dsln(5,joc) &
                   -wk(3,i)*dsln(6,joc)
            wk(3,j)=wk(3,j)-wk(1,i)*dsln(7,joc)-wk(2,i)*dsln(8,joc) &
                   -wk(3,i)*dsln(9,joc)
  205    continue
  206    continue
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 200
         do 210 k=ks,ke
            j=colno(k)
            wk(1,j)=wk(1,j) &
                   -wk(1,i)*zln(1,k)-wk(2,i)*zln(2,k)-wk(3,i)*zln(3,k)
            wk(2,j)=wk(2,j) &
                   -wk(1,i)*zln(4,k)-wk(2,i)*zln(5,k)-wk(3,i)*zln(6,k)
            wk(3,j)=wk(3,j) &
                   -wk(1,i)*zln(7,k)-wk(2,i)*zln(8,k)-wk(3,i)*zln(9,k)
  210    continue
  200 continue
!c permutaion
      do 300 i=1,neqns
         b(1,iperm(i))=wk(1,i)
         b(2,iperm(i))=wk(2,i)
         b(3,iperm(i))=wk(3,i)
  300 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nusol6(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns, &
                      nstop)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),iperm(*)
!GP: DEBUG 13May04  wk(3 ---> wk(6, b(3 ---> b(6, diag(6 ---> diag(21,
!GP: DEBUG 13May04  zln(9 ---> zln(36, dsln(9 ---> dsln(36
!      double precision zln(9,*),diag(6,*),b(3,*),wk(3,*),dsln(9,*)
      double precision zln(36,*),diag(21,*),b(6,*),wk(6,*),dsln(36,*)
!c forward
      do 10 i=1,neqns
         wk(1,i)=b(1,iperm(i))
         wk(2,i)=b(2,iperm(i))
         wk(3,i)=b(3,iperm(i))
         wk(4,i)=b(4,iperm(i))
         wk(5,i)=b(5,iperm(i))
         wk(6,i)=b(6,iperm(i))
   10 continue
      joc=1
      do 100 i=1,neqns
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 110
         call s6pdot(wk(1,i),wk,zln,colno,ks,ke)
  110    continue
         if(i.le.nstop) goto 100
!c        call d6sdot(wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
         call dxsdot(6,wk(1,i),wk(1,nstop),dsln(1,joc),i-nstop)
         joc=joc+i-nstop
  100 continue
      do 120 i=1,neqns
         wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
         wk(3,i)=wk(3,i)-wk(1,i)*diag(4,i)-wk(2,i)*diag(5,i)
         wk(4,i)=wk(4,i)-wk(1,i)*diag(7,i)-wk(2,i)*diag(8,i) &
                -wk(3,i)*diag(9,i)
         wk(5,i)=wk(5,i)-wk(1,i)*diag(11,i)-wk(2,i)*diag(12,i) &
                -wk(3,i)*diag(13,i)-wk(4,i)*diag(14,i)
         wk(6,i)=wk(6,i)-wk(1,i)*diag(16,i)-wk(2,i)*diag(17,i) &
             -wk(3,i)*diag(18,i)-wk(4,i)*diag(19,i)-wk(6,i)*diag(20,i)
         wk(1,i)=wk(1,i)*diag(1,i)
         wk(2,i)=wk(2,i)*diag(3,i)
         wk(3,i)=wk(3,i)*diag(6,i)
         wk(4,i)=wk(4,i)*diag(10,i)
         wk(5,i)=wk(5,i)*diag(15,i)
         wk(6,i)=wk(6,i)*diag(21,i)
         wk(5,i)=wk(5,i)-wk(6,i)*diag(20,i)
         wk(4,i)=wk(4,i)-wk(6,i)*diag(19,i)-wk(5,i)*diag(14,i)
         wk(3,i)=wk(3,i)-wk(6,i)*diag(18,i)-wk(5,i)*diag(13,i) &
                -wk(4,i)*diag(9,i)
         wk(2,i)=wk(2,i)-wk(6,i)*diag(17,i)-wk(5,i)*diag(12,i) &
                -wk(4,i)*diag(8,i)-wk(3,i)*diag(5,i)
         wk(1,i)=wk(1,i)-wk(6,i)*diag(16,i)-wk(5,i)*diag(11,i) &
                -wk(4,i)*diag(7,i)-wk(3,i)*diag(4,i)-wk(2,i)*diag(2,i)
  120 continue
!c back ward
      do 200 i=neqns,1,-1
         if(i.lt.nstop) goto 206
         do 205 j=i-1,nstop,-1
            joc=joc-1
            wk(1,j)=wk(1,j)-wk(1,i)*dsln(1,joc)-wk(2,i)*dsln(2,joc) &
                   -wk(3,i)*dsln(3,joc)-wk(4,i)*dsln(4,joc) &
                   -wk(5,i)*dsln(5,joc)-wk(6,i)*dsln(6,joc)
            wk(2,j)=wk(2,j)-wk(1,i)*dsln(7,joc)-wk(2,i)*dsln(8,joc) &
                   -wk(3,i)*dsln(9,joc)-wk(4,i)*dsln(10,joc) &
                   -wk(5,i)*dsln(11,joc)-wk(6,i)*dsln(12,joc)
            wk(3,j)=wk(3,j)-wk(1,i)*dsln(13,joc)-wk(2,i)*dsln(14,joc) &
                   -wk(3,i)*dsln(15,joc)-wk(4,i)*dsln(16,joc) &
                   -wk(5,i)*dsln(17,joc)-wk(6,i)*dsln(18,joc)
            wk(4,j)=wk(4,j)-wk(1,i)*dsln(19,joc)-wk(2,i)*dsln(20,joc) &
                   -wk(3,i)*dsln(21,joc)-wk(4,i)*dsln(22,joc) &
                   -wk(5,i)*dsln(23,joc)-wk(6,i)*dsln(24,joc)
            wk(5,j)=wk(5,j)-wk(1,i)*dsln(25,joc)-wk(2,i)*dsln(26,joc) &
                   -wk(3,i)*dsln(27,joc)-wk(4,i)*dsln(28,joc) &
                   -wk(5,i)*dsln(29,joc)-wk(6,i)*dsln(30,joc)
            wk(6,j)=wk(6,j)-wk(1,i)*dsln(31,joc)-wk(2,i)*dsln(32,joc) &
                   -wk(3,i)*dsln(33,joc)-wk(4,i)*dsln(34,joc) &
                   -wk(5,i)*dsln(35,joc)-wk(6,i)*dsln(36,joc)
  205    continue
  206    continue
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 200
         do 210 k=ks,ke
            j=colno(k)
            wk(1,j)=wk(1,j)-wk(1,i)*zln(1,joc)-wk(2,i)*zln(2,joc) &
                   -wk(3,i)*zln(3,joc)-wk(4,i)*zln(4,joc) &
                   -wk(5,i)*zln(5,joc)-wk(6,i)*zln(6,joc)
            wk(2,j)=wk(2,j)-wk(1,i)*zln(7,joc)-wk(2,i)*zln(8,joc) &
                   -wk(3,i)*zln(9,joc)-wk(4,i)*zln(10,joc) &
                   -wk(5,i)*zln(11,joc)-wk(6,i)*zln(12,joc)
            wk(3,j)=wk(3,j)-wk(1,i)*zln(13,joc)-wk(2,i)*zln(14,joc) &
                   -wk(3,i)*zln(15,joc)-wk(4,i)*zln(16,joc) &
                   -wk(5,i)*zln(17,joc)-wk(6,i)*zln(18,joc)
            wk(4,j)=wk(4,j)-wk(1,i)*zln(19,joc)-wk(2,i)*zln(20,joc) &
                   -wk(3,i)*zln(21,joc)-wk(4,i)*zln(22,joc) &
                   -wk(5,i)*zln(23,joc)-wk(6,i)*zln(24,joc)
            wk(5,j)=wk(5,j)-wk(1,i)*zln(25,joc)-wk(2,i)*zln(26,joc) &
                   -wk(3,i)*zln(27,joc)-wk(4,i)*zln(28,joc) &
                   -wk(5,i)*zln(29,joc)-wk(6,i)*zln(30,joc)
            wk(6,j)=wk(6,j)-wk(1,i)*zln(31,joc)-wk(2,i)*zln(32,joc) &
                   -wk(3,i)*zln(33,joc)-wk(4,i)*zln(34,joc) &
                   -wk(5,i)*zln(35,joc)-wk(6,i)*zln(36,joc)
  210    continue
  200 continue
!c permutaion
      do 300 i=1,neqns
         b(1,iperm(i))=wk(1,i)
         b(2,iperm(i))=wk(2,i)
         b(3,iperm(i))=wk(3,i)
         b(4,iperm(i))=wk(4,i)
         b(5,iperm(i))=wk(5,i)
         b(6,iperm(i))=wk(6,i)
  300 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine nusolx(xlnzr,colno,dsln,zln,diag,iperm,b,wk,neqns, &
                      nstop,ndeg,ndegl)
!c
!c
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),iperm(*)
      double precision zln(ndeg,ndeg,*),diag(ndegl,*),b(ndeg,*),&
                       wk(ndeg,*),dsln(ndeg,ndeg,*)
!c forward
      do 10 l=1,ndeg
      do 11 i=1,neqns
         wk(l,i)=b(l,iperm(i))
   11 continue
   10 continue
      joc=1
      do 100 i=1,neqns
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 110
         call sxpdot(ndeg,wk(1,i),wk,zln,colno,ks,ke)
  110    continue
         if(i.le.nstop) goto 100
         call dxsdot(ndeg,wk(1,i),wk(1,nstop),dsln(1,1,joc),i-nstop)
         joc=joc+i-nstop
  100 continue
      do 150 i=1,neqns
         locd=0
         do 120 m=1,ndeg-1
            locd=locd+m
            loc1=locd+m
            do 121 n=m+1,ndeg
               wk(n,i)=wk(n,i)-wk(m,i)*diag(loc1,i)
               loc1=loc1+n
  121    continue
  120    continue
         locd=0
         do 130 m=1,ndeg
            locd=locd+m
            wk(m,i)=wk(m,i)*diag(locd,i)
  130    continue
         do 140 n=ndeg,2,-1
         locd=locd-1
         do 141 m=n-1,1,-1
            wk(m,i)=wk(m,i)-wk(n,i)*diag(locd,i)
            locd=locd-1
  141    continue
  140    continue
  150 continue
!c back ward
      do 200 i=neqns,1,-1
         if(i.lt.nstop) goto 208
         do 205 j=i-1,nstop,-1
            joc=joc-1
            do 206 m=1,ndeg
            do 207 n=1,ndeg
            wk(m,j)=wk(m,j)-wk(n,i)*dsln(n,m,joc)
  207    continue
  206    continue
  205    continue
  208    continue
         ks=xlnzr(i)
         ke=xlnzr(i+1)-1
         if(ke.lt.ks) goto 200
         do 210 k=ks,ke
            j=colno(k)
            do 211 m=1,ndeg
            do 212 n=1,ndeg
               wk(m,j)=wk(m,j)-wk(n,i)*zln(n,m,k)
  212    continue
  211    continue
  210    continue
  200 continue
!c permutaion
      do 300 l=1,ndeg
      do 301 i=1,neqns
         b(l,iperm(i))=wk(l,i)
  301 continue
  300 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine posord(parent,btree,invp,iperm,pordr,nch, &
                       neqns,iw,qarent,mch)
      integer  parent(*),pordr(*),btree(2,*),nch(*),invp(*),iperm(*),&
              iw(*),qarent(*),mch(0:neqns+1)
      common /debug/ idbg1
!c
      do 5 i=1,neqns
         mch(i)=0
         pordr(i)=0
    5 continue
      l=1
      locc=neqns+1
   10 continue
      joc=locc
      locc=btree(1,joc)
      if(locc.ne.0) goto 10
      locp=qarent(joc)
      mch(locp)=mch(locp)+1
   20 continue
      pordr(joc)=l
      if(l.ge.neqns) goto 1000
      l=l+1
      locc=btree(2,joc)
      if(locc.ne.0) goto 10
      joc=qarent(joc)
      locp=qarent(joc)
      mch(locp)=mch(locp)+mch(joc)+1
      goto 20
 1000 continue
      do 100 i=1,neqns
         ipinv=pordr(invp(i))
         invp(i)=ipinv
         iperm(ipinv)=i
         iw(pordr(i))=i
  100 continue
      do 110 i=1,neqns
         invpos=iw(i)
         nch(i)=mch(invpos)
         ii=qarent(invpos)
         if(ii.gt.0.and.ii.le.neqns) then
            parent(i)=pordr(ii)
         else
            parent(i)=qarent(invpos)
         endif
  110 continue
      if(idbg1.ne.0) write(6,6020)
      if(idbg1.ne.0) write(6,6000) (pordr(i),i=1,neqns)
      if(idbg1.ne.0) write(6,6030)
      if(idbg1.ne.0) write(6,6050)
      if(idbg1.ne.0) write(6,6000) (parent(i),i=1,neqns)
      if(idbg1.ne.0) write(6,6000) (invp(i),i=1,neqns)
      if(idbg1.ne.0) write(6,6040)
      if(idbg1.ne.0) write(6,6000) (iperm(i),i=1,neqns)
      if(idbg1.ne.0) write(6,6010)
      if(idbg1.ne.0) write(6,6000) (nch(i),i=1,neqns)
 6000 format(10i6)
 6010 format(' nch')
 6020 format(' post order')
 6030 format(/' invp ')
 6040 format(/' iperm ')
 6050 format(/' parent')
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine pre_gnclno(parent,pordr,xleaf,leaf,xlnzr, &
                        neqns,nstop,lncol,ir)
!c
!c
      integer parent(*),pordr(*),xleaf(*),leaf(*),xlnzr(*)
      common /debug/ idbg1
!c
      nc=0
      ir=0
      l=1
      do 100 i=1,neqns
         xlnzr(i)=l
         ks=xleaf(i)
         ke=xleaf(i+1)-1
         if(ke.lt.ks) goto 100
         nxleaf=leaf(ks)
         do 110 k=ks,ke-1
            j=nxleaf
            nxleaf=leaf(k+1)
  105       continue
            if(j.ge.nxleaf) goto 110
            if(j.ge.nstop) then
               goto 100
            endif
            l=l+1
            j=parent(j)
            goto 105
  110    continue
         j=leaf(ke)
  115    continue
         if(j.ge.nstop) goto 100
         if(j.ge.i.or.j.eq.0) goto 100
         l=l+1
         j=parent(j)
         goto 115
  100 continue
      xlnzr(neqns+1)=l
      lncol=l-1
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine prt(ip,n)
      dimension ip(n)
      write(6,6000) (ip(i),i=1,n)
 6000 format(10(2x,i4))
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine ptime(cputim)
      use hecmw_util
      double precision cputim,elaptime
      real x(2)
!c machine dependent
!c**********************************************************************
!c cpu time by hour
!c     cputim=etime(x)
!c     cputim=x(1)
      cputim=hecmw_Wtime()
!c**********************************************************************
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine qmdmrg(xadj,adjncy,deg,qsize,qlink,marker, &
                       deg0,nhdsze,nbrhd,rchset,ovrlp)
      implicit none
      integer adjncy(*),deg(*),qsize(*),qlink(*),marker(*),rchset(*), &
             nbrhd(*),ovrlp(*)
      integer xadj(*),deg0,deg1,head,inhd,iov,irch,&
             j,jstrt,jstop,link,lnode,mark,mrgsze,&
             nabor,nhdsze,node,novrlp,rchsze,root
!c
      if(nhdsze.le.0) return
      do 100 inhd=1,nhdsze
         root=nbrhd(inhd)
         marker(root)=0
  100 continue
      do 1400 inhd=1,nhdsze
         root=nbrhd(inhd)
         marker(root)=-1
         rchsze=0
         novrlp=0
         deg1=0
  200    jstrt=xadj(root)
         jstop=xadj(root+1)-1
         do 600 j=jstrt,jstop
            nabor=adjncy(j)
            root=-nabor
!           if(nabor) 200,700,300
            if(nabor.lt.0)goto 200
            if(nabor.eq.0)goto 700
! 300       mark=marker(nabor)
            mark=marker(nabor)
!           if(mark)600,400,500
            if(mark.lt.0)goto 600
            if(mark.gt.0)goto 500
! 400       rchsze=rchsze+1
            rchsze=rchsze+1
            rchset(rchsze)=nabor
            deg1=deg1+qsize(nabor)
            marker(nabor)=1
            goto 600
  500       if(mark.gt.1) goto 600
            novrlp=novrlp+1
            ovrlp(novrlp)=nabor
            marker(nabor)=2
  600    continue
  700    head=0
         mrgsze=0
         do 1100 iov=1,novrlp
            node=ovrlp(iov)
            jstrt=xadj(node)
            jstop=xadj(node+1)-1
            do 800 j=jstrt,jstop
               nabor=adjncy(j)
               if(marker(nabor).ne.0) goto 800
               marker(node)=1
               goto 1100
  800       continue
            mrgsze=mrgsze+qsize(node)
            marker(node)=-1
            lnode=node
  900       link=qlink(lnode)
            if(link.le.0) goto 1000
            lnode=link
            goto 900
 1000       qlink(lnode)=head
            head=node
 1100    continue
         if(head.le.0) goto 1200
         qsize(head)=mrgsze
         deg(head)=deg0+deg1-1
         marker(head)=2
 1200    root=nbrhd(inhd)
         marker(root)=0
         if(rchsze.le.0) goto 1400
         do 1300 irch=1,rchsze
            node=rchset(irch)
            marker(node)=0
 1300    continue
 1400 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine qmdot(root,xadj,adjncy,marker,rchsze,rchset,nbrhd)
!c
!c
      integer adjncy(*),marker(*),rchset(*),nbrhd(*)
      integer xadj(*),inhd,irch,j,jstrt,jstop,link, &
              nabor,node,rchsze,root
!c
      irch=0
      inhd=0
      node=root
  100 jstrt=xadj(node)
      jstop=xadj(node+1)-2
      if(jstop.lt.jstrt) goto 300
      do 200 j=jstrt,jstop
         irch=irch+1
         adjncy(j)=rchset(irch)
         if(irch.ge.rchsze) goto 400
  200 continue
  300 link=adjncy(jstop+1)
      node=-link
      if(link.lt.0) goto 100
      inhd=inhd+1
      node=nbrhd(inhd)
      adjncy(jstop+1)=-node
      goto 100
  400 adjncy(j+1)=0
      do 600 irch=1,rchsze
         node=rchset(irch)
         if(marker(node).lt.0) goto 600
         jstrt=xadj(node)
         jstop=xadj(node+1)-1
         do 500 j=jstrt,jstop
            nabor=adjncy(j)
            if(marker(nabor).ge.0) goto 500
            adjncy(j)=root
            goto 600
  500    continue
  600 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!

      subroutine qmdrch(root,xadj,adjncy,deg,marker, &
                        rchsze,rchset,nhdsze,nbrhd)
      integer adjncy(*),deg(*),marker(*),rchset(*),nbrhd(*)
      integer xadj(*),i,istrt,istop,j,jstrt,jstop, &
              nabor,nhdsze,node,rchsze,root
!c
      nhdsze=0
      rchsze=0
      istrt=xadj(root)
      istop=xadj(root+1)-1
      if(istop.lt.istrt) return
      do 600 i=istrt,istop
         nabor=adjncy(i)
         if(nabor.eq.0) return
         if(marker(nabor).ne.0) goto 600
         if(deg(nabor).lt.0) goto 200
         rchsze=rchsze+1
         rchset(rchsze)=nabor
         marker(nabor)=1
         goto 600
  200    marker(nabor)=-1
         nhdsze=nhdsze+1
         nbrhd(nhdsze)=nabor
  300    jstrt=xadj(nabor)
         jstop=xadj(nabor+1)-1
         do 500 j=jstrt,jstop
            node=adjncy(j)
            nabor=-node
!           if(node) 300,600,400
            if(node.lt.0) goto 300
            if(node.eq.0) goto 600
! 400       if(marker(node).ne.0) goto 500
            if(marker(node).ne.0) goto 500
            rchsze=rchsze+1
            rchset(rchsze)=node
            marker(node)=1
  500    continue
  600 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine qmdupd(xadj,adjncy,nlist,list,deg,qsize,qlink, &
                        marker,rchset,nbrhd)
      integer adjncy(*),list(*),deg(*),marker(*), &
              rchset(*),nbrhd(*),qsize(*),qlink(*)
      integer xadj(*),deg0,deg1,il,inhd,inode,irch, &
              j,jstrt,jstop,mark,nabor,nhdsze,nlist, &
              node,rchsze
!c
      if(nlist.le.0) return
      deg0=0
      nhdsze=0
      do 200 il=1,nlist
         node=list(il)
         deg0=deg0+qsize(node)
         jstrt=xadj(node)
         jstop=xadj(node+1)-1
         do 100 j=jstrt,jstop
            nabor=adjncy(j)
            if(marker(nabor).ne.0.or. &
               deg(nabor).ge.0) goto 100
            marker(nabor)=-1
            nhdsze=nhdsze+1
            nbrhd(nhdsze)=nabor
  100    continue
  200 continue
!c
      if(nhdsze.gt.0) &
         call qmdmrg(xadj,adjncy,deg,qsize,qlink,marker,deg0,nhdsze, &
                     nbrhd,rchset,nbrhd(nhdsze+1))
      do 600 il=1,nlist
         node=list(il)
         mark=marker(node)
         if(mark.gt.1.or.mark.lt.0) goto 600
         call qmdrch(node,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze, &
                     nbrhd)
         deg1=deg0
         if(rchsze.le.0) goto 400
         do 300 irch=1,rchsze
            inode=rchset(irch)
            deg1=deg1+qsize(inode)
            marker(inode)=0
  300    continue
  400    deg(node)=deg1-1
         if(nhdsze.le.0) goto 600
         do 500 inhd=1,nhdsze
            inode=nbrhd(inhd)
            marker(inode)=0
  500    continue
  600 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine qqsort(iw,ik)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer iw(*)
      common /debug/ idbg1
!c
!c----------------------------------------------------------------------
!c     sort in increasing order up to i
!c
!c     iw   array
!c     ik   number of input/output
!c     i    deal with numbers less than this numberi
!c
!c----------------------------------------------------------------------
!c
      if(ik.le.1) return
      do 100 l=1,ik-1
         do 110 m=l+1,ik
            if(iw(l).lt.iw(m)) goto 110
            itemp=iw(l)
            iw(l)=iw(m)
            iw(m)=itemp
  110    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine rotate(xadj,adjncy,invp,iperm,parent,btree,izz,neqns, &
                       anc,adjt,irr)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer  xadj(*),adjncy(*),parent(*),btree(2,*),invp(*),iperm(*), &
              anc(*),adjt(*)
      common /debug/ idbg1
!c----------------------------------------------------------------------
!c     irr return code irr=0 node izz is not a bottom node
!c                     irr=1          is a bottom node then rotation is
!c                                    performed
!c
!c----------------------------------------------------------------------
      if(izz.eq.0) then
          irr=0
          return
      endif
      izzz=invp(izz)
      if(btree(1,izzz).ne.0) then
          irr=0
!c         return
      endif
      irr=1
!c
!c  ancestors of izzz
!c
      nanc=0
      joc=izzz
  100 continue
      nanc=nanc+1
      anc(nanc)=joc
      joc=parent(joc)
      if(joc.ne.0) goto 100
!c
!c  to find the eligible node from ancestors of izz
!c
!c     adjt = Adj(Tree(y))
      l=1
  200 continue
      do 210 i=1,neqns
         adjt(i)=0
  210 continue
      locc=anc(l)
  220 continue
      joc=locc
      locc=btree(1,joc)
      if(locc.ne.0) goto 220
  230 continue
      do 240 k=xadj(iperm(joc)),xadj(iperm(joc)+1)-1
         adjt(invp(adjncy(k)))=1
  240 continue
      if(joc.ge.anc(l)) goto 250
      locc=btree(2,joc)
      if(locc.ne.0) goto 220
      joc=parent(joc)
      goto 230
  250 continue
      do 260 ll=l+1,nanc
         if(adjt(anc(ll)).eq.0) then
            l=l+1
            goto 200
         endif
  260 continue
      if(l.eq.1) goto 500
       
!c
!c  anc(l-1) is the eligible node
!c
!c (1) number the node not in Ancestor(iy)
      iy=anc(l-1)
      do 300 i=1,neqns
         adjt(i)=0
  300 continue
      do 310 ll=l,nanc
         adjt(anc(ll))=1
  310 continue
      k=0
      do 320 ll=1,neqns
         if(adjt(ll).eq.0) then
            k=k+1
            invp(iperm(ll))=k
         endif
  320 continue
!c (2) followed by nodes in Ancestor(iy)-Adj(T(iy))
      do 340 i=1,neqns
         adjt(i)=0
  340 continue
      locc=iy
  350 continue
      joc=locc
      locc=btree(1,joc)
      if(locc.ne.0) goto 350
  360 continue
      do 370 kk=xadj(iperm(joc)),xadj(iperm(joc)+1)-1
         adjt(invp(adjncy(kk)))=1
  370 continue
      if(joc.ge.iy) goto 380
      locc=btree(2,joc)
      if(locc.ne.0) goto 350
      joc=parent(joc)
      goto 360
  380 continue
      do 390 ll=l,nanc
         if(adjt(anc(ll)).eq.0) then
            k=k+1
            invp(iperm(anc(ll)))=k
         endif
  390 continue
!c (3) and finally number the node in Adj(t(iy))
      do 400 ll=l,nanc
         if(adjt(anc(ll)).ne.0) then
            k=k+1
            invp(iperm(anc(ll)))=k
         endif
  400 continue
      goto 600
!c
!c izz can be numbered last
!c
  500 continue
      k=0
      do 510 i=1,neqns
         if(i.eq.izzz) goto 510
         k=k+1
         invp(iperm(i))=k
  510 continue
      invp(iperm(izzz))=neqns
!c
!c set iperm
!c
  600 continue
      do 610 i=1,neqns
         iperm(invp(i))=i
  610 continue
      if(idbg1.ne.0) write(6,6000) (invp(i),i=1,neqns)
 6000 format(10i6)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s2pdot(bi,b,zln,colno,ks,ke)
      implicit double precision(a-h,o-z)
      integer colno(*)
      double precision zln(4,*),b(2,*),bi(2)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      do 100 jj=ks,ke
         j=colno(jj)
         bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(3,jj)*b(2,j)
         bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(4,jj)*b(2,j)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s2um(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),nch(*),par(*)
      dimension zln(4,*),diag(3,*),temp(4,*),indx(*)
      dimension s(4),zz(4),t(3)
      common /mchdpn/ rmax,rmin,epsm,lratio
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         t(1)=0.0d0
         t(2)=0.0d0
         t(3)=0.0d0
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            do 210 l=1,4
               s(l)=0.0d0
               zz(l)=zln(l,k)
  210       continue
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(3,j)*zln(3,jj)
                  zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(4,j)*zln(3,jj)
                  zz(3)=zz(3)-temp(1,j)*zln(2,jj)-temp(3,j)*zln(4,jj)
                  zz(4)=zz(4)-temp(2,j)*zln(2,jj)-temp(4,j)*zln(4,jj)
               endif
  310       continue
            call inv22(zln(1,k),zz,diag(1,jc))
            do 220 l=1,4
               temp(l,jc)=zz(l)
  220       continue
            t(1)=t(1)+zz(1)*zln(1,k)+zz(3)*zln(3,k)
            t(2)=t(2)+zz(1)*zln(2,k)+zz(3)*zln(4,k)
            t(3)=t(3)+zz(2)*zln(2,k)+zz(4)*zln(4,k)
  200    continue
         diag(1,ic)=diag(1,ic)-t(1)
         diag(2,ic)=diag(2,ic)-t(2)
         diag(3,ic)=diag(3,ic)-t(3)
         call inv2(diag(1,ic),ir)
         nch(ic)=-1
         kk=par(ic)
         nch(kk)=nch(kk)-1
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s2um1(ic,xlnzr,colno,zln,diag,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),par(*)
      dimension zln(4,*),diag(3,*),temp(4,*),indx(*)
      dimension s(4),zz(4)
      common /mchdpn/ rmax,rmin,epsm,lratio
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         do 100 l=1,4
            s(l)=0.0d0
  100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(3,j)*zln(3,jj)
                  s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(4,j)*zln(3,jj)
                  s(3)=s(3)+temp(1,j)*zln(2,jj)+temp(3,j)*zln(4,jj)
                  s(4)=s(4)+temp(2,j)*zln(2,jj)+temp(4,j)*zln(4,jj)
               endif
  310       continue
            do 320 l=1,4
              temp(l,jc)=zln(l,k)-s(l)
              zln(l,k)=temp(l,jc)
              s(l)=0.0d0
  320       continue
  200    continue
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s2um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*)
      dimension zln(4,*),diag(3,*),temp(4,*),indx(*),dsln(4,*)
      joc=0
      do 100 ic=nstop,neqns
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)-1
         do 110 k=ks,ke
            jj=colno(k)
            temp(1,jj)=zln(1,k)
            temp(2,jj)=zln(2,k)
            temp(3,jj)=zln(3,k)
            temp(4,jj)=zln(4,k)
!c           call inv22(zln(1,k),temp(1,jj),diag(1,jj))
            zln(3,k)=temp(3,jj)-temp(1,jj)*diag(2,jj)
            zln(1,k)=temp(1,jj)*diag(1,jj)
            zln(3,k)=zln(3,k)*diag(3,jj)
            zln(1,k)=zln(1,k)-zln(3,k)*diag(2,jj)
!c
            zln(4,k)=temp(4,jj)-temp(2,jj)*diag(2,jj)
            zln(2,k)=temp(2,jj)*diag(1,jj)
            zln(4,k)=zln(4,k)*diag(3,jj)
            zln(2,k)=zln(2,k)-zln(4,k)*diag(2,jj)
!c
            diag(1,ic)=diag(1,ic) &
                       -(temp(1,jj)*zln(1,k)+temp(3,jj)*zln(3,k))
            diag(2,ic)=diag(2,ic) &
                       -(temp(1,jj)*zln(2,k)+temp(3,jj)*zln(4,k))
            diag(3,ic)=diag(3,ic) &
                       -(temp(2,jj)*zln(2,k)+temp(4,jj)*zln(4,k))
            indx(jj)=ic
  110    continue
         do 120 jc=nstop,ic-1
            joc=joc+1
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  dsln(1,joc)=dsln(1,joc) &
                           -(temp(1,j)*zln(1,jj)+temp(3,j)*zln(3,jj))
                  dsln(2,joc)=dsln(2,joc) &
                           -(temp(2,j)*zln(1,jj)+temp(4,j)*zln(3,jj))
                  dsln(3,joc)=dsln(3,joc) &
                           -(temp(1,j)*zln(2,jj)+temp(3,j)*zln(4,jj))
                  dsln(4,joc)=dsln(4,joc) &
                           -(temp(2,j)*zln(2,jj)+temp(4,j)*zln(4,jj))
               endif
  220       continue
  120    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s2um3(n,dsln,diag,indx,temp)
      implicit double precision(a-h,o-z)
      dimension dsln(4,*),diag(3,*),indx(*),temp(4,*),t(4)
      if(n.le.0) goto 1000
      indx(1)=0
      joc=1
      call inv2(diag(1,1),ir)
      do 100 i=2,n
         indx(i)=joc
         do 110 j=1,i-1
            call d2dot(t,dsln(1,indx(i)),dsln(1,indx(j)),j-1)
            dsln(1,joc)=dsln(1,joc)-t(1)
            dsln(2,joc)=dsln(2,joc)-t(2)
            dsln(3,joc)=dsln(3,joc)-t(3)
            dsln(4,joc)=dsln(4,joc)-t(4)
            joc=joc+1
  110    continue
         call v2prod(dsln(1,indx(i)),diag,temp,i-1)
         call d2dot(t,temp,dsln(1,indx(i)),i-1)
         diag(1,i)=diag(1,i)-t(1)
         diag(2,i)=diag(2,i)-t(2)
         diag(3,i)=diag(3,i)-t(4)
         call vcopy(temp,dsln(1,indx(i)),4*(i-1))
         call inv2(diag(1,i),ir)
  100 continue
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s3pdot(bi,b,zln,colno,ks,ke)
      implicit double precision(a-h,o-z)
      integer colno(*)
      double precision zln(9,*),b(3,*),bi(3)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      do 100 jj=ks,ke
         j=colno(jj)
         bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(4,jj)*b(2,j)-zln(7,jj)*b(3,j)
         bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(5,jj)*b(2,j)-zln(8,jj)*b(3,j)
         bi(3)=bi(3)-zln(3,jj)*b(1,j)-zln(6,jj)*b(2,j)-zln(9,jj)*b(3,j)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s3um(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),nch(*),par(*)
      dimension zln(9,*),diag(6,*),temp(9,*),indx(*)
      dimension zz(9),t(6)
      common /mchdpn/ rmax,rmin,epsm,lratio
!c     do ic=2,nstop-1
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
!c$dir max_trips(6)
         do 100 l=1,6
            t(l)=0.0d0
  100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
!c$dir max_trips(9)
            do 210 l=1,9
               zz(l)=zln(l,k)
  210       continue
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(4,j)*zln(4,jj) &
                             -temp(7,j)*zln(7,jj)
                  zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(5,j)*zln(4,jj) &
                             -temp(8,j)*zln(7,jj)
                  zz(3)=zz(3)-temp(3,j)*zln(1,jj)-temp(6,j)*zln(4,jj) &
                             -temp(9,j)*zln(7,jj)
                  zz(4)=zz(4)-temp(1,j)*zln(2,jj)-temp(4,j)*zln(5,jj) &
                             -temp(7,j)*zln(8,jj)
                  zz(5)=zz(5)-temp(2,j)*zln(2,jj)-temp(5,j)*zln(5,jj) &
                             -temp(8,j)*zln(8,jj)
                  zz(6)=zz(6)-temp(3,j)*zln(2,jj)-temp(6,j)*zln(5,jj) &
                             -temp(9,j)*zln(8,jj)
                  zz(7)=zz(7)-temp(1,j)*zln(3,jj)-temp(4,j)*zln(6,jj) &
                             -temp(7,j)*zln(9,jj)
                  zz(8)=zz(8)-temp(2,j)*zln(3,jj)-temp(5,j)*zln(6,jj) &
                             -temp(8,j)*zln(9,jj)
                  zz(9)=zz(9)-temp(3,j)*zln(3,jj)-temp(6,j)*zln(6,jj) &
                             -temp(9,j)*zln(9,jj)
               endif
  310       continue
            call inv33(zln(1,k),zz,diag(1,jc))
!c$dir max_trips(9)
            do 220 l=1,9
               temp(l,jc)=zz(l)
  220       continue
            t(1)=t(1)+zz(1)*zln(1,k)+zz(4)*zln(4,k)+zz(7)*zln(7,k)
            t(2)=t(2)+zz(1)*zln(2,k)+zz(4)*zln(5,k)+zz(7)*zln(8,k)
            t(3)=t(3)+zz(2)*zln(2,k)+zz(5)*zln(5,k)+zz(8)*zln(8,k)
            t(4)=t(4)+zz(1)*zln(3,k)+zz(4)*zln(6,k)+zz(7)*zln(9,k)
            t(5)=t(5)+zz(2)*zln(3,k)+zz(5)*zln(6,k)+zz(8)*zln(9,k)
            t(6)=t(6)+zz(3)*zln(3,k)+zz(6)*zln(6,k)+zz(9)*zln(9,k)
  200    continue
!c$dir max_trips(6)
         do 320 l=1,6
            diag(l,ic)=diag(l,ic)-t(l)
  320    continue
         call inv3(diag(1,ic),ir)
         nch(ic)=-1
         kk=par(ic)
         nch(kk)=nch(kk)-1
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s3um1(ic,xlnzr,colno,zln,diag,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),par(*)
      dimension zln(9,*),diag(6,*),temp(9,*),indx(*)
      dimension s(9),zz(9)
      common /mchdpn/ rmax,rmin,epsm,lratio
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
!c$dir max_trip(9)
         do 100 l=1,9
            s(l)=0.0d0
  100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(4,j)*zln(4,jj) &
                             +temp(7,j)*zln(7,jj)
                  s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(5,j)*zln(4,jj) &
                             +temp(8,j)*zln(7,jj)
                  s(3)=s(3)+temp(3,j)*zln(1,jj)+temp(6,j)*zln(4,jj) &
                             +temp(9,j)*zln(7,jj)
                  s(4)=s(4)+temp(1,j)*zln(2,jj)+temp(4,j)*zln(5,jj) &
                             +temp(7,j)*zln(8,jj)
                  s(5)=s(5)+temp(2,j)*zln(2,jj)+temp(5,j)*zln(5,jj) &
                             +temp(8,j)*zln(8,jj)
                  s(6)=s(6)+temp(3,j)*zln(2,jj)+temp(6,j)*zln(5,jj) &
                             +temp(9,j)*zln(8,jj)
                  s(7)=s(7)+temp(1,j)*zln(3,jj)+temp(4,j)*zln(6,jj) &
                             +temp(7,j)*zln(9,jj)
                  s(8)=s(8)+temp(2,j)*zln(3,jj)+temp(5,j)*zln(6,jj) &
                             +temp(8,j)*zln(9,jj)
                  s(9)=s(9)+temp(3,j)*zln(3,jj)+temp(6,j)*zln(6,jj) &
                             +temp(9,j)*zln(9,jj)
               endif
  310       continue
!c$dir max_trip(9)
            do 320 l=1,9
              temp(l,jc)=zln(l,k)-s(l)
              zln(l,k)=temp(l,jc)
              s(l)=0.0d0
  320       continue
  200    continue
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s3um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*)
      dimension zln(9,*),diag(6,*),temp(neqns,9),dsln(9,*),indx(*)
      joc=0
      do 100 ic=nstop,neqns
!c        do 105 m=1,9
!c           do 105 jj=1,nstop
!c              temp(jj,m)=0.0d0
!c 105    continue
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)-1
         do 110 k=ks,ke
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
  110    continue
         do 111 k=ks,ke
            jj=colno(k)
      zln(4,k)=temp(jj,4)-temp(jj,1)*diag(2,jj)
      zln(7,k)=temp(jj,7)-temp(jj,1)*diag(4,jj)-zln(4,k)*diag(5,jj)
      zln(1,k)=temp(jj,1)*diag(1,jj)
      zln(4,k)=zln(4,k)*diag(3,jj)
      zln(7,k)=zln(7,k)*diag(6,jj)
      zln(4,k)=zln(4,k)-zln(7,k)*diag(5,jj)
      zln(1,k)=zln(1,k)-zln(4,k)*diag(2,jj)-zln(7,k)*diag(4,jj)
!c
      zln(5,k)=temp(jj,5)-temp(jj,2)*diag(2,jj)
      zln(8,k)=temp(jj,8)-temp(jj,2)*diag(4,jj)-zln(5,k)*diag(5,jj)
      zln(2,k)=temp(jj,2)*diag(1,jj)
      zln(5,k)=zln(5,k)*diag(3,jj)
      zln(8,k)=zln(8,k)*diag(6,jj)
      zln(5,k)=zln(5,k)-zln(8,k)*diag(5,jj)
      zln(2,k)=zln(2,k)-zln(5,k)*diag(2,jj)-zln(8,k)*diag(4,jj)
!c
      zln(6,k)=temp(jj,6)-temp(jj,3)*diag(2,jj)
      zln(9,k)=temp(jj,9)-temp(jj,3)*diag(4,jj)-zln(6,k)*diag(5,jj)
      zln(3,k)=temp(jj,3)*diag(1,jj)
      zln(6,k)=zln(6,k)*diag(3,jj)
      zln(9,k)=zln(9,k)*diag(6,jj)
      zln(6,k)=zln(6,k)-zln(9,k)*diag(5,jj)
      zln(3,k)=zln(3,k)-zln(6,k)*diag(2,jj)-zln(9,k)*diag(4,jj)
!c        write(60,6000) k,(zln(llll,k),llll=1,9)
!c6000    format(i6,3d20.10/6x,3d20.10/6x,3d20.10)
  111    continue
!c
         do 112 k=ks,ke
            jj=colno(k)
            diag(1,ic)=diag(1,ic)-temp(jj,1)*zln(1,k) &
                            -temp(jj,4)*zln(4,k)-temp(jj,7)*zln(7,k)
            diag(2,ic)=diag(2,ic)-temp(jj,1)*zln(2,k) &
                            -temp(jj,4)*zln(5,k)-temp(jj,7)*zln(8,k)
            diag(3,ic)=diag(3,ic)-temp(jj,2)*zln(2,k) &
                            -temp(jj,5)*zln(5,k)-temp(jj,8)*zln(8,k)
            diag(4,ic)=diag(4,ic)-temp(jj,1)*zln(3,k) &
                            -temp(jj,4)*zln(6,k)-temp(jj,7)*zln(9,k)
            diag(5,ic)=diag(5,ic)-temp(jj,2)*zln(3,k) &
                            -temp(jj,5)*zln(6,k)-temp(jj,8)*zln(9,k)
            diag(6,ic)=diag(6,ic)-temp(jj,3)*zln(3,k) &
                            -temp(jj,6)*zln(6,k)-temp(jj,9)*zln(9,k)
  112    continue
         do 120 jc=nstop,ic-1
            joc=joc+1
            j1=xlnzr(jc)
            j2=xlnzr(jc+1)
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  dsln(1,joc)=dsln(1,joc)-temp(j,1)*zln(1,jj) &
                             -temp(j,4)*zln(4,jj)-temp(j,7)*zln(7,jj)
                  dsln(2,joc)=dsln(2,joc)-temp(j,2)*zln(1,jj) &
                             -temp(j,5)*zln(4,jj)-temp(j,8)*zln(7,jj)
                  dsln(3,joc)=dsln(3,joc)-temp(j,3)*zln(1,jj) &
                             -temp(j,6)*zln(4,jj)-temp(j,9)*zln(7,jj)
                  dsln(4,joc)=dsln(4,joc)-temp(j,1)*zln(2,jj) &
                             -temp(j,4)*zln(5,jj)-temp(j,7)*zln(8,jj)
                  dsln(5,joc)=dsln(5,joc)-temp(j,2)*zln(2,jj) &
                             -temp(j,5)*zln(5,jj)-temp(j,8)*zln(8,jj)
                  dsln(6,joc)=dsln(6,joc)-temp(j,3)*zln(2,jj) &
                             -temp(j,6)*zln(5,jj)-temp(j,9)*zln(8,jj)
                  dsln(7,joc)=dsln(7,joc)-temp(j,1)*zln(3,jj) &
                             -temp(j,4)*zln(6,jj)-temp(j,7)*zln(9,jj)
                  dsln(8,joc)=dsln(8,joc)-temp(j,2)*zln(3,jj) &
                             -temp(j,5)*zln(6,jj)-temp(j,8)*zln(9,jj)
                  dsln(9,joc)=dsln(9,joc)-temp(j,3)*zln(3,jj) &
                             -temp(j,6)*zln(6,jj)-temp(j,9)*zln(9,jj)
               endif
  220       continue
  120    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s3um3(n,dsln,diag,indx,temp)
      implicit double precision(a-h,o-z)
      dimension dsln(9,*),diag(6,*),indx(*),temp(9,*),t(9)
      if(n.le.0) goto 1000
      indx(1)=0
      joc=1
      call inv3(diag(1,1),ir)
      do 100 i=2,n
         indx(i)=joc
         do 110 j=1,i-1
            call d3dot(t,dsln(1,indx(i)),dsln(1,indx(j)),j-1)
!c$dir max_trips(9)
!c           do 111 l=1,9
!c              dsln(l,joc)=dsln(l,joc)-t(l)
!c 111       continue
            dsln(:,joc)=dsln(:,joc)-t(:)
            joc=joc+1
  110    continue
         call v3prod(dsln(1,indx(i)),diag,temp,i-1)
         call d3dotl(t,temp,dsln(1,indx(i)),i-1)
!c$dir max_trips(6)
!c        do 112 l=1,6
!c           diag(l,i)=diag(l,i)-t(l)
!c 112    continue
         diag(:,i)=diag(:,i)-t(1:6)
         call vcopy(temp,dsln(1,indx(i)),9*(i-1))
         call inv3(diag(1,i),ir)
  100 continue
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s6pdot(bi,b,zln,colno,ks,ke)
      implicit double precision(a-h,o-z)
      integer colno(*)
      double precision zln(36,*),b(6,*),bi(6)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      do 100 jj=ks,ke
         j=colno(jj)
         bi(1)=bi(1) &
                  -zln(1,jj)*b(1,j)-zln(7,jj)*b(2,j)-zln(13,jj)*b(3,j) &
                  -zln(19,jj)*b(4,j)-zln(25,jj)*b(5,j)-zln(31,jj)*b(6,j)
         bi(2)=bi(2) &
                  -zln(2,jj)*b(1,j)-zln(8,jj)*b(2,j)-zln(14,jj)*b(3,j) &
                  -zln(20,jj)*b(4,j)-zln(26,jj)*b(5,j)-zln(32,jj)*b(6,j)
         bi(3)=bi(3) &
                  -zln(3,jj)*b(1,j)-zln(9,jj)*b(2,j)-zln(15,jj)*b(3,j) &
                  -zln(21,jj)*b(4,j)-zln(27,jj)*b(5,j)-zln(33,jj)*b(6,j)
         bi(4)=bi(4) &
                  -zln(4,jj)*b(1,j)-zln(10,jj)*b(2,j)-zln(16,jj)*b(3,j) &
                  -zln(22,jj)*b(4,j)-zln(28,jj)*b(5,j)-zln(34,jj)*b(6,j)
         bi(5)=bi(5) &
                  -zln(5,jj)*b(1,j)-zln(11,jj)*b(2,j)-zln(17,jj)*b(3,j) &
                  -zln(23,jj)*b(4,j)-zln(29,jj)*b(5,j)-zln(35,jj)*b(6,j)
         bi(6)=bi(6) &
                  -zln(6,jj)*b(1,j)-zln(12,jj)*b(2,j)-zln(18,jj)*b(3,j) &
                  -zln(25,jj)*b(4,j)-zln(30,jj)*b(5,j)-zln(36,jj)*b(6,j)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s6um(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),nch(*),par(*)
      dimension zln(36,*),diag(21,*),temp(36,*),indx(*)
      dimension zz(36),t(21)
      common /mchdpn/ rmax,rmin,epsm,lratio
!c     do ic=2,nstop-1
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         do 100 l=1,21
            t(l)=0.0d0
  100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            do 210 l=1,36
               zz(l)=zln(l,k)
  210       continue
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
              zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(7,j)*zln(7,jj) &
                     -temp(13,j)*zln(13,jj)-temp(19,j)*zln(19,jj) &
                     -temp(25,j)*zln(25,jj)-temp(31,j)*zln(31,jj)
              zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(8,j)*zln(7,jj) &
                     -temp(14,j)*zln(13,jj)-temp(20,j)*zln(19,jj) &
                     -temp(26,j)*zln(25,jj)-temp(32,j)*zln(31,jj)
              zz(3)=zz(3)-temp(3,j)*zln(1,jj)-temp(9,j)*zln(7,jj) &
                     -temp(15,j)*zln(13,jj)-temp(21,j)*zln(19,jj) &
                     -temp(27,j)*zln(25,jj)-temp(33,j)*zln(31,jj)
              zz(4)=zz(4)-temp(4,j)*zln(1,jj)-temp(10,j)*zln(7,jj) &
                     -temp(16,j)*zln(13,jj)-temp(22,j)*zln(19,jj) &
                     -temp(28,j)*zln(25,jj)-temp(34,j)*zln(31,jj)
              zz(5)=zz(5)-temp(5,j)*zln(1,jj)-temp(11,j)*zln(7,jj) &
                     -temp(17,j)*zln(13,jj)-temp(23,j)*zln(19,jj) &
                     -temp(29,j)*zln(25,jj)-temp(35,j)*zln(31,jj)
              zz(6)=zz(6)-temp(6,j)*zln(1,jj)-temp(12,j)*zln(7,jj) &
                     -temp(18,j)*zln(13,jj)-temp(24,j)*zln(19,jj) &
                     -temp(30,j)*zln(25,jj)-temp(36,j)*zln(31,jj)
              zz(7)=zz(7)-temp(1,j)*zln(2,jj)-temp(7,j)*zln(8,jj) &
                     -temp(13,j)*zln(14,jj)-temp(19,j)*zln(20,jj) &
                     -temp(25,j)*zln(26,jj)-temp(31,j)*zln(32,jj)
              zz(8)=zz(8)-temp(2,j)*zln(2,jj)-temp(8,j)*zln(8,jj) &
                     -temp(14,j)*zln(14,jj)-temp(20,j)*zln(20,jj) &
                     -temp(26,j)*zln(26,jj)-temp(32,j)*zln(32,jj)
              zz(9)=zz(9)-temp(3,j)*zln(2,jj)-temp(9,j)*zln(8,jj) &
                     -temp(15,j)*zln(14,jj)-temp(21,j)*zln(20,jj) &
                     -temp(27,j)*zln(26,jj)-temp(33,j)*zln(32,jj)
              zz(10)=zz(10)-temp(4,j)*zln(2,jj)-temp(10,j)*zln(8,jj) &
                     -temp(16,j)*zln(14,jj)-temp(22,j)*zln(20,jj) &
                     -temp(28,j)*zln(26,jj)-temp(34,j)*zln(32,jj)
              zz(11)=zz(11)-temp(5,j)*zln(2,jj)-temp(11,j)*zln(8,jj) &
                     -temp(17,j)*zln(14,jj)-temp(23,j)*zln(20,jj) &
                     -temp(29,j)*zln(26,jj)-temp(35,j)*zln(32,jj)
              zz(12)=zz(12)-temp(6,j)*zln(2,jj)-temp(12,j)*zln(8,jj) &
                     -temp(18,j)*zln(14,jj)-temp(24,j)*zln(20,jj) &
                     -temp(30,j)*zln(26,jj)-temp(36,j)*zln(32,jj)
              zz(13)=zz(13)-temp(1,j)*zln(3,jj)-temp(7,j)*zln(9,jj) &
                     -temp(13,j)*zln(15,jj)-temp(19,j)*zln(21,jj) &
                     -temp(25,j)*zln(27,jj)-temp(31,j)*zln(33,jj)
              zz(14)=zz(14)-temp(2,j)*zln(3,jj)-temp(8,j)*zln(9,jj) &
                     -temp(14,j)*zln(15,jj)-temp(20,j)*zln(21,jj) &
                     -temp(26,j)*zln(27,jj)-temp(32,j)*zln(33,jj)
              zz(15)=zz(15)-temp(3,j)*zln(3,jj)-temp(9,j)*zln(9,jj) &
                     -temp(15,j)*zln(15,jj)-temp(21,j)*zln(21,jj) &
                     -temp(27,j)*zln(27,jj)-temp(33,j)*zln(33,jj)
              zz(16)=zz(16)-temp(4,j)*zln(3,jj)-temp(10,j)*zln(9,jj) &
                     -temp(16,j)*zln(15,jj)-temp(22,j)*zln(21,jj) &
                     -temp(28,j)*zln(27,jj)-temp(34,j)*zln(33,jj)
              zz(17)=zz(17)-temp(5,j)*zln(3,jj)-temp(11,j)*zln(9,jj) &
                     -temp(17,j)*zln(15,jj)-temp(23,j)*zln(21,jj) &
                     -temp(29,j)*zln(27,jj)-temp(35,j)*zln(33,jj)
              zz(18)=zz(18)-temp(6,j)*zln(3,jj)-temp(12,j)*zln(9,jj) &
                     -temp(18,j)*zln(15,jj)-temp(24,j)*zln(21,jj) &
                     -temp(30,j)*zln(27,jj)-temp(36,j)*zln(33,jj)
              zz(19)=zz(19)-temp(1,j)*zln(4,jj)-temp(7,j)*zln(10,jj) &
                     -temp(13,j)*zln(16,jj)-temp(19,j)*zln(22,jj) &
                     -temp(25,j)*zln(28,jj)-temp(31,j)*zln(34,jj)
              zz(20)=zz(20)-temp(2,j)*zln(4,jj)-temp(8,j)*zln(10,jj) &
                     -temp(14,j)*zln(16,jj)-temp(20,j)*zln(22,jj) &
                     -temp(26,j)*zln(28,jj)-temp(32,j)*zln(34,jj)
              zz(21)=zz(21)-temp(3,j)*zln(4,jj)-temp(9,j)*zln(10,jj) &
                     -temp(15,j)*zln(16,jj)-temp(21,j)*zln(22,jj) &
                     -temp(27,j)*zln(28,jj)-temp(33,j)*zln(34,jj)
              zz(22)=zz(22)-temp(4,j)*zln(4,jj)-temp(10,j)*zln(10,jj) &
                     -temp(16,j)*zln(16,jj)-temp(22,j)*zln(22,jj) &
                     -temp(28,j)*zln(28,jj)-temp(34,j)*zln(34,jj)
              zz(23)=zz(23)-temp(5,j)*zln(4,jj)-temp(11,j)*zln(10,jj) &
                     -temp(17,j)*zln(16,jj)-temp(23,j)*zln(22,jj) &
                     -temp(29,j)*zln(28,jj)-temp(35,j)*zln(34,jj)
              zz(24)=zz(24)-temp(6,j)*zln(4,jj)-temp(12,j)*zln(10,jj) &
                     -temp(18,j)*zln(16,jj)-temp(24,j)*zln(22,jj) &
                     -temp(30,j)*zln(28,jj)-temp(36,j)*zln(34,jj)
              zz(25)=zz(25)-temp(1,j)*zln(5,jj)-temp(7,j)*zln(11,jj) &
                     -temp(13,j)*zln(17,jj)-temp(19,j)*zln(23,jj) &
                     -temp(25,j)*zln(29,jj)-temp(31,j)*zln(35,jj)
              zz(26)=zz(26)-temp(2,j)*zln(5,jj)-temp(8,j)*zln(11,jj) &
                     -temp(14,j)*zln(17,jj)-temp(20,j)*zln(23,jj) &
                     -temp(26,j)*zln(29,jj)-temp(32,j)*zln(35,jj)
              zz(27)=zz(27)-temp(3,j)*zln(5,jj)-temp(9,j)*zln(11,jj) &
                     -temp(15,j)*zln(17,jj)-temp(21,j)*zln(23,jj) &
                     -temp(27,j)*zln(29,jj)-temp(33,j)*zln(35,jj)
              zz(28)=zz(28)-temp(4,j)*zln(5,jj)-temp(10,j)*zln(11,jj) &
                     -temp(16,j)*zln(17,jj)-temp(22,j)*zln(23,jj) &
                     -temp(28,j)*zln(29,jj)-temp(34,j)*zln(35,jj)
              zz(29)=zz(29)-temp(5,j)*zln(5,jj)-temp(11,j)*zln(11,jj) &
                     -temp(17,j)*zln(17,jj)-temp(23,j)*zln(23,jj) &
                     -temp(29,j)*zln(29,jj)-temp(35,j)*zln(35,jj)
              zz(30)=zz(30)-temp(6,j)*zln(5,jj)-temp(12,j)*zln(11,jj) &
                     -temp(18,j)*zln(17,jj)-temp(24,j)*zln(23,jj) &
                     -temp(30,j)*zln(29,jj)-temp(36,j)*zln(35,jj)
              zz(31)=zz(31)-temp(1,j)*zln(6,jj)-temp(7,j)*zln(12,jj) &
                     -temp(13,j)*zln(18,jj)-temp(19,j)*zln(24,jj) &
                     -temp(25,j)*zln(30,jj)-temp(31,j)*zln(36,jj)
              zz(32)=zz(32)-temp(2,j)*zln(6,jj)-temp(8,j)*zln(12,jj) &
                     -temp(14,j)*zln(18,jj)-temp(20,j)*zln(24,jj) &
                     -temp(26,j)*zln(30,jj)-temp(32,j)*zln(36,jj)
              zz(33)=zz(33)-temp(3,j)*zln(6,jj)-temp(9,j)*zln(12,jj) &
                     -temp(15,j)*zln(18,jj)-temp(21,j)*zln(24,jj) &
                     -temp(27,j)*zln(30,jj)-temp(33,j)*zln(36,jj)
              zz(34)=zz(34)-temp(4,j)*zln(6,jj)-temp(10,j)*zln(12,jj) &
                     -temp(16,j)*zln(18,jj)-temp(22,j)*zln(24,jj) &
                     -temp(28,j)*zln(30,jj)-temp(34,j)*zln(36,jj)
              zz(35)=zz(35)-temp(5,j)*zln(6,jj)-temp(11,j)*zln(12,jj) &
                     -temp(17,j)*zln(18,jj)-temp(23,j)*zln(24,jj) &
                     -temp(29,j)*zln(30,jj)-temp(35,j)*zln(36,jj)
              zz(36)=zz(36)-temp(6,j)*zln(6,jj)-temp(12,j)*zln(12,jj) &
                     -temp(18,j)*zln(18,jj)-temp(24,j)*zln(24,jj) &
                     -temp(30,j)*zln(30,jj)-temp(36,j)*zln(36,jj)
               endif
  310       continue
            call inv66(zln(1,k),zz,diag(1,jc))
!c$dir max_trips(9)
            do 220 l=1,36
               temp(l,jc)=zz(l)
  220       continue
            t(1)=t(1)+zz(1)*zln(1,k)+zz(7)*zln(7,k) &
                       +zz(13)*zln(13,k)+zz(19)*zln(19,k) &
                       +zz(25)*zln(25,k)+zz(31)*zln(31,k)
            t(2)=t(2)+zz(1)*zln(2,k)+zz(7)*zln(8,k) &
                       +zz(13)*zln(14,k)+zz(19)*zln(20,k) &
                       +zz(25)*zln(26,k)+zz(31)*zln(32,k)
            t(3)=t(3)+zz(2)*zln(2,k)+zz(8)*zln(8,k) &
                       +zz(14)*zln(14,k)+zz(20)*zln(20,k) &
                       +zz(26)*zln(26,k)+zz(32)*zln(32,k)
            t(4)=t(4)+zz(1)*zln(3,k)+zz(7)*zln(9,k) &
                       +zz(13)*zln(15,k)+zz(19)*zln(21,k) &
                       +zz(25)*zln(27,k)+zz(31)*zln(33,k)
            t(5)=t(5)+zz(2)*zln(3,k)+zz(8)*zln(9,k) &
                       +zz(14)*zln(15,k)+zz(20)*zln(21,k) &
                       +zz(26)*zln(27,k)+zz(32)*zln(33,k)
            t(6)=t(6)+zz(3)*zln(3,k)+zz(9)*zln(9,k) &
                       +zz(15)*zln(15,k)+zz(21)*zln(21,k) &
                       +zz(27)*zln(27,k)+zz(33)*zln(33,k)
            t(7)=t(7)+zz(1)*zln(4,k)+zz(7)*zln(10,k) &
                       +zz(13)*zln(16,k)+zz(19)*zln(22,k) &
                       +zz(25)*zln(28,k)+zz(31)*zln(34,k)
            t(8)=t(8)+zz(2)*zln(4,k)+zz(8)*zln(10,k) &
                       +zz(14)*zln(16,k)+zz(20)*zln(22,k) &
                       +zz(26)*zln(28,k)+zz(32)*zln(34,k)
            t(9)=t(9)+zz(3)*zln(4,k)+zz(9)*zln(10,k) &
                       +zz(15)*zln(16,k)+zz(21)*zln(22,k) &
                       +zz(27)*zln(28,k)+zz(33)*zln(34,k)
            t(10)=t(10)+zz(4)*zln(4,k)+zz(10)*zln(10,k) &
                       +zz(16)*zln(16,k)+zz(22)*zln(22,k) &
                       +zz(28)*zln(28,k)+zz(34)*zln(34,k)
            t(11)=t(11)+zz(1)*zln(5,k)+zz(7)*zln(11,k) &
                       +zz(13)*zln(17,k)+zz(19)*zln(23,k) &
                       +zz(25)*zln(29,k)+zz(31)*zln(35,k)
            t(12)=t(12)+zz(2)*zln(5,k)+zz(8)*zln(11,k) &
                       +zz(14)*zln(17,k)+zz(20)*zln(23,k) &
                       +zz(26)*zln(29,k)+zz(32)*zln(35,k)
            t(13)=t(13)+zz(3)*zln(5,k)+zz(9)*zln(11,k) &
                       +zz(15)*zln(17,k)+zz(21)*zln(23,k) &
                       +zz(27)*zln(29,k)+zz(33)*zln(35,k)
            t(14)=t(14)+zz(4)*zln(5,k)+zz(10)*zln(11,k) &
                       +zz(16)*zln(17,k)+zz(22)*zln(23,k) &
                       +zz(28)*zln(29,k)+zz(34)*zln(35,k)
            t(15)=t(15)+zz(5)*zln(5,k)+zz(11)*zln(11,k)&
                       +zz(17)*zln(17,k)+zz(23)*zln(23,k) &
                       +zz(29)*zln(29,k)+zz(35)*zln(35,k)
            t(16)=t(16)+zz(1)*zln(6,k)+zz(7)*zln(12,k) &
                       +zz(13)*zln(18,k)+zz(19)*zln(24,k) &
                       +zz(25)*zln(30,k)+zz(31)*zln(36,k)
            t(17)=t(17)+zz(2)*zln(6,k)+zz(8)*zln(12,k) &
                       +zz(14)*zln(18,k)+zz(20)*zln(24,k) &
                       +zz(26)*zln(30,k)+zz(32)*zln(36,k)
            t(18)=t(18)+zz(3)*zln(6,k)+zz(9)*zln(12,k) & 
                       +zz(15)*zln(18,k)+zz(21)*zln(24,k) &
                       +zz(27)*zln(30,k)+zz(33)*zln(36,k)
            t(19)=t(19)+zz(4)*zln(6,k)+zz(10)*zln(12,k) &
                       +zz(16)*zln(18,k)+zz(22)*zln(24,k) &
                       +zz(28)*zln(30,k)+zz(34)*zln(36,k)
            t(20)=t(20)+zz(5)*zln(6,k)+zz(11)*zln(12,k) &
                       +zz(17)*zln(18,k)+zz(23)*zln(24,k) &
                       +zz(29)*zln(30,k)+zz(35)*zln(36,k)
            t(21)=t(21)+zz(6)*zln(6,k)+zz(12)*zln(12,k) &
                       +zz(18)*zln(18,k)+zz(24)*zln(24,k) &
                       +zz(30)*zln(30,k)+zz(36)*zln(36,k)
  200    continue
!c$dir max_trips(6)
         do 320 l=1,21
            diag(l,ic)=diag(l,ic)-t(l)
  320    continue
         call inv6(diag(1,ic),ir)
         nch(ic)=-1
         kk=par(ic)
         nch(kk)=nch(kk)-1
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s6um1(ic,xlnzr,colno,zln,diag,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),par(*)
      dimension zln(9,*),diag(6,*),temp(9,*),indx(*)
      dimension s(9),zz(9)
      common /mchdpn/ rmax,rmin,epsm,lratio
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
!c$dir max_trip(9)
         do 100 l=1,9
            s(l)=0.0d0
  100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(4,j)*zln(4,jj) &
                             +temp(7,j)*zln(7,jj)
                  s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(5,j)*zln(4,jj) &
                             +temp(8,j)*zln(7,jj)
                  s(3)=s(3)+temp(3,j)*zln(1,jj)+temp(6,j)*zln(4,jj) &
                             +temp(9,j)*zln(7,jj)
                  s(4)=s(4)+temp(1,j)*zln(2,jj)+temp(4,j)*zln(5,jj) &
                             +temp(7,j)*zln(8,jj)
                  s(5)=s(5)+temp(2,j)*zln(2,jj)+temp(5,j)*zln(5,jj) &
                             +temp(8,j)*zln(8,jj)
                  s(6)=s(6)+temp(3,j)*zln(2,jj)+temp(6,j)*zln(5,jj) &
                             +temp(9,j)*zln(8,jj)
                  s(7)=s(7)+temp(1,j)*zln(3,jj)+temp(4,j)*zln(6,jj) &
                             +temp(7,j)*zln(9,jj)
                  s(8)=s(8)+temp(2,j)*zln(3,jj)+temp(5,j)*zln(6,jj) &
                             +temp(8,j)*zln(9,jj)
                  s(9)=s(9)+temp(3,j)*zln(3,jj)+temp(6,j)*zln(6,jj) &
                             +temp(9,j)*zln(9,jj)
               endif
  310       continue
!c$dir max_trip(9)
            do 320 l=1,9
              temp(l,jc)=zln(l,k)-s(l)
              zln(l,k)=temp(l,jc)
              s(l)=0.0d0
  320       continue
  200    continue
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s6um2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*)
      dimension zln(36,*),diag(21,*),temp(36,neqns),dsln(36,*),indx(*)
      joc=0
      do 100 ic=nstop,neqns
         do 105 m=1,36
            do 106 jj=1,nstop
               temp(jj,m)=0.0d0
  106    continue
  105    continue
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)-1
         do 110 k=ks,ke
            jj=colno(k)
               temp(:,jj)=zln(:,k)
               indx(jj)=ic
  110    continue
         do 111 k=ks,ke
            jj=colno(k)
            call inv66(zln(1,k),temp,diag(1,jj))
  111    continue
!c
         do 112 k=ks,ke
            jj=colno(k)
            diag(1,ic)=diag(1,ic)-temp(jj,1)*zln(1,k) &
                            -temp(jj,4)*zln(4,k)-temp(jj,7)*zln(7,k)
            diag(2,ic)=diag(2,ic)-temp(jj,1)*zln(2,k) &
                            -temp(jj,4)*zln(5,k)-temp(jj,7)*zln(8,k)
            diag(3,ic)=diag(3,ic)-temp(jj,2)*zln(2,k) &
                            -temp(jj,5)*zln(5,k)-temp(jj,8)*zln(8,k)
            diag(4,ic)=diag(4,ic)-temp(jj,1)*zln(3,k) &
                            -temp(jj,4)*zln(6,k)-temp(jj,7)*zln(9,k)
            diag(5,ic)=diag(5,ic)-temp(jj,2)*zln(3,k) &
                            -temp(jj,5)*zln(6,k)-temp(jj,8)*zln(9,k)
            diag(6,ic)=diag(6,ic)-temp(jj,3)*zln(3,k) &
                            -temp(jj,6)*zln(6,k)-temp(jj,9)*zln(9,k)
  112    continue
         do 120 jc=nstop,ic-1
            joc=joc+1
            j1=xlnzr(jc)
            j2=xlnzr(jc+1)
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  dsln(1,joc)=dsln(1,joc)-temp(j,1)*zln(1,jj) &
                             -temp(j,4)*zln(4,jj)-temp(j,7)*zln(7,jj)
                  dsln(2,joc)=dsln(2,joc)-temp(j,2)*zln(1,jj) &
                             -temp(j,5)*zln(4,jj)-temp(j,8)*zln(7,jj)
                  dsln(3,joc)=dsln(3,joc)-temp(j,3)*zln(1,jj) &
                             -temp(j,6)*zln(4,jj)-temp(j,9)*zln(7,jj)
                  dsln(4,joc)=dsln(4,joc)-temp(j,1)*zln(2,jj) &
                             -temp(j,4)*zln(5,jj)-temp(j,7)*zln(8,jj)
                  dsln(5,joc)=dsln(5,joc)-temp(j,2)*zln(2,jj) &
                             -temp(j,5)*zln(5,jj)-temp(j,8)*zln(8,jj)
                  dsln(6,joc)=dsln(6,joc)-temp(j,3)*zln(2,jj) &
                             -temp(j,6)*zln(5,jj)-temp(j,9)*zln(8,jj)
                  dsln(7,joc)=dsln(7,joc)-temp(j,1)*zln(3,jj) &
                             -temp(j,4)*zln(6,jj)-temp(j,7)*zln(9,jj)
                  dsln(8,joc)=dsln(8,joc)-temp(j,2)*zln(3,jj) &
                             -temp(j,5)*zln(6,jj)-temp(j,8)*zln(9,jj)
                  dsln(9,joc)=dsln(9,joc)-temp(j,3)*zln(3,jj) &
                             -temp(j,6)*zln(6,jj)-temp(j,9)*zln(9,jj)
               endif
  220       continue
  120    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine s6um3(n,dsln,diag,indx,temp)
      implicit double precision(a-h,o-z)
      dimension dsln(9,*),diag(6,*),indx(*),temp(9,*),t(9)
      if(n.le.0) goto 1000
      indx(1)=0
      joc=1
      call inv3(diag(1,1),ir)
      do 100 i=2,n
         indx(i)=joc
         do 110 j=1,i-1
            call d3dot(t,dsln(1,indx(i)),dsln(1,indx(j)),j-1)
!c$dir max_trips(9)
            do 111 l=1,9
               dsln(l,joc)=dsln(l,joc)-t(l)
  111       continue
            joc=joc+1
  110    continue
         call v3prod(dsln(1,indx(i)),diag,temp,i-1)
         call d3dotl(t,temp,dsln(1,indx(i)),i-1)
!c$dir max_trips(6)
         do 112 l=1,6
            diag(l,i)=diag(l,i)-t(l)
  112    continue
         call vcopy(temp,dsln(1,indx(i)),9*(i-1))
         call inv3(diag(1,i),ir)
  100 continue
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      double precision function spdot2(b,zln,colno,ks,ke)
      implicit double precision(a-h,o-z)
      integer colno(*)
      double precision zln(*),b(*)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      s=0.0d0
      do 100 jj=ks,ke
         j=colno(jj)
         s=s+zln(jj)*b(j)
  100 continue
      spdot2=s
      end function
!======================================================================!
!                                                                      !
!======================================================================!


      subroutine staij1(isw,i,j,aij,ir)
!
!
!     USE sp_direct_solver
      implicit double precision(a-h,o-z)
!
!----------------------------------------------------------------------
!
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
!
!----------------------------------------------------------------------
!
      !dimension iv(*),aij(ndeg*ndeg)
      dimension aij(ndeg*ndeg)
      common /debug/ idbg
      common /mchdpn/ rmax,rmin,epsm,lratio
!      DIMENSION diag(neqns*ndeg2l), zln(iv(21)*ndeg2)
!      DIMENSION dsln(len_dsln*ndeg2)
!
      ir=0
         ndeg2=ndeg*ndeg
         ndeg2l=ndeg*(ndeg+1)/2
         if(stage.eq.30) write(6,*) 'warning a matrix was build up '&
                                  //'but never solved.'
      if(stage.eq.10) then
         ALLOCATE(diag(neqns*ndeg2l),STAT = ierror)
         raloc = raloc + neqns*ndeg2l
         ALLOCATE(zln(len_colno*ndeg2),STAT = ierror)
!!!@         CALL hecmw_ALLOCATE(51,'zln','dbl',1,1,len_colno*ndeg2,
!!!@     &                       stat=ierror)
!!!@         IF(ierror.NE.0) STOP "Allocation error zln"
!        iv(27)=iv(10)+len_dsln*ndeg2*lratio
!rmiv
!         iv(27)=iv(8)
         raloc = raloc + len_colno*ndeg2
         ALLOCATE(dsln(len_dsln*ndeg2),STAT = ierror)
!         CALL hecmw_ALLOCATE(imsg,'dsln','dbl',1,1,len_dsln*ndeg2,
!     &                       stat=ierror)
         IF(ierror.NE.0) STOP "Allocation error dsln"
         !CALL hecmw_ALLOCATE(0,'dsln','dbl',1,1,len_dsln*ndeg2,0,0,0,0)
         raloc = raloc + len_dsln*ndeg2
      endif
      if(stage.ne.20) then
!
! for diagonal
! 
          diag = 0.
!         call clear(neqns*ndeg2,diag)
!
! for lower triangle
!
         zln = 0.
!         call clear(iv(21)*ndeg2,zln)
!
! for dense window
!
         dsln = 0.
!         call clear(len_dsln*ndeg2,dsln)
         stage=20
      endif
!         Print *,'********Set Stage 20 *********'
!
      if(ndeg.le.2) then
         call addr0(isw,i,j,aij,invp,xlnzr,colno, &
                 diag,zln,dsln,nstop,ndeg2,ndeg2l,ir)
      elseif(ndeg.eq.3) then
         call addr3(isw,i,j,aij,invp,xlnzr,colno, &
                 diag,zln,dsln,nstop,ir)
      else
         call addrx(isw,i,j,aij,invp,xlnzr,colno,diag, &
                    zln,dsln,nstop,ndeg,ndeg2,ndeg2l,ir)
      endif
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine stiaja(neqns,ia,ja,jcpt,jcolno)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer ia(*),ja(*),jcpt(*),jcolno(*)
      common /debug/idbg
!c
      ia(1)=1
      l=0
      do 100 k=1,neqns
         joc=jcpt(k)
  110    continue
         if(joc.eq.0) goto 120
         ii=jcolno(joc)
         if(ii.eq.k) goto 130
         l=l+1
         ja(l)=ii
  130    continue
         joc=jcpt(joc)
         goto 110
  120    ia(k+1)=l+1
  100 continue
      if(idbg.ne.0) then
         write(6,*) 'ia '
         write(6,60) (ia(i),i=1,neqns)
         write(6,*) 'ja '
         write(6,60) (ja(i),i=1,ia(neqns+1))
      endif
   60 format(10i7)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine stsmat(neqns,nttbr,irow,jcol,jcpt,jcolno)
!c
!c     coded by t.arakawa of RIST on 040510
!c
      integer irow(*),jcol(*),jcpt(*),jcolno(*)
      common /debug/idbg
!c
      do 10 i=1,2*nttbr
         jcpt(i)=0
         jcolno(i)=0
   10 continue
      do 20 i=1,neqns
         jcpt(i)=i+neqns
         jcolno(i+neqns)=i
   20 continue
!c
      k=2*neqns
!c     do 100 l=neqns+1,neqns+nttbr
      do 100 l=1,nttbr
         i=irow(l)
         j=jcol(l)
         if(i.eq.j) goto 100
         joc=jcpt(i)
         locr=i
  110    continue
         if(joc.eq.0) goto 120
         if(jcolno(joc).eq.j) then
            goto 100
         elseif(jcolno(joc).gt.j) then
            goto 130
         endif
         locr=joc
         joc=jcpt(joc)
         goto 110
  120    continue
         k=k+1
         jcpt(locr)=k
         jcolno(k)=j
         goto 150
  130    continue
         k=k+1
         jcpt(locr)=k
         jcpt(k)=joc
         jcolno(k)=j
  150    continue
         joc=jcpt(j)
         locr=j
  160    continue
         if(joc.eq.0) goto 170
         if(jcolno(joc).eq.i) then
            goto 100
         elseif(jcolno(joc).gt.i) then
            goto 180
         endif
         locr=joc
         joc=jcpt(joc)
         goto 160
  170    continue
         k=k+1
         jcpt(locr)=k
         jcolno(k)=i
         goto 100
  180    continue
         k=k+1
         jcpt(locr)=k
         jcpt(k)=joc
         jcolno(k)=i
  100 continue
      if(idbg.ne.0) then
         write(6,*) 'jcolno'
         write(6,60) (jcolno(i),i=1,k)
         write(6,*) 'jcpt'
         write(6,60) (jcpt(i),i=1,k)
   60 format(10i7)
      endif
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sum(ic,xlnzr,colno,zln,diag,nch,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),nch(*),par(*)
      dimension zln(*),diag(*),temp(*),indx(*)
      common /mchdpn/ rmax,rmin,epsm,lratio
      common isem
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         t=0.0d0
!c        do 100 i=1,ic
!c           temp(i)=0.0d0
!c 100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            s=0.0d0
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  s=s+temp(j)*zln(jj)
               endif
  310       continue
!c           j1=xlnzr(jc)
!c           jj=xlnzr(jc+1)-j1
!c           ss=ddoti(jj,zln(j1),colno(j1),temp)
!c           zz=zln(k)-ddoti(jj,zln(j1),colno(j1),temp)
            zz=zln(k)-s
            zln(k)=zz*diag(jc)
            temp(jc)=zz
            t=t+zz*zln(k)
  200    continue
         piv=diag(ic)-t
         if(dabs(piv).gt.rmin) then
            diag(ic)=1.0d0/piv
         endif
    1 continue
         if(isem.eq.1) then
            isem=0
            nch(ic)=-1
            kk=par(ic)
            nch(kk)=nch(kk)-1
            isem=1
         else
            goto 1
         endif
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sum1(ic,xlnzr,colno,zln,diag,par,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),par(*)
      dimension zln(*),diag(*),temp(*),indx(*)
      common /mchdpn/ rmax,rmin,epsm,lratio
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         t=0.0d0
!c        do 100 i=1,ic
!c           temp(i)=0.0d0
!c 100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            s=0.0d0
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  s=s+temp(j)*zln(jj)
               endif
  310       continue
            zz=zln(k)-s
!c           j1=xlnzr(jc)
!c           jj=xlnzr(jc+1)-j1
!c           zz=zln(k)-ddoti(jj,zln(j1),colno(j1),temp)
            zln(k)=zz
            temp(jc)=zz
!c           t=t+zz*zz*diag(jc)
  200    continue
!c        diag(ic)=diag(ic)-t
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*)
      dimension zln(*),diag(*),temp(*),indx(*),dsln(*)
      joc=0
      do 100 ic=nstop,neqns
         do 105 i=1,nstop
            temp(i)=0.0d0
  105    continue
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)-1
         do 110 k=ks,ke
            jj=colno(k)
            temp(jj)=zln(k)
            zln(k)=temp(jj)*diag(jj)
            indx(jj)=ic
            diag(ic)=diag(ic)-temp(jj)*zln(k)
  110    continue
         do 120 jc=nstop,ic-1
            s=0.0d0
            joc=joc+1
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  s=s+temp(j)*zln(jj)
               endif
  220       continue
            if(s.eq.0.0) write(16,*) ic,jc
            dsln(joc)=dsln(joc)-s
!c           j1=xlnzr(jc)
!c           jj=xlnzr(jc+1)-j1
!c           dsln(joc)=dsln(joc)-ddoti(jj,zln(j1),colno(j1),temp)
  120    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sum3(n,dsln,diag,indx,temp)
      implicit double precision(a-h,o-z)
      dimension dsln(*),diag(*),indx(*),temp(*)
      if(n.le.0) goto 1000
      indx(1)=0
      joc=1
      diag(1)=1.0d0/diag(1)
      do 100 i=2,n
         indx(i)=joc
         do 110 j=1,i-1
            dsln(joc)=dsln(joc)-ddot(dsln(indx(i)),dsln(indx(j)),j-1)
            joc=joc+1
  110    continue
         call vprod(dsln(indx(i)),diag,temp,i-1)
         diag(i)=diag(i)-ddot(temp,dsln(indx(i)),i-1)
         call vcopy(temp,dsln(indx(i)),i-1)
         diag(i)=1.0d0/diag(i)
  100 continue
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sxpdot(ndeg,bi,b,zln,colno,ks,ke)
      implicit double precision(a-h,o-z)
      integer colno(*)
      double precision zln(ndeg,ndeg,*),b(ndeg,*),bi(ndeg)
!c
!c----------------------------------------------------------------------
!c
!c      spdot1 performs inner product of sparse vectors
!c
!c
!c      #coded by t.arakawa of RIST on 040510
!c
!c----------------------------------------------------------------------
!c
      do 100 jj=ks,ke
         j=colno(jj)
         do 101 m=1,ndeg
         do 102 n=1,ndeg
            bi(n)=bi(n)-zln(n,m,jj)*b(m,j)
  102 continue
  101 continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sxum(ic,xlnzr,colno,zln,diag,nch,par,temp,indx, &
                      ndeg,ndegl,zz,t)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),nch(*),par(*)
      dimension zln(ndeg,ndeg,*),diag(ndegl,*),temp(ndeg,ndeg,*),indx(*)
      dimension zz(ndeg,ndeg),t(ndegl)
      common /mchdpn/ rmax,rmin,epsm,lratio
      ndeg22=ndeg*ndeg
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         t=0.0
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            zz=zln(:,:,k)
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  do 311 m=1,ndeg,2
                  do 312 n=1,ndeg,2
                  do 313 kk=1,ndeg,2
                     zz(n,m)=zz(n,m)-temp(n,kk,j)*zln(m,kk,jj) &
                                   -temp(n,kk+1,j)*zln(m,kk+1,jj)
                     zz(n,m+1)=zz(n,m+1)-temp(n,kk,j)*zln(m+1,kk,jj)&
                                       -temp(n,kk+1,j)*zln(m+1,kk+1,jj)
                     zz(n+1,m)=zz(n+1,m)-temp(n+1,kk,j)*zln(m,kk,jj) &
                                       -temp(n+1,kk+1,j)*zln(m,kk+1,jj)
                   zz(n+1,m+1)=zz(n+1,m+1)-temp(n+1,kk,j)*zln(m+1,kk,jj)&
                                     -temp(n+1,kk+1,j)*zln(m+1,kk+1,jj)
  313             continue
  312             continue
  311             continue
               endif
  310       continue
            call invxx(zln(1,1,k),zz,diag(1,jc),ndeg)
!c           do 220 m=1,ndeg
!c           do 220 n=1,ndeg
!c              temp(n,m,jc)=zz(n,m)
!c 220       continue
            temp(:,:,jc)=zz
            joc=0
            do 221 n=1,ndeg
            do 222 m=1,n
            joc=joc+1
            do 223 kk=1,ndeg,2
               t(joc)=t(joc) &
              +zz(n,kk)*zln(m,kk,k)+zz(n,kk+1)*zln(m,kk+1,k)
  223       continue
  222       continue
  221       continue
  200    continue
!c        do 320 l=1,ndegl
!c           diag(l,ic)=diag(l,ic)-t(l)
!c 320    continue
         diag(:,ic)=diag(:,ic)-t
         call invx(diag(1,ic),ndeg,ir)
         nch(ic)=-1
         kk=par(ic)
         nch(kk)=nch(kk)-1
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sxum1(ic,xlnzr,colno,zln,diag,par,temp,indx, &
                      ndeg,ndegl,s)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*),par(*)
      dimension zln(ndeg,ndeg,*),diag(ndegl,*),temp(ndeg,ndeg,*),indx(*)
      dimension s(ndeg,ndeg)
      common /mchdpn/ rmax,rmin,epsm,lratio
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)
         do 100 m=1,ndeg
         do 101 n=1,ndeg
            s(n,m)=0.0d0
  101    continue
  100    continue
         do 200 k=ks,ke-1
            jc=colno(k)
            indx(jc)=ic
            do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  do 311 m=1,ndeg
                  do 312 n=1,ndeg
                  do 313 kk=1,ndeg
                     s(n,m)=s(n,m)+temp(n,kk,j)*zln(m,kk,jj)
  313             continue
  312             continue
  311             continue
               endif
  310       continue
            do 320 m=1,ndeg
            do 321 n=1,ndeg
              temp(n,m,jc)=zln(n,m,k)-s(n,m)
              zln(n,m,k)=temp(n,m,jc)
              s(n,m)=0.0d0
  321       continue
  320       continue
  200    continue
         return
         end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sxum2(neqns,nstop,xlnzr,colno,zln,diag,dsln,temp,indx, &
                     ndeg,ndegl)
      implicit double precision(a-h,o-z)
      integer xlnzr(*),colno(*)
      dimension zln(ndeg,ndeg,*),diag(ndegl,*),temp(ndeg,ndeg,*), &
               dsln(ndeg,ndeg,*),indx(*)
      joc=0
      do 100 ic=nstop,neqns
         ks=xlnzr(ic)
         ke=xlnzr(ic+1)-1
         do 110 k=ks,ke
            jj=colno(k)
            do 111 m=1,ndeg
            do 112 n=1,ndeg
               temp(n,m,jj)=zln(n,m,k)
               indx(jj)=ic
  112    continue
  111    continue
  110    continue
         do 120 k=ks,ke
            jj=colno(k)
            call invxx(zln(1,1,k),temp(1,1,jj),diag(1,jj),ndeg)
  120    continue
!c
            locd=0
            do 130 n=1,ndeg
            do 131 m=1,n
            locd=locd+1
            do 132 k=ks,ke
            jj=colno(k)
            do 133 kk=1,ndeg
               diag(locd,ic)=diag(locd,ic)-temp(n,kk,jj)*zln(m,kk,k)
  133       continue
  132       continue
  131       continue
  130       continue
         do 200 jc=nstop,ic-1
            joc=joc+1
            j1=xlnzr(jc)
            j2=xlnzr(jc+1)
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
               j=colno(jj)
               if(indx(j).eq.ic) then
                  do 221 m=1,ndeg
                  do 222 n=1,ndeg
                  do 223 k=1,ndeg
                  dsln(n,m,joc)=dsln(n,m,joc)-temp(n,k,j)*zln(m,k,jj)
  223             continue
  222             continue
  221             continue
               endif
  220       continue
  200    continue
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine sxum3(nn,dsln,diag,indx,temp, &
                     ndeg,ndegl,t)
      implicit double precision(a-h,o-z)
      dimension dsln(ndeg,ndeg,*),diag(ndegl,*),indx(*), &
                temp(ndeg,ndeg,*),t(ndeg,ndeg)
      if(nn.le.0) goto 1000
      indx(1)=0
      joc=1
      call invx(diag(1,1),ndeg,ir)
      do 100 i=2,nn
         indx(i)=joc
         do 110 j=1,i-1
            call dxdot(ndeg,t,dsln(1,1,indx(i)), &
                       dsln(1,1,indx(j)),j-1)
            do 111 m=1,ndeg
            do 112 n=1,ndeg
               dsln(n,m,joc)=dsln(n,m,joc)-t(n,m)
  112       continue
  111       continue
            joc=joc+1
  110    continue
         call vxprod(ndeg,ndegl,dsln(1,1,indx(i)),diag,temp,i-1)
         call dxdotl(ndeg,t,temp,dsln(1,1,indx(i)),i-1)
         locd=0
         do 221 n=1,ndeg
         do 222 m=1,n
         locd=locd+1
            diag(locd,i)=diag(locd,i)-t(n,m)
  222    continue
  221    continue
         call vcopy(temp,dsln(1,1,indx(i)),ndeg*ndeg*(i-1))
         call invx(diag(1,i),ndeg,ir)
  100 continue
 1000 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine v2prod(a,b,c,n)
      implicit double precision(a-h,o-z)
      dimension a(4,n),b(3,n),c(4,n)
      do 100 i=1,n
!c        call inv22(c(1,i),a(1,i),b(1,i))
      c(3,i)=a(3,i)-a(1,i)*b(2,i)
      c(1,i)=a(1,i)*b(1,i)
      c(3,i)=c(3,i)*b(3,i)
      c(1,i)=c(1,i)-c(3,i)*b(2,i)
!c
      c(4,i)=a(4,i)-a(2,i)*b(2,i)
      c(2,i)=a(2,i)*b(1,i)
      c(4,i)=c(4,i)*b(3,i)
      c(2,i)=c(2,i)-c(4,i)*b(2,i)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine v3prod(zln,diag,zz,n)
      implicit double precision(a-h,o-z)
      dimension zln(9,n),diag(6,n),zz(9,n)
      do 100 i=1,n
         zz(4,i)=zln(4,i)-zln(1,i)*diag(2,i)
         zz(7,i)=zln(7,i)-zln(1,i)*diag(4,i)-zz(4,i)*diag(5,i)
         zz(1,i)=zln(1,i)*diag(1,i)
         zz(4,i)=zz(4,i)*diag(3,i)
         zz(7,i)=zz(7,i)*diag(6,i)
         zz(4,i)=zz(4,i)-zz(7,i)*diag(5,i)
         zz(1,i)=zz(1,i)-zz(4,i)*diag(2,i)-zz(7,i)*diag(4,i)
!c
         zz(5,i)=zln(5,i)-zln(2,i)*diag(2,i)
         zz(8,i)=zln(8,i)-zln(2,i)*diag(4,i)-zz(5,i)*diag(5,i)
         zz(2,i)=zln(2,i)*diag(1,i)
         zz(5,i)=zz(5,i)*diag(3,i)
         zz(8,i)=zz(8,i)*diag(6,i)
         zz(5,i)=zz(5,i)-zz(8,i)*diag(5,i)
         zz(2,i)=zz(2,i)-zz(5,i)*diag(2,i)-zz(8,i)*diag(4,i)
!c
         zz(6,i)=zln(6,i)-zln(3,i)*diag(2,i)
         zz(9,i)=zln(9,i)-zln(3,i)*diag(4,i)-zz(6,i)*diag(5,i)
         zz(3,i)=zln(3,i)*diag(1,i)
         zz(6,i)=zz(6,i)*diag(3,i)
         zz(9,i)=zz(9,i)*diag(6,i)
         zz(6,i)=zz(6,i)-zz(9,i)*diag(5,i)
         zz(3,i)=zz(3,i)-zz(6,i)*diag(2,i)-zz(9,i)*diag(4,i)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine v6prod(zln,diag,zz,n)
      implicit double precision(a-h,o-z)
      dimension zln(9,n),diag(6,n),zz(9,n)
      do 100 i=1,n
         zz(4,i)=zln(4,i)-zln(1,i)*diag(2,i)
         zz(7,i)=zln(7,i)-zln(1,i)*diag(4,i)-zz(4,i)*diag(5,i)
         zz(1,i)=zln(1,i)*diag(1,i)
         zz(4,i)=zz(4,i)*diag(3,i)
         zz(7,i)=zz(7,i)*diag(6,i)
         zz(4,i)=zz(4,i)-zz(7,i)*diag(5,i)
         zz(1,i)=zz(1,i)-zz(4,i)*diag(2,i)-zz(7,i)*diag(4,i)
!c
         zz(5,i)=zln(5,i)-zln(2,i)*diag(2,i)
         zz(8,i)=zln(8,i)-zln(2,i)*diag(4,i)-zz(5,i)*diag(5,i)
         zz(2,i)=zln(2,i)*diag(1,i)
         zz(5,i)=zz(5,i)*diag(3,i)
         zz(8,i)=zz(8,i)*diag(6,i)
         zz(5,i)=zz(5,i)-zz(8,i)*diag(5,i)
         zz(2,i)=zz(2,i)-zz(5,i)*diag(2,i)-zz(8,i)*diag(4,i)
!c
         zz(6,i)=zln(6,i)-zln(3,i)*diag(2,i)
         zz(9,i)=zln(9,i)-zln(3,i)*diag(4,i)-zz(6,i)*diag(5,i)
         zz(3,i)=zln(3,i)*diag(1,i)
         zz(6,i)=zz(6,i)*diag(3,i)
         zz(9,i)=zz(9,i)*diag(6,i)
         zz(6,i)=zz(6,i)-zz(9,i)*diag(5,i)
         zz(3,i)=zz(3,i)-zz(6,i)*diag(2,i)-zz(9,i)*diag(4,i)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine vcopy(a,c,n)
      implicit double precision(a-h,o-z)
      dimension a(n),c(n)
!c     do 100 i=1,n
!c        c(i)=a(i)
!c 100 continue
      c=a
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine verif0(neqns,ndeg,nttbr,irow,jcol,val,rhs,x)
      implicit double precision(a-h,o-z)
      dimension irow(*),jcol(*),val(ndeg,ndeg,*),rhs(ndeg,*),x(ndeg,*)
!c
!c----------------------------------------------------------------------
!c
!c     verify the solution(symmetric matrix)
!c
!c----------------------------------------------------------------------
!c
      rel=0.0d0
      do 10 i=1,neqns
         do 11 l=1,ndeg
         rel=rel+dabs(rhs(l,i))
   11 continue
   10 continue
      do 100 k=1,nttbr
         i=irow(k)
         j=jcol(k)
         do 101 l=1,ndeg
         do 102 m=1,ndeg
         rhs(l,i)=rhs(l,i)-val(l,m,k)*x(m,j)
         if(i.ne.j) rhs(l,j)=rhs(l,j)-val(m,l,k)*x(m,i)
  102    continue
  101    continue
  100 continue
      err=0.0d0
      do 200 i=1,neqns
         do 210 l=1,ndeg
         err=err+dabs(rhs(l,i))
  210 continue
  200 continue
      write(6,6000) err,rel,err/rel
!WINDEBUG
!      write(16,6000) err,rel,err/rel
 6000 format(' ***verification***(symmetric)'/ &
             'norm(Ax-b)            =  ',1pd20.10/ &
             'norm(b)               =  ',1pd20.10/ &
             'norm(Ax-b)/norm(b)    =  ',1pd20.10)
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine vprod(a,b,c,n)
      implicit double precision(a-h,o-z)
      dimension a(n),b(n),c(n)
      do 100 i=1,n
         c(i)=a(i)*b(i)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine vxprod(ndeg,ndegl,zln,diag,zz,n)
      implicit double precision(a-h,o-z)
      dimension zln(ndeg*ndeg,n),diag(ndegl,n),zz(ndeg*ndeg,n)
      do 100 i=1,n
         call invxx(zz(1,i),zln(1,i),diag(1,i),ndeg)
  100 continue
      return
      end subroutine
!======================================================================!
!                                                                      !
!======================================================================!
      subroutine zpivot(neqns,neqnsz,nttbr,jcol,irow,zpiv,ir)
      integer jcol(*),irow(*),zpiv(*)
      common /debug/idbg
!c
      ir=0
      do 100 l=1,neqns
         zpiv(l)=1
  100 continue
!c
!c     do 200 l=neqns+1,neqns+nttbr
      do 200 l=1,nttbr
         i=irow(l)
         j=jcol(l)
         if(i.le.0.or.j.le.0) then
            ir=-1
            goto 1000
         elseif(i.gt.neqns.or.j.gt.neqns) then
            ir=1
            goto 1000
         endif
         if(i.eq.j) zpiv(i)=0
  200 continue
      do 310 i=neqns,1,-1
         if(zpiv(i).eq.0) then
         neqnsz=i
         goto 320
         endif
  310 continue
  320 continue
 1000 continue
      if(idbg.ne.0) write(6,60) (zpiv(i),i=1,neqns)
   60 format(20i3)
      return
      end subroutine
      end module hecmw_solver_direct
