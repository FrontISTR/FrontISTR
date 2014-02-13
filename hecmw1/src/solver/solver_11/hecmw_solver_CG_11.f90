!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!
!C***
!C*** module hecmw_solver_CG_11
!C***
!
      module hecmw_solver_CG_11
      contains

!C
!C*** CG
!C
      subroutine hecmw_solve_CG_11                                      &
     &    (N, NP,  NPL, NPU, D,  AL, INL, IAL, AU, INU, IAU,            &
     &     B, X, ALU, RESID, ITER, ERROR, my_rank,                      &
     &     NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                 &
     &                        STACK_EXPORT, NOD_EXPORT,                 &
     &                        SOLVER_COMM , PRECOND, iterPREmax,        &
     &     Tset, Tsol, Tcomm, ITERlog)

      use  hecmw_util
      use  hecmw_solver_SR_11
      use m_hecmw_solve_error
      use m_hecmw_comm_f

      implicit none

      integer(kind=kint )                                 :: N,NP,NPL,NPU
      integer(kind=kint )                                 :: NEIBPETOT
      real   (kind=kreal),                   intent(inout):: RESID
      integer(kind=kint ),                   intent(inout):: ITER
      integer(kind=kint ),                   intent(inout):: ERROR
      integer(kind=kint ),                   intent(in   ):: my_rank

      integer(kind=kint )                  , intent(in)   :: SOLVER_COMM
      integer(kind=kint ), intent(in):: PRECOND, ITERlog

      real   (kind=kreal), intent(inout):: Tset, Tsol, Tcomm

      real   (kind=kreal), dimension(NP )  , intent(inout)::  D
      real   (kind=kreal), dimension(NP )  , intent(inout)::  B
      real   (kind=kreal), dimension(NP )  , intent(inout)::  X
      real   (kind=kreal), dimension(N  )  , intent(inout)::  ALU

      real   (kind=kreal), dimension(NPU)  , intent(inout)::  AU
      real   (kind=kreal), dimension(NPU)  , intent(inout)::  AL

      integer(kind=kint ), dimension(0:NPU)  , intent(in) ::  INU
      integer(kind=kint ), dimension(  NPU)  , intent(in) ::  IAU
      integer(kind=kint ), dimension(0:NPL)  , intent(in) ::  INL
      integer(kind=kint ), dimension(  NPL)  , intent(in) ::  IAL

      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:)
      integer(kind=kint ), pointer :: NOD_IMPORT  (:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:)
      integer(kind=kint ), pointer :: NOD_EXPORT  (:)

      real   (kind=kreal), dimension(:),   allocatable :: WS, WR
      real   (kind=kreal), dimension(:,:), allocatable :: WW

      integer(kind=kint ) :: P, Q, R, Z, ZP, MAXIT, iterPREmax
      real   (kind=kreal) :: TOL

      ! local variables
      integer(kind=kint )::i,j,k,isU,ieU,isL,ieL,iterPRE,inod,indexA,indexB
      integer(kind=kint )::ns,nr
      real   (kind=kreal)::S_TIME,S1_TIME,E_TIME,E1_TIME
      real   (kind=kreal)::WVAL,BNRM20,BNRM2,SW,WV,X1,X2
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,C10,C1,ALPHA,DNRM20,DNRM2

      S_TIME= HECMW_WTIME()
!C
!C-- INIT.
      ERROR= 0

      if( NEIBPETOT > 0 ) then
        ns = STACK_EXPORT(NEIBPETOT)
        nr = STACK_IMPORT(NEIBPETOT)
      else
        ns = 0
        nr = 0
      end if

      allocate (WW(NP,4))
      allocate (WS(ns))
      allocate (WR(nr))

      R = 1
      Z = 2
      Q = 2
      P = 3
      ZP= 4

      MAXIT  = ITER
       TOL   = RESID

      WW= 0.d0
       X= 0.d0
      WR= 0.d0
      WS= 0.d0

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      call hecmw_solve_SEND_RECV_11                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM, my_rank)

!C
!C-- BEGIN calculation

      do j= 1, N
        WVAL= B(j) - D(j) * X(j)
        isU= INU(j-1) + 1
        ieU= INU(j  )
        do i= isU, ieU
          inod= IAU(i)
          WVAL= WVAL - AU(i) * X(inod)
        enddo

        isL= INL(j-1) + 1
        ieL= INL(j  )
        do i= isL, ieL
          inod= IAL(i)
          WVAL= WVAL - AL(i) * X(inod)
        enddo
        WW(j,R)= WVAL
      enddo

      BNRM20 = 0.d0
      do i= 1, N
        BNRM20 = BNRM20 + B(i)**2
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM,SOLVER_COMM)
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      E_TIME= HECMW_WTIME()
      Tset= Tset + E_TIME - S_TIME

      S1_TIME= HECMW_WTIME()
!C===
      do iter= 1, MAXIT
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

!C
!C-- incomplete CHOLESKY

      if (PRECOND.eq.1 .or. PRECOND.eq.2) then

        do i= 1, N
          WW(i,ZP)= WW(i,R)
        enddo

        do i= 1, NP
          WW(i,Z )= 0.d0
        enddo

      do iterPRE= 1, iterPREmax

        do i= 1, N
          isL= INL(i-1) + 1
          ieL= INL(i  )
          WVAL= WW(i,ZP)
          do j= isL, ieL
            inod = IAL(j)
            WVAL=  WVAL -  AL(j) * WW(inod,ZP)
          enddo
          WW(i,ZP)= WVAL * ALU(i)
        enddo

        do i= N, 1, -1
          SW  = 0.0d0
          isU= INU(i-1) + 1
          ieU= INU(i  )
          do j= isU, ieU
            inod = IAU(j)
            SW= SW + AU(j) * WW(inod,ZP)
          enddo
          WW(i,ZP)= WW(i,ZP) - ALU(i) * SW
        enddo

        indexA= ZP
        indexB= Z

        do i= 1, N
          WW(i,indexB)= WW(i,indexB) + WW(i,indexA)
        enddo

        if (iterPRE.eq.iterPREmax) goto 750

      S_TIME= HECMW_WTIME()
      call hecmw_solve_SEND_RECV_11                                     &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,indexB) , SOLVER_COMM,  &
     &       my_rank)
      E_TIME= HECMW_WTIME()
      Tcomm= Tcomm + E_TIME - S_TIME

        do j= 1, N
          WV= WW(j,R) - D(j)*WW(j,indexB)
          do k= INL(j-1)+1, INL(j)
             i= IAL(k)
            WV= WV - AL(k) * WW(i,indexB)
          enddo
          do k= INU(j-1)+1, INU(j)
            i= IAU(k)
            WV= WV - AU(k) * WW(i,indexB)
          enddo
          WW(j,indexA)= WV
        enddo

 750    continue

      enddo

      endif

!C
!C-- diagonal
      if (PRECOND.eq.3) then
      do i= 1, N
        WW(i,Z)=  WW(i,R) * ALU(i)
      enddo
      endif
!C===

!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.0
      do i= 1, N
        RHO0= RHO0 + WW(i,R)*WW(i,Z)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM,SOLVER_COMM)
!C      call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      Tcomm= Tcomm + E_TIME - S_TIME
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then
        do i= 1, N
          WW(i,P)= WW(i,Z)
        enddo

       else

         BETA= RHO / RHO1

         do i= 1, N
           WW(i,P)= WW(i,Z) + BETA*WW(i,P)
         enddo

      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

      S_TIME= HECMW_WTIME()
      call hecmw_solve_SEND_RECV_11                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,P) , SOLVER_COMM,     &
     &     my_rank)
      E_TIME= HECMW_WTIME()
      Tcomm= Tcomm + E_TIME - S_TIME

!C
!C-- BEGIN calculation

      do j= 1, N
        WVAL= D(j) * WW(j,P)

        isU= INU(j-1) + 1
        ieU= INU(j  )

        do i= isU, ieU
          inod= IAU(i)
          WVAL= WVAL + AU(i) * WW(inod,P)
        enddo

        isL= INL(j-1) + 1
        ieL= INL(j  )

        do i= isL, ieL
          inod= IAL(i)
          WVAL= WVAL + AL(i) * WW(inod,P)
        enddo
        WW(j,Q)= WVAL
      enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0

      do i= 1, N
        C10= C10 + WW(i,P)*WW(i,Q)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (C10, C1, HECMW_SUM,SOLVER_COMM)
!C      call MPI_allREDUCE (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      Tcomm= Tcomm + E_TIME - S_TIME

      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===
      DNRM20= 0.d0

      X1= 0.0d0
      X2= 0.0d0

      do i= 1, N
         X(i)  = X (i)   + ALPHA * WW(i,P)
        WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo

      DNRM20 = 0.0
      do i= 1, N
        DNRM20= DNRM20 + WW(i,R)**2
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      Tcomm= Tcomm + E_TIME - S_TIME

        RESID= dsqrt(DNRM2/BNRM2)

        if (my_rank.eq.0.and.ITERlog.eq.1) write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -100

        RHO1 = RHO

      enddo
!C===
      E1_TIME= HECMW_WTIME()
      Tsol= E1_TIME - S1_TIME

!C
!C-- INTERFACE data EXCHANGE

      call hecmw_solve_SEND_RECV_11                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      deallocate (WW, WR, WS)

      close (12)

      end subroutine hecmw_solve_CG_11
      end module     hecmw_solver_CG_11





