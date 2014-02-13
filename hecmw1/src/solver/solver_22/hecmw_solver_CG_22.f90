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

!C
!C***
!C*** module hecmw_solver_CG_22
!C***
!C
      module hecmw_solver_CG_22
      contains
!C
!C*** CG_3
!C
      subroutine hecmw_solve_CG_22                                      &
     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,               &
     &    B, X, ALU, RESID, ITER, ERROR, my_rank,                       &
     &    NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                  &
     &                       STACK_EXPORT, NOD_EXPORT,                  &
     &                       SOLVER_COMM , PRECOND, iterPREmax,         &
     &    Tset, Tsol, Tcomm, ITERlog)

      use hecmw_util
      use  hecmw_solver_SR_22
      use m_hecmw_solve_error
      use m_hecmw_comm_f

      implicit none
!      include  'hecmw_config_f.h'
!      include  'mpif.h'

      integer(kind=kint ), intent(in):: N, NP, NPU, NPL, my_rank
      integer(kind=kint ), intent(in):: NEIBPETOT
      integer(kind=kint ), intent(in):: SOLVER_COMM
      integer(kind=kint ), intent(in):: PRECOND, ITERlog

      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      real(kind=kreal), dimension(2*NP)   , intent(inout):: B, X
      real(kind=kreal), dimension(4*NPL), intent(inout):: AL
      real(kind=kreal), dimension(4*NPU), intent(inout):: AU
      real(kind=kreal), dimension(4*N  ), intent(inout):: ALU
      real(kind=kreal), dimension(4*NP ), intent(inout):: D

      integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
      integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:), NOD_IMPORT(:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:), NOD_EXPORT(:)

      real(kind=kreal), dimension(:),    allocatable       :: WS, WR
      real(kind=kreal), dimension(:,:),  allocatable       :: WW

      integer(kind=kint), parameter ::  R= 1
      integer(kind=kint), parameter ::  Z= 2
      integer(kind=kint), parameter ::  Q= 2
      integer(kind=kint), parameter ::  P= 3
      integer(kind=kint), parameter :: ZP= 4

      integer(kind=kint ) :: MAXIT, iterPREmax

      real   (kind=kreal) :: TOL

      ! local variables
      integer(kind=kint )::i,j,k,isU,ieU,isL,ieL,iterPRE,indexA,indexB,indexC
      integer(kind=kint )::ns,nr
      real   (kind=kreal)::S_TIME,S1_TIME,E_TIME,E1_TIME,START_TIME,END_TIME
      real   (kind=kreal)::BNRM20,BNRM2,X1,X2
      real   (kind=kreal)::RHO,RHO0,RHO1,BETA,C10,C1,ALPHA,DNRM20,DNRM2
      real   (kind=kreal)::COMMtime,COMPtime,WVAL1,WVAL2,SW1,SW2,WV1,WV2

      S_TIME= HECMW_WTIME()
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      ERROR= 0

      COMMtime= 0.d0
      COMPtime= 0.d0

      if( NEIBPETOT > 0 ) then
        ns = STACK_EXPORT(NEIBPETOT)
        nr = STACK_IMPORT(NEIBPETOT)
      else
        ns = 0
        nr = 0
      end if

      allocate (WW(2*NP,4))
      allocate (WS(2*ns), WR(2*nr))

      MAXIT  = ITER
       TOL   = RESID

      X = 0.d0
      WW= 0.d0
      WS= 0.d0
      WR= 0.d0

!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

!C
!C-- BEGIN calculation

      do j= 1, N
        X1= X(2*j-1)
	X2= X(2*j  )
	WVAL1= B(2*j-1) - D(4*j-3)*X1 - D(4*j-2)*X2
	WVAL2= B(2*j  ) - D(4*j-1)*X1 - D(4*j  )*X2
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= X(2*i-1)
          X2= X(2*i  )
          WVAL1= WVAL1 - AL(4*k-3)*X1 - AL(4*k-2)*X2
          WVAL2= WVAL2 - AL(4*k-1)*X1 - AL(4*k  )*X2
        enddo

        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= X(2*i-1)
          X2= X(2*i  )
          WVAL1= WVAL1 - AU(4*k-3)*X1 - AU(4*k-2)*X2
          WVAL2= WVAL2 - AU(4*k-1)*X1 - AU(4*k  )*X2
        enddo
        WW(2*j-1,R)= WVAL1
        WW(2*j  ,R)= WVAL2
      enddo

      BNRM20= 0.d0
      do i= 1, N
        BNRM20= BNRM20 + B(2*i-1)**2 + B(2*i)**2
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      E_TIME= HECMW_WTIME()
      Tset= Tset + E_TIME - S_TIME
!C===

      S1_TIME= HECMW_WTIME()
      do iter= 1, MAXIT
!      call MPI_BARRIER(SOLVER_COMM,ierr)
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

      indexA= ZP
      indexB= Z
      indexC= R

      if (PRECOND.le.2) then
!C
!C== Block SSOR
      do i= 1, N
        WW(2*i-1,indexA)= WW(2*i-1,indexC)
        WW(2*i  ,indexA)= WW(2*i  ,indexC)
      enddo

      do i= 1, NP
        WW(2*i-1,indexB )= 0.d0
        WW(2*i  ,indexB )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax

!C
!C-- FORWARD

        do i= 1, N
          SW1= WW(2*i-1,indexA)
          SW2= WW(2*i  ,indexA)
          isL= INL(i-1)+1
          ieL= INL(i)
         do j= isL, ieL
              k= IAL(j)
             X1= WW(2*k-1,indexA)
             X2= WW(2*k  ,indexA)
            SW1= SW1 - AL(4*j-3)*X1 - AL(4*j-2)*X2
            SW2= SW2 - AL(4*j-1)*X1 - AL(4*j  )*X2
          enddo

          X1= SW1
          X2= SW2
          X2= SW2 - ALU(4*i-1)*X1
          X2= ALU(4*i  )*  X2
          X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
          WW(2*i-1,indexA)= X1
          WW(2*i  ,indexA)= X2
        enddo

!C
!C-- BACKWARD

        do i= N, 1, -1
          isU= INU(i-1) + 1
          ieU= INU(i)
          SW1= 0.d0
          SW2= 0.d0
          do j= isU, ieU
              k= IAU(j)
             X1= WW(2*k-1,indexA)
             X2= WW(2*k  ,indexA)
            SW1= SW1 + AU(4*j-3)*X1 + AU(4*j-2)*X2
            SW2= SW2 + AU(4*j-1)*X1 + AU(4*j  )*X2
          enddo

          X1= SW1
          X2= SW2
          X2= SW2 - ALU(4*i-1)*X1
          X2= ALU(4*i  )*  X2
          X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
          WW(2*i-1,indexA)= WW(2*i-1,indexA) - X1
          WW(2*i  ,indexA)= WW(2*i  ,indexA) - X2
        enddo

!C
!C-- additive Schwartz

        do i= 1, N
          WW(2*i-1,indexB)= WW(2*i-1,indexB) + WW(2*i-1,indexA)
          WW(2*i  ,indexB)= WW(2*i  ,indexB) + WW(2*i  ,indexA)
        enddo


        if (iterPRE.eq.iterPREmax) goto 750

        call hecmw_solve_SEND_RECV_22                                   &
     &     ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,           &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,indexB) ,           &
     &       SOLVER_COMM, my_rank)

        do j= 1, N
           X1= WW(2*j-1,indexB)
           X2= WW(2*j  ,indexB)
          WV1= WW(2*j-1,indexC) - D(4*j-3)*X1-D(4*j-2)*X2
          WV2= WW(2*j  ,indexC) - D(4*j-1)*X1-D(4*j  )*X2
          do k= INL(j-1)+1, INL(j)
              i= IAL(k)
             X1= WW(2*i-1,indexB)
             X2= WW(2*i  ,indexB)
            WV1= WV1 - AL(4*k-3)*X1 - AL(4*k-2)*X2
            WV2= WV2 - AL(4*k-1)*X1 - AL(4*k  )*X2
          enddo
          do k= INU(j-1)+1, INU(j)
              i= IAU(k)
             X1= WW(2*i-1,indexB)
             X2= WW(2*i  ,indexB)
            WV1= WV1 - AU(4*k-3)*X1 - AU(4*k-2)*X2
            WV2= WV2 - AU(4*k-1)*X1 - AU(4*k  )*X2
          enddo
          WW(2*j-1,indexA)= WV1
          WW(2*j  ,indexA)= WV2
        enddo

 750    continue

      enddo
      endif

      if (PRECOND.eq.3) then
!C
!C== Block SCALING

      do i= 1, N
        WW(2*i-1,indexB)= WW(2*i-1,indexC)
        WW(2*i  ,indexB)= WW(2*i  ,indexC)
      enddo

      do i= 1, N
        X1= WW(2*i-1,indexB)
        X2= WW(2*i  ,indexB)
        X2= X2 - ALU(4*i-1)*X1
        X2= ALU(4*i  )*  X2
        X1= ALU(4*i-3)*( X1 - ALU(4*i-2)*X2)
        WW(2*i-1,indexB)= X1
        WW(2*i  ,indexB)= X2
      enddo
      endif
!C===

!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.d0

      do i= 1, N
        RHO0= RHO0 + WW(2*i-1,R)*WW(2*i-1,Z) + WW(2*i  ,R)*WW(2*i  ,Z)
      enddo

      START_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then

        do i= 1, N
          WW(2*i-1,P)= WW(2*i-1,Z)
          WW(2*i  ,P)= WW(2*i  ,Z)
        enddo
       else
         BETA= RHO / RHO1
         do i= 1, N
           WW(2*i-1,P)= WW(2*i-1,Z) + BETA*WW(2*i-1,P)
           WW(2*i  ,P)= WW(2*i  ,Z) + BETA*WW(2*i  ,P)
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
      START_TIME= HECMW_WTIME()
      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(:,P) , SOLVER_COMM,     &
     &     my_rank)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME

      indexA= P
      indexB= Q
!C
!C-- BEGIN calculation

      do j= 1, N
           X1= WW(2*j-1,indexA)
           X2= WW(2*j  ,indexA)
        WVAL1= D(4*j-3)*X1 + D(4*j-2)*X2
        WVAL2= D(4*j-1)*X1 + D(4*j  )*X2
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(2*i-1,indexA)
          X2= WW(2*i  ,indexA)
          WVAL1= WVAL1 + AL(4*k-3)*X1 + AL(4*k-2)*X2
          WVAL2= WVAL2 + AL(4*k-1)*X1 + AL(4*k  )*X2
        enddo
        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(2*i-1,indexA)
          X2= WW(2*i  ,indexA)
          WVAL1= WVAL1 + AU(4*k-3)*X1 + AU(4*k-2)*X2
          WVAL2= WVAL2 + AU(4*k-1)*X1 + AU(4*k  )*X2
        enddo

        WW(2*j-1,indexB)= WVAL1
        WW(2*j  ,indexB)= WVAL2
      enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0
      do i= 1, N
        C10= C10 + WW(2*i-1,P)*WW(2*i-1,Q) + WW(2*i  ,P)*WW(2*i  ,Q)
      enddo

      START_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (C10, C1,  HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME

      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===

      do i= 1, N
         X(2*i-1)  = X (2*i-1)   + ALPHA * WW(2*i-1,P)
         X(2*i  )  = X (2*i  )   + ALPHA * WW(2*i  ,P)
        WW(2*i-1,R)= WW(2*i-1,R) - ALPHA * WW(2*i-1,Q)
        WW(2*i  ,R)= WW(2*i  ,R) - ALPHA * WW(2*i  ,Q)
      enddo

      DNRM20= 0.d0
      do i= 1, N
        DNRM20= DNRM20 + WW(2*i-1,R)**2 + WW(2*i  ,R)**2
      enddo

      START_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      END_TIME= HECMW_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
        if (my_rank.eq.0.and.ITERlog.eq.1) write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
! 1010   format (1pe16.6)
!C#####

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -300

        RHO1 = RHO
      enddo
!C===

!C
!C-- INTERFACE data EXCHANGE
   30 continue
      E1_TIME= HECMW_WTIME()
      COMPtime= E1_TIME - S1_TIME

      Tsol = COMPtime
      Tcomm= COMMtime

      call hecmw_solve_SEND_RECV_22                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      deallocate (WW, WR, WS)

      end subroutine hecmw_solve_CG_22
      end module     hecmw_solver_CG_22

