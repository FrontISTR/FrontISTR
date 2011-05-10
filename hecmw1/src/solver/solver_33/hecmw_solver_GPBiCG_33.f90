!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.1                                                !
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

!C*** 
!C*** module hecmw_solver_GPBiCG_33
!C***
!
      module hecmw_solver_GPBiCG_33
      contains
!C
!C*** hecmw_solve_GPBiCG_33
!C
      subroutine hecmw_solve_GPBiCG_33                                  &
     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,               &
     &    B, X, ALU, RESID, ITER, ERROR, my_rank,                       &
     &    NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,                  &
     &                       STACK_EXPORT, NOD_EXPORT,                  &
     &                       SOLVER_COMM , PRECOND, iterPREmax,         &
     &    Tset, Tsol, Tcomm, ITERlog)     

      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_solver_SR_33

      implicit none

      integer(kind=kint ),                   intent(in   )::  N
      integer(kind=kint ),                   intent(in   )::  NP
      integer(kind=kint ),                   intent(in   )::  NPU, NPL

      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      integer(kind=kint ),                   intent(inout)::  ITER
      integer(kind=kint ),                   intent(inout)::  ERROR
      integer(kind=kint ),                   intent(in   )::  my_rank

      integer(kind=kint )                  , intent(in)   :: SOLVER_COMM

      real   (kind=kreal), dimension(3*NP) , intent(inout)::  B
      real   (kind=kreal), dimension(3*NP) , intent(inout)::  X

      real   (kind=kreal), dimension(9*NPL), intent(inout)::  AL
      real   (kind=kreal), dimension(9*NPU), intent(inout)::  AU
      real   (kind=kreal), dimension(9*N  ), intent(inout):: ALU
      real   (kind=kreal), dimension(9*NP ), intent(inout)::  D

      integer(kind=kint ), dimension(0:NP), intent(in) :: INU, INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
      integer(kind=kint ), dimension(  NPU),intent(in) :: IAU
      integer(kind=kint ), intent(in):: PRECOND, ITERlog

      integer(kind=kint ):: NEIBPETOT, iterPREmax
      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: STACK_IMPORT(:)
      integer(kind=kint ), pointer :: NOD_IMPORT  (:)
      integer(kind=kint ), pointer :: STACK_EXPORT(:)
      integer(kind=kint ), pointer :: NOD_EXPORT  (:)

      real(kind=kreal), dimension(:),    allocatable       ::  WS, WR
      real(kind=kreal), dimension(:,:),  allocatable       ::  WW

      real(kind=kreal), dimension(5) :: C0, CG
      real(kind=kreal), dimension(2) :: EQ

      integer(kind=kint ) :: R, RT, T, T0, TT, P, PT, ZQ
      integer(kind=kint ) :: U, W1, Y, Z, WK, W2, MAXIT
      real   (kind=kreal) :: TOL
      integer(kind=kint ) :: i,j,k,isU,ieU,isL,ieL,iterPRE,indexA,indexB,indexC
      integer(kind=kint ) :: jSR,jSB,ik,kSR,kSB
      real   (kind=kreal) :: S_TIME,S1_TIME,E_TIME,E1_TIME
      real   (kind=kreal) :: BNRM20,BNRM2,X1,X2,X3
      real   (kind=kreal) :: RHO,RHO0,RHO1,BETA,ALPHA,DNRM20,DNRM2,RHO10
      real   (kind=kreal) :: WVAL1,WVAL2,WVAL3,SW1,SW2,SW3
      real   (kind=kreal) :: COMMtime,COMPtime,QSI,ETA,COEF10,COEF1

      S_TIME= HECMW_WTIME()
!C
!C-- INIT.
      ERROR= 0

      allocate (WW(3*NP,14))
      allocate (WS(3*NP))
      allocate (WR(3*NP))

      WW= 0.d0

      MAXIT = ITER
      TOL   = RESID

      COMMtime= 0.d0
      COMPtime= 0.d0

      ERROR= 0

       R= 1
      RT= 2
       T= 3
      TT= 4
      T0= 5
       P= 6
      PT= 7
       U= 8
      W1= 9
       Y=10
       Z=11
      WK=12
      W2=13
      ZQ=14
!C
!C +----------------------+
!C | {r}= {b} - [A]{xini} |
!C +----------------------+
!C===
      indexB= R
      call hecmw_solver_init_resid_33 (indexB)

      do i= 1, N
        WW(3*i-2,RT)= WW(3*i-2,R)
        WW(3*i-1,RT)= WW(3*i-1,R)
        WW(3*i  ,RT)= WW(3*i  ,R)
      enddo
!C==
      BNRM20= 0.d0
      RHO0  = 0.0d0
      do i= 1, N
        BNRM20= BNRM20+B(3*i-2)**2+B(3*i-1)**2+B(3*i)**2
        RHO0  = RHO0 + WW(3*i-2,RT)*WW(3*i-2,R)+WW(3*i-1,RT)*WW(3*i-1,R)&
     &                                         +WW(3*i  ,RT)*WW(3*i  ,R)
      enddo

      call HECMW_allREDUCE_DP1 (BNRM20, BNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      call HECMW_allREDUCE_DP1 (RHO0, RHO, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO0  , RHO,   1, MPI_DOUBLE_PRECISION,       &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)

      if (BNRM2.eq.0.d0) BNRM2= 1.d0      
      E_TIME= HECMW_WTIME()
      Tset= Tset + E_TIME - S_TIME
!C===

!C
!C*************************************************************** ITERATIVE PROC.
!C
      S1_TIME= HECMW_WTIME()
      do iter= 1, MAXIT
!C
!C +----------------+
!C | {r}= [Minv]{r} |
!C +----------------+
!C===
      do j= 1, N
        WW(3*j-2,WK)=  WW(3*j-2,R)
        WW(3*j-1,WK)=  WW(3*j-1,R)
        WW(3*j  ,WK)=  WW(3*j  ,R)
      enddo

      indexA= ZQ
      indexB= R
      indexC= WK
      call hecmw_solver_precond_33 (indexA, indexB, indexC)
!C===

!C
!C +----------------------------------+
!C | {p} = {r} + BETA * ( {p} - {u} ) |
!C +----------------------------------+
!C===
      if (iter.gt.1) then
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R) + BETA*( WW(3*j-2,P)-WW(3*j-2,U))
          WW(3*j-1,P)= WW(3*j-1,R) + BETA*( WW(3*j-1,P)-WW(3*j-1,U))
          WW(3*j  ,P)= WW(3*j  ,R) + BETA*( WW(3*j  ,P)-WW(3*j  ,U))
        enddo
       else
        do j= 1, N
          WW(3*j-2,P)= WW(3*j-2,R)
          WW(3*j-1,P)= WW(3*j-1,R)
          WW(3*j  ,P)= WW(3*j  ,R)
        enddo
      endif
!C===

!C
!C +--------------------------------+
!C | ALPHA= {r_tld}{r}/{r_tld} A{p} |
!C +--------------------------------+
!C===

!C
!C-- calc. {p_tld}= A{p}

      indexA= P
      indexB= PT
      call hecmw_solver_matvec_33  (indexA, indexB)

!C
!C-- calc. ALPHA

      RHO10= 0.d0
      do j= 1, N
        RHO10= RHO10+WW(3*j-2,RT)*WW(3*j-2,PT)+WW(3*j-1,RT)*WW(3*j-1,PT)&
     &                                        +WW(3*j  ,RT)*WW(3*j  ,PT)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (RHO10, RHO1, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (RHO10, RHO1, 1, MPI_DOUBLE_PRECISION,         &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      ALPHA= RHO / RHO1
!C===

!C
!C +------------------------------------------+
!C | {y}= {t} - {r} - ALPHA{w} + ALPHA{p_tld} |
!C | {t}= {r}                  - ALPHA{p_tld} |
!C +------------------------------------------+
!C===
      do j= 1, N
        WW(3*j-2,Y)= WW(3*j-2,T) - WW(3*j-2,WK)                         &
     &                   + ALPHA*(-WW(3*j-2,W1) + WW(3*j-2,PT))
        WW(3*j-1,Y)= WW(3*j-1,T) - WW(3*j-1,WK)                         &
     &                   + ALPHA*(-WW(3*j-1,W1) + WW(3*j-1,PT))
        WW(3*j  ,Y)= WW(3*j  ,T) - WW(3*j  ,WK)                         &
     &                   + ALPHA*(-WW(3*j  ,W1) + WW(3*j  ,PT))

        WW(3*j-2,T)= WW(3*j-2,WK) - ALPHA*WW(3*j-2,PT)
        WW(3*j-1,T)= WW(3*j-1,WK) - ALPHA*WW(3*j-1,PT)
        WW(3*j  ,T)= WW(3*j  ,WK) - ALPHA*WW(3*j  ,PT)
      enddo
!C===

!C
!C +-----------------------+
!C | {t_tld}= [A][Minv]{t} |
!C +-----------------------+
!C===

!C
!C-- calc. {t_tld} and {t0} by [M] inversion
!C         {W2}   = [Minv]{p_tld} 
!C
      indexA= ZQ
      indexB= TT
      indexC= T
      call hecmw_solver_precond_33 (indexA, indexB, indexC)

      indexA= ZQ
      indexB= W2
      indexC= T0
      call hecmw_solver_precond_33 (indexA, indexB, indexC)
      do i= 1, N
        WW(3*i-2,T0)= WW(3*i-2,W2)
        WW(3*i-1,T0)= WW(3*i-1,W2)
        WW(3*i  ,T0)= WW(3*i  ,W2)
      enddo

      indexA= ZQ
      indexB= W2
      indexC= PT
      call hecmw_solver_precond_33 (indexA, indexB, indexC)
!C===

!C
!C-- calc. [A]{t_tld}
      indexA= TT
      indexB= WK
      call hecmw_solver_matvec_33  (indexA, indexB)

      do i= 1, N
        WW(3*i-2,TT)= WW(3*i-2,WK)
        WW(3*i-1,TT)= WW(3*i-1,WK)
        WW(3*i  ,TT)= WW(3*i  ,WK)
      enddo
!C===

!C
!C +-------------------+
!C | calc. QSI and ETA |
!C +-------------------+
!C===
      C0(1)= 0.0d0
      C0(2)= 0.0d0
      C0(3)= 0.0d0
      C0(4)= 0.0d0
      C0(5)= 0.0d0
      CG(1)= 0.0d0
      CG(2)= 0.0d0
      CG(3)= 0.0d0
      CG(4)= 0.0d0
      CG(5)= 0.0d0

      do j= 1, N
        C0(1)= C0(1)+WW(3*j-2, Y)*WW(3*j-2, Y)+WW(3*j-1, Y)*WW(3*j-1, Y)&
     &                                        +WW(3*j  , Y)*WW(3*j  , Y)
        C0(2)= C0(2)+WW(3*j-2,TT)*WW(3*j-2, T)+WW(3*j-1,TT)*WW(3*j-1, T)&
     &                                        +WW(3*j  ,TT)*WW(3*j  , T)
        C0(3)= C0(3)+WW(3*j-2, Y)*WW(3*j-2, T)+WW(3*j-1, Y)*WW(3*j-1, T)&
     &                                        +WW(3*j  , Y)*WW(3*j  , T)
        C0(4)= C0(4)+WW(3*j-2,TT)*WW(3*j-2, Y)+WW(3*j-1,TT)*WW(3*j-1, Y)&
     &                                        +WW(3*j  ,TT)*WW(3*j  , Y)
        C0(5)= C0(5)+WW(3*j-2,TT)*WW(3*j-2,TT)+WW(3*j-1,TT)*WW(3*j-1,TT)&
     &                                        +WW(3*j  ,TT)*WW(3*j  ,TT)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP (C0, CG, 5, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE (C0, CG,  5, MPI_DOUBLE_PRECISION,             &
!C     &                    MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      if (iter.eq.1) then
        EQ(1)= CG(2)/CG(5)
        EQ(2)= 0.d0
       else
        EQ(1)= (CG(1)*CG(2)-CG(3)*CG(4)) / (CG(5)*CG(1)-CG(4)*CG(4))
        EQ(2)= (CG(5)*CG(3)-CG(4)*CG(2)) / (CG(5)*CG(1)-CG(4)*CG(4))
      endif

      QSI= EQ(1)
      ETA= EQ(2)
!C===

!C
!C +----------------------------------------------------------+
!C | {u} = QSI [Minv]{pt} + ETA([Minv]{t0}-[Minv]{r}+BETA*{u} |
!C | {z} = QSI [Minv]{r}  + ETA*{z} - ALPHA*{u}               |
!C +----------------------------------------------------------+
!C===

!C
!C-- compute. {u},{z}

      do j= 1, N
        WW(3*j-2,U)= QSI* WW(3*j-2,W2) +                                &
     &               ETA*(WW(3*j-2,T0) - WW(3*j-2,R) + BETA*WW(3*j-2,U))
        WW(3*j-2,Z)= QSI* WW(3*j-2, R) +                                &
     &               ETA* WW(3*j-2, Z) -              ALPHA*WW(3*j-2,U)
        WW(3*j-1,U)= QSI* WW(3*j-1,W2) +                                &
     &               ETA*(WW(3*j-1,T0) - WW(3*j-1,R) + BETA*WW(3*j-1,U))
        WW(3*j-1,Z)= QSI* WW(3*j-1, R) +                                &
     &               ETA* WW(3*j-1, Z) -              ALPHA*WW(3*j-1,U)
        WW(3*j  ,U)= QSI* WW(3*j  ,W2) +                                &
     &               ETA*(WW(3*j  ,T0) - WW(3*j  ,R) + BETA*WW(3*j  ,U))
        WW(3*j  ,Z)= QSI* WW(3*j  , R) +                                &
     &               ETA* WW(3*j  , Z) -              ALPHA*WW(3*j  ,U)
      enddo
!C===

!C
!C +--------------------+
!C | update {x},{r},{w} |
!C +--------------------+
!C===
      DNRM20= 0.d0
      COEF10= 0.d0

      do j= 1, N
        X (3*j-2)= X(3*j-2) + ALPHA*WW(3*j-2,P) + WW(3*j-2,Z)
        X (3*j-1)= X(3*j-1) + ALPHA*WW(3*j-1,P) + WW(3*j-1,Z)
        X (3*j  )= X(3*j  ) + ALPHA*WW(3*j  ,P) + WW(3*j  ,Z)

        WW(3*j-2,R)= WW(3*j-2,T) - ETA*WW(3*j-2,Y) - QSI*WW(3*j-2,TT)
        WW(3*j-1,R)= WW(3*j-1,T) - ETA*WW(3*j-1,Y) - QSI*WW(3*j-1,TT)
        WW(3*j  ,R)= WW(3*j  ,T) - ETA*WW(3*j  ,Y) - QSI*WW(3*j  ,TT)

        WW(3*j-2,T0)= WW(3*j-2,T)
        WW(3*j-1,T0)= WW(3*j-1,T)
        WW(3*j  ,T0)= WW(3*j  ,T)

        DNRM20= DNRM20 + WW(3*j-2,R)**2+ WW(3*j-1,R)**2+ WW(3*j,R)**2
        COEF10= COEF10 + WW(3*j-2,R)*WW(3*j-2,RT)                       &
     &                 + WW(3*j-1,R)*WW(3*j-1,RT) + WW(3*j,R)*WW(3*j,RT)
      enddo

      S_TIME= HECMW_WTIME()
      call HECMW_allREDUCE_DP1 (DNRM20, DNRM2, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE  (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,      &
!C     &                     MPI_SUM, SOLVER_COMM, ierr)
      call HECMW_allREDUCE_DP1 (COEF10, COEF1, HECMW_SUM, SOLVER_COMM )
!C      call MPI_allREDUCE  (COEF10, COEF1, 1, MPI_DOUBLE_PRECISION,      &
!C     &                     MPI_SUM, SOLVER_COMM, ierr)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      BETA = ALPHA*COEF1 / (QSI*RHO)
      do j= 1, N
        WW(3*j-2,W1)= WW(3*j-2,TT) + BETA*WW(3*j-2,PT)
        WW(3*j-1,W1)= WW(3*j-1,TT) + BETA*WW(3*j-1,PT)
        WW(3*j  ,W1)= WW(3*j  ,TT) + BETA*WW(3*j  ,PT)
      enddo

      RESID= dsqrt(DNRM2/BNRM2)
      RHO  = COEF1

!C##### ITERATION HISTORY
      if (my_rank.eq.0 .and. ITERlog.eq.1)                              &
     &    write (*, 1000) ITER, RESID
 1000   format (i5, 1pe16.6)
!C#####

      if (RESID.le.TOL   ) exit
      if ( ITER.eq.MAXIT ) ERROR= -300
!C===
      enddo

      E1_TIME= HECMW_WTIME()
      COMPtime= E1_TIME - S1_TIME

      Tsol = COMPtime
      Tcomm= COMMtime

!C
!C-- INTERFACE data EXCHANGE
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP, NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,             &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)

      deallocate (WW)
      deallocate (WR)
      deallocate (WS)

      contains

!C
!C***
!C*** hecmw_solver_init_resid_33
!C***
!C

      subroutine hecmw_solver_init_resid_33 (indexB)
      implicit REAL*8 (A-H,O-Z)      
      integer(kind=kint) :: indexB

      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP,  NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,            &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, X , SOLVER_COMM,my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      do j= 1, N
          jSR= 3*j-2
          jSB= 9*j-8
           X1= X(jSR  )
           X2= X(jSR+1)
           X3= X(jSR+2)
        WVAL1= B(jSR  ) - D(jSB  )*X1 - D(jSB+1)*X2 - D(jSB+2)*X3
        WVAL2= B(jSR+1) - D(jSB+3)*X1 - D(jSB+4)*X2 - D(jSB+5)*X3
        WVAL3= B(jSR+2) - D(jSB+6)*X1 - D(jSB+7)*X2 - D(jSB+8)*X3
        do k= INL(j-1)+1, INL(j)
             ik= IAL(k)
            kSR= 3*ik- 2
            kSB= 9* k- 8
             X1= X(kSR  )
             X2= X(kSR+1)
             X3= X(kSR+2)
          WVAL1= WVAL1 - AL(kSB  )*X1 - AL(kSB+1)*X2 - AL(kSB+2)*X3
          WVAL2= WVAL2 - AL(kSB+3)*X1 - AL(kSB+4)*X2 - AL(kSB+5)*X3
          WVAL3= WVAL3 - AL(kSB+6)*X1 - AL(kSB+7)*X2 - AL(kSB+8)*X3
        enddo

        do k= INU(j-1)+1, INU(j)
           ik= IAU(k)
            kSR= 3*ik-2
            kSB= 9* k-8
             X1= X(kSR  )
             X2= X(kSR+1)
             X3= X(kSR+2)
          WVAL1= WVAL1 - AU(kSB  )*X1 - AU(kSB+1)*X2 - AU(kSB+2)*X3
          WVAL2= WVAL2 - AU(kSB+3)*X1 - AU(kSB+4)*X2 - AU(kSB+5)*X3
          WVAL3= WVAL3 - AU(kSB+6)*X1 - AU(kSB+7)*X2 - AU(kSB+8)*X3
        enddo

        WW(jSR  ,indexB)= WVAL1
        WW(jSR+1,indexB)= WVAL2
        WW(jSR+2,indexB)= WVAL3
      enddo

      end subroutine hecmw_solver_init_resid_33

!C
!C***
!C*** hecmw_solver_matvec_33
!C***
!C
      subroutine hecmw_solver_matvec_33 (indexA, indexB)
      implicit REAL*8 (A-H,O-Z)      
      integer(kind=kint) :: indexA,indexB

      S_TIME= HECMW_WTIME()
      call HECMW_SOLVE_SEND_RECV_33                                     &
     &   ( NP,  NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,            &
     &     STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,indexA),              &
     &     SOLVER_COMM,my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      do j= 1, N
        jSR= 3*j-2
        jSB= 9*j-8
           X1= WW(jSR  ,indexA)
           X2= WW(jSR+1,indexA)
           X3= WW(jSR+2,indexA)
        WVAL1= D(jSB  )*X1 + D(jSB+1)*X2 + D(jSB+2)*X3
        WVAL2= D(jSB+3)*X1 + D(jSB+4)*X2 + D(jSB+5)*X3
        WVAL3= D(jSB+6)*X1 + D(jSB+7)*X2 + D(jSB+8)*X3
        do k= INL(j-1)+1, INL(j)
           ik= IAL(k)
          kSR= 3*ik-2
          kSB= 9* k-8
          X1= WW(kSR  ,indexA)
          X2= WW(kSR+1,indexA)
          X3= WW(kSR+2,indexA)
          WVAL1= WVAL1 + AL(kSB  )*X1 + AL(kSB+1)*X2 + AL(kSB+2)*X3
          WVAL2= WVAL2 + AL(kSB+3)*X1 + AL(kSB+4)*X2 + AL(kSB+5)*X3
          WVAL3= WVAL3 + AL(kSB+6)*X1 + AL(kSB+7)*X2 + AL(kSB+8)*X3
        enddo
        do k= INU(j-1)+1, INU(j)
           ik= IAU(k)
          kSR= 3*ik-2
          kSB= 9* k-8
          X1= WW(kSR  ,indexA)
          X2= WW(kSR+1,indexA)
          X3= WW(kSR+2,indexA)
          WVAL1= WVAL1 + AU(kSB  )*X1 + AU(kSB+1)*X2 + AU(kSB+2)*X3
          WVAL2= WVAL2 + AU(kSB+3)*X1 + AU(kSB+4)*X2 + AU(kSB+5)*X3
          WVAL3= WVAL3 + AU(kSB+6)*X1 + AU(kSB+7)*X2 + AU(kSB+8)*X3
        enddo

        WW(jSR  ,indexB)= WVAL1
        WW(jSR+1,indexB)= WVAL2
        WW(jSR+2,indexB)= WVAL3
      enddo

      end subroutine hecmw_solver_matvec_33

!C
!C***
!C*** hecmw_solver_precond_33
!C***
!C
      subroutine hecmw_solver_precond_33 (indexA, indexB, indexC)
      implicit REAL*8 (A-H,O-Z)      
      integer(kind=kint) :: indexA,indexB,indexC

      if (PRECOND.eq.1.or.PRECOND.eq.2) then
!C
!C== Block SSOR

      do j= 1, N
        jSR= 3*j-2
        WW(jSR  ,indexA)= WW(jSR  ,indexC)
        WW(jSR+1,indexA)= WW(jSR+1,indexC)
        WW(jSR+2,indexA)= WW(jSR+2,indexC)
      enddo

      do j= 1, NP
        jSR= 3*j-2
        WW(jSR  ,indexB )= 0.d0
        WW(jSR+1,indexB )= 0.d0
        WW(jSR+2,indexB )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax
        do j= 1+N, NP
          jSR= 3*j-2
          WW(jSR  ,indexA)= 0.d0
          WW(jSR+1,indexA)= 0.d0
          WW(jSR+2,indexA)= 0.d0
        enddo

!C
!C-- FORWARD
        do j= 1, N
          jSR= 3*j-2
          jSB= 9*j-8
          SW1= WW(jSR  ,indexA)
          SW2= WW(jSR+1,indexA)
          SW3= WW(jSR+2,indexA)
          isL= INL(j-1)+1
          ieL= INL(j)
          do k= isL, ieL
             ik= IAL(k)
            kSR= 3*ik-2
            kSB= 9* k-8
             X1= WW(kSR  ,indexA)
             X2= WW(kSR+1,indexA)
             X3= WW(kSR+2,indexA)
            SW1= SW1 - AL(kSB  )*X1 - AL(kSB+1)*X2 - AL(kSB+2)*X3 
            SW2= SW2 - AL(kSB+3)*X1 - AL(kSB+4)*X2 - AL(kSB+5)*X3 
            SW3= SW3 - AL(kSB+6)*X1 - AL(kSB+7)*X2 - AL(kSB+8)*X3 
          enddo

          X1= SW1
          X2= SW2
          X3= SW3

          X2= X2 - ALU(jSB+3)*X1
          X3= X3 - ALU(jSB+6)*X1 - ALU(jSB+7)*X2
          X3= ALU(jSB+8)*  X3
          X2= ALU(jSB+4)*( X2 - ALU(jSB+5)*X3 )
          X1= ALU(jSB  )*( X1 - ALU(jSB+2)*X3 - ALU(jSB+1)*X2)

          WW(jSR  ,indexA)= X1
          WW(jSR+1,indexA)= X2
          WW(jSR+2,indexA)= X3
        enddo

!C
!C-- BACKWARD

        do j= N, 1, -1
          jSR= 3*j-2
          jSB= 9*j-8

          isU= INU(j-1) + 1
          ieU= INU(j) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do k= isU, ieU
             ik= IAU(k)
            kSR= 3*ik-2
            kSB= 9* k-8
             X1= WW(kSR  ,indexA)
             X2= WW(kSR+1,indexA)
             X3= WW(kSR+2,indexA)
            SW1= SW1 + AU(kSB  )*X1 + AU(kSB+1)*X2 + AU(kSB+2)*X3
            SW2= SW2 + AU(kSB+3)*X1 + AU(kSB+4)*X2 + AU(kSB+5)*X3
            SW3= SW3 + AU(kSB+6)*X1 + AU(kSB+7)*X2 + AU(kSB+8)*X3
          enddo

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(jSB+3)*X1
          X3= X3 - ALU(jSB+6)*X1 - ALU(jSB+7)*X2
          X3= ALU(jSB+8)*  X3
          X2= ALU(jSB+4)*( X2 - ALU(jSB+5)*X3 )
          X1= ALU(jSB  )*( X1 - ALU(jSB+2)*X3 - ALU(jSB+1)*X2)
          WW(jSR  ,indexA)=  WW(jSR  ,indexA) - X1
          WW(jSR+1,indexA)=  WW(jSR+1,indexA) - X2
          WW(jSR+2,indexA)=  WW(jSR+2,indexA) - X3
        enddo

!C
!C-- additive Schwartz

      S_TIME= HECMW_WTIME()
        call HECMW_SOLVE_SEND_RECV_33                                   &
     &     ( NP , NEIBPETOT, NEIBPE, STACK_IMPORT, NOD_IMPORT,          &
     &       STACK_EXPORT, NOD_EXPORT, WS, WR, WW(1,indexA) ,           &
     &       SOLVER_COMM,  my_rank)
      E_TIME= HECMW_WTIME()
      COMMtime = COMMtime + E_TIME - S_TIME

      do j= 1, NP
        jSR=  3*j- 2
        WW(jSR  ,indexB)= WW(jSR  ,indexB) + WW(jSR  ,indexA)         
        WW(jSR+1,indexB)= WW(jSR+1,indexB) + WW(jSR+1,indexA)         
        WW(jSR+2,indexB)= WW(jSR+2,indexB) + WW(jSR+2,indexA)         
      enddo

      if (iterPRE.eq.iterPREmax) goto 760

      do j= 1, N
        jSR= 3*j-2
        jSB= 9*j-8
           X1= WW(jSR  ,indexB)
           X2= WW(jSR+1,indexB)
           X3= WW(jSR+2,indexB)
        WVAL1= WW(jSR  ,indexC)                                         &
     &           - D(jSB  )*X1 - D(jSB+1)*X2 - D(jSB+2)*X3           
        WVAL2= WW(jSR+1,indexC)                                         &
     &           - D(jSB+3)*X1 - D(jSB+4)*X2 - D(jSB+5)*X3           
        WVAL3= WW(jSR+2,indexC)                                         &
     &           - D(jSB+6)*X1 - D(jSB+7)*X2 - D(jSB+8)*X3           
        do k= INL(j-1)+1, INL(j)
          ik = IAL(k)
          kSR= 3*ik-2
          kSB= 9* k-8
           X1= WW(kSR  ,indexB)
           X2= WW(kSR+1,indexB)
           X3= WW(kSR+2,indexB)
          WVAL1= WVAL1 - AL(kSB  )*X1 - AL(kSB+1)*X2 - AL(kSB+2)*X3 
          WVAL2= WVAL2 - AL(kSB+3)*X1 - AL(kSB+4)*X2 - AL(kSB+5)*X3 
          WVAL3= WVAL3 - AL(kSB+6)*X1 - AL(kSB+7)*X2 - AL(kSB+8)*X3 
        enddo

        do k= INU(j-1)+1, INU(j)
          ik = IAU(k)
          kSR= 3*ik-2
          kSB= 9* k-8
           X1= WW(kSR  ,indexB)
           X2= WW(kSR+1,indexB)
           X3= WW(kSR+2,indexB)
          WVAL1= WVAL1 - AU(kSB  )*X1 - AU(kSB+1)*X2 - AU(kSB+2)*X3
          WVAL2= WVAL2 - AU(kSB+3)*X1 - AU(kSB+4)*X2 - AU(kSB+5)*X3
          WVAL3= WVAL3 - AU(kSB+6)*X1 - AU(kSB+7)*X2 - AU(kSB+8)*X3
        enddo
        
        WW(jSR  , indexA)= WVAL1
        WW(jSR+1, indexA)= WVAL2
        WW(jSR+2, indexA)= WVAL3
      enddo
 
 760  continue


      enddo
      endif

      if (PRECOND.eq.3) then
!C
!C== Block SCALING

      do j= 1, N
        jSR= 3*j-2
        jSB= 9*j-8

        X1= WW(jSR  ,indexC)
        X2= WW(jSR+1,indexC)
        X3= WW(jSR+2,indexC)

        X2= X2 - ALU(jSB+3)*X1
        X3= X3 - ALU(jSB+6)*X1 - ALU(jSB+7)*X2
        X3= ALU(jSB+8)*  X3
        X2= ALU(jSB+4)*( X2 - ALU(jSB+5)*X3 )
        X1= ALU(jSB  )*( X1 - ALU(jSB+2)*X3 - ALU(jSB+1)*X2)
        WW(jSR  ,indexB)= X1
        WW(jSR+1,indexB)= X2
        WW(jSR+2,indexB)= X3
      enddo
      endif

      end subroutine hecmw_solver_precond_33

      end subroutine  hecmw_solve_GPBiCG_33
      end module     hecmw_solver_GPBiCG_33
