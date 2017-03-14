!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C***
!C*** module hecmw_solver_GMRES_33
!C***
!
      module hecmw_solver_GMRES_33

      public :: hecmw_solve_GMRES_33

      contains
!C
!C*** hecmw_solve_GMRES_33
!C
      subroutine hecmw_solve_GMRES_33( hecMESH,  hecMAT, ITER, RESID, ERROR, &
     &                                    Tset, Tsol, Tcomm )
      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_matrix_misc
      use hecmw_solver_misc
      use hecmw_solver_las_33
      use hecmw_solver_scaling_33
      use hecmw_precond_33
      use hecmw_estimate_condition

      implicit none

      type(hecmwST_local_mesh) :: hecMESH
      type(hecmwST_matrix) :: hecMAT
      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      integer(kind=kint ) :: N, NP, NDOF, NNDOF
      integer(kind=kint ) :: my_rank
      integer(kind=kint ) :: ITERlog, TIMElog
      real(kind=kreal), pointer :: B(:), X(:)

      real(kind=kreal), dimension(:,:),  allocatable :: WW

      integer(kind=kint ) :: MAXIT, NREST

      real   (kind=kreal) :: TOL

      real   (kind=kreal), dimension(:),   allocatable :: SS
      real   (kind=kreal), dimension(:,:), allocatable :: H

      integer(kind=kint ) :: CS, SN

      real   (kind=kreal)   ZERO, ONE
      parameter ( ZERO = 0.0D+0, ONE = 1.0D+0 )

      integer(kind=kint ) :: NRK,i,k,kk,jj,INFO,ik
      integer(kind=kint ) :: IROW
      real   (kind=kreal) :: S_TIME,E_TIME,S1_TIME,E1_TIME
      real   (kind=kreal) :: LDH,LDW,BNRM2,DNRM2,RNORM
      real   (kind=kreal) :: COMMtime,COMPtime, coef,VAL,VCS,VSN,DTEMP,AA,BB,R0,SCALE,RR
      integer(kind=kint ) :: ESTCOND
      real   (kind=kreal) :: t_max,t_min,t_avg,t_sd

      integer(kind=kint), parameter :: R  = 1
      integer(kind=kint), parameter :: ZP = R + 1
      integer(kind=kint), parameter :: ZQ = R + 2
      integer(kind=kint), parameter :: S  = R + 3
      integer(kind=kint), parameter :: W  = S + 1
      integer(kind=kint), parameter :: Y  = W
      integer(kind=kint), parameter :: AV = Y  + 1
      integer(kind=kint), parameter :: V  = AV + 1

      call hecmw_barrier(hecMESH)
      S_TIME= HECMW_WTIME()
!C
!C-- INIT.
      N = hecMAT%N
      NP = hecMAT%NP
      NDOF = hecMAT%NDOF
      NNDOF = N * NDOF
      my_rank = hecMESH%my_rank
      X => hecMAT%X
      B => hecMAT%B

      ITERlog = hecmw_mat_get_iterlog( hecMAT )
      TIMElog = hecmw_mat_get_timelog( hecMAT )
      MAXIT  = hecmw_mat_get_iter( hecMAT )
       TOL   = hecmw_mat_get_resid( hecMAT )
      NREST  = hecmw_mat_get_nrest( hecMAT )
      ESTCOND = hecmw_mat_get_estcond( hecMAT )

      if (NREST >= NDOF*NP-1) NREST = NDOF*NP-2

      ERROR= 0
      NRK= NREST + 7

      allocate (H (NRK,NRK))
      allocate (WW(NDOF*NP,NRK))
      allocate (SS(NRK))

      COMMtime= 0.d0
      COMPtime= 0.d0

      LDH= NREST + 2
      LDW= N

!C
!C-- Store the Givens parameters in matrix H.
      CS= NREST + 1
      SN= CS    + 1

!C
!C-- SCALING
      call hecmw_solver_scaling_fw_33(hecMESH, hecMAT, Tcomm)

!C===
!C +----------------------+
!C | SETUP PRECONDITIONER |
!C +----------------------+
!C===
      call hecmw_precond_33_setup(hecMAT, hecMESH, 0)

!C
!C
!C +--------------------+
!C | {r}= {b} - [A]{x0} |
!C +--------------------+
!C===
      call hecmw_matresid_33(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)
!C===

      call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
      if (BNRM2.eq.0.d0) then
        iter = 0
        MAXIT = 0
        RESID = 0.d0
        X = 0.d0
      endif

      E_TIME= HECMW_WTIME()
      if (TIMElog.eq.2) then
        call hecmw_time_statistics(hecMESH, E_TIME - S_TIME, &
             t_max, t_min, t_avg, t_sd)
        if (hecMESH%my_rank.eq.0) then
          write(*,*) 'Time solver setup'
          write(*,*) '  Max     :',t_max
          write(*,*) '  Min     :',t_min
          write(*,*) '  Avg     :',t_avg
          write(*,*) '  Std Dev :',t_sd
        endif
        Tset = t_max
      else
        Tset = E_TIME - S_TIME
      endif
!C===

      call hecmw_barrier(hecMESH)
      S1_TIME= HECMW_WTIME()
      ITER= 0

      OUTER: do

!C
!C************************************************ GMRES Iteration
!C
        I= 0
!C
!C +---------------+
!C | {v1}= {r}/|r| |
!C +---------------+
!C===
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
        if (DNRM2 == 0.d0) exit ! converged

        RNORM= dsqrt(DNRM2)
        coef= ONE/RNORM
        do ik= 1, NNDOF
          WW(ik,V)= WW(ik,R) * coef
        enddo
!C===

!C
!C +--------------+
!C | {s}= |r|{e1} |
!C +--------------+
!C===
        WW(1 ,S) = RNORM
        do k = 2, NNDOF
          WW(k,S) = ZERO
        enddo
!C===

!C************************************************ GMRES(m) restart
        do I = 1, NREST
        ITER= ITER + 1

!C
!C +-------------------+
!C | {w}= [A][Minv]{v} |
!C +-------------------+
!C===
      call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,V+I-1), WW(:,ZQ), WW(:,ZP), Tcomm)

      call hecmw_matvec_33(hecMESH, hecMAT, WW(:,ZQ), WW(:,W), Tcomm)
!C===

!C
!C +------------------------------+
!C | ORTH. BASIS by GRAMM-SCHMIDT |
!C +------------------------------+
!C   Construct the I-th column of the upper Hessenberg matrix H
!C   using the Gram-Schmidt process on V and W.
!C===
      do K= 1, I
        call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,W), WW(:,V+K-1), VAL, Tcomm)

        do ik= 1, NNDOF
          WW(ik,W)= WW(ik,W) - VAL * WW(ik,V+K-1)
        enddo
        H(K,I)= VAL
      enddo

      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,W), WW(:,W), VAL, Tcomm)
      if (VAL == 0.d0) exit ! converged

      H(I+1,I)= dsqrt(VAL)
      coef= ONE / H(I+1,I)
      do ik= 1, NNDOF
        WW(ik,V+I+1-1)= WW(ik,W) * coef
      enddo
!C===

!C
!C +-----------------+
!C | GIVENS ROTARION |
!C +-----------------+
!C===

!C
!C-- Plane Rotation
      do k = 1, I-1
        VCS= H(k,CS)
        VSN= H(k,SN)
        DTEMP   = VCS*H(k  ,I) + VSN*H(k+1,I)
        H(k+1,I)= VCS*H(k+1,I) - VSN*H(k  ,I)
        H(k  ,I)= DTEMP
      enddo

!C
!C-- Construct Givens Plane Rotation
      AA = H(I  ,I)
      BB = H(I+1,I)
      R0= BB
      if (dabs(AA).gt.dabs(BB)) R0= AA
      SCALE= dabs(AA) + dabs(BB)

      if (SCALE.ne.0.d0) then
        RR= SCALE * dsqrt((AA/SCALE)**2+(BB/SCALE)**2)
        RR= dsign(1.d0,R0)*RR
        H(I,CS)= AA/RR
        H(I,SN)= BB/RR
       else
        H(I,CS)= 1.d0
        H(I,SN)= 0.d0
        RR     = 0.d0
      endif

!C
!C-- Plane Rotation
      VCS= H(I,CS)
      VSN= H(I,SN)
      DTEMP    = VCS*H(I  ,I) + VSN*H(I+1,I)
      H (I+1,I)= VCS*H(I+1,I) - VSN*H(I  ,I)
      H (I  ,I)= DTEMP

      DTEMP    = VCS*WW(I  ,S) + VSN*WW(I+1,S)
      WW(I+1,S)= VCS*WW(I+1,S) - VSN*WW(I  ,S)
      WW(I  ,S)= DTEMP

      RESID = dabs ( WW(I+1,S))/dsqrt(BNRM2)

      if (my_rank.eq.0 .and. ITERlog.eq.1)                              &
     &    write (*, '(2i8, 1pe16.6)') iter,I+1, RESID

      if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
        if (mod(ITER,ESTCOND) == 0) call hecmw_estimate_condition_GMRES(I, H)
      endif

      if ( RESID.le.TOL ) then
!C-- [H]{y}= {s_tld}
         do ik= 1, I
           SS(ik)= WW(ik,S)
         enddo
         IROW= I
         WW(IROW,Y)= SS(IROW) / H(IROW,IROW)

         do kk= IROW-1, 1, -1
           do jj= IROW, kk+1, -1
             SS(kk)= SS(kk) - H(kk,jj)*WW(jj,Y)
           enddo
           WW(kk,Y)= SS(kk) / H(kk,kk)
         enddo

!C-- {x}= {x} + {y}{V}
         do kk= 1, NNDOF
           WW(kk, AV)= 0.d0
         enddo

         jj= IROW
         do jj= 1, IROW
         do kk= 1, NNDOF
           WW(kk,AV)= WW(kk,AV) + WW(jj,Y)*WW(kk,V+jj-1)
         enddo
         enddo

         call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,AV), WW(:,ZQ), WW(:,ZP), Tcomm)

         do kk= 1, NNDOF
           X(kk)= X(kk) + WW(kk,ZQ)
         enddo

         exit OUTER
      endif

      if ( ITER.gt.MAXIT ) then
        ERROR = HECMW_SOLVER_ERROR_NOCONV_MAXIT
        exit OUTER
      end if
      end do
!C===

!C
!C +------------------+
!C | CURRENT SOLUTION |
!C +------------------+
!C===

!C-- [H]{y}= {s_tld}
      do ik= 1, NREST
        SS(ik)= WW(ik,S)
      enddo
      IROW= NREST
      WW(IROW,Y)= SS(IROW) / H(IROW,IROW)

      do kk= IROW-1, 1, -1
        do jj= IROW, kk+1, -1
          SS(kk)= SS(kk) - H(kk,jj)*WW(jj,Y)
        enddo
        WW(kk,Y)= SS(kk) / H(kk,kk)
      enddo

!C-- {x}= {x} + {y}{V}
      do kk= 1, NNDOF
        WW(kk, AV)= 0.d0
      enddo

      jj= IROW
      do jj= 1, IROW
      do kk= 1, NNDOF
        WW(kk,AV)= WW(kk,AV) + WW(jj,Y)*WW(kk,V+jj-1)
      enddo
      enddo

      call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,AV), WW(:,ZQ), WW(:,ZP), Tcomm)

      do kk= 1, NNDOF
        X(kk)= X(kk) + WW(kk,ZQ)
      enddo

!C
!C-- Compute residual vector R, find norm, then check for tolerance.
      call hecmw_matresid_33(hecMESH, hecMAT, X, B, WW(:,R), Tcomm)

      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)

      WW(I+1,S)= dsqrt(DNRM2/BNRM2)
      RESID    = WW( I+1,S )

      if ( RESID.le.TOL )   exit OUTER
      if ( ITER .gt.MAXIT ) then
        ERROR = HECMW_SOLVER_ERROR_NOCONV_MAXIT
        exit OUTER
      end if
!C
!C-- RESTART
      end do OUTER

!C
!C-- iteration FAILED

      if (ERROR == HECMW_SOLVER_ERROR_NOCONV_MAXIT) then
        INFO = ITER

        !C-- [H]{y}= {s_tld}
        do ik= 1, I
          SS(ik)= WW(ik,S)
        enddo
        IROW= I
        WW(IROW,Y)= SS(IROW) / H(IROW,IROW)

        do kk= IROW-1, 1, -1
          do jj= IROW, kk+1, -1
            SS(kk)= SS(kk) - H(kk,jj)*WW(jj,Y)
          enddo
          WW(kk,Y)= SS(kk) / H(kk,kk)
        enddo

        !C-- {x}= {x} + {y}{V}
        do kk= 1, NNDOF
          WW(kk, AV)= 0.d0
        enddo

        jj= IROW
        do jj= 1, IROW
          do kk= 1, NNDOF
            WW(kk,AV)= WW(kk  ,AV) + WW(jj,Y)*WW(kk  ,V+jj-1)
          enddo
        enddo

        call hecmw_precond_33_apply(hecMESH, hecMAT, WW(:,AV), WW(:,ZQ), WW(:,ZP), Tcomm)

        do kk= 1, NNDOF
          X(kk)= X(kk) + WW(kk,ZQ)
        enddo
      end if

      call hecmw_solver_scaling_bk_33(hecMAT)

      if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
        call hecmw_estimate_condition_GMRES(I, H)
      endif
!C
!C-- INTERFACE data EXCHANGE
      S_TIME = HECMW_WTIME()
      call hecmw_update_3_R (hecMESH, X, hecMAT%NP)
      E_TIME = HECMW_WTIME()
      Tcomm = Tcomm + E_TIME - S_TIME

      deallocate (H, WW, SS)
      !call hecmw_precond_33_clear(hecMAT)

      E1_TIME= HECMW_WTIME()
      if (TIMElog.eq.2) then
        call hecmw_time_statistics(hecMESH, E1_TIME - S1_TIME, &
             t_max, t_min, t_avg, t_sd)
        if (hecMESH%my_rank.eq.0) then
          write(*,*) 'Time solver iterations'
          write(*,*) '  Max     :',t_max
          write(*,*) '  Min     :',t_min
          write(*,*) '  Avg     :',t_avg
          write(*,*) '  Std Dev :',t_sd
        endif
        Tsol = t_max
      else
        Tsol = E1_TIME - S1_TIME
      endif

      end subroutine  hecmw_solve_GMRES_33

      end module     hecmw_solver_GMRES_33
