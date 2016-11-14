!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C
!C***
!C*** module hecmw_solver_CG_44
!C***
!C
      module hecmw_solver_CG_44
      contains
!C
!C*** CG_44
!C
      subroutine hecmw_solve_CG_44( hecMESH,  hecMAT, ITER, RESID, ERROR, &
     &                              Tset, Tsol, Tcomm )

      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_matrix_misc
      use hecmw_solver_misc
      use hecmw_solver_las_44
      use hecmw_solver_scaling_44
      use hecmw_precond_44
      use hecmw_jad_type_44

      implicit none

      type(hecmwST_local_mesh) :: hecMESH
      type(hecmwST_matrix) :: hecMAT
      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

      integer(kind=kint ) :: N, NP, NDOF, NNDOF
      integer(kind=kint ) :: my_rank
      integer(kind=kint ) :: ITERlog, TIMElog
      real(kind=kreal), pointer :: B(:), X(:)

      real(kind=kreal), dimension(:,:), allocatable :: WW

      integer(kind=kint), parameter ::  R= 1
      integer(kind=kint), parameter ::  Z= 2
      integer(kind=kint), parameter ::  Q= 2
      integer(kind=kint), parameter ::  P= 3
      integer(kind=kint), parameter :: BT= 1
      integer(kind=kint), parameter ::TATX=2
      integer(kind=kint), parameter :: WK= 4

      integer(kind=kint ) :: MAXIT
      integer(kind=kint ) :: totalmpc

! local variables
      real   (kind=kreal) :: TOL
      integer(kind=kint )::i
      real   (kind=kreal)::S_TIME,S1_TIME,E_TIME,E1_TIME, START_TIME, END_TIME
      real   (kind=kreal)::BNRM2
      real   (kind=kreal)::RHO,RHO1,BETA,C1,ALPHA,DNRM2

      S_TIME= HECMW_WTIME()

!C===
!C +-------+
!C | INIT. |
!C +-------+
!C===
      write (*,'(a)') '### start'
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

      totalmpc = hecMESH%mpc%n_mpc
      call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)
      write (*,'(a)') '### hecmw_allreduceI1_finished'
      ERROR = 0

      allocate (WW(NDOF*NP, 4))
      WW = 0.d0

      call hecmw_mpc_scale_44(hecMESH)
      write (*,'(a)') '### hecmw_mpc_scale_finished'

!C
!C-- SCALING
      call hecmw_solver_scaling_fw_44(hecMESH, hecMAT, Tcomm)
      write (*,'(a)') '### hecmw_scaling_finished'

      IF (hecmw_mat_get_usejad(hecMAT).ne.0) THEN
        call hecmw_JAD_INIT_44(hecMAT)
      ENDIF

!C===
!C +----------------------+
!C | SETUP PRECONDITIONER |
!C +----------------------+
!C===
      call hecmw_precond_44_setup(hecMAT)
      write (*,'(a)') '### hecmw_precond_44_setup_finished'

!C===
!C +----------------------------------------------+
!C | {r0}= [T']({b} - [A]{xg}) - [T'][A][T]{xini} |
!C +----------------------------------------------+
!C===

!C-- {bt}= [T']({b} - [A]{xg})
      if (totalmpc.eq.0) then
        do i = 1, NNDOF
          WW(i,BT) = B(i)
        enddo
       else
        call hecmw_trans_b_44(hecMESH, hecMAT, B, WW(:,BT), WW(:,WK), Tcomm)
        write (*,'(a)') '### hecmw_trans_b_44_finished'
      endif

!C-- compute ||{bt}||
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,BT), WW(:,BT), BNRM2, Tcomm)
      write (*,'(a)') '### hecmw_InnerProduct_R_finished'
      if (BNRM2.eq.0.d0) then
        iter = 0
        MAXIT = 0
        RESID = 0.d0
        X = 0.d0
      endif

!C
!C-- {tatx} = [T'] [A] [T]{x}
      if (totalmpc.eq.0) then
        call hecmw_matvec_44(hecMESH, hecMAT, X, WW(:,TATX), Tcomm)
        write (*,'(a)') '### hecmw_matvec_44_finished'
       else
        call hecmw_TtmatTvec_44(hecMESH, hecMAT, X, WW(:,TATX), WW(:,WK), Tcomm)
        write (*,'(a)') '### hecmw_TtmatTvec_44_finished'
      endif

!C-- {r} = {bt} - {tatx}
      do i = 1, NNDOF
        WW(i,R) = WW(i,BT) - WW(i,TATX)
      enddo

      E_TIME = HECMW_WTIME()
      Tset = Tset + E_TIME - S_TIME

      Tcomm = 0.d0
      S1_TIME = HECMW_WTIME()
!C
!C************************************************* Conjugate Gradient Iteration start
!C
      do iter = 1, MAXIT

!C===
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===
      call hecmw_precond_44_apply(hecMESH, hecMAT, WW(:,R), WW(:,Z), WW(:,WK), Tcomm)
      write (*,'(a)') '### hecmw_precond_44_apply_finished'

!C===
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,Z), RHO, Tcomm)
      write (*,'(a)') '### hecmw_InnerProduct_R_finished'

!C===
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then
        do i = 1, NNDOF
          WW(i,P) = WW(i,Z)
        enddo
       else
         BETA = RHO / RHO1
         do i = 1, NNDOF
           WW(i,P) = WW(i,Z) + BETA*WW(i,P)
         enddo
      endif

!C===
!C +---------------------+
!C | {q}= [T'][A][T] {p} |
!C +---------------------+
!C===
      if (totalmpc.eq.0) then
        call hecmw_matvec_44(hecMESH, hecMAT, WW(:,P), WW(:,Q), Tcomm)
        write (*,'(a)') '### hecmw_matvec_44_finished'
       else
        call hecmw_TtmatTvec_44(hecMESH, hecMAT, WW(:,P), WW(:,Q), WW(:,WK), Tcomm)
        write (*,'(a)') '### hecmw_TtmatTvec_44_finished'
      endif

!C===
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,P), WW(:,Q), C1, Tcomm)
      write (*,'(a)') '### hecmw_InnerProduct_R_finished'

      ALPHA= RHO / C1

!C===
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===
      do i = 1, NNDOF
         X(i)  = X(i)    + ALPHA * WW(i,P)
        WW(i,R)= WW(i,R) - ALPHA * WW(i,Q)
      enddo

      call hecmw_InnerProduct_R(hecMESH, NDOF, WW(:,R), WW(:,R), DNRM2, Tcomm)
      write (*,'(a)') '### hecmw_InnerProduct_R_finished'

      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
      if (my_rank.eq.0.and.ITERLog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID
!C#####

      if ( RESID.le.TOL   ) exit
      if ( ITER .eq.MAXIT ) ERROR = -300

      RHO1 = RHO

      enddo
!C
!C************************************************* Conjugate Gradient Iteration end
!C
      if (totalmpc.ne.0) then
        call hecmw_tback_x_44(hecMESH, X, WW(:,WK), Tcomm)
        write (*,'(a)') '### hecmw_tback_x_44_finished'
      endif

      call hecmw_solver_scaling_bk_44(hecMAT)
      write (*,'(a)') '### hecmw_solver_scaling_bk_44_finished'
!C
!C-- INTERFACE data EXCHANGE
!C
      START_TIME= HECMW_WTIME()
      call hecmw_update_4_R (hecMESH, X, hecMAT%NP)
      write (*,'(a)') '### hecmw_update_4_R_finished'
      END_TIME = HECMW_WTIME()
      Tcomm = Tcomm + END_TIME - START_TIME

      E1_TIME = HECMW_WTIME()
      Tsol = E1_TIME - S1_TIME

      deallocate (WW)
      call hecmw_precond_44_clear(hecMAT)
      write (*,'(a)') '### hecmw_precond_44_clear_finished'

      IF (hecmw_mat_get_usejad(hecMAT).ne.0) THEN
        call hecmw_JAD_FINALIZE_44()
        write (*,'(a)') '### hecmw_JAD_FINALIZE_44_finished'
      ENDIF

      end subroutine hecmw_solve_CG_44
      end module     hecmw_solver_CG_44
