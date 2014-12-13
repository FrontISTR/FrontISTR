!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.5                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                       Kazuya Goto (PExProCS LLC)                     !
!                       Naoki Morita (Univ. of Tokyo)                  !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

      module hecmw_solver_66
      contains
!
!C***
!C*** hecmw_solve_66
!C***
!
        subroutine hecmw_solve_66 (hecMESH, hecMAT)

        use hecmw_util
        use hecmw_solver_CG_66
        use m_hecmw_solve_error
        use m_hecmw_comm_f
        use hecmw_matrix_ass
        use hecmw_matrix_contact
        use hecmw_solver_las_66
        use hecmw_precond_66

        implicit none

        type (hecmwST_matrix)     :: hecMAT
        type (hecmwST_local_mesh) :: hecMESH

        integer(kind=kint) :: ERROR

        integer(kind=kint) :: ITER, METHOD, PRECOND, NSET
        integer(kind=kint) :: iterPREmax, i, BMP

        real(kind=kreal) :: RESID, SIGMA_DIAG, THRESH, FILTER, residual

        integer(kind=kint) :: ITERlog, TIMElog
        real(kind=kreal) :: TIME_setup, TIME_comm, TIME_soltot, TR
        real(kind=kreal) :: time_Ax, time_precond
        real(kind=kreal),    dimension(1) :: RHS
        integer (kind=kint), dimension(1) :: IFLAG

        integer(kind=kint) :: NREST, mtxmarket, mtx33211
        real(kind=kreal)   :: SIGMA,TIME_sol

		real(kind=kreal), allocatable :: vv(:), zz(:)

!C===
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
        ITER      = hecMAT%Iarray(1)
        METHOD    = hecMAT%Iarray(2)
        PRECOND   = hecMAT%Iarray(3)
        NSET      = hecMAT%Iarray(4)
        iterPREmax= hecMAT%Iarray(5)
        NREST     = hecMAT%Iarray(6)
        BMP       = hecMAT%Iarray(36)
        mtxmarket = hecMAT%Iarray(37)
        mtx33211  = hecMAT%Iarray(38)

        ITERlog= hecMAT%Iarray(21)
        TIMElog= hecMAT%Iarray(22)

        TIME_setup  = 0.d0
        TIME_comm   = 0.d0
        TIME_soltot = 0.d0

        RESID     = hecMAT%Rarray(1)
        SIGMA_DIAG= hecMAT%Rarray(2)
        SIGMA     = hecMAT%Rarray(3)

        THRESH= hecMAT%Rarray(4)
        FILTER= hecMAT%Rarray(5)

        if (iterPREmax.lt.0) iterPREmax= 0
        if (iterPREmax.gt.4) iterPREmax= 4

        if (SIGMA_DIAG.lt.0.5d0) SIGMA_DIAG= 0.5d0
        if (SIGMA_DIAG.gt.2.d0) SIGMA_DIAG= 2.d0

        if (SIGMA.lt.0.d0) SIGMA= 0.d0
        if (SIGMA.gt.1.d0) SIGMA= 1.d0

!C===
!C +-------------+
!C | ERROR CHECK |
!C +-------------+
!C===
        ERROR= 0

!C
!C-- ZERO RHS norm
        RHS(1)= 0.d0
        do i= 1, hecMAT%N
          RHS(1)= RHS(1) + hecMAT%B(6*i-5)**2 + hecMAT%B(6*i-4)**2 + hecMAT%B(6*i-3)**2 &
          & + hecMAT%B(6*i-2)**2 + hecMAT%B(6*i-1)**2 + hecMAT%B(6*i  )**2
        enddo
        if (hecMESH%mpc%n_mpc > 0) then
          do i= 1, hecMESH%mpc%n_mpc
            RHS(1)= RHS(1) + hecMESH%mpc%mpc_const(i)**2
          enddo
        endif
        call hecmw_allreduce_R (hecMESH, RHS, 1, hecmw_sum)

        if (RHS(1).eq.0.d0) then
          ERROR= 2002
          call hecmw_solve_error (hecMESH, ERROR)
        endif

!C
!C-- ZERO DIAGONAL component
        IFLAG(1)= 0
        do i= 1, hecMAT%N
          if (dabs(hecMAT%D(36*i-35)).eq.0.d0) IFLAG(1)= 1
          if (dabs(hecMAT%D(36*i-28)).eq.0.d0) IFLAG(1)= 1
          if (dabs(hecMAT%D(36*i-21)).eq.0.d0) IFLAG(1)= 1
          if (dabs(hecMAT%D(36*i-14)).eq.0.d0) IFLAG(1)= 1
          if (dabs(hecMAT%D(36*i-7 )).eq.0.d0) IFLAG(1)= 1
          if (dabs(hecMAT%D(36*i   )).eq.0.d0) IFLAG(1)= 1
        enddo

        call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
        if (IFLAG(1).ne.0 .and. (PRECOND.le.10 .and. iterPREmax.gt.0)) then
          ERROR= 2001
          call hecmw_solve_error (hecMESH, ERROR)
        endif

!C
!C-- INCONSISTENT SOLVER/PRECONDITIONING
        IFLAG(1)= 0
        if (METHOD.le.0 .or. METHOD.ge.5)    IFLAG(1)= 1
        if (PRECOND.le.0 .or. PRECOND.gt.36) IFLAG(1)= 1

        call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
        if (IFLAG(1).ne.0) then
          ERROR= 1001
          call hecmw_solve_error (hecMESH, ERROR)
        endif

        IFLAG(1)= 1
        if (PRECOND.eq. 1) IFLAG(1)= 0
        if (PRECOND.eq. 2) IFLAG(1)= 0
        if (PRECOND.eq. 3) IFLAG(1)= 0
        if (PRECOND.eq.10) IFLAG(1)= 0
        if (PRECOND.eq.11) IFLAG(1)= 0
        if (PRECOND.eq.12) IFLAG(1)= 0

        if (IFLAG(1).ne.0) then
          ERROR= 1001
          call hecmw_solve_error (hecMESH, ERROR)
        endif

        call hecmw_matvec_66_clear_timer()
        call hecmw_precond_66_clear_timer()

!C===
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===

!C
!C-- CG
      if (METHOD.eq.1) then
        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.eq.1)) then
          if (PRECOND.eq.1) &
               write (*,'(a,i3)') '### 6x6 B-IC-CG(0)',iterPREmax
          if (PRECOND.eq.2) &
               write (*,'(a,i3)') '### 6x6 B-SSOR-CG(0)',iterPREmax
          if (PRECOND.eq.3) &
               write (*,'(a,i3)') '### 6x6 B-scale-CG  ',iterPREmax
          if (PRECOND.eq.10) &
               write (*,'(a,i3)') '### 6x6 B-IC-CG(0)',iterPREmax
          if (PRECOND.eq.11) &
               write (*,'(a,i3)') '### 6x6 B-IC-CG(1)',iterPREmax
          if (PRECOND.eq.12) &
               write (*,'(a,i3)') '### 6x6 B-IC-CG(2)',iterPREmax
        endif
        call hecmw_solve_CG_66( hecMESH,  hecMAT, ITER, RESID, ERROR,   &
     &                          TIME_setup, TIME_sol, TIME_comm )
      endif

      if (RESID.gt.hecMAT%Rarray(1)) then
        call hecmw_solve_error (hecMESH, 3001)
      endif

      time_Ax = hecmw_matvec_66_get_timer()
      time_precond = hecmw_precond_66_get_timer()

      if (hecMESH%my_rank.eq.0 .and. TIMElog.eq.1) then
        TR= (TIME_sol-TIME_comm)/(TIME_sol+1.d-24)*100.d0
        write (*,'(/a)')          '### summary of linear solver'
        write (*,'(i10,a, 1pe16.6)')      ITER, ' iterations  ', RESID
        write (*,'(a, 1pe16.6 )') '    set-up time      : ', TIME_setup
        write (*,'(a, 1pe16.6 )') '    solver time      : ', TIME_sol
        write (*,'(a, 1pe16.6 )') '    solver/comm time : ', TIME_comm
        write (*,'(a, 1pe16.6 )') '    solver/matvec    : ', time_Ax
        write (*,'(a, 1pe16.6 )') '    solver/precond   : ', time_precond
        write (*,'(a, 1pe16.6/)') '    work ratio (%)   : ', TR
      endif

      end subroutine hecmw_solve_66
      end module hecmw_solver_66
