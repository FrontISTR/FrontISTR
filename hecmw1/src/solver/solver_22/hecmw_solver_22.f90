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

      module hecmw_solver_22
      contains
!
!C*** 
!C*** hecmw_solver_22
!C***
!
        subroutine hecmw_solve_22 (hecMESH, hecMAT) 

        use hecmw_util
        use hecmw_solver_CG_22
        use hecmw_solver_BLCG_22
        use m_hecmw_solve_error
        use m_hecmw_comm_f
        use hecmw_matrix_ass

        implicit none

        type (hecmwST_matrix)     :: hecMAT
        type (hecmwST_local_mesh) :: hecMESH

        integer(kind=kint) :: ERROR
        real(kind=kreal), dimension(2,2) :: ALU
        real(kind=kreal), dimension(2)   :: PW

        integer(kind=kint) :: ITER, METHOD, PRECOND, NSET
        integer(kind=kint) :: iterPREmax, ii, i, j, k, l

        integer(kind=kint) :: ITERlog, TIMElog
        real(kind=kreal)::TIME_setup,TIME_comm,TIME_soltot,TIME_sol,TR,THRESH
        real(kind=kreal),    dimension(1) :: RHS
        integer (kind=kint), dimension(1) :: IFLAG

        real(kind=kreal)::RESID,SIGMA_DIAG,SIGMA,FILTER,ALO
!C
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
        ITER      = hecMAT%Iarray(1)
        METHOD    = hecMAT%Iarray(2)
        PRECOND   = hecMAT%Iarray(3)
        NSET      = hecMAT%Iarray(4)
        iterPREmax= hecMAT%Iarray(5)

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

        if (iterPREmax.lt.1) iterPREmax= 1
        if (iterPREmax.gt.4) iterPREmax= 4
!C===

!C
!C +-------------+
!C | ERROR CHECK |
!C +-------------+
!C===
        ERROR= 0

!C
!C-- ZERO RHS norm
        RHS(1)= 0.d0
        do i= 1, hecMAT%N
          RHS(1)= RHS(1) + hecMAT%B(2*i-1)**2 + hecMAT%B(2*i)**2
        enddo
        call hecmw_allreduce_R (hecMESH, RHS, 1, hecmw_sum)

        if (RHS(1).eq.0.d0) then
          ERROR= 2002
          call hecmw_solve_error (hecMESH, ERROR)
        endif

!C
!C-- ZERO DIAGONAL component
        IFLAG(1)= 0
        do i= 1, hecMAT%N
          if (dabs(hecMAT%D(4*i-3)).eq.0.d0) IFLAG(1)= 1
          if (dabs(hecMAT%D(4*i  )).eq.0.d0) IFLAG(1)= 1
        enddo

        call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
        if (IFLAG(1).ne.0) then
          ERROR= 2001
          call hecmw_solve_error (hecMESH, ERROR)
        endif

!C
!C-- INCONSISTENT SOLVER/PRECONDITIONING
        IFLAG(1)= 0
        if (METHOD.le.0  .or. METHOD .ge. 2)  IFLAG(1)= 1
        if (PRECOND.le.0 .or. PRECOND.gt.12) IFLAG(1)= 1

        call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
        if (IFLAG(1).ne.0) then
          ERROR= 1001
          call hecmw_solve_error (hecMESH, ERROR)
        endif

        IFLAG(1)= 1
        if (METHOD.eq.1) then
          if (PRECOND.eq. 1) IFLAG(1)= 0
          if (PRECOND.eq. 2) IFLAG(1)= 0
          if (PRECOND.eq. 3) IFLAG(1)= 0
          if (PRECOND.eq.10) IFLAG(1)= 0
          if (PRECOND.eq.11) IFLAG(1)= 0
          if (PRECOND.eq.12) IFLAG(1)= 0
        endif

        if (IFLAG(1).ne.0) then
          ERROR= 1001
          call hecmw_solve_error (hecMESH, ERROR)
        endif
!C===

!C
!C +-----------+
!C | BLOCK LUs |
!C +-----------+
!C===
        if (.not.associated(hecMAT%ALU).and.PRECOND.lt.10) NSET=  0
        if (     associated(hecMAT%ALU).and.PRECOND.lt.10) NSET= -1

        if (NSET.eq.0.and. PRECOND.lt.10) then
          allocate (hecMAT%ALU(4*hecMAT%N))
        endif

        if (NSET.le.0 .and. PRECOND.lt.10) then
          hecMAT%ALU  = 0.d0

          do ii= 1, hecMAT%N
            ALU(1,1)= hecMAT%D(4*ii-3)
            ALU(1,2)= hecMAT%D(4*ii-2)
            ALU(2,1)= hecMAT%D(4*ii-1)
            ALU(2,2)= hecMAT%D(4*ii  )

            do k= 1, 2
               L = k
              ALO= dabs(ALU(L,k))
              do i= k+1, 2
                if (dabs(ALU(i,k)).gt.ALO) then
                   L = i
                  ALO= dabs(ALU(L,k))
                endif
              enddo

              ALU(k,k)= 1.d0/ALU(k,k)
              do i= k+1, 2
                ALU(i,k)= ALU(i,k) * ALU(k,k)
                do j= k+1, 2
                  PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
                enddo
                do j= k+1, 2
                  ALU(i,j)= PW(j)
                enddo
              enddo
            enddo
            hecMAT%ALU(4*ii-3)= ALU(1,1)
            hecMAT%ALU(4*ii-2)= ALU(1,2)
            hecMAT%ALU(4*ii-1)= ALU(2,1)
            hecMAT%ALU(4*ii  )= ALU(2,2)
          enddo
        endif
!C===

!C
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===
!      if (hecMESH%my_rank.eq.0) write (*,'(a)') ' '

      if (METHOD.eq.1 .and. PRECOND.lt.10) then
        ! imposing MPC by penalty
        call hecmw_mat_ass_equation ( hecMESH, hecMAT )

        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.eq.1)) then
          if (PRECOND.eq.1) write (*,'(a)') '### 2x2 B-IC-CG  (0)'
          if (PRECOND.eq.2) write (*,'(a)') '### 2x2 B-SSOR-CG(0)'
          if (PRECOND.eq.3) write (*,'(a)') '### 2x2 B-scale-CG  '
        endif

        call hecmw_solve_CG_22                                          &
     &     ( hecMAT%N, hecMAT%NP, hecMAT%NPL, hecMAT%NPU,               &
     &       hecMAT%D, hecMAT%AL, hecMAT%indexL, hecMAT%itemL,          &
     &                 hecMAT%AU, hecMAT%indexU, hecMAT%itemU,          &
     &       hecMAT%B, hecMAT%X,  hecMAT%ALU, RESID, ITER,              &
     &       ERROR,    hecMESH%my_rank,                                 &
     &       hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                &
     &       hecMESH%import_index, hecMESH%import_item,                 & 
     &       hecMESH%export_index, hecMESH%export_item,                 & 
     &       hecMESH%MPI_COMM, PRECOND, iterPREmax,                     &
     &       TIME_setup, TIME_sol, TIME_comm, ITERlog) 
      endif

      if (METHOD.eq.1 .and. PRECOND.ge.10) then
        ! imposing MPC by penalty
        call hecmw_mat_ass_equation ( hecMESH, hecMAT )

        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.eq.1)) then
          if (PRECOND.eq.10) write (*,'(a)') '### 2x2 B-IC-CG  (0)'
          if (PRECOND.eq.11) write (*,'(a)') '### 2x2 B-IC-CG  (1)'
          if (PRECOND.eq.12) write (*,'(a)') '### 2x2 B-IC-CG  (2)'
        endif

        SIGMA     = 1.d0
        SIGMA_DIAG= 1.d0
        call hecmw_solve_BLCG_22                                        &
     &     ( hecMAT%N, hecMAT%NP, hecMAT%NPL, hecMAT%NPU,               &
     &       hecMAT%D, hecMAT%AL, hecMAT%indexL, hecMAT%itemL,          &
     &                 hecMAT%AU, hecMAT%indexU, hecMAT%itemU,          &
     &       hecMAT%B, hecMAT%X,  RESID, SIGMA, SIGMA_DIAG,             &
     &       ITER, ERROR,  hecMESH%my_rank,                             &
     &       hecMESH%n_neighbor_pe, hecMESH%neighbor_pe,                &
     &       hecMESH%import_index, hecMESH%import_item,                 & 
     &       hecMESH%export_index, hecMESH%export_item,                 & 
     &       hecMESH%MPI_COMM, PRECOND, iterPREmax,                     &
     &       TIME_setup, TIME_sol, TIME_comm, ITERlog) 
      endif

      !hecMAT%ITERactual = ITER
      !hecMAT%RESIDactual= RESID

      if (RESID.gt.hecMAT%Rarray(1)) then
        call hecmw_solve_error (hecMESH, 3001)
      endif

      if (hecMESH%my_rank.eq.0 .and. TIMElog.eq.1) then
        TR= (TIME_sol-TIME_comm)/(TIME_sol+1.d-24)*100.d0
        write (*,'(/a)')          '### summary of linear solver'
        write (*,'(i10,a, 1pe16.6)')      ITER, ' iterations  ', RESID
        write (*,'(a, 1pe16.6 )') '    set-up time      : ', TIME_setup
        write (*,'(a, 1pe16.6 )') '    solver time      : ', TIME_sol
        write (*,'(a, 1pe16.6 )') '    solver/comm time : ', TIME_comm
        write (*,'(a, 1pe16.6/)') '    work ratio (%)   : ', TR
      endif
!C===

      end subroutine hecmw_solve_22
      end module hecmw_solver_22
