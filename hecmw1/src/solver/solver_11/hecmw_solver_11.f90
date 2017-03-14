!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!
!C***
!C*** module hecmw_solver_11
!C***
!

      module hecmw_solver_11
      contains

        subroutine hecmw_solve_11 (hecMESH, hecMAT)

        use hecmw_util
        use hecmw_solver_CG_11
        use m_hecmw_solve_error
        use m_hecmw_comm_f
        use hecmw_matrix_ass
        use hecmw_matrix_dump

        implicit REAL*8 (A-H,O-Z)

        type (hecmwST_matrix)     :: hecMAT
        type (hecmwST_local_mesh) :: hecMESH

        integer(kind=kint) :: ERROR, ICFLAG, PRECOND
        character(len= hecmw_name_len) :: BUF
        data ICFLAG/0/

        integer(kind=kint) :: ITERlog, TIMElog
        real(kind=kreal) :: TIME_setup, TIME_comm, TIME_soltot, TR
        real(kind=kreal),    dimension(1) :: RHS
        integer (kind=kint), dimension(1) :: IFLAG

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
          RHS(1)= RHS(1) + hecMAT%B(i)**2
        enddo
        call hecmw_allreduce_R (hecMESH, RHS, 1, hecmw_sum)

        if (RHS(1).eq.0.d0) then
          ERROR= HECMW_SOLVER_ERROR_ZERO_RHS
          call hecmw_solve_error (hecMESH, ERROR)
          hecMAT%X(:)=0.d0
          return
        endif

!C
!C-- ZERO DIAGONAL component
        IFLAG= 0
        do i= 1, hecMAT%N
          if (dabs(hecMAT%D(i)).eq.0.d0) IFLAG= 1
        enddo

        call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
        if (IFLAG(1).ne.0) then
          ERROR= HECMW_SOLVER_ERROR_ZERO_DIAG
          call hecmw_solve_error (hecMESH, ERROR)
        endif

!C
!C-- INCONSISTENT SOLVER/PRECONDITIONING
        IFLAG= 0
        if (METHOD.le.0  .or. METHOD .ge. 2)  IFLAG(1)= 1
        if (PRECOND.le.0 .or. PRECOND.gt.3)   IFLAG(1)= 1

        call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
        if (IFLAG(1).ne.0) then
          ERROR= HECMW_SOLVER_ERROR_INCONS_PC
          call hecmw_solve_error (hecMESH, ERROR)
        endif

        IFLAG= 1
        if (METHOD.eq.1) then
          if (PRECOND.eq. 1) IFLAG(1)= 0
          if (PRECOND.eq. 2) IFLAG(1)= 0
          if (PRECOND.eq. 3) IFLAG(1)= 0
        endif

        if (IFLAG(1).ne.0) then
          ERROR= HECMW_SOLVER_ERROR_INCONS_PC
          call hecmw_solve_error (hecMESH, ERROR)
        endif
!C===

!C
!C +-----------+
!C | BLOCK LUs |
!C +-----------+
!C===
        if (.not.associated(hecMAT%ALU)) NSET=  0
        if (     associated(hecMAT%ALU)) NSET= -1

        if (NSET.eq.0) then
          allocate (hecMAT%ALU(hecMAT%N))
        endif

        if (NSET.le.0) then
          hecMAT%ALU  = 0.d0

          if (PRECOND.ne.1) then
            do i= 1, hecMAT%N
              hecMAT%ALU(i)= 1.d0/hecMAT%D(i)
            enddo
          endif

          if (PRECOND.eq.1) then
            do i= 1, hecMAT%N
              isL= hecMAT%indexL(i-1) + 1
              ieL= hecMAT%indexL(i)

              W= hecMAT%D(i) * SIGMA_DIAG
              do k= isL, ieL
                SS= 0.d0
                id= hecMAT%itemL(k)

                isU= hecMAT%indexU(id-1) + 1
                ieU= hecMAT%indexU(id)

                do kk= isU, ieU
                  SS= SS + hecMAT%AU(kk) * SIGMA
                enddo
                W= W - hecMAT%AL(k) * SS * hecMAT%ALU(id)
              enddo
              hecMAT%ALU(i)= 1.d0 / W
            enddo
          endif
        endif

!C===

        ! imposing MPC by penalty
        call hecmw_mat_ass_equation ( hecMESH, hecMAT )

        call hecmw_mat_dump(hecMAT, hecMESH)

!C
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===
        call hecmw_solve_CG_11                                          &
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

      !hecMAT%ITERactual = ITER
      !hecMAT%RESIDactual= RESID

      if (ERROR.ne.0) then
        call hecmw_solve_error (hecMESH, ERROR)
      endif

      call hecmw_mat_dump_solution(hecMAT)

      if (hecMESH%my_rank.eq.0 .and. TIMElog.ge.1) then
        TR= (TIME_sol-TIME_comm)/(TIME_sol+1.d-24)*100.d0
        write (*,'(/a)')          '### summary of linear solver'
        write (*,'(i10,a, 1pe16.6)')      ITER, ' iterations  ', RESID
        write (*,'(a, 1pe16.6 )') '    set-up time      : ', TIME_setup
        write (*,'(a, 1pe16.6 )') '    solver time      : ', TIME_sol
        write (*,'(a, 1pe16.6 )') '    solver/comm time : ', TIME_comm
        write (*,'(a, 1pe16.6/)') '    work ratio (%)   : ', TR
      endif

      end subroutine hecmw_solve_11
      end module hecmw_solver_11
