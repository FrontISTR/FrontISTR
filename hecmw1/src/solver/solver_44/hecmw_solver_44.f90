!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_44
contains
!
!C***
!C*** hecmw_solve_44
!C***
!
  subroutine hecmw_solve_44 (hecMESH, hecMAT)

    use hecmw_util
    use hecmw_solver_CG_44
    use hecmw_solver_BiCGSTAB_44
!        use hecmw_solver_GMRES_44
!        use hecmw_solver_GPBiCG_44
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_matrix_ass
    use hecmw_matrix_contact
    use hecmw_solver_las_44
    use hecmw_precond_44
    use hecmw_local_matrix
    use hecmw_matrix_misc
    use hecmw_matrix_dump
    use hecmw_solver_misc

    implicit none

    type (hecmwST_matrix), target :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH

    integer(kind=kint) :: ERROR

    integer(kind=kint) :: ITER, METHOD, PRECOND, NSET
    integer(kind=kint) :: iterPREmax, i

    real(kind=kreal) :: RESID, SIGMA_DIAG, THRESH, FILTER, resid2

    integer(kind=kint) :: ITERlog, TIMElog
    real(kind=kreal) :: TIME_setup, TIME_comm, TIME_sol, TR
    real(kind=kreal) :: time_Ax, time_precond, time_dumm
    real(kind=kreal) :: S_TIME, E_TIME, TIME_mpc_pre, TIME_mpc_post
    real(kind=kreal),    dimension(1) :: RHS
    integer (kind=kint), dimension(1) :: IFLAG

    integer(kind=kint) :: NREST
    real(kind=kreal)   :: SIGMA

    type (hecmwST_matrix), pointer, save :: hecTKT
    logical, save :: first_call = .true.
    integer(kind=kint) :: totalmpc, MPC_METHOD
    real(kind=kreal), allocatable :: Btmp(:)
    real(kind=kreal)::t_max,t_min,t_avg,t_sd

    integer(kind=kint) :: auto_sigma_diag

!C===
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
    ITER      = hecmw_mat_get_iter(hecMAT)
    METHOD    = hecmw_mat_get_method(hecMAT)
    PRECOND   = hecmw_mat_get_precond(hecMAT)
    NSET      = hecmw_mat_get_nset(hecMAT)
    iterPREmax= hecmw_mat_get_iterpremax(hecMAT)
    NREST     = hecmw_mat_get_nrest(hecMAT)

    ITERlog= hecmw_mat_get_iterlog(hecMAT)
    TIMElog= hecmw_mat_get_timelog(hecMAT)

    TIME_setup  = 0.d0
    TIME_comm   = 0.d0
    TIME_sol    = 0.d0

    RESID     = hecmw_mat_get_resid(hecMAT)
    SIGMA_DIAG= hecmw_mat_get_sigma_diag(hecMAT)
    SIGMA     = hecmw_mat_get_sigma(hecMAT)

    THRESH= hecmw_mat_get_thresh(hecMAT)
    FILTER= hecmw_mat_get_filter(hecMAT)

    if (SIGMA_DIAG.lt.0.d0) then
      auto_sigma_diag= 1
      SIGMA_DIAG= 1.d0
    else
      auto_sigma_diag= 0
    endif

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
      RHS(1)= RHS(1) + hecMAT%B(4*i-3)**2 + hecMAT%B(4*i-2)**2 &
          &          + hecMAT%B(4*i-1)**2 + hecMAT%B(4*i  )**2
    enddo
    if (hecMESH%mpc%n_mpc > 0) then
      do i= 1, hecMESH%mpc%n_mpc
        RHS(1)= RHS(1) + hecMESH%mpc%mpc_const(i)**2
      enddo
    endif
    call hecmw_allreduce_R (hecMESH, RHS, 1, hecmw_sum)

    if (RHS(1).eq.0.d0) then
      ERROR= HECMW_SOLVER_ERROR_ZERO_RHS
      call hecmw_solve_error (hecMESH, ERROR)
      hecMAT%X(:)=0.d0
      return
    endif

!C
!C-- ZERO DIAGONAL component
    IFLAG(1)= 0
    do i= 1, hecMAT%N
      if (dabs(hecMAT%D(16*i-15)).eq.0.d0) IFLAG(1)= 1
      if (dabs(hecMAT%D(16*i-10)).eq.0.d0) IFLAG(1)= 1
      if (dabs(hecMAT%D(16*i- 5)).eq.0.d0) IFLAG(1)= 1
      if (dabs(hecMAT%D(16*i   )).eq.0.d0) IFLAG(1)= 1
    enddo

    call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
    if (IFLAG(1).ne.0 .and. (PRECOND.le.10 .and. iterPREmax.gt.0)) then
      ERROR= HECMW_SOLVER_ERROR_ZERO_DIAG
      call hecmw_solve_error (hecMESH, ERROR)
    endif

!C
!C-- INCONSISTENT SOLVER/PRECONDITIONING
    IFLAG(1)= 0
    if (METHOD.le.0 .or. METHOD.ge.7)    IFLAG(1)= 1
    if (PRECOND.le.0 .or. PRECOND.gt.21) IFLAG(1)= 1

    call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
    if (IFLAG(1).ne.0) then
      ERROR= HECMW_SOLVER_ERROR_INCONS_PC
      call hecmw_solve_error (hecMESH, ERROR)
    endif

    IFLAG(1)= 1
    if (PRECOND.eq. 1) IFLAG(1)= 0
    if (PRECOND.eq. 2) IFLAG(1)= 0
    if (PRECOND.eq. 3) IFLAG(1)= 0
    if (PRECOND.eq. 5) IFLAG(1)= 0
    if (PRECOND.eq. 6) IFLAG(1)= 0
    if (PRECOND.eq.10) IFLAG(1)= 0
    if (PRECOND.eq.11) IFLAG(1)= 0
    if (PRECOND.eq.12) IFLAG(1)= 0
    if (PRECOND.eq.20) IFLAG(1)= 0
    if (PRECOND.eq.21) IFLAG(1)= 0

    if (IFLAG(1).ne.0) then
      ERROR= HECMW_SOLVER_ERROR_INCONS_PC
      call hecmw_solve_error (hecMESH, ERROR)
    endif

    !C===
    !C +-------------+
    !C | MPC Preproc |
    !C +-------------+
    !C===
    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)

    S_TIME= HECMW_WTIME()

    if (totalmpc > 0) then
      call hecmw_mpc_scale(hecMESH)

      MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)
      if (MPC_METHOD < 1 .or. 3 < MPC_METHOD) MPC_METHOD = 3

      if (MPC_METHOD == 1) then  ! penalty
        !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: Penalty"
        call hecmw_mat_ass_equation ( hecMESH, hecMAT )
        hecTKT => hecMAT
      else if (MPC_METHOD == 2) then  ! MPCCG
        !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: MPC-CG"
        call hecmw_matvec_44_set_mpcmatvec_flg (.true.)
        allocate(Btmp(hecMAT%NP * hecMAT%NDOF))
        do i=1,hecMAT%NP * hecMAT%NDOF
          Btmp(i) = hecMAT%B(i)
        enddo
        call hecmw_trans_b_44(hecMESH, hecMAT, Btmp, hecMAT%B, time_dumm)
        hecTKT => hecMAT
      else if (MPC_METHOD == 3) then  ! elimination
        !if (hecMESH%my_rank.eq.0) write(0,*) "MPC Method: Elimination"
        if (hecMAT%Iarray(97)>=1) then
          if (first_call) then
            allocate(hecTKT)
            first_call = .false.
          else
            call hecmw_mat_finalize(hecTKT)
          endif
          call hecmw_mat_init(hecTKT)
          call hecmw_trimatmul_TtKT_mpc(hecMESH, hecMAT, hecTKT)
        endif
        call hecmw_trans_b_44(hecMESH, hecMAT, hecMAT%B, hecTKT%B, time_dumm)
      endif
    else
      hecTKT => hecMAT
    endif

    !C
    !C-- RECYCLE SETTING OF PRECONDITIONER
    call hecmw_mat_recycle_precond_setting(hecMAT)

    E_TIME= HECMW_WTIME()
    if (TIMElog.eq.2) then
      call hecmw_time_statistics(hecMESH, E_TIME - S_TIME, &
           t_max, t_min, t_avg, t_sd)
      if (hecMESH%my_rank.eq.0) then
        write(*,*) 'Time MPC pre'
        write(*,*) '  Max     :',t_max
        write(*,*) '  Min     :',t_min
        write(*,*) '  Avg     :',t_avg
        write(*,*) '  Std Dev :',t_sd
      endif
      TIME_mpc_pre = t_max
    else
      TIME_mpc_pre = E_TIME - S_TIME
    endif

    ! exchange diagonal elements of overlap region
    !call hecmw_mat_diag_sr_44(hecMESH, hecTKT)

    call hecmw_mat_dump(hecTKT, hecMESH)

    !call hecmw_matvec_44_set_async(hecTKT)

!C===
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===

    call hecmw_matvec_44_clear_timer()
    call hecmw_precond_44_clear_timer()
!C
!C-- CG
      if (METHOD.eq.1) then
        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
          if (PRECOND.eq.1) &
               write (*,'(a,i3)') '### 4x4 B-IC-CG(0)',iterPREmax
          if (PRECOND.eq.2) &
               write (*,'(a,i3)') '### 4x4 B-SSOR-CG(0)',iterPREmax
          if (PRECOND.eq.3) &
               write (*,'(a,i3)') '### 4x4 B-scale-CG  ',iterPREmax
          if (PRECOND.eq.10) &
               write (*,'(a,i3)') '### 4x4 B-IC-CG(0)',iterPREmax
          if (PRECOND.eq.11) &
               write (*,'(a,i3)') '### 4x4 B-IC-CG(1)',iterPREmax
          if (PRECOND.eq.12) &
               write (*,'(a,i3)') '### 4x4 B-IC-CG(2)',iterPREmax
        endif
        call hecmw_solve_CG_44( hecMESH,  hecTKT, ITER, RESID, ERROR,   &
     &                          TIME_setup, TIME_sol, TIME_comm )
      endif

!C
!C-- BiCGSTAB
      if (METHOD.eq.2) then
        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
          if (PRECOND.eq.1) &
               write (*,'(a,i3)') '### 4x4 B-ILU-BiCGSTAB(0)',iterPREmax
          if (PRECOND.eq.2) &
               write (*,'(a,i3)') '### 4x4 B-SSOR-BiCGSTAB(0)',iterPREmax
          if (PRECOND.eq.3) &
               write (*,'(a,i3)') '### 4x4 B-scale-BiCGSTAB  ',iterPREmax
          if (PRECOND.eq.10) &
               write (*,'(a,i3)') '### 4x4 B-IlU-BiCGSTAB(0)',iterPREmax
          if (PRECOND.eq.11) &
               write (*,'(a,i3)') '### 4x4 B-IlU-BiCGSTAB(1)',iterPREmax
          if (PRECOND.eq.12) &
               write (*,'(a,i3)') '### 4x4 B-IlU-BiCGSTAB(2)',iterPREmax
        endif

        call hecmw_solve_BiCGSTAB_44( hecMESH,  hecTKT, ITER, RESID, ERROR, &
     &                                TIME_setup, TIME_sol, TIME_comm )
      endif

!C
!C-- GMRES
!      if (METHOD.eq.3) then
!        ! imposing MPC by penalty
!        call hecmw_mat_ass_equation ( hecMESH, hecTKT )
!
!        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
!          if (PRECOND.eq.1) &
!               write (*,'(a,i3)') '### 4x4 B-ILU-GMRES(0)',iterPREmax
!          if (PRECOND.eq.2) &
!               write (*,'(a,i3)') '### 4x4 B-SSOR-GMRES(0)',iterPREmax
!          if (PRECOND.eq.3) &
!               write (*,'(a,i3)') '### 4x4 B-scale-GMRES  ',iterPREmax
!          if (PRECOND.eq.10) &
!               write (*,'(a,i3)') '### 4x4 B-ILU-GMRES(0)',iterPREmax
!          if (PRECOND.eq.11) &
!               write (*,'(a,i3)') '### 4x4 B-ILU-GMRES(1)',iterPREmax
!          if (PRECOND.eq.12) &
!               write (*,'(a,i3)') '### 4x4 B-ILU-GMRES(2)',iterPREmax
!        endif
!
!        call hecmw_solve_GMRES_44( hecMESH,  hecTKT, ITER, RESID, ERROR, &
!     &                             TIME_setup, TIME_sol, TIME_comm )
!      endif

!C
!C-- GPBiCG
!      if (METHOD.eq.4) then
!        ! imposing MPC by penalty
!        call hecmw_mat_ass_equation ( hecMESH, hecTKT )
!
!        if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
!          if (PRECOND.eq.1) &
!               write (*,'(a,i3)') '### 4x4 B-ILU-GPBiCG(0)',iterPREmax
!          if (PRECOND.eq.2) &
!               write (*,'(a,i3)') '### 4x4 B-SSOR-GPBiCG(0)',iterPREmax
!          if (PRECOND.eq.3) &
!               write (*,'(a,i3)') '### 4x4 B-scale-GPBiCG  ',iterPREmax
!          if (PRECOND.eq.10) &
!               write (*,'(a)') '### 4x4 B-ILU-GPBiCG(0)',iterPREmax
!          if (PRECOND.eq.11) &
!               write (*,'(a)') '### 4x4 B-ILU-GPBiCG(1)',iterPREmax
!          if (PRECOND.eq.12) &
!               write (*,'(a)') '### 4x4 B-ILU-GPBiCG(2)',iterPREmax
!        endif
!
!        call hecmw_solve_GPBiCG_44( hecMESH,  hecTKT, ITER, RESID, ERROR, &
!     &                              TIME_setup, TIME_sol, TIME_comm )
!      endif

    if (ERROR.ne.0) then
      call hecmw_solve_error (hecMESH, ERROR)
    endif

    !hecTKT => hecMAT%hecTKT
    resid2=hecmw_rel_resid_L2_44(hecMESH, hecTKT)
    if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
      write(*,"(a,1pe12.5)")'### Relative residual =', resid2
      ! if( resid2 >= 1.0d-8) then
      !   write(*,"(a)")'### Relative residual exceeded 1.0d-8---Iterative Solver### '
      ! endif
    endif

    call hecmw_mat_dump_solution(hecTKT)

    !call hecmw_matvec_44_unset_async

    !C===
    !C +--------------+
    !C | MPC Postproc |
    !C +--------------+
    !C===
    call hecmw_barrier(hecMESH)
    S_TIME= HECMW_WTIME()

    if (totalmpc > 0) then
      if (MPC_METHOD == 1) then  ! penalty
        ! do nothing
      else if (MPC_METHOD == 2) then  ! MPCCG
        call hecmw_tback_x_44(hecMESH, hecTKT%X, time_dumm)
        do i=1,hecMAT%NP * hecMAT%NDOF
          hecMAT%B(i) = Btmp(i)
        enddo
        deallocate(Btmp)
      else if (MPC_METHOD == 3) then  ! elimination
        call hecmw_tback_x_44(hecMESH, hecTKT%X, time_dumm)
        do i=1,hecMAT%NP * hecMAT%NDOF
          hecMAT%X(i)=hecTKT%X(i)
        enddo
        hecMAT%Iarray(:)=hecTKT%Iarray(:)
        hecMAT%Rarray(:)=hecTKT%Rarray(:)
        ! call hecmw_mat_finalize(hecTKT)
        ! deallocate(hecTKT)
      endif
    endif

    call hecmw_barrier(hecMESH)
    E_TIME= HECMW_WTIME()
    TIME_mpc_post = E_TIME - S_TIME

    time_Ax = hecmw_matvec_44_get_timer()
    time_precond = hecmw_precond_44_get_timer()

    if (hecMESH%my_rank.eq.0 .and. TIMElog.ge.1) then
      TR= (TIME_sol-TIME_comm)/(TIME_sol+1.d-24)*100.d0
      write (*,'(/a)')          '### summary of linear solver'
      write (*,'(i10,a, 1pe16.6)')      ITER, ' iterations  ', RESID
      write (*,'(a, 1pe16.6 )') '    set-up time      : ', TIME_setup
      write (*,'(a, 1pe16.6 )') '    solver time      : ', TIME_sol
      write (*,'(a, 1pe16.6 )') '    solver/comm time : ', TIME_comm
      write (*,'(a, 1pe16.6 )') '    solver/matvec    : ', time_Ax
      write (*,'(a, 1pe16.6 )') '    solver/precond   : ', time_precond
      if (ITER > 0) &
      write (*,'(a, 1pe16.6 )') '    solver/1 iter    : ', TIME_sol / ITER
      if (totalmpc > 0) then
        write (*,'(a, 1pe16.6 )') '    MPC pre          : ', TIME_mpc_pre
        write (*,'(a, 1pe16.6 )') '    MPC post         : ', TIME_mpc_post
      endif
      write (*,'(a, 1pe16.6/)') '    work ratio (%)   : ', TR
    endif

  end subroutine hecmw_solve_44
end module hecmw_solver_44
