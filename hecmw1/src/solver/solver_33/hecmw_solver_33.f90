!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_33
contains
  !
  !C***
  !C*** hecmw_solve_33
  !C***
  !
  subroutine hecmw_solve_33 (hecMESH, hecMAT)

    use hecmw_util
    use hecmw_solver_CG_33
    use hecmw_solver_BiCGSTAB_33
    use hecmw_solver_GMRES_33
    use hecmw_solver_GPBiCG_33
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_solver_las_33
    use hecmw_precond_33
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
    real(kind=kreal) :: time_Ax, time_precond
    real(kind=kreal),    dimension(1) :: RHS
    integer (kind=kint), dimension(1) :: IFLAG

    integer(kind=kint) :: NREST
    real(kind=kreal)   :: SIGMA

    integer(kind=kint) :: totalmpc, MPC_METHOD
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
      RHS(1)= RHS(1) + hecMAT%B(3*i-2)**2 + hecMAT%B(3*i-1)**2       &
           &                   + hecMAT%B(3*i  )**2
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
    IFLAG(1)= 0
    do i= 1, hecMAT%N
      if (dabs(hecMAT%D(9*i-8)).eq.0.d0) IFLAG(1)= 1
      if (dabs(hecMAT%D(9*i-4)).eq.0.d0) IFLAG(1)= 1
      if (dabs(hecMAT%D(9*i  )).eq.0.d0) IFLAG(1)= 1
    enddo

    call hecmw_allreduce_I (hecMESH, IFLAG, 1, hecmw_sum)
    if (IFLAG(1).ne.0 .and. (PRECOND.lt.10 .and. iterPREmax.gt.0)) then
      ERROR= HECMW_SOLVER_ERROR_ZERO_DIAG
      call hecmw_solve_error (hecMESH, ERROR)
    endif

    !C
    !C-- INCONSISTENT SOLVER/PRECONDITIONING
    IFLAG(1)= 0
    if (METHOD.le.0 .or. METHOD.ge.5)    IFLAG(1)= 1
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
    if (PRECOND.eq.10) IFLAG(1)= 0
    if (PRECOND.eq.11) IFLAG(1)= 0
    if (PRECOND.eq.12) IFLAG(1)= 0
    if (PRECOND.eq.20) IFLAG(1)= 0
    if (PRECOND.eq.21) IFLAG(1)= 0

    if (IFLAG(1).ne.0) then
      ERROR= HECMW_SOLVER_ERROR_INCONS_PC
      call hecmw_solve_error (hecMESH, ERROR)
    endif

    !C
    !C-- IN CASE OF MPC-CG
    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)
    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)
    if (totalmpc > 0 .and. MPC_METHOD == 2) then
      call hecmw_matvec_33_set_mpcmatvec_flg (.true.)
    endif

    !C
    !C-- RECYCLE SETTING OF PRECONDITIONER
    call hecmw_mat_recycle_precond_setting(hecMAT)

    ! exchange diagonal elements of overlap region
    !call hecmw_mat_diag_sr_33(hecMESH, hecMAT)

    call hecmw_mat_dump(hecMAT, hecMESH)

    call hecmw_matvec_33_set_async(hecMAT)

    !C===
    !C +------------------+
    !C | ITERATIVE solver |
    !C +------------------+
    !C===

    do

    if (auto_sigma_diag.eq.1) call hecmw_mat_set_sigma_diag(hecMAT, SIGMA_DIAG)

    call hecmw_matvec_33_clear_timer()
    call hecmw_precond_33_clear_timer()

    !C
    !C-- CG
    if (METHOD.eq.1) then
      if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
        if (PRECOND.le.2) &
             write (*,'(a,i3)') '### 3x3 B-SSOR-CG(0)',iterPREmax
        if (PRECOND.eq.3) &
             write (*,'(a,i3)') '### 3x3 B-scale-CG  ',iterPREmax
        if (PRECOND.eq.10) &
             write (*,'(a,i3)') '### 3x3 B-IC-CG(0)',iterPREmax
        if (PRECOND.eq.11) &
             write (*,'(a,i3)') '### 3x3 B-IC-CG(1)',iterPREmax
        if (PRECOND.eq.12) &
             write (*,'(a,i3)') '### 3x3 B-IC-CG(2)',iterPREmax
        if (PRECOND.eq.5) &
             write (*,'(a,i3)') '### 3x3 ML-CG  ',iterPREmax
      endif
      call hecmw_solve_CG_33( hecMESH,  hecMAT, ITER, RESID, ERROR,   &
           &                          TIME_setup, TIME_sol, TIME_comm )
    endif

    !C
    !C-- BiCGSTAB
    if (METHOD.eq.2) then
      if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
        if (PRECOND.le.2) &
             write (*,'(a,i3)') '### 3x3 B-SSOR-BiCGSTAB(0)',iterPREmax
        if (PRECOND.eq.3) &
             write (*,'(a,i3)') '### 3x3 B-scale-BiCGSTAB  ',iterPREmax
        if (PRECOND.eq.10) &
             write (*,'(a,i3)') '### 3x3 B-IlU-BiCGSTAB(0)',iterPREmax
        if (PRECOND.eq.11) &
             write (*,'(a,i3)') '### 3x3 B-IlU-BiCGSTAB(1)',iterPREmax
        if (PRECOND.eq.12) &
             write (*,'(a,i3)') '### 3x3 B-IlU-BiCGSTAB(2)',iterPREmax
        if (PRECOND.eq.5) &
             write (*,'(a,i3)') '### 3x3 ML-BiCGSTAB  ',iterPREmax
      endif

      call hecmw_solve_BiCGSTAB_33( hecMESH,  hecMAT, ITER, RESID, ERROR, &
           &                                TIME_setup, TIME_sol, TIME_comm )
    endif

    !C
    !C-- GMRES
    if (METHOD.eq.3) then
      if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
        if (PRECOND.le.2) &
             write (*,'(a,i3)') '### 3x3 B-SSOR-GMRES(0)',iterPREmax
        if (PRECOND.eq.3) &
             write (*,'(a,i3)') '### 3x3 B-scale-GMRES  ',iterPREmax
        if (PRECOND.eq.10) &
             write (*,'(a,i3)') '### 3x3 B-ILU-GMRES(0)',iterPREmax
        if (PRECOND.eq.11) &
             write (*,'(a,i3)') '### 3x3 B-ILU-GMRES(1)',iterPREmax
        if (PRECOND.eq.12) &
             write (*,'(a,i3)') '### 3x3 B-ILU-GMRES(2)',iterPREmax
        if (PRECOND.eq.5) &
             write (*,'(a,i3)') '### 3x3 ML-GMRES  ',iterPREmax
      endif

      call hecmw_solve_GMRES_33( hecMESH,  hecMAT, ITER, RESID, ERROR, &
           &                             TIME_setup, TIME_sol, TIME_comm )
    endif

    !C
    !C-- GPBiCG
    if (METHOD.eq.4) then
      if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
        if (PRECOND.le.2) &
             write (*,'(a,i3)') '### 3x3 B-SSOR-GPBiCG(0)',iterPREmax
        if (PRECOND.eq.3) &
             write (*,'(a,i3)') '### 3x3 B-scale-GPBiCG  ',iterPREmax
        if (PRECOND.eq.10) &
             write (*,'(a,i3)') '### 3x3 B-ILU-GPBiCG(0)',iterPREmax
        if (PRECOND.eq.11) &
             write (*,'(a,i3)') '### 3x3 B-ILU-GPBiCG(1)',iterPREmax
        if (PRECOND.eq.12) &
             write (*,'(a,i3)') '### 3x3 B-ILU-GPBiCG(2)',iterPREmax
        if (PRECOND.eq.5) &
             write (*,'(a,i3)') '### 3x3 ML-GPBiCG  ',iterPREmax
      endif

      call hecmw_solve_GPBiCG_33( hecMESH,  hecMAT, ITER, RESID, ERROR, &
           &                              TIME_setup, TIME_sol, TIME_comm )
    endif

    if ((ERROR.eq.HECMW_SOLVER_ERROR_DIVERGE_PC .or. &
         ERROR.eq.HECMW_SOLVER_ERROR_DIVERGE_MAT) .and. &
         (PRECOND.ge.10 .and. PRECOND.lt.20) .and. &
         auto_sigma_diag.eq.1 .and. &
         SIGMA_DIAG.lt.2.d0) then
      SIGMA_DIAG = SIGMA_DIAG + 0.1
      if (hecMESH%my_rank.eq.0) write(*,*) 'Increasing SIGMA_DIAG to', SIGMA_DIAG
    else
      if (auto_sigma_diag.eq.1) call hecmw_mat_set_sigma_diag(hecMAT, -1.d0)
      exit
    endif

    enddo

    if (ERROR.ne.0) then
      call hecmw_solve_error (hecMESH, ERROR)
    endif

    resid2=hecmw_rel_resid_L2_33(hecMESH,hecMAT)
    if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
      write(*,"(a,1pe12.5)")'### Relative residual =', resid2
      ! if( resid2 >= 1.0d-8) then
      !   write(*,"(a)")'### Relative residual exceeded 1.0d-8---Iterative Solver### '
      ! endif
    endif

    call hecmw_mat_dump_solution(hecMAT)

    call hecmw_matvec_33_unset_async

    !C
    !C-- IN CASE OF MPC-CG
    if (totalmpc > 0 .and. MPC_METHOD == 2) then
      call hecmw_matvec_33_set_mpcmatvec_flg (.false.)
    endif

    time_Ax = hecmw_matvec_33_get_timer()
    time_precond = hecmw_precond_33_get_timer()

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
      write (*,'(a, 1pe16.6/)') '    work ratio (%)   : ', TR
    endif

  end subroutine hecmw_solve_33
end module hecmw_solver_33
