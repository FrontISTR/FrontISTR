!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_iterative
contains
  !
  !C***
  !C*** hecmw_solve_nn
  !C***
  !
  subroutine hecmw_solve_iterative (hecMESH, hecMAT)

    use hecmw_util
    use hecmw_solver_CG
    use hecmw_solver_BiCGSTAB
    use hecmw_solver_GMRES
    use hecmw_solver_GPBiCG
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_solver_las
    use hecmw_precond
    use hecmw_matrix_misc
    use hecmw_matrix_dump
    use hecmw_solver_misc

    implicit none

    type (hecmwST_matrix), target :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH

    integer(kind=kint) :: error
    integer(kind=kint) :: ITER, METHOD, PRECOND, NSET, METHOD2
    integer(kind=kint) :: iterPREmax, i, j
    integer(kind=kint) :: ITERlog, TIMElog
    real(kind=kreal)   :: RESID, SIGMA_DIAG, THRESH, FILTER, resid2
    real(kind=kreal)   :: TIME_setup, TIME_comm, TIME_sol, TR
    real(kind=kreal)   :: time_Ax, time_precond

    integer(kind=kint) :: NREST
    real(kind=kreal)   :: SIGMA

    integer(kind=kint) :: totalmpc, MPC_METHOD
    integer(kind=kint) :: auto_sigma_diag

    !C PARAMETERs
    ITER      = hecmw_mat_get_iter(hecMAT)
    METHOD    = hecmw_mat_get_method(hecMAT)
    METHOD2   = hecmw_mat_get_method2(hecMAT)
    PRECOND   = hecmw_mat_get_precond(hecMAT)
    NSET      = hecmw_mat_get_nset(hecMAT)
    iterPREmax= hecmw_mat_get_iterpremax(hecMAT)
    NREST     = hecmw_mat_get_nrest(hecMAT)
    ITERlog   = hecmw_mat_get_iterlog(hecMAT)
    TIMElog   = hecmw_mat_get_timelog(hecMAT)
    TIME_setup= 0.d0
    TIME_comm = 0.d0
    TIME_sol  = 0.d0
    RESID     = hecmw_mat_get_resid(hecMAT)
    SIGMA_DIAG= hecmw_mat_get_sigma_diag(hecMAT)
    SIGMA     = hecmw_mat_get_sigma(hecMAT)
    THRESH    = hecmw_mat_get_thresh(hecMAT)
    FILTER    = hecmw_mat_get_filter(hecMAT)
    if (SIGMA_DIAG.lt.0.d0) then
      auto_sigma_diag= 1
      SIGMA_DIAG= 1.d0
    else
      auto_sigma_diag= 0
    endif

    !C ERROR CHECK
    call hecmw_solve_check_zerodiag(hecMESH, hecMAT) !C-- ZERO DIAGONAL component

    !C-- IN CASE OF MPC-CG
    totalmpc = hecMESH%mpc%n_mpc
    call hecmw_allreduce_I1 (hecMESH, totalmpc, hecmw_sum)
    MPC_METHOD = hecmw_mat_get_mpc_method(hecMAT)
    if (totalmpc > 0 .and. MPC_METHOD == 2) then
      call hecmw_mat_set_flag_mpcmatvec(hecMAT, 1)
    endif

    !C-- RECYCLE SETTING OF PRECONDITIONER
    call hecmw_mat_recycle_precond_setting(hecMAT)

    ! exchange diagonal elements of overlap region
    call hecmw_mat_dump(hecMAT, hecMESH)
    call hecmw_matvec_set_async(hecMAT)

    !C ITERATIVE solver
    error=0
    !! Auto Sigma_diag loop
    do
      call hecmw_mat_set_flag_converged(hecMAT, 0)
      call hecmw_mat_set_flag_diverged(hecMAT, 0)
      if (auto_sigma_diag.eq.1) call hecmw_mat_set_sigma_diag(hecMAT, SIGMA_DIAG)

      call hecmw_matvec_clear_timer()
      call hecmw_precond_clear_timer()
      call hecmw_solve_iterative_printmsg(hecMESH,hecMAT,METHOD)

      select case(METHOD)
        case (1)  !--CG
          hecMAT%symmetric = .true.
          call hecmw_solve_CG( hecMESH, hecMAT, ITER, RESID, error, TIME_setup, TIME_sol, TIME_comm )
        case (2)  !--BiCGSTAB
          hecMAT%symmetric = .false.
          call hecmw_solve_BiCGSTAB( hecMESH,hecMAT, ITER, RESID, error,TIME_setup, TIME_sol, TIME_comm )
        case (3)  !--GMRES
          hecMAT%symmetric = .false.
          call hecmw_solve_GMRES( hecMESH,hecMAT, ITER, RESID, error, TIME_setup, TIME_sol, TIME_comm )
        case (4)  !--GPBiCG
          hecMAT%symmetric = .false.
          call hecmw_solve_GPBiCG( hecMESH,hecMAT, ITER, RESID, error, TIME_setup, TIME_sol, TIME_comm )
        case default
          error = HECMW_SOLVER_ERROR_INCONS_PC  !!未定義なMETHOD!!
          call hecmw_solve_error (hecMESH, error)
      end select

      if (error==HECMW_SOLVER_ERROR_DIVERGE_PC .or. error==HECMW_SOLVER_ERROR_DIVERGE_MAT) then
        call hecmw_mat_set_flag_diverged(hecMAT, 1)
        if ((PRECOND>=10 .and. PRECOND<20) .and. auto_sigma_diag==1 .and. SIGMA_DIAG<2.d0) then
          SIGMA_DIAG = SIGMA_DIAG + 0.1
          if (hecMESH%my_rank.eq.0) write(*,*) 'Increasing SIGMA_DIAG to', SIGMA_DIAG
          cycle
        elseif (METHOD==1 .and. METHOD2>1) then
          if (auto_sigma_diag.eq.1) SIGMA_DIAG = 1.0
          METHOD = METHOD2
          cycle
        endif
      endif

      if (auto_sigma_diag.eq.1) call hecmw_mat_set_sigma_diag(hecMAT, -1.d0)
      exit
    enddo

    if (error.ne.0) then
      call hecmw_solve_error (hecMESH, error)
    endif

    resid2=hecmw_rel_resid_L2(hecMESH,hecMAT)
    if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
      write(*,"(a,1pe12.5)")'### Relative residual =', resid2
    endif
    if (resid2 < hecmw_mat_get_resid(hecMAT)) call hecmw_mat_set_flag_converged(hecMAT, 1)

    call hecmw_mat_dump_solution(hecMAT)
    call hecmw_matvec_unset_async

    !C-- IN CASE OF MPC-CG
    if (totalmpc > 0 .and. MPC_METHOD == 2) then
      call hecmw_mat_set_flag_mpcmatvec(hecMAT, 0)
    endif

    time_Ax = hecmw_matvec_get_timer()
    time_precond = hecmw_precond_get_timer()

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

  end subroutine hecmw_solve_iterative

  subroutine hecmw_solve_check_zerodiag (hecMESH, hecMAT)
    use hecmw_util
    use hecmw_matrix_misc
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix), target :: hecMAT
    integer (kind=kint) :: IFLAG
    integer (kind=kint)::PRECOND,iterPREmax,i,j,error
    PRECOND   = hecmw_mat_get_precond(hecMAT)
    iterPREmax= hecmw_mat_get_iterpremax(hecMAT)
    !C
    !C-- ZERO DIAGONAL component
    error= 0
    do i= 1, hecMAT%N
      do j = 1, hecMAT%NDOF
        if (dabs(hecMAT%D(hecMAT%NDOF*hecMAT%NDOF*(i-1)+(j-1)*(hecMAT%NDOF+1)+1)).eq.0.d0) then
          error=HECMW_SOLVER_ERROR_ZERO_DIAG
        end if
      end do
    enddo

    call hecmw_allreduce_I1 (hecMESH, error, hecmw_max)
    if (error.ne.0 .and. (PRECOND.lt.10 .and. iterPREmax.gt.0)) then
      call hecmw_solve_error (hecMESH, error)
    endif

  end subroutine hecmw_solve_check_zerodiag

  function hecmw_solve_check_zerorhs (hecMESH, hecMAT)
    use hecmw_util
    use hecmw_matrix_misc
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix), target :: hecMAT
    real(kind=kreal),    dimension(1) :: RHS
    integer (kind=kint)::PRECOND,iterPREmax,i,j,error
    logical :: hecmw_solve_check_zerorhs

    PRECOND   = hecmw_mat_get_precond(hecMAT)
    iterPREmax= hecmw_mat_get_iterpremax(hecMAT)
    !C
    !C-- ZERO RHS norm
    error= 0
    hecmw_solve_check_zerorhs = .false.

    RHS(1)= 0.d0
    do i= 1, hecMAT%N
      do j = 1, hecMAT%NDOF
        RHS(1)=RHS(1) + hecMAT%B(hecMAT%NDOF*(i-1)+j)**2
      end do
    enddo
    if (hecMESH%mpc%n_mpc > 0) then
      do i= 1, hecMESH%mpc%n_mpc
        RHS(1)= RHS(1) + hecMESH%mpc%mpc_const(i)**2
      enddo
    endif
    call hecmw_allreduce_R (hecMESH, RHS, 1, hecmw_sum)

    if (RHS(1).eq.0.d0) then
      error= HECMW_SOLVER_ERROR_ZERO_RHS
      call hecmw_solve_error (hecMESH, error)
      hecMAT%X(:)=0.d0
      hecmw_solve_check_zerorhs = .true.
    endif

  end function hecmw_solve_check_zerorhs

  subroutine hecmw_solve_iterative_printmsg (hecMESH, hecMAT, METHOD)
    use hecmw_util
    use hecmw_solver_misc
    use hecmw_matrix_misc

    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix), target :: hecMAT
    integer(kind=kint) :: METHOD
    integer(kind=kint) :: ITER, PRECOND, NSET, iterPREmax, NREST
    integer(kind=kint) :: ITERlog, TIMElog

    character(len=30) :: msg_precond
    character(len=30) :: msg_method

    ITER      = hecmw_mat_get_iter(hecMAT)
    ! METHOD    = hecmw_mat_get_method(hecMAT)
    PRECOND   = hecmw_mat_get_precond(hecMAT)
    NSET      = hecmw_mat_get_nset(hecMAT)
    iterPREmax= hecmw_mat_get_iterpremax(hecMAT)
    NREST     = hecmw_mat_get_nrest(hecMAT)
    ITERlog= hecmw_mat_get_iterlog(hecMAT)
    TIMElog= hecmw_mat_get_timelog(hecMAT)

    select case(METHOD)
      case (1)  !--CG
        msg_method="CG"
      case (2)  !--BiCGSTAB
        msg_method="BiCGSTAB"
      case (3)  !--GMRES
        msg_method="GMRES"
      case (4)  !--GPBiCG
        msg_method="GPBiCG"
      case default
        msg_method="Unlabeled"
    end select
    select case(PRECOND)
      case (1,2)
        msg_precond="SSOR"
      case (3)
        msg_precond="DIAG"
      case (5)
        msg_precond="ML"
      case (7)
        msg_precond="DirectMUMPS"
      case (10, 11, 12)
        write(msg_precond,"(a,i0,a)") "ILU(",PRECOND-10,")"
      case (20)
        msg_precond="SAINV"
      case (21)
        msg_precond="RIF"
      case default
        msg_precond="Unlabeled"
    end select
    if (hecMESH%my_rank.eq.0 .and. (ITERlog.eq.1 .or. TIMElog.ge.1)) then
      write (*,'(a,i0,a,i0,a,a,a,a,a,i0)') '### ',hecMAT%NDOF,'x',hecMAT%NDOF,' BLOCK ', &
        &   trim(msg_method),", ",trim(msg_precond),", ", iterPREmax
    end if
  end subroutine hecmw_solve_iterative_printmsg

end module hecmw_solver_iterative
