!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C***
!C*** module hecmw_solver_GMRESREN
!C***
!
module hecmw_solver_GMRESREN

  public :: hecmw_solve_GMRESREN

contains
  !C
  !C*** hecmw_solve_GMRESREN
  !C
  subroutine hecmw_solve_GMRESREN( hecMESH,  hecMAT, ITER, RESID, error, &
      &                                    Tset, Tsol, Tcomm )
    use hecmw_util
    use m_hecmw_solve_error
    use m_hecmw_comm_f
    use hecmw_matrix_misc
    use hecmw_solver_misc
    use hecmw_solver_las
    use hecmw_solver_scaling
    use hecmw_precond
    use hecmw_estimate_condition

    implicit none

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint ), intent(inout):: ITER, error
    real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm

    integer(kind=kint ) :: N, NP, NDOF, NNDOF
    integer(kind=kint ) :: my_rank
    integer(kind=kint ) :: ITERlog, TIMElog
    real(kind=kreal), pointer :: B(:), X(:)

    real(kind=kreal), dimension(:)  ,  allocatable :: vecR,workPC
    real(kind=kreal), dimension(:,:),  allocatable :: u,c,uin,cin,sBFGS,yBFGS,xi,eta

    integer(kind=kint ) :: MAXIT, NREST

    real   (kind=kreal) :: TOL

    real   (kind=kreal)   ZERO, ONE
    parameter ( ZERO = 0.0D+0, ONE = 1.0D+0 )

    integer(kind=kint ) :: NRK,i,k,kk,jj,INFO,ik,iOrth
    integer(kind=kint ) :: IROW
    real   (kind=kreal) :: S_TIME,E_TIME,S1_TIME,E1_TIME
    real   (kind=kreal) :: LDH,LDW,BNRM2,DNRM2,RNORM
    real   (kind=kreal) :: COMMtime,COMPtime, coef,val,VCS,VSN,DTEMP,AA,BB,R0,scale,RR
    integer(kind=kint ) :: ESTCOND
    real   (kind=kreal) :: t_max,t_min,t_avg,t_sd
    real   (kind=kreal) :: alpha,beta
   

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

    error= 0

    allocate (vecR(NDOF*NP))
    allocate (workPC(NDOF*NP))
    allocate (u  (NDOF*NP,NREST))
    allocate (c  (NDOF*NP,NREST))
    allocate (uin(NDOF*NP,NREST))
    allocate (cin(NDOF*NP,NREST))
    allocate (xi (NDOF*NP,NREST))
    allocate (eta(NDOF*NP,NREST))

    COMMtime= 0.d0
    COMPtime= 0.d0

    !C
    !C-- SCALING
    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

    !C===
    !C +----------------------+
    !C | SETUP PRECONDITIONER |
    !C +----------------------+
    !C===
    call hecmw_precond_setup(hecMAT, hecMESH, 0)


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

      call hecmw_matresid(hecMESH, hecMAT, X, B, vecR, Tcomm)
      do I = 1, NREST
        ITER= ITER + 1

        !C Compute  xi(1) = (I-AM^-1)r
        !C Compute eta(1) = M^-1r
        call hecmw_precond_apply(hecMESH, hecMAT, vecR, eta(:,1), workPC, Tcomm)
        call hecmw_matvec(hecMESH, hecMAT, eta(:,1), xi(:,1), Tcomm)
        call hecmw_xpay_R(NNDOF, -1.0d0, vecR, xi(:,1))

        do iOrth = 1, I-1
           !C alpha = c_{i}^T xi_{i}
           call hecmw_InnerProduct_R(hecMESH, NDOF, c(:,iOrth), xi(:,iOrth), alpha, Tcomm)

           !C  xi(i+1) =  xi(i) - alpha * c(i)
           !C eta(i+1) = eta(i) + alpha * u(i)
           call hecmw_axpyz_R(NNDOF, -alpha, c(:,iOrth),  xi(:,iOrth),  xi(:,iOrth+1))
           call hecmw_axpyz_R(NNDOF,  alpha, u(:,iOrth), eta(:,iOrth), eta(:,iOrth+1))
        enddo

        !C Solve M*r = uin(:,1)
        call hecmw_precond_apply(hecMESH, hecMAT, xi(:,I), uin(:,1), workPC, Tcomm)
        !C cin(:,1) = A*uin(:,1)
        call hecmw_matvec(hecMESH, hecMAT, uin(:,1), cin(:,1), Tcomm)

        do iOrth = 1, I-1
           !C c_{i}^T cin_{i}
           call hecmw_InnerProduct_R(hecMESH, NDOF, c(:,iOrth), cin(:,iOrth), beta, Tcomm)

           call hecmw_axpyz_R(NNDOF, -coef, c(:,iOrth), cin(:,iOrth), cin(:,iOrth+1))
           call hecmw_axpyz_R(NNDOF, -coef, u(:,iOrth), uin(:,iOrth), uin(:,iOrth+1))
        enddo
        call hecmw_InnerProduct_R(hecMESH, NDOF, cin(:,I), cin(:,I), coef, Tcomm)
        coef = 1.0d0 / dsqrt(coef)
        call hecmw_axpby_R(NNDOF, coef, 0.0d0, cin(:,I), c(:,I))
        call hecmw_axpby_R(NNDOF, coef, 0.0d0, uin(:,I), u(:,I))

        call hecmw_InnerProduct_R(hecMESH, NDOF, c(:,I), xi(:,I), coef, Tcomm)
        call hecmw_axpy_R (NNDOF,  coef,   u(:,I), x)
        call hecmw_axpy_R (NNDOF, 1.0d0, eta(:,I), x)
        call hecmw_axpyz_R(NNDOF, -coef,   c(:,I), xi(:,I), vecR)

        call hecmw_InnerProduct_R(hecMESH, NDOF, vecR, vecR, DNRM2, Tcomm)

        RESID= dsqrt(DNRM2/BNRM2)

        !C##### ITERATION HISTORY
        if (my_rank.eq.0.and.ITERLog.eq.1) write (*,'(i7, 1pe16.6)') ITER, RESID
        !C#####

        if ( RESID.le.TOL )   exit OUTER
        if ( ITER.gt.MAXIT ) then
          error = HECMW_SOLVER_ERROR_NOCONV_MAXIT
          exit OUTER
        end if
      end do

    end do OUTER

    call hecmw_solver_scaling_bk(hecMAT)

    !C
    !C-- INTERFACE data EXCHANGE
    S_TIME = HECMW_WTIME()
    !call hecmw_update_m_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)
    call hecmw_update_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)

    E_TIME = HECMW_WTIME()
    Tcomm = Tcomm + E_TIME - S_TIME

    !deallocate (H, WW, SS)
    deallocate (vecR)
    deallocate (workPC)
    deallocate (u  )
    deallocate (c  )
    deallocate (uin)
    deallocate (cin)
    call hecmw_precond_clear(hecMAT)

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

  end subroutine  hecmw_solve_GMRESREN

end module     hecmw_solver_GMRESREN
