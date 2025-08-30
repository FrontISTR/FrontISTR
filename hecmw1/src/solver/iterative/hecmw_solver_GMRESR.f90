!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

!C***
!C*** module hecmw_solver_GMRESR
!C***
!
module hecmw_solver_GMRESR

  public :: hecmw_solve_GMRESR

contains
  !C
  !C*** hecmw_solve_GMRESR
  !C
  subroutine hecmw_solve_GMRESR( hecMESH,  hecMAT, ITER, RESID, error, &
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

    integer(kind=kint),dimension(:), allocatable :: idxBFGS
    real(kind=kreal), dimension(:)  ,  allocatable :: vecR,workPC,tmpVecBFGS,rho,alpha
    real(kind=kreal), dimension(:,:),  allocatable :: u,c,uin,cin,sBFGS,yBFGS

    integer(kind=kint ) :: MAXIT, NREST,NBFGS

    real   (kind=kreal) :: TOL

    real   (kind=kreal)   ZERO, ONE
    parameter ( ZERO = 0.0D+0, ONE = 1.0D+0 )

    integer(kind=kint ) :: NRK,i,k,kk,jj,INFO,ik,iOrth,idx,tmpIdx,iBFGS
    integer(kind=kint ) :: IROW
    real   (kind=kreal) :: S_TIME,E_TIME,S1_TIME,E1_TIME
    real   (kind=kreal) :: LDH,LDW,BNRM2,DNRM2,RNORM
    real   (kind=kreal) :: COMMtime,COMPtime, coef,coef2,val,VCS,VSN,DTEMP,AA,BB,R0,scale,RR
    integer(kind=kint ) :: ESTCOND
    real   (kind=kreal) :: t_max,t_min,t_avg,t_sd

    !! matrix integration for OpenACC
#ifdef _OPENACC
    integer(kind=kint), allocatable :: indexA(:), itemA(:)
    real(kind=kreal), allocatable   :: A(:)
    integer(kind=kint) :: nn, NPA, pre, pp, jS, jE, j
#endif

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
    NBFGS  = hecmw_mat_get_nbfgs( hecMAT )
    ESTCOND = hecmw_mat_get_estcond( hecMAT )

    error= 0

    allocate (vecR(NDOF*NP))
    allocate (workPC(NDOF*NP))
    allocate (u  (NDOF*NP,NREST))
    allocate (c  (NDOF*NP,NREST))
    allocate (uin(NDOF*NP,NREST))
    allocate (cin(NDOF*NP,NREST))

    if(NBFGS>0)then
       allocate (tmpVecBFGS(NDOF*NP))
       allocate (sBFGS(NDOF*NP,NBFGS))
       allocate (yBFGS(NDOF*NP,NBFGS))
       allocate (idxBFGS(NBFGS))
       allocate (rho(NBFGS))
       allocate (alpha(NBFGS))
    endif

    do idx = 1, NBFGS
       idxBFGS(idx) = idx
    enddo

    COMMtime= 0.d0
    COMPtime= 0.d0

    !C
    !C-- SCALING
    call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)

    !C
    !C-- matrix integration for OpenACC
#ifdef _OPENACC
    nn = NDOF * NDOF
    NPA = NP + hecMAT%NPL + hecMAT%NPU
    allocate (indexA(0:NP), itemA(NPA), A(nn * NPA))
    indexA(0) = 0

    pre = 0
    pp = 0
    !$acc parallel loop private(i, j, k, pre, pp, jS, jE)
    do i = 1, NP
      indexA(i) = i + hecMAT%indexL(i) + hecMAT%indexU(i)

      pre = i - 1 + hecMAT%indexU(i - 1)
      jS= hecMAT%indexL(i - 1) + 1
      jE= hecMAT%indexL(i  )
      do j = jS, jE
        pp = pre + j
        itemA(pp) = hecMAT%itemL(j)
        do k = -nn+1, 0
          A(nn * pp + k) = hecMAT%AL(nn * j + k)
        enddo
      enddo

      pp = i + hecMAT%indexU(i - 1) + hecMAT%indexL(i)
      itemA(pp) = i
      do k = -nn+1, 0
        A(nn * pp + k) = hecMAT%D(nn * i + k)
      enddo

      pre = i + hecMAT%indexL(i)
      jS= hecMAT%indexU(i - 1) + 1
      jE= hecMAT%indexU(i  )
      do j = jS, jE
        pp = pre + j
        itemA(pp) = hecMAT%itemU(j)
        do k = -nn+1, 0
          A(nn * pp + k) = hecMAT%AU(nn * j + k)
        enddo
      enddo
    enddo
    !$acc end parallel
#endif

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
    iBFGS=0
    OUTER: do

#ifdef _OPENACC
      call hecmw_matresid_A(hecMESH, hecMAT, indexA, itemA, A, X, B, vecR, Tcomm)
#else
      call hecmw_matresid(hecMESH, hecMAT, X, B, vecR, Tcomm)
#endif
      do I = 1, NREST
        ITER= ITER + 1

        !C Solve M*r = uin(:,1)
        if (iBFGS == 0)then
#ifdef _OPENACC
           call hecmw_precond_apply_A(hecMESH, hecMAT, indexA, itemA, A, vecR, uin(:,1), workPC, Tcomm)
#else
           call hecmw_precond_apply(hecMESH, hecMAT, vecR, uin(:,1), workPC, Tcomm)
#endif
        else
           !$acc kernels
           !$acc loop independent
           do kk= 1, NNDOF
             tmpVecBFGS(kk)= vecR(kk)
           enddo
           !$acc end kernels
           do k = 1,iBFGS
              idx = idxBFGS(k)
              call hecmw_InnerProduct_R(hecMESH, NDOF, sBFGS(:,idx), yBFGS(:,idx), coef, Tcomm)
              rho(k) = 1.0d0 / coef
              call hecmw_InnerProduct_R(hecMESH, NDOF, sBFGS(:,idx), tmpVecBFGS, coef2, Tcomm)
              alpha(k) = rho(k)*coef2
              !$acc kernels
              !$acc loop independent
              do kk= 1, NNDOF
                tmpVecBFGS(kk)= tmpVecBFGS(kk) - alpha(k)*yBFGS(kk,idx)
              enddo
              !$acc end kernels
           enddo
#ifdef _OPENACC
           call hecmw_precond_apply_A(hecMESH, hecMAT, indexA, itemA, A, tmpVecBFGS, uin(:,1), workPC, Tcomm)
#else
           call hecmw_precond_apply(hecMESH, hecMAT, tmpVecBFGS, uin(:,1), workPC, Tcomm)
#endif
           do k = iBFGS,1,-1
              idx = idxBFGS(k)
              call hecmw_InnerProduct_R(hecMESH, NDOF, yBFGS(:,idx), uin(:,1), coef, Tcomm)
              coef2 = rho(k) * coef
              !$acc kernels
              !$acc loop independent
              do kk= 1, NNDOF
                uin(kk,1)= uin(kk,1) + (alpha(k)-coef2)*sBFGS(kk,idx)
              enddo
              !$acc end kernels
           enddo
        endif
        !C cin(:,1) = A*uin(:,1)
#ifdef _OPENACC
        call hecmw_matvec_A(hecMESH, hecMAT, indexA, itemA, A, uin(:,1), cin(:,1), Tcomm)
#else
        call hecmw_matvec(hecMESH, hecMAT, uin(:,1), cin(:,1), Tcomm)
#endif

        do iOrth = 1, I-1
           !C c_{i}^T cin_{i}
           call hecmw_InnerProduct_R(hecMESH, NDOF, c(:,iOrth), cin(:,iOrth), coef, Tcomm)

           !$acc kernels
           !$acc loop independent
           do kk= 1, NNDOF
             cin(kk,iOrth+1)= cin(kk,iOrth) - coef * c(kk,iOrth)
             uin(kk,iOrth+1)= uin(kk,iOrth) - coef * u(kk,iOrth)
           enddo
           !$acc end kernels
        enddo
        call hecmw_InnerProduct_R(hecMESH, NDOF, cin(:,I), cin(:,I), coef, Tcomm)
        coef = 1.0d0 / dsqrt(coef)
        !$acc kernels
        !$acc loop independent
        do kk= 1, NNDOF
          c(kk,I)= coef * cin(kk,I)
          u(kk,I)= coef * uin(kk,I)
        enddo
        !$acc end kernels

        call hecmw_InnerProduct_R(hecMESH, NDOF, c(:,I), vecR, coef, Tcomm)
        !$acc kernels
        !$acc loop independent
        do kk= 1, NNDOF
             x(kk)=    x(kk) + coef*u(kk,I)
          vecR(kk)= vecR(kk) - coef*c(kk,I)
        enddo
        !$acc end kernels

        if (NBFGS > 0)then
           iBFGS = iBFGS + 1
           if (iBFGS == NBFGS+1)then
             tmpIdx = idxBFGS(1)
             do kk = 1, NBFGS-1
                idxBFGS(kk) = idxBFGS(kk+1)
             enddo
             idxBFGS(NBFGS) = tmpIdx
             iBFGS = iBFGS - 1
           endif
           !$acc kernels
           !$acc loop independent
           do kk= 1, NNDOF
             yBFGS(kk,idxBFGS(iBFGS))= coef*c(kk,I)
             sBFGS(kk,idxBFGS(iBFGS))= coef*u(kk,I)
           enddo
           !$acc end kernels
        endif


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

  end subroutine  hecmw_solve_GMRESR

end module     hecmw_solver_GMRESR
