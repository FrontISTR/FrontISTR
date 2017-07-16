!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_solver_GMRES
public :: hecmw_solve_GMRES

contains
subroutine hecmw_solve_GMRES( hecMESH, hecMAT, ITER, RESID, ERROR, Tset, Tsol, Tcomm )

  use hecmw_util
  use m_hecmw_solve_error
  use hecmw_matrix_misc
  use hecmw_solver_misc
  use hecmw_solver_las
  use hecmw_solver_scaling
  use hecmw_precond
  use hecmw_estimate_condition


  implicit none
  type(hecmwST_local_mesh) :: hecMESH
  type(hecmwST_matrix) :: hecMAT
  integer(kind=kint ), intent(inout):: ITER, ERROR
  real   (kind=kreal), intent(inout):: RESID, Tset, Tsol, Tcomm
  integer(kind=kint ) :: N, NP, NDOF, NNDOF, my_rank, ITERlog, TIMElog, MAXIT, ESTCOND
  real   (kind=kreal) :: TOL, S1_TIME, S2_TIME
  real   (kind=kreal), pointer :: B(:), X(:)
  
  integer(kind=kint), parameter :: iR=1, iZP=2, iZQ=3, iW=4, iY=4, iAV=5
  real   (kind=kreal), pointer :: R(:), ZP(:), ZQ(:), W(:), Y(:), AV(:)
  real   (kind=kreal), allocatable,target :: V(:,:), VCOS(:), VSIN(:)
  real   (kind=kreal), allocatable,target :: SS(:), H(:,:), WW(:,:),S(:) 
  integer(kind=kint ) :: NREST, i, j, k
  real   (kind=kreal) :: BNRM2, DNRM2, RNORM
  real   (kind=kreal) :: coef,VAL,DTEMP,AA,BB,R0,SCALE,RR

  call hecmw_barrier(hecMESH)
  S1_TIME= HECMW_WTIME()

  N = hecMAT%N
  NP = hecMAT%NP
  NDOF = hecMAT%NDOF
  NNDOF = N * NDOF
  my_rank = hecMESH%my_rank
  X => hecMAT%X
  B => hecMAT%B

  ITERlog = hecmw_mat_get_iterlog( hecMAT )
  TIMElog = hecmw_mat_get_timelog( hecMAT )
  MAXIT   = hecmw_mat_get_iter   ( hecMAT )
  TOL     = hecmw_mat_get_resid  ( hecMAT )
  ESTCOND = hecmw_mat_get_estcond( hecMAT )
  NREST   = hecmw_mat_get_nrest  ( hecMAT )
  if (NREST >= NDOF*NP-1) NREST = NDOF*NP-2

  ERROR= 0
  allocate (H (NREST+1,NREST));
  allocate (V (NDOF*NP,NREST+1))  
  allocate (S (NREST+1))
  allocate (SS(NREST));
  allocate (VCOS(NREST));
  allocate (VSIN(NREST)); 
  allocate (WW(NDOF*NP,5 ))
  WW=0.0d0
  R  => WW(:, iR);
  ZP => WW(:,iZP);
  ZQ => WW(:,iZQ);
  W  => WW(:, iW);
  Y  => WW(:, iY);
  AV => WW(:,iAV);

  call hecmw_solver_scaling_fw(hecMESH, hecMAT, Tcomm)
  call hecmw_precond_setup(hecMAT, hecMESH, 0)

  !C | {r}= {b} - [A]{x0} |
  call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
  call hecmw_InnerProduct_R(hecMESH, NDOF, B, B, BNRM2, Tcomm)
  if (BNRM2.eq.0.d0) then
    iter = 0
    MAXIT = 0
    RESID = 0.d0
    X = 0.d0
  endif

  Tset = HECMW_WTIME() - S1_TIME
  if (TIMElog.eq.2) then
    write(*,*) 'Time solver setup'
    call hecmw_print_time_statistics(hecMESH, Tset)
  endif

  call hecmw_barrier(hecMESH)
  S1_TIME= HECMW_WTIME()
  ITER= 0

  OUTER: do
  !C************************************************ GMRES Iteration
    I= 0
    !C | {v1}= {r}/|r| |
    call hecmw_InnerProduct_R(hecMESH, NDOF, R, R, DNRM2, Tcomm)
    if (DNRM2 == 0.d0) exit ! converged
    RNORM= dsqrt(DNRM2)
    coef= 1.0d0/RNORM
    do j= 1, NNDOF
      V(j,1)= WW(j,iR) * coef
    enddo
    !C | {s}= |r|{e1} |
    S(1) = RNORM
    do j = 2, NREST+1
      S(j) = 0.0d0
    enddo

    !C************************************************ GMRES(m) restart
    do I = 1, NREST
    ITER= ITER + 1

  !C | {w}= [A][Minv]{v} |
      call hecmw_precond_apply(hecMESH, hecMAT, V(:,I), ZQ, ZP, Tcomm)
      call hecmw_matvec(hecMESH, hecMAT, ZQ, W, Tcomm)

  !C ORTH. BASIS by GRAMM-SCHMIDT |
  !C   Construct the I-th column of the upper Hessenberg matrix H
  !C   using the Gram-Schmidt process on V and W.
      do j= 1, I
        call hecmw_InnerProduct_R(hecMESH, NDOF, W, V(:,j), VAL, Tcomm)
        do k= 1, NNDOF
          WW(k,iW)= WW(k,iW) - VAL * V(k,j)
        enddo
        H(j,I)= VAL
      enddo

      call hecmw_InnerProduct_R(hecMESH, NDOF, W, W, VAL, Tcomm)
      if (VAL == 0.d0) exit ! converged

      H(I+1,I)= dsqrt(VAL)
      coef= 1.0d0 / H(I+1,I)
      do j= 1, NNDOF
        V(j,I+1)= WW(j,iW) * coef
      enddo

  !C GIVENS ROTARION 
  !C-- Plane Rotation
      do k = 1, I-1
        DTEMP   =  VCOS(k)*H(k  ,I) + VSIN(k)*H(k+1,I)
        H(k+1,I)=  VCOS(k)*H(k+1,I) - VSIN(k)*H(k  ,I)
        H(k  ,I)= DTEMP
      enddo

  !C-- Construct Givens Plane Rotation
      AA = H(I  ,I)
      BB = H(I+1,I)
      R0= BB
      if (dabs(AA).gt.dabs(BB)) R0= AA
      SCALE= dabs(AA) + dabs(BB)

      if (SCALE.ne.0.d0) then
        RR= dsign(1.d0,R0)* SCALE * dsqrt((AA/SCALE)**2+(BB/SCALE)**2)
        VCOS(I)= AA/RR
        VSIN(I)= BB/RR
       else
        VCOS(I)= 1.d0
        VSIN(I)= 0.d0
      endif

      !C-- Plane Rotation
      DTEMP    = VCOS(I)*H(I  ,I) + VSIN(I)*H(I+1,I)
      H (I+1,I)= VCOS(I)*H(I+1,I) - VSIN(I)*H(I  ,I)
      H (I  ,I)= DTEMP

      DTEMP    = VCOS(I)*S(I) + VSIN(I)*S(I+1)
      S(I+1)= VCOS(I)*S(I+1) - VSIN(I)*S(I)
      S(I)= DTEMP
      RESID = dabs ( S(I+1))/dsqrt(BNRM2)

      if (my_rank.eq.0 .and. ITERlog.eq.1) write (*, '(2i8, 1pe16.6)') iter,I+1, RESID

      if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
        if (mod(ITER,ESTCOND) == 0) call hecmw_estimate_condition_GMRES(I, H)
      endif

      if ( RESID .le. TOL ) then
      !C-- [H]{y}= {s_tld}
         do j= 1, I
           SS(j)= S(j)
         enddo
         WW(I,iY)= S(I) / H(I,I)

         do j= I-1, 1, -1
           do k= I, j+1, -1
             SS(j)= SS(j) - H(j,k)*WW(k,iY)
           enddo
           WW(j,iY)= SS(j) / H(j,j)
         enddo
      !C-- {x}= {x} + {y}{V}
         do j= 1, NNDOF
           WW(j,iAV)= 0.d0
         enddo
         do j= 1, I
           do k= 1, NNDOF
             WW(k,iAV)= WW(k,iAV) + WW(j,iY)*V(k,j)
           enddo
         enddo
         call hecmw_precond_apply(hecMESH, hecMAT, AV, ZQ, ZP, Tcomm)
         do j= 1, NNDOF
           X(j)= X(j) + WW(j,iZQ)
         enddo
         exit OUTER
      endif
      if ( ITER.gt.MAXIT ) then
        ERROR = HECMW_SOLVER_ERROR_NOCONV_MAXIT
        exit OUTER
      end if
    end do

  !C CURRENT SOLUTION |
    !C-- [H]{y}= {s_tld}
    do j= 1, NREST
      SS(j)= S(j)
    enddo
    WW(NREST,iY)= SS(NREST) / H(NREST,NREST)

    do j= NREST-1, 1, -1
      do k= NREST, j+1, -1
        SS(j)= SS(j) - H(j,k)*WW(k,iY)
      enddo
      WW(j,iY)= SS(j) / H(j,j)
    enddo

    !C-- {x}= {x} + {y}{V}
    do j= 1, NNDOF
      WW(j,iAV)= 0.d0
    enddo

    do j= 1, NREST
      do k= 1, NNDOF
        WW(k,iAV)= WW(k,iAV) + WW(j,iY)*V(k,j)
      enddo
    enddo

    call hecmw_precond_apply(hecMESH, hecMAT, AV, ZQ, ZP, Tcomm)

    do j= 1, NNDOF
      X(j)= X(j) + ZQ(j)
    enddo

    !C-- Compute residual vector R, find norm, then check for tolerance.
    call hecmw_matresid(hecMESH, hecMAT, X, B, R, Tcomm)
    call hecmw_InnerProduct_R(hecMESH, NDOF, R, R, DNRM2, Tcomm)

    RESID= dsqrt(DNRM2/BNRM2)
    !S(I+1) = RESID

    if ( RESID.le.TOL )   exit OUTER
    if ( ITER .gt.MAXIT ) then
      ERROR = HECMW_SOLVER_ERROR_NOCONV_MAXIT
      exit OUTER
    end if

  !C-- RESTART
  end do OUTER


  !C-- iteration FAILED
  if (ERROR == HECMW_SOLVER_ERROR_NOCONV_MAXIT) then

    !C-- [H]{y}= {s_tld}
    do j= 1, I
      SS(j)= S(j)
    enddo
    WW(I,iY)= SS(I) / H(I,I)

    do j= I-1, 1, -1
      do k= I, j+1, -1
        SS(j)= SS(j) - H(j,k)*WW(k,iY)
      enddo
      WW(j,iY)= SS(j) / H(j,j)
    enddo

    !C-- {x}= {x} + {y}{V}
    do j= 1, NNDOF
      WW(j,iAV)= 0.d0
    enddo

    do j= 1, I
      do k= 1, NNDOF
        WW(k,iAV)= WW(k,iAV) + WW(j,iY)*V(k  ,j)
      enddo
    enddo

    call hecmw_precond_apply(hecMESH, hecMAT, AV, ZQ, ZP, Tcomm)

    do j= 1, NNDOF
      X(j)= X(j) + WW(j,iZQ)
    enddo
  end if

  call hecmw_solver_scaling_bk(hecMAT)

  if (ESTCOND /= 0 .and. hecMESH%my_rank == 0) then
    call hecmw_estimate_condition_GMRES(I, H)
  endif

  !C-- INTERFACE data EXCHANGE
  S2_TIME = HECMW_WTIME()
  call hecmw_update_m_R (hecMESH, X, hecMAT%NP, hecMAT%NDOF)
  Tcomm = Tcomm + HECMW_WTIME() - S2_TIME

  deallocate (H, SS, S)
  deallocate (V,VCOS,VSIN)
  deallocate (WW)

  !call hecmw_precond_clear(hecMAT)

  Tsol = HECMW_WTIME() - S1_TIME
  if (TIMElog.eq.2) then
    write(*,*) 'Time solver iterations'
    call hecmw_print_time_statistics(hecMESH, Tsol)
  endif

end subroutine  hecmw_solve_GMRES
end module     hecmw_solver_GMRES
