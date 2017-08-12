!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Lanczos iteration calculation
module m_fstr_EIG_lanczos
  contains

!> SOLVE EIGENVALUE PROBLEM
  subroutine fstr_solve_lanczos(hecMESH, hecMAT, fstrSOLID, fstrEIG)
    use m_fstr
    use hecmw_util
    use m_eigen_lib
    use m_fstr_EIG_lanczos_util
    use m_fstr_EIG_tridiag

    implicit none

    type fstr_eigen_vec
      real(kind=kreal), pointer :: q(:) => null()
    end type fstr_eigen_vec

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix    ) :: hecMAT
    type(fstr_solid        ) :: fstrSOLID
    type(fstr_eigen        ) :: fstrEIG
    type(fstr_tri_diag     ) :: Tri
    type(fstr_eigen_vec), pointer :: Q(:)


    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF
    integer(kind=kint) :: iter, maxiter, nget, ierr
    integer(kind=kint) :: i, j, k, in, jn, kn
    integer(kind=kint) :: ig, ig0, is0, ie0, its0, ite0
    real(kind=kreal) :: t1, t2, tolerance
    real(kind=kreal) :: alpha, beta

    integer(kind=kint), allocatable :: iparm(:)
    real(kind=kreal),   allocatable :: s(:), t(:), p(:)

    real(kind=kreal), pointer :: eigvec(:,:)
    real(kind=kreal), pointer :: eigval(:)







    integer(kind=kint) :: kk, jjiter, ppc, jiter, iiter, kiter, it
    integer(kind=kint) :: IOUT,IREOR,itype,iS,iE,ic_type,icel,jS,nn , ii, iii, ik, in1, in2, in3, nstep, istep
    real(kind=kreal)   :: prechk
    real(kind=kreal)   :: PRECHK1,PRECHK2,cchk0,cchk1,cchk,CERR
    real(kind=kreal)   :: aalf, tmp, tmp2, gm, gm2, r1, r2, r3, r4, r5, r6


      REAL(KIND=KREAL), allocatable ::  LLDIAG(:), LNDIAG(:), LSUB(:)
      REAL(KIND=KREAL), allocatable ::  LZMAT(:,:), LNZMAT(:,:)










    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    allocate(fstrEIG%filter(NPNDOF))
    fstrEIG%filter = 1.0d0

    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      iS0  = hecMESH%node_group%grp_index(ig-1) + 1
      iE0  = hecMESH%node_group%grp_index(ig  )
      it   = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      itS0 = (it - mod(it,10))/10
      itE0 = mod(it,10)

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        do i = itS0,itE0
          fstrEIG%filter((in-1)*NDOF+i) = 0.0d0
        enddo
      enddo
    enddo

    call hecmw_update_m_R(hecMESH, fstrEIG%filter, NP, NDOF)

    in = 0
    do i = 1, NNDOF
      if(fstrEIG%filter(i) == 1.0d0) in = in + 1
    enddo
    call hecmw_allreduce_I1(hecMESH, in, hecmw_sum)

    if(in < fstrEIG%maxiter)then
      if(myrank == 0)then
        WRITE(IMSG,*) '*-------------------------------------------*'
        WRITE(IMSG,*) '  WARNING: maxiter exceeds system matrix size. '
        WRITE(IMSG,*) '  Resetting maxiter to system matrix size.'
        WRITE(IMSG,*) '*-------------------------------------------*'
      endif
      fstrEIG%maxiter = in
    endif

    if(in < fstrEIG%nget)then
      fstrEIG%nget = in
    endif

    nget      = fstrEIG%nget
    maxiter   = fstrEIG%maxiter
    tolerance = fstrEIG%tolerance

    allocate( Q(0:maxiter)                    )
    allocate( Q(0)%q(NPNDOF)                  )
    allocate( Q(1)%q(NPNDOF)                  )
    allocate( fstrEIG%eigval(maxiter)         )
    allocate( fstrEIG%eigvec(NPNDOF, maxiter) )
    allocate( Tri%alpha(maxiter+2)            )
    allocate( Tri%beta (maxiter+2)            )
    allocate( iparm(maxiter)                  )
    allocate( t(NPNDOF) )
    allocate( s(NPNDOF)    )
    allocate( p(NPNDOF)    )

    eigval => fstrEIG%eigval
    eigvec => fstrEIG%eigvec
    eigval = 0.0d0
    eigvec = 0.0d0
    t   = 0.0d0
    p      = 0.0d0
    s      = 0.0d0
    Q(0)%q = 0.0d0
    Q(1)%q = 0.0d0
    Tri%alpha = 0.0d0
    Tri%beta  = 0.0d0

    call lanczos_set_initial_value(hecMESH, hecMAT, fstrEIG, eigvec, p, Q(1)%q, Tri%beta(1))

    hecMAT%Iarray(98) = 1   !Assmebly complete
    hecMAT%Iarray(97) = 1   !Need numerical factorization

    if(myrank == 0)then
      write(IMSG,*)
      write(IMSG,*) ' *****   STAGE Begin Lanczos loop     **'
    endif

    do iter=1, maxiter-1
      !> q = A^{-1} p
      do i=1,NPNDOF
        hecMAT%B(i) = p(i)
      enddo

      call solve_LINEQ( hecMESH, hecMAT )

      allocate(Q(iter+1)%q(NPNDOF))

      do i=1, NPNDOF
        t(i) = hecMAT%X(i) * fstrEIG%filter(i)
      enddo

      !> t = t - beta * q_{i-1}
      !> alpha = p * t
      do i=1, NPNDOF
        t(i) = t(i) - Tri%beta(iter) * Q(iter-1)%q(i)
      enddo

      alpha = 0.0d0
      do i=1, NNDOF
        alpha = alpha + p(i) * t(i)
      enddo
      call hecmw_allreduce_R1(hecMESH, alpha, hecmw_sum)
      Tri%alpha(iter) = alpha

      !> t = t - alpha * q_i
      do i=1, NPNDOF
        t(i) = t(i) - Tri%alpha(iter) * Q(iter)%q(i)
      enddo

      !> re-orthogonalization
      s = 0.0d0

      do i=1, NPNDOF
        s(i) = fstrEIG%mass(i) * t(i)
      enddo

      do j=0, iter
        t1 = 0.0d0
        do i=1, NNDOF
          t1 = t1 + Q(j)%q(i) * s(i)
        enddo
        call hecmw_allreduce_R1(hecMESH, t1, hecmw_sum)
        do i=1, NPNDOF
          t(i) = t(i) - t1 * Q(j)%q(i)
        enddo
      enddo

      !> beta = || {t}^t [M] {t} ||_2
      do i=1, NPNDOF
        s(i) = fstrEIG%mass(i) * t(i)
      enddo

      beta = 0.0d0
      do i=1, NNDOF
        beta = beta + s(i) * t(i)
      enddo
      call hecmw_allreduce_R1(hecMESH, beta, hecmw_sum)
      Tri%beta(iter+1) = dsqrt(beta)

      !> p = s / beta
      !> q = t / beta
      beta = 1.0d0/Tri%beta(iter+1)
      do i=1, NPNDOF
        p(i)           = s(i)    * beta
        Q(iter+1)%q(i) = t(i) * beta
      enddo

      fstrEIG%iter = iter

      if(Tri%beta(iter+1) <= tolerance .and. nget <= iter)then
        write(IDBG,*) '*=====Desired convergence was obtained =====*'
        exit
      endif
    enddo

    iter = fstrEIG%iter

      ALLOCATE( LLDIAG(iter)        )
      ALLOCATE( LNDIAG(iter)        ) !Unordered
      ALLOCATE( LSUB(iter)          )
      ALLOCATE( LZMAT(iter,iter)  )
      ALLOCATE( LNZMAT(iter,iter) ) !Unordered

      DO Jiter = 1,iter
        DO Iiter = 1,iter
          LZMAT(Jiter,Iiter)  = 0.0D0
          LNZMAT(Jiter,Iiter) = 0.0D0 !Unordered
        enddo
      enddo

      DO Iiter = 1,iter
        LLDIAG(Iiter) = Tri%alpha(Iiter)
        LZMAT(Iiter,Iiter) = 1.0D0
      enddo

      LSUB(1) = 0.0
      Iiter   = 0
      DO Iiter = 2,iter
        LSUB(Iiter)  = Tri%beta(Iiter)
      enddo

      call TRIDIAG(iter,iter,LLDIAG,LNDIAG,LSUB,LZMAT,LNZMAT,ierr) !Unordered

      DO Iiter = 1, iter
        IF( LLDIAG(Iiter).NE.0.0D0 ) THEN
          eigval(Iiter) = 1.0D0/LLDIAG(Iiter) + fstrEIG%sigma
        ENDIF
      enddo
      call EVSORT(eigval,iparm,iter)

      cchk  = 0.0
      Kiter = nget+2            !Extra 2 values as a safety feature
      IF(iter .LT. Kiter) Kiter = iter

      DO JJiter = 1,Kiter
        Jiter = iparm(JJiter)
        call VECPRO1(prechk1,LZMAT(1,Jiter),LZMAT(1,Jiter),iter)

        IF(prechk1 .GT. 0.) THEN
          prechk1 = SQRT(prechk1)
        ELSE
          prechk1 = 1.
        ENDIF

        cchk1 = (Tri%beta(Jiter))*(LZMAT(iter,Jiter)/prechk)
        IF(cchk .LT. ABS(cchk1)) THEN
          cchk = ABS(cchk1)
          iiter = jiter
          ppc = prechk1
        ENDIF
      enddo

    DO Iiter = 1, iter
      IF(LLDIAG(Iiter).NE.0.0D0) THEN
        eigval(Iiter) = 1.0D0/LLDIAG(Iiter) + fstrEIG%sigma
      ENDIF
    enddo
    call evsort(eigval,iparm,iter)

    eigvec = 0.0
    k = nget
    IF(k .GT. iter) k = iter
    DO kk = 1,k
      kiter = iparm(kk)
      DO jiter = 1,iter
        DO iiter =1,NPNDOF
          eigvec(iiter,kk) = eigvec(iiter,kk) + Q(jiter)%q(iiter)*LZMAT(jiter,kiter)
        enddo
      enddo
    enddo

    do j=1,nget
      prechk1 = maxval(eigvec(:,j))
      call hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)
      if(prechk1 /= 0.0d0)then
        do i = 1, NNDOF
          eigvec(i,j) = eigvec(i,j)/prechk1
        enddo
      endif
    enddo

    if( allocated(LLDIAG) ) DEALLOCATE(LLDIAG)
    if( allocated(LNDIAG) ) DEALLOCATE(LNDIAG)
    if( allocated(LSUB) )   DEALLOCATE(LSUB)
    if( allocated(LZMAT) )  DEALLOCATE(LZMAT)
    if( allocated(LNZMAT) ) DEALLOCATE(LNZMAT)





    do i=0, iter
      if( associated(Q(i)%q) ) deallocate(Q(i)%q)
    enddo

    t2 = hecmw_Wtime()

    if(myrank == 0)then
      write(IMSG,*)
      write(IMSG,*) ' *     STAGE Output and postprocessing    **'
      write(idbg,'(a,f10.2)') 'Lanczos loop (sec) :', T2 - T1
    endif

  end subroutine fstr_solve_lanczos

end module m_fstr_EIG_lanczos
