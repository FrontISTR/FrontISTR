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

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_solid)         :: fstrSOLID
    type(fstr_eigen)         :: fstrEIG
    type(fstr_tri_diag)      :: Tri
    type(fstr_eigen_vec), pointer :: Q(:)
    integer(kind=kint) :: N, NP, NDOF, NNDOF, NPNDOF
    integer(kind=kint) :: iter, maxiter, nget, ierr
    integer(kind=kint) :: i, j, k, in, jn, kn, ik, it
    integer(kind=kint) :: ig, ig0, is0, ie0, its0, ite0
    real(kind=kreal)   :: t1, t2, tolerance
    real(kind=kreal)   :: alpha, beta, beta0
    real(kind=kreal), allocatable :: s(:), t(:), p(:)
    logical :: is_converge

    N      = hecMAT%N
    NP     = hecMAT%NP
    NDOF   = hecMESH%n_dof
    NNDOF  = N *NDOF
    NPNDOF = NP*NDOF

    allocate(fstrEIG%filter(NPNDOF))
    fstrEIG%filter = 1.0d0

    jn = 0
    do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
      ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      iS0  = hecMESH%node_group%grp_index(ig-1) + 1
      iE0  = hecMESH%node_group%grp_index(ig  )
      it   = fstrSOLID%BOUNDARY_ngrp_type(ig0)
      itS0 = (it - mod(it,10))/10
      itE0 = mod(it,10)

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        if(NDOF < itE0) itE0 = NDOF
        do i = itS0, itE0
          jn = jn + 1
          fstrEIG%filter((in-1)*NDOF+i) = 0.0d0
        enddo
      enddo
    enddo

    call hecmw_allreduce_I1(hecMESH, jn, hecmw_sum)
    if(jn == 0)then
      fstrEIG%is_free = .true.
      if(myrank == 0)then
        write(*,*) '** free modal analysis: shift factor = 0.1'
      endif
    endif

    call hecmw_update_m_R(hecMESH, fstrEIG%filter, NP, NDOF)

    in = 0
    do i = 1, NNDOF
      if(fstrEIG%filter(i) == 1.0d0) in = in + 1
    enddo
    call hecmw_allreduce_I1(hecMESH, in, hecmw_sum)

    fstrEIG%maxiter = fstrEIG%maxiter + 1
    if(in < fstrEIG%maxiter)then
      if(myrank == 0)then
        write(IMSG,*) '** changed maxiter to system matrix size.'
      endif
      fstrEIG%maxiter = in
    endif

    if(in < fstrEIG%nget)then
      fstrEIG%nget = in
    endif

    maxiter   = fstrEIG%maxiter

    allocate(Q(0:maxiter))
    allocate(Q(0)%q(NPNDOF))
    allocate(Q(1)%q(NPNDOF))
    allocate(fstrEIG%eigval(maxiter))
    allocate(fstrEIG%eigvec(NPNDOF, maxiter))
    allocate(Tri%alpha(maxiter))
    allocate(Tri%beta(maxiter))
    allocate(t(NPNDOF))
    allocate(s(NPNDOF))
    allocate(p(NPNDOF))

    fstrEIG%eigval = 0.0d0
    fstrEIG%eigvec = 0.0d0
    t      = 0.0d0
    p      = 0.0d0
    s      = 0.0d0
    Q(0)%q = 0.0d0
    Q(1)%q = 0.0d0
    Tri%alpha = 0.0d0
    Tri%beta  = 0.0d0
    hecMAT%X  = 0.0d0

    call lanczos_set_initial_value(hecMESH, hecMAT, fstrEIG, fstrEIG%eigvec, p, Q(1)%q, Tri%beta(1))

    hecMAT%Iarray(98) = 1 !Assmebly complete
    hecMAT%Iarray(97) = 1 !Need numerical factorization

    if(myrank == 0)then
      write(IMSG,*)
      write(IMSG,*) ' *****   STAGE Begin Lanczos loop     **'
    endif

    do iter = 1, maxiter-1
      !> q = A^{-1} p
      do i = 1, NPNDOF
        hecMAT%B(i) = p(i)
      enddo

      call solve_LINEQ(hecMESH, hecMAT)

      allocate(Q(iter+1)%q(NPNDOF))

      do i = 1, NPNDOF
        t(i) = hecMAT%X(i) * fstrEIG%filter(i)
      enddo

      !> t = t - beta * q_{i-1}
      !> alpha = p * t
      do i = 1, NPNDOF
        t(i) = t(i) - Tri%beta(iter) * Q(iter-1)%q(i)
      enddo

      alpha = 0.0d0
      do i = 1, NNDOF
        alpha = alpha + p(i) * t(i)
      enddo
      call hecmw_allreduce_R1(hecMESH, alpha, hecmw_sum)
      Tri%alpha(iter) = alpha

      !> t = t - alpha * q_i
      do i = 1, NPNDOF
        t(i) = t(i) - Tri%alpha(iter) * Q(iter)%q(i)
      enddo

      !> re-orthogonalization
      s = 0.0d0

      do i = 1, NPNDOF
        s(i) = fstrEIG%mass(i) * t(i)
      enddo

      do j = 0, iter
        t1 = 0.0d0
        do i = 1, NNDOF
          t1 = t1 + Q(j)%q(i) * s(i)
        enddo
        call hecmw_allreduce_R1(hecMESH, t1, hecmw_sum)
        do i = 1, NPNDOF
          t(i) = t(i) - t1 * Q(j)%q(i)
        enddo
      enddo

      !> beta = || {t}^t [M] {t} ||_2
      do i = 1, NPNDOF
        s(i) = fstrEIG%mass(i) * t(i)
      enddo

      beta = 0.0d0
      do i = 1, NNDOF
        beta = beta + s(i) * t(i)
      enddo
      call hecmw_allreduce_R1(hecMESH, beta, hecmw_sum)
      Tri%beta(iter+1) = dsqrt(beta)

      !> p = s / beta
      !> q = t / beta
      beta = 1.0d0/Tri%beta(iter+1)
      do i = 1, NPNDOF
        p(i)           = s(i) * beta
        Q(iter+1)%q(i) = t(i) * beta
      enddo

      fstrEIG%iter = iter
      if(iter == 1) beta0 = Tri%beta(iter+1)

      call tridiag(hecMESH, hecMAT, fstrEIG, Q, Tri, iter, is_converge)

      if(is_converge) exit
    enddo

    do i = 0, iter
      if(associated(Q(i)%q)) deallocate(Q(i)%q)
    enddo
    deallocate(Tri%alpha)
    deallocate(Tri%beta)
    deallocate(t)
    deallocate(s)
    deallocate(p)

    t2 = hecmw_Wtime()

    if(myrank == 0)then
      write(IMSG,*)
      write(IMSG,*) ' *     STAGE Output and postprocessing    **'
      write(idbg,'(a,f10.2)') 'Lanczos loop (sec) :', T2 - T1
    endif

  end subroutine fstr_solve_lanczos

end module m_fstr_EIG_lanczos
