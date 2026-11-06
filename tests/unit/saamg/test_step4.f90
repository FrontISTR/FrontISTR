!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 4 (tentative QR P-hat, smoothing, coarse B)
!
! Verifies:
!   * matrix primitives: from_triplets (dup summation), transpose, to_dense;
!   * tentative prolongator identities  P-hat * B_coarse = B_fine  and
!     P-hat^T P-hat = I  (random nb=3 block B, and the scalar grid B=ones);
!   * smoothed prolongator  P = (I - omega D^{-1} A) P-hat  matches a dense
!     reference to machine precision.
!-------------------------------------------------------------------------------
program test_step4
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_smoother
  use hecmw_precond_SAAMG_aggregate
  use hecmw_precond_SAAMG_prolongation
  implicit none

  integer(kind=kint) :: npass, nfail
  npass = 0; nfail = 0

  call test_primitives(npass, nfail)
  call test_tentative_block(npass, nfail)
  call test_tentative_grid(npass, nfail)
  call test_smoothing(npass, nfail)

  write(*,'(/,a,i0,a,i0,a)') '=== Step 4: ', npass, ' passed, ', nfail, ' failed ==='
  if (nfail > 0) stop 1

contains

  subroutine check(name, cond, npass, nfail)
    character(len=*),   intent(in)    :: name
    logical,            intent(in)    :: cond
    integer(kind=kint), intent(inout) :: npass, nfail
    if (cond) then
      npass = npass + 1; write(*,'(a,a)') '  [PASS] ', name
    else
      nfail = nfail + 1; write(*,'(a,a)') '  [FAIL] ', name
    end if
  end subroutine check

  subroutine build_grid(nx, ny, A)
    integer(kind=kint), intent(in)  :: nx, ny
    type(hecmwST_saamg_bcsr), intent(out) :: A
    real(kind=kreal), allocatable :: dense(:,:)
    integer(kind=kint) :: ix, iy, id, n
    n = nx*ny
    allocate(dense(n,n)); dense = 0.0d0
    do iy = 1, ny
      do ix = 1, nx
        id = (iy-1)*nx + ix
        dense(id,id) = 4.0d0
        if (ix > 1)  dense(id, id-1)  = -1.0d0
        if (ix < nx) dense(id, id+1)  = -1.0d0
        if (iy > 1)  dense(id, id-nx) = -1.0d0
        if (iy < ny) dense(id, id+nx) = -1.0d0
      end do
    end do
    call hecmw_saamg_bcsr_from_dense(dense, n, 1, A)
    deallocate(dense)
  end subroutine build_grid

  !-----------------------------------------------------------------------------
  subroutine test_primitives(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: A, At
    integer(kind=kint) :: ti(4), tj(4)
    real(kind=kreal) :: tv(4), d22(2,2), d32(3,2), d23(2,3)
    real(kind=kreal) :: ref22(2,2), ref23(2,3)

    ! from_triplets with a duplicate (1,1)
    ti = [1,1,2,1]; tj = [1,2,2,1]; tv = [1.0d0,2.0d0,3.0d0,0.5d0]
    call hecmw_saamg_bcsr_from_triplets(2, 2, 1, ti, tj, tv, 4, A)
    call hecmw_saamg_bcsr_to_dense(A, d22)
    ref22 = reshape([1.5d0,0.0d0, 2.0d0,3.0d0],[2,2])  ! col-major: (1,1)=1.5,(2,1)=0,(1,2)=2,(2,2)=3
    call check('from_triplets sums dups', maxval(abs(d22-ref22)) < 1.0d-14, npass, nfail)
    call hecmw_saamg_bcsr_free(A)

    ! transpose of a 3x2 matrix
    ti = [1,2,3,2]; tj = [1,1,2,2]; tv = [1.0d0,2.0d0,3.0d0,4.0d0]
    call hecmw_saamg_bcsr_from_triplets(3, 2, 1, ti, tj, tv, 4, A)
    call hecmw_saamg_bcsr_to_dense(A, d32)
    call hecmw_saamg_bcsr_transpose(A, At)
    call check('transpose dims', At%n == 2 .and. At%ncol == 3, npass, nfail)
    call hecmw_saamg_bcsr_to_dense(At, d23)
    ref23 = transpose(d32)
    call check('transpose values', maxval(abs(d23-ref23)) < 1.0d-14, npass, nfail)
    call hecmw_saamg_bcsr_free(A); call hecmw_saamg_bcsr_free(At)
  end subroutine test_primitives

  !-----------------------------------------------------------------------------
  subroutine test_tentative_block(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: nnode = 6, nb = 3, m = 4
    integer(kind=kint), parameter :: n = nnode*nb, naggr = 2
    type(hecmwST_saamg_bcsr) :: phat
    real(kind=kreal), allocatable :: bcoarse(:,:)
    real(kind=kreal) :: bfine(n,m), pd(n,naggr*m), pb(n,m), ptp(naggr*m,naggr*m)
    integer(kind=kint) :: aggr(nnode), i, j
    integer(kind=8) :: seed

    ! deterministic pseudo-random B (generically full column rank)
    seed = 12345_8
    do j = 1, m
      do i = 1, n
        seed = seed*6364136223846793005_8 + 1442695040888963407_8
        bfine(i,j) = real(seed, kreal)/9.223372036854776d18
      end do
    end do
    aggr = [1,1,1,2,2,2]

    call hecmw_saamg_tentative(bfine, n, m, nb, aggr, naggr, phat, bcoarse)
    call check('phat dims', phat%n==n .and. phat%ncol==naggr*m, npass, nfail)

    call hecmw_saamg_bcsr_to_dense(phat, pd)
    pb  = matmul(pd, bcoarse)
    call check('P-hat B_coarse = B_fine (nb=3)', maxval(abs(pb-bfine)) < 1.0d-11, &
         npass, nfail)
    ptp = matmul(transpose(pd), pd)
    do i = 1, naggr*m
      ptp(i,i) = ptp(i,i) - 1.0d0
    end do
    call check('P-hat^T P-hat = I (nb=3)', maxval(abs(ptp)) < 1.0d-12, npass, nfail)

    call hecmw_saamg_bcsr_free(phat)
    deallocate(bcoarse)
  end subroutine test_tentative_block

  !-----------------------------------------------------------------------------
  subroutine test_tentative_grid(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: nx = 10, ny = 10, m = 1
    type(hecmwST_saamg_bcsr) :: A, phat
    type(hecmwST_saamg_nodegraph) :: g
    integer(kind=kint), allocatable :: aggr(:)
    real(kind=kreal), allocatable :: bcoarse(:,:), bfine(:,:), pd(:,:), pb(:,:), ptp(:,:)
    integer(kind=kint) :: n, naggr, i

    call build_grid(nx, ny, A)
    n = A%n
    call hecmw_saamg_build_nodegraph(A, g)
    call hecmw_saamg_aggregate(g, 3, 9, aggr, naggr)

    allocate(bfine(n,m)); bfine = 1.0d0          ! constant near-kernel
    call hecmw_saamg_tentative(bfine, n, m, 1, aggr, naggr, phat, bcoarse)

    allocate(pd(n, naggr*m), pb(n,m), ptp(naggr*m, naggr*m))
    call hecmw_saamg_bcsr_to_dense(phat, pd)
    pb = matmul(pd, bcoarse)
    call check('P-hat B_coarse = B_fine (grid)', maxval(abs(pb-bfine)) < 1.0d-12, &
         npass, nfail)
    ptp = matmul(transpose(pd), pd)
    do i = 1, naggr*m
      ptp(i,i) = ptp(i,i) - 1.0d0
    end do
    call check('P-hat^T P-hat = I (grid)', maxval(abs(ptp)) < 1.0d-12, npass, nfail)

    deallocate(bfine, bcoarse, pd, pb, ptp, aggr)
    call hecmw_saamg_nodegraph_free(g); call hecmw_saamg_bcsr_free(phat)
    call hecmw_saamg_bcsr_free(A)
  end subroutine test_tentative_grid

  !-----------------------------------------------------------------------------
  subroutine test_smoothing(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: nx = 10, ny = 10, m = 1
    type(hecmwST_saamg_bcsr) :: A, phat, p
    type(hecmwST_saamg_nodegraph) :: g
    type(hecmwST_saamg_blockdiag) :: D
    integer(kind=kint), allocatable :: aggr(:)
    real(kind=kreal), allocatable :: bcoarse(:,:), bfine(:,:)
    real(kind=kreal), allocatable :: ad(:,:), pdhat(:,:), pdref(:,:), pd(:,:), dinvA(:,:)
    integer(kind=kint) :: n, naggr, nc
    real(kind=kreal) :: omega, est

    call build_grid(nx, ny, A)
    n = A%n
    call hecmw_saamg_build_nodegraph(A, g)
    call hecmw_saamg_aggregate(g, 3, 9, aggr, naggr)
    call hecmw_saamg_blockdiag_setup(A, D)
    est = hecmw_saamg_lanczos_lambda_max(A, D, 10)
    omega = (4.0d0/3.0d0) / (1.1d0*est)

    allocate(bfine(n,m)); bfine = 1.0d0
    call hecmw_saamg_tentative(bfine, n, m, 1, aggr, naggr, phat, bcoarse)
    call hecmw_saamg_smooth_prolongator(A, D, omega, phat, p)
    nc = phat%ncol

    ! dense reference: P = (I - omega D^{-1} A) P-hat,  D = 4 I here
    allocate(ad(n,n), pdhat(n,nc), pdref(n,nc), pd(n,nc), dinvA(n,n))
    call hecmw_saamg_bcsr_to_dense(A, ad)
    call hecmw_saamg_bcsr_to_dense(phat, pdhat)
    dinvA = ad / 4.0d0                    ! D = diag(4)
    pdref = pdhat - omega * matmul(dinvA, pdhat)
    call hecmw_saamg_bcsr_to_dense(p, pd)

    call check('smoothed P matches dense ref', maxval(abs(pd-pdref)) < 1.0d-11, &
         npass, nfail)
    write(*,'(a,es10.3,a,i0,a,i0)') '   omega=', omega, ', nnz(phat)=', phat%nnz, &
         ', nnz(P)=', p%nnz

    deallocate(bfine, bcoarse, ad, pdhat, pdref, pd, dinvA, aggr)
    call hecmw_saamg_blockdiag_free(D); call hecmw_saamg_nodegraph_free(g)
    call hecmw_saamg_bcsr_free(A); call hecmw_saamg_bcsr_free(phat)
    call hecmw_saamg_bcsr_free(p)
  end subroutine test_smoothing

end program test_step4
