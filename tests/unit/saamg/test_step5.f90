!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 5 (SpGEMM + Galerkin product)
!
! Verifies:
!   * SpGEMM C = A*B against dense matmul (rectangular);
!   * Galerkin coarse operator Ac = P^T A P against a dense reference, its
!     symmetry, and its dimensions, on the scalar grid.
!-------------------------------------------------------------------------------
program test_step5
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_smoother
  use hecmw_precond_SAAMG_aggregate
  use hecmw_precond_SAAMG_prolongation
  implicit none

  integer(kind=kint) :: npass, nfail
  npass = 0; nfail = 0

  call test_spgemm(npass, nfail)
  call test_galerkin(npass, nfail)

  write(*,'(/,a,i0,a,i0,a)') '=== Step 5: ', npass, ' passed, ', nfail, ' failed ==='
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

  !> Build a CSR (nr x nc) holding all entries of a dense array.
  subroutine full_bcsr(d, nr, nc, A)
    real(kind=kreal),      intent(in)  :: d(:,:)
    integer(kind=kint),    intent(in)  :: nr, nc
    type(hecmwST_saamg_bcsr), intent(out) :: A
    integer(kind=kint), allocatable :: ti(:), tj(:)
    real(kind=kreal),   allocatable :: tv(:)
    integer(kind=kint) :: i, j, t
    allocate(ti(nr*nc), tj(nr*nc), tv(nr*nc))
    t = 0
    do i = 1, nr
      do j = 1, nc
        t = t + 1; ti(t) = i; tj(t) = j; tv(t) = d(i,j)
      end do
    end do
    call hecmw_saamg_bcsr_from_triplets(nr, nc, 1, ti, tj, tv, nr*nc, A)
    deallocate(ti, tj, tv)
  end subroutine full_bcsr

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
  subroutine test_spgemm(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: A, B, C
    real(kind=kreal) :: da(3,4), db(4,2), dc(3,2), ref(3,2)
    integer(kind=kint) :: i, j
    integer(kind=8) :: seed
    seed = 999_8
    do j = 1, 4; do i = 1, 3
      seed = seed*6364136223846793005_8 + 1442695040888963407_8
      da(i,j) = real(seed,kreal)/9.223372036854776d18
    end do; end do
    do j = 1, 2; do i = 1, 4
      seed = seed*6364136223846793005_8 + 1442695040888963407_8
      db(i,j) = real(seed,kreal)/9.223372036854776d18
    end do; end do
    call full_bcsr(da, 3, 4, A)
    call full_bcsr(db, 4, 2, B)
    call hecmw_saamg_spgemm(A, B, C)
    call hecmw_saamg_bcsr_to_dense(C, dc)
    ref = matmul(da, db)
    call check('SpGEMM dims', C%n==3 .and. C%ncol==2, npass, nfail)
    call check('SpGEMM vs dense matmul', maxval(abs(dc-ref)) < 1.0d-13, npass, nfail)
    call hecmw_saamg_bcsr_free(A); call hecmw_saamg_bcsr_free(B)
    call hecmw_saamg_bcsr_free(C)
  end subroutine test_spgemm

  !-----------------------------------------------------------------------------
  subroutine test_galerkin(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: nx = 12, ny = 12, m = 1
    type(hecmwST_saamg_bcsr) :: A, phat, p, Ac
    type(hecmwST_saamg_nodegraph) :: g
    type(hecmwST_saamg_blockdiag) :: D
    integer(kind=kint), allocatable :: aggr(:)
    real(kind=kreal), allocatable :: bcoarse(:,:), bfine(:,:)
    real(kind=kreal), allocatable :: ad(:,:), pd(:,:), acd(:,:), acref(:,:)
    integer(kind=kint) :: n, naggr, nc
    real(kind=kreal) :: omega, est, asym

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
    call hecmw_saamg_galerkin_local(A, p, m, Ac)
    nc = naggr*m

    call check('Ac dims',  Ac%n==nc .and. Ac%ncol==nc, npass, nfail)
    call check('Ac block size = m', Ac%nb==m, npass, nfail)

    ! dense reference Ac = P^T A P
    allocate(ad(n,n), pd(n,nc), acd(nc,nc), acref(nc,nc))
    call hecmw_saamg_bcsr_to_dense(A, ad)
    call hecmw_saamg_bcsr_to_dense(p, pd)
    acref = matmul(transpose(pd), matmul(ad, pd))
    call hecmw_saamg_bcsr_to_dense(Ac, acd)
    call check('Galerkin vs dense P^T A P', maxval(abs(acd-acref)) < 1.0d-11, &
         npass, nfail)

    asym = hecmw_saamg_is_symmetric(Ac)
    call check('Ac symmetric', asym < 1.0d-12, npass, nfail)

    write(*,'(a,i0,a,i0,a,i0,a,f5.2)') '   n=', n, ' -> nc=', nc, &
         ', nnz(A)=', A%nnz, ' op.complexity Ac/A=', real(Ac%nnz,kreal)/real(A%nnz,kreal)

    deallocate(bfine, bcoarse, ad, pd, acd, acref, aggr)
    call hecmw_saamg_blockdiag_free(D); call hecmw_saamg_nodegraph_free(g)
    call hecmw_saamg_bcsr_free(A);    call hecmw_saamg_bcsr_free(phat)
    call hecmw_saamg_bcsr_free(p);    call hecmw_saamg_bcsr_free(Ac)
  end subroutine test_galerkin

end program test_step5
