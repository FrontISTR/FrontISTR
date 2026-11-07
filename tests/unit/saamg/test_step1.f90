!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 1 (internal CSR + generic matvec)
!
! Built against the SA-AMG core modules only, without HEC-MW.  Verifies the
! point-CSR matrix type, generic matvec, the symmetry helper, and exercises a
! tiny unpreconditioned CG so later steps have a known-good base to build on.
!-------------------------------------------------------------------------------
program test_step1
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  implicit none

  integer(kind=kint), parameter :: n = 10
  type(hecmwST_saamg_bcsr) :: A
  real(kind=kreal) :: dense(n,n)
  real(kind=kreal) :: x(n), y(n), b(n)
  real(kind=kreal) :: rel, resnrm
  integer(kind=kint) :: i, iter
  integer(kind=kint) :: npass, nfail

  npass = 0; nfail = 0

  ! 1D Laplacian, tridiag(-1, 2, -1), SPD
  dense = 0.0d0
  do i = 1, n
    dense(i,i) = 2.0d0
  end do
  do i = 2, n
    dense(i,i-1) = -1.0d0
    dense(i-1,i) = -1.0d0
  end do
  call hecmw_saamg_bcsr_from_dense(dense, n, 1, A)

  ! --- structure (1x1 blocks: browptr is the row pointer, nnzb == nnz) ---
  call check('nnz count', A%nnz == (3*n - 2), npass, nfail)
  call check('browptr(1)=1', A%browptr(1) == 1, npass, nfail)
  call check('browptr(n+1)=nnzb+1', A%browptr(n+1) == A%nnzb+1, npass, nfail)

  ! --- matvec : A * ones = [1,0,...,0,1] ---
  x = 1.0d0
  call hecmw_saamg_matvec(A, x, y)
  call check('matvec ones first', abs(y(1) - 1.0d0) < 1.0d-14, npass, nfail)
  call check('matvec ones last',  abs(y(n) - 1.0d0) < 1.0d-14, npass, nfail)
  call check('matvec ones inner', maxval(abs(y(2:n-1))) < 1.0d-14, npass, nfail)

  ! --- matvec against full dense reference on a random-ish vector ---
  do i = 1, n
    x(i) = sin(real(i,kreal))
  end do
  call hecmw_saamg_matvec(A, x, y)
  call check('matvec vs dense', maxval(abs(y - matmul(dense, x))) < 1.0d-13, &
       npass, nfail)

  ! --- symmetry helper ---
  rel = hecmw_saamg_is_symmetric(A)
  call check('symmetric matrix -> rel~0', rel < 1.0d-14, npass, nfail)

  ! --- tiny unpreconditioned CG : solve A x = b , b = ones ---
  b = 1.0d0
  call cg_plain(A, b, x, iter, resnrm)
  call check('CG converged', resnrm < 1.0d-10, npass, nfail)
  write(*,'(a,i0,a,es10.3)') '   CG iters = ', iter, ', final rel.res = ', resnrm

  call hecmw_saamg_bcsr_free(A)

  write(*,'(/,a,i0,a,i0,a)') '=== Step 1: ', npass, ' passed, ', nfail, ' failed ==='
  if (nfail > 0) stop 1

contains

  subroutine check(name, cond, npass, nfail)
    character(len=*),   intent(in)    :: name
    logical,            intent(in)    :: cond
    integer(kind=kint), intent(inout) :: npass, nfail
    if (cond) then
      npass = npass + 1
      write(*,'(a,a)') '  [PASS] ', name
    else
      nfail = nfail + 1
      write(*,'(a,a)') '  [FAIL] ', name
    end if
  end subroutine check

  !> Plain (unpreconditioned) CG, for the checks of this test only.
  subroutine cg_plain(A, b, x, iter, resnrm)
    type(hecmwST_saamg_bcsr),     intent(in)  :: A
    real(kind=kreal),   intent(in)  :: b(:)
    real(kind=kreal),   intent(out) :: x(:)
    integer(kind=kint), intent(out) :: iter
    real(kind=kreal),   intent(out) :: resnrm
    real(kind=kreal), allocatable :: r(:), p(:), q(:)
    real(kind=kreal) :: rho, rho1, alpha, beta, bnrm, pq
    integer(kind=kint) :: m
    m = A%n
    allocate(r(m), p(m), q(m))
    x = 0.0d0
    r = b
    bnrm = sqrt(dot_product(b, b))
    rho1 = 1.0d0; p = 0.0d0
    do iter = 1, 1000
      rho = dot_product(r, r)
      resnrm = sqrt(rho) / bnrm
      if (resnrm < 1.0d-12) exit
      if (iter == 1) then
        p = r
      else
        beta = rho / rho1
        p = r + beta * p
      end if
      call hecmw_saamg_matvec(A, p, q)
      pq = dot_product(p, q)
      alpha = rho / pq
      x = x + alpha * p
      r = r - alpha * q
      rho1 = rho
    end do
    deallocate(r, p, q)
  end subroutine cg_plain

end program test_step1
