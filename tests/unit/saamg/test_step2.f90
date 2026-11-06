!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 2 (block-diagonal D, Lanczos lambda_max, Chebyshev)
!
! Verifies, against LAPACK dsyev references:
!   * node-block diagonal D : D-matvec, and D(D^{-1}x)=x for nb=3 blocks;
!   * Lanczos : estimate <= true lambda_max and covered by the 1.1 safety;
!   * Chebyshev smoother : strong damping of a high-frequency eigenmode, weak
!     damping of a low-frequency one (i.e. it behaves like a smoother).
!-------------------------------------------------------------------------------
program test_step2
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_smoother
  implicit none

  integer(kind=kint) :: npass, nfail
  npass = 0; nfail = 0

  call test_blockdiag(npass, nfail)
  call test_lanczos(npass, nfail)
  call test_chebyshev(npass, nfail)

  write(*,'(/,a,i0,a,i0,a)') '=== Step 2: ', npass, ' passed, ', nfail, ' failed ==='
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

  !> 1D Laplacian tridiag(-1,2,-1), size n, nb=1.
  subroutine build_laplacian(n, A, dense)
    integer(kind=kint), intent(in)  :: n
    type(hecmwST_saamg_bcsr),     intent(out) :: A
    real(kind=kreal),   intent(out) :: dense(n,n)
    integer(kind=kint) :: i
    dense = 0.0d0
    do i = 1, n
      dense(i,i) = 2.0d0
    end do
    do i = 2, n
      dense(i,i-1) = -1.0d0
      dense(i-1,i) = -1.0d0
    end do
    call hecmw_saamg_bcsr_from_dense(dense, n, 1, A)
  end subroutine build_laplacian

  !> Symmetric eigen-decomposition (ascending) via LAPACK dsyev.
  subroutine sym_eig(M, n, w, V)
    integer(kind=kint), intent(in)  :: n
    real(kind=kreal),   intent(in)  :: M(n,n)
    real(kind=kreal),   intent(out) :: w(n), V(n,n)
    real(kind=kreal), allocatable :: work(:)
    integer(kind=kint) :: lwork, info
    external :: dsyev
    V = M
    lwork = 8*n
    allocate(work(lwork))
    call dsyev('V', 'U', n, V, n, w, work, lwork, info)
    if (info /= 0) then; write(*,*) 'dsyev info=', info; stop 1; end if
    deallocate(work)
  end subroutine sym_eig

  !-----------------------------------------------------------------------------
  subroutine test_blockdiag(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: A
    type(hecmwST_saamg_blockdiag) :: D
    real(kind=kreal) :: dense(6,6), x(6), y(6), z(6), bd(6,6)
    integer(kind=kint) :: i

    ! two SPD 3x3 diagonal blocks + symmetric off-diagonal coupling
    dense = 0.0d0
    dense(1:3,1:3) = reshape([4.d0,1.d0,0.d0, 1.d0,3.d0,1.d0, 0.d0,1.d0,2.d0],[3,3])
    dense(4:6,4:6) = reshape([5.d0,0.d0,1.d0, 0.d0,4.d0,0.d0, 1.d0,0.d0,3.d0],[3,3])
    dense(1,4) = 0.5d0; dense(4,1) = 0.5d0
    dense(3,6) = -0.3d0; dense(6,3) = -0.3d0
    call hecmw_saamg_bcsr_from_dense(dense, 6, 3, A)

    ! block-diagonal reference (zero the off-block entries)
    bd = dense
    bd(1:3,4:6) = 0.0d0; bd(4:6,1:3) = 0.0d0

    call hecmw_saamg_blockdiag_setup(A, D)

    do i = 1, 6
      x(i) = real(i,kreal) - 3.5d0
    end do
    call hecmw_saamg_blockdiag_matvec(D, x, y)
    call check('D-matvec vs block-diag dense', &
         maxval(abs(y - matmul(bd, x))) < 1.0d-13, npass, nfail)

    ! D ( D^{-1} x ) = x
    call hecmw_saamg_blockdiag_apply(D, x, z)   ! z = D^{-1} x
    call hecmw_saamg_blockdiag_matvec(D, z, y)  ! y = D z
    call check('D D^{-1} x = x (nb=3)', maxval(abs(y - x)) < 1.0d-12, npass, nfail)

    call hecmw_saamg_blockdiag_free(D)
    call hecmw_saamg_bcsr_free(A)
  end subroutine test_blockdiag

  !-----------------------------------------------------------------------------
  subroutine test_lanczos(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: n = 30
    type(hecmwST_saamg_bcsr) :: A
    type(hecmwST_saamg_blockdiag) :: D
    real(kind=kreal) :: dense(n,n), M(n,n), w(n), V(n,n)
    real(kind=kreal) :: est, true_max

    call build_laplacian(n, A, dense)
    call hecmw_saamg_blockdiag_setup(A, D)        ! D = 2 I
    M = 0.5d0 * dense                       ! D^{-1} A
    call sym_eig(M, n, w, V)
    true_max = w(n)

    est = hecmw_saamg_lanczos_lambda_max(A, D, 10)
    write(*,'(a,es12.5,a,es12.5)') '   lanczos est = ', est, ', true lmax = ', true_max
    call check('lanczos est <= true (from below)', est <= true_max*(1.d0+1.d-9), npass, nfail)
    call check('1.1*est covers true',            1.1d0*est >= true_max,        npass, nfail)
    call check('lanczos est within 20% of true',   est >= 0.8d0*true_max,        npass, nfail)

    call hecmw_saamg_blockdiag_free(D)
    call hecmw_saamg_bcsr_free(A)
  end subroutine test_lanczos

  !-----------------------------------------------------------------------------
  subroutine test_chebyshev(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: n = 30
    type(hecmwST_saamg_bcsr) :: A
    type(hecmwST_saamg_blockdiag) :: D
    type(hecmwST_saamg_chebyshev) :: cheb
    real(kind=kreal) :: dense(n,n), M(n,n), w(n), V(n,n)
    real(kind=kreal) :: est, lambda_hat
    real(kind=kreal) :: xhi(n), xlo(n), zero(n)
    real(kind=kreal) :: damp_hi, damp_lo

    call build_laplacian(n, A, dense)
    call hecmw_saamg_blockdiag_setup(A, D)
    M = 0.5d0 * dense
    call sym_eig(M, n, w, V)

    est = hecmw_saamg_lanczos_lambda_max(A, D, 10)
    lambda_hat = 1.1d0 * est                      ! shared safety-scaled estimate
    call hecmw_saamg_cheb_setup(lambda_hat, 2, 20.0d0, cheb)

    zero = 0.0d0
    ! high-frequency eigenmode (largest eigenvalue, in the Chebyshev band)
    xhi = V(:,n)
    call hecmw_saamg_cheb_apply_local(A, D, cheb, zero, xhi, .false.)
    damp_hi = sqrt(dot_product(xhi,xhi))          ! ||V(:,n)||=1 initially

    ! low-frequency eigenmode (smallest eigenvalue, below the band)
    xlo = V(:,1)
    call hecmw_saamg_cheb_apply_local(A, D, cheb, zero, xlo, .false.)
    damp_lo = sqrt(dot_product(xlo,xlo))

    write(*,'(a,es10.3,a,es10.3)') '   cheb damp high = ', damp_hi, &
         ', damp low = ', damp_lo
    call check('Chebyshev damps high mode strongly', damp_hi < 0.2d0, npass, nfail)
    call check('Chebyshev keeps low mode',           damp_lo > 0.6d0, npass, nfail)
    call check('high damped >> low',  damp_lo > 5.0d0*damp_hi,        npass, nfail)

    call hecmw_saamg_blockdiag_free(D)
    call hecmw_saamg_bcsr_free(A)
  end subroutine test_chebyshev

end program test_step2
