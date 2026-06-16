!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : smoother building blocks
!!
!! Provides, for one level:
!!   * the node-block diagonal preconditioner D (each nb x nb block inverted
!!     explicitly by the BLAS-free block_inverse; nb=1/2/3 use closed forms)
!!     -- handles indefinite blocks (u-p) too; the original blocks are kept so
!!     that the action y = D x is also available;
!!   * a D-inner-product Lanczos estimate of lambda_max(D^{-1} A) (largest Ritz
!!     value of the Krylov tridiagonal), with a fixed-seed start for
!!     reproducibility (callers apply the safety factor);
!!   * the Chebyshev smoother (fixed polynomial of D^{-1}A), whose action is
!!     matvec + block-diagonal solves only (no sequential dependence -> MPI-safe).
!!
!! lambda_max is computed once per level and shared by the prolongator smoothing
!! (omega) and the Chebyshev interval.
module hecmw_precond_SAAMG_smoother
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_comm, only: hecmw_saamg_abort, hecmw_saamg_check_alloc
  implicit none

  private
  public :: hecmwST_saamg_blockdiag
  public :: hecmw_saamg_blockdiag_setup
  public :: hecmw_saamg_blockdiag_apply
  public :: hecmw_saamg_blockdiag_apply_bcsr
  public :: hecmw_saamg_blockdiag_apply_bcsr_inplace
  public :: hecmw_saamg_blockdiag_matvec
  public :: hecmw_saamg_blockdiag_free
  public :: hecmw_saamg_lanczos_lambda_max
  public :: hecmw_saamg_tridiag_max_eig
  public :: hecmwST_saamg_chebyshev
  public :: hecmw_saamg_cheb_setup
  public :: hecmw_saamg_cheb_apply_local

  !> Node-block diagonal preconditioner: nb x nb blocks, one per node.
  type hecmwST_saamg_blockdiag
    integer(kind=kint) :: nnode = 0   !< number of nodes
    integer(kind=kint) :: nb    = 0   !< dofs per node
    real(kind=kreal),   allocatable :: dblk(:)  !< nb*nb*nnode, col-major (original D)
    real(kind=kreal),   allocatable :: inv(:)   !< nb*nb*nnode, col-major D_i^{-1} of each block
  end type hecmwST_saamg_blockdiag

  !> Chebyshev smoother coefficients for one level (interval [lmin, lmax]).
  type hecmwST_saamg_chebyshev
    integer(kind=kint) :: deg   = 2     !< polynomial degree (sweeps)
    real(kind=kreal)   :: theta = 0.0d0 !< interval centre
    real(kind=kreal)   :: delta = 0.0d0 !< interval half-width
    real(kind=kreal)   :: sigma = 0.0d0 !< theta/delta
  end type hecmwST_saamg_chebyshev

contains

  !-----------------------------------------------------------------------------
  ! Node-block diagonal D
  !-----------------------------------------------------------------------------

  !> Extract the nb x nb diagonal blocks of A, keep them, and store the explicit
  !! inverse D_i^{-1} of each via the BLAS-free block_inverse.  Reads the block-CSR
  !! storage directly: the diagonal block of node inode is the stored block with
  !! bcol == inode (square operator, nb == mb), copied verbatim.  The inverse is
  !! applied later as a small dense matmul (no per-block LAPACK call) AND computed
  !! here without LAPACK, so the whole block-diagonal setup is BLAS-free: no nested
  !! BLAS threading under OpenMP, and offloadable to OpenACC (no host-LAPACK call).
  !! Blocks are tiny (nb = 1/3/6) and well-conditioned (nodal stiffness diagonal).
  subroutine hecmw_saamg_blockdiag_setup(A, D)
    implicit none
    type(hecmwST_saamg_bcsr),      intent(in)  :: A
    type(hecmwST_saamg_blockdiag), intent(out) :: D
    integer(kind=kint) :: nb, nnode, inode, off, t, info, astat

    nb = A%nb
    if (nb <= 0 .or. A%mb /= nb .or. mod(A%n, nb) /= 0) then
      write(*,*) 'hecmw_saamg_blockdiag_setup: A not square-blocked (nb,mb)=', nb, A%mb
      call hecmw_saamg_abort('blockdiag_setup: operator not square-blocked')
    end if
    nnode = A%nbrow

    call hecmw_saamg_blockdiag_free(D)
    D%nb = nb; D%nnode = nnode
    allocate(D%dblk(nb*nb*nnode), D%inv(nb*nb*nnode), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'blockdiag_setup (dblk/inv)')
    D%dblk = 0.0d0

    ! each node's block is inverted independently with the BLAS-free block_inverse
    ! (no per-thread workspace, no LAPACK -> no nested BLAS threading, OpenACC-ready)
    !$omp parallel do default(none) shared(A, D, nb, nnode) private(inode, off, t, info)
    do inode = 1, nnode
      off = (inode-1)*nb*nb
      do t = A%browptr(inode), A%browptr(inode+1)-1
        if (A%bcol(t) == inode) then              ! the diagonal block
          D%dblk(off+1:off+nb*nb) = A%bval((t-1)*nb*nb+1:t*nb*nb)
          exit
        end if
      end do
      call block_inverse(nb, D%dblk(off+1:off+nb*nb), D%inv(off+1:off+nb*nb), info)
      if (info /= 0) then
        write(*,'(a,i0)') 'hecmw_saamg_blockdiag_setup: singular diagonal block at node ', inode
        call hecmw_saamg_abort('blockdiag_setup: singular diagonal block')
      end if
    end do
    !$omp end parallel do
  end subroutine hecmw_saamg_blockdiag_setup

  !> Explicit inverse of a small nb x nb block (col-major) WITHOUT LAPACK, so the
  !! block-diagonal setup stays BLAS-free (no nested threading; OpenACC-offloadable
  !! -- this routine is intended to become an `!$acc routine seq`).  nb = 1/2/3 use
  !! branch-free analytic formulas (the dominant heat/plane/solid cases; ideal on a
  !! GPU with no warp divergence); other nb use partial-pivoted Gauss-Jordan, which
  !! stays robust for indefinite (u-p contact) blocks.  info = 0 ok, 1 singular.
  !! For 3x3 the cofactor inverse is as accurate as pivoted LU (the pivoting only
  !! matters for larger systems), so no symmetry/SPD assumption is needed.
  subroutine block_inverse(nb, a, ainv, info)
    implicit none
    integer(kind=kint), intent(in)  :: nb
    real(kind=kreal),   intent(in)  :: a(nb*nb)        ! col-major a(i,j) = a((j-1)*nb+i)
    real(kind=kreal),   intent(out) :: ainv(nb*nb)     ! col-major
    integer(kind=kint), intent(out) :: info
    real(kind=kreal) :: det, idet, piv, f, aug(nb, 2*nb)
    integer(kind=kint) :: i, j, k, p
    info = 0
    select case (nb)
    case (1)
      if (a(1) == 0.0d0) then; info = 1; return; end if
      ainv(1) = 1.0d0 / a(1)
    case (2)
      det = a(1)*a(4) - a(2)*a(3)
      if (det == 0.0d0) then; info = 1; return; end if
      idet = 1.0d0 / det
      ainv(1) =  a(4)*idet; ainv(2) = -a(2)*idet
      ainv(3) = -a(3)*idet; ainv(4) =  a(1)*idet
    case (3)
      det = a(1)*(a(5)*a(9)-a(8)*a(6)) - a(4)*(a(2)*a(9)-a(8)*a(3)) + a(7)*(a(2)*a(6)-a(5)*a(3))
      if (det == 0.0d0) then; info = 1; return; end if
      idet = 1.0d0 / det
      ainv(1) = (a(5)*a(9)-a(8)*a(6))*idet
      ainv(2) = (a(8)*a(3)-a(2)*a(9))*idet
      ainv(3) = (a(2)*a(6)-a(5)*a(3))*idet
      ainv(4) = (a(7)*a(6)-a(4)*a(9))*idet
      ainv(5) = (a(1)*a(9)-a(7)*a(3))*idet
      ainv(6) = (a(4)*a(3)-a(1)*a(6))*idet
      ainv(7) = (a(4)*a(8)-a(7)*a(5))*idet
      ainv(8) = (a(7)*a(2)-a(1)*a(8))*idet
      ainv(9) = (a(1)*a(5)-a(4)*a(2))*idet
    case default
      ! partial-pivoted Gauss-Jordan on the augmented [a | I] -> [I | a^{-1}]
      do j = 1, nb
        do i = 1, nb
          aug(i, j)    = a((j-1)*nb + i)
          aug(i, nb+j) = 0.0d0
        end do
        aug(j, nb+j) = 1.0d0
      end do
      do k = 1, nb
        p = k                                          ! partial pivot
        do i = k+1, nb
          if (abs(aug(i,k)) > abs(aug(p,k))) p = i
        end do
        if (aug(p,k) == 0.0d0) then; info = 1; return; end if
        if (p /= k) then
          do j = 1, 2*nb
            f = aug(k,j); aug(k,j) = aug(p,j); aug(p,j) = f
          end do
        end if
        piv = aug(k,k)
        do j = 1, 2*nb
          aug(k,j) = aug(k,j) / piv
        end do
        do i = 1, nb
          if (i == k) cycle
          f = aug(i,k)
          if (f == 0.0d0) cycle
          do j = 1, 2*nb
            aug(i,j) = aug(i,j) - f*aug(k,j)
          end do
        end do
      end do
      do j = 1, nb
        do i = 1, nb
          ainv((j-1)*nb + i) = aug(i, nb+j)
        end do
      end do
    end select
  end subroutine block_inverse

  !> Apply y = D^{-1} x (block-wise small dense matmul against the stored inverse).
  subroutine hecmw_saamg_blockdiag_apply(D, x, y)
    implicit none
    type(hecmwST_saamg_blockdiag), intent(in)  :: D
    real(kind=kreal),      intent(in)  :: x(:)
    real(kind=kreal),      intent(out) :: y(:)
    integer(kind=kint) :: inode, base, ii, jj, off, nb
    real(kind=kreal) :: s

    nb = D%nb
    ! GPU (OpenACC) vs CPU (OpenMP) via #ifdef _OPENACC (mutually exclusive; safe with
    ! both -acc and -mp).  Managed memory => no data clauses; inner nb x nb seq.
#ifdef _OPENACC
    !$acc parallel loop gang vector private(base, off, ii, jj, s)
#else
    !$omp parallel do default(none) shared(D, x, y, nb) private(inode, base, off, ii, jj, s)
#endif
    do inode = 1, D%nnode
      base = (inode-1)*nb
      off  = (inode-1)*nb*nb
      !$acc loop seq
      do ii = 1, nb                              ! y_i = sum_j Dinv(i,j) x_j  (col-major)
        s = 0.0d0
        !$acc loop seq
        do jj = 1, nb
          s = s + D%inv(off + (jj-1)*nb + ii) * x(base+jj)
        end do
        y(base+ii) = s
      end do
    end do
#ifndef _OPENACC
    !$omp end parallel do
#endif
  end subroutine hecmw_saamg_blockdiag_apply

  !> Apply DC = D^{-1} C block-row-wise to a block-CSR matrix C (reusing the stored
  !! block LU).  Each node block's nb rows are solved together by one multi-RHS
  !! triangular solve over the union of their columns, so a column present in any
  !! row of the block receives all nb output entries.  Used for the SpGEMM-based
  !! prolongator smoothing P = P-hat - omega D^{-1}(A P-hat).
  subroutine hecmw_saamg_blockdiag_apply_bcsr(D, C, DC)
    implicit none
    type(hecmwST_saamg_blockdiag), intent(in)  :: D
    type(hecmwST_saamg_bcsr),      intent(in)  :: C
    type(hecmwST_saamg_bcsr),      intent(out) :: DC
    integer(kind=kint) :: nb, mb, inode, off, t, b0, ii, jj, cc, astat
    real(kind=kreal) :: s

    nb = D%nb; mb = C%mb
    if (C%nb /= nb .or. C%nbrow /= D%nnode) then
      write(*,*) 'hecmw_saamg_blockdiag_apply_bcsr: block-row mismatch', C%nb, nb, &
           C%nbrow, D%nnode
      call hecmw_saamg_abort('blockdiag_apply_bcsr: block-row mismatch')
    end if

    ! DC = D^{-1} C : same block structure as C, each block DC_IK = D_I^{-1} C_IK
    ! (nb x nb inverse times the nb x mb block, small dense matmul)
    call hecmw_saamg_bcsr_free(DC)
    DC%n = C%n; DC%ncol = C%ncol; DC%nb = nb; DC%mb = mb
    DC%nbrow = C%nbrow; DC%nbcol = C%nbcol; DC%nnzb = C%nnzb
    allocate(DC%browptr(C%nbrow+1), DC%bcol(C%nnzb), DC%bval(nb*mb*C%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'blockdiag_apply_bcsr (DC)')
    DC%browptr = C%browptr; DC%bcol = C%bcol

    !$omp parallel do default(none) shared(C, DC, D, nb, mb) private(inode, off, t, b0, cc, ii, jj, s)
    do inode = 1, C%nbrow
      off = (inode-1)*nb*nb
      do t = C%browptr(inode), C%browptr(inode+1)-1
        b0 = (t-1)*nb*mb
        do cc = 1, mb                            ! DC(i,cc) = sum_j Dinv(i,j) C(j,cc)
          do ii = 1, nb
            s = 0.0d0
            do jj = 1, nb
              s = s + D%inv(off + (jj-1)*nb + ii) * C%bval(b0 + (cc-1)*nb + jj)
            end do
            DC%bval(b0 + (cc-1)*nb + ii) = s
          end do
        end do
      end do
    end do
    !$omp end parallel do
    ! block-only result: the sole consumer (smoothed prolongator) reads bval
  end subroutine hecmw_saamg_blockdiag_apply_bcsr

  !> Apply C <- D^{-1} C in place (each nb x mb block solved against the node's LU).
  !! Same as blockdiag_apply_bcsr but overwrites C, avoiding a full second operator
  !! (important for the large level-1 A P-hat product).
  subroutine hecmw_saamg_blockdiag_apply_bcsr_inplace(D, C)
    implicit none
    type(hecmwST_saamg_blockdiag), intent(in)    :: D
    type(hecmwST_saamg_bcsr),      intent(inout) :: C
    integer(kind=kint) :: nb, mb, inode, off, t, b0, ii, jj, cc
    real(kind=kreal) :: tmp(D%nb), s
    nb = D%nb; mb = C%mb
    if (C%nb /= nb .or. C%nbrow /= D%nnode) then
      write(*,*) 'hecmw_saamg_blockdiag_apply_bcsr_inplace: block-row mismatch', &
           C%nb, nb, C%nbrow, D%nnode
      call hecmw_saamg_abort('blockdiag_apply_bcsr_inplace: block-row mismatch')
    end if
    !$omp parallel do default(none) shared(C, D, nb, mb) private(inode, off, t, b0, cc, ii, jj, s, tmp)
    do inode = 1, C%nbrow
      off = (inode-1)*nb*nb
      do t = C%browptr(inode), C%browptr(inode+1)-1
        b0 = (t-1)*nb*mb
        do cc = 1, mb                            ! col cc <- Dinv * col cc (temp: read before overwrite)
          do jj = 1, nb
            tmp(jj) = C%bval(b0 + (cc-1)*nb + jj)
          end do
          do ii = 1, nb
            s = 0.0d0
            do jj = 1, nb
              s = s + D%inv(off + (jj-1)*nb + ii) * tmp(jj)
            end do
            C%bval(b0 + (cc-1)*nb + ii) = s
          end do
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine hecmw_saamg_blockdiag_apply_bcsr_inplace

  !> Apply y = D x (block-wise dense matvec using the original blocks).
  subroutine hecmw_saamg_blockdiag_matvec(D, x, y)
    implicit none
    type(hecmwST_saamg_blockdiag), intent(in)  :: D
    real(kind=kreal),      intent(in)  :: x(:)
    real(kind=kreal),      intent(out) :: y(:)
    integer(kind=kint) :: inode, base, ii, jj, off
    real(kind=kreal) :: s

    !$omp parallel do default(none) shared(D, x, y) private(inode, base, off, ii, jj, s)
    do inode = 1, D%nnode
      base = (inode-1)*D%nb
      off  = (inode-1)*D%nb*D%nb
      do ii = 1, D%nb
        s = 0.0d0
        do jj = 1, D%nb
          s = s + D%dblk(off + (jj-1)*D%nb + ii) * x(base+jj)
        end do
        y(base+ii) = s
      end do
    end do
    !$omp end parallel do
  end subroutine hecmw_saamg_blockdiag_matvec

  subroutine hecmw_saamg_blockdiag_free(D)
    implicit none
    type(hecmwST_saamg_blockdiag), intent(inout) :: D
    if (allocated(D%dblk)) deallocate(D%dblk)
    if (allocated(D%inv))  deallocate(D%inv)
    D%nnode = 0; D%nb = 0
  end subroutine hecmw_saamg_blockdiag_free

  !-----------------------------------------------------------------------------
  ! Lanczos estimate of lambda_max(D^{-1} A)
  !-----------------------------------------------------------------------------

  !> Largest eigenvalue of a symmetric tridiagonal matrix (diagonal diag(1:m),
  !! sub/super-diagonal offd(1:m-1)) via LAPACK dsterf.  Used by the Lanczos
  !! lambda_max estimators.  Falls back to the largest diagonal entry (a Rayleigh
  !! quotient, hence a lower bound) if dsterf fails to converge.
  function hecmw_saamg_tridiag_max_eig(diag, offd, m) result(lmax)
    implicit none
    real(kind=kreal),   intent(in) :: diag(:), offd(:)
    integer(kind=kint), intent(in) :: m
    real(kind=kreal) :: lmax
    real(kind=kreal), allocatable :: d(:), e(:)
    integer(kind=kint) :: info
#ifdef HECMW_WITH_LAPACK
    external :: dsterf
    if (m <= 1) then; lmax = diag(1); return; end if
    allocate(d(m), e(m-1))
    d(1:m) = diag(1:m); e(1:m-1) = offd(1:m-1)
    call dsterf(m, d, e, info)                 ! eigenvalues ascending in d
    if (info == 0) then; lmax = d(m); else; lmax = maxval(diag(1:m)); end if
    deallocate(d, e)
#else
    ! LAPACK-free fallback (largest diagonal = Rayleigh lower bound).  Unreached in
    ! practice: the backend aborts before setup when LAPACK is absent.
    lmax = maxval(diag(1:m))
#endif
  end function hecmw_saamg_tridiag_max_eig

  !> Estimate lambda_max(D^{-1}A) by D-inner-product Lanczos.  M = D^{-1}A is
  !! self-adjoint in <u,v>_D = u^T D v, so a niter-step Lanczos builds a Krylov
  !! tridiagonal T whose largest eigenvalue (Ritz value) approximates lambda_max.
  !! Lanczos converges to the extreme eigenvalues far faster than power iteration,
  !! so a short run gives a reliable estimate even on large/ill-separated spectra
  !! (where power iteration underestimates and the Chebyshev smoother turns
  !! indefinite).  Returns the RAW estimate; callers apply the safety factor.
  function hecmw_saamg_lanczos_lambda_max(A, D, niter) result(lambda)
    implicit none
    type(hecmwST_saamg_bcsr),      intent(in) :: A
    type(hecmwST_saamg_blockdiag), intent(in) :: D
    integer(kind=kint),    intent(in) :: niter
    real(kind=kreal) :: lambda
    real(kind=kreal), allocatable :: q(:), qp(:), w(:), aq(:), dv(:), alpha(:), beta(:)
    real(kind=kreal) :: qDq, wDw, bprev
    integer(kind=kint) :: j, i, m
    integer(kind=8) :: seed

    allocate(q(A%n), qp(A%n), w(A%n), aq(A%n), dv(A%n), alpha(niter), beta(niter))

    ! fixed-seed deterministic start (tiny LCG) mapped to (-1,1)
    seed = 88172645463325252_8
    do i = 1, A%n
      seed = seed*6364136223846793005_8 + 1442695040888963407_8
      q(i) = real(seed, kreal)/9.223372036854776d18
    end do
    ! D-normalize the start vector: <q,q>_D = q^T D q
    call hecmw_saamg_blockdiag_matvec(D, q, dv); qDq = dot_product(q, dv)
    if (qDq <= 0.0d0) then; lambda = 0.0d0; deallocate(q,qp,w,aq,dv,alpha,beta); return; end if
    q = q / sqrt(qDq); qp = 0.0d0; bprev = 0.0d0; m = 0

    do j = 1, niter
      call hecmw_saamg_matvec(A, q, aq)             ! aq = A q
      call hecmw_saamg_blockdiag_apply(D, aq, w)    ! w  = D^{-1} A q = M q
      alpha(j) = dot_product(aq, q)                 ! <M q, q>_D = (A q)^T q
      w = w - alpha(j)*q - bprev*qp                 ! orthogonalize (3-term recurrence)
      m = j
      call hecmw_saamg_blockdiag_matvec(D, w, dv); wDw = dot_product(w, dv)
      if (wDw <= 0.0d0) exit                        ! invariant subspace -> exact
      beta(j) = sqrt(wDw)
      if (j == niter) exit
      qp = q; bprev = beta(j); q = w / beta(j)
    end do

    lambda = hecmw_saamg_tridiag_max_eig(alpha, beta, m)
    deallocate(q, qp, w, aq, dv, alpha, beta)
  end function hecmw_saamg_lanczos_lambda_max

  !-----------------------------------------------------------------------------
  ! Chebyshev smoother
  !-----------------------------------------------------------------------------

  !> Set Chebyshev interval [lmax/alpha, lmax] with lmax = 1.1 * lambda_hat.
  !! lambda_hat is the (already safety-scaled) shared estimate.
  subroutine hecmw_saamg_cheb_setup(lambda_hat, deg, alpha, cheb)
    implicit none
    real(kind=kreal),      intent(in)  :: lambda_hat
    integer(kind=kint),    intent(in)  :: deg
    real(kind=kreal),      intent(in)  :: alpha
    type(hecmwST_saamg_chebyshev), intent(out) :: cheb
    real(kind=kreal) :: lmax, lmin
    lmax = 1.1d0 * lambda_hat
    lmin = lmax / alpha
    cheb%deg   = deg
    cheb%theta = 0.5d0*(lmax + lmin)
    cheb%delta = 0.5d0*(lmax - lmin)
    cheb%sigma = cheb%theta / cheb%delta
  end subroutine hecmw_saamg_cheb_setup

  !> Apply the Chebyshev smoother: x <- x + p_k(D^{-1}A)(b - A x).
  !! zero_init = .true. assumes x=0 on entry (skips the initial matvec).
  subroutine hecmw_saamg_cheb_apply_local(A, D, cheb, b, x, zero_init)
    implicit none
    type(hecmwST_saamg_bcsr),      intent(in)    :: A
    type(hecmwST_saamg_blockdiag), intent(in)    :: D
    type(hecmwST_saamg_chebyshev), intent(in)    :: cheb
    real(kind=kreal),      intent(in)    :: b(:)
    real(kind=kreal),      intent(inout) :: x(:)
    logical,               intent(in)    :: zero_init
    real(kind=kreal), allocatable :: r(:), p(:), dr(:), ap(:)
    real(kind=kreal) :: rho, rho2, c1, c2
    integer(kind=kint) :: it, n

    n = A%n
    allocate(r(n), p(n), dr(n), ap(n))

    if (zero_init) then
      x = 0.0d0
      r = b
    else
      call hecmw_saamg_matvec(A, x, ap)
      r = b - ap
    end if

    call hecmw_saamg_blockdiag_apply(D, r, dr)
    p   = dr / cheb%theta
    rho = 1.0d0 / cheb%sigma

    do it = 1, cheb%deg
      x = x + p
      call hecmw_saamg_matvec(A, p, ap)
      r = r - ap
      rho2 = 1.0d0 / (2.0d0*cheb%sigma - rho)
      call hecmw_saamg_blockdiag_apply(D, r, dr)
      c1 = rho2 * rho
      c2 = 2.0d0 * rho2 / cheb%delta
      p  = c1 * p + c2 * dr
      rho = rho2
    end do

    deallocate(r, p, dr, ap)
  end subroutine hecmw_saamg_cheb_apply_local

end module hecmw_precond_SAAMG_smoother
