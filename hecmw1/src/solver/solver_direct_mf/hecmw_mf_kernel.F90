!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Dense kernels of the multifrontal solver, symmetric (LDLt) and unsymmetric (LU) mode.
!>
!> A tile or panel is a column major block a(lda,*); only the part on or below the diagonal of
!> its leading square is referenced. Factored columns hold the unit lower triangle of L below
!> the diagonal and D on the diagonal. D is block diagonal: ptype(j) is 1 for a 1x1 pivot, 2 and
!> 3 for the first and second column of a 2x2 pivot whose off-diagonal entry is dsub(j) of the
!> first column (the L entry at that position is 0). The BLAS3 path is used when
!> HECMW_WITH_LAPACK is defined, otherwise plain loops.
!>
!> In LU mode the part above the diagonal is held transposed in a second block b of the same
!> shape, b(i,j) = A(j,i) for i > j, so that the row of U of a pivot is a column of b. Factored
!> columns then hold L below the diagonal of a, U on the diagonal of a and in the columns of b.
module hecmw_mf_kernel
  use hecmw_util
  implicit none

  private
  public :: hecmw_mf_kernel_ldlt
  public :: hecmw_mf_kernel_trsm
  public :: hecmw_mf_kernel_panel_nopiv
  public :: hecmw_mf_kernel_panel_piv
  public :: hecmw_mf_kernel_scale
  public :: hecmw_mf_kernel_gemm
  public :: hecmw_mf_kernel_trsv
  public :: hecmw_mf_kernel_trsv_t
  public :: hecmw_mf_kernel_dsolve
  public :: hecmw_mf_kernel_gemv
  public :: hecmw_mf_kernel_gemv_t
  public :: hecmw_mf_kernel_lu
  public :: hecmw_mf_kernel_trsm_rt
  public :: hecmw_mf_kernel_panel_lu_nopiv
  public :: hecmw_mf_kernel_panel_lu_piv
  public :: hecmw_mf_kernel_usolve
  public :: hecmw_mf_kernel_blr_available
  public :: hecmw_mf_kernel_compress
  public :: hecmw_mf_kernel_mult_nn
  public :: hecmw_mf_kernel_mult_tn
  public :: hecmw_mf_kernel_mult_tv
  public :: hecmw_mf_kernel_scale_rows

contains

  !> In place LDLt of the lower triangle of a(n,n) without pivoting; info is the index of the
  !> first pivot whose magnitude is not above zero, 0 on success.
  subroutine hecmw_mf_kernel_ldlt(n, lda, a, zero, info)
    implicit none
    integer(kind=kint), intent(in) :: n, lda
    real(kind=kreal), intent(inout) :: a(lda,*)
    real(kind=kreal), intent(in) :: zero
    integer(kind=kint), intent(out) :: info
    integer(kind=kint) :: i, j, k
    real(kind=kreal) :: d, l

    info = 0
    do k = 1, n
      d = a(k,k)
      if (.not. (abs(d) > zero)) then
        info = k
        return
      endif
      do j = k+1, n
        l = a(j,k) / d
        do i = j, n
          a(i,j) = a(i,j) - a(i,k)*l
        enddo
      enddo
      do i = k+1, n
        a(i,k) = a(i,k) / d
      enddo
    enddo
  end subroutine hecmw_mf_kernel_ldlt

  !> a(m,n) <- a * L^-T * D^-1 with L, D taken from the factored tile l(n,n) (1x1 pivots only).
  subroutine hecmw_mf_kernel_trsm(m, n, ldl, l, lda, a)
    implicit none
    integer(kind=kint), intent(in) :: m, n, ldl, lda
    real(kind=kreal), intent(in) :: l(ldl,*)
    real(kind=kreal), intent(inout) :: a(lda,*)
    integer(kind=kint) :: i, j
#ifdef HECMW_WITH_LAPACK
    external :: dtrsm

    call dtrsm('R', 'L', 'T', 'U', m, n, 1.0d0, l, ldl, a, lda)
    do j = 1, n
      do i = 1, m
        a(i,j) = a(i,j) / l(j,j)
      enddo
    enddo
#else
    integer(kind=kint) :: k

    ! the updates need the unit solve of the earlier columns, so D is divided out afterwards
    do j = 1, n
      do k = 1, j-1
        do i = 1, m
          a(i,j) = a(i,j) - a(i,k)*l(j,k)
        enddo
      enddo
    enddo
    do j = 1, n
      do i = 1, m
        a(i,j) = a(i,j) / l(j,j)
      enddo
    enddo
#endif
  end subroutine hecmw_mf_kernel_trsm

  !> Factor the panel a(m,n) (leading n x n block LDLt, rows below by the panel solve) without
  !> pivoting and check the threshold u afterwards: info is 0 when every pivot magnitude is
  !> above zero and every entry of L is bounded by 1/u, otherwise the first column that fails.
  subroutine hecmw_mf_kernel_panel_nopiv(m, n, lda, a, u, zero, info)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda
    real(kind=kreal), intent(inout) :: a(lda,*)
    real(kind=kreal), intent(in) :: u, zero
    integer(kind=kint), intent(out) :: info
    integer(kind=kint) :: i, j
    real(kind=kreal) :: uinv

    call hecmw_mf_kernel_ldlt(n, lda, a, zero, info)
    if (info /= 0) return
    if (m > n) call hecmw_mf_kernel_trsm(m-n, n, lda, a, lda, a(n+1,1))
    uinv = 1.0d0 / u
    do j = 1, n
      do i = j+1, m
        if (.not. (abs(a(i,j)) <= uinv)) then
          info = j
          return
        endif
      enddo
    enddo
  end subroutine hecmw_mf_kernel_panel_nopiv

  !> Factor the panel a(m,n) with 1x1 and 2x2 pivots chosen among its n columns. A column is
  !> accepted as a 1x1 pivot when it satisfies the threshold u against the whole column; else a
  !> partner is searched, first among the columns of the same block (blk), then among the
  !> others, and the Bunch-Kaufman rule with alpha decides between a 1x1 pivot on the partner
  !> and a 2x2 pivot, the latter accepted when the entries of L stay bounded by 1/u. Columns
  !> without an acceptable pivot are moved to the end of the panel (allow_delay) or stop the
  !> factorization (info = column). On return the first npiv columns are factored, perm(j) is
  !> the original column at position j, and blk is permuted alike. wk holds 2*m words.
  subroutine hecmw_mf_kernel_panel_piv(m, n, lda, a, blk, u, alpha, zero, allow_delay, wk, &
      npiv, perm, ptype, dsub, nswap, n2x2, info)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda
    real(kind=kreal), intent(inout) :: a(lda,*)
    integer(kind=kint), intent(inout) :: blk(n)
    real(kind=kreal), intent(in) :: u, alpha, zero
    logical, intent(in) :: allow_delay
    real(kind=kreal), intent(inout) :: wk(2*m)
    integer(kind=kint), intent(out) :: npiv, perm(n), ptype(n), nswap, n2x2, info
    real(kind=kreal), intent(out) :: dsub(n)
    integer(kind=kint) :: p, nrem, i, r, pass
    real(kind=kreal) :: lam, app, best, uinv
    logical :: accepted, same

    uinv = 1.0d0 / u
    do i = 1, n
      perm(i) = i
      ptype(i) = 0
      dsub(i) = 0.0d0
    enddo
    p = 1
    nrem = n
    nswap = 0
    n2x2 = 0
    info = 0
    do while (p <= nrem)
      lam = 0.0d0
      do i = p+1, m
        lam = max(lam, abs(a(i,p)))
      enddo
      app = abs(a(p,p))
      if (app > zero .and. app >= u*lam) then
        call pivot_1x1()
        cycle
      endif
      accepted = .false.
      if (lam > zero) then
        do pass = 1, 2
          r = 0
          best = zero
          do i = p+1, nrem
            same = (blk(i) == blk(p))
            if (same .neqv. (pass == 1)) cycle
            if (abs(a(i,p)) > best) then
              best = abs(a(i,p))
              r = i
            endif
          enddo
          if (r == 0) cycle
          call try_partner(r, accepted)
          if (accepted) exit
        enddo
      endif
      if (accepted) cycle
      if (.not. allow_delay) then
        info = p
        return
      endif
      if (nrem > p) call swap(p, nrem)
      nrem = nrem - 1
    enddo
    npiv = p - 1

  contains

    !> Bunch-Kaufman choice between a 1x1 pivot on r and the 2x2 pivot (p,r).
    subroutine try_partner(r, accepted)
      integer(kind=kint), intent(in) :: r
      logical, intent(out) :: accepted
      integer(kind=kint) :: i
      real(kind=kreal) :: sig, arr, aa, bb, cc, det, cp, cr, g1, g2

      accepted = .false.
      sig = 0.0d0
      do i = p, r-1
        sig = max(sig, abs(a(r,i)))
      enddo
      do i = r+1, m
        sig = max(sig, abs(a(i,r)))
      enddo
      arr = abs(a(r,r))
      if (arr > zero .and. arr >= alpha*sig) then
        call swap(p, r)
        nswap = nswap + 1
        call pivot_1x1()
        accepted = .true.
        return
      endif
      if (r /= p+1) call swap(p+1, r)
      aa = a(p,p)
      bb = a(p+1,p)
      cc = a(p+1,p+1)
      det = aa*cc - bb*bb
      if (abs(det) > zero*zero) then
        cp = 0.0d0
        cr = 0.0d0
        do i = p+2, m
          cp = max(cp, abs(a(i,p)))
          cr = max(cr, abs(a(i,p+1)))
        enddo
        g1 = (abs(cc)*cp + abs(bb)*cr) / abs(det)
        g2 = (abs(bb)*cp + abs(aa)*cr) / abs(det)
        if (g1 <= uinv .and. g2 <= uinv) then
          if (r /= p+1) nswap = nswap + 1
          call pivot_2x2()
          n2x2 = n2x2 + 1
          accepted = .true.
          return
        endif
      endif
      if (r /= p+1) call swap(p+1, r)
    end subroutine try_partner

    !> Eliminate column p as a 1x1 pivot; the delayed columns beyond nrem are not updated.
    subroutine pivot_1x1()
      integer(kind=kint) :: i, j
      real(kind=kreal) :: d, l

      d = a(p,p)
      do j = p+1, nrem
        l = a(j,p) / d
        do i = j, m
          a(i,j) = a(i,j) - a(i,p)*l
        enddo
      enddo
      do i = p+1, m
        a(i,p) = a(i,p) / d
      enddo
      ptype(p) = 1
      p = p + 1
    end subroutine pivot_1x1

    !> Eliminate columns p, p+1 as a 2x2 pivot E; L = A E^-1 and the update uses A E^-1 A^T.
    subroutine pivot_2x2()
      integer(kind=kint) :: i, j
      real(kind=kreal) :: aa, bb, cc, det

      aa = a(p,p)
      bb = a(p+1,p)
      cc = a(p+1,p+1)
      det = aa*cc - bb*bb
      do i = p+2, m
        wk(i) = (cc*a(i,p) - bb*a(i,p+1)) / det
        wk(m+i) = (aa*a(i,p+1) - bb*a(i,p)) / det
      enddo
      do j = p+2, nrem
        do i = j, m
          a(i,j) = a(i,j) - wk(i)*a(j,p) - wk(m+i)*a(j,p+1)
        enddo
      enddo
      do i = p+2, m
        a(i,p) = wk(i)
        a(i,p+1) = wk(m+i)
      enddo
      a(p+1,p) = 0.0d0
      ptype(p) = 2
      ptype(p+1) = 3
      dsub(p) = bb
      p = p + 2
    end subroutine pivot_2x2

    !> Symmetric exchange of panel positions x < y (both not yet eliminated).
    subroutine swap(x, y)
      integer(kind=kint), intent(in) :: x, y
      integer(kind=kint) :: i, t
      real(kind=kreal) :: v

      v = a(x,x)
      a(x,x) = a(y,y)
      a(y,y) = v
      do i = 1, x-1
        v = a(x,i)
        a(x,i) = a(y,i)
        a(y,i) = v
      enddo
      do i = x+1, y-1
        v = a(i,x)
        a(i,x) = a(y,i)
        a(y,i) = v
      enddo
      do i = y+1, m
        v = a(i,x)
        a(i,x) = a(i,y)
        a(i,y) = v
      enddo
      t = perm(x)
      perm(x) = perm(y)
      perm(y) = t
      t = blk(x)
      blk(x) = blk(y)
      blk(y) = t
    end subroutine swap

  end subroutine hecmw_mf_kernel_panel_piv

  !> w(m,n) <- a * D with the block diagonal D of the factored tile l(n,n).
  subroutine hecmw_mf_kernel_scale(m, n, ldl, l, ptype, dsub, lda, a, w)
    implicit none
    integer(kind=kint), intent(in) :: m, n, ldl, lda
    real(kind=kreal), intent(in) :: l(ldl,*)
    integer(kind=kint), intent(in) :: ptype(n)
    real(kind=kreal), intent(in) :: dsub(n)
    real(kind=kreal), intent(in) :: a(lda,*)
    real(kind=kreal), intent(out) :: w(m,n)
    integer(kind=kint) :: i, j
    real(kind=kreal) :: d11, d21, d22

    do j = 1, n
      if (ptype(j) == 1) then
        do i = 1, m
          w(i,j) = a(i,j) * l(j,j)
        enddo
      else if (ptype(j) == 2) then
        d11 = l(j,j)
        d21 = dsub(j)
        d22 = l(j+1,j+1)
        do i = 1, m
          w(i,j) = a(i,j)*d11 + a(i,j+1)*d21
          w(i,j+1) = a(i,j)*d21 + a(i,j+1)*d22
        enddo
      endif
    enddo
  end subroutine hecmw_mf_kernel_scale

  !> c(m,n) <- c - w(m,k) * b(n,k)^T
  subroutine hecmw_mf_kernel_gemm(m, n, k, ldw, w, ldb, b, ldc, c)
    implicit none
    integer(kind=kint), intent(in) :: m, n, k, ldw, ldb, ldc
    real(kind=kreal), intent(in) :: w(ldw,*)
    real(kind=kreal), intent(in) :: b(ldb,*)
    real(kind=kreal), intent(inout) :: c(ldc,*)
#ifdef HECMW_WITH_LAPACK
    external :: dgemm

    call dgemm('N', 'T', m, n, k, -1.0d0, w, ldw, b, ldb, 1.0d0, c, ldc)
#else
    integer(kind=kint) :: i, j, p

    do j = 1, n
      do p = 1, k
        do i = 1, m
          c(i,j) = c(i,j) - w(i,p)*b(j,p)
        enddo
      enddo
    enddo
#endif
  end subroutine hecmw_mf_kernel_gemm

  !> x <- L^-1 x with the unit lower triangle of the tile l(n,n).
  subroutine hecmw_mf_kernel_trsv(n, ldl, l, x)
    implicit none
    integer(kind=kint), intent(in) :: n, ldl
    real(kind=kreal), intent(in) :: l(ldl,*)
    real(kind=kreal), intent(inout) :: x(n)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = j+1, n
        x(i) = x(i) - l(i,j)*x(j)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_trsv

  !> x <- L^-T x with the unit lower triangle of the tile l(n,n).
  subroutine hecmw_mf_kernel_trsv_t(n, ldl, l, x)
    implicit none
    integer(kind=kint), intent(in) :: n, ldl
    real(kind=kreal), intent(in) :: l(ldl,*)
    real(kind=kreal), intent(inout) :: x(n)
    integer(kind=kint) :: i, j

    do j = n, 1, -1
      do i = j+1, n
        x(j) = x(j) - l(i,j)*x(i)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_trsv_t

  !> x <- D^-1 x with the block diagonal D of the factored tile l(n,n).
  subroutine hecmw_mf_kernel_dsolve(n, ldl, l, ptype, dsub, x)
    implicit none
    integer(kind=kint), intent(in) :: n, ldl
    real(kind=kreal), intent(in) :: l(ldl,*)
    integer(kind=kint), intent(in) :: ptype(n)
    real(kind=kreal), intent(in) :: dsub(n)
    real(kind=kreal), intent(inout) :: x(n)
    integer(kind=kint) :: j
    real(kind=kreal) :: d11, d21, d22, det, x1, x2

    do j = 1, n
      if (ptype(j) == 1) then
        x(j) = x(j) / l(j,j)
      else if (ptype(j) == 2) then
        d11 = l(j,j)
        d21 = dsub(j)
        d22 = l(j+1,j+1)
        det = d11*d22 - d21*d21
        x1 = x(j)
        x2 = x(j+1)
        x(j) = (d22*x1 - d21*x2) / det
        x(j+1) = (d11*x2 - d21*x1) / det
      endif
    enddo
  end subroutine hecmw_mf_kernel_dsolve

  !> y(m) <- y - a(m,n) * x(n)
  subroutine hecmw_mf_kernel_gemv(m, n, lda, a, x, y)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda
    real(kind=kreal), intent(in) :: a(lda,*)
    real(kind=kreal), intent(in) :: x(n)
    real(kind=kreal), intent(inout) :: y(m)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = 1, m
        y(i) = y(i) - a(i,j)*x(j)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_gemv

  !> y(n) <- y - a(m,n)^T * x(m)
  subroutine hecmw_mf_kernel_gemv_t(m, n, lda, a, x, y)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda
    real(kind=kreal), intent(in) :: a(lda,*)
    real(kind=kreal), intent(in) :: x(m)
    real(kind=kreal), intent(inout) :: y(n)
    integer(kind=kint) :: i, j

    do j = 1, n
      do i = 1, m
        y(j) = y(j) - a(i,j)*x(i)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_gemv_t

  !> In place LU of the n x n block held in a (diagonal and below) and b (above, transposed)
  !> without pivoting; info is the index of the first pivot whose magnitude is not above zero,
  !> 0 on success.
  subroutine hecmw_mf_kernel_lu(n, lda, a, ldb, b, zero, info)
    implicit none
    integer(kind=kint), intent(in) :: n, lda, ldb
    real(kind=kreal), intent(inout) :: a(lda,*), b(ldb,*)
    real(kind=kreal), intent(in) :: zero
    integer(kind=kint), intent(out) :: info
    integer(kind=kint) :: i, j, k
    real(kind=kreal) :: d, l, u

    info = 0
    do k = 1, n
      d = a(k,k)
      if (.not. (abs(d) > zero)) then
        info = k
        return
      endif
      do i = k+1, n
        a(i,k) = a(i,k) / d
      enddo
      do j = k+1, n
        u = b(j,k)
        do i = j, n
          a(i,j) = a(i,j) - a(i,k)*u
        enddo
        l = a(j,k)
        do i = j+1, n
          b(i,j) = b(i,j) - b(i,k)*l
        enddo
      enddo
    enddo
  end subroutine hecmw_mf_kernel_lu

  !> a(m,n) <- a * T^-T with the lower triangle of t(n,n), unit diagonal when unit.
  subroutine hecmw_mf_kernel_trsm_rt(m, n, ldt, t, unit, lda, a)
    implicit none
    integer(kind=kint), intent(in) :: m, n, ldt, lda
    real(kind=kreal), intent(in) :: t(ldt,*)
    logical, intent(in) :: unit
    real(kind=kreal), intent(inout) :: a(lda,*)
#ifdef HECMW_WITH_LAPACK
    external :: dtrsm

    if (unit) then
      call dtrsm('R', 'L', 'T', 'U', m, n, 1.0d0, t, ldt, a, lda)
    else
      call dtrsm('R', 'L', 'T', 'N', m, n, 1.0d0, t, ldt, a, lda)
    endif
#else
    integer(kind=kint) :: i, j, k

    do j = 1, n
      do k = 1, j-1
        do i = 1, m
          a(i,j) = a(i,j) - a(i,k)*t(j,k)
        enddo
      enddo
      if (.not. unit) then
        do i = 1, m
          a(i,j) = a(i,j) / t(j,j)
        enddo
      endif
    enddo
#endif
  end subroutine hecmw_mf_kernel_trsm_rt

  !> LU of the panel a(m,n), b(m,n) (leading n x n block, rows below by the panel solves) without
  !> pivoting and the threshold check afterwards, as panel_nopiv. The diagonal of b receives
  !> the diagonal of U.
  subroutine hecmw_mf_kernel_panel_lu_nopiv(m, n, lda, a, ldb, b, u, zero, info)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda, ldb
    real(kind=kreal), intent(inout) :: a(lda,*), b(ldb,*)
    real(kind=kreal), intent(in) :: u, zero
    integer(kind=kint), intent(out) :: info
    integer(kind=kint) :: i, j
    real(kind=kreal) :: uinv

    call hecmw_mf_kernel_lu(n, lda, a, ldb, b, zero, info)
    if (info /= 0) return
    if (m > n) then
      do j = 1, n
        b(j,j) = a(j,j)
      enddo
      call hecmw_mf_kernel_trsm_rt(m-n, n, ldb, b, .false., lda, a(n+1,1))
      call hecmw_mf_kernel_trsm_rt(m-n, n, lda, a, .true., ldb, b(n+1,1))
    endif
    uinv = 1.0d0 / u
    do j = 1, n
      do i = j+1, m
        if (.not. (abs(a(i,j)) <= uinv)) then
          info = j
          return
        endif
      enddo
    enddo
  end subroutine hecmw_mf_kernel_panel_lu_nopiv

  !> LU of the panel a(m,n), b(m,n) with threshold partial pivoting among its n rows: column p
  !> takes the row of largest magnitude among the rows of the same block (blkr against blk(p)),
  !> else among the other panel rows, when it reaches u times the largest magnitude of the whole
  !> column. Rows are exchanged, columns are not. Columns without an acceptable pivot are moved
  !> with their rows to the end of the panel (allow_delay) or stop the factorization (info =
  !> column). On return the first npiv columns are factored, permc(j) and permr(j) are the
  !> original column and row at position j, and blk, blkr are permuted alike.
  subroutine hecmw_mf_kernel_panel_lu_piv(m, n, lda, a, ldb, b, blk, blkr, u, zero, allow_delay, &
      npiv, permc, permr, nswap, info)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda, ldb
    real(kind=kreal), intent(inout) :: a(lda,*), b(ldb,*)
    integer(kind=kint), intent(inout) :: blk(n), blkr(n)
    real(kind=kreal), intent(in) :: u, zero
    logical, intent(in) :: allow_delay
    integer(kind=kint), intent(out) :: npiv, permc(n), permr(n), nswap, info
    integer(kind=kint) :: p, nrem, i, r, pass
    real(kind=kreal) :: lam, best
    logical :: same

    do i = 1, n
      permc(i) = i
      permr(i) = i
    enddo
    p = 1
    nrem = n
    nswap = 0
    info = 0
    do while (p <= nrem)
      lam = 0.0d0
      do i = p, m
        lam = max(lam, abs(a(i,p)))
      enddo
      r = 0
      if (lam > zero) then
        do pass = 1, 2
          best = 0.0d0
          do i = p, nrem
            same = (blkr(i) == blk(p))
            if (same .neqv. (pass == 1)) cycle
            if (abs(a(i,p)) > best) then
              best = abs(a(i,p))
              r = i
            endif
          enddo
          if (r /= 0) then
            if (best >= u*lam) exit
            r = 0
          endif
        enddo
      endif
      if (r == 0) then
        if (.not. allow_delay) then
          info = p
          return
        endif
        if (nrem > p) call swap_sym(p, nrem)
        nrem = nrem - 1
        cycle
      endif
      if (r /= p) then
        call swap_row(p, r)
        nswap = nswap + 1
      endif
      call pivot()
    enddo
    npiv = p - 1

  contains

    !> Eliminate column p; the delayed columns and rows beyond nrem are not updated.
    subroutine pivot()
      integer(kind=kint) :: i, j
      real(kind=kreal) :: d, l, uu

      d = a(p,p)
      do i = p+1, m
        a(i,p) = a(i,p) / d
      enddo
      do j = p+1, nrem
        uu = b(j,p)
        do i = j, m
          a(i,j) = a(i,j) - a(i,p)*uu
        enddo
        l = a(j,p)
        do i = j+1, m
          b(i,j) = b(i,j) - b(i,p)*l
        enddo
      enddo
      p = p + 1
    end subroutine pivot

    !> Exchange of the rows x < y of the panel (both not yet eliminated).
    subroutine swap_row(x, y)
      integer(kind=kint), intent(in) :: x, y
      integer(kind=kint) :: c, t
      real(kind=kreal) :: v

      do c = 1, x
        v = a(x,c)
        a(x,c) = a(y,c)
        a(y,c) = v
      enddo
      do c = x+1, y-1
        v = b(c,x)
        b(c,x) = a(y,c)
        a(y,c) = v
      enddo
      v = b(y,x)
      b(y,x) = a(y,y)
      a(y,y) = v
      do c = y+1, m
        v = b(c,x)
        b(c,x) = b(c,y)
        b(c,y) = v
      enddo
      t = permr(x)
      permr(x) = permr(y)
      permr(y) = t
      t = blkr(x)
      blkr(x) = blkr(y)
      blkr(y) = t
    end subroutine swap_row

    !> Symmetric exchange of panel positions x < y (both not yet eliminated).
    subroutine swap_sym(x, y)
      integer(kind=kint), intent(in) :: x, y
      integer(kind=kint) :: i, t
      real(kind=kreal) :: v

      v = a(x,x)
      a(x,x) = a(y,y)
      a(y,y) = v
      do i = 1, x-1
        v = a(x,i)
        a(x,i) = a(y,i)
        a(y,i) = v
        v = b(x,i)
        b(x,i) = b(y,i)
        b(y,i) = v
      enddo
      do i = x+1, y-1
        v = a(i,x)
        a(i,x) = b(y,i)
        b(y,i) = v
        v = a(y,i)
        a(y,i) = b(i,x)
        b(i,x) = v
      enddo
      v = a(y,x)
      a(y,x) = b(y,x)
      b(y,x) = v
      do i = y+1, m
        v = a(i,x)
        a(i,x) = a(i,y)
        a(i,y) = v
        v = b(i,x)
        b(i,x) = b(i,y)
        b(i,y) = v
      enddo
      t = permc(x)
      permc(x) = permc(y)
      permc(y) = t
      t = permr(x)
      permr(x) = permr(y)
      permr(y) = t
      t = blk(x)
      blk(x) = blk(y)
      blk(y) = t
      t = blkr(x)
      blkr(x) = blkr(y)
      blkr(y) = t
    end subroutine swap_sym

  end subroutine hecmw_mf_kernel_panel_lu_piv

  !> x <- U^-1 x with the diagonal of U on the diagonal of l(n,n) and its strictly upper part
  !> transposed in the strictly lower part of ut(n,n).
  subroutine hecmw_mf_kernel_usolve(n, ldl, l, ldu, ut, x)
    implicit none
    integer(kind=kint), intent(in) :: n, ldl, ldu
    real(kind=kreal), intent(in) :: l(ldl,*), ut(ldu,*)
    real(kind=kreal), intent(inout) :: x(n)
    integer(kind=kint) :: i, j

    do j = n, 1, -1
      do i = j+1, n
        x(j) = x(j) - ut(i,j)*x(i)
      enddo
      x(j) = x(j) / l(j,j)
    enddo
  end subroutine hecmw_mf_kernel_usolve

  !> Whether the BLR compression is available (LAPACK build).
  function hecmw_mf_kernel_blr_available() result(avail)
    implicit none
    logical :: avail

#ifdef HECMW_WITH_LAPACK
    avail = .true.
#else
    avail = .false.
#endif
  end function hecmw_mf_kernel_blr_available

  !> Truncated rank revealing QR of a(m,n) with the truncation |R(k,k)| <= eps |R(1,1)|:
  !> on success the leading rank columns of a hold U and v(n,rank) holds V with a ~ U V^T.
  !> rank is -1, with a destroyed, when the low rank form would not take fewer than m*n words
  !> (the caller keeps its full rank copy) or without LAPACK.
  !> The column pivoted Householder QR stops right at the truncation or at the word bound, so
  !> the cost is O(m*n*rank) where dgeqp3 pays the full decomposition. The pivot column is
  !> selected by downdated norms with the dgeqp3 recomputation guard, but the truncation is
  !> decided on the freshly computed norm so a stale estimate cannot end the sweep early.
  !> rcap tightens the word bound: the sweep gives up (rank = -1) beyond min(rlim, rcap).
  !> The pivoted sweep with the downdated column norms is adapted from LAPACK's dlaqp2 (with
  !> the norm recomputation guard of LAPACK Working Note 176) and the explicit formation of
  !> the Q columns follows dorg2r. LAPACK is distributed under the modified BSD license:
  !>
  !> Copyright (c) 1992-2025 The University of Tennessee and The University of Tennessee
  !>                         Research Foundation. All rights reserved.
  !> Copyright (c) 2000-2025 The University of California Berkeley. All rights reserved.
  !> Copyright (c) 2006-2025 The University of Colorado Denver. All rights reserved.
  !>
  !> Redistribution and use in source and binary forms, with or without modification, are
  !> permitted provided that the following conditions are met:
  !> - Redistributions of source code must retain the above copyright notice, this list of
  !>   conditions and the following disclaimer.
  !> - Redistributions in binary form must reproduce the above copyright notice, this list of
  !>   conditions and the following disclaimer listed in this license in the documentation
  !>   and/or other materials provided with the distribution.
  !> - Neither the name of the copyright holders nor the names of its contributors may be used
  !>   to endorse or promote products derived from this software without specific prior
  !>   written permission.
  !> The copyright holders provide no reassurances that the source code provided does not
  !> infringe any patent, copyright, or any other intellectual property rights of third
  !> parties. The copyright holders disclaim any liability to any recipient for claims brought
  !> against recipient by any third party for infringement of that parties intellectual
  !> property rights.
  !> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
  !> EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  !> MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
  !> THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  !> SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
  !> OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  !> INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  !> LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  !> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  subroutine hecmw_mf_kernel_compress(m, n, lda, a, eps, ldv, v, rank, rcap)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda, ldv
    real(kind=kreal), intent(inout) :: a(lda,*)
    real(kind=kreal), intent(in) :: eps
    real(kind=kreal), intent(out) :: v(ldv,*)
    integer(kind=kint), intent(out) :: rank
    integer(kind=kint), intent(in), optional :: rcap
#ifdef HECMW_WITH_LAPACK
    integer(kind=kint), allocatable :: jpvt(:)
    real(kind=kreal), allocatable :: tau(:), vn1(:), vn2(:), wk(:)
    integer(kind=kint) :: i, j, k, p, r, rlim, itmp
    real(kind=kreal) :: r11, rkk, akk, tol3z, tmp, tmp2, dtmp
    real(kind=kreal), external :: dnrm2, dlamch
    external :: dlarfg, dgemv, dger

    rank = -1
    ! largest rank whose low rank form takes fewer than m*n words; always below min(m,n)
    rlim = int((int(m, 8)*n - 1)/(int(m, 8) + n), kind=kint)
    if (present(rcap)) rlim = min(rlim, max(rcap, 0))
    tol3z = sqrt(dlamch('Epsilon'))
    allocate(jpvt(n), tau(rlim + 1), vn1(n), vn2(n), wk(n))
    do j = 1, n
      jpvt(j) = j
      vn1(j) = dnrm2(m, a(1,j), 1)
      vn2(j) = vn1(j)
    enddo
    r11 = 0.0d0
    r = 0
    do k = 1, rlim + 1
      p = k
      do j = k + 1, n
        if (vn1(j) > vn1(p)) p = j
      enddo
      if (p /= k) then
        do i = 1, m
          dtmp = a(i,k)
          a(i,k) = a(i,p)
          a(i,p) = dtmp
        enddo
        itmp = jpvt(k)
        jpvt(k) = jpvt(p)
        jpvt(p) = itmp
        dtmp = vn1(k)
        vn1(k) = vn1(p)
        vn1(p) = dtmp
        dtmp = vn2(k)
        vn2(k) = vn2(p)
        vn2(p) = dtmp
      endif
      rkk = dnrm2(m - k + 1, a(k,k), 1)
      if (k == 1) r11 = rkk
      if (.not. (rkk > eps*r11)) then
        r = k - 1
        exit
      endif
      if (k > rlim) then
        deallocate(jpvt, tau, vn1, vn2, wk)
        return
      endif
      call dlarfg(m - k + 1, a(k,k), a(k+1,k), 1, tau(k))
      akk = a(k,k)
      a(k,k) = 1.0d0
      call dgemv('T', m - k + 1, n - k, 1.0d0, a(k,k+1), lda, a(k,k), 1, 0.0d0, wk, 1)
      call dger(m - k + 1, n - k, -tau(k), a(k,k), 1, wk, 1, a(k,k+1), lda)
      a(k,k) = akk
      do j = k + 1, n
        if (vn1(j) /= 0.0d0) then
          tmp = max(1.0d0 - (abs(a(k,j))/vn1(j))**2, 0.0d0)
          tmp2 = tmp*(vn1(j)/vn2(j))**2
          if (tmp2 <= tol3z) then
            vn1(j) = dnrm2(m - k, a(k+1,j), 1)
            vn2(j) = vn1(j)
          else
            vn1(j) = vn1(j)*sqrt(tmp)
          endif
        endif
      enddo
    enddo
    do j = 1, n
      do k = 1, r
        if (k <= j) then
          v(jpvt(j), k) = a(k, j)
        else
          v(jpvt(j), k) = 0.0d0
        endif
      enddo
    enddo
    do k = r, 1, -1
      if (k < r) then
        a(k,k) = 1.0d0
        call dgemv('T', m - k + 1, r - k, 1.0d0, a(k,k+1), lda, a(k,k), 1, 0.0d0, wk, 1)
        call dger(m - k + 1, r - k, -tau(k), a(k,k), 1, wk, 1, a(k,k+1), lda)
      endif
      do i = k + 1, m
        a(i,k) = -tau(k)*a(i,k)
      enddo
      a(k,k) = 1.0d0 - tau(k)
      do i = 1, k - 1
        a(i,k) = 0.0d0
      enddo
    enddo
    rank = r
    deallocate(jpvt, tau, vn1, vn2, wk)
#else
    rank = -1
#endif
  end subroutine hecmw_mf_kernel_compress

  !> c(m,n) <- a(m,k) * b(k,n)
  subroutine hecmw_mf_kernel_mult_nn(m, n, k, lda, a, ldb, b, ldc, c)
    implicit none
    integer(kind=kint), intent(in) :: m, n, k, lda, ldb, ldc
    real(kind=kreal), intent(in) :: a(lda,*), b(ldb,*)
    real(kind=kreal), intent(out) :: c(ldc,*)
#ifdef HECMW_WITH_LAPACK
    external :: dgemm

    call dgemm('N', 'N', m, n, k, 1.0d0, a, lda, b, ldb, 0.0d0, c, ldc)
#else
    integer(kind=kint) :: i, j, p

    do j = 1, n
      do i = 1, m
        c(i,j) = 0.0d0
      enddo
      do p = 1, k
        do i = 1, m
          c(i,j) = c(i,j) + a(i,p)*b(p,j)
        enddo
      enddo
    enddo
#endif
  end subroutine hecmw_mf_kernel_mult_nn

  !> c(m,n) <- a(k,m)^T * b(k,n)
  subroutine hecmw_mf_kernel_mult_tn(m, n, k, lda, a, ldb, b, ldc, c)
    implicit none
    integer(kind=kint), intent(in) :: m, n, k, lda, ldb, ldc
    real(kind=kreal), intent(in) :: a(lda,*), b(ldb,*)
    real(kind=kreal), intent(out) :: c(ldc,*)
#ifdef HECMW_WITH_LAPACK
    external :: dgemm

    call dgemm('T', 'N', m, n, k, 1.0d0, a, lda, b, ldb, 0.0d0, c, ldc)
#else
    integer(kind=kint) :: i, j, p

    do j = 1, n
      do i = 1, m
        c(i,j) = 0.0d0
        do p = 1, k
          c(i,j) = c(i,j) + a(p,i)*b(p,j)
        enddo
      enddo
    enddo
#endif
  end subroutine hecmw_mf_kernel_mult_tn

  !> y(n) <- a(m,n)^T * x(m)
  subroutine hecmw_mf_kernel_mult_tv(m, n, lda, a, x, y)
    implicit none
    integer(kind=kint), intent(in) :: m, n, lda
    real(kind=kreal), intent(in) :: a(lda,*)
    real(kind=kreal), intent(in) :: x(m)
    real(kind=kreal), intent(out) :: y(n)
    integer(kind=kint) :: i, j

    do j = 1, n
      y(j) = 0.0d0
      do i = 1, m
        y(j) = y(j) + a(i,j)*x(i)
      enddo
    enddo
  end subroutine hecmw_mf_kernel_mult_tv

  !> w(n,r) <- D * v(n,r) with the block diagonal D of the factored tile l(n,n).
  subroutine hecmw_mf_kernel_scale_rows(n, r, ldl, l, ptype, dsub, ldv, v, ldw, w)
    implicit none
    integer(kind=kint), intent(in) :: n, r, ldl, ldv, ldw
    real(kind=kreal), intent(in) :: l(ldl,*)
    integer(kind=kint), intent(in) :: ptype(n)
    real(kind=kreal), intent(in) :: dsub(n)
    real(kind=kreal), intent(in) :: v(ldv,*)
    real(kind=kreal), intent(out) :: w(ldw,*)
    integer(kind=kint) :: q, c
    real(kind=kreal) :: d11, d21, d22

    q = 1
    do while (q <= n)
      if (ptype(q) == 1) then
        do c = 1, r
          w(q,c) = v(q,c) * l(q,q)
        enddo
        q = q + 1
      else
        d11 = l(q,q)
        d21 = dsub(q)
        d22 = l(q+1,q+1)
        do c = 1, r
          w(q,c) = v(q,c)*d11 + v(q+1,c)*d21
          w(q+1,c) = v(q,c)*d21 + v(q+1,c)*d22
        enddo
        q = q + 2
      endif
    enddo
  end subroutine hecmw_mf_kernel_scale_rows

end module hecmw_mf_kernel
