!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Dense kernels of the multifrontal solver, symmetric (LDLt) mode.
!>
!> A tile or panel is a column major block a(lda,*); only the part on or below the diagonal of
!> its leading square is referenced. Factored columns hold the unit lower triangle of L below
!> the diagonal and D on the diagonal. D is block diagonal: ptype(j) is 1 for a 1x1 pivot, 2 and
!> 3 for the first and second column of a 2x2 pivot whose off-diagonal entry is dsub(j) of the
!> first column (the L entry at that position is 0). The BLAS3 path is used when
!> HECMW_WITH_LAPACK is defined, otherwise plain loops.
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

    do j = 1, n
      do k = 1, j-1
        do i = 1, m
          a(i,j) = a(i,j) - a(i,k)*l(j,k)
        enddo
      enddo
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

end module hecmw_mf_kernel
