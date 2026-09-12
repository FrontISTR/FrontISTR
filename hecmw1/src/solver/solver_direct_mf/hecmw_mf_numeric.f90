!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Numeric factorization (tiled multifrontal LDLt without pivoting) and the triangular
!>        solves, on the read-only structure of hecmw_mf_symbolic.
!>
!> A front of nrow DOFs is partitioned into tiles at node boundaries; the ncol pivot DOFs
!> (the own columns of the supernode) end on a tile boundary. Only the tiles on or below the
!> diagonal are stored, tile column by tile column, each tile column major. The tile columns of
!> the pivot part form the factor panel [L11; L21] kept in lval, the remaining ones form the
!> contribution block that is pushed on the stack for the parent.
module hecmw_mf_numeric
  use hecmw_util
  use hecmw_mf_symbolic
  use hecmw_mf_kernel
  implicit none

  private
  public :: hecmwST_mf_factor
  public :: hecmw_mf_numeric_init
  public :: hecmw_mf_numeric_factor
  public :: hecmw_mf_numeric_solve
  public :: hecmw_mf_numeric_print
  public :: hecmw_mf_numeric_finalize

  type hecmwST_mf_factor
    integer(kind=kint) :: tile = 0
    integer(kind=kint) :: nnode = 0
    integer(kind=kint) :: nsuper = 0
    integer(kind=kint) :: ndof_tot = 0
    integer(kind=kint), allocatable :: cdofptr(:)   !< column k owns the permuted DOFs cdofptr(k):cdofptr(k+1)-1
    integer(kind=kint), allocatable :: pdof(:)      !< pdof(c): original DOF of permuted DOF c
    integer(kind=kint), allocatable :: tptr(:)      !< tile boundaries of supernode s are the DOF offsets
    integer(kind=kint), allocatable :: tbnd(:)      !< tbnd(tptr(s):tptr(s+1)-1), from 0 to nrow
    integer(kind=kint), allocatable :: ntc(:)       !< tiles covering the pivot columns of supernode s
    integer(kind=8), allocatable :: lptr(:)         !< panel of supernode s is lval(lptr(s):lptr(s+1)-1)
    integer(kind=8), allocatable :: cbsize(:)       !< words of the contribution block of supernode s
    integer(kind=kint), allocatable :: chead(:)     !< children of s: chead(s), cnext(...), in decreasing
    integer(kind=kint), allocatable :: cnext(:)     !< order, which is the order they are popped
    real(kind=kreal), allocatable :: lval(:)
    real(kind=kreal), allocatable :: sval(:)        !< stack of contribution blocks
    real(kind=kreal), allocatable :: fval(:)        !< work space of the front being assembled
    real(kind=kreal), allocatable :: wval(:)        !< work tile of the scaled panel L*D
    integer(kind=8) :: factor_words = 0
    integer(kind=8) :: factor_nnz = 0               !< panel entries on or below the diagonal
    integer(kind=8) :: stack_peak = 0
    integer(kind=8) :: front_words = 0
    integer(kind=kint) :: max_tiles = 0
    integer(kind=kint) :: max_rows = 0
    integer(kind=kint) :: max_tile_dim = 0
    logical :: factored = .false.
  end type hecmwST_mf_factor

contains

  !> Tile partitions, factor and stack sizes, and the work space; tile is the target tile size in
  !> DOFs (a tile ends at the last node boundary that keeps it within the target).
  subroutine hecmw_mf_numeric_init(sym, tile, fct)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    integer(kind=kint), intent(in) :: tile
    type(hecmwST_mf_factor), intent(out) :: fct
    integer(kind=8), allocatable :: coloff(:)
    integer(kind=kint), allocatable :: odofptr(:)
    integer(kind=kint) :: n, ns, s, k, i, d, c, nown, nrow_nodes, nt, ntc, nrow, ncol, pass, cur, acc, nd
    integer(kind=8) :: top, hj

    n = sym%nnode
    ns = sym%nsuper
    fct%tile = tile
    fct%nnode = n
    fct%nsuper = ns

    allocate(fct%cdofptr(n+1), odofptr(n+1))
    fct%cdofptr(1) = 1
    odofptr(1) = 1
    do k = 1, n
      fct%cdofptr(k+1) = fct%cdofptr(k) + sym%ndof(k)
      odofptr(k+1) = odofptr(k) + sym%ndof(sym%invp(k))
    enddo
    fct%ndof_tot = fct%cdofptr(n+1) - 1
    allocate(fct%pdof(fct%ndof_tot))
    do k = 1, n
      do d = 1, sym%ndof(k)
        fct%pdof(fct%cdofptr(k)+d-1) = odofptr(sym%perm(k)) + d - 1
      enddo
    enddo
    deallocate(odofptr)

    allocate(fct%tptr(ns+1), fct%ntc(ns))
    fct%max_rows = 0
    fct%max_tile_dim = 0
    do pass = 1, 2
      fct%tptr(1) = 1
      do s = 1, ns
        nown = sym%sptr(s+1) - sym%sptr(s)
        nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
        fct%max_rows = max(fct%max_rows, nrow_nodes)
        nt = 0
        cur = 0
        acc = 0
        ntc = 0
        if (pass == 2) fct%tbnd(fct%tptr(s)) = 0
        do i = 1, nrow_nodes
          nd = sym%ndof(sym%rlist(sym%rptr(s)+i-1))
          if (acc > cur .and. (acc + nd - cur > tile .or. i == nown + 1)) then
            nt = nt + 1
            if (pass == 2) fct%tbnd(fct%tptr(s)+nt) = acc
            fct%max_tile_dim = max(fct%max_tile_dim, acc - cur)
            cur = acc
          endif
          if (i == nown + 1) ntc = nt
          acc = acc + nd
        enddo
        nt = nt + 1
        if (pass == 2) fct%tbnd(fct%tptr(s)+nt) = acc
        fct%max_tile_dim = max(fct%max_tile_dim, acc - cur)
        if (nown == nrow_nodes) ntc = nt
        fct%ntc(s) = ntc
        fct%tptr(s+1) = fct%tptr(s) + nt + 1
        fct%max_tiles = max(fct%max_tiles, nt)
      enddo
      if (pass == 1) allocate(fct%tbnd(fct%tptr(ns+1)-1))
    enddo

    allocate(fct%lptr(ns+1), fct%cbsize(ns), coloff(fct%max_tiles+1))
    fct%lptr(1) = 1
    fct%front_words = 0
    fct%factor_nnz = 0
    do s = 1, ns
      call mf_layout(fct, s, nt, ntc, nrow, ncol, coloff)
      fct%lptr(s+1) = fct%lptr(s) + coloff(ntc+1)
      fct%cbsize(s) = coloff(nt+1) - coloff(ntc+1)
      fct%front_words = max(fct%front_words, coloff(nt+1))
      fct%factor_nnz = fct%factor_nnz + coloff(ntc+1)
      do k = 1, ntc
        hj = fct%tbnd(fct%tptr(s)+k) - fct%tbnd(fct%tptr(s)+k-1)
        fct%factor_nnz = fct%factor_nnz - hj*(hj-1)/2
      enddo
    enddo
    fct%factor_words = fct%lptr(ns+1) - 1
    deallocate(coloff)

    allocate(fct%chead(ns), fct%cnext(ns))
    fct%chead(1:ns) = 0
    do s = 1, ns
      c = sym%sparent(s)
      if (c == 0) cycle
      fct%cnext(s) = fct%chead(c)
      fct%chead(c) = s
    enddo
    top = 0
    fct%stack_peak = 0
    do s = 1, ns
      fct%stack_peak = max(fct%stack_peak, top)
      c = fct%chead(s)
      do while (c /= 0)
        top = top - fct%cbsize(c)
        c = fct%cnext(c)
      enddo
      top = top + fct%cbsize(s)
      fct%stack_peak = max(fct%stack_peak, top)
    enddo

    allocate(fct%lval(fct%factor_words))
    allocate(fct%sval(max(fct%stack_peak, 1_8)))
    allocate(fct%fval(fct%front_words))
    allocate(fct%wval(int(fct%max_tile_dim, 8)**2))
    fct%factored = .false.
  end subroutine hecmw_mf_numeric_init

  subroutine hecmw_mf_numeric_finalize(fct)
    implicit none
    type(hecmwST_mf_factor), intent(inout) :: fct

    fct%nsuper = 0
    fct%factored = .false.
    if (allocated(fct%cdofptr)) deallocate(fct%cdofptr)
    if (allocated(fct%pdof)) deallocate(fct%pdof)
    if (allocated(fct%tptr)) deallocate(fct%tptr)
    if (allocated(fct%tbnd)) deallocate(fct%tbnd)
    if (allocated(fct%ntc)) deallocate(fct%ntc)
    if (allocated(fct%lptr)) deallocate(fct%lptr)
    if (allocated(fct%cbsize)) deallocate(fct%cbsize)
    if (allocated(fct%chead)) deallocate(fct%chead)
    if (allocated(fct%cnext)) deallocate(fct%cnext)
    if (allocated(fct%lval)) deallocate(fct%lval)
    if (allocated(fct%sval)) deallocate(fct%sval)
    if (allocated(fct%fval)) deallocate(fct%fval)
    if (allocated(fct%wval)) deallocate(fct%wval)
  end subroutine hecmw_mf_numeric_finalize

  !> Tile grid of supernode s: nt tiles of which ntc cover the ncol pivot DOFs; coloff(j) is the
  !> word offset of tile column j, coloff(nt+1) the size of the whole lower tile grid.
  subroutine mf_layout(fct, s, nt, ntc, nrow, ncol, coloff)
    implicit none
    type(hecmwST_mf_factor), intent(in) :: fct
    integer(kind=kint), intent(in) :: s
    integer(kind=kint), intent(out) :: nt, ntc, nrow, ncol
    integer(kind=8), intent(out) :: coloff(:)
    integer(kind=kint) :: j, t0

    t0 = fct%tptr(s)
    nt = fct%tptr(s+1) - t0 - 1
    ntc = fct%ntc(s)
    nrow = fct%tbnd(t0+nt)
    ncol = fct%tbnd(t0+ntc)
    coloff(1) = 0
    do j = 1, nt
      coloff(j+1) = coloff(j) + int(nrow - fct%tbnd(t0+j-1), 8)*(fct%tbnd(t0+j) - fct%tbnd(t0+j-1))
    enddo
  end subroutine mf_layout

  !> Word offset (0-based) of tile (i,j), i >= j, in a lower tile grid with boundaries tb(0:nt).
  function mf_tile_off(tb, coloff, i, j) result(off)
    implicit none
    integer(kind=kint), intent(in) :: tb(0:)
    integer(kind=8), intent(in) :: coloff(:)
    integer(kind=kint), intent(in) :: i, j
    integer(kind=8) :: off

    off = coloff(j) + int(tb(i-1) - tb(j-1), 8)*(tb(j) - tb(j-1))
  end function mf_tile_off

  !> Word index (1-based) of entry (r,c), r >= c, both 0-based DOFs of the grid.
  function mf_elem(tb, coloff, dtile, r, c) result(idx)
    implicit none
    integer(kind=kint), intent(in) :: tb(0:)
    integer(kind=8), intent(in) :: coloff(:)
    integer(kind=kint), intent(in) :: dtile(0:)
    integer(kind=kint), intent(in) :: r, c
    integer(kind=8) :: idx
    integer(kind=kint) :: ti, tj

    ti = dtile(r)
    tj = dtile(c)
    idx = mf_tile_off(tb, coloff, ti, tj) + int(c - tb(tj-1), 8)*(tb(ti) - tb(ti-1)) + (r - tb(ti-1)) + 1
  end function mf_elem

  !> Factor the symmetric matrix of hecMAT (lower part referenced) in the order of sym. ierr is
  !> 0 on success, the permuted DOF of the first pivot that is not positive, or -1 when the block
  !> size of hecMAT does not match the structure.
  subroutine hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(out) :: ierr
    integer(kind=8), allocatable :: coloff(:), ccoloff(:)
    integer(kind=kint), allocatable :: rdof(:), pos(:), dtile(:), cmapdof(:), ctb(:)
    integer(kind=kint) :: ns, s, c, nt, ntc, nrow, ncol, nown, nrow_nodes, i, j, k, l, t0, nd, a, b
    integer(kind=kint) :: j0, ki, kk, coff, roff, cnt, cntc, cnrow, cncol, cnown, m, r, cc, hi, wj, hk, info
    integer(kind=8) :: top, base, o, oik, ojk, oij, okk
    real(kind=kreal) :: v

    ierr = 0
    fct%factored = .false.
    ns = fct%nsuper
    do k = 1, sym%nnode
      if (sym%ndof(k) /= hecMAT%NDOF) then
        ierr = -1
        return
      endif
    enddo
    nd = hecMAT%NDOF
    allocate(coloff(fct%max_tiles+1), ccoloff(fct%max_tiles+1), ctb(0:fct%max_tiles))
    allocate(rdof(fct%max_rows+1), pos(sym%nnode), dtile(0:sym%max_front-1), cmapdof(0:sym%max_front-1))
    top = 0

    do s = 1, ns
      nown = sym%sptr(s+1) - sym%sptr(s)
      nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
      rdof(1) = 0
      do i = 1, nrow_nodes
        k = sym%rlist(sym%rptr(s)+i-1)
        rdof(i+1) = rdof(i) + sym%ndof(k)
        pos(k) = i
      enddo
      call mf_layout(fct, s, nt, ntc, nrow, ncol, coloff)
      t0 = fct%tptr(s)
      do j = 1, nt
        dtile(fct%tbnd(t0+j-1):fct%tbnd(t0+j)-1) = j
      enddo
      fct%fval(1:coloff(nt+1)) = 0.0d0

      ! scatter the lower part of the permuted matrix; every neighbor above the diagonal in the
      ! original numbering supplies the transposed block
      do i = 1, nown
        k = sym%sptr(s) + i - 1
        j0 = sym%perm(k)
        coff = rdof(i)
        base = int(j0-1, 8)*nd*nd
        do b = 1, nd
          do a = b, nd
            o = mf_elem(fct%tbnd(t0:t0+nt), coloff, dtile, coff+a-1, coff+b-1)
            fct%fval(o) = fct%fval(o) + hecMAT%D(base + (a-1)*nd + b)
          enddo
        enddo
        do kk = hecMAT%indexL(j0-1)+1, hecMAT%indexL(j0)
          ki = sym%invp(hecMAT%itemL(kk))
          if (ki <= k) cycle
          roff = rdof(pos(ki))
          base = int(kk-1, 8)*nd*nd
          do b = 1, nd
            do a = 1, nd
              o = mf_elem(fct%tbnd(t0:t0+nt), coloff, dtile, roff+a-1, coff+b-1)
              fct%fval(o) = fct%fval(o) + hecMAT%AL(base + (b-1)*nd + a)
            enddo
          enddo
        enddo
        do kk = hecMAT%indexU(j0-1)+1, hecMAT%indexU(j0)
          ki = sym%invp(hecMAT%itemU(kk))
          if (ki <= k) cycle
          roff = rdof(pos(ki))
          base = int(kk-1, 8)*nd*nd
          do b = 1, nd
            do a = 1, nd
              o = mf_elem(fct%tbnd(t0:t0+nt), coloff, dtile, roff+a-1, coff+b-1)
              fct%fval(o) = fct%fval(o) + hecMAT%AU(base + (b-1)*nd + a)
            enddo
          enddo
        enddo
      enddo

      ! extend-add of the children popped from the stack
      c = fct%chead(s)
      do while (c /= 0)
        top = top - fct%cbsize(c)
        cnown = sym%sptr(c+1) - sym%sptr(c)
        call mf_layout(fct, c, cnt, cntc, cnrow, cncol, ccoloff)
        do j = 0, cnt - cntc
          ctb(j) = fct%tbnd(fct%tptr(c)+cntc+j) - cncol
        enddo
        ccoloff(1) = 0
        do j = 1, cnt - cntc
          ccoloff(j+1) = ccoloff(j) + int(ctb(cnt-cntc) - ctb(j-1), 8)*(ctb(j) - ctb(j-1))
        enddo
        m = 0
        do l = sym%rptr(c) + cnown, sym%rptr(c+1) - 1
          k = sym%cmap(sym%cmap_ptr(c) + (l - sym%rptr(c) - cnown))
          do a = 1, sym%ndof(sym%rlist(l))
            cmapdof(m) = rdof(k) + a - 1
            m = m + 1
          enddo
        enddo
        do j = 1, cnt - cntc
          do i = j, cnt - cntc
            base = top + mf_tile_off(ctb, ccoloff, i, j)
            hi = ctb(i) - ctb(i-1)
            do cc = ctb(j-1), ctb(j) - 1
              do r = max(cc, ctb(i-1)), ctb(i) - 1
                v = fct%sval(base + int(cc - ctb(j-1), 8)*hi + (r - ctb(i-1)) + 1)
                o = mf_elem(fct%tbnd(t0:t0+nt), coloff, dtile, cmapdof(r), cmapdof(cc))
                fct%fval(o) = fct%fval(o) + v
              enddo
            enddo
          enddo
        enddo
        c = fct%cnext(c)
      enddo

      ! right looking tiled LDLt of the pivot tile columns, updating the whole lower grid
      do k = 1, ntc
        hk = fct%tbnd(t0+k) - fct%tbnd(t0+k-1)
        okk = mf_tile_off(fct%tbnd(t0:t0+nt), coloff, k, k)
        call hecmw_mf_kernel_ldlt(hk, fct%fval(okk+1), info)
        if (info > 0) then
          ierr = fct%cdofptr(sym%sptr(s)) + fct%tbnd(t0+k-1) + info - 1
          deallocate(coloff, ccoloff, ctb, rdof, pos, dtile, cmapdof)
          return
        endif
        do i = k+1, nt
          hi = fct%tbnd(t0+i) - fct%tbnd(t0+i-1)
          oik = mf_tile_off(fct%tbnd(t0:t0+nt), coloff, i, k)
          call hecmw_mf_kernel_trsm(hi, hk, fct%fval(okk+1), fct%fval(oik+1))
        enddo
        do i = k+1, nt
          hi = fct%tbnd(t0+i) - fct%tbnd(t0+i-1)
          oik = mf_tile_off(fct%tbnd(t0:t0+nt), coloff, i, k)
          call hecmw_mf_kernel_scale(hi, hk, fct%fval(okk+1), fct%fval(oik+1), fct%wval)
          do j = k+1, i
            wj = fct%tbnd(t0+j) - fct%tbnd(t0+j-1)
            ojk = mf_tile_off(fct%tbnd(t0:t0+nt), coloff, j, k)
            oij = mf_tile_off(fct%tbnd(t0:t0+nt), coloff, i, j)
            call hecmw_mf_kernel_gemm(hi, wj, hk, fct%wval, fct%fval(ojk+1), fct%fval(oij+1))
          enddo
        enddo
      enddo

      fct%lval(fct%lptr(s):fct%lptr(s+1)-1) = fct%fval(1:coloff(ntc+1))
      fct%sval(top+1:top+fct%cbsize(s)) = fct%fval(coloff(ntc+1)+1:coloff(nt+1))
      top = top + fct%cbsize(s)
    enddo

    fct%factored = .true.
    deallocate(coloff, ccoloff, ctb, rdof, pos, dtile, cmapdof)
  end subroutine hecmw_mf_numeric_factor

  !> x = A^-1 b with the factor; b and x are in the original DOF numbering.
  subroutine hecmw_mf_numeric_solve(sym, fct, b, x)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(in) :: fct
    real(kind=kreal), intent(in) :: b(:)
    real(kind=kreal), intent(out) :: x(:)
    real(kind=kreal), allocatable :: y(:), v(:)
    integer(kind=8), allocatable :: coloff(:)
    integer(kind=kint) :: ns, s, nt, ntc, nrow, ncol, t0, i, k, hi, hk, d
    integer(kind=8) :: okk, oik

    ns = fct%nsuper
    allocate(y(fct%ndof_tot), v(sym%max_front), coloff(fct%max_tiles+1))
    do i = 1, fct%ndof_tot
      y(i) = b(fct%pdof(i))
    enddo

    do s = 1, ns
      call mf_layout(fct, s, nt, ntc, nrow, ncol, coloff)
      t0 = fct%tptr(s)
      call mf_gather(sym, fct, s, nrow, y, v)
      do k = 1, ntc
        hk = fct%tbnd(t0+k) - fct%tbnd(t0+k-1)
        okk = fct%lptr(s) + mf_tile_off(fct%tbnd(t0:t0+nt), coloff, k, k)
        call hecmw_mf_kernel_trsv(hk, fct%lval(okk), v(fct%tbnd(t0+k-1)+1))
        do i = k+1, nt
          hi = fct%tbnd(t0+i) - fct%tbnd(t0+i-1)
          oik = fct%lptr(s) + mf_tile_off(fct%tbnd(t0:t0+nt), coloff, i, k)
          call hecmw_mf_kernel_gemv(hi, hk, fct%lval(oik), v(fct%tbnd(t0+k-1)+1), v(fct%tbnd(t0+i-1)+1))
        enddo
        do d = 1, hk
          v(fct%tbnd(t0+k-1)+d) = v(fct%tbnd(t0+k-1)+d) / fct%lval(okk + int(d-1, 8)*hk + d - 1)
        enddo
      enddo
      call mf_scatter(sym, fct, s, nrow, v, y)
    enddo

    do s = ns, 1, -1
      call mf_layout(fct, s, nt, ntc, nrow, ncol, coloff)
      t0 = fct%tptr(s)
      call mf_gather(sym, fct, s, nrow, y, v)
      do k = ntc, 1, -1
        hk = fct%tbnd(t0+k) - fct%tbnd(t0+k-1)
        do i = k+1, nt
          hi = fct%tbnd(t0+i) - fct%tbnd(t0+i-1)
          oik = fct%lptr(s) + mf_tile_off(fct%tbnd(t0:t0+nt), coloff, i, k)
          call hecmw_mf_kernel_gemv_t(hi, hk, fct%lval(oik), v(fct%tbnd(t0+i-1)+1), v(fct%tbnd(t0+k-1)+1))
        enddo
        okk = fct%lptr(s) + mf_tile_off(fct%tbnd(t0:t0+nt), coloff, k, k)
        call hecmw_mf_kernel_trsv_t(hk, fct%lval(okk), v(fct%tbnd(t0+k-1)+1))
      enddo
      call mf_scatter(sym, fct, s, ncol, v, y)
    enddo

    do i = 1, fct%ndof_tot
      x(fct%pdof(i)) = y(i)
    enddo
    deallocate(y, v, coloff)
  end subroutine hecmw_mf_numeric_solve

  !> v(1:n) <- the first n front DOFs of supernode s taken from y (permuted DOF numbering).
  subroutine mf_gather(sym, fct, s, n, y, v)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(in) :: fct
    integer(kind=kint), intent(in) :: s, n
    real(kind=kreal), intent(in) :: y(:)
    real(kind=kreal), intent(out) :: v(:)
    integer(kind=kint) :: l, k, d, m

    m = 0
    do l = sym%rptr(s), sym%rptr(s+1) - 1
      k = sym%rlist(l)
      do d = fct%cdofptr(k), fct%cdofptr(k+1) - 1
        if (m == n) return
        m = m + 1
        v(m) = y(d)
      enddo
    enddo
  end subroutine mf_gather

  subroutine mf_scatter(sym, fct, s, n, v, y)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(in) :: fct
    integer(kind=kint), intent(in) :: s, n
    real(kind=kreal), intent(in) :: v(:)
    real(kind=kreal), intent(inout) :: y(:)
    integer(kind=kint) :: l, k, d, m

    m = 0
    do l = sym%rptr(s), sym%rptr(s+1) - 1
      k = sym%rlist(l)
      do d = fct%cdofptr(k), fct%cdofptr(k+1) - 1
        if (m == n) return
        m = m + 1
        y(d) = v(m)
      enddo
    enddo
  end subroutine mf_scatter

  subroutine hecmw_mf_numeric_print(fct)
    implicit none
    type(hecmwST_mf_factor), intent(in) :: fct

    write(*,'(a,i0,a,i0,a,i0)') '[DIRECTmf]: tile size = ', fct%tile, ', tiles (max per front) = ', fct%max_tiles, &
      ', tile dim (max) = ', fct%max_tile_dim
    write(*,'(a,i0,a,f10.3,a)') '[DIRECTmf]: factor words = ', fct%factor_words, &
      ' (', real(fct%factor_words, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
    write(*,'(a,i0,a,i0,a,f10.3,a)') '[DIRECTmf]: stack peak words = ', fct%stack_peak, ', front words = ', &
      fct%front_words, ' (', real(fct%stack_peak + fct%front_words, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
  end subroutine hecmw_mf_numeric_print

end module hecmw_mf_numeric
