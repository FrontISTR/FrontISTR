!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Numeric factorization (tiled multifrontal LDLt with threshold pivoting and delayed
!>        pivots) and the triangular solves, on the read-only structure of hecmw_mf_symbolic.
!>
!> The positions of a front are the own DOFs of the supernode, the DOFs delayed by its children
!> (together the ncol fully summed positions) and the DOFs of the contribution rows. The front
!> is partitioned into tiles, cut at node boundaries where possible, and only the tiles on or
!> below the diagonal are stored, tile column by tile column, each tile column major. Pivots are
!> eliminated in position order; a fully summed column without an acceptable pivot is exchanged
!> to the end of the fully summed part and joins the contribution block, so that the factor
!> panel of a front is the leading npiv columns of its tile grid. Front sizes, tile partitions
!> and the factor layout therefore depend on the values and are rebuilt by every factorization.
!> The factor panel, the position metadata and the contribution block of a supernode are stored
!> per supernode; a contribution block lives from the factorization of its supernode until the
!> extend-add of the parent consumes it.
!>
!> A matrix whose values are found unsymmetric (or whose symmetric flag is off) is factored in
!> LU mode: the part of the front above the diagonal is held transposed in a second grid of the
!> same layout (fvalu, and uval for the stored U panels), pivots are chosen by threshold partial
!> pivoting among the fully summed rows, and the row DOF of a position (frow) may then differ
!> from its column DOF (fsdof). Delayed positions carry both to the parent.
module hecmw_mf_numeric
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_mf_symbolic
  use hecmw_mf_dist
  use hecmw_mf_kernel
  !$ use omp_lib
  implicit none

  private
  public :: hecmwST_mf_factor
  public :: hecmw_mf_numeric_init
  public :: hecmw_mf_numeric_factor
  public :: hecmw_mf_numeric_factor_mpi
  public :: hecmw_mf_numeric_solve
  public :: hecmw_mf_numeric_solve_mpi
  public :: hecmw_mf_numeric_print
  public :: hecmw_mf_numeric_finalize

  !> default pivot threshold: entries of L are bounded by 1/u
  real(kind=kreal), parameter :: MF_PIVOT_U = 0.01d0
  !> Bunch-Kaufman constant
  real(kind=kreal), parameter :: MF_PIVOT_ALPHA = (1.0d0 + sqrt(17.0d0)) / 8.0d0
  !> default zero pivot fraction of the largest entry of the matrix
  real(kind=kreal), parameter :: MF_PIVOT_ZERO = 1.0d-14
  !> default BLR truncation threshold
  real(kind=kreal), parameter :: MF_BLR_EPS = 1.0d-8
  !> the matrix is factored in LU mode above this asymmetry, max|A_ij - A_ji| / max|A_ij|
  real(kind=kreal), parameter :: MF_ASYM_TOL = 1.0d-12
  !> a front with at least this many rows parallelizes its tiles across the team
  integer(kind=kint), parameter :: MF_PAR_ROWS = 512

  !> per-supernode part of the factorization; cval holds the contribution block from the
  !> factorization of the supernode until the extend-add of the parent frees it, dval holds
  !> the forward-solve contribution the same way.
  !> With BLR a stored panel tile may be low rank: lval then holds U (hi x r) followed by
  !> V (tw x r) at bptr with the tile ~ U V^T, a full rank tile keeps the plain hi x tw
  !> layout, and brank(-1 = full rank) tells them apart; uval has its own bptru/branku
  type mf_snode
    integer(kind=kint) :: ncol = 0                  !< fully summed positions
    integer(kind=kint) :: npiv = 0                  !< pivots eliminated (positions 1:npiv)
    integer(kind=kint) :: nt = 0                    !< tiles of the front
    integer(kind=kint) :: ntc = 0                   !< tiles covering the fully summed positions
    integer(kind=8) :: cbsize = 0                   !< words of the contribution block
    integer(kind=kint), allocatable :: tbnd(:)      !< tile boundaries, tbnd(1:nt+1) from 0 to nrow
    integer(kind=kint), allocatable :: fsdof(:)     !< permuted DOF of the column at a position
    integer(kind=kint), allocatable :: frow(:)      !< permuted DOF of the row at a position
    integer(kind=kint), allocatable :: ptype(:)     !< pivot type at a position (see hecmw_mf_kernel)
    real(kind=kreal), allocatable :: dsub(:)        !< off-diagonal entry of a 2x2 pivot at its first position
    integer(kind=8), allocatable :: bptr(:)         !< BLR: word offset of panel tile (i,k) in lval at
                                                    !< mf_bidx(nt,k,i), bptr(ntile+1) the total words
    integer(kind=kint), allocatable :: brank(:)     !< BLR: rank of a panel tile, -1 = full rank
    integer(kind=8), allocatable :: bptru(:)        !< BLR: the same for uval (LU mode)
    integer(kind=kint), allocatable :: branku(:)
    real(kind=kreal), allocatable :: lval(:)        !< factor panel, the leading npiv tile grid columns
    real(kind=kreal), allocatable :: uval(:)        !< U panel of LU mode, same layout as lval
    real(kind=kreal), allocatable :: cval(:)        !< contribution block (lower face, then upper in LU mode)
    real(kind=kreal), allocatable :: dval(:)        !< forward-solve contribution of the rows beyond npiv
  end type mf_snode

  type hecmwST_mf_factor
    integer(kind=kint) :: tile = 0
    integer(kind=kint) :: nnode = 0
    integer(kind=kint) :: nsuper = 0
    integer(kind=kint) :: ndof_tot = 0
    integer(kind=kint), allocatable :: cdofptr(:)   !< column k owns the permuted DOFs cdofptr(k):cdofptr(k+1)-1
    integer(kind=kint), allocatable :: pdof(:)      !< pdof(c): original DOF of permuted DOF c
    integer(kind=kint), allocatable :: chead(:)     !< children of s: chead(s), cnext(...), in decreasing
    integer(kind=kint), allocatable :: cnext(:)     !< order, which is the order they are popped
    integer(kind=kint), allocatable :: cptr(:)      !< the same children in ascending order,
    integer(kind=kint), allocatable :: clist(:)     !< clist(cptr(s):cptr(s+1)-1), for the task dependences
    integer(kind=kint), allocatable :: mirror(:)    !< AU entry holding the transpose of AL entry k
    integer(kind=kint), allocatable :: mirroru(:)   !< AL entry holding the transpose of AU entry k
    ! numeric options; the driver sets them from the solver option lines before every
    ! factorization, 0 (or below) selecting the built-in default, and the factorization
    ! writes the effective values back
    integer(kind=kint) :: mode = 0                  !< 0 auto, 1 force LDLt, 2 force LU
    logical :: blr = .false.                        !< compress the factor panels (BLR)
    real(kind=kreal) :: eps = 0.0d0                 !< BLR truncation threshold
    real(kind=kreal) :: pivot_u = 0.0d0             !< pivot threshold: entries of L bounded by 1/u
    real(kind=kreal) :: pivot_zero = 0.0d0          !< zero pivot fraction of max|A|
    ! layout of the last factorization
    logical :: lu = .false.                         !< LU mode (else LDLt)
    real(kind=kreal) :: asym = 0.0d0                !< max|A_ij - A_ji| / max|A_ij| of the last matrix
    type(mf_snode), allocatable :: sn(:)            !< per-supernode factor storage
    ! estimates from the symbolic structure (no delayed pivots)
    integer(kind=8) :: factor_words = 0
    integer(kind=8) :: factor_nnz = 0               !< panel entries on or below the diagonal
    integer(kind=8) :: stack_peak = 0
    integer(kind=8) :: front_words = 0
    integer(kind=kint) :: max_tiles = 0
    integer(kind=kint) :: max_rows = 0
    integer(kind=kint) :: max_tile_dim = 0
    ! results of the last factorization
    integer(kind=8) :: factor_words_act = 0
    integer(kind=8) :: stack_peak_act = 0
    integer(kind=8) :: front_words_act = 0
    integer(kind=8) :: live_cb = 0                  !< contribution block words currently held
    integer(kind=8) :: live_front = 0               !< front work words currently held
    integer(kind=8) :: front_peak = 0               !< peak of live_front over the factorization
    integer(kind=kint) :: n_pos = 0                 !< inertia: positive, negative eigenvalues of D
    integer(kind=kint) :: n_neg = 0
    integer(kind=kint) :: n_2x2 = 0
    integer(kind=kint) :: n_swap = 0
    integer(kind=kint) :: n_delay = 0               !< delay events (a DOF delayed twice counts twice)
    integer(kind=kint) :: max_growth = 0            !< largest number of delayed DOFs received by a front
    integer(kind=8) :: blr_words_fr = 0             !< words the stored panels would take full rank
    integer(kind=8) :: blr_tiles = 0                !< panel tiles offered to the compression
    integer(kind=8) :: blr_tiles_lr = 0             !< panel tiles kept low rank
    integer(kind=8) :: blr_rank_sum = 0
    integer(kind=kint) :: blr_rank_max = 0
    logical :: factored = .false.
  end type hecmwST_mf_factor

  !> tile grid of the front being processed; tb(0:nt) are the position offsets of the tiles,
  !> the first ntc tiles cover the ncol fully summed positions
  type mf_grid
    integer(kind=kint) :: nt = 0
    integer(kind=kint) :: ntc = 0
    integer(kind=kint) :: nrow = 0
    integer(kind=kint) :: ncol = 0
    integer(kind=kint), allocatable :: tb(:)
    integer(kind=8), allocatable :: coloff(:)
    integer(kind=kint), allocatable :: dtile(:)
  end type mf_grid

  !> work space of the front being processed, one instance per concurrently processed front
  type mf_work
    type(mf_grid) :: g
    real(kind=kreal), allocatable :: fval(:)        !< front being assembled
    real(kind=kreal), allocatable :: fvalu(:)       !< its part above the diagonal, transposed (LU mode)
    real(kind=kreal), allocatable :: pval(:)        !< dense work panel of the pivot search
    real(kind=kreal), allocatable :: pvalu(:)       !< its transposed upper part (LU mode)
    real(kind=kreal), allocatable :: wval(:)        !< work tile of the scaled panel L*D
    real(kind=kreal), allocatable :: wk(:)
    integer(kind=kint), allocatable :: rowoff(:)
    integer(kind=kint), allocatable :: pos(:)
    integer(kind=kint), allocatable :: blk(:)
    integer(kind=kint), allocatable :: blkr(:)
    integer(kind=kint), allocatable :: cmapdof(:)
    integer(kind=kint), allocatable :: ctb(:)
    integer(kind=8), allocatable :: ccoloff(:)
    real(kind=kreal), allocatable :: bval(:)        !< BLR: U, V (and D*V in LDLt mode) of the
    integer(kind=8), allocatable :: boff(:)         !< compressed panel tiles, their slot offsets
    integer(kind=kint), allocatable :: brk(:)       !< and ranks (-1 = kept full rank)
    real(kind=kreal), allocatable :: bvalu(:)       !< BLR: the same for the upper grid (LU mode)
    integer(kind=8), allocatable :: boffu(:)
    integer(kind=kint), allocatable :: brku(:)
  end type mf_work

  !> send buffers of one fan-in message that must outlive the isend until the waitall
  type mf_sendbox
    integer(kind=kint), allocatable :: hdr(:)
    integer(kind=kint), allocatable :: tl(:)
    real(kind=kreal), allocatable :: rv(:)
  end type mf_sendbox

contains

  !> Estimates of the factor and stack sizes from the symbolic structure and the work space;
  !> tile is the target tile size in DOFs.
  subroutine hecmw_mf_numeric_init(sym, tile, fct)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    integer(kind=kint), intent(in) :: tile
    type(hecmwST_mf_factor), intent(out) :: fct
    type(mf_grid) :: g
    integer(kind=8), allocatable :: cbsize(:)
    integer(kind=kint), allocatable :: odofptr(:), wptr(:)
    integer(kind=kint) :: n, ns, s, k, d, c, j, hj
    integer(kind=8) :: top

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

    allocate(fct%chead(ns), fct%cnext(ns))
    fct%chead(1:ns) = 0
    do s = 1, ns
      c = sym%sparent(s)
      if (c == 0) cycle
      fct%cnext(s) = fct%chead(c)
      fct%chead(c) = s
    enddo
    allocate(fct%cptr(ns+1))
    fct%cptr(1:ns+1) = 0
    do s = 1, ns
      c = sym%sparent(s)
      if (c /= 0) fct%cptr(c+1) = fct%cptr(c+1) + 1
    enddo
    fct%cptr(1) = 1
    do s = 1, ns
      fct%cptr(s+1) = fct%cptr(s) + fct%cptr(s+1)
    enddo
    allocate(fct%clist(max(fct%cptr(ns+1)-1, 1)), wptr(ns))
    wptr(1:ns) = fct%cptr(1:ns)
    do s = 1, ns
      c = sym%sparent(s)
      if (c /= 0) then
        fct%clist(wptr(c)) = s
        wptr(c) = wptr(c) + 1
      endif
    enddo
    deallocate(wptr)

    allocate(cbsize(ns))
    fct%factor_words = 0
    fct%factor_nnz = 0
    fct%front_words = 0
    fct%max_tiles = 0
    fct%max_rows = 0
    fct%max_tile_dim = tile
    do s = 1, ns
      call mf_partition(sym, fct, s, 0, g)
      fct%factor_words = fct%factor_words + g%coloff(g%ntc+1)
      cbsize(s) = g%coloff(g%nt+1) - g%coloff(g%ntc+1)
      fct%front_words = max(fct%front_words, g%coloff(g%nt+1))
      fct%factor_nnz = fct%factor_nnz + g%coloff(g%ntc+1)
      do k = 1, g%ntc
        hj = g%tb(k) - g%tb(k-1)
        fct%factor_nnz = fct%factor_nnz - int(hj, 8)*(hj-1)/2
      enddo
      fct%max_tiles = max(fct%max_tiles, g%nt)
      fct%max_rows = max(fct%max_rows, g%nrow)
      do j = 1, g%nt
        fct%max_tile_dim = max(fct%max_tile_dim, g%tb(j) - g%tb(j-1))
      enddo
    enddo
    top = 0
    fct%stack_peak = 0
    do s = 1, ns
      fct%stack_peak = max(fct%stack_peak, top)
      c = fct%chead(s)
      do while (c /= 0)
        top = top - cbsize(c)
        c = fct%cnext(c)
      enddo
      top = top + cbsize(s)
      fct%stack_peak = max(fct%stack_peak, top)
    enddo
    deallocate(cbsize)

    allocate(fct%sn(ns))
    fct%factored = .false.
  end subroutine hecmw_mf_numeric_init

  subroutine hecmw_mf_numeric_finalize(fct)
    implicit none
    type(hecmwST_mf_factor), intent(inout) :: fct

    fct%nsuper = 0
    fct%factored = .false.
    if (allocated(fct%cdofptr)) deallocate(fct%cdofptr)
    if (allocated(fct%pdof)) deallocate(fct%pdof)
    if (allocated(fct%chead)) deallocate(fct%chead)
    if (allocated(fct%cnext)) deallocate(fct%cnext)
    if (allocated(fct%cptr)) deallocate(fct%cptr)
    if (allocated(fct%clist)) deallocate(fct%clist)
    if (allocated(fct%mirror)) deallocate(fct%mirror)
    if (allocated(fct%mirroru)) deallocate(fct%mirroru)
    if (allocated(fct%sn)) deallocate(fct%sn)
  end subroutine hecmw_mf_numeric_finalize

  subroutine mf_grow_r(a, need)
    implicit none
    real(kind=kreal), allocatable, intent(inout) :: a(:)
    integer(kind=8), intent(in) :: need
    real(kind=kreal), allocatable :: t(:)
    integer(kind=8) :: n

    if (allocated(a)) then
      if (size(a, kind=8) >= need) return
      n = max(need, size(a, kind=8) + size(a, kind=8)/2)
      allocate(t(n))
      t(1:size(a, kind=8)) = a(:)
      call move_alloc(t, a)
    else
      allocate(a(max(need, 1_8)))
    endif
  end subroutine mf_grow_r

  subroutine mf_grow_i(a, need)
    implicit none
    integer(kind=kint), allocatable, intent(inout) :: a(:)
    integer(kind=kint), intent(in) :: need
    integer(kind=kint), allocatable :: t(:)
    integer(kind=kint) :: n

    if (allocated(a)) then
      if (size(a) >= need) return
      n = max(need, size(a) + size(a)/2)
      allocate(t(n))
      t(1:size(a)) = a(:)
      call move_alloc(t, a)
    else
      allocate(a(max(need, 1)))
    endif
  end subroutine mf_grow_i

  subroutine mf_grow_i8(a, need)
    implicit none
    integer(kind=8), allocatable, intent(inout) :: a(:)
    integer(kind=kint), intent(in) :: need
    integer(kind=8), allocatable :: t(:)
    integer(kind=kint) :: n

    if (allocated(a)) then
      if (size(a) >= need) return
      n = max(need, size(a) + size(a)/2)
      allocate(t(n))
      t(1:size(a)) = a(:)
      call move_alloc(t, a)
    else
      allocate(a(max(need, 1)))
    endif
  end subroutine mf_grow_i8

  !> Tile grid of supernode s whose fully summed part holds ndel delayed DOFs after the own
  !> DOFs: a tile ends at the last node boundary that keeps it within the target size, the own
  !> DOFs and the fully summed part end on tile boundaries, delayed DOFs are cut in tile sized
  !> chunks. coloff(j) is the word offset of tile column j, coloff(nt+1) the size of the grid.
  subroutine mf_partition(sym, fct, s, ndel, g)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(in) :: fct
    integer(kind=kint), intent(in) :: s, ndel
    type(mf_grid), intent(inout) :: g
    integer(kind=kint) :: nown, nrow_nodes, ncol0, ncb, nt, i, j, nd, cur, acc, tile, pass

    tile = fct%tile
    nown = sym%sptr(s+1) - sym%sptr(s)
    nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
    ncol0 = fct%cdofptr(sym%sptr(s+1)) - fct%cdofptr(sym%sptr(s))
    ncb = 0
    do i = nown + 1, nrow_nodes
      ncb = ncb + sym%ndof(sym%rlist(sym%rptr(s)+i-1))
    enddo
    g%ncol = ncol0 + ndel
    g%nrow = g%ncol + ncb
    do pass = 1, 2
      nt = 0
      cur = 0
      acc = 0
      do i = 1, nrow_nodes
        if (i == nown + 1) then
          do j = 1, ndel, tile
            nt = nt + 1
            cur = acc + j - 1
            if (pass == 2) g%tb(nt) = cur
          enddo
          acc = acc + ndel
          nt = nt + 1
          if (pass == 2) g%tb(nt) = acc
          cur = acc
          g%ntc = nt
        endif
        nd = sym%ndof(sym%rlist(sym%rptr(s)+i-1))
        if (acc > cur .and. acc + nd - cur > tile) then
          nt = nt + 1
          if (pass == 2) g%tb(nt) = acc
          cur = acc
        endif
        acc = acc + nd
      enddo
      if (nown == nrow_nodes) then
        do j = 1, ndel, tile
          nt = nt + 1
          cur = acc + j - 1
          if (pass == 2) g%tb(nt) = cur
        enddo
        acc = acc + ndel
      endif
      nt = nt + 1
      if (pass == 2) g%tb(nt) = acc
      if (nown == nrow_nodes) g%ntc = nt
      if (pass == 1) then
        if (allocated(g%tb)) then
          if (size(g%tb) < nt + 1) deallocate(g%tb, g%coloff)
        endif
        if (.not. allocated(g%tb)) allocate(g%tb(0:nt), g%coloff(nt+1))
        if (allocated(g%dtile)) then
          if (size(g%dtile) < g%nrow) deallocate(g%dtile)
        endif
        if (.not. allocated(g%dtile)) allocate(g%dtile(0:max(g%nrow, 1)-1))
        g%tb(0) = 0
      endif
    enddo
    g%nt = nt
    g%coloff(1) = 0
    do j = 1, nt
      g%coloff(j+1) = g%coloff(j) + int(g%nrow - g%tb(j-1), 8)*(g%tb(j) - g%tb(j-1))
      g%dtile(g%tb(j-1):g%tb(j)-1) = j
    enddo
  end subroutine mf_partition

  !> Word offset (0-based) of tile (i,j), i >= j.
  function mf_off(g, i, j) result(off)
    implicit none
    type(mf_grid), intent(in) :: g
    integer(kind=kint), intent(in) :: i, j
    integer(kind=8) :: off

    off = g%coloff(j) + int(g%tb(i-1) - g%tb(j-1), 8)*(g%tb(j) - g%tb(j-1))
  end function mf_off

  !> Word index (1-based) of the entry at positions (r,c), r >= c, both 1-based.
  function mf_idx(g, r, c) result(idx)
    implicit none
    type(mf_grid), intent(in) :: g
    integer(kind=kint), intent(in) :: r, c
    integer(kind=8) :: idx
    integer(kind=kint) :: ti, tj

    ti = g%dtile(r-1)
    tj = g%dtile(c-1)
    idx = mf_off(g, ti, tj) + int(c - 1 - g%tb(tj-1), 8)*(g%tb(ti) - g%tb(ti-1)) + (r - 1 - g%tb(ti-1)) + 1
  end function mf_idx

  !> Index of the stored panel tile (i,k), i >= k, in the per column enumeration of a front
  !> with nt tiles (every tile of a column, then the next column).
  function mf_bidx(nt, k, i) result(idx)
    implicit none
    integer(kind=kint), intent(in) :: nt, k, i
    integer(kind=kint) :: idx

    idx = (k-1)*(nt+1) - (k-1)*k/2 + (i - k + 1)
  end function mf_bidx

  !> c(m,n) <- c - A * B^T where A (m x k) is fa or the low rank ua(m,ra) va(k,ra)^T and
  !> B (n x k) is fb or ub(n,rb) vb(k,rb)^T. A rank of -1 selects the full rank operand,
  !> a rank of 0 makes the product zero.
  subroutine mf_update_ab(m, n, k, ra, fa, ua, va, rb, fb, ub, vb, c)
    implicit none
    integer(kind=kint), intent(in) :: m, n, k, ra, rb
    real(kind=kreal), intent(in) :: fa(*), ua(*), va(*), fb(*), ub(*), vb(*)
    real(kind=kreal), intent(inout) :: c(*)
    real(kind=kreal), allocatable :: s(:), t(:)

    if (ra == 0 .or. rb == 0) return
    if (ra < 0 .and. rb < 0) then
      call hecmw_mf_kernel_gemm(m, n, k, m, fa, n, fb, m, c)
    else if (ra >= 0 .and. rb < 0) then
      allocate(t(int(n, 8)*ra))
      call hecmw_mf_kernel_mult_nn(n, ra, k, n, fb, k, va, n, t)
      call hecmw_mf_kernel_gemm(m, n, ra, m, ua, n, t, m, c)
      deallocate(t)
    else if (ra < 0) then
      allocate(t(int(m, 8)*rb))
      call hecmw_mf_kernel_mult_nn(m, rb, k, m, fa, k, vb, m, t)
      call hecmw_mf_kernel_gemm(m, n, rb, m, t, n, ub, m, c)
      deallocate(t)
    else
      allocate(s(int(ra, 8)*rb), t(int(m, 8)*rb))
      call hecmw_mf_kernel_mult_tn(ra, rb, k, k, va, k, vb, ra, s)
      call hecmw_mf_kernel_mult_nn(m, rb, ra, m, ua, ra, s, m, t)
      call hecmw_mf_kernel_gemm(m, n, rb, m, t, n, ub, m, c)
      deallocate(s, t)
    endif
  end subroutine mf_update_ab

  !> vec(1:nrow-i0+1) <- entries (i,x) of the front for i = i0..nrow (get) or the reverse (put).
  subroutine mf_col_copy(g, fval, x, i0, vec, put)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(inout) :: fval(:)
    integer(kind=kint), intent(in) :: x, i0
    real(kind=kreal), intent(inout) :: vec(*)
    logical, intent(in) :: put
    integer(kind=kint) :: it, kx, r0, r1, h
    integer(kind=8) :: base

    if (i0 > g%nrow) return
    kx = g%dtile(x-1)
    do it = g%dtile(i0-1), g%nt
      h = g%tb(it) - g%tb(it-1)
      r0 = max(i0, g%tb(it-1) + 1)
      r1 = g%tb(it)
      base = mf_off(g, it, kx) + int(x - 1 - g%tb(kx-1), 8)*h + (r0 - 1 - g%tb(it-1))
      if (put) then
        fval(base+1:base+r1-r0+1) = vec(r0-i0+1:r1-i0+1)
      else
        vec(r0-i0+1:r1-i0+1) = fval(base+1:base+r1-r0+1)
      endif
    enddo
  end subroutine mf_col_copy

  !> vec, holding the entries (i,x) for i = x..nrow, receives the update of the pivots at
  !> positions q1..q2 (L and D taken from the front): vec <- vec - L(:,q) D L(x,q)^T.
  subroutine mf_col_update(g, fval, ptype, dsub, x, q1, q2, vec, t)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(in) :: fval(:)
    integer(kind=kint), intent(in) :: ptype(:)
    real(kind=kreal), intent(in) :: dsub(:)
    integer(kind=kint), intent(in) :: x, q1, q2
    real(kind=kreal), intent(inout) :: vec(*)
    real(kind=kreal), intent(inout) :: t(*)
    integer(kind=kint) :: q
    real(kind=kreal) :: l1, l2, d11, d21, d22

    q = q1
    do while (q <= q2)
      if (ptype(q) == 1) then
        t(q-q1+1) = fval(mf_idx(g, x, q)) * fval(mf_idx(g, q, q))
        q = q + 1
      else
        l1 = fval(mf_idx(g, x, q))
        l2 = fval(mf_idx(g, x, q+1))
        d11 = fval(mf_idx(g, q, q))
        d21 = dsub(q)
        d22 = fval(mf_idx(g, q+1, q+1))
        t(q-q1+1) = l1*d11 + l2*d21
        t(q-q1+2) = l1*d21 + l2*d22
        q = q + 2
      endif
    enddo
    call mf_col_axpy(g, fval, x, q1, q2, t, vec)
  end subroutine mf_col_update

  !> vec(1:nrow-i0+1), holding the entries (i,x) of the grid gval for i = i0..nrow, receives
  !> vec <- vec - gval(:,q) t(q) for the columns q = q1..q2.
  subroutine mf_col_axpy(g, gval, i0, q1, q2, t, vec)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(in) :: gval(:)
    integer(kind=kint), intent(in) :: i0, q1, q2
    real(kind=kreal), intent(in) :: t(*)
    real(kind=kreal), intent(inout) :: vec(*)
    integer(kind=kint) :: q, it, kq, r0, r1, h
    integer(kind=8) :: b

    if (i0 > g%nrow) return
    do q = q1, q2
      kq = g%dtile(q-1)
      do it = g%dtile(i0-1), g%nt
        h = g%tb(it) - g%tb(it-1)
        r0 = max(i0, g%tb(it-1) + 1)
        r1 = g%tb(it)
        b = mf_off(g, it, kq) + int(q - 1 - g%tb(kq-1), 8)*h + (r0 - 1 - g%tb(it-1))
        vec(r0-i0+1:r1-i0+1) = vec(r0-i0+1:r1-i0+1) - gval(b+1:b+r1-r0+1) * t(q-q1+1)
      enddo
    enddo
  end subroutine mf_col_axpy

  !> t(1:q2-q1+1) <- the entries (x,q) of the grid gval for the columns q = q1..q2 < x.
  subroutine mf_row_get(g, gval, x, q1, q2, t)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(in) :: gval(:)
    integer(kind=kint), intent(in) :: x, q1, q2
    real(kind=kreal), intent(out) :: t(*)
    integer(kind=kint) :: q

    do q = q1, q2
      t(q-q1+1) = gval(mf_idx(g, x, q))
    enddo
  end subroutine mf_row_get

  !> Symmetric exchange of the fully summed positions x < y of the front.
  subroutine mf_swap(g, fval, x, y)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(inout) :: fval(:)
    integer(kind=kint), intent(in) :: x, y
    integer(kind=kint) :: i
    integer(kind=8) :: a, b
    real(kind=kreal) :: v

    a = mf_idx(g, x, x)
    b = mf_idx(g, y, y)
    v = fval(a)
    fval(a) = fval(b)
    fval(b) = v
    do i = 1, x-1
      a = mf_idx(g, x, i)
      b = mf_idx(g, y, i)
      v = fval(a)
      fval(a) = fval(b)
      fval(b) = v
    enddo
    do i = x+1, y-1
      a = mf_idx(g, i, x)
      b = mf_idx(g, y, i)
      v = fval(a)
      fval(a) = fval(b)
      fval(b) = v
    enddo
    do i = y+1, g%nrow
      a = mf_idx(g, i, x)
      b = mf_idx(g, i, y)
      v = fval(a)
      fval(a) = fval(b)
      fval(b) = v
    enddo
  end subroutine mf_swap

  !> Exchange of the rows x < y (fully summed positions) of the LU front fval, fvalu.
  subroutine mf_swap_row(g, fval, fvalu, x, y)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(inout) :: fval(:), fvalu(:)
    integer(kind=kint), intent(in) :: x, y
    integer(kind=kint) :: c
    integer(kind=8) :: a, b
    real(kind=kreal) :: v

    do c = 1, x
      a = mf_idx(g, x, c)
      b = mf_idx(g, y, c)
      v = fval(a)
      fval(a) = fval(b)
      fval(b) = v
    enddo
    do c = x+1, y-1
      a = mf_idx(g, c, x)
      b = mf_idx(g, y, c)
      v = fvalu(a)
      fvalu(a) = fval(b)
      fval(b) = v
    enddo
    a = mf_idx(g, y, x)
    b = mf_idx(g, y, y)
    v = fvalu(a)
    fvalu(a) = fval(b)
    fval(b) = v
    do c = y+1, g%nrow
      a = mf_idx(g, c, x)
      b = mf_idx(g, c, y)
      v = fvalu(a)
      fvalu(a) = fvalu(b)
      fvalu(b) = v
    enddo
  end subroutine mf_swap_row

  !> Symmetric exchange of the fully summed positions x < y of the LU front fval, fvalu.
  subroutine mf_swap_lu(g, fval, fvalu, x, y)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(inout) :: fval(:), fvalu(:)
    integer(kind=kint), intent(in) :: x, y
    integer(kind=kint) :: i
    integer(kind=8) :: a, b
    real(kind=kreal) :: v

    a = mf_idx(g, x, x)
    b = mf_idx(g, y, y)
    v = fval(a)
    fval(a) = fval(b)
    fval(b) = v
    do i = 1, x-1
      a = mf_idx(g, x, i)
      b = mf_idx(g, y, i)
      v = fval(a)
      fval(a) = fval(b)
      fval(b) = v
      v = fvalu(a)
      fvalu(a) = fvalu(b)
      fvalu(b) = v
    enddo
    do i = x+1, y-1
      a = mf_idx(g, i, x)
      b = mf_idx(g, y, i)
      v = fval(a)
      fval(a) = fvalu(b)
      fvalu(b) = v
      v = fval(b)
      fval(b) = fvalu(a)
      fvalu(a) = v
    enddo
    a = mf_idx(g, y, x)
    v = fval(a)
    fval(a) = fvalu(a)
    fvalu(a) = v
    do i = y+1, g%nrow
      a = mf_idx(g, i, x)
      b = mf_idx(g, i, y)
      v = fval(a)
      fval(a) = fval(b)
      fval(b) = v
      v = fvalu(a)
      fvalu(a) = fvalu(b)
      fvalu(b) = v
    enddo
  end subroutine mf_swap_lu

  !> Pair every entry of AL with the entry of AU at the transposed position (structure only,
  !> built once per structure); ierr is -2 when the structure is not symmetric.
  subroutine mf_mirror(hecMAT, fct, ierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(out) :: ierr
    integer(kind=kint), allocatable :: ptr(:)
    integer(kind=kint) :: i, j, kk, ku, l

    ierr = 0
    if (allocated(fct%mirror)) then
      if (size(fct%mirror) == hecMAT%NPL .and. size(fct%mirroru) == hecMAT%NPU) return
      deallocate(fct%mirror, fct%mirroru)
    endif
    if (hecMAT%NPL /= hecMAT%NPU) then
      ierr = -2
      return
    endif
    allocate(fct%mirror(max(hecMAT%NPL, 1)), fct%mirroru(max(hecMAT%NPU, 1)), ptr(hecMAT%NP))
    fct%mirroru(:) = 0
    ptr(1:hecMAT%NP) = hecMAT%indexU(0:hecMAT%NP-1) + 1
    do j = 1, hecMAT%NP
      do kk = hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i = hecMAT%itemL(kk)
        ku = ptr(i)
        if (ku <= hecMAT%indexU(i)) then
          if (hecMAT%itemU(ku) /= j) ku = 0
        else
          ku = 0
        endif
        if (ku == 0) then
          do l = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
            if (hecMAT%itemU(l) == j) ku = l
          enddo
          if (ku == 0) then
            ierr = -2
            exit
          endif
        else
          ptr(i) = ku + 1
        endif
        fct%mirror(kk) = ku
        fct%mirroru(ku) = kk
      enddo
      if (ierr /= 0) exit
    enddo
    if (ierr == 0 .and. hecMAT%NPU > 0) then
      if (minval(fct%mirroru(1:hecMAT%NPU)) == 0) ierr = -2
    endif
    deallocate(ptr)
    if (ierr /= 0) deallocate(fct%mirror, fct%mirroru)
  end subroutine mf_mirror

  !> Asymmetry of the values, max|A_ij - A_ji| relative to amax, and the resulting mode.
  subroutine mf_asymmetry(hecMAT, fct, amax)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_factor), intent(inout) :: fct
    real(kind=kreal), intent(in) :: amax
    integer(kind=kint) :: nd, i, kk, a, b
    integer(kind=8) :: base, base2
    real(kind=kreal) :: asym

    nd = hecMAT%NDOF
    asym = 0.0d0
    do i = 1, hecMAT%NP
      base = int(i-1, 8)*nd*nd
      do a = 1, nd
        do b = a+1, nd
          asym = max(asym, abs(hecMAT%D(base + (a-1)*nd + b) - hecMAT%D(base + (b-1)*nd + a)))
        enddo
      enddo
    enddo
    do kk = 1, hecMAT%NPL
      base = int(kk-1, 8)*nd*nd
      base2 = int(fct%mirror(kk)-1, 8)*nd*nd
      do a = 1, nd
        do b = 1, nd
          asym = max(asym, abs(hecMAT%AL(base + (a-1)*nd + b) - hecMAT%AU(base2 + (b-1)*nd + a)))
        enddo
      enddo
    enddo
    fct%asym = 0.0d0
    if (amax > 0.0d0) fct%asym = asym / amax
    fct%lu = (.not. hecMAT%symmetric) .or. (fct%asym > MF_ASYM_TOL)
    if (fct%mode == 1) fct%lu = .false.
    if (fct%mode == 2) fct%lu = .true.
  end subroutine mf_asymmetry

  !> Push the contribution block of the grid gval (positions beyond npiv) at sval(top+1:): the
  !> cndel delayed positions as one tile followed by the contribution tiles, cnb tiles of
  !> boundaries ctb and word offsets ccoloff, kbeg being 1 when there is a delayed tile.
  subroutine mf_store_cb(g, gval, npiv, cndel, cnb, kbeg, ctb, ccoloff, sval, top, par)
    implicit none
    type(mf_grid), intent(in) :: g
    real(kind=kreal), intent(in) :: gval(:)
    integer(kind=kint), intent(in) :: npiv, cndel, cnb, kbeg
    integer(kind=kint), intent(in) :: ctb(0:)
    integer(kind=8), intent(in) :: ccoloff(:)
    real(kind=kreal), intent(inout) :: sval(:)
    integer(kind=8), intent(in) :: top
    logical, intent(in) :: par
    integer(kind=kint) :: i, j, hi, cc, r, m
    integer(kind=8) :: base, o

    do j = 1, min(cnb, kbeg)
      do i = j, cnb
        hi = ctb(i) - ctb(i-1)
        base = top + int(ctb(i-1), 8)*cndel
        do cc = 1, cndel
          do r = max(cc, ctb(i-1) + 1), ctb(i)
            if (r <= cndel) then
              m = npiv + r
            else
              m = g%ncol + r - cndel
            endif
            sval(base + int(cc-1, 8)*hi + (r - ctb(i-1))) = gval(mf_idx(g, m, npiv + cc))
          enddo
        enddo
      enddo
    enddo
    !$omp taskloop default(shared) private(i, hi, m, base, o) grainsize(1) if(par)
    do j = g%ntc + 1, g%nt
      do i = j, g%nt
        hi = g%tb(i) - g%tb(i-1)
        m = g%tb(j) - g%tb(j-1)
        base = top + ccoloff(kbeg+j-g%ntc) + int(ctb(kbeg+i-g%ntc-1) - ctb(kbeg+j-g%ntc-1), 8)*m
        o = mf_off(g, i, j)
        sval(base+1:base+int(hi, 8)*m) = gval(o+1:o+int(hi, 8)*m)
      enddo
    enddo
    !$omp end taskloop
  end subroutine mf_store_cb

  !> Factor the matrix of hecMAT in the order of sym, in LDLt mode (lower part referenced) or
  !> in LU mode when the values are unsymmetric or the symmetric flag is off. ierr is 0 on
  !> success, the permuted DOF of a zero pivot (singular matrix), -1 when the block size of
  !> hecMAT does not match the structure, or -2 when the structure is not symmetric.
  subroutine hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(out) :: ierr
    type(mf_work), allocatable :: wrks(:)
    integer(kind=kint), allocatable :: left(:)
    integer(kind=kint) :: s, k, nd, nthr
    real(kind=kreal) :: amax, zero

    ierr = 0
    fct%factored = .false.
    do k = 1, sym%nnode
      if (sym%ndof(k) /= hecMAT%NDOF) then
        ierr = -1
        return
      endif
    enddo
    nd = hecMAT%NDOF
    amax = maxval(abs(hecMAT%D(1:hecMAT%NP*nd*nd)))
    if (hecMAT%NPL > 0) amax = max(amax, maxval(abs(hecMAT%AL(1:hecMAT%NPL*nd*nd))))
    if (hecMAT%NPU > 0) amax = max(amax, maxval(abs(hecMAT%AU(1:hecMAT%NPU*nd*nd))))
    if (.not. (fct%pivot_u > 0.0d0)) fct%pivot_u = MF_PIVOT_U
    if (.not. (fct%pivot_zero > 0.0d0)) fct%pivot_zero = MF_PIVOT_ZERO
    if (fct%blr) then
      if (.not. hecmw_mf_kernel_blr_available()) fct%blr = .false.
      if (.not. (fct%eps > 0.0d0)) fct%eps = MF_BLR_EPS
    endif
    zero = fct%pivot_zero * amax
    call mf_mirror(hecMAT, fct, ierr)
    if (ierr /= 0) return
    call mf_asymmetry(hecMAT, fct, amax)

    fct%n_pos = 0
    fct%n_neg = 0
    fct%n_2x2 = 0
    fct%n_swap = 0
    fct%n_delay = 0
    fct%max_growth = 0
    fct%factor_words_act = 0
    fct%stack_peak_act = 0
    fct%front_words_act = 0
    fct%live_cb = 0
    fct%live_front = 0
    fct%front_peak = 0
    fct%blr_words_fr = 0
    fct%blr_tiles = 0
    fct%blr_tiles_lr = 0
    fct%blr_rank_sum = 0
    fct%blr_rank_max = 0

    nthr = 1
    !$ nthr = omp_get_max_threads()
    allocate(wrks(0:nthr-1), left(fct%nsuper))
    do s = 1, fct%nsuper
      left(s) = fct%cptr(s+1) - fct%cptr(s)
    enddo
    !$omp parallel default(shared)
    !$omp single
    do s = 1, fct%nsuper
      ! spawn the true leaves only: a climbing task may drive left(s) of an inner node to
      ! zero while this loop is still running
      if (fct%cptr(s+1) == fct%cptr(s)) then
        !$omp task default(shared) firstprivate(s)
        call mf_super_task(hecMAT, sym, fct, wrks, left, s, zero, ierr)
        !$omp end task
      endif
    enddo
    !$omp end single
    !$omp end parallel
    deallocate(wrks, left)
    if (ierr /= 0) return
    fct%factored = .true.
  end subroutine hecmw_mf_numeric_factor

  !> Task body starting at leaf s: factor the front (skipped once an earlier task failed),
  !> release the parent, and continue with the parent in the same task when this was its last
  !> unfinished child. The dependences are counted explicitly (left) because the
  !> depend(iterator(...)) form does not order the tasks with gfortran 13, and the climb stays
  !> within one task so that no reference outlives its scope. The work space of the executing
  !> thread is safe to use because the task holds no task scheduling point.
  subroutine mf_super_task(hecMAT, sym, fct, wrks, left, s, zero, gierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    type(mf_work), intent(inout) :: wrks(0:)
    integer(kind=kint), intent(inout) :: left(:)
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(in) :: zero
    integer(kind=kint), intent(inout) :: gierr
    integer(kind=kint) :: ierr, cur, tid, p, n, ss, c, l, nr, nteam

    tid = 0
    nteam = 1
    !$ tid = omp_get_thread_num()
    !$ nteam = omp_get_num_threads()
    ss = s
    do
      !$omp atomic read
      cur = gierr
      if (cur == 0) then
        nr = 0
        c = fct%chead(ss)
        do while (c /= 0)
          nr = nr + fct%sn(c)%ncol - fct%sn(c)%npiv
          c = fct%cnext(c)
        enddo
        do l = sym%rptr(ss), sym%rptr(ss+1) - 1
          nr = nr + sym%ndof(sym%rlist(l))
        enddo
        if (nr >= MF_PAR_ROWS .and. nteam > 1) then
          ! a large front suspends at its task loops, so it may not borrow the thread's work
          ! space, which another task on this thread could then reuse
          block
            type(mf_work) :: lwrk
            call mf_super_factor(hecMAT, sym, fct, fct%sn(ss), lwrk, ss, zero, ierr)
          end block
        else
          call mf_super_factor(hecMAT, sym, fct, fct%sn(ss), wrks(tid), ss, zero, ierr)
        endif
        if (ierr /= 0) then
          !$omp atomic write
          gierr = ierr
        endif
      endif
      p = sym%sparent(ss)
      if (p == 0) return
      ! critical, not atomic: its flush semantics make the child's writes visible to the
      ! thread that continues with the parent
      !$omp critical (mf_tree)
      left(p) = left(p) - 1
      n = left(p)
      !$omp end critical (mf_tree)
      if (n /= 0) return
      ss = p
    enddo
  end subroutine mf_super_task

  !> Factor the replicated matrix with the supernodes distributed by map: every rank first
  !> factors its subtrees with the task-parallel code, then all ranks walk the upper fronts
  !> in ascending order, the owner of a front receiving the contribution blocks of the
  !> children owned elsewhere so that the extend-add consumes them exactly like local ones
  !> (fan-in). MPI calls stay on the master thread outside the parallel regions because
  !> hecmw initializes MPI without a threading level. Sends are isend so that a rank blocked
  !> in a receive has already issued the sends of its earlier fronts, which makes the
  !> ascending walk deadlock free. Tags encode 3*supernode+kind, relying on a tag space
  !> larger than the MPI minimum of 32767, as every mainstream MPI provides.
  subroutine hecmw_mf_numeric_factor_mpi(hecMAT, sym, map, fct, ierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(out) :: ierr
    type(mf_work), allocatable :: wrks(:)
    type(mf_sendbox), allocatable :: box(:)
    integer(kind=kint), allocatable :: left(:), reqs(:), stats(:,:)
    integer(kind=kint) :: s, k, nd, nthr, nsend, nreq, ib, iu, j, c, gierr, ierr2, iw(1)
    real(kind=kreal) :: amax, zero

    ierr = 0
    fct%factored = .false.
    do k = 1, sym%nnode
      if (sym%ndof(k) /= hecMAT%NDOF) then
        ierr = -1
        return
      endif
    enddo
    nd = hecMAT%NDOF
    amax = maxval(abs(hecMAT%D(1:hecMAT%NP*nd*nd)))
    if (hecMAT%NPL > 0) amax = max(amax, maxval(abs(hecMAT%AL(1:hecMAT%NPL*nd*nd))))
    if (hecMAT%NPU > 0) amax = max(amax, maxval(abs(hecMAT%AU(1:hecMAT%NPU*nd*nd))))
    if (.not. (fct%pivot_u > 0.0d0)) fct%pivot_u = MF_PIVOT_U
    if (.not. (fct%pivot_zero > 0.0d0)) fct%pivot_zero = MF_PIVOT_ZERO
    if (fct%blr) then
      if (.not. hecmw_mf_kernel_blr_available()) fct%blr = .false.
      if (.not. (fct%eps > 0.0d0)) fct%eps = MF_BLR_EPS
    endif
    zero = fct%pivot_zero * amax
    call mf_mirror(hecMAT, fct, ierr)
    if (ierr /= 0) return
    call mf_asymmetry(hecMAT, fct, amax)

    fct%n_pos = 0
    fct%n_neg = 0
    fct%n_2x2 = 0
    fct%n_swap = 0
    fct%n_delay = 0
    fct%max_growth = 0
    fct%factor_words_act = 0
    fct%stack_peak_act = 0
    fct%front_words_act = 0
    fct%live_cb = 0
    fct%live_front = 0
    fct%front_peak = 0
    fct%blr_words_fr = 0
    fct%blr_tiles = 0
    fct%blr_tiles_lr = 0
    fct%blr_rank_sum = 0
    fct%blr_rank_max = 0

    nthr = 1
    !$ nthr = omp_get_max_threads()
    allocate(wrks(0:nthr-1), left(fct%nsuper))
    do s = 1, fct%nsuper
      left(s) = fct%cptr(s+1) - fct%cptr(s)
    enddo
    gierr = 0
    !$omp parallel default(shared)
    !$omp single
    do s = 1, fct%nsuper
      if (fct%cptr(s+1) == fct%cptr(s) .and. map%owner(s) == map%myrank .and. .not. map%upper(s)) then
        !$omp task default(shared) firstprivate(s)
        call mf_super_task_mpi(hecMAT, sym, map, fct, wrks, left, s, zero, gierr)
        !$omp end task
      endif
    enddo
    !$omp end single
    !$omp end parallel
    deallocate(wrks, left)
    iw(1) = gierr
    call hecmw_allreduce_I_comm(iw, 1, hecmw_max, map%comm)
    gierr = iw(1)
    if (gierr /= 0) then
      ierr = gierr
      return
    endif

    nsend = 0
    do s = 1, fct%nsuper
      if (map%owner(s) == map%myrank .and. sym%sparent(s) /= 0) then
        if (map%owner(sym%sparent(s)) /= map%myrank) nsend = nsend + 1
      endif
    enddo
    allocate(box(max(nsend, 1)), reqs(3*max(nsend, 1)), stats(HECMW_STATUS_SIZE, 3*max(nsend, 1)))
    nreq = 0
    ib = 0
    ! the subtree roots feed the fan-in first, then the upper fronts as the walk factors them
    do s = 1, fct%nsuper
      if (map%owner(s) /= map%myrank .or. map%upper(s) .or. sym%sparent(s) == 0) cycle
      if (map%upper(sym%sparent(s)) .and. map%owner(sym%sparent(s)) /= map%myrank) call send_cb(s)
    enddo
    do iu = 1, map%nupper
      s = map%uplist(iu)
      if (map%owner(s) /= map%myrank) cycle
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        c = fct%clist(j)
        if (map%owner(c) /= map%myrank) call recv_cb(c, map%owner(c))
      enddo
      if (gierr == 0) then
        !$omp parallel default(shared)
        !$omp single
        block
          type(mf_work) :: lwrk
          call mf_super_factor(hecMAT, sym, fct, fct%sn(s), lwrk, s, zero, ierr2)
        end block
        !$omp end single
        !$omp end parallel
        if (ierr2 /= 0) gierr = ierr2
      endif
      if (sym%sparent(s) /= 0) then
        if (map%owner(sym%sparent(s)) /= map%myrank) call send_cb(s)
      endif
    enddo
    call hecmw_waitall(nreq, reqs, stats)
    ! the sent contribution blocks were consumed remotely; free the local copies
    do s = 1, fct%nsuper
      if (map%owner(s) == map%myrank .and. sym%sparent(s) /= 0) then
        if (map%owner(sym%sparent(s)) /= map%myrank .and. allocated(fct%sn(s)%cval)) then
          fct%live_cb = fct%live_cb - fct%sn(s)%cbsize
          deallocate(fct%sn(s)%cval)
        endif
      endif
    enddo
    deallocate(box, reqs, stats)
    iw(1) = gierr
    call hecmw_allreduce_I_comm(iw, 1, hecmw_max, map%comm)
    gierr = iw(1)
    ierr = gierr
    if (ierr /= 0) return
    fct%factored = .true.

  contains

    !> isend the contribution block of s0 to the owner of its parent: the header (with the
    !> local error state), the tile boundaries with the delayed column and row DOFs, and the
    !> block values; on an error only the header travels
    subroutine send_cb(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: dst, m, nt1, cndel

      dst = map%owner(sym%sparent(s0))
      ib = ib + 1
      allocate(box(ib)%hdr(5))
      box(ib)%hdr(1) = gierr
      if (gierr == 0) then
        box(ib)%hdr(2) = fct%sn(s0)%ncol
        box(ib)%hdr(3) = fct%sn(s0)%npiv
        box(ib)%hdr(4) = fct%sn(s0)%nt
        box(ib)%hdr(5) = fct%sn(s0)%ntc
      else
        box(ib)%hdr(2:5) = 0
      endif
      nreq = nreq + 1
      call hecmw_isend_int(box(ib)%hdr, 5, dst, 3*s0, map%comm, reqs(nreq))
      if (gierr /= 0) return
      nt1 = fct%sn(s0)%nt + 1
      cndel = fct%sn(s0)%ncol - fct%sn(s0)%npiv
      m = nt1 + 2*cndel
      allocate(box(ib)%tl(m))
      box(ib)%tl(1:nt1) = fct%sn(s0)%tbnd(1:nt1)
      box(ib)%tl(nt1+1:nt1+cndel) = fct%sn(s0)%fsdof(fct%sn(s0)%npiv+1:fct%sn(s0)%ncol)
      box(ib)%tl(nt1+cndel+1:nt1+2*cndel) = fct%sn(s0)%frow(fct%sn(s0)%npiv+1:fct%sn(s0)%ncol)
      nreq = nreq + 1
      call hecmw_isend_int(box(ib)%tl, m, dst, 3*s0+1, map%comm, reqs(nreq))
      if (fct%sn(s0)%cbsize > 0) then
        nreq = nreq + 1
        call hecmw_isend_r(fct%sn(s0)%cval, int(fct%sn(s0)%cbsize, kind=kint), dst, 3*s0+2, map%comm, reqs(nreq))
      endif
    end subroutine send_cb

    !> receive the contribution block of child c0 into its snode, in the shape a local
    !> factorization would have left: the extend-add and the forward solve of the parent
    !> then consume it unchanged
    subroutine recv_cb(c0, src)
      integer(kind=kint), intent(in) :: c0, src
      integer(kind=kint) :: hdr(5), stat(HECMW_STATUS_SIZE), m, nt1, cndel
      integer(kind=kint), allocatable :: tl(:)

      call hecmw_recv_int(hdr, 5, src, 3*c0, map%comm, stat)
      if (hdr(1) /= 0) then
        gierr = hdr(1)
        return
      endif
      fct%sn(c0)%ncol = hdr(2)
      fct%sn(c0)%npiv = hdr(3)
      fct%sn(c0)%nt = hdr(4)
      fct%sn(c0)%ntc = hdr(5)
      nt1 = hdr(4) + 1
      cndel = hdr(2) - hdr(3)
      call mf_grow_i(fct%sn(c0)%tbnd, nt1)
      call mf_grow_i(fct%sn(c0)%fsdof, max(hdr(2), 1))
      call mf_grow_i(fct%sn(c0)%frow, max(hdr(2), 1))
      m = nt1 + 2*cndel
      allocate(tl(m))
      call hecmw_recv_int(tl, m, src, 3*c0+1, map%comm, stat)
      fct%sn(c0)%tbnd(1:nt1) = tl(1:nt1)
      fct%sn(c0)%fsdof(hdr(3)+1:hdr(2)) = tl(nt1+1:nt1+cndel)
      fct%sn(c0)%frow(hdr(3)+1:hdr(2)) = tl(nt1+cndel+1:nt1+2*cndel)
      deallocate(tl)
      fct%sn(c0)%cbsize = mf_cb_words(fct%sn(c0))
      if (fct%lu) fct%sn(c0)%cbsize = 2*fct%sn(c0)%cbsize
      if (allocated(fct%sn(c0)%cval)) deallocate(fct%sn(c0)%cval)
      if (fct%sn(c0)%cbsize > 0) then
        allocate(fct%sn(c0)%cval(fct%sn(c0)%cbsize))
        call hecmw_recv_r(fct%sn(c0)%cval, int(fct%sn(c0)%cbsize, kind=kint), src, 3*c0+2, map%comm, stat)
        fct%live_cb = fct%live_cb + fct%sn(c0)%cbsize
        fct%stack_peak_act = max(fct%stack_peak_act, fct%live_cb)
      endif
    end subroutine recv_cb

  end subroutine hecmw_mf_numeric_factor_mpi

  !> mf_super_task limited to the subtrees of the executing rank: the climb stops below an
  !> upper front, which the sequential fan-in stage factors.
  subroutine mf_super_task_mpi(hecMAT, sym, map, fct, wrks, left, s, zero, gierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    type(mf_work), intent(inout) :: wrks(0:)
    integer(kind=kint), intent(inout) :: left(:)
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(in) :: zero
    integer(kind=kint), intent(inout) :: gierr
    integer(kind=kint) :: ierr, cur, tid, p, n, ss, c, l, nr, nteam

    tid = 0
    nteam = 1
    !$ tid = omp_get_thread_num()
    !$ nteam = omp_get_num_threads()
    ss = s
    do
      !$omp atomic read
      cur = gierr
      if (cur == 0) then
        nr = 0
        c = fct%chead(ss)
        do while (c /= 0)
          nr = nr + fct%sn(c)%ncol - fct%sn(c)%npiv
          c = fct%cnext(c)
        enddo
        do l = sym%rptr(ss), sym%rptr(ss+1) - 1
          nr = nr + sym%ndof(sym%rlist(l))
        enddo
        if (nr >= MF_PAR_ROWS .and. nteam > 1) then
          ! a large front suspends at its task loops, so it may not borrow the thread's work
          ! space, which another task on this thread could then reuse
          block
            type(mf_work) :: lwrk
            call mf_super_factor(hecMAT, sym, fct, fct%sn(ss), lwrk, ss, zero, ierr)
          end block
        else
          call mf_super_factor(hecMAT, sym, fct, fct%sn(ss), wrks(tid), ss, zero, ierr)
        endif
        if (ierr /= 0) then
          !$omp atomic write
          gierr = ierr
        endif
      endif
      p = sym%sparent(ss)
      if (p == 0) return
      if (map%upper(p)) return
      ! critical, not atomic: its flush semantics make the child's writes visible to the
      ! thread that continues with the parent
      !$omp critical (mf_tree)
      left(p) = left(p) - 1
      n = left(p)
      !$omp end critical (mf_tree)
      if (n /= 0) return
      ss = p
    enddo
  end subroutine mf_super_task_mpi

  !> Words of one face of the stored contribution block of a supernode, from its tile
  !> boundaries; mirrors the layout of mf_store_cb.
  function mf_cb_words(sn) result(w)
    implicit none
    type(mf_snode), intent(in) :: sn
    integer(kind=8) :: w
    integer(kind=kint) :: cndel, cnrow, j, prev, b

    cndel = sn%ncol - sn%npiv
    cnrow = cndel + sn%tbnd(sn%nt+1) - sn%ncol
    w = 0
    prev = 0
    if (cndel > 0) then
      w = w + int(cnrow, 8)*cndel
      prev = cndel
    endif
    do j = sn%ntc + 1, sn%nt
      b = cndel + sn%tbnd(j+1) - sn%ncol
      w = w + int(cnrow - prev, 8)*(b - prev)
      prev = b
    enddo
  end function mf_cb_words

  !> Assemble and factor the front of supernode s: scatter the matrix entries, extend-add the
  !> contribution blocks of the children (freed here), factor the fully summed part and store
  !> the factor panel, the position metadata and the contribution block in sn.
  subroutine mf_super_factor(hecMAT, sym, fct, sn, wrk, s, zero, ierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    type(mf_snode), intent(inout) :: sn
    type(mf_work), intent(inout) :: wrk
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(in) :: zero
    integer(kind=kint), intent(out) :: ierr
    integer(kind=kint) :: nown, nrow_nodes, ncol0, ndel, c, i, j, k, l, nd, a, b, m, r, cc
    integer(kind=kint) :: j0, ki, kk, coff, roff, cnown, cnb, cndel, cnrow, hi, cnt, cntc, kbeg, nkc
    integer(kind=kint) :: nswap, n2x2, npos, nneg, ntl, nlr, rmax
    integer(kind=8) :: o, base, base2, half, pw, fw, rsum, pwa, pwu
    integer(kind=8), allocatable :: pcb(:)
    real(kind=kreal) :: v
    logical :: par

    ierr = 0
    nd = hecMAT%NDOF
    nown = sym%sptr(s+1) - sym%sptr(s)
    nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
    ncol0 = fct%cdofptr(sym%sptr(s+1)) - fct%cdofptr(sym%sptr(s))
    ndel = 0
    c = fct%chead(s)
    do while (c /= 0)
      ndel = ndel + fct%sn(c)%ncol - fct%sn(c)%npiv
      c = fct%cnext(c)
    enddo
    !$omp critical (mf_stats)
    fct%max_growth = max(fct%max_growth, ndel)
    !$omp end critical (mf_stats)
    call mf_partition(sym, fct, s, ndel, wrk%g)
    call mf_grow_i(wrk%rowoff, nrow_nodes + 1)
    call mf_grow_i(wrk%pos, sym%nnode)
    if (allocated(wrk%ccoloff)) then
      if (size(wrk%ccoloff) < wrk%g%nt + 2) deallocate(wrk%ccoloff, wrk%ctb)
    endif
    if (.not. allocated(wrk%ccoloff)) allocate(wrk%ccoloff(wrk%g%nt+2), wrk%ctb(0:wrk%g%nt+1))
    sn%ncol = wrk%g%ncol
    sn%nt = wrk%g%nt
    sn%ntc = wrk%g%ntc
    call mf_grow_i(sn%tbnd, wrk%g%nt + 1)
    sn%tbnd(1:wrk%g%nt+1) = wrk%g%tb(0:wrk%g%nt)
    call mf_grow_i(sn%fsdof, wrk%g%ncol)
    call mf_grow_i(sn%frow, wrk%g%ncol)
    call mf_grow_i(sn%ptype, wrk%g%ncol)
    call mf_grow_r(sn%dsub, int(wrk%g%ncol, 8))
    call mf_grow_i(wrk%blk, wrk%g%ncol)
    call mf_grow_i(wrk%blkr, wrk%g%ncol)
    m = 0
    do k = sym%sptr(s), sym%sptr(s+1) - 1
      do i = fct%cdofptr(k), fct%cdofptr(k+1) - 1
        m = m + 1
        sn%fsdof(m) = i
        sn%frow(m) = i
        wrk%blk(m) = k
      enddo
    enddo
    c = fct%chead(s)
    do while (c /= 0)
      do j = fct%sn(c)%npiv + 1, fct%sn(c)%ncol
        m = m + 1
        sn%fsdof(m) = fct%sn(c)%fsdof(j)
        sn%frow(m) = fct%sn(c)%frow(j)
        wrk%blk(m) = 0
      enddo
      c = fct%cnext(c)
    enddo
    wrk%blkr(1:wrk%g%ncol) = wrk%blk(1:wrk%g%ncol)

    ! row offsets (0-based positions) of the rlist nodes: own nodes, then the contribution rows
    ! after the delayed DOFs
    wrk%rowoff(1) = 0
    do i = 1, nrow_nodes
      k = sym%rlist(sym%rptr(s)+i-1)
      wrk%rowoff(i+1) = wrk%rowoff(i) + sym%ndof(k)
      if (i == nown) wrk%rowoff(i+1) = wrk%rowoff(i+1) + ndel
      wrk%pos(k) = i
    enddo
    par = .false.
    !$ par = wrk%g%nrow >= MF_PAR_ROWS .and. omp_get_num_threads() > 1
    call mf_grow_r(wrk%fval, wrk%g%coloff(wrk%g%nt+1))
    fw = wrk%g%coloff(wrk%g%nt+1)
    if (fct%lu) fw = 2*fw
    !$omp critical (mf_stats)
    fct%front_words_act = max(fct%front_words_act, wrk%g%coloff(wrk%g%nt+1))
    fct%live_front = fct%live_front + fw
    fct%front_peak = max(fct%front_peak, fct%live_front)
    !$omp end critical (mf_stats)
    if (fct%lu) call mf_grow_r(wrk%fvalu, wrk%g%coloff(wrk%g%nt+1))
    !$omp taskloop default(shared) if(par)
    do j = 1, wrk%g%nt
      wrk%fval(wrk%g%coloff(j)+1:wrk%g%coloff(j+1)) = 0.0d0
      if (fct%lu) wrk%fvalu(wrk%g%coloff(j)+1:wrk%g%coloff(j+1)) = 0.0d0
    enddo
    !$omp end taskloop

    ! scatter the permuted matrix. LDLt: the lower part, every neighbor above the diagonal in
    ! the original numbering supplying the transposed block. LU: both parts, the block of a
    ! neighbor and its mirror going to the upper and the lower grid
    !$omp taskloop default(shared) private(k, j0, coff, base, base2, o, a, b, kk, ki, roff) &
    !$omp&  grainsize(8) if(par)
    do i = 1, nown
      k = sym%sptr(s) + i - 1
      j0 = sym%perm(k)
      coff = wrk%rowoff(i)
      base = int(j0-1, 8)*nd*nd
      do b = 1, nd
        do a = b, nd
          o = mf_idx(wrk%g, coff+a, coff+b)
          wrk%fval(o) = wrk%fval(o) + hecMAT%D(base + (a-1)*nd + b)
        enddo
        if (fct%lu) then
          do a = 1, b-1
            o = mf_idx(wrk%g, coff+b, coff+a)
            wrk%fvalu(o) = wrk%fvalu(o) + hecMAT%D(base + (a-1)*nd + b)
          enddo
        endif
      enddo
      do kk = hecMAT%indexL(j0-1)+1, hecMAT%indexL(j0)
        ki = sym%invp(hecMAT%itemL(kk))
        if (ki <= k) cycle
        roff = wrk%rowoff(wrk%pos(ki))
        base = int(kk-1, 8)*nd*nd
        if (.not. fct%lu) then
          do b = 1, nd
            do a = 1, nd
              o = mf_idx(wrk%g, roff+a, coff+b)
              wrk%fval(o) = wrk%fval(o) + hecMAT%AL(base + (b-1)*nd + a)
            enddo
          enddo
        else
          base2 = int(fct%mirror(kk)-1, 8)*nd*nd
          do a = 1, nd
            do b = 1, nd
              o = mf_idx(wrk%g, roff+b, coff+a)
              wrk%fvalu(o) = wrk%fvalu(o) + hecMAT%AL(base + (a-1)*nd + b)
              wrk%fval(o) = wrk%fval(o) + hecMAT%AU(base2 + (b-1)*nd + a)
            enddo
          enddo
        endif
      enddo
      do kk = hecMAT%indexU(j0-1)+1, hecMAT%indexU(j0)
        ki = sym%invp(hecMAT%itemU(kk))
        if (ki <= k) cycle
        roff = wrk%rowoff(wrk%pos(ki))
        base = int(kk-1, 8)*nd*nd
        if (.not. fct%lu) then
          do b = 1, nd
            do a = 1, nd
              o = mf_idx(wrk%g, roff+a, coff+b)
              wrk%fval(o) = wrk%fval(o) + hecMAT%AU(base + (b-1)*nd + a)
            enddo
          enddo
        else
          base2 = int(fct%mirroru(kk)-1, 8)*nd*nd
          do a = 1, nd
            do b = 1, nd
              o = mf_idx(wrk%g, roff+b, coff+a)
              wrk%fvalu(o) = wrk%fvalu(o) + hecMAT%AU(base + (a-1)*nd + b)
              wrk%fval(o) = wrk%fval(o) + hecMAT%AL(base2 + (b-1)*nd + a)
            enddo
          enddo
        endif
      enddo
    enddo
    !$omp end taskloop

    ! extend-add of the children; the rows of a child's contribution block are its delayed
    ! DOFs followed by its contribution rows, and the block is freed once consumed
    m = ncol0
    c = fct%chead(s)
    do while (c /= 0)
      cnown = sym%sptr(c+1) - sym%sptr(c)
      cndel = fct%sn(c)%ncol - fct%sn(c)%npiv
      cnt = fct%sn(c)%nt
      cntc = fct%sn(c)%ntc
      if (size(wrk%ccoloff) < cnt + 2) then
        deallocate(wrk%ccoloff, wrk%ctb)
        allocate(wrk%ccoloff(cnt+2), wrk%ctb(0:cnt+1))
      endif
      cnb = 0
      wrk%ctb(0) = 0
      if (cndel > 0) then
        cnb = 1
        wrk%ctb(1) = cndel
      endif
      do j = cntc + 1, cnt
        cnb = cnb + 1
        wrk%ctb(cnb) = cndel + fct%sn(c)%tbnd(j+1) - fct%sn(c)%ncol
      enddo
      cnrow = wrk%ctb(cnb)
      wrk%ccoloff(1) = 0
      do j = 1, cnb
        wrk%ccoloff(j+1) = wrk%ccoloff(j) + int(cnrow - wrk%ctb(j-1), 8)*(wrk%ctb(j) - wrk%ctb(j-1))
      enddo
      half = wrk%ccoloff(cnb+1)
      call mf_grow_i(wrk%cmapdof, cnrow)
      do r = 1, cndel
        wrk%cmapdof(r) = m + r
      enddo
      m = m + cndel
      r = cndel
      do l = sym%rptr(c) + cnown, sym%rptr(c+1) - 1
        k = sym%cmap(sym%cmap_ptr(c) + (l - sym%rptr(c) - cnown))
        do a = 1, sym%ndof(sym%rlist(l))
          r = r + 1
          wrk%cmapdof(r) = wrk%rowoff(k) + a
        enddo
      enddo
      !$omp taskloop default(shared) private(i, base, hi, cc, r, o, v) grainsize(1) if(par)
      do j = 1, cnb
        do i = j, cnb
          base = wrk%ccoloff(j) + int(wrk%ctb(i-1) - wrk%ctb(j-1), 8)*(wrk%ctb(j) - wrk%ctb(j-1))
          hi = wrk%ctb(i) - wrk%ctb(i-1)
          do cc = wrk%ctb(j-1) + 1, wrk%ctb(j)
            do r = max(cc, wrk%ctb(i-1) + 1), wrk%ctb(i)
              o = base + int(cc - 1 - wrk%ctb(j-1), 8)*hi + (r - 1 - wrk%ctb(i-1)) + 1
              v = fct%sn(c)%cval(o)
              if (.not. fct%lu) then
                o = mf_idx(wrk%g, max(wrk%cmapdof(r), wrk%cmapdof(cc)), min(wrk%cmapdof(r), wrk%cmapdof(cc)))
                wrk%fval(o) = wrk%fval(o) + v
              else
                call add_lu(v, wrk%cmapdof(r), wrk%cmapdof(cc))
                if (r > cc) call add_lu(fct%sn(c)%cval(o + half), wrk%cmapdof(cc), wrk%cmapdof(r))
              endif
            enddo
          enddo
        enddo
      enddo
      !$omp end taskloop
      !$omp critical (mf_stats)
      fct%live_cb = fct%live_cb - fct%sn(c)%cbsize
      !$omp end critical (mf_stats)
      deallocate(fct%sn(c)%cval)
      c = fct%cnext(c)
    enddo

    nswap = 0
    n2x2 = 0
    npos = 0
    nneg = 0
    ntl = 0
    nlr = 0
    rsum = 0
    rmax = 0
    if (fct%lu) then
      call mf_factor_front_lu(sn, wrk%g, wrk%fval, wrk%fvalu, wrk%pval, wrk%pvalu, wrk%wk, wrk%blk, wrk%blkr, &
        sym%sparent(s) == 0, fct%pivot_u, zero, fct%blr, fct%eps, wrk%bval, wrk%boff, wrk%brk, &
        wrk%bvalu, wrk%boffu, wrk%brku, nswap, ntl, nlr, rsum, rmax, ierr)
    else
      call mf_factor_front(sn, wrk%g, wrk%fval, wrk%pval, wrk%wval, wrk%wk, wrk%blk, &
        sym%sparent(s) == 0, fct%pivot_u, zero, fct%blr, fct%eps, wrk%bval, wrk%boff, wrk%brk, &
        nswap, n2x2, npos, nneg, ntl, nlr, rsum, rmax, ierr)
      sn%frow(1:wrk%g%ncol) = sn%fsdof(1:wrk%g%ncol)
    endif
    if (ierr /= 0) then
      !$omp critical (mf_stats)
      fct%live_front = fct%live_front - fw
      !$omp end critical (mf_stats)
      return
    endif
    !$omp critical (mf_stats)
    fct%n_swap = fct%n_swap + nswap
    fct%n_2x2 = fct%n_2x2 + n2x2
    fct%n_pos = fct%n_pos + npos
    fct%n_neg = fct%n_neg + nneg
    fct%n_delay = fct%n_delay + wrk%g%ncol - sn%npiv
    fct%blr_tiles = fct%blr_tiles + ntl
    fct%blr_tiles_lr = fct%blr_tiles_lr + nlr
    fct%blr_rank_sum = fct%blr_rank_sum + rsum
    fct%blr_rank_max = max(fct%blr_rank_max, rmax)
    !$omp end critical (mf_stats)

    ! the factor panel: the leading npiv columns of the tile grid, tile by tile
    pw = 0
    do k = 1, wrk%g%ntc
      m = min(wrk%g%tb(k), sn%npiv) - wrk%g%tb(k-1)
      if (m <= 0) exit
      pw = pw + int(wrk%g%nrow - wrk%g%tb(k-1), 8)*m
    enddo
    allocate(pcb(0:wrk%g%ntc))
    nkc = 0
    pcb(0) = 1
    do k = 1, wrk%g%ntc
      m = min(wrk%g%tb(k), sn%npiv) - wrk%g%tb(k-1)
      if (m <= 0) exit
      nkc = k
      pcb(k) = pcb(k-1) + int(wrk%g%nrow - wrk%g%tb(k-1), 8)*m
    enddo
    if (.not. fct%blr) then
      call mf_grow_r(sn%lval, pw)
      if (fct%lu) call mf_grow_r(sn%uval, pw)
      !$omp taskloop default(shared) private(m, i, hi, o, base) grainsize(1) if(par)
      do k = 1, nkc
        m = min(wrk%g%tb(k), sn%npiv) - wrk%g%tb(k-1)
        base = pcb(k-1)
        do i = k, wrk%g%nt
          hi = wrk%g%tb(i) - wrk%g%tb(i-1)
          o = mf_off(wrk%g, i, k)
          sn%lval(base:base+int(hi, 8)*m-1) = wrk%fval(o+1:o+int(hi, 8)*m)
          if (fct%lu) sn%uval(base:base+int(hi, 8)*m-1) = wrk%fvalu(o+1:o+int(hi, 8)*m)
          base = base + int(hi, 8)*m
        enddo
      enddo
      !$omp end taskloop
      !$omp critical (mf_stats)
      fct%factor_words_act = fct%factor_words_act + pw
      !$omp end critical (mf_stats)
    else
      call mf_store_blr(wrk%g, sn%npiv, nkc, wrk%fval, wrk%bval, wrk%boff, wrk%brk, &
        sn%brank, sn%bptr, sn%lval, pwa, par)
      if (fct%lu) then
        call mf_store_blr(wrk%g, sn%npiv, nkc, wrk%fvalu, wrk%bvalu, wrk%boffu, wrk%brku, &
          sn%branku, sn%bptru, sn%uval, pwu, par)
        pwa = pwa + pwu
      endif
      !$omp critical (mf_stats)
      fct%factor_words_act = fct%factor_words_act + pwa
      fct%blr_words_fr = fct%blr_words_fr + pw
      if (fct%lu) fct%blr_words_fr = fct%blr_words_fr + pw
      !$omp end critical (mf_stats)
    endif

    ! the contribution block: the delayed DOFs as one tile followed by the contribution tiles;
    ! in LU mode the upper grid follows the lower one
    cndel = wrk%g%ncol - sn%npiv
    cnb = 0
    wrk%ctb(0) = 0
    if (cndel > 0) then
      cnb = 1
      wrk%ctb(1) = cndel
    endif
    kbeg = cnb
    do j = wrk%g%ntc + 1, wrk%g%nt
      cnb = cnb + 1
      wrk%ctb(cnb) = cndel + wrk%g%tb(j) - wrk%g%ncol
    enddo
    cnrow = wrk%ctb(cnb)
    wrk%ccoloff(1) = 0
    do j = 1, cnb
      wrk%ccoloff(j+1) = wrk%ccoloff(j) + int(cnrow - wrk%ctb(j-1), 8)*(wrk%ctb(j) - wrk%ctb(j-1))
    enddo
    half = wrk%ccoloff(cnb+1)
    sn%cbsize = half
    if (fct%lu) sn%cbsize = 2*half
    if (allocated(sn%cval)) deallocate(sn%cval)
    if (sn%cbsize > 0) then
      allocate(sn%cval(sn%cbsize))
      call mf_store_cb(wrk%g, wrk%fval, sn%npiv, cndel, cnb, kbeg, wrk%ctb, wrk%ccoloff, sn%cval, 0_8, par)
      if (fct%lu) call mf_store_cb(wrk%g, wrk%fvalu, sn%npiv, cndel, cnb, kbeg, wrk%ctb, wrk%ccoloff, sn%cval, half, par)
      !$omp critical (mf_stats)
      fct%live_cb = fct%live_cb + sn%cbsize
      fct%stack_peak_act = max(fct%stack_peak_act, fct%live_cb)
      !$omp end critical (mf_stats)
    endif
    !$omp critical (mf_stats)
    fct%live_front = fct%live_front - fw
    !$omp end critical (mf_stats)

  contains

    !> add v to the entry (pr,pc) of the LU front
    subroutine add_lu(v, pr, pc)
      real(kind=kreal), intent(in) :: v
      integer(kind=kint), intent(in) :: pr, pc
      integer(kind=8) :: o

      if (pr >= pc) then
        o = mf_idx(wrk%g, pr, pc)
        wrk%fval(o) = wrk%fval(o) + v
      else
        o = mf_idx(wrk%g, pc, pr)
        wrk%fvalu(o) = wrk%fvalu(o) + v
      endif
    end subroutine add_lu

  end subroutine mf_super_factor

  !> BLR panel store of one grid: build the tile offsets and ranks of the stored panel from
  !> the compression results, pack the full rank tiles from the grid and U, V of the
  !> compressed tiles from the compression slots; pwa returns the stored words.
  subroutine mf_store_blr(g, npiv, nkc, gval, bval, boff, brk, brank, bptr, dst, pwa, par)
    implicit none
    type(mf_grid), intent(in) :: g
    integer(kind=kint), intent(in) :: npiv, nkc
    real(kind=kreal), intent(in) :: gval(:), bval(:)
    integer(kind=8), intent(in) :: boff(:)
    integer(kind=kint), intent(in) :: brk(:)
    integer(kind=kint), allocatable, intent(inout) :: brank(:)
    integer(kind=8), allocatable, intent(inout) :: bptr(:)
    real(kind=kreal), allocatable, intent(inout) :: dst(:)
    integer(kind=8), intent(out) :: pwa
    logical, intent(in) :: par
    integer(kind=kint) :: ntile, k, i, m, hi, r, idx
    integer(kind=8) :: off, base, ob, o

    ntile = 0
    do k = 1, nkc
      ntile = ntile + g%nt - k + 1
    enddo
    call mf_grow_i(brank, ntile)
    call mf_grow_i8(bptr, ntile + 1)
    off = 0
    do k = 1, nkc
      m = min(g%tb(k), npiv) - g%tb(k-1)
      do i = k, g%nt
        idx = mf_bidx(g%nt, k, i)
        hi = g%tb(i) - g%tb(i-1)
        r = -1
        if (i > g%ntc) r = brk(idx)
        brank(idx) = r
        bptr(idx) = off
        if (r < 0) then
          off = off + int(hi, 8)*m
        else
          off = off + int(r, 8)*(hi + m)
        endif
      enddo
    enddo
    bptr(ntile+1) = off
    pwa = off
    call mf_grow_r(dst, pwa)
    !$omp taskloop default(shared) private(m, i, idx, hi, r, base, ob, o) grainsize(1) if(par)
    do k = 1, nkc
      m = min(g%tb(k), npiv) - g%tb(k-1)
      do i = k, g%nt
        idx = mf_bidx(g%nt, k, i)
        hi = g%tb(i) - g%tb(i-1)
        r = brank(idx)
        base = bptr(idx)
        if (r < 0) then
          o = mf_off(g, i, k)
          dst(base+1:base+int(hi, 8)*m) = gval(o+1:o+int(hi, 8)*m)
        else if (r > 0) then
          ob = boff(idx)
          dst(base+1:base+int(hi, 8)*r) = bval(ob+1:ob+int(hi, 8)*r)
          dst(base+int(hi, 8)*r+1:base+int(hi, 8)*r+int(m, 8)*r) = &
            bval(ob+int(hi, 8)*m+1:ob+int(hi, 8)*m+int(m, 8)*r)
        endif
      enddo
    enddo
    !$omp end taskloop
  end subroutine mf_store_blr

  !> Factor the fully summed part of the assembled front of supernode s, tile column by tile
  !> column. A tile column is first factored without pivoting in the dense work panel and
  !> accepted when the threshold holds; otherwise the panel is refilled and factored with
  !> pivot search, its delayed columns are exchanged with the last unprocessed fully summed
  !> columns and the tile column is refilled until it is complete or no column is left. The
  !> pivots of a tile column then update the trailing tiles. At a root the remaining delayed
  !> columns are factored in a final panel where any column may serve as partner.
  subroutine mf_factor_front(sn, g, fval, pval, wval, wk, blk, isroot, u, zero, blr, eps, &
      bval, boff, brk, nswap, n2x2, npos, nneg, ntl, nlr, rsum, rmax, ierr)
    implicit none
    type(mf_snode), intent(inout) :: sn
    type(mf_grid), intent(in) :: g
    real(kind=kreal), allocatable, intent(inout) :: fval(:)
    real(kind=kreal), allocatable, intent(inout) :: pval(:), wval(:), wk(:)
    integer(kind=kint), intent(inout) :: blk(:)
    logical, intent(in) :: isroot, blr
    real(kind=kreal), intent(in) :: u, zero, eps
    real(kind=kreal), allocatable, intent(inout) :: bval(:)
    integer(kind=8), allocatable, intent(inout) :: boff(:)
    integer(kind=kint), allocatable, intent(inout) :: brk(:)
    integer(kind=kint), intent(out) :: nswap, n2x2, npos, nneg, ntl, nlr, rmax, ierr
    integer(kind=8), intent(out) :: rsum
    integer(kind=kint), allocatable :: perm(:), itmp(:), ip(:), jp(:)
    integer(kind=kint) :: npiv, nfs, k, p0, pa, pb, w, m, jj, np, nsw, n22, info, x, i, j, hi, wj, pk, pt
    integer(kind=kint) :: np2, l2, mn, r, ra, rb
    integer(kind=8) :: oik, ojk, oij, okk, btop, ob, oa, oav, ob2, obv, iw
    real(kind=kreal) :: d11, d21, d22, det
    logical :: par

    ierr = 0
    nswap = 0
    n2x2 = 0
    npos = 0
    nneg = 0
    ntl = 0
    nlr = 0
    rsum = 0
    rmax = 0
    btop = 0
    if (blr) then
      call mf_grow_i(brk, mf_bidx(g%nt, g%ntc, g%nt))
      call mf_grow_i8(boff, mf_bidx(g%nt, g%ntc, g%nt))
      call mf_grow_r(bval, 1_8)
      brk(1:mf_bidx(g%nt, g%ntc, g%nt)) = -1
    endif
    par = .false.
    !$ par = g%nrow >= MF_PAR_ROWS .and. omp_get_num_threads() > 1
    sn%ptype(1:g%ncol) = 0
    sn%dsub(1:g%ncol) = 0.0d0
    npiv = 0
    nfs = g%ncol
    allocate(perm(g%ncol), itmp(g%ncol), ip(g%nt*(g%nt+1)/2), jp(g%nt*(g%nt+1)/2))
    call mf_grow_r(pval, int(g%nrow, 8))
    call mf_grow_r(wk, int(2*g%nrow, 8))
    do k = 1, g%ntc
      p0 = g%tb(k-1)
      if (p0 >= nfs) exit
      pa = p0 + 1
      do while (pa <= min(g%tb(k), nfs))
        pb = min(g%tb(k), nfs)
        w = pb - pa + 1
        m = g%nrow - pa + 1
        call mf_grow_r(pval, int(m, 8)*w)
        call mf_grow_r(wk, int(2*m, 8))
        call fill_panel(npiv > p0)
        call hecmw_mf_kernel_panel_nopiv(m, w, m, pval, u, zero, info)
        if (info == 0) then
          sn%ptype(pa:pb) = 1
          do jj = 1, w
            call mf_col_copy(g, fval, pa+jj-1, pa+jj-1, pval(int(jj-1, 8)*m + jj), .true.)
          enddo
          npiv = npiv + w
          pa = pb + 1
          cycle
        endif
        call fill_panel(npiv > p0)
        call hecmw_mf_kernel_panel_piv(m, w, m, pval, blk(pa:pb), u, MF_PIVOT_ALPHA, zero, &
          .true., wk, np, perm, sn%ptype(pa:pb), sn%dsub(pa:pb), nsw, n22, info)
        call permute_rows(pa, w)
        do jj = 1, np
          call mf_col_copy(g, fval, pa+jj-1, pa+jj-1, pval(int(jj-1, 8)*m + jj), .true.)
        enddo
        npiv = npiv + np
        nswap = nswap + nsw
        n2x2 = n2x2 + n22
        x = pa + np
        do while (x <= pb)
          if (nfs > pb) then
            call mf_swap(g, fval, x, nfs)
            call swap_pos(x, nfs)
            nfs = nfs - 1
            x = x + 1
          else
            nfs = x - 1
            exit
          endif
        enddo
        pa = pa + np
      enddo

      ! the pivots of tile column k update the trailing tiles and the delayed columns left in
      ! tile column k; a parallel front builds the scaled panel L*D for all trailing row tiles
      ! first so that the tile updates are independent
      pk = npiv - p0
      if (blr .and. pk > 0) then
        ! compress the contribution-row tiles of the column before the updates use them; the
        ! fully summed rows stay full rank because later pivot exchanges still permute them.
        ! The grid keeps the full rank values (the delayed column updates read them) and the
        ! slots are reserved up front so that the layout is schedule independent
        do i = g%ntc+1, g%nt
          hi = g%tb(i) - g%tb(i-1)
          boff(mf_bidx(g%nt, k, i)) = btop
          btop = btop + int(hi, 8)*pk + 2_8*int(pk, 8)*min(hi, pk)
        enddo
        call mf_grow_r(bval, btop)
        !$omp taskloop default(shared) private(hi, ob, oik, r) grainsize(1) if(par)
        do i = g%ntc+1, g%nt
          hi = g%tb(i) - g%tb(i-1)
          ob = boff(mf_bidx(g%nt, k, i))
          oik = mf_off(g, i, k)
          bval(ob+1:ob+int(hi, 8)*pk) = fval(oik+1:oik+int(hi, 8)*pk)
          call hecmw_mf_kernel_compress(hi, pk, hi, bval(ob+1), eps, pk, bval(ob+int(hi, 8)*pk+1), r)
          brk(mf_bidx(g%nt, k, i)) = r
        enddo
        !$omp end taskloop
        do i = g%ntc+1, g%nt
          ntl = ntl + 1
          r = brk(mf_bidx(g%nt, k, i))
          if (r >= 0) then
            nlr = nlr + 1
            rsum = rsum + r
            rmax = max(rmax, r)
          endif
        enddo
      endif
      if (pk > 0) then
        okk = mf_off(g, k, k)
        if (blr) then
          ! scaled panels: L*D of the full rank row tiles into wval, D*V of the compressed
          ! ones into their slot, then the (i,j) tile updates in full/low rank combinations
          call mf_grow_r(wval, int(g%nrow - g%tb(k), 8)*pk)
          !$omp taskloop default(shared) private(hi, mn, ob, oik, r) grainsize(1) if(par)
          do i = k+1, g%nt
            hi = g%tb(i) - g%tb(i-1)
            oik = mf_off(g, i, k)
            r = -1
            if (i > g%ntc) r = brk(mf_bidx(g%nt, k, i))
            if (r < 0) then
              call hecmw_mf_kernel_scale(hi, pk, g%tb(k) - g%tb(k-1), fval(okk+1), sn%ptype(p0+1:npiv), &
                sn%dsub(p0+1:npiv), hi, fval(oik+1), wval(int(g%tb(i-1) - g%tb(k), 8)*pk + 1))
            else if (r > 0) then
              mn = min(hi, pk)
              ob = boff(mf_bidx(g%nt, k, i))
              call hecmw_mf_kernel_scale_rows(pk, r, g%tb(k) - g%tb(k-1), fval(okk+1), sn%ptype(p0+1:npiv), &
                sn%dsub(p0+1:npiv), pk, bval(ob+int(hi, 8)*pk+1), pk, bval(ob+int(hi, 8)*pk+int(pk, 8)*mn+1))
            endif
          enddo
          !$omp end taskloop
          np2 = 0
          do i = k+1, g%nt
            do j = k+1, i
              np2 = np2 + 1
              ip(np2) = i
              jp(np2) = j
            enddo
          enddo
          !$omp taskloop default(shared) private(i, j, hi, wj, ra, rb, oa, oav, ob2, obv, oij, iw) &
          !$omp&  grainsize(1) if(par)
          do l2 = 1, np2
            i = ip(l2)
            j = jp(l2)
            hi = g%tb(i) - g%tb(i-1)
            wj = g%tb(j) - g%tb(j-1)
            oij = mf_off(g, i, j)
            iw = int(g%tb(i-1) - g%tb(k), 8)*pk
            ra = -1
            oa = 0
            oav = 0
            if (i > g%ntc) ra = brk(mf_bidx(g%nt, k, i))
            if (ra >= 0) then
              oa = boff(mf_bidx(g%nt, k, i))
              oav = oa + int(hi, 8)*pk + int(pk, 8)*min(hi, pk)
            endif
            rb = -1
            ob2 = 0
            obv = 0
            if (j > g%ntc) rb = brk(mf_bidx(g%nt, k, j))
            if (rb >= 0) then
              ob2 = boff(mf_bidx(g%nt, k, j))
              obv = ob2 + int(wj, 8)*pk
            endif
            call mf_update_ab(hi, wj, pk, ra, wval(iw+1), bval(oa+1), bval(oav+1), &
              rb, fval(mf_off(g, j, k)+1), bval(ob2+1), bval(obv+1), fval(oij+1))
          enddo
          !$omp end taskloop
        else if (par) then
          call mf_grow_r(wval, int(g%nrow - g%tb(k), 8)*pk)
          !$omp taskloop default(shared) private(hi, oik) grainsize(1)
          do i = k+1, g%nt
            hi = g%tb(i) - g%tb(i-1)
            oik = mf_off(g, i, k)
            call hecmw_mf_kernel_scale(hi, pk, g%tb(k) - g%tb(k-1), fval(okk+1), sn%ptype(p0+1:npiv), &
              sn%dsub(p0+1:npiv), hi, fval(oik+1), wval(int(g%tb(i-1) - g%tb(k), 8)*pk + 1))
          enddo
          !$omp end taskloop
          np2 = 0
          do i = k+1, g%nt
            do j = k+1, i
              np2 = np2 + 1
              ip(np2) = i
              jp(np2) = j
            enddo
          enddo
          !$omp taskloop default(shared) private(i, j, hi, wj, ojk, oij) grainsize(1)
          do l2 = 1, np2
            i = ip(l2)
            j = jp(l2)
            hi = g%tb(i) - g%tb(i-1)
            wj = g%tb(j) - g%tb(j-1)
            ojk = mf_off(g, j, k)
            oij = mf_off(g, i, j)
            call hecmw_mf_kernel_gemm(hi, wj, pk, hi, wval(int(g%tb(i-1) - g%tb(k), 8)*pk + 1), wj, &
              fval(ojk+1), hi, fval(oij+1))
          enddo
          !$omp end taskloop
        else
          do i = k+1, g%nt
            hi = g%tb(i) - g%tb(i-1)
            oik = mf_off(g, i, k)
            call mf_grow_r(wval, int(hi, 8)*pk)
            call hecmw_mf_kernel_scale(hi, pk, g%tb(k) - g%tb(k-1), fval(okk+1), sn%ptype(p0+1:npiv), &
              sn%dsub(p0+1:npiv), hi, fval(oik+1), wval)
            do j = k+1, i
              wj = g%tb(j) - g%tb(j-1)
              ojk = mf_off(g, j, k)
              oij = mf_off(g, i, j)
              call hecmw_mf_kernel_gemm(hi, wj, pk, hi, wval, wj, fval(ojk+1), hi, fval(oij+1))
            enddo
          enddo
        endif
        do x = npiv + 1, g%tb(k)
          call mf_col_copy(g, fval, x, x, pval, .false.)
          call mf_col_update(g, fval, sn%ptype, sn%dsub, x, p0+1, npiv, pval, wk)
          call mf_col_copy(g, fval, x, x, pval, .true.)
        enddo
      endif
      if (nfs <= g%tb(k)) exit
    enddo

    if (isroot .and. npiv < g%ncol) then
      pa = npiv + 1
      pb = g%ncol
      w = pb - pa + 1
      m = g%nrow - pa + 1
      call mf_grow_r(pval, int(m, 8)*w)
      call mf_grow_r(wk, int(2*m, 8))
      call fill_panel(.false.)
      call hecmw_mf_kernel_panel_piv(m, w, m, pval, blk(pa:pb), u, MF_PIVOT_ALPHA, zero, &
        .false., wk, np, perm, sn%ptype(pa:pb), sn%dsub(pa:pb), nsw, n22, info)
      call permute_rows(pa, w)
      nswap = nswap + nsw
      n2x2 = n2x2 + n22
      if (info /= 0) then
        ierr = sn%fsdof(pa + info - 1)
        deallocate(perm, itmp, ip, jp)
        return
      endif
      do jj = 1, w
        call mf_col_copy(g, fval, pa+jj-1, pa+jj-1, pval(int(jj-1, 8)*m + jj), .true.)
      enddo
      npiv = g%ncol
    endif
    sn%npiv = npiv

    x = 1
    do while (x <= npiv)
      pt = sn%ptype(x)
      if (pt == 1) then
        if (fval(mf_idx(g, x, x)) > 0.0d0) then
          npos = npos + 1
        else
          nneg = nneg + 1
        endif
        x = x + 1
      else
        d11 = fval(mf_idx(g, x, x))
        d21 = sn%dsub(x)
        d22 = fval(mf_idx(g, x+1, x+1))
        det = d11*d22 - d21*d21
        if (det < 0.0d0) then
          npos = npos + 1
          nneg = nneg + 1
        else if (d11 + d22 > 0.0d0) then
          npos = npos + 2
        else
          nneg = nneg + 2
        endif
        x = x + 2
      endif
    enddo
    deallocate(perm, itmp, ip, jp)

  contains

    !> pval(m,w) <- the columns pa..pb of the front (rows pa..nrow); with upd they receive the
    !> pivots p0+1..npiv of the current tile column, which the front does not carry yet
    subroutine fill_panel(upd)
      logical, intent(in) :: upd
      integer(kind=kint) :: jj

      do jj = 1, w
        call mf_col_copy(g, fval, pa+jj-1, pa+jj-1, pval(int(jj-1, 8)*m + jj), .false.)
        if (upd) &
          call mf_col_update(g, fval, sn%ptype, sn%dsub, pa+jj-1, p0+1, npiv, &
            pval(int(jj-1, 8)*m + jj), wk)
      enddo
    end subroutine fill_panel

    !> apply the panel permutation perm(1:w) to the positions pa..pa+w-1 of the front as
    !> symmetric exchanges (the kernel permuted the panel and the blocks, not the front)
    subroutine permute_rows(pa, w)
      integer(kind=kint), intent(in) :: pa, w
      integer(kind=kint) :: jj, t, q

      do jj = 1, w
        itmp(jj) = jj
      enddo
      do jj = 1, w
        if (itmp(jj) == perm(jj)) cycle
        do t = jj + 1, w
          if (itmp(t) == perm(jj)) exit
        enddo
        call mf_swap(g, fval, pa+jj-1, pa+t-1)
        q = sn%fsdof(pa+jj-1)
        sn%fsdof(pa+jj-1) = sn%fsdof(pa+t-1)
        sn%fsdof(pa+t-1) = q
        itmp(t) = itmp(jj)
        itmp(jj) = perm(jj)
      enddo
    end subroutine permute_rows

    subroutine swap_pos(x, y)
      integer(kind=kint), intent(in) :: x, y
      integer(kind=kint) :: t

      t = sn%fsdof(x)
      sn%fsdof(x) = sn%fsdof(y)
      sn%fsdof(y) = t
      t = blk(x)
      blk(x) = blk(y)
      blk(y) = t
    end subroutine swap_pos

  end subroutine mf_factor_front

  !> LU counterpart of mf_factor_front on the grids fval, fvalu and the panel pair pval, pvalu:
  !> pivots are chosen by threshold partial pivoting among the rows of the panel, the row
  !> exchanges are applied to the front and to frow, and the row of U of a pivot (a column of
  !> fvalu) is complete at write back, so that a column left behind receives the update of the
  !> tile column's pivots as two vector operations, one on each grid.
  subroutine mf_factor_front_lu(sn, g, fval, fvalu, pval, pvalu, wk, blk, blkr, isroot, u, zero, &
      blr, eps, bval, boff, brk, bvalu, boffu, brku, nswap, ntl, nlr, rsum, rmax, ierr)
    implicit none
    type(mf_snode), intent(inout) :: sn
    type(mf_grid), intent(in) :: g
    real(kind=kreal), allocatable, intent(inout) :: fval(:), fvalu(:)
    real(kind=kreal), allocatable, intent(inout) :: pval(:), pvalu(:), wk(:)
    integer(kind=kint), intent(inout) :: blk(:), blkr(:)
    logical, intent(in) :: isroot, blr
    real(kind=kreal), intent(in) :: u, zero, eps
    real(kind=kreal), allocatable, intent(inout) :: bval(:), bvalu(:)
    integer(kind=8), allocatable, intent(inout) :: boff(:), boffu(:)
    integer(kind=kint), allocatable, intent(inout) :: brk(:), brku(:)
    integer(kind=kint), intent(out) :: nswap, ntl, nlr, rmax, ierr
    integer(kind=8), intent(out) :: rsum
    integer(kind=kint), allocatable :: permc(:), permr(:), itmp(:), ip(:), jp(:)
    integer(kind=kint) :: npiv, nfs, k, p0, pa, pb, w, m, np, nsw, info, x, i, j, hi, wj, pk
    integer(kind=kint) :: np2, l2, r, rla, rlb, rua, rub
    integer(kind=8) :: oik, ojk, oij, btop, btopu, ob, ola, olav, olb, olbv, oua, ouav, oub, oubv
    logical :: par

    ierr = 0
    nswap = 0
    ntl = 0
    nlr = 0
    rsum = 0
    rmax = 0
    btop = 0
    btopu = 0
    if (blr) then
      call mf_grow_i(brk, mf_bidx(g%nt, g%ntc, g%nt))
      call mf_grow_i(brku, mf_bidx(g%nt, g%ntc, g%nt))
      call mf_grow_i8(boff, mf_bidx(g%nt, g%ntc, g%nt))
      call mf_grow_i8(boffu, mf_bidx(g%nt, g%ntc, g%nt))
      call mf_grow_r(bval, 1_8)
      call mf_grow_r(bvalu, 1_8)
      brk(1:mf_bidx(g%nt, g%ntc, g%nt)) = -1
      brku(1:mf_bidx(g%nt, g%ntc, g%nt)) = -1
    endif
    par = .false.
    !$ par = g%nrow >= MF_PAR_ROWS .and. omp_get_num_threads() > 1
    sn%ptype(1:g%ncol) = 1
    sn%dsub(1:g%ncol) = 0.0d0
    npiv = 0
    nfs = g%ncol
    allocate(permc(g%ncol), permr(g%ncol), itmp(g%ncol), ip(g%nt*(g%nt+1)/2), jp(g%nt*(g%nt+1)/2))
    call mf_grow_r(pval, int(g%nrow, 8))
    call mf_grow_r(pvalu, int(g%nrow, 8))
    call mf_grow_r(wk, int(2*g%nrow, 8))
    do k = 1, g%ntc
      p0 = g%tb(k-1)
      if (p0 >= nfs) exit
      pa = p0 + 1
      do while (pa <= min(g%tb(k), nfs))
        pb = min(g%tb(k), nfs)
        w = pb - pa + 1
        m = g%nrow - pa + 1
        call mf_grow_r(pval, int(m, 8)*w)
        call mf_grow_r(pvalu, int(m, 8)*w)
        call fill_panel(npiv > p0)
        call hecmw_mf_kernel_panel_lu_nopiv(m, w, m, pval, m, pvalu, u, zero, info)
        if (info == 0) then
          call write_back(w)
          npiv = npiv + w
          pa = pb + 1
          cycle
        endif
        call fill_panel(npiv > p0)
        call hecmw_mf_kernel_panel_lu_piv(m, w, m, pval, m, pvalu, blk(pa:pb), blkr(pa:pb), u, &
          zero, .true., np, permc, permr, nsw, info)
        call permute_front(pa, w)
        call write_back(np)
        npiv = npiv + np
        nswap = nswap + nsw
        x = pa + np
        do while (x <= pb)
          if (nfs > pb) then
            call mf_swap_lu(g, fval, fvalu, x, nfs)
            call swap_pos(x, nfs)
            nfs = nfs - 1
            x = x + 1
          else
            nfs = x - 1
            exit
          endif
        enddo
        pa = pa + np
      enddo

      ! the pivots of tile column k update the trailing tiles of both grids and the delayed
      ! columns left in tile column k
      pk = npiv - p0
      if (blr .and. pk > 0) then
        ! compress the contribution-row tiles of both grids, as in the LDLt mode
        do i = g%ntc+1, g%nt
          hi = g%tb(i) - g%tb(i-1)
          boff(mf_bidx(g%nt, k, i)) = btop
          btop = btop + int(hi, 8)*pk + int(pk, 8)*min(hi, pk)
          boffu(mf_bidx(g%nt, k, i)) = btopu
          btopu = btopu + int(hi, 8)*pk + int(pk, 8)*min(hi, pk)
        enddo
        call mf_grow_r(bval, btop)
        call mf_grow_r(bvalu, btopu)
        !$omp taskloop default(shared) private(hi, ob, oik, r) grainsize(1) if(par)
        do i = g%ntc+1, g%nt
          hi = g%tb(i) - g%tb(i-1)
          oik = mf_off(g, i, k)
          ob = boff(mf_bidx(g%nt, k, i))
          bval(ob+1:ob+int(hi, 8)*pk) = fval(oik+1:oik+int(hi, 8)*pk)
          call hecmw_mf_kernel_compress(hi, pk, hi, bval(ob+1), eps, pk, bval(ob+int(hi, 8)*pk+1), r)
          brk(mf_bidx(g%nt, k, i)) = r
          ob = boffu(mf_bidx(g%nt, k, i))
          bvalu(ob+1:ob+int(hi, 8)*pk) = fvalu(oik+1:oik+int(hi, 8)*pk)
          call hecmw_mf_kernel_compress(hi, pk, hi, bvalu(ob+1), eps, pk, bvalu(ob+int(hi, 8)*pk+1), r)
          brku(mf_bidx(g%nt, k, i)) = r
        enddo
        !$omp end taskloop
        do i = g%ntc+1, g%nt
          ntl = ntl + 2
          r = brk(mf_bidx(g%nt, k, i))
          if (r >= 0) then
            nlr = nlr + 1
            rsum = rsum + r
            rmax = max(rmax, r)
          endif
          r = brku(mf_bidx(g%nt, k, i))
          if (r >= 0) then
            nlr = nlr + 1
            rsum = rsum + r
            rmax = max(rmax, r)
          endif
        enddo
      endif
      if (pk > 0) then
        np2 = 0
        do i = k+1, g%nt
          do j = k+1, i
            np2 = np2 + 1
            ip(np2) = i
            jp(np2) = j
          enddo
        enddo
        if (blr) then
          !$omp taskloop default(shared) private(i, j, hi, wj, oik, ojk, oij, rla, rlb, rua, rub, &
          !$omp&  ola, olav, olb, olbv, oua, ouav, oub, oubv) grainsize(1) if(par)
          do l2 = 1, np2
            i = ip(l2)
            j = jp(l2)
            hi = g%tb(i) - g%tb(i-1)
            wj = g%tb(j) - g%tb(j-1)
            oik = mf_off(g, i, k)
            ojk = mf_off(g, j, k)
            oij = mf_off(g, i, j)
            rla = -1
            ola = 0
            olav = 0
            rua = -1
            oua = 0
            ouav = 0
            if (i > g%ntc) then
              rla = brk(mf_bidx(g%nt, k, i))
              rua = brku(mf_bidx(g%nt, k, i))
              if (rla >= 0) then
                ola = boff(mf_bidx(g%nt, k, i))
                olav = ola + int(hi, 8)*pk
              endif
              if (rua >= 0) then
                oua = boffu(mf_bidx(g%nt, k, i))
                ouav = oua + int(hi, 8)*pk
              endif
            endif
            rlb = -1
            olb = 0
            olbv = 0
            rub = -1
            oub = 0
            oubv = 0
            if (j > g%ntc) then
              rlb = brk(mf_bidx(g%nt, k, j))
              rub = brku(mf_bidx(g%nt, k, j))
              if (rlb >= 0) then
                olb = boff(mf_bidx(g%nt, k, j))
                olbv = olb + int(wj, 8)*pk
              endif
              if (rub >= 0) then
                oub = boffu(mf_bidx(g%nt, k, j))
                oubv = oub + int(wj, 8)*pk
              endif
            endif
            call mf_update_ab(hi, wj, pk, rla, fval(oik+1), bval(ola+1), bval(olav+1), &
              rub, fvalu(ojk+1), bvalu(oub+1), bvalu(oubv+1), fval(oij+1))
            call mf_update_ab(hi, wj, pk, rua, fvalu(oik+1), bvalu(oua+1), bvalu(ouav+1), &
              rlb, fval(ojk+1), bval(olb+1), bval(olbv+1), fvalu(oij+1))
          enddo
          !$omp end taskloop
        else
          !$omp taskloop default(shared) private(i, j, hi, wj, oik, ojk, oij) grainsize(1) if(par)
          do l2 = 1, np2
            i = ip(l2)
            j = jp(l2)
            hi = g%tb(i) - g%tb(i-1)
            wj = g%tb(j) - g%tb(j-1)
            oik = mf_off(g, i, k)
            ojk = mf_off(g, j, k)
            oij = mf_off(g, i, j)
            call hecmw_mf_kernel_gemm(hi, wj, pk, hi, fval(oik+1), wj, fvalu(ojk+1), hi, fval(oij+1))
            call hecmw_mf_kernel_gemm(hi, wj, pk, hi, fvalu(oik+1), wj, fval(ojk+1), hi, fvalu(oij+1))
          enddo
          !$omp end taskloop
        endif
        do x = npiv + 1, g%tb(k)
          call update_pos(x, p0 + 1, npiv)
        enddo
      endif
      if (nfs <= g%tb(k)) exit
    enddo

    if (isroot .and. npiv < g%ncol) then
      pa = npiv + 1
      pb = g%ncol
      w = pb - pa + 1
      m = g%nrow - pa + 1
      call mf_grow_r(pval, int(m, 8)*w)
      call mf_grow_r(pvalu, int(m, 8)*w)
      call fill_panel(.false.)
      call hecmw_mf_kernel_panel_lu_piv(m, w, m, pval, m, pvalu, blk(pa:pb), blkr(pa:pb), u, &
        zero, .false., np, permc, permr, nsw, info)
      call permute_front(pa, w)
      nswap = nswap + nsw
      if (info /= 0) then
        ierr = sn%fsdof(pa + info - 1)
        deallocate(permc, permr, itmp, ip, jp)
        return
      endif
      call write_back(w)
      npiv = g%ncol
    endif
    sn%npiv = npiv
    deallocate(permc, permr, itmp, ip, jp)

  contains

    !> pval(m,w) <- the columns pa..pb of the front (rows pa..nrow), pvalu(m,w) <- the rows
    !> pa..pb (columns beyond the diagonal), transposed; with upd they receive the pivots
    !> p0+1..npiv of the current tile column, which the front does not carry yet
    subroutine fill_panel(upd)
      logical, intent(in) :: upd
      integer(kind=kint) :: jj, xx

      do jj = 1, w
        xx = pa + jj - 1
        call mf_col_copy(g, fval, xx, xx, pval(int(jj-1, 8)*m + jj), .false.)
        if (upd) then
          call mf_row_get(g, fvalu, xx, p0+1, npiv, wk)
          call mf_col_axpy(g, fval, xx, p0+1, npiv, wk, pval(int(jj-1, 8)*m + jj))
        endif
        if (xx < g%nrow) then
          call mf_col_copy(g, fvalu, xx, xx+1, pvalu(int(jj-1, 8)*m + jj + 1), .false.)
          if (upd) then
            call mf_row_get(g, fval, xx, p0+1, npiv, wk)
            call mf_col_axpy(g, fvalu, xx+1, p0+1, npiv, wk, pvalu(int(jj-1, 8)*m + jj + 1))
          endif
        endif
      enddo
    end subroutine fill_panel

    !> the first n panel columns (L, U diagonal) and rows (U) back to the front
    subroutine write_back(n)
      integer(kind=kint), intent(in) :: n
      integer(kind=kint) :: jj, xx

      do jj = 1, n
        xx = pa + jj - 1
        call mf_col_copy(g, fval, xx, xx, pval(int(jj-1, 8)*m + jj), .true.)
        if (xx < g%nrow) call mf_col_copy(g, fvalu, xx, xx+1, pvalu(int(jj-1, 8)*m + jj + 1), .true.)
      enddo
    end subroutine write_back

    !> column x (lower grid) and row x (upper grid) of the front receive the pivots q1..q2
    subroutine update_pos(x, q1, q2)
      integer(kind=kint), intent(in) :: x, q1, q2

      call mf_col_copy(g, fval, x, x, pval, .false.)
      call mf_row_get(g, fvalu, x, q1, q2, wk)
      call mf_col_axpy(g, fval, x, q1, q2, wk, pval)
      call mf_col_copy(g, fval, x, x, pval, .true.)
      if (x < g%nrow) then
        call mf_col_copy(g, fvalu, x, x+1, pvalu, .false.)
        call mf_row_get(g, fval, x, q1, q2, wk)
        call mf_col_axpy(g, fvalu, x+1, q1, q2, wk, pvalu)
        call mf_col_copy(g, fvalu, x, x+1, pvalu, .true.)
      endif
    end subroutine update_pos

    !> apply the panel permutations to the positions pa..pa+w-1 of the front: permc as
    !> symmetric exchanges, then the rows from that order to permr
    subroutine permute_front(pa, w)
      integer(kind=kint), intent(in) :: pa, w
      integer(kind=kint) :: jj, t, q

      do jj = 1, w
        itmp(jj) = jj
      enddo
      do jj = 1, w
        if (itmp(jj) == permc(jj)) cycle
        do t = jj + 1, w
          if (itmp(t) == permc(jj)) exit
        enddo
        call mf_swap_lu(g, fval, fvalu, pa+jj-1, pa+t-1)
        call swap_pos(pa+jj-1, pa+t-1)
        itmp(t) = itmp(jj)
        itmp(jj) = permc(jj)
      enddo
      itmp(1:w) = permc(1:w)
      do jj = 1, w
        if (itmp(jj) == permr(jj)) cycle
        do t = jj + 1, w
          if (itmp(t) == permr(jj)) exit
        enddo
        call mf_swap_row(g, fval, fvalu, pa+jj-1, pa+t-1)
        q = sn%frow(pa+jj-1)
        sn%frow(pa+jj-1) = sn%frow(pa+t-1)
        sn%frow(pa+t-1) = q
        q = blkr(pa+jj-1)
        blkr(pa+jj-1) = blkr(pa+t-1)
        blkr(pa+t-1) = q
        itmp(t) = itmp(jj)
        itmp(jj) = permr(jj)
      enddo
    end subroutine permute_front

    subroutine swap_pos(x, y)
      integer(kind=kint), intent(in) :: x, y
      integer(kind=kint) :: t

      t = sn%fsdof(x)
      sn%fsdof(x) = sn%fsdof(y)
      sn%fsdof(y) = t
      t = sn%frow(x)
      sn%frow(x) = sn%frow(y)
      sn%frow(y) = t
      t = blk(x)
      blk(x) = blk(y)
      blk(y) = t
      t = blkr(x)
      blkr(x) = blkr(y)
      blkr(y) = t
    end subroutine swap_pos

  end subroutine mf_factor_front_lu

  !> Tile grid of the stored panel of a supernode: rows as factored, tile columns cut at npiv.
  subroutine mf_stored_grid(sn, g, tw)
    implicit none
    type(mf_snode), intent(in) :: sn
    type(mf_grid), intent(inout) :: g
    integer(kind=kint), intent(out) :: tw(:)
    integer(kind=kint) :: j

    g%nt = sn%nt
    g%ntc = sn%ntc
    g%ncol = sn%ncol
    if (allocated(g%tb)) then
      if (size(g%tb) < g%nt + 1) deallocate(g%tb, g%coloff)
    endif
    if (.not. allocated(g%tb)) allocate(g%tb(0:g%nt), g%coloff(g%nt+1))
    g%tb(0:g%nt) = sn%tbnd(1:g%nt+1)
    g%nrow = g%tb(g%nt)
    g%coloff(1) = 0
    do j = 1, g%ntc
      tw(j) = max(min(g%tb(j), sn%npiv) - g%tb(j-1), 0)
      g%coloff(j+1) = g%coloff(j) + int(g%nrow - g%tb(j-1), 8)*tw(j)
    enddo
  end subroutine mf_stored_grid

  !> x = A^-1 b with the factor; b and x are in the original DOF numbering. The forward solve
  !> works on the row DOFs (z, indexed by frow), the backward solve on the column DOFs (xx).
  !> A front reads the b entries of its own DOFs only; the forward result of the rows beyond
  !> its pivots is kept per supernode (dval) and added by the parent, mirroring the extend-add.
  subroutine hecmw_mf_numeric_solve(sym, fct, b, x)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    real(kind=kreal), intent(in) :: b(:)
    real(kind=kreal), intent(out) :: x(:)
    real(kind=kreal), allocatable :: z(:), xx(:)
    integer(kind=kint), allocatable :: map(:), left(:)
    integer(kind=kint) :: ns, s, i

    ns = fct%nsuper
    allocate(z(fct%ndof_tot), xx(fct%ndof_tot), map(fct%ndof_tot), left(ns))
    !$omp parallel do
    do i = 1, fct%ndof_tot
      z(i) = b(fct%pdof(i))
    enddo
    !$omp end parallel do

    do s = 1, ns
      left(s) = fct%cptr(s+1) - fct%cptr(s)
    enddo
    !$omp parallel default(shared)
    !$omp single
    do s = 1, ns
      if (fct%cptr(s+1) == fct%cptr(s)) then
        !$omp task default(shared) firstprivate(s)
        call mf_fwd_task(sym, fct, left, s, z, map)
        !$omp end task
      endif
    enddo
    !$omp end single
    !$omp end parallel

    !$omp parallel default(shared)
    !$omp single
    do s = 1, ns
      if (sym%sparent(s) == 0) then
        !$omp task default(shared) firstprivate(s)
        call bwd_task(s)
        !$omp end task
      endif
    enddo
    !$omp end single
    !$omp end parallel

    !$omp parallel do
    do i = 1, fct%ndof_tot
      x(fct%pdof(i)) = xx(i)
    enddo
    !$omp end parallel do
    deallocate(z, xx, map, left)

  contains

    !> backward-solve task of one supernode: solve the front, then spawn the children, which
    !> depend on their parent alone; host association keeps every reference alive until the
    !> parallel region ends
    recursive subroutine bwd_task(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: j, c

      call mf_super_bwd(sym, fct, fct%sn(s0), s0, z, xx)
      do j = fct%cptr(s0), fct%cptr(s0+1) - 1
        c = fct%clist(j)
        !$omp task default(shared) firstprivate(c)
        call bwd_task(c)
        !$omp end task
      enddo
    end subroutine bwd_task

  end subroutine hecmw_mf_numeric_solve

  !> Forward-solve task starting at leaf s, climbing to released parents the way the
  !> factorization tasks do.
  subroutine mf_fwd_task(sym, fct, left, s, z, map)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(inout) :: left(:)
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(inout) :: z(:)
    integer(kind=kint), intent(inout) :: map(:)
    integer(kind=kint) :: p, n, ss

    ss = s
    do
      call mf_super_fwd(sym, fct, fct%sn(ss), ss, z, map)
      p = sym%sparent(ss)
      if (p == 0) return
      ! critical, not atomic: its flush semantics make the child's writes visible to the
      ! thread that continues with the parent
      !$omp critical (mf_tree)
      left(p) = left(p) - 1
      n = left(p)
      !$omp end critical (mf_tree)
      if (n /= 0) return
      ss = p
    enddo
  end subroutine mf_fwd_task

  !> x = A^-1 b on the supernodes distributed by map; b and x are the replicated global
  !> vectors in the original DOF numbering. The forward solve mirrors the fan-in of the
  !> factorization (subtree tasks, then the upper fronts ascending, the dval contributions
  !> of remote children received in place); the backward solve descends the upper fronts and
  !> hands every remote child the solution values its rows reference, a set both sides
  !> derive from the symbolic structure. The solved pivot values of the owned supernodes are
  !> summed over the ranks at the end, which is exact and order independent because every
  !> DOF is eliminated on exactly one rank.
  subroutine hecmw_mf_numeric_solve_mpi(sym, map, fct, b, x)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    real(kind=kreal), intent(in) :: b(:)
    real(kind=kreal), intent(out) :: x(:)
    type(mf_sendbox), allocatable :: box(:)
    real(kind=kreal), allocatable :: z(:), xx(:)
    integer(kind=kint), allocatable :: dmap(:), left(:), reqs(:), stats(:,:)
    integer(kind=kint) :: ns, s, i, iu, j, c, p, m, nreq, ib, nsend, ndn, stat(HECMW_STATUS_SIZE)

    ns = fct%nsuper
    allocate(z(fct%ndof_tot), xx(fct%ndof_tot), dmap(fct%ndof_tot), left(ns))
    !$omp parallel do
    do i = 1, fct%ndof_tot
      z(i) = b(fct%pdof(i))
    enddo
    !$omp end parallel do

    do s = 1, ns
      left(s) = fct%cptr(s+1) - fct%cptr(s)
    enddo
    !$omp parallel default(shared)
    !$omp single
    do s = 1, ns
      if (fct%cptr(s+1) == fct%cptr(s) .and. map%owner(s) == map%myrank .and. .not. map%upper(s)) then
        !$omp task default(shared) firstprivate(s)
        call mf_fwd_task_mpi(sym, map, fct, left, s, z, dmap)
        !$omp end task
      endif
    enddo
    !$omp end single
    !$omp end parallel

    ! the forward stage sends one message per owned supernode with a remote parent, the
    ! backward stage one per remote child of an owned upper front; the counts differ per rank
    nsend = 0
    do s = 1, ns
      if (map%owner(s) == map%myrank .and. sym%sparent(s) /= 0) then
        if (map%owner(sym%sparent(s)) /= map%myrank) nsend = nsend + 1
      endif
    enddo
    ndn = 0
    do iu = 1, map%nupper
      s = map%uplist(iu)
      if (map%owner(s) /= map%myrank) cycle
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        if (map%owner(fct%clist(j)) /= map%myrank) ndn = ndn + 1
      enddo
    enddo
    allocate(box(max(ndn, 1)), reqs(max(nsend, ndn, 1)), stats(HECMW_STATUS_SIZE, max(nsend, ndn, 1)))
    nreq = 0
    do s = 1, ns
      if (map%owner(s) /= map%myrank .or. map%upper(s) .or. sym%sparent(s) == 0) cycle
      if (map%upper(sym%sparent(s)) .and. map%owner(sym%sparent(s)) /= map%myrank) call send_dval(s)
    enddo
    do iu = 1, map%nupper
      s = map%uplist(iu)
      if (map%owner(s) /= map%myrank) cycle
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        c = fct%clist(j)
        if (map%owner(c) == map%myrank) cycle
        m = fct%sn(c)%tbnd(fct%sn(c)%nt+1) - fct%sn(c)%npiv
        if (m > 0) then
          if (allocated(fct%sn(c)%dval)) deallocate(fct%sn(c)%dval)
          allocate(fct%sn(c)%dval(m))
          call hecmw_recv_r(fct%sn(c)%dval, m, map%owner(c), 3*c, map%comm, stat)
        endif
      enddo
      call mf_super_fwd(sym, fct, fct%sn(s), s, z, dmap)
      if (sym%sparent(s) /= 0) then
        if (map%owner(sym%sparent(s)) /= map%myrank) call send_dval(s)
      endif
    enddo
    call hecmw_waitall(nreq, reqs, stats)
    ! the sent forward contributions were consumed remotely; free the local copies
    do s = 1, ns
      if (map%owner(s) == map%myrank .and. sym%sparent(s) /= 0) then
        if (map%owner(sym%sparent(s)) /= map%myrank .and. allocated(fct%sn(s)%dval)) deallocate(fct%sn(s)%dval)
      endif
    enddo

    nreq = 0
    ib = 0
    do iu = map%nupper, 1, -1
      s = map%uplist(iu)
      if (map%owner(s) /= map%myrank) cycle
      p = sym%sparent(s)
      if (p /= 0) then
        if (map%owner(p) /= map%myrank) call recv_xx(s)
      endif
      call mf_super_bwd(sym, fct, fct%sn(s), s, z, xx)
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        c = fct%clist(j)
        if (map%owner(c) /= map%myrank) call send_xx(c)
      enddo
    enddo
    do s = 1, ns
      if (map%owner(s) /= map%myrank .or. map%upper(s) .or. sym%sparent(s) == 0) cycle
      if (map%upper(sym%sparent(s)) .and. map%owner(sym%sparent(s)) /= map%myrank) call recv_xx(s)
    enddo
    !$omp parallel default(shared)
    !$omp single
    do s = 1, ns
      if (map%owner(s) /= map%myrank .or. map%upper(s)) cycle
      p = sym%sparent(s)
      if (p /= 0) then
        if (.not. map%upper(p)) cycle
      endif
      !$omp task default(shared) firstprivate(s)
      call bwd_task_mpi(s)
      !$omp end task
    enddo
    !$omp end single
    !$omp end parallel
    call hecmw_waitall(nreq, reqs, stats)

    z(1:fct%ndof_tot) = 0.0d0
    do s = 1, ns
      if (map%owner(s) /= map%myrank) cycle
      do i = 1, fct%sn(s)%npiv
        z(fct%pdof(fct%sn(s)%fsdof(i))) = xx(fct%sn(s)%fsdof(i))
      enddo
    enddo
    x(1:fct%ndof_tot) = z(1:fct%ndof_tot)
    call hecmw_allreduce_R_comm(x, fct%ndof_tot, hecmw_sum, map%comm)
    deallocate(z, xx, dmap, left, box, reqs, stats)

  contains

    !> isend the forward contribution of s0 to the owner of its parent
    subroutine send_dval(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: mm

      mm = fct%sn(s0)%tbnd(fct%sn(s0)%nt+1) - fct%sn(s0)%npiv
      if (mm <= 0) return
      nreq = nreq + 1
      call hecmw_isend_r(fct%sn(s0)%dval, mm, map%owner(sym%sparent(s0)), 3*s0, map%comm, reqs(nreq))
    end subroutine send_dval

    !> permuted DOFs whose solution the rows of supernode c0 beyond its pivots reference:
    !> the delayed column DOFs, then the column DOFs of the contribution row nodes
    subroutine xx_dofs(c0, mm, list)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint), intent(out) :: mm
      integer(kind=kint), intent(out), allocatable :: list(:)
      integer(kind=kint) :: cnown, d, k, l, idx

      cnown = sym%sptr(c0+1) - sym%sptr(c0)
      mm = fct%sn(c0)%ncol - fct%sn(c0)%npiv
      do l = sym%rptr(c0) + cnown, sym%rptr(c0+1) - 1
        mm = mm + sym%ndof(sym%rlist(l))
      enddo
      allocate(list(max(mm, 1)))
      idx = 0
      do d = fct%sn(c0)%npiv + 1, fct%sn(c0)%ncol
        idx = idx + 1
        list(idx) = fct%sn(c0)%fsdof(d)
      enddo
      do l = sym%rptr(c0) + cnown, sym%rptr(c0+1) - 1
        k = sym%rlist(l)
        do d = fct%cdofptr(k), fct%cdofptr(k+1) - 1
          idx = idx + 1
          list(idx) = d
        enddo
      enddo
    end subroutine xx_dofs

    !> isend the referenced solution values to the owner of child c0
    subroutine send_xx(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint), allocatable :: list(:)
      integer(kind=kint) :: mm, i0

      call xx_dofs(c0, mm, list)
      if (mm <= 0) return
      ib = ib + 1
      allocate(box(ib)%rv(mm))
      do i0 = 1, mm
        box(ib)%rv(i0) = xx(list(i0))
      enddo
      deallocate(list)
      nreq = nreq + 1
      call hecmw_isend_r(box(ib)%rv, mm, map%owner(c0), 3*c0+1, map%comm, reqs(nreq))
    end subroutine send_xx

    !> receive the referenced solution values of s0 from the owner of its parent
    subroutine recv_xx(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint), allocatable :: list(:)
      real(kind=kreal), allocatable :: v(:)
      integer(kind=kint) :: mm, i0

      call xx_dofs(s0, mm, list)
      if (mm <= 0) return
      allocate(v(mm))
      call hecmw_recv_r(v, mm, map%owner(sym%sparent(s0)), 3*s0+1, map%comm, stat)
      do i0 = 1, mm
        xx(list(i0)) = v(i0)
      enddo
      deallocate(v, list)
    end subroutine recv_xx

    !> backward-solve task over an owned subtree, every descendant local by construction
    recursive subroutine bwd_task_mpi(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: jj, cc

      call mf_super_bwd(sym, fct, fct%sn(s0), s0, z, xx)
      do jj = fct%cptr(s0), fct%cptr(s0+1) - 1
        cc = fct%clist(jj)
        !$omp task default(shared) firstprivate(cc)
        call bwd_task_mpi(cc)
        !$omp end task
      enddo
    end subroutine bwd_task_mpi

  end subroutine hecmw_mf_numeric_solve_mpi

  !> mf_fwd_task limited to the subtrees of the executing rank, the climb stopping below an
  !> upper front the way the factorization tasks do.
  subroutine mf_fwd_task_mpi(sym, map, fct, left, s, z, dmap)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(inout) :: left(:)
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(inout) :: z(:)
    integer(kind=kint), intent(inout) :: dmap(:)
    integer(kind=kint) :: p, n, ss

    ss = s
    do
      call mf_super_fwd(sym, fct, fct%sn(ss), ss, z, dmap)
      p = sym%sparent(ss)
      if (p == 0) return
      if (map%upper(p)) return
      ! critical, not atomic: its flush semantics make the child's writes visible to the
      ! thread that continues with the parent
      !$omp critical (mf_tree)
      left(p) = left(p) - 1
      n = left(p)
      !$omp end critical (mf_tree)
      if (n /= 0) return
      ss = p
    enddo
  end subroutine mf_fwd_task_mpi

  !> Forward solve of the front of supernode s: gather the b entries of the own DOFs from z,
  !> add the forward contributions of the children (freed here), apply the panel, write the
  !> solution of the pivot rows to z and keep the rows beyond the pivots in dval for the parent.
  !> Pivoting permutes the fully summed positions, so a child row landing there is located
  !> through map (row DOF to position); the map entries of concurrently processed fronts are
  !> disjoint because their subtrees share no fully summed DOFs. The contribution rows keep
  !> the assembly order and are located by the rlist offsets.
  subroutine mf_super_fwd(sym, fct, sn, s, z, map)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(inout) :: fct
    type(mf_snode), intent(inout) :: sn
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(inout) :: z(:)
    integer(kind=kint), intent(inout) :: map(:)
    type(mf_grid) :: g
    real(kind=kreal), allocatable :: v(:), t(:)
    integer(kind=kint), allocatable :: tw(:), rowoff(:)
    integer(kind=kint) :: nown, nrow_nodes, ncol0, ndel, d0, d1, d, c, i, j, k, l, a, r, rr
    integer(kind=kint) :: cnown, cndel, hk, hi
    integer(kind=8) :: okk, oik

    allocate(tw(sn%nt))
    call mf_stored_grid(sn, g, tw)
    allocate(v(g%nrow))
    d0 = fct%cdofptr(sym%sptr(s))
    d1 = fct%cdofptr(sym%sptr(s+1)) - 1
    do d = 1, sn%ncol
      map(sn%frow(d)) = d
      if (sn%frow(d) >= d0 .and. sn%frow(d) <= d1) then
        v(d) = z(sn%frow(d))
      else
        v(d) = 0.0d0
      endif
    enddo
    v(sn%ncol+1:g%nrow) = 0.0d0

    nown = sym%sptr(s+1) - sym%sptr(s)
    nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
    ncol0 = d1 - d0 + 1
    ndel = sn%ncol - ncol0
    allocate(rowoff(nrow_nodes+1))
    rowoff(1) = 0
    do i = 1, nrow_nodes
      k = sym%rlist(sym%rptr(s)+i-1)
      rowoff(i+1) = rowoff(i) + sym%ndof(k)
      if (i == nown) rowoff(i+1) = rowoff(i+1) + ndel
    enddo
    c = fct%chead(s)
    do while (c /= 0)
      cnown = sym%sptr(c+1) - sym%sptr(c)
      cndel = fct%sn(c)%ncol - fct%sn(c)%npiv
      do j = 1, cndel
        d = fct%sn(c)%frow(fct%sn(c)%npiv + j)
        v(map(d)) = v(map(d)) + fct%sn(c)%dval(j)
      enddo
      r = cndel
      do l = sym%rptr(c) + cnown, sym%rptr(c+1) - 1
        k = sym%cmap(sym%cmap_ptr(c) + (l - sym%rptr(c) - cnown))
        do a = 1, sym%ndof(sym%rlist(l))
          r = r + 1
          if (k <= nown) then
            d = fct%cdofptr(sym%rlist(sym%rptr(s)+k-1)) + a - 1
            v(map(d)) = v(map(d)) + fct%sn(c)%dval(r)
          else
            v(rowoff(k)+a) = v(rowoff(k)+a) + fct%sn(c)%dval(r)
          endif
        enddo
      enddo
      deallocate(fct%sn(c)%dval)
      c = fct%cnext(c)
    enddo
    deallocate(rowoff)

    if (fct%blr) allocate(t(maxval(g%tb(1:g%nt) - g%tb(0:g%nt-1))))
    do k = 1, g%ntc
      if (tw(k) == 0) exit
      hk = g%tb(k) - g%tb(k-1)
      okk = 1 + g%coloff(k)
      if (fct%blr) okk = 1 + sn%bptr(mf_bidx(g%nt, k, k))
      call hecmw_mf_kernel_trsv(tw(k), hk, sn%lval(okk), v(g%tb(k-1)+1))
      if (hk > tw(k)) call hecmw_mf_kernel_gemv(hk - tw(k), tw(k), hk, sn%lval(okk + tw(k)), v(g%tb(k-1)+1), &
        v(g%tb(k-1)+tw(k)+1))
      do i = k+1, g%nt
        hi = g%tb(i) - g%tb(i-1)
        rr = -1
        if (fct%blr) then
          rr = sn%brank(mf_bidx(g%nt, k, i))
          oik = 1 + sn%bptr(mf_bidx(g%nt, k, i))
        else
          oik = 1 + g%coloff(k) + int(g%tb(i-1) - g%tb(k-1), 8)*tw(k)
        endif
        if (rr == 0) cycle
        if (rr < 0) then
          call hecmw_mf_kernel_gemv(hi, tw(k), hi, sn%lval(oik), v(g%tb(k-1)+1), v(g%tb(i-1)+1))
        else
          call hecmw_mf_kernel_mult_tv(tw(k), rr, tw(k), sn%lval(oik + int(hi, 8)*rr), v(g%tb(k-1)+1), t)
          call hecmw_mf_kernel_gemv(hi, rr, hi, sn%lval(oik), t, v(g%tb(i-1)+1))
        endif
      enddo
      if (.not. fct%lu) call hecmw_mf_kernel_dsolve(tw(k), hk, sn%lval(okk), sn%ptype(g%tb(k-1)+1:), &
        sn%dsub(g%tb(k-1)+1:), v(g%tb(k-1)+1))
    enddo
    if (allocated(t)) deallocate(t)

    do d = 1, sn%npiv
      z(sn%frow(d)) = v(d)
    enddo
    if (g%nrow > sn%npiv) then
      allocate(sn%dval(g%nrow - sn%npiv))
      sn%dval(1:g%nrow-sn%npiv) = v(sn%npiv+1:g%nrow)
    endif
    deallocate(v, tw)
  end subroutine mf_super_fwd

  !> Backward solve of the front of supernode s: the pivot rows from z, the rows beyond from
  !> the already solved column DOFs in xx, the pivot solutions back to xx.
  subroutine mf_super_bwd(sym, fct, sn, s, z, xx)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_factor), intent(in) :: fct
    type(mf_snode), intent(in) :: sn
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(in) :: z(:)
    real(kind=kreal), intent(inout) :: xx(:)
    type(mf_grid) :: g
    real(kind=kreal), allocatable :: v(:), t(:)
    integer(kind=kint), allocatable :: tw(:)
    integer(kind=kint) :: nown, d, i, k, l, m, hk, hi, rr
    integer(kind=8) :: okk, okku, oik

    allocate(tw(sn%nt))
    call mf_stored_grid(sn, g, tw)
    allocate(v(g%nrow))
    do d = 1, sn%npiv
      v(d) = z(sn%frow(d))
    enddo
    do d = sn%npiv + 1, sn%ncol
      v(d) = xx(sn%fsdof(d))
    enddo
    nown = sym%sptr(s+1) - sym%sptr(s)
    m = sn%ncol
    do l = sym%rptr(s) + nown, sym%rptr(s+1) - 1
      k = sym%rlist(l)
      do d = fct%cdofptr(k), fct%cdofptr(k+1) - 1
        m = m + 1
        v(m) = xx(d)
      enddo
    enddo

    if (fct%blr) allocate(t(maxval(g%tb(1:g%nt) - g%tb(0:g%nt-1))))
    do k = g%ntc, 1, -1
      if (tw(k) == 0) cycle
      hk = g%tb(k) - g%tb(k-1)
      okk = 1 + g%coloff(k)
      okku = okk
      if (fct%blr) then
        okk = 1 + sn%bptr(mf_bidx(g%nt, k, k))
        if (fct%lu) okku = 1 + sn%bptru(mf_bidx(g%nt, k, k))
      endif
      if (fct%lu) then
        do i = k+1, g%nt
          hi = g%tb(i) - g%tb(i-1)
          rr = -1
          if (fct%blr) then
            rr = sn%branku(mf_bidx(g%nt, k, i))
            oik = 1 + sn%bptru(mf_bidx(g%nt, k, i))
          else
            oik = 1 + g%coloff(k) + int(g%tb(i-1) - g%tb(k-1), 8)*tw(k)
          endif
          if (rr == 0) cycle
          if (rr < 0) then
            call hecmw_mf_kernel_gemv_t(hi, tw(k), hi, sn%uval(oik), v(g%tb(i-1)+1), v(g%tb(k-1)+1))
          else
            call hecmw_mf_kernel_mult_tv(hi, rr, hi, sn%uval(oik), v(g%tb(i-1)+1), t)
            call hecmw_mf_kernel_gemv(tw(k), rr, tw(k), sn%uval(oik + int(hi, 8)*rr), t, v(g%tb(k-1)+1))
          endif
        enddo
        if (hk > tw(k)) call hecmw_mf_kernel_gemv_t(hk - tw(k), tw(k), hk, sn%uval(okku + tw(k)), &
          v(g%tb(k-1)+tw(k)+1), v(g%tb(k-1)+1))
        call hecmw_mf_kernel_usolve(tw(k), hk, sn%lval(okk), hk, sn%uval(okku), v(g%tb(k-1)+1))
      else
        do i = k+1, g%nt
          hi = g%tb(i) - g%tb(i-1)
          rr = -1
          if (fct%blr) then
            rr = sn%brank(mf_bidx(g%nt, k, i))
            oik = 1 + sn%bptr(mf_bidx(g%nt, k, i))
          else
            oik = 1 + g%coloff(k) + int(g%tb(i-1) - g%tb(k-1), 8)*tw(k)
          endif
          if (rr == 0) cycle
          if (rr < 0) then
            call hecmw_mf_kernel_gemv_t(hi, tw(k), hi, sn%lval(oik), v(g%tb(i-1)+1), v(g%tb(k-1)+1))
          else
            call hecmw_mf_kernel_mult_tv(hi, rr, hi, sn%lval(oik), v(g%tb(i-1)+1), t)
            call hecmw_mf_kernel_gemv(tw(k), rr, tw(k), sn%lval(oik + int(hi, 8)*rr), t, v(g%tb(k-1)+1))
          endif
        enddo
        if (hk > tw(k)) call hecmw_mf_kernel_gemv_t(hk - tw(k), tw(k), hk, sn%lval(okk + tw(k)), &
          v(g%tb(k-1)+tw(k)+1), v(g%tb(k-1)+1))
        call hecmw_mf_kernel_trsv_t(tw(k), hk, sn%lval(okk), v(g%tb(k-1)+1))
      endif
    enddo
    if (allocated(t)) deallocate(t)

    do d = 1, sn%npiv
      xx(sn%fsdof(d)) = v(d)
    enddo
    deallocate(v, tw)
  end subroutine mf_super_bwd

  subroutine hecmw_mf_numeric_print(fct)
    implicit none
    type(hecmwST_mf_factor), intent(in) :: fct

    write(*,'(a,i0,a,i0,a,i0)') '[DIRECTmf]: tile size = ', fct%tile, ', tiles (max per front) = ', fct%max_tiles, &
      ', tile dim (max) = ', fct%max_tile_dim
    write(*,'(a,i0,a,f10.3,a)') '[DIRECTmf]: factor words (estimate) = ', fct%factor_words, &
      ' (', real(fct%factor_words, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
    write(*,'(a,i0,a,i0,a,f10.3,a)') '[DIRECTmf]: stack peak words = ', fct%stack_peak, ', front words = ', &
      fct%front_words, ' (', real(fct%stack_peak + fct%front_words, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
    write(*,'(a,i0,a,i0,a,i0,a,f10.3,a)') '[DIRECTmf]: LU mode estimate: factor words = ', 2*fct%factor_words, &
      ', stack peak words = ', 2*fct%stack_peak, ', front words = ', 2*fct%front_words, ' (', &
      real(2*(fct%factor_words + fct%stack_peak + fct%front_words), kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
  end subroutine hecmw_mf_numeric_print

end module hecmw_mf_numeric
