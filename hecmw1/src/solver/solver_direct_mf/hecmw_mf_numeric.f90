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
  public :: hecmw_mf_numeric_front_words
  public :: hecmw_mf_numeric_finalize

  !> default pivot threshold: entries of L are bounded by 1/u
  real(kind=kreal), parameter :: MF_PIVOT_U = 0.01d0
  !> Bunch-Kaufman constant
  real(kind=kreal), parameter :: MF_PIVOT_ALPHA = (1.0d0 + sqrt(17.0d0)) / 8.0d0
  !> default zero pivot fraction of the largest entry of the matrix
  real(kind=kreal), parameter :: MF_PIVOT_ZERO = 1.0d-14
  !> default BLR truncation threshold
  real(kind=kreal), parameter :: MF_BLR_EPS = 1.0d-8
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
    integer(kind=8) :: pwords = 0                   !< stored panel words held on this rank
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
    integer(kind=kint) :: mode = 0                  !< 0 follow the symmetric flag, 1 force LDLt, 2 force LU
    logical :: scan = .false.                       !< scan the numerical asymmetry (log output only)
    logical :: blr = .false.                        !< compress the factor panels (BLR)
    real(kind=kreal) :: eps = 0.0d0                 !< BLR truncation threshold
    real(kind=kreal) :: pivot_u = 0.0d0             !< pivot threshold: entries of L bounded by 1/u
    real(kind=kreal) :: pivot_zero = 0.0d0          !< zero pivot fraction of max|A|
    ! layout of the last factorization
    logical :: lu = .false.                         !< LU mode (else LDLt)
    real(kind=kreal) :: asym = 0.0d0                !< max|A_ij - A_ji| / max|A_ij| of the last matrix (0 unless scanned)
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

  !> pool of pending isends: every send buffered in a box until the waitall of the phase;
  !> one box may back several isends (a multicast), so the requests are counted apart
  type mf_pool
    integer(kind=kint) :: n = 0
    integer(kind=kint) :: nreq = 0
    integer(kind=kint), allocatable :: reqs(:)
    type(mf_sendbox), allocatable :: box(:)
  end type mf_pool

contains

  !> next free box of the pool, growing it as needed
  function mf_pool_slot(pool) result(ib)
    implicit none
    type(mf_pool), intent(inout) :: pool
    integer(kind=kint) :: ib
    type(mf_sendbox), allocatable :: t(:)
    integer(kind=kint) :: i, m

    pool%n = pool%n + 1
    ib = pool%n
    if (.not. allocated(pool%box)) then
      allocate(pool%box(64))
    else if (pool%n > size(pool%box)) then
      m = 2*size(pool%box)
      allocate(t(m))
      do i = 1, size(pool%box)
        if (allocated(pool%box(i)%hdr)) call move_alloc(pool%box(i)%hdr, t(i)%hdr)
        if (allocated(pool%box(i)%tl)) call move_alloc(pool%box(i)%tl, t(i)%tl)
        if (allocated(pool%box(i)%rv)) call move_alloc(pool%box(i)%rv, t(i)%rv)
      enddo
      call move_alloc(t, pool%box)
    endif
  end function mf_pool_slot

  !> record one pending request
  subroutine mf_pool_req(pool, rq)
    implicit none
    type(mf_pool), intent(inout) :: pool
    integer(kind=kint), intent(in) :: rq

    pool%nreq = pool%nreq + 1
    call mf_grow_i(pool%reqs, max(pool%nreq, 64))
    pool%reqs(pool%nreq) = rq
  end subroutine mf_pool_req

  !> wait for every pending isend of the pool and release the boxes
  subroutine mf_pool_wait(pool)
    implicit none
    type(mf_pool), intent(inout) :: pool
    integer(kind=kint), allocatable :: stats(:,:)
    integer(kind=kint) :: i

    if (pool%nreq > 0) then
      allocate(stats(HECMW_STATUS_SIZE, pool%nreq))
      call hecmw_waitall(pool%nreq, pool%reqs, stats)
      deallocate(stats)
    endif
    do i = 1, pool%n
      if (allocated(pool%box(i)%hdr)) deallocate(pool%box(i)%hdr)
      if (allocated(pool%box(i)%tl)) deallocate(pool%box(i)%tl)
      if (allocated(pool%box(i)%rv)) deallocate(pool%box(i)%rv)
    enddo
    pool%n = 0
    pool%nreq = 0
  end subroutine mf_pool_wait

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

  !> The fully summed part of a front as a grid of its own: same tile columns, rows cut at
  !> ncol. The sequential position helpers then work unchanged on the master part of a row
  !> distributed front.
  subroutine mf_fs_grid(g, gm)
    implicit none
    type(mf_grid), intent(in) :: g
    type(mf_grid), intent(inout) :: gm
    integer(kind=kint) :: j

    gm%nt = g%ntc
    gm%ntc = g%ntc
    gm%ncol = g%ncol
    gm%nrow = g%ncol
    if (allocated(gm%tb)) then
      if (size(gm%tb) < gm%nt + 1) deallocate(gm%tb, gm%coloff)
    endif
    if (.not. allocated(gm%tb)) allocate(gm%tb(0:gm%nt), gm%coloff(gm%nt+1))
    if (allocated(gm%dtile)) then
      if (size(gm%dtile) < gm%nrow) deallocate(gm%dtile)
    endif
    if (.not. allocated(gm%dtile)) allocate(gm%dtile(0:max(gm%nrow, 1)-1))
    gm%tb(0:gm%nt) = g%tb(0:g%ntc)
    gm%coloff(1) = 0
    do j = 1, gm%nt
      gm%coloff(j+1) = gm%coloff(j) + int(gm%nrow - gm%tb(j-1), 8)*(gm%tb(j) - gm%tb(j-1))
      gm%dtile(gm%tb(j-1):gm%tb(j)-1) = j
    enddo
  end subroutine mf_fs_grid

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

  !> The mode from the option and the symmetric flag, and, only when the scan option is on,
  !> the asymmetry of the values, max|A_ij - A_ji| relative to amax (log output only).
  subroutine mf_asymmetry(hecMAT, fct, amax)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_factor), intent(inout) :: fct
    real(kind=kreal), intent(in) :: amax
    integer(kind=kint) :: nd, i, kk, a, b
    integer(kind=8) :: base, base2
    real(kind=kreal) :: asym

    fct%lu = (fct%mode == 2) .or. (fct%mode == 0 .and. .not. hecMAT%symmetric)
    fct%asym = 0.0d0
    if (.not. fct%scan) return
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
    if (amax > 0.0d0) fct%asym = asym / amax
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
  !> in LU mode when the symmetric flag is off; the mode option overrides either way. ierr is
  !> 0 on success, the permuted DOF of a zero pivot (singular matrix), -1 when the block size
  !> of hecMAT does not match the structure, or -2 when the structure is not symmetric.
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
  !> factors its subtrees with the task-parallel code, then every rank of a rank set walks
  !> the upper fronts in ascending order. The assembly, the factorization and the storage
  !> of an upper front are distributed over its rank set by contribution row tiles (1D row
  !> distribution), the fully summed part staying with the master (the owner), which keeps
  !> the pivot sequence identical to the sequential factorization. MPI calls stay on the
  !> master thread outside the parallel regions because hecmw initializes MPI without a
  !> threading level. Sends are isend and are consumed within the phase of their front, so
  !> a rank blocked in a receive has already issued every send an earlier front needs,
  !> which makes the ascending walk deadlock free. Tags encode 16*supernode+kind, relying
  !> on a tag space larger than the MPI minimum of 32767, as every mainstream MPI provides.
  subroutine hecmw_mf_numeric_factor_mpi(hecMAT, sym, map, fct, ierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(out) :: ierr
    type(mf_work), allocatable :: wrks(:)
    integer(kind=kint), allocatable :: left(:)
    integer(kind=kint) :: s, k, nd, nthr, iu, gierr, iw(1)
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

    ! no barrier between the stages: a rank done with its subtrees walks straight into the
    ! upper fronts, and a subtree error travels in the child headers, uniformly skipping
    ! every front above it until the final allreduce settles ierr
    do iu = 1, map%nupper
      s = map%uplist(iu)
      if (map%myrank >= map%rbeg(s) .and. map%myrank < map%rbeg(s) + map%rcnt(s)) then
        call mf_super_factor_1d(hecMAT, sym, map, fct, s, zero, gierr)
      endif
    enddo
    iw(1) = gierr
    call hecmw_allreduce_I_comm(iw, 1, hecmw_max, map%comm)
    gierr = iw(1)
    ierr = gierr
    if (ierr /= 0) return
    fct%factored = .true.
  end subroutine hecmw_mf_numeric_factor_mpi

  !> One upper front under the 1D row distribution, executed by every rank of its rank set
  !> in lockstep: the fully summed part of the front lives on the master (the owner) as a
  !> grid of its own, the contribution row tiles live on their owning ranks as row bands
  !> (rows of a tile by all columns up to its diagonal), and a rank assembles, updates and
  !> stores what it holds. The extend-add routes every child contribution entry directly to
  !> the rank holding its front row, both sides deriving the partition from the symbolic
  !> structure and the child headers. The pivot sequence stays identical to the sequential
  !> factorization: the panel of a tile column is attempted without pivoting (the master
  !> factoring the fully summed rows, the tile owners solving and checking their rows, one
  !> max reduction), and on a rejection the whole panel is gathered to the master, which
  !> runs the sequential pivoting kernel and scatters the result. Every message of a front
  !> is consumed within its phase, keeping the ascending walk deadlock free.
  subroutine mf_super_factor_1d(hecMAT, sym, map, fct, s, zero, gierr)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    integer(kind=kint), intent(in) :: s
    real(kind=kreal), intent(in) :: zero
    integer(kind=kint), intent(inout) :: gierr

    type(mf_pool) :: pool
    type(mf_grid) :: g, gm
    integer(kind=kint), allocatable :: ctb(:), towner(:), wrank(:), crow2t(:), tprow(:)
    integer(kind=8), allocatable :: boff(:)
    real(kind=kreal), allocatable, target :: band(:), bandu(:), fvm(:), fvmu(:)
    integer(kind=kint), allocatable :: rowoff(:), pos(:), blk(:), blkr(:), mdel(:)
    integer(kind=kint), allocatable :: perm(:), permr(:), itmp(:), swaps(:)
    integer(kind=kint), allocatable :: cctb(:), cown(:), pmap(:), pair_i(:), pair_j(:)
    integer(kind=kint), allocatable :: cctb2(:), cown2(:), ptv(:)
    real(kind=kreal), allocatable :: pv(:), pvu(:), wk(:), tmat(:), tmatu(:), dgb(:), dgbu(:)
    real(kind=kreal), allocatable :: pfull(:), pfullu(:), cmax(:), pband(:), pbandu(:)
    real(kind=kreal), allocatable, target :: fsl(:), fslu(:), wscm(:), wscb(:), dum(:)
    real(kind=kreal), allocatable :: tvec(:), tvecu(:), dsv(:), dgv(:), dd(:)
    real(kind=kreal), allocatable, target :: bval(:), bvalu(:)
    integer(kind=8), allocatable :: bslot(:), bslotu(:)
    integer(kind=kint), allocatable :: brk(:), brku(:), xrk(:), xrku(:), xw(:)
    integer(kind=8), allocatable :: xoff(:), xoffu(:)
    type rbuf_t
      real(kind=kreal), allocatable :: v(:)
    end type rbuf_t
    type(rbuf_t), allocatable, target :: xbuf(:)
    integer(kind=kint) :: ctli(8), stat(HECMW_STATUS_SIZE)
    integer(kind=kint) :: nd, ncol0, ndel, ncol, nrow, ncb, ncbt, nwk, nmine, myrows
    integer(kind=kint) :: master, me, i, j, k, c, t, m, nkc
    integer(kind=kint) :: npiv, nfs, p0, pa, pb, w, np, pk, pk0, mfs, info, nsw, n22, cndel
    integer(kind=kint) :: nown, nrow_nodes, hi, npair, req, cnb, kbegc, cnrow, wme, ferr
    integer(kind=8) :: bandw, fw, btop, btopu
    real(kind=kreal) :: uinv
    logical :: ismaster, hasband, fsok, acc
    integer(kind=kint) :: nswap_f, n2x2_f, npos_f, nneg_f, ntl_f, nlr_f, rmax_f
    integer(kind=8) :: rsum_f

    nd = hecMAT%NDOF
    me = map%myrank
    master = map%owner(s)
    ismaster = (master == me)
    nown = sym%sptr(s+1) - sym%sptr(s)
    nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
    ncol0 = fct%cdofptr(sym%sptr(s+1)) - fct%cdofptr(sym%sptr(s))
    uinv = 1.0d0 / fct%pivot_u

    ! (1) headers of the children: the rank holding a child multicasts, everyone else
    ! receives, so that all ranks know the delayed growth and the front geometry
    ferr = 0
    do j = fct%cptr(s), fct%cptr(s+1) - 1
      c = fct%clist(j)
      if (map%owner(c) == me) call send_hdr(c)
    enddo
    do j = fct%cptr(s), fct%cptr(s+1) - 1
      c = fct%clist(j)
      if (map%owner(c) /= me) call recv_hdr(c)
    enddo
    if (ferr /= 0) then
      ! a child header carries an error: the whole rank set saw the same headers, so it
      ! skips the front uniformly (no value traffic follows, the message skeleton stays
      ! deterministic), leaving a flagged minimal header for the owner to forward to the
      ! parent front
      if (gierr == 0) gierr = ferr
      fct%sn(s)%ncol = 0
      fct%sn(s)%npiv = 0
      fct%sn(s)%nt = 0
      fct%sn(s)%ntc = 0
      call mf_grow_i(fct%sn(s)%tbnd, 1)
      call mf_grow_i(fct%sn(s)%fsdof, 1)
      call mf_grow_i(fct%sn(s)%frow, 1)
      fct%sn(s)%tbnd(1) = 0
      call mf_pool_wait(pool)
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        c = fct%clist(j)
        if (allocated(fct%sn(c)%cval)) then
          fct%live_cb = fct%live_cb - fct%sn(c)%cbsize
          deallocate(fct%sn(c)%cval)
        endif
      enddo
      return
    endif

    ! (2) geometry: the full grid, the fully summed grid of the master, the contribution
    ! row tiles with their owners, and the row offsets of the assembly
    ndel = 0
    allocate(mdel(max(fct%cptr(s+1) - fct%cptr(s), 1)))
    c = fct%chead(s)
    m = ncol0
    do while (c /= 0)
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        if (fct%clist(j) == c) mdel(j - fct%cptr(s) + 1) = m
      enddo
      m = m + fct%sn(c)%ncol - fct%sn(c)%npiv
      c = fct%cnext(c)
    enddo
    ndel = m - ncol0
    if (ismaster) fct%max_growth = max(fct%max_growth, ndel)
    call mf_partition(sym, fct, s, ndel, g)
    ncol = g%ncol
    nrow = g%nrow
    ncb = nrow - ncol
    call hecmw_mf_dist_rowtiles(sym, map, fct%tile, s, ndel, ncbt, ctb, towner)
    call mf_fs_grid(g, gm)
    allocate(crow2t(max(ncb, 1)))
    do t = 1, ncbt
      crow2t(ctb(t-1)+1:ctb(t)) = t
    enddo
    ! worker ranks ascending (contiguous tile blocks keep towner monotone), my tiles
    allocate(wrank(max(ncbt, 1)), tprow(0:max(ncbt, 1)), boff(max(ncbt, 1)))
    nwk = 0
    do t = 1, ncbt
      if (nwk == 0) then
        nwk = 1
        wrank(1) = towner(t)
      else if (towner(t) /= wrank(nwk)) then
        nwk = nwk + 1
        wrank(nwk) = towner(t)
      endif
    enddo
    nmine = 0
    myrows = 0
    bandw = 0
    tprow(0) = 0
    do t = 1, ncbt
      tprow(t) = myrows
      if (towner(t) == me) then
        nmine = nmine + 1
        boff(t) = bandw
        hi = ctb(t) - ctb(t-1)
        bandw = bandw + int(hi, 8)*(ncol + ctb(t))
        myrows = myrows + hi
      endif
    enddo
    hasband = (nmine > 0)
    allocate(rowoff(nrow_nodes + 1), pos(sym%nnode))
    rowoff(1) = 0
    do i = 1, nrow_nodes
      k = sym%rlist(sym%rptr(s)+i-1)
      rowoff(i+1) = rowoff(i) + sym%ndof(k)
      if (i == nown) rowoff(i+1) = rowoff(i+1) + ndel
      pos(k) = i
    enddo

    ! (3) held storage: the master grid and my band, zeroed
    fw = 0
    if (ismaster) then
      allocate(fvm(gm%coloff(gm%nt+1)))
      if (fct%lu) allocate(fvmu(gm%coloff(gm%nt+1)))
      fw = fw + gm%coloff(gm%nt+1)
    endif
    if (hasband) then
      allocate(band(bandw))
      if (fct%lu) allocate(bandu(bandw))
      fw = fw + bandw
    endif
    if (fct%lu) fw = 2*fw
    fct%front_words_act = max(fct%front_words_act, fw)
    fct%live_front = fct%live_front + fw
    fct%front_peak = max(fct%front_peak, fct%live_front)
    !$omp parallel default(shared)
    !$omp sections
    !$omp section
    if (ismaster) fvm(:) = 0.0d0
    !$omp section
    if (ismaster .and. fct%lu) fvmu(:) = 0.0d0
    !$omp section
    if (hasband) band(:) = 0.0d0
    !$omp section
    if (hasband .and. fct%lu) bandu(:) = 0.0d0
    !$omp end sections
    !$omp end parallel

    ! (4) fully summed metadata on the master, as the sequential assembly builds it
    call fct_sn_meta()

    ! (5) scatter of the replicated matrix into what I hold
    call scatter_mine()

    ! (6) extend-add: slices out for the children I hold, then the slices of every child
    ! applied in chead order, the way the sequential extend-add consumes the blocks
    allocate(cctb(0:g%nt+2), cown(g%nt+2), pmap(max(nrow, 1)))
    do j = fct%cptr(s), fct%cptr(s+1) - 1
      c = fct%clist(j)
      call ea_setup(c, j - fct%cptr(s) + 1)
      call ea_send(c)
    enddo
    c = fct%chead(s)
    do while (c /= 0)
      do j = fct%cptr(s), fct%cptr(s+1) - 1
        if (fct%clist(j) == c) call ea_setup(c, j - fct%cptr(s) + 1)
      enddo
      call ea_recv_apply(c)
      c = fct%cnext(c)
    enddo

    ! (7) factorization of the fully summed part and (8) storage of what I hold
    nswap_f = 0
    n2x2_f = 0
    npos_f = 0
    nneg_f = 0
    ntl_f = 0
    nlr_f = 0
    rsum_f = 0
    rmax_f = 0
    npiv = 0
    if (ncbt == 0) then
      ! no contribution rows (a root front): the master factors the whole front with the
      ! sequential code; all tiles are fully summed, so no tile is ever compressed
      if (ismaster) call factor_root()
    else if (ismaster .or. hasband) then
      call factor_1d()
    endif
    if (gierr == 0 .and. (ismaster .or. hasband)) call store_1d()
    if (ismaster .and. gierr == 0) then
      fct%n_swap = fct%n_swap + nswap_f
      fct%n_2x2 = fct%n_2x2 + n2x2_f
      fct%n_delay = fct%n_delay + ncol - npiv
    endif
    fct%blr_tiles = fct%blr_tiles + ntl_f
    fct%blr_tiles_lr = fct%blr_tiles_lr + nlr_f
    fct%blr_rank_sum = fct%blr_rank_sum + rsum_f
    fct%blr_rank_max = max(fct%blr_rank_max, rmax_f)
    fct%n_pos = fct%n_pos + npos_f
    fct%n_neg = fct%n_neg + nneg_f

    call mf_pool_wait(pool)
    ! the contribution blocks of the children were consumed here and remotely
    do j = fct%cptr(s), fct%cptr(s+1) - 1
      c = fct%clist(j)
      if (allocated(fct%sn(c)%cval)) then
        fct%live_cb = fct%live_cb - fct%sn(c)%cbsize
        deallocate(fct%sn(c)%cval)
      endif
    enddo
    fct%live_front = fct%live_front - fw

  contains

    !> multicast the header of child c (which I hold) to the other ranks of the set
    subroutine send_hdr(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint) :: ib, dst, nt1, cnd, mm

      nt1 = fct%sn(c0)%nt + 1
      cnd = fct%sn(c0)%ncol - fct%sn(c0)%npiv
      mm = nt1 + 2*cnd
      ib = mf_pool_slot(pool)
      allocate(pool%box(ib)%hdr(5), pool%box(ib)%tl(mm))
      ! my local error rides in the header slot: the receivers and I then skip the front
      ! by the same rule, off the headers alone
      pool%box(ib)%hdr(1) = gierr
      if (gierr /= 0) ferr = gierr
      pool%box(ib)%hdr(2) = fct%sn(c0)%ncol
      pool%box(ib)%hdr(3) = fct%sn(c0)%npiv
      pool%box(ib)%hdr(4) = fct%sn(c0)%nt
      pool%box(ib)%hdr(5) = fct%sn(c0)%ntc
      pool%box(ib)%tl(1:nt1) = fct%sn(c0)%tbnd(1:nt1)
      pool%box(ib)%tl(nt1+1:nt1+cnd) = fct%sn(c0)%fsdof(fct%sn(c0)%npiv+1:fct%sn(c0)%ncol)
      pool%box(ib)%tl(nt1+cnd+1:nt1+2*cnd) = fct%sn(c0)%frow(fct%sn(c0)%npiv+1:fct%sn(c0)%ncol)
      do dst = map%rbeg(s), map%rbeg(s) + map%rcnt(s) - 1
        if (dst == me) cycle
        req = 0
        call hecmw_isend_int(pool%box(ib)%hdr, 5, dst, 16*c0, map%comm, req)
        call mf_pool_req(pool, req)
        call hecmw_isend_int(pool%box(ib)%tl, mm, dst, 16*c0+1, map%comm, req)
        call mf_pool_req(pool, req)
      enddo
    end subroutine send_hdr

    !> receive the header of child c into its snode
    subroutine recv_hdr(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint) :: hdr(5), nt1, cnd, mm
      integer(kind=kint), allocatable :: tl(:)

      call hecmw_recv_int(hdr, 5, map%owner(c0), 16*c0, map%comm, stat)
      if (hdr(1) /= 0) then
        gierr = hdr(1)
        ferr = hdr(1)
      endif
      fct%sn(c0)%ncol = hdr(2)
      fct%sn(c0)%npiv = hdr(3)
      fct%sn(c0)%nt = hdr(4)
      fct%sn(c0)%ntc = hdr(5)
      nt1 = hdr(4) + 1
      cnd = hdr(2) - hdr(3)
      call mf_grow_i(fct%sn(c0)%tbnd, nt1)
      call mf_grow_i(fct%sn(c0)%fsdof, max(hdr(2), 1))
      call mf_grow_i(fct%sn(c0)%frow, max(hdr(2), 1))
      mm = nt1 + 2*cnd
      allocate(tl(mm))
      call hecmw_recv_int(tl, mm, map%owner(c0), 16*c0+1, map%comm, stat)
      fct%sn(c0)%tbnd(1:nt1) = tl(1:nt1)
      fct%sn(c0)%fsdof(hdr(3)+1:hdr(2)) = tl(nt1+1:nt1+cnd)
      fct%sn(c0)%frow(hdr(3)+1:hdr(2)) = tl(nt1+cnd+1:nt1+2*cnd)
      deallocate(tl)
    end subroutine recv_hdr

    !> fully summed position metadata on the master (positions, DOFs, node blocks), as the
    !> sequential assembly initializes them; every rank builds blk for the root fallback
    subroutine fct_sn_meta()
      integer(kind=kint) :: mm, kk, ii, cc, jj0

      allocate(blk(max(ncol, 1)), blkr(max(ncol, 1)))
      fct%sn(s)%ncol = ncol
      fct%sn(s)%nt = g%nt
      fct%sn(s)%ntc = g%ntc
      call mf_grow_i(fct%sn(s)%tbnd, g%nt + 1)
      fct%sn(s)%tbnd(1:g%nt+1) = g%tb(0:g%nt)
      call mf_grow_i(fct%sn(s)%fsdof, max(ncol, 1))
      call mf_grow_i(fct%sn(s)%frow, max(ncol, 1))
      call mf_grow_i(fct%sn(s)%ptype, max(ncol, 1))
      call mf_grow_r(fct%sn(s)%dsub, int(max(ncol, 1), 8))
      mm = 0
      do kk = sym%sptr(s), sym%sptr(s+1) - 1
        do ii = fct%cdofptr(kk), fct%cdofptr(kk+1) - 1
          mm = mm + 1
          fct%sn(s)%fsdof(mm) = ii
          fct%sn(s)%frow(mm) = ii
          blk(mm) = kk
        enddo
      enddo
      cc = fct%chead(s)
      do while (cc /= 0)
        do jj0 = fct%sn(cc)%npiv + 1, fct%sn(cc)%ncol
          mm = mm + 1
          fct%sn(s)%fsdof(mm) = fct%sn(cc)%fsdof(jj0)
          fct%sn(s)%frow(mm) = fct%sn(cc)%frow(jj0)
          blk(mm) = 0
        enddo
        cc = fct%cnext(cc)
      enddo
      blkr(1:ncol) = blk(1:ncol)
    end subroutine fct_sn_meta

    !> owning rank of a front grid row
    function drank(row) result(dr)
      integer(kind=kint), intent(in) :: row
      integer(kind=kint) :: dr
      if (row <= ncol) then
        dr = master
      else
        dr = towner(crow2t(row - ncol))
      endif
    end function drank

    !> word index of the entry (row, col) of my band, row in my tile t0
    function bandix(t0, row, col) result(ix)
      integer(kind=kint), intent(in) :: t0, row, col
      integer(kind=8) :: ix
      ix = boff(t0) + int(col - 1, 8)*(ctb(t0) - ctb(t0-1)) + (row - ncol - ctb(t0-1))
    end function bandix

    !> scatter the permuted matrix entries whose front row I hold, following the
    !> sequential scatter (columns of the own nodes, lower part in LDLt, both in LU)
    subroutine scatter_mine()
      integer(kind=kint) :: ii, kk, j0, coff, aa, bb, ki, roff, kk2
      integer(kind=8) :: bs, bs2

      !$omp parallel do default(shared) private(kk, j0, coff, bs, bs2, aa, bb, kk2, ki, roff) schedule(dynamic, 4)
      do ii = 1, nown
        kk = sym%sptr(s) + ii - 1
        j0 = sym%perm(kk)
        coff = rowoff(ii)
        bs = int(j0-1, 8)*nd*nd
        do bb = 1, nd
          do aa = bb, nd
            call put1(coff+aa, coff+bb, hecMAT%D(bs + (aa-1)*nd + bb), .false.)
          enddo
          if (fct%lu) then
            do aa = 1, bb-1
              call put1(coff+bb, coff+aa, hecMAT%D(bs + (aa-1)*nd + bb), .true.)
            enddo
          endif
        enddo
        do kk2 = hecMAT%indexL(j0-1)+1, hecMAT%indexL(j0)
          ki = sym%invp(hecMAT%itemL(kk2))
          if (ki <= kk) cycle
          roff = rowoff(pos(ki))
          bs = int(kk2-1, 8)*nd*nd
          if (.not. fct%lu) then
            do bb = 1, nd
              do aa = 1, nd
                call put1(roff+aa, coff+bb, hecMAT%AL(bs + (bb-1)*nd + aa), .false.)
              enddo
            enddo
          else
            bs2 = int(fct%mirror(kk2)-1, 8)*nd*nd
            do aa = 1, nd
              do bb = 1, nd
                call put1(roff+bb, coff+aa, hecMAT%AL(bs + (aa-1)*nd + bb), .true.)
                call put1(roff+bb, coff+aa, hecMAT%AU(bs2 + (bb-1)*nd + aa), .false.)
              enddo
            enddo
          endif
        enddo
        do kk2 = hecMAT%indexU(j0-1)+1, hecMAT%indexU(j0)
          ki = sym%invp(hecMAT%itemU(kk2))
          if (ki <= kk) cycle
          roff = rowoff(pos(ki))
          bs = int(kk2-1, 8)*nd*nd
          if (.not. fct%lu) then
            do bb = 1, nd
              do aa = 1, nd
                call put1(roff+aa, coff+bb, hecMAT%AU(bs + (bb-1)*nd + aa), .false.)
              enddo
            enddo
          else
            bs2 = int(fct%mirroru(kk2)-1, 8)*nd*nd
            do aa = 1, nd
              do bb = 1, nd
                call put1(roff+bb, coff+aa, hecMAT%AU(bs + (aa-1)*nd + bb), .true.)
                call put1(roff+bb, coff+aa, hecMAT%AL(bs2 + (bb-1)*nd + aa), .false.)
              enddo
            enddo
          endif
        enddo
      enddo
      !$omp end parallel do
    end subroutine scatter_mine

    !> add v0 at front position (pr, pc) if I hold its row; upface selects the upper grid
    !> of the LU mode the way the sequential add_lu does
    subroutine put1(pr, pc, v0, upface)
      integer(kind=kint), intent(in) :: pr, pc
      real(kind=kreal), intent(in) :: v0
      logical, intent(in) :: upface
      integer(kind=kint) :: rr, cc, t0
      integer(kind=8) :: ix
      logical :: up

      rr = max(pr, pc)
      cc = min(pr, pc)
      up = fct%lu .and. (upface .neqv. (pr < pc))
      if (rr <= ncol) then
        if (.not. ismaster) return
        ix = mf_idx(gm, rr, cc)
        if (up) then
          fvmu(ix) = fvmu(ix) + v0
        else
          fvm(ix) = fvm(ix) + v0
        endif
      else
        t0 = crow2t(rr - ncol)
        if (towner(t0) /= me) return
        ix = bandix(t0, rr, cc)
        if (up) then
          bandu(ix) = bandu(ix) + v0
        else
          band(ix) = band(ix) + v0
        endif
      endif
    end subroutine put1

    !> block structure and row map of child c0 (clist ordinal jc): the contribution block
    !> boundaries, the owner of every block row, and the parent front position of every
    !> child contribution row
    subroutine ea_setup(c0, jc)
      integer(kind=kint), intent(in) :: c0, jc
      integer(kind=kint) :: jj0, rr, ll, kk, aa, ntc0, nde, nt2

      cndel = fct%sn(c0)%ncol - fct%sn(c0)%npiv
      cnb = 0
      cctb(0) = 0
      if (cndel > 0) then
        cnb = 1
        cctb(1) = cndel
        cown(1) = map%owner(c0)
      endif
      kbegc = cnb
      ntc0 = fct%sn(c0)%ntc
      do jj0 = ntc0 + 1, fct%sn(c0)%nt
        cnb = cnb + 1
        cctb(cnb) = cndel + fct%sn(c0)%tbnd(jj0+1) - fct%sn(c0)%ncol
      enddo
      if (map%upper(c0)) then
        nde = fct%sn(c0)%ncol - (fct%cdofptr(sym%sptr(c0+1)) - fct%cdofptr(sym%sptr(c0)))
        call hecmw_mf_dist_rowtiles(sym, map, fct%tile, c0, nde, nt2, cctb2, cown2)
        do jj0 = kbegc + 1, cnb
          cown(jj0) = cown2(jj0 - kbegc)
        enddo
      else
        do jj0 = kbegc + 1, cnb
          cown(jj0) = map%owner(c0)
        enddo
      endif
      cnrow = cctb(cnb)
      ! parent position of every child contribution row
      do rr = 1, cndel
        pmap(rr) = mdel(jc) + rr
      enddo
      rr = cndel
      do ll = sym%rptr(c0) + (sym%sptr(c0+1) - sym%sptr(c0)), sym%rptr(c0+1) - 1
        kk = sym%cmap(sym%cmap_ptr(c0) + (ll - sym%rptr(c0) - (sym%sptr(c0+1) - sym%sptr(c0))))
        do aa = 1, sym%ndof(sym%rlist(ll))
          rr = rr + 1
          pmap(rr) = rowoff(kk) + aa
        enddo
      enddo
    end subroutine ea_setup

    !> pack and isend my slices of the contribution block of child c0, one message per
    !> destination rank holding rows of the parent front
    subroutine ea_send(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint) :: dst, nlo, nup, ib

      if (.not. allocated(fct%sn(c0)%cval)) return
      do dst = map%rbeg(s), map%rbeg(s) + map%rcnt(s) - 1
        if (dst == me) cycle
        call ea_count(c0, me, dst, nlo, nup)
        if (nlo == 0) cycle
        ib = mf_pool_slot(pool)
        allocate(pool%box(ib)%rv(nlo + nup))
        call ea_pack(c0, me, dst, pool%box(ib)%rv, nlo)
        req = 0
        call hecmw_isend_r(pool%box(ib)%rv, nlo + nup, dst, 16*c0+2, map%comm, req)
        call mf_pool_req(pool, req)
      enddo
    end subroutine ea_send

    !> receive and apply the slices of child c0 from every rank holding parts of it; my
    !> own part is applied straight from the local storage
    subroutine ea_recv_apply(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint) :: is0, src, nlo, nup, nsrc, jj0
      integer(kind=kint) :: srcs(cnb + 1)
      real(kind=kreal), allocatable :: rb(:)

      ! distinct source ranks in block order
      nsrc = 0
      do jj0 = 1, cnb
        do is0 = 1, nsrc
          if (srcs(is0) == cown(jj0)) exit
        enddo
        if (is0 > nsrc) then
          nsrc = nsrc + 1
          srcs(nsrc) = cown(jj0)
        endif
      enddo
      do is0 = 1, nsrc
        src = srcs(is0)
        call ea_count(c0, src, me, nlo, nup)
        if (nlo == 0) cycle
        if (src == me) then
          call ea_apply_local(c0)
        else
          allocate(rb(nlo + nup))
          call hecmw_recv_r(rb, nlo + nup, src, 16*c0+2, map%comm, stat)
          call ea_apply(c0, src, rb, nlo)
          deallocate(rb)
        endif
      enddo
    end subroutine ea_recv_apply

    !> entries of child c0 held by src and destined to dst
    subroutine ea_count(c0, src, dst, nlo, nup)
      integer(kind=kint), intent(in) :: c0, src, dst
      integer(kind=kint), intent(out) :: nlo, nup
      integer(kind=kint) :: jj0, ii0, cc, rr, hr

      nlo = 0
      nup = 0
      do jj0 = 1, cnb
        do ii0 = jj0, cnb
          if (cown(ii0) /= src) cycle
          do cc = cctb(jj0-1) + 1, cctb(jj0)
            do rr = max(cc, cctb(ii0-1) + 1), cctb(ii0)
              hr = max(pmap(rr), pmap(cc))
              if (drank(hr) /= dst) cycle
              nlo = nlo + 1
              if (fct%lu .and. rr > cc) nup = nup + 1
            enddo
          enddo
        enddo
      enddo
    end subroutine ea_count

    !> pack my entries of child c0 destined to dst: the lower face values first, then the
    !> upper face values of the off-diagonal entries in the same order (LU)
    subroutine ea_pack(c0, src, dst, buf, nlo)
      integer(kind=kint), intent(in) :: c0, src, dst, nlo
      real(kind=kreal), intent(inout) :: buf(:)
      integer(kind=kint) :: jj0, ii0, cc, rr, hr, hh, il, iu
      integer(kind=8) :: sb, halfc, off

      halfc = ea_src_words(c0, src)
      il = 0
      iu = 0
      sb = 0
      do jj0 = 1, cnb
        do ii0 = jj0, cnb
          if (cown(ii0) /= src) cycle
          hh = cctb(ii0) - cctb(ii0-1)
          do cc = cctb(jj0-1) + 1, cctb(jj0)
            do rr = max(cc, cctb(ii0-1) + 1), cctb(ii0)
              hr = max(pmap(rr), pmap(cc))
              if (drank(hr) /= dst) cycle
              off = sb + int(cc - cctb(jj0-1) - 1, 8)*hh + (rr - cctb(ii0-1))
              il = il + 1
              buf(il) = fct%sn(c0)%cval(off)
              if (fct%lu .and. rr > cc) then
                iu = iu + 1
                buf(nlo + iu) = fct%sn(c0)%cval(halfc + off)
              endif
            enddo
          enddo
          sb = sb + int(cctb(jj0) - cctb(jj0-1), 8)*hh
        enddo
      enddo
    end subroutine ea_pack

    !> apply the received slice of child c0 from src into what I hold
    subroutine ea_apply(c0, src, buf, nlo)
      integer(kind=kint), intent(in) :: c0, src, nlo
      real(kind=kreal), intent(in) :: buf(:)
      integer(kind=kint) :: jj0, ii0, cc, rr, hr, il, iu

      il = 0
      iu = 0
      do jj0 = 1, cnb
        do ii0 = jj0, cnb
          if (cown(ii0) /= src) cycle
          do cc = cctb(jj0-1) + 1, cctb(jj0)
            do rr = max(cc, cctb(ii0-1) + 1), cctb(ii0)
              hr = max(pmap(rr), pmap(cc))
              if (drank(hr) /= me) cycle
              il = il + 1
              call put1(pmap(rr), pmap(cc), buf(il), .false.)
              if (fct%lu .and. rr > cc) then
                iu = iu + 1
                call put1(pmap(cc), pmap(rr), buf(nlo + iu), .false.)
              endif
            enddo
          enddo
        enddo
      enddo
    end subroutine ea_apply

    !> apply my own part of child c0 straight from its stored block
    subroutine ea_apply_local(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint) :: jj0, ii0, cc, rr, hr, hh
      integer(kind=8) :: sb, halfc, off

      halfc = ea_src_words(c0, me)
      sb = 0
      do jj0 = 1, cnb
        do ii0 = jj0, cnb
          if (cown(ii0) /= me) cycle
          hh = cctb(ii0) - cctb(ii0-1)
          do cc = cctb(jj0-1) + 1, cctb(jj0)
            do rr = max(cc, cctb(ii0-1) + 1), cctb(ii0)
              hr = max(pmap(rr), pmap(cc))
              if (drank(hr) /= me) cycle
              off = sb + int(cc - cctb(jj0-1) - 1, 8)*hh + (rr - cctb(ii0-1))
              call put1(pmap(rr), pmap(cc), fct%sn(c0)%cval(off), .false.)
              if (fct%lu .and. rr > cc) call put1(pmap(cc), pmap(rr), fct%sn(c0)%cval(halfc + off), .false.)
            enddo
          enddo
          sb = sb + int(cctb(jj0) - cctb(jj0-1), 8)*hh
        enddo
      enddo
    end subroutine ea_apply_local

    !> words of one face of the blocks of child c0 held by src
    function ea_src_words(c0, src) result(ww)
      integer(kind=kint), intent(in) :: c0, src
      integer(kind=8) :: ww
      integer(kind=kint) :: jj0, ii0

      ww = 0
      do jj0 = 1, cnb
        do ii0 = jj0, cnb
          if (cown(ii0) /= src) cycle
          ww = ww + int(cctb(ii0) - cctb(ii0-1), 8)*(cctb(jj0) - cctb(jj0-1))
        enddo
      enddo
    end function ea_src_words

    !> a root front has no contribution rows: the master factors it whole with the
    !> sequential code (the fully summed grid is the full grid there)
    subroutine factor_root()
      type(mf_work) :: lw
      integer(kind=kint) :: ie

      call mf_grow_i(lw%blk, max(ncol, 1))
      call mf_grow_i(lw%blkr, max(ncol, 1))
      lw%blk(1:ncol) = blk(1:ncol)
      lw%blkr(1:ncol) = blkr(1:ncol)
      ie = 0
      !$omp parallel default(shared)
      !$omp single
      if (fct%lu) then
        call mf_factor_front_lu(fct%sn(s), gm, fvm, fvmu, lw%pval, lw%pvalu, lw%wk, lw%blk, lw%blkr, &
          sym%sparent(s) == 0, fct%pivot_u, zero, fct%blr, fct%eps, lw%bval, lw%boff, lw%brk, &
          lw%bvalu, lw%boffu, lw%brku, nswap_f, ntl_f, nlr_f, rsum_f, rmax_f, ie)
      else
        call mf_factor_front(fct%sn(s), gm, fvm, lw%pval, lw%wval, lw%wk, lw%blk, &
          sym%sparent(s) == 0, fct%pivot_u, zero, fct%blr, fct%eps, lw%bval, lw%boff, lw%brk, &
          nswap_f, n2x2_f, npos_f, nneg_f, ntl_f, nlr_f, rsum_f, rmax_f, ie)
        fct%sn(s)%frow(1:ncol) = fct%sn(s)%fsdof(1:ncol)
      endif
      !$omp end single
      !$omp end parallel
      if (ie /= 0) gierr = ie
      npiv = fct%sn(s)%npiv
    end subroutine factor_root

    !> the tile column protocol of the distributed fully summed factorization
    subroutine factor_1d()
      integer(kind=kint) :: kk, iw2, mm

      wme = 0
      do iw2 = 1, nwk
        if (wrank(iw2) == me) wme = iw2
      enddo
      if (fct%blr) then
        mm = mf_bidx(g%nt, g%ntc, g%nt)
        allocate(brk(mm), bslot(mm))
        brk(1:mm) = -1
        if (fct%lu) then
          allocate(brku(mm), bslotu(mm))
          brku(1:mm) = -1
        endif
      endif
      btop = 0
      btopu = 0
      allocate(perm(max(ncol, 1)), permr(max(ncol, 1)), itmp(4*max(ncol, 1)), swaps(2*max(ncol, 1)))
      allocate(cmax(max(ncol, 1)), ptv(max(ncol, 1)), dsv(max(ncol, 1)), dgv(max(ncol, 1)), dum(1))
      allocate(pair_i(g%nt*(g%nt+1)/2), pair_j(g%nt*(g%nt+1)/2))
      call mf_grow_r(wk, int(2*nrow, 8))
      call mf_grow_r(pv, int(max(ncol, 1), 8))
      if (.not. fct%lu) then
        ! the upper grid buffers stay tiny but allocated, so they can pass as arguments
        call mf_grow_r(pvu, 1_8)
        call mf_grow_r(pbandu, 1_8)
        call mf_grow_r(pfullu, 1_8)
        call mf_grow_r(tmatu, 1_8)
        call mf_grow_r(dgbu, 1_8)
        call mf_grow_r(fslu, 1_8)
        call mf_grow_r(tvecu, 1_8)
      endif
      if (ismaster) then
        if (fct%lu) then
          fct%sn(s)%ptype(1:ncol) = 1
        else
          fct%sn(s)%ptype(1:ncol) = 0
        endif
        fct%sn(s)%dsub(1:ncol) = 0.0d0
      endif
      npiv = 0
      nfs = ncol
      do kk = 1, g%ntc
        p0 = g%tb(kk-1)
        if (p0 >= nfs) exit
        pa = p0 + 1
        do while (pa <= min(g%tb(kk), nfs))
          pb = min(g%tb(kk), nfs)
          w = pb - pa + 1
          pk0 = npiv - p0
          call panel_attempt()
          if (acc) then
            npiv = npiv + w
            pa = pb + 1
          else
            call panel_piv_1d()
            npiv = npiv + np
            pa = pa + np
          endif
        enddo
        pk = npiv - p0
        call postcol(kk)
        if (nfs <= g%tb(kk)) exit
      enddo
      if (ismaster) then
        if (.not. fct%lu) then
          fct%sn(s)%frow(1:ncol) = fct%sn(s)%fsdof(1:ncol)
          call inertia_1d()
        endif
      endif
      fct%sn(s)%npiv = npiv
    end subroutine factor_1d

    !> update vectors of the pivots already eliminated in this tile column, for the panel
    !> columns pa..pb (what the sequential fill_panel computes on the fly)
    subroutine build_tmat()
      integer(kind=kint) :: jj2, xx, q
      real(kind=kreal) :: l1, l2, d11, d21, d22

      if (pk0 <= 0) return
      call mf_grow_r(tmat, int(w, 8)*pk0)
      if (fct%lu) call mf_grow_r(tmatu, int(w, 8)*pk0)
      do jj2 = 1, w
        xx = pa + jj2 - 1
        if (.not. fct%lu) then
          q = p0 + 1
          do while (q <= npiv)
            if (fct%sn(s)%ptype(q) == 1) then
              tmat((jj2-1)*pk0 + (q-p0)) = fvm(mf_idx(gm, xx, q)) * fvm(mf_idx(gm, q, q))
              q = q + 1
            else
              l1 = fvm(mf_idx(gm, xx, q))
              l2 = fvm(mf_idx(gm, xx, q+1))
              d11 = fvm(mf_idx(gm, q, q))
              d21 = fct%sn(s)%dsub(q)
              d22 = fvm(mf_idx(gm, q+1, q+1))
              tmat((jj2-1)*pk0 + (q-p0)) = l1*d11 + l2*d21
              tmat((jj2-1)*pk0 + (q-p0+1)) = l1*d21 + l2*d22
              q = q + 2
            endif
          enddo
        else
          call mf_row_get(gm, fvmu, xx, p0+1, npiv, tmat((jj2-1)*pk0+1))
          call mf_row_get(gm, fvm, xx, p0+1, npiv, tmatu((jj2-1)*pk0+1))
        endif
      enddo
    end subroutine build_tmat

    !> the fully summed rows of the panel columns pa..pb into pv (and pvu)
    subroutine fill_master_panel()
      integer(kind=kint) :: jj2, xx

      call mf_grow_r(pv, int(mfs, 8)*w)
      if (fct%lu) call mf_grow_r(pvu, int(mfs, 8)*w)
      do jj2 = 1, w
        xx = pa + jj2 - 1
        call mf_col_copy(gm, fvm, xx, xx, pv(int(jj2-1, 8)*mfs + jj2), .false.)
        if (pk0 > 0) call mf_col_axpy(gm, fvm, xx, p0+1, npiv, tmat((jj2-1)*pk0+1), pv(int(jj2-1, 8)*mfs + jj2))
        if (fct%lu .and. xx < ncol) then
          call mf_col_copy(gm, fvmu, xx, xx+1, pvu(int(jj2-1, 8)*mfs + jj2 + 1), .false.)
          if (pk0 > 0) &
            call mf_col_axpy(gm, fvmu, xx+1, p0+1, npiv, tmatu((jj2-1)*pk0+1), pvu(int(jj2-1, 8)*mfs + jj2 + 1))
        endif
      enddo
    end subroutine fill_master_panel

    !> my band rows of the panel columns pa..pb into pband (and pbandu), with the tile
    !> column updates; solveit applies the panel solve and flags the columns whose
    !> threshold fails (cmax 1.0), as the sequential check does over these rows
    subroutine fill_band_panel(solveit)
      logical, intent(in) :: solveit
      integer(kind=kint) :: jj2, t2, h2, r0, i2
      integer(kind=8) :: bx, px

      call mf_grow_r(pband, int(myrows, 8)*w)
      if (fct%lu) call mf_grow_r(pbandu, int(myrows, 8)*w)
      do jj2 = 1, w
        do t2 = 1, ncbt
          if (towner(t2) /= me) cycle
          h2 = ctb(t2) - ctb(t2-1)
          r0 = tprow(t2)
          bx = boff(t2) + int(pa + jj2 - 2, 8)*h2
          px = int(jj2-1, 8)*myrows + r0
          pband(px+1:px+h2) = band(bx+1:bx+h2)
          if (fct%lu) pbandu(px+1:px+h2) = bandu(bx+1:bx+h2)
          if (pk0 > 0) then
            call hecmw_mf_kernel_gemv(h2, pk0, h2, band(boff(t2) + int(p0, 8)*h2 + 1), &
              tmat((jj2-1)*pk0+1), pband(px+1))
            if (fct%lu) call hecmw_mf_kernel_gemv(h2, pk0, h2, bandu(boff(t2) + int(p0, 8)*h2 + 1), &
              tmatu((jj2-1)*pk0+1), pbandu(px+1))
          endif
        enddo
      enddo
      if (solveit) then
        if (.not. fct%lu) then
          call hecmw_mf_kernel_trsm(myrows, w, w, dgb, myrows, pband)
        else
          call hecmw_mf_kernel_trsm_rt(myrows, w, w, dgbu, .false., myrows, pband)
          call hecmw_mf_kernel_trsm_rt(myrows, w, w, dgb, .true., myrows, pbandu)
        endif
        do jj2 = 1, w
          cmax(jj2) = 0.0d0
          px = int(jj2-1, 8)*myrows
          do i2 = 1, myrows
            if (.not. (abs(pband(px+i2)) <= uinv)) cmax(jj2) = 1.0d0
          enddo
        enddo
      endif
    end subroutine fill_band_panel

    !> the solved panel columns pa..pb back into my band
    subroutine band_writeback()
      integer(kind=kint) :: jj2, t2, h2, r0
      integer(kind=8) :: bx, px

      do jj2 = 1, w
        do t2 = 1, ncbt
          if (towner(t2) /= me) cycle
          h2 = ctb(t2) - ctb(t2-1)
          r0 = tprow(t2)
          bx = boff(t2) + int(pa + jj2 - 2, 8)*h2
          px = int(jj2-1, 8)*myrows + r0
          band(bx+1:bx+h2) = pband(px+1:px+h2)
          if (fct%lu) bandu(bx+1:bx+h2) = pbandu(px+1:px+h2)
        enddo
      enddo
    end subroutine band_writeback

    !> exchange the columns x0 and y0 of my band (bandu when ub)
    subroutine band_swap_cols(x0, y0, ub)
      integer(kind=kint), intent(in) :: x0, y0
      logical, intent(in) :: ub
      integer(kind=kint) :: t2, h2, i2
      integer(kind=8) :: ax, bx
      real(kind=kreal) :: v0

      do t2 = 1, ncbt
        if (towner(t2) /= me) cycle
        h2 = ctb(t2) - ctb(t2-1)
        ax = boff(t2) + int(x0-1, 8)*h2
        bx = boff(t2) + int(y0-1, 8)*h2
        if (.not. ub) then
          do i2 = 1, h2
            v0 = band(ax+i2)
            band(ax+i2) = band(bx+i2)
            band(bx+i2) = v0
          enddo
        else
          do i2 = 1, h2
            v0 = bandu(ax+i2)
            bandu(ax+i2) = bandu(bx+i2)
            bandu(bx+i2) = v0
          enddo
        endif
      enddo
    end subroutine band_swap_cols

    !> isend n0 control ints to every worker except me
    subroutine ctl_send_i(vals, n0)
      integer(kind=kint), intent(in) :: vals(:)
      integer(kind=kint), intent(in) :: n0
      integer(kind=kint) :: ib, iw2

      ib = mf_pool_slot(pool)
      allocate(pool%box(ib)%hdr(n0))
      pool%box(ib)%hdr(1:n0) = vals(1:n0)
      do iw2 = 1, nwk
        if (wrank(iw2) == me) cycle
        req = 0
        call hecmw_isend_int(pool%box(ib)%hdr, n0, wrank(iw2), 16*s+3, map%comm, req)
        call mf_pool_req(pool, req)
      enddo
    end subroutine ctl_send_i

    !> isend the boxed reals of slot ib to every worker except me
    subroutine ctl_send_box(ib, n0)
      integer(kind=kint), intent(in) :: ib, n0
      integer(kind=kint) :: iw2

      do iw2 = 1, nwk
        if (wrank(iw2) == me) cycle
        req = 0
        call hecmw_isend_r(pool%box(ib)%rv, n0, wrank(iw2), 16*s+4, map%comm, req)
        call mf_pool_req(pool, req)
      enddo
    end subroutine ctl_send_box

    !> panel of the positions pa..pb attempted without pivoting: the master factors the
    !> fully summed rows and broadcasts the update vectors and the factored diagonal
    !> blocks, every tile owner solves and checks its rows, and one flag reduction over
    !> the rank set decides the acceptance, reproducing the sequential accept/reject
    subroutine panel_attempt()
      real(kind=kreal), allocatable :: rcm(:)
      integer(kind=kint) :: jj2, iw2, mm, ib, fac
      integer(kind=8) :: ox

      fac = 1
      if (fct%lu) fac = 2
      mfs = ncol - pa + 1
      fsok = .false.
      call mf_grow_r(dgb, int(w, 8)*w)
      if (fct%lu) call mf_grow_r(dgbu, int(w, 8)*w)
      if (ismaster) then
        call build_tmat()
        call fill_master_panel()
        if (.not. fct%lu) then
          call hecmw_mf_kernel_panel_nopiv(mfs, w, mfs, pv, fct%pivot_u, zero, info)
        else
          call hecmw_mf_kernel_panel_lu_nopiv(mfs, w, mfs, pv, mfs, pvu, fct%pivot_u, zero, info)
          if (info == 0 .and. mfs == w) then
            do jj2 = 1, w
              pvu(int(jj2-1, 8)*mfs + jj2) = pv(int(jj2-1, 8)*mfs + jj2)
            enddo
          endif
        endif
        fsok = (info == 0)
        if (fsok) then
          do jj2 = 1, w
            dgb((jj2-1)*w+1:jj2*w) = pv(int(jj2-1, 8)*mfs + 1 : int(jj2-1, 8)*mfs + w)
            if (fct%lu) dgbu((jj2-1)*w+1:jj2*w) = pvu(int(jj2-1, 8)*mfs + 1 : int(jj2-1, 8)*mfs + w)
          enddo
        endif
        ctli(1:8) = 0
        ctli(1) = 1
        ctli(2) = pa
        ctli(3) = pb
        ctli(4) = npiv
        ctli(5) = merge(1, 0, fsok)
        call ctl_send_i(ctli, 8)
        mm = fac*w*pk0
        if (fsok) mm = mm + fac*w*w
        if (mm > 0 .and. nwk > merge(1, 0, wme > 0)) then
          ib = mf_pool_slot(pool)
          allocate(pool%box(ib)%rv(mm))
          ox = 0
          if (pk0 > 0) then
            pool%box(ib)%rv(ox+1:ox+int(w, 8)*pk0) = tmat(1:int(w, 8)*pk0)
            ox = ox + int(w, 8)*pk0
            if (fct%lu) then
              pool%box(ib)%rv(ox+1:ox+int(w, 8)*pk0) = tmatu(1:int(w, 8)*pk0)
              ox = ox + int(w, 8)*pk0
            endif
          endif
          if (fsok) then
            pool%box(ib)%rv(ox+1:ox+int(w, 8)*w) = dgb(1:int(w, 8)*w)
            ox = ox + int(w, 8)*w
            if (fct%lu) then
              pool%box(ib)%rv(ox+1:ox+int(w, 8)*w) = dgbu(1:int(w, 8)*w)
            endif
          endif
          call ctl_send_box(ib, mm)
        endif
      else
        call hecmw_recv_int(ctli, 8, master, 16*s+3, map%comm, stat)
        pa = ctli(2)
        pb = ctli(3)
        w = pb - pa + 1
        pk0 = npiv - p0
        fsok = (ctli(5) == 1)
        mm = fac*w*pk0
        if (fsok) mm = mm + fac*w*w
        if (mm > 0) then
          allocate(rcm(mm))
          call hecmw_recv_r(rcm, mm, master, 16*s+4, map%comm, stat)
          ox = 0
          if (pk0 > 0) then
            call mf_grow_r(tmat, int(w, 8)*pk0)
            tmat(1:int(w, 8)*pk0) = rcm(ox+1:ox+int(w, 8)*pk0)
            ox = ox + int(w, 8)*pk0
            if (fct%lu) then
              call mf_grow_r(tmatu, int(w, 8)*pk0)
              tmatu(1:int(w, 8)*pk0) = rcm(ox+1:ox+int(w, 8)*pk0)
              ox = ox + int(w, 8)*pk0
            endif
          endif
          if (fsok) then
            dgb(1:int(w, 8)*w) = rcm(ox+1:ox+int(w, 8)*w)
            ox = ox + int(w, 8)*w
            if (fct%lu) dgbu(1:int(w, 8)*w) = rcm(ox+1:ox+int(w, 8)*w)
          endif
          deallocate(rcm)
        endif
      endif
      ! the threshold over the distributed rows: fail flags of my rows, one reduction
      acc = .false.
      if (fsok) then
        if (hasband) then
          call fill_band_panel(.true.)
        else
          cmax(1:w) = 0.0d0
        endif
        if (ismaster) then
          acc = .true.
          do jj2 = 1, w
            if (cmax(jj2) /= 0.0d0) acc = .false.
          enddo
          allocate(rcm(w))
          do iw2 = 1, nwk
            if (wrank(iw2) == me) cycle
            call hecmw_recv_r(rcm, w, wrank(iw2), 16*s+5, map%comm, stat)
            do jj2 = 1, w
              if (rcm(jj2) /= 0.0d0) acc = .false.
            enddo
          enddo
          deallocate(rcm)
          ctli(1:8) = 0
          ctli(1) = 2
          ctli(2) = merge(1, 0, acc)
          call ctl_send_i(ctli, 8)
        else
          ib = mf_pool_slot(pool)
          allocate(pool%box(ib)%rv(w))
          pool%box(ib)%rv(1:w) = cmax(1:w)
          req = 0
          call hecmw_isend_r(pool%box(ib)%rv, w, master, 16*s+5, map%comm, req)
          call mf_pool_req(pool, req)
          call hecmw_recv_int(ctli, 8, master, 16*s+3, map%comm, stat)
          acc = (ctli(2) == 1)
        endif
      endif
      if (acc) then
        if (ismaster) then
          do jj2 = 1, w
            call mf_col_copy(gm, fvm, pa+jj2-1, pa+jj2-1, pv(int(jj2-1, 8)*mfs + jj2), .true.)
            if (fct%lu .and. pa+jj2-1 < ncol) &
              call mf_col_copy(gm, fvmu, pa+jj2-1, pa+jj2, pvu(int(jj2-1, 8)*mfs + jj2 + 1), .true.)
          enddo
          if (.not. fct%lu) fct%sn(s)%ptype(pa:pb) = 1
        endif
        if (hasband) call band_writeback()
      endif
    end subroutine panel_attempt

    !> the rejected panel gathered whole to the master, which runs the sequential pivoting
    !> kernel, permutes its fully summed grid and scatters the eliminated columns and the
    !> exchanges back to the tile owners
    subroutine panel_piv_1d()
      real(kind=kreal), allocatable :: rw(:)
      integer(kind=kint) :: jj2, iw2, t2, h2, x2, q2, mm, ib, mrows, fs1
      integer(kind=8) :: px, qx

      mrows = nrow - pa + 1
      fs1 = ncol - pa + 1
      if (hasband) call fill_band_panel(.false.)
      if (ismaster) then
        call mf_grow_r(pfull, int(mrows, 8)*w)
        if (fct%lu) call mf_grow_r(pfullu, int(mrows, 8)*w)
        call fill_master_panel()
        do jj2 = 1, w
          px = int(jj2-1, 8)*mrows
          pfull(px+jj2:px+fs1) = pv(int(jj2-1, 8)*mfs + jj2 : int(jj2-1, 8)*mfs + fs1)
          if (fct%lu) then
            if (pa+jj2-1 < ncol) pfullu(px+jj2+1:px+fs1) = pvu(int(jj2-1, 8)*mfs + jj2 + 1 : int(jj2-1, 8)*mfs + fs1)
          endif
        enddo
        if (hasband) call gather_rows(pband, pbandu, me)
        do iw2 = 1, nwk
          if (wrank(iw2) == me) cycle
          mm = wrows(wrank(iw2))
          allocate(rw(merge(2, 1, fct%lu)*mm*w))
          call hecmw_recv_r(rw, merge(2, 1, fct%lu)*mm*w, wrank(iw2), 16*s+5, map%comm, stat)
          call gather_rows_buf(rw, mm, wrank(iw2))
          deallocate(rw)
        enddo
        if (.not. fct%lu) then
          call hecmw_mf_kernel_panel_piv(mrows, w, mrows, pfull, blk(pa:pb), fct%pivot_u, MF_PIVOT_ALPHA, &
            zero, .true., wk, np, perm, fct%sn(s)%ptype(pa:pb), fct%sn(s)%dsub(pa:pb), nsw, n22, info)
          n2x2_f = n2x2_f + n22
        else
          call hecmw_mf_kernel_panel_lu_piv(mrows, w, mrows, pfull, mrows, pfullu, blk(pa:pb), blkr(pa:pb), &
            fct%pivot_u, zero, .true., np, perm, permr, nsw, info)
        endif
        nswap_f = nswap_f + nsw
        call master_perm_dance()
        do jj2 = 1, np
          px = int(jj2-1, 8)*mrows
          call mf_col_copy(gm, fvm, pa+jj2-1, pa+jj2-1, pfull(px + jj2), .true.)
          if (fct%lu .and. pa+jj2-1 < ncol) &
            call mf_col_copy(gm, fvmu, pa+jj2-1, pa+jj2, pfullu(px + jj2 + 1), .true.)
        enddo
        npair = 0
        x2 = pa + np
        do while (x2 <= pb)
          if (nfs > pb) then
            if (.not. fct%lu) then
              call mf_swap(gm, fvm, x2, nfs)
            else
              call mf_swap_lu(gm, fvm, fvmu, x2, nfs)
            endif
            call swap_meta(x2, nfs)
            npair = npair + 1
            swaps(2*npair-1) = x2
            swaps(2*npair) = nfs
            nfs = nfs - 1
            x2 = x2 + 1
          else
            nfs = x2 - 1
            exit
          endif
        enddo
        ctli(1:8) = 0
        ctli(1) = 4
        ctli(2) = np
        ctli(3) = npair
        ctli(4) = nfs
        call ctl_send_i(ctli, 8)
        mm = merge(2, 1, fct%lu)*w + 2*npair
        ib = mf_pool_slot(pool)
        allocate(pool%box(ib)%hdr(mm))
        pool%box(ib)%hdr(1:w) = perm(1:w)
        q2 = w
        if (fct%lu) then
          pool%box(ib)%hdr(q2+1:q2+w) = permr(1:w)
          q2 = q2 + w
        endif
        pool%box(ib)%hdr(q2+1:q2+2*npair) = swaps(1:2*npair)
        do iw2 = 1, nwk
          if (wrank(iw2) == me) cycle
          req = 0
          call hecmw_isend_int(pool%box(ib)%hdr, mm, wrank(iw2), 16*s+3, map%comm, req)
          call mf_pool_req(pool, req)
        enddo
        if (np > 0) then
          do iw2 = 1, nwk
            if (wrank(iw2) == me) cycle
            mm = wrows(wrank(iw2))
            if (mm == 0) cycle
            ib = mf_pool_slot(pool)
            allocate(pool%box(ib)%rv(merge(2, 1, fct%lu)*int(mm, 8)*np))
            call scatter_rows_buf(pool%box(ib)%rv, mm, wrank(iw2))
            req = 0
            call hecmw_isend_r(pool%box(ib)%rv, merge(2, 1, fct%lu)*mm*np, wrank(iw2), 16*s+4, map%comm, req)
            call mf_pool_req(pool, req)
          enddo
        endif
        if (hasband) then
          call mf_grow_r(pband, int(myrows, 8)*max(np, 1))
          if (fct%lu) call mf_grow_r(pbandu, int(myrows, 8)*max(np, 1))
          do jj2 = 1, np
            do t2 = 1, ncbt
              if (towner(t2) /= me) cycle
              h2 = ctb(t2) - ctb(t2-1)
              px = int(jj2-1, 8)*myrows + tprow(t2)
              qx = int(jj2-1, 8)*mrows + fs1 + ctb(t2-1)
              pband(px+1:px+h2) = pfull(qx+1:qx+h2)
              if (fct%lu) pbandu(px+1:px+h2) = pfullu(qx+1:qx+h2)
            enddo
          enddo
          call band_apply_piv()
        endif
      else
        ! send my raw rows, then apply the master's result to my band
        ib = mf_pool_slot(pool)
        mm = merge(2, 1, fct%lu)*myrows*w
        allocate(pool%box(ib)%rv(mm))
        pool%box(ib)%rv(1:int(myrows, 8)*w) = pband(1:int(myrows, 8)*w)
        if (fct%lu) pool%box(ib)%rv(int(myrows, 8)*w+1:mm) = pbandu(1:int(myrows, 8)*w)
        req = 0
        call hecmw_isend_r(pool%box(ib)%rv, mm, master, 16*s+5, map%comm, req)
        call mf_pool_req(pool, req)
        call hecmw_recv_int(ctli, 8, master, 16*s+3, map%comm, stat)
        np = ctli(2)
        npair = ctli(3)
        mm = merge(2, 1, fct%lu)*w + 2*npair
        call hecmw_recv_int(itmp, mm, master, 16*s+3, map%comm, stat)
        perm(1:w) = itmp(1:w)
        q2 = w
        if (fct%lu) then
          permr(1:w) = itmp(q2+1:q2+w)
          q2 = q2 + w
        endif
        swaps(1:2*npair) = itmp(q2+1:q2+2*npair)
        if (np > 0 .and. myrows > 0) then
          mm = merge(2, 1, fct%lu)*myrows*np
          allocate(rw(mm))
          call hecmw_recv_r(rw, mm, master, 16*s+4, map%comm, stat)
          call mf_grow_r(pband, int(myrows, 8)*max(np, 1))
          pband(1:int(myrows, 8)*np) = rw(1:int(myrows, 8)*np)
          if (fct%lu) then
            call mf_grow_r(pbandu, int(myrows, 8)*max(np, 1))
            pbandu(1:int(myrows, 8)*np) = rw(int(myrows, 8)*np+1:mm)
          endif
          deallocate(rw)
        endif
        call band_apply_piv()
        nfs = ctli(4)
      endif
    end subroutine panel_piv_1d

    !> rows of my band from pband into the gathered panel (the master's own part)
    subroutine gather_rows(pb0, pb0u, rk)
      real(kind=kreal), intent(in) :: pb0(:), pb0u(:)
      integer(kind=kint), intent(in) :: rk
      integer(kind=kint) :: jj2, t2, h2, r0, mrows, fs1
      integer(kind=8) :: px, qx

      mrows = nrow - pa + 1
      fs1 = ncol - pa + 1
      r0 = 0
      do t2 = 1, ncbt
        if (towner(t2) /= rk) cycle
        h2 = ctb(t2) - ctb(t2-1)
        do jj2 = 1, w
          px = int(jj2-1, 8)*myrows + r0
          qx = int(jj2-1, 8)*mrows + fs1 + ctb(t2-1)
          pfull(qx+1:qx+h2) = pb0(px+1:px+h2)
          if (fct%lu) pfullu(qx+1:qx+h2) = pb0u(px+1:px+h2)
        enddo
        r0 = r0 + h2
      enddo
    end subroutine gather_rows

    !> rows of the worker rk from its message into the gathered panel
    subroutine gather_rows_buf(rw, nr0, rk)
      real(kind=kreal), intent(in) :: rw(:)
      integer(kind=kint), intent(in) :: nr0, rk
      integer(kind=kint) :: jj2, t2, h2, r0, mrows, fs1
      integer(kind=8) :: px, qx

      mrows = nrow - pa + 1
      fs1 = ncol - pa + 1
      r0 = 0
      do t2 = 1, ncbt
        if (towner(t2) /= rk) cycle
        h2 = ctb(t2) - ctb(t2-1)
        do jj2 = 1, w
          px = int(jj2-1, 8)*nr0 + r0
          qx = int(jj2-1, 8)*mrows + fs1 + ctb(t2-1)
          pfull(qx+1:qx+h2) = rw(px+1:px+h2)
          if (fct%lu) pfullu(qx+1:qx+h2) = rw(int(nr0, 8)*w + px + 1 : int(nr0, 8)*w + px + h2)
        enddo
        r0 = r0 + h2
      enddo
    end subroutine gather_rows_buf

    !> rows of the worker rk of the eliminated panel columns into its scatter message
    subroutine scatter_rows_buf(rw, nr0, rk)
      real(kind=kreal), intent(inout) :: rw(:)
      integer(kind=kint), intent(in) :: nr0, rk
      integer(kind=kint) :: jj2, t2, h2, r0, mrows, fs1
      integer(kind=8) :: px, qx

      mrows = nrow - pa + 1
      fs1 = ncol - pa + 1
      r0 = 0
      do t2 = 1, ncbt
        if (towner(t2) /= rk) cycle
        h2 = ctb(t2) - ctb(t2-1)
        do jj2 = 1, np
          px = int(jj2-1, 8)*nr0 + r0
          qx = int(jj2-1, 8)*mrows + fs1 + ctb(t2-1)
          rw(px+1:px+h2) = pfull(qx+1:qx+h2)
          if (fct%lu) rw(int(nr0, 8)*np + px + 1 : int(nr0, 8)*np + px + h2) = pfullu(qx+1:qx+h2)
        enddo
        r0 = r0 + h2
      enddo
    end subroutine scatter_rows_buf

    !> rows of the band held by rank rk
    function wrows(rk) result(nr0)
      integer(kind=kint), intent(in) :: rk
      integer(kind=kint) :: nr0, t2

      nr0 = 0
      do t2 = 1, ncbt
        if (towner(t2) == rk) nr0 = nr0 + ctb(t2) - ctb(t2-1)
      enddo
    end function wrows

    !> the panel permutation applied to the master's fully summed grid and the position
    !> metadata, as the sequential permute_rows / permute_front do
    subroutine master_perm_dance()
      integer(kind=kint) :: jj2, t2, q2

      do jj2 = 1, w
        itmp(jj2) = jj2
      enddo
      do jj2 = 1, w
        if (itmp(jj2) == perm(jj2)) cycle
        do t2 = jj2 + 1, w
          if (itmp(t2) == perm(jj2)) exit
        enddo
        if (.not. fct%lu) then
          call mf_swap(gm, fvm, pa+jj2-1, pa+t2-1)
          q2 = fct%sn(s)%fsdof(pa+jj2-1)
          fct%sn(s)%fsdof(pa+jj2-1) = fct%sn(s)%fsdof(pa+t2-1)
          fct%sn(s)%fsdof(pa+t2-1) = q2
        else
          call mf_swap_lu(gm, fvm, fvmu, pa+jj2-1, pa+t2-1)
          call swap_meta(pa+jj2-1, pa+t2-1)
        endif
        itmp(t2) = itmp(jj2)
        itmp(jj2) = perm(jj2)
      enddo
      if (fct%lu) then
        itmp(1:w) = perm(1:w)
        do jj2 = 1, w
          if (itmp(jj2) == permr(jj2)) cycle
          do t2 = jj2 + 1, w
            if (itmp(t2) == permr(jj2)) exit
          enddo
          call mf_swap_row(gm, fvm, fvmu, pa+jj2-1, pa+t2-1)
          q2 = fct%sn(s)%frow(pa+jj2-1)
          fct%sn(s)%frow(pa+jj2-1) = fct%sn(s)%frow(pa+t2-1)
          fct%sn(s)%frow(pa+t2-1) = q2
          q2 = blkr(pa+jj2-1)
          blkr(pa+jj2-1) = blkr(pa+t2-1)
          blkr(pa+t2-1) = q2
          itmp(t2) = itmp(jj2)
          itmp(jj2) = permr(jj2)
        enddo
      endif
    end subroutine master_perm_dance

    !> exchange of the position metadata x0 <-> y0 (fsdof and blk; frow and blkr in LU)
    subroutine swap_meta(x0, y0)
      integer(kind=kint), intent(in) :: x0, y0
      integer(kind=kint) :: q2

      q2 = fct%sn(s)%fsdof(x0)
      fct%sn(s)%fsdof(x0) = fct%sn(s)%fsdof(y0)
      fct%sn(s)%fsdof(y0) = q2
      q2 = blk(x0)
      blk(x0) = blk(y0)
      blk(y0) = q2
      if (fct%lu) then
        q2 = fct%sn(s)%frow(x0)
        fct%sn(s)%frow(x0) = fct%sn(s)%frow(y0)
        fct%sn(s)%frow(y0) = q2
        q2 = blkr(x0)
        blkr(x0) = blkr(y0)
        blkr(y0) = q2
      endif
    end subroutine swap_meta

    !> the panel permutation and the delayed exchanges applied to my band columns, then
    !> the eliminated columns overwritten with the scattered rows
    subroutine band_apply_piv()
      integer(kind=kint) :: jj2, t2, h2, i2
      integer(kind=8) :: px, bx

      do jj2 = 1, w
        itmp(jj2) = jj2
      enddo
      do jj2 = 1, w
        if (itmp(jj2) == perm(jj2)) cycle
        do t2 = jj2 + 1, w
          if (itmp(t2) == perm(jj2)) exit
        enddo
        call band_swap_cols(pa+jj2-1, pa+t2-1, .false.)
        if (fct%lu) call band_swap_cols(pa+jj2-1, pa+t2-1, .true.)
        itmp(t2) = itmp(jj2)
        itmp(jj2) = perm(jj2)
      enddo
      if (fct%lu) then
        itmp(1:w) = perm(1:w)
        do jj2 = 1, w
          if (itmp(jj2) == permr(jj2)) cycle
          do t2 = jj2 + 1, w
            if (itmp(t2) == permr(jj2)) exit
          enddo
          call band_swap_cols(pa+jj2-1, pa+t2-1, .true.)
          itmp(t2) = itmp(jj2)
          itmp(jj2) = permr(jj2)
        enddo
      endif
      do jj2 = 1, np
        do t2 = 1, ncbt
          if (towner(t2) /= me) cycle
          h2 = ctb(t2) - ctb(t2-1)
          px = int(jj2-1, 8)*myrows + tprow(t2)
          bx = boff(t2) + int(pa + jj2 - 2, 8)*h2
          band(bx+1:bx+h2) = pband(px+1:px+h2)
          if (fct%lu) bandu(bx+1:bx+h2) = pbandu(px+1:px+h2)
        enddo
      enddo
      do i2 = 1, npair
        call band_swap_cols(swaps(2*i2-1), swaps(2*i2), .false.)
        if (fct%lu) call band_swap_cols(swaps(2*i2-1), swaps(2*i2), .true.)
      enddo
    end subroutine band_apply_piv

    !> after the pivots of tile column kk: the master hands the tile owners the pivot
    !> scaling data, its fully summed rows of the eliminated columns and the update
    !> vectors of the columns left delayed in the column; then every rank compresses,
    !> exchanges and updates what it owns, in the arithmetic of the sequential trailing
    !> update
    subroutine postcol(kk)
      integer(kind=kint), intent(in) :: kk
      real(kind=kreal), allocatable :: rcm(:)
      integer(kind=kint) :: q2, x2, mm, ib, nfsr, ndl, fac
      integer(kind=8) :: ox

      fac = 1
      if (fct%lu) fac = 2
      nfsr = ncol - g%tb(kk)
      ndl = g%tb(kk) - npiv
      if (ismaster) then
        ctli(1:8) = 0
        ctli(1) = 5
        ctli(2) = npiv
        ctli(3) = nfs
        call ctl_send_i(ctli, 8)
        if (pk > 0) then
          if (.not. fct%lu) then
            do q2 = p0 + 1, npiv
              ptv(q2-p0) = fct%sn(s)%ptype(q2)
              dsv(q2-p0) = fct%sn(s)%dsub(q2)
              dgv(q2-p0) = fvm(mf_idx(gm, q2, q2))
            enddo
          endif
          ! the fully summed rows of the eliminated columns, tile block by tile block, so
          !  a tile is a contiguous block on both sides
          call mf_grow_r(fsl, max(int(nfsr, 8)*pk, 1_8))
          if (fct%lu) call mf_grow_r(fslu, max(int(nfsr, 8)*pk, 1_8))
          do q2 = kk + 1, g%ntc
            x2 = g%tb(q2) - g%tb(q2-1)
            fsl(fsl_off(q2, kk)+1:fsl_off(q2, kk)+int(x2, 8)*pk) = &
              fvm(mf_off(gm, q2, kk)+1:mf_off(gm, q2, kk)+int(x2, 8)*pk)
            if (fct%lu) fslu(fsl_off(q2, kk)+1:fsl_off(q2, kk)+int(x2, 8)*pk) = &
              fvmu(mf_off(gm, q2, kk)+1:mf_off(gm, q2, kk)+int(x2, 8)*pk)
          enddo
          call mf_grow_r(tvec, max(int(ndl, 8)*pk, 1_8))
          if (fct%lu) call mf_grow_r(tvecu, max(int(ndl, 8)*pk, 1_8))
          do x2 = npiv + 1, g%tb(kk)
            if (.not. fct%lu) then
              call build_tcol(x2, tvec(int(x2-npiv-1, 8)*pk + 1))
            else
              call mf_row_get(gm, fvmu, x2, p0+1, npiv, tvec(int(x2-npiv-1, 8)*pk + 1))
              call mf_row_get(gm, fvm, x2, p0+1, npiv, tvecu(int(x2-npiv-1, 8)*pk + 1))
            endif
          enddo
          if (.not. fct%lu) call ctl_send_i(ptv, pk)
          mm = merge(0, 2*pk, fct%lu) + fac*nfsr*pk + fac*ndl*pk
          if (mm > 0 .and. nwk > merge(1, 0, wme > 0)) then
            ib = mf_pool_slot(pool)
            allocate(pool%box(ib)%rv(mm))
            ox = 0
            if (.not. fct%lu) then
              pool%box(ib)%rv(1:pk) = dsv(1:pk)
              pool%box(ib)%rv(pk+1:2*pk) = dgv(1:pk)
              ox = 2*pk
            endif
            pool%box(ib)%rv(ox+1:ox+int(nfsr, 8)*pk) = fsl(1:int(nfsr, 8)*pk)
            ox = ox + int(nfsr, 8)*pk
            if (fct%lu) then
              pool%box(ib)%rv(ox+1:ox+int(nfsr, 8)*pk) = fslu(1:int(nfsr, 8)*pk)
              ox = ox + int(nfsr, 8)*pk
            endif
            pool%box(ib)%rv(ox+1:ox+int(ndl, 8)*pk) = tvec(1:int(ndl, 8)*pk)
            ox = ox + int(ndl, 8)*pk
            if (fct%lu) pool%box(ib)%rv(ox+1:ox+int(ndl, 8)*pk) = tvecu(1:int(ndl, 8)*pk)
            call ctl_send_box(ib, mm)
          endif
        endif
      else
        call hecmw_recv_int(ctli, 8, master, 16*s+3, map%comm, stat)
        nfs = ctli(3)
        if (pk > 0) then
          if (.not. fct%lu) call hecmw_recv_int(ptv, pk, master, 16*s+3, map%comm, stat)
          mm = merge(0, 2*pk, fct%lu) + fac*nfsr*pk + fac*ndl*pk
          if (mm > 0) then
            allocate(rcm(mm))
            call hecmw_recv_r(rcm, mm, master, 16*s+4, map%comm, stat)
            ox = 0
            if (.not. fct%lu) then
              dsv(1:pk) = rcm(1:pk)
              dgv(1:pk) = rcm(pk+1:2*pk)
              ox = 2*pk
            endif
            call mf_grow_r(fsl, max(int(nfsr, 8)*pk, 1_8))
            fsl(1:int(nfsr, 8)*pk) = rcm(ox+1:ox+int(nfsr, 8)*pk)
            ox = ox + int(nfsr, 8)*pk
            if (fct%lu) then
              call mf_grow_r(fslu, max(int(nfsr, 8)*pk, 1_8))
              fslu(1:int(nfsr, 8)*pk) = rcm(ox+1:ox+int(nfsr, 8)*pk)
              ox = ox + int(nfsr, 8)*pk
            endif
            call mf_grow_r(tvec, max(int(ndl, 8)*pk, 1_8))
            tvec(1:int(ndl, 8)*pk) = rcm(ox+1:ox+int(ndl, 8)*pk)
            ox = ox + int(ndl, 8)*pk
            if (fct%lu) then
              call mf_grow_r(tvecu, max(int(ndl, 8)*pk, 1_8))
              tvecu(1:int(ndl, 8)*pk) = rcm(ox+1:ox+int(ndl, 8)*pk)
            endif
            deallocate(rcm)
          endif
        endif
      endif
      if (pk <= 0) return
      if (fct%blr .and. hasband) call compress_band(kk)
      if (hasband) call exchange_band(kk)
      call update_col(kk)
    end subroutine postcol

    !> word offset of the fully summed tile j2 in the tile blocked fsl buffer of column kk
    function fsl_off(j2, kk) result(ox)
      integer(kind=kint), intent(in) :: j2, kk
      integer(kind=8) :: ox

      ox = int(g%tb(j2-1) - g%tb(kk), 8)*pk
    end function fsl_off

    !> update vector of the delayed or panel column x2 for the pivots p0+1..npiv (the t of
    !> the sequential mf_col_update)
    subroutine build_tcol(x2, tv)
      integer(kind=kint), intent(in) :: x2
      real(kind=kreal), intent(out) :: tv(*)
      integer(kind=kint) :: q2
      real(kind=kreal) :: l1, l2, d11, d21, d22

      q2 = p0 + 1
      do while (q2 <= npiv)
        if (fct%sn(s)%ptype(q2) == 1) then
          tv(q2-p0) = fvm(mf_idx(gm, x2, q2)) * fvm(mf_idx(gm, q2, q2))
          q2 = q2 + 1
        else
          l1 = fvm(mf_idx(gm, x2, q2))
          l2 = fvm(mf_idx(gm, x2, q2+1))
          d11 = fvm(mf_idx(gm, q2, q2))
          d21 = fct%sn(s)%dsub(q2)
          d22 = fvm(mf_idx(gm, q2+1, q2+1))
          tv(q2-p0) = l1*d11 + l2*d21
          tv(q2-p0+1) = l1*d21 + l2*d22
          q2 = q2 + 2
        endif
      enddo
    end subroutine build_tcol

    !> compress the contribution row tiles I own in tile column kk, the slots reserved up
    !> front so the layout is schedule independent, as the sequential compression does
    subroutine compress_band(kk)
      integer(kind=kint), intent(in) :: kk
      integer(kind=kint) :: t2, h2, r2, idx2

      do t2 = 1, ncbt
        if (towner(t2) /= me) cycle
        h2 = ctb(t2) - ctb(t2-1)
        idx2 = mf_bidx(g%nt, kk, g%ntc + t2)
        bslot(idx2) = btop
        if (fct%lu) then
          btop = btop + int(h2, 8)*pk + int(pk, 8)*min(h2, pk)
          bslotu(idx2) = btopu
          btopu = btopu + int(h2, 8)*pk + int(pk, 8)*min(h2, pk)
        else
          btop = btop + int(h2, 8)*pk + 2_8*int(pk, 8)*min(h2, pk)
        endif
      enddo
      call mf_grow_r(bval, max(btop, 1_8))
      if (fct%lu) call mf_grow_r(bvalu, max(btopu, 1_8))
      !$omp parallel do default(shared) private(t2, h2, r2, idx2) schedule(dynamic, 1)
      do t2 = 1, ncbt
        if (towner(t2) /= me) cycle
        h2 = ctb(t2) - ctb(t2-1)
        idx2 = mf_bidx(g%nt, kk, g%ntc + t2)
        bval(bslot(idx2)+1:bslot(idx2)+int(h2, 8)*pk) = &
          band(boff(t2)+int(p0, 8)*h2+1:boff(t2)+int(p0, 8)*h2+int(h2, 8)*pk)
        call hecmw_mf_kernel_compress(h2, pk, h2, bval(bslot(idx2)+1), fct%eps, pk, &
          bval(bslot(idx2)+int(h2, 8)*pk+1), r2)
        brk(idx2) = r2
        if (fct%lu) then
          bvalu(bslotu(idx2)+1:bslotu(idx2)+int(h2, 8)*pk) = &
            bandu(boff(t2)+int(p0, 8)*h2+1:boff(t2)+int(p0, 8)*h2+int(h2, 8)*pk)
          call hecmw_mf_kernel_compress(h2, pk, h2, bvalu(bslotu(idx2)+1), fct%eps, pk, &
            bvalu(bslotu(idx2)+int(h2, 8)*pk+1), r2)
          brku(idx2) = r2
        endif
      enddo
      !$omp end parallel do
      do t2 = 1, ncbt
        if (towner(t2) /= me) cycle
        idx2 = mf_bidx(g%nt, kk, g%ntc + t2)
        ntl_f = ntl_f + merge(2, 1, fct%lu)
        if (brk(idx2) >= 0) then
          nlr_f = nlr_f + 1
          rsum_f = rsum_f + brk(idx2)
          rmax_f = max(rmax_f, brk(idx2))
        endif
        if (fct%lu) then
          if (brku(idx2) >= 0) then
            nlr_f = nlr_f + 1
            rsum_f = rsum_f + brku(idx2)
            rmax_f = max(rmax_f, brku(idx2))
          endif
        endif
      enddo
    end subroutine compress_band

    !> exchange the tile column kk data of the owned tiles among the tile owners: my
    !> tiles go to every later owner (an update (i,j) only reads tiles j <= i), the
    !> earlier owners' tiles are received and indexed for the update phase
    subroutine exchange_band(kk)
      integer(kind=kint), intent(in) :: kk
      integer(kind=kint), allocatable :: irk(:)
      integer(kind=kint) :: t2, h2, r2, ru2, iw2, ib, nti
      integer(kind=8) :: sz, ox

      if (.not. allocated(xw)) then
        allocate(xw(max(ncbt, 1)), xoff(max(ncbt, 1)), xrk(max(ncbt, 1)))
        allocate(xoffu(max(ncbt, 1)), xrku(max(ncbt, 1)), xbuf(max(nwk, 1)))
      endif
      ! pack my tiles once, then isend to every later owner
      if (wme < nwk) then
        sz = 0
        nti = 0
        do t2 = 1, ncbt
          if (towner(t2) /= me) cycle
          h2 = ctb(t2) - ctb(t2-1)
          nti = nti + 1
          sz = sz + xtile_words(t2, kk, .false.)
          if (fct%lu) sz = sz + xtile_words(t2, kk, .true.)
        enddo
        ib = mf_pool_slot(pool)
        allocate(pool%box(ib)%rv(max(sz, 1_8)))
        if (fct%blr) allocate(pool%box(ib)%hdr(max(merge(2, 1, fct%lu)*nti, 1)))
        ox = 0
        nti = 0
        do t2 = 1, ncbt
          if (towner(t2) /= me) cycle
          call xtile_pack(t2, kk, .false., pool%box(ib)%rv, ox)
          if (fct%blr) then
            nti = nti + 1
            pool%box(ib)%hdr(merge(2, 1, fct%lu)*(nti-1)+1) = xrank_of(t2, kk, .false.)
            if (fct%lu) pool%box(ib)%hdr(merge(2, 1, fct%lu)*(nti-1)+2) = xrank_of(t2, kk, .true.)
          endif
          if (fct%lu) call xtile_pack(t2, kk, .true., pool%box(ib)%rv, ox)
        enddo
        do iw2 = wme + 1, nwk
          if (fct%blr) then
            req = 0
            call hecmw_isend_int(pool%box(ib)%hdr, merge(2, 1, fct%lu)*max(nti, 1), wrank(iw2), 16*s+6, map%comm, req)
            call mf_pool_req(pool, req)
          endif
          req = 0
          call hecmw_isend_r(pool%box(ib)%rv, int(max(sz, 1_8), kind=kint), wrank(iw2), 16*s+7, map%comm, req)
          call mf_pool_req(pool, req)
        enddo
      endif
      ! receive from the earlier owners and index every tile below my range
      do iw2 = 1, wme - 1
        nti = 0
        do t2 = 1, ncbt
          if (towner(t2) == wrank(iw2)) nti = nti + 1
        enddo
        allocate(irk(merge(2, 1, fct%lu)*max(nti, 1)))
        if (fct%blr) then
          call hecmw_recv_int(irk, merge(2, 1, fct%lu)*max(nti, 1), wrank(iw2), 16*s+6, map%comm, stat)
        else
          irk(:) = -1
        endif
        sz = 0
        nti = 0
        do t2 = 1, ncbt
          if (towner(t2) /= wrank(iw2)) cycle
          h2 = ctb(t2) - ctb(t2-1)
          nti = nti + 1
          r2 = irk(merge(2, 1, fct%lu)*(nti-1)+1)
          xw(t2) = iw2
          xrk(t2) = r2
          xoff(t2) = sz
          if (r2 < 0) then
            sz = sz + int(h2, 8)*pk
          else
            sz = sz + int(h2, 8)*r2 + int(pk, 8)*r2
          endif
          if (fct%lu) then
            ru2 = irk(merge(2, 1, fct%lu)*(nti-1)+2)
            xrku(t2) = ru2
            xoffu(t2) = sz
            if (ru2 < 0) then
              sz = sz + int(h2, 8)*pk
            else
              sz = sz + int(h2, 8)*ru2 + int(pk, 8)*ru2
            endif
          endif
        enddo
        deallocate(irk)
        if (allocated(xbuf(iw2)%v)) then
          if (size(xbuf(iw2)%v, kind=8) < sz) deallocate(xbuf(iw2)%v)
        endif
        if (.not. allocated(xbuf(iw2)%v)) allocate(xbuf(iw2)%v(max(sz, 1_8)))
        call hecmw_recv_r(xbuf(iw2)%v, int(max(sz, 1_8), kind=kint), wrank(iw2), 16*s+7, map%comm, stat)
      enddo
    end subroutine exchange_band

    !> stored rank of my tile t2 in column kk (-1 = full rank)
    function xrank_of(t2, kk, uface) result(r2)
      integer(kind=kint), intent(in) :: t2, kk
      logical, intent(in) :: uface
      integer(kind=kint) :: r2

      r2 = -1
      if (fct%blr) then
        if (uface) then
          r2 = brku(mf_bidx(g%nt, kk, g%ntc + t2))
        else
          r2 = brk(mf_bidx(g%nt, kk, g%ntc + t2))
        endif
      endif
    end function xrank_of

    !> words of the exchanged data of my tile t2 in column kk
    function xtile_words(t2, kk, uface) result(ww)
      integer(kind=kint), intent(in) :: t2, kk
      logical, intent(in) :: uface
      integer(kind=8) :: ww
      integer(kind=kint) :: h2, r2

      h2 = ctb(t2) - ctb(t2-1)
      r2 = xrank_of(t2, kk, uface)
      if (r2 < 0) then
        ww = int(h2, 8)*pk
      else
        ww = int(h2, 8)*r2 + int(pk, 8)*r2
      endif
    end function xtile_words

    !> the exchanged data of my tile t2 (the eliminated columns, or U and V) into buf
    subroutine xtile_pack(t2, kk, uface, buf, ox)
      integer(kind=kint), intent(in) :: t2, kk
      logical, intent(in) :: uface
      real(kind=kreal), intent(inout) :: buf(:)
      integer(kind=8), intent(inout) :: ox
      integer(kind=kint) :: h2, r2, idx2
      integer(kind=8) :: bx

      h2 = ctb(t2) - ctb(t2-1)
      r2 = xrank_of(t2, kk, uface)
      if (r2 < 0) then
        bx = boff(t2) + int(p0, 8)*h2
        if (uface) then
          buf(ox+1:ox+int(h2, 8)*pk) = bandu(bx+1:bx+int(h2, 8)*pk)
        else
          buf(ox+1:ox+int(h2, 8)*pk) = band(bx+1:bx+int(h2, 8)*pk)
        endif
        ox = ox + int(h2, 8)*pk
      else
        idx2 = mf_bidx(g%nt, kk, g%ntc + t2)
        if (uface) then
          bx = bslotu(idx2)
          buf(ox+1:ox+int(h2, 8)*r2) = bvalu(bx+1:bx+int(h2, 8)*r2)
          buf(ox+int(h2, 8)*r2+1:ox+int(h2, 8)*r2+int(pk, 8)*r2) = &
            bvalu(bx+int(h2, 8)*pk+1:bx+int(h2, 8)*pk+int(pk, 8)*r2)
        else
          bx = bslot(idx2)
          buf(ox+1:ox+int(h2, 8)*r2) = bval(bx+1:bx+int(h2, 8)*r2)
          buf(ox+int(h2, 8)*r2+1:ox+int(h2, 8)*r2+int(pk, 8)*r2) = &
            bval(bx+int(h2, 8)*pk+1:bx+int(h2, 8)*pk+int(pk, 8)*r2)
        endif
        ox = ox + int(h2, 8)*r2 + int(pk, 8)*r2
      endif
    end subroutine xtile_pack

    !> the trailing updates of tile column kk on what I hold: the master updates its fully
    !> summed tile pairs, every tile owner updates its band tiles and the delayed columns,
    !> with the operand kinds and per tile arithmetic of the sequential trailing update
    subroutine update_col(kk)
      integer(kind=kint), intent(in) :: kk
      integer(kind=kint) :: t2, h2, i2, j2, l2, q2, x2, mn2, r2
      integer(kind=8) :: okk2, ox

      ! scaled panels of the tiles I own (LDLt)
      if (.not. fct%lu) then
        if (ismaster .and. g%ntc > kk) then
          call mf_grow_r(wscm, max(int(ncol - g%tb(kk), 8)*pk, 1_8))
          okk2 = mf_off(gm, kk, kk)
          !$omp parallel do default(shared) private(i2, h2) schedule(dynamic, 1)
          do i2 = kk + 1, g%ntc
            h2 = g%tb(i2) - g%tb(i2-1)
            call hecmw_mf_kernel_scale(h2, pk, g%tb(kk) - g%tb(kk-1), fvm(okk2+1), &
              fct%sn(s)%ptype(p0+1:npiv), fct%sn(s)%dsub(p0+1:npiv), h2, fvm(mf_off(gm, i2, kk)+1), &
              wscm(int(g%tb(i2-1) - g%tb(kk), 8)*pk + 1))
          enddo
          !$omp end parallel do
        endif
        if (hasband) then
          call mf_grow_r(dd, int(pk, 8)*pk)
          dd(1:int(pk, 8)*pk) = 0.0d0
          do q2 = 1, pk
            dd(int(q2-1, 8)*pk + q2) = dgv(q2)
          enddo
          call mf_grow_r(wscb, max(int(myrows, 8)*pk, 1_8))
          !$omp parallel do default(shared) private(t2, h2, r2, mn2, ox) schedule(dynamic, 1)
          do t2 = 1, ncbt
            if (towner(t2) /= me) cycle
            h2 = ctb(t2) - ctb(t2-1)
            r2 = xrank_of(t2, kk, .false.)
            if (r2 < 0) then
              call hecmw_mf_kernel_scale(h2, pk, pk, dd, ptv, dsv, h2, &
                band(boff(t2) + int(p0, 8)*h2 + 1), wscb(int(tprow(t2), 8)*pk + 1))
            else if (r2 > 0) then
              mn2 = min(h2, pk)
              ox = bslot(mf_bidx(g%nt, kk, g%ntc + t2))
              call hecmw_mf_kernel_scale_rows(pk, r2, pk, dd, ptv, dsv, pk, &
                bval(ox + int(h2, 8)*pk + 1), pk, bval(ox + int(h2, 8)*pk + int(pk, 8)*mn2 + 1))
            endif
          enddo
          !$omp end parallel do
        endif
      endif
      ! the tile pairs I update
      npair = 0
      if (ismaster) then
        do i2 = kk + 1, g%ntc
          do j2 = kk + 1, i2
            npair = npair + 1
            pair_i(npair) = i2
            pair_j(npair) = j2
          enddo
        enddo
      endif
      do t2 = 1, ncbt
        if (towner(t2) /= me) cycle
        do j2 = kk + 1, g%ntc + t2
          npair = npair + 1
          pair_i(npair) = g%ntc + t2
          pair_j(npair) = j2
        enddo
      enddo
      !$omp parallel do default(shared) private(l2) schedule(dynamic, 1)
      do l2 = 1, npair
        call update_pair(kk, pair_i(l2), pair_j(l2))
      enddo
      !$omp end parallel do
      ! the columns left delayed in the tile column receive the pivots as vector updates
      if (ismaster) then
        do x2 = npiv + 1, g%tb(kk)
          if (.not. fct%lu) then
            call mf_col_copy(gm, fvm, x2, x2, pv, .false.)
            call mf_col_update(gm, fvm, fct%sn(s)%ptype, fct%sn(s)%dsub, x2, p0+1, npiv, pv, wk)
            call mf_col_copy(gm, fvm, x2, x2, pv, .true.)
          else
            call mf_col_copy(gm, fvm, x2, x2, pv, .false.)
            call mf_row_get(gm, fvmu, x2, p0+1, npiv, wk)
            call mf_col_axpy(gm, fvm, x2, p0+1, npiv, wk, pv)
            call mf_col_copy(gm, fvm, x2, x2, pv, .true.)
            if (x2 < ncol) then
              call mf_col_copy(gm, fvmu, x2, x2+1, pv, .false.)
              call mf_row_get(gm, fvm, x2, p0+1, npiv, wk)
              call mf_col_axpy(gm, fvmu, x2+1, p0+1, npiv, wk, pv)
              call mf_col_copy(gm, fvmu, x2, x2+1, pv, .true.)
            endif
          endif
        enddo
      endif
      if (hasband) then
        !$omp parallel do default(shared) private(t2, h2, x2) schedule(dynamic, 1)
        do t2 = 1, ncbt
          if (towner(t2) /= me) cycle
          h2 = ctb(t2) - ctb(t2-1)
          do x2 = npiv + 1, g%tb(kk)
            call hecmw_mf_kernel_gemv(h2, pk, h2, band(boff(t2) + int(p0, 8)*h2 + 1), &
              tvec(int(x2-npiv-1, 8)*pk + 1), band(boff(t2) + int(x2-1, 8)*h2 + 1))
            if (fct%lu) call hecmw_mf_kernel_gemv(h2, pk, h2, bandu(boff(t2) + int(p0, 8)*h2 + 1), &
              tvecu(int(x2-npiv-1, 8)*pk + 1), bandu(boff(t2) + int(x2-1, 8)*h2 + 1))
          enddo
        enddo
        !$omp end parallel do
      endif
    end subroutine update_col

    !> one trailing tile update (i2, j2) of tile column kk, in the operand combination of
    !> the sequential mf_update_ab
    subroutine update_pair(kk, i2, j2)
      integer(kind=kint), intent(in) :: kk, i2, j2
      real(kind=kreal), pointer, contiguous :: fa(:), ua(:), va(:), fb(:), ub(:), vb(:)
      integer(kind=kint) :: hi2, wj2, ra2, rb2, t2, rua2, rlb2
      integer(kind=8) :: ox, oc

      hi2 = g%tb(i2) - g%tb(i2-1)
      wj2 = g%tb(j2) - g%tb(j2-1)
      if (i2 <= g%ntc) then
        ! a fully summed pair of the master: both sides full rank
        if (.not. fct%lu) then
          call hecmw_mf_kernel_gemm(hi2, wj2, pk, hi2, wscm(int(g%tb(i2-1) - g%tb(kk), 8)*pk + 1), &
            wj2, fvm(mf_off(gm, j2, kk)+1), hi2, fvm(mf_off(gm, i2, j2)+1))
        else
          call hecmw_mf_kernel_gemm(hi2, wj2, pk, hi2, fvm(mf_off(gm, i2, kk)+1), &
            wj2, fvmu(mf_off(gm, j2, kk)+1), hi2, fvm(mf_off(gm, i2, j2)+1))
          call hecmw_mf_kernel_gemm(hi2, wj2, pk, hi2, fvmu(mf_off(gm, i2, kk)+1), &
            wj2, fvm(mf_off(gm, j2, kk)+1), hi2, fvmu(mf_off(gm, i2, j2)+1))
        endif
        return
      endif
      t2 = i2 - g%ntc
      oc = boff(t2) + int(hi2, 8)*g%tb(j2-1)
      if (.not. fct%lu) then
        ra2 = xrank_of(t2, kk, .false.)
        if (ra2 < 0) then
          fa => wscb(int(tprow(t2), 8)*pk + 1 : int(tprow(t2), 8)*pk + int(hi2, 8)*pk)
          ua => dum(1:1)
          va => dum(1:1)
        else
          ox = bslot(mf_bidx(g%nt, kk, i2))
          fa => dum(1:1)
          ua => bval(ox+1 : ox+int(hi2, 8)*max(ra2, 1))
          va => bval(ox + int(hi2, 8)*pk + int(pk, 8)*min(hi2, pk) + 1 : &
                     ox + int(hi2, 8)*pk + int(pk, 8)*min(hi2, pk) + int(pk, 8)*max(ra2, 1))
        endif
        call bside(kk, j2, .false., rb2, fb, ub, vb)
        call mf_update_ab(hi2, wj2, pk, ra2, fa, ua, va, rb2, fb, ub, vb, band(oc+1:oc+int(hi2, 8)*wj2))
      else
        ra2 = xrank_of(t2, kk, .false.)
        rua2 = xrank_of(t2, kk, .true.)
        call aside_lu(kk, t2, hi2, .false., ra2, fa, ua, va)
        call bside(kk, j2, .true., rb2, fb, ub, vb)
        call mf_update_ab(hi2, wj2, pk, ra2, fa, ua, va, rb2, fb, ub, vb, band(oc+1:oc+int(hi2, 8)*wj2))
        call aside_lu(kk, t2, hi2, .true., rua2, fa, ua, va)
        call bside(kk, j2, .false., rlb2, fb, ub, vb)
        call mf_update_ab(hi2, wj2, pk, rua2, fa, ua, va, rlb2, fb, ub, vb, bandu(oc+1:oc+int(hi2, 8)*wj2))
      endif
    end subroutine update_pair

    !> A side operands of my tile t2 in the LU mode (raw block, or U and V)
    subroutine aside_lu(kk, t2, hi2, uface, ra2, fa, ua, va)
      integer(kind=kint), intent(in) :: kk, t2, hi2
      logical, intent(in) :: uface
      integer(kind=kint), intent(in) :: ra2
      real(kind=kreal), pointer, contiguous, intent(out) :: fa(:), ua(:), va(:)
      integer(kind=8) :: ox

      if (ra2 < 0) then
        if (uface) then
          fa => bandu(boff(t2) + int(p0, 8)*hi2 + 1 : boff(t2) + int(p0, 8)*hi2 + int(hi2, 8)*pk)
        else
          fa => band(boff(t2) + int(p0, 8)*hi2 + 1 : boff(t2) + int(p0, 8)*hi2 + int(hi2, 8)*pk)
        endif
        ua => dum(1:1)
        va => dum(1:1)
      else
        fa => dum(1:1)
        if (uface) then
          ox = bslotu(mf_bidx(g%nt, kk, g%ntc + t2))
          ua => bvalu(ox+1 : ox+int(hi2, 8)*max(ra2, 1))
          va => bvalu(ox + int(hi2, 8)*pk + 1 : ox + int(hi2, 8)*pk + int(pk, 8)*max(ra2, 1))
        else
          ox = bslot(mf_bidx(g%nt, kk, g%ntc + t2))
          ua => bval(ox+1 : ox+int(hi2, 8)*max(ra2, 1))
          va => bval(ox + int(hi2, 8)*pk + 1 : ox + int(hi2, 8)*pk + int(pk, 8)*max(ra2, 1))
        endif
      endif
    end subroutine aside_lu

    !> B side operands of tile j2 of column kk: the master's fully summed rows, my own
    !> tile, or the exchanged data of another owner; uface selects the upper grid data
    subroutine bside(kk, j2, uface, rb2, fb, ub, vb)
      integer(kind=kint), intent(in) :: kk, j2
      logical, intent(in) :: uface
      integer(kind=kint), intent(out) :: rb2
      real(kind=kreal), pointer, contiguous, intent(out) :: fb(:), ub(:), vb(:)
      integer(kind=kint) :: tj2, hj2, iwx
      integer(kind=8) :: ox

      hj2 = g%tb(j2) - g%tb(j2-1)
      rb2 = -1
      fb => dum(1:1)
      ub => dum(1:1)
      vb => dum(1:1)
      if (j2 <= g%ntc) then
        if (ismaster) then
          if (uface) then
            fb => fvmu(mf_off(gm, j2, kk)+1 : mf_off(gm, j2, kk)+int(hj2, 8)*pk)
          else
            fb => fvm(mf_off(gm, j2, kk)+1 : mf_off(gm, j2, kk)+int(hj2, 8)*pk)
          endif
        else
          if (uface) then
            fb => fslu(fsl_off(j2, kk)+1 : fsl_off(j2, kk)+int(hj2, 8)*pk)
          else
            fb => fsl(fsl_off(j2, kk)+1 : fsl_off(j2, kk)+int(hj2, 8)*pk)
          endif
        endif
        return
      endif
      tj2 = j2 - g%ntc
      if (towner(tj2) == me) then
        rb2 = xrank_of(tj2, kk, uface)
        if (rb2 < 0) then
          if (uface) then
            fb => bandu(boff(tj2) + int(p0, 8)*hj2 + 1 : boff(tj2) + int(p0, 8)*hj2 + int(hj2, 8)*pk)
          else
            fb => band(boff(tj2) + int(p0, 8)*hj2 + 1 : boff(tj2) + int(p0, 8)*hj2 + int(hj2, 8)*pk)
          endif
        else
          if (uface) then
            ox = bslotu(mf_bidx(g%nt, kk, j2))
            ub => bvalu(ox+1 : ox+int(hj2, 8)*max(rb2, 1))
            vb => bvalu(ox + int(hj2, 8)*pk + 1 : ox + int(hj2, 8)*pk + int(pk, 8)*max(rb2, 1))
          else
            ox = bslot(mf_bidx(g%nt, kk, j2))
            ub => bval(ox+1 : ox+int(hj2, 8)*max(rb2, 1))
            if (fct%lu) then
              vb => bval(ox + int(hj2, 8)*pk + 1 : ox + int(hj2, 8)*pk + int(pk, 8)*max(rb2, 1))
            else
              vb => bval(ox + int(hj2, 8)*pk + 1 : ox + int(hj2, 8)*pk + int(pk, 8)*max(rb2, 1))
            endif
          endif
        endif
      else
        iwx = xw(tj2)
        if (uface) then
          rb2 = xrku(tj2)
          ox = xoffu(tj2)
        else
          rb2 = xrk(tj2)
          ox = xoff(tj2)
        endif
        if (rb2 < 0) then
          fb => xbuf(iwx)%v(ox+1 : ox+int(hj2, 8)*pk)
        else
          ub => xbuf(iwx)%v(ox+1 : ox+int(hj2, 8)*max(rb2, 1))
          vb => xbuf(iwx)%v(ox + int(hj2, 8)*max(rb2, 0) + 1 : ox + int(hj2, 8)*max(rb2, 0) + int(pk, 8)*max(rb2, 1))
        endif
      endif
    end subroutine bside

    !> inertia of the eliminated pivots from the master grid, as the sequential count
    subroutine inertia_1d()
      integer(kind=kint) :: x2
      real(kind=kreal) :: d11, d21, d22, det

      x2 = 1
      do while (x2 <= npiv)
        if (fct%sn(s)%ptype(x2) == 1) then
          if (fvm(mf_idx(gm, x2, x2)) > 0.0d0) then
            npos_f = npos_f + 1
          else
            nneg_f = nneg_f + 1
          endif
          x2 = x2 + 1
        else
          d11 = fvm(mf_idx(gm, x2, x2))
          d21 = fct%sn(s)%dsub(x2)
          d22 = fvm(mf_idx(gm, x2+1, x2+1))
          det = d11*d22 - d21*d21
          if (det < 0.0d0) then
            npos_f = npos_f + 1
            nneg_f = nneg_f + 1
          else if (d11 + d22 > 0.0d0) then
            npos_f = npos_f + 2
          else
            nneg_f = nneg_f + 2
          endif
          x2 = x2 + 2
        endif
      enddo
    end subroutine inertia_1d

    !> store what I hold: the master its fully summed panel tiles, every tile owner its
    !> band tiles (full rank or compressed), all addressed through bptr; then the owned
    !> blocks of the contribution block in the canonical block order the extend-add of
    !> the parent traverses
    subroutine store_1d()
      integer(kind=kint) :: kk, t2, h2, mw, idx2, r2, i2, j2, cc2, rr2, colx, cnbP, kbegP, cndelP
      integer(kind=8) :: off, offu, pw0, bx, sb, halfP
      integer(kind=kint), allocatable :: ctbP(:)

      nkc = 0
      do kk = 1, g%ntc
        if (min(g%tb(kk), npiv) - g%tb(kk-1) <= 0) exit
        nkc = kk
      enddo
      if (nkc > 0) then
        idx2 = mf_bidx(g%nt, nkc, g%nt)
        call mf_grow_i(fct%sn(s)%brank, idx2)
        call mf_grow_i8(fct%sn(s)%bptr, idx2 + 1)
        if (fct%lu) then
          call mf_grow_i(fct%sn(s)%branku, idx2)
          call mf_grow_i8(fct%sn(s)%bptru, idx2 + 1)
        endif
      endif
      off = 0
      offu = 0
      pw0 = 0
      do kk = 1, nkc
        mw = min(g%tb(kk), npiv) - g%tb(kk-1)
        do i2 = kk, g%nt
          if (.not. tile_is_mine(i2)) cycle
          h2 = g%tb(i2) - g%tb(i2-1)
          idx2 = mf_bidx(g%nt, kk, i2)
          pw0 = pw0 + int(h2, 8)*mw
          r2 = -1
          if (fct%blr .and. i2 > g%ntc) r2 = brk(idx2)
          fct%sn(s)%brank(idx2) = r2
          fct%sn(s)%bptr(idx2) = off
          if (r2 < 0) then
            off = off + int(h2, 8)*mw
          else
            off = off + int(r2, 8)*(h2 + mw)
          endif
          if (fct%lu) then
            r2 = -1
            if (fct%blr .and. i2 > g%ntc) r2 = brku(idx2)
            fct%sn(s)%branku(idx2) = r2
            fct%sn(s)%bptru(idx2) = offu
            if (r2 < 0) then
              offu = offu + int(h2, 8)*mw
            else
              offu = offu + int(r2, 8)*(h2 + mw)
            endif
          endif
        enddo
      enddo
      call mf_grow_r(fct%sn(s)%lval, max(off, 1_8))
      if (fct%lu) call mf_grow_r(fct%sn(s)%uval, max(offu, 1_8))
      do kk = 1, nkc
        mw = min(g%tb(kk), npiv) - g%tb(kk-1)
        do i2 = kk, g%nt
          if (.not. tile_is_mine(i2)) cycle
          call store_tile(kk, i2, mw)
        enddo
      enddo
      fct%sn(s)%pwords = off + offu
      ! the sequential accounting counts the U panel only with BLR, where its stored
      ! words differ from the L panel
      fct%factor_words_act = fct%factor_words_act + off
      if (fct%lu .and. fct%blr) fct%factor_words_act = fct%factor_words_act + offu
      if (fct%blr) then
        fct%blr_words_fr = fct%blr_words_fr + merge(2, 1, fct%lu)*pw0
      endif

      ! the contribution block, canonical block order (see the extend-add traversal)
      cndelP = ncol - npiv
      allocate(ctbP(0:g%nt - g%ntc + 2))
      cnbP = 0
      ctbP(0) = 0
      if (cndelP > 0) then
        cnbP = 1
        ctbP(1) = cndelP
      endif
      kbegP = cnbP
      do j2 = g%ntc + 1, g%nt
        cnbP = cnbP + 1
        ctbP(cnbP) = cndelP + g%tb(j2) - ncol
      enddo
      halfP = 0
      do j2 = 1, cnbP
        do i2 = j2, cnbP
          if (.not. cbblk_is_mine(i2, kbegP)) cycle
          halfP = halfP + int(ctbP(i2) - ctbP(i2-1), 8)*(ctbP(j2) - ctbP(j2-1))
        enddo
      enddo
      fct%sn(s)%cbsize = merge(2, 1, fct%lu)*halfP
      if (allocated(fct%sn(s)%cval)) deallocate(fct%sn(s)%cval)
      if (fct%sn(s)%cbsize > 0) then
        allocate(fct%sn(s)%cval(fct%sn(s)%cbsize))
        sb = 0
        do j2 = 1, cnbP
          do i2 = j2, cnbP
            if (.not. cbblk_is_mine(i2, kbegP)) cycle
            h2 = ctbP(i2) - ctbP(i2-1)
            do cc2 = ctbP(j2-1) + 1, ctbP(j2)
              colx = cb_col(cc2, cndelP)
              if (i2 <= kbegP) then
                ! the delayed square of the master, rows on or below the diagonal
                do rr2 = max(cc2, ctbP(i2-1) + 1), ctbP(i2)
                  fct%sn(s)%cval(sb + int(cc2 - ctbP(j2-1) - 1, 8)*h2 + (rr2 - ctbP(i2-1))) = &
                    fvm(mf_idx(gm, npiv + rr2, npiv + cc2))
                  if (fct%lu) fct%sn(s)%cval(halfP + sb + int(cc2 - ctbP(j2-1) - 1, 8)*h2 + (rr2 - ctbP(i2-1))) = &
                    fvmu(mf_idx(gm, npiv + rr2, npiv + cc2))
                enddo
              else
                t2 = i2 - kbegP
                bx = boff(t2) + int(colx - 1, 8)*h2
                fct%sn(s)%cval(sb + int(cc2 - ctbP(j2-1) - 1, 8)*h2 + 1 : &
                               sb + int(cc2 - ctbP(j2-1) - 1, 8)*h2 + h2) = band(bx+1:bx+h2)
                if (fct%lu) fct%sn(s)%cval(halfP + sb + int(cc2 - ctbP(j2-1) - 1, 8)*h2 + 1 : &
                                           halfP + sb + int(cc2 - ctbP(j2-1) - 1, 8)*h2 + h2) = bandu(bx+1:bx+h2)
              endif
            enddo
            sb = sb + int(ctbP(j2) - ctbP(j2-1), 8)*h2
          enddo
        enddo
        fct%live_cb = fct%live_cb + fct%sn(s)%cbsize
        fct%stack_peak_act = max(fct%stack_peak_act, fct%live_cb)
      endif
      deallocate(ctbP)
    end subroutine store_1d

    !> whether front tile i2 is stored on this rank
    function tile_is_mine(i2) result(mine0)
      integer(kind=kint), intent(in) :: i2
      logical :: mine0

      if (i2 <= g%ntc) then
        mine0 = ismaster
      else
        mine0 = (towner(i2 - g%ntc) == me)
      endif
    end function tile_is_mine

    !> whether contribution block row i2 (canonical numbering) is stored on this rank
    function cbblk_is_mine(i2, kbegP) result(mine0)
      integer(kind=kint), intent(in) :: i2, kbegP
      logical :: mine0

      if (i2 <= kbegP) then
        mine0 = ismaster
      else
        mine0 = (towner(i2 - kbegP) == me)
      endif
    end function cbblk_is_mine

    !> front column of contribution block column cc2 (delayed columns first)
    function cb_col(cc2, cndelP) result(colx)
      integer(kind=kint), intent(in) :: cc2, cndelP
      integer(kind=kint) :: colx

      if (cc2 <= cndelP) then
        colx = npiv + cc2
      else
        colx = ncol + (cc2 - cndelP)
      endif
    end function cb_col

    !> one stored panel tile (i2, kk): the master from its grid, a tile owner from its
    !> band or its compression slot
    subroutine store_tile(kk, i2, mw)
      integer(kind=kint), intent(in) :: kk, i2, mw
      integer(kind=kint) :: h2, r2, t2, idx2
      integer(kind=8) :: bx, dx

      h2 = g%tb(i2) - g%tb(i2-1)
      idx2 = mf_bidx(g%nt, kk, i2)
      dx = fct%sn(s)%bptr(idx2)
      if (i2 <= g%ntc) then
        bx = mf_off(gm, i2, kk)
        fct%sn(s)%lval(dx+1:dx+int(h2, 8)*mw) = fvm(bx+1:bx+int(h2, 8)*mw)
        if (fct%lu) then
          dx = fct%sn(s)%bptru(idx2)
          fct%sn(s)%uval(dx+1:dx+int(h2, 8)*mw) = fvmu(bx+1:bx+int(h2, 8)*mw)
        endif
        return
      endif
      t2 = i2 - g%ntc
      r2 = fct%sn(s)%brank(idx2)
      if (r2 < 0) then
        bx = boff(t2) + int(g%tb(kk-1), 8)*h2
        fct%sn(s)%lval(dx+1:dx+int(h2, 8)*mw) = band(bx+1:bx+int(h2, 8)*mw)
      else
        bx = bslot(idx2)
        fct%sn(s)%lval(dx+1:dx+int(h2, 8)*r2) = bval(bx+1:bx+int(h2, 8)*r2)
        fct%sn(s)%lval(dx+int(h2, 8)*r2+1:dx+int(h2, 8)*r2+int(mw, 8)*r2) = &
          bval(bx+int(h2, 8)*mw+1:bx+int(h2, 8)*mw+int(mw, 8)*r2)
      endif
      if (fct%lu) then
        dx = fct%sn(s)%bptru(idx2)
        r2 = fct%sn(s)%branku(idx2)
        if (r2 < 0) then
          bx = boff(t2) + int(g%tb(kk-1), 8)*h2
          fct%sn(s)%uval(dx+1:dx+int(h2, 8)*mw) = bandu(bx+1:bx+int(h2, 8)*mw)
        else
          bx = bslotu(idx2)
          fct%sn(s)%uval(dx+1:dx+int(h2, 8)*r2) = bvalu(bx+1:bx+int(h2, 8)*r2)
          fct%sn(s)%uval(dx+int(h2, 8)*r2+1:dx+int(h2, 8)*r2+int(mw, 8)*r2) = &
            bvalu(bx+int(h2, 8)*mw+1:bx+int(h2, 8)*mw+int(mw, 8)*r2)
        endif
      endif
    end subroutine store_tile

  end subroutine mf_super_factor_1d

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
      sn%pwords = pw
      if (fct%lu) sn%pwords = 2*pw
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
      sn%pwords = pwa
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
  !> vectors in the original DOF numbering. The solves of the upper fronts follow the 1D
  !> row distribution of the factorization: the master runs the fully summed sweeps, every
  !> tile owner applies its stored rows, and the forward contributions travel as row
  !> slices to the ranks holding the parent rows. In the backward substitution the running
  !> right hand side of a tile column travels along the tile owners in ascending order (a
  !> token chain), so the reduction keeps the sequential fold order and the solution stays
  !> bitwise identical to the sequential solve. The solved pivot values of the owned
  !> supernodes are summed over the ranks at the end, which is exact and order independent
  !> because every DOF is eliminated on exactly one rank.
  subroutine hecmw_mf_numeric_solve_mpi(sym, map, fct, b, x)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_factor), intent(inout) :: fct
    real(kind=kreal), intent(in) :: b(:)
    real(kind=kreal), intent(out) :: x(:)
    type(mf_pool) :: pool, gpool
    real(kind=kreal), allocatable :: z(:), xx(:)
    real(kind=kreal), allocatable :: v(:), vb(:), tok(:), pivs(:), tbuf(:)
    integer(kind=kint), allocatable :: dmap(:), left(:)
    integer(kind=kint), allocatable :: ctb(:), towner(:), crow2t(:), tprow(:), wrank(:), rdof(:), rowoff(:)
    integer(kind=kint), allocatable :: cctb2(:), cown2(:), pmapr(:), pdofm(:), crown(:), mdel(:)
    integer(kind=kint) :: ns, s, i, iu, j, c, p, stat(HECMW_STATUS_SIZE), req
    integer(kind=kint) :: ncol, npv, ntc0, nt0, nrow0, ncb0, ncbt, nwk, wme, myrows, master, ndel0
    integer(kind=kint) :: cnrow0, cndel0, cncbt
    logical :: ismaster, hasband

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

    do iu = 1, map%nupper
      s = map%uplist(iu)
      if (map%myrank >= map%rbeg(s) .and. map%myrank < map%rbeg(s) + map%rcnt(s)) call fwd_1d(s)
    enddo

    do iu = map%nupper, 1, -1
      s = map%uplist(iu)
      if (map%myrank >= map%rbeg(s) .and. map%myrank < map%rbeg(s) + map%rcnt(s)) call bwd_1d(s)
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
    call mf_pool_wait(gpool)

    z(1:fct%ndof_tot) = 0.0d0
    do s = 1, ns
      if (map%owner(s) /= map%myrank) cycle
      do i = 1, fct%sn(s)%npiv
        z(fct%pdof(fct%sn(s)%fsdof(i))) = xx(fct%sn(s)%fsdof(i))
      enddo
    enddo
    x(1:fct%ndof_tot) = z(1:fct%ndof_tot)
    call hecmw_allreduce_R_comm(x, fct%ndof_tot, hecmw_sum, map%comm)
    deallocate(z, xx, dmap, left)

  contains

    !> row geometry of the upper front s from its stored metadata and the mapping
    subroutine geo_1d(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: t2, i2, k2, d2, nown2, a2

      ncol = fct%sn(s0)%ncol
      npv = fct%sn(s0)%npiv
      ntc0 = fct%sn(s0)%ntc
      nt0 = fct%sn(s0)%nt
      nrow0 = fct%sn(s0)%tbnd(nt0+1)
      ncb0 = nrow0 - ncol
      master = map%owner(s0)
      ismaster = (master == map%myrank)
      ndel0 = ncol - (fct%cdofptr(sym%sptr(s0+1)) - fct%cdofptr(sym%sptr(s0)))
      call hecmw_mf_dist_rowtiles(sym, map, fct%tile, s0, ndel0, ncbt, ctb, towner)
      if (allocated(crow2t)) then
        if (size(crow2t) < ncb0) deallocate(crow2t)
      endif
      if (.not. allocated(crow2t)) allocate(crow2t(max(ncb0, 1)))
      do t2 = 1, ncbt
        crow2t(ctb(t2-1)+1:ctb(t2)) = t2
      enddo
      if (allocated(wrank)) then
        if (size(wrank) < ncbt + 1) deallocate(wrank, tprow)
      endif
      if (.not. allocated(wrank)) allocate(wrank(max(ncbt, 1)), tprow(0:max(ncbt, 1)))
      nwk = 0
      wme = 0
      myrows = 0
      tprow(0) = 0
      do t2 = 1, ncbt
        if (nwk == 0) then
          nwk = 1
          wrank(1) = towner(t2)
        else if (towner(t2) /= wrank(nwk)) then
          nwk = nwk + 1
          wrank(nwk) = towner(t2)
        endif
        if (wrank(nwk) == map%myrank) wme = nwk
        tprow(t2) = myrows
        if (towner(t2) == map%myrank) myrows = myrows + ctb(t2) - ctb(t2-1)
      enddo
      hasband = (myrows > 0)
      nown2 = sym%sptr(s0+1) - sym%sptr(s0)
      if (allocated(rowoff)) then
        if (size(rowoff) < sym%rptr(s0+1) - sym%rptr(s0) + 1) deallocate(rowoff)
      endif
      if (.not. allocated(rowoff)) allocate(rowoff(sym%rptr(s0+1) - sym%rptr(s0) + 1))
      rowoff(1) = 0
      do i2 = 1, sym%rptr(s0+1) - sym%rptr(s0)
        k2 = sym%rlist(sym%rptr(s0)+i2-1)
        rowoff(i2+1) = rowoff(i2) + sym%ndof(k2)
        if (i2 == nown2) rowoff(i2+1) = rowoff(i2+1) + ndel0
      enddo
      ! the column DOFs of the contribution rows
      if (allocated(rdof)) then
        if (size(rdof) < ncb0) deallocate(rdof)
      endif
      if (.not. allocated(rdof)) allocate(rdof(max(ncb0, 1)))
      i2 = 0
      do k2 = sym%rptr(s0) + nown2, sym%rptr(s0+1) - 1
        do a2 = 1, sym%ndof(sym%rlist(k2))
          i2 = i2 + 1
          rdof(i2) = fct%cdofptr(sym%rlist(k2)) + a2 - 1
        enddo
      enddo
      ! the delayed position offset of every child (chead order)
      if (allocated(mdel)) then
        if (size(mdel) < fct%cptr(s0+1) - fct%cptr(s0)) deallocate(mdel)
      endif
      if (.not. allocated(mdel)) allocate(mdel(max(fct%cptr(s0+1) - fct%cptr(s0), 1)))
      d2 = fct%cdofptr(sym%sptr(s0+1)) - fct%cdofptr(sym%sptr(s0))
      c = fct%chead(s0)
      do while (c /= 0)
        do j = fct%cptr(s0), fct%cptr(s0+1) - 1
          if (fct%clist(j) == c) mdel(j - fct%cptr(s0) + 1) = d2
        enddo
        d2 = d2 + fct%sn(c)%ncol - fct%sn(c)%npiv
        c = fct%cnext(c)
      enddo
    end subroutine geo_1d

    !> owning rank of a front row of the current front
    function drank_s(row) result(dr)
      integer(kind=kint), intent(in) :: row
      integer(kind=kint) :: dr

      if (row <= ncol) then
        dr = master
      else
        dr = towner(crow2t(row - ncol))
      endif
    end function drank_s

    !> row map of child c0 (clist ordinal jc): the owner of every contribution row, its
    !> parent front position and, for a fully summed destination, the permuted DOF the
    !> master assembles into
    subroutine crows_setup(s0, c0, jc)
      integer(kind=kint), intent(in) :: s0, c0, jc
      integer(kind=kint) :: rr, ll, kk2, aa, nde, nown2, cnown2

      cndel0 = fct%sn(c0)%ncol - fct%sn(c0)%npiv
      cnrow0 = cndel0 + fct%sn(c0)%tbnd(fct%sn(c0)%nt+1) - fct%sn(c0)%ncol
      if (allocated(pmapr)) then
        if (size(pmapr) < cnrow0) deallocate(pmapr, pdofm, crown)
      endif
      if (.not. allocated(pmapr)) allocate(pmapr(max(cnrow0, 1)), pdofm(max(cnrow0, 1)), crown(max(cnrow0, 1)))
      cncbt = 0
      if (map%upper(c0)) then
        nde = fct%sn(c0)%ncol - (fct%cdofptr(sym%sptr(c0+1)) - fct%cdofptr(sym%sptr(c0)))
        call hecmw_mf_dist_rowtiles(sym, map, fct%tile, c0, nde, cncbt, cctb2, cown2)
        do rr = 1, cnrow0 - cndel0
          do ll = 1, cncbt
            if (rr <= cctb2(ll)) exit
          enddo
          crown(cndel0 + rr) = cown2(ll)
        enddo
      else
        crown(cndel0+1:cnrow0) = map%owner(c0)
      endif
      crown(1:cndel0) = map%owner(c0)
      nown2 = sym%sptr(s0+1) - sym%sptr(s0)
      cnown2 = sym%sptr(c0+1) - sym%sptr(c0)
      do rr = 1, cndel0
        pmapr(rr) = mdel(jc) + rr
        pdofm(rr) = fct%sn(c0)%frow(fct%sn(c0)%npiv + rr)
      enddo
      rr = cndel0
      do ll = sym%rptr(c0) + cnown2, sym%rptr(c0+1) - 1
        kk2 = sym%cmap(sym%cmap_ptr(c0) + (ll - sym%rptr(c0) - cnown2))
        do aa = 1, sym%ndof(sym%rlist(ll))
          rr = rr + 1
          pmapr(rr) = rowoff(kk2) + aa
          if (kk2 <= nown2) then
            pdofm(rr) = fct%cdofptr(sym%rlist(sym%rptr(s0)+kk2-1)) + aa - 1
          else
            pdofm(rr) = 0
          endif
        enddo
      enddo
    end subroutine crows_setup

    !> slice my rows of the forward contribution of child c0 to the ranks holding the
    !> parent rows
    subroutine dv_send(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint) :: dst, rr, nsl, il, ib

      if (.not. allocated(fct%sn(c0)%dval)) return
      do dst = map%rbeg(s), map%rbeg(s) + map%rcnt(s) - 1
        if (dst == map%myrank) cycle
        nsl = 0
        do rr = 1, cnrow0
          if (crown(rr) == map%myrank .and. drank_s(pmapr(rr)) == dst) nsl = nsl + 1
        enddo
        if (nsl == 0) cycle
        ib = mf_pool_slot(pool)
        allocate(pool%box(ib)%rv(nsl))
        il = 0
        do rr = 1, cnrow0
          if (crown(rr) /= map%myrank) cycle
          if (drank_s(pmapr(rr)) /= dst) cycle
          il = il + 1
          pool%box(ib)%rv(il) = fct%sn(c0)%dval(dv_local(c0, rr))
        enddo
        req = 0
        call hecmw_isend_r(pool%box(ib)%rv, nsl, dst, 16*c0+8, map%comm, req)
        call mf_pool_req(pool, req)
      enddo
    end subroutine dv_send

    !> local dval index of child contribution row rr on this rank
    function dv_local(c0, rr) result(il)
      integer(kind=kint), intent(in) :: c0, rr
      integer(kind=kint) :: il, t2, ll

      if (.not. map%upper(c0)) then
        il = rr
      else if (rr <= cndel0) then
        il = rr
      else
        ! my compact row numbering over my child tiles, after the delayed rows when I am
        ! also the master of the child
        do t2 = 1, cncbt
          if (rr - cndel0 <= cctb2(t2)) exit
        enddo
        il = 0
        if (map%owner(c0) == map%myrank) il = cndel0
        do ll = 1, t2 - 1
          if (cown2(ll) == map%myrank) il = il + cctb2(ll) - cctb2(ll-1)
        enddo
        il = il + (rr - cndel0 - cctb2(t2-1))
      endif
    end function dv_local

    !> receive and add the forward contributions of child c0 (my own part straight from
    !> the local dval), the fully summed destinations through the permuted position map
    subroutine dv_apply(s0, c0)
      integer(kind=kint), intent(in) :: s0, c0
      real(kind=kreal), allocatable :: rb(:)
      integer(kind=kint) :: is0, src, nsl, il, rr, nsrc
      integer(kind=kint) :: srcs(cnrow0 + 1)

      nsrc = 0
      do rr = 1, cnrow0
        do is0 = 1, nsrc
          if (srcs(is0) == crown(rr)) exit
        enddo
        if (is0 > nsrc) then
          nsrc = nsrc + 1
          srcs(nsrc) = crown(rr)
        endif
      enddo
      do is0 = 1, nsrc
        src = srcs(is0)
        nsl = 0
        do rr = 1, cnrow0
          if (crown(rr) == src .and. drank_s(pmapr(rr)) == map%myrank) nsl = nsl + 1
        enddo
        if (nsl == 0) cycle
        if (src == map%myrank) then
          do rr = 1, cnrow0
            if (crown(rr) /= map%myrank) cycle
            if (drank_s(pmapr(rr)) /= map%myrank) cycle
            call dv_add(rr, fct%sn(c0)%dval(dv_local(c0, rr)))
          enddo
        else
          allocate(rb(nsl))
          call hecmw_recv_r(rb, nsl, src, 16*c0+8, map%comm, stat)
          il = 0
          do rr = 1, cnrow0
            if (crown(rr) /= src) cycle
            if (drank_s(pmapr(rr)) /= map%myrank) cycle
            il = il + 1
            call dv_add(rr, rb(il))
          enddo
          deallocate(rb)
        endif
      enddo
    end subroutine dv_apply

    !> add one forward contribution at child row rr into what I hold
    subroutine dv_add(rr, val)
      integer(kind=kint), intent(in) :: rr
      real(kind=kreal), intent(in) :: val
      integer(kind=kint) :: pr, t2

      pr = pmapr(rr)
      if (pr <= ncol) then
        v(dmap(pdofm(rr))) = v(dmap(pdofm(rr))) + val
      else
        t2 = crow2t(pr - ncol)
        vb(tprow(t2) + (pr - ncol - ctb(t2-1))) = vb(tprow(t2) + (pr - ncol - ctb(t2-1))) + val
      endif
    end subroutine dv_add

    !> forward solve of the upper front s0: the children's contributions arrive as row
    !> slices, the master runs the fully summed sweep and hands the tile owners the post
    !> substitution pivot values, and every owner applies its stored rows to its band rows
    subroutine fwd_1d(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: d2, d0, d1, kk, tw2, hk2, t2, h2, r2, i2, iw2, ofs, idx2, jc, ib, ndv
      integer(kind=8) :: okk2, oik2

      call geo_1d(s0)
      call mf_grow_r(v, int(max(ncol, 1), 8))
      call mf_grow_r(vb, int(max(myrows, 1), 8))
      call mf_grow_r(tbuf, int(max(ncol, nrow0 - ncol, 1), 8))
      if (ismaster) then
        d0 = fct%cdofptr(sym%sptr(s0))
        d1 = fct%cdofptr(sym%sptr(s0+1)) - 1
        do d2 = 1, ncol
          dmap(fct%sn(s0)%frow(d2)) = d2
          if (fct%sn(s0)%frow(d2) >= d0 .and. fct%sn(s0)%frow(d2) <= d1) then
            v(d2) = z(fct%sn(s0)%frow(d2))
          else
            v(d2) = 0.0d0
          endif
        enddo
      endif
      vb(1:max(myrows, 1)) = 0.0d0
      do j = fct%cptr(s0), fct%cptr(s0+1) - 1
        c = fct%clist(j)
        call crows_setup(s0, c, j - fct%cptr(s0) + 1)
        call dv_send(c)
      enddo
      c = fct%chead(s0)
      do while (c /= 0)
        jc = 0
        do j = fct%cptr(s0), fct%cptr(s0+1) - 1
          if (fct%clist(j) == c) jc = j - fct%cptr(s0) + 1
        enddo
        call crows_setup(s0, c, jc)
        call dv_apply(s0, c)
        c = fct%cnext(c)
      enddo
      call mf_grow_r(pivs, int(max(npv, 1), 8))
      if (ismaster) then
        ofs = 0
        do kk = 1, ntc0
          tw2 = min(fct%sn(s0)%tbnd(kk+1), npv) - fct%sn(s0)%tbnd(kk)
          if (tw2 <= 0) exit
          hk2 = fct%sn(s0)%tbnd(kk+1) - fct%sn(s0)%tbnd(kk)
          okk2 = 1 + fct%sn(s0)%bptr(mf_bidx(nt0, kk, kk))
          call hecmw_mf_kernel_trsv(tw2, hk2, fct%sn(s0)%lval(okk2), v(fct%sn(s0)%tbnd(kk)+1))
          if (hk2 > tw2) call hecmw_mf_kernel_gemv(hk2 - tw2, tw2, hk2, fct%sn(s0)%lval(okk2 + tw2), &
            v(fct%sn(s0)%tbnd(kk)+1), v(fct%sn(s0)%tbnd(kk)+tw2+1))
          pivs(ofs+1:ofs+tw2) = v(fct%sn(s0)%tbnd(kk)+1:fct%sn(s0)%tbnd(kk)+tw2)
          ofs = ofs + tw2
          do i2 = kk + 1, ntc0
            h2 = fct%sn(s0)%tbnd(i2+1) - fct%sn(s0)%tbnd(i2)
            oik2 = 1 + fct%sn(s0)%bptr(mf_bidx(nt0, kk, i2))
            call hecmw_mf_kernel_gemv(h2, tw2, h2, fct%sn(s0)%lval(oik2), v(fct%sn(s0)%tbnd(kk)+1), &
              v(fct%sn(s0)%tbnd(i2)+1))
          enddo
          if (.not. fct%lu) call hecmw_mf_kernel_dsolve(tw2, hk2, fct%sn(s0)%lval(okk2), &
            fct%sn(s0)%ptype(fct%sn(s0)%tbnd(kk)+1:), fct%sn(s0)%dsub(fct%sn(s0)%tbnd(kk)+1:), &
            v(fct%sn(s0)%tbnd(kk)+1))
        enddo
        do d2 = 1, npv
          z(fct%sn(s0)%frow(d2)) = v(d2)
        enddo
        if (npv > 0) then
          ib = mf_pool_slot(pool)
          allocate(pool%box(ib)%rv(npv))
          pool%box(ib)%rv(1:npv) = pivs(1:npv)
          do iw2 = 1, nwk
            if (wrank(iw2) == map%myrank) cycle
            req = 0
            call hecmw_isend_r(pool%box(ib)%rv, npv, wrank(iw2), 16*s0+9, map%comm, req)
            call mf_pool_req(pool, req)
          enddo
        endif
      else if (hasband .and. npv > 0) then
        call hecmw_recv_r(pivs, npv, master, 16*s0+9, map%comm, stat)
      endif
      if (hasband) then
        ofs = 0
        do kk = 1, ntc0
          tw2 = min(fct%sn(s0)%tbnd(kk+1), npv) - fct%sn(s0)%tbnd(kk)
          if (tw2 <= 0) exit
          do t2 = 1, ncbt
            if (towner(t2) /= map%myrank) cycle
            h2 = ctb(t2) - ctb(t2-1)
            idx2 = mf_bidx(nt0, kk, ntc0 + t2)
            r2 = fct%sn(s0)%brank(idx2)
            oik2 = 1 + fct%sn(s0)%bptr(idx2)
            if (r2 == 0) cycle
            if (r2 < 0) then
              call hecmw_mf_kernel_gemv(h2, tw2, h2, fct%sn(s0)%lval(oik2), pivs(ofs+1), vb(tprow(t2)+1))
            else
              call hecmw_mf_kernel_mult_tv(tw2, r2, tw2, fct%sn(s0)%lval(oik2 + int(h2, 8)*r2), pivs(ofs+1), tbuf)
              call hecmw_mf_kernel_gemv(h2, r2, h2, fct%sn(s0)%lval(oik2), tbuf, vb(tprow(t2)+1))
            endif
          enddo
          ofs = ofs + tw2
        enddo
      endif
      ! keep my rows of the forward contribution: the delayed rows on the master, the
      ! band rows on their owners
      ndv = myrows
      if (ismaster) ndv = ndv + (ncol - npv)
      if (allocated(fct%sn(s0)%dval)) deallocate(fct%sn(s0)%dval)
      if (ndv > 0) then
        allocate(fct%sn(s0)%dval(ndv))
        ofs = 0
        if (ismaster) then
          fct%sn(s0)%dval(1:ncol-npv) = v(npv+1:ncol)
          ofs = ncol - npv
        endif
        if (myrows > 0) fct%sn(s0)%dval(ofs+1:ofs+myrows) = vb(1:myrows)
      endif
      call mf_pool_wait(pool)
      do j = fct%cptr(s0), fct%cptr(s0+1) - 1
        c = fct%clist(j)
        if (allocated(fct%sn(c)%dval)) deallocate(fct%sn(c)%dval)
      enddo
    end subroutine fwd_1d

    !> backward solve of the upper front s0: the master receives the referenced solution
    !> values from the parent and hands every tile owner the values of its rows; the
    !> right hand side of a tile column then travels along the tile owners in ascending
    !> order, each subtracting its stored rows, which keeps the sequential fold order
    subroutine bwd_1d(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: d2, kk, tw2, hk2, t2, h2, i2, iw2, ib, nr2, il
      integer(kind=8) :: okk2, okku2, oik2

      call geo_1d(s0)
      ! a rank without the master role or tiles has no part in the backward solve
      if (.not. (ismaster .or. hasband)) return
      call mf_grow_r(v, int(max(ncol, 1), 8))
      call mf_grow_r(vb, int(max(myrows, 1), 8))
      call mf_grow_r(tok, int(max(ncol, 1), 8))
      call mf_grow_r(tbuf, int(max(ncol, nrow0 - ncol, 1), 8))
      if (ismaster) then
        p = sym%sparent(s0)
        if (p /= 0) then
          if (map%owner(p) /= map%myrank) call recv_xx(s0)
        endif
        do iw2 = 1, nwk
          if (wrank(iw2) == map%myrank) cycle
          nr2 = 0
          do t2 = 1, ncbt
            if (towner(t2) == wrank(iw2)) nr2 = nr2 + ctb(t2) - ctb(t2-1)
          enddo
          if (nr2 == 0) cycle
          ib = mf_pool_slot(pool)
          allocate(pool%box(ib)%rv(nr2))
          il = 0
          do t2 = 1, ncbt
            if (towner(t2) /= wrank(iw2)) cycle
            do d2 = ctb(t2-1) + 1, ctb(t2)
              il = il + 1
              pool%box(ib)%rv(il) = xx(rdof(d2))
            enddo
          enddo
          req = 0
          call hecmw_isend_r(pool%box(ib)%rv, nr2, wrank(iw2), 16*s0+10, map%comm, req)
          call mf_pool_req(pool, req)
        enddo
        do d2 = 1, npv
          v(d2) = z(fct%sn(s0)%frow(d2))
        enddo
        do d2 = npv + 1, ncol
          v(d2) = xx(fct%sn(s0)%fsdof(d2))
        enddo
      endif
      if (hasband) then
        if (ismaster) then
          do t2 = 1, ncbt
            if (towner(t2) /= map%myrank) cycle
            do d2 = ctb(t2-1) + 1, ctb(t2)
              vb(tprow(t2) + (d2 - ctb(t2-1))) = xx(rdof(d2))
            enddo
          enddo
        else
          call hecmw_recv_r(vb, myrows, master, 16*s0+10, map%comm, stat)
        endif
      endif
      do kk = ntc0, 1, -1
        tw2 = min(fct%sn(s0)%tbnd(kk+1), npv) - fct%sn(s0)%tbnd(kk)
        if (tw2 <= 0) cycle
        hk2 = fct%sn(s0)%tbnd(kk+1) - fct%sn(s0)%tbnd(kk)
        if (ismaster) then
          tok(1:tw2) = v(fct%sn(s0)%tbnd(kk)+1:fct%sn(s0)%tbnd(kk)+tw2)
          do i2 = kk + 1, ntc0
            h2 = fct%sn(s0)%tbnd(i2+1) - fct%sn(s0)%tbnd(i2)
            if (fct%lu) then
              oik2 = 1 + fct%sn(s0)%bptru(mf_bidx(nt0, kk, i2))
              call hecmw_mf_kernel_gemv_t(h2, tw2, h2, fct%sn(s0)%uval(oik2), v(fct%sn(s0)%tbnd(i2)+1), tok)
            else
              oik2 = 1 + fct%sn(s0)%bptr(mf_bidx(nt0, kk, i2))
              call hecmw_mf_kernel_gemv_t(h2, tw2, h2, fct%sn(s0)%lval(oik2), v(fct%sn(s0)%tbnd(i2)+1), tok)
            endif
          enddo
        endif
        call token_chain(s0, kk, tw2)
        if (ismaster) then
          okk2 = 1 + fct%sn(s0)%bptr(mf_bidx(nt0, kk, kk))
          if (fct%lu) then
            okku2 = 1 + fct%sn(s0)%bptru(mf_bidx(nt0, kk, kk))
            if (hk2 > tw2) call hecmw_mf_kernel_gemv_t(hk2 - tw2, tw2, hk2, fct%sn(s0)%uval(okku2 + tw2), &
              v(fct%sn(s0)%tbnd(kk)+tw2+1), tok)
            call hecmw_mf_kernel_usolve(tw2, hk2, fct%sn(s0)%lval(okk2), hk2, fct%sn(s0)%uval(okku2), tok)
          else
            if (hk2 > tw2) call hecmw_mf_kernel_gemv_t(hk2 - tw2, tw2, hk2, fct%sn(s0)%lval(okk2 + tw2), &
              v(fct%sn(s0)%tbnd(kk)+tw2+1), tok)
            call hecmw_mf_kernel_trsv_t(tw2, hk2, fct%sn(s0)%lval(okk2), tok)
          endif
          v(fct%sn(s0)%tbnd(kk)+1:fct%sn(s0)%tbnd(kk)+tw2) = tok(1:tw2)
        endif
      enddo
      if (ismaster) then
        do d2 = 1, npv
          xx(fct%sn(s0)%fsdof(d2)) = v(d2)
        enddo
        do j = fct%cptr(s0), fct%cptr(s0+1) - 1
          c = fct%clist(j)
          if (map%owner(c) /= map%myrank) call send_xx(c)
        enddo
      endif
      call mf_pool_wait(pool)
    end subroutine bwd_1d

    !> the running right hand side of tile column kk along the tile owners in ascending
    !> order; each owner subtracts its stored rows in place, keeping the sequential fold
    subroutine token_chain(s0, kk, tw2)
      integer(kind=kint), intent(in) :: s0, kk, tw2
      integer(kind=kint) :: iw2, prv, nxt, ib

      if (nwk == 0) return
      if (ismaster) then
        prv = -1
        do iw2 = 1, nwk
          if (wrank(iw2) == map%myrank) then
            if (prv >= 0) then
              call hecmw_recv_r(tok, tw2, prv, 16*s0+11, map%comm, stat)
              prv = -1
            endif
            call chain_apply(s0, kk, tw2)
          else
            if (prv < 0) then
              ib = mf_pool_slot(pool)
              allocate(pool%box(ib)%rv(tw2))
              pool%box(ib)%rv(1:tw2) = tok(1:tw2)
              req = 0
              call hecmw_isend_r(pool%box(ib)%rv, tw2, wrank(iw2), 16*s0+11, map%comm, req)
              call mf_pool_req(pool, req)
            endif
            prv = wrank(iw2)
          endif
        enddo
        if (prv >= 0) call hecmw_recv_r(tok, tw2, prv, 16*s0+11, map%comm, stat)
      else
        if (wme == 1) then
          prv = master
        else
          prv = wrank(wme-1)
        endif
        call hecmw_recv_r(tok, tw2, prv, 16*s0+11, map%comm, stat)
        call chain_apply(s0, kk, tw2)
        if (wme == nwk) then
          nxt = master
        else
          nxt = wrank(wme+1)
        endif
        ib = mf_pool_slot(pool)
        allocate(pool%box(ib)%rv(tw2))
        pool%box(ib)%rv(1:tw2) = tok(1:tw2)
        req = 0
        call hecmw_isend_r(pool%box(ib)%rv, tw2, nxt, 16*s0+11, map%comm, req)
        call mf_pool_req(pool, req)
      endif
    end subroutine token_chain

    !> subtract my stored rows of tile column kk from the running right hand side
    subroutine chain_apply(s0, kk, tw2)
      integer(kind=kint), intent(in) :: s0, kk, tw2
      integer(kind=kint) :: t2, h2, r2, idx2
      integer(kind=8) :: oik2

      do t2 = 1, ncbt
        if (towner(t2) /= map%myrank) cycle
        h2 = ctb(t2) - ctb(t2-1)
        idx2 = mf_bidx(nt0, kk, ntc0 + t2)
        if (fct%lu) then
          r2 = fct%sn(s0)%branku(idx2)
          oik2 = 1 + fct%sn(s0)%bptru(idx2)
          if (r2 == 0) cycle
          if (r2 < 0) then
            call hecmw_mf_kernel_gemv_t(h2, tw2, h2, fct%sn(s0)%uval(oik2), vb(tprow(t2)+1), tok)
          else
            call hecmw_mf_kernel_mult_tv(h2, r2, h2, fct%sn(s0)%uval(oik2), vb(tprow(t2)+1), tbuf)
            call hecmw_mf_kernel_gemv(tw2, r2, tw2, fct%sn(s0)%uval(oik2 + int(h2, 8)*r2), tbuf, tok)
          endif
        else
          r2 = fct%sn(s0)%brank(idx2)
          oik2 = 1 + fct%sn(s0)%bptr(idx2)
          if (r2 == 0) cycle
          if (r2 < 0) then
            call hecmw_mf_kernel_gemv_t(h2, tw2, h2, fct%sn(s0)%lval(oik2), vb(tprow(t2)+1), tok)
          else
            call hecmw_mf_kernel_mult_tv(h2, r2, h2, fct%sn(s0)%lval(oik2), vb(tprow(t2)+1), tbuf)
            call hecmw_mf_kernel_gemv(tw2, r2, tw2, fct%sn(s0)%lval(oik2 + int(h2, 8)*r2), tbuf, tok)
          endif
        endif
      enddo
    end subroutine chain_apply

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

    !> isend the referenced solution values to the owner of child c0; the send stays
    !> pending until the end of the backward stage, when the owner has received it
    subroutine send_xx(c0)
      integer(kind=kint), intent(in) :: c0
      integer(kind=kint), allocatable :: list(:)
      integer(kind=kint) :: mm, i0, ib

      call xx_dofs(c0, mm, list)
      if (mm <= 0) return
      ib = mf_pool_slot(gpool)
      allocate(gpool%box(ib)%rv(mm))
      do i0 = 1, mm
        gpool%box(ib)%rv(i0) = xx(list(i0))
      enddo
      deallocate(list)
      req = 0
      call hecmw_isend_r(gpool%box(ib)%rv, mm, map%owner(c0), 16*c0+12, map%comm, req)
      call mf_pool_req(gpool, req)
    end subroutine send_xx

    !> receive the referenced solution values of s0 from the owner of its parent
    subroutine recv_xx(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint), allocatable :: list(:)
      real(kind=kreal), allocatable :: rv(:)
      integer(kind=kint) :: mm, i0

      call xx_dofs(s0, mm, list)
      if (mm <= 0) return
      allocate(rv(mm))
      call hecmw_recv_r(rv, mm, map%owner(sym%sparent(s0)), 16*s0+12, map%comm, stat)
      do i0 = 1, mm
        xx(list(i0)) = rv(i0)
      enddo
      deallocate(rv, list)
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

  !> Stored panel words of supernode s held on this rank (test support).
  function hecmw_mf_numeric_front_words(fct, s) result(w)
    implicit none
    type(hecmwST_mf_factor), intent(in) :: fct
    integer(kind=kint), intent(in) :: s
    integer(kind=8) :: w

    w = fct%sn(s)%pwords
  end function hecmw_mf_numeric_front_words

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
