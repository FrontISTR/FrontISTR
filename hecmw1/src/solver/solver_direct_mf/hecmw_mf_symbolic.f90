!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Symbolic factorization (elimination tree, column counts, supernodes, row structure)
!>        for the multifrontal direct solver, on the node level graph.
!>
!> References:
!>  Liu, "The role of elimination trees in sparse factorization", SIAM J. Matrix Anal. Appl. 11 (1990)
!>  Gilbert, Ng, Peyton, "An efficient algorithm to compute row and column counts for sparse
!>    Cholesky factorization", SIAM J. Matrix Anal. Appl. 15 (1994)
!>  Ashcraft, Grimes, "The influence of relaxed supernode partitions on the multifrontal method",
!>    ACM TOMS 15 (1989)
module hecmw_mf_symbolic
  use hecmw_util
  use hecmw_mf_graph
  implicit none

  private
  public :: hecmwST_mf_symbolic
  public :: hecmw_mf_symbolic_build
  public :: hecmw_mf_symbolic_cluster
  public :: hecmw_mf_symbolic_check
  public :: hecmw_mf_symbolic_print
  public :: hecmw_mf_symbolic_finalize

  !> All indices are 1-based. "Column" means a node of the graph in the eliminated order,
  !> which is a topological order of the elimination tree (parent(k) > k); every supernode
  !> owns a contiguous range of columns and supernodes are numbered in postorder of the
  !> supernodal tree. Without relaxed merges the column order is a postorder as well.
  !> The numeric stage reads this structure only; it never modifies it.
  type hecmwST_mf_symbolic
    integer(kind=kint) :: nnode = 0
    integer(kind=kint) :: nsuper = 0
    integer(kind=kint) :: nmerge = 0                !< supernode merges done by relaxed amalgamation
    integer(kind=kint), allocatable :: perm(:)      !< perm(k): original node eliminated at column k
    integer(kind=kint), allocatable :: invp(:)      !< invp(i): column of original node i
    integer(kind=kint), allocatable :: ndof(:)      !< ndof(k): DOFs of column k
    integer(kind=kint), allocatable :: parent(:)    !< elimination tree over columns, 0 at roots
    integer(kind=kint), allocatable :: colcnt(:)    !< nodes in L(:,k) including k (exact, no relaxation)
    integer(kind=kint), allocatable :: sptr(:)      !< supernode s owns columns sptr(s):sptr(s+1)-1
    integer(kind=kint), allocatable :: sparent(:)   !< supernodal elimination tree, 0 at roots
    integer(kind=kint), allocatable :: rptr(:)      !< rows of supernode s are rlist(rptr(s):rptr(s+1)-1)
    integer(kind=kint), allocatable :: rlist(:)     !< strictly ascending; the own columns come first,
                                                    !< then the rows of the contribution block
    integer(kind=kint), allocatable :: cmap_ptr(:)  !< cmap(cmap_ptr(c):cmap_ptr(c+1)-1) maps the
    integer(kind=kint), allocatable :: cmap(:)      !< contribution rows of child c, in rlist order, to
                                                    !< positions (1-based) within the parent's rlist
    integer(kind=kint) :: max_front = 0             !< largest frontal matrix size in DOFs
    integer(kind=8) :: factor_nnz = 0               !< DOF entries of L including the diagonal blocks
    integer(kind=8) :: flops = 0                    !< sum over pivots of (front rows remaining)^2
  end type hecmwST_mf_symbolic

contains

  !> perm(k) is the original node placed at position k by a fill reducing ordering; any
  !> topological reordering of the elimination tree gives the same factor structure, so the
  !> stored perm differs from the input by postordering and amalgamation only.
  !> relax_nodes: a child supernode is merged into its parent while the merged column count
  !> (in nodes) does not exceed this value; 0 keeps the fundamental supernodes.
  subroutine hecmw_mf_symbolic_build(graph, perm, relax_nodes, sym)
    implicit none
    type(hecmwST_mf_graph), intent(in) :: graph
    integer(kind=kint), intent(in) :: perm(:)
    integer(kind=kint), intent(in) :: relax_nodes
    type(hecmwST_mf_symbolic), intent(out) :: sym
    integer(kind=kint), allocatable :: post(:)
    integer(kind=kint) :: n, k

    n = graph%nnode
    sym%nnode = n
    allocate(sym%perm(n), sym%invp(n), sym%ndof(n), sym%parent(n), sym%colcnt(n))
    allocate(post(n))
    sym%perm(1:n) = perm(1:n)
    do k = 1, n
      sym%invp(sym%perm(k)) = k
    enddo

    call mf_etree(graph, sym%perm, sym%invp, sym%parent)
    call mf_postorder(n, sym%parent, post)
    call mf_relabel(n, post, sym%perm, sym%invp, sym%parent)
    call mf_colcounts(graph, sym%perm, sym%invp, sym%parent, sym%colcnt)
    call mf_supernodes(n, sym%parent, sym%colcnt, relax_nodes, sym%nsuper, sym%sptr, sym%sparent, post, sym%nmerge)
    call mf_relabel(n, post, sym%perm, sym%invp, sym%parent, sym%colcnt)
    do k = 1, n
      sym%ndof(k) = graph%ndof(sym%perm(k))
    enddo
    call mf_rowstruct(graph, sym)
    call mf_estimate(sym)
    deallocate(post)
  end subroutine hecmw_mf_symbolic_build

  subroutine hecmw_mf_symbolic_finalize(sym)
    implicit none
    type(hecmwST_mf_symbolic), intent(inout) :: sym

    sym%nnode = 0
    sym%nsuper = 0
    if (allocated(sym%perm)) deallocate(sym%perm)
    if (allocated(sym%invp)) deallocate(sym%invp)
    if (allocated(sym%ndof)) deallocate(sym%ndof)
    if (allocated(sym%parent)) deallocate(sym%parent)
    if (allocated(sym%colcnt)) deallocate(sym%colcnt)
    if (allocated(sym%sptr)) deallocate(sym%sptr)
    if (allocated(sym%sparent)) deallocate(sym%sparent)
    if (allocated(sym%rptr)) deallocate(sym%rptr)
    if (allocated(sym%rlist)) deallocate(sym%rlist)
    if (allocated(sym%cmap_ptr)) deallocate(sym%cmap_ptr)
    if (allocated(sym%cmap)) deallocate(sym%cmap)
  end subroutine hecmw_mf_symbolic_finalize

  !> Elimination tree of the permuted matrix (Liu, The role of elimination trees in sparse
  !> factorization, 1990): column k adopts as its children the roots of every subtree that
  !> row k touches below the diagonal. The forest built so far is kept as disjoint sets
  !> whose representative is the subtree root; the root search halves the paths it walks.
  subroutine mf_etree(graph, perm, invp, parent)
    implicit none
    type(hecmwST_mf_graph), intent(in) :: graph
    integer(kind=kint), intent(in) :: perm(:)
    integer(kind=kint), intent(in) :: invp(:)
    integer(kind=kint), intent(out) :: parent(:)
    integer(kind=kint), allocatable :: setp(:)
    integer(kind=kint) :: n, k, u, r, p

    n = graph%nnode
    allocate(setp(n))
    parent(1:n) = 0
    do k = 1, n
      setp(k) = k
    enddo
    do k = 1, n
      do p = graph%xadj(perm(k)), graph%xadj(perm(k)+1)-1
        u = invp(graph%adjncy(p))
        if (u >= k) cycle
        r = u
        do while (setp(r) /= r)
          setp(r) = setp(setp(r))
          r = setp(r)
        enddo
        if (r /= k) then
          parent(r) = k
          setp(r) = k
        endif
      enddo
    enddo
    deallocate(setp)
  end subroutine mf_etree

  !> post(k): the column visited k-th in a depth first postorder (children in increasing order).
  subroutine mf_postorder(n, parent, post)
    implicit none
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(in) :: parent(:)
    integer(kind=kint), intent(out) :: post(:)
    integer(kind=kint), allocatable :: cptr(:), clist(:), cur(:), stk(:)
    integer(kind=kint) :: j, k, p, v, c, top

    ! the children of v are clist(cptr(v-1)+1 : cptr(v)), ascending because the fill
    ! walks j upward; cur(v) is the next child of v the walk has not descended into
    allocate(cptr(0:n), clist(n), cur(n), stk(n))
    cptr(0:n) = 0
    do j = 1, n
      p = parent(j)
      if (p /= 0) cptr(p) = cptr(p) + 1
    enddo
    do j = 1, n
      cptr(j) = cptr(j-1) + cptr(j)
    enddo
    cur(1:n) = 0
    do j = 1, n
      p = parent(j)
      if (p == 0) cycle
      cur(p) = cur(p) + 1
      clist(cptr(p-1) + cur(p)) = j
    enddo
    k = 0
    do j = 1, n
      if (parent(j) /= 0) cycle
      top = 1
      stk(1) = j
      cur(j) = cptr(j-1)
      do while (top > 0)
        v = stk(top)
        if (cur(v) < cptr(v)) then
          cur(v) = cur(v) + 1
          c = clist(cur(v))
          top = top + 1
          stk(top) = c
          cur(c) = cptr(c-1)
        else
          top = top - 1
          k = k + 1
          post(k) = v
        endif
      enddo
    enddo
    deallocate(cptr, clist, cur, stk)
  end subroutine mf_postorder

  !> Renumber columns so that old column post(k) becomes column k.
  subroutine mf_relabel(n, post, perm, invp, parent, colcnt)
    implicit none
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(in) :: post(:)
    integer(kind=kint), intent(inout) :: perm(:)
    integer(kind=kint), intent(inout) :: invp(:)
    integer(kind=kint), intent(inout) :: parent(:)
    integer(kind=kint), intent(inout), optional :: colcnt(:)
    integer(kind=kint), allocatable :: newpos(:), work(:)
    integer(kind=kint) :: k

    allocate(newpos(n), work(n))
    do k = 1, n
      newpos(post(k)) = k
    enddo
    do k = 1, n
      work(k) = perm(post(k))
    enddo
    perm(1:n) = work(1:n)
    do k = 1, n
      invp(perm(k)) = k
    enddo
    do k = 1, n
      if (parent(post(k)) == 0) then
        work(k) = 0
      else
        work(k) = newpos(parent(post(k)))
      endif
    enddo
    parent(1:n) = work(1:n)
    if (present(colcnt)) then
      do k = 1, n
        work(k) = colcnt(post(k))
      enddo
      colcnt(1:n) = work(1:n)
    endif
    deallocate(newpos, work)
  end subroutine mf_relabel

  !> Column counts of L (Gilbert, Ng and Peyton, An efficient algorithm to compute row and
  !> column counts for sparse Cholesky factorization, 1994); parent must be postordered.
  !> A row subtree is counted through its skeleton: an entry (i, j) of the strict lower
  !> part contributes a new path to count j exactly when j is a leaf of row i's subtree,
  !> and the overlap of two successive leaves is taken back at their least common
  !> ancestor. The counts of a column then sum up the tree.
  subroutine mf_colcounts(graph, perm, invp, parent, colcnt)
    implicit none
    type(hecmwST_mf_graph), intent(in) :: graph
    integer(kind=kint), intent(in) :: perm(:)
    integer(kind=kint), intent(in) :: invp(:)
    integer(kind=kint), intent(in) :: parent(:)
    integer(kind=kint), intent(out) :: colcnt(:)
    integer(kind=kint), allocatable :: fdesc(:), maxfd(:), lastleaf(:), setp(:)
    integer(kind=kint) :: n, j, i, p, q, w

    n = graph%nnode
    allocate(fdesc(n), maxfd(n), lastleaf(n), setp(n))
    ! with a postordered tree the subtree of j is the contiguous range fdesc(j)..j, so
    ! the first descendant is the minimum over the children (available when j is reached)
    do j = 1, n
      fdesc(j) = j
    enddo
    do j = 1, n
      p = parent(j)
      if (p /= 0) fdesc(p) = min(fdesc(p), fdesc(j))
    enddo
    do j = 1, n
      if (fdesc(j) == j) then
        colcnt(j) = 1
      else
        colcnt(j) = 0
      endif
      maxfd(j) = 0
      lastleaf(j) = 0
      setp(j) = j
    enddo
    do j = 1, n
      p = parent(j)
      ! the parent sees j's whole subtree except j itself through the summing pass
      if (p /= 0) colcnt(p) = colcnt(p) - 1
      do i = graph%xadj(perm(j)), graph%xadj(perm(j)+1)-1
        q = invp(graph%adjncy(i))
        if (q <= j) cycle
        ! j is a leaf of row q's skeleton when no earlier leaf reached into j's subtree
        if (fdesc(j) <= maxfd(q)) cycle
        maxfd(q) = fdesc(j)
        colcnt(j) = colcnt(j) + 1
        if (lastleaf(q) /= 0) then
          w = lca(lastleaf(q))
          colcnt(w) = colcnt(w) - 1
        endif
        lastleaf(q) = j
      enddo
      ! j's set joins the parent: the representative of a processed column is the lowest
      ! unprocessed ancestor, which is the least common ancestor with any later leaf
      if (p /= 0) setp(j) = p
    enddo
    do j = 1, n
      if (parent(j) /= 0) colcnt(parent(j)) = colcnt(parent(j)) + colcnt(j)
    enddo
    deallocate(fdesc, maxfd, lastleaf, setp)

  contains

    !> set representative of v, halving the path it walks
    integer(kind=kint) function lca(v)
      integer(kind=kint), intent(in) :: v
      integer(kind=kint) :: r
      r = v
      do while (setp(r) /= r)
        setp(r) = setp(setp(r))
        r = setp(r)
      enddo
      lca = r
    end function lca

  end subroutine mf_colcounts

  !> Fundamental supernodes, relaxed amalgamation, and the column renumbering post(:) that
  !> makes every merged supernode contiguous and the supernodes postordered.
  subroutine mf_supernodes(n, parent, colcnt, relax_nodes, nsuper, sptr, sparent, post, nmerge)
    implicit none
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(in) :: parent(:)
    integer(kind=kint), intent(in) :: colcnt(:)
    integer(kind=kint), intent(in) :: relax_nodes
    integer(kind=kint), intent(out) :: nsuper
    integer(kind=kint), allocatable, intent(out) :: sptr(:)
    integer(kind=kint), allocatable, intent(out) :: sparent(:)
    integer(kind=kint), intent(out) :: post(:)
    integer(kind=kint), intent(out) :: nmerge
    integer(kind=kint), allocatable :: nchild(:), superof(:), fptr(:), fparent(:), rep(:), size(:)
    integer(kind=kint), allocatable :: mid(:), mcols_ptr(:), mcols(:), mparent(:), morder(:), newid(:)
    integer(kind=kint) :: nfund, k, s, p, m, nm, l, c

    allocate(nchild(n), superof(n))
    nchild(1:n) = 0
    do k = 1, n
      if (parent(k) /= 0) nchild(parent(k)) = nchild(parent(k)) + 1
    enddo
    nfund = 0
    do k = 1, n
      if (k > 1) then
        if (parent(k-1) == k .and. colcnt(k-1) == colcnt(k) + 1 .and. nchild(k) == 1) then
          superof(k) = nfund
          cycle
        endif
      endif
      nfund = nfund + 1
      superof(k) = nfund
    enddo
    allocate(fptr(nfund+1), fparent(nfund), rep(nfund), size(nfund))
    fptr(1) = 1
    do k = 1, n
      fptr(superof(k)+1) = k + 1
    enddo
    do s = 1, nfund
      p = parent(fptr(s+1)-1)
      if (p == 0) then
        fparent(s) = 0
      else
        fparent(s) = superof(p)
      endif
      rep(s) = s
      size(s) = fptr(s+1) - fptr(s)
    enddo

    ! greedy bottom-up merge; children are processed before their parent by the postorder
    nmerge = 0
    do s = 1, nfund
      p = fparent(s)
      if (p == 0) cycle
      if (size(s) + size(p) <= relax_nodes) then
        rep(s) = p
        size(p) = size(p) + size(s)
        nmerge = nmerge + 1
      endif
    enddo

    ! merged supernodes are the fixed points of rep; mid numbers them in increasing order
    allocate(mid(nfund))
    nm = 0
    do s = 1, nfund
      if (rep(s) == s) then
        nm = nm + 1
        mid(s) = nm
      else
        mid(s) = 0
      endif
    enddo
    do s = 1, nfund
      p = s
      do while (rep(p) /= p)
        p = rep(p)
      enddo
      mid(s) = mid(p)
    enddo
    allocate(mparent(nm), mcols_ptr(nm+1), mcols(n), morder(nm))
    mcols_ptr(1:nm+1) = 0
    do k = 1, n
      m = mid(superof(k))
      mcols_ptr(m+1) = mcols_ptr(m+1) + 1
    enddo
    mcols_ptr(1) = 1
    do m = 1, nm
      mcols_ptr(m+1) = mcols_ptr(m) + mcols_ptr(m+1)
    enddo
    do k = 1, n
      m = mid(superof(k))
      mcols(mcols_ptr(m)) = k
      mcols_ptr(m) = mcols_ptr(m) + 1
    enddo
    do m = nm, 2, -1
      mcols_ptr(m) = mcols_ptr(m-1)
    enddo
    mcols_ptr(1) = 1
    do s = 1, nfund
      if (rep(s) /= s) cycle
      if (fparent(s) == 0) then
        mparent(mid(s)) = 0
      else
        mparent(mid(s)) = mid(fparent(s))
      endif
    enddo
    call mf_postorder(nm, mparent, morder)

    nsuper = nm
    allocate(sptr(nsuper+1), sparent(nsuper), newid(nm))
    do m = 1, nm
      newid(morder(m)) = m
    enddo
    l = 0
    sptr(1) = 1
    do m = 1, nm
      c = morder(m)
      do k = mcols_ptr(c), mcols_ptr(c+1)-1
        l = l + 1
        post(l) = mcols(k)
      enddo
      sptr(m+1) = l + 1
      if (mparent(c) == 0) then
        sparent(m) = 0
      else
        sparent(m) = newid(mparent(c))
      endif
    enddo
    deallocate(nchild, superof, fptr, fparent, rep, size, mid, mcols_ptr, mcols, mparent, morder, newid)
  end subroutine mf_supernodes

  !> Row structure of each supernode by merging the graph rows with the children's structure,
  !> and the child-to-parent row maps used by extend-add.
  subroutine mf_rowstruct(graph, sym)
    implicit none
    type(hecmwST_mf_graph), intent(in) :: graph
    type(hecmwST_mf_symbolic), intent(inout) :: sym
    integer(kind=kint), allocatable :: mark(:), head(:), next(:), buf(:), pos(:)
    integer(kind=kint) :: n, ns, s, c, k, fc, lc, p, i, cnt, cap, l, ncol

    n = sym%nnode
    ns = sym%nsuper
    allocate(mark(n), head(ns), next(ns), buf(n), pos(n))
    allocate(sym%rptr(ns+1), sym%cmap_ptr(ns+1))
    mark(1:n) = 0
    head(1:ns) = 0
    do s = ns, 1, -1
      p = sym%sparent(s)
      if (p == 0) cycle
      next(s) = head(p)
      head(p) = s
    enddo
    cap = 0
    do k = 1, n
      cap = cap + sym%colcnt(k)
    enddo
    allocate(sym%rlist(cap))
    sym%rptr(1) = 1
    do s = 1, ns
      fc = sym%sptr(s)
      lc = sym%sptr(s+1) - 1
      cnt = 0
      do k = fc, lc
        cnt = cnt + 1
        buf(cnt) = k
        mark(k) = s
      enddo
      do k = fc, lc
        do p = graph%xadj(sym%perm(k)), graph%xadj(sym%perm(k)+1)-1
          i = sym%invp(graph%adjncy(p))
          if (i <= lc) cycle
          if (mark(i) == s) cycle
          mark(i) = s
          cnt = cnt + 1
          buf(cnt) = i
        enddo
      enddo
      c = head(s)
      do while (c /= 0)
        ncol = sym%sptr(c+1) - sym%sptr(c)
        do l = sym%rptr(c) + ncol, sym%rptr(c+1) - 1
          i = sym%rlist(l)
          if (i <= lc) cycle
          if (mark(i) == s) cycle
          mark(i) = s
          cnt = cnt + 1
          buf(cnt) = i
        enddo
        c = next(c)
      enddo
      call mf_sort(buf(lc-fc+2:cnt), cnt-(lc-fc+1))
      if (sym%rptr(s) + cnt - 1 > cap) then
        cap = max(2*cap, sym%rptr(s) + cnt - 1)
        call mf_grow(sym%rlist, cap)
      endif
      sym%rlist(sym%rptr(s):sym%rptr(s)+cnt-1) = buf(1:cnt)
      sym%rptr(s+1) = sym%rptr(s) + cnt
    enddo

    sym%cmap_ptr(1) = 1
    do c = 1, ns
      ncol = sym%sptr(c+1) - sym%sptr(c)
      sym%cmap_ptr(c+1) = sym%cmap_ptr(c) + (sym%rptr(c+1) - sym%rptr(c) - ncol)
    enddo
    allocate(sym%cmap(sym%cmap_ptr(ns+1)-1))
    do s = 1, ns
      if (head(s) == 0) cycle
      do l = sym%rptr(s), sym%rptr(s+1) - 1
        pos(sym%rlist(l)) = l - sym%rptr(s) + 1
      enddo
      c = head(s)
      do while (c /= 0)
        ncol = sym%sptr(c+1) - sym%sptr(c)
        k = sym%cmap_ptr(c)
        do l = sym%rptr(c) + ncol, sym%rptr(c+1) - 1
          sym%cmap(k) = pos(sym%rlist(l))
          k = k + 1
        enddo
        c = next(c)
      enddo
    enddo
    deallocate(mark, head, next, buf, pos)
  end subroutine mf_rowstruct

  subroutine mf_grow(a, cap)
    implicit none
    integer(kind=kint), allocatable, intent(inout) :: a(:)
    integer(kind=kint), intent(in) :: cap
    integer(kind=kint), allocatable :: tmp(:)

    allocate(tmp(cap))
    tmp(1:size(a)) = a(1:size(a))
    call move_alloc(tmp, a)
  end subroutine mf_grow

  !> In-place heap sort of a(1:m) in ascending order.
  subroutine mf_sort(a, m)
    implicit none
    integer(kind=kint), intent(inout) :: a(:)
    integer(kind=kint), intent(in) :: m
    integer(kind=kint) :: i, j, c, t

    do i = m/2, 1, -1
      j = i
      t = a(j)
      do
        c = 2*j
        if (c > m) exit
        if (c < m) then
          if (a(c+1) > a(c)) c = c + 1
        endif
        if (a(c) <= t) exit
        a(j) = a(c)
        j = c
      enddo
      a(j) = t
    enddo
    do i = m, 2, -1
      t = a(i)
      a(i) = a(1)
      j = 1
      do
        c = 2*j
        if (c > i-1) exit
        if (c < i-1) then
          if (a(c+1) > a(c)) c = c + 1
        endif
        if (a(c) <= t) exit
        a(j) = a(c)
        j = c
      enddo
      a(j) = t
    enddo
  end subroutine mf_sort

  !> Factor size and work in DOF units. Each supernode is a dense front of nrow DOF rows with
  !> ncol pivot DOFs; pivot k (0-based) stores nrow-k entries and costs (nrow-k)^2 multiply-adds
  !> on the full front (an LDLt kernel performs about half of that).
  subroutine mf_estimate(sym)
    implicit none
    type(hecmwST_mf_symbolic), intent(inout) :: sym
    integer(kind=kint) :: s, l, k
    integer(kind=8) :: nrow, ncol, a, b

    sym%factor_nnz = 0
    sym%flops = 0
    sym%max_front = 0
    do s = 1, sym%nsuper
      ncol = 0
      do k = sym%sptr(s), sym%sptr(s+1) - 1
        ncol = ncol + sym%ndof(k)
      enddo
      nrow = 0
      do l = sym%rptr(s), sym%rptr(s+1) - 1
        nrow = nrow + sym%ndof(sym%rlist(l))
      enddo
      sym%max_front = max(sym%max_front, int(nrow, kind=kint))
      sym%factor_nnz = sym%factor_nnz + ncol*nrow - ncol*(ncol-1)/2
      a = nrow
      b = nrow - ncol
      sym%flops = sym%flops + (a*(a+1)*(2*a+1) - b*(b+1)*(2*b+1))/6
    enddo
  end subroutine mf_estimate

  !> Consistency checks of a built structure against its graph; every violation is reported on
  !> unit 6 and counted in nerr.
  subroutine hecmw_mf_symbolic_check(graph, sym, nerr)
    implicit none
    type(hecmwST_mf_graph), intent(in) :: graph
    type(hecmwST_mf_symbolic), intent(in) :: sym
    integer(kind=kint), intent(out) :: nerr
    integer(kind=kint), allocatable :: fdesc(:), dsize(:), superof(:), mark(:)
    integer(kind=kint) :: n, ns, k, s, p, c, l, i, fc, lc, ncol, nbad
    integer(kind=8) :: nnz_node, nnz_cnt
    type(hecmwST_mf_symbolic) :: tmp

    nerr = 0
    n = sym%nnode
    ns = sym%nsuper
    allocate(fdesc(n), dsize(n), superof(n), mark(n))

    nbad = 0
    do k = 1, n
      if (sym%perm(k) < 1 .or. sym%perm(k) > n) then
        nbad = nbad + 1
      else if (sym%invp(sym%perm(k)) /= k) then
        nbad = nbad + 1
      endif
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: perm/invp are not inverse, violations = ', nbad
    endif

    nbad = 0
    do k = 1, n
      fdesc(k) = k
      dsize(k) = 0
    enddo
    do k = 1, n
      p = sym%parent(k)
      if (p == 0) cycle
      if (p <= k .or. p > n) then
        nbad = nbad + 1
        cycle
      endif
      fdesc(p) = min(fdesc(p), fdesc(k))
      dsize(p) = dsize(p) + (k - fdesc(k) + 1)
    enddo
    if (sym%nmerge == 0) then
      do k = 1, n
        if (dsize(k) /= k - fdesc(k)) nbad = nbad + 1
      enddo
    endif
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: parent is not a topologically ordered tree, violations = ', nbad
    endif

    nbad = 0
    if (sym%sptr(1) /= 1 .or. sym%sptr(ns+1) /= n+1) nbad = nbad + 1
    do s = 1, ns
      if (sym%sptr(s+1) <= sym%sptr(s)) nbad = nbad + 1
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a)') '[DIRECTmf] check: sptr does not partition the columns'
      deallocate(fdesc, dsize, superof, mark)
      return
    endif
    do s = 1, ns
      superof(sym%sptr(s):sym%sptr(s+1)-1) = s
    enddo
    do s = 1, ns
      lc = sym%sptr(s+1) - 1
      do k = sym%sptr(s), lc - 1
        if (sym%parent(k) == 0) then
          nbad = nbad + 1
        else if (superof(sym%parent(k)) /= s) then
          nbad = nbad + 1
        endif
      enddo
      if (sym%parent(lc) == 0) then
        if (sym%sparent(s) /= 0) nbad = nbad + 1
      else
        if (sym%sparent(s) /= superof(sym%parent(lc))) nbad = nbad + 1
        if (sym%sparent(s) <= s) nbad = nbad + 1
      endif
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: supernodes inconsistent with the elimination tree, violations = ', nbad
    endif

    nbad = 0
    do s = 1, ns
      fc = sym%sptr(s)
      lc = sym%sptr(s+1) - 1
      ncol = lc - fc + 1
      if (sym%rptr(s+1) - sym%rptr(s) < ncol) then
        nbad = nbad + 1
        cycle
      endif
      do l = sym%rptr(s), sym%rptr(s+1) - 1
        i = sym%rlist(l)
        if (l - sym%rptr(s) < ncol) then
          if (i /= fc + (l - sym%rptr(s))) nbad = nbad + 1
        else
          if (i <= lc .or. i > n) nbad = nbad + 1
          if (l > sym%rptr(s)) then
            if (i <= sym%rlist(l-1)) nbad = nbad + 1
          endif
        endif
      enddo
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: rlist not ascending or not led by the own columns, violations = ', nbad
    endif

    ! every entry of A below the diagonal lies in the structure of L
    nbad = 0
    mark(1:n) = 0
    do s = 1, ns
      do l = sym%rptr(s), sym%rptr(s+1) - 1
        mark(sym%rlist(l)) = s
      enddo
      do k = sym%sptr(s), sym%sptr(s+1) - 1
        do p = graph%xadj(sym%perm(k)), graph%xadj(sym%perm(k)+1)-1
          i = sym%invp(graph%adjncy(p))
          if (i > k .and. mark(i) /= s) nbad = nbad + 1
        enddo
      enddo
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: graph edges missing from rlist, violations = ', nbad
    endif

    ! colcnt(k) equals the rows of the supernode from k down when no relaxation was applied
    nbad = 0
    nnz_node = 0
    nnz_cnt = 0
    do s = 1, ns
      fc = sym%sptr(s)
      lc = sym%sptr(s+1) - 1
      do k = fc, lc
        c = (sym%rptr(s+1) - sym%rptr(s)) - (k - fc)
        nnz_node = nnz_node + c
        nnz_cnt = nnz_cnt + sym%colcnt(k)
        if (sym%colcnt(k) > c) nbad = nbad + 1
        if (sym%nmerge == 0 .and. sym%colcnt(k) /= c) nbad = nbad + 1
      enddo
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: colcnt inconsistent with rlist, violations = ', nbad
    endif
    if (nnz_cnt > nnz_node .or. (sym%nmerge == 0 .and. nnz_cnt /= nnz_node)) then
      nerr = nerr + 1
      write(*,'(a,i0,a,i0)') '[DIRECTmf] check: sum(colcnt) = ', nnz_cnt, ' vs rlist node entries = ', nnz_node
    endif

    nbad = 0
    do c = 1, ns
      p = sym%sparent(c)
      ncol = sym%sptr(c+1) - sym%sptr(c)
      if (sym%cmap_ptr(c+1) - sym%cmap_ptr(c) /= sym%rptr(c+1) - sym%rptr(c) - ncol) then
        nbad = nbad + 1
        cycle
      endif
      if (p == 0) cycle
      k = sym%cmap_ptr(c)
      do l = sym%rptr(c) + ncol, sym%rptr(c+1) - 1
        i = sym%cmap(k)
        if (i < 1 .or. i > sym%rptr(p+1) - sym%rptr(p)) then
          nbad = nbad + 1
        else if (sym%rlist(sym%rptr(p) + i - 1) /= sym%rlist(l)) then
          nbad = nbad + 1
        endif
        k = k + 1
      enddo
    enddo
    if (nbad > 0) then
      nerr = nerr + 1
      write(*,'(a,i0)') '[DIRECTmf] check: cmap does not map onto the parent rlist, violations = ', nbad
    endif

    tmp%nnode = n
    tmp%nsuper = ns
    allocate(tmp%ndof(n), tmp%sptr(ns+1), tmp%rptr(ns+1), tmp%rlist(sym%rptr(ns+1)-1))
    tmp%ndof(1:n) = sym%ndof(1:n)
    tmp%sptr(1:ns+1) = sym%sptr(1:ns+1)
    tmp%rptr(1:ns+1) = sym%rptr(1:ns+1)
    tmp%rlist(1:sym%rptr(ns+1)-1) = sym%rlist(1:sym%rptr(ns+1)-1)
    call mf_estimate(tmp)
    if (tmp%factor_nnz /= sym%factor_nnz .or. tmp%flops /= sym%flops .or. tmp%max_front /= sym%max_front) then
      nerr = nerr + 1
      write(*,'(a)') '[DIRECTmf] check: stored estimates differ from recomputation'
    endif
    call hecmw_mf_symbolic_finalize(tmp)
    deallocate(fdesc, dsize, superof, mark)
  end subroutine hecmw_mf_symbolic_check

  subroutine hecmw_mf_symbolic_print(sym)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    integer(kind=kint) :: k, ndof_min, ndof_max
    integer(kind=8) :: ndof_tot

    ndof_tot = 0
    ndof_min = huge(1_kint)
    ndof_max = 0
    do k = 1, sym%nnode
      ndof_tot = ndof_tot + sym%ndof(k)
      ndof_min = min(ndof_min, sym%ndof(k))
      ndof_max = max(ndof_max, sym%ndof(k))
    enddo
    write(*,'(a,i0,a,i0,a,i0,a,i0,a)') '[DIRECTmf]: nodes = ', sym%nnode, ', dofs = ', ndof_tot, &
      ' (ndof ', ndof_min, '..', ndof_max, ')'
    write(*,'(a,i0,a,i0,a,i0)') '[DIRECTmf]: supernodes = ', sym%nsuper, ', relaxed merges = ', sym%nmerge, &
      ', max front (dof) = ', sym%max_front
    write(*,'(a,i0,a,f10.3,a)') '[DIRECTmf]: factor nnz (dof) = ', sym%factor_nnz, &
      ' (', real(sym%factor_nnz, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
    write(*,'(a,i0,a,f12.3,a)') '[DIRECTmf]: factorization flops = ', sym%flops, &
      ' (', real(sym%flops, kind=kreal)/1.0d9, ' GFLOP)'
  end subroutine hecmw_mf_symbolic_print

  !> Refined permutation for the BLR compression: the own nodes of every supernode are
  !> reordered by a BFS on the induced subgraph, restarted per connected component from an
  !> unvisited node of least internal degree, so that geometrically close separator nodes
  !> stay in neighboring columns and the panel tiles keep low ranks. The order across
  !> supernodes is unchanged (any order within a supernode is a valid elimination order);
  !> the caller rebuilds the symbolic structure with the returned permutation.
  subroutine hecmw_mf_symbolic_cluster(graph, sym, perm)
    implicit none
    type(hecmwST_mf_graph), intent(in) :: graph
    type(hecmwST_mf_symbolic), intent(in) :: sym
    integer(kind=kint), intent(out) :: perm(:)
    integer(kind=kint), allocatable :: own(:), deg(:), queue(:)
    integer(kind=kint) :: s, k, i, j, p, pos, head, tail, seed, mindeg

    allocate(own(graph%nnode), deg(graph%nnode), queue(graph%nnode))
    do s = 1, sym%nsuper
      do k = sym%sptr(s), sym%sptr(s+1) - 1
        own(sym%perm(k)) = s
      enddo
    enddo
    pos = 0
    do s = 1, sym%nsuper
      do k = sym%sptr(s), sym%sptr(s+1) - 1
        i = sym%perm(k)
        deg(i) = 0
        do p = graph%xadj(i), graph%xadj(i+1) - 1
          if (own(graph%adjncy(p)) == s) deg(i) = deg(i) + 1
        enddo
      enddo
      do
        seed = 0
        mindeg = huge(seed)
        do k = sym%sptr(s), sym%sptr(s+1) - 1
          i = sym%perm(k)
          if (own(i) == s .and. deg(i) < mindeg) then
            mindeg = deg(i)
            seed = i
          endif
        enddo
        if (seed == 0) exit
        own(seed) = -s
        queue(1) = seed
        head = 1
        tail = 1
        do while (head <= tail)
          i = queue(head)
          head = head + 1
          pos = pos + 1
          perm(pos) = i
          do p = graph%xadj(i), graph%xadj(i+1) - 1
            j = graph%adjncy(p)
            if (own(j) == s) then
              own(j) = -s
              tail = tail + 1
              queue(tail) = j
            endif
          enddo
        enddo
      enddo
    enddo
    deallocate(own, deg, queue)
  end subroutine hecmw_mf_symbolic_cluster

end module hecmw_mf_symbolic
