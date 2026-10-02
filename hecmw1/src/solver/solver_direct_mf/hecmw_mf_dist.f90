!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief MPI distribution layer of the multifrontal solver: the supernode-to-rank mapping
!>        and the broadcast of the symbolic structure. Every rank computes the mapping from
!>        the (identical) symbolic structure alone, without communication.
module hecmw_mf_dist
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_mf_symbolic
  implicit none

  private
  public :: hecmwST_mf_map
  public :: hecmw_mf_dist_map_build
  public :: hecmw_mf_dist_map_finalize
  public :: hecmw_mf_dist_rowtiles
  public :: hecmw_mf_dist_fsgrid
  public :: hecmw_mf_dist_fsrank
  public :: hecmw_mf_dist_symbolic_bcast
  public :: hecmwST_mf_gmat
  public :: hecmw_mf_dist_gmat_build
  public :: hecmw_mf_dist_gmat_part
  public :: hecmw_mf_dist_gmat_vals
  public :: hecmw_mf_dist_gmat_extract
  public :: hecmw_mf_dist_gmat_words
  public :: hecmw_mf_dist_gmat_finalize
  public :: hecmw_mf_dist_gather_vec

  !> Subtree-to-subcube mapping over the supernodal tree: every supernode carries the rank
  !> set rbeg/rcnt that processes it; a supernode whose rank set spans more than one rank
  !> is an upper front, whose fully summed tiles are block cyclic on the 2D process grid of
  !> the front (hecmw_mf_dist_fsgrid, the master coordinating and holding the metadata)
  !> while the contribution row tiles are distributed over the rank set (1D row
  !> distribution). Below a single-rank supernode the whole subtree belongs to that rank.
  type hecmwST_mf_map
    integer(kind=kint) :: nprocs = 1
    integer(kind=kint) :: myrank = 0
    integer(kind=kint) :: comm = 0
    integer(kind=kint) :: nupper = 0
    integer(kind=kint), allocatable :: owner(:)   !< owner(s): rank that factors the fully summed part of s
    logical, allocatable :: upper(:)              !< the rank set of s spans more than one rank
    integer(kind=kint), allocatable :: uplist(:)  !< upper supernodes, ascending (children first)
    integer(kind=kint), allocatable :: rbeg(:)    !< first rank of the rank set of s
    integer(kind=kint), allocatable :: rcnt(:)    !< ranks in the rank set of s
  end type hecmwST_mf_map

  !> The part of the global matrix this rank holds, assembled from the internal rows of
  !> the distributed matrix (complete in the overlapped HEC-MW assembly). The global node
  !> numbering concatenates the internal nodes rank by rank, the way the sparse matrix
  !> interface of the external direct solvers numbers them. The structure is gathered once
  !> (hecmw_mf_dist_gmat_build), the full profile living on rank 0 only until the ordering
  !> has read it; the partial profile and the value routing are then derived from the
  !> mapping (hecmw_mf_dist_gmat_part), and the values move to their holding ranks before
  !> every factorization (hecmw_mf_dist_gmat_vals). The block stream of a rank lists, row
  !> by row, the diagonal block followed by the lower then upper neighbor blocks in the
  !> local item order, and both the send and the receive lists enumerate it in that order.
  type hecmwST_mf_gmat
    type(hecmwST_matrix) :: mat                   !< the partial matrix (full index arrays, retained blocks)
    integer(kind=kint) :: nblk = 0                !< global structure stream length in nd*nd blocks
    integer(kind=kint), allocatable :: nn(:)      !< internal nodes of rank r at nn(r+1)
    integer(kind=kint), allocatable :: ndisp(:)   !< global node offset of rank r at ndisp(r+1)
    integer(kind=kint), allocatable :: vblk(:)    !< structure stream blocks per rank
    integer(kind=kint), allocatable :: stream(:)  !< gathered structure stream, freed by gmat_part
    integer(kind=kint), allocatable :: xptr(:)    !< per global row the missing transposed columns (0:ng), freed by gmat_part
    integer(kind=kint), allocatable :: xcol(:)    !< their column ids, ascending within a row
    integer(kind=kint), allocatable :: gid(:)     !< user node id of every global row, for messages
    integer(kind=kint), allocatable :: scnt(:)    !< blocks sent to rank r at scnt(r+1)
    integer(kind=kint), allocatable :: ssel(:)    !< sent blocks as local stream indices, grouped by destination
    integer(kind=kint), allocatable :: rcnt(:)    !< blocks received from rank r at rcnt(r+1)
    integer(kind=kint), allocatable :: rdst(:)    !< received block slot: 1..NPL in AL, then NPU in AU, then row in D
  end type hecmwST_mf_gmat

  !> A subtree heavier than this multiple of the mean rank load of its front's range is
  !> promoted to an upper front instead of packed onto a single rank. Lowering it refines
  !> the packing granularity at the cost of more upper fronts (protocol overhead).
  real(kind=kreal), parameter :: MF_MAP_PROMOTE = 1.2d0

  !> the promotion threshold of a subtree at least MF_MAP_PROMOTE_ABS heavy (in the weight
  !> units of wt, flops plus MF_MAP_WPF per factor word): the parallelism gained by
  !> refining such a subtree outweighs the per-front protocol cost, so the packing is
  !> refined more aggressively there, while light subtrees keep the conservative threshold
  real(kind=kreal), parameter :: MF_MAP_PROMOTE_HEAVY = 0.5d0
  integer(kind=8), parameter :: MF_MAP_PROMOTE_ABS = 500000000000_8

  !> flops charged per factor panel word in the subtree weights (the memory-bound cost of
  !> assembling and storing a word, relative to one flop of elimination)
  integer(kind=8), parameter :: MF_MAP_WPF = 2500_8

  !> fully summed tiles per rank of the 2D process grid of a front: the grid takes one rank
  !> per this many tiles, so a small fully summed part stays on few ranks (a 1 x 1 grid
  !> keeps it on the master) and the protocol overhead of spreading it is gated by size
  integer(kind=kint), parameter :: MF_2D_FSTPR = 8

contains

  !> Build the mapping from the symbolic structure. The roots share all ranks (several
  !> roots split them proportionally to the subtree weights, factorization flops
  !> estimates); an upper front assigns the subtrees below it by promotion and packing: a
  !> subtree heavier than the promotion threshold becomes an upper front on the same rank
  !> range and hands its children down to the pool, and the remaining subtrees are packed
  !> greedily onto the least loaded rank of the range (LPT). The result is deterministic
  !> and identical on every rank. The owner of an upper front is the owner of its heaviest
  !> child subtree, keeping the fan-in communication local to the subcube.
  subroutine hecmw_mf_dist_map_build(sym, map)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(out) :: map
    integer(kind=8), allocatable :: wt(:)
    integer(kind=kint), allocatable :: cptr(:), clist(:), wptr(:)
    integer(kind=kint) :: ns, s, p, c, i, nroot, nup, cbest
    integer(kind=8) :: w, wbest, b

    ns = sym%nsuper
    map%nprocs = hecmw_comm_get_size()
    map%myrank = hecmw_comm_get_rank()
    map%comm = hecmw_comm_get_comm()
    allocate(map%owner(ns), map%upper(ns))

    ! subtree weights: the factorization flops estimate of the supernode (the formula of
    ! mf_estimate) plus the factor panel words scaled by MF_MAP_WPF, accumulated bottom-up
    ! (the ascending numbering places children before parents); the words term charges the
    ! memory-bound assembly and store work of the many small fronts, which the flops alone
    ! underrate
    allocate(wt(ns))
    do s = 1, ns
      w = 0
      do i = sym%rptr(s), sym%rptr(s+1) - 1
        w = w + sym%ndof(sym%rlist(i))
      enddo
      b = w - sym_cdofcount(s)
      wt(s) = (w*(w+1)*(2*w+1) - b*(b+1)*(2*b+1))/6 + MF_MAP_WPF*w*sym_cdofcount(s)
    enddo
    do s = 1, ns
      p = sym%sparent(s)
      if (p /= 0) wt(p) = wt(p) + wt(s)
    enddo

    ! children in ascending order
    allocate(cptr(ns+1), clist(max(ns, 1)), wptr(ns))
    cptr(1:ns+1) = 0
    do s = 1, ns
      p = sym%sparent(s)
      if (p /= 0) cptr(p+1) = cptr(p+1) + 1
    enddo
    cptr(1) = 1
    do s = 1, ns
      cptr(s+1) = cptr(s) + cptr(s+1)
    enddo
    wptr(1:ns) = cptr(1:ns)
    do s = 1, ns
      p = sym%sparent(s)
      if (p /= 0) then
        clist(wptr(p)) = s
        wptr(p) = wptr(p) + 1
      endif
    enddo

    ! rank ranges: the roots share all ranks, an upper front promotes and packs the
    ! subtrees below it within its range; parents carry larger indices, so a descending
    ! sweep sets every range before it is read, and a front promoted higher up finds its
    ! children already assigned
    allocate(map%rbeg(ns), map%rcnt(ns))
    map%rcnt(1:ns) = -1
    nroot = 0
    do s = 1, ns
      if (sym%sparent(s) == 0) nroot = nroot + 1
    enddo
    i = 0
    do s = 1, ns
      if (sym%sparent(s) == 0) then
        i = i + 1
        wptr(i) = s
      endif
    enddo
    call mf_map_split(wt, wptr(1:nroot), 0, map%nprocs, map%rbeg, map%rcnt)
    do s = ns, 1, -1
      if (map%rcnt(s) < 0) then
        p = sym%sparent(s)
        map%rbeg(s) = map%rbeg(p)
        map%rcnt(s) = 1
      endif
      if (map%rcnt(s) > 1) call mf_map_pack(s)
    enddo

    ! owners bottom-up: a single-rank supernode is owned by its rank, an upper front by the
    ! owner of its heaviest child subtree (smallest child index on ties)
    nup = 0
    do s = 1, ns
      if (map%rcnt(s) <= 1) then
        map%owner(s) = map%rbeg(s)
        map%upper(s) = .false.
      else
        cbest = 0
        wbest = -1
        do i = cptr(s), cptr(s+1) - 1
          c = clist(i)
          if (wt(c) > wbest) then
            wbest = wt(c)
            cbest = c
          endif
        enddo
        if (cbest > 0) then
          map%owner(s) = map%owner(cbest)
        else
          map%owner(s) = map%rbeg(s)
        endif
        map%upper(s) = .true.
        nup = nup + 1
      endif
    enddo
    map%nupper = nup
    allocate(map%uplist(max(nup, 1)))
    nup = 0
    do s = 1, ns
      if (map%upper(s)) then
        nup = nup + 1
        map%uplist(nup) = s
      endif
    enddo
    deallocate(wt, cptr, clist, wptr)

  contains

    !> DOFs of the own columns of supernode s
    function sym_cdofcount(s0) result(nc)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint) :: nc, k
      nc = 0
      do k = sym%sptr(s0), sym%sptr(s0+1) - 1
        nc = nc + sym%ndof(k)
      enddo
    end function sym_cdofcount

    !> Assign the yet unassigned subtrees under upper front s0 to the ranks of its range.
    !> A pooled subtree heavier than the promotion threshold times the mean rank load of
    !> the pool (MF_MAP_PROMOTE, or MF_MAP_PROMOTE_HEAVY for a subtree at least
    !> MF_MAP_PROMOTE_ABS heavy) is promoted to an upper front on the same range and
    !> replaced by its children; its own front work leaves the pool as the distribution
    !> spreads it over the whole range, so the threshold shrinks and the sweep repeats
    !> until stable. The remaining subtrees go heaviest first onto the least loaded rank
    !> (LPT).
    subroutine mf_map_pack(s0)
      integer(kind=kint), intent(in) :: s0
      integer(kind=kint), allocatable :: pool(:)
      integer(kind=8), allocatable :: load(:)
      integer(kind=kint) :: npool, i0, j0, c0, r0, nr0, rmin
      integer(kind=8) :: tot
      real(kind=kreal) :: th
      logical :: grew

      r0 = map%rbeg(s0)
      nr0 = map%rcnt(s0)
      allocate(pool(ns), load(0:nr0-1))
      npool = 0
      tot = 0
      do i0 = cptr(s0), cptr(s0+1) - 1
        c0 = clist(i0)
        if (map%rcnt(c0) < 0) then
          npool = npool + 1
          pool(npool) = c0
          tot = tot + wt(c0)
        endif
      enddo
      grew = .true.
      do while (grew)
        grew = .false.
        i0 = 1
        do while (i0 <= npool)
          c0 = pool(i0)
          th = MF_MAP_PROMOTE
          if (wt(c0) >= MF_MAP_PROMOTE_ABS) th = MF_MAP_PROMOTE_HEAVY
          if (real(wt(c0), kind=kreal)*nr0 > th*real(tot, kind=kreal)) then
            map%rbeg(c0) = r0
            map%rcnt(c0) = nr0
            tot = tot - wt(c0)
            pool(i0) = pool(npool)
            npool = npool - 1
            do j0 = cptr(c0), cptr(c0+1) - 1
              npool = npool + 1
              pool(npool) = clist(j0)
              tot = tot + wt(clist(j0))
            enddo
            grew = .true.
          else
            i0 = i0 + 1
          endif
        enddo
      enddo
      call mf_sort_lpt(pool, npool)
      load(0:nr0-1) = 0
      do i0 = 1, npool
        rmin = 0
        do j0 = 1, nr0 - 1
          if (load(j0) < load(rmin)) rmin = j0
        enddo
        c0 = pool(i0)
        map%rbeg(c0) = r0 + rmin
        map%rcnt(c0) = 1
        ! a zero weight still loads its rank one unit, so tiny subtrees spread out
        load(rmin) = load(rmin) + max(wt(c0), 1_8)
      enddo
      deallocate(pool, load)
    end subroutine mf_map_pack

    !> heapsort of list(1:n0) into the LPT order: descending subtree weight, ties by
    !> ascending supernode index
    subroutine mf_sort_lpt(list, n0)
      integer(kind=kint), intent(inout) :: list(:)
      integer(kind=kint), intent(in) :: n0
      integer(kind=kint) :: i0, j0, k0, t0

      do i0 = n0/2, 1, -1
        j0 = i0
        do
          k0 = 2*j0
          if (k0 > n0) exit
          if (k0 < n0) then
            if (lpt_after(list(k0+1), list(k0))) k0 = k0 + 1
          endif
          if (.not. lpt_after(list(k0), list(j0))) exit
          t0 = list(j0)
          list(j0) = list(k0)
          list(k0) = t0
          j0 = k0
        enddo
      enddo
      do i0 = n0, 2, -1
        t0 = list(1)
        list(1) = list(i0)
        list(i0) = t0
        j0 = 1
        do
          k0 = 2*j0
          if (k0 > i0 - 1) exit
          if (k0 < i0 - 1) then
            if (lpt_after(list(k0+1), list(k0))) k0 = k0 + 1
          endif
          if (.not. lpt_after(list(k0), list(j0))) exit
          t0 = list(j0)
          list(j0) = list(k0)
          list(k0) = t0
          j0 = k0
        enddo
      enddo
    end subroutine mf_sort_lpt

    !> supernode a0 follows b0 in the LPT order (lighter subtree, larger index on ties)
    logical function lpt_after(a0, b0)
      integer(kind=kint), intent(in) :: a0, b0
      lpt_after = (wt(a0) < wt(b0)) .or. (wt(a0) == wt(b0) .and. a0 > b0)
    end function lpt_after

  end subroutine hecmw_mf_dist_map_build

  !> Split the rank range [r0, r0+np) among the given supernodes proportionally to their
  !> subtree weights, by cumulative integer arithmetic (deterministic). A zero or tiny
  !> weight still gets one rank, which may be shared with a neighbor.
  subroutine mf_map_split(wt, list, r0, np, rbeg, rcnt)
    implicit none
    integer(kind=8), intent(in) :: wt(:)
    integer(kind=kint), intent(in) :: list(:)
    integer(kind=kint), intent(in) :: r0, np
    integer(kind=kint), intent(inout) :: rbeg(:), rcnt(:)
    integer(kind=8) :: ctot, c
    integer(kind=kint) :: i, s, b0, b1

    ctot = 0
    do i = 1, size(list)
      ctot = ctot + wt(list(i))
    enddo
    c = 0
    do i = 1, size(list)
      s = list(i)
      if (ctot > 0) then
        ! rounded boundaries: a floor would push both halves of an even split onto the
        ! first rank whenever the cut misses an integer
        b0 = int((2_8*np*c + ctot)/(2_8*ctot), kind=kint)
        c = c + wt(s)
        b1 = int((2_8*np*c + ctot)/(2_8*ctot), kind=kint)
      else
        b0 = 0
        b1 = 0
      endif
      if (b0 > np - 1) b0 = np - 1
      if (b1 <= b0) b1 = b0 + 1
      rbeg(s) = r0 + b0
      rcnt(s) = min(b1, np) - b0
    enddo
  end subroutine mf_map_split

  subroutine hecmw_mf_dist_map_finalize(map)
    implicit none
    type(hecmwST_mf_map), intent(inout) :: map

    map%nprocs = 1
    map%nupper = 0
    if (allocated(map%owner)) deallocate(map%owner)
    if (allocated(map%upper)) deallocate(map%upper)
    if (allocated(map%uplist)) deallocate(map%uplist)
    if (allocated(map%rbeg)) deallocate(map%rbeg)
    if (allocated(map%rcnt)) deallocate(map%rcnt)
  end subroutine hecmw_mf_dist_map_finalize

  !> Contribution row tiles of upper front s and their owning ranks. The boundaries ctb
  !> (0:ncbt) are offsets from the end of the fully summed part; they reproduce the front
  !> partition of the numeric stage, whose cuts beyond the fully summed part shift with the
  !> delayed growth ndel without changing, so senders and receivers of a contribution slice
  !> derive the same tiles before the front exists. Contiguous tile blocks go to the ranks
  !> of s proportionally to the trailing update weight rows x band width of a tile:
  !> contiguous, unlike a cyclic assignment, keeps the owner changes along an ascending
  !> tile walk (the backward solve token chain) at one per rank.
  subroutine hecmw_mf_dist_rowtiles(sym, map, tile, s, ndel, ncbt, ctb, towner)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    integer(kind=kint), intent(in) :: tile, s, ndel
    integer(kind=kint), intent(out) :: ncbt
    integer(kind=kint), allocatable, intent(inout) :: ctb(:)
    integer(kind=kint), allocatable, intent(inout) :: towner(:)
    integer(kind=kint) :: nown, nrow_nodes, ncol, i, k, nd, rel, cur, nr, t, b
    integer(kind=8) :: ctot, c, w

    nown = sym%sptr(s+1) - sym%sptr(s)
    nrow_nodes = sym%rptr(s+1) - sym%rptr(s)
    ncol = ndel
    do k = sym%sptr(s), sym%sptr(s+1) - 1
      ncol = ncol + sym%ndof(k)
    enddo
    if (allocated(ctb)) then
      if (size(ctb) < nrow_nodes - nown + 1) deallocate(ctb, towner)
    endif
    if (.not. allocated(ctb)) allocate(ctb(0:max(nrow_nodes - nown, 1)), towner(max(nrow_nodes - nown, 1)))
    ncbt = 0
    ctb(0) = 0
    rel = 0
    cur = 0
    do i = nown + 1, nrow_nodes
      nd = sym%ndof(sym%rlist(sym%rptr(s)+i-1))
      if (rel > cur .and. rel + nd - cur > tile) then
        ncbt = ncbt + 1
        ctb(ncbt) = rel
        cur = rel
      endif
      rel = rel + nd
    enddo
    if (rel > cur .or. ncbt == 0) then
      if (rel > 0) then
        ncbt = ncbt + 1
        ctb(ncbt) = rel
      endif
    endif
    if (ncbt == 0) return

    nr = map%rcnt(s)
    ctot = 0
    do t = 1, ncbt
      ctot = ctot + int(ctb(t) - ctb(t-1), 8)*(ncol + ctb(t))
    enddo
    c = 0
    do t = 1, ncbt
      w = int(ctb(t) - ctb(t-1), 8)*(ncol + ctb(t))
      b = int(((2_8*c + w)*nr)/(2_8*ctot), kind=kint)
      if (b > nr - 1) b = nr - 1
      towner(t) = map%rbeg(s) + b
      c = c + w
    enddo
  end subroutine hecmw_mf_dist_rowtiles

  !> 2D process grid pr x pc of the fully summed part of upper front s, from its fully
  !> summed tile count ntc (known to every rank of s after the child headers). The grid
  !> grows with the tile count at MF_2D_FSTPR tiles per rank up to the rank set of s and
  !> stays near square with pr >= pc; a 1 x 1 grid keeps the fully summed part on the
  !> master, reproducing the undistributed layout.
  subroutine hecmw_mf_dist_fsgrid(map, s, ntc, pr, pc)
    implicit none
    type(hecmwST_mf_map), intent(in) :: map
    integer(kind=kint), intent(in) :: s, ntc
    integer(kind=kint), intent(out) :: pr, pc
    integer(kind=kint) :: ng

    ng = (ntc*(ntc+1)/2) / MF_2D_FSTPR
    if (ng > map%rcnt(s)) ng = map%rcnt(s)
    if (ng < 1) ng = 1
    pc = 1
    do while ((pc+1)*(pc+1) <= ng)
      pc = pc + 1
    enddo
    pr = ng / pc
  end subroutine hecmw_mf_dist_fsgrid

  !> Owning rank of fully summed tile (i, j), j <= i, of upper front s on its pr x pc
  !> grid: block cyclic over the tile coordinates, the grid ranks being the ranks of s in
  !> ascending order with the master moved to the front (so grid position 0 is the master).
  function hecmw_mf_dist_fsrank(map, s, pr, pc, i, j) result(r)
    implicit none
    type(hecmwST_mf_map), intent(in) :: map
    integer(kind=kint), intent(in) :: s, pr, pc, i, j
    integer(kind=kint) :: r, gp

    gp = mod(i-1, pr)*pc + mod(j-1, pc)
    if (gp == 0) then
      r = map%owner(s)
    else
      r = map%rbeg(s) + gp - 1
      if (r >= map%owner(s)) r = r + 1
    endif
  end function hecmw_mf_dist_fsrank

  !> Broadcast the symbolic structure built on the root rank to all ranks. The int8
  !> estimates (factor_nnz, flops) are not transferred: nothing reads them off the root.
  subroutine hecmw_mf_dist_symbolic_bcast(sym, root)
    implicit none
    type(hecmwST_mf_symbolic), intent(inout) :: sym
    integer(kind=kint), intent(in) :: root
    integer(kind=kint) :: hdr(6)
    integer(kind=kint) :: n, ns, nrl, ncm, comm

    if (hecmw_comm_get_size() == 1) return
    comm = hecmw_comm_get_comm()
    if (hecmw_comm_get_rank() == root) then
      hdr(1) = sym%nnode
      hdr(2) = sym%nsuper
      hdr(3) = sym%nmerge
      hdr(4) = sym%max_front
      hdr(5) = sym%rptr(sym%nsuper+1) - 1
      hdr(6) = sym%cmap_ptr(sym%nsuper+1) - 1
    endif
    call hecmw_bcast_I_comm(hdr, 6, root, comm)
    n = hdr(1)
    ns = hdr(2)
    nrl = hdr(5)
    ncm = hdr(6)
    if (hecmw_comm_get_rank() /= root) then
      call hecmw_mf_symbolic_finalize(sym)
      sym%nnode = n
      sym%nsuper = ns
      sym%nmerge = hdr(3)
      sym%max_front = hdr(4)
      allocate(sym%perm(n), sym%invp(n), sym%ndof(n), sym%parent(n), sym%colcnt(n))
      allocate(sym%sptr(ns+1), sym%sparent(ns), sym%rptr(ns+1), sym%rlist(nrl))
      allocate(sym%cmap_ptr(ns+1), sym%cmap(max(ncm, 1)))
    endif
    call hecmw_bcast_I_comm(sym%perm, n, root, comm)
    call hecmw_bcast_I_comm(sym%invp, n, root, comm)
    call hecmw_bcast_I_comm(sym%ndof, n, root, comm)
    call hecmw_bcast_I_comm(sym%parent, n, root, comm)
    call hecmw_bcast_I_comm(sym%colcnt, n, root, comm)
    call hecmw_bcast_I_comm(sym%sptr, ns + 1, root, comm)
    call hecmw_bcast_I_comm(sym%sparent, ns, root, comm)
    call hecmw_bcast_I_comm(sym%rptr, ns + 1, root, comm)
    call hecmw_bcast_I_comm(sym%rlist, nrl, root, comm)
    call hecmw_bcast_I_comm(sym%cmap_ptr, ns + 1, root, comm)
    if (ncm > 0) call hecmw_bcast_I_comm(sym%cmap, ncm, root, comm)
  end subroutine hecmw_mf_dist_symbolic_bcast

  !> Gather the nonzero structure of the internal rows of every rank. The stream carries,
  !> per row, the neighbor count and the global column ids, and is kept in gmat until
  !> hecmw_mf_dist_gmat_part derives the partial profile and the value routing from it;
  !> only rank 0 builds the full matrix profile (structure only, the columns of a row
  !> sorted ascending), which the ordering and the symbolic stage read.
  subroutine hecmw_mf_dist_gmat_build(hecMESH, hecMAT, gmat)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_gmat), intent(inout) :: gmat
    integer(kind=kint), allocatable :: ibuf(:), ilen(:), idisp(:), gc(:)
    integer(kind=kint) :: np, me, comm, nd, n, ng, i, j, k, l, m, r, ptr, ncols, grow, maxcols
    integer(kind=kint) :: nl, nu, g

    np = hecmw_comm_get_size()
    me = hecmw_comm_get_rank()
    comm = hecmw_comm_get_comm()
    nd = hecMAT%NDOF
    n = hecMAT%N
    call hecmw_mf_dist_gmat_finalize(gmat)
    allocate(gmat%nn(np), gmat%ndisp(np+1), gmat%vblk(np), ilen(np), idisp(np))
    call hecmw_allgather_int_1(n, gmat%nn, comm)
    gmat%ndisp(1) = 0
    do r = 1, np
      gmat%ndisp(r+1) = gmat%ndisp(r) + gmat%nn(r)
    enddo
    ng = gmat%ndisp(np+1)
    allocate(gmat%gid(ng))
    call hecmw_allgatherv_int(hecMESH%global_node_ID, n, gmat%gid, gmat%nn, gmat%ndisp, comm)

    ! the int stream: per internal row its column count and the global column ids
    m = n + (hecMAT%indexL(n) - hecMAT%indexL(0)) + (hecMAT%indexU(n) - hecMAT%indexU(0))
    allocate(ibuf(max(m, 1)))
    ptr = 0
    do i = 1, n
      ptr = ptr + 1
      ibuf(ptr) = (hecMAT%indexL(i) - hecMAT%indexL(i-1)) + (hecMAT%indexU(i) - hecMAT%indexU(i-1))
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        ptr = ptr + 1
        ibuf(ptr) = mf_global_node(hecMESH, gmat, me, hecMAT%itemL(k))
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        ptr = ptr + 1
        ibuf(ptr) = mf_global_node(hecMESH, gmat, me, hecMAT%itemU(k))
      enddo
    enddo
    call hecmw_allgather_int_1(m, ilen, comm)
    idisp(1) = 0
    do r = 2, np
      idisp(r) = idisp(r-1) + ilen(r-1)
    enddo
    allocate(gmat%stream(idisp(np) + ilen(np)))
    call hecmw_allgatherv_int(ibuf, m, gmat%stream, ilen, idisp, comm)
    deallocate(ibuf)
    ! per row one diagonal block plus the neighbor blocks, so the block count of a rank
    ! equals its int stream length
    gmat%vblk(1:np) = ilen(1:np)
    gmat%nblk = idisp(np) + ilen(np)
    call mf_gmat_closure(gmat, ng)
    gmat%mat%N = ng
    gmat%mat%NP = ng
    gmat%mat%NDOF = nd
    nullify(gmat%mat%indexL, gmat%mat%indexU, gmat%mat%itemL, gmat%mat%itemU)
    nullify(gmat%mat%D, gmat%mat%AL, gmat%mat%AU)
    nullify(gmat%mat%B, gmat%mat%X, gmat%mat%A, gmat%mat%indexA, gmat%mat%itemA)
    deallocate(ilen, idisp)
    if (me /= 0) return

    ! the full profile of rank 0: count the lower/upper split per global row, then sort
    ! the columns of a row ascending
    allocate(gmat%mat%indexL(0:ng), gmat%mat%indexU(0:ng))
    gmat%mat%indexL(0:ng) = 0
    gmat%mat%indexU(0:ng) = 0
    maxcols = 0
    ptr = 0
    do r = 1, np
      do i = 1, gmat%nn(r)
        grow = gmat%ndisp(r) + i
        ncols = gmat%stream(ptr+1)
        maxcols = max(maxcols, ncols + gmat%xptr(grow) - gmat%xptr(grow-1))
        do j = 1, ncols
          g = gmat%stream(ptr+1+j)
          if (g < grow) then
            gmat%mat%indexL(grow) = gmat%mat%indexL(grow) + 1
          else
            gmat%mat%indexU(grow) = gmat%mat%indexU(grow) + 1
          endif
        enddo
        do k = gmat%xptr(grow-1)+1, gmat%xptr(grow)
          if (gmat%xcol(k) < grow) then
            gmat%mat%indexL(grow) = gmat%mat%indexL(grow) + 1
          else
            gmat%mat%indexU(grow) = gmat%mat%indexU(grow) + 1
          endif
        enddo
        ptr = ptr + 1 + ncols
      enddo
    enddo
    do i = 1, ng
      gmat%mat%indexL(i) = gmat%mat%indexL(i-1) + gmat%mat%indexL(i)
      gmat%mat%indexU(i) = gmat%mat%indexU(i-1) + gmat%mat%indexU(i)
    enddo
    gmat%mat%NPL = gmat%mat%indexL(ng)
    gmat%mat%NPU = gmat%mat%indexU(ng)
    allocate(gmat%mat%itemL(max(gmat%mat%NPL, 1)), gmat%mat%itemU(max(gmat%mat%NPU, 1)))
    allocate(gc(max(maxcols, 1)))
    ptr = 0
    do r = 1, np
      do i = 1, gmat%nn(r)
        grow = gmat%ndisp(r) + i
        ncols = gmat%stream(ptr+1)
        do j = 1, ncols
          gc(j) = gmat%stream(ptr+1+j)
        enddo
        ptr = ptr + 1 + ncols
        do k = gmat%xptr(grow-1)+1, gmat%xptr(grow)
          ncols = ncols + 1
          gc(ncols) = gmat%xcol(k)
        enddo
        ! insertion sort by the global column id (unique within a row)
        do j = 2, ncols
          g = gc(j)
          l = j - 1
          do while (l >= 1)
            if (gc(l) <= g) exit
            gc(l+1) = gc(l)
            l = l - 1
          enddo
          gc(l+1) = g
        enddo
        nl = 0
        nu = 0
        do j = 1, ncols
          if (gc(j) < grow) then
            nl = nl + 1
            gmat%mat%itemL(gmat%mat%indexL(grow-1) + nl) = gc(j)
          else
            nu = nu + 1
            gmat%mat%itemU(gmat%mat%indexU(grow-1) + nu) = gc(j)
          endif
        enddo
      enddo
    enddo
    deallocate(gc)
  end subroutine hecmw_mf_dist_gmat_build

  !> Transposed positions absent from the gathered structure, per global row (xptr/xcol).
  !> An assembly can touch the row of a node another rank owns only one-sidedly (the fill
  !> of the contact elimination does), so the union of the internal rows is not always
  !> structurally symmetric. The profile builders add these positions to keep the factored
  !> structure symmetric; the value exchange never routes a block there, so they stay zero.
  subroutine mf_gmat_closure(gmat, ng)
    implicit none
    type(hecmwST_mf_gmat), intent(inout) :: gmat
    integer(kind=kint), intent(in) :: ng
    integer(kind=kint), allocatable :: rp(:), cols(:), fil(:)
    integer(kind=kint) :: ptr, i, j, k, g, ncols, nmiss

    ! the stream lists the rows in global order, so the CSR copy is two straight walks
    allocate(rp(0:ng), fil(ng))
    rp(0:ng) = 0
    ptr = 0
    do i = 1, ng
      ncols = gmat%stream(ptr+1)
      rp(i) = rp(i-1) + ncols
      ptr = ptr + 1 + ncols
    enddo
    allocate(cols(max(rp(ng), 1)))
    ptr = 0
    do i = 1, ng
      ncols = gmat%stream(ptr+1)
      do j = 1, ncols
        cols(rp(i-1)+j) = gmat%stream(ptr+1+j)
      enddo
      ptr = ptr + 1 + ncols
    enddo
    ! insertion sort within a row (the rows are short)
    do i = 1, ng
      do j = rp(i-1)+2, rp(i)
        g = cols(j)
        k = j - 1
        do while (k >= rp(i-1)+1)
          if (cols(k) <= g) exit
          cols(k+1) = cols(k)
          k = k - 1
        enddo
        cols(k+1) = g
      enddo
    enddo
    allocate(gmat%xptr(0:ng))
    gmat%xptr(0:ng) = 0
    do i = 1, ng
      do j = rp(i-1)+1, rp(i)
        g = cols(j)
        if (.not. found(g, i)) gmat%xptr(g) = gmat%xptr(g) + 1
      enddo
    enddo
    do i = 1, ng
      gmat%xptr(i) = gmat%xptr(i-1) + gmat%xptr(i)
    enddo
    nmiss = gmat%xptr(ng)
    allocate(gmat%xcol(max(nmiss, 1)))
    if (nmiss > 0) then
      ! the outer row ascends, so every xcol row collects its columns ascending
      fil(1:ng) = 0
      do i = 1, ng
        do j = rp(i-1)+1, rp(i)
          g = cols(j)
          if (.not. found(g, i)) then
            fil(g) = fil(g) + 1
            gmat%xcol(gmat%xptr(g-1) + fil(g)) = i
          endif
        enddo
      enddo
    endif
    deallocate(rp, cols, fil)

  contains

    !> row r0 of the sorted structure contains column c0
    logical function found(r0, c0)
      integer(kind=kint), intent(in) :: r0, c0
      integer(kind=kint) :: lo0, hi0, mid0
      found = .false.
      lo0 = rp(r0-1) + 1
      hi0 = rp(r0)
      do while (lo0 <= hi0)
        mid0 = (lo0 + hi0) / 2
        if (cols(mid0) == c0) then
          found = .true.
          return
        else if (cols(mid0) < c0) then
          lo0 = mid0 + 1
        else
          hi0 = mid0 - 1
        endif
      enddo
    end function found

  end subroutine mf_gmat_closure

  !> Global node id of a local node: internal nodes by the rank offset, external nodes
  !> through the (owner rank, local id) pair of node_ID, as the sparse matrix interface
  !> of the external direct solvers converts them.
  function mf_global_node(hecMESH, gmat, me, j) result(g)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_mf_gmat), intent(in) :: gmat
    integer(kind=kint), intent(in) :: me, j
    integer(kind=kint) :: g

    if (j <= hecMESH%nn_internal) then
      g = gmat%ndisp(me+1) + j
    else
      g = gmat%ndisp(hecMESH%node_ID(2*j) + 1) + hecMESH%node_ID(2*j-1)
    endif
  end function mf_global_node

  !> Derive from the mapping the partial profile of this rank and the value routing, and
  !> drop the gathered structure stream (rank 0 also drops the full profile the ordering
  !> has read). A block pair belongs to the supernode of its earlier column and is held by
  !> the rank set of that supernode, a delay-proof superset of the per-entry owners; the
  !> index arrays keep the full length, the retained columns of a row are sorted ascending
  !> like the full profile, and D stays full length with the dropped blocks zeroed. The
  !> send and receive lists both enumerate the blocks of a source in its stream order, so
  !> the value transfer needs no structure traffic.
  subroutine hecmw_mf_dist_gmat_part(sym, map, gmat)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_mf_gmat), intent(inout) :: gmat
    integer(kind=kint), allocatable :: c2s(:), gc(:), bb(:), sptr(:)
    integer(kind=kint) :: np, me, ng, nd2, i, j, l, m, r, g, b, ptr, ncols, grow, maxcols
    integer(kind=kint) :: ka, nl, nu, nkeep, nrecv, pos, r0, nr0, d

    np = map%nprocs
    me = map%myrank
    ng = gmat%ndisp(np+1)
    nd2 = gmat%mat%NDOF * gmat%mat%NDOF
    if (associated(gmat%mat%indexL)) then
      deallocate(gmat%mat%indexL, gmat%mat%indexU, gmat%mat%itemL, gmat%mat%itemU)
    endif
    allocate(c2s(sym%nnode))
    call mf_gmat_c2s(sym, c2s)

    ! count the retained lower/upper split per row and the received blocks per source
    allocate(gmat%rcnt(np), gmat%mat%indexL(0:ng), gmat%mat%indexU(0:ng))
    gmat%rcnt(1:np) = 0
    gmat%mat%indexL(0:ng) = 0
    gmat%mat%indexU(0:ng) = 0
    maxcols = 0
    ptr = 0
    do r = 1, np
      do i = 1, gmat%nn(r)
        grow = gmat%ndisp(r) + i
        ka = sym%invp(grow)
        ncols = gmat%stream(ptr+1)
        maxcols = max(maxcols, ncols + gmat%xptr(grow) - gmat%xptr(grow-1))
        if (mine(ka, ka)) gmat%rcnt(r) = gmat%rcnt(r) + 1
        do j = 1, ncols
          g = gmat%stream(ptr+1+j)
          if (mine(ka, sym%invp(g))) then
            gmat%rcnt(r) = gmat%rcnt(r) + 1
            if (g < grow) then
              gmat%mat%indexL(grow) = gmat%mat%indexL(grow) + 1
            else
              gmat%mat%indexU(grow) = gmat%mat%indexU(grow) + 1
            endif
          endif
        enddo
        ! the closure positions are retained like their mirrored partners but receive no block
        do j = gmat%xptr(grow-1)+1, gmat%xptr(grow)
          g = gmat%xcol(j)
          if (mine(ka, sym%invp(g))) then
            if (g < grow) then
              gmat%mat%indexL(grow) = gmat%mat%indexL(grow) + 1
            else
              gmat%mat%indexU(grow) = gmat%mat%indexU(grow) + 1
            endif
          endif
        enddo
        ptr = ptr + 1 + ncols
      enddo
    enddo
    do i = 1, ng
      gmat%mat%indexL(i) = gmat%mat%indexL(i-1) + gmat%mat%indexL(i)
      gmat%mat%indexU(i) = gmat%mat%indexU(i-1) + gmat%mat%indexU(i)
    enddo
    gmat%mat%NPL = gmat%mat%indexL(ng)
    gmat%mat%NPU = gmat%mat%indexU(ng)
    allocate(gmat%mat%itemL(max(gmat%mat%NPL, 1)), gmat%mat%itemU(max(gmat%mat%NPU, 1)))
    allocate(gmat%mat%D(int(ng, 8)*nd2))
    allocate(gmat%mat%AL(max(int(gmat%mat%NPL, 8)*nd2, 1_8)))
    allocate(gmat%mat%AU(max(int(gmat%mat%NPU, 8)*nd2, 1_8)))
    ! the closure positions are never overwritten by the value exchange, so they must be zero
    gmat%mat%D(:) = 0.0d0
    gmat%mat%AL(:) = 0.0d0
    gmat%mat%AU(:) = 0.0d0

    ! second sweep: sort the retained columns of a row ascending and record the landing
    ! slot of every received block, following the stream order of its source
    nrecv = 0
    do r = 1, np
      nrecv = nrecv + gmat%rcnt(r)
    enddo
    allocate(gmat%rdst(max(nrecv, 1)), gc(max(maxcols, 1)), bb(max(maxcols, 1)))
    ptr = 0
    b = 0
    do r = 1, np
      do i = 1, gmat%nn(r)
        grow = gmat%ndisp(r) + i
        ka = sym%invp(grow)
        ncols = gmat%stream(ptr+1)
        if (mine(ka, ka)) then
          b = b + 1
          gmat%rdst(b) = gmat%mat%NPL + gmat%mat%NPU + grow
        endif
        nkeep = 0
        do j = 1, ncols
          g = gmat%stream(ptr+1+j)
          if (mine(ka, sym%invp(g))) then
            nkeep = nkeep + 1
            b = b + 1
            gc(nkeep) = g
            bb(nkeep) = b
          endif
        enddo
        ! retained closure positions take an item slot but no landing slot (bb = 0)
        do j = gmat%xptr(grow-1)+1, gmat%xptr(grow)
          g = gmat%xcol(j)
          if (mine(ka, sym%invp(g))) then
            nkeep = nkeep + 1
            gc(nkeep) = g
            bb(nkeep) = 0
          endif
        enddo
        ptr = ptr + 1 + ncols
        ! insertion sort by the global column id (unique within a row)
        do j = 2, nkeep
          g = gc(j)
          m = bb(j)
          l = j - 1
          do while (l >= 1)
            if (gc(l) <= g) exit
            gc(l+1) = gc(l)
            bb(l+1) = bb(l)
            l = l - 1
          enddo
          gc(l+1) = g
          bb(l+1) = m
        enddo
        nl = 0
        nu = 0
        do j = 1, nkeep
          if (gc(j) < grow) then
            nl = nl + 1
            pos = gmat%mat%indexL(grow-1) + nl
            gmat%mat%itemL(pos) = gc(j)
            if (bb(j) > 0) gmat%rdst(bb(j)) = pos
          else
            nu = nu + 1
            pos = gmat%mat%indexU(grow-1) + nu
            gmat%mat%itemU(pos) = gc(j)
            if (bb(j) > 0) gmat%rdst(bb(j)) = gmat%mat%NPL + pos
          endif
        enddo
      enddo
    enddo

    ! send side: the blocks of my stream segment by destination rank set, kept in the
    ! stream order within a destination
    allocate(gmat%scnt(np), sptr(np))
    gmat%scnt(1:np) = 0
    ptr = 0
    do r = 1, me
      ptr = ptr + gmat%vblk(r)
    enddo
    b = ptr
    do i = 1, gmat%nn(me+1)
      grow = gmat%ndisp(me+1) + i
      ka = sym%invp(grow)
      ncols = gmat%stream(ptr+1)
      call destrange(ka, ka, r0, nr0)
      gmat%scnt(r0+1:r0+nr0) = gmat%scnt(r0+1:r0+nr0) + 1
      do j = 1, ncols
        call destrange(ka, sym%invp(gmat%stream(ptr+1+j)), r0, nr0)
        gmat%scnt(r0+1:r0+nr0) = gmat%scnt(r0+1:r0+nr0) + 1
      enddo
      ptr = ptr + 1 + ncols
    enddo
    sptr(1) = 0
    do r = 2, np
      sptr(r) = sptr(r-1) + gmat%scnt(r-1)
    enddo
    i = sptr(np) + gmat%scnt(np)
    allocate(gmat%ssel(max(i, 1)))
    ptr = b
    b = 0
    do i = 1, gmat%nn(me+1)
      grow = gmat%ndisp(me+1) + i
      ka = sym%invp(grow)
      ncols = gmat%stream(ptr+1)
      b = b + 1
      call destrange(ka, ka, r0, nr0)
      do d = r0 + 1, r0 + nr0
        sptr(d) = sptr(d) + 1
        gmat%ssel(sptr(d)) = b
      enddo
      do j = 1, ncols
        b = b + 1
        call destrange(ka, sym%invp(gmat%stream(ptr+1+j)), r0, nr0)
        do d = r0 + 1, r0 + nr0
          sptr(d) = sptr(d) + 1
          gmat%ssel(sptr(d)) = b
        enddo
      enddo
      ptr = ptr + 1 + ncols
    enddo
    deallocate(gmat%stream, gmat%xptr, gmat%xcol, c2s, gc, bb, sptr)

  contains

    !> this rank holds the pair at permuted column positions (ka0, kb0)
    logical function mine(ka0, kb0)
      integer(kind=kint), intent(in) :: ka0, kb0
      integer(kind=kint) :: s0
      s0 = c2s(min(ka0, kb0))
      mine = me >= map%rbeg(s0) .and. me < map%rbeg(s0) + map%rcnt(s0)
    end function mine

    !> rank range holding the pair at permuted column positions (ka0, kb0)
    subroutine destrange(ka0, kb0, r1, nr1)
      integer(kind=kint), intent(in) :: ka0, kb0
      integer(kind=kint), intent(out) :: r1, nr1
      integer(kind=kint) :: s0
      s0 = c2s(min(ka0, kb0))
      r1 = map%rbeg(s0)
      nr1 = map%rcnt(s0)
    end subroutine destrange

  end subroutine hecmw_mf_dist_gmat_part

  !> Move the values of the internal rows to their holding ranks through the routing of
  !> hecmw_mf_dist_gmat_part (built once per structure): the local block stream is packed
  !> by destination, exchanged all to all and scattered into the partial arrays. Every
  !> retained slot is written on every call.
  subroutine hecmw_mf_dist_gmat_vals(hecMAT, gmat)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_gmat), intent(inout) :: gmat
    real(kind=kreal), allocatable :: vbuf(:), sbuf(:), rbuf(:)
    integer(kind=kint), allocatable :: scs(:), sdisp(:), rcs(:), rdisp(:)
    integer(kind=kint) :: np, me, comm, nd2, n, i, k, r, t, nsend, nrecv
    integer(kind=8) :: ptr, b, base

    np = hecmw_comm_get_size()
    me = hecmw_comm_get_rank()
    comm = hecmw_comm_get_comm()
    gmat%mat%symmetric = hecMAT%symmetric
    nd2 = hecMAT%NDOF * hecMAT%NDOF
    n = hecMAT%N
    allocate(vbuf(max(int(gmat%vblk(me+1), 8)*nd2, 1_8)))
    ptr = 0
    do i = 1, n
      vbuf(ptr+1:ptr+nd2) = hecMAT%D(int(i-1, 8)*nd2+1:int(i, 8)*nd2)
      ptr = ptr + nd2
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        vbuf(ptr+1:ptr+nd2) = hecMAT%AL(int(k-1, 8)*nd2+1:int(k, 8)*nd2)
        ptr = ptr + nd2
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        vbuf(ptr+1:ptr+nd2) = hecMAT%AU(int(k-1, 8)*nd2+1:int(k, 8)*nd2)
        ptr = ptr + nd2
      enddo
    enddo
    allocate(scs(np), sdisp(np), rcs(np), rdisp(np))
    nsend = 0
    nrecv = 0
    do r = 1, np
      scs(r) = gmat%scnt(r) * nd2
      rcs(r) = gmat%rcnt(r) * nd2
      sdisp(r) = nsend * nd2
      rdisp(r) = nrecv * nd2
      nsend = nsend + gmat%scnt(r)
      nrecv = nrecv + gmat%rcnt(r)
    enddo
    allocate(sbuf(max(int(nsend, 8)*nd2, 1_8)), rbuf(max(int(nrecv, 8)*nd2, 1_8)))
    do i = 1, nsend
      b = int(gmat%ssel(i) - 1, 8)*nd2
      base = int(i-1, 8)*nd2
      sbuf(base+1:base+nd2) = vbuf(b+1:b+nd2)
    enddo
    call hecmw_alltoallv_real(sbuf, scs, sdisp, rbuf, rcs, rdisp, comm)
    do b = 1, nrecv
      t = gmat%rdst(b)
      base = (b-1)*nd2
      if (t <= gmat%mat%NPL) then
        gmat%mat%AL(int(t-1, 8)*nd2+1:int(t, 8)*nd2) = rbuf(base+1:base+nd2)
      else if (t <= gmat%mat%NPL + gmat%mat%NPU) then
        t = t - gmat%mat%NPL
        gmat%mat%AU(int(t-1, 8)*nd2+1:int(t, 8)*nd2) = rbuf(base+1:base+nd2)
      else
        t = t - gmat%mat%NPL - gmat%mat%NPU
        gmat%mat%D(int(t-1, 8)*nd2+1:int(t, 8)*nd2) = rbuf(base+1:base+nd2)
      endif
    enddo
    deallocate(vbuf, sbuf, rbuf, scs, sdisp, rcs, rdisp)
  end subroutine hecmw_mf_dist_gmat_vals

  !> Words this rank holds for the matrix (the values with the integer structure and
  !> routing at two integers per word) and the words the replicated design held per rank,
  !> for the memory log.
  subroutine hecmw_mf_dist_gmat_words(gmat, wpart, wrepl)
    implicit none
    type(hecmwST_mf_gmat), intent(in) :: gmat
    integer(kind=8), intent(out) :: wpart, wrepl
    integer(kind=8) :: ng, ni, nd2

    ng = size(gmat%mat%indexL, kind=8) - 1
    nd2 = int(gmat%mat%NDOF, 8)**2
    wpart = size(gmat%mat%D, kind=8) + size(gmat%mat%AL, kind=8) + size(gmat%mat%AU, kind=8)
    ni = 2*(ng+1) + size(gmat%mat%itemL, kind=8) + size(gmat%mat%itemU, kind=8) &
      + size(gmat%ssel, kind=8) + size(gmat%rdst, kind=8) &
      + size(gmat%scnt, kind=8) + size(gmat%rcnt, kind=8) + 3*size(gmat%nn, kind=8) + 1
    wpart = wpart + (ni + 1)/2
    wrepl = int(gmat%nblk, 8)*nd2 + (int(gmat%nblk, 8) - ng + 2*(ng+1) + int(gmat%nblk, 8) + 1)/2
  end subroutine hecmw_mf_dist_gmat_words

  !> Supernode of every permuted column position, from the own column ranges of sym.
  subroutine mf_gmat_c2s(sym, c2s)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    integer(kind=kint), intent(out) :: c2s(:)
    integer(kind=kint) :: s, k

    do s = 1, sym%nsuper
      do k = sym%sptr(s), sym%sptr(s+1) - 1
        c2s(k) = s
      enddo
    enddo
  end subroutine mf_gmat_c2s

  !> Extract from the replicated matrix the part this rank reads under the distribution:
  !> the block pairs of the own columns of its supernodes. A pair belongs to the supernode
  !> of its earlier column in the elimination order; a single-rank supernode keeps its
  !> pairs on its rank, an upper front on its whole rank set, a superset of the per-entry
  !> owners that stays valid under any delayed growth. The index arrays keep the full
  !> length and the retained blocks keep their original order within a row, so the
  !> factorization reads bitwise the same values; D stays full length with the dropped
  !> blocks zeroed, an O(n) vector-class array indexed in place by the shared assembly.
  subroutine hecmw_mf_dist_gmat_extract(full, sym, map, part)
    implicit none
    type(hecmwST_matrix), intent(in) :: full
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(in) :: map
    type(hecmwST_matrix), intent(inout) :: part
    integer(kind=kint), allocatable :: c2s(:)
    integer(kind=kint) :: me, nd2, ng, i, k, ka, nl, nu
    integer(kind=8) :: bs, bd

    me = map%myrank
    nd2 = full%NDOF * full%NDOF
    ng = full%NP
    part%N = full%N
    part%NP = full%NP
    part%NDOF = full%NDOF
    part%symmetric = full%symmetric
    allocate(c2s(sym%nnode))
    call mf_gmat_c2s(sym, c2s)
    allocate(part%indexL(0:ng), part%indexU(0:ng))
    part%indexL(0) = 0
    part%indexU(0) = 0
    do i = 1, ng
      ka = sym%invp(i)
      nl = 0
      do k = full%indexL(i-1)+1, full%indexL(i)
        if (mine(ka, sym%invp(full%itemL(k)))) nl = nl + 1
      enddo
      nu = 0
      do k = full%indexU(i-1)+1, full%indexU(i)
        if (mine(ka, sym%invp(full%itemU(k)))) nu = nu + 1
      enddo
      part%indexL(i) = part%indexL(i-1) + nl
      part%indexU(i) = part%indexU(i-1) + nu
    enddo
    part%NPL = part%indexL(ng)
    part%NPU = part%indexU(ng)
    allocate(part%itemL(max(part%NPL, 1)), part%itemU(max(part%NPU, 1)))
    allocate(part%D(int(ng, 8)*nd2))
    allocate(part%AL(max(int(part%NPL, 8)*nd2, 1_8)))
    allocate(part%AU(max(int(part%NPU, 8)*nd2, 1_8)))
    nullify(part%B, part%X, part%A, part%indexA, part%itemA)
    nl = 0
    nu = 0
    do i = 1, ng
      ka = sym%invp(i)
      bd = int(i-1, 8)*nd2
      if (mine(ka, ka)) then
        part%D(bd+1:bd+nd2) = full%D(bd+1:bd+nd2)
      else
        part%D(bd+1:bd+nd2) = 0.0d0
      endif
      do k = full%indexL(i-1)+1, full%indexL(i)
        if (.not. mine(ka, sym%invp(full%itemL(k)))) cycle
        nl = nl + 1
        part%itemL(nl) = full%itemL(k)
        bs = int(k-1, 8)*nd2
        bd = int(nl-1, 8)*nd2
        part%AL(bd+1:bd+nd2) = full%AL(bs+1:bs+nd2)
      enddo
      do k = full%indexU(i-1)+1, full%indexU(i)
        if (.not. mine(ka, sym%invp(full%itemU(k)))) cycle
        nu = nu + 1
        part%itemU(nu) = full%itemU(k)
        bs = int(k-1, 8)*nd2
        bd = int(nu-1, 8)*nd2
        part%AU(bd+1:bd+nd2) = full%AU(bs+1:bs+nd2)
      enddo
    enddo
    deallocate(c2s)

  contains

    !> this rank holds the pair at permuted column positions (ka0, kb0)
    logical function mine(ka0, kb0)
      integer(kind=kint), intent(in) :: ka0, kb0
      integer(kind=kint) :: s0
      s0 = c2s(min(ka0, kb0))
      mine = me >= map%rbeg(s0) .and. me < map%rbeg(s0) + map%rcnt(s0)
    end function mine

  end subroutine hecmw_mf_dist_gmat_extract

  !> Gather the internal parts of a distributed nodal vector into the replicated global
  !> vector, in the global node numbering.
  subroutine hecmw_mf_dist_gather_vec(gmat, v, gv)
    implicit none
    type(hecmwST_mf_gmat), intent(in) :: gmat
    real(kind=kreal), intent(in) :: v(:)
    real(kind=kreal), intent(out) :: gv(:)
    integer(kind=kint), allocatable :: vlen(:), vdisp(:)
    integer(kind=kint) :: np, me, comm, nd, r

    np = hecmw_comm_get_size()
    me = hecmw_comm_get_rank()
    comm = hecmw_comm_get_comm()
    nd = gmat%mat%NDOF
    allocate(vlen(np), vdisp(np))
    do r = 1, np
      vlen(r) = gmat%nn(r) * nd
      vdisp(r) = gmat%ndisp(r) * nd
    enddo
    call hecmw_allgatherv_real(v, vlen(me+1), gv, vlen, vdisp, comm)
    deallocate(vlen, vdisp)
  end subroutine hecmw_mf_dist_gather_vec

  subroutine hecmw_mf_dist_gmat_finalize(gmat)
    implicit none
    type(hecmwST_mf_gmat), intent(inout) :: gmat

    gmat%nblk = 0
    if (allocated(gmat%nn)) then
      deallocate(gmat%nn, gmat%ndisp, gmat%vblk)
      if (allocated(gmat%stream)) deallocate(gmat%stream)
      if (allocated(gmat%xptr)) deallocate(gmat%xptr, gmat%xcol)
      if (allocated(gmat%gid)) deallocate(gmat%gid)
      if (allocated(gmat%scnt)) deallocate(gmat%scnt, gmat%ssel, gmat%rcnt, gmat%rdst)
      if (associated(gmat%mat%indexL)) then
        deallocate(gmat%mat%indexL, gmat%mat%indexU, gmat%mat%itemL, gmat%mat%itemU)
      endif
      if (associated(gmat%mat%D)) deallocate(gmat%mat%D, gmat%mat%AL, gmat%mat%AU)
      nullify(gmat%mat%indexL, gmat%mat%indexU, gmat%mat%itemL, gmat%mat%itemU)
      nullify(gmat%mat%D, gmat%mat%AL, gmat%mat%AU)
    endif
  end subroutine hecmw_mf_dist_gmat_finalize

end module hecmw_mf_dist
