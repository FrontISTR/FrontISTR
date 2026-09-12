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
  public :: hecmw_mf_dist_symbolic_bcast
  public :: hecmwST_mf_gmat
  public :: hecmw_mf_dist_gmat_build
  public :: hecmw_mf_dist_gmat_vals
  public :: hecmw_mf_dist_gmat_finalize
  public :: hecmw_mf_dist_gather_vec

  !> Subtree-to-subcube mapping over the supernodal tree (proportional mapping): every
  !> supernode carries the rank set rbeg/rcnt that processes it; a supernode whose rank set
  !> spans more than one rank is an upper front, whose fully summed part the owner (master)
  !> factors while the contribution row tiles are distributed over the rank set (1D row
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

  !> The global matrix replicated on every rank, assembled from the internal rows of the
  !> distributed matrix (complete in the overlapped HEC-MW assembly). The global node
  !> numbering concatenates the internal nodes rank by rank, the way the sparse matrix
  !> interface of the external direct solvers numbers them. The structure is gathered once;
  !> the values are regathered by hecmw_mf_dist_gmat_vals before every factorization, the
  !> gathered block stream landing in its slots through the dst map built alongside the
  !> structure.
  type hecmwST_mf_gmat
    type(hecmwST_matrix) :: mat
    integer(kind=kint) :: nblk = 0                !< gathered value stream length in nd*nd blocks
    integer(kind=kint), allocatable :: nn(:)      !< internal nodes of rank r at nn(r+1)
    integer(kind=kint), allocatable :: ndisp(:)   !< global node offset of rank r at ndisp(r+1)
    integer(kind=kint), allocatable :: vblk(:)    !< value stream blocks per rank
    integer(kind=kint), allocatable :: dst(:)     !< block slot: 1..NPL in AL, then NPU in AU, then N in D
  end type hecmwST_mf_gmat

contains

  !> Build the mapping from the symbolic structure. The rank set of the virtual root is all
  !> ranks; the set of an upper front is split among its children in ascending child order,
  !> proportionally to the subtree weights (factor word estimates), so the result is
  !> deterministic and identical on every rank. The owner of an upper front is the owner of
  !> its heaviest child subtree, keeping the fan-in communication local to the subcube.
  subroutine hecmw_mf_dist_map_build(sym, map)
    implicit none
    type(hecmwST_mf_symbolic), intent(in) :: sym
    type(hecmwST_mf_map), intent(out) :: map
    integer(kind=8), allocatable :: wt(:)
    integer(kind=kint), allocatable :: cptr(:), clist(:), wptr(:)
    integer(kind=kint) :: ns, s, p, c, i, nroot, nup, cbest
    integer(kind=8) :: w, wbest

    ns = sym%nsuper
    map%nprocs = hecmw_comm_get_size()
    map%myrank = hecmw_comm_get_rank()
    map%comm = hecmw_comm_get_comm()
    allocate(map%owner(ns), map%upper(ns))

    ! subtree weights: factor panel words of the supernode, accumulated bottom-up (the
    ! ascending numbering places children before parents)
    allocate(wt(ns))
    do s = 1, ns
      w = 0
      do i = sym%rptr(s), sym%rptr(s+1) - 1
        w = w + sym%ndof(sym%rlist(i))
      enddo
      wt(s) = w * int(sym_cdofcount(s), 8)
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

    ! rank ranges: the roots share all ranks, an upper front splits its range among its
    ! children; parents carry larger indices, so a descending sweep sets every range
    ! before it is read
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
      if (map%rcnt(s) > 1) call mf_map_split(wt, clist(cptr(s):cptr(s+1)-1), map%rbeg(s), map%rcnt(s), map%rbeg, map%rcnt)
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

  !> Gather the nonzero structure of the internal rows of every rank and build the global
  !> matrix profile with the dst slot map of the value stream. The stream carries, row by
  !> row, the diagonal block followed by the lower then upper neighbor blocks in the local
  !> item order; a block lands in the global lower or upper part by comparing the global
  !> ids, the columns of a row sorted ascending.
  subroutine hecmw_mf_dist_gmat_build(hecMESH, hecMAT, gmat)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_gmat), intent(inout) :: gmat
    integer(kind=kint), allocatable :: ibuf(:), irbuf(:), ilen(:), idisp(:), gc(:), bb(:)
    integer(kind=kint) :: np, me, comm, nd, n, ng, i, j, k, l, m, r, b, ptr, ncols, grow, maxcols
    integer(kind=kint) :: nl, nu, g, pos

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
    allocate(irbuf(idisp(np) + ilen(np)))
    call hecmw_allgatherv_int(ibuf, m, irbuf, ilen, idisp, comm)
    deallocate(ibuf)
    ! per row one diagonal block plus the neighbor blocks, so the block count of a rank
    ! equals its int stream length
    gmat%vblk(1:np) = ilen(1:np)
    gmat%nblk = idisp(np) + ilen(np)

    ! count the lower/upper split per global row
    gmat%mat%N = ng
    gmat%mat%NP = ng
    gmat%mat%NDOF = nd
    allocate(gmat%mat%indexL(0:ng), gmat%mat%indexU(0:ng))
    gmat%mat%indexL(0:ng) = 0
    gmat%mat%indexU(0:ng) = 0
    maxcols = 0
    ptr = 0
    do r = 1, np
      do i = 1, gmat%nn(r)
        grow = gmat%ndisp(r) + i
        ncols = irbuf(ptr+1)
        maxcols = max(maxcols, ncols)
        do j = 1, ncols
          g = irbuf(ptr+1+j)
          if (g < grow) then
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
    allocate(gmat%mat%D(int(ng, 8)*nd*nd))
    allocate(gmat%mat%AL(max(int(gmat%mat%NPL, 8)*nd*nd, 1_8)))
    allocate(gmat%mat%AU(max(int(gmat%mat%NPU, 8)*nd*nd, 1_8)))
    nullify(gmat%mat%B, gmat%mat%X, gmat%mat%A, gmat%mat%indexA, gmat%mat%itemA)

    ! second sweep: sort the columns of a row ascending and record the slot of every block
    allocate(gmat%dst(gmat%nblk), gc(max(maxcols, 1)), bb(max(maxcols, 1)))
    ptr = 0
    b = 0
    do r = 1, np
      do i = 1, gmat%nn(r)
        grow = gmat%ndisp(r) + i
        ncols = irbuf(ptr+1)
        b = b + 1
        gmat%dst(b) = gmat%mat%NPL + gmat%mat%NPU + grow
        do j = 1, ncols
          gc(j) = irbuf(ptr+1+j)
          bb(j) = b + j
        enddo
        b = b + ncols
        ptr = ptr + 1 + ncols
        ! insertion sort by the global column id (unique within a row)
        do j = 2, ncols
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
        do j = 1, ncols
          if (gc(j) < grow) then
            nl = nl + 1
            pos = gmat%mat%indexL(grow-1) + nl
            gmat%mat%itemL(pos) = gc(j)
            gmat%dst(bb(j)) = pos
          else
            nu = nu + 1
            pos = gmat%mat%indexU(grow-1) + nu
            gmat%mat%itemU(pos) = gc(j)
            gmat%dst(bb(j)) = gmat%mat%NPL + pos
          endif
        enddo
      enddo
    enddo
    deallocate(irbuf, ilen, idisp, gc, bb)
  end subroutine hecmw_mf_dist_gmat_build

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

  !> Regather the values of the internal rows into the global matrix (the structure and the
  !> dst map come from hecmw_mf_dist_gmat_build).
  subroutine hecmw_mf_dist_gmat_vals(hecMAT, gmat)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_gmat), intent(inout) :: gmat
    real(kind=kreal), allocatable :: vbuf(:), vrbuf(:)
    integer(kind=kint), allocatable :: vlen(:), vdisp(:)
    integer(kind=kint) :: np, me, comm, nd2, n, i, k, r, t
    integer(kind=8) :: ptr, b, base

    np = hecmw_comm_get_size()
    me = hecmw_comm_get_rank()
    comm = hecmw_comm_get_comm()
    nd2 = hecMAT%NDOF * hecMAT%NDOF
    n = hecMAT%N
    allocate(vlen(np), vdisp(np))
    do r = 1, np
      vlen(r) = gmat%vblk(r) * nd2
    enddo
    vdisp(1) = 0
    do r = 2, np
      vdisp(r) = vdisp(r-1) + vlen(r-1)
    enddo
    allocate(vbuf(max(vlen(me+1), 1)), vrbuf(int(gmat%nblk, 8)*nd2))
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
    call hecmw_allgatherv_real(vbuf, vlen(me+1), vrbuf, vlen, vdisp, comm)
    do b = 1, gmat%nblk
      t = gmat%dst(b)
      base = (b-1)*nd2
      if (t <= gmat%mat%NPL) then
        gmat%mat%AL(int(t-1, 8)*nd2+1:int(t, 8)*nd2) = vrbuf(base+1:base+nd2)
      else if (t <= gmat%mat%NPL + gmat%mat%NPU) then
        t = t - gmat%mat%NPL
        gmat%mat%AU(int(t-1, 8)*nd2+1:int(t, 8)*nd2) = vrbuf(base+1:base+nd2)
      else
        t = t - gmat%mat%NPL - gmat%mat%NPU
        gmat%mat%D(int(t-1, 8)*nd2+1:int(t, 8)*nd2) = vrbuf(base+1:base+nd2)
      endif
    enddo
    deallocate(vbuf, vrbuf, vlen, vdisp)
  end subroutine hecmw_mf_dist_gmat_vals

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
      deallocate(gmat%nn, gmat%ndisp, gmat%vblk, gmat%dst)
      deallocate(gmat%mat%indexL, gmat%mat%indexU, gmat%mat%itemL, gmat%mat%itemU)
      deallocate(gmat%mat%D, gmat%mat%AL, gmat%mat%AU)
    endif
  end subroutine hecmw_mf_dist_gmat_finalize

end module hecmw_mf_dist
