!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : node graph + aggregation
!!
!! Builds the node-block strong-coupling graph (theta=0 default: a block nonzero
!! is an edge) and runs the Vanek 3-phase greedy aggregation, followed by the
!! "no small aggregates" forced merge (every aggregate ends with >= min_size
!! nodes).  All processing is node-based: dofs of one node are never split.
!!
!! Output: aggr(inode) = aggregate id in 1..naggr, or 0 for an excluded node
!! (no off-diagonal coupling -> Dirichlet/isolated; no coarse representative).
module hecmw_precond_SAAMG_aggregate
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  implicit none

  private
  public :: hecmwST_saamg_nodegraph
  public :: hecmw_saamg_build_nodegraph
  public :: hecmw_saamg_nodegraph_free
  public :: hecmw_saamg_aggregate
  public :: hecmw_saamg_write_vtk

  !> Node-level adjacency graph (CSR, 1-based; no self loops).
  type hecmwST_saamg_nodegraph
    integer(kind=kint) :: n = 0                 !< number of nodes
    integer(kind=kint), allocatable :: index(:) !< size n+1, neighbor ranges
    integer(kind=kint), allocatable :: item(:)  !< neighbor node ids
  end type hecmwST_saamg_nodegraph

contains

  !> Build the node graph from the block-CSR operator A.
  !! With theta (optional) > 0, apply a strength-of-connection filter: an
  !! off-diagonal node-block (I,J) is an edge only if its Frobenius norm is
  !! strong relative to the diagonal blocks,
  !!     ||A_IJ||_F >= theta * sqrt(||A_II||_F * ||A_JJ||_F).
  !! theta absent or <= 0 reproduces the structural graph (every nonzero block is
  !! an edge), matching the original theta=0 behavior.
  subroutine hecmw_saamg_build_nodegraph(A, g, theta)
    implicit none
    type(hecmwST_saamg_bcsr),      intent(in)  :: A
    type(hecmwST_saamg_nodegraph), intent(out) :: g
    real(kind=kreal), optional,  intent(in)  :: theta
    integer(kind=kint) :: nb, nnode, inode, ii, k, jnode, deg, pos, t, nt
    integer(kind=kint), allocatable :: mark(:), touched(:)
    real(kind=kreal),   allocatable :: dnorm(:), acc(:)
    real(kind=kreal) :: th2
    logical :: filt

    nb = A%nb
    nnode = A%nbrow
    call hecmw_saamg_nodegraph_free(g)
    g%n = nnode
    filt = present(theta)
    if (filt) filt = (theta > 0.0d0)
    th2 = 0.0d0
    if (filt) th2 = theta*theta

    allocate(g%index(nnode+1), mark(nnode), touched(nnode), acc(nnode), dnorm(nnode))
    mark = 0; acc = 0.0d0; dnorm = 0.0d0

    ! diagonal block Frobenius norms (only needed for the strength filter):
    ! sum of squares of the diagonal block's entries, read directly from bval
    if (filt) then
      do inode = 1, nnode
        do k = A%browptr(inode), A%browptr(inode+1)-1
          if (A%bcol(k) == inode) then
            do ii = 1, nb*A%mb
              dnorm(inode) = dnorm(inode) + A%bval((k-1)*nb*A%mb+ii)**2
            end do
            exit
          end if
        end do
        dnorm(inode) = sqrt(dnorm(inode))
      end do
    end if

    ! pass 1: degree count (off-diagonal blocks surviving the strength filter)
    g%index = 0
    do inode = 1, nnode
      call gather_blocks(A, nb, inode, mark, touched, acc, nt)
      deg = 0
      do t = 1, nt
        jnode = touched(t)
        if (.not. filt) then
          deg = deg + 1
        else if (acc(jnode) >= th2 * dnorm(inode) * dnorm(jnode)) then
          deg = deg + 1
        end if
      end do
      g%index(inode+1) = deg
    end do
    do inode = 1, nnode
      g%index(inode+1) = g%index(inode+1) + g%index(inode)
    end do
    allocate(g%item(g%index(nnode+1)))

    ! pass 2: fill the surviving neighbor ids (same gather order as pass 1)
    mark = 0
    do inode = 1, nnode
      call gather_blocks(A, nb, inode, mark, touched, acc, nt)
      pos = g%index(inode)
      do t = 1, nt
        jnode = touched(t)
        if (.not. filt) then
          pos = pos + 1; g%item(pos) = jnode
        else if (acc(jnode) >= th2 * dnorm(inode) * dnorm(jnode)) then
          pos = pos + 1; g%item(pos) = jnode
        end if
      end do
    end do
    deallocate(mark, touched, acc, dnorm)
  end subroutine hecmw_saamg_build_nodegraph

  !> Collect node inode's off-diagonal neighbor blocks into touched(1:nt) (block
  !! order) and store each block's sum of squared entries in acc(jnode) (=
  !! ||A_IJ||_F^2), reading the block storage directly.  Block columns are unique
  !! per block row, so no dedup is needed; mark is stamped for the caller's reuse.
  subroutine gather_blocks(A, nb, inode, mark, touched, acc, nt)
    type(hecmwST_saamg_bcsr), intent(in) :: A
    integer(kind=kint),    intent(in)    :: nb, inode
    integer(kind=kint),    intent(inout) :: mark(:), touched(:)
    real(kind=kreal),      intent(inout) :: acc(:)
    integer(kind=kint),    intent(out)   :: nt
    integer(kind=kint) :: t, jnode, nrow, e, b0
    real(kind=kreal)   :: bn
    ! Only internal node columns form graph edges: nrow = A%nbrow is the number of
    ! (internal) row-nodes; halo columns (jnode > nrow, present when A is the
    ! rectangular distributed operator) are dropped -> uncoupled aggregation.
    ! For a square (sequential) A this skip never triggers.
    nrow = A%nbrow
    nt = 0
    do t = A%browptr(inode), A%browptr(inode+1)-1
      jnode = A%bcol(t)
      if (jnode == inode .or. jnode > nrow) cycle
      b0 = (t-1)*nb*A%mb; bn = 0.0d0
      do e = 1, nb*A%mb
        bn = bn + A%bval(b0+e)*A%bval(b0+e)
      end do
      mark(jnode) = inode; nt = nt + 1; touched(nt) = jnode; acc(jnode) = bn
    end do
  end subroutine gather_blocks

  subroutine hecmw_saamg_nodegraph_free(g)
    implicit none
    type(hecmwST_saamg_nodegraph), intent(inout) :: g
    if (allocated(g%index)) deallocate(g%index)
    if (allocated(g%item))  deallocate(g%item)
    g%n = 0
  end subroutine hecmw_saamg_nodegraph_free

  !> Write a legacy-VTK point cloud colored by aggregate id (ParaView diagnostic).
  !! coord is (x,y,z) interleaved, length 3*nnode; aggr(i) in 0..naggr (0=excluded).
  subroutine hecmw_saamg_write_vtk(fname, coord, nnode, aggr)
    implicit none
    character(len=*),   intent(in) :: fname
    real(kind=kreal),   intent(in) :: coord(:)
    integer(kind=kint), intent(in) :: nnode
    integer(kind=kint), intent(in) :: aggr(:)
    integer(kind=kint) :: iu, i

    open(newunit=iu, file=fname, status='replace', action='write')
    write(iu,'(a)') '# vtk DataFile Version 3.0'
    write(iu,'(a)') 'SA-AMG aggregates'
    write(iu,'(a)') 'ASCII'
    write(iu,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(iu,'(a,i0,a)') 'POINTS ', nnode, ' double'
    do i = 1, nnode
      write(iu,'(es16.8,1x,es16.8,1x,es16.8)') &
           coord(3*i-2), coord(3*i-1), coord(3*i)
    end do
    write(iu,'(a,i0,1x,i0)') 'CELLS ', nnode, 2*nnode
    do i = 1, nnode
      write(iu,'(a,i0)') '1 ', i-1
    end do
    write(iu,'(a,i0)') 'CELL_TYPES ', nnode
    do i = 1, nnode
      write(iu,'(a)') '1'
    end do
    write(iu,'(a,i0)') 'POINT_DATA ', nnode
    write(iu,'(a)') 'SCALARS agg_id int 1'
    write(iu,'(a)') 'LOOKUP_TABLE default'
    do i = 1, nnode
      write(iu,'(i0)') aggr(i)
    end do
    close(iu)
  end subroutine hecmw_saamg_write_vtk

  !> Vanek 3-phase aggregation + min-size forced merge.
  !! aggr(1:n): 1..naggr for aggregated nodes, 0 for excluded nodes.
  subroutine hecmw_saamg_aggregate(g, min_size, max_size, aggr, naggr)
    implicit none
    type(hecmwST_saamg_nodegraph), intent(in)  :: g
    integer(kind=kint),          intent(in)  :: min_size, max_size
    integer(kind=kint), allocatable, intent(out) :: aggr(:)
    integer(kind=kint),          intent(out) :: naggr
    integer(kind=kint), allocatable :: state(:), state1(:), sz(:)
    integer(kind=kint) :: n, i, k

    n = g%n
    allocate(state(n), state1(n), sz(n))
    state = 0

    ! excluded: no neighbors
    do i = 1, n
      if (g%index(i+1) == g%index(i)) state(i) = -1
    end do

    ! --- phase 1: cores of diameter <= 2 ---
    ! Seed at a node whose whole neighborhood is still free, then take the seed
    ! plus all its neighbors.  Phase 1 is intentionally left uncapped: the
    ! resulting moderately-large cores (tens of nodes on dense meshes) coarsen
    ! aggressively and keep the smoothed Galerkin operators sparse (low operator
    ! complexity).  Oversized aggregates are prevented in phase 3, not here.
    naggr = 0
    do i = 1, n
      if (state(i) /= 0) cycle
      if (.not. all_neighbors_free(g, state, i)) cycle
      naggr = naggr + 1
      state(i) = naggr
      do k = g%index(i)+1, g%index(i+1)
        state(g%item(k)) = naggr
      end do
    end do
    state1 = state   ! snapshot of phase-1 membership

    sz = 0
    do i = 1, n
      if (state(i) > 0) sz(state(i)) = sz(state(i)) + 1
    end do

    ! --- phase 2: attach leftovers to strongest phase-1 aggregate ---
    call phase2_attach(g, state, state1, sz, max_size)

    ! --- phase 3: remaining unassigned -> bounded BFS aggregates (<= max_size) ---
    call phase3_components(g, state, sz, naggr, max_size)

    ! --- forced merge of aggregates smaller than min_size ---
    call merge_small(g, state, sz, naggr, min_size)

    ! --- compact aggregate ids; excluded -> 0 ---
    call compact_ids(state, n, naggr)

    allocate(aggr(n))
    do i = 1, n
      if (state(i) < 0) then
        aggr(i) = 0
      else
        aggr(i) = state(i)
      end if
    end do

    deallocate(state, state1, sz)
  end subroutine hecmw_saamg_aggregate

  !-----------------------------------------------------------------------------
  ! helpers (private)
  !-----------------------------------------------------------------------------

  logical function all_neighbors_free(g, state, i) result(ok)
    type(hecmwST_saamg_nodegraph), intent(in) :: g
    integer(kind=kint),          intent(in) :: state(:), i
    integer(kind=kint) :: k
    ok = .true.
    do k = g%index(i)+1, g%index(i+1)
      if (state(g%item(k)) /= 0) then
        ok = .false.; return
      end if
    end do
  end function all_neighbors_free

  !> Attach each still-unassigned node to the neighboring phase-1 aggregate to
  !! which it is most strongly connected (here: most edges), respecting max_size.
  subroutine phase2_attach(g, state, state1, sz, max_size)
    type(hecmwST_saamg_nodegraph), intent(in)    :: g
    integer(kind=kint),          intent(inout) :: state(:), sz(:)
    integer(kind=kint),          intent(in)    :: state1(:), max_size
    integer(kind=kint), allocatable :: str(:), touched(:)
    integer(kind=kint) :: n, i, k, a, nt, t, best, beststr
    n = g%n
    allocate(str(n), touched(n)); str = 0
    do i = 1, n
      if (state(i) /= 0) cycle
      nt = 0
      do k = g%index(i)+1, g%index(i+1)
        a = state1(g%item(k))
        if (a > 0) then
          if (str(a) == 0) then; nt = nt + 1; touched(nt) = a; end if
          str(a) = str(a) + 1
        end if
      end do
      best = 0; beststr = 0
      do t = 1, nt
        a = touched(t)
        if (sz(a) < max_size .and. str(a) > beststr) then
          beststr = str(a); best = a
        end if
      end do
      if (best > 0) then
        state(i) = best; sz(best) = sz(best) + 1
      end if
      do t = 1, nt
        str(touched(t)) = 0
      end do
    end do
    deallocate(str, touched)
  end subroutine phase2_attach

  !> Group remaining unassigned nodes into NEW size-bounded aggregates.
  !! A connected blob of leftovers is chopped into BFS aggregates of at most
  !! max_size nodes each (rather than one aggregate per connected component),
  !! so a large leftover region cannot collapse into a single giant aggregate.
  subroutine phase3_components(g, state, sz, naggr, max_size)
    type(hecmwST_saamg_nodegraph), intent(in)    :: g
    integer(kind=kint),          intent(inout) :: state(:), sz(:), naggr
    integer(kind=kint),          intent(in)    :: max_size
    integer(kind=kint), allocatable :: queue(:)
    integer(kind=kint) :: n, i, head, tail, u, k, v, cnt
    n = g%n
    allocate(queue(n))
    do i = 1, n
      if (state(i) /= 0) cycle
      naggr = naggr + 1
      head = 1; tail = 1; queue(1) = i; state(i) = naggr; cnt = 1
      do while (head <= tail .and. cnt < max_size)
        u = queue(head); head = head + 1
        do k = g%index(u)+1, g%index(u+1)
          v = g%item(k)
          if (state(v) == 0) then
            state(v) = naggr; cnt = cnt + 1
            tail = tail + 1; queue(tail) = v
            if (cnt >= max_size) exit
          end if
        end do
      end do
      sz(naggr) = cnt
    end do
    deallocate(queue)
  end subroutine phase3_components

  !> Forced merge: dissolve every aggregate with size < min_size into the
  !! neighboring aggregate it is most strongly connected to (size cap ignored).
  !! An aggregate with no edge to any other aggregate is left as-is (warned).
  subroutine merge_small(g, state, sz, naggr, min_size)
    type(hecmwST_saamg_nodegraph), intent(in)    :: g
    integer(kind=kint),          intent(inout) :: state(:), sz(:)
    integer(kind=kint),          intent(in)    :: naggr, min_size
    integer(kind=kint), allocatable :: str(:), touched(:), stuck(:)
    integer(kind=kint) :: n, i, k, a, b, small, nt, t, best, beststr
    n = g%n
    allocate(str(naggr), touched(naggr), stuck(naggr))
    str = 0; stuck = 0

    do
      ! pick a non-stuck aggregate below min_size
      small = 0
      do a = 1, naggr
        if (sz(a) > 0 .and. sz(a) < min_size .and. stuck(a) == 0) then
          small = a; exit
        end if
      end do
      if (small == 0) exit

      ! tally connection strength to other aggregates
      nt = 0
      do i = 1, n
        if (state(i) /= small) cycle
        do k = g%index(i)+1, g%index(i+1)
          b = state(g%item(k))
          if (b > 0 .and. b /= small) then
            if (str(b) == 0) then; nt = nt + 1; touched(nt) = b; end if
            str(b) = str(b) + 1
          end if
        end do
      end do
      best = 0; beststr = 0
      do t = 1, nt
        b = touched(t)
        if (str(b) > beststr) then; beststr = str(b); best = b; end if
        str(b) = 0
      end do

      if (best == 0) then
        stuck(small) = 1            ! isolated small cluster: cannot merge
        write(*,'(a,i0,a)') 'hecmw_saamg merge_small: aggregate ', small, &
             ' has no neighbor aggregate; left below min_size'
        cycle
      end if
      ! merge small -> best
      do i = 1, n
        if (state(i) == small) state(i) = best
      end do
      sz(best) = sz(best) + sz(small)
      sz(small) = 0
    end do

    deallocate(str, touched, stuck)
  end subroutine merge_small

  !> Renumber aggregates with sz>0 to 1..naggr_new (in place on state).
  subroutine compact_ids(state, n, naggr)
    integer(kind=kint), intent(inout) :: state(:), naggr
    integer(kind=kint), intent(in)    :: n
    integer(kind=kint), allocatable :: map(:)
    integer(kind=kint) :: i, a, cnt
    allocate(map(naggr)); map = 0
    cnt = 0
    do i = 1, n
      a = state(i)
      if (a > 0) then
        if (map(a) == 0) then; cnt = cnt + 1; map(a) = cnt; end if
        state(i) = map(a)
      end if
    end do
    naggr = cnt
    deallocate(map)
  end subroutine compact_ids

end module hecmw_precond_SAAMG_aggregate
