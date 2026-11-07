!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : Step 3 (node graph + Vanek aggregation)
!
! Builds a 2D 5-point grid (plus one isolated node) and checks the structural
! invariants the later QR / Galerkin steps rely on:
!   * isolated node is excluded (aggr = 0);
!   * every other node belongs to exactly one aggregate (partition);
!   * every aggregate has >= min_size nodes (forced merge);
!   * every aggregate is connected in the node graph;
!   * real coarsening occurs (naggr << nnode).
! Also checks node-graph construction for an nb=3 block matrix.
!-------------------------------------------------------------------------------
program test_step3
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_aggregate
  implicit none

  integer(kind=kint) :: npass, nfail
  npass = 0; nfail = 0

  call test_nodegraph_block(npass, nfail)
  call test_grid_aggregation(npass, nfail)

  write(*,'(/,a,i0,a,i0,a)') '=== Step 3: ', npass, ' passed, ', nfail, ' failed ==='
  if (nfail > 0) stop 1

contains

  subroutine check(name, cond, npass, nfail)
    character(len=*),   intent(in)    :: name
    logical,            intent(in)    :: cond
    integer(kind=kint), intent(inout) :: npass, nfail
    if (cond) then
      npass = npass + 1; write(*,'(a,a)') '  [PASS] ', name
    else
      nfail = nfail + 1; write(*,'(a,a)') '  [FAIL] ', name
    end if
  end subroutine check

  !> 2D Nx*Ny grid, 5-point stencil, nb=1; optionally append one isolated node.
  subroutine build_grid(nx, ny, isolated, A, nnode)
    integer(kind=kint), intent(in)  :: nx, ny
    logical,            intent(in)  :: isolated
    type(hecmwST_saamg_bcsr), intent(out) :: A
    integer(kind=kint), intent(out) :: nnode
    real(kind=kreal), allocatable :: dense(:,:)
    integer(kind=kint) :: ix, iy, id, n
    nnode = nx*ny
    n = nnode; if (isolated) n = n + 1
    allocate(dense(n,n)); dense = 0.0d0
    do iy = 1, ny
      do ix = 1, nx
        id = (iy-1)*nx + ix
        dense(id,id) = 4.0d0
        if (ix > 1)  dense(id, id-1)  = -1.0d0
        if (ix < nx) dense(id, id+1)  = -1.0d0
        if (iy > 1)  dense(id, id-nx) = -1.0d0
        if (iy < ny) dense(id, id+nx) = -1.0d0
      end do
    end do
    if (isolated) dense(n,n) = 1.0d0   ! diagonal only -> excluded
    call hecmw_saamg_bcsr_from_dense(dense, n, 1, A)
    deallocate(dense)
  end subroutine build_grid

  !-----------------------------------------------------------------------------
  subroutine test_nodegraph_block(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    type(hecmwST_saamg_bcsr) :: A
    type(hecmwST_saamg_nodegraph) :: g
    real(kind=kreal) :: dense(6,6)
    ! two nb=3 nodes coupled by one off-diagonal entry
    dense = 0.0d0
    dense(1:3,1:3) = reshape([4.d0,1.d0,0.d0, 1.d0,3.d0,1.d0, 0.d0,1.d0,2.d0],[3,3])
    dense(4:6,4:6) = reshape([5.d0,0.d0,1.d0, 0.d0,4.d0,0.d0, 1.d0,0.d0,3.d0],[3,3])
    dense(1,4) = 0.5d0; dense(4,1) = 0.5d0
    call hecmw_saamg_bcsr_from_dense(dense, 6, 3, A)
    call hecmw_saamg_build_nodegraph(A, g)
    call check('nb=3 collapses to 2 nodes', g%n == 2, npass, nfail)
    call check('node 1 degree 1', g%index(2)-g%index(1) == 1, npass, nfail)
    call check('node 2 degree 1', g%index(3)-g%index(2) == 1, npass, nfail)
    call check('edge 1->2', g%item(g%index(1)+1) == 2, npass, nfail)
    call hecmw_saamg_nodegraph_free(g)
    call hecmw_saamg_bcsr_free(A)
  end subroutine test_nodegraph_block

  !-----------------------------------------------------------------------------
  subroutine test_grid_aggregation(npass, nfail)
    integer(kind=kint), intent(inout) :: npass, nfail
    integer(kind=kint), parameter :: nx = 12, ny = 12, min_size = 3, max_size = 9
    type(hecmwST_saamg_bcsr) :: A
    type(hecmwST_saamg_nodegraph) :: g
    integer(kind=kint), allocatable :: aggr(:), sz(:)
    integer(kind=kint) :: nnode, ntot, naggr, i, ia, nexcl, nassigned
    logical :: ok_partition, ok_minsize, ok_connected, ok_excluded

    call build_grid(nx, ny, .true., A, nnode)
    ntot = A%n                       ! nnode + 1 isolated
    call hecmw_saamg_build_nodegraph(A, g)
    call hecmw_saamg_aggregate(g, min_size, max_size, aggr, naggr)

    ! isolated node (last) excluded
    ok_excluded = (aggr(ntot) == 0)

    ! partition: each grid node in 1..naggr exactly once; count sizes
    allocate(sz(naggr)); sz = 0
    nexcl = 0; nassigned = 0; ok_partition = .true.
    do i = 1, ntot
      ia = aggr(i)
      if (ia == 0) then
        nexcl = nexcl + 1
      else if (ia >= 1 .and. ia <= naggr) then
        sz(ia) = sz(ia) + 1; nassigned = nassigned + 1
      else
        ok_partition = .false.
      end if
    end do
    if (nassigned /= ntot - nexcl) ok_partition = .false.

    ! min-size
    ok_minsize = .true.
    do ia = 1, naggr
      if (sz(ia) < min_size) ok_minsize = .false.
    end do

    ! connectivity of each aggregate
    ok_connected = aggregates_connected(g, aggr, naggr)

    write(*,'(a,i0,a,i0,a,i0,a,f5.2)') '   nodes=', ntot, ' excluded=', nexcl, &
         ' naggr=', naggr, ' coarsening=', real(naggr,kreal)/real(ntot-nexcl,kreal)
    write(*,'(a,i0,a,i0)') '   aggregate size min=', minval(sz), ' max=', maxval(sz)

    call check('isolated node excluded',        ok_excluded,  npass, nfail)
    call check('exactly one isolated excluded',  nexcl == 1,   npass, nfail)
    call check('valid partition / full coverage', ok_partition, npass, nfail)
    call check('all aggregates >= min_size',     ok_minsize,   npass, nfail)
    call check('all aggregates connected',       ok_connected, npass, nfail)
    call check('real coarsening (3x+)',  naggr*3 <= ntot-nexcl, npass, nfail)

    deallocate(sz, aggr)
    call hecmw_saamg_nodegraph_free(g)
    call hecmw_saamg_bcsr_free(A)
  end subroutine test_grid_aggregation

  !> Check every aggregate forms a connected subgraph (BFS within aggregate).
  logical function aggregates_connected(g, aggr, naggr) result(ok)
    type(hecmwST_saamg_nodegraph), intent(in) :: g
    integer(kind=kint),          intent(in) :: aggr(:), naggr
    integer(kind=kint), allocatable :: seen(:), stack(:), start(:), cnt(:)
    integer(kind=kint) :: n, i, a, top, u, k, v, reached
    n = g%n
    allocate(seen(n), stack(n), start(naggr), cnt(naggr))
    seen = 0; start = 0; cnt = 0
    do i = 1, n
      a = aggr(i)
      if (a > 0) then
        if (start(a) == 0) start(a) = i
        cnt(a) = cnt(a) + 1
      end if
    end do
    ok = .true.
    do a = 1, naggr
      if (start(a) == 0) cycle
      reached = 0
      top = 1; stack(1) = start(a); seen(start(a)) = a
      do while (top > 0)
        u = stack(top); top = top - 1
        reached = reached + 1
        do k = g%index(u)+1, g%index(u+1)
          v = g%item(k)
          if (aggr(v) == a .and. seen(v) /= a) then
            seen(v) = a; top = top + 1; stack(top) = v
          end if
        end do
      end do
      if (reached /= cnt(a)) then; ok = .false.; exit; end if
    end do
    deallocate(seen, stack, start, cnt)
  end function aggregates_connected

end program test_step3
