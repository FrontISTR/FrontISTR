!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : lightweight comm table
!!
!! A per-level communication table, decoupled from hecmwST_local_mesh, so that
!! coarse levels (which have no mesh) carry their own halo-exchange descriptor.
!! The halo exchange mirrors hecmw_solver_SR (HECMW_SOLVE_SEND_RECV) but is kept
!! self-contained (raw MPI guarded by HECMW_SERIAL) so the SA-AMG core stays
!! FrontISTR-independent and still compiles in the stand-alone gfortran harness
!! (HECMW_SERIAL -> the exchange is a no-op).
!!
!! Vector convention: a level vector has length nb*nnode; the internal part is
!! 1..nb*nint, the halo part nb*nint+1..nb*nnode is filled by comm_update from
!! the owning ranks.  Reductions/norms must use the internal part only.
module hecmw_precond_SAAMG_comm
  use hecmw_precond_SAAMG_util
#ifndef HECMW_SERIAL
  ! HEC-MW convention: go through the comm wrappers, not raw MPI.  Guarded so the
  ! standalone (-DHECMW_SERIAL) keeps the self-contained no-op stubs and needs no
  ! HEC-MW comm layer.
  use hecmw_util,     only: HECMW_STATUS_SIZE, hecmw_sum, hecmw_max, &
       hecmw_abort, hecmw_comm_get_comm
  use m_hecmw_comm_f, only: hecmw_isend_r, hecmw_irecv_r, hecmw_isend_int, &
       hecmw_irecv_int, hecmw_waitall, hecmw_comm_size, hecmw_allgather_int_1, &
       hecmw_alltoall_int, hecmw_allgatherv_int, hecmw_allgatherv_real, &
       hecmw_alltoallv_int, hecmw_alltoallv_real, hecmw_allreduce_I_comm, &
       hecmw_allreduce_R_comm
#ifdef _OPENACC
  ! GPU halo exchange: explicit device buffers + host_data use_device (CUDA-aware MPI),
  ! mirroring hecmw_solver_SR.F90.  acc_malloc/acc_map_data/acc_free + c_devptr.
  use openacc
#endif
#endif
  implicit none

  private
  public :: hecmwST_saamg_comm
  public :: hecmw_saamg_comm_update
  public :: hecmw_saamg_comm_reverse_add
  public :: hecmw_saamg_comm_update_i
  public :: hecmw_saamg_comm_update_var
  public :: hecmw_saamg_comm_update_var_i
  public :: hecmw_saamg_comm_allgather_int
  public :: hecmw_saamg_comm_exchange_neighbor_int
  public :: hecmw_saamg_comm_allreduce_max_r
  public :: hecmw_saamg_comm_allreduce_sum_r
  public :: hecmw_saamg_comm_allreduce_max_int
  public :: hecmw_saamg_comm_allreduce_sum_int
  public :: hecmw_saamg_comm_allgatherv_triplets
  public :: hecmw_saamg_comm_allgatherv_real
  public :: hecmw_saamg_comm_alltoall_int
  public :: hecmw_saamg_comm_alltoallv_int
  public :: hecmw_saamg_comm_alltoallv_triplets
  public :: hecmw_saamg_comm_alltoallv_real
  public :: hecmw_saamg_comm_size
  public :: hecmw_saamg_comm_copy
  public :: hecmw_saamg_comm_free
  public :: hecmw_saamg_abort
  public :: hecmw_saamg_check_alloc
  public :: hecmw_saamg_lapack_available
  public :: hecmw_saamg_require_lapack
  public :: hecmw_saamg_comm_init_serial

  !> Lightweight communication table (mirrors the SEND_RECV argument layout).
  type hecmwST_saamg_comm
    integer(kind=kint) :: comm    = 0   !< MPI communicator (passed in; never MPI_COMM_WORLD directly)
    integer(kind=kint) :: my_rank = 0
    integer(kind=kint) :: nint    = 0   !< # internal nodes at this level
    integer(kind=kint) :: nnode   = 0   !< # total nodes (internal + halo)
    integer(kind=kint) :: nb      = 0   !< dofs per node at this level (block size)
    integer(kind=kint) :: n_neighbor = 0
    integer(kind=kint), allocatable :: neighbor(:)       !< neighbor ranks, size n_neighbor
    integer(kind=kint), allocatable :: import_index(:)   !< 0:n_neighbor, halo nodes to receive
    integer(kind=kint), allocatable :: import_item(:)    !< local halo node ids (nint+1..nnode)
    integer(kind=kint), allocatable :: export_index(:)   !< 0:n_neighbor, internal nodes to send
    integer(kind=kint), allocatable :: export_item(:)    !< local internal node ids (1..nint)
    integer(kind=kint), allocatable :: gnode(:)          !< global node id per local node (partition-invariant seeds)
  end type hecmwST_saamg_comm

contains

  !> Fatal-error termination with a clear message.  Under MPI this is a COLLECTIVE
  !! abort (hecmw_abort -> MPI_Abort) so a setup failure on one rank tears the job
  !! down cleanly instead of a bare `stop` that leaves the other ranks hanging in
  !! the next collective.  In the stand-alone (HECMW_SERIAL) harness it is a plain
  !! stop.  Used by the SA-AMG core in place of bare `stop` on unrecoverable errors.
  subroutine hecmw_saamg_abort(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') '#### SA-AMG fatal error: '//trim(msg)
#ifndef HECMW_SERIAL
    call hecmw_abort(hecmw_comm_get_comm())
#else
    stop 1
#endif
  end subroutine hecmw_saamg_abort

  !> Whether this build links LAPACK (SA-AMG's setup needs it: dense coarsest
  !! factorization, per-aggregate QR, Lanczos tridiagonal eigenvalue).  The
  !! HECMW_WITH_LAPACK macro is set by the build (--with-lapack / -DWITH_LAPACK).
  logical function hecmw_saamg_lapack_available()
#ifdef HECMW_WITH_LAPACK
    hecmw_saamg_lapack_available = .true.
#else
    hecmw_saamg_lapack_available = .false.
#endif
  end function hecmw_saamg_lapack_available

  !> Abort with a uniform message when a LAPACK-only code path is reached in a
  !! build without LAPACK.  Defensive: the backend checks availability up front
  !! and aborts before setup, so this normally never fires.
  subroutine hecmw_saamg_require_lapack(site)
    character(len=*), intent(in) :: site
    call hecmw_saamg_abort('SA-AMG (PRECOND=22) requires a LAPACK-enabled build ['//trim(site)// &
         ']; rebuild with --with-lapack (setup.sh) or -DWITH_LAPACK=ON (cmake), or select another preconditioner')
  end subroutine hecmw_saamg_require_lapack

  !> Report a failed allocation (stat /= 0) with a clear message and a collective
  !! abort, instead of a bare runtime stop.  Used on the size-scaling allocations
  !! so an out-of-memory / size-overflow failure is diagnosable.  (Under default
  !! Linux overcommit a true OOM arrives as a SIGKILL before allocate returns, so
  !! this catches the catchable cases: size overflow, address-space exhaustion,
  !! and cgroup/ulimit-bound failures.)
  subroutine hecmw_saamg_check_alloc(ier, what)
    integer(kind=kint), intent(in) :: ier
    character(len=*),   intent(in) :: what
    if (ier /= 0) call hecmw_saamg_abort('memory allocation failed: '//trim(what))
  end subroutine hecmw_saamg_check_alloc

  !> Build a trivial single-rank communication table for `nnode` nodes of block
  !! size `nb`: no neighbors, internal == all nodes, global id = local id.  This
  !! mirrors what comm_from_mesh produces on one rank (n_neighbor==0 leaves the
  !! neighbor/import/export arrays unallocated -- the exchange routines short-circuit
  !! on n_neighbor==0).  Lets the sequential wrapper drive the distributed core.
  subroutine hecmw_saamg_comm_init_serial(cmt, nnode, nb)
    type(hecmwST_saamg_comm), intent(out) :: cmt
    integer(kind=kint),     intent(in)  :: nnode, nb
    integer(kind=kint) :: i
    call hecmw_saamg_comm_free(cmt)
#ifndef HECMW_SERIAL
    cmt%comm = hecmw_comm_get_comm()   ! valid single-process communicator for global reductions
#else
    cmt%comm = 0
#endif
    cmt%my_rank    = 0
    cmt%nint       = nnode
    cmt%nnode      = nnode
    cmt%nb         = nb
    cmt%n_neighbor = 0
    allocate(cmt%gnode(nnode))
    do i = 1, nnode
      cmt%gnode(i) = i
    end do
  end subroutine hecmw_saamg_comm_init_serial

  !> Halo exchange of a block vector X (length nb*nnode): fill the halo region
  !! nb*nint+1..nb*nnode from the owning ranks.  No-op when serial / no neighbors.
  subroutine hecmw_saamg_comm_update(cmt, nb, X)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(in)    :: nb
    real(kind=kreal),       intent(inout) :: X(:)
#ifndef HECMW_SERIAL
    real(kind=kreal),   allocatable :: WS(:), WR(:)
    integer(kind=kint), allocatable :: req1(:), req2(:)
    integer(kind=kint), allocatable :: sta1(:,:), sta2(:,:)
    integer(kind=kint) :: neib, istart, inum, k, kk, ii, nreq1, nreq2, ns, nr
#ifdef _OPENACC
    type(c_devptr) :: WS_dev, WR_dev
#endif

    if (cmt%n_neighbor == 0) return
    ns = cmt%export_index(cmt%n_neighbor)
    nr = cmt%import_index(cmt%n_neighbor)
    allocate(WS(nb*ns), WR(nb*nr))
    allocate(req1(cmt%n_neighbor), req2(cmt%n_neighbor))
    allocate(sta1(HECMW_STATUS_SIZE, cmt%n_neighbor), sta2(HECMW_STATUS_SIZE, cmt%n_neighbor))
#ifdef _OPENACC
    ! Map WS/WR to explicit device buffers so pack/unpack run on device and
    ! host_data use_device passes a true device pointer to CUDA-aware MPI (SR.F90 idiom).
    WS_dev = acc_malloc(kreal * nb * ns)
    WR_dev = acc_malloc(kreal * nb * nr)
    call acc_map_data(WS, WS_dev, kreal * nb * ns)
    call acc_map_data(WR, WR_dev, kreal * nb * nr)
#endif

    ! pack + post sends
    nreq1 = 0
    do neib = 1, cmt%n_neighbor
      istart = cmt%export_index(neib-1)
      inum   = cmt%export_index(neib) - istart
      if (inum == 0) cycle
#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent collapse(2)
      do k = istart+1, istart+inum
        do kk = 1, nb
          ii = cmt%export_item(k)
          WS(nb*(k-1)+kk) = X(nb*(ii-1)+kk)
        end do
      end do
      !$acc end kernels
      nreq1 = nreq1 + 1
      !$acc host_data use_device(WS)
      call hecmw_isend_r(WS(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 0, cmt%comm, req1(nreq1))
      !$acc end host_data
#else
      do k = istart+1, istart+inum
        ii = cmt%export_item(k)
        do kk = 1, nb
          WS(nb*(k-1)+kk) = X(nb*(ii-1)+kk)
        end do
      end do
      nreq1 = nreq1 + 1
      call hecmw_isend_r(WS(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 0, cmt%comm, req1(nreq1))
#endif
    end do

    ! post receives
    nreq2 = 0
    do neib = 1, cmt%n_neighbor
      istart = cmt%import_index(neib-1)
      inum   = cmt%import_index(neib) - istart
      if (inum == 0) cycle
      nreq2 = nreq2 + 1
#ifdef _OPENACC
      !$acc host_data use_device(WR)
      call hecmw_irecv_r(WR(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 0, cmt%comm, req2(nreq2))
      !$acc end host_data
#else
      call hecmw_irecv_r(WR(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 0, cmt%comm, req2(nreq2))
#endif
    end do

    call hecmw_waitall(nreq2, req2, sta2)

    ! unpack into halo region
    do neib = 1, cmt%n_neighbor
      istart = cmt%import_index(neib-1)
      inum   = cmt%import_index(neib) - istart
#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent collapse(2)
      do k = istart+1, istart+inum
        do kk = 1, nb
          ii = cmt%import_item(k)
          X(nb*(ii-1)+kk) = WR(nb*(k-1)+kk)
        end do
      end do
      !$acc end kernels
#else
      do k = istart+1, istart+inum
        ii = cmt%import_item(k)
        do kk = 1, nb
          X(nb*(ii-1)+kk) = WR(nb*(k-1)+kk)
        end do
      end do
#endif
    end do

    call hecmw_waitall(nreq1, req1, sta1)
#ifdef _OPENACC
    call acc_unmap_data(WS)
    call acc_unmap_data(WR)
    call acc_free(WS_dev)
    call acc_free(WR_dev)
#endif
    deallocate(WS, WR, req1, req2, sta1, sta2)
#else
    ! serial: no halo, nothing to exchange
    if (nb < 0) X(1) = X(1)   ! silence unused-arg warnings without effect
#endif
  end subroutine hecmw_saamg_comm_update

  !> Reverse halo exchange with accumulation: send each halo (import) node's value
  !! to its owner and ADD it into the owner's (export) node.  The transpose of
  !! comm_update; used to assemble restriction contributions (P^T res) that land on
  !! coarse rows owned by neighbors.  After this, halo entries of X are stale (only
  !! the internal/owned entries are meaningful).  No-op when serial / no neighbors.
  subroutine hecmw_saamg_comm_reverse_add(cmt, nb, X)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(in)    :: nb
    real(kind=kreal),       intent(inout) :: X(:)
#ifndef HECMW_SERIAL
    real(kind=kreal),   allocatable :: WS(:), WR(:)
    integer(kind=kint), allocatable :: req1(:), req2(:), sta1(:,:), sta2(:,:)
    integer(kind=kint) :: neib, istart, inum, k, kk, ii, nreq1, nreq2, ns, nr
#ifdef _OPENACC
    type(c_devptr) :: WS_dev, WR_dev
#endif

    if (cmt%n_neighbor == 0) return
    ns = cmt%import_index(cmt%n_neighbor)   ! sending halo (import) values
    nr = cmt%export_index(cmt%n_neighbor)   ! receiving into owned (export) nodes
    allocate(WS(nb*ns), WR(nb*nr))
    allocate(req1(cmt%n_neighbor), req2(cmt%n_neighbor))
    allocate(sta1(HECMW_STATUS_SIZE, cmt%n_neighbor), sta2(HECMW_STATUS_SIZE, cmt%n_neighbor))
#ifdef _OPENACC
    WS_dev = acc_malloc(kreal * nb * ns)
    WR_dev = acc_malloc(kreal * nb * nr)
    call acc_map_data(WS, WS_dev, kreal * nb * ns)
    call acc_map_data(WR, WR_dev, kreal * nb * nr)
#endif

    nreq1 = 0
    do neib = 1, cmt%n_neighbor
      istart = cmt%import_index(neib-1)
      inum   = cmt%import_index(neib) - istart
      if (inum == 0) cycle
#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent collapse(2)
      do k = istart+1, istart+inum
        do kk = 1, nb
          ii = cmt%import_item(k)
          WS(nb*(k-1)+kk) = X(nb*(ii-1)+kk)
        end do
      end do
      !$acc end kernels
      nreq1 = nreq1 + 1
      !$acc host_data use_device(WS)
      call hecmw_isend_r(WS(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 2, cmt%comm, req1(nreq1))
      !$acc end host_data
#else
      do k = istart+1, istart+inum
        ii = cmt%import_item(k)
        do kk = 1, nb
          WS(nb*(k-1)+kk) = X(nb*(ii-1)+kk)
        end do
      end do
      nreq1 = nreq1 + 1
      call hecmw_isend_r(WS(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 2, cmt%comm, req1(nreq1))
#endif
    end do
    nreq2 = 0
    do neib = 1, cmt%n_neighbor
      istart = cmt%export_index(neib-1)
      inum   = cmt%export_index(neib) - istart
      if (inum == 0) cycle
      nreq2 = nreq2 + 1
#ifdef _OPENACC
      !$acc host_data use_device(WR)
      call hecmw_irecv_r(WR(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 2, cmt%comm, req2(nreq2))
      !$acc end host_data
#else
      call hecmw_irecv_r(WR(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 2, cmt%comm, req2(nreq2))
#endif
    end do
    call hecmw_waitall(nreq2, req2, sta2)
    ! unpack: accumulate halo contributions into owned nodes.  Per-neighbor export_item
    ! values are distinct, so collapse(2) within one neighbor has no += race; different
    ! neighbors accumulate in separate (sequential) kernels.
    do neib = 1, cmt%n_neighbor
      istart = cmt%export_index(neib-1)
      inum   = cmt%export_index(neib) - istart
#ifdef _OPENACC
      !$acc kernels
      !$acc loop independent collapse(2)
      do k = istart+1, istart+inum
        do kk = 1, nb
          ii = cmt%export_item(k)
          X(nb*(ii-1)+kk) = X(nb*(ii-1)+kk) + WR(nb*(k-1)+kk)
        end do
      end do
      !$acc end kernels
#else
      do k = istart+1, istart+inum
        ii = cmt%export_item(k)
        do kk = 1, nb
          X(nb*(ii-1)+kk) = X(nb*(ii-1)+kk) + WR(nb*(k-1)+kk)
        end do
      end do
#endif
    end do
    call hecmw_waitall(nreq1, req1, sta1)
#ifdef _OPENACC
    call acc_unmap_data(WS)
    call acc_unmap_data(WR)
    call acc_free(WS_dev)
    call acc_free(WR_dev)
#endif
    deallocate(WS, WR, req1, req2, sta1, sta2)
#else
    if (nb < 0) X(1) = X(1)
#endif
  end subroutine hecmw_saamg_comm_reverse_add

  !> Integer halo exchange of a block vector IX (length nb*cmt%nnode): fill the
  !! halo region from owners.  Used to propagate coarse-node global ids and
  !! aggregate ids across rank boundaries.  No-op when serial / no neighbors.
  subroutine hecmw_saamg_comm_update_i(cmt, nb, IX)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(in)    :: nb
    integer(kind=kint),     intent(inout) :: IX(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: WS(:), WR(:)
    integer(kind=kint), allocatable :: req1(:), req2(:)
    integer(kind=kint), allocatable :: sta1(:,:), sta2(:,:)
    integer(kind=kint) :: neib, istart, inum, k, kk, ii, nreq1, nreq2, ns, nr

    if (cmt%n_neighbor == 0) return
    ns = cmt%export_index(cmt%n_neighbor)
    nr = cmt%import_index(cmt%n_neighbor)
    allocate(WS(nb*ns), WR(nb*nr))
    allocate(req1(cmt%n_neighbor), req2(cmt%n_neighbor))
    allocate(sta1(HECMW_STATUS_SIZE, cmt%n_neighbor), sta2(HECMW_STATUS_SIZE, cmt%n_neighbor))

    nreq1 = 0
    do neib = 1, cmt%n_neighbor
      istart = cmt%export_index(neib-1)
      inum   = cmt%export_index(neib) - istart
      if (inum == 0) cycle
      do k = istart+1, istart+inum
        ii = cmt%export_item(k)
        do kk = 1, nb
          WS(nb*(k-1)+kk) = IX(nb*(ii-1)+kk)
        end do
      end do
      nreq1 = nreq1 + 1
      call hecmw_isend_int(WS(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 0, cmt%comm, req1(nreq1))
    end do
    nreq2 = 0
    do neib = 1, cmt%n_neighbor
      istart = cmt%import_index(neib-1)
      inum   = cmt%import_index(neib) - istart
      if (inum == 0) cycle
      nreq2 = nreq2 + 1
      call hecmw_irecv_int(WR(nb*istart+1:nb*istart+nb*inum), nb*inum, &
           cmt%neighbor(neib), 0, cmt%comm, req2(nreq2))
    end do
    call hecmw_waitall(nreq2, req2, sta2)
    do neib = 1, cmt%n_neighbor
      istart = cmt%import_index(neib-1)
      inum   = cmt%import_index(neib) - istart
      do k = istart+1, istart+inum
        ii = cmt%import_item(k)
        do kk = 1, nb
          IX(nb*(ii-1)+kk) = WR(nb*(k-1)+kk)
        end do
      end do
    end do
    call hecmw_waitall(nreq1, req1, sta1)
    deallocate(WS, WR, req1, req2, sta1, sta2)
#else
    if (nb < 0) IX(1) = IX(1)
#endif
  end subroutine hecmw_saamg_comm_update_i

  !> Variable-length real halo exchange.  Node i carries cnt(i) blocks of pb
  !! reals each, packed contiguously in X starting at block-offset off(i)
  !! (0-based); X has length off(nnode+1)*pb.  Internal nodes' data and cnt
  !! are filled by the caller; halo (import) nodes' cnt must already hold the
  !! owner's count (exchange it first with comm_update_i, payload 1), and this
  !! routine fills their X region.  Unlike comm_update there is no per-node
  !! padding to a global maximum: the messages carry exactly sum(cnt) blocks.
  subroutine hecmw_saamg_comm_update_var(cmt, pb, cnt, off, X)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(in)    :: pb
    integer(kind=kint),     intent(in)    :: cnt(:)
    integer(kind=kint),     intent(in)    :: off(:)
    real(kind=kreal),       intent(inout) :: X(:)
#ifndef HECMW_SERIAL
    real(kind=kreal),   allocatable :: WS(:), WR(:)
    integer(kind=kint), allocatable :: req1(:), req2(:)
    integer(kind=kint), allocatable :: sta1(:,:), sta2(:,:)
    integer(kind=kint), allocatable :: sdisp(:), rdisp(:)
    integer(kind=kint) :: neib, k, ii, c, o, wp, nreq1, nreq2, b, nblk, astat

    if (cmt%n_neighbor == 0) return
    allocate(sdisp(0:cmt%n_neighbor), rdisp(0:cmt%n_neighbor))
    sdisp(0) = 0; rdisp(0) = 0
    do neib = 1, cmt%n_neighbor
      nblk = 0
      do k = cmt%export_index(neib-1)+1, cmt%export_index(neib)
        nblk = nblk + cnt(cmt%export_item(k))
      end do
      sdisp(neib) = sdisp(neib-1) + nblk
      nblk = 0
      do k = cmt%import_index(neib-1)+1, cmt%import_index(neib)
        nblk = nblk + cnt(cmt%import_item(k))
      end do
      rdisp(neib) = rdisp(neib-1) + nblk
    end do
    allocate(WS(sdisp(cmt%n_neighbor)*pb), WR(rdisp(cmt%n_neighbor)*pb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'comm_update_var (WS/WR)')
    allocate(req1(cmt%n_neighbor), req2(cmt%n_neighbor))
    allocate(sta1(HECMW_STATUS_SIZE, cmt%n_neighbor), sta2(HECMW_STATUS_SIZE, cmt%n_neighbor))

    ! pack + post sends
    nreq1 = 0
    do neib = 1, cmt%n_neighbor
      if (sdisp(neib) == sdisp(neib-1)) cycle
      wp = sdisp(neib-1)
      do k = cmt%export_index(neib-1)+1, cmt%export_index(neib)
        ii = cmt%export_item(k); c = cnt(ii); o = off(ii)
        do b = 1, c*pb
          WS(wp*pb+b) = X(o*pb+b)
        end do
        wp = wp + c
      end do
      nreq1 = nreq1 + 1
      call hecmw_isend_r(WS(sdisp(neib-1)*pb+1:sdisp(neib)*pb), &
           (sdisp(neib)-sdisp(neib-1))*pb, cmt%neighbor(neib), 0, cmt%comm, req1(nreq1))
    end do

    ! post receives
    nreq2 = 0
    do neib = 1, cmt%n_neighbor
      if (rdisp(neib) == rdisp(neib-1)) cycle
      nreq2 = nreq2 + 1
      call hecmw_irecv_r(WR(rdisp(neib-1)*pb+1:rdisp(neib)*pb), &
           (rdisp(neib)-rdisp(neib-1))*pb, cmt%neighbor(neib), 0, cmt%comm, req2(nreq2))
    end do

    call hecmw_waitall(nreq2, req2, sta2)

    ! unpack into halo region
    do neib = 1, cmt%n_neighbor
      wp = rdisp(neib-1)
      do k = cmt%import_index(neib-1)+1, cmt%import_index(neib)
        ii = cmt%import_item(k); c = cnt(ii); o = off(ii)
        do b = 1, c*pb
          X(o*pb+b) = WR(wp*pb+b)
        end do
        wp = wp + c
      end do
    end do

    call hecmw_waitall(nreq1, req1, sta1)
    deallocate(WS, WR, req1, req2, sta1, sta2, sdisp, rdisp)
#else
    if (pb < 0) X(1) = X(1)   ! silence unused-arg warnings without effect
#endif
  end subroutine hecmw_saamg_comm_update_var

  !> Variable-length integer halo exchange (1 int per block).  Companion to
  !! comm_update_var for the per-block column ids: node i carries cnt(i) ints
  !! at IX block-offset off(i).
  subroutine hecmw_saamg_comm_update_var_i(cmt, cnt, off, IX)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(in)    :: cnt(:)
    integer(kind=kint),     intent(in)    :: off(:)
    integer(kind=kint),     intent(inout) :: IX(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: WS(:), WR(:)
    integer(kind=kint), allocatable :: req1(:), req2(:)
    integer(kind=kint), allocatable :: sta1(:,:), sta2(:,:)
    integer(kind=kint), allocatable :: sdisp(:), rdisp(:)
    integer(kind=kint) :: neib, k, ii, c, o, wp, nreq1, nreq2, b, nblk, astat

    if (cmt%n_neighbor == 0) return
    allocate(sdisp(0:cmt%n_neighbor), rdisp(0:cmt%n_neighbor))
    sdisp(0) = 0; rdisp(0) = 0
    do neib = 1, cmt%n_neighbor
      nblk = 0
      do k = cmt%export_index(neib-1)+1, cmt%export_index(neib)
        nblk = nblk + cnt(cmt%export_item(k))
      end do
      sdisp(neib) = sdisp(neib-1) + nblk
      nblk = 0
      do k = cmt%import_index(neib-1)+1, cmt%import_index(neib)
        nblk = nblk + cnt(cmt%import_item(k))
      end do
      rdisp(neib) = rdisp(neib-1) + nblk
    end do
    allocate(WS(sdisp(cmt%n_neighbor)), WR(rdisp(cmt%n_neighbor)), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'comm_update_var_i (WS/WR)')
    allocate(req1(cmt%n_neighbor), req2(cmt%n_neighbor))
    allocate(sta1(HECMW_STATUS_SIZE, cmt%n_neighbor), sta2(HECMW_STATUS_SIZE, cmt%n_neighbor))

    nreq1 = 0
    do neib = 1, cmt%n_neighbor
      if (sdisp(neib) == sdisp(neib-1)) cycle
      wp = sdisp(neib-1)
      do k = cmt%export_index(neib-1)+1, cmt%export_index(neib)
        ii = cmt%export_item(k); c = cnt(ii); o = off(ii)
        do b = 1, c
          WS(wp+b) = IX(o+b)
        end do
        wp = wp + c
      end do
      nreq1 = nreq1 + 1
      call hecmw_isend_int(WS(sdisp(neib-1)+1:sdisp(neib)), &
           sdisp(neib)-sdisp(neib-1), cmt%neighbor(neib), 0, cmt%comm, req1(nreq1))
    end do

    nreq2 = 0
    do neib = 1, cmt%n_neighbor
      if (rdisp(neib) == rdisp(neib-1)) cycle
      nreq2 = nreq2 + 1
      call hecmw_irecv_int(WR(rdisp(neib-1)+1:rdisp(neib)), &
           rdisp(neib)-rdisp(neib-1), cmt%neighbor(neib), 0, cmt%comm, req2(nreq2))
    end do

    call hecmw_waitall(nreq2, req2, sta2)

    do neib = 1, cmt%n_neighbor
      wp = rdisp(neib-1)
      do k = cmt%import_index(neib-1)+1, cmt%import_index(neib)
        ii = cmt%import_item(k); c = cnt(ii); o = off(ii)
        do b = 1, c
          IX(o+b) = WR(wp+b)
        end do
        wp = wp + c
      end do
    end do

    call hecmw_waitall(nreq1, req1, sta1)
    deallocate(WS, WR, req1, req2, sta1, sta2, sdisp, rdisp)
#else
    if (size(cnt) < 0) IX(1) = IX(1)   ! silence unused-arg warnings without effect
#endif
  end subroutine hecmw_saamg_comm_update_var_i

  !> Number of ranks in this level's communicator (1 when serial).
  function hecmw_saamg_comm_size(cmt) result(np)
    type(hecmwST_saamg_comm), intent(in) :: cmt
    integer(kind=kint) :: np
#ifndef HECMW_SERIAL
    call hecmw_comm_size(cmt%comm, np)
#else
    np = 1
    if (cmt%my_rank < 0) np = 1
#endif
  end function hecmw_saamg_comm_size

  !> Allgather one integer from every rank into rbuf (size must be #ranks).
  subroutine hecmw_saamg_comm_allgather_int(cmt, sval, rbuf)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: sval
    integer(kind=kint),     intent(out) :: rbuf(:)
#ifndef HECMW_SERIAL
    call hecmw_allgather_int_1(sval, rbuf, cmt%comm)
#else
    rbuf(1) = sval
#endif
  end subroutine hecmw_saamg_comm_allgather_int

  !> Exchange one integer with each neighbor: sendvals(k) is sent to neighbor(k),
  !! recvvals(k) is received from neighbor(k).  Used to check comm-table symmetry.
  subroutine hecmw_saamg_comm_exchange_neighbor_int(cmt, sendvals, recvvals)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: sendvals(:)
    integer(kind=kint),     intent(out) :: recvvals(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: req(:), sta(:,:)
    integer(kind=kint) :: neib, nreq
    if (cmt%n_neighbor == 0) return
    allocate(req(2*cmt%n_neighbor), sta(HECMW_STATUS_SIZE, 2*cmt%n_neighbor))
    nreq = 0
    do neib = 1, cmt%n_neighbor
      nreq = nreq + 1
      call hecmw_irecv_int(recvvals(neib:neib), 1, cmt%neighbor(neib), 1, cmt%comm, req(nreq))
    end do
    do neib = 1, cmt%n_neighbor
      nreq = nreq + 1
      call hecmw_isend_int(sendvals(neib:neib), 1, cmt%neighbor(neib), 1, cmt%comm, req(nreq))
    end do
    call hecmw_waitall(nreq, req, sta)
    deallocate(req, sta)
#else
    if (size(sendvals) > 0) recvvals = sendvals
#endif
  end subroutine hecmw_saamg_comm_exchange_neighbor_int

  !> Global max-reduction of a scalar over the communicator (no-op when serial).
  subroutine hecmw_saamg_comm_allreduce_max_r(cmt, s)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    real(kind=kreal),       intent(inout) :: s
#ifndef HECMW_SERIAL
    real(kind=kreal) :: buf(1)
    buf(1) = s; call hecmw_allreduce_R_comm(buf, 1, hecmw_max, cmt%comm); s = buf(1)
#else
    if (cmt%my_rank < 0) s = s
#endif
  end subroutine hecmw_saamg_comm_allreduce_max_r

  !> Global max-reduction of an integer over the communicator (no-op when serial).
  subroutine hecmw_saamg_comm_allreduce_max_int(cmt, n)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(inout) :: n
#ifndef HECMW_SERIAL
    integer(kind=kint) :: buf(1)
    buf(1) = n; call hecmw_allreduce_I_comm(buf, 1, hecmw_max, cmt%comm); n = buf(1)
#else
    if (cmt%my_rank < 0) n = n
#endif
  end subroutine hecmw_saamg_comm_allreduce_max_int

  !> Global sum-reduction of an integer over the communicator (no-op when serial).
  subroutine hecmw_saamg_comm_allreduce_sum_int(cmt, n)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    integer(kind=kint),     intent(inout) :: n
#ifndef HECMW_SERIAL
    integer(kind=kint) :: buf(1)
    buf(1) = n; call hecmw_allreduce_I_comm(buf, 1, hecmw_sum, cmt%comm); n = buf(1)
#else
    if (cmt%my_rank < 0) n = n
#endif
  end subroutine hecmw_saamg_comm_allreduce_sum_int

  !> Global sum-reduction of a scalar over the communicator (no-op when serial).
  subroutine hecmw_saamg_comm_allreduce_sum_r(cmt, s)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    real(kind=kreal),       intent(inout) :: s
#ifndef HECMW_SERIAL
    real(kind=kreal) :: buf(1)
    buf(1) = s; call hecmw_allreduce_R_comm(buf, 1, hecmw_sum, cmt%comm); s = buf(1)
#else
    if (cmt%my_rank < 0) s = s
#endif
  end subroutine hecmw_saamg_comm_allreduce_sum_r

  !> Allgatherv a triplet stream (ti,tj,tv)[1:nloc] from every rank into the
  !! globally-concatenated (gti,gtj,gtv)[1:ntot] (allocated here).  Used to gather
  !! every rank's partial coarse-operator contributions for redundant assembly.
  subroutine hecmw_saamg_comm_allgatherv_triplets(cmt, nloc, ti, tj, tv, ntot, gti, gtj, gtv)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: nloc
    integer(kind=kint),     intent(in)  :: ti(:), tj(:)
    real(kind=kreal),       intent(in)  :: tv(:)
    integer(kind=kint),              intent(out) :: ntot
    integer(kind=kint), allocatable, intent(out) :: gti(:), gtj(:)
    real(kind=kreal),   allocatable, intent(out) :: gtv(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: counts(:), displs(:)
    integer(kind=kint) :: nprocs, p
    call hecmw_comm_size(cmt%comm, nprocs)
    allocate(counts(nprocs), displs(nprocs))
    call hecmw_allgather_int_1(nloc, counts, cmt%comm)
    displs(1) = 0
    do p = 2, nprocs
      displs(p) = displs(p-1) + counts(p-1)
    end do
    ntot = displs(nprocs) + counts(nprocs)
    allocate(gti(max(ntot,1)), gtj(max(ntot,1)), gtv(max(ntot,1)))
    call hecmw_allgatherv_int(ti, nloc, gti, counts, displs, cmt%comm)
    call hecmw_allgatherv_int(tj, nloc, gtj, counts, displs, cmt%comm)
    call hecmw_allgatherv_real(tv, nloc, gtv, counts, displs, cmt%comm)
    deallocate(counts, displs)
#else
    ntot = nloc
    allocate(gti(max(ntot,1)), gtj(max(ntot,1)), gtv(max(ntot,1)))
    gti(1:nloc) = ti(1:nloc); gtj(1:nloc) = tj(1:nloc); gtv(1:nloc) = tv(1:nloc)
#endif
  end subroutine hecmw_saamg_comm_allgatherv_triplets

  !> Allgatherv a real vector vloc[1:nloc] from every rank into gv[1:ntot]
  !! (allocated here).  Used to gather the distributed coarse RHS / near-kernel.
  subroutine hecmw_saamg_comm_allgatherv_real(cmt, nloc, vloc, ntot, gv)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: nloc
    real(kind=kreal),       intent(in)  :: vloc(:)
    integer(kind=kint),              intent(out) :: ntot
    real(kind=kreal),   allocatable, intent(out) :: gv(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: counts(:), displs(:)
    integer(kind=kint) :: nprocs, p
    call hecmw_comm_size(cmt%comm, nprocs)
    allocate(counts(nprocs), displs(nprocs))
    call hecmw_allgather_int_1(nloc, counts, cmt%comm)
    displs(1) = 0
    do p = 2, nprocs
      displs(p) = displs(p-1) + counts(p-1)
    end do
    ntot = displs(nprocs) + counts(nprocs)
    allocate(gv(max(ntot,1)))
    call hecmw_allgatherv_real(vloc, nloc, gv, counts, displs, cmt%comm)
    deallocate(counts, displs)
#else
    ntot = nloc
    allocate(gv(max(ntot,1)))
    gv(1:nloc) = vloc(1:nloc)
#endif
  end subroutine hecmw_saamg_comm_allgatherv_real

  !> MPI_Alltoall of one integer per rank: scnt(r) (count this rank sends to rank r-1)
  !! -> rcnt(r) (count this rank receives from rank r-1).  Both size = #ranks.
  subroutine hecmw_saamg_comm_alltoall_int(cmt, scnt, rcnt)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: scnt(:)
    integer(kind=kint),     intent(out) :: rcnt(:)
#ifndef HECMW_SERIAL
    call hecmw_alltoall_int(scnt, 1, rcnt, 1, cmt%comm)
#else
    rcnt(1) = scnt(1)
#endif
  end subroutine hecmw_saamg_comm_alltoall_int

  !> MPI_Alltoallv of integers.  scnt(#ranks) = send counts per destination rank;
  !! sbuf must already be ordered by destination rank.  Returns ntot received and
  !! the concatenated rbuf (allocated here, ordered by source rank).
  subroutine hecmw_saamg_comm_alltoallv_int(cmt, scnt, sbuf, ntot, rbuf)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: scnt(:), sbuf(:)
    integer(kind=kint),              intent(out) :: ntot
    integer(kind=kint), allocatable, intent(out) :: rbuf(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: rcnt(:), sdis(:), rdis(:)
    integer(kind=kint) :: nprocs, p
    call hecmw_comm_size(cmt%comm, nprocs)
    allocate(rcnt(nprocs), sdis(nprocs), rdis(nprocs))
    call hecmw_alltoall_int(scnt, 1, rcnt, 1, cmt%comm)
    sdis(1) = 0; rdis(1) = 0
    do p = 2, nprocs
      sdis(p) = sdis(p-1) + scnt(p-1); rdis(p) = rdis(p-1) + rcnt(p-1)
    end do
    ntot = rdis(nprocs) + rcnt(nprocs)
    allocate(rbuf(max(ntot,1)))
    call hecmw_alltoallv_int(sbuf, scnt, sdis, rbuf, rcnt, rdis, cmt%comm)
    deallocate(rcnt, sdis, rdis)
#else
    ntot = scnt(1); allocate(rbuf(max(ntot,1)))
    if (ntot > 0) rbuf(1:ntot) = sbuf(1:ntot)
#endif
  end subroutine hecmw_saamg_comm_alltoallv_int

  !> MPI_Alltoallv of triplets (2 int + 1 real).  scnt(#ranks) = send counts per
  !! destination; si/sj/sv ordered by destination rank.  Returns ntot received and
  !! ri/rj/rv (allocated here).  Routes sparse-matrix partial contributions to the
  !! owning ranks (possibly non-neighbors) for distributed Galerkin assembly.
  subroutine hecmw_saamg_comm_alltoallv_triplets(cmt, scnt, si, sj, sv, ntot, ri, rj, rv)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: scnt(:), si(:), sj(:)
    real(kind=kreal),       intent(in)  :: sv(:)
    integer(kind=kint),              intent(out) :: ntot
    integer(kind=kint), allocatable, intent(out) :: ri(:), rj(:)
    real(kind=kreal),   allocatable, intent(out) :: rv(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: rcnt(:), sdis(:), rdis(:)
    integer(kind=kint) :: nprocs, p
    call hecmw_comm_size(cmt%comm, nprocs)
    allocate(rcnt(nprocs), sdis(nprocs), rdis(nprocs))
    call hecmw_alltoall_int(scnt, 1, rcnt, 1, cmt%comm)
    sdis(1) = 0; rdis(1) = 0
    do p = 2, nprocs
      sdis(p) = sdis(p-1) + scnt(p-1); rdis(p) = rdis(p-1) + rcnt(p-1)
    end do
    ntot = rdis(nprocs) + rcnt(nprocs)
    allocate(ri(max(ntot,1)), rj(max(ntot,1)), rv(max(ntot,1)))
    call hecmw_alltoallv_int(si, scnt, sdis, ri, rcnt, rdis, cmt%comm)
    call hecmw_alltoallv_int(sj, scnt, sdis, rj, rcnt, rdis, cmt%comm)
    call hecmw_alltoallv_real(sv, scnt, sdis, rv, rcnt, rdis, cmt%comm)
    deallocate(rcnt, sdis, rdis)
#else
    ntot = scnt(1); allocate(ri(max(ntot,1)), rj(max(ntot,1)), rv(max(ntot,1)))
    if (ntot > 0) then; ri(1:ntot) = si(1:ntot); rj(1:ntot) = sj(1:ntot); rv(1:ntot) = sv(1:ntot); end if
#endif
  end subroutine hecmw_saamg_comm_alltoallv_triplets

  !> MPI_Alltoallv of reals only, reusing a fixed routing (same scnt as a prior
  !! _alltoallv_triplets call).  For numeric-only Galerkin refresh: the triplet
  !! POSITIONS are unchanged, so only the values are re-routed.  sval ordered by
  !! destination rank; returns ntot and rval in the identical receive order.
  subroutine hecmw_saamg_comm_alltoallv_real(cmt, scnt, sval, ntot, rval)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: scnt(:)
    real(kind=kreal),       intent(in)  :: sval(:)
    integer(kind=kint),              intent(out) :: ntot
    real(kind=kreal),   allocatable, intent(out) :: rval(:)
#ifndef HECMW_SERIAL
    integer(kind=kint), allocatable :: rcnt(:), sdis(:), rdis(:)
    integer(kind=kint) :: nprocs, p
    call hecmw_comm_size(cmt%comm, nprocs)
    allocate(rcnt(nprocs), sdis(nprocs), rdis(nprocs))
    call hecmw_alltoall_int(scnt, 1, rcnt, 1, cmt%comm)
    sdis(1) = 0; rdis(1) = 0
    do p = 2, nprocs
      sdis(p) = sdis(p-1) + scnt(p-1); rdis(p) = rdis(p-1) + rcnt(p-1)
    end do
    ntot = rdis(nprocs) + rcnt(nprocs)
    allocate(rval(max(ntot,1)))
    call hecmw_alltoallv_real(sval, scnt, sdis, rval, rcnt, rdis, cmt%comm)
    deallocate(rcnt, sdis, rdis)
#else
    ntot = scnt(1); allocate(rval(max(ntot,1)))
    if (ntot > 0) rval(1:ntot) = sval(1:ntot)
#endif
  end subroutine hecmw_saamg_comm_alltoallv_real

  !> Deep copy a communication table:  dst = src.
  subroutine hecmw_saamg_comm_copy(src, dst)
    type(hecmwST_saamg_comm), intent(in)  :: src
    type(hecmwST_saamg_comm), intent(out) :: dst
    integer(kind=kint) :: nnb
    call hecmw_saamg_comm_free(dst)
    dst%comm = src%comm; dst%my_rank = src%my_rank
    dst%nint = src%nint; dst%nnode = src%nnode; dst%nb = src%nb
    dst%n_neighbor = src%n_neighbor
    if (allocated(src%gnode)) then
      allocate(dst%gnode(size(src%gnode))); dst%gnode = src%gnode
    end if
    nnb = src%n_neighbor
    if (nnb == 0) return
    allocate(dst%neighbor(nnb), dst%import_index(0:nnb), dst%export_index(0:nnb))
    dst%neighbor     = src%neighbor
    dst%import_index = src%import_index
    dst%export_index = src%export_index
    allocate(dst%import_item(size(src%import_item)), dst%export_item(size(src%export_item)))
    dst%import_item = src%import_item
    dst%export_item = src%export_item
  end subroutine hecmw_saamg_comm_copy

  subroutine hecmw_saamg_comm_free(cmt)
    type(hecmwST_saamg_comm), intent(inout) :: cmt
    if (allocated(cmt%neighbor))     deallocate(cmt%neighbor)
    if (allocated(cmt%import_index)) deallocate(cmt%import_index)
    if (allocated(cmt%import_item))  deallocate(cmt%import_item)
    if (allocated(cmt%export_index)) deallocate(cmt%export_index)
    if (allocated(cmt%export_item))  deallocate(cmt%export_item)
    if (allocated(cmt%gnode))        deallocate(cmt%gnode)
    cmt%n_neighbor = 0; cmt%nint = 0; cmt%nnode = 0
  end subroutine hecmw_saamg_comm_free

end module hecmw_precond_SAAMG_comm
