!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : distributed coarsening
!!
!! MPI-distributed construction of one coarse level.  Aggregation is "uncoupled":
!! each rank aggregates only its internal nodes (halo neighbors are dropped by the
!! node-graph builder), so an aggregate never crosses a rank boundary and is owned
!! entirely by one rank.  Coarse nodes are globally numbered via allgather offsets;
!! the coarse communication table is derived directly from the fine import/export
!! lists plus the exchanged global aggregate ids -- no extra communication.
!!
!! F4b-1 scope: aggregation + coarse global numbering + coarse comm table.  The
!! distributed prolongator (F4b-2) and Galerkin operator (F4c) build on the maps
!! produced here (aggr, agg_loc, the coarse comm table).
module hecmw_precond_SAAMG_core
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_comm
  use hecmw_precond_SAAMG_matrix,    only: hecmwST_saamg_bcsr, hecmw_saamg_bcsr_free, &
       hecmw_saamg_bcsr_from_triplets, &
       hecmw_saamg_matvec, hecmw_saamg_matvec_d, &
       hecmw_saamg_spgemm, hecmw_saamg_bcsr_transpose, hecmw_saamg_bcsr_copy, &
       hecmw_saamg_bcsr_move, hecmw_saamg_triple_blk, hecmw_saamg_transpose_blk_rows
  use hecmw_precond_SAAMG_smoother,  only: hecmwST_saamg_blockdiag, hecmw_saamg_blockdiag_setup, &
       hecmw_saamg_blockdiag_apply, hecmw_saamg_blockdiag_matvec, hecmw_saamg_blockdiag_free, &
       hecmwST_saamg_chebyshev, hecmw_saamg_cheb_setup, hecmw_saamg_tridiag_max_eig
  use hecmw_precond_SAAMG_aggregate, only: hecmwST_saamg_nodegraph, hecmw_saamg_build_nodegraph, &
       hecmw_saamg_aggregate, hecmw_saamg_nodegraph_free
  use hecmw_precond_SAAMG_prolongation, only: hecmw_saamg_tentative, hecmw_saamg_smooth_prolongator
  use hecmw_precond_SAAMG_coarse,    only: hecmwST_saamg_coarse, hecmw_saamg_coarse_setup, &
       hecmw_saamg_coarse_solve, hecmw_saamg_coarse_free
  use hecmw_precond_SAAMG_coarse_mumps, only: hecmwST_saamg_cmumps, hecmw_saamg_cmumps_available, &
       hecmw_saamg_cmumps_setup, hecmw_saamg_cmumps_refresh, hecmw_saamg_cmumps_solve, &
       hecmw_saamg_cmumps_free
  use hecmw_precond_SAAMG_param,     only: hecmwST_saamg_params
  implicit none

  private
  public :: hecmw_saamg_coarsen_struct
  public :: hecmw_saamg_verify_commtable
  public :: hecmw_saamg_tentative_ext
  public :: hecmw_saamg_lambda_max
  public :: hecmw_saamg_prolongator
  public :: hecmw_saamg_verify_prows
  public :: hecmw_saamg_galerkin_global
  public :: hecmw_saamg_galerkin
  public :: hecmw_saamg_verify_galerkin
  public :: hecmw_saamg_nc_global
  public :: hecmwST_saamg_hier
  public :: hecmw_saamg_setup
  public :: hecmw_saamg_refresh
  public :: hecmw_saamg_apply
  public :: hecmw_saamg_free

  abstract interface
    !> Optional fast fine-level matvec hook: y = A x for the finest operator.  The
    !! FrontISTR backend binds this to HEC-MW's native blocked matvec on hecMAT (so
    !! the hot smoother avoids the scalar-CSR matvec); standalone leaves it null and
    !! the scalar matvec_d is used.  x carries internal+halo (halo refreshed inside);
    !! y returns internal rows (1:n).
    subroutine hecmw_saamg_matvec_if(x, y)
      import :: kreal
      real(kind=kreal), intent(inout) :: x(:)
      real(kind=kreal), intent(out)   :: y(:)
    end subroutine hecmw_saamg_matvec_if
  end interface

  !> A fully-distributed multilevel SA-AMG hierarchy (S2/S3): every level keeps a
  !! distributed rectangular operator + comm table + Chebyshev smoother, and each
  !! coarse operator is assembled distributed via row-owner routing (hecmw_saamg_galerkin).
  !! Only the tiny COARSEST level is solved globally -- by distributed MUMPS (no
  !! gather) or, for the dense fallback, gathered redundantly and factored on every
  !! rank.
  !> Cached routing for numeric-only Galerkin refresh (S4).  The coarse operator's
  !! SPARSITY is fixed across a Newton refresh (aggregation / P-hat / A pattern all
  !! reused), so only the values move.  This records the structural maps so a refresh
  !! skips the operator-comm-table handshake, column localization and triplet sort:
  !!   scnt(nprocs)  send counts per destination rank (same as the build)
  !!   sperm(nt)     triplet t -> its slot in the send buffer (bucket-by-owner order)
  !!   acc_map(nrecv) received triplet k -> its slot in A_next%bval (accumulate)
  type hecmwST_saamg_gplan
    logical :: ready = .false.
    integer(kind=kint) :: nt = 0, nrecv = 0
    integer(kind=kint), allocatable :: scnt(:), sperm(:), acc_map(:)
  end type hecmwST_saamg_gplan

  !> One distributed level of the hierarchy.  The operator A is rectangular
  !! (nint_l*nb_l rows x nnode_l*nb_l cols) with comm table cmt; P maps this level
  !! to the next (nint_l*nb_l x ncnode_l*m).  Coarse levels have nb_l = m.
  type hecmwST_saamg_dlevel
    type(hecmwST_saamg_bcsr)      :: A       !< this level's operator (rect)
    type(hecmwST_saamg_comm)      :: cmt     !< operator comm table (used for matvec/restrict/prolong)
    type(hecmwST_saamg_blockdiag) :: D
    type(hecmwST_saamg_chebyshev) :: cheb
    type(hecmwST_saamg_bcsr)      :: P       !< prolongator to the next level
    type(hecmwST_saamg_bcsr)      :: Pt      !< P^T
    type(hecmwST_saamg_bcsr)      :: phat_ext!< tentative prolongator (kept for refresh)
    type(hecmwST_saamg_comm)      :: cmt_c   !< 1-ring coarse table from coarsen_struct (refresh)
    type(hecmwST_saamg_gplan)     :: gplan   !< cached Galerkin routing for numeric refresh (S4)
    integer(kind=kint) :: nb_l = 0, naggr_local = 0, ncnode = 0, my_off = 0
    integer(kind=kint), allocatable :: chalo_gid(:)  !< global id per P halo coarse node
    real(kind=kreal),   allocatable :: x(:), rhs(:), res(:), tmp(:), rc(:)  !< V-cycle work
    real(kind=kreal),   allocatable :: cheb_r(:), cheb_p(:), cheb_dr(:), cheb_ap(:)  !< Chebyshev scratch (hoisted, allocated once)
    logical :: is_coarsest = .false.
  end type hecmwST_saamg_dlevel

  !> The fully-distributed multilevel hierarchy: lev(1)=finest .. lev(nlevel).  The
  !! coarsest operator is gathered redundantly (small) and factored by dense LDL^T.
  type hecmwST_saamg_hier
    integer(kind=kint) :: nlevel = 0, m = 0, nc_coarsest = 0
    type(hecmwST_saamg_params) :: prm
    type(hecmwST_saamg_dlevel), allocatable :: lev(:)
    type(hecmwST_saamg_coarse) :: coarse                !< coarsest dense LDL^T (fallback)
    type(hecmwST_saamg_cmumps) :: cmumps                !< coarsest distributed MUMPS (F7)
    logical :: coarse_is_mumps = .false.              !< which coarsest backend is in use
    logical :: coarse_is_smoother = .false.           !< coarsest = Chebyshev smoother sweeps (no direct solve)
    logical :: symmetric = .true.                     !< F8: non-sym => SYM=0/LU coarsest
    real(kind=kreal), allocatable :: rcg(:), zcg(:)   !< coarsest gather work (dense path)
    integer(kind=kint), allocatable :: aggr_fine(:)   !< finest aggregation, global id per internal node (dump_vtk)
    procedure(hecmw_saamg_matvec_if), pointer, nopass :: fine_matvec => null()  !< fast finest matvec (HEC-MW), else scalar
    logical :: hook_disabled = .false.  !< temporarily bypass fine_matvec (refresh self-check)
  end type hecmwST_saamg_hier


contains

  !> Uncoupled aggregation + coarse communication table for one distributed level.
  !! Inputs: A = rectangular distributed operator (nint*nb rows x nnode*nb cols),
  !! cmt = this level's comm table, m = near-kernel size (coarse block size).
  !! Outputs:
  !!   aggr(1:nint)     local aggregate id per internal fine node (0 = excluded)
  !!   naggr_local      # aggregates (coarse nodes) owned by this rank
  !!   agg_loc(1:nnode) localized coarse-node id per fine node: internal fine ->
  !!                    1..naggr_local, halo fine -> naggr_local+1..ncnode (the
  !!                    coarse halo nodes owned by neighbors), 0 = excluded
  !!   ncnode           total coarse nodes referenced locally (internal + halo)
  !!   cmt_c            coarse communication table (nint=naggr_local, nnode=ncnode)
  subroutine hecmw_saamg_coarsen_struct(A, cmt, min_size, max_size, theta, m, &
       aggr, naggr_local, agg_loc, ncnode, cmt_c, my_coarse_off, chalo_gid)
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: min_size, max_size, m
    real(kind=kreal),       intent(in)  :: theta
    integer(kind=kint), allocatable, intent(out) :: aggr(:)
    integer(kind=kint),              intent(out) :: naggr_local
    integer(kind=kint), allocatable, intent(out) :: agg_loc(:)
    integer(kind=kint),              intent(out) :: ncnode
    type(hecmwST_saamg_comm),          intent(out) :: cmt_c
    !> global coarse-node offset of this rank, and the global coarse id of each
    !! halo coarse node (local id naggr_local+t -> chalo_gid(t)).  Optional: the
    !! F4b checks do not need them; the F4c Galerkin does (local<->global cols).
    integer(kind=kint), optional,              intent(out) :: my_coarse_off
    integer(kind=kint), optional, allocatable, intent(out) :: chalo_gid(:)

    type(hecmwST_saamg_nodegraph) :: g
    integer(kind=kint), allocatable :: agg_global(:), counts(:), displs(:), halo_gid(:), cmark(:)
    integer(kind=kint) :: nb, nint, nnode, nprocs, my_off, i, r, gg, t, found, nhc
    integer(kind=kint) :: neib, nnb, h, e, c, ks, ke, cnt, pos

    nb    = A%nb
    nint  = A%n / nb
    nnode = A%ncol / nb

    ! 1. uncoupled aggregation (build_nodegraph drops halo neighbors automatically)
    call hecmw_saamg_build_nodegraph(A, g, theta)
    call hecmw_saamg_aggregate(g, min_size, max_size, aggr, naggr_local)
    call hecmw_saamg_nodegraph_free(g)

    ! 2. coarse global numbering: offset = sum of naggr_local over lower ranks
    nprocs = hecmw_saamg_comm_size(cmt)
    allocate(counts(nprocs), displs(nprocs+1))
    call hecmw_saamg_comm_allgather_int(cmt, naggr_local, counts)
    displs(1) = 0
    do r = 1, nprocs
      displs(r+1) = displs(r) + counts(r)
    end do
    my_off = displs(cmt%my_rank+1)

    ! 3. global coarse id per fine node; halo filled from owners
    allocate(agg_global(nnode)); agg_global = 0
    do i = 1, nint
      if (aggr(i) > 0) agg_global(i) = my_off + aggr(i)
    end do
    call hecmw_saamg_comm_update_i(cmt, 1, agg_global)

    ! 4. localize: internal coarse -> 1..naggr_local, remote coarse -> naggr_local+1..ncnode
    allocate(agg_loc(nnode)); agg_loc = 0
    allocate(halo_gid(nnode)); nhc = 0
    do i = 1, nint
      agg_loc(i) = aggr(i)
    end do
    do i = nint+1, nnode
      gg = agg_global(i)
      if (gg == 0) cycle
      found = 0
      do t = 1, nhc
        if (halo_gid(t) == gg) then; found = t; exit; end if
      end do
      if (found == 0) then
        nhc = nhc + 1; halo_gid(nhc) = gg; found = nhc
      end if
      agg_loc(i) = naggr_local + found
    end do
    ncnode = naggr_local + nhc

    if (present(my_coarse_off)) my_coarse_off = my_off
    if (present(chalo_gid)) then
      allocate(chalo_gid(max(nhc,1)))
      if (nhc > 0) chalo_gid(1:nhc) = halo_gid(1:nhc)
    end if

    ! 5. coarse comm table (derived from fine import/export + agg maps)
    nnb = cmt%n_neighbor
    call hecmw_saamg_comm_free(cmt_c)
    cmt_c%comm    = cmt%comm
    cmt_c%my_rank = cmt%my_rank
    cmt_c%nint    = naggr_local
    cmt_c%nnode   = ncnode
    cmt_c%nb      = m
    cmt_c%n_neighbor = nnb
    if (nnb == 0) then
      deallocate(agg_global, counts, displs, halo_gid)
      return
    end if
    allocate(cmt_c%neighbor(nnb))
    allocate(cmt_c%import_index(0:nnb), cmt_c%export_index(0:nnb))
    cmt_c%neighbor(1:nnb) = cmt%neighbor(1:nnb)
    cmt_c%import_index(0) = 0; cmt_c%export_index(0) = 0
    allocate(cmark(ncnode)); cmark = 0

    ! coarse import: distinct halo coarse ids per fine-import neighbor (count, fill)
    do neib = 1, nnb
      ks = cmt%import_index(neib-1)+1; ke = cmt%import_index(neib)
      cnt = 0
      do h = ks, ke
        c = agg_loc(cmt%import_item(h))
        if (c > 0 .and. cmark(c) /= neib) then; cmark(c) = neib; cnt = cnt + 1; end if
      end do
      cmt_c%import_index(neib) = cmt_c%import_index(neib-1) + cnt
    end do
    allocate(cmt_c%import_item(cmt_c%import_index(nnb)))
    cmark = 0
    do neib = 1, nnb
      ks = cmt%import_index(neib-1)+1; ke = cmt%import_index(neib)
      pos = cmt_c%import_index(neib-1)
      do h = ks, ke
        c = agg_loc(cmt%import_item(h))
        if (c > 0 .and. cmark(c) /= neib) then
          cmark(c) = neib; pos = pos + 1; cmt_c%import_item(pos) = c
        end if
      end do
    end do

    ! coarse export: distinct local coarse ids per fine-export neighbor (count, fill)
    cmark = 0
    do neib = 1, nnb
      ks = cmt%export_index(neib-1)+1; ke = cmt%export_index(neib)
      cnt = 0
      do e = ks, ke
        c = aggr(cmt%export_item(e))
        if (c > 0 .and. cmark(c) /= neib) then; cmark(c) = neib; cnt = cnt + 1; end if
      end do
      cmt_c%export_index(neib) = cmt_c%export_index(neib-1) + cnt
    end do
    allocate(cmt_c%export_item(cmt_c%export_index(nnb)))
    cmark = 0
    do neib = 1, nnb
      ks = cmt%export_index(neib-1)+1; ke = cmt%export_index(neib)
      pos = cmt_c%export_index(neib-1)
      do e = ks, ke
        c = aggr(cmt%export_item(e))
        if (c > 0 .and. cmark(c) /= neib) then
          cmark(c) = neib; pos = pos + 1; cmt_c%export_item(pos) = c
        end if
      end do
    end do

    deallocate(agg_global, counts, displs, halo_gid, cmark)
  end subroutine hecmw_saamg_coarsen_struct

  !> Distributed tentative prolongator with halo-extended rows.
  !! Builds P-hat over this rank's internal aggregates (per-aggregate QR), then
  !! exchanges each fine node's nb x m P-hat block so that halo fine rows carry
  !! the owner's block.  The result phat_ext (nnode*nb x ncnode*m, columns in the
  !! localized coarse numbering) is what the smoother's A*P-hat needs, since A has
  !! halo columns.  bcoarse (naggr_local*m x m) is the coarse near-kernel for the
  !! rows this rank owns.
  subroutine hecmw_saamg_tentative_ext(cmt, nb, m, bfine_int, aggr, naggr_local, &
       agg_loc, ncnode, phat_ext, bcoarse)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: nb, m, naggr_local, ncnode
    real(kind=kreal),       intent(in)  :: bfine_int(:,:)     ! (nint*nb, m)
    integer(kind=kint),     intent(in)  :: aggr(:), agg_loc(:)
    type(hecmwST_saamg_bcsr), intent(out) :: phat_ext
    real(kind=kreal), allocatable, intent(out) :: bcoarse(:,:)
    type(hecmwST_saamg_bcsr) :: phat_int
    real(kind=kreal),   allocatable :: pvec(:), tv(:)
    integer(kind=kint), allocatable :: ti(:), tj(:)
    integer(kind=kint) :: nint, nnode, i, r, c, t, ci, cbase, nt, maxt
    real(kind=kreal) :: v

    nint  = cmt%nint
    nnode = cmt%nnode

    ! tentative on internal aggregates (columns 1..naggr_local*m)
    call hecmw_saamg_tentative(bfine_int, nint*nb, m, nb, aggr, naggr_local, phat_int, bcoarse)

    ! pack each internal node's nb x m tentative block into pvec, then fill halo
    ! blocks from owners.  An assigned node has exactly one block (its aggregate),
    ! and pvec's (c-1)*nb+r layout IS the column-major block layout -> direct copy.
    allocate(pvec(nb*m*nnode)); pvec = 0.0d0
    do i = 1, nint
      if (aggr(i) <= 0) cycle
      do t = phat_int%browptr(i), phat_int%browptr(i+1)-1
        pvec((i-1)*nb*m + 1 : i*nb*m) = phat_int%bval((t-1)*nb*m + 1 : t*nb*m)
      end do
    end do
    call hecmw_saamg_comm_update(cmt, nb*m, pvec)

    ! assemble phat_ext (nnode*nb x ncnode*m) with localized coarse columns
    maxt = nnode*nb*m
    allocate(ti(maxt), tj(maxt), tv(maxt)); nt = 0
    do i = 1, nnode
      ci = agg_loc(i)
      if (ci <= 0) cycle
      cbase = (ci-1)*m
      do r = 1, nb
        do c = 1, m
          v = pvec((i-1)*nb*m + (c-1)*nb + r)
          if (v /= 0.0d0) then
            nt = nt + 1; ti(nt) = (i-1)*nb + r; tj(nt) = cbase + c; tv(nt) = v
          end if
        end do
      end do
    end do
    call hecmw_saamg_bcsr_from_triplets(nnode*nb, ncnode*m, nb, ti, tj, tv, nt, phat_ext, mb=m)

    deallocate(pvec, ti, tj, tv)
    call hecmw_saamg_bcsr_free(phat_int)
  end subroutine hecmw_saamg_tentative_ext

  !> Distributed D-inner-product Lanczos estimate of lambda_max(D^{-1} A) using the
  !! halo-aware matvec and globally-reduced (allreduce) inner products, so all ranks
  !! agree on one lambda_max (required for a globally consistent prolongator-smoothing
  !! omega).
  function hecmw_saamg_lambda_max(A, cmt, D, niter) result(lambda)
    type(hecmwST_saamg_bcsr),      intent(in) :: A
    type(hecmwST_saamg_comm),      intent(in) :: cmt
    type(hecmwST_saamg_blockdiag), intent(in) :: D
    integer(kind=kint),          intent(in) :: niter
    real(kind=kreal) :: lambda
    real(kind=kreal), allocatable :: q(:), qp(:), w(:), aq(:), dv(:), alpha(:), beta(:)
    real(kind=kreal) :: qDq, wDw, bprev, al, be
    integer(kind=kint) :: nb, nint, nnode, ni, j, i, node, dof, gid, m

    nb = A%nb; nint = A%n / nb; nnode = A%ncol / nb; ni = nint*nb
    allocate(q(nnode*nb), qp(ni), w(ni), aq(ni), dv(ni), alpha(niter), beta(niter))
    q = 0.0d0
    ! partition-invariant seed: a deterministic function of the global node id, so
    ! lambda_max is identical regardless of how the mesh is partitioned (otherwise
    ! a partition-dependent estimate can make the Chebyshev smoother indefinite)
    do i = 1, ni
      node = (i-1)/nb + 1; dof = mod(i-1, nb) + 1
      if (allocated(cmt%gnode)) then; gid = cmt%gnode(node); else; gid = node; end if
      q(i) = sin(0.1d0*real(gid, kreal) + 0.37d0*real(dof, kreal) + 0.5d0)
    end do
    ! D-inner-product Lanczos for lambda_max(D^{-1}A): M = D^{-1}A is self-adjoint in
    ! <u,v>_D = u^T D v.  Build the Krylov tridiagonal T and take its largest eigen-
    ! value; this converges to lambda_max far faster than power iteration, so a short
    ! run is accurate even on large/ill-separated spectra (where power iteration
    ! underestimates and the Chebyshev smoother turns indefinite).
    call hecmw_saamg_blockdiag_matvec(D, q(1:ni), dv); qDq = dot_product(q(1:ni), dv)
    call hecmw_saamg_comm_allreduce_sum_r(cmt, qDq)
    if (qDq <= 0.0d0) then; lambda = 0.0d0; deallocate(q,qp,w,aq,dv,alpha,beta); return; end if
    q(1:ni) = q(1:ni) / sqrt(qDq); qp = 0.0d0; bprev = 0.0d0; m = 0

    do j = 1, niter
      call hecmw_saamg_matvec_d(cmt, A, q, aq)          ! aq = A q (q halo refreshed)
      call hecmw_saamg_blockdiag_apply(D, aq, w)        ! w  = D^{-1} A q = M q
      al = dot_product(aq, q(1:ni)); call hecmw_saamg_comm_allreduce_sum_r(cmt, al)
      alpha(j) = al                                     ! <M q, q>_D = (A q)^T q
      w(1:ni) = w(1:ni) - al*q(1:ni) - bprev*qp(1:ni)   ! 3-term recurrence
      m = j
      call hecmw_saamg_blockdiag_matvec(D, w, dv); wDw = dot_product(w, dv)
      call hecmw_saamg_comm_allreduce_sum_r(cmt, wDw)
      if (wDw <= 0.0d0) exit                            ! invariant subspace -> exact
      be = sqrt(wDw); beta(j) = be
      if (j == niter) exit
      qp(1:ni) = q(1:ni); bprev = be; q(1:ni) = w(1:ni) / be
    end do

    lambda = hecmw_saamg_tridiag_max_eig(alpha, beta, m)
    deallocate(q, qp, w, aq, dv, alpha, beta)
  end function hecmw_saamg_lambda_max

  !> Distributed smoothed prolongator P = (I - omega D^{-1} A) P-hat_ext.
  !! Builds the halo-extended tentative P-hat, the block-diagonal D of A, the
  !! shared lambda_max -> omega, then reuses the sequential smooth_prolongator
  !! (which, given the rectangular A with A%n = nint*nb rows, emits only internal
  !! fine rows of P while the SpGEMM A*P-hat_ext correctly sums over halo columns).
  subroutine hecmw_saamg_prolongator(A, cmt, m, safety, lanczos_iter, bfine_int, &
       aggr, naggr_local, agg_loc, ncnode, P, bcoarse, omega)
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: m, lanczos_iter, naggr_local, ncnode
    real(kind=kreal),       intent(in)  :: safety, bfine_int(:,:)
    integer(kind=kint),     intent(in)  :: aggr(:), agg_loc(:)
    type(hecmwST_saamg_bcsr), intent(out) :: P
    real(kind=kreal), allocatable, intent(out) :: bcoarse(:,:)
    real(kind=kreal),       intent(out) :: omega
    type(hecmwST_saamg_bcsr)      :: phat_ext
    type(hecmwST_saamg_blockdiag) :: D
    real(kind=kreal) :: lambda_hat

    call hecmw_saamg_tentative_ext(cmt, A%nb, m, bfine_int, aggr, naggr_local, &
         agg_loc, ncnode, phat_ext, bcoarse)
    call hecmw_saamg_blockdiag_setup(A, D)
    lambda_hat = safety * hecmw_saamg_lambda_max(A, cmt, D, lanczos_iter)
    omega = (4.0d0/3.0d0) / lambda_hat
    call hecmw_saamg_smooth_prolongator(A, D, omega, phat_ext, P)
    call hecmw_saamg_bcsr_free(phat_ext)
    call hecmw_saamg_blockdiag_free(D)
  end subroutine hecmw_saamg_prolongator

  !> Total number of coarse nodes across all ranks (= sum of naggr_local).
  function hecmw_saamg_nc_global(cmt, naggr_local) result(nc)
    type(hecmwST_saamg_comm), intent(in) :: cmt
    integer(kind=kint),     intent(in) :: naggr_local
    integer(kind=kint) :: nc, nprocs
    integer(kind=kint), allocatable :: counts(:)
    nprocs = hecmw_saamg_comm_size(cmt)
    allocate(counts(nprocs))
    call hecmw_saamg_comm_allgather_int(cmt, naggr_local, counts)
    nc = sum(counts(1:nprocs))
    deallocate(counts)
  end function hecmw_saamg_nc_global

  !> Distributed Galerkin via redundant assembly: each rank forms its partial
  !! contribution Ac_loc = (Pext_internal)^T (A Pext) in GLOBAL coarse-dof ids,
  !! then all partials are Allgatherv'd and summed (duplicates added) into the
  !! full global coarse operator Ac_global, identical on every rank.  This avoids
  !! distributed coarse comm tables and non-neighbor routing: the 2-ring coarse
  !! couplings enter through Pext's halo rows (columns of A*Pext) and the sum over
  !! every rank's partial reconstructs each global entry exactly once.
  !> This rank's partial coarse-operator contributions Ac_loc=(Pext_internal)^T(A Pext)
  !! as GLOBAL coarse-dof triplets (gti,gtj,gtv)[1:nt].  Shared by the redundant gather
  !! (galerkin_global, coarsest) and the distributed assembly (galerkin).  Rows are
  !! EXTENDED coarse nodes (some owned by other ranks) -> caller gathers or routes by
  !! row owner.  (Pext/C/Ac_loc construction unchanged from the original routine.)
  subroutine ac_global_triplets(A, cmt, P, naggr_local, my_off, chalo_gid, m, gti, gtj, gtv, nt, mem_level)
    type(hecmwST_saamg_bcsr), intent(in)  :: A, P
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: naggr_local, my_off, m
    integer(kind=kint),     intent(in)  :: chalo_gid(:)
    integer(kind=kint), allocatable, intent(out) :: gti(:), gtj(:)
    real(kind=kreal),   allocatable, intent(out) :: gtv(:)
    integer(kind=kint),              intent(out) :: nt
    integer(kind=kint),     intent(in)  :: mem_level   !< >=2 prints the [mem] peak probes (Pext/Ac_loc/buffers)
    type(hecmwST_saamg_bcsr) :: Pext, Pt, Ac_loc
    integer(kind=kint), allocatable :: gcolbuf(:), gext(:), extbuf(:), bcnt(:), boff(:)
    real(kind=kreal),   allocatable :: gvalbuf(:)
    integer(kind=kint) :: nb, nint, nnode, ntot, bs, i, t, base, s, cn, gcn, pos
    integer(kind=kint) :: nexth, gnode, ext, found, ncext, ib, ecn, r, cc, erow, ecol, b0, astat
    real(kind=kreal) :: v

    nb = P%nb; nint = cmt%nint; nnode = cmt%nnode; bs = nb*m

    ! --- exchange P's block rows: one nb x m block per (fine node, coarse node) ---
    ! VARIABLE-LENGTH exchange: each node carries exactly its own block count (no
    ! padding to a global maximum).  First exchange the per-node block counts, then
    ! the column ids and values packed contiguously per node.  This avoids the
    ! maxe*nnode padding (which dominated level-1 setup memory: most nodes have far
    ! fewer than the global-max blocks).
    allocate(bcnt(nnode), boff(nnode+1))
    bcnt = 0
    do i = 1, nint
      bcnt(i) = P%browptr(i+1) - P%browptr(i)
    end do
    call hecmw_saamg_comm_update_i(cmt, 1, bcnt)   ! halo nodes get owner's count
    boff(1) = 0
    do i = 1, nnode
      boff(i+1) = boff(i) + bcnt(i)
    end do
    ntot = boff(nnode+1)
    if (ntot < 1) ntot = 1
    allocate(gcolbuf(ntot), gvalbuf(ntot*bs), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'galerkin P-row exchange (gcolbuf/gvalbuf)')
    do i = 1, nint                                 ! pack internal P rows, no padding
      base = boff(i)
      do t = P%browptr(i), P%browptr(i+1)-1
        s = t - P%browptr(i) + 1
        cn = P%bcol(t)
        if (cn <= naggr_local) then; gcn = my_off + cn; else; gcn = chalo_gid(cn - naggr_local); end if
        gcolbuf(base + s) = gcn
        gvalbuf((base+s-1)*bs + 1 : (base+s)*bs) = P%bval((t-1)*bs+1 : t*bs)
      end do
    end do
    call hecmw_saamg_comm_update_var_i(cmt, bcnt, boff, gcolbuf)
    call hecmw_saamg_comm_update_var(cmt, bs, bcnt, boff, gvalbuf)

    ! --- localize received coarse nodes and build Pext DIRECTLY from gvalbuf ---
    ! (no intermediate block-triplet buffer: P-rows are already grouped by node, so
    !  Pext's block storage is filled in place -- avoids holding ~2x Pext on level 1).
    allocate(gext(ntot), extbuf(ntot), stat=astat); nexth = 0
    call hecmw_saamg_check_alloc(astat, 'galerkin Pext localize (gext/extbuf)')
    allocate(Pext%browptr(nnode+1)); Pext%browptr(1) = 1
    do i = 1, nnode                              ! pass 1: localize coarse nodes
      base = boff(i)
      do s = 1, bcnt(i)
        gnode = gcolbuf(base+s)
        if (gnode > my_off .and. gnode <= my_off+naggr_local) then
          ext = gnode - my_off
        else
          found = 0
          do t = 1, nexth
            if (gext(t) == gnode) then; found = t; exit; end if
          end do
          if (found == 0) then; nexth = nexth + 1; gext(nexth) = gnode; found = nexth; end if
          ext = naggr_local + found
        end if
        extbuf(base+s) = ext
      end do
      Pext%browptr(i+1) = Pext%browptr(i) + bcnt(i)
    end do
    ncext = naggr_local + nexth
    Pext%n = nnode*nb; Pext%ncol = ncext*m; Pext%nb = nb; Pext%mb = m
    Pext%nbrow = nnode; Pext%nbcol = ncext
    Pext%nnzb = Pext%browptr(nnode+1) - 1; Pext%nnz = bs*Pext%nnzb
    allocate(Pext%bcol(Pext%nnzb), Pext%bval(bs*Pext%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'galerkin Pext (bcol/bval)')
    do i = 1, nnode                              ! pass 2: fill bcol/bval from gvalbuf
      base = boff(i); pos = Pext%browptr(i)
      do s = 1, bcnt(i)
        Pext%bcol(pos) = extbuf(base+s)
        Pext%bval((pos-1)*bs+1 : pos*bs) = gvalbuf((base+s-1)*bs + 1 : (base+s)*bs)
        pos = pos + 1
      end do
    end do
    deallocate(gcolbuf, gvalbuf, extbuf, bcnt, boff)
    if (mem_level >= 2) call report_mem(cmt, '  galerkin: after Pext  ')

    ! FUSED Galerkin: Ac_loc = Pint^T A Pext computed directly (no A*Pext intermediate
    ! -- the large level-1 product that dominates setup memory).  Pt = (Pext's internal
    ! rows)^T = Pint^T, formed without an explicit Pint copy.
    call hecmw_saamg_transpose_blk_rows(Pext, nint, Pt)
    call hecmw_saamg_triple_blk(Pt, A, Pext, Ac_loc)
    if (mem_level >= 2) call report_mem(cmt, '  galerkin: after Ac_loc')
    call hecmw_saamg_bcsr_free(Pt); call hecmw_saamg_bcsr_free(Pext)

    ! --- emit Ac_loc as GLOBAL scalar triplets (block read; skip exact zeros) ---
    allocate(gti(Ac_loc%nnzb*m*m), gtj(Ac_loc%nnzb*m*m), gtv(Ac_loc%nnzb*m*m), stat=astat); nt = 0
    call hecmw_saamg_check_alloc(astat, 'galerkin coarse triplets (gti/gtj/gtv)')
    do ib = 1, Ac_loc%nbrow
      erow = e2g(ib)
      do t = Ac_loc%browptr(ib), Ac_loc%browptr(ib+1)-1
        ecn = Ac_loc%bcol(t); ecol = e2g(ecn); b0 = (t-1)*m*m
        ! emit ALL m*m entries (no exact-zero skip): the triplet SEQUENCE must be
        ! purely structural so the refresh routing plan (sperm/acc_map, built once)
        ! stays valid when only values change.
        do cc = 1, m
          do r = 1, m
            v = Ac_loc%bval(b0 + (cc-1)*m + r)
            nt = nt + 1
            gti(nt) = (erow-1)*m + r; gtj(nt) = (ecol-1)*m + cc; gtv(nt) = v
          end do
        end do
      end do
    end do
    call hecmw_saamg_bcsr_free(Ac_loc); deallocate(gext)
  contains
    integer(kind=kint) function e2g(e) result(g)
      integer(kind=kint), intent(in) :: e
      if (e <= naggr_local) then; g = my_off + e; else; g = gext(e - naggr_local); end if
    end function e2g
  end subroutine ac_global_triplets

  !> Redundant Galerkin: gather every rank's partial triplets and sum -> the full
  !! global coarse operator (identical on every rank).  Used ONLY at the coarsest
  !! level now (the distributed levels use hecmw_saamg_galerkin).
  subroutine hecmw_saamg_galerkin_global(A, cmt, P, naggr_local, ncnode, &
       my_off, chalo_gid, m, Ac_global)
    type(hecmwST_saamg_bcsr), intent(in)  :: A, P
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: naggr_local, ncnode, my_off, m
    integer(kind=kint),     intent(in)  :: chalo_gid(:)
    type(hecmwST_saamg_bcsr), intent(out) :: Ac_global
    integer(kind=kint), allocatable :: gti(:), gtj(:), agti(:), agtj(:)
    real(kind=kreal),   allocatable :: gtv(:), agtv(:)
    integer(kind=kint) :: nt, ntot, nc_global, idummy
    idummy = ncnode
    call ac_global_triplets(A, cmt, P, naggr_local, my_off, chalo_gid, m, gti, gtj, gtv, nt, 0)
    nc_global = hecmw_saamg_nc_global(cmt, naggr_local)
    call hecmw_saamg_comm_allgatherv_triplets(cmt, nt, gti, gtj, gtv, ntot, agti, agtj, agtv)
    call hecmw_saamg_bcsr_from_triplets(nc_global*m, nc_global*m, m, agti, agtj, agtv, ntot, Ac_global)
    deallocate(gti, gtj, gtv, agti, agtj, agtv)
  end subroutine hecmw_saamg_galerkin_global

  !> Distributed Galerkin: assemble THIS rank's owned rows of the coarse operator
  !! A_next (rectangular naggr_local*m x ncnode_op*m) by routing each rank's partial
  !! triplets to the OWNING rank (MPI_Alltoallv; owners may be non-neighbors), summing
  !! duplicates; and build the coarse OPERATOR comm table cmt_op from A_next's column
  !! sparsity.  Each global Ac(I,J) is summed exactly once on its owner.
  subroutine hecmw_saamg_galerkin(A, cmt, P, naggr_local, ncnode, my_off, chalo_gid, m, &
       A_next, cmt_op, gplan, mem_level)
    type(hecmwST_saamg_bcsr), intent(in)  :: A, P
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: naggr_local, ncnode, my_off, m
    integer(kind=kint),     intent(in)  :: chalo_gid(:)
    type(hecmwST_saamg_bcsr), intent(out) :: A_next
    type(hecmwST_saamg_comm), intent(out) :: cmt_op
    type(hecmwST_saamg_gplan), optional, intent(out) :: gplan  !< S4: record refresh routing
    integer(kind=kint), optional, intent(in) :: mem_level    !< >=2 prints [mem] peak probes (default off)
    integer(kind=kint), allocatable :: gti(:), gtj(:), si(:), sj(:), ri(:), rj(:)
    real(kind=kreal),   allocatable :: gtv(:), sv(:), rv(:)
    integer(kind=kint), allocatable :: coff(:), scnt(:), spos(:), lti(:), ltj(:), chalo_op(:)
    real(kind=kreal),   allocatable :: ltv(:)
    integer(kind=kint), allocatable :: sperm(:)
    integer(kind=kint) :: nt, nprocs, jp, t, gnode, r, nrecv, grow, gcol, gcnode, gcdof
    integer(kind=kint) :: lcnode, found, nhc, ncnode_op, mlv

    mlv = 0; if (present(mem_level)) mlv = mem_level
    call ac_global_triplets(A, cmt, P, naggr_local, my_off, chalo_gid, m, gti, gtj, gtv, nt, mlv)

    nprocs = hecmw_saamg_comm_size(cmt)
    allocate(coff(0:nprocs))
    call coarse_offsets(cmt, naggr_local, nprocs, coff)

    ! route triplets to row owners (bucket by owner, Alltoallv)
    allocate(scnt(nprocs)); scnt = 0
    do t = 1, nt
      gnode = (gti(t)-1)/m + 1; r = owner_of(gnode, coff, nprocs); scnt(r+1) = scnt(r+1) + 1
    end do
    allocate(si(max(nt,1)), sj(max(nt,1)), sv(max(nt,1)), spos(nprocs))
    if (present(gplan)) allocate(sperm(max(nt,1)))
    spos(1) = 0
    do jp = 2, nprocs; spos(jp) = spos(jp-1) + scnt(jp-1); end do
    do t = 1, nt
      gnode = (gti(t)-1)/m + 1; r = owner_of(gnode, coff, nprocs)
      spos(r+1) = spos(r+1) + 1
      si(spos(r+1)) = gti(t); sj(spos(r+1)) = gtj(t); sv(spos(r+1)) = gtv(t)
      if (present(gplan)) sperm(t) = spos(r+1)
    end do
    deallocate(gti, gtj, gtv)
    if (mlv >= 2) call report_mem(cmt, '  galerkin: send-buf si ')
    call hecmw_saamg_comm_alltoallv_triplets(cmt, scnt, si, sj, sv, nrecv, ri, rj, rv)
    deallocate(si, sj, sv, spos)
    if (mlv >= 2) call report_mem(cmt, '  galerkin: recv-buf ri ')

    ! localize columns -> owned (1..naggr_local) + halo coarse (naggr_local+1..).
    ! PRE-SEED the halo set with P's own coarse columns (the 1-ring set, in P's
    ! order) so that cmt_op's first ncnode nodes == P's coarse nodes with identical
    ! local ids; any extra 2-ring nodes from the operator are appended.  Then P's
    ! columns are directly valid in cmt_op (no re-localization needed in the V-cycle).
    allocate(lti(max(nrecv,1)), ltj(max(nrecv,1)), ltv(max(nrecv,1)), chalo_op(max(nrecv+ncnode,1)))
    nhc = ncnode - naggr_local
    do jp = 1, nhc; chalo_op(jp) = chalo_gid(jp); end do
    do t = 1, nrecv
      grow = ri(t); gcol = rj(t)
      gcnode = (gcol-1)/m + 1; gcdof = mod(gcol-1, m) + 1
      if (gcnode > my_off .and. gcnode <= my_off+naggr_local) then
        lcnode = gcnode - my_off
      else
        found = 0
        do jp = 1, nhc
          if (chalo_op(jp) == gcnode) then; found = jp; exit; end if
        end do
        if (found == 0) then; nhc = nhc + 1; chalo_op(nhc) = gcnode; found = nhc; end if
        lcnode = naggr_local + found
      end if
      lti(t) = grow - my_off*m                 ! local owned row dof (1..naggr_local*m)
      ltj(t) = (lcnode-1)*m + gcdof; ltv(t) = rv(t)
    end do
    deallocate(ri, rj, rv)
    ncnode_op = naggr_local + nhc
    call hecmw_saamg_bcsr_from_triplets(naggr_local*m, ncnode_op*m, m, lti, ltj, ltv, nrecv, A_next)

    ! S4: record the refresh routing (scnt + send permutation + value-slot map).
    ! acc_map(k) = the value slot in A_next for received triplet k (lti,ltj).
    if (present(gplan)) then
      gplan%nt = nt; gplan%nrecv = nrecv
      call move_alloc(scnt, gplan%scnt)
      call move_alloc(sperm, gplan%sperm)
      allocate(gplan%acc_map(max(nrecv,1)))
      call build_accmap(A_next, m, lti, ltj, nrecv, gplan%acc_map)
      gplan%ready = .true.
    else
      deallocate(scnt)
    end if
    deallocate(lti, ltj, ltv)

    call build_operator_commtable(cmt, naggr_local, my_off, m, ncnode_op, chalo_op, nhc, coff, nprocs, cmt_op)
    deallocate(chalo_op, coff)
  end subroutine hecmw_saamg_galerkin

  !> Numeric-only Galerkin refresh (S4): recompute the coarse operator VALUES into
  !! the existing A_next pattern using the cached routing plan.  Skips the operator
  !! comm-table handshake, column localization and triplet sort -- only the values
  !! are recomputed (SpGEMM) and re-routed (values-only Alltoallv).  A_next must be
  !! the pattern produced by the matching hecmw_saamg_galerkin call.
  subroutine hecmw_saamg_galerkin_refresh(A, cmt, P, naggr_local, my_off, chalo_gid, m, &
       gplan, A_next)
    type(hecmwST_saamg_bcsr),  intent(in)    :: A, P
    type(hecmwST_saamg_comm),  intent(in)    :: cmt
    integer(kind=kint),      intent(in)    :: naggr_local, my_off, m
    integer(kind=kint),      intent(in)    :: chalo_gid(:)
    type(hecmwST_saamg_gplan), intent(in)    :: gplan
    type(hecmwST_saamg_bcsr),  intent(inout) :: A_next  !< pattern reused, values overwritten
    integer(kind=kint), allocatable :: gti(:), gtj(:)
    real(kind=kreal),   allocatable :: gtv(:), sv(:), rv(:)
    integer(kind=kint) :: nt, t, nrecv

    call ac_global_triplets(A, cmt, P, naggr_local, my_off, chalo_gid, m, gti, gtj, gtv, nt, 0)
    deallocate(gti, gtj)                       ! positions unchanged -> not needed
    ! reorder values into the cached send-buffer order, route values only
    allocate(sv(max(nt,1)))
    do t = 1, nt; sv(gplan%sperm(t)) = gtv(t); end do
    deallocate(gtv)
    call hecmw_saamg_comm_alltoallv_real(cmt, gplan%scnt, sv, nrecv, rv)
    deallocate(sv)
    ! scatter-accumulate into the fixed A_next block pattern (bval slots via acc_map)
    A_next%bval(:) = 0.0d0
    do t = 1, nrecv
      A_next%bval(gplan%acc_map(t)) = A_next%bval(gplan%acc_map(t)) + rv(t)
    end do
    deallocate(rv)
  end subroutine hecmw_saamg_galerkin_refresh

  !> Free a Galerkin refresh plan.
  subroutine hecmw_saamg_gplan_free(gp)
    type(hecmwST_saamg_gplan), intent(inout) :: gp
    if (allocated(gp%scnt))    deallocate(gp%scnt)
    if (allocated(gp%sperm))   deallocate(gp%sperm)
    if (allocated(gp%acc_map)) deallocate(gp%acc_map)
    gp%ready = .false.; gp%nt = 0; gp%nrecv = 0
  end subroutine hecmw_saamg_gplan_free

  !> Build acc_map(k) = the slot in A%bval where received triplet (lti(k),ltj(k))
  !! lands: locate the block (block row ib, block col jb) by linear search (coarse
  !! block rows are short), then the entry (r,c) within that m x m block.  Valid for
  !! any refresh as long as the triplet positions (lti,ltj order) are unchanged.
  subroutine build_accmap(A, m, lti, ltj, nrecv, acc_map)
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    integer(kind=kint),    intent(in)  :: m, lti(:), ltj(:), nrecv
    integer(kind=kint),    intent(out) :: acc_map(:)
    integer(kind=kint) :: t, ib, jb, r, cc, k
    do t = 1, nrecv
      ib = (lti(t)-1)/m + 1; r  = mod(lti(t)-1, m) + 1
      jb = (ltj(t)-1)/m + 1; cc = mod(ltj(t)-1, m) + 1
      acc_map(t) = 0
      do k = A%browptr(ib), A%browptr(ib+1)-1
        if (A%bcol(k) == jb) then
          acc_map(t) = (k-1)*m*m + (cc-1)*m + r; exit
        end if
      end do
    end do
  end subroutine build_accmap

  !> S1 self-check: build the coarse operator BOTH ways -- distributed
  !! (hecmw_saamg_galerkin, with row-routing + operator comm table) and
  !! redundant (hecmw_saamg_galerkin_global, Allgatherv).  Gather the
  !! distributed one and confirm it equals the redundant one (probe by matvec on a
  !! seeded global coarse vector).  Also check the operator comm table is symmetric.
  subroutine hecmw_saamg_verify_galerkin(A, cmt, P, naggr_local, ncnode, my_off, &
       chalo_gid, m, reldiff, sym_ok)
    type(hecmwST_saamg_bcsr), intent(in)  :: A, P
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: naggr_local, ncnode, my_off, m
    integer(kind=kint),     intent(in)  :: chalo_gid(:)
    real(kind=kreal),       intent(out) :: reldiff
    logical,                intent(out) :: sym_ok
    type(hecmwST_saamg_bcsr) :: Ac_global, A_next, Ac_dist
    type(hecmwST_saamg_comm) :: cmt_op
    integer(kind=kint), allocatable :: gti(:), gtj(:), agti(:), agtj(:)
    real(kind=kreal),   allocatable :: gtv(:), agtv(:), v(:), y1(:), y2(:)
    integer(kind=kint) :: nc_global, nt, ntot, ib, t, gnode, lcn, gdof, c, d, r, b0, ndim, okint
    real(kind=kreal) :: dmax, ymax, vv

    call hecmw_saamg_galerkin_global(A, cmt, P, naggr_local, ncnode, my_off, chalo_gid, m, Ac_global)
    call hecmw_saamg_galerkin(A, cmt, P, naggr_local, ncnode, my_off, chalo_gid, m, A_next, cmt_op)
    call hecmw_saamg_verify_commtable(cmt_op, sym_ok)
    okint = 0; if (.not. sym_ok) okint = 1
    call hecmw_saamg_comm_allreduce_max_int(cmt, okint); sym_ok = (okint == 0)

    ! gather the distributed A_next (owned rows, global ids) into a global square
    nc_global = hecmw_saamg_nc_global(cmt, naggr_local)
    allocate(gti(max(A_next%nnzb*m*m,1)), gtj(max(A_next%nnzb*m*m,1)), gtv(max(A_next%nnzb*m*m,1))); nt = 0
    do ib = 1, A_next%nbrow
      do t = A_next%browptr(ib), A_next%browptr(ib+1)-1
        lcn = A_next%bcol(t); gnode = cmt_op%gnode(lcn); b0 = (t-1)*m*m
        do gdof = 1, m
          do r = 1, m
            vv = A_next%bval(b0 + (gdof-1)*m + r)
            if (vv == 0.0d0) cycle
            nt = nt + 1
            gti(nt) = my_off*m + (ib-1)*m + r; gtj(nt) = (gnode-1)*m + gdof; gtv(nt) = vv
          end do
        end do
      end do
    end do
    call hecmw_saamg_comm_allgatherv_triplets(cmt, nt, gti, gtj, gtv, ntot, agti, agtj, agtv)
    call hecmw_saamg_bcsr_from_triplets(nc_global*m, nc_global*m, m, agti, agtj, agtv, ntot, Ac_dist)

    ! probe equality via matvec on a seeded vector (both are redundant squares)
    ndim = nc_global*m
    allocate(v(ndim), y1(ndim), y2(ndim))
    do c = 1, nc_global
      do d = 1, m; v((c-1)*m+d) = sin(0.21d0*real(c,kreal) + 0.13d0*real(d,kreal)); end do
    end do
    call hecmw_saamg_matvec(Ac_dist,   v, y1)
    call hecmw_saamg_matvec(Ac_global, v, y2)
    dmax = 0.0d0; ymax = 0.0d0
    do c = 1, ndim; dmax = max(dmax, abs(y1(c)-y2(c))); ymax = max(ymax, abs(y2(c))); end do
    reldiff = dmax / max(ymax, 1.0d-300)

    deallocate(gti, gtj, gtv, agti, agtj, agtv, v, y1, y2)
    call hecmw_saamg_bcsr_free(Ac_global); call hecmw_saamg_bcsr_free(A_next)
    call hecmw_saamg_bcsr_free(Ac_dist);   call hecmw_saamg_comm_free(cmt_op)
  end subroutine hecmw_saamg_verify_galerkin

  !> coff(0:nprocs): prefix sum of naggr_local over ranks (coarse-node offsets).
  subroutine coarse_offsets(cmt, naggr_local, nprocs, coff)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: naggr_local, nprocs
    integer(kind=kint),     intent(out) :: coff(0:nprocs)
    integer(kind=kint), allocatable :: counts(:)
    integer(kind=kint) :: r
    allocate(counts(nprocs))
    call hecmw_saamg_comm_allgather_int(cmt, naggr_local, counts)
    coff(0) = 0
    do r = 1, nprocs; coff(r) = coff(r-1) + counts(r); end do
    deallocate(counts)
  end subroutine coarse_offsets

  !> Owner rank (0-based) of global coarse node gnode: coff(r) < gnode <= coff(r+1).
  integer(kind=kint) function owner_of(gnode, coff, nprocs) result(r)
    integer(kind=kint), intent(in) :: gnode, nprocs
    integer(kind=kint), intent(in) :: coff(0:)
    integer(kind=kint) :: lo, hi, mid
    lo = 0; hi = nprocs-1
    do while (lo < hi)
      mid = (lo+hi)/2
      if (gnode <= coff(mid+1)) then; hi = mid; else; lo = mid+1; end if
    end do
    r = lo
  end function owner_of

  !> Build the coarse-level operator comm table from A_next's halo coarse columns
  !! (chalo_op(1:nhc) = global node id per halo coarse node, local id naggr_local+t).
  !! import = halo coarse grouped by owner; export = dual via a symmetric count/id
  !! handshake (Alltoall counts then Alltoallv requested global ids).  Owners may be
  !! non-neighbors.  import_item[k] <-> owner's export_item[k] match by send order.
  subroutine build_operator_commtable(cmt, naggr_local, my_off, m, ncnode_op, chalo_op, nhc, &
       coff, nprocs, cmt_op)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: naggr_local, my_off, m, ncnode_op, nhc, nprocs
    integer(kind=kint),     intent(in)  :: chalo_op(:), coff(0:)
    type(hecmwST_saamg_comm), intent(out) :: cmt_op
    integer(kind=kint), allocatable :: want(:), recvwant(:), neigh(:), rank2nb(:)
    integer(kind=kint), allocatable :: scnt(:), spos(:), ipos(:), sbuf(:), rbuf(:), rdis(:)
    integer(kind=kint) :: t, r, nnb, p, ntot, i, src, kk, base, nb2

    allocate(want(nprocs), recvwant(nprocs)); want = 0
    do t = 1, nhc
      r = owner_of(chalo_op(t), coff, nprocs); want(r+1) = want(r+1) + 1
    end do
    call hecmw_saamg_comm_alltoall_int(cmt, want, recvwant)

    call hecmw_saamg_comm_free(cmt_op)
    cmt_op%comm = cmt%comm; cmt_op%my_rank = cmt%my_rank
    cmt_op%nint = naggr_local; cmt_op%nnode = ncnode_op; cmt_op%nb = m
    allocate(cmt_op%gnode(ncnode_op))
    do i = 1, naggr_local; cmt_op%gnode(i) = my_off + i; end do
    do t = 1, nhc; cmt_op%gnode(naggr_local + t) = chalo_op(t); end do

    nnb = 0
    do r = 0, nprocs-1
      if (r == cmt%my_rank) cycle
      if (want(r+1) > 0 .or. recvwant(r+1) > 0) nnb = nnb + 1
    end do
    cmt_op%n_neighbor = nnb
    if (nnb == 0) then; deallocate(want, recvwant); return; end if

    allocate(neigh(nnb), rank2nb(0:nprocs-1)); rank2nb = 0
    nb2 = 0
    do r = 0, nprocs-1
      if (r == cmt%my_rank) cycle
      if (want(r+1) > 0 .or. recvwant(r+1) > 0) then
        nb2 = nb2 + 1; neigh(nb2) = r; rank2nb(r) = nb2
      end if
    end do
    allocate(cmt_op%neighbor(nnb)); cmt_op%neighbor = neigh

    ! import: halo coarse local ids grouped by neighbor; also pack the global ids
    ! I want (bucketed by owner) to send -> owners build their export from it.
    allocate(cmt_op%import_index(0:nnb)); cmt_op%import_index(0) = 0
    allocate(cmt_op%export_index(0:nnb)); cmt_op%export_index(0) = 0
    do p = 1, nnb
      cmt_op%import_index(p) = cmt_op%import_index(p-1) + want(neigh(p)+1)
      cmt_op%export_index(p) = cmt_op%export_index(p-1) + recvwant(neigh(p)+1)
    end do
    allocate(cmt_op%import_item(cmt_op%import_index(nnb)))
    allocate(cmt_op%export_item(cmt_op%export_index(nnb)))

    allocate(scnt(nprocs), spos(nprocs), ipos(nnb), sbuf(max(nhc,1)))
    scnt = 0
    do t = 1, nhc
      r = owner_of(chalo_op(t), coff, nprocs); scnt(r+1) = scnt(r+1) + 1
    end do
    spos(1) = 0
    do p = 2, nprocs; spos(p) = spos(p-1) + scnt(p-1); end do
    do p = 1, nnb; ipos(p) = cmt_op%import_index(p-1); end do
    do t = 1, nhc
      r = owner_of(chalo_op(t), coff, nprocs)
      spos(r+1) = spos(r+1) + 1; sbuf(spos(r+1)) = chalo_op(t)        ! global id wanted from r
      ipos(rank2nb(r)) = ipos(rank2nb(r)) + 1
      cmt_op%import_item(ipos(rank2nb(r))) = naggr_local + t          ! my local halo coarse id
    end do
    call hecmw_saamg_comm_alltoallv_int(cmt, scnt, sbuf, ntot, rbuf)  ! rbuf = ids others want from me

    ! rbuf ordered by source rank (per recvwant); map global id -> my local owned id
    allocate(rdis(nprocs)); rdis(1) = 0
    do p = 2, nprocs; rdis(p) = rdis(p-1) + recvwant(p-1); end do
    do src = 0, nprocs-1
      if (recvwant(src+1) == 0) cycle
      base = cmt_op%export_index(rank2nb(src)-1)
      do kk = 1, recvwant(src+1)
        cmt_op%export_item(base + kk) = rbuf(rdis(src+1) + kk) - my_off   ! local owned coarse id
      end do
    end do

    deallocate(want, recvwant, neigh, rank2nb, scnt, spos, ipos, sbuf, rbuf, rdis)
  end subroutine build_operator_commtable

  !> Distributed Chebyshev smoother: x <- x + p_k(D^{-1}A)(b - A x), using the
  !! halo-aware matvec.  b is the internal RHS (length nint*nb); x has halo room
  !! (length nnode*nb) and only its internal part is the meaningful result.
  !! Optional r_out returns the final residual b - A x (length nint*nb), maintained
  !! incrementally by the iteration, so the caller can restrict it without a fresh
  !! matvec (fused V-cycle: drops one fine matvec per cycle).
  subroutine hecmw_saamg_cheb_apply(cmt, A, D, cheb, b, x, zero_init, r, p, dr, ap, mv, r_out)
    type(hecmwST_saamg_comm),      intent(in)    :: cmt
    type(hecmwST_saamg_bcsr),      intent(in)    :: A
    type(hecmwST_saamg_blockdiag), intent(in)    :: D
    type(hecmwST_saamg_chebyshev), intent(in)    :: cheb
    real(kind=kreal),            intent(in)    :: b(:)
    real(kind=kreal),            intent(inout) :: x(:)
    logical,                     intent(in)    :: zero_init
    real(kind=kreal),            intent(inout) :: r(:), p(:), dr(:), ap(:)  !< hoisted scratch (sizes A%n, A%ncol, A%n, A%n)
    procedure(hecmw_saamg_matvec_if), optional :: mv  !< fast finest matvec; else scalar matvec_d
    real(kind=kreal),            optional, intent(out) :: r_out(:)  !< final residual b - A x
    real(kind=kreal) :: rho, rho2, c1, c2
    integer(kind=kint) :: n, nt, it
    logical :: use_mv

    use_mv = present(mv)
    n = A%n; nt = A%ncol
    ! Vector ops wrapped in !$acc kernels so the whole sweep stays device-resident
    ! between the device matvec and D^-1 (no host round-trip under managed memory).
    if (zero_init) then
      !$acc kernels
      x(1:n) = 0.0d0
      r(1:n) = b(1:n)
      !$acc end kernels
    else
      if (use_mv) then; call mv(x, ap); else; call hecmw_saamg_matvec_d(cmt, A, x, ap); end if
      !$acc kernels
      r(1:n) = b(1:n) - ap(1:n)
      !$acc end kernels
    end if
    call hecmw_saamg_blockdiag_apply(D, r, dr)
    !$acc kernels
    p = 0.0d0; p(1:n) = dr(1:n) / cheb%theta
    !$acc end kernels
    rho = 1.0d0 / cheb%sigma
    do it = 1, cheb%deg
      !$acc kernels
      x(1:n) = x(1:n) + p(1:n)
      !$acc end kernels
      ! The matvec/residual/search-direction below only feed the NEXT sweep (and,
      ! on the final sweep, r_out).  On the last sweep skip the residual matvec
      ! entirely when the caller does not want r_out (post-smoother) -- x is already
      ! the final smoothed result, so this is bit-identical and saves one fine matvec.
      if (it == cheb%deg .and. .not. present(r_out)) exit
      if (use_mv) then; call mv(p, ap); else; call hecmw_saamg_matvec_d(cmt, A, p, ap); end if
      !$acc kernels
      r(1:n) = r(1:n) - ap(1:n)
      !$acc end kernels
      if (it == cheb%deg) exit                  ! r_out wanted: residual done, no next-sweep p needed
      rho2 = 1.0d0 / (2.0d0*cheb%sigma - rho)
      call hecmw_saamg_blockdiag_apply(D, r, dr)
      c1 = rho2 * rho
      c2 = 2.0d0 * rho2 / cheb%delta
      !$acc kernels
      p(1:n) = c1 * p(1:n) + c2 * dr(1:n)
      !$acc end kernels
      rho = rho2
    end do
    if (present(r_out)) then
      !$acc kernels
      r_out(1:n) = r(1:n)   ! residual b - A x (incrementally maintained)
      !$acc end kernels
    end if
  end subroutine hecmw_saamg_cheb_apply

  !> Memory diagnostic: report resident-set size (VmRSS) summed and maxed over the
  !! ranks at a labelled point.  Enabled by the !SOLVER LOGLEVEL>=2; lets users see
  !! the setup memory high-water mark (e.g. to diagnose OOM at scale).  Reads
  !! /proc/self/status (Linux); on other platforms it reports "unavailable".
  subroutine report_mem(cmt, label)
    type(hecmwST_saamg_comm), intent(in) :: cmt
    character(len=*),       intent(in) :: label
    integer(kind=kint) :: iu, ios, kb
    character(len=128) :: line
    real(kind=kreal) :: tot, mx
    kb = 0
    open(newunit=iu, file='/proc/self/status', action='read', iostat=ios)
    if (ios == 0) then
      do
        read(iu, '(a)', iostat=ios) line; if (ios /= 0) exit
        if (line(1:6) == 'VmRSS:') then; read(line(7:), *) kb; exit; end if
      end do
      close(iu)
    end if
    tot = real(kb,kreal)/1024.0d0; mx = tot
    call hecmw_saamg_comm_allreduce_sum_r(cmt, tot)   ! collective: run on every rank
    call hecmw_saamg_comm_allreduce_max_r(cmt, mx)
    if (cmt%my_rank == 0) then
      if (tot > 0.0d0) then       ! VmRSS read succeeded (Linux)
        write(*,'(a,a,a,f9.1,a,f9.1,a)') '#### [mem] ', label, &
             ': total=', tot/1024.0d0, ' GB  max=', mx, ' MB/rank'
      else                        ! /proc/self/status unavailable (non-Linux)
        write(*,'(a,a,a)') '#### [mem] ', label, ': VmRSS unavailable (non-Linux platform)'
      end if
    end if
  end subroutine report_mem

  !> Wall-clock reading in seconds (system_clock), used for the [time] probes and
  !! the overall setup timer.  Monotonic; differences give elapsed wall time.
  function saamg_wtime() result(t)
    real(kind=kreal)   :: t
    integer(kind=8)    :: cc, rate
    call system_clock(cc, rate)
    t = real(cc,kreal) / real(rate,kreal)
  end function saamg_wtime

  !> Timing diagnostic: report the wall time of a labelled setup section at level
  !! lev.  Enabled by !SOLVER LOGLEVEL>=2 (same gate as the [mem] probes); the
  !! reported value is the max over ranks (the slowest rank gates a collective).
  subroutine report_time(cmt, lev, label, dt)
    type(hecmwST_saamg_comm), intent(in) :: cmt
    integer(kind=kint),     intent(in) :: lev
    character(len=*),       intent(in) :: label   !< section name (caller pads for alignment)
    real(kind=kreal),       intent(in) :: dt      !< local elapsed seconds
    real(kind=kreal) :: mx
    mx = dt
    call hecmw_saamg_comm_allreduce_max_r(cmt, mx)   ! collective: run on every rank
    if (cmt%my_rank == 0) then
      if (lev > 0) then   ! per-level section
        write(*,'(a,i0,1x,a,a,f9.3,a)') '#### [time] L', lev, label, ': ', mx, ' s'
      else                ! whole-setup total (lev=0)
        write(*,'(a,a,a,f9.3,a)') '#### [time] ', label, ': ', mx, ' s'
      end if
    end if
  end subroutine report_time

  !> Build a fully-distributed multilevel SA-AMG hierarchy.  bfine_int is the
  !! near-kernel on internal nodes, built with a GLOBALLY consistent reference
  !! (zero / global centroid) so owner and halo agree on shared nodes.
  subroutine hecmw_saamg_setup(A, cmt, bfine_int, m, prm, dh, fine_matvec, move_in)
    type(hecmwST_saamg_bcsr),      intent(inout) :: A   !< moved into the hierarchy if move_in
    type(hecmwST_saamg_comm),      intent(in)  :: cmt
    real(kind=kreal),            intent(in)  :: bfine_int(:,:)
    integer(kind=kint),          intent(in)  :: m
    type(hecmwST_saamg_params),    intent(in)  :: prm
    type(hecmwST_saamg_hier), intent(out) :: dh
    procedure(hecmw_saamg_matvec_if), optional :: fine_matvec  !< fast finest matvec (HEC-MW)
    logical, optional,           intent(in)  :: move_in  !< take ownership of A (move, no copy)
    real(kind=kreal), allocatable :: B(:,:), Bnext(:,:)
    real(kind=kreal) :: omega
    integer(kind=kint) :: l, nglob, ni, nnzg, astat
    logical :: stalled, take_ownership
    real(kind=kreal) :: t0, t5, ts

    t0 = saamg_wtime()
    dh%m = m; dh%prm = prm; dh%symmetric = prm%symmetric
    if (present(fine_matvec)) dh%fine_matvec => fine_matvec
    allocate(dh%lev(prm%max_level))
    if (prm%loglevel >= 2) call report_mem(cmt, 'setup entry            ')
    ! NB: Fortran .and. is NOT short-circuit, so `present(move_in) .and. move_in`
    ! may read move_in when absent (NULL deref).  Test present() first, separately.
    take_ownership = .false.
    if (present(move_in)) take_ownership = move_in
    if (take_ownership) then   ! take ownership: avoid holding A twice
      call hecmw_saamg_bcsr_move(A, dh%lev(1)%A)
    else
      call hecmw_saamg_bcsr_copy(A, dh%lev(1)%A)
    end if
    call hecmw_saamg_comm_copy(cmt, dh%lev(1)%cmt)
    ! NB: use dh%lev(1)%A (NOT A) from here -- A is emptied when move_in is set.
    dh%lev(1)%nb_l = dh%lev(1)%A%nb
    allocate(B(dh%lev(1)%A%n, m), stat=astat); call hecmw_saamg_check_alloc(astat, 'setup near-kernel B     ')
    B = bfine_int(1:dh%lev(1)%A%n, 1:m)
    if (prm%loglevel >= 2) call report_mem(cmt, 'after A0 block-CSR copy ')

    l = 1
    do
      ts = saamg_wtime()
      call level_smoother(dh%lev(l), prm, omega)
      if (prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'smoother      ', saamg_wtime()-ts)
      nglob = dh%lev(l)%A%n
      call hecmw_saamg_comm_allreduce_sum_int(dh%lev(l)%cmt, nglob)
      ! coarsest when the global size is small (but always coarsen >=1 so the
      ! coarsest is a coarse level with contiguous coarse offsets), or max depth.
      if ((nglob <= prm%coarse_size .and. l > 1) .or. l == prm%max_level) then
        ts = saamg_wtime()
        call make_coarsest(dh, l)
        if (prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'coarsest solve', saamg_wtime()-ts)
        exit
      end if
      call build_coarse_level(dh, l, B, m, omega, Bnext, stalled)
      if (stalled) then
        ts = saamg_wtime()
        call make_coarsest(dh, l)
        if (prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'coarsest solve', saamg_wtime()-ts)
        exit
      end if
      call move_alloc(Bnext, B)
      l = l + 1
    end do
    dh%nlevel = l
    deallocate(B)

    ! per-level V-cycle work vectors (x has halo room; rc spans the next level)
    do l = 1, dh%nlevel
      ni = dh%lev(l)%A%n
      allocate(dh%lev(l)%x(dh%lev(l)%cmt%nnode * dh%lev(l)%nb_l), stat=astat)
      call hecmw_saamg_check_alloc(astat, 'V-cycle work x')
      allocate(dh%lev(l)%rhs(ni), dh%lev(l)%res(ni), dh%lev(l)%tmp(ni), stat=astat)
      call hecmw_saamg_check_alloc(astat, 'V-cycle work rhs/res/tmp')
      ! Chebyshev scratch hoisted here (allocated once, reused every apply) so the
      ! per-apply path does no allocation -- avoids managed-malloc churn on GPU.
      allocate(dh%lev(l)%cheb_r(ni), dh%lev(l)%cheb_dr(ni), dh%lev(l)%cheb_ap(ni), &
               dh%lev(l)%cheb_p(dh%lev(l)%A%ncol), stat=astat)
      call hecmw_saamg_check_alloc(astat, 'V-cycle cheb scratch')
      if (l < dh%nlevel) then
        allocate(dh%lev(l)%rc(dh%lev(l+1)%cmt%nnode * m), stat=astat)
        call hecmw_saamg_check_alloc(astat, 'V-cycle work rc')
      end if
    end do

    call warn_if_indefinite(dh)
    t5 = saamg_wtime()
    if (prm%timelog >= 2) call report_time(dh%lev(1)%cmt, 0, 'setup total   ', t5-t0)
    if (prm%verbose) then    ! NB: collectives below run on ALL ranks; only the write is rank-0
      if (dh%lev(1)%cmt%my_rank == 0) &
        write(*,'(a,i0,a,i0,a,f7.3)') '#### SA-AMG setup: nlevel=', dh%nlevel, &
             '  cycle(gamma)=', dh%prm%ncycle, '  tot(s)=', t5-t0
      do l = 1, dh%nlevel
        nglob = dh%lev(l)%A%n; call hecmw_saamg_comm_allreduce_sum_int(dh%lev(l)%cmt, nglob)
        nnzg  = dh%lev(l)%A%nnz; call hecmw_saamg_comm_allreduce_sum_int(dh%lev(l)%cmt, nnzg)
        if (dh%lev(1)%cmt%my_rank == 0) write(*,'(a,i0,a,i0,a,i0)') &
             '####   level ', l, ': global dof=', nglob, '  global nnz=', nnzg
      end do
      if (dh%lev(1)%cmt%my_rank == 0) then
        if (dh%coarse_is_smoother) then
          write(*,'(a,i0,a,i0,a)') '####   coarsest solver: smoother (Chebyshev deg ', &
               dh%prm%cheb_deg, ', ', dh%nc_coarsest, ' dof)'
        else if (dh%coarse_is_mumps) then
          write(*,'(a,a,a,i0,a)') '####   coarsest solver: distributed MUMPS ', &
               merge('SYM=2', 'SYM=0', dh%symmetric), ' (', dh%nc_coarsest, ' dof)'
        else
          write(*,'(a,a,a,i0,a)') '####   coarsest solver: dense ', &
               merge('LDL^T', 'LU   ', dh%symmetric), ' (', dh%nc_coarsest, ' dof)'
        end if
      end if
    end if
  end subroutine hecmw_saamg_setup

  !> Build the level-l smoother: block-diagonal D, distributed lambda_max,
  !! Chebyshev coefficients; return the prolongator-smoothing omega.
  subroutine level_smoother(lev, prm, omega)
    type(hecmwST_saamg_dlevel), intent(inout) :: lev
    type(hecmwST_saamg_params), intent(in)    :: prm
    real(kind=kreal),         intent(out)   :: omega
    real(kind=kreal) :: lambda, lambda_hat
    call hecmw_saamg_blockdiag_free(lev%D)
    call hecmw_saamg_blockdiag_setup(lev%A, lev%D)
    lambda = hecmw_saamg_lambda_max(lev%A, lev%cmt, lev%D, prm%lanczos_iter)
    lambda_hat = prm%safety * lambda
    call hecmw_saamg_cheb_setup(lambda_hat, prm%cheb_deg, prm%cheb_alpha, lev%cheb)
    omega = (4.0d0/3.0d0) / lambda_hat
  end subroutine level_smoother

  !> Coarsen level l -> level l+1 (distributed): uncoupled aggregation + tentative
  !! P-hat + smoothed P + distributed Galerkin.  Bnext = owned coarse near-kernel.
  !! stalled = the coarsening barely reduced the size (treat l as the coarsest).
  subroutine build_coarse_level(dh, l, B, m, omega, Bnext, stalled)
    type(hecmwST_saamg_hier), intent(inout) :: dh
    integer(kind=kint),          intent(in)    :: l, m
    real(kind=kreal),            intent(in)    :: B(:,:), omega
    real(kind=kreal), allocatable, intent(out) :: Bnext(:,:)
    logical,                     intent(out)   :: stalled
    integer(kind=kint), allocatable :: aggr(:), agg_loc(:)
    real(kind=kreal),   allocatable :: bcoarse(:,:)
    integer(kind=kint) :: na_glob, nf_glob, i, nint
    real(kind=kreal) :: ts
    ts = saamg_wtime()
    call hecmw_saamg_coarsen_struct(dh%lev(l)%A, dh%lev(l)%cmt, dh%prm%min_size, &
         dh%prm%max_size, dh%prm%theta, m, aggr, dh%lev(l)%naggr_local, agg_loc, &
         dh%lev(l)%ncnode, dh%lev(l)%cmt_c, dh%lev(l)%my_off, dh%lev(l)%chalo_gid)
    ! capture the finest aggregation (global coarse-node id per internal node) for VTK
    if (l == 1 .and. dh%prm%dump_vtk) then
      nint = dh%lev(1)%A%n / dh%lev(1)%A%nb
      if (allocated(dh%aggr_fine)) deallocate(dh%aggr_fine)
      allocate(dh%aggr_fine(nint))
      do i = 1, nint
        if (aggr(i) > 0) then
          dh%aggr_fine(i) = dh%lev(1)%my_off + aggr(i)
        else
          dh%aggr_fine(i) = 0
        end if
      end do
    end if
    na_glob = dh%lev(l)%naggr_local
    call hecmw_saamg_comm_allreduce_sum_int(dh%lev(l)%cmt, na_glob)
    nf_glob = dh%lev(l)%A%n / dh%lev(l)%A%nb
    call hecmw_saamg_comm_allreduce_sum_int(dh%lev(l)%cmt, nf_glob)
    if (real(na_glob,kreal) > dh%prm%stagn * real(nf_glob,kreal)) then
      stalled = .true.; deallocate(aggr, agg_loc); return
    end if
    stalled = .false.
    if (dh%prm%loglevel >= 2) call report_mem(dh%lev(l)%cmt, 'L  after coarsen_struct')
    if (dh%prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'coarsen_struct', saamg_wtime()-ts)
    ts = saamg_wtime()
    call hecmw_saamg_tentative_ext(dh%lev(l)%cmt, dh%lev(l)%A%nb, m, B, aggr, &
         dh%lev(l)%naggr_local, agg_loc, dh%lev(l)%ncnode, dh%lev(l)%phat_ext, bcoarse)
    if (dh%prm%loglevel >= 2) call report_mem(dh%lev(l)%cmt, 'L  after tentative_ext ')
    if (dh%prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'tentative_ext ', saamg_wtime()-ts)
    ts = saamg_wtime()
    call hecmw_saamg_smooth_prolongator(dh%lev(l)%A, dh%lev(l)%D, omega, dh%lev(l)%phat_ext, dh%lev(l)%P)
    if (dh%prm%loglevel >= 2) call report_mem(dh%lev(l)%cmt, 'L  after smooth_prolong')
    if (dh%prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'smooth_prolong', saamg_wtime()-ts)
    ! Galerkin first (it builds its own internal P^T); defer dh%Pt (only the V-cycle
    ! apply needs it) so it is not held during the memory-heavy A*Pext / P^T*C products.
    ts = saamg_wtime()
    call hecmw_saamg_galerkin(dh%lev(l)%A, dh%lev(l)%cmt, dh%lev(l)%P, dh%lev(l)%naggr_local, &
         dh%lev(l)%ncnode, dh%lev(l)%my_off, dh%lev(l)%chalo_gid, m, dh%lev(l+1)%A, dh%lev(l+1)%cmt, &
         dh%lev(l)%gplan, mem_level=dh%prm%loglevel)
    if (dh%prm%loglevel >= 2) call report_mem(dh%lev(l)%cmt, 'L  after galerkin      ')
    if (dh%prm%timelog >= 2) call report_time(dh%lev(l)%cmt, l, 'galerkin      ', saamg_wtime()-ts)
    call hecmw_saamg_bcsr_transpose(dh%lev(l)%P, dh%lev(l)%Pt)
    dh%lev(l+1)%nb_l   = m
    dh%lev(l+1)%my_off = dh%lev(l)%my_off    ! owned coarse global offset for level l+1
    call move_alloc(bcoarse, Bnext)
    deallocate(aggr, agg_loc)
  end subroutine build_coarse_level

  !> Mark level l the coarsest and factor its operator.  With MUMPS (F7) the
  !! distributed coarse operator is factored in parallel WITHOUT a gather; otherwise
  !! it is gathered redundantly (tiny) and factored by dense LDL^T.
  subroutine make_coarsest(dh, l)
    type(hecmwST_saamg_hier), intent(inout) :: dh
    integer(kind=kint),          intent(in)    :: l
    type(hecmwST_saamg_bcsr) :: Ac
    integer(kind=kint) :: nglob_dof, nown
    logical :: use_mumps
    dh%lev(l)%is_coarsest = .true.
    ! smoother-as-coarse-solve (ML CoarseSolver=Smoother): no direct factorization --
    ! reuse the level's Chebyshev smoother (D/cheb already built by level_smoother).
    ! Fully local+halo, so it scales like the smoother (no gather / host solve / factor).
    ! Accurate only when the coarsest is small; coarse_size default (100) keeps it so.
    if (dh%prm%coarsest_solver == 1) then
      dh%coarse_is_smoother = .true.
      dh%coarse_is_mumps    = .false.
      nown = dh%lev(l)%A%n / dh%m
      dh%nc_coarsest = hecmw_saamg_nc_global(dh%lev(l)%cmt, nown) * dh%m
      return
    end if
    select case (dh%prm%coarsest_solver)
    case (2);    use_mumps = .false.                        ! force dense LDL^T
    case (3);    use_mumps = .true.                         ! force MUMPS (aborts if not built)
    case default; use_mumps = hecmw_saamg_cmumps_available() ! auto: MUMPS if available
    end select
    if (use_mumps) then
      dh%coarse_is_mumps = .true.
      ! the coarsest level is never coarsened, so its naggr_local is unset; the
      ! owned coarse-node count is A%n/m and its global offset is my_off (inherited).
      nown = dh%lev(l)%A%n / dh%m
      dh%nc_coarsest = hecmw_saamg_nc_global(dh%lev(l)%cmt, nown) * dh%m
      call hecmw_saamg_cmumps_setup(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%my_off, dh%m, &
           nown, dh%symmetric, dh%cmumps)
    else
      dh%coarse_is_mumps = .false.
      call gather_level_operator(dh%lev(l), Ac, nglob_dof)
      dh%nc_coarsest = nglob_dof
      call hecmw_saamg_coarse_setup(Ac, dh%symmetric, dh%coarse)
      call hecmw_saamg_bcsr_free(Ac)
      allocate(dh%rcg(nglob_dof), dh%zcg(nglob_dof))
    end if
  end subroutine make_coarsest

  !> Gather a distributed level operator into the full global square operator (in
  !! global node/dof numbering via the level comm table's gnode), for the coarsest
  !! dense solve.  nglob_dof = global dof count.
  subroutine gather_level_operator(lev, Ac, nglob_dof)
    type(hecmwST_saamg_dlevel), intent(in)  :: lev
    type(hecmwST_saamg_bcsr),   intent(out) :: Ac
    integer(kind=kint),       intent(out) :: nglob_dof
    integer(kind=kint), allocatable :: gti(:), gtj(:), agti(:), agtj(:)
    real(kind=kreal),   allocatable :: gtv(:), agtv(:)
    integer(kind=kint) :: nb_l, nglob_nodes, ib, t, k, cn, rdof, cdof, nt, ntot, grn, gcn, b0
    real(kind=kreal) :: v
    nb_l = lev%A%nb
    nglob_nodes = lev%A%nbrow
    call hecmw_saamg_comm_allreduce_sum_int(lev%cmt, nglob_nodes)
    nglob_dof = nglob_nodes * nb_l
    allocate(gti(max(lev%A%nnzb*nb_l*nb_l,1)), gtj(max(lev%A%nnzb*nb_l*nb_l,1)), &
             gtv(max(lev%A%nnzb*nb_l*nb_l,1))); nt = 0
    do ib = 1, lev%A%nbrow
      grn = lev%cmt%gnode(ib)
      do k = lev%A%browptr(ib), lev%A%browptr(ib+1)-1
        cn = lev%A%bcol(k); gcn = lev%cmt%gnode(cn); b0 = (k-1)*nb_l*nb_l
        do cdof = 1, nb_l
          do rdof = 1, nb_l
            v = lev%A%bval(b0 + (cdof-1)*nb_l + rdof)
            if (v == 0.0d0) cycle
            nt = nt + 1
            gti(nt) = (grn-1)*nb_l + rdof
            gtj(nt) = (gcn-1)*nb_l + cdof
            gtv(nt) = v
          end do
        end do
      end do
    end do
    call hecmw_saamg_comm_allgatherv_triplets(lev%cmt, nt, gti, gtj, gtv, ntot, agti, agtj, agtv)
    call hecmw_saamg_bcsr_from_triplets(nglob_dof, nglob_dof, nb_l, agti, agtj, agtv, ntot, Ac)
    deallocate(gti, gtj, gtv, agti, agtj, agtv)
  end subroutine gather_level_operator

  !> Numeric-only refresh (Newton): reuse each level's aggregation / comm table /
  !! tentative P-hat AND the cached Galerkin routing plan (S4), recompute only
  !! D / lambda / smoothed P / coarse operator VALUES.  No comm-table rebuild:
  !! the coarse operator pattern and every level's cmt are untouched, so only the
  !! values move (SpGEMM + values-only Alltoallv).
  subroutine hecmw_saamg_refresh(dh, A_new)
    type(hecmwST_saamg_hier), intent(inout) :: dh
    type(hecmwST_saamg_bcsr), intent(in)    :: A_new
    type(hecmwST_saamg_bcsr) :: Ac
    real(kind=kreal) :: omega
    integer(kind=kint) :: l, nglob_dof
    call hecmw_saamg_bcsr_copy(A_new, dh%lev(1)%A)
    do l = 1, dh%nlevel
      if (dh%lev(l)%is_coarsest) then
        if (dh%coarse_is_smoother) then
          ! rebuild D / lambda_max / Chebyshev on the refreshed coarsest operator
          call level_smoother(dh%lev(l), dh%prm, omega)
        else if (dh%coarse_is_mumps) then
          call hecmw_saamg_cmumps_refresh(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%my_off, &
               dh%m, dh%lev(l)%A%n / dh%m, dh%cmumps)
        else
          call gather_level_operator(dh%lev(l), Ac, nglob_dof)
          call hecmw_saamg_coarse_setup(Ac, dh%symmetric, dh%coarse)
          call hecmw_saamg_bcsr_free(Ac)
        end if
        exit
      end if
      call level_smoother(dh%lev(l), dh%prm, omega)
      call hecmw_saamg_bcsr_free(dh%lev(l)%P)
      call hecmw_saamg_smooth_prolongator(dh%lev(l)%A, dh%lev(l)%D, omega, dh%lev(l)%phat_ext, dh%lev(l)%P)
      call hecmw_saamg_bcsr_free(dh%lev(l)%Pt)
      call hecmw_saamg_bcsr_transpose(dh%lev(l)%P, dh%lev(l)%Pt)
      ! recompute A_{l+1} values into its existing pattern via the cached plan
      call hecmw_saamg_galerkin_refresh(dh%lev(l)%A, dh%lev(l)%cmt, dh%lev(l)%P, &
           dh%lev(l)%naggr_local, dh%lev(l)%my_off, dh%lev(l)%chalo_gid, dh%m, &
           dh%lev(l)%gplan, dh%lev(l+1)%A)
    end do
    call warn_if_indefinite(dh)
  end subroutine hecmw_saamg_refresh

  !> Warn (once, on rank 0) when the coarsest operator is indefinite, which means
  !! the input matrix is not SPD -- outside SA-AMG's design assumption.
  subroutine warn_if_indefinite(dh)
    type(hecmwST_saamg_hier), intent(in) :: dh
    integer(kind=kint) :: nneg
    if (dh%coarse_is_smoother) return   ! no factorization -> no inertia to inspect
    nneg = dh%coarse%n_neg
    if (dh%coarse_is_mumps) nneg = dh%cmumps%n_neg
    if (dh%lev(1)%cmt%my_rank /= 0 .or. nneg <= 0) return
    write(*,'(a,i0,a)') ' #### SA-AMG WARNING: coarse operator has ', &
         nneg, ' negative eigenvalue(s) -- the matrix appears non-SPD.'
    write(*,'(a)') ' #### SA-AMG assumes a symmetric positive-definite operator; CG may fail'
    write(*,'(a)') ' #### with "indefinite preconditioner".  Non-SPD systems (u-p, friction'
    write(*,'(a)') ' #### contact, some smoothed-strain elements) need an indefinite solver.'
  end subroutine warn_if_indefinite

  !> Apply one distributed multilevel V-cycle:  z = M^{-1} r  (internal, length nint*nb).
  subroutine hecmw_saamg_apply(dh, r, z)
    type(hecmwST_saamg_hier), intent(inout) :: dh
    real(kind=kreal),            intent(in)    :: r(:)
    real(kind=kreal),            intent(out)   :: z(:)
    integer(kind=kint) :: n
    n = dh%lev(1)%A%n
    !$acc kernels
    dh%lev(1)%rhs(1:n) = r(1:n)
    !$acc end kernels
    call mg_cycle(dh, 1)
    !$acc kernels
    z(1:n) = dh%lev(1)%x(1:n)
    !$acc end kernels
  end subroutine hecmw_saamg_apply

  !> One multigrid gamma-cycle at distributed level l (gamma=prm%ncycle: 1=V, 2=W):
  !! lev%x <- M_l^{-1} lev%rhs.  Restriction uses the operator comm table of level
  !! l+1 (reverse-accumulate to owners); prolongation forward-broadcasts the coarse
  !! halo with the same table.  W-cycle re-visits the coarse level with init_zero=.false.
  recursive subroutine mg_cycle(dh, l, init_zero)
    type(hecmwST_saamg_hier), intent(inout) :: dh
    integer(kind=kint),          intent(in)    :: l
    logical, optional,           intent(in)    :: init_zero  !< .false. => warm start (W-cycle re-visit)
    real(kind=kreal), allocatable :: gv(:)
    integer(kind=kint) :: ni, nc, m, ntot, off, i, rn, rdof, nbl, j, gamma
    logical :: use_fmv, zstart

    m = dh%m; ni = dh%lev(l)%A%n; nbl = dh%lev(l)%nb_l
    zstart = .true.; if (present(init_zero)) zstart = init_zero
    ! the finest operator equals hecMAT -> use HEC-MW's native blocked matvec there
    use_fmv = (l == 1 .and. associated(dh%fine_matvec) .and. .not. dh%hook_disabled)
    if (dh%lev(l)%is_coarsest) then
      if (dh%coarse_is_smoother) then
        ! coarsest "solve" = Chebyshev smoother sweeps (zero initial guess), same as a
        ! pre-smooth; fully local+halo, no gather/factor -- scales to high rank counts.
        call hecmw_saamg_cheb_apply(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%D, dh%lev(l)%cheb, &
             dh%lev(l)%rhs, dh%lev(l)%x, .true., &
             dh%lev(l)%cheb_r, dh%lev(l)%cheb_p, dh%lev(l)%cheb_dr, dh%lev(l)%cheb_ap)
      else if (dh%coarse_is_mumps) then
        ! distributed sparse direct: gather RHS to host, solve, scatter owned slice
        call hecmw_saamg_cmumps_solve(dh%lev(l)%cmt, dh%lev(l)%rhs, ni, dh%cmumps, dh%lev(l)%x)
      else
        ! redundant dense: gather full RHS, dense LDL^T, take this rank's owned slice
        call hecmw_saamg_comm_allgatherv_real(dh%lev(l)%cmt, ni, dh%lev(l)%rhs, ntot, gv)
        dh%rcg(1:dh%nc_coarsest) = gv(1:dh%nc_coarsest); deallocate(gv)
        call hecmw_saamg_coarse_solve(dh%coarse, dh%rcg, dh%zcg)
        do i = 1, ni
          rn = (i-1)/nbl + 1; rdof = mod(i-1, nbl) + 1
          dh%lev(l)%x(i) = dh%zcg((dh%lev(l)%cmt%gnode(rn)-1)*nbl + rdof)
        end do
      end if
      return
    end if
    ! pre-smoothing (zero initial guess).  The smoother maintains res = rhs - A x
    ! incrementally, so it is returned in res and restricted directly -- this fuses
    ! away the separate residual matvec (one fewer fine matvec per V-cycle).
    if (use_fmv) then
      call hecmw_saamg_cheb_apply(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%D, dh%lev(l)%cheb, &
           dh%lev(l)%rhs, dh%lev(l)%x, zstart, &
           dh%lev(l)%cheb_r, dh%lev(l)%cheb_p, dh%lev(l)%cheb_dr, dh%lev(l)%cheb_ap, &
           dh%fine_matvec, r_out=dh%lev(l)%res)
    else
      call hecmw_saamg_cheb_apply(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%D, dh%lev(l)%cheb, &
           dh%lev(l)%rhs, dh%lev(l)%x, zstart, &
           dh%lev(l)%cheb_r, dh%lev(l)%cheb_p, dh%lev(l)%cheb_dr, dh%lev(l)%cheb_ap, &
           r_out=dh%lev(l)%res)
    end if
    ! restrict: rc = P^T res (in P's coarse cols), accumulate halo coarse to owners
    !$acc kernels
    dh%lev(l)%rc(:) = 0.0d0
    !$acc end kernels
    call hecmw_saamg_matvec(dh%lev(l)%Pt, dh%lev(l)%res, dh%lev(l)%rc)
    call hecmw_saamg_comm_reverse_add(dh%lev(l+1)%cmt, m, dh%lev(l)%rc)
    nc = dh%lev(l+1)%A%n
    !$acc kernels
    dh%lev(l+1)%rhs(1:nc) = dh%lev(l)%rc(1:nc)
    !$acc end kernels
    ! coarse-grid correction, repeated gamma times per level (V=1, W=2).  The coarse
    ! RHS rhs_{l+1} stays FIXED across passes; each pass warm-starts x_{l+1} from the
    ! previous iterate (init_zero=.false. for j>1).  Directly above the coarsest level
    ! gamma is forced to 1: the coarsest is an exact solve of that fixed RHS, so a
    ! second visit would redo identical work.  M^{-1} stays a fixed linear map of r
    ! (top level always zero-start), so CG/BiCGSTAB remain valid.
    gamma = dh%prm%ncycle
    if (dh%lev(l+1)%is_coarsest) gamma = 1
    do j = 1, gamma
      call mg_cycle(dh, l+1, init_zero=(j==1))
    end do
    ! prolongate: broadcast coarse halo, x += P x_{l+1}
    call hecmw_saamg_comm_update(dh%lev(l+1)%cmt, m, dh%lev(l+1)%x)
    call hecmw_saamg_matvec(dh%lev(l)%P, dh%lev(l+1)%x, dh%lev(l)%tmp)
    !$acc kernels
    dh%lev(l)%x(1:ni) = dh%lev(l)%x(1:ni) + dh%lev(l)%tmp(1:ni)
    !$acc end kernels
    ! post-smoothing
    if (use_fmv) then
      call hecmw_saamg_cheb_apply(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%D, dh%lev(l)%cheb, &
           dh%lev(l)%rhs, dh%lev(l)%x, .false., &
           dh%lev(l)%cheb_r, dh%lev(l)%cheb_p, dh%lev(l)%cheb_dr, dh%lev(l)%cheb_ap, &
           dh%fine_matvec)
    else
      call hecmw_saamg_cheb_apply(dh%lev(l)%cmt, dh%lev(l)%A, dh%lev(l)%D, dh%lev(l)%cheb, &
           dh%lev(l)%rhs, dh%lev(l)%x, .false., &
           dh%lev(l)%cheb_r, dh%lev(l)%cheb_p, dh%lev(l)%cheb_dr, dh%lev(l)%cheb_ap)
    end if
  end subroutine mg_cycle

  subroutine hecmw_saamg_free(dh)
    type(hecmwST_saamg_hier), intent(inout) :: dh
    integer(kind=kint) :: l
    if (allocated(dh%lev)) then
      do l = 1, dh%nlevel
        call hecmw_saamg_bcsr_free(dh%lev(l)%A)
        call hecmw_saamg_comm_free(dh%lev(l)%cmt)
        call hecmw_saamg_blockdiag_free(dh%lev(l)%D)
        call hecmw_saamg_bcsr_free(dh%lev(l)%P)
        call hecmw_saamg_bcsr_free(dh%lev(l)%Pt)
        call hecmw_saamg_bcsr_free(dh%lev(l)%phat_ext)
        call hecmw_saamg_comm_free(dh%lev(l)%cmt_c)
        call hecmw_saamg_gplan_free(dh%lev(l)%gplan)
        if (allocated(dh%lev(l)%chalo_gid)) deallocate(dh%lev(l)%chalo_gid)
        if (allocated(dh%lev(l)%x))   deallocate(dh%lev(l)%x)
        if (allocated(dh%lev(l)%rhs)) deallocate(dh%lev(l)%rhs)
        if (allocated(dh%lev(l)%res)) deallocate(dh%lev(l)%res)
        if (allocated(dh%lev(l)%tmp)) deallocate(dh%lev(l)%tmp)
        if (allocated(dh%lev(l)%rc))  deallocate(dh%lev(l)%rc)
        if (allocated(dh%lev(l)%cheb_r))  deallocate(dh%lev(l)%cheb_r)
        if (allocated(dh%lev(l)%cheb_p))  deallocate(dh%lev(l)%cheb_p)
        if (allocated(dh%lev(l)%cheb_dr)) deallocate(dh%lev(l)%cheb_dr)
        if (allocated(dh%lev(l)%cheb_ap)) deallocate(dh%lev(l)%cheb_ap)
      end do
      deallocate(dh%lev)
    end if
    call hecmw_saamg_coarse_free(dh%coarse)
    call hecmw_saamg_cmumps_free(dh%cmumps)
    dh%coarse_is_mumps = .false.
    if (allocated(dh%rcg)) deallocate(dh%rcg)
    if (allocated(dh%zcg)) deallocate(dh%zcg)
    if (allocated(dh%aggr_fine)) deallocate(dh%aggr_fine)
    dh%nlevel = 0
  end subroutine hecmw_saamg_free


  !> Halo-exchange the smoothed prolongator's rows (F4c Galerkin needs P rows for
  !! halo fine nodes, since A has halo columns).  P columns are converted to GLOBAL
  !! coarse-dof ids before sending (each rank localizes differently); the owner's
  !! row arrives padded to maxe entries per fine row.  This routine VERIFIES the
  !! row exchange by a round-trip checksum: the owner also sends, as a per-node
  !! scalar, sum(P_row^2); the importer compares it to the checksum of the row it
  !! received.  mism (global max) must be ~0.  The packed buffers are the basis of
  !! the F4c distributed Galerkin (built next).
  subroutine hecmw_saamg_verify_prows(cmt, m, P, naggr_local, my_off, chalo_gid, mism)
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: m, naggr_local, my_off
    type(hecmwST_saamg_bcsr), intent(in)  :: P
    integer(kind=kint),     intent(in)  :: chalo_gid(:)
    real(kind=kreal),       intent(out) :: mism
    integer(kind=kint), allocatable :: gcolbuf(:)
    real(kind=kreal),   allocatable :: gvalbuf(:), csum(:)
    integer(kind=kint) :: nb, nint, nnode, maxe, bs, i, t, cnt, cn, gcn, off, voff, e
    real(kind=kreal) :: crowsum

    nb = P%nb; nint = cmt%nint; nnode = cmt%nnode; bs = nb*m
    maxe = 0
    do i = 1, nint
      maxe = max(maxe, P%browptr(i+1) - P%browptr(i))
    end do
    call hecmw_saamg_comm_allreduce_max_int(cmt, maxe)
    if (maxe < 1) maxe = 1

    allocate(gcolbuf(maxe*nnode), gvalbuf(maxe*bs*nnode))
    gcolbuf = 0; gvalbuf = 0.0d0
    do i = 1, nint
      off = (i-1)*maxe; voff = (i-1)*maxe*bs; cnt = 0
      do t = P%browptr(i), P%browptr(i+1)-1
        cnt = cnt + 1
        cn = P%bcol(t)
        if (cn <= naggr_local) then; gcn = my_off + cn; else; gcn = chalo_gid(cn - naggr_local); end if
        gcolbuf(off + cnt) = gcn
        gvalbuf(voff + (cnt-1)*bs + 1 : voff + cnt*bs) = P%bval((t-1)*bs+1 : t*bs)
      end do
    end do
    call hecmw_saamg_comm_update_i(cmt, maxe, gcolbuf)    ! halo rows: global coarse nodes
    call hecmw_saamg_comm_update(cmt, maxe*bs, gvalbuf)   ! halo rows: blocks

    ! owner checksum sum(P_row^2) per node (= sum of block entries^2), exchanged as a scalar
    allocate(csum(nnode)); csum = 0.0d0
    do i = 1, nint
      do t = P%browptr(i), P%browptr(i+1)-1
        do e = 1, bs
          csum(i) = csum(i) + P%bval((t-1)*bs+e)**2
        end do
      end do
    end do
    call hecmw_saamg_comm_update(cmt, 1, csum)

    mism = 0.0d0
    do i = nint+1, nnode
      crowsum = 0.0d0; voff = (i-1)*maxe*bs
      do e = 1, maxe*bs
        crowsum = crowsum + gvalbuf(voff+e)**2
      end do
      mism = max(mism, abs(crowsum - csum(i)))
    end do
    call hecmw_saamg_comm_allreduce_max_r(cmt, mism)
    deallocate(gcolbuf, gvalbuf, csum)
  end subroutine hecmw_saamg_verify_prows

  !> Check coarse comm-table symmetry: for each neighbor, the number of coarse
  !! nodes it exports to me must equal the number I import from it.  Returns the
  !! per-rank result in ok (collective consistency is checked by the caller).
  subroutine hecmw_saamg_verify_commtable(cmt_c, ok)
    type(hecmwST_saamg_comm), intent(in)  :: cmt_c
    logical,                intent(out) :: ok
    integer(kind=kint), allocatable :: sendv(:), recvv(:)
    integer(kind=kint) :: neib, nnb, mycnt

    ok = .true.
    nnb = cmt_c%n_neighbor
    if (nnb == 0) return
    allocate(sendv(nnb), recvv(nnb)); recvv = 0
    do neib = 1, nnb
      sendv(neib) = cmt_c%export_index(neib) - cmt_c%export_index(neib-1)
    end do
    call hecmw_saamg_comm_exchange_neighbor_int(cmt_c, sendv, recvv)
    do neib = 1, nnb
      mycnt = cmt_c%import_index(neib) - cmt_c%import_index(neib-1)
      if (recvv(neib) /= mycnt) ok = .false.   ! neighbor's export to me /= my import from it
    end do
    deallocate(sendv, recvv)
  end subroutine hecmw_saamg_verify_commtable

end module hecmw_precond_SAAMG_core
