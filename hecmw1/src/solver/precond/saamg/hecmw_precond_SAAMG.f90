!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : FrontISTR backend (id 22)
!!
!! Top-level preconditioner backend registered as PRECOND=22.  Mirrors the
!! setup / apply / clear lifecycle of the other backends (DIAG, ML, ...): the
!! hierarchy is held in module state and rebuilt only when the symbolic/numeric
!! setup flags (hecMAT%Iarray(98)/(97)) request it.  apply works in place on the
!! single working vector ZP (ZP <- M^{-1} ZP), one multigrid (V/W) cycle per call.
module hecmw_precond_SAAMG
  use hecmw_util
  use hecmw_matrix_misc, only: hecmw_mat_get_solver_opt, hecmw_mat_get_solver_ropt, &
       hecmw_mat_get_loglevel, hecmw_mat_get_timelog
  use hecmw_precond_SAAMG_adapt
  use hecmw_precond_SAAMG_comm,      only: hecmwST_saamg_comm, hecmw_saamg_comm_free, &
       hecmw_saamg_comm_allreduce_sum_r, hecmw_saamg_comm_allreduce_sum_int, &
       hecmw_saamg_check_alloc, hecmw_saamg_lapack_available
  use hecmw_precond_SAAMG_core,      only: hecmwST_saamg_hier, hecmw_saamg_setup, &
       hecmw_saamg_refresh, hecmw_saamg_apply, hecmw_saamg_free
  use hecmw_precond_SAAMG_coarse_mumps, only: hecmw_saamg_cmumps_available
  use hecmw_precond_SAAMG_aggregate, only: hecmw_saamg_write_vtk
  use hecmw_precond_SAAMG_matrix,    only: hecmwST_saamg_bcsr, hecmw_saamg_bcsr_free
  use hecmw_precond_SAAMG_param,     only: hecmwST_saamg_params
  use hecmw_precond_rbm,             only: hecmw_precond_rbm_from_mesh
  use hecmw_precond_SAAMG_verify,    only: hecmw_saamg_verify_matvec, hecmw_saamg_verify_coarsen, &
       hecmw_saamg_verify_prolong, hecmw_saamg_verify_smoothp, hecmw_saamg_verify_refresh
  use hecmw_solver_las,              only: hecmw_matvec
  use hecmw_mat_id
  implicit none

  private
  public :: hecmw_precond_SAAMG_setup
  public :: hecmw_precond_SAAMG_apply
  public :: hecmw_precond_SAAMG_clear

  type(hecmwST_saamg_hier), save :: g_dh    !< fully-distributed multilevel hierarchy
  type(hecmwST_saamg_comm),      save :: g_cmt   !< finest-level comm table (for MPI; empty on 1 rank)
  logical,            save :: INITIALIZED = .false.
  integer(kind=kint), save :: g_n = 0          !< finest internal dof count (= N*NDOF)
  integer(kind=kint), save :: g_id = 0         !< hecmw_mat_id handle for the fine-level matvec hook
  ! Size part of the structure signature of the matrix the hierarchy was built
  ! from.  The recycling flags (Iarray(98)/(97)) are a caller-side CONTRACT:
  ! structure changes must be reported via Iarray(98).  The one known caller whose
  ! structure can change under a values-only flag -- the contact-elimination path,
  ! whose reduced system depends on a value-driven slave-dof choice -- now detects
  ! this itself and raises Iarray(98) at the source (notify_structure_change in
  ! solve_LINEQ_contact_elim; full pattern hash there).  This O(1) size check is
  ! kept as a last-resort safety net against OTHER contract violations: a stale
  ! hierarchy applied to a matrix of different size corrupts memory (observed on
  ! mobile_case np=8 before the fixes).  On mismatch we fall through to a full
  ! rebuild (notice under LOGLEVEL>=1).
  integer(kind=kint), save :: g_sig(4) = 0     !< [N, NP, NPL, NPU] at build time

contains

  subroutine hecmw_precond_SAAMG_setup(hecMAT, hecMESH, sym)
    implicit none
    type(hecmwST_matrix),     intent(inout) :: hecMAT
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    integer(kind=kint),       intent(in)    :: sym
    type(hecmwST_saamg_bcsr)    :: A
    type(hecmwST_saamg_params) :: prm
    real(kind=kreal), allocatable :: B(:,:)
    integer(kind=kint) :: ndof, nint, iopt(10), m, astat, myrank, sig_bad
    real(kind=kreal)   :: ropt(10)
    logical            :: will_use_mumps

    ! SA-AMG's setup needs LAPACK (dense coarsest factorization, per-aggregate QR,
    ! Lanczos tridiagonal eigenvalue).  LAPACK is optional in FrontISTR, so fail
    ! early with a clear message rather than at an obscure link/runtime point.
    if (.not. hecmw_saamg_lapack_available()) then
      if (hecmw_comm_get_rank() == 0) write(*,'(a)') &
        '#### SA-AMG (PRECOND=22) requires a LAPACK-enabled build. Rebuild with '// &
        '--with-lapack (setup.sh) or -DWITH_LAPACK=ON (cmake), or select another preconditioner.'
      call hecmw_abort(hecmw_comm_get_comm())
    end if

    ! Recycle-flag branches below are guarded by the structure signature (g_sig):
    ! the flags are a caller-side contract that some callers break (see g_sig
    ! comment).  The verdict MUST be identical on every rank -- reuse / refresh /
    ! rebuild all execute different collective patterns, so a per-rank decision
    ! deadlocks (one rank's contact state can change while another's does not).
    ! Hence the local check is allreduced before branching.
    if (INITIALIZED .and. hecMAT%Iarray(98) == 0) then
      sig_bad = 0
      if (.not. sig_sizes_match(hecMAT)) sig_bad = 1
      call hecmw_saamg_comm_allreduce_sum_int(g_cmt, sig_bad)
      if (sig_bad > 0) then
        if (hecmw_mat_get_loglevel(hecMAT) >= 1 .and. hecmw_comm_get_rank() == 0) write(*,'(a)') &
          '#### SA-AMG: matrix structure changed under a values-only recycle flag -- full rebuild'
      else if (hecMAT%Iarray(97) == 0) then
        ! nothing changed -> reuse the hierarchy as-is (also covers precond recycling)
        call saamg_register_mat(hecMAT, hecMESH)   ! keep the matvec hook pointing at the current matrix
        return
      else
        ! structure unchanged, values changed -> numeric-only refresh (Newton reuse):
        ! reuse aggregation / coarse comm table / tentative P-hat, recompute the rest.
        call saamg_register_mat(hecMAT, hecMESH)   ! refresh the fine-level matvec hook to the new values
        call hecmw_saamg_from_hecmat(hecMAT, A, include_halo=.true.)
        call hecmw_saamg_refresh(g_dh, A)
        call hecmw_saamg_bcsr_free(A)
        hecMAT%Iarray(97) = 0
        return
      end if
    end if

    ! structure changed -> full rebuild
    if (INITIALIZED) call hecmw_precond_SAAMG_clear(hecMAT%NDOF)
    call saamg_register_mat(hecMAT, hecMESH)     ! register for the fine-level matvec hook (HEC-MW matvec)

    ndof = hecMAT%NDOF
    ! near-kernel size m by problem type: 1=heat(const), 2=plane(2D RBM, 3 modes),
    ! 3=solid / 6=shell (3D RBM, 6 modes).
    select case (ndof)
    case (1);    m = 1
    case (2);    m = 3
    case (3, 6); m = 6
    case default
      write(*,'(a)') '#### SA-AMG (PRECOND=22): only NDOF=1 (heat) / 2 (plane) / 3 (solid) / 6 (shell) supported'
      call hecmw_abort(hecmw_comm_get_comm())
    end select
    nint = hecMAT%N

    ! map control-file options; 0 = keep default.  Integer line (Iarray 41:50) and
    ! real line (Rarray 41:50).  Integer slots 1-7 mirror the ML (PRECOND=5) layout
    ! so ML users can reuse the same option line (see fstr_ctrl_common precond==22
    ! and hecmw_ML_wrapper.c, which reads opt[0..6] = slots 1-7); slots 8-10 are
    ! outside what ML reads and hold SA-AMG-specific knobs.  The real line is not
    ! read by ML at all, so it carries the remaining knobs (including a few integer
    ! ones, rounded on read).
    !   1 = coarsest solver  (ML CoarseSolver: 0=auto, 1=Smoother, 2=dense, 3=MUMPS)
    !   2 = smoother type     (ML SmootherType: 0/1=Chebyshev; 2/3 N/A -> Chebyshev)
    !   3 = cycle             (ML MGType: 0=default(=W)/1=V, 2=W; 3=FullV N/A -> W)
    !   4 = max levels        (ML MaxLevels)
    !   5 = RESERVED          (ML CoarsenScheme; SA-AMG is always uncoupled -> warn
    !                          and ignore.  Deliberately kept unused: an ML deck
    !                          carries 1..5 here, which used to be silently taken as
    !                          max_size and wrecked the hierarchy)
    !   6 = Chebyshev degree  (ML NumSweeps)
    !   7 = coarse size       (ML MaxCoarseSize)
    !   8 = max aggregate size (SA-AMG-specific: the aggressive-coarsening lever)
    !   9 = galerkin_lowmem   (Ac=P^T A P: 0 = default (2-stage everywhere, faster),
    !                          > 0 = low-memory (fuse the finest level only, same
    !                          peak as fully fused but faster; deep levels 2-stage))
    !  10 = RESERVED          (free slot for future SA-AMG-specific knobs)
    ! Real line: 1 = theta (strength threshold), 2 = cheb_alpha, 3 = safety,
    !            4 = taper_k (coarsening taper, 3-1: 0 = default 100, > 0 = use as K,
    !                         < 0 = disable the taper / legacy behavior)
    !            5 = agg_order (aggregation seed-scan ordering: 0 = default (BFS),
    !                         < 0 = natural node order / legacy, 1..4 = explicit mode
    !                         (1=bfs, 2=gid-hash, 3=mindeg, 4=maxdeg))
    !            6 = min aggregate size (0 = default 3)
    !            7 = verify (0 = off, > 0 = on) / 8 = dump_vtk (0 = off, > 0 = on)
    call hecmw_mat_get_solver_opt(hecMAT, iopt)
    myrank = hecmw_comm_get_rank()
    ! slot 1: coarsest solver.  The external encoding (ML CoarseSolver-compatible)
    ! and the internal prm%coarsest_solver share the SAME numbering, so no
    ! translation is needed: 0=auto, 1=smoother, 2=dense, 3=MUMPS.
    select case (iopt(1))
    case (0, 1, 2, 3) ; prm%coarsest_solver = iopt(1)
    case default ; if (myrank == 0) write(*,'(a,i0,a)') &
        '#### SA-AMG: invalid coarse solver ', iopt(1), ' (ignored) -- using auto'
    end select
    ! slot 2: smoother type -- SA-AMG only implements Chebyshev
    select case (iopt(2))
    case (0, 1) ; ! Chebyshev (only supported smoother)
    case default ; if (myrank == 0) write(*,'(a,i0,a)') &
        '#### SA-AMG: smoother type ', iopt(2), ' not supported -- using Chebyshev'
    end select
    ! slot 3: multigrid cycle order gamma (ML MGType-compatible: 1=V, 2=W).
    ! 0 = keep default = W (prm%ncycle=2): never slower than V and markedly fewer
    ! iterations on deep hierarchies; nlevel=2 degenerates to V (is_coarsest guard).
    select case (iopt(3))
    case (0)    ; ! default -> keep prm%ncycle (=2, W-cycle)
    case (1)    ; prm%ncycle = 1        ! V-cycle (explicit)
    case (2)    ; prm%ncycle = 2        ! W-cycle
    case (3)    ; prm%ncycle = 2 ; if (myrank == 0) write(*,'(a)') &
        '#### SA-AMG: Full-V cycle not supported -- using W-cycle'
    case default ; if (myrank == 0) write(*,'(a,i0,a)') &
        '#### SA-AMG: invalid cycle ', iopt(3), ' (ignored) -- using default (W)'
    end select
    if (iopt(4) > 0) prm%max_level = iopt(4)   ! ML MaxLevels
    ! slot 5: reserved.  ML puts CoarsenScheme here (1..5 = UncoupledMIS/METIS/
    ! ParMETIS/Zoltan/DD); SA-AMG always uses uncoupled aggregation, so the value is
    ! ignored -- but say so, otherwise a reused ML option line looks like it selects
    ! a coarsening scheme.
    if (iopt(5) /= 0 .and. myrank == 0) then
      write(*,'(a,i0,a)') '#### SA-AMG: int slot5 (ML CoarsenScheme) = ', iopt(5), &
           ' is ignored -- SA-AMG always uses uncoupled aggregation'
      write(*,'(a)')      '####         (max aggregate size is int slot8)'
    end if
    if (iopt(6) > 0) prm%cheb_deg  = iopt(6)   ! ML NumSweeps (Chebyshev degree)
    if (iopt(7) > 0) prm%coarse_size = iopt(7) ! ML MaxCoarseSize
    if (iopt(8) > 0) prm%max_size  = iopt(8)   ! SA-AMG-specific (aggressive coarsening lever)
    if (iopt(9) > 0) prm%galerkin_lowmem = .true.  ! SA-AMG-specific (fuse the finest level only)
    ! verbose has no iopt slot: enable it via !SOLVER LOGLEVEL>=1 (independent of
    ! TIMELOG; LOGLEVEL unset = -1 leaves verbose off).
    prm%loglevel = hecmw_mat_get_loglevel(hecMAT)   ! independent LOGLEVEL (-1 if unset); >=2 -> [mem] probes
    prm%verbose  = (prm%loglevel >= 1)
    prm%timelog  = hecmw_mat_get_timelog(hecMAT)    ! !SOLVER TIMELOG (2=VERBOSE -> per-section [time] probes)
    call hecmw_mat_get_solver_ropt(hecMAT, ropt)
    if (ropt(1) > 0.0d0) prm%theta      = ropt(1)
    if (ropt(2) > 0.0d0) prm%cheb_alpha = ropt(2)
    if (ropt(3) > 0.0d0) prm%safety     = ropt(3)
    ! real slot 4: coarsening-taper K (an integer carried on the real line: the int
    ! line is reserved for the ML-compatible and high-traffic knobs, and the taper is
    ! rarely touched).  0 = default (prm%taper_k = 100), > 0 = use as K,
    ! < 0 = disable the taper (legacy pre-taper coarsening).
    ! (int(x+0.5), not the intrinsic nint(): a local variable `nint` shadows it.)
    if (ropt(4) > 0.0d0) prm%taper_k = int(ropt(4) + 0.5d0, kind=kint)
    if (ropt(4) < 0.0d0) prm%taper_k = 0
    ! real slot 5: aggregation seed-scan ordering.  0 = default (BFS), < 0 = natural
    ! (legacy), 1..4 = explicit mode.
    if (ropt(5) > 0.0d0 .and. ropt(5) < 4.5d0) prm%agg_order = int(ropt(5) + 0.5d0, kind=kint)
    if (ropt(5) < 0.0d0) prm%agg_order = 0
    ! real slot 6: min aggregate size (integer carried on the real line; rarely
    ! touched -- it only sets the forced-merge threshold for leftover aggregates).
    if (ropt(6) > 0.0d0) prm%min_size = int(ropt(6) + 0.5d0, kind=kint)
    ! real slots 7/8: diagnostics (0 = off, > 0 = on), alongside LOGLEVEL/TIMELOG.
    prm%verify   = (ropt(7) > 0.0d0)
    prm%dump_vtk = (ropt(8) > 0.0d0)   ! write finest aggregation to saamg_agg.<rank>.vtk
    ! symmetric (CG, sym=1) vs non-symmetric (BiCGSTAB/GMRES, sym=0): the latter uses a
    ! general (SYM=0 / LU) coarsest so the actual non-symmetric coarse is factored exactly.
    prm%symmetric = (sym == 1)

    ! coarse_size auto-default by coarsest solver (only when the user left it unset):
    ! dense LDL^T is O(N^3) so keep the coarsest small (100 dof); a distributed sparse
    ! MUMPS coarsest is cheap at tens of thousands of dof -> use a much larger default
    ! so the hierarchy stays shallow (matches the ML+MUMPS convention of 50000).
    if (iopt(7) == 0) then
      will_use_mumps = (prm%coarsest_solver == 3) .or. &
           (prm%coarsest_solver == 0 .and. hecmw_saamg_cmumps_available())
      if (will_use_mumps) prm%coarse_size = 50000
    end if

    ! finest-level communication table (empty on 1 rank -> sequential behavior)
    call hecmw_saamg_comm_from_mesh(hecMESH, ndof, g_cmt)
    ! F4a: verify the distributed (halo-aware) matvec is consistent across ranks
    if (prm%verify) call hecmw_saamg_verify_matvec(hecMAT, hecMESH, g_cmt, ndof)
    ! F4b: verify uncoupled aggregation + coarse comm-table construction
    if (prm%verify) call hecmw_saamg_verify_coarsen(hecMAT, g_cmt, m, prm)
    ! F4b-2: verify halo-extended tentative P-hat (P-hat * B_coarse = B_fine, all rows)
    if (prm%verify) call hecmw_saamg_verify_prolong(hecMAT, hecMESH, g_cmt, ndof, m, prm)
    ! F4c-1: verify the distributed smoothed prolongator (global ||P B_c - B_f|| rank-independent)
    if (prm%verify) call hecmw_saamg_verify_smoothp(hecMAT, hecMESH, g_cmt, ndof, m, prm)

    ! fully-distributed multilevel hierarchy: finest operator keeps halo columns
    call hecmw_saamg_from_hecmat(hecMAT, A, include_halo=.true.)
    allocate(B(nint*ndof, m), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'near-kernel B (rigid-body modes)')
    ! shared rigid-body near-kernel (with ML)
    call hecmw_precond_rbm_from_mesh(hecMESH, ndof, B)
    ! move A into the hierarchy (no redundant copy) unless verify still needs it
    call hecmw_saamg_setup(A, g_cmt, B, m, prm, g_dh, saamg_fine_matvec, &
         move_in=(.not. prm%verify))
    ! verify: a numeric refresh reproduces (a) the full build for the SAME A and
    ! (b) a fresh full build for a value-perturbed A' (same pattern) -- S4 routing
    if (prm%verify) call hecmw_saamg_verify_refresh(g_dh, A, g_cmt, B, m, prm, nint*ndof)
    if (prm%dump_vtk) call dump_agg_vtk(g_dh, hecMESH, ndof, nint)
    call hecmw_saamg_bcsr_free(A)
    deallocate(B)

    g_n = nint*ndof
    g_sig = (/ hecMAT%N, hecMAT%NP, hecMAT%NPL, hecMAT%NPU /)
    INITIALIZED = .true.
    hecMAT%Iarray(98) = 0
    hecMAT%Iarray(97) = 0
  end subroutine hecmw_precond_SAAMG_setup

  !> O(1) part of the structure signature: matrix/halo sizes.
  logical function sig_sizes_match(hecMAT) result(ok)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    ok = (g_sig(1) == hecMAT%N .and. g_sig(2) == hecMAT%NP .and. &
          g_sig(3) == hecMAT%NPL .and. g_sig(4) == hecMAT%NPU)
  end function sig_sizes_match


  !> Dump the finest-level aggregation to a per-rank VTK point cloud (ParaView).
  !! Each rank writes its internal nodes (coords from the mesh) colored by the GLOBAL
  !! coarse-node id, so aggregate colors are consistent across the partition.  Loaded
  !! together in ParaView the pieces tile the whole domain (saamg_agg.<rank>.vtk).
  subroutine dump_agg_vtk(dh, hecMESH, ndof, nint)
    implicit none
    type(hecmwST_saamg_hier), intent(in) :: dh
    type(hecmwST_local_mesh),    intent(in) :: hecMESH
    integer(kind=kint),          intent(in) :: ndof, nint
    character(len=128) :: fname
    integer(kind=kint) :: rank, idummy
    idummy = ndof
    if (.not. allocated(dh%aggr_fine)) return
    rank = hecmw_comm_get_rank()
    write(fname,'(a,i0,a)') 'saamg_agg.', rank, '.vtk'
    ! hecMESH%node is 3*n_node interleaved; internal nodes are 1..nint
    call hecmw_saamg_write_vtk(trim(fname), hecMESH%node, nint, dh%aggr_fine)
    if (rank == 0) write(*,'(a)') &
         '#### SA-AMG: wrote finest aggregation to saamg_agg.<rank>.vtk (ParaView)'
  end subroutine dump_agg_vtk

  subroutine hecmw_precond_SAAMG_apply(ZP, NDOF)
    implicit none
    real(kind=kreal),   intent(inout) :: ZP(:)
    integer(kind=kint), intent(in)    :: NDOF
    real(kind=kreal), allocatable :: r(:), z(:)
    allocate(r(g_n), z(g_n))
    r(1:g_n) = ZP(1:g_n)
    z(1:g_n) = 0.0d0
    call hecmw_saamg_apply(g_dh, r, z)   ! z = M^{-1} r (one distributed V-cycle)
    ZP(1:g_n) = z(1:g_n)
    deallocate(r, z)
  end subroutine hecmw_precond_SAAMG_apply

  subroutine hecmw_precond_SAAMG_clear(NDOF)
    implicit none
    integer(kind=kint), intent(in) :: NDOF
    if (INITIALIZED) then
      call hecmw_saamg_free(g_dh)
      call hecmw_saamg_comm_free(g_cmt)
      INITIALIZED = .false.
      g_n = 0
      g_sig = 0
    end if
    if (g_id > 0) then; call hecmw_mat_id_clear(g_id); g_id = 0; end if
  end subroutine hecmw_precond_SAAMG_clear

  !> (Re)register the current (hecMAT, hecMESH) in the hecmw_mat_id registry so the
  !! fine-level matvec hook can look them up at apply time.  Releases any previous
  !! handle first (no leak across Newton refresh / recycling).
  subroutine saamg_register_mat(hecMAT, hecMESH)
    implicit none
    type(hecmwST_matrix),     intent(in), target :: hecMAT
    type(hecmwST_local_mesh), intent(in), target :: hecMESH
    if (g_id > 0) call hecmw_mat_id_clear(g_id)
    call hecmw_mat_id_set(hecMAT, hecMESH, g_id)
  end subroutine saamg_register_mat

  !> Fine-level matvec hook bound into the distributed V-cycle: y = A x via HEC-MW's
  !! native blocked matvec on the registered hecMAT (the finest SA-AMG operator IS
  !! hecMAT).  Avoids the scalar-CSR matvec in the hot Chebyshev smoother.  x carries
  !! internal+halo (halo refreshed inside hecmw_matvec); y returns internal rows.
  subroutine saamg_fine_matvec(x, y)
    implicit none
    real(kind=kreal), intent(inout) :: x(:)
    real(kind=kreal), intent(out)   :: y(:)
    type(hecmwST_matrix),     pointer :: mat
    type(hecmwST_local_mesh), pointer :: mesh
    call hecmw_mat_id_get(g_id, mat, mesh)
    call hecmw_matvec(mesh, mat, x, y)
  end subroutine saamg_fine_matvec

end module hecmw_precond_SAAMG
