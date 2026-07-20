!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : tunable parameters
!!
!! Low-level module holding only the parameter struct, so the distributed core
!! (hecmw_precond_SAAMG_core) and the backend can depend on it from a low level
!! without a cyclic dependency.
module hecmw_precond_SAAMG_param
  use hecmw_precond_SAAMG_util
  implicit none

  private
  public :: hecmwST_saamg_params

  !> Tunable parameters (defaults match the implementation plan, ch.4).
  type hecmwST_saamg_params
    integer(kind=kint) :: cheb_deg   = 2      !< Chebyshev degree
    real(kind=kreal)   :: cheb_alpha = 20.0d0 !< eigenvalue ratio (interval lower end)
    integer(kind=kint) :: min_size   = 3      !< min aggregate size
    integer(kind=kint) :: max_size   = 96     !< max aggregate size (phase 2).  With block-CSR setup
                                              !! and Lanczos lambda_max, aggressive coarsening lowers
                                              !! BOTH setup and solve on large meshes (hinge 252k: -16%
                                              !! total; Mold 3M: 250->180s, iters 205->194) and is
                                              !! harmless on small/contact (coarse saturates).  128 is
                                              !! marginally faster on very large meshes; 96 keeps a
                                              !! safety margin (~64 nodes/agg vs ML UncoupledMIS ~35).
    real(kind=kreal)   :: theta      = 0.0d0  !< strength-of-connection threshold (0 = structural graph)
    integer(kind=kint) :: taper_k    = 100    !< coarsening taper: at levels >= 2 cap the effective
                                              !! max_size at (global nodes)/taper_k, i.e. keep at least
                                              !! ~taper_k aggregates per level so the deep levels shrink
                                              !! gradually instead of collapsing into a few aggregates
                                              !! (weak coarse correction, worst at low rank counts --
                                              !! uncoupled aggregation floors the coarsest at ~1 aggregate
                                              !! per rank).  Measured: Mold 3M np=1 W-cycle 740 -> 169
                                              !! iters, total 3.7x; high np neutral.  0 = off (legacy).
    integer(kind=kint) :: agg_order  = 1      !< phase-1 seed-scan ordering of the aggregation
                                              !! (the seeds form a distance-2 MIS; the ordering decides
                                              !! WHICH admissible seed is taken first and thus the
                                              !! aggregate shapes): 0 = natural node order (legacy),
                                              !! 1 = BFS/graph order (default; contiguous tiling ->
                                              !! uniform aggregates.  Measured: never worse than natural
                                              !! except Mold np=1/4, and much better on badly numbered
                                              !! meshes: 3M np=8 W 284->162, 33M 230->172, v6 36->28
                                              !! beating ML, contact total iters -36%),
                                              !! 2 = gid-hash, 3 = min-degree, 4 = max-degree (experimental)
    logical            :: galerkin_lowmem = .false. !< Galerkin Ac = P^T A P memory/speed mode.
                                              !! .false. (default) = 2-stage at EVERY level: materialize
                                              !!   C = A*Pext then P^T*C.  Clean SpGEMM kernels + reuse ->
                                              !!   fastest setup (Mold 3M np=8: 10.4 -> 5.0 s, iters
                                              !!   unchanged), but C is the setup memory high-water mark
                                              !!   (level 1: +1.4 GB peak on 3M).
                                              !! .true. (low-memory) = FUSE THE FINEST LEVEL ONLY (l==1),
                                              !!   2-stage on the deep levels.  The finest level's C=A*P is
                                              !!   the memory peak; fusing just that level removes the peak
                                              !!   (same peak as fusing all levels -- deep-level C is ~1/30
                                              !!   and never the high-water mark), while the cheap deep
                                              !!   levels keep the fast 2-stage kernels.  Use when
                                              !!   memory-bound at scale.  (Fusing ALL levels is strictly
                                              !!   worse: same peak, slower -- so it is not offered.)
    integer(kind=kint) :: lanczos_iter = 10   !< Lanczos iterations for lambda_max(D^{-1}A)
    real(kind=kreal)   :: safety     = 1.1d0  !< eigenvalue safety factor
    integer(kind=kint) :: max_level  = 20     !< hierarchy depth cap
    integer(kind=kint) :: coarse_size = 100   !< coarsest size threshold (stop coarsening)
    integer(kind=kint) :: coarsest_solver = 0 !< 0=auto (MUMPS if built else dense), 1=smoother sweeps, 2=dense, 3=MUMPS (same encoding as !SOLVER slot1, ML CoarseSolver-compatible)
    integer(kind=kint) :: ncycle     = 2      !< multigrid cycle order gamma per level (1=V, 2=W); default W (never slower than V; far fewer iters on deep hierarchies). nlevel=2 -> W==V (is_coarsest guard)
    real(kind=kreal)   :: stagn      = 0.7d0  !< coarsening stagnation guard (n_c/n_f)
    logical            :: symmetric  = .true. !< matrix assumed symmetric (CG); .false. => non-sym coarsest (SYM=0/LU)
    logical            :: verbose    = .false.!< print hierarchy diagnostics/timing at setup (LOGLEVEL>=1)
    integer(kind=kint) :: loglevel   = -1     !< diagnostic verbosity from !SOLVER LOGLEVEL (-1=unset); >=2 adds [mem] probes
    integer(kind=kint) :: timelog    = 0      !< !SOLVER TIMELOG (0=NO,1=YES,2=VERBOSE); ==2 emits per-section [time] probes
    logical            :: verify     = .false.!< run per-level self-checks at setup
    logical            :: dump_vtk   = .false.!< dump finest aggregation to VTK (ParaView)
  end type hecmwST_saamg_params

end module hecmw_precond_SAAMG_param
