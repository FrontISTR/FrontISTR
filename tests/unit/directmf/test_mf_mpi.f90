!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> Unit checks of the MPI distributed factorization and solve of hecmw_mf_numeric
!> (subtree-to-subcube mapping, fan-in upper fronts). Run under mpirun with any process
!> count; every rank builds the identical replicated matrix, the symbolic stage runs on
!> rank 0 and is broadcast, mirroring the driver flow. The distributed factorization is
!> fed the partial matrix of hecmw_mf_dist_gmat_extract, so the checks cover the partial
!> holding the driver path builds. Links against the hecmw library and LAPACK; see
!> CMakeLists.txt.
!>
!> For each case:
!>   (a) the relative residual of the serial run is below 1e-8 and the one of the distributed
!>       run matches the single process level
!>   (b) a repeated distributed factorization and solve is bitwise reproducible
!>   (c) the returned solution is bitwise identical on every rank
!>   (d) the per-rank factor words sum to the serial factor words
!>   (e) the solution is compared bitwise with the serial (np=1 equivalent) run on the same
!>       structure: required at np=1, reported as bitwise_serial=T/F at np>1 (the cross-np
!>       agreement is a design goal, checked through this common reference)
!>   (g) the stored panels of the upper fronts sum to the serial ones and, when a front has
!>       contribution rows, are held by at least two ranks (the 1D row distribution works)
!>   (f) with BLR, iterative refinement with the distributed factor reaches 1e-10
!>   (h) the zerodiag cases (every diagonal block zero, so delayed pivots climb to the root
!>       and the distributed exchanges with delayed columns and the root final panel are
!>       exercised) fail unless delayed pivots occurred
!>   (i) the per-rank inertia counts sum to the serial ones
!> and separately:
!>   (r1), (r2) the BLR rank reuse (run_reuse_case)
!>   (j) the BLR compression of a tile holding a NaN or an Inf gives up
!> A process grid with two columns (hecmw_mf_dist_fsgrid) needs 4 or more processes and a front
!> of 8 or more fully summed tiles; the zd_grid2d_d1_t2, zd_grid3d_d1_t3 and sad_grid2d_d3_t3 cases and
!> their LU counterparts update tiles on it.
program test_mf_mpi
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_mf_dist
  use hecmw_mf_numeric
  use hecmw_mf_kernel, only: hecmw_mf_kernel_compress
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  implicit none

  integer(kind=kint) :: nfail, ntest, me, np, ierr
  integer(kind=8) :: seed
  real(kind=kreal), parameter :: irtol = 1.0d-10

  call MPI_Init(ierr)
  call hecmw_init_nocontrol
  me = hecmw_comm_get_rank()
  np = hecmw_comm_get_size()
  nfail = 0
  ntest = 0
  seed = 31415

  call run_case('spd_grid3d_d3_metis_t16',   'randspd',   3, 7, 6, 5, 3, 5, 16, 0, .false., 0.0d0)
  call run_case('spd_grid3d_d3_metis_t64',   'randspd',   3, 9, 8, 7, 3, 5, 64, 0, .false., 0.0d0)
  call run_case('sad_grid2d_d3_relax_t16',   'smoothsad', 2, 12, 10, 1, 3, 5, 16, 2, .false., 0.0d0)
  call run_case('lu_grid2d_d3_metis_t16',    'randoff',   2, 12, 10, 1, 3, 0, 16, 0, .false., 0.0d0)
  call run_case('blr_grid3d_d3_metis_t16',   'smoothspd', 3, 7, 6, 5, 3, 5, 16, 0, .true., 1.0d-8)
  call run_case('blrlu_grid3d_d3_metis_t16', 'smoothlu',  3, 7, 6, 5, 3, 5, 16, 0, .true., 1.0d-8)
  call run_case('zd_grid2d_d1_t8',    'zerodiag',   2, 16, 14, 1, 1, 5, 8, 0, .false., 0.0d0)
  call run_case('zd_grid2d_d3_t16',   'zerodiag',   2, 12, 10, 1, 3, 5, 16, 0, .false., 0.0d0)
  call run_case('zdlu_grid2d_d1_t8',  'zerodiaglu', 2, 16, 14, 1, 1, 5, 8, 0, .false., 0.0d0)
  call run_case('zdlu_grid2d_d3_t16', 'zerodiaglu', 2, 12, 10, 1, 3, 5, 16, 0, .false., 0.0d0)
  call run_case('zd_grid2d_d1_t2',    'zerodiag',   2, 16, 14, 1, 1, 5, 2, 0, .false., 0.0d0)
  call run_case('zd_grid2d_d1_t3',    'zerodiag',   2, 16, 14, 1, 1, 5, 3, 0, .false., 0.0d0)
  call run_case('zd_grid2d_d3_t3',    'zerodiag',   2, 12, 10, 1, 3, 5, 3, 0, .false., 0.0d0)
  call run_case('zdlu_grid2d_d3_t3',  'zerodiaglu', 2, 12, 10, 1, 3, 5, 3, 0, .false., 0.0d0)
  call run_case('zd_grid2d_d3_t5',    'zerodiag',   2, 12, 10, 1, 3, 5, 5, 0, .false., 0.0d0)
  call run_case('zd_grid2d_d3_t7',    'zerodiag',   2, 12, 10, 1, 3, 5, 7, 0, .false., 0.0d0)
  call run_case('zd_grid3d_d1_t3',    'zerodiag',   3, 7, 6, 5, 1, 5, 3, 0, .false., 0.0d0)
  call run_case('zd_grid3d_d3_t5',    'zerodiag',   3, 7, 6, 5, 3, 5, 5, 0, .false., 0.0d0)
  call run_case('sad_grid2d_d3_t3',   'smoothsad',  2, 12, 10, 1, 3, 5, 3, 2, .false., 0.0d0)
  call run_case('sad_grid2d_d3_t5',   'smoothsad',  2, 12, 10, 1, 3, 5, 5, 2, .false., 0.0d0)
  call run_case('zdblr_grid2d_d3_t3', 'zerodiag',   2, 12, 10, 1, 3, 5, 3, 0, .true., 1.0d-8)
  call run_case('zdblrlu_grid2d_d3_t3', 'zerodiaglu', 2, 12, 10, 1, 3, 5, 3, 0, .true., 1.0d-8)
  call run_case('zdblr_grid3d_d3_t5', 'zerodiag',   3, 7, 6, 5, 3, 5, 5, 0, .true., 1.0d-8)
  call run_case('zdlu_grid2d_d1_t2',  'zerodiaglu', 2, 16, 14, 1, 1, 5, 2, 0, .false., 0.0d0)
  call run_case('zdlu_grid3d_d1_t3',  'zerodiaglu', 3, 7, 6, 5, 1, 5, 3, 0, .false., 0.0d0)
  call run_case('nsd_grid2d_d3_t3',   'smoothnsd',  2, 12, 10, 1, 3, 5, 3, 2, .false., 0.0d0)
  call run_reuse_case('reuse_grid3d_d3_t16', 3, 7, 6, 5, 3, 5, 16)
  call run_compress_nonfinite('compress_nonfinite')

  if (me == 0) write(*,'(a,i0,a,i0,a,i0)') 'test_mf_mpi: np=', np, ' failed ', nfail, ' of ', ntest
  call MPI_Finalize(ierr)
  if (nfail > 0) stop 1

contains

  subroutine run_case(name, kind, dim, nx, ny, nz, ndof, relax, tile, mstride, blr, eps)
    character(len=*), intent(in) :: name, kind
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, ndof, relax, tile, mstride
    logical, intent(in) :: blr
    real(kind=kreal), intent(in) :: eps
    type(hecmwST_mf_graph) :: graph
    type(hecmwST_mf_symbolic) :: sym
    type(hecmwST_mf_map) :: map
    type(hecmwST_mf_factor) :: fct0, fct
    type(hecmwST_matrix) :: hecMAT, pmat
    integer(kind=kint), allocatable :: xadj(:), adjncy(:), perm(:), invp(:), nd(:)
    logical, allocatable :: ismult(:)
    real(kind=kreal), allocatable :: r(:), x0(:), x1(:), x2(:), dd(:)
    integer(kind=kint) :: n, it, ss, iw(1), iref(2), iw2(2)
    integer(kind=8) :: wref, wupref, wup
    logical :: hascb
    real(kind=kreal) :: resid0, resid, rir, bnrm, limit, w1(1)
    logical :: ok, bser
    character(len=8) :: sbser

    ntest = ntest + 1
    ok = .true.
    call make_graph(dim, nx, ny, nz, mstride, n, xadj, adjncy, ismult)
    allocate(nd(n), perm(n), invp(n))
    nd(1:n) = ndof
    call hecmw_mf_graph_set(graph, n, nd, xadj, adjncy)
    call make_matrix(n, ndof, xadj, adjncy, ismult, kind, hecMAT)

    ! symbolic on rank 0, broadcast (the driver flow)
    if (me == 0) then
      call hecmw_ordering_gen(n, n + (xadj(n+1)-1)/2, xadj, adjncy, perm, invp, 0, 0)
      call hecmw_mf_symbolic_build(graph, perm, relax, sym)
      if (blr) then
        call hecmw_mf_symbolic_cluster(graph, sym, perm)
        call hecmw_mf_symbolic_finalize(sym)
        call hecmw_mf_symbolic_build(graph, perm, relax, sym)
      endif
    endif
    call hecmw_mf_dist_symbolic_bcast(sym, 0)
    call hecmw_mf_dist_map_build(sym, map)
    call hecmw_mf_dist_gmat_extract(hecMAT, sym, map, pmat)

    allocate(x0(n*ndof), x1(n*ndof), x2(n*ndof), dd(n*ndof), r(n*ndof))
    bnrm = sqrt(dot_product(hecMAT%B, hecMAT%B))

    ! serial reference, run identically on every rank (equals the np=1 result by design)
    call hecmw_mf_numeric_init(graph, sym, tile, blr, merge(MF_BLR_ETA, 0, blr), fct0)
    fct0%blr = blr
    fct0%eps = eps
    call hecmw_mf_numeric_factor(hecMAT, sym, fct0, ierr)
    if (.not. reduced(ierr == 0)) then
      if (me == 0) write(*,'(a,a,a,i0)') '  ', name, ': serial factorization failed, ierr = ', ierr
      ok = .false.
    endif
    resid0 = 0.0d0
    wref = 0
    if (ok) then
      call hecmw_mf_numeric_solve(sym, fct0, hecMAT%B, x0)
      call residual(hecMAT, x0, r)
      resid0 = sqrt(dot_product(r, r)) / bnrm
      ! The limits of (a) follow this residual, so a wrong serial solve must not loosen them.
      if (.not. reduced(resid0 <= 1.0d-8)) then
        w1(1) = resid0
        call hecmw_allreduce_R_comm(w1, 1, hecmw_max, hecmw_comm_get_comm())
        if (me == 0) write(*,'(a,a,a,1pe12.4)') '  ', name, ': relative residual of the serial run = ', w1(1)
        ok = .false.
      endif
      wref = fct0%factor_words_act
      iref(1) = fct0%n_pos
      iref(2) = fct0%n_neg
      wupref = 0
      do it = 1, map%nupper
        wupref = wupref + hecmw_mf_numeric_front_words(fct0, map%uplist(it))
      enddo
    endif
    call hecmw_mf_numeric_finalize(fct0)

    ! distributed factorization and solve
    call hecmw_mf_numeric_init(graph, sym, tile, blr, merge(MF_BLR_ETA, 0, blr), fct)
    fct%blr = blr
    fct%eps = eps
    resid = 0.0d0
    if (ok) then
      call hecmw_mf_numeric_factor_mpi(pmat, sym, map, fct, ierr)
      if (ierr /= 0) then
        if (me == 0) write(*,'(a,a,a,i0)') '  ', name, ': distributed factorization failed, ierr = ', ierr
        ok = .false.
      endif
    endif
    if (ok) then
      call hecmw_mf_numeric_solve_mpi(sym, map, fct, hecMAT%B, x1)
      ! (a) residual at the single process level
      call residual(hecMAT, x1, r)
      resid = sqrt(dot_product(r, r)) / bnrm
      limit = 1.0d-9
      if (blr) limit = max(100.0d0*eps, 10.0d0*resid0)
      if (kind(1:8) == 'zerodiag') limit = max(1.0d-8, 100.0d0*resid0)
      if (.not. reduced(resid <= limit)) then
        w1(1) = resid
        call hecmw_allreduce_R_comm(w1, 1, hecmw_max, hecmw_comm_get_comm())
        if (me == 0) write(*,'(a,a,a,1pe12.4)') '  ', name, ': relative residual = ', w1(1)
        ok = .false.
      endif
    endif
    if (ok) then
      ! (b) repeated distributed run is bitwise reproducible
      call hecmw_mf_numeric_factor_mpi(pmat, sym, map, fct, ierr)
      call hecmw_mf_numeric_solve_mpi(sym, map, fct, hecMAT%B, x2)
      if (.not. reduced(all(x2 == x1) .and. ierr == 0)) then
        if (me == 0) write(*,'(a,a,a)') '  ', name, ': repeated distributed run is not bitwise reproducible'
        ok = .false.
      endif
    endif
    if (ok) then
      ! (c) the solution is replicated bitwise across the ranks
      x2(1:n*ndof) = x1(1:n*ndof)
      call hecmw_allreduce_R_comm(x2, n*ndof, hecmw_max, hecmw_comm_get_comm())
      if (.not. reduced(all(x2 == x1))) then
        if (me == 0) write(*,'(a,a,a)') '  ', name, ': solution differs between ranks'
        ok = .false.
      endif
    endif
    if (ok) then
      ! (d) the per-rank factor words sum to the serial factor words
      w1(1) = real(fct%factor_words_act, kind=kreal)
      call hecmw_allreduce_R_comm(w1, 1, hecmw_sum, hecmw_comm_get_comm())
      if (nint(w1(1), kind=8) /= wref) then
        if (me == 0) write(*,'(a,a,a,i0,a,i0)') '  ', name, ': factor words sum ', nint(w1(1), kind=8), &
          ' vs serial ', wref
        ok = .false.
      endif
    endif
    if (ok) then
      ! (i) the per-rank inertia counts sum to the serial ones
      iw2(1) = fct%n_pos
      iw2(2) = fct%n_neg
      call hecmw_allreduce_I_comm(iw2, 2, hecmw_sum, hecmw_comm_get_comm())
      if (any(iw2 /= iref)) then
        if (me == 0) write(*,'(a,a,a,i0,a,i0,a,i0,a,i0)') '  ', name, ': inertia (+,-) = ', iw2(1), ',', iw2(2), &
          ' vs serial ', iref(1), ',', iref(2)
        ok = .false.
      endif
    endif
    if (ok) then
      ! (g) the upper front panels sum to the serial ones and are spread over the ranks
      wup = 0
      hascb = .false.
      do it = 1, map%nupper
        ss = map%uplist(it)
        wup = wup + hecmw_mf_numeric_front_words(fct, ss)
        if (sym%rptr(ss+1) - sym%rptr(ss) > sym%sptr(ss+1) - sym%sptr(ss)) hascb = .true.
      enddo
      w1(1) = real(wup, kind=kreal)
      call hecmw_allreduce_R_comm(w1, 1, hecmw_sum, hecmw_comm_get_comm())
      iw(1) = merge(1, 0, wup > 0)
      call hecmw_allreduce_I_comm(iw, 1, hecmw_sum, hecmw_comm_get_comm())
      if (nint(w1(1), kind=8) /= wupref .or. (np > 1 .and. hascb .and. iw(1) < 2)) then
        if (me == 0) write(*,'(a,a,a,i0,a,i0,a,i0)') '  ', name, ': upper front words ', nint(w1(1), kind=8), &
          ' vs serial ', wupref, ', holding ranks ', iw(1)
        ok = .false.
      endif
    endif
    if (ok .and. kind(1:8) == 'zerodiag') then
      ! (h) the case must exercise delayed pivots
      iw(1) = fct%n_delay
      call hecmw_allreduce_I_comm(iw, 1, hecmw_max, hecmw_comm_get_comm())
      if (iw(1) == 0) then
        if (me == 0) write(*,'(a,a,a)') '  ', name, ': no delayed pivots exercised'
        ok = .false.
      endif
    endif
    bser = .false.
    if (ok) then
      ! (e) bitwise agreement with the serial reference; required at np=1
      bser = reduced(all(x1 == x0))
      if (np == 1 .and. .not. bser) then
        if (me == 0) write(*,'(a,a,a)') '  ', name, ': np=1 run is not bitwise the serial one'
        ok = .false.
      endif
    endif
    rir = 0.0d0
    if (ok .and. blr) then
      ! (f) iterative refinement with the distributed factor
      x2(1:n*ndof) = x1(1:n*ndof)
      do it = 1, 40
        call residual(hecMAT, x2, r)
        rir = sqrt(dot_product(r, r)) / bnrm
        if (reduced(rir <= irtol)) exit
        call hecmw_mf_numeric_solve_mpi(sym, map, fct, r, dd)
        x2 = x2 + dd
      enddo
      if (.not. reduced(rir <= irtol)) then
        if (me == 0) write(*,'(a,a,a,1pe12.4)') '  ', name, ': refinement stalled at residual = ', rir
        ok = .false.
      endif
    endif

    iw(1) = merge(0, 1, ok)
    call hecmw_allreduce_I_comm(iw, 1, hecmw_max, hecmw_comm_get_comm())
    ok = iw(1) == 0
    if (me == 0) then
      sbser = 'F'
      if (bser) sbser = 'T'
      if (ok) then
        write(*,'(a,a,a,i0,a,i0,a,i0,a,i0,a,1pe10.2,a,1pe10.2,a,a)') 'PASS ', name, ': ndof_tot=', n*ndof, &
          ' nsuper=', sym%nsuper, ' upper=', map%nupper, ' delayed=', fct%n_delay, &
          ' resid=', resid, ' serial_resid=', resid0, ' bitwise_serial=', trim(sbser)
      else
        write(*,'(a,a)') 'FAIL ', name
      endif
    endif
    if (.not. ok) nfail = nfail + 1

    call hecmw_mf_numeric_finalize(fct)
    call hecmw_mf_dist_map_finalize(map)
    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    deallocate(hecMAT%D, hecMAT%AL, hecMAT%AU, hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU)
    deallocate(hecMAT%B, hecMAT%X)
    deallocate(pmat%D, pmat%AL, pmat%AU, pmat%indexL, pmat%indexU, pmat%itemL, pmat%itemU)
    deallocate(xadj, adjncy, perm, invp, nd, ismult)
    deallocate(x0, x1, x2, dd, r)
  end subroutine run_case

  !> BLR rank reuse: (r1) a sequential LU factorization with the reuse after an LDLt one on the
  !> same layout, (r2) a distributed BLR factorization with the reuse after one without BLR
  !> keeps the low rank tiles of the first BLR factorization
  subroutine run_reuse_case(name, dim, nx, ny, nz, ndof, relax, tile)
    character(len=*), intent(in) :: name
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, ndof, relax, tile
    type(hecmwST_mf_graph) :: graph
    type(hecmwST_mf_symbolic) :: sym
    type(hecmwST_mf_map) :: map
    type(hecmwST_mf_factor) :: fct
    type(hecmwST_matrix) :: hecMAT, pmat
    integer(kind=kint), allocatable :: xadj(:), adjncy(:), perm(:), invp(:), nd(:)
    logical, allocatable :: ismult(:)
    real(kind=kreal), allocatable :: x(:), r(:)
    integer(kind=kint) :: n, iw(2)
    real(kind=kreal) :: bnrm, resid
    logical :: ok

    ntest = ntest + 1
    ok = .true.
    call make_graph(dim, nx, ny, nz, 0, n, xadj, adjncy, ismult)
    allocate(nd(n), perm(n), invp(n))
    nd(1:n) = ndof
    call hecmw_mf_graph_set(graph, n, nd, xadj, adjncy)
    call make_matrix(n, ndof, xadj, adjncy, ismult, 'smoothspd', hecMAT)
    if (me == 0) then
      call hecmw_ordering_gen(n, n + (xadj(n+1)-1)/2, xadj, adjncy, perm, invp, 0, 0)
      call hecmw_mf_symbolic_build(graph, perm, relax, sym)
      call hecmw_mf_symbolic_cluster(graph, sym, perm)
      call hecmw_mf_symbolic_finalize(sym)
      call hecmw_mf_symbolic_build(graph, perm, relax, sym)
    endif
    call hecmw_mf_dist_symbolic_bcast(sym, 0)
    call hecmw_mf_dist_map_build(sym, map)
    call hecmw_mf_dist_gmat_extract(hecMAT, sym, map, pmat)
    allocate(x(n*ndof), r(n*ndof))
    bnrm = sqrt(dot_product(hecMAT%B, hecMAT%B))

    ! (r2)
    call hecmw_mf_numeric_init(graph, sym, tile, .true., MF_BLR_ETA, fct)
    fct%eps = 1.0d-8
    fct%blr_reuse = .true.
    fct%blr = .true.
    call hecmw_mf_numeric_factor_mpi(pmat, sym, map, fct, ierr)
    iw(1) = int(fct%blr_tiles_lr, kind=kint)
    fct%blr = .false.
    if (ierr == 0) call hecmw_mf_numeric_factor_mpi(pmat, sym, map, fct, ierr)
    fct%blr = .true.
    if (ierr == 0) call hecmw_mf_numeric_factor_mpi(pmat, sym, map, fct, ierr)
    iw(2) = int(fct%blr_tiles_lr, kind=kint)
    call hecmw_allreduce_I_comm(iw, 2, hecmw_sum, hecmw_comm_get_comm())
    if (me == 0) write(*,'(a,a,a,i0,a,i0)') '  ', name, ': low rank tiles first ', iw(1), &
      ', after a factorization without BLR ', iw(2)
    if (.not. (ierr == 0 .and. iw(2) == iw(1))) ok = .false.
    call hecmw_mf_numeric_finalize(fct)

    ! (r1)
    call hecmw_mf_numeric_init(graph, sym, tile, .true., MF_BLR_ETA, fct)
    fct%blr = .true.
    fct%eps = 1.0d-8
    fct%blr_reuse = .true.
    fct%mode = 1
    call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
    if (ierr == 0) then
      fct%mode = 2
      call hecmw_mf_numeric_factor(hecMAT, sym, fct, ierr)
    endif
    resid = huge(1.0d0)
    if (ierr == 0) then
      call hecmw_mf_numeric_solve(sym, fct, hecMAT%B, x)
      call residual(hecMAT, x, r)
      resid = sqrt(dot_product(r, r)) / bnrm
    endif
    if (.not. (ierr == 0 .and. resid <= 1.0d-6)) then
      if (me == 0) write(*,'(a,a,a,i0,a,1pe12.4)') '  ', name, ': LU after LDLt, ierr = ', ierr, &
        ' residual = ', resid
      ok = .false.
    endif
    call hecmw_mf_numeric_finalize(fct)

    iw(1) = merge(0, 1, ok)
    call hecmw_allreduce_I_comm(iw, 1, hecmw_max, hecmw_comm_get_comm())
    ok = iw(1) == 0
    if (me == 0) then
      if (ok) then
        write(*,'(a,a,a,1pe10.2)') 'PASS ', name, ': LU after LDLt resid=', resid
      else
        write(*,'(a,a)') 'FAIL ', name
      endif
    endif
    if (.not. ok) nfail = nfail + 1
    call hecmw_mf_dist_map_finalize(map)
    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    deallocate(hecMAT%D, hecMAT%AL, hecMAT%AU, hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU)
    deallocate(hecMAT%B, hecMAT%X)
    deallocate(pmat%D, pmat%AL, pmat%AU, pmat%indexL, pmat%indexU, pmat%itemL, pmat%itemU)
    deallocate(xadj, adjncy, perm, invp, nd, ismult, x, r)
  end subroutine run_reuse_case

  !> (j) the BLR compression of a tile holding a NaN or an Inf gives up (rank -1)
  subroutine run_compress_nonfinite(name)
    character(len=*), intent(in) :: name
    real(kind=kreal) :: a(12,8), v(8,8), bad(2)
    integer(kind=kint) :: i, j, k, rk(2)
    logical :: ok

    ntest = ntest + 1
    bad(1) = ieee_value(1.0d0, ieee_quiet_nan)
    bad(2) = ieee_value(1.0d0, ieee_positive_inf)
    do k = 1, 2
      do j = 1, 8
        do i = 1, 12
          a(i,j) = 1.0d0 / real(i + j, kind=kreal)
        enddo
      enddo
      a(5,3) = bad(k)
      call hecmw_mf_kernel_compress(12, 8, 12, a, 1.0d-8, 8, v, rk(k))
    enddo
    ok = all(rk == -1)
    if (me == 0) then
      if (ok) then
        write(*,'(a,a)') 'PASS ', name
      else
        write(*,'(a,a,a,i0,a,i0)') 'FAIL ', name, ': rank with NaN ', rk(1), ', with Inf ', rk(2)
      endif
    endif
    if (.not. ok) nfail = nfail + 1
  end subroutine run_compress_nonfinite

  !> a condition that must hold on every rank
  function reduced(cond) result(res)
    logical, intent(in) :: cond
    logical :: res
    integer(kind=kint) :: iw(1)

    iw(1) = merge(1, 0, cond)
    call hecmw_allreduce_I_comm(iw, 1, hecmw_min, hecmw_comm_get_comm())
    res = iw(1) == 1
  end function reduced

  !> Grid graph (chain, 2D 9-point, 3D 27-point); with mstride > 0 every mstride-th grid node
  !> gets a multiplier node attached by a single edge (ismult).
  subroutine make_graph(dim, nx, ny, nz, mstride, n, xadj, adjncy, ismult)
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, mstride
    integer(kind=kint), intent(out) :: n
    integer(kind=kint), allocatable, intent(out) :: xadj(:), adjncy(:)
    logical, allocatable, intent(out) :: ismult(:)
    integer(kind=kint), allocatable :: mult(:)
    integer(kind=kint) :: ix, iy, iz, jx, jy, jz, i, j, cnt, my, mz, pass, ng, nm

    my = ny
    mz = nz
    if (dim < 2) my = 1
    if (dim < 3) mz = 1
    ng = nx*my*mz
    allocate(mult(ng))
    mult(1:ng) = 0
    nm = 0
    if (mstride > 0) then
      do i = 1, ng, mstride
        nm = nm + 1
        mult(i) = ng + nm
      enddo
    endif
    n = ng + nm
    allocate(xadj(n+1), ismult(n))
    ismult(1:ng) = .false.
    ismult(ng+1:n) = .true.
    do pass = 1, 2
      cnt = 0
      xadj(1) = 1
      do iz = 1, mz
        do iy = 1, my
          do ix = 1, nx
            i = ix + nx*((iy-1) + my*(iz-1))
            do jz = max(1,iz-1), min(mz,iz+1)
              do jy = max(1,iy-1), min(my,iy+1)
                do jx = max(1,ix-1), min(nx,ix+1)
                  j = jx + nx*((jy-1) + my*(jz-1))
                  if (j == i) cycle
                  cnt = cnt + 1
                  if (pass == 2) adjncy(cnt) = j
                enddo
              enddo
            enddo
            if (mult(i) > 0) then
              cnt = cnt + 1
              if (pass == 2) adjncy(cnt) = mult(i)
            endif
            xadj(i+1) = cnt + 1
          enddo
        enddo
      enddo
      do i = 1, ng
        if (mult(i) == 0) cycle
        cnt = cnt + 1
        if (pass == 2) adjncy(cnt) = i
        xadj(mult(i)+1) = cnt + 1
      enddo
      if (pass == 1) allocate(adjncy(cnt))
    enddo
    deallocate(mult)
  end subroutine make_graph

  !> Block matrix in the row major block layout of hecmwST_matrix, the kinds of test_mf_blr:
  !> 'smooth' couplings are outer products of per-node vectors (low rank fronts), 'rand'
  !> couplings are random blocks; smoothsad carries zero diagonal multiplier nodes (delays),
  !> smoothlu unsymmetric couplings (LU), randoff a small random diagonal (row swaps).
  subroutine make_matrix(n, ndof, xadj, adjncy, ismult, kind, hecMAT)
    integer(kind=kint), intent(in) :: n, ndof
    integer(kind=kint), intent(in) :: xadj(:), adjncy(:)
    logical, intent(in) :: ismult(:)
    character(len=*), intent(in) :: kind
    type(hecmwST_matrix), intent(out) :: hecMAT
    integer(kind=kint) :: i, j, p, nl, nu, a, b, ku, kl, nd2
    integer(kind=kint), allocatable :: upos(:)
    real(kind=kreal), allocatable :: rowsum(:), uvec(:,:), vvec(:,:)
    real(kind=kreal) :: v
    logical :: smooth, lu

    smooth = (kind(1:6) == 'smooth')
    lu = (kind == 'smoothlu' .or. kind == 'smoothnsd' .or. kind == 'randoff' .or. kind == 'zerodiaglu')
    nd2 = ndof*ndof
    hecMAT%N = n
    hecMAT%NP = n
    hecMAT%NDOF = ndof
    hecMAT%symmetric = .not. lu
    allocate(hecMAT%indexL(0:n), hecMAT%indexU(0:n))
    allocate(hecMAT%itemL(xadj(n+1)-1), hecMAT%itemU(xadj(n+1)-1))
    hecMAT%indexL(0) = 0
    hecMAT%indexU(0) = 0
    nl = 0
    nu = 0
    do i = 1, n
      do p = xadj(i), xadj(i+1)-1
        j = adjncy(p)
        if (j < i) then
          nl = nl + 1
          hecMAT%itemL(nl) = j
        else
          nu = nu + 1
          hecMAT%itemU(nu) = j
        endif
      enddo
      hecMAT%indexL(i) = nl
      hecMAT%indexU(i) = nu
    enddo
    hecMAT%NPL = nl
    hecMAT%NPU = nu
    allocate(hecMAT%D(n*nd2), hecMAT%AL(nl*nd2), hecMAT%AU(nu*nd2), hecMAT%B(n*ndof), hecMAT%X(n*ndof))
    allocate(rowsum(n*ndof), upos(n), uvec(ndof, n), vvec(ndof, n))
    rowsum(:) = 0.0d0
    do i = 1, n
      do a = 1, ndof
        uvec(a, i) = 2.0d0*rand01() - 1.0d0
        vvec(a, i) = 2.0d0*rand01() - 1.0d0
      enddo
    enddo

    do i = 1, n
      do ku = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        j = hecMAT%itemU(ku)
        do b = 1, ndof
          do a = 1, ndof
            if (smooth .and. .not. (ismult(i) .or. ismult(j))) then
              ! rank-1 coupling; multiplier couplings stay full rank (else singular)
              v = uvec(a, i) * vvec(b, j)
            else
              v = 2.0d0*rand01() - 1.0d0
            endif
            hecMAT%AU((ku-1)*nd2 + (a-1)*ndof + b) = v
          enddo
        enddo
      enddo
    enddo
    upos(1:n) = hecMAT%indexU(0:n-1) + 1
    do j = 1, n
      do kl = hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i = hecMAT%itemL(kl)
        ku = upos(i)
        upos(i) = upos(i) + 1
        do a = 1, ndof
          do b = 1, ndof
            hecMAT%AL((kl-1)*nd2 + (a-1)*ndof + b) = hecMAT%AU((ku-1)*nd2 + (b-1)*ndof + a)
          enddo
        enddo
      enddo
    enddo
    if (lu) then
      ! unsymmetric values: an independent outer product (or random blocks) in the lower part
      do j = 1, n
        do kl = hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
          i = hecMAT%itemL(kl)
          do b = 1, ndof
            do a = 1, ndof
              if (smooth .and. .not. (ismult(i) .or. ismult(j))) then
                v = uvec(a, j) * vvec(b, i)
              else
                v = 2.0d0*rand01() - 1.0d0
              endif
              hecMAT%AL((kl-1)*nd2 + (a-1)*ndof + b) = v
            enddo
          enddo
        enddo
      enddo
    endif

    do i = 1, n
      do kl = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        do a = 1, ndof
          do b = 1, ndof
            rowsum((i-1)*ndof + a) = rowsum((i-1)*ndof + a) + abs(hecMAT%AL((kl-1)*nd2 + (a-1)*ndof + b))
          enddo
        enddo
      enddo
      do ku = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        do a = 1, ndof
          do b = 1, ndof
            rowsum((i-1)*ndof + a) = rowsum((i-1)*ndof + a) + abs(hecMAT%AU((ku-1)*nd2 + (a-1)*ndof + b))
          enddo
        enddo
      enddo
      do a = 1, ndof
        do b = a+1, ndof
          v = 2.0d0*rand01() - 1.0d0
          hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) = v
          hecMAT%D((i-1)*nd2 + (b-1)*ndof + a) = v
          rowsum((i-1)*ndof + a) = rowsum((i-1)*ndof + a) + abs(v)
          rowsum((i-1)*ndof + b) = rowsum((i-1)*ndof + b) + abs(v)
        enddo
      enddo
      do a = 1, ndof
        hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = rowsum((i-1)*ndof + a) + 1.0d0 + rand01()
      enddo
      if (kind == 'randoff') then
        do a = 1, ndof
          do b = 1, ndof
            if (a /= b) hecMAT%D((i-1)*nd2 + (a-1)*ndof + b) = 2.0d0*rand01() - 1.0d0
          enddo
          hecMAT%D((i-1)*nd2 + (a-1)*ndof + a) = 0.1d0*(2.0d0*rand01() - 1.0d0)
        enddo
      endif
      if ((kind == 'smoothsad' .or. kind == 'smoothnsd') .and. ismult(i)) then
        hecMAT%D((i-1)*nd2 + 1:i*nd2) = 0.0d0
      endif
      if (kind(1:8) == 'zerodiag') hecMAT%D((i-1)*nd2 + 1:i*nd2) = 0.0d0
    enddo

    do i = 1, n*ndof
      hecMAT%B(i) = 2.0d0*rand01() - 1.0d0
      hecMAT%X(i) = 0.0d0
    enddo
    deallocate(rowsum, upos, uvec, vvec)
  end subroutine make_matrix

  !> r = B - K X with the full (D, AL, AU) block storage.
  subroutine residual(hecMAT, x, r)
    type(hecmwST_matrix), intent(in) :: hecMAT
    real(kind=kreal), intent(in) :: x(:)
    real(kind=kreal), intent(out) :: r(:)
    integer(kind=kint) :: i, j, k, a, b, ndof, nd2

    ndof = hecMAT%NDOF
    nd2 = ndof*ndof
    r(:) = hecMAT%B(:)
    do i = 1, hecMAT%NP
      do a = 1, ndof
        do b = 1, ndof
          r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%D((i-1)*nd2 + (a-1)*ndof + b)*x((i-1)*ndof + b)
        enddo
      enddo
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        j = hecMAT%itemL(k)
        do a = 1, ndof
          do b = 1, ndof
            r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%AL((k-1)*nd2 + (a-1)*ndof + b)*x((j-1)*ndof + b)
          enddo
        enddo
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        j = hecMAT%itemU(k)
        do a = 1, ndof
          do b = 1, ndof
            r((i-1)*ndof + a) = r((i-1)*ndof + a) - hecMAT%AU((k-1)*nd2 + (a-1)*ndof + b)*x((j-1)*ndof + b)
          enddo
        enddo
      enddo
    enddo
  end subroutine residual

  function rand01() result(v)
    real(kind=kreal) :: v

    seed = mod(1103515245_8*seed + 12345_8, 2147483648_8)
    v = real(seed, kind=kreal) / 2147483648.0d0
  end function rand01

end program test_mf_mpi
