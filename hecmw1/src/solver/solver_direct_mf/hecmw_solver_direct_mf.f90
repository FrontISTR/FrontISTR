!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Multifrontal direct solver (METHOD=DIRECTmf): LDLt with threshold pivoting for
!>        symmetric (possibly indefinite) matrices, LU for structurally symmetric matrices
!>        with unsymmetric values. Multi-process runs hold the global matrix partially per
!>        rank, run the symbolic stage on rank 0 and distribute the numeric stage
!>        subtree-to-subcube.
module hecmw_solver_direct_mf
  use hecmw_util
  use m_hecmw_comm_f
  use hecmw_matrix_misc
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_mf_dist
  use hecmw_mf_numeric
  use hecmw_mf_kernel
  use hecmw_solver_las
  use hecmw_solver_misc
  use hecmw_solver_direct
  !$ use omp_lib
  implicit none

  private
  public :: hecmw_solve_direct_mf
  public :: hecmw_mf_compare_builtin

  !> relaxed amalgamation target in columns (DOFs); converted to nodes with the block size
  integer(kind=kint), parameter :: MF_RELAX_COLS = 16
  !> target tile size in DOFs
  integer(kind=kint), parameter :: MF_TILE = 256
  !> default refinement steps of a BLR solve; a negative option turns the refinement off
  integer(kind=kint), parameter :: MF_IR_STEPS = 5
  !> default refinement stopping residual
  real(kind=kreal), parameter :: MF_IR_TOL = 1.0d-8

  type(hecmwST_mf_symbolic), save :: SYM
  type(hecmwST_mf_factor), save :: FCT
  type(hecmwST_mf_map), save :: MAP
  type(hecmwST_mf_gmat), save :: GMAT

contains

  subroutine hecmw_solve_direct_mf(hecMESH, hecMAT, imsg)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: imsg
    type(hecmwST_mf_graph) :: graph
    integer(kind=kint), allocatable :: perm(:), invp(:)
    integer(kind=kint) :: loglevel, ordering, n, nerr, relax, tile, ierr, idof, nthreads, irmax, it, eta
    real(kind=kreal) :: t1, t2, irtol, bnrm, rnrm, tcomm
    real(kind=kreal), allocatable :: rr(:), dd(:)
    logical :: clustered

    loglevel = hecmw_mat_get_loglevel(hecMAT)
    if (loglevel < 0) loglevel = max(hecmw_mat_get_timelog(hecMAT), hecmw_mat_get_iterlog(hecMAT))
    if (hecmw_comm_get_rank() /= 0) loglevel = 0

    if (hecMESH%PETOT > 1) then
      call hecmw_solve_direct_mf_dist(hecMESH, hecMAT, imsg, loglevel)
      return
    endif

    if (hecMAT%Iarray(98) == 1) then
      t1 = hecmw_wtime()
      call hecmw_mf_graph_from_hecmat(hecMAT, graph)
      n = graph%nnode
      allocate(perm(n), invp(n))
      ordering = hecMAT%Iarray(41)
      call hecmw_ordering_gen(n, n + (graph%xadj(n+1)-1)/2, graph%xadj, graph%adjncy, perm, invp, ordering, loglevel)
      relax = hecMAT%Iarray(46)
      if (relax <= 0) relax = MF_RELAX_COLS
      relax = max(1, relax / hecMAT%NDOF)
      call hecmw_mf_symbolic_finalize(SYM)
      call hecmw_mf_symbolic_build(graph, perm, relax, SYM)
      clustered = hecMAT%Iarray(43) /= 0 .and. hecmw_mf_kernel_blr_available()
      if (clustered) then
        call hecmw_mf_symbolic_cluster(graph, SYM, perm)
        call hecmw_mf_symbolic_finalize(SYM)
        call hecmw_mf_symbolic_build(graph, perm, relax, SYM)
      endif
      tile = hecMAT%Iarray(45)
      if (tile <= 0) tile = MF_TILE
      eta = hecMAT%Iarray(47)
      if (eta == 0) eta = MF_BLR_ETA
      if (eta < 0) eta = 0
      call hecmw_mf_numeric_finalize(FCT)
      call hecmw_mf_numeric_init(graph, SYM, tile, clustered, eta, FCT)
      t2 = hecmw_wtime()
      if (loglevel > 0) then
        write(*,'(a,f10.3,a)') '[DIRECTmf]: symbolic fct done (', t2 - t1, ' sec)'
        if (clustered) write(*,'(a)') '[DIRECTmf]: separator nodes clustered for BLR'
        call hecmw_mf_symbolic_print(SYM)
        call hecmw_mf_numeric_print(FCT)
      endif
      if (loglevel > 1) then
        call hecmw_mf_symbolic_check(graph, SYM, nerr)
        write(*,'(a,i0)') '[DIRECTmf]: self check violations = ', nerr
        if (.not. clustered) then
          call hecmw_mf_compare_builtin(hecMESH, hecMAT, graph, ordering, nerr)
          write(*,'(a,i0)') '[DIRECTmf]: mismatches against built-in DIRECT = ', nerr
        endif
      endif
      deallocate(perm, invp)
      call hecmw_mf_graph_finalize(graph)
    endif

    hecMAT%Iarray(98) = 0

    if (hecMAT%Iarray(97) == 1) then
      FCT%mode = hecMAT%Iarray(42)
      FCT%scan = loglevel > 1
      FCT%blr = hecMAT%Iarray(43) /= 0
      FCT%eps = hecMAT%Rarray(41)
      FCT%pivot_u = hecMAT%Rarray(43)
      FCT%pivot_zero = hecMAT%Rarray(44)
      FCT%blr_beta = hecMAT%Rarray(45)
      if (FCT%blr_beta <= 0.0d0) FCT%blr_beta = 1.0d0
      FCT%blr_reuse = hecMAT%Iarray(48) /= 0
      if (FCT%blr .and. .not. hecmw_mf_kernel_blr_available()) then
        FCT%blr = .false.
        if (loglevel > 0) write(*,'(a)') '[DIRECTmf]: BLR disabled (built without LAPACK)'
      endif
      t1 = hecmw_wtime()
      call hecmw_mf_numeric_factor(hecMAT, SYM, FCT, ierr)
      t2 = hecmw_wtime()
      if (ierr /= 0) then
        if (ierr > 0) then
          idof = FCT%pdof(ierr)
          write(imsg,'(a,i0,a,i0,a)') 'ERROR: DIRECTmf: zero pivot at node ', (idof-1)/hecMAT%NDOF + 1, &
            ' dof ', mod(idof-1, hecMAT%NDOF) + 1, ' (matrix is singular)'
          write(*,'(a,i0,a,i0,a)') 'ERROR: DIRECTmf: zero pivot at node ', (idof-1)/hecMAT%NDOF + 1, &
            ' dof ', mod(idof-1, hecMAT%NDOF) + 1, ' (matrix is singular)'
        else if (ierr == -1) then
          write(imsg,*) 'ERROR: DIRECTmf: block size of the matrix does not match the symbolic structure'
          write(*,*) 'ERROR: DIRECTmf: block size of the matrix does not match the symbolic structure'
        else
          write(imsg,*) 'ERROR: DIRECTmf: the nonzero structure of the matrix is not symmetric'
          write(*,*) 'ERROR: DIRECTmf: the nonzero structure of the matrix is not symmetric'
        endif
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      hecMAT%Iarray(97) = 0
      if (loglevel > 0) then
        write(*,'(a,f10.3,a)') '[DIRECTmf]: numeric fct done (', t2 - t1, ' sec)'
        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        write(*,'(a,i0)') '[DIRECTmf]: threads = ', nthreads
        write(*,'(a,1pe9.2,a,1pe9.2)') '[DIRECTmf]: options: pivot u = ', FCT%pivot_u, &
          ', zero pivot factor = ', FCT%pivot_zero
        if (FCT%lu) then
          if (FCT%mode == 2) then
            write(*,'(a)') '[DIRECTmf]: mode = LU (forced by option)'
          else
            write(*,'(a)') '[DIRECTmf]: mode = LU (unsymmetric matrix)'
          endif
          if (loglevel > 1) write(*,'(a,1pe9.2)') '[DIRECTmf]: asymmetry = ', FCT%asym
          write(*,'(a,i0,a,i0,a,i0,a)') '[DIRECTmf]: row swaps = ', FCT%n_swap, ', delayed = ', FCT%n_delay, &
            ' (max front growth = ', FCT%max_growth, ' dofs)'
        else
          if (FCT%mode == 1) then
            write(*,'(a)') '[DIRECTmf]: mode = LDLt (forced by option)'
          else
            write(*,'(a)') '[DIRECTmf]: mode = LDLt'
          endif
          if (loglevel > 1) write(*,'(a,1pe9.2)') '[DIRECTmf]: asymmetry = ', FCT%asym
          write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') '[DIRECTmf]: inertia (+,-) = (', FCT%n_pos, ',', FCT%n_neg, &
            '), 2x2 pivots = ', FCT%n_2x2, ', swaps = ', FCT%n_swap, ', delayed = ', FCT%n_delay, &
            ' (max front growth = ', FCT%max_growth, ' dofs)'
        endif
        write(*,'(a,i0,a,f10.3,a,i0,a,i0)') '[DIRECTmf]: factor words = ', FCT%factor_words_act, ' (', &
          real(FCT%factor_words_act, kind=kreal)*8.0d0/1024.0d0**3, ' GB), stack peak words = ', FCT%stack_peak_act, &
          ', front words = ', FCT%front_words_act
        if (FCT%blr) then
          write(*,'(a,1pe9.2,a,i0,a,0pf6.1,a)') '[DIRECTmf]: BLR: eps = ', FCT%eps, &
            ', factor words full rank = ', FCT%blr_words_fr, ' (', &
            100.0d0*(1.0d0 - real(FCT%factor_words_act, kind=kreal)/max(real(FCT%blr_words_fr, kind=kreal), 1.0d0)), &
            '% saved)'
          write(*,'(a,i0,a,f5.2,a,l1)') '[DIRECTmf]: BLR: skip eta = ', FCT%blr_eta, &
            ', gain cap beta = ', FCT%blr_beta, ', rank reuse = ', FCT%blr_reuse
          write(*,'(a,i0,a,i0,a,i0,a,i0,a,f8.1)') '[DIRECTmf]: BLR: tiles compressed = ', FCT%blr_tiles_lr, ' of ', &
            FCT%blr_tiles, ', skipped = ', FCT%blr_tiles_skip, ', rank max = ', FCT%blr_rank_max, ', rank avg = ', &
            real(FCT%blr_rank_sum, kind=kreal)/max(real(FCT%blr_tiles_lr, kind=kreal), 1.0d0)
        endif
        write(*,'(a,i0,a,f10.3,a)') '[DIRECTmf]: front words held concurrently (peak) = ', FCT%front_peak, ' (', &
          real(FCT%front_peak, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
      endif
    endif

    if (.not. FCT%factored) then
      write(imsg,*) 'ERROR: DIRECTmf: numeric factorization not performed'
      write(*,*) 'ERROR: DIRECTmf: numeric factorization not performed'
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    t1 = hecmw_wtime()
    call hecmw_mf_numeric_solve(SYM, FCT, hecMAT%B, hecMAT%X)
    t2 = hecmw_wtime()
    if (loglevel > 0) write(*,'(a,f10.3,a)') '[DIRECTmf]: solve done (', t2 - t1, ' sec)'

    irmax = hecMAT%Iarray(44)
    if (FCT%blr .and. irmax == 0) irmax = MF_IR_STEPS
    if (irmax > 0) then
      irtol = hecMAT%Rarray(42)
      if (.not. (irtol > 0.0d0)) irtol = MF_IR_TOL
      n = hecMAT%NP * hecMAT%NDOF
      allocate(rr(n), dd(n))
      tcomm = 0.0d0
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnrm, tcomm)
      t1 = hecmw_wtime()
      do it = 1, irmax
        call hecmw_matresid(hecMESH, hecMAT, hecMAT%X, hecMAT%B, rr, tcomm)
        call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, rr, rr, rnrm, tcomm)
        rnrm = sqrt(rnrm / max(bnrm, tiny(bnrm)))
        if (loglevel > 0) write(*,'(a,i0,a,1pe11.4)') '[DIRECTmf]: refinement ', it - 1, ': residual = ', rnrm
        if (rnrm <= irtol) exit
        call hecmw_mf_numeric_solve(SYM, FCT, rr, dd)
        hecMAT%X(1:n) = hecMAT%X(1:n) + dd(1:n)
      enddo
      t2 = hecmw_wtime()
      if (loglevel > 0) write(*,'(a,f10.3,a)') '[DIRECTmf]: refinement done (', t2 - t1, ' sec)'
      deallocate(rr, dd)
    endif
  end subroutine hecmw_solve_direct_mf

  !> Multi-process solve: hecmw_mf_dist gathers the global structure (the full profile on
  !> rank 0 only), rank 0 runs the symbolic stage on it and broadcasts the structure, each
  !> rank keeps the partial matrix its supernodes read with the values transferred by
  !> destination before every factorization, and the numeric stage is distributed
  !> subtree-to-subcube. The built-in comparison of loglevel > 1 is skipped (it reads the
  !> local matrix, which no longer matches the global structure). Pivot and BLR statistics
  !> are reduced over the ranks for the log; the factor stays distributed, each rank holding
  !> the panels of its own supernodes only.
  subroutine hecmw_solve_direct_mf_dist(hecMESH, hecMAT, imsg, loglevel)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: imsg, loglevel
    type(hecmwST_mf_graph) :: graph
    integer(kind=kint), allocatable :: perm(:), invp(:)
    real(kind=kreal), allocatable :: gb(:), gx(:), rr(:), wpr(:), wmr(:)
    integer(kind=kint) :: ordering, n, nerr, relax, tile, ierr, idof, nthreads, irmax, it, i, ofs, eta
    integer(kind=8) :: wrepl
    real(kind=kreal) :: t1, t2, irtol, bnrm, rnrm, tcomm
    logical :: clustered

    clustered = .false.
    if (hecMAT%Iarray(98) == 1) then
      t1 = hecmw_wtime()
      call hecmw_mf_dist_gmat_build(hecMESH, hecMAT, GMAT)
      call hecmw_mf_symbolic_finalize(SYM)
      relax = hecMAT%Iarray(46)
      if (relax <= 0) relax = MF_RELAX_COLS
      relax = max(1, relax / hecMAT%NDOF)
      clustered = hecMAT%Iarray(43) /= 0 .and. hecmw_mf_kernel_blr_available()
      if (hecmw_comm_get_rank() == 0) then
        call hecmw_mf_graph_from_hecmat(GMAT%mat, graph)
        n = graph%nnode
        allocate(perm(n), invp(n))
        ordering = hecMAT%Iarray(41)
        call hecmw_ordering_gen(n, n + (graph%xadj(n+1)-1)/2, graph%xadj, graph%adjncy, perm, invp, ordering, loglevel)
        call hecmw_mf_symbolic_build(graph, perm, relax, SYM)
        if (clustered) then
          call hecmw_mf_symbolic_cluster(graph, SYM, perm)
          call hecmw_mf_symbolic_finalize(SYM)
          call hecmw_mf_symbolic_build(graph, perm, relax, SYM)
        endif
        if (loglevel > 1) then
          call hecmw_mf_symbolic_check(graph, SYM, nerr)
          write(*,'(a,i0)') '[DIRECTmf]: self check violations = ', nerr
        endif
        deallocate(perm, invp)
      endif
      call hecmw_mf_dist_symbolic_bcast(SYM, 0)
      call hecmw_mf_dist_map_finalize(MAP)
      call hecmw_mf_dist_map_build(SYM, MAP)
      call hecmw_mf_dist_gmat_part(SYM, MAP, GMAT)
      allocate(wmr(MAP%nprocs))
      call mf_gather_matwords(GMAT, wmr, wrepl)
      tile = hecMAT%Iarray(45)
      if (tile <= 0) tile = MF_TILE
      eta = hecMAT%Iarray(47)
      if (eta == 0) eta = MF_BLR_ETA
      if (eta < 0) eta = 0
      call hecmw_mf_numeric_finalize(FCT)
      call hecmw_mf_numeric_init(graph, SYM, tile, clustered, eta, FCT)
      call hecmw_mf_graph_finalize(graph)
      if (clustered .and. eta > 0) call hecmw_mf_numeric_admis_bcast(FCT, 0)
      t2 = hecmw_wtime()
      if (loglevel > 0) then
        write(*,'(a,f10.3,a)') '[DIRECTmf]: symbolic fct done (', t2 - t1, ' sec)'
        if (clustered) write(*,'(a)') '[DIRECTmf]: separator nodes clustered for BLR'
        write(*,'(a,i0,a,i0)') '[DIRECTmf]: MPI ranks = ', MAP%nprocs, ', upper fronts = ', MAP%nupper
        write(*,'(a,i0,a,*(i0,1x))') '[DIRECTmf]: matrix words replicated = ', wrepl, ', held per rank = ', &
          (nint(wmr(i), kind=8), i = 1, MAP%nprocs)
        call hecmw_mf_symbolic_print(SYM)
        call hecmw_mf_numeric_print(FCT)
      endif
      deallocate(wmr)
      hecMAT%Iarray(98) = 0
    endif

    if (hecMAT%Iarray(97) == 1) then
      FCT%mode = hecMAT%Iarray(42)
      FCT%scan = loglevel > 1
      FCT%blr = hecMAT%Iarray(43) /= 0
      FCT%eps = hecMAT%Rarray(41)
      FCT%pivot_u = hecMAT%Rarray(43)
      FCT%pivot_zero = hecMAT%Rarray(44)
      FCT%blr_beta = hecMAT%Rarray(45)
      if (FCT%blr_beta <= 0.0d0) FCT%blr_beta = 1.0d0
      FCT%blr_reuse = hecMAT%Iarray(48) /= 0
      if (FCT%blr .and. .not. hecmw_mf_kernel_blr_available()) then
        FCT%blr = .false.
        if (loglevel > 0) write(*,'(a)') '[DIRECTmf]: BLR disabled (built without LAPACK)'
      endif
      call hecmw_mf_dist_gmat_vals(hecMAT, GMAT)
      t1 = hecmw_wtime()
      call hecmw_mf_numeric_factor_mpi(GMAT%mat, SYM, MAP, FCT, ierr)
      t2 = hecmw_wtime()
      if (ierr /= 0) then
        if (hecmw_comm_get_rank() == 0) then
          if (ierr > 0) then
            idof = FCT%pdof(ierr)
            write(imsg,'(a,i0,a,i0,a)') 'ERROR: DIRECTmf: zero pivot at node ', (idof-1)/hecMAT%NDOF + 1, &
              ' dof ', mod(idof-1, hecMAT%NDOF) + 1, ' (matrix is singular)'
            write(*,'(a,i0,a,i0,a)') 'ERROR: DIRECTmf: zero pivot at node ', (idof-1)/hecMAT%NDOF + 1, &
              ' dof ', mod(idof-1, hecMAT%NDOF) + 1, ' (matrix is singular)'
          else if (ierr == -1) then
            write(imsg,*) 'ERROR: DIRECTmf: block size of the matrix does not match the symbolic structure'
            write(*,*) 'ERROR: DIRECTmf: block size of the matrix does not match the symbolic structure'
          else
            write(imsg,*) 'ERROR: DIRECTmf: the nonzero structure of the matrix is not symmetric'
            write(*,*) 'ERROR: DIRECTmf: the nonzero structure of the matrix is not symmetric'
          endif
        endif
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      hecMAT%Iarray(97) = 0
      allocate(wpr(MAP%nprocs))
      call mf_reduce_stats(FCT, wpr)
      if (loglevel > 0) then
        write(*,'(a,f10.3,a)') '[DIRECTmf]: numeric fct done (', t2 - t1, ' sec)'
        nthreads = 1
        !$ nthreads = omp_get_max_threads()
        write(*,'(a,i0,a,i0)') '[DIRECTmf]: ranks = ', MAP%nprocs, ', threads per rank = ', nthreads
        write(*,'(a,1pe9.2,a,1pe9.2)') '[DIRECTmf]: options: pivot u = ', FCT%pivot_u, &
          ', zero pivot factor = ', FCT%pivot_zero
        if (FCT%lu) then
          if (FCT%mode == 2) then
            write(*,'(a)') '[DIRECTmf]: mode = LU (forced by option)'
          else
            write(*,'(a)') '[DIRECTmf]: mode = LU (unsymmetric matrix)'
          endif
          if (loglevel > 1) write(*,'(a,1pe9.2)') '[DIRECTmf]: asymmetry = ', FCT%asym
          write(*,'(a,i0,a,i0,a,i0,a)') '[DIRECTmf]: row swaps = ', FCT%n_swap, ', delayed = ', FCT%n_delay, &
            ' (max front growth = ', FCT%max_growth, ' dofs)'
        else
          if (FCT%mode == 1) then
            write(*,'(a)') '[DIRECTmf]: mode = LDLt (forced by option)'
          else
            write(*,'(a)') '[DIRECTmf]: mode = LDLt'
          endif
          if (loglevel > 1) write(*,'(a,1pe9.2)') '[DIRECTmf]: asymmetry = ', FCT%asym
          write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') '[DIRECTmf]: inertia (+,-) = (', FCT%n_pos, ',', FCT%n_neg, &
            '), 2x2 pivots = ', FCT%n_2x2, ', swaps = ', FCT%n_swap, ', delayed = ', FCT%n_delay, &
            ' (max front growth = ', FCT%max_growth, ' dofs)'
        endif
        write(*,'(a,i0,a,f10.3,a,i0,a,i0)') '[DIRECTmf]: factor words = ', FCT%factor_words_act, ' (', &
          real(FCT%factor_words_act, kind=kreal)*8.0d0/1024.0d0**3, ' GB), stack peak words = ', FCT%stack_peak_act, &
          ', front words = ', FCT%front_words_act
        write(*,'(a,*(i0,1x))') '[DIRECTmf]: factor words per rank = ', (nint(wpr(i), kind=8), i = 1, MAP%nprocs)
        if (FCT%blr) then
          write(*,'(a,1pe9.2,a,i0,a,0pf6.1,a)') '[DIRECTmf]: BLR: eps = ', FCT%eps, &
            ', factor words full rank = ', FCT%blr_words_fr, ' (', &
            100.0d0*(1.0d0 - real(FCT%factor_words_act, kind=kreal)/max(real(FCT%blr_words_fr, kind=kreal), 1.0d0)), &
            '% saved)'
          write(*,'(a,i0,a,f5.2,a,l1)') '[DIRECTmf]: BLR: skip eta = ', FCT%blr_eta, &
            ', gain cap beta = ', FCT%blr_beta, ', rank reuse = ', FCT%blr_reuse
          write(*,'(a,i0,a,i0,a,i0,a,i0,a,f8.1)') '[DIRECTmf]: BLR: tiles compressed = ', FCT%blr_tiles_lr, ' of ', &
            FCT%blr_tiles, ', skipped = ', FCT%blr_tiles_skip, ', rank max = ', FCT%blr_rank_max, ', rank avg = ', &
            real(FCT%blr_rank_sum, kind=kreal)/max(real(FCT%blr_tiles_lr, kind=kreal), 1.0d0)
        endif
        write(*,'(a,i0,a,f10.3,a)') '[DIRECTmf]: front words held concurrently (peak) = ', FCT%front_peak, ' (', &
          real(FCT%front_peak, kind=kreal)*8.0d0/1024.0d0**3, ' GB)'
      endif
      deallocate(wpr)
    endif

    if (.not. FCT%factored) then
      if (hecmw_comm_get_rank() == 0) then
        write(imsg,*) 'ERROR: DIRECTmf: numeric factorization not performed'
        write(*,*) 'ERROR: DIRECTmf: numeric factorization not performed'
      endif
      call hecmw_abort(hecmw_comm_get_comm())
    endif
    allocate(gb(FCT%ndof_tot), gx(FCT%ndof_tot))
    ofs = GMAT%ndisp(hecmw_comm_get_rank()+1) * hecMAT%NDOF
    call hecmw_mf_dist_gather_vec(GMAT, hecMAT%B, gb)
    t1 = hecmw_wtime()
    call hecmw_mf_numeric_solve_mpi(SYM, MAP, FCT, gb, gx)
    t2 = hecmw_wtime()
    if (loglevel > 0) write(*,'(a,f10.3,a)') '[DIRECTmf]: solve done (', t2 - t1, ' sec)'
    do i = 1, hecMAT%N * hecMAT%NDOF
      hecMAT%X(i) = gx(ofs + i)
    enddo

    irmax = hecMAT%Iarray(44)
    if (FCT%blr .and. irmax == 0) irmax = MF_IR_STEPS
    if (irmax > 0) then
      irtol = hecMAT%Rarray(42)
      if (.not. (irtol > 0.0d0)) irtol = MF_IR_TOL
      n = hecMAT%NP * hecMAT%NDOF
      allocate(rr(n))
      tcomm = 0.0d0
      call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, hecMAT%B, hecMAT%B, bnrm, tcomm)
      t1 = hecmw_wtime()
      do it = 1, irmax
        call hecmw_matresid(hecMESH, hecMAT, hecMAT%X, hecMAT%B, rr, tcomm)
        call hecmw_InnerProduct_R(hecMESH, hecMAT%NDOF, rr, rr, rnrm, tcomm)
        rnrm = sqrt(rnrm / max(bnrm, tiny(bnrm)))
        if (loglevel > 0) write(*,'(a,i0,a,1pe11.4)') '[DIRECTmf]: refinement ', it - 1, ': residual = ', rnrm
        if (rnrm <= irtol) exit
        call hecmw_mf_dist_gather_vec(GMAT, rr, gb)
        call hecmw_mf_numeric_solve_mpi(SYM, MAP, FCT, gb, gx)
        do i = 1, hecMAT%N * hecMAT%NDOF
          hecMAT%X(i) = hecMAT%X(i) + gx(ofs + i)
        enddo
      enddo
      t2 = hecmw_wtime()
      if (loglevel > 0) write(*,'(a,f10.3,a)') '[DIRECTmf]: refinement done (', t2 - t1, ' sec)'
      deallocate(rr)
    endif
    call hecmw_update_R(hecMESH, hecMAT%X, hecMAT%NP, hecMAT%NDOF)
    deallocate(gb, gx)
  end subroutine hecmw_solve_direct_mf_dist

  !> Gather the per-rank retained matrix words for the log and return the words the
  !> replicated design held per rank. The int8 words travel as reals, exact below 2^53.
  subroutine mf_gather_matwords(gmat, wmr, wrepl)
    implicit none
    type(hecmwST_mf_gmat), intent(in) :: gmat
    real(kind=kreal), intent(out) :: wmr(:)
    integer(kind=8), intent(out) :: wrepl
    real(kind=kreal) :: w1(1)
    integer(kind=kint), allocatable :: rcs(:), disp(:)
    integer(kind=kint) :: comm, np, r
    integer(kind=8) :: wpart

    comm = hecmw_comm_get_comm()
    np = hecmw_comm_get_size()
    call hecmw_mf_dist_gmat_words(gmat, wpart, wrepl)
    allocate(rcs(np), disp(np))
    do r = 1, np
      rcs(r) = 1
      disp(r) = r - 1
    enddo
    w1(1) = real(wpart, kind=kreal)
    call hecmw_allgatherv_real(w1, 1, wmr, rcs, disp, comm)
    deallocate(rcs, disp)
  end subroutine mf_gather_matwords

  !> Reduce the factorization statistics over the ranks for the log: counts and words are
  !> summed (the per-rank stack and front peaks summing to a bound on the concurrent global
  !> memory), the growth and rank maxima taken; wpr returns the per-rank factor words. The
  !> int8 words travel as reals, exact below 2^53.
  subroutine mf_reduce_stats(fct, wpr)
    implicit none
    type(hecmwST_mf_factor), intent(inout) :: fct
    real(kind=kreal), intent(out) :: wpr(:)
    real(kind=kreal) :: rs(13), rm(3), w1(1)
    integer(kind=kint), allocatable :: rcs(:), disp(:)
    integer(kind=kint) :: comm, np, r

    comm = hecmw_comm_get_comm()
    np = hecmw_comm_get_size()
    allocate(rcs(np), disp(np))
    do r = 1, np
      rcs(r) = 1
      disp(r) = r - 1
    enddo
    w1(1) = real(fct%factor_words_act, kind=kreal)
    call hecmw_allgatherv_real(w1, 1, wpr, rcs, disp, comm)
    deallocate(rcs, disp)
    rs(1) = real(fct%n_pos, kind=kreal)
    rs(2) = real(fct%n_neg, kind=kreal)
    rs(3) = real(fct%n_2x2, kind=kreal)
    rs(4) = real(fct%n_swap, kind=kreal)
    rs(5) = real(fct%n_delay, kind=kreal)
    rs(6) = real(fct%factor_words_act, kind=kreal)
    rs(7) = real(fct%stack_peak_act, kind=kreal)
    rs(8) = real(fct%front_peak, kind=kreal)
    rs(9) = real(fct%blr_words_fr, kind=kreal)
    rs(10) = real(fct%blr_tiles, kind=kreal)
    rs(11) = real(fct%blr_tiles_lr, kind=kreal)
    rs(12) = real(fct%blr_rank_sum, kind=kreal)
    rs(13) = real(fct%blr_tiles_skip, kind=kreal)
    call hecmw_allreduce_R_comm(rs, 13, hecmw_sum, comm)
    rm(1) = real(fct%max_growth, kind=kreal)
    rm(2) = real(fct%front_words_act, kind=kreal)
    rm(3) = real(fct%blr_rank_max, kind=kreal)
    call hecmw_allreduce_R_comm(rm, 3, hecmw_max, comm)
    fct%n_pos = nint(rs(1), kind=kint)
    fct%n_neg = nint(rs(2), kind=kint)
    fct%n_2x2 = nint(rs(3), kind=kint)
    fct%n_swap = nint(rs(4), kind=kint)
    fct%n_delay = nint(rs(5), kind=kint)
    fct%factor_words_act = nint(rs(6), kind=8)
    fct%stack_peak_act = nint(rs(7), kind=8)
    fct%front_peak = nint(rs(8), kind=8)
    fct%blr_words_fr = nint(rs(9), kind=8)
    fct%blr_tiles = nint(rs(10), kind=8)
    fct%blr_tiles_lr = nint(rs(11), kind=8)
    fct%blr_rank_sum = nint(rs(12), kind=8)
    fct%blr_tiles_skip = nint(rs(13), kind=8)
    fct%max_growth = nint(rm(1), kind=kint)
    fct%front_words_act = nint(rm(2), kind=8)
    fct%blr_rank_max = nint(rm(3), kind=kint)
  end subroutine mf_reduce_stats

  !> Compare the elimination tree and the column counts with the symbolic stage of the built-in
  !> direct solver, on the same node graph and the same permutation. Both trees are compared in
  !> the original node labels since the two codes postorder differently.
  subroutine hecmw_mf_compare_builtin(hecMESH, hecMAT, graph, ordering, nerr)
    implicit none
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_graph), intent(in) :: graph
    integer(kind=kint), intent(in) :: ordering
    integer(kind=kint), intent(out) :: nerr
    type(cholesky_factor) :: fct
    type(hecmwST_mf_symbolic) :: ref
    integer(kind=kint), allocatable :: cnt(:)
    integer(kind=kint) :: n, k, l, ir, pb, myk, mp, nbad_tree, nbad_cnt

    nerr = 0
    call SETIJ(hecMESH, hecMAT, fct)
    call MATINI(fct, ordering, 0, ir)
    if (ir /= 0) then
      write(*,'(a,i0)') '[DIRECTmf] compare: built-in MATINI failed, ir = ', ir
      nerr = 1
      return
    endif
    n = fct%NEQns
    call hecmw_mf_symbolic_build(graph, fct%IPErm(1:n), 0, ref)

    nbad_tree = 0
    do k = 1, n
      pb = fct%PARent(k)
      myk = ref%invp(fct%IPErm(k))
      mp = ref%parent(myk)
      if (pb == 0 .or. pb == n + 1) then
        if (mp /= 0) nbad_tree = nbad_tree + 1
      else
        if (mp == 0) then
          nbad_tree = nbad_tree + 1
        else if (ref%perm(mp) /= fct%IPErm(pb)) then
          nbad_tree = nbad_tree + 1
        endif
      endif
    enddo

    allocate(cnt(n))
    cnt(1:n) = 1
    do k = 1, n
      do l = fct%XLNzr(k), fct%XLNzr(k+1) - 1
        cnt(fct%COLno(l)) = cnt(fct%COLno(l)) + 1
      enddo
    enddo
    nbad_cnt = 0
    do k = 1, n
      if (cnt(k) /= ref%colcnt(ref%invp(fct%IPErm(k)))) nbad_cnt = nbad_cnt + 1
    enddo
    if (nbad_tree > 0) write(*,'(a,i0)') '[DIRECTmf] compare: elimination tree mismatches = ', nbad_tree
    if (nbad_cnt > 0) write(*,'(a,i0)') '[DIRECTmf] compare: column count mismatches = ', nbad_cnt
    nerr = nbad_tree + nbad_cnt

    deallocate(cnt)
    call hecmw_mf_symbolic_finalize(ref)
    deallocate(fct%IROw, fct%JCOl, fct%IPErm, fct%INVp, fct%PARent, fct%NCH, fct%XLNzr, fct%COLno)
  end subroutine hecmw_mf_compare_builtin

end module hecmw_solver_direct_mf
