!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Multifrontal direct solver (METHOD=DIRECTmf): sequential LDLt with threshold pivoting
!>        for symmetric (possibly indefinite) matrices, LU for structurally symmetric matrices
!>        with unsymmetric values. Multi-process runs are still delegated to the built-in
!>        parallel direct solver.
module hecmw_solver_direct_mf
  use hecmw_util
  use hecmw_matrix_misc
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_mf_numeric
  use hecmw_solver_direct
  use hecmw_solver_direct_parallel
  implicit none

  private
  public :: hecmw_solve_direct_mf
  public :: hecmw_mf_compare_builtin

  !> relaxed amalgamation target in columns (DOFs); converted to nodes with the block size
  integer(kind=kint), parameter :: MF_RELAX_COLS = 16
  !> target tile size in DOFs
  integer(kind=kint), parameter :: MF_TILE = 256

  type(hecmwST_mf_symbolic), save :: SYM
  type(hecmwST_mf_factor), save :: FCT

contains

  subroutine hecmw_solve_direct_mf(hecMESH, hecMAT, imsg)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: imsg
    type(hecmwST_mf_graph) :: graph
    integer(kind=kint), allocatable :: perm(:), invp(:)
    integer(kind=kint) :: loglevel, ordering, n, nerr, relax, ierr, idof
    real(kind=kreal) :: t1, t2

    loglevel = hecmw_mat_get_loglevel(hecMAT)
    if (loglevel < 0) loglevel = max(hecmw_mat_get_timelog(hecMAT), hecmw_mat_get_iterlog(hecMAT))
    if (hecmw_comm_get_rank() /= 0) loglevel = 0

    if (hecMAT%Iarray(98) == 1) then
      t1 = hecmw_wtime()
      call hecmw_mf_graph_from_hecmat(hecMAT, graph)
      n = graph%nnode
      allocate(perm(n), invp(n))
      ordering = hecMAT%Iarray(41)
      call hecmw_ordering_gen(n, n + (graph%xadj(n+1)-1)/2, graph%xadj, graph%adjncy, perm, invp, ordering, loglevel)
      relax = max(1, MF_RELAX_COLS / hecMAT%NDOF)
      call hecmw_mf_symbolic_finalize(SYM)
      call hecmw_mf_symbolic_build(graph, perm, relax, SYM)
      call hecmw_mf_numeric_finalize(FCT)
      call hecmw_mf_numeric_init(SYM, MF_TILE, FCT)
      t2 = hecmw_wtime()
      if (loglevel > 0) then
        write(*,'(a,f10.3,a)') '[DIRECTmf]: symbolic fct done (', t2 - t1, ' sec)'
        call hecmw_mf_symbolic_print(SYM)
        call hecmw_mf_numeric_print(FCT)
      endif
      if (loglevel > 1) then
        call hecmw_mf_symbolic_check(graph, SYM, nerr)
        write(*,'(a,i0)') '[DIRECTmf]: self check violations = ', nerr
        call hecmw_mf_compare_builtin(hecMESH, hecMAT, graph, ordering, nerr)
        write(*,'(a,i0)') '[DIRECTmf]: mismatches against built-in DIRECT = ', nerr
      endif
      deallocate(perm, invp)
      call hecmw_mf_graph_finalize(graph)
    endif

    if (hecMESH%PETOT > 1) then
      call hecmw_solve_direct_parallel(hecMESH, hecMAT, imsg)
      return
    endif
    hecMAT%Iarray(98) = 0

    if (hecMAT%Iarray(97) == 1) then
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
        if (FCT%lu) then
          write(*,'(a,1pe9.2,a)') '[DIRECTmf]: mode = LU (asymmetry = ', FCT%asym, ')'
          write(*,'(a,i0,a,i0,a,i0,a)') '[DIRECTmf]: row swaps = ', FCT%n_swap, ', delayed = ', FCT%n_delay, &
            ' (max front growth = ', FCT%max_growth, ' dofs)'
        else
          write(*,'(a,1pe9.2,a)') '[DIRECTmf]: mode = LDLt (asymmetry = ', FCT%asym, ')'
          write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') '[DIRECTmf]: inertia (+,-) = (', FCT%n_pos, ',', FCT%n_neg, &
            '), 2x2 pivots = ', FCT%n_2x2, ', swaps = ', FCT%n_swap, ', delayed = ', FCT%n_delay, &
            ' (max front growth = ', FCT%max_growth, ' dofs)'
        endif
        write(*,'(a,i0,a,f10.3,a,i0,a,i0)') '[DIRECTmf]: factor words = ', FCT%factor_words_act, ' (', &
          real(FCT%factor_words_act, kind=kreal)*8.0d0/1024.0d0**3, ' GB), stack peak words = ', FCT%stack_peak_act, &
          ', front words = ', FCT%front_words_act
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
  end subroutine hecmw_solve_direct_mf

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
