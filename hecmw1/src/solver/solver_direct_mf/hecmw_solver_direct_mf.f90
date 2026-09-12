!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Multifrontal direct solver (METHOD=DIRECTmf). Only the symbolic stage exists so far;
!>        the numeric factorization and the solve are delegated to the built-in direct solver.
module hecmw_solver_direct_mf
  use hecmw_util
  use hecmw_matrix_misc
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_solver_direct
  use hecmw_solver_direct_parallel
  implicit none

  private
  public :: hecmw_solve_direct_mf
  public :: hecmw_mf_compare_builtin

  !> relaxed amalgamation target in columns (DOFs); converted to nodes with the block size
  integer(kind=kint), parameter :: MF_RELAX_COLS = 16

  type(hecmwST_mf_symbolic), save :: SYM

contains

  subroutine hecmw_solve_direct_mf(hecMESH, hecMAT, imsg)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type(hecmwST_matrix), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: imsg
    type(hecmwST_mf_graph) :: graph
    integer(kind=kint), allocatable :: perm(:), invp(:)
    integer(kind=kint) :: loglevel, ordering, n, nerr, relax
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
      t2 = hecmw_wtime()
      if (loglevel > 0) then
        write(*,'(a,f10.3,a)') '[DIRECTmf]: symbolic fct done (', t2 - t1, ' sec)'
        call hecmw_mf_symbolic_print(SYM)
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
    else
      call hecmw_solve_direct(hecMESH, hecMAT, imsg)
    endif
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
