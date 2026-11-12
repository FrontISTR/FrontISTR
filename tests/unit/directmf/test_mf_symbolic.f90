!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> Unit checks of hecmw_mf_symbolic on synthetic grid graphs.
!> Links against the hecmw library; see CMakeLists.txt.
!>
!> For each case the structure is verified by hecmw_mf_symbolic_check, and (with relax=0) the
!> elimination tree, column counts and factor size are compared with an independent reference
!> symbolic factorization (George-Liu row merge) and with the built-in DIRECT solver.
program test_mf_symbolic
  use hecmw_util
  use hecmw_ordering
  use hecmw_mf_graph
  use hecmw_mf_symbolic
  use hecmw_solver_direct_mf
  implicit none

  integer(kind=kint) :: nfail, ntest
  integer(kind=8) :: seed

  nfail = 0
  ntest = 0
  seed = 12345

  call run_case('grid2d_9pt_d1_id',    2, 12, 10,  1, 1, 'identity', 0)
  call run_case('grid2d_9pt_d1_rand',  2, 12, 10,  1, 1, 'random',   0)
  call run_case('grid2d_9pt_d3_metis', 2, 20, 17,  1, 3, 'metis',    0)
  call run_case('grid2d_9pt_d3_relax', 2, 20, 17,  1, 3, 'metis',    5)
  call run_case('grid2d_9pt_mix_rand', 2, 15, 13,  1, 0, 'random',   0)
  call run_case('grid2d_9pt_mix_relax',2, 15, 13,  1, 0, 'random',   4)
  call run_case('grid3d_27pt_d1_rev',  3,  7,  6,  5, 1, 'reverse',  0)
  call run_case('grid3d_27pt_d3_metis',3,  9,  8,  7, 3, 'metis',    0)
  call run_case('grid3d_27pt_d3_relax',3,  9,  8,  7, 3, 'metis',    5)
  call run_case('grid3d_27pt_mix_metis',3, 8,  7,  6, 0, 'metis',    0)
  call run_case('chain_d2_id',         1, 50,  1,  1, 2, 'identity', 0)
  call run_case('chain_d2_relax',      1, 50,  1,  1, 2, 'identity', 8)
  call run_case('single_node',         1,  1,  1,  1, 3, 'identity', 0)

  write(*,'(a,i0,a,i0)') 'test_mf_symbolic: failed ', nfail, ' of ', ntest
  if (nfail > 0) stop 1

contains

  !> dim: 1 chain, 2 grid with 9-point (Q4) stencil, 3 grid with 27-point (HEX8) stencil
  !> ndof_mode: >0 constant block size, 0 mixed sizes 1..3
  subroutine run_case(name, dim, nx, ny, nz, ndof_mode, ordering, relax)
    character(len=*), intent(in) :: name, ordering
    integer(kind=kint), intent(in) :: dim, nx, ny, nz, ndof_mode, relax
    type(hecmwST_mf_graph) :: graph
    type(hecmwST_mf_symbolic) :: sym
    integer(kind=kint), allocatable :: ndof(:), xadj(:), adjncy(:), perm(:), invp(:)
    integer(kind=kint), allocatable :: ref_parent(:), ref_colcnt(:)
    integer(kind=8) :: ref_nnz
    integer(kind=kint) :: n, i, k, nerr, nbad, mp
    logical :: ok

    ntest = ntest + 1
    call make_grid(dim, nx, ny, nz, n, xadj, adjncy)
    allocate(ndof(n), perm(n), invp(n))
    do i = 1, n
      if (ndof_mode > 0) then
        ndof(i) = ndof_mode
      else
        ndof(i) = 1 + mod(i, 3)
      endif
    enddo
    call hecmw_mf_graph_set(graph, n, ndof, xadj, adjncy)

    select case (ordering)
    case ('identity')
      do i = 1, n
        perm(i) = i
      enddo
    case ('reverse')
      do i = 1, n
        perm(i) = n + 1 - i
      enddo
    case ('random')
      call shuffle(n, perm)
    case ('metis')
      call hecmw_ordering_gen(n, n + (xadj(n+1)-1)/2, xadj, adjncy, perm, invp, 0, 0)
    end select

    call hecmw_mf_symbolic_build(graph, perm, relax, sym)
    ok = .true.

    call hecmw_mf_symbolic_check(graph, sym, nerr)
    if (nerr /= 0) then
      write(*,'(a,a,a,i0)') '  ', name, ': self check violations = ', nerr
      ok = .false.
    endif

    ! independent reference in the original node labels
    allocate(ref_parent(n), ref_colcnt(n))
    call reference_symbolic(graph, perm, ref_parent, ref_colcnt, ref_nnz)
    nbad = 0
    do i = 1, n
      k = sym%invp(i)
      mp = sym%parent(k)
      if (mp == 0) then
        if (ref_parent(i) /= 0) nbad = nbad + 1
      else
        if (sym%perm(mp) /= ref_parent(i)) nbad = nbad + 1
      endif
      if (sym%colcnt(k) /= ref_colcnt(i)) nbad = nbad + 1
    enddo
    if (nbad /= 0) then
      write(*,'(a,a,a,i0)') '  ', name, ': tree/colcnt mismatches vs reference = ', nbad
      ok = .false.
    endif
    if (relax == 0 .and. sym%factor_nnz /= ref_nnz) then
      write(*,'(a,a,a,i0,a,i0)') '  ', name, ': factor_nnz ', sym%factor_nnz, ' vs reference ', ref_nnz
      ok = .false.
    endif
    if (relax > 0 .and. sym%factor_nnz < ref_nnz) then
      write(*,'(a,a,a)') '  ', name, ': relaxed factor_nnz smaller than exact'
      ok = .false.
    endif

    if (ndof_mode > 0) then
      call compare_builtin(graph, ndof_mode, nerr)
      if (nerr /= 0) then
        write(*,'(a,a,a,i0)') '  ', name, ': mismatches vs built-in DIRECT = ', nerr
        ok = .false.
      endif
    endif

    if (ok) then
      write(*,'(a,a,a,i0,a,i0,a,i0,a,i0)') 'PASS ', name, ': n=', n, ' nsuper=', sym%nsuper, &
        ' merges=', sym%nmerge, ' factor_nnz=', sym%factor_nnz
    else
      write(*,'(a,a)') 'FAIL ', name
      nfail = nfail + 1
    endif

    call hecmw_mf_symbolic_finalize(sym)
    call hecmw_mf_graph_finalize(graph)
    deallocate(ndof, xadj, adjncy, perm, invp, ref_parent, ref_colcnt)
  end subroutine run_case

  subroutine make_grid(dim, nx, ny, nz, n, xadj, adjncy)
    integer(kind=kint), intent(in) :: dim, nx, ny, nz
    integer(kind=kint), intent(out) :: n
    integer(kind=kint), allocatable, intent(out) :: xadj(:), adjncy(:)
    integer(kind=kint) :: ix, iy, iz, jx, jy, jz, i, j, cnt, my, mz, pass

    my = ny
    mz = nz
    if (dim < 2) my = 1
    if (dim < 3) mz = 1
    n = nx*my*mz
    allocate(xadj(n+1))
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
            xadj(i+1) = cnt + 1
          enddo
        enddo
      enddo
      if (pass == 1) allocate(adjncy(cnt))
    enddo
  end subroutine make_grid

  subroutine shuffle(n, perm)
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(out) :: perm(:)
    integer(kind=kint) :: i, j, t

    do i = 1, n
      perm(i) = i
    enddo
    do i = n, 2, -1
      seed = mod(1103515245_8*seed + 12345_8, 2147483648_8)
      j = 1 + int(mod(seed, int(i, 8)), kind=kint)
      t = perm(i)
      perm(i) = perm(j)
      perm(j) = t
    enddo
  end subroutine shuffle

  !> George-Liu symbolic factorization: struct(L(:,j)) = adj_upper(j) U (U_{children c} struct(L(:,c)) \ {c});
  !> parent(c) = min struct(L(:,c)) \ {c}. Results are returned in the original node labels.
  subroutine reference_symbolic(graph, perm, parent, colcnt, nnz)
    type(hecmwST_mf_graph), intent(in) :: graph
    integer(kind=kint), intent(in) :: perm(:)
    integer(kind=kint), intent(out) :: parent(:), colcnt(:)
    integer(kind=8), intent(out) :: nnz
    integer(kind=kint), allocatable :: invp(:), sptr(:), slist(:), head(:), next(:), mark(:), buf(:), par(:)
    integer(kind=kint) :: n, j, k, c, p, i, cnt, cap, minrow
    integer(kind=8) :: rowdof

    n = graph%nnode
    allocate(invp(n), sptr(n+1), head(n), next(n), mark(n), buf(n), par(n))
    do k = 1, n
      invp(perm(k)) = k
    enddo
    cap = graph%xadj(n+1)
    allocate(slist(cap))
    head(1:n) = 0
    mark(1:n) = 0
    sptr(1) = 1
    nnz = 0
    do j = 1, n
      cnt = 0
      mark(j) = j
      do p = graph%xadj(perm(j)), graph%xadj(perm(j)+1)-1
        i = invp(graph%adjncy(p))
        if (i > j .and. mark(i) /= j) then
          mark(i) = j
          cnt = cnt + 1
          buf(cnt) = i
        endif
      enddo
      c = head(j)
      do while (c /= 0)
        do k = sptr(c), sptr(c+1)-1
          i = slist(k)
          if (i > j .and. mark(i) /= j) then
            mark(i) = j
            cnt = cnt + 1
            buf(cnt) = i
          endif
        enddo
        c = next(c)
      enddo
      if (sptr(j) + cnt - 1 > cap) then
        cap = max(2*cap, sptr(j) + cnt - 1)
        call grow(slist, cap)
      endif
      slist(sptr(j):sptr(j)+cnt-1) = buf(1:cnt)
      sptr(j+1) = sptr(j) + cnt
      minrow = 0
      rowdof = 0
      do k = 1, cnt
        if (minrow == 0 .or. buf(k) < minrow) minrow = buf(k)
        rowdof = rowdof + graph%ndof(perm(buf(k)))
      enddo
      par(j) = minrow
      if (minrow /= 0) then
        next(j) = head(minrow)
        head(minrow) = j
      endif
      colcnt(perm(j)) = cnt + 1
      nnz = nnz + int(graph%ndof(perm(j)),8)*rowdof + int(graph%ndof(perm(j)),8)*(graph%ndof(perm(j))+1)/2
    enddo
    do j = 1, n
      if (par(j) == 0) then
        parent(perm(j)) = 0
      else
        parent(perm(j)) = perm(par(j))
      endif
    enddo
    deallocate(invp, sptr, slist, head, next, mark, buf, par)
  end subroutine reference_symbolic

  subroutine grow(a, cap)
    integer(kind=kint), allocatable, intent(inout) :: a(:)
    integer(kind=kint), intent(in) :: cap
    integer(kind=kint), allocatable :: tmp(:)

    allocate(tmp(cap))
    tmp(1:size(a)) = a(1:size(a))
    call move_alloc(tmp, a)
  end subroutine grow

  !> Assemble a hecmwST_matrix pattern from the graph and run hecmw_mf_compare_builtin on it.
  subroutine compare_builtin(graph, ndof, nerr)
    type(hecmwST_mf_graph), intent(in) :: graph
    integer(kind=kint), intent(in) :: ndof
    integer(kind=kint), intent(out) :: nerr
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix) :: hecMAT
    integer(kind=kint) :: n, i, p, j, nl, nu

    n = graph%nnode
    hecMESH%N_DOF = ndof
    hecMAT%N = n
    hecMAT%NP = n
    hecMAT%NDOF = ndof
    allocate(hecMAT%indexL(0:n), hecMAT%indexU(0:n))
    allocate(hecMAT%itemL(graph%xadj(n+1)-1), hecMAT%itemU(graph%xadj(n+1)-1))
    hecMAT%indexL(0) = 0
    hecMAT%indexU(0) = 0
    nl = 0
    nu = 0
    do i = 1, n
      do p = graph%xadj(i), graph%xadj(i+1)-1
        j = graph%adjncy(p)
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
    call hecmw_mf_compare_builtin(hecMESH, hecMAT, graph, 0, nerr)
    deallocate(hecMAT%indexL, hecMAT%indexU, hecMAT%itemL, hecMAT%itemU)
  end subroutine compare_builtin

end program test_mf_symbolic
