!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : redundant sparse DIRECTmf coarsest solver
!!
!! Third direct coarsest-level backend: the same redundant data flow as the dense
!! LDL^T backend (coarse operator gathered to every rank, RHS allgathered per
!! V-cycle), but factored SPARSE by the built-in multifrontal solver (DIRECTmf)
!! through its graph / symbolic / numeric API.  No external library is needed, so
!! a MUMPS-less build gets a sparse coarsest and the dense O(n^3) size limit does
!! not apply; a symmetric-indefinite coarse operator is factored LDL^T with
!! threshold pivoting (inertia reported like the other backends) and a
!! non-symmetric one (friction contact) LU.  A numeric-only refresh re-factors
!! the new values on the retained symbolic structure.
module hecmw_precond_SAAMG_coarse_mf
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_comm, only: hecmw_saamg_abort, hecmw_saamg_check_alloc
  use hecmw_util,        only: hecmwST_matrix
  use hecmw_ordering,    only: hecmw_ordering_gen
  use hecmw_mf_graph,    only: hecmwST_mf_graph, hecmw_mf_graph_from_hecmat, hecmw_mf_graph_finalize
  use hecmw_mf_symbolic, only: hecmwST_mf_symbolic, hecmw_mf_symbolic_build, hecmw_mf_symbolic_finalize
  use hecmw_mf_numeric,  only: hecmwST_mf_factor, hecmw_mf_numeric_init, hecmw_mf_numeric_factor, &
       hecmw_mf_numeric_solve, hecmw_mf_numeric_finalize
  implicit none

  private
  public :: hecmwST_saamg_cmf
  public :: hecmw_saamg_cmf_setup
  public :: hecmw_saamg_cmf_refresh
  public :: hecmw_saamg_cmf_solve
  public :: hecmw_saamg_cmf_free

  !> DIRECTmf driver defaults, mirrored here because this backend calls the
  !! factorization API directly (there is no solver option line on this path)
  integer(kind=kint), parameter :: CMF_RELAX_COLS = 16  !< relaxed amalgamation target in DOFs
  integer(kind=kint), parameter :: CMF_TILE = 256       !< factor tile size in DOFs

  !> Redundant sparse DIRECTmf coarsest solver state: the gathered coarse operator
  !! as a synthetic NDOF=m block matrix plus the persistent symbolic structure and
  !! factor (the symbolic part is reused by the numeric-only refresh).
  type hecmwST_saamg_cmf
    logical            :: ready = .false.
    logical            :: symmetric = .true.  !< .true. = LDL^T (threshold pivoting), .false. = LU
    integer(kind=kint) :: n = 0
    integer(kind=kint) :: n_neg = 0           !< # negative eigenvalues (inertia; LDL^T only)
    type(hecmwST_matrix)      :: mat          !< pointers defined only while ready
    type(hecmwST_mf_symbolic) :: symb
    type(hecmwST_mf_factor)   :: fct
  end type hecmwST_saamg_cmf

contains

  !> Convert the redundantly gathered coarse operator Ac to the block matrix layout,
  !! run the symbolic analysis (graph, fill-reducing ordering, supernode structure)
  !! and factor: symmetric -> LDL^T, else -> LU (mode follows mat%symmetric).
  subroutine hecmw_saamg_cmf_setup(Ac, symmetric, cmf)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)    :: Ac
    logical,                  intent(in)    :: symmetric
    type(hecmwST_saamg_cmf),  intent(inout) :: cmf
    logical :: ok

    call hecmw_saamg_cmf_free(cmf)
    cmf%symmetric = symmetric
    cmf%n = Ac%n
    call cmf_build_pattern(Ac, cmf%mat)
    cmf%mat%symmetric = symmetric
    cmf%ready = .true.
    call cmf_fill_values(Ac, cmf%mat, ok)   ! ok guaranteed: the pattern was built from this Ac
    call cmf_analyze(cmf)
    call cmf_factor(cmf)
  end subroutine hecmw_saamg_cmf_setup

  !> Numeric refresh: reload the values of the re-gathered operator into the stored
  !! pattern and re-run the numeric factorization on the retained symbolic structure.
  !! The gather drops numerically zero entries, so a refreshed operator can (rarely)
  !! carry a block outside the stored pattern; then the symbolic structure is stale
  !! and the backend falls back to a full re-setup.
  subroutine hecmw_saamg_cmf_refresh(Ac, cmf)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)    :: Ac
    type(hecmwST_saamg_cmf),  intent(inout) :: cmf
    logical :: ok

    call cmf_fill_values(Ac, cmf%mat, ok)
    if (.not. ok) then
      call hecmw_saamg_cmf_setup(Ac, cmf%symmetric, cmf)
      return
    end if
    call cmf_factor(cmf)
  end subroutine hecmw_saamg_cmf_refresh

  !> Solve  Ac x = b  with the stored factor (plain-array interface, every rank
  !! solves the full gathered system redundantly, like the dense backend).
  subroutine hecmw_saamg_cmf_solve(cmf, b, x)
    implicit none
    type(hecmwST_saamg_cmf), intent(inout) :: cmf
    real(kind=kreal),        intent(in)    :: b(:)
    real(kind=kreal),        intent(out)   :: x(:)

    call hecmw_mf_numeric_solve(cmf%symb, cmf%fct, b, x)
  end subroutine hecmw_saamg_cmf_solve

  subroutine hecmw_saamg_cmf_free(cmf)
    implicit none
    type(hecmwST_saamg_cmf), intent(inout) :: cmf

    call hecmw_mf_numeric_finalize(cmf%fct)
    call hecmw_mf_symbolic_finalize(cmf%symb)
    ! hecmwST_matrix components are pointers with no default initialization, so they
    ! are touched only when this instance was actually set up (cf. the MUMPS backend)
    if (cmf%ready) then
      deallocate(cmf%mat%D, cmf%mat%AL, cmf%mat%AU)
      deallocate(cmf%mat%indexL, cmf%mat%indexU, cmf%mat%itemL, cmf%mat%itemU)
    end if
    cmf%ready = .false.; cmf%n = 0; cmf%n_neg = 0
  end subroutine hecmw_saamg_cmf_free

  !> Build the block-CSR pattern of the synthetic matrix from Ac's block pattern,
  !! SYMMETRIZED as the union of (i,j) and (j,i): the gather drops numerically zero
  !! entries, so one triangle of a block pair can be missing from Ac, while the
  !! multifrontal solver requires a structurally symmetric input (its mirror maps
  !! pair every AL block with an AU block).  Missing counterparts stay zero blocks.
  subroutine cmf_build_pattern(Ac, mat)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)    :: Ac
    type(hecmwST_matrix),     intent(inout) :: mat
    integer(kind=kint), allocatable :: tptr(:), tcur(:), tlist(:), mark(:), cols(:)
    integer(kind=kint) :: n, m, i, j, t, q, w, nc, nl, nu, astat

    n = Ac%nbrow; m = Ac%nb
    mat%N = n; mat%NP = n; mat%NDOF = m; mat%NPA = 0
    mat%Iarray = 0; mat%Rarray = 0.0d0
    ! transpose block pattern: tlist(tptr(j):tptr(j+1)-1) = rows i with (i,j) stored
    allocate(tptr(n+1), tcur(n), tlist(Ac%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf pattern (transpose)')
    tptr(1:n+1) = 0
    do i = 1, n
      do t = Ac%browptr(i), Ac%browptr(i+1)-1
        j = Ac%bcol(t); tptr(j+1) = tptr(j+1) + 1
      end do
    end do
    tptr(1) = 1
    do j = 1, n
      tptr(j+1) = tptr(j+1) + tptr(j)
      tcur(j) = tptr(j)
    end do
    do i = 1, n
      do t = Ac%browptr(i), Ac%browptr(i+1)-1
        j = Ac%bcol(t); tlist(tcur(j)) = i; tcur(j) = tcur(j) + 1
      end do
    end do

    allocate(mark(n), cols(n), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf pattern (row union)')
    mark(1:n) = 0
    allocate(mat%indexL(0:n), mat%indexU(0:n), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf pattern (indexL/indexU)')
    mat%indexL(0) = 0; mat%indexU(0) = 0
    do i = 1, n                                  ! pass 1: count the L/U unions
      call cmf_row_union(Ac, tptr, tlist, mark, cols, i, nc)
      nl = 0; nu = 0
      do q = 1, nc
        if (cols(q) < i) nl = nl + 1
        if (cols(q) > i) nu = nu + 1
        mark(cols(q)) = 0
      end do
      mat%indexL(i) = mat%indexL(i-1) + nl
      mat%indexU(i) = mat%indexU(i-1) + nu
    end do
    mat%NPL = mat%indexL(n); mat%NPU = mat%indexU(n)
    allocate(mat%itemL(max(mat%NPL,1)), mat%itemU(max(mat%NPU,1)), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf pattern (itemL/itemU)')
    do i = 1, n                                  ! pass 2: fill sorted item lists
      call cmf_row_union(Ac, tptr, tlist, mark, cols, i, nc)
      do q = 2, nc                               ! insertion sort (rows are short)
        j = cols(q); w = q - 1
        do while (w >= 1)
          if (cols(w) <= j) exit
          cols(w+1) = cols(w); w = w - 1
        end do
        cols(w+1) = j
      end do
      nl = mat%indexL(i-1); nu = mat%indexU(i-1)
      do q = 1, nc
        j = cols(q); mark(j) = 0
        if (j < i) then
          nl = nl + 1; mat%itemL(nl) = j
        else if (j > i) then
          nu = nu + 1; mat%itemU(nu) = j
        end if
      end do
    end do
    deallocate(tptr, tcur, tlist, mark, cols)
    allocate(mat%D(n*m*m), mat%AL(max(mat%NPL,1)*m*m), mat%AU(max(mat%NPU,1)*m*m), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf coarse matrix values (D/AL/AU)')
  end subroutine cmf_build_pattern

  !> Collect the unmarked block columns of row i, from Ac's row and its transpose
  !! (the symmetrizing union).  Leaves the collected columns marked; the caller
  !! resets the marks from cols(1:nc).
  subroutine cmf_row_union(Ac, tptr, tlist, mark, cols, i, nc)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)    :: Ac
    integer(kind=kint),       intent(in)    :: tptr(:), tlist(:), i
    integer(kind=kint),       intent(inout) :: mark(:), cols(:)
    integer(kind=kint),       intent(out)   :: nc
    integer(kind=kint) :: t, q, j

    nc = 0
    do t = Ac%browptr(i), Ac%browptr(i+1)-1
      j = Ac%bcol(t)
      if (mark(j) == 0) then
        mark(j) = 1; nc = nc + 1; cols(nc) = j
      end if
    end do
    do q = tptr(i), tptr(i+1)-1
      j = tlist(q)
      if (mark(j) == 0) then
        mark(j) = 1; nc = nc + 1; cols(nc) = j
      end if
    end do
  end subroutine cmf_row_union

  !> Scatter Ac's block values into the stored pattern (zeros first, so union blocks
  !! without a stored counterpart stay zero).  Blocks are transposed on the fly:
  !! bcsr blocks are column major, hecmwST_matrix blocks row major.  ok = .false.
  !! when Ac carries a block outside the stored pattern (the pattern changed).
  subroutine cmf_fill_values(Ac, mat, ok)
    implicit none
    type(hecmwST_saamg_bcsr), intent(in)    :: Ac
    type(hecmwST_matrix),     intent(inout) :: mat
    logical,                  intent(out)   :: ok
    integer(kind=kint), allocatable :: pos(:)
    integer(kind=kint) :: n, m, mm, i, j, t, l, u, r, c, b0, d0, astat

    n = mat%NP; m = mat%NDOF; mm = m*m
    ok = .true.
    allocate(pos(n), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf value scatter (pos)')
    pos(1:n) = 0
    mat%D(1:n*mm) = 0.0d0
    if (mat%NPL > 0) mat%AL(1:mat%NPL*mm) = 0.0d0
    if (mat%NPU > 0) mat%AU(1:mat%NPU*mm) = 0.0d0
    do i = 1, n
      do l = mat%indexL(i-1)+1, mat%indexL(i)
        pos(mat%itemL(l)) = l
      end do
      do u = mat%indexU(i-1)+1, mat%indexU(i)
        pos(mat%itemU(u)) = -u
      end do
      do t = Ac%browptr(i), Ac%browptr(i+1)-1
        j = Ac%bcol(t); b0 = (t-1)*mm
        if (j == i) then
          d0 = (i-1)*mm
          do c = 1, m
            do r = 1, m
              mat%D(d0 + (r-1)*m + c) = Ac%bval(b0 + (c-1)*m + r)
            end do
          end do
        else if (pos(j) > 0) then
          d0 = (pos(j)-1)*mm
          do c = 1, m
            do r = 1, m
              mat%AL(d0 + (r-1)*m + c) = Ac%bval(b0 + (c-1)*m + r)
            end do
          end do
        else if (pos(j) < 0) then
          d0 = (-pos(j)-1)*mm
          do c = 1, m
            do r = 1, m
              mat%AU(d0 + (r-1)*m + c) = Ac%bval(b0 + (c-1)*m + r)
            end do
          end do
        else
          ok = .false.
        end if
      end do
      do l = mat%indexL(i-1)+1, mat%indexL(i)
        pos(mat%itemL(l)) = 0
      end do
      do u = mat%indexU(i-1)+1, mat%indexU(i)
        pos(mat%itemU(u)) = 0
      end do
    end do
    deallocate(pos)
  end subroutine cmf_fill_values

  !> Symbolic stage: node graph, fill-reducing ordering (0 = solver default:
  !! METIS if built, else QMD), supernode structure, factor storage layout.
  subroutine cmf_analyze(cmf)
    implicit none
    type(hecmwST_saamg_cmf), intent(inout) :: cmf
    type(hecmwST_mf_graph) :: graph
    integer(kind=kint), allocatable :: perm(:), invp(:)
    integer(kind=kint) :: n, relax, astat

    call hecmw_mf_graph_from_hecmat(cmf%mat, graph)
    n = graph%nnode
    allocate(perm(n), invp(n), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'cmf ordering (perm/invp)')
    call hecmw_ordering_gen(n, n + (graph%xadj(n+1)-1)/2, graph%xadj, graph%adjncy, perm, invp, 0, 0)
    relax = max(1, CMF_RELAX_COLS / cmf%mat%NDOF)
    call hecmw_mf_symbolic_finalize(cmf%symb)
    call hecmw_mf_symbolic_build(graph, perm, relax, cmf%symb)
    call hecmw_mf_numeric_finalize(cmf%fct)
    call hecmw_mf_numeric_init(graph, cmf%symb, CMF_TILE, .false., 0, cmf%fct)
    deallocate(perm, invp)
    call hecmw_mf_graph_finalize(graph)
  end subroutine cmf_analyze

  !> Numeric factorization with the DIRECTmf defaults (no BLR, default pivot
  !! thresholds); LDL^T / LU follows mat%symmetric.  The inertia of an LDL^T
  !! factorization is reported like the dense and MUMPS backends.
  subroutine cmf_factor(cmf)
    implicit none
    type(hecmwST_saamg_cmf), intent(inout) :: cmf
    integer(kind=kint) :: ierr

    call hecmw_mf_numeric_factor(cmf%mat, cmf%symb, cmf%fct, ierr)
    if (ierr /= 0) then
      if (ierr > 0) then
        write(*,'(a,i0)') 'hecmw_saamg_cmf: DIRECTmf zero pivot at coarse dof ', cmf%fct%pdof(ierr)
        call hecmw_saamg_abort('cmf coarsest: factorization failed (singular coarse operator)')
      else
        write(*,'(a,i0)') 'hecmw_saamg_cmf: DIRECTmf structural error, ierr=', ierr
        call hecmw_saamg_abort('cmf coarsest: coarse matrix does not match the symbolic structure')
      end if
    end if
    if (cmf%fct%lu) then
      cmf%n_neg = 0
    else
      cmf%n_neg = cmf%fct%n_neg
    end if
  end subroutine cmf_factor

end module hecmw_precond_SAAMG_coarse_mf
