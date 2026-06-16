!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : distributed MUMPS coarsest solver
!!
!! Alternative coarsest-level backend (F7): a parallel sparse direct factorization
!! of the DISTRIBUTED coarse operator, with NO redundant gather.  Each rank feeds
!! its owned coarse rows to MUMPS as distributed assembled input (ICNTL(18)=3, lower
!! triangle, SYM=2) over the level communicator; MUMPS factors in parallel.  Per
!! V-cycle the coarse RHS is gathered to the host, solved, and the owned solution
!! slice scattered back.  This removes the last O(P) redundant gather so a much
!! larger (sparse) coarsest may be used, keeping the hierarchy shallow.
!!
!! All MUMPS code is guarded by HECMW_WITH_MUMPS; without it this module compiles
!! to stubs (hecmw_saamg_cmumps_available() = .false.) so the standalone build and
!! MUMPS-less configurations fall back to the dense LDL^T coarsest.
module hecmw_precond_SAAMG_coarse_mumps
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_comm
#ifdef HECMW_WITH_MUMPS
  ! HEC-MW convention: gather/scatter/reduce via wrappers (MUMPS implies MPI).
  use hecmw_util,     only: hecmw_sum
  use m_hecmw_comm_f, only: hecmw_comm_size, hecmw_allgather_int_1, &
       hecmw_allreduce_I_comm, hecmw_gatherv_real, hecmw_scatterv_real
#endif
  implicit none

  private
  public :: hecmwST_saamg_cmumps
  public :: hecmw_saamg_cmumps_available
  public :: hecmw_saamg_cmumps_setup
  public :: hecmw_saamg_cmumps_refresh
  public :: hecmw_saamg_cmumps_solve
  public :: hecmw_saamg_cmumps_free

#ifdef HECMW_WITH_MUMPS
  include 'dmumps_struc.h'
#endif

  !> Distributed MUMPS coarsest solver state.  Holds the persistent factorization
  !! (between setup and the per-V-cycle solves) plus the gather/scatter layout.
  type hecmwST_saamg_cmumps
    logical            :: ready = .false.
    logical            :: symmetric = .true.                   !< SYM=2 (sym) vs SYM=0 (general)
    integer(kind=kint) :: n_global = 0, n_loc = 0, n_neg = 0
    integer(kind=kint), allocatable :: rcounts(:), rdispls(:)  !< host gather/scatter (dof)
#ifdef HECMW_WITH_MUMPS
    type(dmumps_struc) :: id
#endif
  end type hecmwST_saamg_cmumps

contains

  !> .true. only when compiled with MUMPS; the caller uses the dense coarsest else.
  logical function hecmw_saamg_cmumps_available() result(ok)
#ifdef HECMW_WITH_MUMPS
    ok = .true.
#else
    ok = .false.
#endif
  end function hecmw_saamg_cmumps_available

  !> Factor the distributed coarse operator.  Aloc = this rank's owned coarse rows
  !! (n_loc = naggr_local*m rows, columns localized; global column dof via cmt%gnode);
  !! my_off = owned coarse-node offset, m = coarse block size.  Builds the global
  !! lower-triangle distributed COO and runs MUMPS analyze+factor (JOB=4).
  subroutine hecmw_saamg_cmumps_setup(cmt, Aloc, my_off, m, naggr_local, symmetric, cm)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    type(hecmwST_saamg_bcsr), intent(in)    :: Aloc
    integer(kind=kint),     intent(in)    :: my_off, m, naggr_local
    logical,                intent(in)    :: symmetric
    type(hecmwST_saamg_cmumps), intent(inout) :: cm
#ifdef HECMW_WITH_MUMPS
    integer(kind=kint) :: nprocs, p, nloc, nzloc, gbuf(1)
    call hecmw_saamg_cmumps_free(cm)
    cm%symmetric = symmetric
    cm%n_loc = naggr_local * m
    nloc = cm%n_loc
    call hecmw_comm_size(cmt%comm, nprocs)
    gbuf(1) = nloc
    call hecmw_allreduce_I_comm(gbuf, 1, hecmw_sum, cmt%comm)
    cm%n_global = gbuf(1)
    ! host gather/scatter layout (owned dof blocks are contiguous in global order)
    allocate(cm%rcounts(nprocs), cm%rdispls(nprocs))
    call hecmw_allgather_int_1(nloc, cm%rcounts, cmt%comm)
    cm%rdispls(1) = 0
    do p = 2, nprocs; cm%rdispls(p) = cm%rdispls(p-1) + cm%rcounts(p-1); end do

    ! init MUMPS instance over the level communicator
    cm%id%COMM = cmt%comm
    cm%id%PAR  = 1
    if (symmetric) then; cm%id%SYM = 2; else; cm%id%SYM = 0; end if  ! SYM=2 sym-indef / SYM=0 general
    cm%id%JOB  = -1
    call DMUMPS(cm%id)

    ! distributed COO (global numbering): lower triangle for SYM=2, full for SYM=0
    call coo_count(Aloc, cmt, my_off, m, symmetric, nzloc)
    cm%id%NZ_loc = nzloc
    allocate(cm%id%IRN_loc(nzloc), cm%id%JCN_loc(nzloc), cm%id%A_loc(nzloc))
    call coo_fill(Aloc, cmt, my_off, m, symmetric, cm%id%IRN_loc, cm%id%JCN_loc, cm%id%A_loc)
    cm%id%N      = cm%n_global
    cm%id%ICNTL(18) = 3            ! distributed assembled matrix input

    call set_quiet_options(cm%id)

    ! host RHS/solution buffer (reused across solves); other ranks need length 1
    if (cmt%my_rank == 0) then
      allocate(cm%id%RHS(cm%n_global))
    else
      allocate(cm%id%RHS(1))
    end if

    call mumps_run(cm%id, 4, cmt%my_rank)     ! analyze + factor
    if (symmetric) then; cm%n_neg = cm%id%INFOG(12); else; cm%n_neg = 0; end if  ! inertia: SYM=2 only
    cm%ready = .true.
#else
    call cmumps_unavailable()
#endif
  end subroutine hecmw_saamg_cmumps_setup

  !> Numeric refresh: pattern unchanged, only values change.  Refill A_loc and run
  !! MUMPS numeric factorization (JOB=2, reusing the JOB=1 analysis).
  subroutine hecmw_saamg_cmumps_refresh(cmt, Aloc, my_off, m, naggr_local, cm)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    type(hecmwST_saamg_bcsr), intent(in)    :: Aloc
    integer(kind=kint),     intent(in)    :: my_off, m, naggr_local
    type(hecmwST_saamg_cmumps), intent(inout) :: cm
#ifdef HECMW_WITH_MUMPS
    integer(kind=kint) :: idummy
    idummy = naggr_local
    call coo_fill(Aloc, cmt, my_off, m, cm%symmetric, cm%id%IRN_loc, cm%id%JCN_loc, cm%id%A_loc)
    call mumps_run(cm%id, 2, cmt%my_rank)     ! numeric factorization only
    if (cm%symmetric) then; cm%n_neg = cm%id%INFOG(12); else; cm%n_neg = 0; end if
#else
    call cmumps_unavailable()
#endif
  end subroutine hecmw_saamg_cmumps_refresh

  !> Solve A_coarse x = b.  rhs_loc/sol_loc are this rank's owned slices (length
  !! n_loc, local owned order = global order on this rank).  Gathers RHS to the
  !! host, runs MUMPS solve (JOB=3), scatters the owned solution slice back.
  subroutine hecmw_saamg_cmumps_solve(cmt, rhs_loc, n_loc, cm, sol_loc)
    type(hecmwST_saamg_comm), intent(in)    :: cmt
    real(kind=kreal),       intent(in)    :: rhs_loc(:)
    integer(kind=kint),     intent(in)    :: n_loc
    type(hecmwST_saamg_cmumps), intent(inout) :: cm
    real(kind=kreal),       intent(out)   :: sol_loc(:)
#ifdef HECMW_WITH_MUMPS
    call hecmw_gatherv_real(rhs_loc, n_loc, cm%id%RHS, cm%rcounts, cm%rdispls, 0, cmt%comm)
    call mumps_run(cm%id, 3, cmt%my_rank)     ! solve (centralized RHS/solution on host)
    call hecmw_scatterv_real(cm%id%RHS, cm%rcounts, cm%rdispls, sol_loc, n_loc, 0, cmt%comm)
#else
    call cmumps_unavailable()
#endif
  end subroutine hecmw_saamg_cmumps_solve

  subroutine hecmw_saamg_cmumps_free(cm)
    type(hecmwST_saamg_cmumps), intent(inout) :: cm
#ifdef HECMW_WITH_MUMPS
    ! Only touch the dmumps_struc pointers when this instance was actually set up:
    ! for a fresh (never-set-up) instance those pointers are UNDEFINED, not null
    ! (dmumps_struc has no default initialization), so associated() is unsafe.
    if (cm%ready) then
      cm%id%JOB = -2
      call DMUMPS(cm%id)                       ! free MUMPS internal data
      if (associated(cm%id%IRN_loc)) deallocate(cm%id%IRN_loc)
      if (associated(cm%id%JCN_loc)) deallocate(cm%id%JCN_loc)
      if (associated(cm%id%A_loc))   deallocate(cm%id%A_loc)
      if (associated(cm%id%RHS))     deallocate(cm%id%RHS)
    end if
#endif
    if (allocated(cm%rcounts)) deallocate(cm%rcounts)
    if (allocated(cm%rdispls)) deallocate(cm%rdispls)
    cm%ready = .false.; cm%n_global = 0; cm%n_loc = 0; cm%n_neg = 0
  end subroutine hecmw_saamg_cmumps_free

#ifdef HECMW_WITH_MUMPS
  !> Count owned-row lower-triangle entries (global row >= global col).  Reads the
  !! block storage directly; ALL block entries are emitted (explicit zeros included)
  !! so the COO pattern is purely structural -- cmumps_refresh reuses it for a
  !! value-only refactor, which requires the pattern not to depend on values.
  subroutine coo_count(A, cmt, my_off, m, symmetric, nz)
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: my_off, m
    logical,                intent(in)  :: symmetric
    integer(kind=kint),     intent(out) :: nz
    integer(kind=kint) :: ib, t, jb, r, c, grow, gcol
    nz = 0
    do ib = 1, A%nbrow
      do t = A%browptr(ib), A%browptr(ib+1)-1
        jb = A%bcol(t)
        do c = 1, m
          gcol = (cmt%gnode(jb)-1)*m + c
          do r = 1, m
            grow = my_off*m + (ib-1)*m + r
            if (.not. symmetric .or. grow >= gcol) nz = nz + 1
          end do
        end do
      end do
    end do
  end subroutine coo_count

  !> Fill the distributed COO arrays (global numbering): lower triangle (grow>=gcol)
  !! for SYM=2, the full matrix for SYM=0.  Block storage; ALL entries (structural
  !! pattern, zeros included) so cmumps_refresh can reuse the pattern.
  subroutine coo_fill(A, cmt, my_off, m, symmetric, irn, jcn, aval)
    type(hecmwST_saamg_bcsr), intent(in)  :: A
    type(hecmwST_saamg_comm), intent(in)  :: cmt
    integer(kind=kint),     intent(in)  :: my_off, m
    logical,                intent(in)  :: symmetric
    integer(kind=kint),     intent(out) :: irn(:), jcn(:)
    real(kind=kreal),       intent(out) :: aval(:)
    integer(kind=kint) :: ib, t, jb, r, c, grow, gcol, b0, nz
    real(kind=kreal) :: v
    nz = 0
    do ib = 1, A%nbrow
      do t = A%browptr(ib), A%browptr(ib+1)-1
        jb = A%bcol(t); b0 = (t-1)*m*m
        do c = 1, m
          gcol = (cmt%gnode(jb)-1)*m + c
          do r = 1, m
            v = A%bval(b0+(c-1)*m+r)
            grow = my_off*m + (ib-1)*m + r
            if (.not. symmetric .or. grow >= gcol) then
              nz = nz + 1; irn(nz) = grow; jcn(nz) = gcol; aval(nz) = v
            end if
          end do
        end do
      end do
    end do
  end subroutine coo_fill

  !> Quiet output + ordering / robustness options (set after JOB=-1, before JOB=4).
  subroutine set_quiet_options(id)
    type(dmumps_struc), intent(inout) :: id
    id%ICNTL(1) = 6; id%ICNTL(2) = 0; id%ICNTL(3) = 0; id%ICNTL(4) = 0  ! silent
    id%ICNTL(28) = 0                ! ordering: auto seq/par
    id%ICNTL(7)  = 7                ! sequential ordering: auto
    id%ICNTL(29) = 0                ! parallel ordering: auto
    id%ICNTL(14) = 20               ! working-space relaxation (%)
    id%ICNTL(10) = 0                ! no iterative refinement (cheap coarse solve)
    id%ICNTL(22) = 0                ! in-core
  end subroutine set_quiet_options

  !> Run one MUMPS phase with the standard -9 (insufficient memory) relaxation loop.
  subroutine mumps_run(id, job, myrank)
    type(dmumps_struc), intent(inout) :: id
    integer(kind=kint), intent(in)    :: job, myrank
    integer(kind=kint) :: istat
    id%JOB = job
    do
      call DMUMPS(id)
      istat = id%INFOG(1)
      if (istat >= 0) exit
      if (istat == -9 .and. id%ICNTL(14) < 200) then
        id%ICNTL(14) = id%ICNTL(14) + 20
        if (myrank == 0) write(*,'(a,i0)') &
             ' #### SA-AMG MUMPS: increasing relaxation ICNTL(14) to ', id%ICNTL(14)
      else
        if (myrank == 0) write(*,'(a,i0,a,i0)') &
             ' #### SA-AMG MUMPS error: INFOG(1)=', istat, ' INFOG(2)=', id%INFOG(2)
        exit
      end if
    end do
  end subroutine mumps_run
#else
  subroutine cmumps_unavailable()
    call hecmw_saamg_abort('MUMPS coarsest requested but not compiled (HECMW_WITH_MUMPS)')
  end subroutine cmumps_unavailable
#endif

end module hecmw_precond_SAAMG_coarse_mumps
