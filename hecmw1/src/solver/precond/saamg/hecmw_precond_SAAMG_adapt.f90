!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
!> \brief  Smoothed Aggregation AMG preconditioner : FrontISTR adapter
!!
!! Bridges FrontISTR data structures to the mesh-independent SA-AMG core:
!!   * converts the BCRS hecmwST_matrix (separate D / AL / AU node blocks) into
!!     the core's internal block-CSR (include_halo keeps the rectangular couplings
!!     to halo columns for the distributed operator);
!!   * builds the finest-level communication table from the mesh import/export.
!! The rigid-body near-kernel is built in the backend (build_global_rbm) because it
!! needs the global-centroid allreduce for rank consistency.
module hecmw_precond_SAAMG_adapt
  use hecmw_util
  use hecmw_precond_SAAMG_matrix, only: hecmwST_saamg_bcsr, hecmw_saamg_bcsr_free
  use hecmw_precond_SAAMG_comm,   only: hecmwST_saamg_comm, hecmw_saamg_comm_free, &
       hecmw_saamg_check_alloc
  implicit none

  private
  public :: hecmw_saamg_from_hecmat
  public :: hecmw_saamg_comm_from_mesh

contains

  !> Convert hecmwST_matrix (FrontISTR block storage) -> internal block-CSR.
  !! Rows are always internal only (1..hecMAT%N); the result has nint*ndof rows.
  !! include_halo = .false. (default): square nint*ndof matrix, halo columns
  !!   dropped -- the sequential operator.
  !! include_halo = .true.: rectangular nint*ndof x NP*ndof matrix that keeps the
  !!   couplings to halo nodes (columns nint*ndof+1..NP*ndof) -- the distributed
  !!   finest operator, used with a comm table by hecmw_saamg_matvec_d.
  subroutine hecmw_saamg_from_hecmat(hecMAT, A, include_halo)
    implicit none
    type(hecmwST_matrix),  intent(in)  :: hecMAT
    type(hecmwST_saamg_bcsr), intent(out) :: A
    logical, optional,     intent(in)  :: include_halo
    integer(kind=kint) :: ndof, ndof2, nnode, ncolnode, i, ii, kk, j, jn, pos, cnt, b0, astat
    logical :: halo

    halo = .false.
    if (present(include_halo)) halo = include_halo

    ndof  = hecMAT%NDOF
    ndof2 = ndof*ndof
    nnode = hecMAT%N                     ! internal rows
    if (halo) then
      ncolnode = hecMAT%NP               ! columns reference internal + halo
    else
      ncolnode = hecMAT%N
    end if

    ! Build the block-CSR directly from the hecMAT BCRS (diagonal D + lower AL +
    ! upper AU, one nb x nb block per neighbor) WITHOUT an intermediate triplet
    ! buffer -- avoids holding the operator twice (critical for large meshes).
    ! hecMAT stores each block row-major ((ii-1)*ndof+kk); transpose to col-major.
    call hecmw_saamg_bcsr_free(A)
    A%n = nnode*ndof; A%ncol = ncolnode*ndof; A%nb = ndof; A%mb = ndof
    A%nbrow = nnode; A%nbcol = ncolnode
    allocate(A%browptr(nnode+1))

    ! pass 1: count blocks per row (diagonal + kept L + kept U)
    A%browptr(1) = 1
    do i = 1, nnode
      cnt = 1                                       ! diagonal block
      do j = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        if (halo .or. hecMAT%itemL(j) <= nnode) cnt = cnt + 1
      end do
      do j = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        if (halo .or. hecMAT%itemU(j) <= nnode) cnt = cnt + 1
      end do
      A%browptr(i+1) = A%browptr(i) + cnt
    end do
    A%nnzb = A%browptr(nnode+1) - 1
    A%nnz  = ndof2 * A%nnzb
    allocate(A%bcol(A%nnzb), A%bval(ndof2*A%nnzb), stat=astat)
    call hecmw_saamg_check_alloc(astat, 'from_hecmat finest operator (bcol/bval)')

    ! pass 2: fill blocks in place (diagonal first, then L, then U)
    do i = 1, nnode
      pos = A%browptr(i)
      A%bcol(pos) = i; b0 = (pos-1)*ndof2          ! diagonal block
      do ii = 1, ndof
        do kk = 1, ndof
          A%bval(b0 + (kk-1)*ndof+ii) = hecMAT%D(ndof2*(i-1) + (ii-1)*ndof + kk)
        end do
      end do
      pos = pos + 1
      do j = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        jn = hecMAT%itemL(j)
        if (.not. halo .and. jn > nnode) cycle
        A%bcol(pos) = jn; b0 = (pos-1)*ndof2
        do ii = 1, ndof
          do kk = 1, ndof
            A%bval(b0 + (kk-1)*ndof+ii) = hecMAT%AL(ndof2*(j-1) + (ii-1)*ndof + kk)
          end do
        end do
        pos = pos + 1
      end do
      do j = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        jn = hecMAT%itemU(j)
        if (.not. halo .and. jn > nnode) cycle
        A%bcol(pos) = jn; b0 = (pos-1)*ndof2
        do ii = 1, ndof
          do kk = 1, ndof
            A%bval(b0 + (kk-1)*ndof+ii) = hecMAT%AU(ndof2*(j-1) + (ii-1)*ndof + kk)
          end do
        end do
        pos = pos + 1
      end do
    end do
  end subroutine hecmw_saamg_from_hecmat

  !> Build the finest-level communication table from the FrontISTR mesh.
  !! Copies the node-based halo descriptor (neighbor PEs, import/export CSR
  !! lists) verbatim; these operate on nodes with nb = ndof dofs each.
  subroutine hecmw_saamg_comm_from_mesh(hecMESH, ndof, cmt)
    implicit none
    type(hecmwST_local_mesh), intent(in)  :: hecMESH
    integer(kind=kint),       intent(in)  :: ndof
    type(hecmwST_saamg_comm),   intent(out) :: cmt
    integer(kind=kint) :: nnb, nimp, nexp

    call hecmw_saamg_comm_free(cmt)
    cmt%comm       = hecMESH%MPI_COMM
    cmt%my_rank    = hecMESH%my_rank
    cmt%nint       = hecMESH%nn_internal
    cmt%nnode      = hecMESH%n_node
    cmt%nb         = ndof
    cmt%n_neighbor = hecMESH%n_neighbor_pe
    ! global node id per local node (for partition-invariant Lanczos seeds)
    allocate(cmt%gnode(hecMESH%n_node))
    cmt%gnode(1:hecMESH%n_node) = hecMESH%global_node_ID(1:hecMESH%n_node)
    nnb = cmt%n_neighbor
    if (nnb == 0) return

    nimp = hecMESH%import_index(nnb)
    nexp = hecMESH%export_index(nnb)
    allocate(cmt%neighbor(nnb))
    allocate(cmt%import_index(0:nnb), cmt%import_item(nimp))
    allocate(cmt%export_index(0:nnb), cmt%export_item(nexp))
    cmt%neighbor(1:nnb)        = hecMESH%neighbor_pe(1:nnb)
    cmt%import_index(0:nnb)    = hecMESH%import_index(0:nnb)
    cmt%export_index(0:nnb)    = hecMESH%export_index(0:nnb)
    if (nimp > 0) cmt%import_item(1:nimp) = hecMESH%import_item(1:nimp)
    if (nexp > 0) cmt%export_item(1:nexp) = hecMESH%export_item(1:nexp)
  end subroutine hecmw_saamg_comm_from_mesh

end module hecmw_precond_SAAMG_adapt
