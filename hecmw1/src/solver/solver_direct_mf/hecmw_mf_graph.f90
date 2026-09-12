!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> @brief Compressed (node level) graph of a block sparse matrix for the multifrontal solver
module hecmw_mf_graph
  use hecmw_util
  implicit none

  private
  public :: hecmwST_mf_graph
  public :: hecmw_mf_graph_from_hecmat
  public :: hecmw_mf_graph_set
  public :: hecmw_mf_graph_finalize

  !> Node i owns the ndof(i) consecutive DOFs dofptr(i):dofptr(i+1)-1.
  !> xadj/adjncy is the symmetric adjacency without self loops (1-based CSR); the same
  !> layout hecmw_ordering_gen expects.
  type hecmwST_mf_graph
    integer(kind=kint) :: nnode = 0
    integer(kind=kint), allocatable :: ndof(:)
    integer(kind=kint), allocatable :: dofptr(:)
    integer(kind=kint), allocatable :: xadj(:)
    integer(kind=kint), allocatable :: adjncy(:)
  end type hecmwST_mf_graph

contains

  !> Build the node graph from the nonzero pattern of hecMAT (values are not referenced).
  !> All NP nodes are used, as the built-in direct solver does, so the graph covers the
  !> external nodes too.
  subroutine hecmw_mf_graph_from_hecmat(hecMAT, graph)
    implicit none
    type(hecmwST_matrix), intent(in) :: hecMAT
    type(hecmwST_mf_graph), intent(out) :: graph
    integer(kind=kint) :: n, i, k, l

    n = hecMAT%NP
    graph%nnode = n
    allocate(graph%ndof(n))
    allocate(graph%dofptr(n+1))
    allocate(graph%xadj(n+1))
    graph%ndof(1:n) = hecMAT%NDOF
    graph%dofptr(1) = 1
    graph%xadj(1) = 1
    do i = 1, n
      graph%dofptr(i+1) = graph%dofptr(i) + graph%ndof(i)
      graph%xadj(i+1) = graph%xadj(i) &
        + (hecMAT%indexL(i) - hecMAT%indexL(i-1)) + (hecMAT%indexU(i) - hecMAT%indexU(i-1))
    enddo
    allocate(graph%adjncy(graph%xadj(n+1)-1))
    l = 0
    do i = 1, n
      do k = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
        l = l + 1
        graph%adjncy(l) = hecMAT%itemL(k)
      enddo
      do k = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
        l = l + 1
        graph%adjncy(l) = hecMAT%itemU(k)
      enddo
    enddo
  end subroutine hecmw_mf_graph_from_hecmat

  !> Build the node graph from explicit arrays (synthetic graphs in tests).
  subroutine hecmw_mf_graph_set(graph, nnode, ndof, xadj, adjncy)
    implicit none
    type(hecmwST_mf_graph), intent(out) :: graph
    integer(kind=kint), intent(in) :: nnode
    integer(kind=kint), intent(in) :: ndof(:)
    integer(kind=kint), intent(in) :: xadj(:)
    integer(kind=kint), intent(in) :: adjncy(:)
    integer(kind=kint) :: i

    graph%nnode = nnode
    allocate(graph%ndof(nnode))
    allocate(graph%dofptr(nnode+1))
    allocate(graph%xadj(nnode+1))
    allocate(graph%adjncy(xadj(nnode+1)-1))
    graph%ndof(1:nnode) = ndof(1:nnode)
    graph%dofptr(1) = 1
    do i = 1, nnode
      graph%dofptr(i+1) = graph%dofptr(i) + graph%ndof(i)
    enddo
    graph%xadj(1:nnode+1) = xadj(1:nnode+1)
    graph%adjncy(1:xadj(nnode+1)-1) = adjncy(1:xadj(nnode+1)-1)
  end subroutine hecmw_mf_graph_set

  subroutine hecmw_mf_graph_finalize(graph)
    implicit none
    type(hecmwST_mf_graph), intent(inout) :: graph

    graph%nnode = 0
    if (allocated(graph%ndof)) deallocate(graph%ndof)
    if (allocated(graph%dofptr)) deallocate(graph%dofptr)
    if (allocated(graph%xadj)) deallocate(graph%xadj)
    if (allocated(graph%adjncy)) deallocate(graph%adjncy)
  end subroutine hecmw_mf_graph_finalize

end module hecmw_mf_graph
