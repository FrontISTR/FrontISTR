module my_hecmw_util
  implicit none
  include 'mpif.h'

  integer(kind=4), parameter:: kint  = 4
  integer(kind=4), parameter:: kreal = 8

  type hecmwST_matrix
    integer(kind=kint) ::  NP, NPL, NDOF
    real(kind=kreal), pointer :: D(:), B(:), X(:)
    real(kind=kreal), pointer :: AL(:)
    integer(kind=kint), pointer :: indexL(:)
    integer(kind=kint), pointer ::  itemL(:)
    integer(kind=kint ), dimension(100) :: Iarray
    real   (kind=kreal), dimension(100) :: Rarray
  end type hecmwST_matrix

contains

  subroutine hecmw_abort(comm)
    integer(kind=kint) :: comm
    call MPI_ABORT(mpi_comm_world)
  end subroutine hecmw_abort

  function hecmw_comm_get_comm() result(comm)
    integer(kind=kint) :: comm
    comm = mpi_comm_world
  end function hecmw_comm_get_comm

end module my_hecmw_util
