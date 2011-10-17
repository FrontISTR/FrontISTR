module my_hecmw_util_lag
  implicit none

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

end module my_hecmw_util_lag



