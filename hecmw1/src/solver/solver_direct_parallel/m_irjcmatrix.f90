!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module m_irjcmatrix
! irow, jcol type sparse matrix

type irjc_matrix
integer :: neqns ! dimension of matrix. each element have (ndeg,ndeg) small matrix
integer :: nttbr ! number of nonzero elements
integer :: ndeg  ! degree of freedom of each element
integer, pointer :: irow(:), jcol(:) ! location of i'th element
real(8), pointer :: val(:,:) ! concrete value of nonzero elements. small matrix is expanded as ndeg*ndeg.
end type irjc_matrix
end module m_irjcmatrix
