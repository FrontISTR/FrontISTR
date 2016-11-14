!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_irjc_matrix_lag

use my_hecmw_util_lag

! Handle irow, jcol type sparse matrix information.
! Eqch matrix element have ndeg*ndeg small matrix.
! The number of non-zero elements are stored in nttbr.

type irjc_square_matrix
  ! number of equations
  integer(kind=kint) :: neqns

  ! number of non-zero elements
  integer(kind=kint) :: nttbr

  ! degree of freedom of each element
  integer(kind=kint) :: ndeg
  integer, pointer :: irow(:), jcol(:) ! irow(nttbr_t), jcol(nttbr_t) location of i'th element
  real(kind=kreal), pointer :: val(:,:) ! val(ndeg*ndeg, nttbr_t) concrete value of nonzero elements. small matrix is expanded as ndeg*ndeg.
end type irjc_square_matrix

type irjc_mn_matrix !m*n matrix
  ! number of rows.
  integer(kind=kint) :: nrows

  ! number of columns.
  integer(kind=kint) :: ncols

  ! number of non-zero elements
  integer(kind=kint) :: nttbr

  ! degree of freedom of each element
  integer(kind=kint) :: ndeg
  integer(kind=kint), pointer :: irow(:), jcol(:) ! irow(nttbr_t), jcol(nttbr_t) location of i'th element
  real(kind=kreal), pointer :: val(:,:) ! val(ndeg*ndeg, nttbr_t) concrete value of nonzero elements. small matrix is expanded as ndeg*ndeg.
end type irjc_mn_matrix
end module m_irjc_matrix_lag
