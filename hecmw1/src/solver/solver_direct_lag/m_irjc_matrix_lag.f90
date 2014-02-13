!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2008/03/13                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takeshi Kitayama (Univ. of Tokyo)              !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
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
