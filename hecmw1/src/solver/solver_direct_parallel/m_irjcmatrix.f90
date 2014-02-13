!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2007/11/21                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takeshi Kitayama  (Univ. of Tokyo)             !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

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
