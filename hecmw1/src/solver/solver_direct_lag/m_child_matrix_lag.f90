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
module m_child_matrix_lag

use my_hecmw_util_lag
use m_irjc_matrix_lag

type child_matrix

  integer(kind=kint) :: ndeg    ! same for A and C
  integer(kind=kint) :: ista_c  ! begining index of row of C
  integer(kind=kint) :: neqns_t ! total number of equations.

  ! A region
  type(irjc_square_matrix) :: a

  ! C region
  type(irjc_mn_matrix) :: c

end type child_matrix
end module m_child_matrix_lag
