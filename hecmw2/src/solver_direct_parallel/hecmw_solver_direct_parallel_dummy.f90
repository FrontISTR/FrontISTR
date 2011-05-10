!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.1                                                !
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
! dummy module for serial compile
module hecmw_solver_direct_parallel

integer(kind=4),parameter:: kint  = 4
integer(kind=4),parameter:: kreal = 8
type irjc_square_matrix
  integer(kind=kint) :: neqns
  integer(kind=kint) :: nttbr
  integer(kind=kint) :: ndeg 
  integer, pointer :: irow(:), jcol(:)
  real(kind=kreal), pointer :: val(:,:)
end type irjc_square_matrix

public hecmw_solve_direct_parallel

contains
subroutine hecmw_solve_direct_parallel(a0,b)
type (irjc_square_matrix), intent(inout) :: a0
real(kind=kreal), intent(inout) :: b(:)
write(*,*)'hecmw_solver_direct_parallel_dummy'
return
end subroutine hecmw_solve_direct_parallel
end module hecmw_solver_direct_parallel
