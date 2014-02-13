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
! dummy module for serial compile
module hecmw_solver_direct_parallel

use hecmw_util

public hecmw_solve_direct_parallel

contains
subroutine hecmw_solve_direct_parallel(hecMESH, hecMAT, ii)
type (hecmwST_local_mesh) :: hecMESH
type (hecmwST_matrix    ) :: hecMAT
integer(kind=kint)        :: ii 
write(*,*)'hecmw_solver_direct_parallel_dummy'!DEBUG
return
end subroutine hecmw_solve_direct_parallel
end module hecmw_solver_direct_parallel
