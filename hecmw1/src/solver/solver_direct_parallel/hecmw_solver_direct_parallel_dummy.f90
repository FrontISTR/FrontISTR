!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
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
