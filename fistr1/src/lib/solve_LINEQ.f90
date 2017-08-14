!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This program is a HECMW interface to a set of linear iterative and direct
!! solvers. The interface may be called from within a HECMW application, with
!! an appropriate choice of TYPE (iterative, direct), and METHOD (depending
!! on the iterative solver used).
module m_solve_LINEQ
  implicit none

  private
  public :: solve_LINEQ

  contains

  SUBROUTINE solve_LINEQ(hecMESH,hecMAT,imsg)
    USE hecmw
    USE hecmw_solver

    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix    ) :: hecMAT
    INTEGER(kind=kint) imsg
    
    call hecmw_solve(hecMESH,hecMAT,imsg)

    RETURN

  end subroutine solve_LINEQ

end module m_solve_LINEQ
