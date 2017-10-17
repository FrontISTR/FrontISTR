!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module m_solve_LINEQ
  implicit none

  private

  public :: solve_LINEQ

contains

  subroutine solve_LINEQ(hecMESH, hecMAT)
    use hecmw
    use hecmw_solver

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT

    call hecmw_solve(hecMESH, hecMAT)

  end subroutine solve_LINEQ
end module m_solve_LINEQ
