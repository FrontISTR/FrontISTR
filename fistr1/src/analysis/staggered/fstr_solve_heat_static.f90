!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function to control heat analysis
module m_fstr_solve_heat_static

contains

  subroutine fstr_solve_heat_static(hecMESH, hecMAT, fstrSOLID, fstrHEAT, fstrMAT, fstrPARAM, conMAT)
    use m_fstr
    use m_heat_init
    use m_fstr_solve_NLGEOM
    type (hecmwST_local_mesh) :: hecMESH
    type (hecmwST_matrix)     :: hecMAT
    type (fstr_param)         :: fstrPARAM
    type (fstr_solid)         :: fstrSOLID
    type (fstr_heat)          :: fstrHEAT
    type (fstrST_matrix_contact_lagrange) :: fstrMAT
    type (hecmwST_matrix), optional :: conMAT

    call heat_init(hecMESH, fstrHEAT)

    call FSTR_SOLVE_NLGEOM(hecMESH, hecMAT, fstrSOLID, fstrMAT, fstrPARAM, &
      & conMAT=conMAT, fstrHEAT=fstrHEAT)

  end subroutine fstr_solve_heat_static

end module m_fstr_solve_heat_static
