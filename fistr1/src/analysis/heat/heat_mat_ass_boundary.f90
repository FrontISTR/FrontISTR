!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine for all boundary conditions
!! needed in heat anaylsis
module m_heat_mat_ass_boundary
contains
  !C***
  !C*** MAT_ASS_BOUNDARY
  !C***
  subroutine heat_mat_ass_boundary ( hecMESH,hecMAT,hecMATmpc,fstrHEAT,ATIME, BTIME, DTIME )

    use m_fstr
    use m_heat_mat_ass_bc_CFLUX
    use m_heat_mat_ass_bc_DFLUX
    use m_heat_mat_ass_bc_FIXT
    use m_heat_mat_ass_bc_FILM
    use m_heat_mat_ass_bc_RADIATE

    implicit none
    real(kind=kreal) :: ATIME, BTIME, CTIME, DTIME
    type(fstr_heat)          :: fstrHEAT
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_matrix), pointer :: hecMATmpc
    type(hecmwST_local_mesh) :: hecMESH

    CTIME = ATIME + BTIME

    !C
    !C +---------+
    !C | !CFLUX  |
    !C +---------+
    !C===
    call heat_mat_ass_bc_CFLUX ( hecMAT, fstrHEAT, CTIME )
    !C
    !C +---------+
    !C | !DFLUX  |
    !C +---------+
    !C===
    call heat_mat_ass_bc_DFLUX ( hecMESH, hecMAT, fstrHEAT, CTIME, DTIME )
    !C
    !C +--------+
    !C | !FILM  |
    !C +--------+
    !C===
    call heat_mat_ass_bc_FILM ( hecMESH, hecMAT, fstrHEAT, CTIME )
    !C
    !C +-----------+
    !C | !RADIATE  |
    !C +-----------+
    !C===
    call heat_mat_ass_bc_RADIATE ( hecMESH, hecMAT, fstrHEAT, CTIME )
    !C
    !C +------+
    !C | MPC  |
    !C +------+
    !C===
    call hecmw_mpc_mat_ass( hecMESH, hecMAT, hecMATmpc )
    call hecmw_mpc_trans_rhs( hecMESH, hecMAT, hecMATmpc )
    !C
    !C +------------+
    !C | !BOUNDARY  |
    !C +------------+
    !C===
    call heat_mat_ass_bc_FIXT ( hecMATmpc, fstrHEAT, CTIME )
    return

  end subroutine heat_mat_ass_boundary
end module m_heat_mat_ass_boundary
