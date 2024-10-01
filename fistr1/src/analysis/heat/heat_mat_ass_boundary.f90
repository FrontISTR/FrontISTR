!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides a subroutine for all boundary conditions
!! needed in heat analysis
module m_heat_mat_ass_boundary
contains

  subroutine heat_mat_ass_boundary(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, fstrHEAT, next_time, delta_time)
    use m_fstr
    use m_heat_mat_ass_bc_CFLUX
    use m_heat_mat_ass_bc_DFLUX
    use m_heat_mat_ass_bc_FIXT
    use m_heat_mat_ass_bc_FILM
    use m_heat_mat_ass_bc_RADIATE
    implicit none
    type(fstr_heat) :: fstrHEAT
    type(hecmwST_matrix) :: hecMAT
    type(hecmwST_matrix), pointer :: hecMATmpc
    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_local_mesh), pointer :: hecMESHmpc
    real(kind=kreal) :: next_time, delta_time, beta

    beta = fstrHEAT%beta

    !> !CFLUX
    call heat_mat_ass_bc_CFLUX(hecMAT, fstrHEAT, next_time, delta_time, beta)

    !> !DFLUX
    call heat_mat_ass_bc_DFLUX(hecMESH, hecMAT, fstrHEAT, next_time, delta_time, beta)

    !> !FILM
    call heat_mat_ass_bc_FILM(hecMESH, hecMAT, fstrHEAT, next_time, delta_time, beta)

    !> !RADIATE
    call heat_mat_ass_bc_RADIATE(hecMESH, hecMAT, fstrHEAT, next_time, delta_time, beta)

    !> !BOUNDARY
    call heat_mat_ass_bc_FIXT(hecMAT, fstrHEAT, next_time, delta_time, beta)

    !> MPC
    call hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    call hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)

  end subroutine heat_mat_ass_boundary
end module m_heat_mat_ass_boundary
