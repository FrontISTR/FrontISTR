!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function for stationary heat analysis
module m_heat_solve_SS
contains

  subroutine heat_solve_SS ( hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrHEAT,ISTEP )
    use m_fstr
    use m_heat_mat_ass_conductivity
    use m_heat_mat_ass_boundary
    use m_heat_init
    use m_heat_solve_main
    use m_solve_lineq
    use m_heat_io
    implicit none
    integer(kind=kint) :: ISTEP, iterALL, ITM, i, INCR, LMAX, LMIN, inod, ii, bup_n_dof
    real(kind=kreal)   :: CTIME, STIME, val, CHK, TMAX, TMIN, temp
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(fstr_solid)          :: fstrSOLID
    type(hecmwST_result_data) :: fstrRESULT
    type(fstr_param)          :: fstrPARAM
    type(fstr_heat)           :: fstrHEAT
    type(hecmwST_matrix), pointer   :: hecMATmpc
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: label
    character(len=HECMW_NAME_LEN)   :: nameID

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMATmpc)

    STIME = 0.0d0
    ETIME = 0.0d0
    hecMAT%X = 0.0d0
    hecMAT%NDOF = 1
    hecMAT%Iarray(98) = 1 !Assmebly complete
    fstrHEAT%beta = 1.0d0

    call heat_solve_main(hecMESH, hecMAT, hecMATmpc, fstrSOLID, fstrPARAM, fstrHEAT, ISTEP, 0.0d0, 0.0d0)

    call hecmw_mpc_mat_finalize( hecMESH, hecMAT, hecMATmpc )
    call heat_output_log(hecMESH, fstrPARAM, fstrHEAT, ISTEP, 1.0d0)
    call heat_output_result(hecMESH, fstrHEAT, 1, 1)
    call heat_output_visual(hecMESH, fstrRESULT, fstrHEAT, 1, 1)

  end subroutine heat_solve_SS
end module m_heat_solve_SS