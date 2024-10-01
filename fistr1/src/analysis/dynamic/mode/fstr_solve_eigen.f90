!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function to control eigen analysis
module m_fstr_solve_eigen
contains

  !> solve eigenvalue probrem
  subroutine fstr_solve_eigen( hecMESH, hecMAT, fstrEIG, fstrSOLID, &
      & fstrRESULT, fstrPARAM, hecLagMAT)
    use hecmw_util
    use m_fstr
    use m_fstr_StiffMatrix
    use m_fstr_AddBC
    use m_fstr_EIG_setMASS
    use m_fstr_EIG_lanczos
    use m_fstr_EIG_output
    use m_static_lib
    use m_hecmw2fstr_mesh_conv
    use m_fstr_spring

    implicit none

    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(fstr_solid)          :: fstrSOLID
    type(hecmwST_result_data) :: fstrRESULT
    type(fstr_param)          :: fstrPARAM
    type(fstr_eigen)          :: fstrEIG
    type(hecmwST_matrix_lagrange) :: hecLagMAT

    type(hecmwST_local_mesh), pointer :: hecMESHmpc
    type(hecmwST_matrix), pointer :: hecMATmpc
    real(kind=kreal) :: t1, t2

    t1 = hecmw_Wtime()

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    fstrSOLID%dunode = 0.0d0
    call fstr_StiffMatrix(hecMESH, hecMAT, fstrSOLID, 0.0d0, 0.0d0)
    call fstr_AddSPRING(1, hecMESH, hecMAT, fstrSOLID, fstrPARAM)

    call fstr_AddBC(1,  hecMESH, hecMAT, fstrSOLID, fstrPARAM, hecLagMAT, 2)
    call hecmw_mpc_mat_ass(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)
    call hecmw_mpc_trans_rhs(hecMESH, hecMAT, hecMATmpc)

    call setMASS(fstrSOLID, hecMESH, hecMAT, fstrEIG)
    call hecmw_mpc_trans_mass(hecMESH, hecMAT, fstrEIG%mass)

    call fstr_solve_lanczos(hecMESHmpc, hecMATmpc, fstrSOLID, fstrEIG)

    call hecmw_mpc_tback_eigvec(hecMESH, hecMAT, fstrEIG%iter, fstrEIG%eigvec)

    call fstr_eigen_output(hecMESH, hecMAT, fstrEIG)

    call fstr_eigen_make_result(hecMESH, hecMAT, fstrEIG, fstrRESULT)

    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMESHmpc, hecMATmpc)

    t2 = hecmw_Wtime()

    if(myrank == 0)then
      write(IMSG,'("### FSTR_SOLVE_EIGEN FINISHED!")')
      write(*,'("### FSTR_SOLVE_EIGEN FINISHED!")')
    endif

  end subroutine fstr_solve_eigen
end module m_fstr_solve_eigen
