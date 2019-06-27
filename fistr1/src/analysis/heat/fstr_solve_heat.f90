!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function to control heat analysis
module m_fstr_solve_heat
contains

  subroutine fstr_solve_heat(hecMESH, hecMAT, fstrRESULT, fstrPARAM, fstrHEAT)
    use m_fstr
    use m_heat_init
    use m_heat_solve_TRAN
    implicit none
    integer(kind=kint) :: i, in, ISTEP
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(hecmwST_result_data) :: fstrRESULT
    type(fstr_param)          :: fstrPARAM
    type(fstr_heat)           :: fstrHEAT

    integer(kind=kint) :: total_step

    call heat_init(hecMESH, fstrHEAT)
    call heat_init_log(hecMESH)

    total_step = 1
    do ISTEP = 1, fstrHEAT%STEPtot
      fstrHEAT%is_steady = 0
      if(fstrHEAT%STEP_DLTIME(ISTEP) <= 0.0d0) fstrHEAT%is_steady = 1

      if(hecMESH%my_rank == 0)then
        write(IMSG,"(a,i8,a,i8)")"* Current step / Total step: ", ISTEP, "/", fstrHEAT%STEPtot
        write(IMSG,"(a,i8)")"* max iteration at each step: ", fstrPARAM%ITMAX(ISTEP)
        if(fstrHEAT%is_steady == 1)then
          write(IMSG,"(a)")"* Steady state analysis"
        else
          write(IMSG,"(a)")"* Transient anslysis"
        endif
      endif

      call heat_solve_TRAN(hecMESH, hecMAT, fstrRESULT, fstrPARAM, fstrHEAT, ISTEP, total_step)
    enddo

    call heat_finalize(fstrHEAT)
  end subroutine fstr_solve_heat
end module m_fstr_solve_heat
