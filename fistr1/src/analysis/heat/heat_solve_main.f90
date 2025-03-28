!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function for stationary heat analysis
module m_heat_solve_main
contains

  subroutine heat_solve_main(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, &
     & fstrSOLID, fstrPARAM, fstrHEAT, ISTEP, iterALL, next_time, delta_time)
    use m_fstr
    use m_heat_mat_ass_conductivity
    use m_heat_mat_ass_capacity
    use m_heat_mat_ass_boundary
    use m_solve_lineq
    use m_fstr_elemact
    implicit none
    integer(kind=kint) :: i, iterALL, ISTEP, bup_n_dof
    real(kind=kreal)   :: delta_time, next_time
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(fstr_solid)          :: fstrSOLID
    type(fstr_heat)           :: fstrHEAT
    type(fstr_param)          :: fstrPARAM
    type(hecmwST_local_mesh), pointer :: hecMESHmpc
    type(hecmwST_matrix), pointer :: hecMATmpc
    logical :: is_congerged

    if( fstrHEAT%elemact%ELEMACT_egrp_tot>0 ) &
      &  call fstr_update_elemact_heat( hecMESH, fstrHEAT%elemact, next_time, fstrSOLID%elements )

    iterALL = 0
    do
      iterALL = iterALL + 1
      hecMAT%X = 0.0d0

      call heat_mat_ass_conductivity(hecMESH, hecMAT, fstrSOLID, fstrHEAT, fstrHEAT%beta)
      if(fstrHEAT%is_steady == 0) call heat_mat_ass_capacity(hecMESH, hecMAT, fstrSOLID, fstrHEAT, delta_time)
      call heat_mat_ass_boundary(hecMESH, hecMAT, hecMESHmpc, hecMATmpc, &
       & fstrSOLID, fstrHEAT, next_time, delta_time)

      hecMATmpc%Iarray(97) = 1 !Need numerical factorization
      bup_n_dof = hecMESH%n_dof
      hecMESH%n_dof = 1
      call solve_LINEQ(hecMESHmpc, hecMATmpc)
      hecMESH%n_dof = bup_n_dof
      call hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)

      do i = 1, hecMESH%n_node
        fstrHEAT%TEMPC(i) = fstrHEAT%TEMP(i)
        fstrHEAT%TEMP (i) = hecMAT%X(i)
      enddo

      if( fstrHEAT%elemact%ELEMACT_egrp_tot>0 ) &
        &  call fstr_updatedof_elemact( 1, hecMESH, fstrHEAT%elemact, fstrSOLID%elements, fstrHEAT%TEMP0, fstrHEAT%TEMP )

      call heat_check_convergence(hecMESH, fstrHEAT, fstrPARAM, ISTEP, iterALL, is_congerged)

      if(is_congerged .or. fstrHEAT%is_iter_max_limit .or. fstrHEAT%beta == 0.0d0) exit
    enddo
  end subroutine heat_solve_main

  subroutine heat_check_convergence(hecMESH, fstrHEAT, fstrPARAM, ISTEP, iterALL, is_congerged)
    use m_fstr
    implicit none
    integer(kind=kint) :: i, iterALL, iter_max, ISTEP
    real(kind=kreal)   :: val, iter_eps
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_heat)           :: fstrHEAT
    type(fstr_param)          :: fstrPARAM
    logical :: is_congerged

    is_congerged = .false.
    fstrHEAT%is_iter_max_limit = .false.
    iter_max = fstrPARAM%itmax(ISTEP)
    iter_eps = fstrPARAM%eps(ISTEP)

    val = 0.0d0
    do i = 1, hecMESH%nn_internal
      val = val + (fstrHEAT%TEMP(i) - fstrHEAT%TEMPC(i))**2
    enddo
    call hecmw_allREDUCE_R1(hecMESH, val, hecmw_sum)
    val = dsqrt(val)

    if(hecMESH%my_rank == 0)then
      write(*,'(i8,1pe16.6)') iterALL, val
    endif

    if(val < iter_eps)then
      if(hecMESH%my_rank == 0) write(IMSG,*) ' !!! CONVERGENCE ACHIEVED '
      is_congerged = .true.
    endif

    if(iter_max <= iterALL)then
      if(hecMESH%my_rank == 0) write(*,*) ' !!! ITERATION COUNT OVER : MAX = ', iter_max
      fstrHEAT%is_iter_max_limit = .true.
    endif
  end subroutine heat_check_convergence
end module m_heat_solve_main
