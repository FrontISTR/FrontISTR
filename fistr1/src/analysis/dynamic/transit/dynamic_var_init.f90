!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to initialize variables
!> when initial velocity or acceleration boundary conditions are given.
!> attention : just for rigid motion in the initial state.

module m_dynamic_init_variables

  use m_fstr
  use m_dynamic_mat_ass_load
  use m_solve_LINEQ
  use hecmw_solver_las, only: hecmw_matvec

contains

  subroutine dynamic_init_varibles( hecMESH, hecMAT, fstrSOLID, fstrEIG, fstrDYNAMIC, fstrPARAM, mass_matrix )

    implicit none

    type(hecmwST_local_mesh) :: hecMESH
    type(hecmwST_matrix)     :: hecMAT
    type(fstr_eigen)         :: fstrEIG
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC
    type(fstr_param)         :: fstrPARAM
    type(hecmwST_matrix), pointer, optional :: mass_matrix

    integer(kind=kint) :: j, n_internal_dof
    real(kind=kreal), allocatable :: mass_vector(:)
    logical :: use_consistent_mass

    call dynamic_mat_ass_load (1, 0.d0, hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM)

    use_consistent_mass = .false.
    if( present(mass_matrix) ) use_consistent_mass = associated(mass_matrix)
    n_internal_dof = hecMAT%N*hecMAT%NDOF

    if( use_consistent_mass ) then
      allocate(mass_vector(mass_matrix%NP*mass_matrix%NDOF))
      call hecmw_mat_copy_val(mass_matrix, hecMAT)
      hecMAT%X = 0.0d0

      if( fstrSOLID%VELOCITY_type == kbcInitial ) then
        call hecmw_matvec(hecMESH, mass_matrix, fstrDYNAMIC%VEL(:,1), mass_vector)
        hecMAT%B(1:n_internal_dof) = hecMAT%B(1:n_internal_dof) - &
          fstrDYNAMIC%ray_m*mass_vector(1:n_internal_dof)
        call solve_LINEQ(hecMESH, hecMAT)
        fstrDYNAMIC%ACC(1:n_internal_dof,1) = hecMAT%X(1:n_internal_dof)
        call hecmw_update_R(hecMESH, fstrDYNAMIC%ACC(:,1), hecMAT%NP, hecMAT%NDOF)
      elseif( fstrSOLID%ACCELERATION_type == kbcInitial ) then
        call hecmw_matvec(hecMESH, mass_matrix, fstrDYNAMIC%ACC(:,1), mass_vector)
        hecMAT%B(1:n_internal_dof) = hecMAT%B(1:n_internal_dof) - mass_vector(1:n_internal_dof)
        hecMAT%D = fstrDYNAMIC%ray_m*hecMAT%D
        hecMAT%AL = fstrDYNAMIC%ray_m*hecMAT%AL
        hecMAT%AU = fstrDYNAMIC%ray_m*hecMAT%AU
        call solve_LINEQ(hecMESH, hecMAT)
        fstrDYNAMIC%VEL(1:n_internal_dof,1) = hecMAT%X(1:n_internal_dof)
        call hecmw_update_R(hecMESH, fstrDYNAMIC%VEL(:,1), hecMAT%NP, hecMAT%NDOF)
      endif

      deallocate(mass_vector)
      return
    endif

    if( fstrSOLID%VELOCITY_type == kbcInitial ) then
      do j = 1, hecMESH%n_node*hecMESH%n_dof
        fstrDYNAMIC%ACC(j,1)=(hecMAT%B(j)-fstrDYNAMIC%ray_m*fstrEIG%mass(j)*fstrDYNAMIC%VEL(j,1))/&
          fstrEIG%mass(j)
      enddo
    elseif( fstrSOLID%ACCELERATION_type == kbcInitial ) then
      do j = 1, hecMESH%n_node*hecMESH%n_dof
        fstrDYNAMIC%VEL(j,1)=(hecMAT%B(j)-fstrEIG%mass(j)*fstrDYNAMIC%ACC(j,1))/&
          (fstrDYNAMIC%ray_m*fstrEIG%mass(j))
      enddo
    endif

  end subroutine dynamic_init_varibles

end module m_dynamic_init_variables
