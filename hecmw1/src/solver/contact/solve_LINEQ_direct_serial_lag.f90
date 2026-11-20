!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_solve_LINEQ_direct_serial_lag

  use hecmw_util
  use m_set_arrays_directsolver_contact
  use hecmw_solver_direct_serial_lag
  use hecmw_matrix_dump
  use hecmw_matrix_ass
  use hecmw_matrix_misc
  use hecmw_ebc_defer

contains

  subroutine solve_LINEQ_serial_lag_hecmw_init(hecMAT,is_sym)
    implicit none
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    logical                                  :: is_sym         !< symmetry of matrix

    call set_pointersANDindices_directsolver(hecMAT,is_sym)

  end subroutine solve_LINEQ_serial_lag_hecmw_init


  subroutine solve_LINEQ_serial_lag_hecmw(hecMESH,hecMAT,hecEBC)
    implicit none
    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (hecmwST_ebc)                       :: hecEBC         !< prescribed displacements imposed after the MPC processing
    integer (kind=4)                         :: ntdf, ilag_sta

    real(kind=8), allocatable               :: b(:)           !< right-hand side vector
    integer(kind=kint) :: mpc_method, npndof, i

    mpc_method = hecmw_mat_get_mpc_method(hecMAT)
    if (mpc_method /= 1 .and. mpc_method /= 3) then
      mpc_method = 1
      call hecmw_mat_set_mpc_method(hecMAT,mpc_method)
    endif
    if (mpc_method /= 1) then
      write(*,*) 'ERROR: MPCMETHOD other than penalty is not available for DIRECT solver', &
          ' in contact analysis without elimination'
      stop
    endif
    call hecmw_mat_ass_equation(hecMESH, hecMAT)
    call hecmw_mat_ass_equation_rhs(hecMESH, hecMAT)
    call hecmw_ebc_apply(hecMESH, hecMAT, hecEBC)

    call hecmw_mat_dump(hecMAT, hecMESH)

    call set_values_directsolver(hecMAT)

    npndof = hecMAT%NP*hecMAT%NDOF
    ntdf = npndof + hecMAT%lagmat%num_lagrange
    ilag_sta = npndof + 1

    allocate(b(ntdf))
    b(1:npndof) = hecMAT%B(1:npndof)
    do i = 1, hecMAT%lagmat%num_lagrange
      b(npndof+i) = hecMAT%lagmat%B_lagrange(i)
    enddo

    call hecmw_solve_direct_serial_lag(ntdf, ilag_sta, pointers, indices, values, b)

    hecMAT%X(1:npndof) = b(1:npndof)
    do i = 1, hecMAT%lagmat%num_lagrange
      hecMAT%lagmat%X_lagrange(i) = b(npndof+i)
    enddo

    deallocate(b)

    call hecmw_mat_dump_solution(hecMAT)

  end subroutine solve_LINEQ_serial_lag_hecmw


end module m_solve_LINEQ_direct_serial_lag
