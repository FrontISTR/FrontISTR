!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_solve_LINEQ_direct_serial_lag

  use hecmw_util
  use m_set_arrays_directsolver_contact
  use hecmw_solver_direct_serial_lag
  ! use hecmw_solver_

contains

  subroutine solve_LINEQ_serial_lag_hecmw_init(hecMAT,hecLagMAT,is_sym)
    implicit none
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (hecmwST_matrix_lagrange)           :: hecLagMAT        !< type hecmwST_matrix_lagrange
    logical                                  :: is_sym         !< symmetry of matrix

    call set_pointersANDindices_directsolver(hecMAT,hecLagMAT,is_sym)

  end subroutine solve_LINEQ_serial_lag_hecmw_init


  subroutine solve_LINEQ_serial_lag_hecmw(hecMESH,hecMAT,hecLagMAT)
    implicit none
    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (hecmwST_matrix_lagrange)           :: hecLagMAT        !< type hecmwST_matrix_lagrange
    integer (kind=4)                         :: ntdf, ilag_sta
    integer (kind=4)                         :: numNon0
    integer (kind=4)                         :: ierr, nprocs, myrank

    real(kind=8), allocatable               :: b(:)           !< right-hand side vector

    call hecmw_mat_dump(hecMAT, hecMESH)

    call set_values_directsolver(hecMAT,hecLagMAT)

    ntdf = hecMAT%NP*hecMAT%NDOF + hecLagMAT%num_lagrange
    ilag_sta = hecMAT%NP*hecMAT%NDOF + 1
    numNon0 = hecMAT%NPU*hecMAT%NDOF**2+hecMAT%NP*hecMAT%NDOF*(ntdf+1)/2 &
      + (hecLagMAT%numU_lagrange)*hecMAT%NDOF+hecLagMAT%num_lagrange

    allocate(b(size(hecMAT%B)))
    b = hecMAT%B

    call hecmw_solve_direct_serial_lag(ntdf, ilag_sta, numNon0, pointers, indices, values, b)

    hecMAT%X = b

    deallocate(b)

    call hecmw_mat_dump_solution(hecMAT)

  end subroutine solve_LINEQ_serial_lag_hecmw


end module m_solve_LINEQ_direct_serial_lag
