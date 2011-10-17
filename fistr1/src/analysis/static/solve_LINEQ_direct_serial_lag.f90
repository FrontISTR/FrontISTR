module m_solve_LINEQ_direct_serial_lag

  use hecmw_util
  use fstr_matrix_con_contact
  use m_set_arrays_directsolver_contact
  use hecmw_solver_direct_serial_lag
! use hecmw_solver_    

contains

  subroutine solve_LINEQ_serial_lag_hecmw(hecMAT,fstrMAT)
  
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange
    integer (kind=4)                         :: ntdf, ilag_sta
    integer (kind=4)                         :: numNon0
    integer (kind=4)                         :: ierr, nprocs, myrank
    
    real(kind=8), allocatable               :: b(:)           !< right-hand side vector

    ntdf = hecMAT%NP*hecMAT%NDOF + fstrMAT%num_lagrange
    ilag_sta = hecMAT%NP*hecMAT%NDOF + 1
    numNon0 = hecMAT%NPU*hecMAT%NDOF**2+hecMAT%NP*hecMAT%NDOF*(ndof+1)/2 &
            + (fstrMAT%numU_lagrange)*hecMAT%NDOF+fstrMAT%num_lagrange

    allocate(b(size(hecMAT%B)))
    b = hecMAT%B

    call hecmw_solve_direct_serial_lag(ntdf, ilag_sta, numNon0, pointers, indices, values, b)

    hecMAT%X = b

    deallocate(b)

  end subroutine solve_LINEQ_serial_lag_hecmw


end module m_solve_LINEQ_direct_serial_lag
