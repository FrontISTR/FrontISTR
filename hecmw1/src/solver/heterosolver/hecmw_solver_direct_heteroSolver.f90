!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface for HeteroSolver
module hecmw_solver_direct_HeteroSolver
  use hecmw_util
  use m_sparse_matrix
  use m_sparse_matrix_hec
  use m_hecmw_HeteroSolver_wrapper
  use hecmw_matrix_ass
  use hecmw_matrix_dump

  private
  public :: hecmw_solve_direct_HeteroSolver

  logical, save :: INITIALIZED = .false.
  type (sparse_matrix), save :: spMAT

contains

  subroutine hecmw_solve_direct_HeteroSolver(hecMESH,hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: myrank
    real(kind=kreal) :: t1,t2,t3

#ifdef HECMW_WITH_HETEROSOLVER
    write(0,*) 'hecmw_solve_direct_HeteroSolver called'
    call hecmw_mat_dump(hecMAT, hecMESH)

    t1=hecmw_wtime()
    myrank=hecmw_comm_get_rank()

    if (.not. INITIALIZED) then
      spmat_type = SPARSE_MATRIX_TYPE_CSR
      if (hecMAT%symmetric) then
        spmat_symtype = SPARSE_MATRIX_SYMTYPE_SPD
      else
        spmat_symtype = SPARSE_MATRIX_SYMTYPE_ASYM
      end if
      call sparse_matrix_set_type(spMAT, spmat_type, spmat_symtype)
      call sparse_matrix_hec_init_prof(spMAT, hecMAT, hecMESH)
      call sparse_matrix_hec_set_vals(spMAT, hecMAT)
      !call sparse_matrix_dump(spMAT)
      call hecmw_HeteroSolver_wrapper_init(spMAT)
      INITIALIZED = .true.
      hecMAT%Iarray(97) = 1
      hecMAT%Iarray(98) = 1
    endif

    !* Flag to activate symbolic factorization: 1(yes) 0(no)  hecMESH%Iarray(98)
    !* Flag to activate numeric  factorization: 1(yes) 0(no)  hecMESH%Iarray(97)
    if (hecMAT%Iarray(97) .eq. 1) then
      call hecmw_HeteroSolver_wrapper_preprocess(spMAT)
      if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
           write(*,*) ' [HeteroSolver]: Analysis and symbolic factorization completed.'
      hecMAT%Iarray(97) = 0
      hecMAT%Iarray(98) = 1
    endif
    if (hecMAT%Iarray(98) .eq. 1) then
      call hecmw_HeteroSolver_wrapper_factorize(spMAT)
      if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
           write(*,*) ' [HeteroSolver]: Numerical factorization completed.'
      hecMAT%Iarray(98) = 0
    endif
    t2=hecmw_wtime()

    ! SOLUTION
    call sparse_matrix_hec_set_rhs(spMAT, hecMAT)
    call hecmw_HeteroSolver_wrapper_solve(spMAT)
    call sparse_matrix_hec_get_rhs(spMAT, hecMAT)
    if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
         write(*,*) ' [HeteroSolver]: Solution completed.'

    t3=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) then
      write(*,*) 'setup time : ',t2-t1
      write(*,*) 'solve time : ',t3-t2
    endif

    call hecmw_mat_dump_solution(hecMAT)
#else
    stop "HeteroSolver not available"
#endif
  end subroutine hecmw_solve_direct_HeteroSolver
end module hecmw_solver_direct_HeteroSolver
