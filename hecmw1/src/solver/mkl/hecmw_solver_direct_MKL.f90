!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface for Pardiso

module hecmw_solver_direct_MKL
  use hecmw_util
  use m_sparse_matrix
  use m_sparse_matrix_hec
  use m_hecmw_MKL_wrapper
  use hecmw_matrix_ass
  use hecmw_matrix_dump

  private
  public :: hecmw_solve_direct_MKL

  logical, save :: INITIALIZED = .false.
  type (sparse_matrix), save :: spMAT

contains

  subroutine hecmw_solve_direct_MKL(hecMESH,hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: phase_start
    integer(kind=kint) :: istat,myrank
    real(kind=kreal) :: t1,t2

    t1=hecmw_wtime()

    call hecmw_mat_dump(hecMAT, hecMESH)

    myrank=hecmw_comm_get_rank()

    if (INITIALIZED .and. hecMAT%Iarray(98) .eq. 1) then
      call sparse_matrix_finalize(spMAT)
      INITIALIZED = .false.
    endif

    if (.not. INITIALIZED) then
      spmat_type = SPARSE_MATRIX_TYPE_CSR
      if (hecMAT%symmetric) then
        spmat_symtype = SPARSE_MATRIX_SYMTYPE_SYM
      else
        spmat_symtype = SPARSE_MATRIX_SYMTYPE_ASYM
      end if
      call sparse_matrix_set_type(spMAT, spmat_type, spmat_symtype)

      INITIALIZED = .true.
      hecMAT%Iarray(98) = 1
    endif

    !* Flag to activate symbolic factorization: 1(yes) 0(no)  hecMESH%Iarray(98)
    !* Flag to activate numeric  factorization: 1(yes) 0(no)  hecMESH%Iarray(97)

    if (hecMAT%Iarray(98) .eq. 1) then
      ! ANALYSIS and FACTORIZATION
      call sparse_matrix_hec_init_prof(spMAT, hecMAT, hecMESH)
      call sparse_matrix_hec_set_vals(spMAT, hecMAT)
      !call sparse_matrix_dump(spMAT)
      phase_start = 1
      hecMAT%Iarray(98) = 0
      hecMAT%Iarray(97) = 0
    endif
    if (hecMAT%Iarray(97) .eq. 1) then
      ! FACTORIZATION
      call sparse_matrix_hec_set_vals(spMAT, hecMAT)
      !call sparse_matrix_dump(spMAT)
      phase_start = 2
      hecMAT%Iarray(97) = 0
    endif
    call sparse_matrix_hec_set_rhs(spMAT, hecMAT)

    t2=hecmw_wtime()
    if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
       write(*,'(A,f10.3)') ' [Pardiso]: Setup completed.          time(sec)=',t2-t1

    ! SOLVE
    call hecmw_mkl_wrapper(spMAT, phase_start, hecMAT%X, istat)

    call hecmw_mat_dump_solution(hecMAT)

  end subroutine hecmw_solve_direct_MKL

end module hecmw_solver_direct_MKL
