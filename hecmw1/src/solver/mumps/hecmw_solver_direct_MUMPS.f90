!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface for MUMPS
module hecmw_solver_direct_MUMPS
  use hecmw_util
  use m_sparse_matrix
  use m_sparse_matrix_hec
  use m_hecmw_MUMPS_wrapper
  use hecmw_matrix_ass
  use hecmw_matrix_dump

  private
  public :: hecmw_solve_direct_MUMPS

  logical, save :: INITIALIZED = .false.
  type (sparse_matrix), save :: spMAT

contains

  subroutine hecmw_solve_direct_MUMPS(hecMESH,hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: mumps_job
    integer(kind=kint) :: istat,myrank
    real(kind=kreal) :: t1,t2,t3

    call hecmw_mat_dump(hecMAT, hecMESH)

    t1=hecmw_wtime()
    myrank=hecmw_comm_get_rank()

    if (INITIALIZED .and. hecMAT%Iarray(98) .eq. 1) then
      mumps_job=-2
      call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MUMPS returned with error', istat
        stop
      endif
      call sparse_matrix_finalize(spMAT)
      INITIALIZED = .false.
    endif

    if (.not. INITIALIZED) then
      spmat_type = SPARSE_MATRIX_TYPE_COO
      if (hecMAT%symmetric) then
        spmat_symtype = SPARSE_MATRIX_SYMTYPE_SPD
      else
        spmat_symtype = SPARSE_MATRIX_SYMTYPE_ASYM
      end if
      call sparse_matrix_set_type(spMAT, spmat_type, spmat_symtype)
      mumps_job = -1
      call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MUMPS returned with error', istat
        stop
      endif
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
      mumps_job=4
      call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MUMPS returned with error', istat
        stop
      endif
      if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
           write(*,*) ' [MUMPS]: Analysis and Factorization completed.'
      hecMAT%Iarray(98) = 0
      hecMAT%Iarray(97) = 0
    endif
    if (hecMAT%Iarray(97) .eq. 1) then
      ! FACTORIZATION
      call sparse_matrix_hec_set_vals(spMAT, hecMAT)
      !call sparse_matrix_dump(spMAT)
      mumps_job=2
      call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MUMPS returned with error', istat
        stop
      endif
      if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
           write(*,*) ' [MUMPS]: Factorization completed.'
      hecMAT%Iarray(97) = 0
    endif

    t2=hecmw_wtime()

    ! SOLUTION
    call sparse_matrix_hec_set_rhs(spMAT, hecMAT)
    mumps_job=3
    call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
    if (istat < 0) then
      write(*,*) 'ERROR: MUMPS returned with error', istat
      stop
    endif
    call sparse_matrix_hec_get_rhs(spMAT, hecMAT)
    if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
         write(*,*) ' [MUMPS]: Solution completed.'

    t3=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) then
      write(*,*) 'setup time : ',t2-t1
      write(*,*) 'solve time : ',t3-t2
    endif

    call hecmw_mat_dump_solution(hecMAT)

    !call sparse_matrix_finalize(spMAT)
  end subroutine hecmw_solve_direct_MUMPS

end module hecmw_solver_direct_MUMPS
