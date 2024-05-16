!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface of MUMPS for
!! contact problems using Lagrange multiplier.
module m_solve_LINEQ_MUMPS_contact
  use hecmw_util
  use m_sparse_matrix
  use m_sparse_matrix_contact
  use m_hecmw_MUMPS_wrapper
  use hecmw_matrix_dump
  use hecmw_matrix_ass
  use hecmw_matrix_misc

  private
  public :: solve_LINEQ_MUMPS_contact_init
  public :: solve_LINEQ_MUMPS_contact

  logical, save :: INITIALIZED = .false.
  logical, save :: NEED_ANALYSIS = .true.
  type (sparse_matrix), save :: spMAT

contains

  subroutine solve_LINEQ_MUMPS_contact_init(hecMESH,hecMAT,hecLagMAT,is_sym)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(in) :: hecLagMAT !< type hecmwST_matrix_lagrange
    logical, intent(in) :: is_sym
    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: mumps_job
    integer(kind=kint) :: istat

    if (INITIALIZED) then
      mumps_job=-2
      call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MUMPS returned with error', istat
        stop
      endif
      call sparse_matrix_finalize(spMAT)
      INITIALIZED = .false.
    endif

    if (is_sym) then
      !spmat_symtype = SPARSE_MATRIX_SYMTYPE_SPD
      spmat_symtype = SPARSE_MATRIX_SYMTYPE_SYM
    else
      spmat_symtype = SPARSE_MATRIX_SYMTYPE_ASYM
    endif
    spmat_type = SPARSE_MATRIX_TYPE_COO
    call sparse_matrix_set_type(spMAT, spmat_type, spmat_symtype)
    mumps_job=-1
    call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
    if (istat < 0) then
      write(*,*) 'ERROR: MUMPS returned with error', istat
      stop
    endif

    NEED_ANALYSIS = .true.
    INITIALIZED = .true.
  end subroutine solve_LINEQ_MUMPS_contact_init

  subroutine solve_LINEQ_MUMPS_contact(hecMESH,hecMAT,hecLagMAT,istat,conMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (hecmwST_matrix_lagrange), intent(inout) :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint), intent(out) :: istat
    type (hecmwST_matrix), intent(in) :: conMAT

    integer(kind=kint) :: mumps_job, mpc_method

    mpc_method = hecmw_mat_get_mpc_method(hecMAT)
    if (mpc_method < 1 .or. 3 < mpc_method) then
      mpc_method = 1
      call hecmw_mat_set_mpc_method(hecMAT,mpc_method)
    endif
    if (mpc_method /= 1) then
      write(*,*) 'ERROR: MPCMETHOD other than penalty is not available for MUMPS solver', &
          ' in contact analysis without elimination'
      stop
    endif
    call hecmw_mat_ass_equation(hecMESH, hecMAT)
    call hecmw_mat_ass_equation_rhs(hecMESH, hecMAT)

    call hecmw_mat_dump(hecMAT, hecMESH)

    ! ANALYSIS
    if (NEED_ANALYSIS) then
      call sparse_matrix_contact_init_prof(spMAT, hecMAT, hecLagMAT, hecMESH )
      mumps_job=1
      call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MUMPS returned with error', istat
        stop
      endif
      if ( hecmw_comm_get_rank()==0 ) write(*,*) ' [MUMPS]: Analysis completed.'
      NEED_ANALYSIS = .false.
    endif

    ! FACTORIZATION and SOLUTION
    call sparse_matrix_para_contact_set_vals(spMAT, hecMAT, hecLagMAT, conMAT)
    call sparse_matrix_para_contact_set_rhs(spMAT, hecMAT, hecLagMAT, conMAT)
    mumps_job=5
    call hecmw_mumps_wrapper(spMAT, mumps_job, istat)
    if (istat < 0) then
      write(*,*) 'ERROR: MUMPS returned with error', istat
      return
    endif
    call sparse_matrix_contact_get_rhs(spMAT, hecMAT, hecLagMAT)
    if (hecmw_comm_get_rank()==0) write(*,*) ' [MUMPS]: Factorization and Solution completed.'

    call hecmw_mat_dump_solution(hecMAT)

    !call sparse_matrix_finalize(spMAT)
  end subroutine solve_LINEQ_MUMPS_contact

end module m_solve_LINEQ_MUMPS_contact
