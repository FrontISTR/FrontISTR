!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to solve sparse system of
!> \linear equitions using intel MKL direct sparse solver

module m_solve_LINEQ_MKL_contact
  use hecmw_util
  use m_fstr
  use m_sparse_matrix
  use m_sparse_matrix_contact
  use fstr_matrix_con_contact
  use m_hecmw_MKL_wrapper
  use m_hecmw_ClusterMKL_wrapper

  implicit none

  private
  public :: solve_LINEQ_MKL_contact_init
  public :: solve_LINEQ_MKL_contact

  logical, save :: NEED_ANALYSIS = .true.
  type (sparse_matrix), save :: spMAT

contains

  subroutine solve_LINEQ_MKL_contact_init(hecMESH,is_sym)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    logical, intent(in) :: is_sym

    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: i

    call sparse_matrix_finalize(spMAT)

    if (is_sym) then
      spmat_symtype = SPARSE_MATRIX_SYMTYPE_SYM
    else
      spmat_symtype = SPARSE_MATRIX_SYMTYPE_ASYM
    end if
    if(hecMESH%PETOT.GT.1) then
      spmat_type = SPARSE_MATRIX_TYPE_COO
    else
      spmat_type = SPARSE_MATRIX_TYPE_CSR
    endif
    call sparse_matrix_set_type(spMAT, spmat_type, spmat_symtype)

    NEED_ANALYSIS = .true.
  end subroutine solve_LINEQ_MKL_contact_init

  !> \brief This subroutine executes the MKL solver
  subroutine solve_LINEQ_MKL_contact(hecMESH,hecMAT,fstrMAT,istat,conMAT)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    type (fstrST_matrix_contact_lagrange), intent(inout) :: fstrMAT !< type fstrST_matrix_contact_lagrange
    integer(kind=kint), intent(out) :: istat
    type (hecmwST_matrix), intent(in),optional :: conMAT

    integer(kind=kint)  :: phase_start
    real(kind=kreal)    :: t1,t2

    t1=hecmw_wtime()
    call hecmw_mat_dump(hecMAT, hecMESH)

    if (NEED_ANALYSIS) then
      !constrtuct new structure
      call sparse_matrix_contact_init_prof(spMAT, hecMAT, fstrMAT, hecMESH)
    endif

    !  ----  For Parallel Contact with Multi-Partition Domains
    if(paraContactFlag.and.present(conMAT)) then
      call sparse_matrix_para_contact_set_vals(spMAT, hecMAT, fstrMAT, conMAT)
      call sparse_matrix_para_contact_set_rhs(spMAT, hecMAT, fstrMAT, conMAT)
    else
      call sparse_matrix_contact_set_vals(spMAT, hecMAT, fstrMAT)
      !call sparse_matrix_dump(spMAT)
      call sparse_matrix_contact_set_rhs(spMAT, hecMAT, fstrMAT)
    endif

    t2=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) then
      if( hecMESH%PETOT .GT. 1 ) then
         write(*,'(A,f10.3)') ' [Cluster Pardiso]: Setup completed.            time(sec)=',t2-t1
      else
         write(*,'(A,f10.3)') ' [Pardiso]: Setup completed.          time(sec)=',t2-t1
      end if
    endif

    phase_start = 2
    if (NEED_ANALYSIS) then
      phase_start = 1
      NEED_ANALYSIS = .false.
    endif

    ! SOLVE
    if( hecMESH%PETOT.GT.1 ) then
      call hecmw_clustermkl_wrapper(spMAT, phase_start, hecMAT%X, istat)
      call sparse_matrix_contact_get_rhs(spMAT, hecMAT, fstrMAT)
    else
      call hecmw_mkl_wrapper(spMAT, phase_start, hecMAT%X, istat)
    endif

    call hecmw_mat_dump_solution(hecMAT)
  end subroutine solve_LINEQ_MKL_contact

end module m_solve_LINEQ_MKL_contact
