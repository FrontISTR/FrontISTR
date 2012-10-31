!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by K. Goto (VINAS)                                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module provides linear equation solver interface for MUMPS
module m_solve_LINEQ_MUMPS
  use m_fstr
  use m_sparse_matrix
  use m_sparse_matrix_hec
  use m_MUMPS_wrapper

  private
  public :: solve_LINEQ_MUMPS

  logical, save :: INITIALIZED = .false.
  type (sparse_matrix), save :: spMAT

contains

  subroutine solve_LINEQ_MUMPS(hecMESH,hecMAT)
    implicit none
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: mumps_job
    integer(kind=kint) :: istat

    if (INITIALIZED .and. hecMAT%Iarray(98) .eq. 1) then
       mumps_job=-2
       call mumps_wrapper(spMAT, mumps_job, istat)
       if (istat < 0) then
         write(*,*) 'ERROR: MUMPS returned with error', istat
         stop
       endif
       call sparse_matrix_finalize(spMAT)
       INITIALIZED = .false.
    endif

    if (.not. INITIALIZED) then
       spmat_type = SPARSE_MATRIX_TYPE_COO
       spmat_symtype = SPARSE_MATRIX_SYMTYPE_SPD
       call sparse_matrix_set_type(spMAT, spmat_type, spmat_symtype)
       mumps_job = -1
       call mumps_wrapper(spMAT, mumps_job, istat)
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
       call mumps_wrapper(spMAT, mumps_job, istat)
       if (istat < 0) then
         write(*,*) 'ERROR: MUMPS returned with error', istat
         stop
       endif
       if (myrank==0) write(*,*) ' [MUMPS]: Analysis and Factorization completed.'
       hecMAT%Iarray(98) = 0
       hecMAT%Iarray(97) = 0
    endif
    if (hecMAT%Iarray(97) .eq. 1) then
       ! FACTORIZATION
       call sparse_matrix_hec_set_vals(spMAT, hecMAT)
       !call sparse_matrix_dump(spMAT)
       mumps_job=2
       call mumps_wrapper(spMAT, mumps_job, istat)
       if (istat < 0) then
         write(*,*) 'ERROR: MUMPS returned with error', istat
         stop
       endif
       if (myrank==0) write(*,*) ' [MUMPS]: Factorization completed.'
       hecMAT%Iarray(97) = 0
    endif
    ! SOLUTION
    call sparse_matrix_hec_set_rhs(spMAT, hecMAT)
    mumps_job=3
    call mumps_wrapper(spMAT, mumps_job, istat)
    if (istat < 0) then
      write(*,*) 'ERROR: MUMPS returned with error', istat
      stop
    endif
    call sparse_matrix_hec_get_rhs(spMAT, hecMAT)
    if (myrank==0) write(*,*) ' [MUMPS]: Solution completed.'

    !call sparse_matrix_finalize(spMAT)
  end subroutine solve_LINEQ_MUMPS

end module m_solve_LINEQ_MUMPS
