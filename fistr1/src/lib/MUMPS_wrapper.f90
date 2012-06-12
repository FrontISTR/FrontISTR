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
!> This module provides wrapper for parallel sparse direct solver MUMPS
module m_MUMPS_wrapper
  use m_fstr
  use m_sparse_matrix
  include 'dmumps_struc.h'

  private
  public :: mumps_wrapper

  type (dmumps_struc), save :: mumps_par

contains

  subroutine mumps_wrapper(spMAT, job, istat)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    integer(kind=kint), intent(out) :: istat
    integer(kind=kint) :: ierr

    if (spMAT%type /= SPARSE_MATRIX_TYPE_COO) then
       write(*,*) 'ERROR: MUMPS require COO type sparse matrix'
       call hecmw_abort(hecmw_comm_get_comm())
    endif

    if (job==-1) then
       mumps_par%COMM = hecmw_comm_get_comm()
       mumps_par%JOB = -1
       if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) then
          mumps_par%SYM = 1
       else if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) then
          mumps_par%SYM = 2
       else
          mumps_par%SYM = 0
       endif
       mumps_par%PAR = 1
    elseif (job>0) then
       call set_mumps_pointers(mumps_par, spMAT)
       if (job==3 .or. job==5 .or. job==6) then
          if (myrank == 0) then
            allocate(mumps_par%RHS(mumps_par%N), stat=ierr)
            if (ierr /= 0) then
              write(*,*) " Allocation error, mumps_par%RHS"
              call hecmw_abort(hecmw_comm_get_comm())
            endif
          endif
          call sparse_matrix_gather_rhs(spMAT, mumps_par%RHS)
       endif
    endif

    if (spMAT%timelog > 0) then
       mumps_par%ICNTL(1)=6
       mumps_par%ICNTL(2)=0
       mumps_par%ICNTL(3)=6
       mumps_par%ICNTL(4)=1
    else
       mumps_par%ICNTL(1)=6
       mumps_par%ICNTL(2)=0
       mumps_par%ICNTL(3)=0
       mumps_par%ICNTL(4)=1
    endif

    mumps_par%JOB = job
    call DMUMPS(mumps_par)
    istat = mumps_par%INFOG(1)
    if (istat < 0) then
       write(*,*) 'ERROR: MUMPS job=',job
       return
    endif

    if (job==3 .or. job==5 .or. job==6) then
       call sparse_matrix_scatter_rhs(spMAT, mumps_par%RHS)
       if (myrank == 0) deallocate(mumps_par%RHS)
    endif
  end subroutine mumps_wrapper

  subroutine set_mumps_pointers(mumps_par, spMAT)
    implicit none
    type (dmumps_struc), intent(inout) :: mumps_par
    type (sparse_matrix), intent(in) :: spMAT
    mumps_par%N = spMAT%N
    ! Distributed assembled matrix input
    mumps_par%ICNTL(18) = 3
    mumps_par%NZ_loc = spMAT%NZ
    mumps_par%IRN_loc => spMAT%IRN
    mumps_par%JCN_loc => spMAT%JCN
    mumps_par%A_loc => spMAT%A
  end subroutine set_mumps_pointers

end module m_MUMPS_wrapper
