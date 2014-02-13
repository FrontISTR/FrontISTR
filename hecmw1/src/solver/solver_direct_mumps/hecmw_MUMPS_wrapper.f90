!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
!> This module provides wrapper for parallel sparse direct solver MUMPS
module m_hecmw_MUMPS_wrapper
  use hecmw_util
  use m_hecmw_comm_f
  use m_sparse_matrix
  include 'dmumps_struc.h'

  private
  public :: hecmw_mumps_wrapper

  type (dmumps_struc), save :: mumps_par

contains

  subroutine hecmw_mumps_wrapper(spMAT, job, istat)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    integer(kind=kint), intent(out) :: istat
    integer(kind=kint) :: ierr,myrank

    myrank=hecmw_comm_get_rank()

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
          else
            allocate(mumps_par%RHS(1), stat=ierr)
          endif
          if (ierr /= 0) then
            write(*,*) " Allocation error, mumps_par%RHS"
            call hecmw_abort(hecmw_comm_get_comm())
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
       mumps_par%ICNTL(4)=0
    endif

    mumps_par%JOB = job
    do
       call DMUMPS(mumps_par)
       istat = mumps_par%INFOG(1)
       if (istat >= 0) exit
       if (istat == -9 .and. mumps_par%ICNTL(14) < 200) then
          mumps_par%ICNTL(14) = mumps_par%ICNTL(14) + 20
          if (myrank == 0) &
               write(*,*) 'INFO: MUMPS increasing relaxation parameter to', &
               mumps_par%ICNTL(14)
       elseif (istat < 0) then
          if (myrank == 0) then
             write(*,*) 'ERROR: MUMPS job=',job,&
                  ', INFOG(1)=',istat,', INFOG(2)=',mumps_par%INFOG(2)
          endif
          return
       endif
    enddo

    if (job==-1) then
       ! ordering: 0:auto, 1:seq, 2:par
       mumps_par%ICNTL(28)=0
       ! seq ord: 0:AMD, 1:USER, 2:AMF, 3:scotch, 4:pord, 5:metis, 6:QAMD, 7:auto
       mumps_par%ICNTL(7)=7
       ! par ord: 0:auto, 1:ptscotch, 2:parmetis
       mumps_par%ICNTL(29)=0
       ! relaxation parameter
       mumps_par%ICNTL(14)=20
       ! iterative refinement
       mumps_par%ICNTL(10)=3
       mumps_par%CNTL(2)=1.0e-8
       ! Out-Of-Core: 0:IN-CORE only, 1:OOC
       mumps_par%ICNTL(22)=0
    endif
    if (job==3 .or. job==5 .or. job==6) then
       call sparse_matrix_scatter_rhs(spMAT, mumps_par%RHS)
       deallocate(mumps_par%RHS)
    endif
  end subroutine hecmw_mumps_wrapper

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

end module m_hecmw_MUMPS_wrapper
