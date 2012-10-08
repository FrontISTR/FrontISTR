!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.4                                   !
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
  use hecmw_util
  use m_hecmw_comm_f
  use m_sparse_matrix
  include 'dmumps_struc.h'

  private
  public :: mumps_wrapper

  type (dmumps_struc), save :: mumps_par
  real(kreal),public  ::  rhs_b,rhs_x

contains

  subroutine mumps_wrapper(spMAT, job, paraContactFlag, istat)
    implicit none
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    logical, intent(in) :: paraContactFlag
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
          if(paraContactFlag) then
            if(myrank == 0) then
              mumps_par%RHS(:) = mumps_par%RHS(:) + spMAT%rhs_con_sum(:)
!              deallocate(spMAT%rhs_con_sum,stat=ierr)
              rhs_b = dot_product(mumps_par%RHS,mumps_par%RHS)
            endif
            deallocate(spMAT%rhs_con_sum,stat=ierr)
            call hecmw_bcast_R1_comm (rhs_b, 0, mumps_par%COMM)
!            call MPI_BCAST(rhs_b,1,MPI_DOUBLE_PRECISION,0,mumps_par%COMM,ierr)
          else
            if(myrank == 0) then
              rhs_b = dot_product(mumps_par%RHS,mumps_par%RHS)
            endif
          endif
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
    endif
    if (job==3 .or. job==5 .or. job==6) then
       if(myrank == 0) then
         rhs_x = dot_product(mumps_par%RHS,mumps_par%RHS)
       endif
       call hecmw_bcast_R1_comm (rhs_x, 0, mumps_par%COMM)
!       call MPI_BCAST(rhs_x,1,MPI_DOUBLE_PRECISION,0,mumps_par%COMM,ierr)
       call sparse_matrix_scatter_rhs(spMAT, mumps_par%RHS)
       deallocate(mumps_par%RHS)
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
