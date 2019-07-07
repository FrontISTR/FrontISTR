!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface for Pardiso
#ifdef WITH_MKL
include 'mkl_pardiso.f90'
#endif

module m_hecmw_MKL_wrapper
  use hecmw_util
  use m_sparse_matrix

#ifdef WITH_MKL
  use mkl_pardiso
#endif

  implicit none

  private                         ! default
  public :: hecmw_mkl_wrapper     ! only entry point of Parallel Direct Solver is public

  logical, save :: INITIALIZED = .false.
#ifdef WITH_MKL
  type(MKL_PARDISO_HANDLE) :: pt(64)
#endif
  integer maxfct, mnum, mtype, nrhs, msglvl
  integer idum(1), i
  data nrhs /1/, maxfct /1/, mnum /1/

  integer, save :: iparm(64)

contains

  subroutine hecmw_mkl_wrapper(spMAT, phase_start, solx, istat)
    implicit none
    type (sparse_matrix), intent(inout)      :: spMAT
    integer(kind=kint), intent(in)           :: phase_start
    integer(kind=kint), intent(out)          :: istat
    real(kind=kreal), pointer, intent(inout) :: solx(:)

    integer(kind=kint) :: myrank, phase
    real(kind=kreal)   :: t2,t3,t4,t5

#ifdef WITH_MKL

    myrank=hecmw_comm_get_rank()

    if (.not. INITIALIZED) then
      do i=1,64
        pt(i)%dummy = 0
      enddo
      iparm(:) = 0
      iparm(1) = 1 ! no solver default
      iparm(2) = 3 ! fill-in reordering from METIS
      iparm(3) = 1
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(13) = 1 ! maximum weighted matching algorithm is switched-off
      iparm(18) = -1
      iparm(19) = -1
      msglvl = 0 ! print statistical information
      INITIALIZED = .true.
    endif

    if ( phase_start == 1 ) then
      phase = -1
      call pardiso(pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, idum,  &
        nrhs, iparm, msglvl, spMAT%RHS, solx, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat
        return
      endif

      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) then
        mtype = 2
      else if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) then
        mtype = -2
      else
        mtype = 11
      endif

      ! ANALYSIS
      t2=hecmw_wtime()

      phase = 11
      call pardiso(pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, idum,  &
        nrhs, iparm, msglvl, spMAT%RHS, solx, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat
        return
      endif

      t3=hecmw_wtime()
      if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
        write(*,'(A,f10.3)') ' [Pardiso]: Analysis completed.       time(sec)=',t3-t2
    endif

    ! FACTORIZATION
    if ( phase_start .le. 2 ) then
      t3=hecmw_wtime()

      phase = 22
      call pardiso(pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, idum,  &
        nrhs, iparm, msglvl, spMAT%RHS, solx, istat)
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat
        return
      endif

      t4=hecmw_wtime()
      if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
        write(*,'(A,f10.3)') ' [Pardiso]: Factorization completed.  time(sec)=',t4-t3
    endif

    ! SOLUTION
    t4=hecmw_wtime()
    phase = 33
    call pardiso(pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, idum,  &
      nrhs, iparm, msglvl, spMAT%RHS, solx, istat)
    if (istat < 0) then
      write(*,*) 'ERROR: MKL returned with error', istat
      return
    endif

    t5=hecmw_wtime()
    if (myrank==0 .and. (spMAT%iterlog > 0 .or. spMAT%timelog > 0)) &
       write(*,'(A,f10.3)') ' [Pardiso]: Solution completed.       time(sec)=',t5-t4

#else
    stop "MKL Pardiso not available"
#endif
  end subroutine hecmw_mkl_wrapper

end module m_hecmw_MKL_wrapper
