!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides linear equation solver interface for Cluster Pardiso
#ifdef WITH_MKL
include 'mkl_cluster_sparse_solver.f90'
#endif

module m_hecmw_ClusterMKL_wrapper
  use hecmw_util
  use m_sparse_matrix

#ifdef WITH_MKL
  use mkl_cluster_sparse_solver
#endif

  implicit none

  private                         ! default
  public hecmw_clustermkl_wrapper ! only entry point of Parallel Direct Solver is public

  logical, save :: INITIALIZED = .false.
#ifdef WITH_MKL
  type(MKL_CLUSTER_SPARSE_SOLVER_HANDLE) :: pt(64)
#endif
  integer maxfct, mnum, mtype, nrhs, msglvl
  integer idum(1), i
  data nrhs /1/, maxfct /1/, mnum /1/

  integer, save :: iparm(64)

contains

  subroutine hecmw_clustermkl_wrapper(spMAT, phase_start, solx, istat)
    implicit none
    type (sparse_matrix), intent(inout)      :: spMAT
    integer(kind=kint), intent(in)           :: phase_start
    integer(kind=kint), intent(out)          :: istat
    real(kind=kreal), pointer, intent(inout) :: solx(:)

    integer(kind=kint) :: myrank, phase
    real(kind=kreal)   :: t1,t2,t3,t4,t5

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
      iparm(11) = 0 ! use nonsymmetric permutation and scaling MPS
      iparm(13) = 0 ! maximum weighted matching algorithm is switched-off
      iparm(18) = -1
      iparm(19) = -1
      msglvl = 0 ! print statistical information
      INITIALIZED = .true.
    else
      if( phase_start == 1 ) then
        phase = -1
        call cluster_sparse_solver ( &
          pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, &
          idum, nrhs, iparm, msglvl, spMAT%RHS, solx, hecmw_comm_get_comm(), istat )
        if (istat < 0) then
          write(*,*) 'ERROR: MKL returned with error in phase -1', istat
          return
        endif
      end if
    endif

    iparm(40) = 2 ! Input: distributed matrix/rhs/solution format
    iparm(41) = spMAT%DISPLS(myrank+1)+1
    iparm(42) = spMAT%DISPLS(myrank+1)+spMAT%N_COUNTS(myrank+1)

    ! Additional setup
    t1=hecmw_wtime()
    if ( phase_start == 1 ) then
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) then
        mtype = 2
      else if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) then
        mtype = -2
      else
        mtype = 11
      endif
    endif

    call sort_column_ascending_order(spMAT,myrank)

    t2=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) &
      write(*,'(A,f10.3)') ' [Cluster Pardiso]: Additional Setup completed. time(sec)=',t2-t1

    if ( phase_start == 1 ) then
      ! ANALYSIS
      t2=hecmw_wtime()

      phase = 11
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, &
        idum, nrhs, iparm, msglvl, spMAT%RHS, solx, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error in phase 1', istat
        return
      endif

      t3=hecmw_wtime()
      if (myrank==0 .and. spMAT%timelog > 0) &
        write(*,'(A,f10.3)') ' [Cluster Pardiso]: Analysis completed.         time(sec)=',t3-t2
    endif

    ! FACTORIZATION
    if ( phase_start .le. 2 ) then
      t3=hecmw_wtime()

      phase = 22
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, &
        idum, nrhs, iparm, msglvl, spMAT%RHS, solx, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error in phase 2', istat
        stop
      endif
      t4=hecmw_wtime()
      if (myrank==0 .and. spMAT%timelog > 0) &
        write(*,'(A,f10.3)') ' [Cluster Pardiso]: Factorization completed.    time(sec)=',t4-t3
    endif

    t4=hecmw_wtime()

    ! SOLUTION
    phase = 33
    call cluster_sparse_solver ( &
      pt, maxfct, mnum, mtype, phase, spMAT%N, spMAT%A, spMAT%IRN, spMAT%JCN, &
      idum, nrhs, iparm, msglvl, spMAT%RHS, solx, hecmw_comm_get_comm(), istat )
    if (istat < 0) then
      write(*,*) 'ERROR: MKL returned with error in phase 3', istat
      stop
    endif

    ! solution is written to hecMAT%X. you have to edit lib/solve_LINEQ.f90
    ! to use this routine properly.

    t5=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) &
       write(*,'(A,f10.3)') ' [Cluster Pardiso]: Solution completed.         time(sec)=',t5-t4

#else
    stop "MKL Pardiso not available"
#endif
  end subroutine hecmw_clustermkl_wrapper

#ifdef WITH_MKL
  subroutine sort_column_ascending_order(spMAT,myrank)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in)   :: myrank

    integer(kind=kint) :: i, j, is, iE

    do i=2,spMAT%N_loc+1
      if(spMAT%IRN(i)-spMAT%IRN(i-1)>10000) &
        &  write(*,*) 'warning: n may be too large for set_mkl_set_prof0 for row ',i-1
    enddo

    !!$omp parallel private(i, iS, iE)
    !!$omp do
    do i=1,spMAT%N_loc
      is = spMAT%IRN(i)
      iE = spMAT%IRN(i+1)-1
      call sort_column_ascending_order0(iE-is+1,spMAT%JCN(is:iE),spMAT%A(is:iE))
    enddo
    !!$omp end do
    !!$omp end parallel
  end subroutine sort_column_ascending_order

  subroutine sort_column_ascending_order0(n,ja,a)
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(inout) :: ja(:)
    real(kind=kreal), intent(inout) :: a(:)

    integer(kind=kint) :: i, j, work(2,n)
    real(kind=kreal)   :: oa(n)

    do i=1,n
      work(1,i) = ja(i)
      work(2,i) = i
      oa(i) = a(i)
    enddo

    call myqsort(n,work)

    do i=1,n
      ja(i) = work(1,i)
      a(i) = oa(work(2,i))
    enddo

  end subroutine sort_column_ascending_order0

  recursive subroutine myqsort(n,work)
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(inout) :: work(:,:)

    integer(kind=kint) :: i, j, key, work1(2,n), n_low, n_high
    logical :: sorted

    if(n<2) return
    !return if work is already sorted
    sorted = .true.
    do i=1,n-1
      if(work(1,i)>work(1,i+1)) then
        sorted = .false.
        key = (work(1,i)+work(1,i+1))/2
        exit
      endif
    enddo
    if(sorted) return
    if(n<5) then
      call myinssort(n,work)
      return
    endif

    n_low=0
    do i=1,n
      if(work(1,i)<key) then
        n_low=n_low+1
        work1(1:2,n_low) = work(1:2,i)
      endif
    enddo

    n_high=0
    do i=1,n
      if(work(1,i)>=key) then
        n_high=n_high+1
        work1(1:2,n_low+n_high) = work(1:2,i)
      endif
    enddo

    if(n_low>0) call myqsort(n_low,work1(1:2,1:n_low))
    if(n_high>0) call myqsort(n_high,work1(1:2,n_low+1:n_low+n_high))

    do i=1,n
      work(1:2,i) = work1(1:2,i)
    enddo

  end subroutine myqsort

  subroutine myinssort(n,work)
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(inout) :: work(:,:)

    integer(kind=kint) :: i,j,tmp(2)

    do i=2,n
      do j=i,2,-1
        if(work(1,j)<work(1,j-1)) then
          tmp(1:2) = work(1:2,j)
          work(1:2,j) = work(1:2,j-1)
          work(1:2,j-1) = tmp(1:2)
        endif
      enddo
    enddo

  end subroutine myinssort

#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_hecmw_ClusterMKL_wrapper
