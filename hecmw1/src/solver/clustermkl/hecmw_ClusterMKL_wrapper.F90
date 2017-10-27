!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
#ifdef WITH_MKL
include 'mkl_cluster_sparse_solver.f90'
#endif

module m_hecmw_ClusterMKL_wrapper
  use hecmw_util

#ifdef WITH_MKL
  use m_sparse_matrix
  use m_sparse_matrix_hec
  use hecmw_matrix_ass
  use hecmw_matrix_dump
  use mkl_cluster_sparse_solver
#endif

  implicit none

  private                            ! default
  public hecmw_clustermkl_wrapper ! only entry point of Parallel Direct Solver is public

#ifdef WITH_MKL
  logical, save :: INITIALIZED = .false.
  type (sparse_matrix), save :: spMAT
  type(MKL_CLUSTER_SPARSE_SOLVER_HANDLE) :: pt(64)
  integer maxfct, mnum, mtype, phase, nrhs, error, msglvl
  integer idum(1), n, nnz, nloc
  real(kind=kreal)  ddum(1)
  data nrhs /1/, maxfct /1/, mnum /1/
  integer, parameter :: imsg=51

  integer :: iparm(64)
  integer, allocatable, dimension(:) :: ia
  integer, allocatable, dimension(:) :: ja
  real*8, allocatable, dimension(:) :: a
  real*8, allocatable, dimension(:) :: b
  real*8, allocatable, dimension(:) :: x

  integer myrank0
#endif

contains

  subroutine hecmw_clustermkl_wrapper(hecMESH, hecMAT)
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT

#ifdef WITH_MKL
    integer(kind=kint) :: spmat_type
    integer(kind=kint) :: spmat_symtype
    integer(kind=kint) :: istat,myrank
    real(kind=kreal) :: t1,t2,t3,t4,t5

    call hecmw_mat_dump(hecMAT, hecMESH)

    t1=hecmw_wtime()
    myrank=hecmw_comm_get_rank()

    if (INITIALIZED .and. hecMAT%Iarray(98) .eq. 1) then
      phase = -1
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
        idum, nrhs, iparm, msglvl, ddum, ddum, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat
        stop
      endif
      call sparse_matrix_finalize(spMAT)
      if(allocated(a)) deallocate(a)
      if(allocated(ia)) deallocate(ia)
      if(allocated(ja)) deallocate(ja)
      if(allocated(b)) deallocate(b)
      if(allocated(x)) deallocate(x)
      INITIALIZED = .false.
    endif

    if (.not. INITIALIZED) then
      spmat_type = SPARSE_MATRIX_TYPE_CSR
      spmat_symtype = SPARSE_MATRIX_SYMTYPE_SYM
      !spmat_symtype = SPARSE_MATRIX_SYMTYPE_ASYM
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
      call setup_mkl_parameters(spMAT,myrank)
      call setup_mkl_set_val(spMAT,myrank,.false.)
      t2=hecmw_wtime()
      phase = 11
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, n, spMAT%A, spMAT%IRN, spMAT%JCN, &
        idum, nrhs, iparm, msglvl, spMAT%rhs, x, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat, myrank
        stop
      endif
      t3=hecmw_wtime()
      phase = 22
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, n, spMAT%A, spMAT%IRN, spMAT%JCN, &
        idum, nrhs, iparm, msglvl, spMAT%rhs, x, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat
        stop
      endif
      !if (myrank==0) write(*,*) ' [MKL]: Analysis and Factorization completed.'
      hecMAT%Iarray(98) = 0
      hecMAT%Iarray(97) = 0
    endif
    if (hecMAT%Iarray(97) .eq. 1) then
      ! FACTORIZATION
      call sparse_matrix_hec_set_prof(spMAT, hecMAT)
      call sparse_matrix_hec_set_vals(spMAT, hecMAT)
      call setup_mkl_set_val(spMAT,myrank,.true.)
      !call sparse_matrix_dump(spMAT)
      t2=hecmw_wtime()
      t3=hecmw_wtime()
      phase = 22
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, n, spMAT%A, spMAT%IRN, spMAT%JCN, &
        idum, nrhs, iparm, msglvl, spMAT%rhs, x, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        write(*,*) 'ERROR: MKL returned with error', istat
        stop
      endif
      !if (myrank==0) write(*,*) ' [MKL]: Factorization completed.'
      hecMAT%Iarray(97) = 0
    endif

    t4=hecmw_wtime()

    ! SOLUTION
    call sparse_matrix_hec_set_rhs(spMAT, hecMAT)
    phase = 33
    call cluster_sparse_solver ( &
      pt, maxfct, mnum, mtype, phase, n, spMAT%A, spMAT%IRN, spMAT%JCN, &
      idum, nrhs, iparm, msglvl, spMAT%rhs, x, hecmw_comm_get_comm(), istat )
    if (istat < 0) then
      write(*,*) 'ERROR: MKL returned with error', istat
      stop
    endif
    call sparse_matrix_hec_get_rhs(spMAT, hecMAT)
    ! if (myrank==0) write(*,*) ' [MKL]: Solution completed.'
    ! solution is written to hecMAT%X. you have to edit lib/solve_LINEQ.f90
    ! to use this routine properly.

    t5=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) then
      write(imsg,'(A)') '############## MATRIX SOLVER TIME INFORMATION ##############'
      write(imsg,*) ' setup time         : ',t2-t1
      write(imsg,*) ' ordering time      : ',t3-t2
      write(imsg,*) ' factorization time : ',t4-t3
      write(imsg,*) ' solve time         : ',t5-t4
      write(imsg,'(A,i2)') 'number of performed iterative refinement steps',iparm(7)
      write(imsg,'(A,i8)') 'number of positive eigenvalues',iparm(22)
      write(imsg,'(A,i8)') 'number of negative eigenvalues',iparm(23)
      write(imsg,'(A,i8)') 'number of perturbed pivots',iparm(14)
      write(imsg,'(A)') '############ END MATRIX SOLVER TIME INFORMATION ############'
    endif

    call hecmw_mat_dump_solution(hecMAT)

#else
    stop "PARADISO not available"
#endif
  end subroutine hecmw_clustermkl_wrapper

#ifdef WITH_MKL

  subroutine hecmw_MKL_wrapper(spMAT, job, istat)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in) :: job
    integer(kind=kint), intent(inout) :: istat
    integer(kind=kint) :: ierr,myrank

    integer(kind=kint) :: phase

    myrank=hecmw_comm_get_rank()

  end subroutine hecmw_MKL_wrapper

  subroutine setup_mkl_parameters(spMAT,myrank)
    type (sparse_matrix), intent(in) :: spMAT
    integer(kind=kint), intent(in)   :: myrank

    integer(kind=kint) :: i, j

    n = spMAT%N
    nnz = spMAT%NZ
    nloc = spMAT%N_loc
    nrhs = 1
    do i=1,64
      pt(i)%dummy = 0
    enddo
    iparm(:) = 0
    iparm(1) = 1 ! no solver default
    iparm(2) = 3 ! fill-in reordering from METIS
    iparm(6) = 1 ! =0 solution on the first n compoments of x
    iparm(8) = 2 ! numbers of iterative refinement steps
    iparm(10) = 13 ! perturbe the pivot elements with 1E-13
    iparm(11) = 0 ! use nonsymmetric permutation and scaling MPS
    iparm(13) = 0 ! maximum weighted matching algorithm is switched-off
    iparm(27) = 0 ! check matrix(temp)
    iparm(40) = 2 ! Input: distributed matrix/rhs/solution format
    error = 0  ! initialize error flag
    msglvl = 0 ! print statistical information

    if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SPD) then
      mtype = 2
    else if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_SYM) then
      mtype = -2
    else
      mtype = 11
    endif
    iparm(41) = spMAT%DISPLS(myrank+1)+1
    iparm(42) = spMAT%DISPLS(myrank+1)+spMAT%N_COUNTS(myrank+1)

  end subroutine setup_mkl_parameters

  subroutine setup_mkl_set_val(spMAT,myrank,ordered)
    type (sparse_matrix), intent(inout) :: spMAT
    integer(kind=kint), intent(in)   :: myrank
    logical, intent(in)   :: ordered

    integer(kind=kint) :: i, j, is, iE

    do i=2,nloc+1
      if(spMAT%IRN(i)-spMAT%IRN(i-1)>10000) &
        &  write(*,*) 'warning: n may be too large for set_mkl_set_prof0 for row ',i-1
    enddo

    !!$omp parallel private(i, iS, iE)
    !!$omp do
    do i=1,nloc
      is = spMAT%IRN(i)
      iE = spMAT%IRN(i+1)-1
      call set_mkl_set_val0(iE-is+1,spMAT%JCN(is:iE),spMAT%A(is:iE))
    enddo
    !!$omp end do
    !!$omp end parallel
  end subroutine setup_mkl_set_val

  subroutine set_mkl_set_val0(n,ja,a)
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

  end subroutine set_mkl_set_val0

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
