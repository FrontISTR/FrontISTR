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
  use m_hecmw_comm_f
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
  integer(kind=kint), save          :: nn
  integer(kind=kint), pointer, save :: irow(:), jcol(:)
  real(kind=kreal), pointer, save   :: aval(:), rhs(:), solx(:)

  integer(kind=kint), parameter :: debug=0

contains

  subroutine hecmw_clustermkl_wrapper(spMAT, phase_start, solution, istat)
    implicit none
    type (sparse_matrix), intent(inout)      :: spMAT
    integer(kind=kint), intent(in)           :: phase_start
    integer(kind=kint), intent(out)          :: istat
    real(kind=kreal), pointer, intent(inout) :: solution(:)

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
      iparm(8) = 2
      if (spMAT%symtype == SPARSE_MATRIX_SYMTYPE_ASYM) then
        iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      else
        iparm(10) = 8 ! perturbe the pivot elements with 1E-8
      endif
      iparm(11) = 1 ! Enable scaling
      iparm(13) = 1 ! Enable matching
      iparm(18) = 0
      iparm(19) = 0
      msglvl = 0 ! print statistical information
      INITIALIZED = .true.
    else
      if( phase_start == 1 ) then
        phase = -1
        call cluster_sparse_solver ( &
          pt, maxfct, mnum, mtype, phase, nn, aval, irow, jcol, &
          idum, nrhs, iparm, msglvl, rhs, solx, hecmw_comm_get_comm(), istat )
        if (istat < 0) then
          write(*,*) 'ERROR: MKL returned with error in phase -1', istat
          return
        endif
        if( spMAT%type == SPARSE_MATRIX_TYPE_COO ) deallocate(irow,jcol,aval,rhs,solx)
      end if
    endif

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

    if( spMAT%type == SPARSE_MATRIX_TYPE_CSR ) then
      iparm(40) = 2 ! Input: distributed matrix/rhs/solution format
      iparm(41) = spMAT%DISPLS(myrank+1)+1
      iparm(42) = spMAT%DISPLS(myrank+1)+spMAT%N_COUNTS(myrank+1)
      call sort_column_ascending_order(spMAT,myrank)
      nn = spMAT%N
      irow => spMAT%IRN
      jcol => spMAT%JCN
      aval => spMAT%A
      rhs  => spMAT%rhs
      solx => solution
    else if( spMAT%type == SPARSE_MATRIX_TYPE_COO ) then
      iparm(40) = 0 ! Input: centralized input format( mkl don't have distributed nonsymmetric matrix input...)
      !input matrix is gathered to rank0 if matrix is given in COO format
      call export_spMAT_to_CentralizedCRS(spMAT,myrank,nn,irow,jcol,aval,rhs,solx)
    end if

    t2=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) &
      write(*,'(A,f10.3)') ' [Cluster Pardiso]: Additional Setup completed. time(sec)=',t2-t1

    if ( phase_start == 1 ) then
      ! ANALYSIS
      t2=hecmw_wtime()

      phase = 11
      call cluster_sparse_solver ( &
        pt, maxfct, mnum, mtype, phase, nn, aval, irow, jcol, &
        idum, nrhs, iparm, msglvl, rhs, solx, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        if (myrank==0 .and. spMAT%timelog > 0) call print_iparm_paramters()
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
        pt, maxfct, mnum, mtype, phase, nn, aval, irow, jcol, &
        idum, nrhs, iparm, msglvl, rhs, solx, hecmw_comm_get_comm(), istat )
      if (istat < 0) then
        if (myrank==0 .and. spMAT%timelog > 0) call print_iparm_paramters()
        write(*,*) 'ERROR: MKL returned with error in phase 2', istat
        return
      endif
      t4=hecmw_wtime()
      if (myrank==0 .and. spMAT%timelog > 0) &
        write(*,'(A,f10.3)') ' [Cluster Pardiso]: Factorization completed.    time(sec)=',t4-t3
    endif

    t4=hecmw_wtime()

    ! SOLUTION
    phase = 33
    call cluster_sparse_solver ( &
      pt, maxfct, mnum, mtype, phase, nn, aval, irow, jcol, &
      idum, nrhs, iparm, msglvl, rhs, solx, hecmw_comm_get_comm(), istat )
    if (istat < 0) then
      if (myrank==0 .and. spMAT%timelog > 0) call print_iparm_paramters()
      write(*,*) 'ERROR: MKL returned with error in phase 3', istat
      return
    endif

    if( spMAT%type == SPARSE_MATRIX_TYPE_COO ) then !scatter global solution
      call sparse_matrix_scatter_rhs(spMAT, solx)
      do i=1,spMAT%N_loc
        solution(i) = spMAT%rhs(i)
      end do
    endif

    ! solution is written to hecMAT%X. you have to edit lib/solve_LINEQ.f90
    ! to use this routine properly.

    t5=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0) then
      write(*,'(A,f10.3)') ' [Cluster Pardiso]: Solution completed.         time(sec)=',t5-t4
    end if
    if( debug>0 .and. myrank==0 ) call print_iparm_paramters()

#else
    stop "MKL Pardiso not available"
#endif
  end subroutine hecmw_clustermkl_wrapper

#ifdef WITH_MKL
  subroutine print_iparm_paramters()
    write(*,'(A60,I8)') 'Number of iterative refinement steps performed: ',iparm(7)
    write(*,'(A60,I8)') 'Number of perturbed pivots: ',iparm(14)
    write(*,'(A60,I8)') 'Peak memory on symbolic factorization: ',iparm(15)
    write(*,'(A60,I8)') 'Permanent memory on symbolic factorization: ',iparm(16)
    write(*,'(A60,I8)') 'Size of factors/Peak memory on num. fact. and sol: ',iparm(17)
    write(*,'(A60,I8)') 'The number of non-zero elements in the factors: ',iparm(18)
    if( mtype < 11 ) then
      write(*,'(A60,I8)') 'Number of positive eigenvalues: ',iparm(22)
      write(*,'(A60,I8)') 'Number of negative eigenvalues: ',iparm(23)
    end if
  end subroutine

  subroutine export_spMAT_to_CentralizedCRS(spMAT,myrank,n,ia,ja,a,b,x)
    type (sparse_matrix), intent(in)           :: spMAT
    integer(kind=kint), intent(in)             :: myrank
    integer(kind=kint), intent(out)            :: n
    integer(kind=kint), pointer, intent(inout) :: ia(:), ja(:)
    real(kind=kreal), pointer, intent(inout)   :: a(:), b(:), x(:)

    integer(kind=kint) :: i,k
    integer(kind=kint) :: nprocs, ierr, nnz, info
    integer(kind=kint), allocatable :: DISPMAT(:), NCOUNTS(:)
    real(kind=kreal)   :: t1,t2,t3,t4,t5

    t1=hecmw_wtime()

    nprocs = hecmw_comm_get_size()
    allocate(DISPMAT(nprocs),NCOUNTS(nprocs))
    if (nprocs > 1) then
      call HECMW_ALLGATHER_INT_1(spMAT%NZ, NCOUNTS, hecmw_comm_get_comm())
    endif
    n = spMAT%N
    nnz = 0
    DISPMAT(1)=0
    do i=1,nprocs-1
      DISPMAT(i+1)=DISPMAT(i)+NCOUNTS(i)
      nnz = nnz + NCOUNTS(i)
    enddo
    nnz = nnz + NCOUNTS(nprocs)

    if( myrank == 0 ) then
      ierr = 0
      if( .not. associated(ia) ) allocate(ia(nnz+n), stat=ierr)
      if (ierr /= 0) then
        write(*,*) " Allocation error in Setup Cluster MKL, ia",ierr,nnz
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      if( .not. associated(ja) ) allocate(ja(nnz+n), stat=ierr)
      if (ierr /= 0) then
        write(*,*) " Allocation error in Setup Cluster MKL, ja",ierr,nnz
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      if( .not. associated(a) ) allocate(a(nnz+n), stat=ierr)
      if (ierr /= 0) then
        write(*,*) " Allocation error in Setup Cluster MKL, a",ierr,nnz
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      if( .not. associated(b) ) allocate(b(n), stat=ierr)
      if (ierr /= 0) then
        write(*,*) " Allocation error in Setup Cluster MKL, b",ierr,n
        call hecmw_abort(hecmw_comm_get_comm())
      endif
      if( .not. associated(x) ) allocate(x(n), stat=ierr)
      if (ierr /= 0) then
        write(*,*) " Allocation error in Setup Cluster MKL, x",ierr,n
        call hecmw_abort(hecmw_comm_get_comm())
      endif
    else !dummy
      allocate(ia(1),ja(1),a(1),b(1),x(1))
    end if

    t2=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0 .and. debug > 0 ) &
       write(*,'(A,f10.3)') ' [Cluster Pardiso]:   - Allocate Matrix         time(sec)=',t2-t1

    !gather matrix components to rank 0

    call hecmw_gatherv_int(spMAT%IRN, spMAT%NZ, ia, NCOUNTS, DISPMAT, 0, hecmw_comm_get_comm())
    call hecmw_gatherv_int(spMAT%JCN, spMAT%NZ, ja, NCOUNTS, DISPMAT, 0, hecmw_comm_get_comm())
    call hecmw_gatherv_real(spMAT%A, spMAT%NZ, a, NCOUNTS, DISPMAT, 0, hecmw_comm_get_comm())
    if( myrank == 0 ) then !cluster pardiso must be given diagonal components even if they are zero.
      do i=1,n
        ia(nnz+i) = i
        ja(nnz+i) = i
        a(nnz+i) = 0.d0
      end do
      nnz = nnz + n
    endif
    call sparse_matrix_gather_rhs(spMAT, b)

    t3=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0 .and. debug > 0 ) &
       write(*,'(A,f10.3)') ' [Cluster Pardiso]:   - Gather Matrix           time(sec)=',t3-t2

    !convert COO to CRS
    if( myrank==0 ) then
      call coo2csr(n, nnz, ia, ja, a)
      if( debug>0 ) call check_csr(n, nnz, ia, ja, a)
    endif

    t4=hecmw_wtime()
    if (myrank==0 .and. spMAT%timelog > 0 .and. debug > 0 ) &
       write(*,'(A,f10.3)') ' [Cluster Pardiso]:   - Convert Matrix Format   time(sec)=',t4-t3

    deallocate(DISPMAT,NCOUNTS)
  end subroutine

  subroutine coo2csr(n, nnz, ia, ja, a)
    integer(kind=kint), intent(in)             :: n
    integer(kind=kint), intent(inout)          :: nnz
    integer(kind=kint), pointer, intent(inout) :: ia(:)
    integer(kind=kint), pointer, intent(inout) :: ja(:)
    real(kind=kreal), pointer, intent(inout)   :: a(:)

    integer(kind=kint) :: i, k, idx, iS, iE, nnz_new, prev
    integer(kind=kint), allocatable :: counter(:), counter_curr(:), ia0(:), ja0(:)
    real(kind=kreal), allocatable   :: a0(:)

    allocate(counter(0:n),counter_curr(n),ia0(n+1),ja0(nnz),a0(nnz))

    counter(:) = 0
    do i=1,nnz
      idx = ia(i)
      counter(idx) = counter(idx) + 1
    end do
    ia0(:)=0
    do i=1,n
      ia0(i+1) = ia0(i)+counter(i)
    end do

    counter_curr(:) = 0
    do i=1,nnz
      idx = ia(i)
      counter_curr(idx) = counter_curr(idx) + 1
      idx = ia0(idx)+counter_curr(idx)
      ja0(idx) = ja(i)
      a0(idx) = a(i)
    end do

    !$omp parallel private(i, iS, iE)
    !$omp do
    do i=1,n
      iS = ia0(i)+1
      iE = ia0(i+1)
      call sort_column_ascending_order0(iE-iS+1,ja0(iS:iE),a0(iS:iE))
    enddo
    !$omp end do
    !$omp end parallel

    ia(1)=1
    do i=1,n
      ia(i+1) = ia(i)+counter(i)
    end do

    prev = 0
    nnz_new = 0
    ia(1) = 1
    do i=1,n
      do k=ia0(i)+1,ia0(i+1)
        if( ja0(k) == prev ) then
          a(nnz_new) = a(nnz_new) + a0(k)
        else
          nnz_new = nnz_new + 1
          ja(nnz_new) = ja0(k)
          a(nnz_new) = a0(k)
        end if
        prev = ja0(k)
      end do
      ia(i+1) = nnz_new + 1
      prev = 0
    end do

    nnz = nnz_new

    deallocate(counter,counter_curr,ia0,ja0,a0)

  end subroutine

  subroutine check_csr(n, nnz, ia, ja, a)
    integer(kind=kint), intent(in)            :: n
    integer(kind=kint), intent(in)            :: nnz
    integer(kind=kint), pointer, intent(in)   :: ia(:)
    integer(kind=kint), pointer, intent(in)   :: ja(:)
    real(kind=kreal), pointer, intent(in)     :: a(:)

    integer :: i,k

    if( ia(n+1)-1 /= nnz ) then
      write(*,*) "Error in check_csr(1): ia(n+1)-1 /= nnz"
      call hecmw_abort(hecmw_comm_get_comm())
    end if

    do i=1,n
      do k=ia(i),ia(i+1)-1
        if( ja(k) <= 0 ) then
          write(*,*) "Error in check_csr(2): ja(k) <= 0"
          call hecmw_abort(hecmw_comm_get_comm())
        end if
        if( ja(k) > n ) then
          write(*,*) "Error in check_csr(3): ja(k) > n"
          call hecmw_abort(hecmw_comm_get_comm())
        end if
        if( k > ia(i) ) then
          if( ja(k) <= ja(k-1) ) then
            write(*,*) "Error in check_csr(4): ja(k) <= ja(k-1)"
            call hecmw_abort(hecmw_comm_get_comm())
          end if
        end if
      end do
    end do

  end subroutine

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

    integer(kind=kint) :: i, key, work1(2,n), n_low, n_high, next, minidx, maxidx
    logical :: sorted

    if(n<2) return
    !return if work is already sorted
    sorted = .true.
    minidx = work(1,1)
    maxidx = work(1,1)
    do i=1,n-1
      next = work(1,i+1)
      if( work(1,i) > next ) sorted = .false.
      if( next < minidx ) minidx = next
      if( next > maxidx ) maxidx = next
    enddo
    if(sorted) return

    if( n<5 .or. maxidx-minidx < 2) then
      call myinssort(n,work)
      return
    endif

    key = (minidx+maxidx)/2

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
