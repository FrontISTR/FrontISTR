!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
! module for Parallel Direct Solver
module hecmw_solver_direct_parallel

#ifndef HECMW_SERIAL
  use m_child_matrix
  use m_irjc_matrix
  use m_matrix_partition_info
  use m_elap
#endif

  use hecmw_util

#ifndef HECMW_SERIAL
  use hecmw_matrix_ass
  use hecmw_matrix_dump
#endif

  ! access control !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private                            ! default

  public hecmw_solve_direct_parallel ! only entry point of Parallel Direct Solver is public

#ifndef HECMW_SERIAL

  ! internal type definition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type procinfo ! process informations on MPI
    integer(kind=kint) :: myid
    integer(kind=kint) :: imp      ! mother  process id
    logical :: isparent ! true if this process is parent
    logical :: ischild  ! true if this process is child

    integer(kind=kint) :: nchildren ! number of child process
    integer(kind=kint), pointer :: ichildren(:) ! array of children process number

    integer(kind=kint) :: ndiv ! count for matrix division.
  end type procinfo

  type dsinfo ! direct solver information
    integer(kind=kint) :: ndeg   ! dimension of small matrix
    integer(kind=kint) :: neqns  ! number of equations
    integer(kind=kint) :: nstop  ! begining point of C
    integer(kind=kint) :: stage  ! calculation stage
    integer(kind=kint) :: lncol  ! length of col
    integer(kind=kint) :: lndsln ! length of dsln

    integer(kind=kint), pointer :: zpiv(:)    ! in zpivot()
    integer(kind=kint), pointer :: iperm(:)   ! permtation vector
    integer(kind=kint), pointer :: invp(:)    ! inverse permtation of iperm
    integer(kind=kint), pointer :: parent(:)  !
    integer(kind=kint), pointer :: nch(:)     !
    integer(kind=kint), pointer :: xlnzr(:)   ! ia index of whole sparse matrix. (neqns_t + 1)
    integer(kind=kint), pointer :: colno(:)   ! ja index of whole sparse matrix.

    real(kind=kreal), pointer :: diag(:,:)  ! diagonal element
    real(kind=kreal), pointer :: zln(:,:)   ! non diagonal sparse
    real(kind=kreal), pointer :: dsln(:,:)  ! non diagonal dens
  end type dsinfo

  ! internal global variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(procinfo) :: m_pds_procinfo ! initialized in initproc()
  real(kind=kreal), parameter :: rmin = 1.00D-200 ! for inv3() pivot

  integer :: imsg ! output file handler

  integer, parameter :: ilog  = 16 ! according to FSTR
  logical, parameter :: ldbg  = .false.
  integer, parameter :: idbg  = 52 ! according to FSTR
  logical :: lelap = .false.

#endif

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hecmw_solve_direct_parallel(hecMESH, hecMAT, ii)

    implicit none

#ifndef HECMW_SERIAL
    include 'mpif.h'
#endif

    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT
    integer(kind=kint), intent(in) :: ii ! output file handler

#ifndef HECMW_SERIAL

    logical, save :: first_time = .true.

    integer(kind=kint) :: ierr

    ! start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call hecmw_mat_dump(hecMAT, hecMESH)

    imsg=ii ! set message file

    ! set timelog
    if (hecMAT%Iarray(22) .ge. 1) then ! = timelog = kYES (kYES is defined in m_fstr
      lelap = .true.
    end if

    ! set location in process tree
    if (first_time) then
      call initproc()
      hecMAT%Iarray(97) = 0 ! numeric  factorization done flag
      hecMAT%Iarray(98) = 0 ! symbolic factorization done flag
      first_time = .false.
    end if

    if ((hecMAT%Iarray(97) .ne. 0) .or. (hecMAT%Iarray(98) .ne. 0)) then
      write(ilog,*) 'Error: Recalculation of LDU decompose is currently not surported'
      call hecmw_abort( hecmw_comm_get_comm())
    end if

    ! set elap time information
    call initelap(lelap, idbg) !TODO it should be replaced with lelap

    if (m_pds_procinfo%isparent) then
      call elapout('hecmw_solve_direct_parallel: entering sp_direct_parent') !elap
      call sp_direct_parent(hecMESH, hecMAT)
    else if (m_pds_procinfo%ischild) then
      call elapout('hecmw_solve_direct_parallel: entering sp_direct_child') !elap
      call sp_direct_child()
    else
      call elapout('hecmw_solve_direct_parallel: never come here') !elap
      call hecmw_abort( hecmw_comm_get_comm())
    end if

    call MPI_BCAST(hecMAT%b, hecMESH%n_dof*hecMAT%NP, MPI_REAL8, m_pds_procinfo%imp, MPI_COMM_WORLD, ierr)

    call hecmw_mat_dump_solution(hecMAT)

#endif

    return
  end subroutine hecmw_solve_direct_parallel


#ifndef HECMW_SERIAL

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sp_direct_parent(hecMESH, hecMAT)

    implicit none

    include 'mpif.h'

    !I/O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type (hecmwST_local_mesh), intent(inout) :: hecMESH
    type (hecmwST_matrix    ), intent(inout) :: hecMAT


    !internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! given A0 x = b0
    type (irjc_square_matrix) :: a0 ! given left side matrix assembled from sp_matrix
    real(kind=kreal), allocatable :: b(:,:) ! (ndeg,neqns) right hand side value vector of equation.

    ! for divided matrixes
    type(matrix_partition_info) :: pmi
    integer(kind=kint), pointer, save :: iperm_rev(:)
    integer(kind=kint), pointer, save :: iofst_dm(:)

    integer(kind=kint), save               :: neqns_d     ! number of eqns in D matrix
    real(kind=kreal), pointer, save      :: dsln(:,:)   ! non-diagonal elements of dens D matrix
    real(kind=kreal), pointer, save      :: diag(:,:)   ! diagonal elements of dens D matrix
    integer(kind=kint), pointer, save      :: part_all(:)    ! index of corresponding dm of a0 row
    integer(kind=kint), pointer, save      :: iperm_all(:)   ! index in partitioned matrix of a0 row
    type(child_matrix), pointer, save :: dm(:)       !divided matrices

    real(kind=kreal), allocatable :: bd(:,:) ! for right hand side value

    ! internal use
    real(kind=kreal), allocatable :: oldb(:,:)
    logical, save :: nusol_ready = .false.
    integer(kind=kint), save :: ndeg, nndeg, ndegt
    integer(kind=kint), save :: neqns_c, iofst_a2, iofst_c, ndm

    ! misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=kint) :: ierr
    integer(kind=kint) :: i,j,k,l,m,n

    ! for MPI
    integer(kind=kint) :: istatus(MPI_STATUS_SIZE)
    integer(kind=kint) :: ICP
    real(kind=kreal), allocatable :: spdslnval(:,:), diagbuf(:,:), bdbuf(:,:)
    integer(kind=kint), allocatable :: spdslnidx(:)
    integer(kind=kint) :: nspdsln

    ierr=0
    call elapout('sp_direct_parent: entered') !elap


    if (.not. nusol_ready) then

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP01: get a0 from FEM data format hecMESH
      !

      call getA0(hecMESH, hecMAT, a0)
      ndeg=a0%ndeg
      nndeg=ndeg*ndeg
      ndegt = (ndeg+1)*ndeg/2 !triangle element in diag

      call elapout('sp_direct_parent: make a0 done') !elap

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP1X: Parallel LDU decompose of givem A0
      !
      ! STEP11: divide given large matrix A into a1 and a2 ... extend to 2^n division.
      ! C is added to bottom of A1, A2.
      ! D is kept as dsln (off-diagonal), diag (diagonal).

      call elapout('sp_direct_parent: enter matrix partition') !elap
      call matrix_partition_recursive_bisection(a0, m_pds_procinfo%ndiv, pmi)
      call elapout('sp_direct_parent: end matrix partition') !elap
      neqns_d = pmi%neqns_d
      dsln=>pmi%dsln
      diag=>pmi%diag
      part_all=>pmi%part_all
      iperm_all=>pmi%iperm_all
      dm=>pmi%dm
      ndm=pmi%ndm

      ! permtation vector for right hand side value
      allocate(iofst_dm(0:ndm), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      iofst_dm(0)=0
      iofst_dm(1)=neqns_d
      do i=2,ndm
        iofst_dm(i)=iofst_dm(i-1)+dm(i-1)%a%neqns
      end do
      do i=1,a0%neqns
        iperm_all(i)=iperm_all(i)+iofst_dm(part_all(i))
      end do

      allocate(iperm_rev(a0%neqns), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      do i=1, a0%neqns
        iperm_rev(iperm_all(i)) = i
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP12: Send divided left hand side matrixes to child processes.
      !
      ! nstop (separater of A and C) is also send
      ! A and C will be LDU decomposed in child processes.
      ! D region update data will be returned.

      call elapout('sp_direct_parent: send divided matrix to children') !elap
      do i=1,m_pds_procinfo%nchildren
        ICP = m_pds_procinfo%ichildren(i)
        call MPI_SEND(dm(i)%a%ndeg,   1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%a%neqns,  1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%a%nttbr,  1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)

        call MPI_SEND(dm(i)%a%irow, dm(i)%a%nttbr, MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%a%jcol, dm(i)%a%nttbr, MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%a%val,  dm(i)%a%nttbr*dm(i)%a%ndeg*dm(i)%a%ndeg, MPI_REAL8,ICP,1,MPI_COMM_WORLD,istatus,ierr)

        call MPI_SEND(dm(i)%c%ndeg,   1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%c%nttbr,  1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%c%nrows,  1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%c%ncols,  1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)

        call MPI_SEND(dm(i)%c%irow, dm(i)%c%nttbr, MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%c%jcol, dm(i)%c%nttbr, MPI_INTEGER,ICP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(dm(i)%c%val,  dm(i)%c%nttbr*dm(i)%c%ndeg*dm(i)%c%ndeg, MPI_REAL8,ICP,1,MPI_COMM_WORLD,istatus,ierr)
      end do

      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! clean up
      do i=1,m_pds_procinfo%nchildren
        deallocate(dm(i)%a%irow,dm(i)%a%jcol,dm(i)%a%val)
        deallocate(dm(i)%c%irow,dm(i)%c%jcol,dm(i)%c%val)
      end do

      call elapout('sp_direct_parent: end send matrix') !elap

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP13: Receive D region and update D matrix as D' = D - D1' - D2' ...
      !
      ! D is receive as dens matrix, which format is given in s3um2() in serial solver.
      ! to decompose this dens D matrix, use s3um3() on parent.

      call elapout('sp_direct_parent:  receive D matrix element') !elap
      allocate(diagbuf(ndegt, neqns_d), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if

      do k=1,m_pds_procinfo%nchildren
        ICP=m_pds_procinfo%ichildren(k)
        call MPI_RECV(nspdsln,  1,MPI_INTEGER,ICP,1,MPI_COMM_WORLD,istatus,ierr)
        allocate(spdslnidx(nspdsln),spdslnval(nndeg,nspdsln),stat=ierr)
        if(ierr .ne. 0) then
          call errtrp('stop due to allocation error.')
        end if
        call MPI_RECV(spdslnidx,  nspdsln,    MPI_INTEGER,ICP,1,MPI_COMM_WORLD,istatus,ierr)
        call MPI_RECV(spdslnval,  nspdsln*nndeg,MPI_REAL8,ICP,1,MPI_COMM_WORLD,istatus,ierr)
        call MPI_RECV(diagbuf,    neqns_d*ndegt,MPI_REAL8,ICP,1,MPI_COMM_WORLD,istatus,ierr)

        ! off diagonal
        do i=1,nspdsln
          dsln(:,spdslnidx(i)) = dsln(:,spdslnidx(i)) + spdslnval(:,i) ! because of child process dsln is already substructed in s3um2()
        end do

        ! diagonal
        do i=1,neqns_d
          do j=1,ndegt
            diag(j,i) = diag(j,i) + diagbuf(j,i)
          end do
        end do
        deallocate(spdslnidx, spdslnval)
      end do
      deallocate(diagbuf)
      call elapout('sp_direct_parent:  end receive D matrix element') !elap

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP13: ! LDU decompose dens D
      !
      call elapout('sp_direct_parent:  LDU decompose of D. entering nufct0_parent') !elap
      call nufct0_parent(dsln, diag, neqns_d, ndeg)
      call elapout('sp_direct_parent:  exit nufct0_parent') !elap


      nusol_ready = .true.
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP2X: Solve Ax=b0
    !
    ! STEP21: divide right hand side b0 and send it to child in charge.
    !
    ! right hand side vector b is reordered as follows:
    ! b(1:neqns_d)=for parent process
    ! b(neqns_d+1:number of b in dm(1)) for dm(1)
    ! b(end of b1:number of b in dm(2)) for dm(2)
    ! ..
    ! b(end of b_n-1:neqns_t) for dm(ndm)

    ! set right hand side vector (b)
    allocate(b(ndeg,a0%neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    do i=1,a0%neqns
      do j=1,ndeg
        b(j,i)=hecMAT%b(ndeg*(i-1)+j)
      end do
    end do

    ! for verify
    allocate(oldb(ndeg,a0%neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    oldb=b

    call reovec(b, iperm_all)
    do i=1,m_pds_procinfo%nchildren
      ICP=m_pds_procinfo%ichildren(i)
      call MPI_SEND(b(1,iofst_dm(i)+1), dm(i)%ndeg*dm(i)%a%neqns, MPI_REAL8, ICP, 1,MPI_COMM_WORLD, ierr)
    end do

    allocate(bd(ndeg,neqns_d), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    bd(:,1:neqns_d) = b(:,1:neqns_d)

    call elapout('sp_direct_parent:  end send b') !elap

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP22 forward substitution for D region. and get Y
    !
    ! b1, b2... are sended to child processes.
    ! bd is substituted to D in locally, and results from child processes
    ! (C1-Y1 substitute, C2-Y2 substitute...) are receive from child processes.
    ! these value also add for Yd
    call elapout('sp_direct_parent:  begin receive bd') !elap
    allocate(bdbuf(ndeg,neqns_d), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    bdbuf=0
    do k=1,m_pds_procinfo%nchildren
      ICP=m_pds_procinfo%ichildren(k)
      call MPI_RECV(bdbuf, ndeg*neqns_d, MPI_REAL8, ICP, 1,MPI_COMM_WORLD, istatus, ierr)
      do i=1,neqns_d
        do j=1,ndeg
          bd(j,i) = bd(j,i) + bdbuf(j,i)
        end do
      end do
    end do

    call elapout('sp_direct_parent:  end receive bd') !elap

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP22 solve Ax=b for dens matrix, using updated bd
    !
    call elapout('sp_direct_parent:  begin solve Ax_d=b_d') !elap
    call nusol0_parent(dsln, diag, bd, neqns_d, ndeg)
    call elapout('sp_direct_parent:  end solve Ax_d=b_d') !elap


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP23 send Xd to children
    !
    call elapout('sp_direct_parent:  begin send Xd') !elap
    call MPI_BCAST(bd, ndeg*neqns_d, MPI_REAL8, m_pds_procinfo%imp, MPI_COMM_WORLD, ierr)
    call elapout('sp_direct_parent:  end send Xd') !elap

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP24 Receive results Xi from children.
    !
    call elapout('sp_direct_parent:  begin receive X') !elap
    do k=1,m_pds_procinfo%nchildren
      ICP=m_pds_procinfo%ichildren(k)
      call MPI_RECV(b(1,iofst_dm(k)+1), dm(k)%ndeg*dm(k)%a%neqns, MPI_REAL8, ICP, 1,MPI_COMM_WORLD, istatus, ierr)
    end do
    b(:,1:neqns_d)=bd(:,:) ! set xd
    call elapout('sp_direct_parent:  end receive X') !elap

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP25 restore final result
    !
    ! permtation divided matrix => initial order
    call elapout('sp_direct_parent:  begin permutate X') !elap
    call reovec(b, iperm_rev)
    call elapout('sp_direct_parent:  end permutate X') !elap

    ! verify result
    call verif0(a0%neqns, ndeg, a0%nttbr, a0%irow, a0%jcol, a0%val, oldb, b) !verify result oldb will be broken.

    ! set result to FEM data
    do i=1,a0%neqns
      do j=1,ndeg
        hecMAT%b(ndeg*(i-1)+j)=b(j,i)
      end do
    end do

    call elapout('sp_direct_parent: end solve Ax=b') !elap
    deallocate(b, bd, bdbuf, oldb)
    return
  end subroutine sp_direct_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine getA0(hecMESH,hecMAT, a0)

    implicit none

    type (hecmwST_local_mesh), intent(in)  :: hecMESH
    type (hecmwST_matrix    ), intent(in)  :: hecMAT
    type (irjc_square_matrix), intent(out) :: a0

    integer(kind=kint) :: i,j,k,l,ierr,numnp,ndof,kk,ntotal
    integer(kind=kint) :: iiS,iiE,kki,kkj,ndof2

    NUMNP = hecMAT%NP
    NDOF  = hecMESH%n_dof
    ntotal = NUMNP*NDOF

    !*NUFACT variables
    a0%neqns = NUMNP
    a0%ndeg  = NDOF
    a0%nttbr = hecMAT%NP+hecMAT%NPL !+hecMAT%NPU if unsymmetric

    !*Allocations
    allocate(a0%irow(a0%nttbr),stat=ierr)
    allocate(a0%jcol(a0%nttbr),stat=ierr)
    allocate(a0%val(a0%ndeg*a0%ndeg, a0%nttbr),stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    kk = 0
    ndof2 = NDOF*NDOF
    do j= 1, NUMNP
      !*Diagonal
      kk = kk + 1
      a0%irow(kk) = j
      a0%jcol(kk) = j
      call vlcpy(a0%val(:,kk),hecMAT%D(ndof2*(j-1)+1:ndof2*j),ndof)
      !*Lower
      do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
        i= hecMAT%itemL(k)
        kk = kk + 1
        a0%irow(kk) = j
        a0%jcol(kk) = i
        call vlcpy(a0%val(:,kk),hecMAT%AL(ndof2*(k-1)+1:ndof2*k),ndof)
      enddo
    enddo

    return

  contains

    subroutine vlcpy(a,b,n)
      implicit none
      real(kind=kreal),   intent(out) :: a(:)
      real(kind=kreal),   intent(in)  :: b(:)
      integer(kind=kint), intent(in)  :: n

      integer(kind=kint) :: i,j

      do i = 1,n
        do j = 1,n
          a((j-1)*n+i) = b((i-1)*n+j) !transpose
        end do
      end do
      return
    end subroutine vlcpy

  end subroutine getA0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sp_direct_child()
    ! parallel direct solver child process

    implicit none

    include 'mpif.h'

    type (child_matrix)  :: cm

    real(kind=kreal), allocatable :: b(:,:) ! (ndeg, neqns)
    type (dsinfo) :: dsi
    logical, save :: nusol_ready = .false.

    ! for MPI
    integer(kind=kint) :: istatus(MPI_STATUS_SIZE)
    integer(kind=kint) :: imp, ierr
    integer(kind=kint) :: i,j,k,l,m,n

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call elapout('sp_direct_child: entered') !elap

    imp=m_pds_procinfo%imp

    if (.not. nusol_ready) then

      call elapout('sp_direct_child: waiting matrix from parent via MPI') !elap
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP01: get a,c from parent
      ! C matrix is placed below A.
      !
      call MPI_RECV(cm%a%ndeg,   1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%a%neqns,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%a%nttbr,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      allocate(cm%a%irow(cm%a%nttbr), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      allocate(cm%a%jcol(cm%a%nttbr), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      allocate(cm%a%val(cm%a%ndeg*cm%a%ndeg, cm%a%nttbr), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      call MPI_RECV(cm%a%irow, cm%a%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%a%jcol, cm%a%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%a%val,  cm%a%nttbr*cm%a%ndeg*cm%a%ndeg, MPI_REAL8,IMP,1,MPI_COMM_WORLD,istatus,ierr)


      call MPI_RECV(cm%c%ndeg,   1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%c%nttbr,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%c%nrows,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%c%ncols,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      allocate(cm%c%irow(cm%c%nttbr), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      allocate(cm%c%jcol(cm%c%nttbr), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      allocate(cm%c%val(cm%c%ndeg*cm%c%ndeg, cm%c%nttbr), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      call MPI_RECV(cm%c%irow, cm%c%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%c%jcol, cm%c%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(cm%c%val,  cm%c%nttbr*cm%c%ndeg*cm%c%ndeg, MPI_REAL8,IMP,1,MPI_COMM_WORLD,istatus,ierr)

      cm%ndeg    = cm%a%ndeg
      cm%ista_c  = cm%a%neqns+1
      cm%neqns_t = cm%a%neqns + cm%c%nrows
      call elapout('sp_direct_child: end get matrix from parent via MPI') !elap
      !call MPI_BARRIER(MPI_COMM_WORLD, ierr)



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! STEP1x: LDU decompose for given matrix
      !

      ! set up dsi for allocate matrix array for fill-in
      call elapout('sp_direct_child: entering matini_para') !elap
      call matini_para(cm, dsi, ierr)
      call elapout('sp_direct_child: exit matini_para') !elap

      call elapout('sp_direct_child: entering staij1') !elap
      ! set real8 value
      do i=1,cm%a%nttbr
        call staij1(0, cm%a%irow(i), cm%a%jcol(i), cm%a%val(:,i), dsi, ierr)
      end do
      do i=1,cm%c%nttbr
        !  call staij1(0, cm%c%irow(i)+cm%a%neqns, dsi%iperm(cm%c%jcol(i)), cm%c%val(:,i), dsi, ierr)
        call staij1(0, cm%c%irow(i)+cm%a%neqns, cm%c%jcol(i), cm%c%val(:,i), dsi, ierr)
      end do
      call elapout('sp_direct_child: end staij1') !elap

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! following STEP12-15 will be done in nufct0_child()
      !
      ! STEP12: LDU decompose of A (1..nstop-1)
      ! STEP13: LDU decompose of C (nstop..neqnsA+neqnsd)
      ! STEP14: update D region.
      ! STEP15: send D region to parent
      call elapout('sp_direct_child: entering nufct0_child') !elap
      call nufct0_child(dsi, ierr)
      call elapout('sp_direct_child: exit nufct0_child') !elap

      nusol_ready = .true.
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP2x: solve ax=b by using forward an backword substitution
    !
    ! STEP21: receive b from parent
    !
    allocate(b(cm%ndeg, cm%neqns_t), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    ! wait for right hand side vector
    call MPI_RECV(b, cm%ndeg*cm%a%neqns, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, istatus, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! following STEP22-24 is done in nusol0-child()
    !
    ! STEP22: forward substitution for A
    ! STEP23: forward substitution for C and send it (yi) to parent
    ! STEP24: divide with diagonal matrix
    ! STEP25: receive xd from parent and do backword substitution
    call elapout('sp_direct_child: enter nusol0_child') !elap
    call nusol0_child(b, dsi, ierr)
    call elapout('sp_direct_child: exit nusol0_child') !elap

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP26: send final result to parent
    !
    call elapout('sp_direct_child: begin send result to parent') !elap
    call MPI_SEND(b, cm%ndeg*cm%a%neqns, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)
    call elapout('sp_direct_child: end send result to parent') !elap

    return

  end subroutine sp_direct_child


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initproc()
    implicit none

    include 'mpif.h'

    integer(kind=kint) :: npe, myid
    integer(kind=kint) :: ndiv
    integer(kind=kint) :: ierr
    integer(kind=kint) :: i,j,k,l

    m_pds_procinfo%isparent=.false.
    m_pds_procinfo%ischild=.false.

    ! get process number and number of whole processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, npe, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

    m_pds_procinfo%myid = myid

    ndiv=0
    do
      if (2**(ndiv + 1) .gt. npe) then
        exit
      end if
      ndiv = ndiv + 1
    end do
    m_pds_procinfo%ndiv = ndiv

    if (npe .ne. 2**ndiv + 1) then
      write(ilog,*) 'Error: please use 2**n+1 (3,5,9,17...) processes for parallel direct solver.'
      write(6,*)    'Error: please use 2**n+1 (3,5,9,17...) processes for parallel direct solver.'
      call hecmw_abort( hecmw_comm_get_comm())
      stop
    end if

    if (myid.eq.0) then
      write(idbg,*)'parent process.'
      m_pds_procinfo%isparent=.true.
      m_pds_procinfo%nchildren=2**ndiv
      allocate(m_pds_procinfo%ichildren(m_pds_procinfo%nchildren))
      do i=1, 2**ndiv
        m_pds_procinfo%ichildren(i)=i
      end do
    else
      write(idbg,*)'child process.'
      m_pds_procinfo%ischild=.true.
      m_pds_procinfo%imp=0
    end if

    return
  end subroutine initproc

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine errtrp(mes)
    character(*) mes
    write(6,*)  'Error in : process ', m_pds_procinfo%myid
    write(ilog,*) mes

    call hecmw_abort( hecmw_comm_get_comm())
    stop
  end subroutine errtrp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine matini_para(cm,dsi,ir)

    !----------------------------------------------------------------------
    !
    !     matini initializes storage for sparse matrix solver.
    !     this routine is used for both symmetric matrices
    !     and must be called once at the beginning
    !
    !    (i)
    !        neqns     number of unknowns
    !        nttbr     number of non0s, pattern of non-zero elements are
    !                  given like following.
    !                  nonz(A)={(i,j);i=irow(l),j=jcol(l); 1<= l <= nttbr}
    !        irow
    !        jcol      to define non-zero pattern
    !        lenv      length of the array v (iv)
    !
    !    (o)
    !        dsi       matrix informations
    !        ir        return code
    !                              =0    normal
    !                              =-1   non positive index
    !                              =1    too big index
    !                              =10   insufficient storage
    !
    !
    !        stage   10  after initialization
    !                20  building up matrix
    !                30  after LU decomposition
    !                40  after solving
    !
    !         # coded by t.arakawa
    !         # reviced by t.kitayama of Univ. Tokyo on 20071120
    !
    !----------------------------------------------------------------------

    implicit none

    type(child_matrix), intent(in)  :: cm
    type(dsinfo),       intent(out) :: dsi
    integer(kind=kint), intent(out) :: ir

    integer(kind=kint), pointer :: irow_a(:), jcol_a(:)
    integer(kind=kint), pointer :: irow_c(:), jcol_c(:)

    integer(kind=kint), pointer :: ia(:)     ! in stiaja() neqns+2
    integer(kind=kint), pointer :: ja(:)     ! in stiaja() 2*nttbr
    integer(kind=kint), pointer :: jcpt(:)   ! in stsmat() 2*nttbr
    integer(kind=kint), pointer :: jcolno(:) ! in stsmat() 2*nttbr

    integer(kind=kint), pointer :: iperm_a(:)
    integer(kind=kint), pointer :: invp_a(:)

    integer(kind=kint), pointer :: xlnzr_a(:)
    integer(kind=kint), pointer :: colno_a(:)

    integer(kind=kint), pointer :: xlnzr_c(:)
    integer(kind=kint), pointer :: colno_c(:)


    integer(kind=kint), pointer :: adjncy(:) ! in genqmd() 2*nttbr
    integer(kind=kint), pointer :: qlink(:)  ! in genqmd() neqne+2
    integer(kind=kint), pointer :: qsize(:)  ! in genqmd() neqne+2
    integer(kind=kint), pointer :: nbrhd(:)  ! in genqmd() neqne+2
    integer(kind=kint), pointer :: rchset(:) ! in genqmd() neqne+2

    integer(kind=kint), pointer :: cstr(:)

    integer(kind=kint), pointer :: adjt(:)   ! in rotate() neqne+2
    integer(kind=kint), pointer :: anc(:)    ! in rotate() neqne+2

    integer(kind=kint), pointer :: lwk3arr(:)
    integer(kind=kint), pointer :: lwk2arr(:)
    integer(kind=kint), pointer :: lwk1arr(:)
    integer(kind=kint), pointer :: lbtreearr(:,:) ! genbtq() (2,neqns+1)
    integer(kind=kint), pointer :: lleafarr(:)
    integer(kind=kint), pointer :: lxleafarr(:)
    integer(kind=kint), pointer :: ladparr(:)
    integer(kind=kint), pointer :: lpordrarr(:)

    integer(kind=kint) :: neqns_a, nttbr_a, neqns_a1, nstop, neqns_t, neqns_d, nttbr_c, ndeg
    integer(kind=kint) :: lncol_a, lncol_c
    integer(kind=kint) :: neqnsz, nofsub, izz, izz0, lnleaf ! dummy variables
    integer(kind=kint) :: ir1
    integer(kind=kint) :: i, j, k , ipass, ks, ke, ierr

    ndeg    = cm%ndeg

    neqns_t =  cm%neqns_t
    neqns_d =  cm%c%nrows

    neqns_a =  cm%a%neqns
    nttbr_a =  cm%a%nttbr
    irow_a  => cm%a%irow
    jcol_a  => cm%a%jcol

    nttbr_c =  cm%c%nttbr
    irow_c  => cm%c%irow
    jcol_c  => cm%c%jcol

    dsi%neqns=neqns_t ! because direct solver treat A + C as one matrix.
    dsi%ndeg=ndeg

    neqns_a1=neqns_a+2
    ir=0
    ierr=0
    izz0 = 0.0d0
    !
    !  set z pivot
    !
    allocate(dsi%zpiv(neqns_a), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call zpivot(neqns_a,neqnsz,nttbr_a,jcol_a,irow_a,dsi%zpiv,ir1)
    if(ir1.ne.0) then
      ir=ir1
      goto 1000
    endif
    !
    !  build jcpt,jcolno
    !
    allocate(jcpt(2*nttbr_a), jcolno(2*nttbr_a), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call stsmat(neqns_a,nttbr_a,irow_a,jcol_a,jcpt,jcolno)
    !
    !  build ia,ja
    !
    allocate(ia(neqns_a1), ja(2*nttbr_a), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call stiaja(neqns_a, neqns_a,ia,ja,jcpt,jcolno)
    !
    !  get permutation vector iperm,invp
    !

    ! setup identity permtation for C matrix
    allocate(iperm_a(neqns_a), invp_a(neqns_a), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call idntty(neqns_a,invp_a,iperm_a)

    ! reorder A matrix
    allocate(adjncy(2*nttbr_a),qlink(neqns_a1),qsize(neqns_a1),nbrhd(neqns_a1),rchset(neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    allocate(lwk2arr(neqns_a1),lwk1arr(neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call genqmd(neqns_a,ia,ja,iperm_a,invp_a,lwk1arr,lwk2arr,rchset,nbrhd,qsize,qlink,nofsub,adjncy)
    deallocate(adjncy, qlink, qsize, nbrhd, rchset)

    ! set ia, ja
    call stiaja(neqns_a, neqns_a, ia, ja, jcpt,jcolno)

    !   build up the parent vector parent vector will be saved in
    !   work2 for a while
    allocate(cstr(neqns_a1),adjt(neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    10 continue
    call genpaq(ia,ja,invp_a,iperm_a,lwk2arr,neqns_a,cstr)

    !   build up the binary tree
    allocate (lbtreearr(2,neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call genbtq(ia, ja, invp_a, iperm_a,lwk2arr,lbtreearr,dsi%zpiv,izz,neqns_a)

    !   rotate the binary tree to avoid a zero pivot
    if(izz.eq.0) goto 20
    if(izz0.eq.0) izz0=izz
    if(izz0.ne.izz) goto 30
    call rotate(ia, ja, invp_a, iperm_a, lwk2arr,lbtreearr,izz,neqns_a,anc,adjt,ir1)
    goto 10
    30 continue
    call bringu(dsi%zpiv,iperm_a, invp_a, lwk2arr,izz,neqns_a,ir1)
    goto 10

    !   post ordering
    20 continue
    allocate(lwk3arr(0:neqns_a1),lpordrarr(neqns_a1),dsi%parent(neqns_a1), dsi%nch(neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call posord(dsi%parent,lbtreearr,invp_a,iperm_a,lpordrarr,dsi%nch,neqns_a,lwk1arr,lwk2arr,lwk3arr)

    !   generate skelton graph
    allocate(lleafarr(nttbr_a),lxleafarr(neqns_a1),ladparr(neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call gnleaf(ia, ja, invp_a, iperm_a, lpordrarr,dsi%nch,ladparr,lxleafarr,lleafarr,neqns_a,lnleaf)

    !   build up xlnzr,colno  (this is the symbolic fct.)
    nstop = cm%ista_c
    call countclno(dsi%parent, lxleafarr, lleafarr, neqns_a, nstop, lncol_a, ir1) ! only for A
    allocate(colno_a(lncol_a),xlnzr_a(neqns_a1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    call gnclno(dsi%parent,lpordrarr,lxleafarr,lleafarr,xlnzr_a, colno_a, neqns_a, nstop,lncol_a,ir1) ! only for A
    if(ir1.ne.0) then
      ir=10
      goto 1000
    endif

    ! do symbolic LDU decomposition for C region.
    call lduDecomposeC(xlnzr_a,colno_a,invp_a,iperm_a, ndeg, nttbr_c, irow_c, &
      jcol_c, cm%c%ncols, cm%c%nrows, xlnzr_c, colno_c, lncol_c)

    ! set calculated informations to dsi.
    allocate(dsi%xlnzr(neqns_t + 1), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    dsi%xlnzr(1:neqns_a)=xlnzr_a(:)
    dsi%xlnzr(neqns_a+1:neqns_t+1)=xlnzr_c(:)+xlnzr_a(neqns_a+1)-1

    dsi%lncol=lncol_a + lncol_c
    allocate(dsi%colno(lncol_a + lncol_c), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    dsi%colno(1:lncol_a)=colno_a(:)
    dsi%colno(lncol_a+1:lncol_a+lncol_c)=colno_c(:)

    allocate(dsi%invp(neqns_t), dsi%iperm(neqns_t), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    dsi%invp(1:neqns_a)=invp_a(1:neqns_a)
    dsi%iperm(1:neqns_a)=iperm_a(1:neqns_a)
    do i=neqns_a+1,neqns_t
      dsi%invp(i)=i
      dsi%iperm(i)=i
    end do

    deallocate(xlnzr_a, colno_a, xlnzr_c, colno_c, invp_a, iperm_a)

    dsi%nstop=nstop
    dsi%stage=10
    1000 continue

    return
  end subroutine matini_para

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufct0_child(dsi,ir)

    implicit none
    type(dsinfo),       intent(inout) :: dsi
    integer(kind=kint), intent(out)   :: ir
    !
    !     this performs Cholesky factorization
    !

    if(dsi%stage.ne.20) then
      ir=40
      goto 1000
    else
      ir=0
    endif
    if(dsi%ndeg.eq.1) then
      call nufct1_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%neqns,dsi%parent,dsi%nch,dsi%nstop,ir)
    else if(dsi%ndeg.eq.2) then
      call nufct2_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%neqns,dsi%parent,dsi%nch,dsi%nstop,ir)
    else if(dsi%ndeg.eq.3) then
      call nufct3_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%neqns,dsi%parent,dsi%nch,dsi%nstop,ir)
      !      else if(dsi%ndeg.eq.6) then !TODO implement it
      !        call nufct6_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%neqns,dsi%parent,dsi%nch,dsi%nstop,ir)
    else
      call nufctx_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%neqns,dsi%parent,dsi%nch,dsi%nstop,dsi%ndeg,ir)
    end if

    dsi%stage=30
    1000 continue
    return
  end subroutine nufct0_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufct1_child(xlnzr,colno,zln,diag,neqns,parent,nch,nstop,ir)

    implicit none
    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),parent(:)
    integer(kind=kint), intent(in)  :: neqns, nstop, ir
    integer(kind=kint), intent(out) :: nch(:)
    real(kind=kreal),   intent(out) :: zln(:,:),diag(:,:) !zln(1,:), diag(1,:)

    integer(kind=kint) :: neqns_c
    integer(kind=kint) :: i,j,k,l, ic,ierr,imp
    integer(kind=kint)          :: nspdsln
    integer(kind=kint), pointer :: spdslnidx(:)
    real(kind=kreal),   pointer :: spdslnval(:,:)
    include 'mpif.h'

    !----------------------------------------------------------------------
    !
    !     nufct1 performs cholesky factorization in row order for ndeg=1
    !
    !     (i) xlnzr,colno,zln,diag
    !         symbolicaly factorized
    !
    !     (o) zln,diag,dsln
    !
    !         #coded by t.arakawa
    !
    !----------------------------------------------------------------------

    !
    ! phase I
    ! LDU decompose of A (1..nstop-1)
    !
    diag(1,1)=1.0d0/diag(1,1)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    do 100 ic=2,nstop-1
      call sum(ic,xlnzr,colno,zln(1,:),diag(1,:),nch,parent,neqns)
      100 continue
      !
      ! phase II
      ! LDU decompose of C (nstop..neqnsA+neqnsd)
      !
      do 200 ic=nstop,neqns
        call sum1(ic,xlnzr,colno,zln(1,:),diag(1,:),parent,neqns)
        200 continue
        !
        ! phase III
        ! Update D region.
        !

        ! clear dummy diagonal value for D region
        do i=nstop,neqns
          diag(:,i)=0.0
        end do

        neqns_c = neqns - nstop + 1
        call sum2_child(neqns,nstop,xlnzr,colno,zln(1,:),diag(1,:),spdslnidx,spdslnval,nspdsln)
        ! send D region to parent
        imp = m_pds_procinfo%imp
        call MPI_SEND(nspdsln, 1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(spdslnidx,  nspdsln,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(spdslnval,  nspdsln,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(diag(1,nstop),  neqns_c,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
        return
  end subroutine nufct1_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufct2_child(xlnzr, colno, zln, diag, neqns, parent, nch, nstop, ir)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:), colno(:), parent(:)
    integer(kind=kint), intent(out) :: nch(:)
    real(kind=kreal),   intent(out) :: zln(:,:), diag(:,:) !zln(6,*), diag(3,*)
    integer(kind=kint), intent(in)  :: neqns, nstop
    integer(kind=kint), intent(out) :: ir

    integer(kind=kint) :: i,j,k,l, ic, imp, ierr, neqns_c
    integer(kind=kint) :: nspdsln
    integer(kind=kint), pointer :: spdslnidx(:)
    real(kind=kreal),   pointer :: spdslnval(:,:)
    include 'mpif.h'

    !----------------------------------------------------------------------
    !
    !     nufct2 performs cholesky factorization in row order for ndeg=2
    !
    !     (i) xlnzr,colno,zln,diag
    !         symbolicaly factorized
    !
    !     (o) zln,diag,dsln
    !
    !         #coded by t.arakawa
    !
    !----------------------------------------------------------------------


    ! For parallel calculation, factorization for A and C region
    ! will be done.
    !
    ! Creation for D region also done.
    !
    ! Factorization for D region is omitted.

    !
    ! phase I
    ! LDU decompose of A (1..nstop-1)
    !
    call elapout('nufct2_child: begin phase I LDU decompose of A') !elap
    if(nstop.gt.1) call inv2(diag(:,1),ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    do ic=2,nstop-1
      call s2um(ic,xlnzr,colno,zln,diag,nch,parent,neqns)
    end do
    !
    ! phase II
    ! LDU decompose of C (nstop..neqnsA+neqnsd)
    !
    call elapout('nufct2_child: begin phase II LDU decompose of C') !elap
    do ic=nstop,neqns
      call s2um1(ic,xlnzr,colno,zln,diag,nch,parent,neqns)
    end do

    !
    ! phase III
    ! Update D region.
    !

    call elapout('nufct2_child: begin phase III update D region') !elap

    ! clear dummy diagonal value for D region ! currently dummy value setting was not done
    do i=nstop,neqns
      diag(:,i)=0.0
    end do

    neqns_c = neqns - nstop + 1
    call s2um2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln)

    ! send D region to parent
    imp = m_pds_procinfo%imp
    call MPI_SEND(nspdsln, 1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spdslnidx,  nspdsln,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spdslnval,  nspdsln*4,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(diag(1,nstop),  neqns_c*3,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)


    return

  end subroutine nufct2_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufct3_child(xlnzr,colno,zln,diag,neqns,parent,nch,nstop,ir)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),parent(:)
    integer(kind=kint), intent(out) :: nch(:)
    real(kind=kreal),   intent(out) :: zln(:,:),diag(:,:) !zln(9,*),diag(6,*)
    integer(kind=kint), intent(in)  :: neqns, nstop, ir

    integer(kind=kint) :: i,j,k,l, ic, imp, ierr, neqns_c
    integer(kind=kint)          :: nspdsln
    integer(kind=kint), pointer :: spdslnidx(:)
    real(kind=kreal), pointer :: spdslnval(:,:)
    include 'mpif.h'
    !
    !----------------------------------------------------------------------
    !
    !     nufct3 performs cholesky factorization in row order for ndeg=3
    !
    !     (i) xlnzr,colno,zln,diag
    !         symbolicaly factorized
    !
    !     (o) zln,diag,dsln
    !
    !         #coded by t.arakawa
    !
    !----------------------------------------------------------------------

    !
    ! For parallel calculation, factorization for A and C region
    ! will be done.
    !
    ! Creation for D region also done.
    !
    ! Factorization for D region is omitted.
    !

    !
    ! phase I
    ! LDU decompose of A (1..nstop-1)
    !
    call elapout('nufct3_child: begin phase I LDU decompose of A') !elap
    if(nstop.gt.1) call inv3(diag(:,1),ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    do 100 ic=2,nstop-1
      call s3um(ic,xlnzr,colno,zln,diag,nch,parent,neqns)
      100 continue
      !
      ! phase II
      ! LDU decompose of C (nstop..neqnsA+neqnsd)
      !
      call elapout('nufct3_child: begin phase II LDU decompose of C') !elap
      do 200 ic=nstop,neqns
        call s3um1(ic,xlnzr,colno,zln,diag,nch,parent,neqns)
        200 continue

        !
        ! phase III
        ! Update D region.
        !

        call elapout('nufct3_child: begin phase III update D region') !elap

        ! clear dummy diagonal value for D region ! currently dummy value setting was not done
        do i=nstop,neqns
          diag(:,i)=0.0
        end do

        neqns_c = neqns - nstop + 1
        call s3um2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln)

        ! send D region to parent
        imp = m_pds_procinfo%imp
        call MPI_SEND(nspdsln, 1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(spdslnidx,  nspdsln,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(spdslnval,  nspdsln*9,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(diag(1,nstop),  neqns_c*6,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)

        return
  end subroutine nufct3_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !subroutine nufct6_child(xlnzr, colno, zln, diag, neqns, parent, nch, nstop, ir)
  !TODO caution! currently this routine is not implemented because of s6um1() is not implemented.
  !
  !implicit none
  !
  !integer :: xlnzr(:), colno(:), parent(:), nch(:)
  !real(8) :: zln(:,:), diag(:,:) !zln(36,*), diag(21,*)
  !
  !integer :: neqns, nstop, ir, imp, ierr
  !
  !integer :: neqns_c
  !integer :: l, ic
  !integer :: i
  !
  !real(8), allocatable :: dsln(:,:)
  !include 'mpif.h'
  !
  !!----------------------------------------------------------------------
  !!
  !!     nufct6 performs cholesky factorization in row order for ndeg=3
  !!
  !!     (i) xlnzr,colno,zln,diag
  !!         symbolicaly factorized
  !!
  !!     (o) zln,diag,dsln
  !!
  !!         #coded by t.arakawa
  !!
  !!----------------------------------------------------------------------
  !
  !
  !! For parallel calculation, factorization for A and C region
  !! will be done.
  !!
  !! Creation for D region also done.
  !!
  !! Factorization for D region is omitted.
  !!
  !
  !!
  !! phase I
  !! LDU decompose of A (1..nstop-1)
  !!
  !call elapout('nufct6_child: begin phase I LDU decompose of A') !elap
  !if(nstop.gt.1) call inv6(diag(:,1),ir)
  !l=parent(1)
  !nch(l)=nch(l)-1
  !nch(1)=-1
  !do ic=2,nstop-1
  !  call s6um(ic,xlnzr,colno,zln,diag,nch,parent,neqns)
  !end do
  !
  !!
  !! phase II
  !! LDU decompose of C (nstop..neqnsA+neqnsd)
  !!
  !call elapout('nufct6_child: begin phase II LDU decompose of C') !elap
  !do ic=nstop,neqns
  !  call s6um1(ic,xlnzr,colno,zln,diag,nch,parent,neqns)
  !end do
  !end subroutine nufct6_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufctx_child(xlnzr,colno,zln,diag,neqns,parent,nch,nstop,ndeg,ir)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),parent(:)
    integer(kind=kint), intent(out) :: nch(:)
    real(kind=kreal),   intent(out) :: zln(:,:),diag(:,:) !zln(9,*),diag(6,*)
    integer(kind=kint), intent(in)  :: neqns, nstop, ndeg
    integer(kind=kint), intent(out) :: ir

    integer(kind=kint) :: i,j,k,l, ic, neqns_c, ndeg2, ndegl, imp, ierr
    integer(kind=kint)          :: nspdsln
    integer(kind=kint), pointer :: spdslnidx(:)
    real(kind=kreal),   pointer :: spdslnval(:,:)
    include 'mpif.h'
    !----------------------------------------------------------------------
    !
    !     nufctx performs cholesky factorization in row order for every ndeg
    !
    !     (i) xlnzr,colno,zln,diag
    !         symbolicaly factorized
    !
    !     (o) zln,diag,dsln
    !
    !         #coded by t.arakawa
    !
    !----------------------------------------------------------------------

    !
    ! For parallel calculation, factorization for A and C region
    ! will be done.
    !
    ! Creation for D region also done.
    !
    ! Factorization for D region is omitted.
    !

    ndeg2=ndeg*ndeg
    ndegl=(ndeg+1)*ndeg/2
    !
    ! phase I
    ! LDU decompose of A (1..nstop-1)
    if(nstop.gt.1) call invx(diag,ndeg,ir)
    l=parent(1)
    nch(l)=nch(l)-1
    nch(1)=-1
    do ic=2,nstop-1
      call sxum(ic,xlnzr,colno,zln,diag,nch,parent,neqns,ndeg,ndegl)
    end do

    !
    ! phase II
    ! LDU decompose of C (nstop..neqnsA+neqnsd)
    !
    do ic=nstop,neqns
      call sxum1(ic,xlnzr,colno,zln,diag,nch,parent,neqns,ndeg,ndegl)
    end do

    !
    ! phase III
    ! Update D region.
    !
    neqns_c = neqns - nstop + 1
    call sxum2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln,ndeg,ndegl)

    ! send D region to parent
    imp = m_pds_procinfo%imp
    call MPI_SEND(nspdsln, 1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spdslnidx,  nspdsln,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(spdslnval,  nspdsln*ndeg2,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
    call MPI_SEND(diag(1,nstop),  neqns_c*ndegl,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
    return
  end subroutine nufctx_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol0_child(b,dsi,ir)

    !----------------------------------------------------------------------
    !
    !     this performs forward elimination and backward substitution
    !
    !     (i/o)
    !           b        on entry     right hand side vector
    !                    on exit      solution vector
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------

    implicit none

    real(kind=kreal),   intent(inout) :: b(:,:)
    type(dsinfo),       intent(inout) :: dsi
    integer(kind=kint), intent(out)   :: ir

    integer(kind=kint) :: neqns, nstop, ndeg

    if(dsi%stage.ne.30 .and. dsi%stage.ne.40) then
      ir=50
      goto 1000
    else
      ir=0
    end if
    if(dsi%ndeg.eq.1) then
      call nusol1_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%iperm,b,dsi%neqns,dsi%nstop)
    else if(dsi%ndeg .eq. 2) then
      call nusol2_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%iperm,b,dsi%neqns,dsi%nstop)
    else if(dsi%ndeg .eq. 3) then
      call nusol3_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%iperm,b,dsi%neqns,dsi%nstop)!org
      !      else if(dsi%ndeg .eq. 6) then !TODO implement it
      !        call nusol6_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%iperm,b,dsi%neqns,dsi%nstop)
    else
      call nusolx_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%iperm,b,dsi%neqns,dsi%nstop,dsi%ndeg)
    endif
    dsi%stage=40
    1000 continue
    return
  end subroutine nusol0_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol1_child(xlnzr,colno,zln,diag,iperm,b,neqns,nstop)

    implicit none

    integer(kind=kint), intent(in)    :: xlnzr(:), colno(:), iperm(:)
    real(kind=kreal),   intent(in)    :: zln(:,:), diag(:,:)
    real(kind=kreal),   intent(inout) :: b(:,:)
    integer(kind=kint), intent(in)    :: neqns, nstop

    integer(kind=kint) :: neqns_a, neqns_c
    integer(kind=kint) :: k, ks, ke, i, j, imp, ierr
    real(kind=kreal), allocatable :: wk(:), wk_d(:)

    include 'mpif.h'

    ! forward

    !     now nstop is begining point of C
    neqns_a = nstop - 1
    neqns_c = neqns - nstop + 1

    allocate(wk(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk = 0

    do 10 i=1,neqns_a
      wk(i)=b(1,iperm(i))
      10 continue

      ! STEP22: forward substitution for A
      do 100 i=1,neqns_a
        ks=xlnzr(i)
        ke=xlnzr(i+1)-1
        if(ke.lt.ks) goto 110
        wk(i)=wk(i)-spdot2(wk,zln(1,:),colno,ks,ke)
        110    continue
        100 continue

        ! STEP23: forward substitution for C and send it (yi) to parent
        allocate(wk_d(nstop:neqns), stat=ierr)
        if(ierr .ne. 0) then
          call errtrp('stop due to allocation error.')
        end if
        wk_d=0

        do 101 i=nstop,neqns
          ks=xlnzr(i)
          ke=xlnzr(i+1)-1
          if(ke.lt.ks) goto 111
          wk_d(i)=wk_d(i)-spdot2(wk,zln(1,:),colno,ks,ke)
          111   continue
          101 continue
          imp = m_pds_procinfo%imp
          call MPI_SEND(wk_d, neqns_c, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)

          ! STEP24: divide with diagonal matrix
          do 120 i=1,neqns
            wk(i)=wk(i)*diag(1,i)
            120 continue

            ! STEP25: receive xd from parent and do backword substitution
            call MPI_BCAST(wk_d, neqns_c, MPI_REAL8, IMP, MPI_COMM_WORLD, ierr)

            wk(nstop:neqns)=wk_d(nstop:neqns)

            do 200 i=neqns,1,-1
              ks=xlnzr(i)
              ke=xlnzr(i+1)-1
              if(ke.lt.ks) goto 200

              do 210 k=ks,ke
                j=colno(k)
                wk(j)=wk(j)-wk(i)*zln(1,k)
                210    continue
                200 continue

                ! permutaion
                do 300 i=1,neqns
                  b(1,iperm(i))=wk(i)
                  300 continue
                  return

  end subroutine nusol1_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol2_child(xlnzr, colno, zln, diag, iperm, b, neqns, nstop)
    ! perform forward substitution for sparse matrix,
    ! send bd to parent and receive xd from parent,
    ! backword substitution using xd ad send final result x to parent.

    use m_elap

    implicit none

    integer(kind=kint), intent(in)    :: xlnzr(:), colno(:), iperm(:)
    real(kind=kreal),   intent(in)    :: zln(:,:), diag(:,:) !zln(4,*), diag(3,*), b(2,*)
    real(kind=kreal),   intent(inout) :: b(:,:) ! b(2,*)

    real(kind=kreal), allocatable :: wk(:,:), wk_d(:,:)
    integer(kind=kint) :: neqns_c, neqns_a, nstop, neqns, ks, ke
    integer(kind=kint) :: i, j, k, l, imp, ierr

    include 'mpif.h'

    !now nstop is begining point of C
    neqns_a = nstop - 1
    neqns_c = neqns - nstop + 1

    allocate(wk(2,neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk = 0
    do i=1,neqns_a
      wk(1,i) = b(1,iperm(i))
      wk(2,i) = b(2,iperm(i))
    end do

    ! STEP22: forward substitution for A
    call elapout('nusol2_child: begin forward substitution for A') !elap
    do i=1, neqns_a
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      if(ke.ge.ks) then
        call s2pdot(wk(:,i),wk,zln,colno,ks,ke)
      end if
    end do

    ! STEP23: forward substitution for C and send it (yi) to parent
    call elapout('nusol2_child: begin forward substitution for C') !elap
    allocate(wk_d(2,nstop:neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk_d=0

    do i=nstop,neqns
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      if(ke.ge.ks) then
        call s2pdot(wk_d(:,i),wk,zln,colno,ks,ke)
      end if
    end do

    call elapout('nusol2_child: wait to send wk_d') !elap
    imp = m_pds_procinfo%imp
    call MPI_SEND(wk_d, 2*neqns_c, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)

    ! STEP24: divide with diagonal matrix
    do i=1,neqns_a
      wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
      wk(1,i)=wk(1,i)*diag(1,i)
      wk(2,i)=wk(2,i)*diag(3,i)
      wk(1,i)=wk(1,i)-wk(2,i)*diag(2,i)
    end do

    ! STEP25: receive xd from parent and do backword substitution
    call elapout('nusol2_child: wait until receive wk_d') !elap
    call MPI_BCAST(wk_d, 2*neqns_c, MPI_REAL8, IMP, MPI_COMM_WORLD, ierr)
    call elapout('nusol2_child: end receive wk_d') !elap
    call elapout('nusol2_child: begin backword substitution') !elap

    wk(:,nstop:neqns)=wk_d(:,nstop:neqns)
    do i=neqns,1,-1
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      if(ke.ge.ks) then
        do k=ks,ke
          j=colno(k)
          wk(1,j)=wk(1,j)-wk(1,i)*zln(1,k)-wk(2,i)*zln(2,k)
          wk(2,j)=wk(2,j)-wk(1,i)*zln(3,k)-wk(2,i)*zln(4,k)
        end do
      end if
    end do
    call elapout('nusol2_child: end backword substitution') !elap

    ! permutaion
    do i=1,neqns_a
      b(1,iperm(i))=wk(1,i)
      b(2,iperm(i))=wk(2,i)
    end do

    call elapout('nusol2_child: end') !elap
    return

  end subroutine nusol2_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol3_child(xlnzr,colno,zln,diag,iperm,b,neqns, nstop)

    ! perform forward substitution for sparse matrix,
    ! send bd to parent and receive xd from parent,
    ! backword substitution using xd ad send final result x to parent.

    use m_elap

    implicit none

    integer(kind=kint), intent(in)    :: xlnzr(:), colno(:), iperm(:)
    real(kind=kreal),   intent(in)    :: zln(:,:), diag(:,:) !zln(9,*),diag(6,*)
    real(kind=kreal),   intent(inout) :: b(:,:) ! b(3,*)

    real(kind=kreal), allocatable :: wk(:,:), wk_d(:,:)
    integer(kind=kint) :: neqns_c, neqns_a, nstop, neqns, ks, ke
    integer(kind=kint) :: i, j, k, l, imp, ierr

    include 'mpif.h'

    !     now nstop is begining point of C

    neqns_a = nstop - 1
    neqns_c = neqns - nstop + 1

    call elapout('nusol3_child: entered') !elap

    allocate(wk(3,neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk = 0
    do 10 i=1,neqns_a
      wk(1,i)=b(1,iperm(i))
      wk(2,i)=b(2,iperm(i))
      wk(3,i)=b(3,iperm(i))
      10 continue


      ! STEP22: forward substitution for A
      call elapout('nusol3_child: begin forward substitution for A') !elap

      do 100 i=1, neqns_a
        ks=xlnzr(i)
        ke=xlnzr(i+1)-1
        if(ke.lt.ks) goto 110
        call s3pdot(wk(:,i),wk,zln,colno,ks,ke)
        110   continue
        100 continue


        ! STEP23: forward substitution for C and send it (yi) to parent
        call elapout('nusol3_child: begin forward substitution for C') !elap
        allocate(wk_d(3,nstop:neqns), stat=ierr)
        if(ierr .ne. 0) then
          call errtrp('stop due to allocation error.')
        end if
        wk_d=0

        do 101 i=nstop,neqns
          ks=xlnzr(i)
          ke=xlnzr(i+1)-1
          if(ke.lt.ks) goto 111
          call s3pdot(wk_d(:,i),wk,zln,colno,ks,ke)
          111   continue
          101 continue

          call elapout('nusol3_child: wait to send wk_d') !elap
          imp = m_pds_procinfo%imp
          call MPI_SEND(wk_d, 3*neqns_c, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)

          ! STEP24: divide with diagonal matrix
          call elapout('nusol3_child: divide with diagonal matrix') !elap
          do 120 i=1,neqns_a
            wk(2,i)=wk(2,i)-wk(1,i)*diag(2,i)
            wk(3,i)=wk(3,i)-wk(1,i)*diag(4,i)-wk(2,i)*diag(5,i)
            wk(1,i)=wk(1,i)*diag(1,i)
            wk(2,i)=wk(2,i)*diag(3,i)
            wk(3,i)=wk(3,i)*diag(6,i)
            wk(2,i)=wk(2,i)-wk(3,i)*diag(5,i)
            wk(1,i)=wk(1,i)-wk(2,i)*diag(2,i)-wk(3,i)*diag(4,i)
            120 continue

            ! STEP25: receive xd from parent and do backword substitution
            call elapout('nusol3_child: wait until receive wk_d') !elap
            call MPI_BCAST(wk_d, 3*neqns_c, MPI_REAL8, IMP, MPI_COMM_WORLD, ierr)
            call elapout('nusol3_child: end receive wk_d') !elap
            call elapout('nusol3_child: begin backword substitution') !elap

            wk(:,nstop:neqns)=wk_d(:,nstop:neqns)

            do 200 i=neqns,1,-1
              ks=xlnzr(i)
              ke=xlnzr(i+1)-1
              if(ke.lt.ks) goto 200

              do 210 k=ks,ke
                j=colno(k)

                wk(1,j)=wk(1,j)-wk(1,i)*zln(1,k)-wk(2,i)*zln(2,k)-wk(3,i)*zln(3,k)
                wk(2,j)=wk(2,j)-wk(1,i)*zln(4,k)-wk(2,i)*zln(5,k)-wk(3,i)*zln(6,k)
                wk(3,j)=wk(3,j)-wk(1,i)*zln(7,k)-wk(2,i)*zln(8,k)-wk(3,i)*zln(9,k)
                210    continue
                200 continue
                call elapout('nusol3_child: end backword substitution') !elap


                ! permutaion
                do 300 i=1,neqns_a
                  b(1,iperm(i))=wk(1,i)
                  b(2,iperm(i))=wk(2,i)
                  b(3,iperm(i))=wk(3,i)
                  300 continue

                  call elapout('nusol3_child: end') !elap
                  return
  end subroutine nusol3_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nusolx_child(xlnzr, colno, zln, diag, iperm, b, neqns, nstop, ndeg)

    implicit none

    integer(kind=kint), intent(in)    :: xlnzr(:), colno(:), iperm(:)
    real(kind=kreal),   intent(in)    :: zln(:,:), diag(:,:)
    real(kind=kreal),   intent(inout) :: b(:,:)
    integer(kind=kint), intent(in)    :: neqns, nstop, ndeg

    real(kind=kreal), allocatable :: wk(:,:), wk_d(:,:)
    integer(kind=kint) :: neqns_c, neqns_a, ks, ke, locd, loc1
    integer(kind=kint) :: i, j, k, l, m, n, imp, ierr

    include 'mpif.h'

    !now nstop is begining point of C
    neqns_a = nstop - 1
    neqns_c = neqns - nstop + 1

    allocate(wk(ndeg,neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk = 0
    do i=1,neqns_a
      wk(1,i)=b(1,iperm(i))
      wk(2,i)=b(2,iperm(i))
      wk(3,i)=b(3,iperm(i))
    end do

    ! STEP22: forward substitution for A
    do i=1,neqns_a
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      if(ke.ge.ks) then
        call sxpdot(ndeg,wk(1,i),wk,zln,colno,ks,ke)
      end if
    end do

    ! STEP23: forward substitution for C and send it (yi) to parent
    allocate(wk_d(ndeg,nstop:neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk_d=0
    do i=nstop,neqns
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      if(ke.ge.ks) then
        call sxpdot(ndeg,wk_d(:,i),wk,zln,colno,ks,ke)
      end if
    end do

    imp = m_pds_procinfo%imp
    call MPI_SEND(wk_d, ndeg*neqns_c, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)


    ! STEP24: divide with diagonal matrix
    do i=1,neqns_a
      locd=0
      do m=1,ndeg-1
        locd=locd+m
        loc1=locd+m
        do n=m+1,ndeg
          wk(n,i)=wk(n,i)-wk(m,i)*diag(loc1,i)
          loc1=loc1+n
        end do
      end do

      locd=0
      do m=1,ndeg
        locd=locd+m
        wk(m,i)=wk(m,i)*diag(locd,i)
      end do

      do n=ndeg,2,-1
        locd=locd-1
        do m=n-1,1,-1
          wk(m,i)=wk(m,i)-wk(n,i)*diag(locd,i)
          locd=locd-1
        end do
      end do
    end do

    call MPI_BCAST(wk_d, ndeg*neqns_c, MPI_REAL8, IMP, MPI_COMM_WORLD, ierr)
    wk(:,nstop:neqns)=wk_d(:,nstop:neqns)

    ! back ward
    do i=neqns,1,-1
      ks=xlnzr(i)
      ke=xlnzr(i+1)-1
      if(ke.ge.ks) then
        do k=ks,ke
          j=colno(k)
          do m=1,ndeg
            do n=1,ndeg
              wk(m,j)=wk(m,j)-wk(n,i)*zln(n+(m-1)*ndeg,k)
            end do
          end do
        end do
      end if
    end do

    ! permutaion
    do l=1,ndeg
      do i=1,neqns_a
        b(l,iperm(i))=wk(l,i)
      end do
    end do

    return
  end subroutine nusolx_child
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufct0_parent(dsln, diag, neqns, ndeg)
    ! select LDU decomposer for dens matrix according to ndeg

    implicit none

    real(kind=kreal),   intent(inout) :: dsln(:,:)
    real(kind=kreal),   intent(inout) :: diag(:,:)
    integer(kind=kint), intent(in)    :: neqns, ndeg

    integer(kind=kint) :: ndegl

    if (ndeg .eq. 1) then
      call sum3(neqns, dsln(1,:), diag(1,:))
    else if (ndeg .eq. 3) then
      call s3um3(neqns, dsln, diag)
    else
      ndegl = (ndeg+1)*ndeg/2
      call sxum3(neqns, dsln, diag, ndeg, ndegl)
    end if

    return
  end subroutine nufct0_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol0_parent(dsln, diag, b, neqns, ndeg)
    ! select solvers according to ndeg

    implicit none

    real(kind=kreal),   intent(in)    :: dsln(:,:)
    real(kind=kreal),   intent(in)    :: diag(:,:)
    real(kind=kreal),   intent(inout) :: b(:,:)

    integer(kind=kint), intent(in)    :: neqns, ndeg

    if (ndeg .eq. 1) then
      call nusol1_parent(dsln(1,:), diag(1,:), b(1,:), neqns)
    else if (ndeg .eq. 3) then
      call nusol3_parent(dsln, diag, b, neqns)
    else
      call nusolx_parent(dsln, diag, b, neqns, ndeg)
    end if

    return
  end subroutine nusol0_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol1_parent(dsln, diag, b, neqns)
    ! solve Ax=b for dens matrix with ndeg=1
    ! require dsln, diag is already LDU decomposed.
    ! currently not tested 20071121

    implicit none

    real(kind=kreal),   intent(in)    :: dsln(:) !((neqns+1)*neqns/2)
    real(kind=kreal),   intent(in)    :: diag(:) !(neqns)
    real(kind=kreal),   intent(inout) :: b(:)    !(3,neqns)
    integer(kind=kint), intent(in)    :: neqns

    integer(kind=kint) :: i,j,k,l,loc

    ! forward substitution
    do i=2,neqns
      k=(i-1)*(i-2)/2 + 1 ! first element of i'th row.
      b(i)=b(i)-dot_product(b(1:i-1),dsln(k:k+i-2))
    end do

    ! divide by D (because of diag is already inverted (1/Dii))
    b(:)=b(:)*diag(:)

    ! Backword substitution.
    ! Substitute Zi into D and get Xd results.
    loc=(neqns-1)*neqns/2
    do i=neqns,1,-1
      do j=i-1,1,-1
        b(j)=b(j)-b(i)*dsln(loc)
        loc=loc-1
      end do
    end do

    return
  end subroutine nusol1_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol2_parent(dsln, diag, b, neqns)
    ! solve Ax=b for dens matrix with ndeg=3
    ! require dsln, diag is already LDU decomposed.

    implicit none

    real(kind=kreal),   intent(in)    :: dsln(:,:) !(4, (neqns+1)*newns/2)
    real(kind=kreal),   intent(in)    :: diag(:,:) !(3,neqns)
    real(kind=kreal),   intent(inout) :: b(:,:) !(2,neqns)
    integer(kind=kint), intent(in)    :: neqns

    integer(kind=kint) :: i,j,k,l,loc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP22 forward substitution
    !
    do i=2,neqns
      k=(i-1)*(i-2)/2 + 1 ! first element of i'th row.
      call d2sdot(b(:,i),b,dsln(:, k:k+i-2),i-1)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP23 divide Yd by diagonal elment of D and get Zi=Yi/Di
    !
    do i=1,neqns
      b(2,i)=b(2,i)-b(1,i)*diag(2,i)
      b(1,i)=b(1,i)*diag(1,i)
      b(2,i)=b(2,i)*diag(3,i)
      b(1,i)=b(1,i)-b(2,i)*diag(2,i)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP24 Backword substitution.
    ! Substitute Zi into D and get Xd results.
    !
    loc=(neqns-1)*neqns/2
    do i=neqns,1,-1
      do j=i-1,1,-1
        b(1,j)=b(1,j)-b(1,i)*dsln(1,loc)-b(2,i)*dsln(2,loc)
        b(2,j)=b(2,j)-b(1,i)*dsln(3,loc)-b(2,i)*dsln(4,loc)
        loc=loc-1
      end do
    end do

    return

  end subroutine nusol2_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusol3_parent(dsln, diag, b, neqns)
    ! solve Ax=b for dens matrix with ndeg=3
    ! require dsln, diag is already LDU decomposed.

    implicit none

    real(kind=kreal),   intent(in)    :: dsln(:,:) !(9,(neqns+1)*neqns/2)
    real(kind=kreal),   intent(in)    :: diag(:,:) !(6,neqns)
    real(kind=kreal),   intent(inout) :: b(:,:)    !(3,neqns)
    integer(kind=kint), intent(in)    :: neqns

    integer(kind=kint) :: i,j,k,l,loc

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP22 forward substitution
    !
    do i=2,neqns
      k=(i-1)*(i-2)/2 + 1 ! first element of i'th row.
      call d3sdot(b(:,i),b,dsln(:, k:k+i-2),i-1)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP23 divide Yd by diagonal elment of D and get Zi=Yi/Di
    !
    do i=1,neqns
      b(2,i)=b(2,i)-b(1,i)*diag(2,i)
      b(3,i)=b(3,i)-b(1,i)*diag(4,i)-b(2,i)*diag(5,i)
      b(1,i)=b(1,i)*diag(1,i)
      b(2,i)=b(2,i)*diag(3,i)
      b(3,i)=b(3,i)*diag(6,i)
      b(2,i)=b(2,i)-b(3,i)*diag(5,i)
      b(1,i)=b(1,i)-b(2,i)*diag(2,i)-b(3,i)*diag(4,i)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP24 Backword substitution.
    ! Substitute Zi into D and get Xd results.
    !
    loc=(neqns-1)*neqns/2
    do i=neqns,1,-1
      do j=i-1,1,-1
        b(1,j)=b(1,j)-b(1,i)*dsln(1,loc)-b(2,i)*dsln(2,loc)-b(3,i)*dsln(3,loc)
        b(2,j)=b(2,j)-b(1,i)*dsln(4,loc)-b(2,i)*dsln(5,loc)-b(3,i)*dsln(6,loc)
        b(3,j)=b(3,j)-b(1,i)*dsln(7,loc)-b(2,i)*dsln(8,loc)-b(3,i)*dsln(9,loc)
        loc=loc-1
      end do
    end do

    return
  end subroutine nusol3_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nusolx_parent (dsln,diag,b,neqns,ndeg)

    implicit none

    real(kind=kreal),   intent(in)    :: diag(:,:), dsln(:,:)
    real(kind=kreal),   intent(inout) :: b(:,:)
    integer(kind=kint), intent(in)    :: neqns, ndeg

    integer(kind=kint) :: i,j,k,l,m,n,loc, locd, loc1

    ! forward
    do i=2,neqns
      k=(i-1)*(i-2)/2 + 1 ! first element of i'th row.
      call dxsdot(ndeg,b(:,i),b,dsln(:,k:k+i-2),i-1)
    end do

    ! divide
    do i=1,neqns
      locd=0
      do m=1,ndeg-1
        locd=locd+m
        loc1=locd+m
        do n=m+1,ndeg
          b(n,i)=b(n,i)-b(m,i)*diag(loc1,i)
          loc1=loc1+n
        end do
      end do

      locd=0
      do m=1,ndeg
        locd=locd+m
        b(m,i)=b(m,i)*diag(locd,i)
      end do

      do n=ndeg,2,-1
        locd=locd-1
        do m=n-1,1,-1
          b(m,i)=b(m,i)-b(n,i)*diag(locd,i)
          locd=locd-1
        end do
      end do
    end do
    ! back ward
    loc=(neqns-1)*neqns/2
    do i=neqns,1,-1
      do j=i-1,1,-1
        do m=1,ndeg
          do n=1,ndeg
            b(m,j)=b(m,j)-b(n,i)*dsln((m-1)*ndeg+n,loc)
          end do
        end do
        loc=loc-1
      end do
    end do
    return
  end subroutine nusolx_parent

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s3um2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln)

    implicit none

    integer(kind=kint), intent(in)    :: neqns, nstop
    integer(kind=kint), intent(in)    :: xlnzr(:),colno(:)
    real(kind=kreal),   intent(inout) :: zln(:,:),diag(:,:) !zln(9,*),diag(6,*),dsln(9,*) !
    integer(kind=kint), pointer       :: spdslnidx(:)
    real(kind=kreal),   pointer       :: spdslnval(:,:)
    integer(kind=kint), intent(out)   :: nspdsln

    real(kind=kreal),   allocatable :: temp(:,:)
    integer(kind=kint), allocatable :: indx(:)
    logical :: ftflag
    integer(kind=kint) :: i,j,k,l,m,n, ic,ks,ke,ii,jj,jc,j1,j2,loc,ispdsln, ierr

    allocate(temp(neqns,9),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    nspdsln=0
    do ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do k=ks,ke
        jj=colno(k)
        indx(jj)=ic
      end do
      do jc=nstop,ic-1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        do jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            nspdsln=nspdsln+1
            exit
          endif
        end do
      end do
    end do
    allocate(spdslnidx(nspdsln), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    allocate(spdslnval(9,nspdsln), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    loc=0
    ispdsln=0
    spdslnval=0
    ftflag = .true.
    do 100 ic=nstop,neqns
      !        do 105 m=1,9
      !           do 105 jj=1,nstop
      !              temp(jj,m)=0.0d0
      ! 105    continue
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do 110 k=ks,ke
        jj=colno(k)
        temp(jj,1)=zln(1,k)
        temp(jj,2)=zln(2,k)
        temp(jj,3)=zln(3,k)
        temp(jj,4)=zln(4,k)
        temp(jj,5)=zln(5,k)
        temp(jj,6)=zln(6,k)
        temp(jj,7)=zln(7,k)
        temp(jj,8)=zln(8,k)
        temp(jj,9)=zln(9,k)
        indx(jj)=ic
        110    continue
        do 111 k=ks,ke
          jj=colno(k)
          zln(4,k)=temp(jj,4)-temp(jj,1)*diag(2,jj)
          zln(7,k)=temp(jj,7)-temp(jj,1)*diag(4,jj)-zln(4,k)*diag(5,jj)
          zln(1,k)=temp(jj,1)*diag(1,jj)
          zln(4,k)=zln(4,k)*diag(3,jj)
          zln(7,k)=zln(7,k)*diag(6,jj)
          zln(4,k)=zln(4,k)-zln(7,k)*diag(5,jj)
          zln(1,k)=zln(1,k)-zln(4,k)*diag(2,jj)-zln(7,k)*diag(4,jj)
          !
          zln(5,k)=temp(jj,5)-temp(jj,2)*diag(2,jj)
          zln(8,k)=temp(jj,8)-temp(jj,2)*diag(4,jj)-zln(5,k)*diag(5,jj)
          zln(2,k)=temp(jj,2)*diag(1,jj)
          zln(5,k)=zln(5,k)*diag(3,jj)
          zln(8,k)=zln(8,k)*diag(6,jj)
          zln(5,k)=zln(5,k)-zln(8,k)*diag(5,jj)
          zln(2,k)=zln(2,k)-zln(5,k)*diag(2,jj)-zln(8,k)*diag(4,jj)
          !
          zln(6,k)=temp(jj,6)-temp(jj,3)*diag(2,jj)
          zln(9,k)=temp(jj,9)-temp(jj,3)*diag(4,jj)-zln(6,k)*diag(5,jj)
          zln(3,k)=temp(jj,3)*diag(1,jj)
          zln(6,k)=zln(6,k)*diag(3,jj)
          zln(9,k)=zln(9,k)*diag(6,jj)
          zln(6,k)=zln(6,k)-zln(9,k)*diag(5,jj)
          zln(3,k)=zln(3,k)-zln(6,k)*diag(2,jj)-zln(9,k)*diag(4,jj)
          !        write(60,6000) k,(zln(llll,k),llll=1,9)
          !6000    format(i6,3d20.10/6x,3d20.10/6x,3d20.10)
          111    continue
          !
          do 112 k=ks,ke
            jj=colno(k)
            diag(1,ic)=diag(1,ic)-temp(jj,1)*zln(1,k)-temp(jj,4)*zln(4,k)-temp(jj,7)*zln(7,k)
            diag(2,ic)=diag(2,ic)-temp(jj,1)*zln(2,k)-temp(jj,4)*zln(5,k)-temp(jj,7)*zln(8,k)
            diag(3,ic)=diag(3,ic)-temp(jj,2)*zln(2,k)-temp(jj,5)*zln(5,k)-temp(jj,8)*zln(8,k)
            diag(4,ic)=diag(4,ic)-temp(jj,1)*zln(3,k)-temp(jj,4)*zln(6,k)-temp(jj,7)*zln(9,k)
            diag(5,ic)=diag(5,ic)-temp(jj,2)*zln(3,k)-temp(jj,5)*zln(6,k)-temp(jj,8)*zln(9,k)
            diag(6,ic)=diag(6,ic)-temp(jj,3)*zln(3,k)-temp(jj,6)*zln(6,k)-temp(jj,9)*zln(9,k)
            112    continue
            do 120 jc=nstop,ic-1
              loc=loc+1
              j1=xlnzr(jc)
              j2=xlnzr(jc+1)
              do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
                j=colno(jj)
                if(indx(j).eq.ic) then
                  if (ftflag) then
                    ispdsln=ispdsln+1
                    ftflag=.false.
                  end if
                  spdslnidx(ispdsln)=loc
                  spdslnval(1,ispdsln)=spdslnval(1,ispdsln)-temp(j,1)*zln(1,jj)-temp(j,4)*zln(4,jj)-temp(j,7)*zln(7,jj)
                  spdslnval(2,ispdsln)=spdslnval(2,ispdsln)-temp(j,2)*zln(1,jj)-temp(j,5)*zln(4,jj)-temp(j,8)*zln(7,jj)
                  spdslnval(3,ispdsln)=spdslnval(3,ispdsln)-temp(j,3)*zln(1,jj)-temp(j,6)*zln(4,jj)-temp(j,9)*zln(7,jj)
                  spdslnval(4,ispdsln)=spdslnval(4,ispdsln)-temp(j,1)*zln(2,jj)-temp(j,4)*zln(5,jj)-temp(j,7)*zln(8,jj)
                  spdslnval(5,ispdsln)=spdslnval(5,ispdsln)-temp(j,2)*zln(2,jj)-temp(j,5)*zln(5,jj)-temp(j,8)*zln(8,jj)
                  spdslnval(6,ispdsln)=spdslnval(6,ispdsln)-temp(j,3)*zln(2,jj)-temp(j,6)*zln(5,jj)-temp(j,9)*zln(8,jj)
                  spdslnval(7,ispdsln)=spdslnval(7,ispdsln)-temp(j,1)*zln(3,jj)-temp(j,4)*zln(6,jj)-temp(j,7)*zln(9,jj)
                  spdslnval(8,ispdsln)=spdslnval(8,ispdsln)-temp(j,2)*zln(3,jj)-temp(j,5)*zln(6,jj)-temp(j,8)*zln(9,jj)
                  spdslnval(9,ispdsln)=spdslnval(9,ispdsln)-temp(j,3)*zln(3,jj)-temp(j,6)*zln(6,jj)-temp(j,9)*zln(9,jj)
                endif
                220       continue
                ftflag = .true.
                120    continue
                100 continue
                return
  end subroutine s3um2_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine zpivot(neqns,neqnsz,nttbr,jcol,irow,zpiv,ir)

    implicit none

    integer(kind=kint), intent(in)  :: jcol(:),irow(:)
    integer(kind=kint), intent(out) :: zpiv(:)
    integer(kind=kint), intent(in)  :: neqns,nttbr
    integer(kind=kint), intent(out) :: neqnsz,ir

    integer(kind=kint) :: i,j,k,l

    ir=0
    do 100 l=1,neqns
      zpiv(l)=1
      100 continue

      do 200 l=1,nttbr
        i=irow(l)
        j=jcol(l)
        if(i.le.0.or.j.le.0) then
          ir=-1
          goto 1000
        elseif(i.gt.neqns.or.j.gt.neqns) then
          ir=1
          goto 1000
        endif
        if(i.eq.j) zpiv(i)=0
        200 continue

        do 310 i=neqns,1,-1
          if(zpiv(i).eq.0) then
            neqnsz=i
            goto 320
          endif
          310 continue
          320 continue
          1000 continue
          if(ldbg) write(idbg,*) '# zpivot ########################'
          if(ldbg) write(idbg,60) (zpiv(i),i=1,neqns)
          60 format(20i3)
          return
  end subroutine zpivot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine stsmat(neqns,nttbr,irow,jcol,jcpt,jcolno)

    implicit none

    integer(kind=kint), intent(in)  :: irow(:), jcol(:)
    integer(kind=kint), intent(out) :: jcpt(:), jcolno(:)
    integer(kind=kint), intent(in)  :: neqns, nttbr

    integer(kind=kint) :: i,j,k,l,loc,locr

    do 10 i=1,2*nttbr
      jcpt(i)=0
      jcolno(i)=0
      10 continue
      do 20 i=1,neqns
        jcpt(i)=i+neqns
        jcolno(i+neqns)=i
        20 continue

        k=2*neqns
        do 100 l=1,nttbr
          i=irow(l)
          j=jcol(l)
          if(i.eq.j) goto 100
          loc=jcpt(i)
          locr=i
          110    continue
          if(loc.eq.0) goto 120
          if(jcolno(loc).eq.j) then
            goto 100
          elseif(jcolno(loc).gt.j) then
            goto 130
          endif
          locr=loc
          loc=jcpt(loc)
          goto 110
          120    continue
          k=k+1
          jcpt(locr)=k
          jcolno(k)=j
          goto 150
          130    continue
          k=k+1
          jcpt(locr)=k
          jcpt(k)=loc
          jcolno(k)=j
          150    continue
          loc=jcpt(j)
          locr=j
          160    continue
          if(loc.eq.0) goto 170
          if(jcolno(loc).eq.i) then
            goto 100
          elseif(jcolno(loc).gt.i) then
            goto 180
          endif
          locr=loc
          loc=jcpt(loc)
          goto 160
          170    continue
          k=k+1
          jcpt(locr)=k
          jcolno(k)=i
          goto 100
          180    continue
          k=k+1
          jcpt(locr)=k
          jcpt(k)=loc
          jcolno(k)=i
          100 continue
          if(ldbg) then
            write(idbg,*) 'jcolno'
            write(idbg,60) (jcolno(i),i=1,k)
            write(idbg,*) 'jcpt'
            write(idbg,60) (jcpt(i),i=1,k)
            60 format(10i7)
          endif
          return
  end subroutine stsmat

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine stiaja(neqns,neqnsz,ia,ja,jcpt,jcolno)

    implicit none
    !
    !     coded by t.arakawa
    !
    integer(kind=kint), intent(in)  :: jcpt(:),jcolno(:)
    integer(kind=kint), intent(out) :: ia(:),ja(:)
    integer(kind=kint), intent(in)  :: neqns, neqnsz

    integer(kind=kint) :: i,j,k,l,ii,loc
    !

    ia(1)=1
    l=0
    do 100 k=1,neqns
      loc=jcpt(k)
      110    continue
      if(loc.eq.0) goto 120
      ii=jcolno(loc)
      if(ii.eq.k.or.ii.gt.neqnsz) goto 130
      l=l+1
      ja(l)=ii
      130    continue
      loc=jcpt(loc)
      goto 110
      120    ia(k+1)=l+1
      100 continue
      if(ldbg) then
        write(idbg,*) 'stiaja(): ia '
        write(idbg,60) (ia(i),i=1,neqns+1)
        write(idbg,*) 'stiaja(): ja '
        write(idbg,60) (ja(i),i=1,ia(neqns+1))
      endif
      60 format(10i7)
      return
  end subroutine stiaja

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine idntty(neqns,invp,iperm)

    implicit none

    integer(kind=kint), intent(out) :: invp(:),iperm(:)
    integer(kind=kint), intent(in)  :: neqns

    integer(kind=kint) :: i

    do 100 i=1,neqns
      invp(i)=i
      iperm(i)=i
      100 continue
      return
  end subroutine idntty

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine genqmd(neqns,xadj,adj0,perm,invp,deg,marker,rchset,nbrhd,qsize,qlink,nofsub,adjncy)

    implicit none

    integer(kind=kint), intent(in)  :: adj0(:),xadj(:)
    integer(kind=kint), intent(out) :: rchset(:),nbrhd(:),adjncy(:),perm(:),invp(:),deg(:),marker(:),qsize(:),qlink(:)
    integer(kind=kint), intent(in)  :: neqns
    integer(kind=kint), intent(out) :: nofsub

    integer(kind=kint) :: inode,ip,irch,mindeg,nhdsze,node,np,num,nump1,nxnode,rchsze,search,thresh,ndeg
    integer(kind=kint) :: i,j,k,l

    mindeg=neqns
    nofsub=0
    do 10 i=1,xadj(neqns+1)-1
      adjncy(i)=adj0(i)
      10 continue
      do 100 node=1,neqns
        perm(node)=node
        invp(node)=node
        marker(node)=0
        qsize(node)=1
        qlink(node)=0
        ndeg=xadj(node+1)-xadj(node)
        deg(node)=ndeg
        if(ndeg.lt.mindeg) mindeg=ndeg
        100 continue

        num=0
        200 search=1
        thresh=mindeg
        mindeg=neqns
        300 nump1=num+1
        if(nump1.gt.search) search=nump1
        do 400 j=search,neqns
          node=perm(j)
          if(marker(node).lt.0) goto 400
          ndeg=deg(node)
          if(ndeg.le.thresh) goto 500
          if(ndeg.lt.mindeg) mindeg=ndeg
          400 continue
          goto 200

          500 search=j
          nofsub=nofsub+deg(node)
          marker(node)=1
          call qmdrch(node,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze,nbrhd)
          nxnode=node
          600 num=num+1
          np=invp(nxnode)
          ip=perm(num)
          perm(np)=ip
          invp(ip)=np
          perm(num)=nxnode
          invp(nxnode)=num
          deg(nxnode)=-1
          nxnode=qlink(nxnode)
          if(nxnode.gt.0) goto 600
          if(rchsze.le.0) goto 800
          !
          call qmdupd(xadj,adjncy,rchsze,rchset,deg,qsize,qlink,marker,rchset(rchsze+1:),nbrhd(nhdsze+1:))
          marker(node)=0
          do 700 irch=1,rchsze
            inode=rchset(irch)
            if(marker(inode).lt.0) goto 700
            marker(inode)=0
            ndeg=deg(inode)
            if(ndeg.lt.mindeg) mindeg=ndeg
            if(ndeg.gt.thresh) goto 700
            mindeg=thresh
            thresh=ndeg
            search=invp(inode)
            700 continue
            if(nhdsze.gt.0) call qmdot(node,xadj,adjncy,marker,rchsze,rchset,nbrhd)
            800 if(num.lt.neqns) goto 300
            return
  end subroutine genqmd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine genpaq(xadj,adjncy,invp,iperm,parent,neqns,ancstr)

    implicit none

    integer(kind=kint), intent(in)  :: xadj(:),adjncy(:),invp(:),iperm(:)
    integer(kind=kint), intent(out) :: parent(:),ancstr(:)
    integer(kind=kint), intent(in)  :: neqns

    integer(kind=kint) :: i,j,k,l,ip,it

    do 100 i=1,neqns
      parent(i)=0
      ancstr(i)=0
      ip=iperm(i)
      do 110 k=xadj(ip),xadj(ip+1)-1
        l=invp(adjncy(k))
        if(l.ge.i) goto 110
        112       continue
        if(ancstr(l).eq.0) goto 111
        if(ancstr(l).eq.i) goto 110
        it=ancstr(l)
        ancstr(l)=i
        l=it
        goto 112
        111       continue
        ancstr(l)=i
        parent(l)=i
        110    continue
        100 continue
        do 200 i=1,neqns
          if(parent(i).eq.0) parent(i)=neqns+1
          200 continue
          parent(neqns+1)=0
          if(ldbg) write(idbg,6010)
          if(ldbg) write(idbg,6000) (i,parent(i),i=1,neqns)
          6000 format(2i6)
          6010 format(' parent')
          return
  end subroutine genpaq

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine genbtq(xadj,adjncy,invp,iperm,parent,btree,zpiv,izz,neqns)

    implicit none

    integer(kind=kint), intent(in)  :: xadj(:),adjncy(:),parent(:),invp(:),iperm(:),zpiv(:)
    integer(kind=kint), intent(out) :: btree(:,:) ! btree is (2,:)
    integer(kind=kint), intent(in)  :: neqns
    integer(kind=kint), intent(out) :: izz

    integer(kind=kint) :: i,j,k,l,ip,ib,inext

    do 10 i=1,neqns+1
      btree(1,i)=0
      btree(2,i)=0
      10 continue
      do 100 i=1,neqns+1
        ip=parent(i)
        if(ip.le.0) goto 100
        ib=btree(1,ip)
        if(ib.eq.0) then
          btree(1,ip)=i
        else
          101       continue
          inext=btree(2,ib)
          if(inext.eq.0) then
            btree(2,ib)=i
          else
            ib=inext
            goto 101
          endif
        endif
        100 continue
        !
        ! find zeropivot
        !
        do 200 i=1,neqns
          if(zpiv(i).ne.0) then
            if(btree(1,invp(i)).eq.0) then
              izz=i
              goto 210
            endif
          endif
          200 continue
          izz=0
          210 continue
          if(ldbg) write(idbg,6010)
          if(ldbg) write(idbg,6000) (i,btree(1,i),btree(2,i),i=1,neqns)
          if(ldbg) write(idbg,6020) izz
          !     if(idbg1.ge.2) write(10,6100) neqns
          !     if(idbg1.ge.2) write(10,6100) (btree(1,i),btree(2,i),i=1,neqns)
          6000 format(i6,'(',2i6,')')
          6010 format(' binary tree')
          6020 format(' the first zero pivot is ',i4)
          6100 format(2i8)
          return
  end subroutine genbtq

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rotate(xadj,adjncy,invp,iperm,parent,btree,izz,neqns,anc,adjt,irr)

    implicit none

    integer(kind=kint), intent(in)  :: xadj(:),adjncy(:),parent(:),btree(:,:)
    integer(kind=kint), intent(out) :: anc(:),adjt(:),invp(:),iperm(:)
    integer(kind=kint), intent(in)  :: neqns,izz
    integer(kind=kint), intent(out) :: irr

    integer(kind=kint) :: i,j,k,l,izzz,nanc,loc,locc,ll,kk,iy

    !----------------------------------------------------------------------
    !     irr return code irr=0 node izz is not a bottom node
    !                     irr=1          is a bottom node then rotation is
    !                                    performed
    !
    !----------------------------------------------------------------------
    if(izz.eq.0) then
      irr=0
      return
    endif
    izzz=invp(izz)
    if(btree(1,izzz).ne.0) then
      irr=0
      !         return
    endif
    irr=1
    !
    !  ancestors of izzz
    !
    nanc=0
    loc=izzz
    100 continue
    nanc=nanc+1
    anc(nanc)=loc
    loc=parent(loc)
    if(loc.ne.0) goto 100
    !
    !  to find the eligible node from ancestors of izz
    !
    !     adjt = Adj(Tree(y))
    l=1
    200 continue
    do 210 i=1,neqns
      adjt(i)=0
      210 continue
      locc=anc(l)
      220 continue
      loc=locc
      locc=btree(1,loc)
      if(locc.ne.0) goto 220
      230 continue
      do 240 k=xadj(iperm(loc)),xadj(iperm(loc)+1)-1
        adjt(invp(adjncy(k)))=1
        240 continue
        if(loc.ge.anc(l)) goto 250
        locc=btree(2,loc)
        if(locc.ne.0) goto 220
        loc=parent(loc)
        goto 230
        250 continue
        do 260 ll=l+1,nanc
          if(adjt(anc(ll)).eq.0) then
            l=l+1
            goto 200
          endif
          260 continue
          if(l.eq.1) goto 500

          !
          !  anc(l-1) is the eligible node
          !
          ! (1) number the node not in Ancestor(iy)
          iy=anc(l-1)
          do 300 i=1,neqns
            adjt(i)=0
            300 continue
            do 310 ll=l,nanc
              adjt(anc(ll))=1
              310 continue
              k=0
              do 320 ll=1,neqns
                if(adjt(ll).eq.0) then
                  k=k+1
                  invp(iperm(ll))=k
                endif
                320 continue
                ! (2) followed by nodes in Ancestor(iy)-Adj(T(iy))
                330 continue
                do 340 i=1,neqns
                  adjt(i)=0
                  340 continue
                  locc=iy
                  350 continue
                  loc=locc
                  locc=btree(1,loc)
                  if(locc.ne.0) goto 350
                  360 continue
                  do 370 kk=xadj(iperm(loc)),xadj(iperm(loc)+1)-1
                    adjt(invp(adjncy(kk)))=1
                    370 continue
                    if(loc.ge.iy) goto 380
                    locc=btree(2,loc)
                    if(locc.ne.0) goto 350
                    loc=parent(loc)
                    goto 360
                    380 continue
                    do 390 ll=l,nanc
                      if(adjt(anc(ll)).eq.0) then
                        k=k+1
                        invp(iperm(anc(ll)))=k
                      endif
                      390 continue
                      ! (3) and finally number the node in Adj(t(iy))
                      do 400 ll=l,nanc
                        if(adjt(anc(ll)).ne.0) then
                          k=k+1
                          invp(iperm(anc(ll)))=k
                        endif
                        400 continue
                        goto 600
                        !
                        ! izz can be numbered last
                        !
                        500 continue
                        k=0
                        do 510 i=1,neqns
                          if(i.eq.izzz) goto 510
                          k=k+1
                          invp(iperm(i))=k
                          510 continue
                          invp(iperm(izzz))=neqns
                          !
                          ! set iperm
                          !
                          600 continue
                          do 610 i=1,neqns
                            iperm(invp(i))=i
                            610 continue
                            if(ldbg) write(idbg,6000) (invp(i),i=1,neqns)
                            6000 format(10i6)
                            return
  end subroutine rotate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bringu(zpiv,iperm,invp,parent,izz,neqns,irr)

    implicit none

    integer(kind=kint), intent(in)  :: zpiv(:),parent(:)
    integer(kind=kint), intent(out) :: iperm(:),invp(:)
    integer(kind=kint), intent(in)  :: neqns,izz
    integer(kind=kint), intent(out) :: irr

    integer(kind=kint) :: i,j,k,l,ib0,ib,ibp,izzp

    !----------------------------------------------------------------------
    !
    !      bringu brings up zero pivots from bottom of the elimination tree
    !      to higher nodes
    !
    !      irr = 0     complete
    !          = 1     impossible
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------

    irr=0
    ib0=invp(izz)
    ib=ib0
    100 continue
    if(ib.le.0) goto 1000
    ibp=parent(ib)
    izzp=iperm(ibp)
    if(zpiv(izzp).eq.0) goto 110
    ib=ibp
    goto 100
    110 continue
    invp(izz)=ibp
    invp(izzp)=ib0
    iperm(ibp)=izz
    iperm(ib0)=izzp
    if(ldbg) then
      do 200 i=1,neqns
        if(invp(iperm(i)).ne.i) goto 210
        if(iperm(invp(i)).ne.i) goto 210
        200    continue
        goto 220
        210    continue
        write(20,*) 'permutation error'
        stop
    endif
    220 continue
    return
    1000 continue
    irr=1
  end subroutine bringu

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine posord(parent,btree,invp,iperm,pordr,nch,neqns,iw,qarent,mch)

    implicit none

    integer(kind=kint), intent(in)  :: btree(:,:),qarent(:)
    integer(kind=kint), intent(out) :: pordr(:),invp(:),iperm(:),nch(:),iw(:),parent(:),mch(0:neqns+1)
    integer(kind=kint), intent(in)  :: neqns

    integer(kind=kint) :: i,j,k,l,locc,loc,locp,invpos,ipinv,ii

    do 5 i=1,neqns
      mch(i)=0
      pordr(i)=0
      5 continue
      l=1
      locc=neqns+1
      10 continue
      loc=locc
      locc=btree(1,loc)
      if(locc.ne.0) goto 10
      locp=qarent(loc)
      mch(locp)=mch(locp)+1
      20 continue
      pordr(loc)=l
      if(l.ge.neqns) goto 1000
      l=l+1
      locc=btree(2,loc)
      if(locc.ne.0) goto 10
      loc=qarent(loc)
      locp=qarent(loc)
      mch(locp)=mch(locp)+mch(loc)+1
      goto 20
      1000 continue
      do 100 i=1,neqns
        ipinv=pordr(invp(i))
        invp(i)=ipinv
        iperm(ipinv)=i
        iw(pordr(i))=i
        100 continue
        do 110 i=1,neqns
          invpos=iw(i)
          nch(i)=mch(invpos)
          ii=qarent(invpos)
          if(ii.gt.0.and.ii.le.neqns) then
            parent(i)=pordr(ii)
          else
            parent(i)=qarent(invpos)
          endif
          110 continue
          if(ldbg) write(idbg,6020)
          if(ldbg) write(idbg,6000) (pordr(i),i=1,neqns)
          if(ldbg) write(idbg,6030)
          if(ldbg) write(idbg,6050)
          if(ldbg) write(idbg,6000) (parent(i),i=1,neqns)
          if(ldbg) write(idbg,6000) (invp(i),i=1,neqns)
          if(ldbg) write(idbg,6040)
          if(ldbg) write(idbg,6000) (iperm(i),i=1,neqns)
          if(ldbg) write(idbg,6010)
          if(ldbg) write(idbg,6000) (nch(i),i=1,neqns)
          6000 format(10i6)
          6010 format(' nch')
          6020 format(' post order')
          6030 format(/' invp ')
          6040 format(/' iperm ')
          6050 format(/' parent')
          return
  end subroutine posord

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gnleaf(xadj,adjncy,invp,iperm,pordr,nch,adjncp,xleaf,leaf,neqns,lnleaf)

    implicit none

    integer(kind=kint), intent(in)  :: xadj(:),adjncy(:),pordr(:),nch(:),invp(:),iperm(:)
    integer(kind=kint), intent(out) :: xleaf(:),leaf(:),adjncp(:)
    integer(kind=kint), intent(in)  :: neqns

    integer(kind=kint) i,j,k,l,m,n,ik,istart,ip,iq,lnleaf,lc1,lc

    l=1
    ik=0
    istart=0
    do 100 i=1,neqns
      xleaf(i)=l
      ip=iperm(i)
      do 105 k=xadj(ip),xadj(ip+1)-1
        iq=invp(adjncy(k))
        if(iq.lt.i) then
          ik=ik+1
          adjncp(ik)=iq
        endif
        105    continue
        m=ik-istart
        if(m.eq.0) goto 131
        call qqsort(adjncp(istart+1:),m)
        lc1=adjncp(istart+1)
        if(lc1.ge.i) goto 100
        leaf(l)=lc1
        l=l+1
        do 130 k=istart+2,ik
          lc=adjncp(k)
          !           if(lc.ge.i) goto 125
          if(lc1.lt.lc-nch(lc)) then
            leaf(l)=lc
            l=l+1
          endif
          125       continue
          lc1=lc
          130    continue
          ik=1
          istart=ik
          131    continue
          100 continue
          xleaf(neqns+1)=l
          lnleaf=l-1
          if(ldbg) write(idbg,6020)
          if(ldbg) write(idbg,6000) (xleaf(i),i=1,neqns+1)
          if(ldbg) write(idbg,6010) lnleaf
          if(ldbg) write(idbg,6000) (leaf(i),i=1,lnleaf)
          return
          6000 format(10i6)
          6010 format(' leaf (len = ',i6,')')
          6020 format(' xleaf')
  end subroutine gnleaf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine countclno(parent,xleaf,leaf,neqns,nstop,lncol,ir)

    implicit none

    ! Count total number of non-zero elements
    ! which include fill-in.
    ! A and C region of given sparse matrix will consider.
    ! D region will not consider because of D is treat as
    ! dens matrix.
    !
    integer(kind=kint), intent(in)  :: parent(:),xleaf(:),leaf(:)
    integer(kind=kint), intent(in)  :: neqns, nstop
    integer(kind=kint), intent(out) :: lncol, ir

    integer(kind=kint) :: i,j,k,l,nc,ks,ke,nxleaf

    nc=0
    ir=0
    l=1
    do 100 i=1,neqns
      ks=xleaf(i)
      ke=xleaf(i+1)-1
      if(ke.lt.ks) goto 100
      nxleaf=leaf(ks)
      do 110 k=ks,ke-1
        j=nxleaf
        nxleaf=leaf(k+1)
        105       continue
        if(j.ge.nxleaf) goto 110
        if(j.ge.nstop) then
          goto 100
        endif
        l=l+1
        j=parent(j)
        goto 105
        110    continue
        j=leaf(ke)
        115    continue
        if(j.ge.nstop) goto 100
        if(j.ge.i.or.j.eq.0) goto 100
        l=l+1
        j=parent(j)
        goto 115
        100 continue
        lncol=l-1
        return
  end subroutine countclno

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gnclno(parent,pordr,xleaf,leaf,xlnzr,colno,neqns,nstop,lncol,ir)

    implicit none

    integer(kind=kint), intent(in)  :: parent(:),pordr(:),xleaf(:),leaf(:)
    integer(kind=kint), intent(out) :: colno(:),xlnzr(:)
    integer(kind=kint), intent(in)  :: neqns, nstop
    integer(kind=kint), intent(out) :: lncol,ir

    integer(kind=kint) :: i,j,k,l,nc,ks,ke,nxleaf

    nc=0
    ir=0
    l=1
    do 100 i=1,neqns
      xlnzr(i)=l
      ks=xleaf(i)
      ke=xleaf(i+1)-1
      if(ke.lt.ks) goto 100
      nxleaf=leaf(ks)
      do 110 k=ks,ke-1
        j=nxleaf
        nxleaf=leaf(k+1)
        105       continue
        if(j.ge.nxleaf) goto 110
        if(j.ge.nstop) then
          goto 100
        endif
        colno(l)=j
        l=l+1
        j=parent(j)
        goto 105
        110    continue
        j=leaf(ke)
        115    continue
        if(j.ge.nstop) goto 100
        if(j.ge.i.or.j.eq.0) goto 100
        colno(l)=j
        l=l+1
        j=parent(j)
        goto 115
        100 continue
        xlnzr(neqns+1)=l
        lncol=l-1
        if(ldbg) write(idbg,6010)
        !     if(idbg1.ne.0) write(6,6000) (xlnzr(i),i=1,neqns+1)
        if(ldbg) write(idbg,6020) lncol
        if(ldbg) then
          do 200 k=1,neqns
            write(idbg,6100) k
            write(idbg,6000) (colno(i),i=xlnzr(k),xlnzr(k+1)-1)
            200    continue
        endif
        6000 format(10i4)
        6010 format(' xlnzr')
        6020 format(' colno (lncol =',i10,')')
        6100 format(/' row = ',i6)
        return
  end subroutine gnclno

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine qmdrch(root,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze,nbrhd)

    implicit none

    integer(kind=kint), intent(in)  :: deg(:),xadj(:),adjncy(:)
    integer(kind=kint), intent(out) :: rchset(:),marker(:),nbrhd(:)
    integer(kind=kint), intent(in)  :: root
    integer(kind=kint), intent(out) :: nhdsze,rchsze

    integer(kind=kint) :: i,j,k,l, istrt, istop, jstrt, jstop, nabor, node

    nhdsze=0
    rchsze=0
    istrt=xadj(root)
    istop=xadj(root+1)-1
    if(istop.lt.istrt) return
    do 600 i=istrt,istop
      nabor=adjncy(i)
      if(nabor.eq.0) return
      if(marker(nabor).ne.0) goto 600
      if(deg(nabor).lt.0) goto 200
      rchsze=rchsze+1
      rchset(rchsze)=nabor
      marker(nabor)=1
      goto 600
      200    marker(nabor)=-1
      nhdsze=nhdsze+1
      nbrhd(nhdsze)=nabor
      300    jstrt=xadj(nabor)
      jstop=xadj(nabor+1)-1
      do 500 j=jstrt,jstop
        node=adjncy(j)
        nabor=-node
        if(node) 300,600,400
        400       if(marker(node).ne.0) goto 500
        rchsze=rchsze+1
        rchset(rchsze)=node
        marker(node)=1
        500    continue
        600 continue
        return
  end subroutine qmdrch

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine qmdupd(xadj,adjncy,nlist,list,deg,qsize,qlink,marker,rchset,nbrhd)

    implicit none

    integer(kind=kint), intent(in)  :: adjncy(:),list(:),xadj(:)
    integer(kind=kint), intent(out) :: marker(:),nbrhd(:),rchset(:),deg(:),qsize(:),qlink(:)
    integer(kind=kint), intent(in)  :: nlist

    integer(kind=kint) :: i,j,k,l, deg0,deg1,il,inhd,inode,irch,jstrt,jstop,mark,nabor,nhdsze,node,rchsze

    if(nlist.le.0) return
    deg0=0
    nhdsze=0
    do 200 il=1,nlist
      node=list(il)
      deg0=deg0+qsize(node)
      jstrt=xadj(node)
      jstop=xadj(node+1)-1

      do 100 j=jstrt,jstop
        nabor=adjncy(j)
        if(marker(nabor).ne.0.or.deg(nabor).ge.0) goto 100
        marker(nabor)=-1
        nhdsze=nhdsze+1
        nbrhd(nhdsze)=nabor
        100    continue
        200 continue

        if(nhdsze.gt.0) call qmdmrg(xadj,adjncy,deg,qsize,qlink,marker,deg0,nhdsze,nbrhd,rchset,nbrhd(nhdsze+1:))
        do 600 il=1,nlist
          node=list(il)
          mark=marker(node)
          if(mark.gt.1.or.mark.lt.0) goto 600
          call qmdrch(node,xadj,adjncy,deg,marker,rchsze,rchset,nhdsze,nbrhd)
          deg1=deg0
          if(rchsze.le.0) goto 400
          do 300 irch=1,rchsze
            inode=rchset(irch)
            deg1=deg1+qsize(inode)
            marker(inode)=0
            300    continue
            400    deg(node)=deg1-1
            if(nhdsze.le.0) goto 600
            do 500 inhd=1,nhdsze
              inode=nbrhd(inhd)
              marker(inode)=0
              500    continue
              600 continue
              return
  end subroutine qmdupd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine qmdot(root,xadj,adjncy,marker,rchsze,rchset,nbrhd)

    implicit none

    integer(kind=kint), intent(in)  :: marker(:),rchset(:),nbrhd(:),xadj(:)
    integer(kind=kint), intent(out) :: adjncy(:)
    integer(kind=kint), intent(in)  :: rchsze,root

    integer(kind=kint) :: i,j,k,l,irch,inhd,node,jstrt,jstop,link,nabor

    irch=0
    inhd=0
    node=root
    100 jstrt=xadj(node)
    jstop=xadj(node+1)-2
    if(jstop.lt.jstrt) goto 300
    do 200 j=jstrt,jstop
      irch=irch+1
      adjncy(j)=rchset(irch)
      if(irch.ge.rchsze) goto 400
      200 continue
      300 link=adjncy(jstop+1)
      node=-link
      if(link.lt.0) goto 100
      inhd=inhd+1
      node=nbrhd(inhd)
      adjncy(jstop+1)=-node
      goto 100
      400 adjncy(j+1)=0
      do 600 irch=1,rchsze
        node=rchset(irch)
        if(marker(node).lt.0) goto 600
        jstrt=xadj(node)
        jstop=xadj(node+1)-1
        do 500 j=jstrt,jstop
          nabor=adjncy(j)
          if(marker(nabor).ge.0) goto 500
          adjncy(j)=root
          goto 600
          500    continue
          600 continue
          return
  end subroutine qmdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine qmdmrg(xadj,adjncy,deg,qsize,qlink,marker,deg0,nhdsze,nbrhd,rchset,ovrlp)

    implicit none

    integer(kind=kint), intent(in)  :: adjncy(:),nbrhd(:),xadj(:)
    integer(kind=kint), intent(out) :: deg(:),marker(:),rchset(:),ovrlp(:),qsize(:),qlink(:)
    integer(kind=kint), intent(in)  :: nhdsze

    integer(kind=kint) :: i,j,k,l, deg0,deg1,head,inhd,iov,irch,jstrt,jstop,link,lnode,mark,mrgsze,nabor,node,novrlp,rchsze,root


    if(nhdsze.le.0) return
    do 100 inhd=1,nhdsze
      root=nbrhd(inhd)
      marker(root)=0
      100 continue
      do 1400 inhd=1,nhdsze
        root=nbrhd(inhd)
        marker(root)=-1
        rchsze=0
        novrlp=0
        deg1=0
        200    jstrt=xadj(root)
        jstop=xadj(root+1)-1
        do 600 j=jstrt,jstop
          nabor=adjncy(j)
          root=-nabor
          if(nabor) 200,700,300
          300       mark=marker(nabor)
          if(mark)600,400,500
          400       rchsze=rchsze+1
          rchset(rchsze)=nabor
          deg1=deg1+qsize(nabor)
          marker(nabor)=1
          goto 600
          500       if(mark.gt.1) goto 600
          novrlp=novrlp+1
          ovrlp(novrlp)=nabor
          marker(nabor)=2
          600    continue
          700    head=0
          mrgsze=0
          do 1100 iov=1,novrlp
            node=ovrlp(iov)
            jstrt=xadj(node)
            jstop=xadj(node+1)-1
            do 800 j=jstrt,jstop
              nabor=adjncy(j)
              if(marker(nabor).ne.0) goto 800
              marker(node)=1
              goto 1100
              800       continue
              mrgsze=mrgsze+qsize(node)
              marker(node)=-1
              lnode=node
              900       link=qlink(lnode)
              if(link.le.0) goto 1000
              lnode=link
              goto 900
              1000       qlink(lnode)=head
              head=node
              1100    continue
              if(head.le.0) goto 1200
              qsize(head)=mrgsze
              deg(head)=deg0+deg1-1
              marker(head)=2
              1200    root=nbrhd(inhd)
              marker(root)=0
              if(rchsze.le.0) goto 1400
              do 1300 irch=1,rchsze
                node=rchset(irch)
                marker(node)=0
                1300    continue
                1400 continue
                return
  end subroutine qmdmrg

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lduDecomposeC(xlnzr_a,colno_a,invp_a,iperm_a,ndeg,nttbr_c,irow_c,jcol_c,ncol,nrow,xlnzr_c,colno_c,lncol_c)
    ! find fill-in position in C which placed under A and set it in xlnzr_c, colno_c

    use m_crs_matrix
    implicit none
    !input
    integer(kind=kint), intent(in) :: xlnzr_a(:)
    integer(kind=kint), intent(in) :: colno_a(:)
    integer(kind=kint), intent(in) :: iperm_a(:)
    integer(kind=kint), intent(in) :: invp_a(:)
    integer(kind=kint), intent(in) :: ndeg
    integer(kind=kint), intent(in) :: nttbr_c
    integer(kind=kint), intent(in) :: irow_c(:)
    integer(kind=kint), intent(inout) :: jcol_c(:)
    integer(kind=kint), intent(in) :: ncol
    integer(kind=kint), intent(in) :: nrow

    !output
    integer(kind=kint), pointer :: xlnzr_c(:)
    integer(kind=kint), pointer :: colno_c(:)
    integer(kind=kint), intent(out) :: lncol_c

    ! internal
    integer(kind=kint) :: i,j,k,l,m,n
    integer(kind=kint) :: ks, ke, ipass, ierr
    logical, allocatable :: cnz(:)
    type(crs_matrix) :: crs_c

    !permtate column in C for crs_c
    do i=1,nttbr_c
      jcol_c(i)=invp_a(jcol_c(i))
    end do

    ! make Compact Column Storoge using symbolic information.
    call symbolicirjctocrs(ndeg, nttbr_c, irow_c, jcol_c, ncol, nrow, crs_c)

    ! symbolic LDU factorization for C matrix
    allocate(cnz(ncol), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    do ipass = 1,2
      lncol_c = 0
      do k=1,nrow
        ! set cnz as non-zero pattern of C
        cnz = .false.
        ks = crs_c%ia(k)
        ke = crs_c%ia(k+1)-1
        if (ke .lt. ks) then
          if (ipass .eq. 2) then
            xlnzr_c(k+1)=lncol_c+1
          end if
          cycle ! in case of zero vector, no need to check dot product.
        end if

        do i=ks,ke
          cnz(crs_c%ja(i)) = .true.
        end do

        ! check for non-zero dot product and update cnz for each point of cnz
        do i=2,ncol
          ks = xlnzr_a(i)
          ke = xlnzr_a(i+1)-1
          if (ke .lt. ks) then ! in case of column of A is zero vector.
            cycle
          end if
          do j=ks,ke
            if (cnz(colno_a(j))) then
              cnz(i) = .true.
              exit
            end if
          end do
        end do

        do i=1,ncol
          if (cnz(i)) then
            lncol_c = lncol_c + 1
            if (ipass .eq. 2) then
              colno_c(lncol_c) = i
            end if
          end if
        end do
        if (ipass .eq. 2) then
          xlnzr_c(k+1)=lncol_c + 1
        end if
      end do

      if (ipass .eq. 1) then
        allocate(xlnzr_c(nrow+1),colno_c(lncol_c), stat=ierr)
        if(ierr .ne. 0) then
          call errtrp('stop due to allocation error.')
        end if
        xlnzr_c(1)=1
      end if
    end do

    ! restore order of C column.
    do i=1,nttbr_c
      jcol_c(i)=iperm_a(jcol_c(i))
    end do

    return

  end subroutine lduDecomposeC

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine qqsort(iw,ik)

    implicit none

    integer(kind=kint), intent(out) :: iw(:)
    integer(kind=kint), intent(in)  :: ik

    integer(kind=kint) :: l,m,itemp

    !----------------------------------------------------------------------
    !     sort in increasing order up to i
    !
    !     iw   array
    !     ik   number of input/output
    !     i    deal with numbers less than this numberi
    !
    !----------------------------------------------------------------------

    if(ik.le.1) return
    do 100 l=1,ik-1
      do 110 m=l+1,ik
        if(iw(l).lt.iw(m)) goto 110
        itemp=iw(l)
        iw(l)=iw(m)
        iw(m)=itemp
        110    continue
        100 continue
        200 continue
        return
  end subroutine qqsort


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine staij1(isw,i,j,aij,dsi,ir)

    implicit none

    !----------------------------------------------------------------------
    !
    !      this routine sets an non-zero entry  of the matrix.
    !      (symmetric version)
    !
    !      (i)
    !          isw      =0    set the value
    !                   =1    add the value
    !          i        row entry
    !          j        column entry
    !          aij      value
    !
    !      (o)
    !          iv       communication array
    !
    !        #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    type(dsinfo) :: dsi
    real(kind=kreal),   intent(out)  :: aij(:) ! ndeg*ndeg
    integer(kind=kint), intent(in)   :: isw, i, j
    integer(kind=kint), intent(out)  :: ir

    integer(kind=kint) :: ndeg, neqns, nstop, ndeg2, ndeg2l, ierr
    ndeg=dsi%ndeg
    neqns=dsi%neqns
    nstop=dsi%nstop
    ndeg2=ndeg*ndeg
    ndeg2l=ndeg*(ndeg+1)/2

    ir=0
    ierr=0


    ! array allocation
    if(dsi%stage.ne.20) then
      if(dsi%stage.eq.30) write(ilog,*) 'Warning a matrix was build up but never solved.'
      !
      ! for diagonal
      !
      allocate(dsi%diag(ndeg2l,neqns), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      dsi%diag=0
      !
      ! for lower triangle
      !
      allocate(dsi%zln(ndeg2,dsi%lncol), stat=ierr)
      if(ierr .ne. 0) then
        call errtrp('stop due to allocation error.')
      end if
      dsi%zln=0
      !
      ! for dense window !TODO delete this and corresponding line in addr3()
      !
      !         allocate(dsi%dsln(ndeg2,dsi%lndsln))! because of there is no dense window
      !         dsi%dsln=0

      dsi%stage=20
    endif

    ! set value
    if(ndeg.le.2) then
      call addr0(isw,i,j,aij,dsi%invp,dsi%xlnzr,dsi%colno,dsi%diag,dsi%zln,dsi%dsln,nstop,dsi%ndeg,ir)
    elseif(ndeg.eq.3) then
      call addr3(isw,i,j,aij,dsi%invp,dsi%xlnzr,dsi%colno,dsi%diag,dsi%zln,dsi%dsln,nstop,ir)
    else
      call addrx(isw,i,j,aij,dsi%invp,dsi%xlnzr,dsi%colno,dsi%diag,dsi%zln,dsi%dsln,nstop,ndeg,ndeg2,ndeg2l,ir)
    endif
    1000 continue
    return
  end subroutine staij1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! After here, routines specilized for ndeg = 1
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! LDU decompose of A (1..nstop-1) region
  subroutine sum(ic,xlnzr,colno,zln,diag,nch,par,neqns)

    implicit none

    integer(kind=kint), intent(in)    :: xlnzr(:),colno(:),par(:)
    integer(kind=kint), intent(in)    :: ic, neqns
    real(kind=kreal),   intent(inout) :: zln(:),diag(:)
    integer(kind=kint), intent(out)   :: nch(:)

    real(kind=kreal) :: s, t, zz, piv
    integer(kind=kint) :: ks, ke, kk, k, jc, jj, j, ierr
    integer(kind=kint) :: isem
    real(kind=kreal),allocatable :: temp(:)
    integer(kind=kint),allocatable :: indx(:)
    allocate(temp(neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    isem = 0

    2 continue
    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t=0.0d0
    !        do 100 i=1,ic
    !           temp(i)=0.0d0
    ! 100    continue
    do 200 k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      s=0.0d0
      do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        if(indx(j).eq.ic) then
          s=s+temp(j)*zln(jj)
        endif
        310       continue
        !           j1=xlnzr(jc)
        !           jj=xlnzr(jc+1)-j1
        !           ss=ddoti(jj,zln(j1),colno(j1),temp)
        !           zz=zln(k)-ddoti(jj,zln(j1),colno(j1),temp)
        zz=zln(k)-s
        zln(k)=zz*diag(jc)
        temp(jc)=zz
        t=t+zz*zln(k)
        200    continue
        piv=diag(ic)-t
        if(dabs(piv).gt.rmin) then
          diag(ic)=1.0d0/piv
        endif
        1 continue
        if(isem.eq.1) then
          isem=0
          nch(ic)=-1
          kk=par(ic)
          nch(kk)=nch(kk)-1
          isem=1
        else
          goto 1
        endif
        return
  end subroutine sum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! LDU decompose of C (nstop..neqnsA+neqnsd) region
  subroutine sum1(ic,xlnzr,colno,zln,diag,par,neqns)

    implicit none

    integer(kind=kint), intent(in)    :: xlnzr(:),colno(:),par(:)
    integer(kind=kint), intent(in)    :: ic, neqns
    real(kind=kreal),   intent(inout) :: zln(:),diag(:)

    real(kind=kreal) :: s, t, zz
    integer(kind=kint) :: ks, ke, k, jc, j, jj, ierr
    real(kind=kreal),allocatable :: temp(:)
    integer(kind=kint),allocatable :: indx(:)

    ierr=0

    allocate(temp(neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t=0.0d0
    !        do 100 i=1,ic
    !           temp(i)=0.0d0
    ! 100    continue
    do 200 k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      s=0.0d0
      do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        if(indx(j).eq.ic) then
          s=s+temp(j)*zln(jj)
        endif
        310       continue
        zz=zln(k)-s
        !           j1=xlnzr(jc)
        !           jj=xlnzr(jc+1)-j1
        !           zz=zln(k)-ddoti(jj,zln(j1),colno(j1),temp)
        zln(k)=zz
        temp(jc)=zz
        !           t=t+zz*zz*diag(jc)
        200    continue
        !        diag(ic)=diag(ic)-t
        return
  end subroutine sum1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! LDU decompose and Update D region.
  subroutine sum2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln)

    implicit none

    integer(kind=kint), intent(in)    :: neqns, nstop
    integer(kind=kint), intent(in)    :: xlnzr(:),colno(:)
    real(kind=kreal),   intent(inout) :: zln(:),diag(:)
    integer(kind=kint), pointer       :: spdslnidx(:)
    real(kind=kreal),   pointer       :: spdslnval(:,:)
    integer(kind=kint), intent(out)   :: nspdsln

    real(kind=kreal) :: s, t
    integer(kind=kint) :: ks, ke, kk, k, jc, jj, j, j1,j2
    integer(kind=kint) :: ic, i, loc, ierr
    integer(kind=kint) :: ispdsln
    logical :: ftflag
    real(kind=kreal),allocatable :: temp(:)
    integer(kind=kint),allocatable :: indx(:)
    ierr=0
    allocate(temp(neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    nspdsln=0
    do ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do k=ks,ke
        jj=colno(k)
        indx(jj)=ic
      end do
      do jc=nstop,ic-1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        do jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            nspdsln=nspdsln+1
            exit
          endif
        end do
      end do
    end do
    allocate(spdslnidx(nspdsln),spdslnval(1,nspdsln), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    loc=0
    ispdsln=0
    spdslnval=0
    ftflag = .true.
    do 100 ic=nstop,neqns
      do 105 i=1,nstop
        temp(i)=0.0d0
        105    continue
        ks=xlnzr(ic)
        ke=xlnzr(ic+1)-1
        do 110 k=ks,ke
          jj=colno(k)
          temp(jj)=zln(k)
          zln(k)=temp(jj)*diag(jj)
          indx(jj)=ic
          diag(ic)=diag(ic)-temp(jj)*zln(k)
          110    continue
          do 120 jc=nstop,ic-1
            s=0.0d0
            loc=loc+1
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
              j=colno(jj)
              if(indx(j).eq.ic) then
                if (ftflag) then
                  ispdsln=ispdsln+1
                  ftflag=.false.
                end if
                s=s+temp(j)*zln(jj)
              endif
              220       continue
              spdslnidx(ispdsln)=loc
              spdslnval(1,ispdsln)=spdslnval(1,ispdsln)-s
              ftflag = .true.
              !           j1=xlnzr(jc)
              !           jj=xlnzr(jc+1)-j1
              !           dsln(loc)=dsln(loc)-ddoti(jj,zln(j1),colno(j1),temp)
              120    continue
              100 continue
              return
  end subroutine sum2_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sum3(n,dsln,diag)

    implicit none

    real(kind=kreal),   intent(inout) :: dsln(:),diag(:)
    integer(kind=kint), intent(in)    :: n

    integer(kind=kint) :: i, j, loc, ierr
    real(kind=kreal),allocatable :: temp(:)
    integer(kind=kint),allocatable :: indx(:)
    allocate(temp(n),indx(n), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    if(n.le.0) goto 1000
    indx(1)=0
    loc=1
    diag(1)=1.0d0/diag(1)
    do 100 i=2,n
      indx(i)=loc
      do 110 j=1,i-1
        dsln(loc)=dsln(loc)-dot_product(dsln(indx(i):indx(i)+j-2),dsln(indx(j):indx(j)+j-2))
        loc=loc+1
        110    continue
        temp(1:i-1)=dsln(indx(i):indx(i)+i-2)*diag(1:i-1)
        diag(i)=diag(i)-dot_product(temp(1:i-1),dsln(indx(i):indx(i)+i-2))
        dsln(indx(i):indx(i)+i-2)=temp(1:i-1)
        diag(i)=1.0d0/diag(i)
        100 continue
        1000 continue
        return
  end subroutine sum3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=kreal) function spdot2(b,zln,colno,ks,ke)

    implicit none

    integer(kind=kint), intent(in) :: colno(:)
    integer(kind=kint), intent(in) :: ks,ke
    real(kind=kreal),   intent(in) :: zln(:),b(:)

    integer(kind=kint) :: j,jj
    real(kind=kreal) :: s

    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    s=0.0d0
    do 100 jj=ks,ke
      j=colno(jj)
      s=s+zln(jj)*b(j)
      100 continue
      spdot2=s
  end function spdot2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=kreal) function ddot(a,b,n)

    implicit none

    real(kind=kreal),   intent(in) :: a(n),b(n)
    integer(kind=kint), intent(in) :: n

    real(kind=kreal)   :: s
    integer(kind=kint) :: i

    s=0.0d0
    do 100 i=1,n
      s=s+a(i)*b(i)
      100 continue
      ddot=s
      return
  end function ddot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine addr0(isw,i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ndeg,ir)

    implicit none

    integer(kind=kint), intent(in) :: isw ! 0: renew diag, dsln, zln  other: add to diag, dsln, zln
    integer(kind=kint), intent(in) :: i,j,nstop, ndeg, invp(:),xlnzr(:),colno(:)
    real(kind=kreal),   intent(inout) :: zln(:,:),diag(:,:),dsln(:,:),aij(:)
    integer(kind=kint), intent(out) :: ir

    integer(kind=kint) :: ndeg2, ii, jj, itrans, k, i0, j0, l, ks, ke
    integer(kind=kint), parameter :: idbg=0
    ndeg2=ndeg*ndeg

    ir=0
    ii=invp(i)
    jj=invp(j)
    if(idbg.ne.0) write(idbg,*) 'addr0',ii,jj,aij
    if(ii.eq.jj) then
      if(ndeg2.eq.1) then
        if(isw.eq.0) then
          diag(1,ii)=aij(1)
        else
          diag(1,ii)=diag(1,ii)+aij(1)
        endif
      elseif(ndeg2.eq.4) then
        if(isw.eq.0) then
          diag(1,ii)=aij(1)
          diag(2,ii)=aij(2)
          diag(3,ii)=aij(4)
        else
          diag(1,ii)=diag(1,ii)+aij(1)
          diag(2,ii)=diag(2,ii)+aij(2)
          diag(3,ii)=diag(3,ii)+aij(4)
        endif
      endif
      goto 1000
    endif
    itrans=0
    if(jj.gt.ii) then
      k=jj
      jj=ii
      ii=k
      itrans=1
    endif
    if(jj.ge.nstop) then
      i0=ii-nstop
      j0=jj-nstop+1
      k=i0*(i0-1)/2+j0
      if(ndeg2.eq.1) then
        dsln(1,k)=aij(1)
        goto 1000
      elseif(ndeg2.eq.4) then
        if(itrans.eq.0) then
          do 3 l=1,ndeg2
            dsln(l,k)=aij(l)
            3           continue
            goto 1000
          else
            dsln(1,k)=aij(1)
            dsln(2,k)=aij(3)
            dsln(3,k)=aij(2)
            dsln(4,k)=aij(4)
            goto 1000
        endif
      endif
    endif
    ks=xlnzr(ii)
    ke=xlnzr(ii+1)-1
    do 100 k=ks,ke
      if(colno(k).eq.jj) then
        if(isw.eq.0) then
          if(ndeg2.eq.1) then
            zln(1,k)=aij(1)
          elseif(ndeg2.eq.4) then
            if(itrans.eq.0) then
              do 4 l=1,ndeg2
                zln(l,k)=aij(l)
                4                continue
              else
                zln(1,k)=aij(1)
                zln(2,k)=aij(3)
                zln(3,k)=aij(2)
                zln(4,k)=aij(4)
            endif
          endif
        else
          if(ndeg2.eq.1) then
            zln(1,k)=zln(1,k)+aij(1)
          elseif(ndeg2.eq.4) then
            if(itrans.eq.0) then
              do 5 l=1,ndeg2
                zln(l,k)=zln(l,k)+aij(l)
                5                continue
              else
                zln(1,k)=zln(1,k)+aij(1)
                zln(2,k)=zln(2,k)+aij(3)
                zln(3,k)=zln(3,k)+aij(2)
                zln(4,k)=zln(4,k)+aij(4)
            endif
          endif
        endif
        goto 1000
      endif
      100 continue
      ir=20
      1000 continue
      return
  end subroutine addr0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! After here, routines specilized for ndeg = 3
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s3um(ic,xlnzr,colno,zln,diag,nch,par,neqns)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),par(:)
    integer(kind=kint), intent(out) :: nch(:)
    real(kind=kreal),   intent(out) :: zln(:,:), diag(:,:) ! zln(9,*),diag(6,*),temp(9,*),indx(*)
    integer(kind=kint), intent(in)  :: ic,neqns

    real(kind=kreal),   allocatable :: temp(:,:)
    integer(kind=kint), allocatable :: indx(:)
    real(kind=kreal) :: zz(9),t(6)
    integer(kind=kint) :: i,j,k,l,ks,ke,kk,jc,jj,ir, ierr

    allocate(temp(9,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    !$dir max_trips(6)
    do 100 l=1,6
      t(l)=0.0d0
      100    continue
      do 200 k=ks,ke-1
        jc=colno(k)
        indx(jc)=ic
        !$dir max_trips(9)
        do 210 l=1,9
          zz(l)=zln(l,k)
          210       continue
          do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
            j=colno(jj)
            if(indx(j).eq.ic) then
              zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(4,j)*zln(4,jj)-temp(7,j)*zln(7,jj)
              zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(5,j)*zln(4,jj)-temp(8,j)*zln(7,jj)
              zz(3)=zz(3)-temp(3,j)*zln(1,jj)-temp(6,j)*zln(4,jj)-temp(9,j)*zln(7,jj)
              zz(4)=zz(4)-temp(1,j)*zln(2,jj)-temp(4,j)*zln(5,jj)-temp(7,j)*zln(8,jj)
              zz(5)=zz(5)-temp(2,j)*zln(2,jj)-temp(5,j)*zln(5,jj)-temp(8,j)*zln(8,jj)
              zz(6)=zz(6)-temp(3,j)*zln(2,jj)-temp(6,j)*zln(5,jj)-temp(9,j)*zln(8,jj)
              zz(7)=zz(7)-temp(1,j)*zln(3,jj)-temp(4,j)*zln(6,jj)-temp(7,j)*zln(9,jj)
              zz(8)=zz(8)-temp(2,j)*zln(3,jj)-temp(5,j)*zln(6,jj)-temp(8,j)*zln(9,jj)
              zz(9)=zz(9)-temp(3,j)*zln(3,jj)-temp(6,j)*zln(6,jj)-temp(9,j)*zln(9,jj)
            endif
            310       continue

            call inv33(zln(:,k),zz,diag(:,jc))

            !$dir max_trips(9)
            do 220 l=1,9
              temp(l,jc)=zz(l)
              220       continue
              t(1)=t(1)+zz(1)*zln(1,k)+zz(4)*zln(4,k)+zz(7)*zln(7,k)
              t(2)=t(2)+zz(1)*zln(2,k)+zz(4)*zln(5,k)+zz(7)*zln(8,k)
              t(3)=t(3)+zz(2)*zln(2,k)+zz(5)*zln(5,k)+zz(8)*zln(8,k)
              t(4)=t(4)+zz(1)*zln(3,k)+zz(4)*zln(6,k)+zz(7)*zln(9,k)
              t(5)=t(5)+zz(2)*zln(3,k)+zz(5)*zln(6,k)+zz(8)*zln(9,k)
              t(6)=t(6)+zz(3)*zln(3,k)+zz(6)*zln(6,k)+zz(9)*zln(9,k)
              200    continue
              !$dir max_trips(6)
              do 320 l=1,6
                diag(l,ic)=diag(l,ic)-t(l)
                320    continue

                call inv3(diag(:,ic),ir)
                nch(ic)=-1
                kk=par(ic)
                nch(kk)=nch(kk)-1

                return
  end subroutine s3um

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s3um1(ic,xlnzr,colno,zln,diag,nch,par,neqns)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),nch(:),par(:)
    real(kind=kreal),   intent(in)  :: diag(:,:)! zln(9,*),diag(6,*)
    real(kind=kreal),   intent(out) :: zln(:,:)
    integer(kind=kint), intent(in)  :: ic,neqns

    integer(kind=kint) :: i,j,k,l,ks,ke,jc,jj,ierr
    real(kind=kreal) :: s(9),zz(9)
    real(kind=kreal),allocatable :: temp(:,:)
    integer(kind=kint),allocatable :: indx(:)

    ierr=0
    allocate(temp(9,neqns),indx(neqns), stat=ierr)
    temp = 0.0d0
    indx = 0

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    !$dir max_trip(9)
    do 100 l=1,9
      s(l)=0.0d0
      100    continue
      do 200 k=ks,ke-1
        jc=colno(k)
        indx(jc)=ic

        do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(4,j)*zln(4,jj)+temp(7,j)*zln(7,jj)
            s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(5,j)*zln(4,jj)+temp(8,j)*zln(7,jj)
            s(3)=s(3)+temp(3,j)*zln(1,jj)+temp(6,j)*zln(4,jj)+temp(9,j)*zln(7,jj)
            s(4)=s(4)+temp(1,j)*zln(2,jj)+temp(4,j)*zln(5,jj)+temp(7,j)*zln(8,jj)
            s(5)=s(5)+temp(2,j)*zln(2,jj)+temp(5,j)*zln(5,jj)+temp(8,j)*zln(8,jj)
            s(6)=s(6)+temp(3,j)*zln(2,jj)+temp(6,j)*zln(5,jj)+temp(9,j)*zln(8,jj)
            s(7)=s(7)+temp(1,j)*zln(3,jj)+temp(4,j)*zln(6,jj)+temp(7,j)*zln(9,jj)
            s(8)=s(8)+temp(2,j)*zln(3,jj)+temp(5,j)*zln(6,jj)+temp(8,j)*zln(9,jj)
            s(9)=s(9)+temp(3,j)*zln(3,jj)+temp(6,j)*zln(6,jj)+temp(9,j)*zln(9,jj)
          endif
          310       continue
          !$dir max_trip(9)
          do 320 l=1,9
            temp(l,jc)=zln(l,k)-s(l)
            zln(l,k)=temp(l,jc)
            s(l)=0.0d0
            320       continue

            200    continue

    deallocate(temp,indx)
  end subroutine s3um1


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s3um3(n,dsln,diag)

    implicit none

    real(kind=kreal),   intent(out) :: dsln(:,:),diag(:,:)!dsln(9,*),diag(6,*)
    integer(kind=kint), intent(in)  :: n

    real(kind=kreal)   :: t(9)
    integer(kind=kint) :: i,j,k,l,loc,ir, ierr
    real(kind=kreal),   allocatable :: temp(:,:)
    integer(kind=kint), allocatable :: indx(:)

    allocate(temp(9,n),indx(n), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    if(n.le.0) goto 1000
    indx(1)=1
    loc=1
    call inv3(diag(:,1),ir)
    do 100 i=2,n
      indx(i)=loc
      do 110 j=1,i-1
        call d3dot(t,dsln(:,indx(i):indx(i)+j-2), dsln(:,indx(j):indx(j)+j-2),j-1)
        !$dir max_trips(9)
        !           do 111 l=1,9
        !              dsln(l,loc)=dsln(l,loc)-t(l)
        ! 111       continue
        dsln(:,loc)=dsln(:,loc)-t(:)
        loc=loc+1
        110    continue
        call v3prod(dsln(:,indx(i):indx(i)+i-2), diag,temp,i-1)
        call d3dotl(t,temp,dsln(:,indx(i):indx(i)+i-2),i-1)
        !$dir max_trips(6)
        !        do 112 l=1,6
        !           diag(l,i)=diag(l,i)-t(l)
        ! 112    continue
        diag(:,i)=diag(:,i)-t(1:6)
        dsln(:,indx(i):indx(i)+i-2)=temp(:,1:i-1)
        call inv3(diag(:,i),ir)
        100 continue
        1000 continue
        return
  end subroutine s3um3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine d3sdot(wi,a,b,n)

    implicit none

    real(kind=kreal),   intent(in)  :: a(:,:),b(:,:) !wi(3),a(3,*),b(9,*) !
    real(kind=kreal),   intent(out) :: wi(:)
    integer(kind=kint), intent(in)  :: n

    integer(kind=kint) :: jj
    !
    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    do 100 jj=1,n
      wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(4,jj)-a(3,jj)*b(7,jj)
      wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(5,jj)-a(3,jj)*b(8,jj)
      wi(3)=wi(3)-a(1,jj)*b(3,jj)-a(2,jj)*b(6,jj)-a(3,jj)*b(9,jj)
      100 continue
      return
  end subroutine d3sdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s3pdot(bi,b,zln,colno,ks,ke)

    implicit none

    integer(kind=kint), intent(in)  :: colno(:)
    real(kind=kreal),   intent(in)  :: zln(:,:),b(:,:) !zln(9,*),b(3,*),bi(3) !
    real(kind=kreal),   intent(out) :: bi(:)
    integer(kind=kint), intent(in)  :: ks,ke

    integer(kind=kint) :: j,jj
    !
    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    do 100 jj=ks,ke
      j=colno(jj)
      bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(4,jj)*b(2,j)-zln(7,jj)*b(3,j)
      bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(5,jj)*b(2,j)-zln(8,jj)*b(3,j)
      bi(3)=bi(3)-zln(3,jj)*b(1,j)-zln(6,jj)*b(2,j)-zln(9,jj)*b(3,j)
      100 continue
      return
  end subroutine s3pdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inv33(zln,zz,diag)

    implicit none

    real(kind=kreal), intent(in)  :: zz(9),diag(6)
    real(kind=kreal), intent(out) :: zln(9)

    zln(4)=zz(4)-zz(1)*diag(2)
    zln(7)=zz(7)-zz(1)*diag(4)-zln(4)*diag(5)
    zln(1)=zz(1)*diag(1)
    zln(4)=zln(4)*diag(3)
    zln(7)=zln(7)*diag(6)
    zln(4)=zln(4)-zln(7)*diag(5)
    zln(1)=zln(1)-zln(4)*diag(2)-zln(7)*diag(4)
    !
    zln(5)=zz(5)-zz(2)*diag(2)
    zln(8)=zz(8)-zz(2)*diag(4)-zln(5)*diag(5)
    zln(2)=zz(2)*diag(1)
    zln(5)=zln(5)*diag(3)
    zln(8)=zln(8)*diag(6)
    zln(5)=zln(5)-zln(8)*diag(5)
    zln(2)=zln(2)-zln(5)*diag(2)-zln(8)*diag(4)
    !
    zln(6)=zz(6)-zz(3)*diag(2)
    zln(9)=zz(9)-zz(3)*diag(4)-zln(6)*diag(5)
    zln(3)=zz(3)*diag(1)
    zln(6)=zln(6)*diag(3)
    zln(9)=zln(9)*diag(6)
    zln(6)=zln(6)-zln(9)*diag(5)
    zln(3)=zln(3)-zln(6)*diag(2)-zln(9)*diag(4)
    return
  end subroutine inv33

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inv3(dsln,ir)

    implicit none

    real(kind=kreal) :: dsln(6),t(2)
    integer(kind=kint) :: ir

    ir=0
    if(dabs(dsln(1)).lt.rmin) then
      goto 999
    endif
    dsln(1)=1.0d0/dsln(1)
    t(1)=dsln(2)*dsln(1)
    dsln(3)=dsln(3)-t(1)*dsln(2)
    dsln(2)=t(1)
    if(dabs(dsln(3)).lt.rmin) then
      goto 999
    endif
    dsln(3)=1.0d0/dsln(3)
    t(1)=dsln(4)*dsln(1)
    dsln(5)=dsln(5)-dsln(2)*dsln(4)
    t(2)=dsln(5)*dsln(3)
    dsln(6)=dsln(6)-t(1)*dsln(4)-t(2)*dsln(5)
    dsln(4)=t(1)
    dsln(5)=t(2)
    if(dabs(dsln(6)).lt.rmin) then
      goto 999
    endif
    dsln(6)=1.0d0/dsln(6)
    !      write(6,*) "dsln",dsln(1),dsln(3),dsln(6)
    return
    999 continue
    write(ilog,*) "singular"
    dsln(1)=1.0d0
    dsln(2)=0.0d0
    dsln(3)=1.0d0
    dsln(4)=0.0d0
    dsln(5)=0.0d0
    dsln(6)=1.0d0
    return
  end subroutine inv3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine d3dot(t,a,b,n)
    implicit none

    real(kind=kreal),  intent(in)  :: a(:,:),b(:,:)
    real(kind=kreal),  intent(out) :: t(:)
    integer(kind=kint),intent(in)  :: n

    integer(kind=kint) :: l,jj
    !      double precision a(9,n),b(9,n)
    !      double precision t(9)
    !
    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !      it might be 'DENS' kitayama
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    !$dir max_trips(9)
    do 10 l=1,9
      t(l)=0.0d0
      10 continue
      do 100 jj=1,n
        t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
        t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
        t(3)=t(3)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
        t(4)=t(4)+a(1,jj)*b(2,jj)+a(4,jj)*b(5,jj)+a(7,jj)*b(8,jj)
        t(5)=t(5)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
        t(6)=t(6)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
        t(7)=t(7)+a(1,jj)*b(3,jj)+a(4,jj)*b(6,jj)+a(7,jj)*b(9,jj)
        t(8)=t(8)+a(2,jj)*b(3,jj)+a(5,jj)*b(6,jj)+a(8,jj)*b(9,jj)
        t(9)=t(9)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
        100 continue
        return
  end subroutine d3dot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine v3prod(zln,diag,zz,n)

    implicit none

    real(kind=kreal),   intent(in)  :: zln(:,:),diag(:,:)
    real(kind=kreal),   intent(out) :: zz(:,:)
    integer(kind=kint), intent(in)  :: n

    integer(kind=kint) :: i

    do 100 i=1,n
      zz(4,i)=zln(4,i)-zln(1,i)*diag(2,i)
      zz(7,i)=zln(7,i)-zln(1,i)*diag(4,i)-zz(4,i)*diag(5,i)
      zz(1,i)=zln(1,i)*diag(1,i)
      zz(4,i)=zz(4,i)*diag(3,i)
      zz(7,i)=zz(7,i)*diag(6,i)
      zz(4,i)=zz(4,i)-zz(7,i)*diag(5,i)
      zz(1,i)=zz(1,i)-zz(4,i)*diag(2,i)-zz(7,i)*diag(4,i)
      !
      zz(5,i)=zln(5,i)-zln(2,i)*diag(2,i)
      zz(8,i)=zln(8,i)-zln(2,i)*diag(4,i)-zz(5,i)*diag(5,i)
      zz(2,i)=zln(2,i)*diag(1,i)
      zz(5,i)=zz(5,i)*diag(3,i)
      zz(8,i)=zz(8,i)*diag(6,i)
      zz(5,i)=zz(5,i)-zz(8,i)*diag(5,i)
      zz(2,i)=zz(2,i)-zz(5,i)*diag(2,i)-zz(8,i)*diag(4,i)
      !
      zz(6,i)=zln(6,i)-zln(3,i)*diag(2,i)
      zz(9,i)=zln(9,i)-zln(3,i)*diag(4,i)-zz(6,i)*diag(5,i)
      zz(3,i)=zln(3,i)*diag(1,i)
      zz(6,i)=zz(6,i)*diag(3,i)
      zz(9,i)=zz(9,i)*diag(6,i)
      zz(6,i)=zz(6,i)-zz(9,i)*diag(5,i)
      zz(3,i)=zz(3,i)-zz(6,i)*diag(2,i)-zz(9,i)*diag(4,i)
      100 continue
      return
  end subroutine v3prod
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine d3dotl(t,a,b,n)
    implicit none

    real(kind=kreal),   intent(in)  :: a(:,:),b(:,:)
    real(kind=kreal),   intent(out) :: t(:)
    integer(kind=kint), intent(in)  :: n

    integer(kind=kint) :: l,jj
    !      double precision t(6),a(9,n),b(9,n)
    !
    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    !$dir max_trips(6)
    do 10 l=1,6
      t(l)=0.0d0
      10 continue
      do 100 jj=1,n
        t(1)=t(1)+a(1,jj)*b(1,jj)+a(4,jj)*b(4,jj)+a(7,jj)*b(7,jj)
        t(2)=t(2)+a(2,jj)*b(1,jj)+a(5,jj)*b(4,jj)+a(8,jj)*b(7,jj)
        t(3)=t(3)+a(2,jj)*b(2,jj)+a(5,jj)*b(5,jj)+a(8,jj)*b(8,jj)
        t(4)=t(4)+a(3,jj)*b(1,jj)+a(6,jj)*b(4,jj)+a(9,jj)*b(7,jj)
        t(5)=t(5)+a(3,jj)*b(2,jj)+a(6,jj)*b(5,jj)+a(9,jj)*b(8,jj)
        t(6)=t(6)+a(3,jj)*b(3,jj)+a(6,jj)*b(6,jj)+a(9,jj)*b(9,jj)
        100 continue
        return
  end subroutine d3dotl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine addr3(isw,i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ir)

    implicit none

    integer(kind=kint), intent(in)   :: invp(:),xlnzr(:),colno(:)
    real(kind=kreal),   intent(in)   :: aij(:) !zln(9,*),diag(6,*),dsln(9,*),aij(9)
    real(kind=kreal),   intent(out)  :: zln(:,:),diag(:,:),dsln(:,:) !zln(9,*),diag(6,*),dsln(9,*),aij(9)
    integer(kind=kint), intent(in)   :: isw,i,j,nstop
    integer(kind=kint), intent(out)  :: ir

    integer(kind=kint), parameter :: ndeg2=9
    integer(kind=kint), parameter :: ndeg2l=6
    integer(kind=kint) :: k,l,ii,jj,itrans,i0,j0,ks,ke

    ir=0
    ii=invp(i)
    jj=invp(j)
    if(ldbg) write(idbg,*) 'addr3',ii,jj,aij

    ! diagonal
    if(ii.eq.jj) then
      diag(1,ii)=aij(1)
      diag(2,ii)=aij(2)
      diag(3,ii)=aij(5)
      diag(4,ii)=aij(3)
      diag(5,ii)=aij(6)
      diag(6,ii)=aij(9)
      goto 1000
    endif
    itrans=0
    if(jj.gt.ii) then
      k=jj
      jj=ii
      ii=k
      itrans=1
    endif

    ! D region
    if(jj.ge.nstop) then
      i0=ii-nstop
      j0=jj-nstop+1
      k=i0*(i0-1)/2+j0
      if(itrans.eq.0) then
        do 110 l=1,ndeg2
          dsln(l,k)=aij(l)
          110       continue
          goto 1000
        else
          dsln(1,k)=aij(1)
          dsln(2,k)=aij(4)
          dsln(3,k)=aij(7)
          dsln(4,k)=aij(2)
          dsln(5,k)=aij(5)
          dsln(6,k)=aij(8)
          dsln(7,k)=aij(3)
          dsln(8,k)=aij(6)
          dsln(9,k)=aij(9)
          goto 1000
      endif
    endif

    ! A and C region
    ks=xlnzr(ii)
    ke=xlnzr(ii+1)-1
    do 100 k=ks,ke
      if(colno(k).eq.jj) then
        if(itrans.eq.0) then
          do 120 l=1,ndeg2
            zln(l,k)=aij(l)
            120          continue
          else
            zln(1,k)=aij(1)
            zln(2,k)=aij(4)
            zln(3,k)=aij(7)
            zln(4,k)=aij(2)
            zln(5,k)=aij(5)
            zln(6,k)=aij(8)
            zln(7,k)=aij(3)
            zln(8,k)=aij(6)
            zln(9,k)=aij(9)
        endif
        goto 1000
      endif
      100 continue
      ir=20
      1000 continue
      return
  end subroutine addr3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s2um(ic,xlnzr,colno,zln,diag,nch,par,neqns)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),par(:)
    integer(kind=kint), intent(out) :: nch(:)
    real(kind=kreal),   intent(out) :: zln(:,:),diag(:,:) ! zln(6,*), diag(3,*)
    integer(kind=kint), intent(in)  :: ic,neqns

    integer(kind=kint) :: i,j,k,l,ks,ke,jj,jc,ir,kk, ierr
    real(kind=kreal),   allocatable :: temp(:,:)
    integer(kind=kint), allocatable :: indx(:)
    real(kind=kreal) :: s(4),zz(4),t(3)

    allocate(temp(4,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t(1)=0.0d0
    t(2)=0.0d0
    t(3)=0.0d0
    do 200 k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      do 210 l=1,4
        s(l)=0.0d0
        zz(l)=zln(l,k)
        210       continue
        do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            zz(1)=zz(1)-temp(1,j)*zln(1,jj)-temp(3,j)*zln(3,jj)
            zz(2)=zz(2)-temp(2,j)*zln(1,jj)-temp(4,j)*zln(3,jj)
            zz(3)=zz(3)-temp(1,j)*zln(2,jj)-temp(3,j)*zln(4,jj)
            zz(4)=zz(4)-temp(2,j)*zln(2,jj)-temp(4,j)*zln(4,jj)
          endif
          310       continue
          call inv22(zln(:,k),zz,diag(:,jc))
          do 220 l=1,4
            temp(l,jc)=zz(l)
            220       continue
            t(1)=t(1)+zz(1)*zln(1,k)+zz(3)*zln(3,k)
            t(2)=t(2)+zz(2)*zln(1,k)+zz(4)*zln(3,k)
            t(3)=t(3)+zz(2)*zln(2,k)+zz(4)*zln(4,k)
            200    continue
            diag(1,ic)=diag(1,ic)-t(1)
            diag(2,ic)=diag(2,ic)-t(2)
            diag(3,ic)=diag(3,ic)-t(3)
            call inv2(diag(:,ic),ir)
            nch(ic)=-1
            kk=par(ic)
            nch(kk)=nch(kk)-1
            return
  end subroutine s2um

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s2um1(ic,xlnzr,colno,zln,diag,nch,par,neqns)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(:),colno(:),nch(:),par(:)
    real(kind=kreal),   intent(in)  :: diag(:,:) !zln(4,*),diag(3,*)
    real(kind=kreal),   intent(out) :: zln(:,:)
    integer(kind=kint), intent(in)  :: ic,neqns

    integer(kind=kint) :: i,j,k,l,ks,ke,jc,jj, ierr
    real(kind=kreal) :: s(4),zz(4)
    real(kind=kreal),   allocatable :: temp(:,:)
    integer(kind=kint), allocatable :: indx(:)

    allocate(temp(4,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    do 100 l=1,4
      s(l)=0.0d0
      100    continue
      do 200 k=ks,ke-1
        jc=colno(k)
        indx(jc)=ic
        do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            s(1)=s(1)+temp(1,j)*zln(1,jj)+temp(3,j)*zln(3,jj)
            s(2)=s(2)+temp(2,j)*zln(1,jj)+temp(4,j)*zln(3,jj)
            s(3)=s(3)+temp(1,j)*zln(2,jj)+temp(3,j)*zln(4,jj)
            s(4)=s(4)+temp(2,j)*zln(2,jj)+temp(4,j)*zln(4,jj)
          endif
          310       continue
          do 320 l=1,4
            temp(l,jc)=zln(l,k)-s(l)
            zln(l,k)=temp(l,jc)
            s(l)=0.0d0
            320       continue
            200    continue
            return
  end subroutine s2um1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s2um2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln)

    implicit none

    integer(kind=kint), intent(in) :: neqns, nstop
    integer(kind=kint), intent(in) ::  xlnzr(:),colno(:)
    real(kind=kreal),   intent(inout) :: zln(:,:),diag(:,:)!zln(4,*),diag(3,*)
    integer(kind=kint), pointer     :: spdslnidx(:)
    real(kind=kreal),   pointer     :: spdslnval(:,:)
    integer(kind=kint), intent(out) :: nspdsln

    integer(kind=kint) :: i,j,k,l,m,n, ic,ks,ke,ii,jj,jc,j1,j2,loc,ispdsln, ierr
    real(kind=kreal),   allocatable :: temp(:,:)
    integer(kind=kint), allocatable :: indx(:)
    logical :: ftflag

    allocate(temp(4,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    nspdsln=0
    do ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do k=ks,ke
        jj=colno(k)
        indx(jj)=ic
      end do
      do jc=nstop,ic-1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        do jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            nspdsln=nspdsln+1
            exit
          endif
        end do
      end do
    end do
    allocate(spdslnidx(nspdsln),spdslnval(4,nspdsln), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    loc=0
    ispdsln=0
    spdslnval=0
    ftflag = .true.
    do 100 ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do 110 k=ks,ke
        jj=colno(k)
        temp(1,jj)=zln(1,k)
        temp(2,jj)=zln(2,k)
        temp(3,jj)=zln(3,k)
        temp(4,jj)=zln(4,k)
        !           call inv22(zln(1,k),temp(1,jj),diag(1,jj))
        zln(3,k)=temp(3,jj)-temp(1,jj)*diag(2,jj)
        zln(1,k)=temp(1,jj)*diag(1,jj)
        zln(3,k)=zln(3,k)*diag(3,jj)
        zln(1,k)=zln(1,k)-zln(3,k)*diag(2,jj)
        !
        zln(4,k)=temp(4,jj)-temp(2,jj)*diag(2,jj)
        zln(2,k)=temp(2,jj)*diag(1,jj)
        zln(4,k)=zln(4,k)*diag(3,jj)
        zln(2,k)=zln(2,k)-zln(4,k)*diag(2,jj)
        !
        diag(1,ic)=diag(1,ic)-(temp(1,jj)*zln(1,k)+temp(3,jj)*zln(3,k))
        diag(2,ic)=diag(2,ic)-(temp(1,jj)*zln(2,k)+temp(3,jj)*zln(4,k))
        diag(3,ic)=diag(3,ic)-(temp(2,jj)*zln(2,k)+temp(4,jj)*zln(4,k))
        indx(jj)=ic
        110    continue
        do 120 jc=nstop,ic-1
          loc=loc+1
          do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
            j=colno(jj)
            if(indx(j).eq.ic) then
              if (ftflag) then
                ispdsln=ispdsln+1
                ftflag=.false.
              end if
              spdslnidx(ispdsln)=loc
              spdslnval(1,ispdsln)=spdslnval(1,ispdsln)-(temp(1,j)*zln(1,jj)+temp(3,j)*zln(3,jj))
              spdslnval(2,ispdsln)=spdslnval(2,ispdsln)-(temp(2,j)*zln(1,jj)+temp(4,j)*zln(3,jj))
              spdslnval(3,ispdsln)=spdslnval(3,ispdsln)-(temp(1,j)*zln(2,jj)+temp(3,j)*zln(4,jj))
              spdslnval(4,ispdsln)=spdslnval(4,ispdsln)-(temp(2,j)*zln(2,jj)+temp(4,j)*zln(4,jj))
            endif
            220       continue
            ftflag = .true.
            120    continue
            100 continue
            return
  end subroutine s2um2_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inv22(zln,zz,diag)

    implicit none

    real(kind=kreal), intent(in)  :: zz(4),diag(3)
    real(kind=kreal), intent(out) :: zln(4)

    zln(3)=zz(3)-zz(1)*diag(2)
    zln(1)=zz(1)*diag(1)
    zln(3)=zln(3)*diag(3)
    zln(1)=zln(1)-zln(3)*diag(2)

    zln(4)=zz(4)-zz(2)*diag(2)
    zln(2)=zz(2)*diag(1)
    zln(4)=zln(4)*diag(3)
    zln(2)=zln(2)-zln(4)*diag(2)

    return
  end subroutine inv22

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inv2(dsln,ir)

    implicit none

    real(kind=kreal),   intent(out) :: dsln(3)
    integer(kind=kint), intent(out) :: ir

    real(kind=kreal) :: t

    ir=0
    if(dabs(dsln(1)).lt.rmin) then
      ir=10
      return
    endif
    dsln(1)=1.0d0/dsln(1)
    t=dsln(2)*dsln(1)
    dsln(3)=dsln(3)-t*dsln(2)
    dsln(2)=t
    if(dabs(dsln(3)).lt.rmin) then
      ir=10
      return
    endif
    dsln(3)=1.0d0/dsln(3)
    return
  end subroutine inv2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine d2sdot(wi,a,b,n)

    implicit none

    real(kind=kreal), intent(in)    :: a(:,:),b(:,:)!wi(2),a(2,*),b(4,*)
    real(kind=kreal), intent(inout) :: wi(:)
    integer(kind=kint), intent(in)  :: n

    integer(kind=kint) :: jj
    !
    !----------------------------------------------------------------------
    !
    !      d2sdot performs inner product of dens vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    do 100 jj=1,n
      wi(1)=wi(1)-a(1,jj)*b(1,jj)-a(2,jj)*b(3,jj)
      wi(2)=wi(2)-a(1,jj)*b(2,jj)-a(2,jj)*b(4,jj)
      100 continue
      return
  end subroutine d2sdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine s2pdot(bi,b,zln,colno,ks,ke)

    implicit none

    integer(kind=kint), intent(in)  :: colno(:)
    integer(kind=kint), intent(in)  :: ks,ke
    real(kind=kreal),   intent(in)  :: zln(:,:),b(:,:) !zln(4,*),b(2,*),bi(2)
    real(kind=kreal),   intent(out) :: bi(:) !zln(4,*),b(2,*),bi(2)

    integer(kind=kint) :: jj,j

    !----------------------------------------------------------------------
    !
    !      s2pdot performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------

    do 100 jj=ks,ke
      j=colno(jj)
      bi(1)=bi(1)-zln(1,jj)*b(1,j)-zln(3,jj)*b(2,j)
      bi(2)=bi(2)-zln(2,jj)*b(1,j)-zln(4,jj)*b(2,j)
      100 continue
      return
  end subroutine s2pdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine addrx(isw,i,j,aij,invp,xlnzr,colno,diag,zln,dsln,nstop,ndeg,ndeg2,ndeg2l,ir)

    implicit none

    integer(kind=kint), intent(in)  :: invp(*),xlnzr(*),colno(*)
    real(kind=kreal),   intent(in)  :: aij(ndeg,ndeg)
    real(kind=kreal),   intent(out) :: zln(ndeg,ndeg,*),diag(ndeg2l,*),dsln(ndeg,ndeg,*)
    integer(kind=kint), intent(in)  :: isw,i,j,nstop,ndeg,ndeg2,ndeg2l
    integer(kind=kint), intent(out) :: ir

    integer(kind=kint) :: ii,jj,k,l,m,n,ks,ke,itrans,i0,j0

    ir=0
    ii=invp(i)
    jj=invp(j)
    if(ldbg) write(idbg,*) 'addrx',ii,jj,aij
    if(ii.eq.jj) then
      l=0
      do 100 n=1,ndeg
        do 110 m=1,n
          l=l+1
          diag(l,ii)=aij(n,m)
          110       continue
          100    continue
          goto 1000
    endif
    itrans=0
    if(jj.gt.ii) then
      k=jj
      jj=ii
      ii=k
      itrans=1
    endif
    if(jj.ge.nstop) then
      i0=ii-nstop
      j0=jj-nstop+1
      k=i0*(i0-1)/2+j0
      if(itrans.eq.0) then
        do 120 m=1,ndeg
          do 130 n=1,ndeg
            dsln(n,m,k)=aij(n,m)
            130          continue
            120       continue
            goto 1000
          else
            do 140 m=1,ndeg
              do 150 n=1,ndeg
                dsln(n,m,k)=aij(m,n)
                150          continue
                140       continue
                goto 1000
      endif
    endif
    ks=xlnzr(ii)
    ke=xlnzr(ii+1)-1
    do 200 k=ks,ke
      if(colno(k).eq.jj) then
        if(itrans.eq.0) then
          do 160 m=1,ndeg
            do 170 n=1,ndeg
              zln(n,m,k)=aij(n,m)
              170             continue
              160          continue
            else
              do 180 m=1,ndeg
                do 190 n=1,ndeg
                  zln(n,m,k)=aij(m,n)
                  190             continue
                  180          continue
        endif
        goto 1000
      endif
      200 continue
      ir=20
      1000 continue
      return
  end subroutine addrx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dxdot(ndeg,t,a,b,l)

    implicit none

    real(kind=kreal),   intent(in)  :: a(ndeg,ndeg,*),b(ndeg,ndeg,*)
    real(kind=kreal),   intent(out) :: t(ndeg,ndeg)
    integer(kind=kint), intent(in) :: ndeg,l

    integer(kind=kint) :: k,jj,n,m
    !
    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    do 221 n=1,ndeg
      do 221 m=1,ndeg
        t(n,m)=0.0d0
        do 221 k=1,ndeg
          do 100 jj=1,l
            t(n,m)=t(n,m)+a(n,k,jj)*b(m,k,jj)
            100    continue
            221    continue
            return
  end subroutine dxdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dxdotl(ndeg,t,a,b,l)

    implicit none

    real(kind=kreal),   intent(in)  :: a(ndeg,ndeg,*),b(ndeg,ndeg,*)
    real(kind=kreal),   intent(out) :: t(ndeg,ndeg)
    integer(kind=kint), intent(in)  :: ndeg,l

    integer(kind=kint) :: n,m,jj,k
    !
    !----------------------------------------------------------------------
    !
    !      spdot1 performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    do 221 n=1,ndeg
      do 221 m=1,n
        t(n,m)=0.0d0
        do 221 k=1,ndeg
          do 100 jj=1,l
            t(n,m)=t(n,m)+a(n,k,jj)*b(m,k,jj)
            100    continue
            221    continue
            return
  end subroutine dxdotl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dxsdot(ndeg,wi,a,b,n)

    implicit none

    real(kind=kreal),   intent(in)  :: a(ndeg,*),b(ndeg,ndeg,*)
    real(kind=kreal),   intent(out) :: wi(ndeg)
    integer(kind=kint), intent(in)  :: ndeg, n

    integer(kind=kint) :: jj, k, l
    !
    !----------------------------------------------------------------------
    !
    !      dxsdot performs inner product of dens vectors
    !
    !
    !      #coded by t.arakawa
    !      #reviced by t.kitayama 20071122
    !
    !----------------------------------------------------------------------
    !
    do jj=1,n
      do k=1,ndeg
        do l=1,ndeg
          wi(l)=wi(l)-b(l,k,jj)*a(k,jj)
        end do
      end do
    end do
    return
  end subroutine dxsdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine invx(dsln,ndeg,ir)

    implicit none

    real(kind=kreal),   intent(inout) :: dsln(*)
    integer(kind=kint), intent(in)    :: ndeg
    integer(kind=kint), intent(out)   :: ir

    integer(kind=kint) :: i,j,k,l,ld,l0,k0,ll
    real(kind=kreal) :: tem,t

    ir=0
    l=1
    dsln(1)=1.0d0/dsln(1)
    do 100 i=2,ndeg
      ld=0
      l0=l
      do 110 j=1,i-1
        l=l+1
        do 120 k=1,j-1
          ld=ld+1
          dsln(l)=dsln(l)-dsln(l0+k)*dsln(ld)
          120       continue
          ld=ld+1
          110    continue
          t=0.0d0
          k0=0
          ll=0
          do 130 k=l-i+2,l
            ll=ll+1
            k0=k0+ll
            tem=dsln(k)*dsln(k0)
            t=t+tem*dsln(k)
            dsln(k)=tem
            130    continue
            l=l+1
            dsln(l)=dsln(l)-t
            dsln(l)=1.0d0/dsln(l)
            100 continue
            return
  end subroutine invx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine invxx(zln,zz,diag,ndeg)

    implicit none

    real(kind=kreal),   intent(in)  :: zz(ndeg,ndeg),diag(*)
    real(kind=kreal),   intent(out) :: zln(ndeg,ndeg)
    integer(kind=kint), intent(in) :: ndeg

    integer(kind=kint) :: i,j,k,l,m,n,loc,loc1

    zln=zz
    do 100 l=1,ndeg
      loc=0
      do 120 m=1,ndeg
        loc=loc+m
        loc1=loc+m
        do 120 n=m+1,ndeg
          zln(l,n)=zln(l,n)-zln(l,m)*diag(loc1)
          loc1=loc1+n
          120    continue
          loc=0
          do 130 m=1,ndeg
            loc=loc+m
            zln(l,m)=zln(l,m)*diag(loc)
            130    continue
            do 140 n=ndeg,1,-1
              loc=loc-1
              do 140 m=n-1,1,-1
                zln(l,m)=zln(l,m)-zln(l,n)*diag(loc)
                loc=loc-1
                140    continue
                100 continue
                return
  end subroutine invxx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sxpdot(ndeg,bi,b,zln,colno,ks,ke)

    implicit none

    integer(kind=kint), intent(in)  :: colno(*)
    real(kind=kreal),   intent(in)  :: zln(ndeg,ndeg,*),b(ndeg,*)
    real(kind=kreal),   intent(out) :: bi(ndeg)
    integer(kind=kint),intent(in)   :: ndeg,ks,ke

    integer(kind=kint) :: j,jj,m,n
    !
    !----------------------------------------------------------------------
    !
    !      sxpdot performs inner product of sparse vectors
    !
    !
    !      #coded by t.arakawa
    !
    !----------------------------------------------------------------------
    !
    do 100 jj=ks,ke
      j=colno(jj)
      do 100 m=1,ndeg
        do 100 n=1,ndeg
          bi(n)=bi(n)-zln(n,m,jj)*b(m,j)
          100 continue
          return
  end subroutine sxpdot

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sxum(ic,xlnzr,colno,zln,diag,nch,par,neqns,ndeg,ndegl)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(*),colno(*),par(*)
    integer(kind=kint), intent(out) :: nch(*)
    real(kind=kreal),   intent(out) :: zln(ndeg,ndeg,*),diag(ndegl,*)
    integer(kind=kint), intent(in)  :: ic,neqns,ndeg,ndegl

    real(kind=kreal) :: zz(ndeg,ndeg),t(ndegl)
    integer(kind=kint) :: i,j,k,l,m,n,ndeg22,ks,ke,jc,loc,jj,kk,ir, ierr
    real(kind=kreal),allocatable :: temp(:,:,:)
    integer(kind=kint),allocatable :: indx(:)

    ndeg22=ndeg*ndeg
    allocate(temp(ndeg,ndeg,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    t=0.0

    do 200 k=ks,ke-1
      jc=colno(k)
      indx(jc)=ic
      zz=zln(:,:,k)
      do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
        j=colno(jj)
        if(indx(j).eq.ic) then
          do 311 m=1,ndeg
            do 311 n=1,ndeg
              do 311 kk=1,ndeg
                zz(n,m)=zz(n,m)-temp(n,kk,j)*zln(m,kk,jj)

                311             continue
        endif
        310       continue
        call invxx(zln(1,1,k),zz,diag(1,jc),ndeg)
        temp(:,:,jc)=zz
        loc=0
        do 221 n=1,ndeg
          do 221 m=1,n
            loc=loc+1
            do 221 kk=1,ndeg
              t(loc)=t(loc)+zz(n,kk)*zln(m,kk,k)
              221       continue
              200    continue
              diag(:,ic)=diag(:,ic)-t
              call invx(diag(1,ic),ndeg,ir)
              nch(ic)=-1
              kk=par(ic)
              nch(kk)=nch(kk)-1
              return
  end subroutine sxum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sxum1(ic,xlnzr,colno,zln,diag,nch,par,neqns,ndeg,ndegl)

    implicit none

    integer(kind=kint), intent(in)  :: xlnzr(*),colno(*),nch(*),par(*)
    real(kind=kreal),   intent(in)  :: diag(ndegl,*)
    real(kind=kreal),   intent(out) :: zln(ndeg,ndeg,*)
    integer(kind=kint), intent(in)  :: ic,neqns,ndeg,ndegl

    real(kind=kreal)   :: s(ndeg,ndeg)
    integer(kind=kint) :: i,j,k,l,m,n,ks,ke,jc,jj,kk, ierr
    real(kind=kreal),allocatable :: temp(:,:,:)
    integer(kind=kint),allocatable :: indx(:)

    allocate(temp(ndeg,ndeg,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    ks=xlnzr(ic)
    ke=xlnzr(ic+1)
    do 100 m=1,ndeg
      do 100 n=1,ndeg
        s(n,m)=0.0d0
        100    continue
        do 200 k=ks,ke-1
          jc=colno(k)
          indx(jc)=ic
          do 310 jj=xlnzr(jc),xlnzr(jc+1)-1
            j=colno(jj)
            if(indx(j).eq.ic) then
              do 311 m=1,ndeg
                do 311 n=1,ndeg
                  do 311 kk=1,ndeg
                    s(n,m)=s(n,m)+temp(n,kk,j)*zln(m,kk,jj)
                    311             continue
            endif
            310       continue
            do 320 m=1,ndeg
              do 320 n=1,ndeg
                temp(n,m,jc)=zln(n,m,k)-s(n,m)
                zln(n,m,k)=temp(n,m,jc)
                s(n,m)=0.0d0
                320       continue
                200    continue
                return
  end subroutine sxum1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sxum2_child(neqns,nstop,xlnzr,colno,zln,diag,spdslnidx,spdslnval,nspdsln,ndeg,ndegl)

    implicit none

    integer(kind=kint), intent(in) :: neqns, nstop
    integer(kind=kint), intent(in) :: xlnzr(*),colno(*)
    real(kind=kreal),   intent(inout) :: zln(ndeg,ndeg,*),diag(ndegl,*)
    integer(kind=kint), pointer     :: spdslnidx(:)
    real(kind=kreal),   pointer     :: spdslnval(:,:)
    integer(kind=kint), intent(out) :: nspdsln
    integer(kind=kint), intent(in)  :: ndeg, ndegl

    integer(kind=kint) :: i,j,k,l,m,n, ic,ks,ke,ii,jj,jc,j1,j2,loc,locd,kk, ierr
    integer(kind=kint) :: ispdsln
    real(kind=kreal),   allocatable :: temp(:,:,:)
    integer(kind=kint), allocatable :: indx(:)
    logical :: ftflag

    allocate(temp(ndeg,ndeg,neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    nspdsln=0
    do ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do k=ks,ke
        jj=colno(k)
        indx(jj)=ic
      end do
      do jc=nstop,ic-1
        j1=xlnzr(jc)
        j2=xlnzr(jc+1)
        do jj=xlnzr(jc),xlnzr(jc+1)-1
          j=colno(jj)
          if(indx(j).eq.ic) then
            nspdsln=nspdsln+1
            exit
          endif
        end do
      end do
    end do
    allocate(spdslnidx(nspdsln),spdslnval(ndeg*ndeg,nspdsln), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    loc=0
    ispdsln=0
    spdslnval=0
    ftflag = .true.
    do 100 ic=nstop,neqns
      ks=xlnzr(ic)
      ke=xlnzr(ic+1)-1
      do 110 k=ks,ke
        jj=colno(k)
        do 110 m=1,ndeg
          do 110 n=1,ndeg
            temp(n,m,jj)=zln(n,m,k)
            indx(jj)=ic
            110    continue
            do 111 k=ks,ke
              jj=colno(k)
              call invxx(zln(1,1,k),temp(1,1,jj),diag(1,jj),ndeg)
              111    continue
              !
              locd=0
              do 112 n=1,ndeg
                do 112 m=1,n
                  locd=locd+1
                  do 112 k=ks,ke
                    jj=colno(k)
                    do 112 kk=1,ndeg
                      diag(locd,ic)=diag(locd,ic)-temp(n,kk,jj)*zln(m,kk,k)
                      112       continue
                      do 120 jc=nstop,ic-1
                        loc=loc+1
                        j1=xlnzr(jc)
                        j2=xlnzr(jc+1)
                        do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
                          j=colno(jj)
                          if(indx(j).eq.ic) then
                            if (ftflag) then
                              ispdsln=ispdsln+1
                              ftflag=.false.
                            end if
                            spdslnidx(ispdsln)=loc
                            do 221 m=1,ndeg
                              do 221 n=1,ndeg
                                do 221 k=1,ndeg
                                  spdslnval(ndeg*(m-1)+n,ispdsln)=spdslnval(ndeg*(m-1)+n,ispdsln)-temp(n,k,j)*zln(m,k,jj)
                                  221             continue
                          endif
                          220       continue
                          ftflag = .true.
                          120    continue
                          100 continue
                          return
  end subroutine sxum2_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sxum3(neqns,dsln,diag,ndeg,ndegl)

    implicit none

    real(kind=kreal),   intent(inout):: dsln(ndeg,ndeg,*),diag(ndegl,*)
    integer(kind=kint), intent(in) :: neqns, ndeg, ndegl

    integer(kind=kint) :: loc, locd, ir, i,j,n,m, ierr
    integer(kind=kint), allocatable :: indx(:)
    real(kind=kreal),   allocatable :: temp(:,:,:)
    real(kind=kreal),   allocatable :: t(:,:)

    allocate(indx(neqns),temp(ndeg,ndeg,neqns),t(ndeg,ndeg), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    if(neqns.le.0) goto 1000
    indx(1)=1 ! it will work...
    loc=1
    call invx(diag(1,1),ndeg,ir)
    do 100 i=2,neqns
      indx(i)=loc
      do 110 j=1,i-1
        call dxdot(ndeg,t,dsln(1,1,indx(i)),dsln(1,1,indx(j)),j-1)
        do 111 m=1,ndeg
          do 111 n=1,ndeg
            dsln(n,m,loc)=dsln(n,m,loc)-t(n,m)
            111       continue
            loc=loc+1
            110    continue
            call vxprod(ndeg,ndegl,dsln(1,1,indx(i)),diag,temp,i-1)
            call dxdotl(ndeg,t,temp,dsln(1,1,indx(i)),i-1)
            locd=0
            do 221 n=1,ndeg
              do 221 m=1,n
                locd=locd+1
                diag(locd,i)=diag(locd,i)-t(n,m)
                221    continue
                call vcopy(temp,dsln(1,1,indx(i)),ndeg*ndeg*(i-1))
                call invx(diag(1,i),ndeg,ir)
                100 continue
                1000 continue
                return
  end subroutine sxum3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine vcopy(a,c,n)
    implicit none

    integer(kind=kint) :: n
    real(kind=kreal)   :: a(n),c(n)
    !     do 100 i=1,n
    !        c(i)=a(i)
    ! 100 continue
    c=a
    return
  end subroutine vcopy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine verif0(neqns,ndeg,nttbr,irow,jcol,val,rhs,x)

    implicit none

    integer(kind=kint), intent(in) :: irow(*),jcol(*)
    integer(kind=kint), intent(in) :: neqns,ndeg,nttbr
    real(kind=kreal),   intent(in) :: val(ndeg,ndeg,*),x(ndeg,*)
    real(kind=kreal),   intent(out) :: rhs(ndeg,*)

    integer(kind=kint) :: i,j,k,l,m
    real(kind=kreal) :: rel,err
    !
    !----------------------------------------------------------------------
    !
    !     verify the solution(symmetric matrix)
    !
    !----------------------------------------------------------------------
    !
    rel=0.0d0
    do 10 i=1,neqns
      do 10 l=1,ndeg
        rel=rel+dabs(rhs(l,i))
        10 continue
        do 100 k=1,nttbr
          i=irow(k)
          j=jcol(k)
          do 101 l=1,ndeg
            do 102 m=1,ndeg
              rhs(l,i)=rhs(l,i)-val(l,m,k)*x(m,j)
              if(i.ne.j) rhs(l,j)=rhs(l,j)-val(m,l,k)*x(m,i)
              102    continue
              101    continue
              100 continue
              err=0.0d0
              do 200 i=1,neqns
                do 200 l=1,ndeg
                  err=err+dabs(rhs(l,i))
                  200 continue
                  if (m_pds_procinfo%myid .eq. 0) then
                    write(imsg,6000) err,rel,err/rel
                  end if
                  6000 format(' ***verification***(symmetric)'/&
                    &       'norm(Ax-b)            =  ',1pd20.10/&
                    &       'norm(b)               =  ',1pd20.10/&
                    &       'norm(Ax-b)/norm(b)    =  ',1pd20.10)
                  6010 format(1p4d15.7)
                  return
  end subroutine verif0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vxprod(ndeg,ndegl,zln,diag,zz,n)

    implicit none

    real(kind=kreal),   intent(in)  :: zln(ndeg*ndeg,n),diag(ndegl,n)
    real(kind=kreal),   intent(out) :: zz(ndeg*ndeg,n)
    integer(kind=kint), intent(in)  :: ndeg,ndegl,n

    integer(kind=kint) :: i

    do 100 i=1,n
      call invxx(zz(1,i),zln(1,i),diag(1,i),ndeg)
      100 continue
      return
  end subroutine vxprod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif

end module hecmw_solver_direct_parallel
