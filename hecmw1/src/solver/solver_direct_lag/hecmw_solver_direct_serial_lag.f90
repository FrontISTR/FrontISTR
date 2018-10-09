!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
! module for Parallel Direct Solver
module hecmw_solver_direct_serial_lag

  use m_child_matrix_lag
  use m_irjc_matrix_lag
  use my_hecmw_util_lag

  ! access control !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private                            ! default

  public hecmw_solve_direct_serial_lag ! only entry point of Parallel Direct Solver is public

  ! internal type definition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

  real(kind=kreal), parameter :: rmin = 1.00D-200 ! for inv3() pivot

  integer :: imsg ! output file handler

  integer, parameter :: ilog  = 16 ! according to FSTR
  logical, parameter :: ldbg  = .false.
  !integer, parameter :: idbg  = 52 ! according to FSTR
  integer :: idbg  = 10 ! debug output fort.*
  logical :: lelap = .false. ! debug out

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hecmw_solve_direct_serial_lag(nrows, ilag_sta, nttbr, pointers, indices, values, b)
    ! wrapper for parallel direct solver hecmw_solve_direct_parallel_internal()

    implicit none


    ! arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! input
    ! CRS style matrix
    integer(kind=kint), intent(in) :: nrows ! number of rows of whole square matrix, including lagrange elements
    integer(kind=kint), intent(in) :: ilag_sta ! row index of start point of lagrange elements
    integer(kind=kint), intent(in) :: nttbr    ! number of none zero elements of whole square matrix. include lagrange elements and both lower, upper matrix
    integer(kind=kint), intent(in) :: pointers(:) ! whole matrix in CRS format. size is (0:nrows)
    integer(kind=kint), intent(in) :: indices(:)  ! whole matrix in CRS format. size is (nttbr)
    real(kind=kreal),   intent(in) :: values(:)   ! whole matrix in CRS format. size is (nttbr)

    !input/output
    real(kind=kreal),   intent(inout) :: b(:) ! right hand side vector. result will return via this array.


    ! internal valuables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! stiffness matrix
    type(irjc_square_matrix), target :: a0

    ! lagrange elements
    type(irjc_mn_matrix), target :: lag

    ! degree of freedom is fixed as 1
    integer(kind=kint), parameter :: ndeg_prm = 1

    ! for ASTOM matrix
    type(irjc_square_matrix), target :: a0tmp
    type(irjc_mn_matrix), target :: lagtmp

    ! misc
    integer(kind=kint) :: i,j,k,l,ii,jj,kk,ll

    ! change CRS style matrix to irow jcol style !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! set up a0
    a0%neqns = ilag_sta - 1   ! non lagrange region only
    a0%ndeg  = ndeg_prm

    ! count ordinally A0 elements in upper triangle matrix from ASTOM
    a0%nttbr = 0
    do i= 1, nrows
      do l= pointers(i), pointers(i+1)-1
        j= indices(l)
        if (j .le. a0%neqns) then
          a0%nttbr = a0%nttbr+1
        end if
      enddo
    enddo
    !write (idbg,*) 'a0%nttbr', a0%nttbr
    allocate( a0%irow(a0%nttbr), a0%jcol(a0%nttbr), a0%val(ndeg_prm,a0%nttbr))
    allocate( a0tmp%irow(a0%nttbr), a0tmp%jcol(a0%nttbr), a0tmp%val(ndeg_prm,a0%nttbr))

    kk = 0
    do i= 1, nrows
      do l= pointers(i), pointers(i+1)-1
        j= indices(l)
        if (j .le. a0%neqns) then
          !*Lower
          kk = kk + 1
          a0tmp%irow(kk) = i
          a0tmp%jcol(kk) = j
          a0tmp%val(1,kk)=values(l)
        end if
      enddo
    enddo

    ! exchange irow and jcol, reordering
    kk=0
    do i=1,a0%neqns
      do k=1,a0%nttbr
        if (a0tmp%jcol(k) == i) then
          kk=kk+1
          a0%irow(kk)=a0tmp%jcol(k)
          a0%jcol(kk)=a0tmp%irow(k)
          a0%val(1,kk)=a0tmp%val(1,k)
        end if
      end do
    end do

    ! set up lagrange elements !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lag%nrows = nrows - ilag_sta + 1    ! number of rows of lagrange region
    lag%ncols = ilag_sta - 1            ! none lagrange region only
    lag%ndeg  = ndeg_prm

    ! count and allocate lagrange matrix
    lag%nttbr = 0
    do i= 1, nrows  ! loop for whole matrix
      do l= pointers(i), pointers(i+1)-1
        j= indices(l)
        if ((j .ge. ilag_sta) .and. (i .lt. ilag_sta) ) then
          lag%nttbr = lag%nttbr+1
        end if
      enddo
    enddo
    !write (idbg,*) 'lag%nttbr', lag%nttbr
    allocate( lag%irow(lag%nttbr), lag%jcol(lag%nttbr), lag%val(ndeg_prm,lag%nttbr))
    allocate( lagtmp%irow(lag%nttbr), lagtmp%jcol(lag%nttbr), lagtmp%val(ndeg_prm,lag%nttbr))

    kk = 0
    do i= 1, nrows
      do l= pointers(i), pointers(i+1)-1
        j= indices(l)
        if ((j .ge. ilag_sta) .and. (i .lt. ilag_sta) ) then
          kk = kk + 1
          lagtmp%irow(kk) = i
          lagtmp%jcol(kk) = j - ilag_sta +1
          lagtmp%val(1,kk)=values(l)
        end if
      enddo
    enddo

    ! exchange irow and jcol, reordering
    kk=0
    do i=1,lag%nrows
      do k=1,lag%nttbr
        if (lagtmp%jcol(k) == i) then
          kk=kk+1
          lag%irow(kk)=lagtmp%jcol(k)
          lag%jcol(kk)=lagtmp%irow(k)
          lag%val(1,kk)=lagtmp%val(1,k)
        end if
      end do
    end do
    call hecmw_solve_direct_serial_lag_in(a0,lag, b)

    return
  end subroutine hecmw_solve_direct_serial_lag

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hecmw_solve_direct_serial_lag_in(a0,lag, b)

    implicit none


    type (irjc_square_matrix), intent(inout) :: a0 ! given left side matrix assembled from sp_matrix
    type (irjc_mn_matrix), intent(inout) :: lag ! lagrange elements
    real(kind=kreal), intent(inout) :: b(:) ! (a0%neqns) right hand side value vector include both original stifness matrix and followed by lagrange right hand side value.

    logical, save :: first_time = .true.

    integer(kind=kint) :: ierr


    ! start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    imsg=99 ! set message file

    call sp_direct_parent(a0,lag,b)

    return
  end subroutine hecmw_solve_direct_serial_lag_in

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sp_direct_parent(a0, lag, b_in)

    implicit none


    !I/O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! given A0 x = b0
    type (irjc_square_matrix), intent(inout) :: a0 ! given left side matrix assembled from sp_matrix
    type (irjc_mn_matrix), intent(inout) :: lag ! lagrange elements
    real(kind=kreal), intent(inout) :: b_in(:) ! (ndeg,neqns) right hand side value vector of equation.

    !internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! for sp_direct_child
    type (child_matrix)  :: cm  ! irow, jcol matrix
    type (dsinfo) :: dsi        ! direct solver info


    integer(kind=kint)               :: neqns_a     ! number of eqns in A matrix
    integer(kind=kint)               :: neqns_lag     ! number of eqns in D matrix
    real(kind=kreal), pointer      :: dsln(:,:)   ! non-diagonal elements of dens D matrix
    real(kind=kreal), pointer      :: diag(:,:)   ! diagonal elements of dens D matrix
    type(child_matrix), pointer :: dm(:)       !divided matrices

    real(kind=kreal), allocatable :: bd(:,:) ! for right hand side value

    ! internal use
    real(kind=kreal), allocatable :: b(:,:)
    real(kind=kreal), allocatable :: oldb(:,:)
    real(kind=kreal), allocatable :: b_a(:)
    real(kind=kreal), allocatable :: b_lag(:)



    logical, save :: nusol_ready = .false.
    integer(kind=kint), save :: ndeg, nndeg, ndegt
    integer(kind=kint), save :: neqns_c, iofst_a2, iofst_c, ndm

    ! misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(kind=kint) :: ierr
    integer(kind=kint) :: i,j,k,l,m,n


    real(kind=kreal), pointer :: spdslnval(:,:),  bdbuf(:,:)
    integer(kind=kint), pointer :: spdslnidx(:)
    integer(kind=kint) :: nspdsln


    !! temporaly
    integer(kind=kint), pointer :: iperm_all_inc_lag(:), part_all_inc_lag(:), iperm_rev_inc_lag(:)
    integer(kind=kint) :: child_lag_nrows, child_lag_ncols, child_lag_nttbr, offset_irow
    integer(kind=kint) :: ii, jj
    real(kind=kreal), pointer :: dsln_lag(:,:)
    real(kind=kreal), pointer :: diag_lag(:,:)
    real(kind=kreal), allocatable :: wk(:), wk_d(:)
    integer(kind=kint) :: ks, ke


    ierr=0


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP01: get a0 from FEM data format hecMAT
    !

    ndeg=a0%ndeg
    nndeg=ndeg*ndeg
    ndegt = (ndeg+1)*ndeg/2 !triangle element in diag


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP1X: LDU decompose of givem A0
    !
    ! STEP11: build up matrix.
    ! A is given a0
    ! C is lagrange region
    ! D is kept as dsln (off-diagonal), diag (diagonal).


    neqns_a   = a0%neqns
    neqns_lag = lag%nrows
    allocate(dsln_lag(1,neqns_lag*(neqns_lag - 1)/2)) ! size of lagrange
    allocate(diag_lag(1,neqns_lag))                 ! size of lagrange
    dsln_lag=0.0d0 ! bcause of lagrange
    diag_lag=0.0d0 ! bcause of lagrange

    ! set matrix for LDU decomposition


    ! A matrix
    cm%a%ndeg  = a0%ndeg
    cm%a%neqns = a0%neqns
    cm%a%nttbr = a0%nttbr
    allocate(cm%a%irow(cm%a%nttbr), cm%a%jcol(cm%a%nttbr), cm%a%val(1, cm%a%nttbr))
    cm%a%irow(:) = a0%irow(:)
    cm%a%jcol(:) = a0%jcol(:)
    cm%a%val(1,:) = a0%val(1,:)

    ! C matrix (lagrange region)
    cm%c%ndeg = lag%ndeg
    cm%c%nttbr = lag%nttbr
    cm%c%nrows = lag%nrows
    cm%c%ncols = lag%ncols
    allocate(cm%c%irow(cm%c%nttbr), cm%c%jcol(cm%c%nttbr), cm%c%val(1, cm%c%nttbr))
    cm%c%irow(:) = lag%irow(:)
    cm%c%jcol(:) = lag%jcol(:)
    cm%c%val(1,:) = lag%val(1,:)

    cm%ndeg    = cm%a%ndeg
    cm%ista_c  = cm%a%neqns+1
    cm%neqns_t = cm%a%neqns + cm%c%nrows


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP1x: LDU decompose for given matrix
    !

    ! set up dsi for allocate matrix array for fill-in
    call matini_para(cm, dsi, ierr)

    ! set real8 value
    do i=1,cm%a%nttbr
      call staij1(0, cm%a%irow(i), cm%a%jcol(i), cm%a%val(:,i), dsi, ierr)
    end do
    do i=1,cm%c%nttbr
      !  call staij1(0, cm%c%irow(i)+cm%a%neqns, dsi%iperm(cm%c%jcol(i)), cm%c%val(:,i), dsi, ierr)
      call staij1(0, cm%c%irow(i)+cm%a%neqns, cm%c%jcol(i), cm%c%val(:,i), dsi, ierr)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! following STEP12-15 will be done in nufct0_child()
    ! and return D region
    !
    ! STEP12: LDU decompose of A (1..nstop-1)
    ! STEP13: LDU decompose of C (nstop..neqnsA+neqnsd)
    ! STEP14: update D region.

    call nufct0_child(dsi, ierr, nspdsln, spdslnidx, spdslnval, diag_lag)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP13: Receive D region and update D matrix as D' = D - D1' - D2' ...
    !
    ! D is receive as dens matrix, which format is given in s3um2() in serial solver.
    ! to decompose this dens D matrix, use s3um3() on parent.

    ! off diagonal
    do i=1,nspdsln
      dsln_lag(:,spdslnidx(i)) = dsln_lag(:,spdslnidx(i)) + spdslnval(:,i) ! because of child process dsln is already substructed in s3um2()
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP13: ! LDU decompose dens D
    !
    call nufct0_parent(dsln_lag, diag_lag, neqns_lag, ndeg)

    nusol_ready = .true.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP2X: Solve Ax=b0
    !

    ! forward substitution for A region

    ! set right hand side vector (b)
    allocate(b(ndeg,a0%neqns+lag%nrows), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    do i=1,a0%neqns+lag%nrows
      b(1,i)=b_in(i) !for ndeg=1
    end do

    ! for verify
    allocate(oldb(ndeg,a0%neqns+lag%nrows), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    do i=1, a0%neqns
      oldb(ndeg,i)=b(ndeg,i)
    end do
    do i=a0%neqns+1, a0%neqns+lag%nrows
      oldb(ndeg,i)=b(ndeg,i)
    end do

    allocate(b_a(neqns_a))
    do i=1,neqns_a
      b_a(i)=b_in(i)
    end do

    allocate(b_lag(neqns_lag))
    do i=1, neqns_lag
      b_lag(i)=b_in(neqns_a + i)
    end do

    ! old nusol0_child !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP22: forward substitution for A
    allocate(wk(neqns_a+neqns_lag), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk = 0

    do i=1,neqns_a
      wk(i)=b_a(dsi%iperm(i)) ! it sholud be permtated
    end do

    ! STEP22: forward substitution for A
    do i=1,neqns_a
      ks=dsi%xlnzr(i)
      ke=dsi%xlnzr(i+1)-1
      if(ke.lt.ks) then ! logic inverted
        cycle
      end if
      wk(i)=wk(i)-spdot2(wk,dsi%zln(1,:),dsi%colno,ks,ke)
    end do

    ! STEP23: forward substitution for C and send it (yi) to parent
    allocate(wk_d(dsi%nstop:dsi%neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if
    wk_d=0

    do i=dsi%nstop,dsi%neqns
      ks=dsi%xlnzr(i)
      ke=dsi%xlnzr(i+1)-1
      if(ke.lt.ks) then ! logic inverted
        cycle
      end if
      wk_d(i)=wk_d(i)-spdot2(wk,dsi%zln(1,:),dsi%colno,ks,ke)
    end do
    ! Now wk_d is given. it should be used in parent.


    ! end old nusol0_child !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP22 forward substitution for D region. and get Y
    !
    ! b1, b2... are sended to child processes.
    ! bd is substituted to D in locally, and results from child processes
    ! (C1-Y1 substitute, C2-Y2 substitute...) are receive from child processes.
    ! these value also add for Yd

    ! update b_lag
    do i=1, neqns_lag
      b_lag(i)=b_lag(i) + wk_d(neqns_a + i)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP22 solve Ax=b for dens matrix, using updated bd
    !
    call nusol1_parent(dsln_lag(1,:), diag_lag(1,:), b_lag, neqns_lag)
    !now b_lag is correct answer


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP23 solve A region
    !

    ! prepare Z
    do i=1,neqns_a
      wk(i)=wk(i)*dsi%diag(1,i)
    end do

    ! D region is already solved
    do i=1,neqns_lag
      wk(neqns_a + i)=b_lag(i)
    end do

    do i=dsi%neqns,1,-1
      ks=dsi%xlnzr(i)
      ke=dsi%xlnzr(i+1)-1
      if(ke.lt.ks) then
        cycle
      end if
      do k=ks,ke
        j=dsi%colno(k)
        wk(j)=wk(j)-wk(i)*dsi%zln(1,k)
      end do
    end do

    do i=1, neqns_a
      b(1,dsi%iperm(i))=wk(i)
    end do
    do i=1, neqns_lag
      b(1,neqns_a +i )=b_lag(i)
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! STEP25 restore final result
    !

    ! verify result
    call verif0(ndeg, a0%neqns, a0%nttbr, a0%irow, a0%jcol, a0%val, lag%nrows, lag%nttbr, lag%irow, lag%jcol, lag%val, oldb, b) !verify result oldb will be broken.

    do i=1,a0%neqns+lag%nrows
      b_in(i)=b(1,i)
    end do


    deallocate(spdslnidx, spdslnval)
    deallocate(dsln_lag, diag_lag)

    deallocate(b, oldb, b_a, b_lag)

    deallocate(cm%a%irow, cm%a%jcol, cm%a%val)
    deallocate(cm%c%irow, cm%c%jcol, cm%c%val)
    deallocate(dsi%zpiv)
    deallocate(dsi%iperm)
    deallocate(dsi%invp)
    deallocate(dsi%parent)
    deallocate(dsi%nch)
    deallocate(dsi%xlnzr)
    deallocate(dsi%colno)
    deallocate(dsi%diag)
    deallocate(dsi%zln)
    !deallocate(dsi%dsln) ! not allocated

    return
  end subroutine sp_direct_parent


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine errtrp(mes)
    character(*) mes
    write(ilog,*) mes

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
    izz0=0
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

    do i=1,lncol_a
      dsi%colno(i)=colno_a(i)
    end do
    do i=1, lncol_c
      dsi%colno(lncol_a + i)=colno_c(i)
    end do

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

  subroutine nufct0_child(dsi,ir, nspdsln, spdslnidx, spdslnval, diag_lag)

    implicit none
    type(dsinfo),       intent(inout) :: dsi
    integer(kind=kint), intent(out)   :: ir

    integer(kind=kint), intent(inout) :: nspdsln
    real(kind=kreal), pointer :: spdslnval(:,:), bdbuf(:,:)
    integer(kind=kint), pointer :: spdslnidx(:)

    real(kind=kreal), intent(inout) :: diag_lag(:,:)

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
      call nufct1_child(dsi%xlnzr,dsi%colno,dsi%zln,dsi%diag,dsi%neqns,dsi%parent,dsi%nch,dsi%nstop,ir, &
        nspdsln, spdslnidx, spdslnval, diag_lag)
    else if(dsi%ndeg.eq.2) then
      write(idbg,*) 'ndeg=1 only'
      stop
    else if(dsi%ndeg.eq.3) then
      write(idbg,*) 'ndeg=1 only'
      stop
    else if(dsi%ndeg.eq.6) then
      write(idbg,*) 'ndeg=1 only'
      stop
    else
      write(idbg,*) 'ndeg=1 only'
      stop
    end if

    dsi%stage=30
    1000 continue
    return
  end subroutine nufct0_child

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nufct1_child(xlnzr,colno,zln,diag,neqns,parent,nch,nstop,ir, nspdsln, spdslnidx, spdslnval, diag_lag)

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
    real(kind=kreal),   intent(out) :: diag_lag(:,:) ! diag(1,:)

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
        !      imp = m_pds_procinfo%imp
        !      call MPI_SEND(nspdsln, 1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
        !      call MPI_SEND(spdslnidx,  nspdsln,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,ierr)
        !      call MPI_SEND(spdslnval,  nspdsln,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)
        !      call MPI_SEND(diag(1,nstop),  neqns_c,MPI_REAL8,IMP,1,MPI_COMM_WORLD,ierr)

        do i=1, neqns_c
          diag_lag(1,i) = diag(1,nstop+i-1) ! lagrange region
        end do

        return
  end subroutine nufct1_child



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
      write(idbg,*) 'ndeg=1 only'
      stop
    else
      write(idbg,*) 'ndeg=1 only'
      stop
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
      write(idbg,*) 'ndeg=1 only'
      stop
    else
      write(idbg,*) 'ndeg=1 only'
      stop
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

    use m_crs_matrix_lag
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
          cycle ! in case of zero vector, no need to check dot product. not cycle?
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
      write(idbg,*) 'ndeg=1 only'
      stop
    else
      write(idbg,*) 'ndeg=1 only'
      stop
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
        !         if(isem.eq.1) then !DBG
        isem=0
        nch(ic)=-1
        kk=par(ic)
        nch(kk)=nch(kk)-1
        isem=1
        !         else  !DBG
        !            goto 1 !DBG
        !         endi !DBG
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
    integer(kind=kint) :: i

    ierr=0

    allocate(temp(neqns),indx(neqns), stat=ierr)
    if(ierr .ne. 0) then
      call errtrp('stop due to allocation error.')
    end if

    do i=1,neqns
      temp(i)=0
    end do

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
    temp=0

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
            loc=loc+1
            do 220 jj=xlnzr(jc),xlnzr(jc+1)-1
              j=colno(jj)
              if(indx(j).eq.ic) then
                if (ftflag) then
                  ispdsln=ispdsln+1
                  ftflag=.false.
                end if
                spdslnidx(ispdsln)=loc
                spdslnval(1,ispdsln)=spdslnval(1,ispdsln)-temp(j)*zln(jj)
              endif
              220       continue
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

  subroutine verif0(ndeg,neqns_a0,nttbr_a0,irow_a0,jcol_a0,val_a0,neqns_l,nttbr_l,irow_l,jcol_l,val_l,rhs,x)

    implicit none

    integer(kind=kint), intent(in) :: ndeg

    integer(kind=kint), intent(in) :: irow_a0(:),jcol_a0(:)
    integer(kind=kint), intent(in) :: neqns_a0,nttbr_a0
    real(kind=kreal),   intent(in) :: val_a0(:,:)
    integer(kind=kint), intent(in) :: irow_l(:),jcol_l(:)
    integer(kind=kint), intent(in) :: neqns_l,nttbr_l
    real(kind=kreal),   intent(in) :: val_l(:,:)

    real(kind=kreal),   intent(in) :: x(:,:)
    real(kind=kreal),   intent(out) :: rhs(:,:)

    integer(kind=kint) :: i,j,k,l,m
    real(kind=kreal) :: rel,err
    !
    !----------------------------------------------------------------------
    !
    !     verify the solution(symmetric matrix)
    !
    !     ndeg=1 only
    !
    !     include lagrange elements as L region.
    !
    !     A0        |        +
    !       A0      |        |
    !         A0    |        | neqns_a0
    !           A0  |        |
    !             A0|        +
    !     ----------+---     -
    !               |0       +
    !               | 0      | neqns_l
    !     laglange  |  0     +
    !
    !
    !----------------------------------------------------------------------
    !
    rel=0.0d0
    do 10 i=1,neqns_a0+neqns_l
      do 10 l=1,ndeg
        rel=rel+dabs(rhs(l,i))
        10 continue

        ! A0 region
        do k=1,nttbr_a0
          i=irow_a0(k)
          j=jcol_a0(k)
          do l=1,ndeg
            do m=1,ndeg
              rhs(l,i)=rhs(l,i)-val_a0(1,k)*x(m,j)
              if(i.ne.j) rhs(l,j)=rhs(l,j)-val_a0(1,k)*x(m,i)
            end do
          end do
        end do

        ! lagrange region
        do k=1,nttbr_l
          i=irow_l(k)+neqns_a0
          j=jcol_l(k)
          do l=1,ndeg
            do m=1,ndeg
              rhs(l,i)=rhs(l,i)-val_l(1,k)*x(m,j)
              if(i.ne.j) rhs(l,j)=rhs(l,j)-val_l(1,k)*x(m,i)
            end do
          end do
        end do

        err=0.0d0
        do 200 i=1,neqns_a0 + neqns_l
          do 200 l=1,ndeg
            err=err+dabs(rhs(l,i))
            200 continue

            write(imsg,6000) err,rel,err/rel
            6000 format(' ***verification***(symmetric)'/&
              &       'norm(Ax-b)            =  ',1pd20.10/&
              &       'norm(b)               =  ',1pd20.10/&
              &       'norm(Ax-b)/norm(b)    =  ',1pd20.10)
            6010 format(1p4d15.7)
            return
  end subroutine verif0

end module hecmw_solver_direct_serial_lag
