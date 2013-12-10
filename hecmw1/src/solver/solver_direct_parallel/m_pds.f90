!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2007/11/21                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takeshi Kitayama  (Univ. of Tokyo)             !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

! module for Parallel Direct Solver
module m_pds

use m_cclsmatrix

implicit none

private

public sp_LINEQ ! entry point of Parallel Direct Solver


! required global informations for solver !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for trunk and leaf solve
integer :: m_pds_neqns   ! size of given left side sparse matrix.
integer :: m_pds_ndeg    ! degree of freedum of each element in sparse matrix.

! for trunksolve
integer :: m_pds_neqns1  ! size of divided sparse matrix a1.
integer :: m_pds_neqns2  ! size of divided sparse matrix a2.
integer :: m_pds_neqnsd  ! number of columns of connection sparse matrix c. number is given as ndim is NOT multiplied.
integer :: m_pds_nd      ! size of connection dens matrix d. number is given as ndim is already multiplied.
real(8), pointer :: m_pds_d(:,:)                 ! Connection dens matrix for trunksolve. Size is (nd,nd)
type(ccls_matrix), pointer :: m_pds_c1, m_pds_c2 ! Connection matrix. size is (neqns1, neqnsd), (neqns2, neqnsd)
integer, pointer :: m_pds_part(:)                ! Right side vector partition information. Size is (neqns)
logical :: m_pds_isready_dc = .false.            ! set true after d is updated correctly.

! for leafsolve
integer, parameter :: m_pds_lenv=80000000     ! allocate v as size of lenv in leaf process.
real(8), pointer   :: m_pds_v(:) ! Communication vector for serial direct solver. 
logical :: m_pds_isready_v = .false.          ! set true after v is setted correctly via nufct0().

! process position informations on MPI processes binary tree. !!!!!!!!!!!!!!!!!!!!!!!!!!

integer :: npe     ! total number of process
integer :: myid    ! my process id
integer :: imp     ! mother  process id
integer :: ip1     ! left  side child process id
integer :: ip2     ! right side child process id
logical :: isroot  ! true if this process is root  of binary tree
logical :: istrunk ! true if this process is trunk of binary tree (not root, not leaf)
logical :: isleaf  ! true if this process is leaf  of binary tree


contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sp_LINEQ(SP_MAT)

use m_elap
implicit none


character(132),parameter::version='sp_LINEQ: $Revision: 1.2 $ $Date: 2007/11/21 02:13:45 $'

integer :: ierr
type (sp_matrix) :: sp_mat

! start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call elapinitmpi()
!write(20,*) 'sp_LINEQ: entered'; call flush(20)!DEBUG

! set location in process tree
call setbt()

if (isroot) then
!  write(20,*) 'sp_LINEQ: entering sp_direct_root'; call flush(20)!DEBUG
  call sp_direct_root(sp_mat)
else if (istrunk) then
!  write(20,*) 'sp_LINEQ: entering sp_direct_trunk'; call flush(20)!DEBUG
  call sp_direct_trunk()
else if (isleaf) then
!write(20,*) 'sp_LINEQ: entering sp_direct_leaf'; call flush(20)!DEBUG
  call sp_direct_leaf()
else
!  write(20,*) 'sp_LINEQ: never come here'; call flush(20)
  stop
end if

return
end subroutine sp_LINEQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sp_direct_root(sp_mat)

use m_irjcmatrix
use m_cclsmatrix
use m_elap

implicit none

character(132), parameter :: version='sp_direct_root: $Revision: 1.2 $ $Date: 2007/11/21 02:13:45 $'

include 'mpif.h'

!I/O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (sp_matrix) :: sp_mat ! Interface to FEM code


!internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! given A0 x = b0
type (irjc_matrix), target :: a0 ! given left side matrix assembled from sp_matrix
real(8), allocatable, dimension(:) :: r ! right hand side value vector of equation.
integer :: ll, loc
integer :: jj, ii
integer :: nprof, nprof2

integer :: ndeg

! for matrix dividing and restore valuables
integer, allocatable, dimension(:) :: iperm ! relation of index of large matrix and small matrix

! for divided matrixes (a1, a2)
type (irjc_matrix) :: a1, a2
type (ccls_matrix), target :: c1, c2

! for connection region
real(8), allocatable, dimension(:,:) :: d1, d2  ! to update D

! misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: nndeg
integer :: ierr
integer :: i,j,k,l,m,n

! for MPI
integer, dimension(MPI_STATUS_SIZE) :: istatus
logical :: do_calc

! for elaps time
logical, parameter :: elap=.true.
real(8) :: t00s, t00e
real(8) :: t10s, t10e
real(8) :: t20s, t20e, t21s, t21e
real(8) :: t30s, t30e, t31s, t31e
real(8) :: t40s, t40e
real(8) :: t50s, t50e
real(8) :: t60s, t60e

!test
integer, allocatable, dimension(:),   target :: part 
real(8), allocatable, dimension(:,:), target :: d
integer :: neqns1, neqns2, neqnsd, nd
real(8), allocatable, dimension(:) :: oldb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP01: get a0 from FEM data format SP_MAT
!

!call ptime(curt);write(20,'(a,1f15.5)') 'sp_direct_root: begin ', curt - epocht; call flush(20)!DEBUG
call elapout('sp_direct_root: begin')

call ptime(t10s)!ELAP

a0%neqns=SP_MAT%n
a0%ndeg=SP_MAT%ndeg

ndeg=a0%ndeg
nndeg=ndeg*ndeg

ll=0
nprof=SP_MAT%ipoi(SP_MAT%ncol+1)-1
nprof2=nprof/2+a0%neqns
allocate(a0%irow(nprof/2+a0%neqns),a0%jcol(nprof/2+a0%neqns),a0%val(nndeg,nprof/2+a0%neqns))
do j=1,a0%neqns
  jj=SP_MAT%iperm(j)
  ll=ll+1
  a0%irow(ll)=j
  a0%jcol(ll)=j
  a0%val(1:nndeg,ll)=SP_MAT%d(jj,1:nndeg)
  loc=SP_MAT%istat(jj)
  do while(loc.ne.0)
    ii=SP_MAT%irowno(loc)
    i=SP_MAT%invp(ii)
    if(i.gt.j) then
      ll=ll+1
      a0%irow(ll)=i
      a0%jcol(ll)=j
      a0%val(1:nndeg,ll)=SP_MAT%elm(loc,1:nndeg)
    end if
    loc=SP_MAT%irpt(loc)
  end do
end do
a0%nttbr=ll

! set right hand side vector (b)
allocate(r(a0%neqns*ndeg))
do i=1,a0%neqns
  do j=1,ndeg
    r(ndeg*(i-1)+j)=SP_MAT%b(j,i)
  end do
end do

call ptime(t10e)!ELAP
call elapout('sp_direct_root: get matrix information done ')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP02: divide given large matrix A into a1 and a2
!
call elapout('sp_direct_root: begin divide matrix ')
call ptime(t20s)!ELAP

allocate(part(a0%neqns))
allocate(iperm(a0%neqns))
part=0
iperm=0
call mkpart(a0%neqns, a0%nttbr, ndeg, a0%irow, a0%jcol, neqns1, neqns2, neqnsd, part, iperm) ! set matrix partition information

nd = neqnsd*ndeg

allocate(d(nd,nd))
d=0
call divmat(a0, part, iperm, neqns1, neqns2, neqnsd, ndeg, a1, a2, c1, c2, d) 
call ptime(t20e)!ELAP
call elapout('sp_direct_root: end divide matrix')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP03: Send divided left hand side matrixes to child process recursively.
!         To prepare solver.
!

call ptime(t30s)!ELAP
call elapout('sp_direct_root: begin send a1 ')
!send a1
call MPI_SEND(a1%neqns,          1, MPI_INTEGER,IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%nttbr,          1, MPI_INTEGER,IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%ndeg,           1, MPI_INTEGER,IP1,1,MPI_COMM_WORLD,ierr)

call MPI_SEND(a1%irow,      a1%nttbr, MPI_INTEGER,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%jcol,      a1%nttbr, MPI_INTEGER,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%val, a1%nttbr*nndeg, MPI_REAL8,  IP1,1,MPI_COMM_WORLD,ierr)
call elapout('sp_direct_root: end send a1 ')

call elapout('sp_direct_root: begin send a2 ')
!send a2
call MPI_SEND(a2%neqns,          1, MPI_INTEGER,IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%nttbr,          1, MPI_INTEGER,IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%ndeg,           1, MPI_INTEGER,IP2,1,MPI_COMM_WORLD,ierr)

call MPI_SEND(a2%irow,      a2%nttbr, MPI_INTEGER,  IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%jcol,      a2%nttbr, MPI_INTEGER,  IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%val, a2%nttbr*nndeg, MPI_REAL8,  IP2,1,MPI_COMM_WORLD,ierr)
call elapout('sp_direct_root: end send a2 ')
call ptime(t30e)!ELAP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP04: Set up trunksolver.
!         Send C matrixes to child process calcCtAC().
!         Get result D1, D2 and update D as D'= D -D1 -D2
!         LDU decompose D' and prepare dense solve.
!

call ptime(t40s)!ELAP
call elapout('sp_direct_root: begin send c1 ')
call MPI_SEND(c1%neqns,            1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ncol,             1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%nttbr,            1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ndeg,             1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ia,       c1%ncol+1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ja,        c1%nttbr, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%val, c1%nttbr*nndeg, MPI_REAL8,   IP1,1,MPI_COMM_WORLD,ierr)
call elapout('sp_direct_root: end send c1 ')

call elapout('sp_direct_root: begin send c2')
call MPI_SEND(c2%neqns,            1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ncol,             1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%nttbr,            1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ndeg,             1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ia,       c2%ncol+1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ja,        c2%nttbr, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%val, c2%nttbr*nndeg, MPI_REAL8,   IP2,1,MPI_COMM_WORLD,ierr)
call elapout('sp_direct_root: end send c2')

allocate(d1(nd,nd), d2(nd,nd))
d1=0
d2=0

call elapout('sp_direct_root: wait until receive D1, D2')
call MPI_RECV(d1, nd*nd, MPI_REAL8, IP1, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(d2, nd*nd, MPI_REAL8, IP2, 1,MPI_COMM_WORLD, istatus, ierr)
call elapout('sp_direct_root: receive d1, d2 from children')


call elapout('sp_direct_root: modify D and LDU decompose it.')
d=d-d1
d=d-d2

call densldu(d,nd)

deallocate(d1, d2)
call elapout('sp_direct_root: LDU decompose for denssol done.')

! set solver information
m_pds_neqns  = a0%neqns
m_pds_ndeg   = a0%ndeg
m_pds_neqns1 = neqns1
m_pds_neqns2 = neqns2
m_pds_neqnsd = neqnsd
m_pds_nd     = nd

m_pds_part => part
m_pds_c1   => c1
m_pds_c2   => c2
m_pds_d    => d
m_pds_isready_dc=.true. ! on m_pds

call ptime(t40e)!ELAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP05: Solve Ax=b0 using child processes.
!         Results is given as "r"
!

call ptime(t50s)!ELAP

allocate(oldb(a0%neqns*ndeg))
oldb=r
call elapout('sp_direct_root: entering trunksolve')
call trunksolve(r)
call elapout('sp_direct_root: end trunksolve')

call verif0(a0%neqns, a0%ndeg, a0%nttbr, a0%irow, a0%jcol, a0%val, oldb, r) !verify result oldb will be broken.
deallocate(oldb)

do i=1,a0%neqns
  do j=1,ndeg
    SP_MAT%x(j,i)=r(ndeg*(i-1)+j)
  end do
end do

call ptime(t50e)!ELAP
! send stop flag
call elapout('sp_direct_root: send stop flag to children')
do_calc=.false.
call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP2,1,MPI_COMM_WORLD,ierr)

call elapout('sp_direct_root: end')

!  write(20,'(a)')          '# profile sp_direct_root ########################################'
!  write(20,*) 'Elaptime (sec)'
!    write(20,'(a,1f15.5)') 'STEP01: Make matrix from SP_MAT:           ', t10e  - t10s
!    write(20,'(a,1f15.5)') 'STEP02: Divide matrix:                     ', t20e  - t20s
!    write(20,'(a,1f15.5)') 'STEP03: Send divided matrix to children:   ', t30e  - t30s
!    write(20,'(a,1f15.5)') 'STEP04: Set up trunksolver:                ', t40e  - t40s
!    write(20,'(a,1f15.5)') 'STEP05: Solve Ax=b:                        ', t50e  - t50s
!  write(20,'(a)')          '# End profile data######################################'

return 
end subroutine sp_direct_root

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sp_direct_trunk()

use m_irjcmatrix
use m_cclsmatrix
use m_elap

implicit none

character(132), parameter :: version='sp_direct_trunk: $Revision: 1.2 $ $Date: 2007/11/21 02:13:45 $'

include 'mpif.h'

!internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! given A0 x = b0
type (irjc_matrix), target :: a0 ! given left side matrix assembled from sp_matrix

integer :: ndeg

! for matrix dividing and restore valuables
integer, allocatable, dimension(:) :: iperm ! relation of index of large matrix and small matrix

! for divided matrixes (a1, a2, c1, c2)
type (irjc_matrix) :: a1, a2
type (ccls_matrix), target :: c1, c2

! for connection region
real(8), allocatable, dimension(:,:) :: d1, d2  ! to update D

! misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: nndeg
integer :: ierr
integer :: i,j,k,l,m,n

! for MPI
integer, dimension(MPI_STATUS_SIZE) :: istatus

! for elaps time
logical, parameter :: elap=.true.
real(8) :: t00s, t00e
real(8) :: t10s, t10e
real(8) :: t20s, t20e, t21s, t21e
real(8) :: t30s, t30e, t31s, t31e
real(8) :: t40s, t40e
real(8) :: t50s, t50e
real(8) :: t60s, t60e

!test
integer, allocatable, dimension(:),   target :: part 
real(8), allocatable, dimension(:,:), target :: d
integer :: neqns1, neqns2, neqnsd, nd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP01 receive a0
!

call elapout('sp_direct_trunk: begin')
call elapout('sp_direct_trunk: waiting matrix via MPI')

call MPI_RECV(a0%neqns,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a0%nttbr,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a0%ndeg,   1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
allocate(a0%irow(a0%nttbr))
allocate(a0%jcol(a0%nttbr))
allocate(a0%val(a0%ndeg*a0%ndeg,a0%nttbr))

call elapout('sp_direct_trunk: begen get matrix via MPI')
call MPI_RECV(a0%irow, a0%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a0%jcol, a0%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a0%val,  a0%nttbr*a0%ndeg*a0%ndeg, MPI_REAL8,IMP,1,MPI_COMM_WORLD,istatus,ierr)

call elapout('sp_direct_trunk: end get matrix via MPI')

ndeg=a0%ndeg
nndeg=ndeg*ndeg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP02: divide given large matrix A into a1 and a2
!
call elapout('sp_direct_trunk: begin divide matrix')

allocate(part(a0%neqns))
allocate(iperm(a0%neqns))
part=0
iperm=0
call mkpart(a0%neqns, a0%nttbr, ndeg, a0%irow, a0%jcol, neqns1, neqns2, neqnsd, part, iperm) ! set matrix partition information

nd = neqnsd*ndeg

allocate(d(nd,nd))
d=0
call divmat(a0, part, iperm, neqns1, neqns2, neqnsd, ndeg, a1, a2, c1, c2, d) 
call elapout('sp_direct_trunk: end divide matrix')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP03: Send divided left hand side matrixes to child process recursively.
!         To prepare solver.
!

call elapout('sp_direct_trunk: begin send a1 ')
!send a1
call MPI_SEND(a1%neqns,          1, MPI_INTEGER,IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%nttbr,          1, MPI_INTEGER,IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%ndeg,           1, MPI_INTEGER,IP1,1,MPI_COMM_WORLD,ierr)

call MPI_SEND(a1%irow,      a1%nttbr, MPI_INTEGER,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%jcol,      a1%nttbr, MPI_INTEGER,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a1%val, a1%nttbr*nndeg, MPI_REAL8,  IP1,1,MPI_COMM_WORLD,ierr)
call elapout('sp_direct_trunk: end send a1 ')

call elapout('sp_direct_trunk: begin send a2 ')
!send a2
call MPI_SEND(a2%neqns,          1, MPI_INTEGER,IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%nttbr,          1, MPI_INTEGER,IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%ndeg,           1, MPI_INTEGER,IP2,1,MPI_COMM_WORLD,ierr)

call MPI_SEND(a2%irow,      a2%nttbr, MPI_INTEGER,  IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%jcol,      a2%nttbr, MPI_INTEGER,  IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(a2%val, a2%nttbr*nndeg, MPI_REAL8,  IP2,1,MPI_COMM_WORLD,ierr)
call elapout('sp_direct_trunk: end send a2')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP04: Set up trunksolver.
!         Send C matrixes to child process calcCtAC().
!         Get result D1, D2 and update D as D'= D -D1 -D2
!         LDU decompose D' and prepare dense solve.
!

call elapout('sp_direct_trunk: send c1')
call MPI_SEND(c1%neqns,            1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ncol,             1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%nttbr,            1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ndeg,             1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ia,       c1%ncol+1, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%ja,        c1%nttbr, MPI_INTEGER, IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c1%val, c1%nttbr*nndeg, MPI_REAL8,   IP1,1,MPI_COMM_WORLD,ierr)

call elapout('sp_direct_trunk: send c2')
call MPI_SEND(c2%neqns,            1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ncol,             1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%nttbr,            1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ndeg,             1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ia,       c2%ncol+1, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%ja,        c2%nttbr, MPI_INTEGER, IP2,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(c2%val, c2%nttbr*nndeg, MPI_REAL8,   IP2,1,MPI_COMM_WORLD,ierr)

allocate(d1(nd,nd), d2(nd,nd))
d1=0
d2=0

call elapout('sp_direct_trunk: wait until receive D1, D2')
call MPI_RECV(d1, nd*nd, MPI_REAL8, IP1, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(d2, nd*nd, MPI_REAL8, IP2, 1,MPI_COMM_WORLD, istatus, ierr)
call elapout('sp_direct_trunk: receive d1, d2 from children')


call elapout('sp_direct_trunk: modify D and LDU decompose it.')
d=d-d1
d=d-d2

call densldu(d,nd)

deallocate(d1, d2)
call elapout('sp_direct_trunk: LDU decompose for denssol done.')

! set solver information
m_pds_neqns  = a0%neqns
m_pds_ndeg   = a0%ndeg
m_pds_neqns1 = neqns1
m_pds_neqns2 = neqns2
m_pds_neqnsd = neqnsd
m_pds_nd     = nd

m_pds_part => part
m_pds_c1   => c1
m_pds_c2   => c2
m_pds_d    => d
m_pds_isready_dc=.true. ! on m_pds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP05: Receive C matrixes and return d' to parent
!

call elapout('sp_direct_trunk: begin calcCtAC()')
call calcCtAC()
call elapout('sp_direct_trunk: end calcCtAC()')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP06: Start solver to wait right hand side vector
!

call elapout('sp_direct_trunk: begin recsolve()')
call recsolve()
call elapout('sp_direct_trunk: end recsolve()')

return 

end subroutine sp_direct_trunk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sp_direct_leaf()

use m_irjcmatrix
use m_elap

implicit none

include 'mpif.h'

!internal variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (irjc_matrix) :: a

! for MPI
integer, dimension(MPI_STATUS_SIZE) :: istatus
integer :: ierr

!misc
integer :: i,j,k,l

! for elaps time
real(8) :: t00s, t00e
real(8) :: t10s, t10e
real(8) :: t20s, t20e, t21s, t21e
real(8) :: t30s, t30e, t31s, t31e
real(8) :: t40s, t40e
real(8) :: t50s, t50e
real(8) :: t60s, t60e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP01 receive a0
!

call elapout('sp_direct_leaf: begin')
call elapout('sp_direct_leaf: waiting matrix via MPI')

call MPI_RECV(a%neqns,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a%nttbr,  1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a%ndeg,   1,MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
allocate(a%irow(a%nttbr))
allocate(a%jcol(a%nttbr))
allocate(a%val(a%ndeg*a%ndeg,a%nttbr))

call elapout('sp_direct_leaf: begin receive matrix')
call ptime(t10s)!ELAP
call MPI_RECV(a%irow, a%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a%jcol, a%nttbr, MPI_INTEGER,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call MPI_RECV(a%val,  a%nttbr*a%ndeg*a%ndeg, MPI_REAL8,IMP,1,MPI_COMM_WORLD,istatus,ierr)
call ptime(t10e)!ELAP
call elapout('sp_direct_leaf: end get matrix via MPI')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP02: Prepare serial solver
!

call ptime(t20s)!ELAP
allocate(m_pds_v(m_pds_lenv)) ! communication matrix hold in m_pds

call elapout('sp_direct_leaf: begin matini()')
call matini(a%neqns, a%nttbr, a%irow, a%jcol, m_pds_lenv, m_pds_v, ierr)
if (ierr .ne. 0) call errtrp('sp_direct_leaf: matini')
call elapout('sp_direct_leaf: end matini()')

call elapout('sp_direct_leaf: begin staij1()')
do i=1,a%nttbr
  call staij1(0, a%irow(i), a%jcol(i), a%val(1,i), m_pds_v, a%ndeg, ierr)
  if (ierr .ne. 0) call errtrp('sp_direct_leaf: staij1')
end do
call elapout('sp_direct_leaf: end staij1()')

call elapout('sp_direct_leaf: begin nufct0()')
call nufct0(m_pds_v, ierr)
if (ierr .ne. 0) call errtrp('sp_direct_leaf: nufct0')
call elapout('sp_direct_leaf: end nufct0()')

! set solver information in module m_pds
m_pds_neqns = a%neqns
m_pds_ndeg  = a%ndeg
m_pds_isready_v=.true.

call ptime(t20e)!ELAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP03: Receive C matrixes and return d' to parent
!

call elapout('sp_direct_leaf: entering calcCtAC')
call ptime(t30s)!ELAP
call calcCtAC()
call ptime(t30e)!ELAP
call elapout('sp_direct_leaf: exit calcCtAC')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! STEP04: Start solver to wait right hand side vector
!

call elapout('sp_direct_leaf: entering recsolve')
call ptime(t40s)!ELAP
call recsolve()
call ptime(t40e)!ELAP
call elapout('sp_direct_leaf: exit recsolve')


!  write(20,'(a)')          '# profile sp_direct_leaf ########################################'
!  write(20,*) 'Elaptime (sec)'
!    write(20,'(a,1f15.5)') 'STEP01: Recive a from parent:              ', t10e  - t10s
!    write(20,'(a,1f15.5)') 'STEP02: Prepare serial solver:             ', t20e  - t20s
!    write(20,'(a,1f15.5)') 'STEP03: calcCtAC:                          ', t30e  - t30s
!    write(20,'(a,1f15.5)') 'STEP04: recsolve:                          ', t40e  - t40s
!  write(20,'(a)')          '# End profile data######################################'

return

end subroutine sp_direct_leaf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine recsolve()
! Recursive solver for Ax=b
! Before call this routine, trunk or leaf solver setting must finished in m_pds.
! Wait right hand side vector b via MPI_RECV.
! Solve equation using leafsolve() or trunksolve().
! Send result to parents and do again.
implicit none
include 'mpif.h'

real(8), dimension(m_pds_ndeg*m_pds_neqns) :: b
logical :: do_calc

! for MPI
integer, dimension(MPI_STATUS_SIZE) :: istatus
integer :: ierr

do
!write(20,*)'recsolve: do_calc flag'; call flush(20)!DEBUG
  call MPI_RECV(do_calc, 1, MPI_LOGICAL, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
  if (.false.==do_calc) then
    if (.true. == istrunk) then
      do_calc=.false.
      call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP1,1,MPI_COMM_WORLD,ierr)
      call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP2,1,MPI_COMM_WORLD,ierr)
    end if
    return
  end if

!write(20,*)'recsolve: waiting right hand side vector'; call flush(20)!DEBUG
!write(20,*)'recsolve: m_pds_ndeg',m_pds_ndeg; call flush(20)!DEBUG
!write(20,*)'recsolve: m_pds_neqns',m_pds_neqns; call flush(20)!DEBUG
!write(20,*)'recsolve: total MPI length',m_pds_ndeg*m_pds_neqns; call flush(20)!DEBUG

  call MPI_RECV(b, m_pds_ndeg*m_pds_neqns, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, istatus, ierr)

  if (.true. == istrunk) then
    call trunksolve(b)
  else if(.true. == isleaf) then
    call leafsolve(b)
  else
    call errtrp('recsolve: never come here')
  end if
  call MPI_SEND(b, m_pds_ndeg*m_pds_neqns, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)
end do

return

end subroutine recsolve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine trunksolve(b)
! Do recursive calculation to solve A x =b0
! Before start this routine, m_pds_d, c1 and c2 must be setted on module m_pds correctly.
!
! Get right hand side vector b0 and divide it to b1, b2, bd.
! Send b1, b2 to child process and receive result of A1 xab1 = b1.
! Update right hand side vector as bd' = bd - C1^t xab1 - C2^t xab2
! Using LDU decomposed D' on module m_pds, solve D'x_d= bd' myself.
! Calc v1=C1 x_d, v2=C2 x_d.
! Send v1, v2 to child process and receive result of A1 w1 = v1.
! Calc x1=-w1+xab1, x2=-w2+xab2. 
! Reorder x1, x2, xd and return it as result
use m_cclsmatrix
implicit none
include 'mpif.h'

! I/O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8), dimension(:), intent(inout) :: b

! internal
real(8), allocatable, dimension(:) :: b1, x1, v1, w1, xab1, vtmp1 
real(8), allocatable, dimension(:) :: b2, x2, v2, w2, xab2, vtmp2
real(8), allocatable, dimension(:) :: bd, xd
type (ccls_matrix) :: c1, c2

integer :: neqns, ndeg, neqns1, neqns2, neqnsd, nd
logical :: do_calc

! misc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,k,l,m,n

! for MPI
integer, dimension(MPI_STATUS_SIZE) :: istatus
integer :: ierr

! start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.false.==m_pds_isready_dc) then
  call errtrp('trunksolve is not ready.')
  stop
end if

ndeg   = m_pds_ndeg
neqns  = m_pds_neqns
neqns1 = m_pds_neqns1
neqns2 = m_pds_neqns2
neqnsd = m_pds_neqnsd
nd     = m_pds_nd

c1 = m_pds_c1
c2 = m_pds_c2

allocate(b1(neqns1*ndeg), x1(neqns1*ndeg), v1(neqns1*ndeg), w1(neqns1*ndeg), xab1(neqns1*ndeg), vtmp1(neqns1*ndeg))
allocate(b2(neqns2*ndeg), x2(neqns2*ndeg), v2(neqns2*ndeg), w2(neqns2*ndeg), xab2(neqns2*ndeg), vtmp2(neqns2*ndeg))
allocate(bd(nd), xd(nd))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Update right hand side vector as bd' = bd - C1^t xab1 - C2^t xab2
!

!write(20,*) 'trunksolve: begin modify bd'; call flush(20)!DEBUG

call divvec(neqns, ndeg, m_pds_part, b, b1, b2, bd)

!write(20,*) 'trunksolve: send right hand side vector to children'; call flush(20)!DEBUG
! send to children and get result for A1 xab1 = b1, A2 xab2 = b2
do_calc=.true.
call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP2,1,MPI_COMM_WORLD,ierr)

call MPI_SEND(b1,   ndeg*neqns1, MPI_REAL8,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(b2,   ndeg*neqns2, MPI_REAL8,  IP2,1,MPI_COMM_WORLD,ierr)

call MPI_RECV(xab1, ndeg*neqns1, MPI_REAL8, IP1, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(xab2, ndeg*neqns2, MPI_REAL8, IP2, 1,MPI_COMM_WORLD, istatus, ierr)

do k=1,nd
  call m_cclsmatrix_getvec(c1, k, vtmp1)
  do i=1,neqns1*ndeg
    bd(k) = bd(k) -(vtmp1(i) * xab1(i))
  end do
  call m_cclsmatrix_getvec(c2, k, vtmp2)
  do i=1,neqns2*ndeg
    bd(k) = bd(k) -(vtmp2(i) * xab2(i))
  end do
end do
!write(20,*) 'trunksolve: end modify bd'; call flush(20)!DEBUG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calc D'x=bd'
! D' is already LDU decomposed.
!

call denssolve(m_pds_d, nd, bd)

xd=bd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calc v1=C1 x_d, v2=C2 x_d.
!

v1=0
do k=1, nd
  call m_cclsmatrix_getvec(c1, k, vtmp1)
  do i=1, neqns1*ndeg
    v1(i) = v1(i) + vtmp1(i)*xd(k)
  end do
end do

v2=0
do k=1, nd
  call m_cclsmatrix_getvec(c2, k, vtmp2)
  do i=1, neqns2*ndeg
    v2(i) = v2(i) + vtmp2(i)*xd(k)
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! calc A1 w1 = v1, a2 w2 = v2
!
do_calc=.true.
call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(do_calc,   1, MPI_LOGICAL,  IP2,1,MPI_COMM_WORLD,ierr)

call MPI_SEND(v1,   ndeg*neqns1, MPI_REAL8,  IP1,1,MPI_COMM_WORLD,ierr)
call MPI_SEND(v2,   ndeg*neqns2, MPI_REAL8,  IP2,1,MPI_COMM_WORLD,ierr)

call MPI_RECV(w1, ndeg*neqns1, MPI_REAL8, IP1, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(w2, ndeg*neqns2, MPI_REAL8, IP2, 1,MPI_COMM_WORLD, istatus, ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Calc x1=-w1+xab1, x2=-w2+xab2. 
!

x1 = -w1 + xab1
x2 = -w2 + xab2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! reorder x1, x2, xd and get finel result (x1, x2, xd) to b
!

call reovec(neqns, ndeg, m_pds_part, b, x1, x2, xd)

return

end subroutine trunksolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine leafsolve(b)

implicit none

real(8), dimension(:), intent(inout) :: b
integer ierr

! start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (.false.==m_pds_isready_v) then
  call errtrp('leafsolve is not ready.')
end if

call nusol0(b, m_pds_v, ierr)
if (ierr .ne. 0) call errtrp('leafsolve: nusol0')
return

end subroutine leafsolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calcCtAC()

use m_cclsmatrix
use m_elap
implicit none
include 'mpif.h'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (ccls_matrix) :: pc                     ! parent C 
real(8), allocatable, dimension(:,:) :: pd   ! parent D' for update
integer :: npd                               ! size of parent dens D'

real(8), allocatable, dimension(:)   :: vtmp, vpc
integer :: ndeg, nndeg
integer :: jcol, iofset, idx
integer :: i,j,k,l

! for MPI
integer, dimension(MPI_STATUS_SIZE) :: istatus
integer :: ierr

! for elaps time
real(8) :: t10s, t10e, t20s, t20e, t30s, t30e

! start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! receive parent C
call elapout('calcCtAC: waiting c matrix via MPI')
call MPI_RECV(pc%neqns, 1, MPI_INTEGER, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(pc%ncol,  1, MPI_INTEGER, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(pc%nttbr, 1, MPI_INTEGER, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(pc%ndeg,  1, MPI_INTEGER, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
ndeg=pc%ndeg
nndeg=ndeg*ndeg
npd=pc%ncol*pc%ndeg
call cclsallocation(pc)

call elapout('calcCtAC: begin receive C matrix')
call ptime(t10s)!ELAP
call MPI_RECV(pc%ia,        pc%ncol+1, MPI_INTEGER, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(pc%ja,         pc%nttbr, MPI_INTEGER, IMP, 1,MPI_COMM_WORLD, istatus, ierr)
call MPI_RECV(pc%val,  pc%nttbr*nndeg, MPI_REAL8,   IMP, 1,MPI_COMM_WORLD, istatus, ierr)
call ptime(t10e)!ELAP
call elapout('calcCtAC: end receive C matrix')

allocate(vtmp(pc%neqns*pc%ndeg))
allocate(vpc(pc%neqns*pc%ndeg))
allocate(pd(npd,npd))

! Calc Ct A^-1 C
! Take care C is belong to parent.
call elapout('calcCtAC: start solve CtAC')
call ptime(t20s)!ELAP
pd=0
do i=1, npd
  call m_cclsmatrix_getvec(pc, i, vtmp)

!write(20,*) 'before solve vtmp',vtmp!DEBUG

  if (.true. == istrunk) then
    call trunksolve(vtmp)
  else if(.true. == isleaf) then
    call leafsolve(vtmp)
  else
    call errtrp('calcCtAC: never come here')
  end if

!write(20,*) 'after solve vtmp',vtmp!DEBUG

  do j=1, npd
    jcol = (j+ndeg-1) / ndeg       ! column number in sparse matrix
    iofset = mod(j+ndeg-1, ndeg)   ! offset in val. 0offset
    do k=pc%ia(jcol),pc%ia(jcol+1)-1
      idx=pc%ja(k)                 ! row number in sparse matrix
      do l=1,ndeg
        pd(j,i)=pd(j,i)+pc%val(ndeg*iofset + l,k) * vtmp(ndeg*(idx-1)+l)
      end do
    end do
  end do
end do

!write(20,*)'pd',pd!DEBUG

call ptime(t20e)!ELAP
call elapout('calcCtAC: end solve CtAC')

call elapout('calcCtAC: start send CtAC')
call ptime(t30s)!ELAP
call MPI_SEND(pd, npd*npd, MPI_REAL8, IMP, 1,MPI_COMM_WORLD, ierr)
call ptime(t30e)!ELAP
deallocate(vtmp,pd,vpc)
call elapout('calcCtAC: end send CtAC')

!  write(20,'(a)')          '# profile calcCtAC ########################################'
!  write(20,*) 'Elaptime (sec)'
!  write(20,'(a,1f15.5)') 'STEP01: Recive C from parent:              ', t10e  - t10s
!  write(20,'(a,1f15.5)') 'STEP02: calc CtAC:                         ', t20e  - t20s
!  write(20,'(a,1f15.5)') 'STEP03: send CtAC to parent:               ', t30e  - t30s
!  write(20,'(a)')          '# End profile data######################################'

return

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mkpart(neqns, nttbr, ndeg, irow, jcol, neqns1, neqns2, neqnsd, part, iperm)
! Make partitioning information for given matrix 
! which specified by neqns, nttbr, ndeg, irow, jcol.
! 
! number of partitioned array a1, a2, d (neqns1, neqns2, neqnsd) and
! partition information array part, iperm will return.

implicit none

integer, intent(in) :: neqns, nttbr, ndeg
integer, dimension(:), intent(in) :: irow
integer, dimension(:), intent(in) :: jcol

integer, intent(out) :: neqns1, neqns2 ! dimension of each divided matrix (sparse)
integer, intent(out) :: neqnsd ! dimension of join matrix D (ndeg is not multiplied)
integer, dimension(:), intent(out) :: part  ! belonging of each element in large matrix
integer, dimension(:), intent(out) :: iperm ! relation of index of large matrix and small matrix

integer :: ipart1(neqns),ipart2(neqns),isp(neqns)
integer :: loc1, loc2, loc3
integer :: i,j,k,l,m,n

! make partitioning array
call bi_part_directive(neqns, nttbr, irow, jcol, neqns1, neqns2, neqnsd)
call get_part_result(neqns1, ipart1, neqns2, ipart2, neqnsd, isp)

do i=1,neqns1
   part(ipart1(i))=1
enddo
do i=1,neqns2
   part(ipart2(i))=2
enddo
do i=1,neqnsd
   part(isp(i))=3
enddo
loc1=0
loc2=0
loc3=0
do i=1,neqns
   if(part(i).eq.1) then
     loc1=loc1+1
     iperm(i)=loc1
   endif
   if(part(i).eq.2) then
     loc2=loc2+1
     iperm(i)=loc2
   endif
   if(part(i).eq.3) then
     loc3=loc3+1
     iperm(i)=loc3
   endif
enddo

return
end subroutine mkpart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine divmat(a0, part, iperm, neqns1, neqns2, neqnsd, ndeg, a1, a2, c1, c2, d)

use m_irjcmatrix
use m_cclsmatrix

implicit none

! Divide given matrix a0 to a1 and a2 
! accortind to partitioning informatin array "part" and "iperm" which set by mkpart()
! this subroutine set following valuables.
!
! a1, a2 : divided small matrix
! c1, c2 : connection region
! d  : connecting dens matrix

! I/O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (irjc_matrix), intent(in) :: a0
integer, dimension(:), intent(in) :: part  ! belonging of each element in large matrix
integer, dimension(:), intent(in) :: iperm ! relation of index of large matrix and small matrix
integer, intent(in) :: neqns1 ! size of a1
integer, intent(in) :: neqns2 ! size of a2
integer, intent(in) :: neqnsd ! size of D (ndeg is not multiplied. so actual size of d is neqnsd*ndeg)
integer, intent(in) :: ndeg   ! degree of freedum of each element in sparse matrix

type (irjc_matrix), intent(out) :: a1, a2
type (ccls_matrix), intent(out) :: c1, c2
real(8), dimension(:,:), intent(out) :: d

! internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, allocatable :: jstat1(:), irpt1(:), irowno1(:)
integer, allocatable :: jstat2(:), irpt2(:), irowno2(:)
integer :: nttbr1, nttbr2, ntt1, ntt2

integer :: ipass, nndeg
integer :: i,j,k,l,m,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(jstat1(neqnsd),    jstat2(neqnsd))
allocate(irpt1(a0%nttbr),   irpt2(a0%nttbr))
allocate(irowno1(a0%nttbr), irowno2(a0%nttbr))

jstat1=0
jstat2=0
irpt1=0
irpt2=0
irowno1=0
irowno2=0
ntt1=0
ntt2=0

nndeg=ndeg*ndeg

!write(20,*) 'divmat: begin divide matrix'; call flush(20)!DEBUG
do ipass=1,2
  nttbr1=0
  nttbr2=0
  do l=1,a0%nttbr
    i=a0%irow(l)
    j=a0%jcol(l)

    if(part(i).eq.1.and.part(j).eq.1) then
      nttbr1=nttbr1+1
      if(ipass.eq.2) then
        a1%irow(nttbr1)=iperm(i)
        a1%jcol(nttbr1)=iperm(j)
        do k=1,nndeg
          a1%val(k,nttbr1)=a0%val(k,l)
        end do
      end if
    end if

    if(part(i).eq.2.and.part(j).eq.2) then
      nttbr2=nttbr2+1
      if(ipass.eq.2) then
        a2%irow(nttbr2)=iperm(i)
        a2%jcol(nttbr2)=iperm(j)
        do k=1,nndeg
          a2%val(k,nttbr2)=a0%val(k,l)
        end do
      end if
    end if

    if(part(i).eq.1.and.part(j).eq.3) then
      if(ipass.eq.1)then
        call reserv(iperm(i),iperm(j),jstat1,irpt1,irowno1,ntt1)
      endif
      if(ipass.eq.2)then
        call stval(c1,iperm(i),iperm(j),a0%val(1,l),0)
      endif
    end if

    if(part(i).eq.2.and.part(j).eq.3) then
      if(ipass.eq.1)then
        call reserv(iperm(i),iperm(j),jstat2,irpt2,irowno2,ntt2)
      endif
      if(ipass.eq.2)then
        call stval(c2,iperm(i),iperm(j),a0%val(1,l),0)
      endif
    end if

    if(part(i).eq.3.and.part(j).eq.1) then
      if(ipass.eq.1)then
        call reserv(iperm(j),iperm(i),jstat1,irpt1,irowno1,ntt1)
      endif
      if(ipass.eq.2)then
        call stval(c1,iperm(j),iperm(i),a0%val(1,l),1)
      endif
    end if

    if(part(i).eq.3.and.part(j).eq.2) then
      if(ipass.eq.1)then
        call reserv(iperm(j),iperm(i),jstat2,irpt2,irowno2,ntt2)
      endif
      if(ipass.eq.2)then
        call stval(c2,iperm(j),iperm(i),a0%val(1,l),1)
      endif
    end if

    if(part(i).eq.3.and.part(j).eq.3) then
      if(ipass.eq.2) then
        do m=1,ndeg
          do k=1,ndeg
            d(ndeg*(iperm(i)-1)+k,ndeg*(iperm(j)-1)+m)=a0%val(ndeg*(m-1)+k,l)
            if(i.ne.j) then
              d(ndeg*(iperm(j)-1)+k,ndeg*(iperm(i)-1)+m)=a0%val(ndeg*(k-1)+m,l)
            end if
          end do
        end do
      end if
    end if

    if(part(i).eq.1.and.part(j).eq.2)then
      write(20,*) "divmat: wrong",i,j
      stop
    endif
    if(part(i).eq.2.and.part(j).eq.1)then
      write(20,*) "divmat: wrong",i,j
      stop
    endif
  end do

  if(ipass.eq.1)then
    a1%neqns=neqns1
    a1%ndeg=ndeg
    a1%nttbr=nttbr1

    a2%neqns=neqns2
    a2%ndeg=ndeg
    a2%nttbr=nttbr2

    allocate(a1%irow(nttbr1), a1%jcol(nttbr1), a1%val(ndeg*ndeg,nttbr1))
    allocate(a2%irow(nttbr2), a2%jcol(nttbr2), a2%val(ndeg*ndeg,nttbr2))

    c1%neqns=neqns1
    c1%ncol=neqnsd
    c1%nttbr=ntt1
    c1%ndeg=ndeg
    allocate(c1%ia(c1%ncol+1))
    allocate(c1%ja(c1%nttbr))
    allocate(c1%val(ndeg*ndeg,c1%nttbr))

    c2%neqns=neqns2
    c2%ncol=neqnsd
    c2%nttbr=ntt2
    c2%ndeg=ndeg
    allocate(c2%ia(c2%ncol+1))
    allocate(c2%ja(c2%nttbr))
    allocate(c2%val(ndeg*ndeg,c2%nttbr))


!    call cclsallocation(c1)
!    call cclsallocation(c2)
  
    call stiajac(c1,jstat1,irpt1,irowno1)
    call stiajac(c2,jstat2,irpt2,irowno2)
    deallocate(jstat1,irpt1,irowno1)
    deallocate(jstat2,irpt2,irowno2)
  endif
end do


!write(20,*) 'divmat: end divide matrix'; call flush(20)!DEBUG

return

end subroutine divmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine divvec(neqns, ndeg, part, r, r1, r2, rd)

implicit none

integer, intent(in) :: neqns, ndeg
real(8), dimension(:), intent(in)  :: r    ! original right hand side vector
integer, dimension(:), intent(in)  :: part ! belonging of each element in large matrix

real(8), dimension(:), intent(out) :: r1
real(8), dimension(:), intent(out) :: r2
real(8), dimension(:), intent(out) :: rd

integer :: idx1, idx2, idxd, idxr
integer :: i

idx1=1
idx2=1
idxd=1
idxr=1
do i=1,neqns
   if(part(i).eq.1) then
     r1(idx1:idx1+ndeg-1)=r(idxr:idxr+ndeg-1)
     idx1=idx1+ndeg
     idxr=idxr+ndeg
   endif
   if(part(i).eq.2) then
     r2(idx2:idx2+ndeg-1)=r(idxr:idxr+ndeg-1)
     idx2=idx2+ndeg
     idxr=idxr+ndeg
   endif
   if(part(i).eq.3) then
     rd(idxd:idxd+ndeg-1)=r(idxr:idxr+ndeg-1)
     idxd=idxd+ndeg
     idxr=idxr+ndeg
   endif
end do

end subroutine divvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine reovec(neqns, ndeg, part, r, r1, r2, rd)

implicit none

integer, intent(in) :: neqns, ndeg
integer, dimension(:),       intent(in)  :: part ! belonging of each element in large matrix
real(8), dimension(:), intent(in)  :: r1
real(8), dimension(:), intent(in)  :: r2
real(8), dimension(:), intent(in)  :: rd

real(8), dimension(:), intent(out) :: r    ! reordered right hand side vector

integer :: idx1, idx2, idxd, idxr
integer :: i

idx1=1
idx2=1
idxd=1
idxr=1

do i=1,neqns
   if(part(i).eq.1) then
     r(idxr:idxr+ndeg-1)=r1(idx1:idx1+ndeg-1)
     idx1=idx1+ndeg
     idxr=idxr+ndeg
   endif
   if(part(i).eq.2) then
     r(idxr:idxr+ndeg-1)=r2(idx2:idx2+ndeg-1)
     idx2=idx2+ndeg
     idxr=idxr+ndeg
   endif
   if(part(i).eq.3) then
     r(idxr:idxr+ndeg-1)=rd(idxd:idxd+ndeg-1)
     idxd=idxd+ndeg
     idxr=idxr+ndeg
   endif
enddo

end subroutine reovec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine densldu(a,n)
! do LDU decomposition for given matrix a.
! after this routine a is returned as below:
!
! DUUUUU
! LDUUUU
! LLDUUU
! LLLDUU
! LLLLDU
! LLLLLD
!
! in nortation of
!
! a = l*d*u
!
! l = 1000000 d = D000000 u = 1UUUUUU
!     L100000     0D00000     01UUUUU
!     LL10000     00D0000     001UUUU
!     LLL1000     000D000     0001UUU
!     LLLL100     0000D00     00001UU
!     LLLLL10     00000D0     000001U
!     LLLLLL1     000000D     0000001

!$Id: m_pds.f90,v 1.2 2007/11/21 02:13:45 kitayama Exp $
implicit none
integer :: i,j,k,l,m,n
real(8), dimension(n,n) :: a


!write(*,*) 'in ldu: a'
!do i=1,n
!  do j=1,n
!write(*,*) a(i,j)
!  end do
!end do

do j=2, n
  a(1,j) = a(1,j)/a(1,1)
end do

do i=2,n
  do k=1,i-1
    a(i,i) = a(i,i) - a(i,k)*a(k,i)
  end do

  do j=i+1,n
    do k=1,i-1
      a(j,i)=a(j,i)-a(j,k)*a(k,i)
    end do

    do k=1,i-1
      a(i,j)=a(i,j)-a(i,k)*a(k,j)
    end do
    a(i,j) = a(i,j)/a(i,i)
  end do
end do

! treat a as L D U
do j=1,n
  do i=j+1, n
    a(i,j) = a(i,j)/a(j,j)
  end do
end do

!write(*,*) "in ldu: LDU decomposed a"
!do i=1,n
!  do j=1,n
!write(*,*) a(i,j)
!  end do
!end do


end subroutine densldu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine denssolve(a,ndim,r)
!$Id: m_pds.f90,v 1.2 2007/11/21 02:13:45 kitayama Exp $
! solve Ax=b for multiple right hand side vector r=b1, b2, b3...
! result is return as r
integer :: ndim ! dimension of matrix a
real(8), dimension(:,:) :: a ! coffecient matrix which already decomposed as LDU by ldu()
real(8), dimension(:)   :: r   ! right hand side vectors

!is it needed?
!real(8), dimension(ndim, nrhs) :: x !TMP

integer :: i,j,k,l,m,n

! solve Ly=b (y=DUx)
do i=1,ndim
  do j=1,i-1
    r(i)=r(i)-a(i,j)*r(j)
  end do
end do

!solve Dz=y (z=Ux)
do i=1, ndim
  r(i)=r(i)/a(i,i)
end do

!solve Ux=z
do i=ndim-1, 1, -1
  do j=i+1, ndim
    r(i)=r(i)-a(i,j)*r(j)
  end do
end do

end subroutine denssolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setbt()
implicit none

include 'mpif.h'

integer :: ierr

isroot=.false.
istrunk=.false.
isleaf=.false.

! get process number and number of whole processes
call MPI_COMM_SIZE(MPI_COMM_WORLD, npe, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

if ((mod(npe,2)==0) .and. (myid == npe-1)) then
!  write(20,*) 'I am not a member of tree.'
  call fin()
  stop
end if

if (myid==0) then
  isroot=.true.
  ip1=myid*2+1
  ip2=myid*2+2
!  write(20,*)'I am root.'
else
  imp=(myid-1)/2
end if

! trunk or leaf
if (myid*2+1 < npe-1) then
  istrunk=.true.
  ip1=myid*2+1
  ip2=myid*2+2
!  write(20,*)'I am trunk.'
else
  isleaf=.true.
!  write(20,*)'I am leaf.'
end if
return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine errtrp(mes)
character(*) mes
write(6,*)  'Error in : process ', myid
write(20,*) mes

call finalize()
stop
end subroutine errtrp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_pds
