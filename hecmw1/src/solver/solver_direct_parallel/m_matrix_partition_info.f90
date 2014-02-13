!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2008/03/13                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takeshi Kitayama (Univ. of Tokyo)              !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
module m_matrix_partition_info

use hecmw_util
use m_irjc_matrix
use m_child_matrix

private !default access control

public matrix_partition_info
public matrix_partition_recursive_bisection
public reovec

! Handle matrix partitioning informations.
type matrix_partition_info
  integer(kind=kint)          :: ndm          ! number of dm
  type(child_matrix), pointer :: dm(:)        ! partitioned matrix. 2^n number
  integer(kind=kint)          :: neqns_d      ! number of eqns in D matrix
  real(kind=kreal), pointer   :: dsln(:,:)    ! non-diagonal elements of dens D matrix which kept in proc0
  real(kind=kreal), pointer   :: diag(:,:)    ! diagonal elements of dens D matrix
  integer(kind=kint), pointer :: part_all(:)  ! index of corresponding dm of a0 row
  integer(kind=kint), pointer :: iperm_all(:) ! index in partitioned matrix of a0 row
end type

! node of bisection partition tree
type matrix_partition_node
  integer(kind=kint) :: id    ! position in array
  integer(kind=kint) :: depth ! depth in tree. root is 0.
  integer(kind=kint) :: idm   ! corresponding divided matrix no. root and tree node have 0. leaf node have 1..2^n

  integer(kind=kint) :: neqns_a, nttbr_a ! size and number of non-zero element in a0.
  integer(kind=kint) :: neqns_d ! size of D. 
  integer(kind=kint), pointer :: irow(:), jcol(:)
  integer(kind=kint), pointer :: idx_g_a(:), idx_g_d(:) ! global indices of each point
end type

integer(kind=kint) :: maxdepth   ! root is 0.
integer(kind=kint) :: iend_trunk ! position of last trunk node.
integer(kind=kint) :: ista_leaf  ! position of first leaf node.
integer(kind=kint) :: iend_leaf  ! position of last leaf node.

type(matrix_partition_node), pointer :: tree_array(:)

integer(kind=kint), parameter :: ilog = 16

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_partition_recursive_bisection(a0, ndiv, pmi)
implicit none
type(irjc_square_matrix), intent(inout) :: a0
integer(kind=kint), intent(in) :: ndiv
type(matrix_partition_info), intent(out) :: pmi

type(matrix_partition_node), pointer :: mdroot ! root of partition tree

integer(kind=kint) :: i, ndeg, nndeg, ndegt


! make matrix partition tree using METIS library.
! leaf of tree is according to divided child matrix.
! trunk of tree is according to D region (parent process).
write(6,*)  'METIS 4.0 library is used for matrix partitioning.'
write(6,*)  'Copyright 1997, Regents of the University of Minnesota'
write(6,*)  'http://glaros.dtc.umn.edu/gkhome/metis/metis/overview'
pmi%ndm=2**ndiv
call matrix_partition_tree_initialize(ndiv)
call matrix_partition_tree_get_root(mdroot)
mdroot%neqns_a=a0%neqns
mdroot%nttbr_a=a0%nttbr
allocate(mdroot%idx_g_a(a0%neqns))
do i=1,a0%neqns
  mdroot%idx_g_a(i)=i
end do
allocate(mdroot%irow(size(a0%irow)))
allocate(mdroot%jcol(size(a0%jcol)))
mdroot%irow=a0%irow
mdroot%jcol=a0%jcol
call make_matrix_partition_tree(mdroot)

! get global partition array from matrix partition tree.
allocate(pmi%part_all(a0%neqns), pmi%iperm_all(a0%neqns))
call mkpart(pmi%part_all, pmi%iperm_all, pmi%neqns_d)

ndeg=a0%ndeg
nndeg = ndeg*ndeg
ndegt = (ndeg+1)*ndeg/2  ! triangle number. 6 for ndeg=3, 3 for ndeg=2
allocate(pmi%dsln(nndeg, pmi%neqns_d*(pmi%neqns_d - 1) / 2))
allocate(pmi%diag(ndegt, pmi%neqns_d))
allocate(pmi%dm(2**ndiv))
pmi%dsln=0
pmi%diag=0
call divmat(a0, pmi%part_all, pmi%iperm_all, pmi%dm, pmi%dsln, pmi%diag, pmi%neqns_d)

end subroutine matrix_partition_recursive_bisection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_partition_tree_initialize(ndiv)
! allocate 2^ndiv matrix_partition_node and set belongings in processes for each node.

implicit none
integer(kind=kint), intent(in) :: ndiv ! number of bisection

type(matrix_partition_node), pointer :: this

integer(kind=kint) :: tsize, i


! allocate matrix divide tree
tsize=1
do i=1,ndiv
  tsize=tsize+2**i
end do
allocate(tree_array(tsize))
iend_trunk = tsize - 2**ndiv
ista_leaf  = iend_trunk+1
iend_leaf  = tsize
maxdepth = ndiv
do i=1,tsize
  if (i .le. iend_trunk) then
    this=>tree_array(i)
    this%id=i
    this%idm=0
    this%depth=depthInBinaryTree(this%id)
  else if (i .ge. ista_leaf) then
    this=>tree_array(i)
    this%id=i
    this%idm=i - iend_trunk
    this%depth=maxdepth
  end if
end do
end subroutine matrix_partition_tree_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_partition_tree_get_root(root)
implicit none
type(matrix_partition_node),pointer :: root
root=>tree_array(1)
end subroutine matrix_partition_tree_get_root

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_partition_tree_get_left_child(this, left)
implicit none
type(matrix_partition_node), pointer :: this
type(matrix_partition_node), pointer :: left
integer(kind=kint) :: id_this, id_left
id_this = this%id
id_left = id_this*2
left => tree_array(id_left)
end subroutine matrix_partition_tree_get_left_child

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_partition_tree_get_right_child(this, right)
implicit none
type(matrix_partition_node), pointer :: this
type(matrix_partition_node), pointer :: right
integer(kind=kint) :: id_this, id_right
id_this = this%id
id_right = id_this*2+1
right => tree_array(id_right)
end subroutine matrix_partition_tree_get_right_child

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine make_matrix_partition_tree(this)
! make matrix partition tree according to graph structure of given matrix via matrix_partition_node.
implicit none
type(matrix_partition_node), pointer :: this, left, right
integer(kind=kint), allocatable :: ip_left(:), ip_right(:), ip_d(:), part(:), iperm(:)
integer(kind=kint) :: loc1, loc2, loc3, ipass
integer(kind=kint) :: i,j,k,l

allocate(ip_left(this%neqns_a), ip_right(this%neqns_a), ip_d(this%neqns_a)) ! local permtation index used by metis.
allocate(part(this%neqns_a), iperm(this%neqns_a))


if (this%depth .eq. maxdepth) then ! leaf node
  return
end if

call matrix_partition_tree_get_left_child(this, left)
call matrix_partition_tree_get_right_child(this, right)

call bi_part_directive(this%neqns_a, this%nttbr_a, this%irow, this%jcol, left%neqns_a, right%neqns_a, this%neqns_d)

if ((left%neqns_a + right%neqns_a + this%neqns_d) .ne. this%neqns_a) then
  write(6,*) 'Error: Fail to matrix partitioning in parallel direct solver.'
  write(6,*) 'Matrix is too small.'
  write(ilog,*) 'Error: Fail to matrix partitioning in parallel direct solver.'
  write(ilog,*) 'Matrix is too small.'
  call hecmw_abort( hecmw_comm_get_comm())
end if

call get_part_result(left%neqns_a, ip_left, right%neqns_a, ip_right, this%neqns_d, ip_d)

! set belongings
do i=1,left%neqns_a  ! left child
  part(ip_left(i))=1
end do
do i=1,right%neqns_a ! right child
  part(ip_right(i))=2
end do
do i=1,this%neqns_d  ! remain in this node
  part(ip_d(i))=3
end do

! set global indices of each point in terms of local array indices.
allocate(left%idx_g_a(left%neqns_a))
allocate(right%idx_g_a(right%neqns_a))
allocate(this%idx_g_d(this%neqns_d))

loc1=0
loc2=0
loc3=0
do i=1,this%neqns_a
  if(part(i) .eq. 1) then
    loc1=loc1+1
    iperm(i)=loc1
    left%idx_g_a(loc1)=this%idx_g_a(i)
  end if
  if(part(i) .eq. 2) then
    loc2=loc2+1
    iperm(i)=loc2
    right%idx_g_a(loc2)=this%idx_g_a(i)
  end if
  if(part(i) .eq. 3) then
    loc3=loc3+1
    iperm(i)=loc3
    this%idx_g_d(loc3)=this%idx_g_a(i)
  end if
end do

! devide a0 to a_left and a_right
do ipass=1,2
  left%nttbr_a  = 0
  right%nttbr_a = 0

  do l=1,this%nttbr_a
    i=this%irow(l)
    j=this%jcol(l)

    if ((part(i) .eq. 1) .and. (part(j) .eq. 1)) then !left
      left%nttbr_a=left%nttbr_a + 1
      if(ipass .eq. 2) then
        left%irow(left%nttbr_a)=iperm(i)
        left%jcol(left%nttbr_a)=iperm(j)
      end if
      cycle
    end if

    if ((part(i) .eq. 2) .and. (part(j) .eq. 2)) then !right
      right%nttbr_a=right%nttbr_a + 1
      if(ipass .eq. 2) then
        right%irow(right%nttbr_a)=iperm(i)
        right%jcol(right%nttbr_a)=iperm(j)
      end if
      cycle
    end if
  end do

  if (ipass .eq. 1) then
    allocate(left%irow(left%nttbr_a),   left%jcol(left%nttbr_a))
    allocate(right%irow(right%nttbr_a), right%jcol(right%nttbr_a))
  end if
end do

deallocate(ip_left, ip_right, ip_d, iperm, part, this%irow, this%jcol)

call make_matrix_partition_tree(left)
call make_matrix_partition_tree(right)
return
end subroutine make_matrix_partition_tree

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mkpart(part_g, iperm_g, neqns_d)
! make global partition array from matrix partition tree.
implicit none
integer(kind=kint), intent(out) :: part_g(:)
integer(kind=kint), intent(out) :: iperm_g(:)
integer(kind=kint), intent(out) :: neqns_d

type(matrix_partition_node), pointer :: node
integer(kind=kint) :: ncount
integer(kind=kint) :: i,j



! for leaf
do i=iend_leaf, ista_leaf, -1
  node=>tree_array(i)
  ncount=0         ! node local numbering
  do j=1,node%neqns_a
    ncount=ncount+1
    part_g(node%idx_g_a(j))=node%idm ! index of divided matrix. 1..2**ndiv
    iperm_g(node%idx_g_a(j))=ncount  ! divided matrix local numbering
  end do
  deallocate(node%idx_g_a)
end do

! for trunk
ncount=0 !global numbering
do i=iend_trunk, 1, -1
  node => tree_array(i)
  do j=1,node%neqns_d
    ncount=ncount+1
    part_g(node%idx_g_d(j))=node%idm ! parent D. 0
    iperm_g(node%idx_g_d(j))=ncount  ! global numbering
  end do
  deallocate(node%idx_g_d)
end do

neqns_d = ncount ! size of D
deallocate(tree_array) ! tree_array is not used after here.
return
end subroutine mkpart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine divmat(a0, part_g, iperm_g, dm, dsln, diag, neqns_d)
implicit none
! divide a0 to n children and D region on parent, according to the global permtation information.

type(irjc_square_matrix), intent(in) :: a0
integer(kind=kint), intent(in) :: part_g(:), iperm_g(:)

type(child_matrix), intent(out), target :: dm(:)
real(kind=kreal), intent(out) :: dsln(:,:), diag(:,:) ! for dens D

integer(kind=kint), intent(in) :: neqns_d

! internal
type(child_matrix), pointer :: dmc
integer(kind=kint) :: ndeg, ipass, ipos, itmp
integer(kind=kint) :: i,j,k,l,m,n, ii, jj

! set neqns for each divided matrix 
ndeg = a0%ndeg
do i=1,size(dm)
  dm(i)%a%neqns = 0
end do

do i=1,a0%neqns
  k=part_g(i)
  if (k .ne. 0 ) then
    dm(k)%a%neqns = dm(k)%a%neqns + 1
  end if
end do

do i=1,size(dm)
  dm(i)%ndeg    = ndeg
  dm(i)%a%ndeg  = ndeg
  dm(i)%c%ndeg  = ndeg
  dm(i)%c%nrows = neqns_d
  dm(i)%c%ncols = dm(i)%a%neqns
  dm(i)%ista_c  = dm(i)%a%neqns + 1
end do

! divide matrix element to child matrix
do ipass=1,2
  do i=1,size(dm)
    dm(i)%a%nttbr = 0
    dm(i)%c%nttbr = 0
  end do

  do l=1,a0%nttbr
    i=a0%irow(l)
    j=a0%jcol(l)
  
    if (part_g(i) .eq. part_g(j)) then ! same matrix
      if ((part_g(i) .eq. 0) .or. (part_g(j) .eq. 0)) then       ! D
        if (iperm_g(i) .eq. iperm_g(j)) then                     ! diag in D
          ipos=1
          do m=1,ndeg
            do n=1,m
              diag(ipos, iperm_g(i))=a0%val(ndeg*(m-1)+n,l)
              ipos=ipos+1
            end do
          end do
  
        else if (iperm_g(i) .gt. iperm_g(j)) then                ! dsln in D
  
          !count index in dsln for specified i,j.
          ii = iperm_g(i)
          jj = iperm_g(j)
  
          k = (ii-1)*(ii-2)/2 + jj
  
          dsln(:,k)=a0%val(:,l)
  
        else                                                     ! dsln in D with inverse
          !count index in dsln for specified i,j.
          jj = iperm_g(i)
          ii = iperm_g(j)
  
          k = (ii-1)*(ii-2)/2 + jj
  
          do m=1,ndeg
            do n=1,ndeg
              dsln((m-1)*ndeg+n,k)=a0%val(m+(n-1)*ndeg,l)
            end do
          end do

        end if
      else !An
        dmc=>dm(part_g(i))
        dmc%a%nttbr = dmc%a%nttbr + 1
        if(ipass.eq.2) then
          dmc%a%irow(dmc%a%nttbr)  = iperm_g(i)
          dmc%a%jcol(dmc%a%nttbr)  = iperm_g(j)
          dmc%a%val(:,dmc%a%nttbr) = a0%val(:,l)
        end if
      end if
      cycle
    end if
  
    if (part_g(i) .eq. 0) then !C run right in ndm part_g(j)
      dmc=>dm(part_g(j))
      dmc%c%nttbr = dmc%c%nttbr + 1
      if(ipass.eq.2) then
        dmc%c%irow(dmc%c%nttbr)=iperm_g(i)
        dmc%c%jcol(dmc%c%nttbr)=iperm_g(j)
        dmc%c%val(:,dmc%c%nttbr)=a0%val(:,l)
      end if
      cycle
    end if
  
    if (part_g(j) .eq. 0) then !C run down in ndm part_g(i). Invert val
      dmc=>dm(part_g(i))
      dmc%c%nttbr = dmc%c%nttbr + 1
      if(ipass.eq.2) then
        dmc%c%irow(dmc%c%nttbr)=iperm_g(j)
        dmc%c%jcol(dmc%c%nttbr)=iperm_g(i)
        do m=1,ndeg
          do n=1,ndeg
            dmc%c%val(n+ndeg*(m-1),dmc%c%nttbr)=a0%val(ndeg*(n-1)+m,l)
          end do
        end do
      end if
      cycle
    end if
  
    stop !never come here
  end do
  
  if (ipass .eq. 1) then
    do i=1,size(dm)
      allocate(dm(i)%a%irow(dm(i)%a%nttbr))
      allocate(dm(i)%a%jcol(dm(i)%a%nttbr))
      allocate(dm(i)%a%val(ndeg*ndeg,dm(i)%a%nttbr))

      allocate(dm(i)%c%irow(dm(i)%c%nttbr))
      allocate(dm(i)%c%jcol(dm(i)%c%nttbr))
      allocate(dm(i)%c%val(ndeg*ndeg,dm(i)%c%nttbr))
    end do
  end if
end do

return
end subroutine divmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine reovec(r, iperm)
! reorder vector as given iperm

implicit none

integer(kind=kint), dimension(:),   intent(in)    :: iperm ! permtation vector
real(kind=kreal), dimension(:,:), intent(inout) :: r     ! reordered result will be return

real(kind=kreal), allocatable, dimension(:,:) :: tmp
integer(kind=kint) :: ndeg, neqns
integer(kind=kint), dimension(2) :: ishape
integer(kind=kint) :: i,j,k,l

ishape = shape(r)
ndeg   = ishape(1)
neqns  = ishape(2)
allocate(tmp(ndeg, neqns))


do i=1, neqns
  do j=1, ndeg
    tmp(j,iperm(i))=r(j,i)
  end do
end do

r=tmp

deallocate(tmp)

return
end subroutine reovec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function depthInBinaryTree(n)
! return depth of given integer "n" in 1-ofset binary tree. root have 0 depth.
implicit none
integer(kind=kint) :: depthInBinaryTree
integer(kind=kint), intent(in)  :: n
integer(kind=kint) :: i
i=0
do while (n .ge. 2**(i+1))
  i=i+1
end do
depthInBinaryTree = i ! root is 0
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_matrix_partition_info
