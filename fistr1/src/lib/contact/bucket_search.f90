!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides bucket-search functionality
!>  It provides definition of bucket info and its access routines
module bucket_search
  use hecmw
  implicit none

  private
  public :: bucketDB
  public :: bucketDB_init
  public :: bucketDB_finalize
  public :: bucketDB_setup
  public :: bucketDB_getBucketID
  public :: bucketDB_registerPre
  public :: bucketDB_allocate
  public :: bucketDB_register
  public :: bucketDB_getNumCand
  public :: bucketDB_getCand

  integer(kind=kint), parameter :: DEBUG = 0

  !> Structure for a single bucket (private to this module)
  type bucket
    integer(kind=kint) :: n                             !< number of members stored in this bucket
    integer(kind=kint) :: n_max                         !< number of maximum members that can be stored in this bucket
    integer(kind=kint), pointer :: member(:) => null()  !< array pointer of members
  end type bucket

  !> Structure for bucket search
  type bucketDB
    private
    real(kind=kreal) :: x_min(3)                           !< min coordinate of rectangular area
    real(kind=kreal) :: x_max(3)                           !< max coordinate of rectangular area
    real(kind=kreal) :: d(3)                               !< bucket sizes in x,y,z direction
    integer(kind=kint) :: ndiv(3)                          !< divisions in x,y,z direction
    type(bucket), pointer :: buckets(:,:,:) => null()      !< 3D array of buckets
    integer(kind=kint) :: n_tot                            !< total number of members
    integer(kind=kint), pointer :: member_all(:) => null() !< array of members (referenced by each bucket)
  end type bucketDB

contains

  !> Assertion routine for debugging
  subroutine assert(cond, mesg)
    implicit none
    logical, intent(in) :: cond  !< condition statement that should be true
    character(len=*) :: mesg     !< error message when the condition is false
    if (DEBUG > 0) then
      if (.not. cond) then
        write(0,*) 'ASSERTION FAILED: ',mesg
        call hecmw_abort( hecmw_comm_get_comm() )
      endif
    endif
  end subroutine assert

!!!
!!! routines for type(bucket)
!!!

  !> Initializer
  subroutine bucket_init(bkt)
    implicit none
    type(bucket), intent(inout) :: bkt  !< bucket
    bkt%n = 0
    bkt%n_max = 0
    nullify(bkt%member)
  end subroutine bucket_init

  !> Finalizer
  subroutine bucket_finalize(bkt)
    implicit none
    type(bucket), intent(inout) :: bkt  !< bucket
    !if (bkt%n > 0) deallocate(bkt%member)
    nullify(bkt%member)
    bkt%n = 0
    bkt%n_max = 0
  end subroutine bucket_finalize

  !> Just increment count before actual registration of member
  subroutine bucket_incr_count(bkt)
    implicit none
    type(bucket), intent(inout) :: bkt  !< bucket
    !$omp atomic
    bkt%n = bkt%n + 1
  end subroutine bucket_incr_count

  !> Assign memory to member array pointer
  subroutine bucket_assign(bkt, mem)
    implicit none
    type(bucket), intent(inout) :: bkt     !< bucket
    integer(kind=kint), pointer :: mem(:)  !< pointer to integer array
    bkt%member => mem
    bkt%n_max = bkt%n
    bkt%n = 0
  end subroutine bucket_assign

  !> Register member
  subroutine bucket_register(bkt, sid)
    implicit none
    type(bucket), intent(inout) :: bkt     !< bucket
    integer(kind=kint), intent(in) :: sid  !< member ID
    integer(kind=kint) :: idx
    !$omp atomic capture
    bkt%n = bkt%n + 1
    idx = bkt%n
    !$omp end atomic
    call assert(idx <= bkt%n_max, 'bucket_register: too many members')
    bkt%member(idx) = sid
  end subroutine bucket_register

  !> Get number of members
  function bucket_get_n(bkt)
    implicit none
    integer(kind=kint) :: bucket_get_n  !< number of members
    type(bucket), intent(in) :: bkt     !< bucket
    bucket_get_n = bkt%n
  end function bucket_get_n

  !> Get members
  subroutine bucket_get_member(bkt, n, memb)
    implicit none
    type(bucket), intent(in) :: bkt             !< bucket
    integer(kind=kint), intent(in) :: n         !< number of members
    integer(kind=kint), intent(out) :: memb(n)  !< array to store members
    call assert(n == bkt%n, 'bucket_get_member: wrong n')
    memb(1:n) = bkt%member(1:n)
  end subroutine bucket_get_member

!!!
!!! routines for type(bucketDB)
!!!

  !> Initializer
  subroutine bucketDB_init(bktdb)
    implicit none
    type(bucketDB), intent(inout) :: bktdb  !< bucket info
    bktdb%x_min(:) = 0.d0
    bktdb%x_max(:) = 0.d0
    bktdb%d(:) = 0.d0
    bktdb%ndiv(:) = 0
    nullify(bktdb%buckets)
    bktdb%n_tot = 0
    nullify(bktdb%member_all)
  end subroutine bucketDB_init

  !> Finalizer
  subroutine bucketDB_finalize(bktdb)
    implicit none
    type(bucketDB), intent(inout) :: bktdb  !< bucket info
    integer(kind=kint) :: i, j, k
    if (bktdb%n_tot > 0) then
      deallocate(bktdb%member_all)
      bktdb%n_tot = 0
    endif
    if (any(bktdb%ndiv == 0)) then
      bktdb%ndiv(:) = 0
      return
    endif
    do k = 1, bktdb%ndiv(3)
      do j = 1, bktdb%ndiv(2)
        do i = 1, bktdb%ndiv(1)
          call bucket_finalize(bktdb%buckets(i,j,k))
        enddo
      enddo
    enddo
    deallocate(bktdb%buckets)
    bktdb%ndiv(:) = 0
  end subroutine bucketDB_finalize

  !> Setup basic info of buckets
  subroutine bucketDB_setup(bktdb, x_min, x_max, dmin, n_tot)
    implicit none
    type(bucketDB), intent(inout) :: bktdb    !< bucket info
    real(kind=kreal), intent(in) :: x_min(3)  !< min coordinate of rectangle area covered by buckets
    real(kind=kreal), intent(in) :: x_max(3)  !< max coordinate of rectangle area covered by buckets
    real(kind=kreal), intent(in) :: dmin      !< minimal size of bucket
    integer(kind=kint), intent(in) :: n_tot   !< total number of members to be stored
    real(kind=kreal) :: xrange(3)
    integer(kind=kint) :: i, j, k
    real(kind=kreal), parameter :: EPS = 1.d-6
    if (DEBUG >= 1) write(0,*) 'DEBUG: bucketDB_setup', x_min, x_max, dmin, n_tot
    if (associated(bktdb%buckets)) deallocate(bktdb%buckets)
    bktdb%x_min(:) = x_min(:)
    bktdb%x_max(:) = x_max(:)
    xrange(:) = x_max(:) - x_min(:)
    call assert(all(xrange > 0.d0), 'bucketDB_setup: invalid x_min, x_max')
    do i = 1, 3
      bktdb%ndiv(i) = max(floor(xrange(i) / dmin), 1)
      bktdb%d(i) = xrange(i) / bktdb%ndiv(i) * (1.d0 + EPS)
    enddo
    if (DEBUG >= 1) write(0,*) 'DEBUG: bucketDB_setup: ndiv, d: ', bktdb%ndiv, bktdb%d
    call assert(all(bktdb%d > 0.d0), 'bucketDB_setup: invalid bktdb%d')
    allocate(bktdb%buckets(bktdb%ndiv(1), bktdb%ndiv(2), bktdb%ndiv(3)))
    do k = 1, bktdb%ndiv(3)
      do j = 1, bktdb%ndiv(2)
        do i = 1, bktdb%ndiv(1)
          call bucket_init(bktdb%buckets(i,j,k))
        enddo
      enddo
    enddo
    if (bktdb%n_tot /= n_tot) then
      if (associated(bktdb%member_all)) deallocate(bktdb%member_all)
      allocate(bktdb%member_all(n_tot))
      bktdb%n_tot = n_tot
    endif
  end subroutine bucketDB_setup

  !> Encode 3D address of bucket into integer ID
  function encode_bid(bktdb, baddr)
    implicit none
    integer(kind=kint) :: encode_bid            !< bucket ID
    type(bucketDB), intent(in) :: bktdb         !< bucket info
    integer(kind=kint), intent(in) :: baddr(3)  !< 3D bucket address of bucket
    if (any(baddr <= 0) .or. any(baddr > bktdb%ndiv)) then
      encode_bid = -1
    else
      encode_bid = &
           (baddr(3)-1) * bktdb%ndiv(1) * bktdb%ndiv(2) + (baddr(2)-1) * bktdb%ndiv(1) + baddr(1)
    endif
  end function encode_bid

  !> Decode integer ID of bucket into 3D address
  function decode_bid(bktdb, bid)
    implicit none
    integer(kind=kint) :: decode_bid(3)    !< 3D address of bucket
    type(bucketDB), intent(in) :: bktdb    !< bucket info
    integer(kind=kint), intent(in) :: bid  !< bucket ID
    call assert(bid <= bktdb%ndiv(1)*bktdb%ndiv(2)*bktdb%ndiv(3), 'decode_bid: out of range')
    if (bid < 0) then
      decode_bid(:) = -1
    else
      decode_bid(1) = mod(bid-1, bktdb%ndiv(1)) + 1
      decode_bid(2) = mod((bid-1)/bktdb%ndiv(1), bktdb%ndiv(2)) + 1
      decode_bid(3) = (bid-1)/(bktdb%ndiv(1) * bktdb%ndiv(2)) + 1
      call assert(encode_bid(bktdb, decode_bid) == bid, 'decode_bid')
    endif
  end function decode_bid

  !> Get bucket ID that includes given point
  function bucketDB_getBucketID(bktdb, x)
    implicit none
    integer(kind=kint) :: bucketDB_getBucketID  !< bucket ID
    type(bucketDB), intent(in) :: bktdb         !< bucket info
    real(kind=kreal), intent(in) :: x(3)        !< coordinate of point
    integer(kind=kint) :: baddr(3)
    integer(kind=kint) :: i
    if (bktdb%n_tot == 0) then
      bucketDB_getBucketID = -1
      return
    endif
    do i = 1, 3
      call assert(bktdb%d(i) > 0.d0, 'bucketDB_getBucketID: bktdb%d(i) is zero')
      baddr(i) = floor((x(i) - bktdb%x_min(i)) / bktdb%d(i)) + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_getBucketID: ',x,baddr
    bucketDB_getBucketID = encode_bid(bktdb, baddr)
  end function bucketDB_getBucketID

  !> Pre-register for just counting members to be actually registered
  !! Bucket ID has to be obtained with bucketDB_getBucketID
  subroutine bucketDB_registerPre(bktdb, bid)
    implicit none
    type(bucketDB), intent(inout) :: bktdb  !< bucket info
    integer(kind=kint), intent(in) :: bid   !< bucket ID
    integer(kind=kint) :: baddr(3)
    baddr = decode_bid(bktdb, bid)
    call assert(all(baddr > 0) .and. all(baddr <= bktdb%ndiv), 'bucketDB_register_pre: block ID out of range')
    call bucket_incr_count(bktdb%buckets(baddr(1),baddr(2),baddr(3)))
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_registerPre: ', baddr
  end subroutine bucketDB_registerPre

  !> Allocate memory before actually registering members
  !! Before allocating memory, bucketDB_registerPre has to be called for all members to be registered
  subroutine bucketDB_allocate(bktdb)
    implicit none
    type(bucketDB), intent(inout) :: bktdb  !< bucket info
    integer(kind=kint) :: i, j, k, count, n
    integer(kind=kint), pointer :: pmemb(:)
    count = 0
    do k = 1, bktdb%ndiv(3)
      do j = 1, bktdb%ndiv(2)
        do i = 1, bktdb%ndiv(1)
          !call bucket_allocate(bktdb%buckets(i,j,k))
          n = bucket_get_n(bktdb%buckets(i,j,k))
          pmemb => bktdb%member_all(count+1:count+n)
          call bucket_assign(bktdb%buckets(i,j,k), pmemb)
          count = count + n
        enddo
      enddo
    enddo
  end subroutine bucketDB_allocate

  !> Register member
  !! Before actually register, bucketDB_allocate has to be called
  subroutine bucketDB_register(bktdb, bid, sid)
    implicit none
    type(bucketDB), intent(inout) :: bktdb  !< bucket info
    integer(kind=kint), intent(in) :: bid   !< bucket ID
    integer(kind=kint), intent(in) :: sid   !< member ID
    integer(kind=kint) :: baddr(3)
    baddr = decode_bid(bktdb, bid)
    call assert(all(baddr > 0) .and. all(baddr <= bktdb%ndiv), 'bucketDB_register: block ID our of range')
    call bucket_register(bktdb%buckets(baddr(1),baddr(2),baddr(3)), sid)
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_register: ', baddr, sid
  end subroutine bucketDB_register

  !> Get number of candidates within neighboring buckets of a given bucket
  !! Bucket ID has to be obtained with bucketDB_getBucketID
  function bucketDB_getNumCand(bktdb, bid)
    implicit none
    integer(kind=kint) :: bucketDB_getNumCand  !< number of candidates
    type(bucketDB), intent(in) :: bktdb        !< bucket info
    integer(kind=kint), intent(in) :: bid      !< bucket ID
    integer(kind=kint) :: baddr(3), ncand, i, j, k, is, ie, js, je, ks, ke
    if (bid < 0) then
      bucketDB_getNumCand = 0
      return
    endif
    baddr = decode_bid(bktdb, bid)
    ncand = 0
    is = max(baddr(1)-1, 1)
    ie = min(baddr(1)+1, bktdb%ndiv(1))
    js = max(baddr(2)-1, 1)
    je = min(baddr(2)+1, bktdb%ndiv(2))
    ks = max(baddr(3)-1, 1)
    ke = min(baddr(3)+1, bktdb%ndiv(3))
    do k = ks, ke
      do j = js, je
        do i = is, ie
          ncand = ncand + bucket_get_n(bktdb%buckets(i,j,k))
        enddo
      enddo
    enddo
    bucketDB_getNumCand = ncand
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_getNumCand: ',ncand
  end function bucketDB_getNumCand

  !> Get candidates within neighboring buckets of a given bucket
  !! Number of candidates has to be obtained with bucketDB_getNumCand beforehand
  subroutine bucketDB_getCand(bktdb, bid, ncand, cand)
    implicit none
    type(bucketDB), intent(in) :: bktdb                     !< bucket info
    integer(kind=kint), intent(in) :: bid                   !< bucket ID
    integer(kind=kint), intent(in) :: ncand                 !< number of candidates
    integer(kind=kint), intent(out), target :: cand(ncand)  !< array to store candidates
    integer(kind=kint) :: baddr(3), i, j, k, n, cnt, is, ie, js, je, ks, ke
    integer(kind=kint), pointer :: pcand(:)
    if (bid < 0 .or. ncand == 0) return
    baddr = decode_bid(bktdb, bid)
    cnt = 0
    is = max(baddr(1)-1, 1)
    ie = min(baddr(1)+1, bktdb%ndiv(1))
    js = max(baddr(2)-1, 1)
    je = min(baddr(2)+1, bktdb%ndiv(2))
    ks = max(baddr(3)-1, 1)
    ke = min(baddr(3)+1, bktdb%ndiv(3))
    do k = ks, ke
      do j = js, je
        do i = is, ie
          n = bucket_get_n(bktdb%buckets(i,j,k))
          pcand => cand(cnt+1:cnt+n)
          call bucket_get_member(bktdb%buckets(i,j,k), n, pcand)
          cnt = cnt + n
          call assert(cnt <= ncand, 'bucketDB_get_cand: array overflow')
        enddo
      enddo
    enddo
    call assert(cnt == ncand, 'bucketDB_get_cand: count mismatch')
    if (DEBUG >= 3) write(0,*) '    DEBUG: bucketDB_getCand: ',cand
  end subroutine bucketDB_getCand

end module bucket_search
