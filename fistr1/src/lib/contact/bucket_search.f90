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

  type bucket
    integer(kind=kint) :: n, n_max
    integer(kind=kint), pointer :: member(:) => null()
  end type bucket

  type bucketDB
    private
    real(kind=kreal) :: x_min(3)
    real(kind=kreal) :: x_max(3)
    real(kind=kreal) :: d(3)
    integer(kind=kint) :: ndiv(3)
    type(bucket), pointer :: buckets(:,:,:) => null()
    integer(kind=kint) :: n_max
    integer(kind=kint), pointer :: member_all(:) => null()
  end type bucketDB

contains

  subroutine assert(cond, mesg)
    implicit none
    logical, intent(in) :: cond
    character(len=*) :: mesg
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

  subroutine bucket_init(bkt)
    implicit none
    type(bucket), intent(inout) :: bkt
    bkt%n = 0
    bkt%n_max = 0
    nullify(bkt%member)
  end subroutine bucket_init

  subroutine bucket_finalize(bkt)
    implicit none
    type(bucket), intent(inout) :: bkt
    !if (bkt%n > 0) deallocate(bkt%member)
    nullify(bkt%member)
    bkt%n = 0
    bkt%n_max = 0
  end subroutine bucket_finalize

  subroutine bucket_incr_count(bkt)
    implicit none
    type(bucket), intent(inout) :: bkt
    !$omp atomic
    bkt%n = bkt%n + 1
  end subroutine bucket_incr_count

  ! subroutine bucket_allocate(bkt)
  !   implicit none
  !   type(bucket), intent(inout) :: bkt
  !   if (bkt%n /= bkt%n_max) then
  !     if (associated(bkt%member)) deallocate(bkt%member)
  !   endif
  !   if (bkt%n > 0) allocate(bkt%member(bkt%n))
  !   bkt%n_max = bkt%n
  !   bkt%n = 0
  ! end subroutine bucket_allocate

  subroutine bucket_assign(bkt, memb)
    implicit none
    type(bucket), intent(inout) :: bkt
    integer(kind=kint), pointer :: memb(:)
    bkt%member => memb
    bkt%n_max = bkt%n
    bkt%n = 0
  end subroutine bucket_assign

  subroutine bucket_register(bkt, sid)
    implicit none
    type(bucket), intent(inout) :: bkt
    integer(kind=kint), intent(in) :: sid
    !$omp atomic
    bkt%n = bkt%n + 1
    call assert(bkt%n <= bkt%n_max, 'bucket_register: too many members')
    bkt%member(bkt%n) = sid
  end subroutine bucket_register

  function bucket_get_n(bkt)
    implicit none
    integer(kind=kint) :: bucket_get_n
    type(bucket), intent(in) :: bkt
    bucket_get_n = bkt%n
  end function bucket_get_n

  subroutine bucket_get_member(bkt, n, memb)
    implicit none
    type(bucket), intent(in) :: bkt
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(out) :: memb(n)
    call assert(n == bkt%n, 'bucket_get_member: wrong n')
    memb(1:n) = bkt%member(1:n)
  end subroutine bucket_get_member

!!!
!!! routines for type(bucketDB)
!!!

  subroutine bucketDB_init(bktdb)
    implicit none
    type(bucketDB), intent(inout) :: bktdb
    bktdb%x_min(:) = 0.d0
    bktdb%x_max(:) = 0.d0
    bktdb%d(:) = 0.d0
    bktdb%ndiv(:) = 0
    nullify(bktdb%buckets)
    bktdb%n_max = 0
    nullify(bktdb%member_all)
  end subroutine bucketDB_init

  subroutine bucketDB_finalize(bktdb)
    implicit none
    type(bucketDB), intent(inout) :: bktdb
    integer(kind=kint) :: i, j, k
    if (bktdb%n_max > 0) then
      deallocate(bktdb%member_all)
      bktdb%n_max = 0
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

  subroutine bucketDB_setup(bktdb, x_min, x_max, dmin, n_max)
    implicit none
    type(bucketDB), intent(inout) :: bktdb
    real(kind=kreal), intent(in) :: x_min(3)
    real(kind=kreal), intent(in) :: x_max(3)
    real(kind=kreal), intent(in) :: dmin
    integer(kind=kint), intent(in) :: n_max
    real(kind=kreal) :: xrange(3)
    integer(kind=kint) :: i, j, k
    real(kind=kreal), parameter :: EPS = 1.d-6
    if (DEBUG >= 1) write(0,*) 'DEBUG: bucketDB_setup', x_min, x_max, dmin, n_max
    if (associated(bktdb%buckets)) deallocate(bktdb%buckets)
    bktdb%x_min(:) = x_min(:)
    bktdb%x_max(:) = x_max(:)
    xrange(:) = x_max(:) - x_min(:)
    call assert(all(xrange > 0.d0), 'bucketDB_setup: invalid x_min, x_max')
    do i = 1, 3
      bktdb%ndiv(i) = max(floor(xrange(i) / dmin), 1)
      bktdb%d(i) = xrange(i) / bktdb%ndiv(i) * (1.d0 + EPS)
    enddo
    allocate(bktdb%buckets(bktdb%ndiv(1), bktdb%ndiv(2), bktdb%ndiv(3)))
    do k = 1, bktdb%ndiv(3)
      do j = 1, bktdb%ndiv(2)
        do i = 1, bktdb%ndiv(1)
          call bucket_init(bktdb%buckets(i,j,k))
        enddo
      enddo
    enddo
    if (bktdb%n_max /= n_max) then
      if (associated(bktdb%member_all)) deallocate(bktdb%member_all)
      allocate(bktdb%member_all(n_max))
      bktdb%n_max = n_max
    endif
  end subroutine bucketDB_setup

  function encode_bid(bktdb, bid)
    implicit none
    integer(kind=kint) :: encode_bid
    type(bucketDB), intent(in) :: bktdb
    integer(kind=kint), intent(in) :: bid(3)
    if (any(bid <= 0) .or. any(bid > bktdb%ndiv)) then
      encode_bid = -1
    else
      encode_bid = &
           (bid(3)-1) * bktdb%ndiv(1) * bktdb%ndiv(2) + (bid(2)-1) * bktdb%ndiv(1) + bid(1)
    endif
  end function encode_bid

  function decode_bid(bktdb, bidenc)
    implicit none
    integer(kind=kint) :: decode_bid(3)
    type(bucketDB), intent(in) :: bktdb
    integer(kind=kint), intent(in) :: bidenc
    call assert(bidenc <= bktdb%ndiv(1)*bktdb%ndiv(2)*bktdb%ndiv(3), 'decode_bid: out of range')
    if (bidenc < 0) then
      decode_bid(:) = -1
    else
      decode_bid(1) = mod(bidenc-1, bktdb%ndiv(1)) + 1
      decode_bid(2) = mod((bidenc-1)/bktdb%ndiv(1), bktdb%ndiv(2)) + 1
      decode_bid(3) = (bidenc-1)/(bktdb%ndiv(1) * bktdb%ndiv(2)) + 1
      call assert(encode_bid(bktdb, decode_bid) == bidenc, 'decode_bid')
    endif
  end function decode_bid

  function bucketDB_getBucketID(bktdb, x)
    implicit none
    integer(kind=kint) :: bucketDB_getBucketID
    type(bucketDB), intent(in) :: bktdb
    real(kind=kreal), intent(in) :: x(3)
    integer(kind=kint) :: bid(3)
    integer(kind=kint) :: i
    do i = 1, 3
      bid(i) = floor((x(i) - bktdb%x_min(i)) / bktdb%d(i)) + 1
    enddo
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_getBucketID: ',x,bid
    bucketDB_getBucketID = encode_bid(bktdb, bid)
  end function bucketDB_getBucketID

  subroutine bucketDB_registerPre(bktdb, bidenc)
    implicit none
    type(bucketDB), intent(inout) :: bktdb
    integer(kind=kint), intent(in) :: bidenc
    integer(kind=kint) :: bid(3)
    bid = decode_bid(bktdb, bidenc)
    call assert(all(bid > 0) .and. all(bid <= bktdb%ndiv), 'bucketDB_register_pre: block ID our of range')
    call bucket_incr_count(bktdb%buckets(bid(1),bid(2),bid(3)))
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_registerPre: ', bid
  end subroutine bucketDB_registerPre

  subroutine bucketDB_allocate(bktdb)
    implicit none
    type(bucketDB), intent(inout) :: bktdb
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

  subroutine bucketDB_register(bktdb, bidenc, sid)
    implicit none
    type(bucketDB), intent(inout) :: bktdb
    integer(kind=kint), intent(in) :: bidenc
    integer(kind=kint), intent(in) :: sid
    integer(kind=kint) :: bid(3)
    bid = decode_bid(bktdb, bidenc)
    call assert(all(bid > 0) .and. all(bid <= bktdb%ndiv), 'bucketDB_register: block ID our of range')
    call bucket_register(bktdb%buckets(bid(1),bid(2),bid(3)), sid)
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_register: ', bid, sid
  end subroutine bucketDB_register

  function bucketDB_getNumCand(bktdb, bidenc)
    implicit none
    integer(kind=kint) :: bucketDB_getNumCand
    type(bucketDB), intent(in) :: bktdb
    integer(kind=kint), intent(in) :: bidenc
    integer(kind=kint) :: bid(3), ncand, i, j, k
    bid = decode_bid(bktdb, bidenc)
    ncand = 0
    do k = bid(3)-1, bid(3)+1
      if (k <= 0 .or. k > bktdb%ndiv(3)) cycle
      do j = bid(2)-1, bid(2)+1
        if (j <= 0 .or. j > bktdb%ndiv(2)) cycle
        do i = bid(1)-1, bid(1)+1
          if (i <= 0 .or. i > bktdb%ndiv(1)) cycle
          ncand = ncand + bucket_get_n(bktdb%buckets(i,j,k))
        enddo
      enddo
    enddo
    bucketDB_getNumCand = ncand
    if (DEBUG >= 2) write(0,*) '  DEBUG: bucketDB_getNumCand: ',ncand
  end function bucketDB_getNumCand

  subroutine bucketDB_getCand(bktdb, bidenc, ncand, cand)
    implicit none
    type(bucketDB), intent(in) :: bktdb
    integer(kind=kint), intent(in) :: bidenc
    integer(kind=kint), intent(in) :: ncand
    integer(kind=kint), intent(out), target :: cand(ncand)
    integer(kind=kint) :: bid(3), i, j, k, n, cnt
    integer(kind=kint), pointer :: pcand(:)
    bid = decode_bid(bktdb, bidenc)
    cnt = 0
    do k = bid(3)-1, bid(3)+1
      if (k <= 0 .or. k > bktdb%ndiv(3)) cycle
      do j = bid(2)-1, bid(2)+1
        if (j <= 0 .or. j > bktdb%ndiv(2)) cycle
        do i = bid(1)-1, bid(1)+1
          if (i <= 0 .or. i > bktdb%ndiv(1)) cycle
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
