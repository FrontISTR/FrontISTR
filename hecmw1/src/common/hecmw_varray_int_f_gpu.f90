!-------------------------------------------------------------------------------
! Copyright (c) 2022 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief varray int

module hecmw_varray_int_gpu
  use hecmw_util, only: kint
  implicit none

  private
  public :: hecmwST_varray_int_gpu_container
  public :: hecmwST_varray_int_gpu
  public :: HECMW_varray_int_gpu_container_initialize
  public :: HECMW_varray_int_gpu_container_finalize
  public :: HECMW_varray_int_gpu_initialize
  public :: HECMW_varray_int_gpu_initialize_all
  public :: HECMW_varray_int_gpu_finalize
  public :: HECMW_varray_int_gpu_finalize_all
  public :: HECMW_varray_int_gpu_clear
  public :: HECMW_varray_int_gpu_add
  public :: HECMW_varray_int_gpu_add_if_not_exits
  public :: HECMW_varray_int_gpu_insert
  public :: HECMW_varray_int_gpu_insert_if_not_exists
  public :: HECMW_varray_int_gpu_expand
  public :: HECMW_varray_int_gpu_print
  public :: HECMW_varray_int_gpu_print_all
  public :: HECMW_varray_int_gpu_find
  public :: HECMW_varray_int_gpu_get_nitem
  public :: HECMW_varray_int_gpu_get_item
  public :: HECMW_varray_int_gpu_get_item_all

  type hecmwST_varray_int_gpu_buffer
    private
    integer :: cap
    integer :: size
    integer(kind=kint), pointer :: buffer(:) => null()
  end type hecmwST_varray_int_gpu_buffer

  type hecmwST_varray_int_gpu_container
    private
    type(hecmwST_varray_int_gpu_buffer) :: buffer(50)
    integer, pointer :: copy_buffer(:) => null()
  end type hecmwST_varray_int_gpu_container

  type hecmwST_varray_int_gpu
    private
    integer(kind=kint) :: nitem
    integer(kind=kint), pointer :: items(:)
  end type hecmwST_varray_int_gpu

contains

  subroutine HECMW_varray_int_gpu_container_initialize(container)
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container

    integer :: i

    do i = 1, size(container%buffer)
      container%buffer(i)%cap = 0
      container%buffer(i)%size = 0
    enddo

    container%buffer(1)%cap = 64
    allocate(container%buffer(1)%buffer(container%buffer(1)%cap))

    allocate(container%copy_buffer(1024))
  end subroutine HECMW_varray_int_gpu_container_initialize

  subroutine HECMW_varray_int_gpu_container_finalize(container)
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container

    integer :: i

    do i = 1, size(container%buffer)
      if(container%buffer(i)%cap > 0) then
        deallocate(container%buffer(i)%buffer)
      endif
    enddo

    deallocate(container%copy_buffer)
  end subroutine HECMW_varray_int_gpu_container_finalize

  subroutine HECMW_varray_int_gpu_container_get_list(container, alloc_size, list)
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    integer, intent(in) :: alloc_size
    integer(kind=kint), pointer, intent(inout) :: list(:)

    integer :: pre
    integer :: i
    integer :: j
    integer :: p

    do i = 1, size(container%buffer)
      if(container%buffer(i)%cap == 0) then
        container%buffer(i)%cap = pre * 2
        allocate(container%buffer(i)%buffer(container%buffer(i)%cap))
      endif

      if(container%buffer(i)%cap - container%buffer(i)%size > alloc_size) then
        list => container%buffer(i)%buffer(container%buffer(i)%size+1 : container%buffer(i)%size+alloc_size)
        container%buffer(i)%size = container%buffer(i)%size + alloc_size
        exit
      endif

      pre = container%buffer(i)%cap
    end do
  end subroutine HECMW_varray_int_gpu_container_get_list

  subroutine HECMW_varray_int_gpu_container_release_list(container, list)
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    integer(kind=kint), pointer, intent(inout) :: list(:)

    !-- do nothing
  end subroutine HECMW_varray_int_gpu_container_release_list

  subroutine HECMW_varray_int_gpu_initialize( container, ilist, n_init_in )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist
    integer(kind=kint), intent(in) :: n_init_in

    integer(kind=kint) :: n_init

    n_init = max(n_init_in,1)

    ilist%nitem = 0
    call HECMW_varray_int_gpu_container_get_list(container, n_init, ilist%items)
    ilist%items(:) = 0
  end subroutine HECMW_varray_int_gpu_initialize

  subroutine HECMW_varray_int_gpu_initialize_all( container, ilists, nlists, n_init_in )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), allocatable, intent(inout) :: ilists(:)
    integer(kind=kint), intent(in) :: nlists
    integer(kind=kint), intent(in) :: n_init_in

    integer(kind=kint) :: i

    allocate(ilists(nlists))
    do i=1,size(ilists)
      call HECMW_varray_int_gpu_initialize( container, ilists(i), n_init_in )
    end do
  end subroutine HECMW_varray_int_gpu_initialize_all

  subroutine HECMW_varray_int_gpu_finalize( container, ilist )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist

    ilist%nitem = 0
    call HECMW_varray_int_gpu_container_release_list(container, ilist%items)
  end subroutine HECMW_varray_int_gpu_finalize

  subroutine HECMW_varray_int_gpu_finalize_all( container, ilists )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), allocatable, intent(inout) :: ilists(:)

    integer(kind=kint) :: i

    do i=1,size(ilists)
      ilists(i)%nitem = 0
      call HECMW_varray_int_gpu_container_release_list(container, ilists(i)%items)
    end do
    deallocate(ilists)
  end subroutine HECMW_varray_int_gpu_finalize_all

  subroutine HECMW_varray_int_gpu_clear( ilist )
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist

    ilist%nitem = 0
    ilist%items(:) = 0
  end subroutine HECMW_varray_int_gpu_clear

  subroutine HECMW_varray_int_gpu_enlarge( container, ilist )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist

    integer(kind=kint) :: i, size_ilist
    integer(kind=kint), pointer :: tmp(:)

    size_ilist = size(ilist%items)

    if(size_ilist <= size(container%copy_buffer)) then
      tmp => container%copy_buffer(:)
    else
      allocate(tmp(size_ilist))
    endif

    do i=1,ilist%nitem
      tmp(i) = ilist%items(i)
    end do

    call HECMW_varray_int_gpu_container_release_list(container, ilist%items)
    call HECMW_varray_int_gpu_container_get_list(container, size_ilist*2, ilist%items)

    ilist%items(:) = 0
    do i=1,ilist%nitem
      ilist%items(i) = tmp(i)
    end do

    if(size_ilist <= size(container%copy_buffer)) then
      nullify(tmp)
    else
      deallocate(tmp)
    endif
  end subroutine HECMW_varray_int_gpu_enlarge

  subroutine HECMW_varray_int_gpu_add( container, ilist, ival )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist
    integer(kind=kint), intent(in) :: ival

    if( ilist%nitem == size(ilist%items) ) &
      &  call HECMW_varray_int_gpu_enlarge( container, ilist )

    ilist%nitem = ilist%nitem + 1
    ilist%items(ilist%nitem) = ival
  end subroutine HECMW_varray_int_gpu_add

  subroutine HECMW_varray_int_gpu_add_if_not_exits( container, ilist, ival )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist
    integer(kind=kint), intent(in) :: ival

    integer(kind=kint) :: i

    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) return
    end do

    call HECMW_varray_int_gpu_add( container, ilist, ival  )
  end subroutine HECMW_varray_int_gpu_add_if_not_exits

  subroutine HECMW_varray_int_gpu_insert( container, ilist, ival )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist
    integer(kind=kint), intent(in) :: ival

    integer(kind=kint) :: i, cval, tmpval

    cval = ival
    do i=1,ilist%nitem
      if( ilist%items(i) < cval ) cycle
      tmpval = ilist%items(i)
      ilist%items(i) = cval
      cval = tmpval
    end do

    call HECMW_varray_int_gpu_add( container, ilist, cval )
  end subroutine HECMW_varray_int_gpu_insert

  subroutine HECMW_varray_int_gpu_insert_if_not_exists( container, ilist, ival )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist
    integer(kind=kint), intent(in) :: ival

    integer(kind=kint) :: i

    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) return
    end do

    call HECMW_varray_int_gpu_insert( container, ilist, ival )
  end subroutine HECMW_varray_int_gpu_insert_if_not_exists

  subroutine HECMW_varray_int_gpu_expand( container, ilist, n, vals )
    type(hecmwST_varray_int_gpu_container), intent(inout) :: container
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist
    integer(kind=kint), intent(in) :: n
    integer(kind=kint), intent(in) :: vals(:)

    do while( ilist%nitem + n > size(ilist%items) )
      call HECMW_varray_int_gpu_enlarge( container, ilist )
    end do

    ilist%items(ilist%nitem+1:ilist%nitem+n) = vals(1:n)
    ilist%nitem = ilist%nitem + n
  end subroutine HECMW_varray_int_gpu_expand

  subroutine HECMW_varray_int_gpu_print( ilist )
    type( hecmwST_varray_int_gpu ), intent(inout) :: ilist

    write(*,*) 'n, maxn:', ilist%nitem, size(ilist%items)
    write(*,*) 'items:', ilist%items(1:ilist%nitem)
  end subroutine HECMW_varray_int_gpu_print

  subroutine HECMW_varray_int_gpu_print_all( ilists )
    type( hecmwST_varray_int_gpu ), allocatable, intent(inout) :: ilists(:)

    integer(kind=kint) :: i

    do i=1,size(ilists)
      write(*,*) "i, n, maxn: ", i, ilists(i)%nitem, size(ilists(i)%items)
      write(*,*) 'items:', ilists(i)%items(1:ilists(i)%nitem)
    end do
  end subroutine HECMW_varray_int_gpu_print_all

  integer(kind=kint) function HECMW_varray_int_gpu_find( ilist, ival )
    type( hecmwST_varray_int_gpu ), intent(in) :: ilist
    integer(kind=kint), intent(in)          :: ival

    integer(kind=kint) :: i

    HECMW_varray_int_gpu_find = -1
    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) HECMW_varray_int_gpu_find = i
    end do
  end function HECMW_varray_int_gpu_find

  integer(kind=kint) function HECMW_varray_int_gpu_get_nitem( ilist )
    type( hecmwST_varray_int_gpu ), intent(in) :: ilist

    HECMW_varray_int_gpu_get_nitem = ilist%nitem
  end function HECMW_varray_int_gpu_get_nitem

  integer(kind=kint) function HECMW_varray_int_gpu_get_item( ilist, n )
    type( hecmwST_varray_int_gpu ), intent(in) :: ilist
    integer(kind=kint), intent(in) :: n

    HECMW_varray_int_gpu_get_item = ilist%items(n)
  end function HECMW_varray_int_gpu_get_item

  subroutine HECMW_varray_int_gpu_get_item_all( ilist, array )
    type( hecmwST_varray_int_gpu ), intent(in) :: ilist
    integer(kind=kint), intent(inout) :: array(:)

    array(1:ilist%nitem) = ilist%items(1:ilist%nitem)
  end subroutine HECMW_varray_int_gpu_get_item_all

end module hecmw_varray_int_gpu
