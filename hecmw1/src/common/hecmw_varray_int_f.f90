!-------------------------------------------------------------------------------
! Copyright (c) 2022 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief varray int

module hecmw_varray_int
  implicit none

  private
  public :: hecmwST_varray_int
  public :: HECMW_varray_int_initialize
  public :: HECMW_varray_int_initialize_all
  public :: HECMW_varray_int_finalize
  public :: HECMW_varray_int_finalize_all
  public :: HECMW_varray_int_clear
  public :: HECMW_varray_int_add
  public :: HECMW_varray_int_add_if_not_exits
  public :: HECMW_varray_int_insert
  public :: HECMW_varray_int_insert_if_not_exists
  public :: HECMW_varray_int_print
  public :: HECMW_varray_int_print_all
  public :: HECMW_varray_int_find
  public :: HECMW_varray_int_get_nitem
  public :: HECMW_varray_int_get_item

  type hecmwST_varray_int
    private
    integer :: nitem
    integer, allocatable :: items(:)
  end type hecmwST_varray_int

contains

  subroutine HECMW_varray_int_initialize( ilist, n_init_in )
    type( hecmwST_varray_int ), intent(inout) :: ilist
    integer, intent(in) :: n_init_in

    integer :: n_init

    n_init = max(n_init_in,1)

    ilist%nitem = 0
    allocate( ilist%items(n_init) )
    ilist%items(:) = 0
  end subroutine HECMW_varray_int_initialize

  subroutine HECMW_varray_int_initialize_all( ilists, nlists, n_init_in )
    type( hecmwST_varray_int ), allocatable, intent(inout) :: ilists(:)
    integer, intent(in) :: nlists
    integer, intent(in) :: n_init_in

    integer :: i

    allocate(ilists(nlists))
    do i=1,size(ilists)
      call HECMW_varray_int_initialize( ilists(i), n_init_in )
    end do
  end subroutine HECMW_varray_int_initialize_all

  subroutine HECMW_varray_int_finalize( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    ilist%nitem = 0
    deallocate( ilist%items )
  end subroutine HECMW_varray_int_finalize

  subroutine HECMW_varray_int_finalize_all( ilists )
    type( hecmwST_varray_int ), allocatable, intent(inout) :: ilists(:)

    integer :: i

    do i=1,size(ilists)
      ilists(i)%nitem = 0
      deallocate( ilists(i)%items )
    end do
    deallocate(ilists)
  end subroutine HECMW_varray_int_finalize_all

  subroutine HECMW_varray_int_clear( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    ilist%nitem = 0
    ilist%items(:) = 0
  end subroutine HECMW_varray_int_clear

  subroutine HECMW_varray_int_enlarge( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    integer :: i, size_ilist
    integer, allocatable :: tmp(:)

    size_ilist = size(ilist%items)
    allocate(tmp(size_ilist))
    do i=1,ilist%nitem
      tmp(i) = ilist%items(i)
    end do

    deallocate(ilist%items)
    allocate(ilist%items(size_ilist*2))

    ilist%items(:) = 0
    do i=1,ilist%nitem
      ilist%items(i) = tmp(i)
    end do
  end subroutine HECMW_varray_int_enlarge

  subroutine HECMW_varray_int_add( ilist, ival )
    type( hecmwST_varray_int ), intent(inout) :: ilist
    integer, intent(in) :: ival

    if( ilist%nitem == size(ilist%items) ) &
      &  call HECMW_varray_int_enlarge( ilist )

    ilist%nitem = ilist%nitem + 1
    ilist%items(ilist%nitem) = ival
  end subroutine HECMW_varray_int_add

  subroutine HECMW_varray_int_add_if_not_exits( ilist, ival )
    type( hecmwST_varray_int ), intent(inout) :: ilist
    integer, intent(in) :: ival

    integer :: i

    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) return
    end do

    call HECMW_varray_int_add( ilist, ival  )
  end subroutine HECMW_varray_int_add_if_not_exits

  subroutine HECMW_varray_int_insert( ilist, ival )
    type( hecmwST_varray_int ), intent(inout) :: ilist
    integer, intent(in) :: ival

    integer :: i, cval, tmpval

    cval = ival
    do i=1,ilist%nitem
      if( ilist%items(i) < cval ) cycle
      tmpval = ilist%items(i)
      ilist%items(i) = cval
      cval = tmpval
    end do

    call HECMW_varray_int_add( ilist, cval )
  end subroutine HECMW_varray_int_insert

  subroutine HECMW_varray_int_insert_if_not_exists( ilist, ival )
    type( hecmwST_varray_int ), intent(inout) :: ilist
    integer, intent(in) :: ival

    integer :: i

    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) return
    end do

    call HECMW_varray_int_insert( ilist, ival )
  end subroutine HECMW_varray_int_insert_if_not_exists

  subroutine HECMW_varray_int_expand( ilist, n, vals )
    type( hecmwST_varray_int ), intent(inout) :: ilist
    integer, intent(in) :: n
    integer, intent(in) :: vals(:)

    do while( ilist%nitem + n > size(ilist%items) )
      call HECMW_varray_int_enlarge( ilist )
    end do

    ilist%items(ilist%nitem+1:ilist%nitem+n) = vals(1:n)
    ilist%nitem = ilist%nitem + n
  end subroutine HECMW_varray_int_expand

  subroutine HECMW_varray_int_print( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    write(*,*) 'n, maxn:', ilist%nitem, size(ilist%items)
    write(*,*) 'items:', ilist%items(1:ilist%nitem)
  end subroutine HECMW_varray_int_print

  subroutine HECMW_varray_int_print_all( ilists )
    type( hecmwST_varray_int ), allocatable, intent(inout) :: ilists(:)

    integer :: i

    do i=1,size(ilists)
      write(*,*) "i, n, maxn: ", i, ilists(i)%nitem, size(ilists(i)%items)
      write(*,*) 'items:', ilists(i)%items(1:ilists(i)%nitem)
    end do
  end subroutine HECMW_varray_int_print_all

  integer function HECMW_varray_int_find( ilist, ival )
    type( hecmwST_varray_int ), intent(in) :: ilist
    integer, intent(in)          :: ival

    integer :: i

    HECMW_varray_int_find = -1
    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) HECMW_varray_int_find = i
    end do
  end function HECMW_varray_int_find

  integer function HECMW_varray_int_get_nitem( ilist )
    type( hecmwST_varray_int), intent(in) :: ilist

    HECMW_varray_int_get_nitem = ilist%nitem
  end function HECMW_varray_int_get_nitem

  integer function HECMW_varray_int_get_item( ilist, n )
    type( hecmwST_varray_int), intent(in) :: ilist
    integer, intent(in) :: n

    HECMW_varray_int_get_item = ilist%items(n)
  end function HECMW_varray_int_get_item

end module hecmw_varray_int
