!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief varray int

module hecmw_varray_int
  implicit none

  type hecmwST_varray_int
    integer :: nitem
    integer, allocatable :: items(:)
  end type hecmwST_varray_int

  public

contains

  subroutine HECMW_varray_int_initialize( n_init_in, ilist )
    integer, intent(in) :: n_init_in
    type( hecmwST_varray_int ), intent(inout) :: ilist

    integer :: n_init

    n_init = max(n_init_in,1)

    ilist%nitem = 0
    allocate( ilist%items(n_init) )
    ilist%items(:) = 0
  end subroutine

  subroutine HECMW_varray_int_initialize_all( nlists, n_init_in, ilists )
    integer, intent(in) :: nlists
    integer, intent(in) :: n_init_in
    type( hecmwST_varray_int ), allocatable, intent(inout) :: ilists(:)

    integer :: i

    allocate(ilists(nlists))
    do i=1,size(ilists)
      call HECMW_varray_int_initialize( n_init_in, ilists(i) )
    end do
  end subroutine

  subroutine HECMW_varray_int_finalize( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    ilist%nitem = 0
    deallocate( ilist%items )
  end subroutine

  subroutine HECMW_varray_int_finalize_all( ilists )
    type( hecmwST_varray_int ), allocatable, intent(inout) :: ilists(:)

    integer :: i

    do i=1,size(ilists)
      ilists(i)%nitem = 0
      deallocate( ilists(i)%items )
    end do
    deallocate(ilists)
  end subroutine

  subroutine HECMW_varray_int_clear( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    ilist%nitem = 0
    ilist%items(:) = 0
  end subroutine

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
  end subroutine

  subroutine HECMW_varray_int_add( ival, ilist )
    integer, intent(in) :: ival
    type( hecmwST_varray_int ), intent(inout) :: ilist

    if( ilist%nitem == size(ilist%items) ) &
      &  call HECMW_varray_int_enlarge( ilist )

    ilist%nitem = ilist%nitem + 1
    ilist%items(ilist%nitem) = ival
  end subroutine

  subroutine HECMW_varray_int_add_if_not_exits( ival, ilist )
    integer, intent(in) :: ival
    type( hecmwST_varray_int ), intent(inout) :: ilist

    integer :: i

    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) return
    end do

    call HECMW_varray_int_add( ival, ilist )
  end subroutine

  subroutine HECMW_varray_int_insert( ival, ilist )
    integer, intent(in) :: ival
    type( hecmwST_varray_int ), intent(inout) :: ilist

    integer :: i, cval, tmpval

    cval = ival
    do i=1,ilist%nitem
      if( ilist%items(i) < cval ) cycle
      tmpval = ilist%items(i)
      ilist%items(i) = cval
      cval = tmpval
    end do

    call HECMW_varray_int_add( cval, ilist )
  end subroutine

  subroutine HECMW_varray_int_insert_if_not_exists( ival, ilist )
    integer, intent(in) :: ival
    type( hecmwST_varray_int ), intent(inout) :: ilist

    integer :: i

    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) return
    end do

    call HECMW_varray_int_insert( ival, ilist )
  end subroutine

  subroutine HECMW_varray_int_expand( n, vals, ilist )
    integer, intent(in) :: n
    integer, intent(in) :: vals(:)
    type( hecmwST_varray_int ), intent(inout) :: ilist

    do while( ilist%nitem + n > size(ilist%items) )
      call HECMW_varray_int_enlarge( ilist )
    end do

    ilist%items(ilist%nitem+1:ilist%nitem+n) = vals(1:n)
    ilist%nitem = ilist%nitem + n
  end subroutine

  subroutine HECMW_varray_int_print( ilist )
    type( hecmwST_varray_int ), intent(inout) :: ilist

    write(*,*) 'n, maxn:', ilist%nitem, size(ilist%items)
    write(*,*) 'items:', ilist%items(1:ilist%nitem)
  end subroutine

  subroutine HECMW_varray_int_print_all( ilists )
    type( hecmwST_varray_int ), allocatable, intent(inout) :: ilists(:)

    integer :: i

    do i=1,size(ilists)
      write(*,*) "i, n, maxn: ", i, ilists(i)%nitem, size(ilists(i)%items)
      write(*,*) 'items:', ilists(i)%items(1:ilists(i)%nitem)
    end do
  end subroutine

  integer function HECMW_varray_int_find( ival, ilist )
    integer, intent(in)          :: ival
    type( hecmwST_varray_int ), intent(in) :: ilist

    integer :: i

    HECMW_varray_int_find = -1
    do i=1,ilist%nitem
      if( ilist%items(i) == ival ) HECMW_varray_int_find = i
    end do
  end function

end module hecmw_varray_int
