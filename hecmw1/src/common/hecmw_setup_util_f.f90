!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains auxiliary functions in calculation setup

module hecmw_setup_util
  use hecmw

  !> container of character array pointer, because of gfortran's bug
  type hecmw_str_arr
    character(len=HECMW_NAME_LEN), pointer :: s(:)
  end type hecmw_str_arr

  !> private parameter and pointers to group members in hecMESH
  integer(kind=kint),private :: grp_type     ! 1:node_grp, 2:elem_grp, 3:surf_grp
  integer(kind=kint),pointer,private :: n_grp
  integer(kind=kint),pointer,private :: grp_index(:)
  integer(kind=kint),pointer,private :: grp_item(:)
  type(hecmw_str_arr),private :: grp_name

  ! private subroutines ------------
  private :: set_group_pointers
!  private :: append_single_group

contains

  function hecmw_str2index( s, x )
    implicit none
    logical hecmw_str2index
    character(*) :: s
    integer :: i, n, a, i0,i9, m, x, b
    logical :: fg

    hecmw_str2index = .false.
    i0 = iachar('0')
    i9 = iachar('9')
    n = len_trim(s)
    x = 0
    b = 1
    fg = .true.
    do i=n,1,-1
      fg = .false.
      a = iachar(s(i:i))
      if( a < i0 .or. a > i9 ) return
      m = a-i0
      x = x + b * m
      b = b*10
    end do
    hecmw_str2index = .true.
  end function hecmw_str2index

  subroutine hecmw_strupr( s )
    implicit none
    character(*) :: s
    integer :: i, n, a

    n = len_trim(s)
    do i = 1, n
      a = iachar(s(i:i))
      if( a >= iachar('a') .and. a <= iachar('z')) then
        s(i:i) = achar(a - 32)
      end if
    end do
  end subroutine hecmw_strupr

  function hecmw_streqr( s1, s2 )
    implicit none
    character(*) :: s1, s2
    logical :: hecmw_streqr
    integer :: i, n, a1, a2

    hecmw_streqr = .false.
    n = len_trim(s1)
    if( n /= len_trim(s2)) return
    call hecmw_strupr(s1)
    call hecmw_strupr(s2)
    do i = 1, n
      a1 = iachar(s1(i:i))
      a2 = iachar(s2(i:i))
      if( a1 /= a2 ) then
        return
      end if
    end do
    hecmw_streqr = .true.
  end function hecmw_streqr

  subroutine hecmw_setup_util_err_stop( msg )
    implicit none
    character(*) :: msg

    write(*,*) msg
!    write(imsg,*) msg
    call hecmw_abort( hecmw_comm_get_comm())
  end subroutine hecmw_setup_util_err_stop

  ! grp_type_name : 'node_grp', 'elem_grp' or 'surf_grp'

  subroutine set_group_pointers( hecMESH, grp_type_name )
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*) :: grp_type_name

    if( grp_type_name == 'node_grp' ) then
      grp_type = 1
      n_grp      => hecMESH%node_group%n_grp
      grp_name%s => hecMESH%node_group%grp_name
      grp_index  => hecMESH%node_group%grp_index
      grp_item   => hecMESH%node_group%grp_item
    else if( grp_type_name == 'elem_grp' ) then
      grp_type = 2
      n_grp      => hecMESH%elem_group%n_grp
      grp_name%s => hecMESH%elem_group%grp_name
      grp_index  => hecMESH%elem_group%grp_index
      grp_item   => hecMESH%elem_group%grp_item
    else if( grp_type_name == 'surf_grp' ) then
      grp_type = 3
      n_grp      => hecMESH%surf_group%n_grp
      grp_name%s => hecMESH%surf_group%grp_name
      grp_index  => hecMESH%surf_group%grp_index
      grp_item   => hecMESH%surf_group%grp_item
    else
      stop 'assert in set_group_pointers'
    end if
  end subroutine set_group_pointers

  subroutine backset_group_pointers( hecMESH, grp_type_name )
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*) :: grp_type_name

    if( grp_type_name == 'node_grp' ) then
      grp_type = 1
      hecMESH%node_group%grp_name => grp_name%s
      hecMESH%node_group%grp_index => grp_index
      hecMESH%node_group%grp_item => grp_item
    else if( grp_type_name == 'elem_grp' ) then
      grp_type = 2
      hecMESH%elem_group%grp_name => grp_name%s
      hecMESH%elem_group%grp_index => grp_index
      hecMESH%elem_group%grp_item => grp_item
    else if( grp_type_name == 'surf_grp' ) then
      grp_type = 3
      hecMESH%surf_group%grp_name => grp_name%s
      hecMESH%surf_group%grp_index => grp_index
      hecMESH%surf_group%grp_item => grp_item
    else
      stop 'assert in backset_group_pointers'
    end if
  end subroutine backset_group_pointers

  function node_global_to_local( hecMESH, list, n )
    implicit none
    type (hecmwST_local_mesh), target :: hecMESH
    integer(kind=kint) :: list(:)
    integer(kind=kint) :: n, i, j, cache
    logical:: fg
    integer(kind=kint):: node_global_to_local

    node_global_to_local = 0
    cache = 1
    aa:do j=1, n
      fg = .false.

      do i=cache, hecMESH%n_node
        if( hecMESH%global_node_ID(i) == list(j)) then
          list(j) = i
          cache = i+1
          fg = .true.
          node_global_to_local = node_global_to_local +1
          cycle aa
        endif
      enddo

      do i=1, cache-1
        if( hecMESH%global_node_ID(i) == list(j)) then
          list(j) = i
          cache = i+1
          fg = .true.
          node_global_to_local = node_global_to_local +1
          cycle aa
        endif
      enddo

      if( .not. fg ) then
        list(j) = -1 ! not exist node
      endif
    enddo aa
  end function node_global_to_local

  function elem_global_to_local( hecMESH, list, n )
    implicit none
    type (hecmwST_local_mesh), target :: hecMESH
    integer(kind=kint), pointer :: list(:)
    integer(kind=kint) :: n, i, j
    logical :: fg
    integer(kind=kint) :: elem_global_to_local

    elem_global_to_local = 0
    do j=1, n
      fg = .false.
      do i=1, hecMESH%n_elem
        if( hecMESH%global_elem_ID(i) == list(j)) then
          list(j) = i
          fg = .true.
          elem_global_to_local = elem_global_to_local+1
          exit
        endif
      end do
      if( .not. fg ) then
        list(j) = -1
      endif
    end do
  end function elem_global_to_local

  function append_single_group( hecMESH, grp_type_name, no_count, no_list )
    implicit none
    type (hecmwST_local_mesh), target :: hecMESH
    character(len=*) :: grp_type_name
    integer(kind=kint) :: no_count
    integer(kind=kint),pointer :: no_list(:)
    integer(kind=kint):: append_single_group
    integer(kind=kint) :: old_grp_number, new_grp_number
    integer(kind=kint) :: old_item_number, new_item_number
    integer(kind=kint) :: i,j,k, exist_n
    integer(kind=kint), save :: grp_count = 1
    character(50) :: grp_name_s

    exist_n = 0
    call set_group_pointers( hecMESH, grp_type_name )
    if( grp_type_name == 'node_grp') then
      exist_n = node_global_to_local( hecMESH, no_list, no_count )
    else if( grp_type_name == 'elem_grp') then
      exist_n = elem_global_to_local( hecMESH, no_list, no_count )
    endif

    old_grp_number = n_grp
    new_grp_number = old_grp_number + no_count

    old_item_number = grp_index(n_grp)
    new_item_number = old_item_number + exist_n

    call hecmw_expand_name_array( grp_name, old_grp_number, new_grp_number )
    call hecmw_expand_index_array( grp_index, old_grp_number + 1, new_grp_number+1)
    call hecmw_expand_integer_array( grp_item, old_item_number, new_item_number )

    n_grp = new_grp_number

    j = old_grp_number + 1
    k = old_item_number + 1
    do i = 1, no_count
      write( grp_name_s, '(a,i0,a,i0)') 'FSTR_', grp_count, '_', i
      grp_name%s(j) = grp_name_s
      if( no_list(i) >= 0) then
        grp_item(k) = no_list(i)
        grp_index(j) = grp_index(j-1)+1
        k = k + 1
      else
        grp_index(j) = grp_index(j-1)
      endif
      j = j + 1
    end do
    grp_count = grp_count + 1
    call backset_group_pointers( hecMESH, grp_type_name )
    append_single_group = exist_n
  end function append_single_group

  subroutine append_new_group(hecMESH, grp_type_name, name, count, list, grp_id)
    implicit none
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    character(len=*), intent(in) :: grp_type_name
    character(len=HECMW_NAME_LEN), intent(in) :: name
    integer(kind=kint), intent(in) :: count
    integer(kind=kint), intent(in) :: list(:)
    integer(kind=kint), intent(out) :: grp_id
    integer(kind=kint) :: id, old_grp_number, new_grp_number, old_item_number, new_item_number, k

    call set_group_pointers( hecMESH, grp_type_name )
    do id = 1, n_grp
      if( hecmw_streqr(grp_name%s(id), name) ) then
        write(*,*) '### Error: Group already exists: ', name
        stop
      endif
    enddo

    old_grp_number = n_grp
    new_grp_number = old_grp_number + 1

    old_item_number = grp_index(n_grp)
    new_item_number = old_item_number + count

    call hecmw_expand_name_array( grp_name, old_grp_number, new_grp_number )
    call hecmw_expand_index_array( grp_index, old_grp_number + 1, new_grp_number + 1)
    call hecmw_expand_integer_array( grp_item, old_item_number, new_item_number )

    n_grp = new_grp_number
    grp_id = new_grp_number
    grp_name%s(grp_id) = name
    do k = 1, count
      grp_item(old_item_number + k) = list(k)
    enddo
    grp_index(grp_id) = grp_index(grp_id-1) + count
    call backset_group_pointers( hecMESH, grp_type_name )
  end subroutine append_new_group

  !------------------------------------------------------------------------------
  ! JP-0
  ! grp_type_name : 'node_grp', 'elem_grp' or 'surf_grp'
  ! name : group name
  ! return : number of member in specified group

  function get_grp_member_n( hecMESH, grp_type_name, name )
    implicit none
    integer(kind=kint) :: get_grp_member_n
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*) :: grp_type_name
    character(len=*) :: name
    integer(kind=kint) :: i

    call set_group_pointers( hecMESH, grp_type_name )

    do i = 1, n_grp
      if( hecmw_streqr(grp_name%s(i),name)) then
        get_grp_member_n = grp_index(i) - grp_index(i-1)
        return
      end if
    end do
    get_grp_member_n = 0
    return
  end function get_grp_member_n

  !------------------------------------------------------------------------------
  ! JP-1
  ! grp_type_name : 'node_grp', 'elem_grp' or 'surf_grp'
  ! name : group name
  ! return : number of member in specified group

  function get_grp_id( hecMESH, grp_type_name, name )
    implicit none
    integer(kind=kint) :: get_grp_id
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*) :: grp_type_name
    character(len=*) :: name
    integer(kind=kint) :: i

    call set_group_pointers( hecMESH, grp_type_name )

    do i = 1, n_grp
      if( hecmw_streqr(grp_name%s(i), name)) then
        get_grp_id = i
        return
      end if
    end do
    get_grp_id = 0
    return
  end function get_grp_id

  !------------------------------------------------------------------------------
  ! JP-2
  ! grp_type_name : 'node_grp', 'elem_grp' or 'surf_grp'
  ! name : group name
  ! member1 : id list for node or element
  ! member2 : id list for surface ( only 'surf_grp' specified )
  ! return : number of member in specified group

  function get_grp_member( hecMESH, grp_type_name, name, member1, member2 )
    implicit none
    integer(kind=kint) :: get_grp_member
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*) :: grp_type_name
    character(len=*) :: name
    integer(kind=kint),pointer :: member1(:)
    integer(kind=kint),pointer, optional :: member2(:)
    integer(kind=kint) :: i, j, k, sn, en

    get_grp_member = -1
    if( grp_type_name == 'surf_grp' .and. (.not. present( member2 ))) then
      stop 'assert in get_grp_member: not present member2 '
    end if

    call set_group_pointers( hecMESH, grp_type_name )

    do i = 1, n_grp
      if( hecmw_streqr(grp_name%s(i), name)) then
        sn = grp_index(i-1) + 1
        en = grp_index(i)
        k = 1
        if( grp_type == 3 ) then ! == surf_grp
          do j = sn, en
            member1(k) = grp_item(2*j-1)
            member2(k) = grp_item(2*j)
            k = k + 1
          end do
        else
          do j = sn, en
            member1(k) = grp_item(j)
            k = k + 1
          end do
        end if
        get_grp_member = en - sn + 1
        return
      end if
    end do
    get_grp_member = 0
    return
  end function get_grp_member

  subroutine node_grp_name_to_id( hecMESH, header_name, n, grp_id_name, grp_ID )
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    character(len=*)          :: header_name
    character(HECMW_NAME_LEN) :: grp_id_name(:)
    integer(kind=kint),pointer :: grp_ID(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: i, id
    character(len=256) :: msg

    do i = 1, n
      grp_ID(i) = -1
      do id = 1, hecMESH%node_group%n_grp
        if( hecmw_streqr(hecMESH%node_group%grp_name(id),grp_id_name(i))) then
          grp_ID(i) = id
          exit
        end if
      end do
      if( grp_ID(i) == -1 ) then
        write(msg,*) '### Error: ', header_name,' : Node group "',&
          grp_id_name(i),'" does not exist.'
        call hecmw_setup_util_err_stop(msg)
      end if
    end do
  end subroutine node_grp_name_to_id

  subroutine elem_grp_name_to_id( hecMESH, header_name, n, grp_id_name, grp_ID )
    implicit none
    type (hecmwST_local_mesh) :: hecMESH
    character(len=*)          :: header_name
    character(HECMW_NAME_LEN) :: grp_id_name(:)
    integer(kind=kint) :: grp_ID(:)
    integer(kind=kint) :: n
    integer(kind=kint) :: i, id
    character(len=256) :: msg

    do i = 1, n
      grp_ID(i) = -1
      do id = 1, hecMESH%elem_group%n_grp
        if (hecmw_streqr(hecMESH%elem_group%grp_name(id), grp_id_name(i))) then
          grp_ID(i) = id
          exit
        end if
      end do
      if( grp_ID(i) == -1 ) then
        write(msg,*) '### Error: ', header_name,' : Element group "',&
          grp_id_name(i),'" does not exist.'
        call hecmw_setup_util_err_stop(msg)
      end if
    end do
  end subroutine elem_grp_name_to_id

  !------------------------------------------------------------------------------
  ! JP-5
  ! JP-6
  !

  subroutine node_grp_name_to_id_ex( hecMESH, header_name, n, grp_id_name, grp_ID )
    implicit none
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*)                 :: header_name
    integer(kind=kint) :: n
    character(len=HECMW_NAME_LEN) :: grp_id_name(:)
    integer(kind=kint) :: grp_ID(:)

    integer(kind=kint) :: i, id
    integer(kind=kint) :: no, no_count, exist_n
    integer(kind=kint),pointer :: no_list(:)
    character(HECMW_NAME_LEN) :: name
    character(len=256) :: msg

    allocate( no_list( n ))
    no_count = 0
    do i = 1, n
      if( hecmw_str2index( grp_id_name(i), no )) then
        no_count = no_count + 1
        no_list(no_count) = no
        grp_ID(i) = hecMESH%node_group%n_grp + no_count
      else
        grp_ID(i) = -1
        do id = 1, hecMESH%node_group%n_grp
          if (hecmw_streqr(hecMESH%node_group%grp_name(id), grp_id_name(i))) then
            grp_ID(i) = id
            exit
          end if
        end do
        if( grp_ID(i) == -1 ) then
          write(msg,*) '### Error: ', header_name,' : Node group "',grp_id_name(i),'" does not exist.'
          call hecmw_setup_util_err_stop(msg)
        end if
      end if
    end do

    if( no_count > 0 ) then
      name = 'node_grp'
      exist_n = append_single_group( hecMESH, name, no_count, no_list )
      !    if( exist_n < no_count ) then
      !      write(*,*) '### Warning: ', header_name, ': following nodes are not exist'
      !      write(imsg,*) '### Warning: ', header_name, ': following nodes are not exist'
      !      do i=1, no_count
      !        if( no_list(i)<0 ) then
      !          write(*,*) -no_list(i)
      !          write(imsg,*) -no_list(i)
      !        end if
      !      end do
      !    end if
    end if

    deallocate( no_list )
  end subroutine node_grp_name_to_id_ex

  subroutine elem_grp_name_to_id_ex( hecMESH, header_name, n, grp_id_name, grp_ID )
    implicit none
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*)                 :: header_name
    integer(kind=kint) :: n
    character(HECMW_NAME_LEN) :: grp_id_name(:)
    integer(kind=kint) :: grp_ID(:)
    integer(kind=kint) :: i, id
    integer(kind=kint) :: no, no_count, exist_n
    integer(kind=kint),pointer :: no_list(:)
    character(HECMW_NAME_LEN) :: name
    character(len=256) :: msg

    allocate( no_list( n ))
    no_count = 0
    do i = 1, n
      if( hecmw_str2index( grp_id_name(i), no )) then
        no_count = no_count + 1
        no_list(no_count) = no
        grp_ID(i) = hecMESH%elem_group%n_grp + no_count
      else
        grp_ID(i) = -1
        do id = 1, hecMESH%elem_group%n_grp
          if (hecmw_streqr(hecMESH%elem_group%grp_name(id), grp_id_name(i))) then
            grp_ID(i) = id
            exit
          end if
        end do
        if( grp_ID(i) == -1 ) then
          write(msg,*) '### Error: ', header_name,' : Element group "',&
            grp_id_name(i),'" does not exist.'
          call hecmw_setup_util_err_stop(msg)
        end if
      end if
    end do

    if( no_count > 0 ) then
      name = 'elem_grp'
      exist_n = append_single_group( hecMESH, name, no_count, no_list )
      if( exist_n < no_count ) then
        write(*,*) '### Warning: ', header_name, ': following elements are not exist'
!        write(imsg,*) '### Warning: ', header_name, ': following elements are not exist'
        do i=1, no_count
          if( no_list(i)<0 ) then
            write(*,*) -no_list(i)
!            write(imsg,*) -no_list(i)
          end if
        end do
      end if
    end if

    deallocate( no_list )
  end subroutine elem_grp_name_to_id_ex

  !------------------------------------------------------------------------------

  subroutine surf_grp_name_to_id_ex( hecMESH, header_name, n, grp_id_name, grp_ID )
    implicit none
    type (hecmwST_local_mesh),target :: hecMESH
    character(len=*)                 :: header_name
    integer(kind=kint) :: n
    character(len=HECMW_NAME_LEN) :: grp_id_name(:)
    integer(kind=kint) :: grp_ID(:)
    integer(kind=kint) :: i, id
    character(len=256) :: msg

    do i = 1, n
      grp_ID(i) = -1
      do id = 1, hecMESH%surf_group%n_grp
        if (hecmw_streqr(hecMESH%surf_group%grp_name(id), grp_id_name(i))) then
          grp_ID(i) = id
          exit
        end if
      end do
      if( grp_ID(i) == -1 ) then
        write(msg,*) '### Error: ', header_name,' : Surface group "',grp_id_name(i),'" does not exist.'
        call hecmw_setup_util_err_stop(msg)
      end if
    end do
  end subroutine surf_grp_name_to_id_ex

  !------------------------------------------------------------------------------
  ! JP-8

  function get_node_grp_member_n( hecMESH, grp_name_array, n )
    implicit none
    integer(kind=kint) :: get_node_grp_member_n
    type (hecmwST_local_mesh), target :: hecMESH
    type(hecmw_str_arr) :: grp_name_array
    integer(kind=kint) :: n
    integer(kind=kint) :: i,j, m

    m = 0;
    do i = 1, n
      call set_group_pointers( hecMESH, grp_name_array%s(i) )
      do j = 1, n_grp
        if( hecmw_streqr(grp_name%s(j), grp_name_array%s(i))) then
          m = m + grp_index(j) - grp_index(j-1)
        end if
      end do
    end do
    get_node_grp_member_n = m
    return
  end function get_node_grp_member_n

  !------------------------------------------------------------------------------

  subroutine hecmw_expand_index_array( array, old_size, new_size )
    implicit none
    integer(kind=kint), pointer :: array(:)
    integer(kind=kint) :: old_size, new_size,i
    integer(kind=kint), pointer :: temp(:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(0:old_size-1))
      do i=0, old_size-1
        temp(i) = array(i)
      end do
      deallocate(array)
      allocate(array(0:new_size-1))
      array = 0
      do i=0, old_size-1
        array(i) = temp(i)
      end do
      deallocate(temp)
    else
      allocate(array(0:new_size-1))
      array = 0
    end if
  end subroutine hecmw_expand_index_array

  subroutine hecmw_expand_char_array( array, old_size, new_size )
    implicit none
    character(len=HECMW_NAME_LEN), pointer :: array(:)
    integer(kind=kint) :: old_size, new_size,i
    character(len=HECMW_NAME_LEN), pointer :: temp(:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(old_size))
      do i=1, old_size
        temp(i) = array(i)
      end do
      deallocate(array)
      allocate(array(new_size))
      array = ''
      do i=1, old_size
        array(i) = temp(i)
      end do
      deallocate(temp)
    else
      allocate(array(new_size))
      array = ''
    end if
  end subroutine hecmw_expand_char_array

  subroutine hecmw_expand_integer_array( array, old_size, new_size )
    implicit none
    integer(kind=kint), pointer :: array(:)
    integer(kind=kint) :: old_size, new_size,i
    integer(kind=kint), pointer :: temp(:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(old_size))
      do i=1, old_size
        temp(i) = array(i)
      end do
      deallocate(array)
      allocate(array(new_size))
      array = 0
      do i=1, old_size
        array(i) = temp(i)
      end do
      deallocate(temp)
    else
      allocate(array(new_size))
      array = 0
    end if
  end subroutine hecmw_expand_integer_array

  subroutine hecmw_expand_real_array( array, old_size, new_size )
    implicit none
    real(kind=kreal), pointer :: array(:)
    integer(kind=kint) :: old_size, new_size, i
    real(kind=kreal), pointer :: temp(:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(old_size))
      do i=1, old_size
        temp(i) = array(i)
      end do
      deallocate(array)
      allocate(array(new_size))
      array = 0
      do i=1, old_size
        array(i) = temp(i)
      end do
      deallocate(temp)
    else
      allocate(array(new_size))
      array = 0
    end if
  end subroutine hecmw_expand_real_array

  ! array( old_size, column ) -> array( new_size, column )
  subroutine hecmw_expand_integer_array2( array, column, old_size, new_size )
    implicit none
    integer(kind=kint), pointer :: array(:,:)
    integer(kind=kint) :: column, old_size, new_size, i,j
    integer(kind=kint), pointer :: temp(:,:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(old_size,column))
      do i=1, old_size
        do j=1,column
          temp(i,j) = array(i,j)
        end do
      end do
      deallocate(array)
      allocate(array(new_size,column))
      array = 0
      do i=1, old_size
        do j=1,column
          array(i,j) = temp(i,j)
        end do
      end do
      deallocate(temp)
    else
      allocate(array(new_size, column))
      array = 0
    end if
  end subroutine hecmw_expand_integer_array2


  ! array( old_size, column ) -> array( new_size, column )

  subroutine hecmw_expand_real_array2( array, column, old_size, new_size )
    implicit none
    real(kind=kreal), pointer :: array(:,:)
    integer(kind=kint) :: column, old_size, new_size, i,j
    real(kind=kreal), pointer :: temp(:,:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array ) ) then
      allocate(temp(old_size,column))
      do i=1, old_size
        do j=1,column
          temp(i,j) = array(i,j)
        end do
      end do
      deallocate(array)
      allocate(array(new_size,column))
      array = 0
      do i=1, old_size
        do j=1,column
          array(i,j) = temp(i,j)
        end do
      end do
      deallocate(temp)
    else
      allocate(array(new_size, column))
      array = 0
    end if
  end subroutine hecmw_expand_real_array2

  subroutine hecmw_expand_name_array( array, old_size, new_size )
    implicit none
    type(hecmw_str_arr) :: array
    integer(kind=kint) :: old_size, new_size, i
    character(len=HECMW_NAME_LEN), pointer :: temp(:)

    if( old_size >= new_size ) then
      return
    end if

    if( associated( array%s ) ) then
      allocate(temp(old_size))
      do i=1, old_size
        temp(i) = array%s(i)
      end do
      deallocate(array%s)
      allocate(array%s(new_size))
      do i=1, old_size
        array%s(i) = temp(i)
      end do
      deallocate(temp)
    else
      allocate(array%s(new_size))
    end if
  end subroutine hecmw_expand_name_array

  subroutine hecmw_delete_index_array( array, old_size, nindex )
    implicit none
    integer(kind=kint), pointer    :: array(:)  !< array to be modified
    integer(kind=kint), intent(in) :: old_size  !< current array size
    integer(kind=kint), intent(in) :: nindex    !< number of items to be deleted
    integer(kind=kint) :: i
    integer(kind=kint), pointer :: temp(:)

    if( old_size < nindex ) then
      return
    end if

    if( old_size == nindex ) then
      deallocate( array )
      return
    endif

    allocate(temp(0:old_size-1))
    do i=0, old_size-nindex-1
      temp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(0:old_size-nindex-1))
    array = 0
    do i=0, old_size-nindex-1
      array(i) = temp(i)
    end do
    deallocate(temp)
  end subroutine hecmw_delete_index_array

  subroutine hecmw_delete_integer_array( array, old_size, nitem )
    implicit none
    integer(kind=kint), pointer    :: array(:)  !< array to be modified
    integer(kind=kint), intent(in) :: old_size  !< current array size
    integer(kind=kint), intent(in) :: nitem     !< number of items to be deleted
    integer(kind=kint) :: i
    integer(kind=kint), pointer :: temp(:)

    if( old_size < nitem ) then
      return
    end if

    if( old_size == nitem ) then
      deallocate( array )
      return
    endif

    allocate(temp(old_size))
    do i=1, old_size-nitem
      temp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(old_size-nitem))
    array = 0
    do i=1, old_size-nitem
      array(i) = temp(i)
    end do
    deallocate(temp)
  end subroutine hecmw_delete_integer_array

  subroutine hecmw_delete_real_array( array, old_size, nitem )
    implicit none
    real(kind=kreal), pointer :: array(:)!< array to be modified
    integer(kind=kint), intent(in) :: old_size  !< current array size
    integer(kind=kint), intent(in) :: nitem     !< number of items to be deleted
    integer(kind=kint) :: i
    real(kind=kreal), pointer :: temp(:)

    if( old_size < nitem ) then
      return
    end if

    if( old_size == nitem ) then
      deallocate( array )
      return
    endif

    allocate(temp(old_size))
    do i=1, old_size-nitem
      temp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(old_size-nitem))
    array = 0
    do i=1, old_size-nitem
      array(i) = temp(i)
    end do
    deallocate(temp)
  end subroutine hecmw_delete_real_array

  !-----------------------------------------------------------------------------!

  subroutine reallocate_integer( array, n )
    implicit none
    integer(kind=kint),pointer :: array(:)
    integer(kind=kint) :: n;

    if( associated( array )) deallocate(array)
    allocate( array(n));
  end subroutine reallocate_integer

  subroutine reallocate_real( array, n )
    implicit none
    real(kind=kreal),pointer :: array(:)
    integer(kind=kint) :: n;

    if( associated( array )) deallocate(array)
    allocate( array(n));
  end subroutine reallocate_real

end module
