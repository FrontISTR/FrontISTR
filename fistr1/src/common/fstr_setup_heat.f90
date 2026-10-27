!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to append boundary conditions of
!! heat analysis to fstr_heat, independently of the input file format.
!! header_name is the card name shown in error messages. A blank name of
!! amplitude means that it is not given.
module m_fstr_setup_heat
  use m_fstr
  use hecmw_setup_util, only &
    : get_grp_member &
    , get_grp_member_n &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array &
    , hecmw_expand_integer_array2 &
    , hecmw_expand_real_array2
  use fstr_setup_util, only &
    : amp_name_to_id &
    , get_local_member_index &
    , get_sorted_local_member_index
  implicit none

contains

  !-----------------------------------------------------------------------------!
  !> Append !FIXTEMP                                                            !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_fixtemp( hecMESH, fstrPARAM, fstrHEAT, header_name, &
      amp_name, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_param), target                  :: fstrPARAM
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0
    do i = 1, n
      !rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      rtc = get_sorted_local_member_index( hecMESH, fstrPARAM, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'node_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-8
    old_size = fstrHEAT%T_FIX_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%T_FIX_node, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%T_FIX_ampl, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%T_FIX_val,  old_size, new_size )
    fstrHEAT%T_FIX_tot = new_size

    head = old_size + 1
    member => fstrHEAT%T_FIX_node(head:)
    id = head
    do i = 1, n
      !rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      rtc = get_sorted_local_member_index( hecMESH, fstrPARAM, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'node_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) then
        member => member( member_n+1 : )
      endif
      do j = 1, member_n
        fstrHEAT%T_FIX_val  (id) = vals(i)
        fstrHEAT%T_FIX_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_fixtemp

  !-----------------------------------------------------------------------------!
  !> Append !CFLUX                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_cflux( hecMESH, fstrHEAT, header_name, amp_name, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0

    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'node_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-9
    old_size = fstrHEAT%Q_NOD_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%Q_NOD_node, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_NOD_ampl, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%Q_NOD_val,  old_size, new_size )
    fstrHEAT%Q_NOD_tot = new_size

    head = old_size + 1
    member => fstrHEAT%Q_NOD_node(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'node_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%Q_NOD_val  (id) = vals(i)
        fstrHEAT%Q_NOD_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_cflux

  !-----------------------------------------------------------------------------!
  !> Append !DFLUX                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_dflux( hecMESH, fstrHEAT, header_name, amp_name, n, grp_names, &
      load_types, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: load_types(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'elem_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-10
    old_size = fstrHEAT%Q_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_elem, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_ampl, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_surf, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%Q_SUF_val,  old_size, new_size )
    fstrHEAT%Q_SUF_tot = new_size

    head = old_size + 1
    member => fstrHEAT%Q_SUF_elem(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'elem_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%Q_SUF_surf (id) = load_types(i)
        fstrHEAT%Q_SUF_val  (id) = vals(i)
        fstrHEAT%Q_SUF_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_dflux

  !-----------------------------------------------------------------------------!
  !> Append !SFLUX                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_sflux( hecMESH, fstrHEAT, header_name, amp_name, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member1(:), member2(:)

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0
    do i = 1, n
      m = m + get_grp_member_n( hecMESH, 'surf_grp', grp_names(i) )
    end do

    if (m == 0) return

    ! JP-11
    old_size = fstrHEAT%Q_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_elem, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_ampl, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_surf, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%Q_SUF_val,  old_size, new_size )
    fstrHEAT%Q_SUF_tot = new_size

    head = old_size + 1
    member1 => fstrHEAT%Q_SUF_elem(head:)
    member2 => fstrHEAT%Q_SUF_surf(head:)
    id = head
    do i = 1, n
      member_n = get_grp_member( hecMESH, 'surf_grp', grp_names(i), member1, member2 )
      if( i<n ) then
        member1 => member1( member_n+1 : )
        member2 => member2( member_n+1 : )
      end if
      do j = 1, member_n
        fstrHEAT%Q_SUF_val  (id) = vals(i)
        fstrHEAT%Q_SUF_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_sflux

  !-----------------------------------------------------------------------------!
  !> Append !FILM                                                               !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_film( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, load_types, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: load_types(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'elem_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-12
    old_size = fstrHEAT%H_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%H_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%H_SUF_val,  2, old_size, new_size )
    fstrHEAT%H_SUF_tot = new_size

    head = old_size + 1
    member => fstrHEAT%H_SUF_elem(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'elem_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%H_SUF_surf (id)   = load_types(i)
        fstrHEAT%H_SUF_val  (id,1) = vals(i)
        fstrHEAT%H_SUF_val  (id,2) = sinks(i)
        fstrHEAT%H_SUF_ampl (id,1) = amp_id1
        fstrHEAT%H_SUF_ampl (id,2) = amp_id2
        id= id + 1
      end do
    end do

  end subroutine fstr_heat_append_film

  !-----------------------------------------------------------------------------!
  !> Append !SFILM                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_sfilm( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member1(:), member2(:)

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      m = m + get_grp_member_n( hecMESH, 'surf_grp', grp_names(i) )
    end do

    if (m == 0) return

    ! JP-13
    old_size = fstrHEAT%H_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%H_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%H_SUF_val,  2, old_size, new_size )
    fstrHEAT%H_SUF_tot = new_size

    head = old_size + 1
    member1 => fstrHEAT%H_SUF_elem(head:)
    member2 => fstrHEAT%H_SUF_surf(head:)
    id = head
    do i = 1, n
      member_n = get_grp_member( hecMESH, 'surf_grp', grp_names(i), member1, member2 )
      if( i<n ) then
        member1 => member1( member_n+1 : )
        member2 => member2( member_n+1 : )
      end if
      do j = 1, member_n
        fstrHEAT%H_SUF_val  (id,1) = vals(i)
        fstrHEAT%H_SUF_val  (id,2) = sinks(i)
        fstrHEAT%H_SUF_ampl (id,1) = amp_id1
        fstrHEAT%H_SUF_ampl (id,2) = amp_id2
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_sfilm

  !-----------------------------------------------------------------------------!
  !> Append !RADIATE                                                            !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_radiate( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, load_types, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: load_types(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'elem_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-14
    old_size = fstrHEAT%R_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%R_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%R_SUF_val,  2, old_size, new_size )
    fstrHEAT%R_SUF_tot = new_size

    head = old_size + 1
    member => fstrHEAT%R_SUF_elem(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'elem_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%R_SUF_surf (id)   = load_types(i)
        fstrHEAT%R_SUF_val  (id,1) = vals(i)
        fstrHEAT%R_SUF_val  (id,2) = sinks(i)
        fstrHEAT%R_SUF_ampl (id,1) = amp_id1
        fstrHEAT%R_SUF_ampl (id,2) = amp_id2
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_radiate

  !-----------------------------------------------------------------------------!
  !> Append !SRADIATE                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_sradiate( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member1(:), member2(:)

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      m = m + get_grp_member_n( hecMESH, 'surf_grp', grp_names(i) )
    end do

    if (m == 0) return

    ! JP-15
    old_size = fstrHEAT%R_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%R_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%R_SUF_val,  2, old_size, new_size )
    fstrHEAT%R_SUF_tot = new_size

    head = old_size + 1
    member1 => fstrHEAT%R_SUF_elem(head:)
    member2 => fstrHEAT%R_SUF_surf(head:)
    id = head
    do i = 1, n
      member_n = get_grp_member( hecMESH, 'surf_grp', grp_names(i), member1, member2 )
      if( i<n ) then
        member1 => member1( member_n+1 : )
        member2 => member2( member_n+1 : )
      end if
      do j = 1, member_n
        fstrHEAT%R_SUF_val  (id,1) = vals(i)
        fstrHEAT%R_SUF_val  (id,2) = sinks(i)
        fstrHEAT%R_SUF_ampl (id,1) = amp_id1
        fstrHEAT%R_SUF_ampl (id,2) = amp_id2
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_sradiate

end module m_fstr_setup_heat
