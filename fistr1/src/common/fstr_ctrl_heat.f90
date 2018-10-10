!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains control file data obtaining functions for heat conductive analysis
module fstr_ctrl_heat
  use m_fstr
  use hecmw

  include 'fstr_ctrl_util_f.inc'

  private :: pc_strupr

contains

  subroutine pc_strupr( s )
    implicit none
    character(*) :: s
    integer :: i, n, a, da

    n = len_trim(s)
    da = iachar('a') - iachar('A')
    do i = 1, n
      a = iachar(s(i:i))
      if( a > iachar('Z')) then
        a = a - da
        s(i:i) = achar(a)
      end if
    end do
  end subroutine pc_strupr



  !> Read in !HEAT
  function fstr_ctrl_get_HEAT( ctrl, dt, etime, dtmin, deltmx, itmax, eps )
    implicit none
    integer(kind=kint) :: ctrl
    real(kind=kreal),pointer :: dt(:)
    real(kind=kreal),pointer :: etime(:)
    real(kind=kreal),pointer :: dtmin(:)
    real(kind=kreal),pointer :: deltmx(:)
    integer(kind=kint),pointer :: itmax(:)
    real(kind=kreal),pointer :: eps(:)
    integer(kind=kint) :: fstr_ctrl_get_HEAT

    integer(kind=kint) :: result

    fstr_ctrl_get_HEAT = -1

    ! JP-7
    if( fstr_ctrl_get_data_array_ex( ctrl, 'rrrrir ', dt, etime, dtmin, deltmx, itmax, eps )/= 0) return

    fstr_ctrl_get_HEAT = 0
  end function fstr_ctrl_get_HEAT

  !> Read in !FIXTEMP
  function fstr_ctrl_get_FIXTEMP( ctrl, amp, node_grp_name, node_grp_name_len, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target :: node_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: node_grp_name_p
    integer(kind=kint) :: node_grp_name_len
    real(kind=kreal), pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_FIXTEMP

    character(len=HECMW_NAME_LEN) :: data_fmt,ss

    fstr_ctrl_get_FIXTEMP = -1

    ! JP-8
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',          0,     'S',   amp    )/= 0) return

    write(ss,*)  node_grp_name_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'r '

    node_grp_name_p => node_grp_name(1)
    fstr_ctrl_get_FIXTEMP = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_grp_name_p, value )

  end function fstr_ctrl_get_FIXTEMP


  !> Read in !CFLUX (heat)
  function fstr_ctrl_get_CFLUX( ctrl, amp, node_grp_name, node_grp_name_len, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target :: node_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: node_grp_name_p
    integer(kind=kint) :: node_grp_name_len
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_CFLUX

    character(len=HECMW_NAME_LEN) :: data_fmt,ss

    fstr_ctrl_get_CFLUX= -1

    ! JP-9
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',          0,     'S',   amp    )/= 0) return

    write(ss,*)  node_grp_name_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'r '

    node_grp_name_p => node_grp_name(1)
    fstr_ctrl_get_CFLUX = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_grp_name_p, value )

  end function fstr_ctrl_get_CFLUX

  !> Read in !DFLUX (heat)
  function fstr_ctrl_get_DFLUX( ctrl, amp, elem_grp_name, elem_grp_name_len, load_type, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target :: elem_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: elem_grp_name_p
    integer(kind=kint) :: elem_grp_name_len
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_DFLUX

    integer(kind=kint), parameter :: type_name_size = 5
    integer(kind=kint) :: i, n
    character(len=HECMW_NAME_LEN) :: data_fmt,s1,s2
    character(len=type_name_size),pointer :: type_name_list(:)
    character(len=type_name_size),pointer :: type_name_list_p
    integer(kind=kint) :: rcode
    integer(kind=kint) :: lid = -1

    fstr_ctrl_get_DFLUX = -1

    ! JP-10
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',          0,     'S',   amp    )/= 0) return

    write(s1,*)  elem_grp_name_len
    write(s2,*)  type_name_size
    write(data_fmt,'(a,a,a,a,a)') 'S',trim(adjustl(s1)),'S',trim(adjustl(s2)),'r '

    n = fstr_ctrl_get_data_line_n(ctrl)
    allocate( type_name_list(n) )

    elem_grp_name_p => elem_grp_name(1)
    type_name_list_p => type_name_list(1)
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, elem_grp_name_p, type_name_list_p, value)

      if( rcode /= 0 ) then
        deallocate( type_name_list )
        return
      end if

      do i = 1, n
        lid = -1;
        call pc_strupr( type_name_list(i) )
        if(      type_name_list(i)(1:2) == 'BF'  ) then; lid = 0
        else if( type_name_list(i)(1:2) == 'S0'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'S1'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'S2'  ) then; lid = 2
        else if( type_name_list(i)(1:2) == 'S3'  ) then; lid = 3
        else if( type_name_list(i)(1:2) == 'S4'  ) then; lid = 4
        else if( type_name_list(i)(1:2) == 'S5'  ) then; lid = 5
        else if( type_name_list(i)(1:2) == 'S6'  ) then; lid = 6
        end if
        if( lid < 0 ) then
          write(ILOG,*) 'Error : !DFLUX : Load  type ',type_name_list(i),' is unknown'
          deallocate( type_name_list )
          return
        end if
        load_type(i) = lid
      end do

      deallocate( type_name_list )
      fstr_ctrl_get_DFLUX = 0
  end function fstr_ctrl_get_DFLUX

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !SFLUX (heat)
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_SFLUX( ctrl, amp, surface_grp_name, surface_grp_name_len, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target :: surface_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: surface_grp_name_p
    integer(kind=kint) :: surface_grp_name_len
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_SFLUX

    character(len=HECMW_NAME_LEN) :: data_fmt,ss

    fstr_ctrl_get_SFLUX = -1

    ! JP-11
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',          0,     'S',   amp    )/= 0) return

    write(ss,*)  surface_grp_name_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'r '

    surface_grp_name_p => surface_grp_name(1)
    fstr_ctrl_get_SFLUX = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, surface_grp_name_p, value )
  end function fstr_ctrl_get_SFLUX


  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !FILM (heat)
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_FILM( ctrl, amp1, amp2, elem_grp_name, elem_grp_name_len, load_type, value, sink)
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp1
    character(len=HECMW_NAME_LEN) :: amp2
    character(len=HECMW_NAME_LEN),target :: elem_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: elem_grp_name_p
    integer(kind=kint) :: elem_grp_name_len
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: sink(:)
    integer(kind=kint) :: fstr_ctrl_get_FILM

    integer(kind=kint),parameter :: type_name_size = 5
    integer(kind=kint) :: i, n
    character(len=HECMW_NAME_LEN) :: data_fmt,s1,s2
    character(len=type_name_size),pointer :: type_name_list(:)
    character(len=type_name_size),pointer :: type_name_list_p
    integer(kind=kint) :: lid
    integer(kind=kint) :: rcode

    fstr_ctrl_get_FILM = -1

    ! JP-12
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP1 ',  '# ',         0,     'S',   amp1   )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP2 ',  '# ',         0,     'S',   amp2   )/= 0) return

    write(s1,*)  elem_grp_name_len
    write(s2,*)  type_name_size
    write(data_fmt,'(a,a,a,a,a)') 'S',trim(adjustl(s1)),'S',trim(adjustl(s2)),'Rr '

    n = fstr_ctrl_get_data_line_n(ctrl)
    allocate( type_name_list(n) )

    elem_grp_name_p => elem_grp_name(1)
    type_name_list_p => type_name_list(1)
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, elem_grp_name_p, type_name_list_p, value, sink)

      if( rcode /= 0 ) then
        deallocate( type_name_list )
        return
      end if

      do i = 1, n
        lid = -1;
        call pc_strupr( type_name_list(i) )
        if(      type_name_list(i)(1:2) == 'F0'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'F1'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'F2'  ) then; lid = 2
        else if( type_name_list(i)(1:2) == 'F3'  ) then; lid = 3
        else if( type_name_list(i)(1:2) == 'F4'  ) then; lid = 4
        else if( type_name_list(i)(1:2) == 'F5'  ) then; lid = 5
        else if( type_name_list(i)(1:2) == 'F6'  ) then; lid = 6
        end if
        if( lid < 0 ) then
          write(ILOG,*) 'Error : !FILM : Load  type ',type_name_list(i),' is unknown'
          deallocate( type_name_list )
          return
        end if
        load_type(i) = lid
      end do

      deallocate( type_name_list )

      fstr_ctrl_get_FILM = 0
  end function fstr_ctrl_get_FILM

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !SFILM  (heat)
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_SFILM( ctrl, amp1, amp2, surface_grp_name, surface_grp_name_len, value, sink)
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp1
    character(len=HECMW_NAME_LEN) :: amp2
    character(len=HECMW_NAME_LEN),target :: surface_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: surface_grp_name_p
    integer(kind=kint) :: surface_grp_name_len
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: sink(:)
    integer(kind=kint) :: fstr_ctrl_get_SFILM

    character(len=HECMW_NAME_LEN) :: data_fmt,ss

    fstr_ctrl_get_SFILM = -1

    ! JP-13
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP1 ',  '# ',          0,     'S',   amp1    )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP2 ',  '# ',          0,     'S',   amp2    )/= 0) return

    write(ss,*)  surface_grp_name_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'Rr '

    surface_grp_name_p => surface_grp_name(1)
    fstr_ctrl_get_SFILM = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, surface_grp_name_p, value, sink )
  end function fstr_ctrl_get_SFILM


  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !RADIATE (heat)
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_RADIATE( ctrl, amp1, amp2, elem_grp_name, elem_grp_name_len, load_type, value, sink)
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp1
    character(len=HECMW_NAME_LEN) :: amp2
    character(len=HECMW_NAME_LEN),target :: elem_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: elem_grp_name_p
    integer(kind=kint) :: elem_grp_name_len
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: sink(:)
    integer(kind=kint) :: fstr_ctrl_get_RADIATE

    integer(kind=kint),parameter :: type_name_size = 5
    integer(kind=kint) :: i, n
    character(len=HECMW_NAME_LEN) :: data_fmt,s1,s2
    character(len=type_name_size),pointer :: type_name_list(:)
    character(len=type_name_size),pointer :: type_name_list_p
    integer(kind=kint) :: lid
    integer(kind=kint) :: rcode

    fstr_ctrl_get_RADIATE = -1

    ! JP-14
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP1 ',  '# ',         0,     'S',   amp1   )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP2 ',  '# ',         0,     'S',   amp2   )/= 0) return

    write(s1,*)  elem_grp_name_len
    write(s2,*)  type_name_size
    write(data_fmt,'(a,a,a,a,a)') 'S',trim(adjustl(s1)),'S',trim(adjustl(s2)),'Rr '

    n = fstr_ctrl_get_data_line_n(ctrl)
    allocate( type_name_list(n) )

    elem_grp_name_p => elem_grp_name(1)
    type_name_list_p => type_name_list(1)
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, elem_grp_name_p ,type_name_list_p, value, sink)

      if( rcode /= 0 ) then
        deallocate( type_name_list )
        return
      end if

      do i = 1, n
        lid = -1;
        call pc_strupr( type_name_list(i) )
        if(      type_name_list(i)(1:2) == 'R0'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'R1'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'R2'  ) then; lid = 2
        else if( type_name_list(i)(1:2) == 'R3'  ) then; lid = 3
        else if( type_name_list(i)(1:2) == 'R4'  ) then; lid = 4
        else if( type_name_list(i)(1:2) == 'R5'  ) then; lid = 5
        else if( type_name_list(i)(1:2) == 'R6'  ) then; lid = 6
        end if
        if( lid < 0 ) then
          write(ILOG,*) 'Error : !RADIATE : Load  type ',type_name_list(i),' is unknown'
          deallocate( type_name_list )
          return
        end if
        load_type(i) = lid
      end do

      deallocate( type_name_list )
      fstr_ctrl_get_RADIATE = 0
  end function fstr_ctrl_get_RADIATE


  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !SRADIATE  (heat)
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_SRADIATE( ctrl, amp1, amp2, surface_grp_name, surface_grp_name_len, value, sink)
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp1
    character(len=HECMW_NAME_LEN) :: amp2
    character(len=HECMW_NAME_LEN),target :: surface_grp_name(:)
    character(len=HECMW_NAME_LEN),pointer:: surface_grp_name_p
    integer(kind=kint) :: surface_grp_name_len
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: sink(:)
    integer(kind=kint) :: fstr_ctrl_get_SRADIATE

    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=HECMW_NAME_LEN) :: s1

    fstr_ctrl_get_SRADIATE = -1

    ! JP-15
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP1 ',  '# ',          0,     'S',   amp1    )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP2 ',  '# ',          0,     'S',   amp2    )/= 0) return

    write(s1,*) surface_grp_name_len;
    write(data_fmt,'(a,a,a)') 'S', trim(adjustl(s1)), 'Rr  '

    surface_grp_name_p => surface_grp_name(1)
    fstr_ctrl_get_SRADIATE = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, surface_grp_name_p, value, sink )

  end function fstr_ctrl_get_SRADIATE

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !WELD_LINE  (heat)
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_WELDLINE( ctrl, hecMESH, grp_name_len, weldline )
    use fstr_setup_util
    implicit none
    integer(kind=kint), intent(in)       :: ctrl
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint), intent(in)       :: grp_name_len
    type(tWeldLine), intent(inout)       :: weldline
    integer(kind=kint) :: fstr_ctrl_get_WELDLINE

    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=HECMW_NAME_LEN) :: s1, grp_id_name(1)
    integer :: grp_id(1)

    fstr_ctrl_get_WELDLINE = -1
    if( fstr_ctrl_get_data_ex( ctrl, 1, 'RRRR ',   weldline%I, weldline%U, weldline%coe, weldline%v )/=0 ) return
    write(s1,*) grp_name_len
    write(data_fmt,'(a,a,a)') 'S', trim(adjustl(s1)), 'IRRRR  '
    if( fstr_ctrl_get_data_ex( ctrl, 2, data_fmt,  grp_id_name, weldline%xyz, weldline%n1, &
      weldline%n2, weldline%distol, weldline%tstart )/=0 ) return
    call elem_grp_name_to_id( hecMESH, 'WELD_LINE ', 1, grp_id_name, grp_id )
    weldline%egrpid = grp_id(1)

    fstr_ctrl_get_WELDLINE = 0
  end function fstr_ctrl_get_WELDLINE

  !* ----------------------------------------------------------------------------------------------- *!
end module fstr_ctrl_heat
