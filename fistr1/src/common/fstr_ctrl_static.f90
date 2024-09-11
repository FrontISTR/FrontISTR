!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains control file data obtaining functions for static analysis
module fstr_ctrl_static
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

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !STATIC
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_STATIC( ctrl, &
      & dtime, etime, itime, eps, restart_nout, &
      & idx_elpl, &
      & iout_list, &
      & sig_y0, h_dash, &
      & nout, nout_monit, node_monit_1, elem_monit_1, intg_monit_1 )
    implicit none
    integer(kind=kint) :: ctrl
    real(kind=kreal)   :: dtime
    real(kind=kreal)   :: etime
    integer(kind=kint) :: itime
    real(kind=kreal)   :: eps
    integer(kind=kint) :: restart_nout
    integer(kind=kint) :: idx_elpl
    real(kind=kreal)   :: sig_y0, h_dash
    integer(kind=kint) :: nout, nout_monit, node_monit_1, elem_monit_1, intg_monit_1
    integer(kind=kint) :: iout_list(6)
    integer(kind=kint) :: fstr_ctrl_get_STATIC

    fstr_ctrl_get_STATIC = -1

    if( fstr_ctrl_get_data_ex( ctrl, 1, 'rriri ', dtime, etime, itime, eps, restart_nout ) /= 0 ) return
    if( fstr_ctrl_get_data_ex( ctrl, 2, 'i ', idx_elpl ) /= 0 ) return
    if( fstr_ctrl_get_data_ex( ctrl, 3, 'iiiiii ', &
      & iout_list(1), iout_list(2), iout_list(3), iout_list(4), iout_list(5), iout_list(6)) /= 0 ) return
    if( fstr_ctrl_get_data_ex( ctrl, 4, 'rr ', sig_y0, h_dash ) /= 0 ) return
    if( fstr_ctrl_get_data_ex( ctrl, 5, 'iiiii ', &
      & nout, nout_monit, node_monit_1, elem_monit_1, intg_monit_1 ) /= 0 ) return

    fstr_ctrl_get_STATIC = 0
  end function fstr_ctrl_get_STATIC

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !BOUNDARY
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_BOUNDARY( ctrl, amp, node_id, node_id_len, dof_ids, dof_ide, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target  :: node_id(:)
    character(len=HECMW_NAME_LEN),pointer :: node_id_p
    integer(kind=kint) :: node_id_len
    integer(kind=kint),pointer :: dof_ids (:)
    integer(kind=kint),pointer :: dof_ide (:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_BOUNDARY

    character(len=HECMW_NAME_LEN) :: data_fmt,ss
    write(ss,*)  node_id_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'IIr '

    fstr_ctrl_get_BOUNDARY = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
    node_id_p => node_id(1)
    fstr_ctrl_get_BOUNDARY = &
      fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_p, dof_ids, dof_ide, value )

  end function fstr_ctrl_get_BOUNDARY


  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !CLOAD
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_CLOAD( ctrl, amp, node_id, node_id_len, dof_id, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target  :: node_id(:)
    character(len=HECMW_NAME_LEN),pointer :: node_id_p
    integer(kind=kint) :: node_id_len
    integer(kind=kint),pointer :: dof_id(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_CLOAD

    character(len=HECMW_NAME_LEN) :: data_fmt,ss
    write(ss,*)  node_id_len
    write( data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'IR '

    fstr_ctrl_get_CLOAD = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
    node_id_p => node_id(1)
    fstr_ctrl_get_CLOAD = &
      fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_p, dof_id, value )

  end function fstr_ctrl_get_CLOAD

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !DLOAD
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_DLOAD( ctrl, amp, follow, element_id, element_id_len, load_type, params )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    integer(kind=kint) :: follow
    character(len=HECMW_NAME_LEN),target :: element_id(:)
    integer(kind=kint) :: element_id_len
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: params(:,:)
    integer(kind=kint) :: fstr_ctrl_get_DLOAD

    character(len=HECMW_NAME_LEN),pointer :: type_name_list(:)
    character(len=HECMW_NAME_LEN),pointer :: type_name_list_p
    character(len=HECMW_NAME_LEN),pointer :: element_id_p

    integer(kind=kint) :: i, n
    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: data_fmt,s1,s2
    integer(kind=kint) :: lid

    fstr_ctrl_get_DLOAD = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
    follow = follow+1
    if( fstr_ctrl_get_param_ex( ctrl, 'FOLLOW ','NO,YES ', 0, 'P', follow ) /= 0) return
    follow = follow-1

    write(s1,*)  element_id_len
    write(s2,*)  HECMW_NAME_LEN
    write( data_fmt, '(a,a,a,a,a)') 'S', trim(adjustl(s1)), 'S', trim(adjustl(s2)),'Rrrrrrr '

    n = fstr_ctrl_get_data_line_n(ctrl)
    allocate( type_name_list(n) )
    !!
    !! for avoiding stack overflow with intel 9 compiler
    !!
    element_id_p => element_id(1)
    type_name_list_p => type_name_list(1)

      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, element_id_p, type_name_list_p, &
        params(0,:), params(1,:), params(2,:), params(3,:), params(4,:),params(5,:), &
        params(6,:) )

      if( rcode /= 0 ) then
        deallocate( type_name_list )
        return
      end if

      do i=1, n
        call pc_strupr( type_name_list(i) )
        lid = -1;
        if(      type_name_list(i)(1:2) == 'BX'  ) then; lid = 1
        else if( type_name_list(i)(1:2) == 'BY'  ) then; lid = 2
        else if( type_name_list(i)(1:2) == 'BZ'  ) then; lid = 3
        else if( type_name_list(i)(1:4) == 'GRAV') then; lid = 4
        else if( type_name_list(i)(1:4) == 'CENT') then; lid = 5
        else if( type_name_list(i)(1:2) == 'PP'  ) then; lid = 10
        else if( type_name_list(i)(1:2) == 'P0'  ) then; lid = 10
        else if( type_name_list(i)(1:2) == 'PX'  ) then
          lid = 10; params(1,:)=1.d0; params(2,:)=0.d0; params(3,:)=0.d0
        else if( type_name_list(i)(1:2) == 'PY'  ) then
          lid = 10; params(1,:)=0.d0; params(2,:)=1.d0; params(3,:)=0.d0
        else if( type_name_list(i)(1:2) == 'PZ'  ) then
          lid = 10; params(1,:)=0.d0; params(2,:)=0.d0; params(3,:)=1.d0
        else if( type_name_list(i)(1:2) == 'P1'  ) then; lid = 10
        else if( type_name_list(i)(1:2) == 'P2'  ) then; lid = 20
        else if( type_name_list(i)(1:2) == 'P3'  ) then; lid = 30
        else if( type_name_list(i)(1:2) == 'P4'  ) then; lid = 40
        else if( type_name_list(i)(1:2) == 'P5'  ) then; lid = 50
        else if( type_name_list(i)(1:2) == 'P6'  ) then; lid = 60
        else if( type_name_list(i)(1:1) == 'S'   ) then; lid = 100
        else
          write(ILOG, *) 'Error : !DLOAD : Load  type ',type_name_list(i), ' is unknown'
          deallocate( type_name_list )
          return
        end if
        load_type(i) = lid
      end do

      deallocate( type_name_list )
      fstr_ctrl_get_DLOAD = 0

  end function fstr_ctrl_get_DLOAD



  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !REFTEMP
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_REFTEMP( ctrl, value )
    implicit none
    integer(kind=kint) :: ctrl
    real(kind=kreal)   :: value
    integer(kind=kint) :: fstr_ctrl_get_REFTEMP,rcode

    rcode = fstr_ctrl_get_data_array_ex( ctrl, 'r ', value )
    fstr_ctrl_get_REFTEMP = rcode

  end function fstr_ctrl_get_REFTEMP

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !TEMPERATURE
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_TEMPERATURE( ctrl, irres, tstep, tintl, rtype, node_id, node_id_len, value )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: irres
    integer(kind=kint) :: tstep
    integer(kind=kint) :: tintl
    integer(kind=kint) :: rtype
    character(len=HECMW_NAME_LEN), target :: node_id(:)
    character(len=HECMW_NAME_LEN), pointer:: node_id_p
    integer(kind=kint) :: node_id_len
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_TEMPERATURE, rcode

    character(len=HECMW_NAME_LEN) :: data_fmt,ss

    irres = 0
    if( fstr_ctrl_get_param_ex( ctrl, 'READRESULT ', '# ', 0, 'I', irres )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'SSTEP ',      '# ', 0, 'I', tstep )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'INTERVAL ',   '# ', 0, 'I', tintl )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'READTYPE ',   'STEP,TIME ', 0, 'P', rtype )/= 0) return
    if( irres > 0 ) then
      fstr_ctrl_get_TEMPERATURE = 0
      return
    endif

    write(ss,*)  node_id_len
    write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'r '

    node_id_p => node_id(1)
    rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_p, value )
    fstr_ctrl_get_TEMPERATURE = rcode

  end function fstr_ctrl_get_TEMPERATURE


  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !SPRING
  !* ----------------------------------------------------------------------------------------------- *!

  function fstr_ctrl_get_SPRING( ctrl, amp, node_id, node_id_len, dof_id, value )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=HECMW_NAME_LEN) :: amp
    character(len=HECMW_NAME_LEN),target  :: node_id(:)
    character(len=HECMW_NAME_LEN),pointer :: node_id_p
    integer(kind=kint) :: node_id_len
    integer(kind=kint),pointer :: dof_id(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: fstr_ctrl_get_SPRING

    character(len=HECMW_NAME_LEN) :: data_fmt,ss
    write(ss,*)  node_id_len
    write( data_fmt, '(a,a,a)') 'S', trim(adjustl(ss)), 'IR '

    fstr_ctrl_get_SPRING = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'AMP ',  '# ',  0, 'S', amp )/= 0) return
    node_id_p => node_id(1)
    fstr_ctrl_get_SPRING = &
      fstr_ctrl_get_data_array_ex( ctrl, data_fmt, node_id_p, dof_id, value )

  end function fstr_ctrl_get_SPRING


  !----------------------------------------------------------------------
  !> Read in !ULOAD
  integer function fstr_ctrl_get_USERLOAD( ctrl )
    use mULoad
    integer(kind=kint), intent(in)    :: ctrl

    character(len=256) :: fname

    fstr_ctrl_get_USERLOAD = -1
    if( fstr_ctrl_get_param_ex( ctrl, 'FILE  ', '# ',           0,   'S',   fname )/=0 ) return
    if( fname=="" ) stop "You must define a file name before read in user-defined material"
    if( ureadload(fname)/=0 ) return

    fstr_ctrl_get_USERLOAD = 0
  end function fstr_ctrl_get_USERLOAD

end module fstr_ctrl_static




