!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to initialize fstr_heat and to read
!! the cards of heat analysis into it. Boundary conditions are appended
!! independently of the input file format.
!! header_name is the card name shown in error messages. A blank name of
!! amplitude means that it is not given.
module m_fstr_heat_setup
  use m_fstr
  use hecmw_setup_util, only &
    : hecmw_streqr &
    , get_grp_member &
    , get_grp_member_n &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array &
    , hecmw_expand_integer_array2 &
    , hecmw_expand_real_array2
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_ctrl_err_stop &
    , reallocate_integer &
    , reallocate_real &
    , amp_name_to_id &
    , get_local_member_index &
    , get_sorted_local_member_index
  use fstr_ctrl_heat
  use fstr_ctrl_util_f
  implicit none

contains

  !> Initial setting of heat analysis
  subroutine fstr_heat_init( fstrHEAT )
    implicit none
    type(fstr_heat) :: fstrHEAT

    fstrHEAT%STEPtot     = 0
    fstrHEAT%MATERIALtot = 0
    fstrHEAT%AMPLITUDEtot= 0
    fstrHEAT%T_FIX_tot   = 0
    fstrHEAT%Q_NOD_tot   = 0
    fstrHEAT%Q_VOL_tot   = 0
    fstrHEAT%Q_SUF_tot   = 0
    fstrHEAT%R_SUF_tot   = 0
    fstrHEAT%H_SUF_tot   = 0
    fstrHEAT%WL_tot      = 0
    fstrHEAT%beta        = -1.0d0
  end subroutine fstr_heat_init

  !-----------------------------------------------------------------------------!
  !> Read in !HEAT                                                             !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_HEAT( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    integer(kind=kint) :: n
    character(len=HECMW_NAME_LEN) :: mName
    integer(kind=kint) :: i

    n = fstr_ctrl_get_data_line_n( ctrl )

    if( n == 0 ) return

    call reallocate_real( P%PARAM%dtime, n)
    call reallocate_real( P%PARAM%etime, n)
    call reallocate_real( P%PARAM%dtmin, n)
    call reallocate_real( P%PARAM%delmax,n)
    call reallocate_integer( P%PARAM%itmax, n)
    call reallocate_real( P%PARAM%eps,   n)
    P%PARAM%analysis_n = n

    P%PARAM%dtime = 0
    P%PARAM%etime = 0
    P%PARAM%dtmin = 0
    P%PARAM%delmax = 0
    P%PARAM%itmax = 20
    P%PARAM%eps = 1.0e-6
    P%PARAM%timepoint_id = 0

    rcode = fstr_ctrl_get_HEAT(   ctrl,        &
      P%PARAM%dtime,     &
      P%PARAM%etime,     &
      P%PARAM%dtmin,     &
      P%PARAM%delmax,    &
      P%PARAM%itmax,     &
      P%PARAM%eps,       &
      mName,             &
      P%HEAT%beta)
    if( rcode /= 0 ) then
      call fstr_ctrl_err_stop
    end if

    if( associated(P%PARAM%timepoints) ) then
      do i=1,size(P%PARAM%timepoints)
        if( hecmw_streqr( P%PARAM%timepoints(i)%name, mName ) ) then
          P%PARAM%timepoint_id = i; exit
        endif
      enddo
    endif

    call reallocate_real( P%HEAT%STEP_DLTIME, n)
    call reallocate_real( P%HEAT%STEP_EETIME, n)
    call reallocate_real( P%HEAT%STEP_DELMIN, n)
    call reallocate_real( P%HEAT%STEP_DELMAX, n)
    P%HEAT%STEPtot = n

    P%HEAT%STEP_DLTIME = P%PARAM%dtime
    P%HEAT%STEP_EETIME = P%PARAM%etime
    P%HEAT%STEP_DELMIN = P%PARAM%dtmin
    P%HEAT%STEP_DELMAX = P%PARAM%delmax
    P%HEAT%timepoint_id = P%PARAM%timepoint_id

  end subroutine fstr_setup_HEAT

  !-----------------------------------------------------------------------------!
  !> Read in !FIXTEMP                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_FIXTEMP( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_FIXTEMP( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_fixtemp( P%MESH, P%PARAM, P%HEAT, '!FIXTEMP', amp, n, grp_id_name, value )

    deallocate( grp_id_name )
    deallocate( value )
  end subroutine fstr_setup_FIXTEMP

  !-----------------------------------------------------------------------------!
  !> Read in !CFLUX                                                            !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_CFLUX( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_CFLUX( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_cflux( P%MESH, P%HEAT, '!CFLUX', amp, n, grp_id_name, value )

    deallocate( grp_id_name )
    deallocate( value )
  end subroutine fstr_setup_CFLUX

  !-----------------------------------------------------------------------------!
  !> Read in !DFLUX                                                            !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_DFLUX( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( load_type(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_DFLUX( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, load_type, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_dflux( P%MESH, P%HEAT, '!DFLUX', amp, n, grp_id_name, load_type, value )

    deallocate( grp_id_name )
    deallocate( load_type )
    deallocate( value )
  end subroutine fstr_setup_DFLUX

  !-----------------------------------------------------------------------------!
  !> Read in !SFLUX                                                            !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_SFLUX( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_SFLUX( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_sflux( P%MESH, P%HEAT, '!SFLUX', amp, n, grp_id_name, value )

    deallocate( grp_id_name )
    deallocate( value )
  end subroutine fstr_setup_SFLUX

  !-----------------------------------------------------------------------------!
  !> Read in !FILM                                                             !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_FILM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( load_type(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '

    rcode = fstr_ctrl_get_FILM( ctrl, amp1, amp2, &
      grp_id_name, HECMW_NAME_LEN, load_type, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_film( P%MESH, P%HEAT, '!FILM', amp1, amp2, &
      n, grp_id_name, load_type, value, shink )

    deallocate( grp_id_name )
    deallocate( load_type )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_FILM

  !-----------------------------------------------------------------------------!
  !> Read in !SFILM                                                            !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_SFILM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '
    rcode = fstr_ctrl_get_SFILM( ctrl, amp1, amp2, &
      grp_id_name, HECMW_NAME_LEN, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_sfilm( P%MESH, P%HEAT, '!SFILM', amp1, amp2, &
      n, grp_id_name, value, shink )

    deallocate( grp_id_name )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_SFILM

  !-----------------------------------------------------------------------------!
  !> Read in !RADIATE                                                          !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_RADIATE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( load_type(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '
    rcode = fstr_ctrl_get_RADIATE( ctrl, amp1, amp2, &
      grp_id_name, HECMW_NAME_LEN, load_type, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_radiate( P%MESH, P%HEAT, '!RADIATE', amp1, amp2, &
      n, grp_id_name, load_type, value, shink )

    deallocate( grp_id_name )
    deallocate( load_type )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_RADIATE

  !-----------------------------------------------------------------------------!
  !> Read in !SRADIATE                                                         !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_SRADIATE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '
    rcode = fstr_ctrl_get_SRADIATE( ctrl, amp1, amp2, grp_id_name, HECMW_NAME_LEN, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_sradiate( P%MESH, P%HEAT, '!SRADIATE', amp1, amp2, &
      n, grp_id_name, value, shink )

    deallocate( grp_id_name )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_SRADIATE

  !-----------------------------------------------------------------------------!
  !> Read in !WELD_LINE                                                         !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_WELD_LINE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_heat), pointer :: fstrHEAT

    hecMESH => P%MESH
    fstrHEAT => P%HEAT

    fstrHEAT%WL_tot = fstrHEAT%WL_tot+1
    if( fstr_ctrl_get_WELDLINE( ctrl, hecMESH, HECMW_NAME_LEN, fstrHEAT%weldline(fstrHEAT%WL_tot) )/=0 ) then
      write(*,*) '### Error: Fail in read in Weld Line definition : ' , fstrHEAT%WL_tot
      write(ILOG,*) '### Error: Fail in read in Weld Line definition : ', fstrHEAT%WL_tot
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_WELD_LINE

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

end module m_fstr_heat_setup
