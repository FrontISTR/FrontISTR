!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to read the cards that go into fstr_param, and the cards that go into no
!! structure (the variables of m_fstr, the amplitudes of the mesh and the state of reading the cards).
module m_fstr_param_setup
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_card_slot &
    , fstr_ctrl_err_stop &
    , append_new_amplitude
  use fstr_ctrl_common
  use fstr_ctrl_static
  use mContactParam
  use m_utilities
  use fstr_ctrl_util_f
  implicit none

contains

  !-----------------------------------------------------------------------------!
  !> Read in !VERSION                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_VERSION( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_data_ex( ctrl, 1, 'i ', P%STATE%version )

  end subroutine fstr_setup_VERSION

  !-----------------------------------------------------------------------------!
  !> Read in !SOLUTION                                                         !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_SOLUTION( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_SOLUTION( ctrl, P%PARAM%solution_type, P%PARAM%nlgeom )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_SOLUTION

  !-----------------------------------------------------------------------------!
  !> Read in !NONLINEAR_SOLVER                                                         !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_NONLINEAR_SOLVER( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_NONLINEAR_SOLVER( ctrl, P%PARAM%nlsolver_method )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_NONLINEAR_SOLVER

  !-----------------------------------------------------------------------------!
  !> Read in !SOLVER                                                           !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_SOLVER( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    if( slot%first > 1 ) then
      write(ILOG,*) '### Error : !SOLVER exists twice in FSTR control file.'
      stop HECMW_EXIT_INPUT
    endif

    !   nier       => svIarray(1)
    !   method     => svIarray(2)
    !   precond    => svIarray(3)
    !   nset       => svIarray(4)
    !   iterpremax => svIarray(5)
    !   nrest      => svIarray(6)
    !   scaling    => svIarray(7)
    !   iterlog    => svIarray(21)
    !   timelog    => svIarray(22)
    !   steplog    => svIarray(23)
    !   loglevel   => svIarray(24)
    !   dumptype   => svIarray(31)
    !   dumpexit   => svIarray(32)
    !   usejad     => svIarray(33)
    !   ncolor_in  => svIarray(34)
    !   mpc_method => svIarray(13)
    !   estcond    => svIarray(14)
    !   contact_elim=> svIarray(15)
    !   recompute_residual => svIarray(16)
    !   method2    => svIarray(8)
    !   recyclepre => svIarray(35)
    !   matvec_impl=> svIarray(36)
    !   precond_impl=> svIarray(37)
    !   solver_opt => svIarray(41:50)
    !   nBFGS      => svIarray(60)

    !   resid      => svRarray(1)
    !   sigma_diag => svRarray(2)
    !   omega      => svRarray(3)
    !   thresh     => svRarray(4)
    !   filter     => svRarray(5)
    !   solver_ropt=> svRarray(41:50)

    rcode = fstr_ctrl_get_SOLVER( ctrl,                      &
      p%MAT%Iarray(2), p%MAT%Iarray(3), p%MAT%Iarray(4), &
      p%MAT%Iarray(21), p%MAT%Iarray(22), p%MAT%Iarray(23), &
      p%MAT%Iarray(1), p%MAT%Iarray(5), p%MAT%Iarray(6), p%MAT%Iarray(60), p%MAT%Iarray(7), &
      p%MAT%Iarray(31), p%MAT%Iarray(32), p%MAT%Iarray(33), p%MAT%Iarray(34), &
      p%MAT%Iarray(13), p%MAT%Iarray(14), p%MAT%Iarray(8), &
      p%MAT%Iarray(35), p%MAT%Iarray(41:50), p%MAT%Iarray(15), &
      p%MAT%Rarray(1), p%MAT%Rarray(2), p%MAT%Rarray(3), &
      p%MAT%Rarray(4), p%MAT%Rarray(5), p%MAT%Rarray(41:50), p%MAT%Iarray(24), &
      p%MAT%Iarray(36), p%MAT%Iarray(37), p%MAT%Iarray(16) )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    if( p%MAT%Iarray(2) <= 100 ) then
      p%MAT%Iarray(99) = 1  ! indirect method
    else
      p%MAT%Iarray(99) = p%MAT%Iarray(2)-99 !2  ! direct method
    end if

  end subroutine fstr_setup_SOLVER

  !* ----------------------------------------------------------------------------------------------- *!
  !> Read in !ORIENTATION
  !* ----------------------------------------------------------------------------------------------- *!

  subroutine fstr_setup_ORIENTATION( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    if( fstr_get_ORIENTATION( ctrl, P%MESH, slot%first, g_LocalCoordSys(slot%first) )/=0 ) then
      write(*,*) '### Error: Fail in read in ORIENTATION definition : ', slot%first
      write(ILOG,*) '### Error: Fail in read in ORIENTATION definition : ', slot%first
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_ORIENTATION

  !> Reads !ORIENTATION into coordsys. Returns 0, or -1 for an error.
  integer function fstr_get_ORIENTATION( ctrl, hecMESH, cnt, coordsys )
    use hecmw_setup_util, only: node_global_to_local
    implicit none
    integer(kind=kint)         :: ctrl
    type( hecmwST_local_mesh ) :: hecMESH
    integer                    :: cnt
    type( tLocalCoordSys )     :: coordsys

    character(len=HECMW_NAME_LEN) :: grp_id_name(1)

    integer :: nid, dtype
    character(len=HECMW_NAME_LEN) :: data_fmt
    real(kind=kreal) :: fdum, xyza(3), xyzb(3), xyzc(3), ff1(3), ff2(3), ff3(3)

    fstr_get_ORIENTATION = -1

    nid = 1
    coordsys%sys_type = 10

    nid = 1
    data_fmt = 'COORDINATES,NODES,LOCAL_NODES '
    if( fstr_ctrl_get_param_ex( ctrl, 'DEFINITION ', data_fmt, 0, 'P', nid )/=0 ) return
    dtype = nid-1
    coordsys%sys_type = coordsys%sys_type + dtype

    if( fstr_ctrl_get_param_ex( ctrl, 'NAME ',  '# ',  1, 'S', grp_id_name(1) )/= 0) return
    coordsys%sys_name = grp_id_name(1)

    if( dtype==0 ) then
      data_fmt = "RRRRRRrrr "
      xyzc(:) = 0.d0
      if( fstr_ctrl_get_data_ex( ctrl, 1, data_fmt, xyza(1), xyza(2),  &
        xyza(3), xyzb(1), xyzb(2), xyzb(3), xyzc(1), xyzc(2), xyzc(3) )/=0 ) return
      if( coordsys%sys_type==10 ) then
        ff1 = xyza-xyzc
        fdum = dsqrt( dot_product(ff1, ff1) )
        if( fdum==0.d0 ) return
        ff1 = ff1/fdum
        ff2 = xyzb-xyzc
        call cross_product(ff1,ff2,ff3)
        coordsys%CoordSys(1,:) = ff1

        fdum = dsqrt( dot_product(ff3, ff3) )
        if( fdum==0.d0 ) return
        coordsys%CoordSys(3,:) = ff3/fdum

        call cross_product(coordsys%CoordSys(3,:), coordsys%CoordSys(1,:), coordsys%CoordSys(2,:) )
      else
        coordsys%CoordSys(1,:) = xyza
        coordsys%CoordSys(2,:) = xyzb
      endif

    else
      coordsys%node_ID(3) = 0   ! global origin
      data_fmt = "IIi "
      if( fstr_ctrl_get_data_ex( ctrl, 1, data_fmt, coordsys%node_ID(1),  &
        coordsys%node_ID(2), coordsys%node_ID(3) )/=0 ) return
      ! DEFINITION=LOCAL_NODES: three local node numbers (1..number of nodes of
      ! the element) defining the coordinate system of each element. All three
      ! are mandatory here, unlike the third node of DEFINITION=NODES.
      if( dtype==2 ) then
        if( maxval(coordsys%node_ID(1:3)) > 10 .or. minval(coordsys%node_ID(1:3)) < 1 ) then
          write(*,*) "!ORIENTATION, DEFINITION=LOCAL_NODES needs three local node numbers (1-10)!"
          write(IDBG,*) "!ORIENTATION, DEFINITION=LOCAL_NODES needs three local node numbers (1-10)!"
          return
        endif
        fstr_get_ORIENTATION = 0
        return
      endif
      if( coordsys%node_ID(3) == 0 ) then
        nid = node_global_to_local( hecMESH, coordsys%node_ID(1:2), 2 )
        if( nid/=0 .and. nid/=2 ) then
          write(*,*) "We cannot define coordinate system using nodes in other CPU!"
          write(IDBG,*) "We cannot define coordinate system using nodes in other CPU!"
          return
        endif
      else
        nid = node_global_to_local( hecMESH, coordsys%node_ID, 3 )
        if( nid/=0 .and. nid/=3 ) then
          write(*,*) "We cannot define coordinate system using nodes in other CPU!"
          write(IDBG,*) "We cannot define coordinate system using nodes in other CPU!"
          return
        endif
      endif
    endif

    fstr_get_ORIENTATION = 0
  end function fstr_get_ORIENTATION

  !> Read in !INITIAL_CONDITION
  subroutine fstr_setup_INITIAL_CONDITION( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    if( fstr_get_INITIAL_CONDITION( ctrl, g_InitialCnd(slot%first), P%MESH )/=0 ) then
       write(*,*) '### Error: Fail in read in INITIAL CONDITION definition : ' ,slot%first
       write(ILOG,*) '### Error: Fail in read in INITIAL CONDITION definition : ', slot%first
       stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_INITIAL_CONDITION

  !> Reads !INITIAL_CONDITION into cond. Returns 0, or nonzero for an error.
  integer(kind=kint) function fstr_get_INITIAL_CONDITION( ctrl, cond, hecMESH )
    implicit none
    integer(kind=kint)        :: ctrl
    type( tInitialCondition ) :: cond
    type(hecmwST_local_mesh)  :: hecMESH
    integer, pointer          :: grp_id(:), dof(:)
    real(kind=kreal), pointer :: temp(:)
    character(len=HECMW_NAME_LEN), pointer :: grp_id_name(:)
    character(len=HECMW_NAME_LEN) :: data_fmt, ss
    integer :: i,j,n, iS, iE, gid, nid, rcode

    fstr_get_INITIAL_CONDITION = -1

    ss = 'TEMPERATURE,VELOCITY,ACCELERATION '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', ss, 1, 'P', nid )
    if( nid==1 ) then
      cond%cond_name = "temperature"
      allocate( cond%intval(hecMESH%n_node) )
      allocate( cond%realval(hecMESH%n_node) )
    elseif( nid==2 ) then
      cond%cond_name = "velocity"
      allocate( cond%intval(hecMESH%n_node) )
      allocate( cond%realval(hecMESH%n_node) )
    elseif( nid==3 ) then
      cond%cond_name = "acceleration"
      allocate( cond%intval(hecMESH%n_node) )
      allocate( cond%realval(hecMESH%n_node) )
    else
      return
    endif

    cond%intval = -1
    cond%realval = 0.d0

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n<=0 ) return
    allocate( temp(n), grp_id_name(n), grp_id(n), dof(n) )
    dof = 0
    write(ss,*)  HECMW_NAME_LEN
    if( nid==1 ) then
      write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'R '
      fstr_get_INITIAL_CONDITION = &
            fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_id_name, temp )
    else
        write(data_fmt,'(a,a,a)') 'S',trim(adjustl(ss)),'IR '
        fstr_get_INITIAL_CONDITION = &
            fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_id_name, dof, temp )
    endif

    if( fstr_get_INITIAL_CONDITION /= 0 ) then
        if( associated(grp_id) )    deallocate( grp_id )
        if( associated(temp) )    deallocate( temp )
        if( associated(dof) )       deallocate( dof )
        if( associated(grp_id_name) )    deallocate( grp_id_name )
        return
    end if

    call node_grp_name_to_id_ex( hecMESH, '!INITIAL CONDITION', n, grp_id_name, grp_id )
    do i=1,n
      gid = grp_id(i)
      iS = hecMESH%node_group%grp_index(gid-1) + 1
      iE = hecMESH%node_group%grp_index(gid  )
      do j=iS, iE
        nid = hecMESH%node_group%grp_item(j)
        cond%realval(nid) = temp(i)
        cond%intval(nid) = dof(i)
      enddo
    enddo

    if( associated(grp_id) )    deallocate( grp_id )
    if( associated(temp) )      deallocate( temp )
    if( associated(dof) )       deallocate( dof )
    if( associated(grp_id_name) )    deallocate( grp_id_name )
end function fstr_get_INITIAL_CONDITION

  !-----------------------------------------------------------------------------!
  !> Read in !WRITE                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_WRITE( ctrl, counter, P )
    implicit none
    integer(kind=kint) :: ctrl
    integer(kind=kint) :: counter
    type(fstr_param_pack) :: P
    integer(kind=kint) :: res, visual, neutral

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_WRITE( ctrl, res, visual, neutral )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( res == 1 ) P%PARAM%fg_result = 1
    if( visual == 1 ) P%PARAM%fg_visual = 1
    if( neutral == 1 ) P%PARAM%fg_neutral = 1

  end subroutine fstr_setup_WRITE

  !-----------------------------------------------------------------------------!
  !> Read in !ECHO                                                            !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_ECHO( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_ECHO( ctrl,        &
      P%PARAM%fg_echo )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_ECHO

  !-----------------------------------------------------------------------------!
  !> Read in !RESTART                                                         !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_RESTART( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, nout
    nout = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'FREQUENCY ', '# ', 0, 'I', nout )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    rcode = fstr_ctrl_get_param_ex( ctrl, 'VERSION ', '# ', 0, 'I', P%PARAM%restart_version )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    P%SOLID%restart_nout= nout
    P%DYN%restart_nout= nout
    P%HEAT%restart_nout= nout

  end subroutine fstr_setup_RESTART

  !-----------------------------------------------------------------------------!
  !> Read in !AMPLITUDE                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_AMPLITUDE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    real(kind=kreal), pointer :: val(:), table(:)
    character(len=HECMW_NAME_LEN) :: name
    integer :: nline, n, type_def, type_time, type_val, rcode

    nline = fstr_ctrl_get_data_line_n( ctrl )
    if( nline<=0 ) return
    allocate( val(nline*4) )
    allocate( table(nline*4) )
    rcode = fstr_ctrl_get_AMPLITUDE( ctrl, nline, name, type_def, type_time, type_val, &
        n, val, table )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call append_new_amplitude( P%MESH%amp, name, type_def, type_time, type_val, n, val, table )

    if( associated(val) ) deallocate( val )
    if( associated(table) ) deallocate( table )
  end subroutine fstr_setup_AMPLITUDE

  !-----------------------------------------------------------------------------!
  !> Read in !REFTEMP                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_REFTEMP( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode

    rcode = fstr_ctrl_get_REFTEMP( ctrl, P%PARAM%ref_temp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_REFTEMP

  !-----------------------------------------------------------------------------!
  !> Read in !MPC                                                              !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_MPC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    !        integer(kind=kint) :: type
    !        integer(kind=kint),pointer :: node1_ptr(:)
    !        integer(kind=kint),pointer :: node2_ptr(:)
    !        integer(kind=kint),pointer :: dof_ptr(:)
    !        integer(kind=kint) :: n, old_size, new_size
    !
    !        rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', 'RIGID ', 1, 'P', type )
    !        if( rcode < 0 ) call fstr_ctrl_err_stop
    !
    !        n = fstr_ctrl_get_data_line_n( ctrl )
    !        if( n == 0 ) return
    !        old_size = P%MPC_RD%nmpc
    !        new_size = old_size + n
    !        P%MPC_RD%nmpc = new_size
    !
    !        call hecmw_expand_integer_array ( P%MPC_RD%node1,  old_size, new_size )
    !        call hecmw_expand_integer_array ( P%MPC_RD%node2,  old_size, new_size )
    !        call hecmw_expand_integer_array ( P%MPC_RD%dof,    old_size, new_size )
    !
    !        node1_ptr => P%MPC_RD%node1(old_size+1:)
    !        node2_ptr => P%MPC_RD%node2(old_size+1:)
    !        dof_ptr   => P%MPC_RD%dof(old_size+1:)
    !
    !        rcode = fstr_ctrl_get_MPC( ctrl, type, node1_ptr, node2_ptr, dof_ptr )
    !        if( rcode /= 0 ) call fstr_ctrl_err_stop
    !
    !        if( node_global_to_local( P%MESH, node1_ptr, n ) /= n ) then
    !                call fstr_setup_util_err_stop( '### Error : not exist node (!MPC)' )
    !        endif
    !        if( node_global_to_local( P%MESH, node2_ptr, n ) /= n ) then
    !                call fstr_setup_util_err_stop( '### Error : not exist node (!MPC)' )
    !        endif

    !   penalty => svRarray(11)
    rcode = fstr_ctrl_get_MPC( ctrl, p%MAT%Rarray(11))
    if( rcode /= 0) call fstr_ctrl_err_stop
  end subroutine fstr_setup_MPC

  !-----------------------------------------------------------------------------!
  !> Read in !CONTACT                                                           !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_CONTACTALGO( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode


    rcode = fstr_ctrl_get_CONTACTALGO( ctrl, P%PARAM%contact_algo, P%PARAM%augiter, &
      P%PARAM%fric_cone_follow, P%PARAM%eps_fric_band )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_CONTACTALGO

  !-----------------------------------------------------------------------------!
  !> Read in !OUTPUT_SSTYPE                                                         !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_OUTPUT_SSTYPE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, nid
    character(len=HECMW_NAME_LEN) :: data_fmt

    data_fmt = 'SOLUTION,MATERIAL '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', data_fmt, 0, 'P', nid )
    OPSSTYPE = nid
    if( rcode /= 0 ) call fstr_ctrl_err_stop

  end subroutine fstr_setup_OUTPUT_SSTYPE

  !-----------------------------------------------------------------------------!
  !> Read in !AUTOINC_PARAM                                                     !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_AUTOINC_PARAM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    fstrPARAM => P%PARAM

    if( fstr_get_AUTOINC( ctrl, fstrPARAM%ainc(slot%first) ) /=0 ) then
      write(*,*) '### Error: Fail in read in AUTOINC_PARAM definition : ' , slot%first
      write(ILOG,*) '### Error: Fail in read in AUTOINC_PARAM definition : ', slot%first
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_AUTOINC_PARAM

  !-----------------------------------------------------------------------------!
  !> Read in !CONVERG_PARAM                                                     !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_CONVERG_PARAM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    fstrPARAM => P%PARAM

    if( fstr_ctrl_get_CONVERGPARAM( ctrl, fstrPARAM%cnvparam(slot%first) ) /=0 ) then
      write(*,*) '### Error: Fail in read in CONVERG_PARAM definition : ' , slot%first
      write(ILOG,*) '### Error: Fail in read in CONVERG_PARAM definition : ', slot%first
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_CONVERG_PARAM

  !-----------------------------------------------------------------------------!
  !> Read in !TIME_POINTS                                                       !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_TIME_POINTS( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    fstrPARAM => P%PARAM

    if( fstr_ctrl_get_TIMEPOINTS( ctrl, fstrPARAM%timepoints(slot%first) )/=0 ) then
      write(*,*) '### Error: Fail in read in TIME_POINTS definition : ' , slot%first
      write(ILOG,*) '### Error: Fail in read in TIME_POINTS definition : ', slot%first
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_TIME_POINTS

  !-----------------------------------------------------------------------------!
  !> Read in !CONTACT_PARAM                                                     !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_CONTACT_PARAM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    fstrPARAM => P%PARAM

    if( fstr_ctrl_get_CONTACTPARAM( ctrl, fstrPARAM%contactparam(slot%first) ) /=0 ) then
      write(*,*) '### Error: Fail in read in CONTACT_PARAM definition : ' , slot%first
      write(ILOG,*) '### Error: Fail in read in CONTACT_PARAM definition : ', slot%first
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_CONTACT_PARAM

  !-----------------------------------------------------------------------------!
  !> Read in !CONTACT_INTERFERENCE                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_CONTACT_INTERFERENCE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    integer(kind=kint) :: c_contact_if, n, i
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID

    n = fstr_ctrl_get_data_line_n( ctrl )
    c_contact_if = slot%first - 1
    if( fstr_ctrl_get_CONTACT_IF( ctrl, n, fstrPARAM%contact_if(slot%first:slot%last) ) /= 0 ) then
      write(*,*) '### Error: Fail in read in CONTACT_INTERFERENCE definition : ' , c_contact_if
      write(ILOG,*) '### Error: Fail in read in CONTACT_INTERFERENCE definition : ', c_contact_if
      stop HECMW_EXIT_INPUT
    endif
    do i=1, n
      if( check_apply_Contact_IF(fstrPARAM%contact_if(c_contact_if+i), fstrSOLID%contacts) /= 0) then
        write(*,*) '### Error:(INTERFERENCE) Inconsistence of contact_pair in CONTACTS: ' , i+c_contact_if
        write(ILOG,*) '### Error:(INTERFERENCE)  Inconsistence of contact_pair in CONTACTS: ', i+c_contact_if
        stop HECMW_EXIT_MODEL
      end if
    end do
  end subroutine fstr_setup_CONTACT_INTERFERENCE

end module m_fstr_param_setup
