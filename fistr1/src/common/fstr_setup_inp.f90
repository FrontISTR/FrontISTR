!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to read in the analysis data of an inp file
!! and do necessary preparation for following calculation, like m_fstr_setup for a cnt file.
!! hecmw1 reads the mesh, the materials and the sections of the same file.
!! The cards of the inp format are read in its syntax, and the cards of FrontISTR's own settings
!! are the cnt cards with '*' as the header character, read by the readers of the cnt cards.
module m_fstr_setup_inp
  use m_fstr
  use hecmw_setup_util, only &
    : hecmw_expand_integer_array &
    , hecmw_expand_real_array
  use fstr_setup_util
  use m_fstr_solid_setup, only &
    : fstr_setup_SECTION &
    , fstr_setup_WRITE_flags &
    , fstr_setup_WRITE_store &
    , fstr_setup_step_links &
    , fstr_solid_append_boundary &
    , fstr_solid_append_cload
  use m_fstr_param_setup, only &
    : fstr_setup_SOLVER &
    , fstr_setup_NONLINEAR_SOLVER &
    , fstr_setup_AUTOINC_PARAM &
    , fstr_setup_CONVERG_PARAM
  use m_fstr_setup_finish
  use fstr_ctrl_util_f
  implicit none
  private
  public :: fstr_setup_inp

contains

  !=============================================================================!
  !> Read in and initialize the analysis data of an inp file                    !
  !=============================================================================!
  subroutine fstr_setup_inp( inp_filename, hecMESH, hecMAT, fstrPARAM,  &
      fstrSOLID, fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL, fstrFREQ )
    character(len=HECMW_FILENAME_LEN) :: inp_filename
    type(hecmwST_local_mesh),target :: hecMESH
    type(hecmwST_matrix),target :: hecMAT
    type(fstr_param),target   :: fstrPARAM
    type(fstr_solid),target   :: fstrSOLID
    type(fstr_eigen),target     :: fstrEIG
    type(fstr_heat),target    :: fstrHEAT
    type(fstr_dynamic),target :: fstrDYNAMIC
    type(fstr_couple),target  :: fstrCPL
    type(fstr_freqanalysis), target :: fstrFREQ

    integer(kind=kint) :: ctrl, rcode
    type(fstr_param_pack) :: P
    type(fstr_card_table) :: T
    integer(kind=kint) :: total(FSTR_CARD_MAX)
    type(fstr_setup_counts) :: counts

    ! packaging
    P%MESH   => hecMESH
    P%MAT    => hecMAT
    P%PARAM  => fstrPARAM
    P%SOLID  => fstrSOLID
    P%EIGEN  => fstrEIG
    P%HEAT   => fstrHEAT
    P%DYN    => fstrDYNAMIC
    P%CPL    => fstrCPL
    P%FREQ   => fstrFREQ

    fstrPARAM%contact_algo = kcaALagrange

    ! Defaults of an inp file: the linear solver is DIRECTmf unless *SOLVER is given,
    ! with the solver type of a direct method as fstr_setup_SOLVER sets it, and the result file is written.
    hecMAT%Iarray(2)  = 106
    hecMAT%Iarray(99) = hecMAT%Iarray(2) - 99
    fstrPARAM%fg_result = kON

    call fstr_setup_inp_cards( T )

    ctrl = fstr_ctrl_open_inp( inp_filename )
    if( ctrl < 0 ) then
      write(*,*) '### Error: Cannot read the inp file : ', trim(inp_filename)
      write(ILOG,*) '### Error: Cannot read the inp file : ', trim(inp_filename)
      call fstr_ctrl_err_stop
    end if

    P%STATE%step = 0
    call fstr_setup_inp_read_cards( ctrl, T, 1, total, P )

    ! ----- the other counts stay 0, since the other cards are not stored in the arrays that fstr_setup_prepare allocates
    counts%steps          = fstr_card_total( T, total, 'STEP' )
    counts%autoinc_params = fstr_card_total( T, total, 'AUTOINC_PARAM' )
    counts%converg_params = fstr_card_total( T, total, 'CONVERG_PARAM' )
    if( counts%steps > 0 .and. fstr_card_total( T, total, 'STATIC' ) /= counts%steps ) &
      call fstr_setup_util_err_stop( '### Error: *STEP : each step needs one *STATIC' )
    call fstr_setup_inp_groups( P, counts%steps )
    call fstr_setup_prepare( P, counts )
    call fstr_setup_inp_formulation( P )

    ! -----
    rcode = fstr_ctrl_rewind( ctrl )
    P%STATE%step = 0
    call fstr_setup_inp_read_cards( ctrl, T, 2, total, P )
    call fstr_setup_inp_steps( P )

    call fstr_setup_finish( P, counts )
    rcode = fstr_ctrl_close( ctrl )

  end subroutine fstr_setup_inp

  !> Adds the cards of an inp file to the table T, by the names without the '*' of the header.
  !> The cards of the inp format are read in the first reading, so that fstr_setup_prepare knows the analysis type,
  !> the steps and all the boundary conditions. The second reading stores the steps that fstr_setup_prepare allocates.
  subroutine fstr_setup_inp_cards( T )
    implicit none
    type(fstr_card_table) :: T

    ! The cards of the inp format. Without *STEP, the conditions go into the default step, as in a cnt file
    ! without !STEP. *STEP also refers by name to *AUTOINC_PARAM and *CONVERG_PARAM read before it.
    call fstr_card_add( T, 'STEP',     first=fstr_setup_inp_STEP, second=fstr_setup_inp_STEP_store )
    call fstr_card_add( T, 'END STEP' )
    call fstr_card_add( T, 'STATIC',   first=fstr_setup_inp_STATIC, second=fstr_setup_inp_STATIC_store )
    call fstr_card_add( T, 'BOUNDARY', first=fstr_setup_inp_BOUNDARY )
    call fstr_card_add( T, 'CLOAD',    first=fstr_setup_inp_CLOAD )

    ! FrontISTR's own cards, read as the cnt cards
    call fstr_card_add( T, 'SOLVER',   first=fstr_setup_SOLVER )
    call fstr_card_add( T, 'NONLINEAR_SOLVER', first=fstr_setup_NONLINEAR_SOLVER )
    call fstr_card_add( T, 'AUTOINC_PARAM',    second=fstr_setup_AUTOINC_PARAM )
    call fstr_card_add( T, 'CONVERG_PARAM',    second=fstr_setup_CONVERG_PARAM )
    call fstr_card_add( T, 'SECTION',  second=fstr_setup_SECTION )
    call fstr_card_add( T, 'WRITE',    first=fstr_setup_inp_WRITE, second=fstr_setup_WRITE_store )

    ! Read by hecmw1
    call fstr_card_add( T, 'HEADING' )
    call fstr_card_add( T, 'NODE' )
    call fstr_card_add( T, 'ELEMENT' )
    call fstr_card_add( T, 'NSET' )
    call fstr_card_add( T, 'ELSET' )
    call fstr_card_add( T, 'MATERIAL' )
    call fstr_card_add( T, 'ELASTIC' )
    call fstr_card_add( T, 'DENSITY' )
    call fstr_card_add( T, 'CONDUCTIVITY' )
    call fstr_card_add( T, 'SPECIFIC HEAT' )
    call fstr_card_add( T, 'SOLID SECTION' )
    call fstr_card_add( T, 'SHELL SECTION' )
    call fstr_card_add( T, 'EQUATION' )
    call fstr_card_add( T, 'AMPLITUDE' )

    ! The controls of the iterations and the times of the output requests, skipped with a warning.
    ! They do not change the results at the end of a step.
    call fstr_card_add( T, 'CONTROLS',       first=fstr_setup_inp_skipped )
    call fstr_card_add( T, 'TIME POINTS',    first=fstr_setup_inp_skipped )

    ! Output requests, skipped with a warning. FrontISTR writes its own result file.
    call fstr_card_add( T, 'NODE PRINT',     first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'EL PRINT',       first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'NODE FILE',      first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'EL FILE',        first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'CONTACT PRINT',  first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'CONTACT FILE',   first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'SECTION PRINT',  first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'OUTPUT',         first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'NODE OUTPUT',    first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'ELEMENT OUTPUT', first=fstr_setup_inp_output_request )
    call fstr_card_add( T, 'CONTACT OUTPUT', first=fstr_setup_inp_output_request )
  end subroutine fstr_setup_inp_cards

  !> Reads the cards of the inp file in order, and calls the reader of each card in the table T for the reading ipass
  !> (1 or 2), as fstr_setup_read_cards does for a cnt file. The first reading stops at a card not in T, since the cards
  !> of the inp format that are not read yet would change the analysis. A parameter of an analysis card that its reader
  !> does not read is skipped with a warning.
  subroutine fstr_setup_inp_read_cards( ctrl, T, ipass, total, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_card_table), intent(in) :: T
    integer(kind=kint), intent(in) :: ipass
    integer(kind=kint), intent(inout) :: total(:)
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, k, comma, next_comma, i, c
    integer(kind=kint) :: pos(size(total))
    character(len=HECMW_NAME_LEN) :: header_name
    character(len=256) :: header, param
    logical :: supported

    pos(:) = 0
    do
      rcode = fstr_ctrl_get_c_h_name_inp( ctrl, header_name, HECMW_NAME_LEN )
      if( ipass == 1 ) then
        select case( trim(header_name) )
        case( '*STEP', '*END STEP', '*STATIC', '*BOUNDARY', '*CLOAD' )
          rcode = fstr_ctrl_get_line( ctrl, fstr_ctrl_get_c_h_pos(ctrl), header, len(header) )
          if( rcode /= 0 ) call fstr_ctrl_err_stop
          comma = index( header, ',' )
          do while( comma > 0 )
            header = adjustl( header(comma+1:) )
            next_comma = index( header, ',' )
            param = header
            if( next_comma > 0 ) param = header(:next_comma-1)
            if( len_trim(param) > 0 ) then
              i = index( param, '=' )
              if( i > 0 ) param = param(:i-1)
              do i = 1, len_trim(param)
                c = iachar( param(i:i) )
                if( c >= iachar('a') .and. c <= iachar('z') ) param(i:i) = achar( c - iachar('a') + iachar('A') )
              enddo
              supported = .false.
              select case( trim(header_name) )
              case( '*STEP' )
                supported = param == 'NLGEOM' .or. param == 'INC' .or. param == 'AUTOINCPARAM' &
                  .or. param == 'CONVERGPARAM'
              case( '*STATIC' )
                supported = param == 'DIRECT'
              case( '*BOUNDARY' )
                supported = param == 'AMPLITUDE' .or. param == 'OP' .or. param == 'FIXED'
              case( '*CLOAD' )
                supported = param == 'AMPLITUDE' .or. param == 'OP'
              end select
              ! FrontISTR cannot shift the time of an amplitude
              if( param == 'TIME DELAY' .and. ( header_name == '*BOUNDARY' .or. header_name == '*CLOAD' ) ) &
                call fstr_setup_inp_err_stop( ctrl, 'TIME DELAY is not supported yet' )
              if( .not. supported .and. P%MESH%my_rank == 0 ) write(*,'(a,a,a,i0,a,a,a)') ' Warning: ', &
                trim(header_name), ' (line:', fstr_ctrl_get_c_h_line_no( ctrl ), ') : unsupported parameter ', &
                trim(param), ' is skipped'
            endif
            comma = next_comma
          enddo
        end select
      endif
      k = fstr_card_find( T, header_name(2:), P%STATE%version )
      if( k > 0 ) then
        call fstr_card_read( ctrl, T, k, ipass, header_name, pos, total, P )
      else if( ipass == 1 ) then
        call fstr_setup_inp_err_stop( ctrl, 'not supported in an inp file yet' )
      endif
      if( fstr_ctrl_seek_next_header( ctrl ) == 0 ) exit
    end do
    if( ipass == 1 ) total(:) = pos(:)
  end subroutine fstr_setup_inp_read_cards

  !> Stops with msg about the card at the current header of ctrl and its line in the inp file
  subroutine fstr_setup_inp_err_stop( ctrl, msg )
    implicit none
    integer(kind=kint) :: ctrl
    character(len=*), intent(in) :: msg

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: header_name
    character(len=HECMW_MSG_LEN) :: s

    rcode = fstr_ctrl_get_c_h_name_inp( ctrl, header_name, HECMW_NAME_LEN )
    write( s, '(a,a,a,i0,a,a)' ) '### Error: ', trim(header_name), ' (line:', fstr_ctrl_get_c_h_line_no( ctrl ), &
      ') : ', msg
    call fstr_setup_util_err_stop( trim(s) )
  end subroutine fstr_setup_inp_err_stop

  !> Stops at AMPLITUDE= of a condition in an analysis of more than one step. See fstr_setup_inp_STEP.
  subroutine fstr_setup_inp_amp_err_stop( ctrl )
    implicit none
    integer(kind=kint) :: ctrl

    call fstr_setup_inp_err_stop( ctrl, 'AMPLITUDE= of a condition with more than one step is not supported yet' )
  end subroutine fstr_setup_inp_amp_err_stop

  !> OP=NEW of a *BOUNDARY or a *CLOAD in a later step is used only on the first card of the name in the step, as in
  !> the inp format. The card is not the first if the groups grpid(1:n) of the conditions read so far have the step.
  subroutine fstr_setup_inp_first_op( ctrl, P, header_name, grpid, n, step, op )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    character(len=*), intent(in) :: header_name
    integer(kind=kint), pointer :: grpid(:)
    integer(kind=kint), intent(in) :: n, step
    integer(kind=kint), intent(inout) :: op

    if( n == 0 ) return
    if( .not. any( grpid(1:n) == step ) ) return
    op = 1
    if( P%MESH%my_rank == 0 ) write(*,'(a,a,a,i0,a,a,a)') ' Warning: ', header_name, ' (line:', &
      fstr_ctrl_get_c_h_line_no( ctrl ), ') : OP=NEW is not used, since OP of the first ', header_name, ' of a step is used'
  end subroutine fstr_setup_inp_first_op

  !-----------------------------------------------------------------------------!
  !> Read in *STEP                                                              !
  !-----------------------------------------------------------------------------!
  !> P%STATE%step is set to this step. NLGEOM given in the first step makes the analysis geometrically nonlinear,
  !> as NLSTATIC of !SOLUTION, and stays in the later steps as in the inp format. Since a FrontISTR analysis is
  !> geometrically nonlinear in all the steps or in none, NLGEOM given first in a later step stops.
  !> AMPLITUDE= of a condition stops with more than one step, since FrontISTR takes the amplitude again in each step
  !> a condition is active, while the inp format keeps the value at the end of the step in the later steps.
  subroutine fstr_setup_inp_STEP( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, is_nlgeom, nlgeom

    P%STATE%step = P%STATE%slot%first
    if( P%STATE%step > 1 ) then
      if( P%SOLID%BOUNDARY_ngrp_tot > 0 ) then
        if( any( P%SOLID%BOUNDARY_ngrp_amp(1:P%SOLID%BOUNDARY_ngrp_tot) > 0 ) ) call fstr_setup_inp_amp_err_stop( ctrl )
      endif
      if( P%SOLID%CLOAD_ngrp_tot > 0 ) then
        if( any( P%SOLID%CLOAD_ngrp_amp(1:P%SOLID%CLOAD_ngrp_tot) > 0 ) ) call fstr_setup_inp_amp_err_stop( ctrl )
      endif
    endif

    is_nlgeom = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'NLGEOM ', '# ', 0, 'E', is_nlgeom )
    if( is_nlgeom == 0 ) return
    nlgeom = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'NLGEOM ', 'YES,NO ', 0, 'P', nlgeom )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( nlgeom == 1 ) then
      if( P%STATE%step > 1 .and. .not. P%PARAM%nlgeom ) &
        call fstr_setup_inp_err_stop( ctrl, 'NLGEOM given first in a later step is not supported yet' )
      P%PARAM%nlgeom = .true.
    else if( P%PARAM%nlgeom .and. P%MESH%my_rank == 0 ) then
      write(*,'(a,i0,a)') ' Warning: *STEP (line:', fstr_ctrl_get_c_h_line_no( ctrl ), &
        ') : NLGEOM=NO is not used, since NLGEOM of an earlier step stays'
    endif
  end subroutine fstr_setup_inp_STEP

  !> Read in *STEP after fstr_setup_prepare and store the step: INC is the upper bound of the automatic increments
  !> (100 by default), and FrontISTR's own AUTOINCPARAM and CONVERGPARAM refer to *AUTOINC_PARAM and *CONVERG_PARAM
  !> as those of !STEP.
  subroutine fstr_setup_inp_STEP_store( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, inc
    character(len=HECMW_NAME_LEN) :: ainc_name, cnv_name

    P%STATE%step = P%STATE%slot%first
    inc = 100
    rcode = fstr_ctrl_get_param_ex( ctrl, 'INC ', '# ', 0, 'I', inc )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( inc < 1 ) call fstr_setup_inp_err_stop( ctrl, 'INC must be positive' )
    ainc_name = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AUTOINCPARAM ', '# ', 0, 'S', ainc_name )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    cnv_name = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'CONVERGPARAM ', '# ', 0, 'S', cnv_name )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    P%SOLID%step_ctrl(P%STATE%step)%num_substep = inc
    call fstr_setup_step_links( P%PARAM, P%SOLID, P%DYN, P%STATE%step, ' ', ainc_name, cnv_name )
  end subroutine fstr_setup_inp_STEP_store

  !-----------------------------------------------------------------------------!
  !> Read in *STATIC                                                            !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_inp_STATIC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    if( P%STATE%step > 0 .and. P%STATE%slot%first /= P%STATE%step ) &
      call fstr_setup_inp_err_stop( ctrl, 'each step needs one *STATIC' )
    P%PARAM%solution_type = kstSTATIC
  end subroutine fstr_setup_inp_STATIC

  !> Read in *STATIC after fstr_setup_prepare and store its time increments in the step. The data line is the initial
  !> increment, the time period of the step, and the minimum and the maximum increments, used in a geometrically
  !> nonlinear analysis: DIRECT fixes the increment to the initial one, and the increments are automatic without it.
  !> A step of a linear analysis is solved at once over the time period 1, whatever the data line is.
  subroutine fstr_setup_inp_STATIC_store( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, is_direct, k
    real(kind=kreal) :: dt, period, dtmin, dtmax

    k = P%STATE%step
    if( k == 0 .or. .not. P%PARAM%nlgeom ) return
    is_direct = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DIRECT ', '# ', 0, 'E', is_direct )
    dt = 1.0d0
    period = 1.0d0
    dtmin = -1.0d0
    dtmax = 1.0d+30
    if( fstr_ctrl_get_data_line_n( ctrl ) > 0 ) then
      rcode = fstr_ctrl_get_data_ex( ctrl, 1, 'rrrr ', dt, period, dtmin, dtmax )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif
    if( dt <= 0.0d0 .or. period <= 0.0d0 ) &
      call fstr_setup_inp_err_stop( ctrl, 'the time increment and the time period must be positive' )
    if( dtmin < 0.0d0 ) dtmin = min( dt, 1.0d-5*period )

    P%SOLID%step_ctrl(k)%initdt = dt
    P%SOLID%step_ctrl(k)%elapsetime = period
    if( is_direct == 1 ) then
      P%SOLID%step_ctrl(k)%inc_type = stepFixedInc
      P%SOLID%step_ctrl(k)%num_substep = max( int( (period + 0.999999999d0*dt)/dt ), 1 )
    else
      P%SOLID%step_ctrl(k)%inc_type = stepAutoInc
      P%SOLID%step_ctrl(k)%mindt = min( dtmin, dt )
      P%SOLID%step_ctrl(k)%maxdt = dtmax
    endif
  end subroutine fstr_setup_inp_STATIC_store

  !-----------------------------------------------------------------------------!
  !> Read in *BOUNDARY                                                          !
  !-----------------------------------------------------------------------------!
  !> A data line is a node or a node set, the first dof, the last dof (the first dof if omitted)
  !> and the value (0 if omitted). The value is the total displacement, which TYPE=TOTAL of !BOUNDARY takes
  !> from the displacement at the start of the step, and FIXED holds the displacement at the end of the
  !> previous step, as an incremental displacement of 0.
  !> Each dof is an entry of its own, whose group (GRPID) is the step it starts in until fstr_setup_inp_groups.
  !> The entries of the earlier steps stop in this step at the same nodes and dof, and all of them with OP=NEW
  !> of the first *BOUNDARY in the step (fstr_setup_inp_first_op).
  !> P%STATE%boundary_end keeps the step each entry stops in.
  subroutine fstr_setup_inp_BOUNDARY( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n, m, i, d, op, is_fixed, istot, step, old_size, e, f
    character(len=HECMW_NAME_LEN) :: amp, rotc_name, data_fmt
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:), names(:)
    integer(kind=kint), allocatable :: dof_ids(:), dof_ide(:), dofs(:)
    real(kind=kreal), allocatable :: vals(:), dof_vals(:)
    logical, allocatable :: mark(:)

    amp = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AMPLITUDE ', '# ', 0, 'S', amp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( amp /= ' ' .and. P%STATE%step > 1 ) call fstr_setup_inp_amp_err_stop( ctrl )
    op = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'OP ', 'MOD,NEW ', 0, 'P', op )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    is_fixed = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'FIXED ', '# ', 0, 'E', is_fixed )
    rotc_name = ' '

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_names(n), dof_ids(n), dof_ide(n), vals(n) )
    dof_ide = 0
    vals = 0.0d0
    if( n > 0 ) then
      write( data_fmt, '(a,i0,a)' ) 'S', HECMW_NAME_LEN, 'Iir '
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_names, dof_ids, dof_ide, vals )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif
    where( dof_ide == 0 ) dof_ide = dof_ids
    istot = 1
    if( is_fixed == 1 ) then
      istot = 0
      vals = 0.0d0
    endif

    m = sum( max( dof_ide - dof_ids + 1, 0 ) )
    allocate( names(m), dofs(m), dof_vals(m) )
    m = 0
    do i = 1, n
      do d = dof_ids(i), dof_ide(i)
        m = m + 1
        names(m) = grp_names(i)
        dofs(m) = d
        dof_vals(m) = vals(i)
      enddo
    enddo

    step = max( P%STATE%step, 1 )
    old_size = P%SOLID%BOUNDARY_ngrp_tot
    if( op == 2 .and. P%STATE%step > 1 ) &
      call fstr_setup_inp_first_op( ctrl, P, '*BOUNDARY', P%SOLID%BOUNDARY_ngrp_GRPID, old_size, step, op )
    if( op == 2 .and. P%STATE%step > 1 ) then
      do f = 1, old_size
        if( P%STATE%boundary_end(f) == 0 .and. P%SOLID%BOUNDARY_ngrp_GRPID(f) < step ) P%STATE%boundary_end(f) = step
      enddo
    endif
    call fstr_solid_append_boundary( P%MESH, P%SOLID, '*BOUNDARY', step, istot, rotc_name, &
      amp, m, names, dofs, dofs, dof_vals )
    call hecmw_expand_integer_array( P%STATE%boundary_end, old_size, P%SOLID%BOUNDARY_ngrp_tot )

    if( P%STATE%step <= 1 ) return
    allocate( mark(P%MESH%n_node) )
    mark = .false.
    do e = old_size + 1, P%SOLID%BOUNDARY_ngrp_tot
      do f = 1, old_size
        if( P%STATE%boundary_end(f) /= 0 .or. P%SOLID%BOUNDARY_ngrp_GRPID(f) >= step ) cycle
        if( P%SOLID%BOUNDARY_ngrp_type(f) /= P%SOLID%BOUNDARY_ngrp_type(e) ) cycle
        if( fstr_setup_inp_same_nodes( P%MESH, P%SOLID%BOUNDARY_ngrp_ID(f), P%SOLID%BOUNDARY_ngrp_ID(e), mark ) ) &
          P%STATE%boundary_end(f) = step
      enddo
    enddo
  end subroutine fstr_setup_inp_BOUNDARY

  !-----------------------------------------------------------------------------!
  !> Read in *CLOAD                                                             !
  !-----------------------------------------------------------------------------!
  !> A data line is a node or a node set, the dof and the value. The group (GRPID) of a load is the step it starts in,
  !> and a load is active to the end, so that FrontISTR keeps its value in the later steps. A load replacing
  !> those of the earlier steps at the same nodes and dof, and OP=NEW of the first *CLOAD in the step, add the load
  !> that cancels them by fstr_setup_inp_cancel_cload. The loads given at the same nodes and dof in a step are summed.
  subroutine fstr_setup_inp_CLOAD( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n, op, step, old_size, e, f
    character(len=HECMW_NAME_LEN) :: amp, rotc_name, data_fmt
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:)
    integer(kind=kint), allocatable :: dofs(:)
    real(kind=kreal), allocatable :: vals(:)
    logical, allocatable :: mark(:)

    amp = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AMPLITUDE ', '# ', 0, 'S', amp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( amp /= ' ' .and. P%STATE%step > 1 ) call fstr_setup_inp_amp_err_stop( ctrl )
    op = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'OP ', 'MOD,NEW ', 0, 'P', op )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    rotc_name = ' '

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_names(n), dofs(n), vals(n) )
    if( n > 0 ) then
      write( data_fmt, '(a,i0,a)' ) 'S', HECMW_NAME_LEN, 'IR '
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_names, dofs, vals )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif

    step = max( P%STATE%step, 1 )
    old_size = P%SOLID%CLOAD_ngrp_tot
    if( P%STATE%step > 1 ) then
      allocate( mark(P%MESH%n_node) )
      mark = .false.
    endif
    if( op == 2 .and. P%STATE%step > 1 ) &
      call fstr_setup_inp_first_op( ctrl, P, '*CLOAD', P%SOLID%CLOAD_ngrp_GRPID, old_size, step, op )
    if( op == 2 .and. P%STATE%step > 1 ) then
      do f = 1, old_size
        if( P%SOLID%CLOAD_ngrp_GRPID(f) >= step ) cycle
        call fstr_setup_inp_cancel_cload( P, f, step, P%SOLID%CLOAD_ngrp_tot + 1, mark )
      enddo
    endif
    old_size = P%SOLID%CLOAD_ngrp_tot
    call fstr_solid_append_cload( P%MESH, P%SOLID, '*CLOAD', step, rotc_name, &
      amp, n, grp_names, dofs, vals )

    if( P%STATE%step <= 1 ) return
    do e = old_size + 1, old_size + n
      call fstr_setup_inp_cancel_cload( P, e, step, e, mark )
    enddo
  end subroutine fstr_setup_inp_CLOAD

  !> Adds the load that cancels, from step on, the loads of the earlier steps at the nodes and the dof of the load
  !> src, as a load of the group step whose value is minus their sum. Nothing is added if a load of step before
  !> the load upto has the same nodes and dof, since the loads of the earlier steps there are cancelled already.
  !> mark is all .false. on entry and on return.
  subroutine fstr_setup_inp_cancel_cload( P, src, step, upto, mark )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: src, step, upto
    logical :: mark(:)

    integer(kind=kint) :: f, n
    real(kind=kreal) :: val, scale

    val = 0.0d0
    scale = 0.0d0
    do f = 1, P%SOLID%CLOAD_ngrp_tot
      if( P%SOLID%CLOAD_ngrp_DOF(f) /= P%SOLID%CLOAD_ngrp_DOF(src) ) cycle
      if( P%SOLID%CLOAD_ngrp_GRPID(f) > step ) cycle
      if( P%SOLID%CLOAD_ngrp_GRPID(f) == step .and. f >= upto ) cycle
      if( .not. fstr_setup_inp_same_nodes( P%MESH, P%SOLID%CLOAD_ngrp_ID(f), P%SOLID%CLOAD_ngrp_ID(src), mark ) ) cycle
      if( P%SOLID%CLOAD_ngrp_GRPID(f) == step ) return
      val = val + P%SOLID%CLOAD_ngrp_val(f)
      scale = scale + abs( P%SOLID%CLOAD_ngrp_val(f) )
    enddo
    if( fstr_setup_inp_cancelled( val, scale ) ) return

    n = P%SOLID%CLOAD_ngrp_tot + 1
    call hecmw_expand_integer_array( P%SOLID%CLOAD_ngrp_GRPID, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%CLOAD_ngrp_ID, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%CLOAD_ngrp_DOF, n-1, n )
    call hecmw_expand_real_array( P%SOLID%CLOAD_ngrp_val, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%CLOAD_ngrp_amp, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%CLOAD_ngrp_rotID, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%CLOAD_ngrp_centerID, n-1, n )
    P%SOLID%CLOAD_ngrp_tot = n
    P%SOLID%CLOAD_ngrp_GRPID(n) = step
    P%SOLID%CLOAD_ngrp_ID(n) = P%SOLID%CLOAD_ngrp_ID(src)
    P%SOLID%CLOAD_ngrp_DOF(n) = P%SOLID%CLOAD_ngrp_DOF(src)
    P%SOLID%CLOAD_ngrp_val(n) = -val
    P%SOLID%CLOAD_ngrp_amp(n) = P%SOLID%CLOAD_ngrp_amp(src)
    P%SOLID%CLOAD_ngrp_rotID(n) = P%SOLID%CLOAD_ngrp_rotID(src)
    P%SOLID%CLOAD_ngrp_centerID(n) = P%SOLID%CLOAD_ngrp_centerID(src)
  end subroutine fstr_setup_inp_cancel_cload

  !-----------------------------------------------------------------------------!
  !> Read in *WRITE                                                             !
  !-----------------------------------------------------------------------------!
  !> The same as !WRITE of a cnt file without VISUAL, since the visualization is set by !VISUAL of a cnt file.
  subroutine fstr_setup_inp_WRITE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    call fstr_setup_WRITE_flags( ctrl, P )
    if( P%PARAM%fg_visual == kON ) call fstr_setup_inp_err_stop( ctrl, 'VISUAL is not supported in an inp file yet' )
  end subroutine fstr_setup_inp_WRITE

  !> An output request of the inp format, skipped with a warning
  subroutine fstr_setup_inp_output_request( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: header_name

    if( P%MESH%my_rank /= 0 ) return
    rcode = fstr_ctrl_get_c_h_name_inp( ctrl, header_name, HECMW_NAME_LEN )
    write(*,'(a,a,a,i0,a)') ' Warning: ', trim(header_name), ' (line:', fstr_ctrl_get_c_h_line_no( ctrl ), &
      ') is skipped. The results are written to the result file of FrontISTR'
  end subroutine fstr_setup_inp_output_request

  !> A card of the inp format skipped with a warning
  subroutine fstr_setup_inp_skipped( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: header_name

    if( P%MESH%my_rank /= 0 ) return
    rcode = fstr_ctrl_get_c_h_name_inp( ctrl, header_name, HECMW_NAME_LEN )
    write(*,'(a,a,a,i0,a)') ' Warning: ', trim(header_name), ' (line:', fstr_ctrl_get_c_h_line_no( ctrl ), &
      ') is skipped'
  end subroutine fstr_setup_inp_skipped

  !> After the first reading of an inp file of nstep steps, stops at a load replacing on some nodes the loads of
  !> the earlier steps given to other nodes, and sets the group (GRPID) of each !BOUNDARY entry to the pair of
  !> the step it starts in and the step it stops in, which fstr_setup_inp_steps takes back.
  !> A load of a step stops at the nodes of an earlier load that is not cancelled in the step, which FrontISTR
  !> cannot do, since a load has a value for all the nodes of its group.
  subroutine fstr_setup_inp_groups( P, nstep )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: nstep

    integer(kind=kint) :: e, f, g, k, s, last
    real(kind=kreal) :: val, scale
    logical, allocatable :: mark(:)
    character(len=HECMW_MSG_LEN) :: msg

    if( nstep == 0 ) then
      if( associated(P%STATE%boundary_end) ) deallocate( P%STATE%boundary_end )
      return
    endif
    allocate( mark(P%MESH%n_node) )
    mark = .false.
    do e = 1, P%SOLID%CLOAD_ngrp_tot
      k = P%SOLID%CLOAD_ngrp_GRPID(e)
      if( k == 1 ) cycle
      do f = 1, P%SOLID%CLOAD_ngrp_tot
        if( P%SOLID%CLOAD_ngrp_GRPID(f) >= k .or. P%SOLID%CLOAD_ngrp_DOF(f) /= P%SOLID%CLOAD_ngrp_DOF(e) ) cycle
        if( fstr_setup_inp_same_nodes( P%MESH, P%SOLID%CLOAD_ngrp_ID(f), P%SOLID%CLOAD_ngrp_ID(e), mark ) ) cycle
        if( .not. fstr_setup_inp_overlap( P%MESH, P%SOLID%CLOAD_ngrp_ID(f), P%SOLID%CLOAD_ngrp_ID(e), mark ) ) cycle
        ! the loads at the nodes of f, unless a load of step k cancels them
        val = 0.0d0
        scale = 0.0d0
        do g = 1, P%SOLID%CLOAD_ngrp_tot
          if( P%SOLID%CLOAD_ngrp_GRPID(g) > k .or. P%SOLID%CLOAD_ngrp_DOF(g) /= P%SOLID%CLOAD_ngrp_DOF(f) ) cycle
          if( .not. fstr_setup_inp_same_nodes( P%MESH, P%SOLID%CLOAD_ngrp_ID(g), P%SOLID%CLOAD_ngrp_ID(f), mark ) ) cycle
          if( P%SOLID%CLOAD_ngrp_GRPID(g) == k ) then
            val = 0.0d0
            exit
          endif
          val = val + P%SOLID%CLOAD_ngrp_val(g)
          scale = scale + abs( P%SOLID%CLOAD_ngrp_val(g) )
        enddo
        if( fstr_setup_inp_cancelled( val, scale ) ) cycle
        write( msg, '(a,i0,a)' ) '### Error: *CLOAD : a load of step ', k, ' is given to some nodes of a load of '// &
          'an earlier step given to another node set, which is not supported yet'
        call fstr_setup_util_err_stop( trim(msg) )
      enddo
    enddo

    do e = 1, P%SOLID%BOUNDARY_ngrp_tot
      s = P%SOLID%BOUNDARY_ngrp_GRPID(e)
      last = P%STATE%boundary_end(e)
      if( last == 0 ) last = nstep + 1
      P%SOLID%BOUNDARY_ngrp_GRPID(e) = fstr_setup_inp_group_id( s, last, nstep )
    enddo
    if( associated(P%STATE%boundary_end) ) deallocate( P%STATE%boundary_end )
  end subroutine fstr_setup_inp_groups

  !> Whether the sum val of loads, whose absolute values sum to scale, is zero up to the round-off of the loads
  !> cancelling each other
  logical function fstr_setup_inp_cancelled( val, scale )
    implicit none
    real(kind=kreal), intent(in) :: val, scale

    fstr_setup_inp_cancelled = abs( val ) <= 1.0d-12*scale
  end function fstr_setup_inp_cancelled

  !> The group (GRPID) of a !BOUNDARY entry active from the step first to the step before last, of nstep steps
  integer(kind=kint) function fstr_setup_inp_group_id( first, last, nstep )
    implicit none
    integer(kind=kint), intent(in) :: first, last, nstep

    fstr_setup_inp_group_id = (last - 1)*nstep + first
  end function fstr_setup_inp_group_id

  !> Makes each step of the inp file take the conditions active in it: the !BOUNDARY entries by the steps in their
  !> groups (fstr_setup_inp_groups), and the loads of the step and the earlier steps.
  subroutine fstr_setup_inp_steps( P )
    implicit none
    type(fstr_param_pack) :: P

    integer(kind=kint) :: nstep, k, e, g, n, i
    integer(kind=kint), allocatable :: groups(:)

    if( .not. associated(P%SOLID%step_ctrl) ) return
    nstep = size(P%SOLID%step_ctrl)
    allocate( groups(max(P%SOLID%BOUNDARY_ngrp_tot, 1)) )
    do k = 1, nstep
      n = 0
      do e = 1, P%SOLID%BOUNDARY_ngrp_tot
        g = P%SOLID%BOUNDARY_ngrp_GRPID(e)
        if( any( groups(1:n) == g ) ) cycle
        if( mod(g-1, nstep) + 1 > k .or. (g-1)/nstep + 1 <= k ) cycle
        n = n + 1
        groups(n) = g
      enddo
      if( n > 0 ) then
        allocate( P%SOLID%step_ctrl(k)%Boundary(n) )
        P%SOLID%step_ctrl(k)%Boundary(:) = groups(1:n)
      endif
      if( P%SOLID%CLOAD_ngrp_tot + P%SOLID%DLOAD_ngrp_tot + P%SOLID%TEMP_ngrp_tot + P%SOLID%SPRING_ngrp_tot > 0 ) then
        allocate( P%SOLID%step_ctrl(k)%Load(k) )
        do i = 1, k
          P%SOLID%step_ctrl(k)%Load(i) = i
        enddo
      endif
    enddo
  end subroutine fstr_setup_inp_steps

  !> Whether the node groups ig1 and ig2 have the same nodes. mark is all .false. on entry and on return.
  logical function fstr_setup_inp_same_nodes( hecMESH, ig1, ig2, mark )
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint), intent(in) :: ig1, ig2
    logical :: mark(:)

    integer(kind=kint) :: is1, ie1, is2, ie2

    fstr_setup_inp_same_nodes = ig1 == ig2
    if( fstr_setup_inp_same_nodes ) return
    is1 = hecMESH%node_group%grp_index(ig1-1) + 1
    ie1 = hecMESH%node_group%grp_index(ig1)
    is2 = hecMESH%node_group%grp_index(ig2-1) + 1
    ie2 = hecMESH%node_group%grp_index(ig2)
    if( ie1 - is1 /= ie2 - is2 ) return
    mark( hecMESH%node_group%grp_item(is1:ie1) ) = .true.
    fstr_setup_inp_same_nodes = all( mark( hecMESH%node_group%grp_item(is2:ie2) ) )
    mark( hecMESH%node_group%grp_item(is1:ie1) ) = .false.
  end function fstr_setup_inp_same_nodes

  !> Whether the node groups ig1 and ig2 have a node in common. mark is all .false. on entry and on return.
  logical function fstr_setup_inp_overlap( hecMESH, ig1, ig2, mark )
    implicit none
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint), intent(in) :: ig1, ig2
    logical :: mark(:)

    integer(kind=kint) :: is1, ie1, is2, ie2

    is1 = hecMESH%node_group%grp_index(ig1-1) + 1
    ie1 = hecMESH%node_group%grp_index(ig1)
    is2 = hecMESH%node_group%grp_index(ig2-1) + 1
    ie2 = hecMESH%node_group%grp_index(ig2)
    mark( hecMESH%node_group%grp_item(is1:ie1) ) = .true.
    fstr_setup_inp_overlap = any( mark( hecMESH%node_group%grp_item(is2:ie2) ) )
    mark( hecMESH%node_group%grp_item(is1:ie1) ) = .false.
  end function fstr_setup_inp_overlap

  !> Sets the formulation of 361 and 362 of the sections by the inp element names, which hecmw1 keeps in sect_opt:
  !> C3D8 is computed with FBAR, C3D8I with IC, C3D20 with FI, and C3D8R and C3D20R with RI. Called after the
  !> defaults by the analysis type are set and before *SECTION is read, which overrides it.
  subroutine fstr_setup_inp_formulation( P )
    implicit none
    type(fstr_param_pack) :: P

    integer(kind=kint) :: i

    do i = 1, P%MESH%section%n_sect
      select case( P%MESH%section%sect_opt(i) )
      case( HECMW_SECT_OPT_SOLID_FI )
        P%SOLID%sections(i)%elemopt361 = kel361FI
        P%SOLID%sections(i)%elemopt362 = kel362FI
      case( HECMW_SECT_OPT_SOLID_IC )
        P%SOLID%sections(i)%elemopt361 = kel361IC
      case( HECMW_SECT_OPT_SOLID_RI )
        P%SOLID%sections(i)%elemopt361 = kel361RI
        P%SOLID%sections(i)%elemopt362 = kel362RI
      case( HECMW_SECT_OPT_SOLID_FBAR )
        P%SOLID%sections(i)%elemopt361 = kel361FBAR
      end select
    enddo
  end subroutine fstr_setup_inp_formulation

end module m_fstr_setup_inp
