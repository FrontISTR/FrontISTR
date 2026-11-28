!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to read in the analysis data of an inp file
!! and do necessary preparation for following calculation, like m_fstr_setup for a cnt file.
!! hecmw1 reads the mesh, the materials and the sections of the same file,
!! and the material constants that hecmw1 does not read are added here to the materials of the mesh.
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
    , fstr_setup_MATERIAL &
    , fstr_setup_EXPANSION_COEFF &
    , fstr_solid_append_boundary &
    , fstr_solid_append_cload &
    , fstr_solid_append_dload &
    , fstr_solid_append_temperature &
    , fstr_expand_dload_array
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

    ! ----- the other counts stay 0: the materials are those of the mesh,
    ! and the other cards are not stored in the arrays that fstr_setup_prepare allocates
    counts%steps          = fstr_card_total( T, total, 'STEP' )
    counts%autoinc_params = fstr_card_total( T, total, 'AUTOINC_PARAM' )
    counts%converg_params = fstr_card_total( T, total, 'CONVERG_PARAM' )
    if( counts%steps > 0 .and. fstr_card_total( T, total, 'analysis' ) /= counts%steps ) &
      call fstr_setup_util_err_stop( '### Error: *STEP : each step needs one analysis card' )
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
  !> the steps and all the boundary conditions. The second reading stores the steps that fstr_setup_prepare allocates,
  !> and the material constants into the materials that fstr_setup_prepare makes from the mesh.
  subroutine fstr_setup_inp_cards( T )
    implicit none
    type(fstr_card_table) :: T

    ! The cards of the inp format. Without *STEP, the conditions go into the default step, as in a cnt file
    ! without !STEP. *STEP also refers by name to *AUTOINC_PARAM and *CONVERG_PARAM read before it.
    ! The analysis cards share a sequence, so that a step without one or with two is found.
    call fstr_card_add( T, 'STEP',     first=fstr_setup_inp_STEP, second=fstr_setup_inp_STEP_store )
    call fstr_card_add( T, 'END STEP' )
    call fstr_card_add( T, 'STATIC',   first=fstr_setup_inp_STATIC, second=fstr_setup_inp_STATIC_store, seq='analysis' )
    call fstr_card_add( T, 'FREQUENCY', first=fstr_setup_inp_FREQUENCY, seq='analysis' )
    call fstr_card_add( T, 'BOUNDARY', first=fstr_setup_inp_BOUNDARY )
    call fstr_card_add( T, 'CLOAD',    first=fstr_setup_inp_CLOAD )
    call fstr_card_add( T, 'DLOAD',    first=fstr_setup_inp_DLOAD )
    call fstr_card_add( T, 'DSLOAD',   first=fstr_setup_inp_DSLOAD )
    call fstr_card_add( T, 'TEMPERATURE', first=fstr_setup_inp_TEMPERATURE )
    call fstr_card_add( T, 'INITIAL CONDITIONS', first=fstr_setup_inp_INITIAL_CONDITIONS )

    ! The material constants that hecmw1 does not read, stored in the material of the last *MATERIAL
    call fstr_card_add( T, 'MATERIAL',  second=fstr_setup_MATERIAL )
    call fstr_card_add( T, 'EXPANSION', second=fstr_setup_inp_EXPANSION )

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
    call fstr_card_add( T, 'SURFACE' )
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
  !> of the inp format that are not read yet would change the analysis. A parameter of a card of the analysis data that
  !> its reader does not read is skipped with a warning.
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
        case( '*STEP', '*END STEP', '*STATIC', '*FREQUENCY', '*BOUNDARY', '*CLOAD', '*DLOAD', '*DSLOAD', &
            '*TEMPERATURE', '*INITIAL CONDITIONS', '*EXPANSION' )
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
              case( '*FREQUENCY' )
                supported = param == 'SOLVER' .or. param == 'STORAGE'
              case( '*BOUNDARY' )
                supported = param == 'AMPLITUDE' .or. param == 'OP' .or. param == 'FIXED'
              case( '*CLOAD', '*DLOAD' )
                supported = param == 'AMPLITUDE' .or. param == 'OP'
              case( '*DSLOAD' )
                supported = param == 'AMPLITUDE'
              case( '*TEMPERATURE' )
                supported = param == 'AMPLITUDE' .or. param == 'USER' .or. param == 'SUBMODEL' .or. param == 'FILE'
              case( '*INITIAL CONDITIONS' )
                supported = param == 'TYPE'
              case( '*EXPANSION' )
                supported = param == 'TYPE' .or. param == 'ZERO'
              end select
              ! FrontISTR cannot shift the time of an amplitude
              if( param == 'TIME DELAY' .and. ( header_name == '*BOUNDARY' .or. header_name == '*CLOAD' &
                .or. header_name == '*DLOAD' .or. header_name == '*DSLOAD' ) ) &
                call fstr_setup_inp_err_stop( ctrl, 'TIME DELAY is not supported yet' )
              if( ( header_name == '*DLOAD' .and. ( param == 'LOAD CASE' .or. param == 'SECTOR' ) ) .or. &
                ( header_name == '*DSLOAD' .and. param == 'SUBMODEL' ) ) &
                call fstr_setup_inp_err_stop( ctrl, trim(param)//' is not supported yet' )
              if( .not. supported ) call fstr_setup_inp_warn( ctrl, P, &
                'unsupported parameter '//trim(param)//' is skipped' )
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

  !> Writes msg about the card at the current header of ctrl and its line in the inp file as a warning
  subroutine fstr_setup_inp_warn( ctrl, P, msg )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    character(len=*), intent(in) :: msg

    integer(kind=kint) :: rcode
    character(len=HECMW_NAME_LEN) :: header_name

    if( P%MESH%my_rank /= 0 ) return
    rcode = fstr_ctrl_get_c_h_name_inp( ctrl, header_name, HECMW_NAME_LEN )
    write(*,'(a,a,a,i0,a,a)') ' Warning: ', trim(header_name), ' (line:', fstr_ctrl_get_c_h_line_no( ctrl ), &
      ') : ', msg
  end subroutine fstr_setup_inp_warn

  !> Stops unless the analysis card at the current header of ctrl is the only one in its step, and each step before
  !> has one. The steps after it are checked by the total of the analysis cards in fstr_setup_inp.
  subroutine fstr_setup_inp_analysis_card( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    if( P%STATE%slot%first > max( P%STATE%step, 1 ) ) &
      call fstr_setup_inp_err_stop( ctrl, 'only one analysis card is allowed in a step' )
    if( P%STATE%slot%first < P%STATE%step ) call fstr_setup_inp_err_stop( ctrl, 'each step needs one analysis card' )
  end subroutine fstr_setup_inp_analysis_card

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
      if( P%SOLID%DLOAD_ngrp_tot > 0 ) then
        if( any( P%SOLID%DLOAD_ngrp_amp(1:P%SOLID%DLOAD_ngrp_tot) > 0 ) ) call fstr_setup_inp_amp_err_stop( ctrl )
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

    call fstr_setup_inp_analysis_card( ctrl, P )
    if( P%PARAM%solution_type == kstEIGEN ) &
      call fstr_setup_inp_err_stop( ctrl, 'a step after a step of *FREQUENCY is not supported yet' )

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
  !> Read in *FREQUENCY                                                         !
  !-----------------------------------------------------------------------------!
  !> The data line is the number of eigenfrequencies, and the lower and the upper value of a frequency range.
  !> The lowest eigenfrequencies are computed whatever the range is, since the range of the inp format only selects
  !> the ones written out. All of them are written to the result file.
  subroutine fstr_setup_inp_FREQUENCY( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, storage
    real(kind=kreal) :: fmin, fmax
    character(len=HECMW_NAME_LEN) :: solver

    call fstr_setup_inp_analysis_card( ctrl, P )
    if( P%STATE%step > 1 ) call fstr_setup_inp_err_stop( ctrl, '*FREQUENCY after the first step is not supported yet' )

    solver = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'SOLVER ', '# ', 0, 'S', solver )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( solver == 'MATRIXSTORAGE' ) then
      call fstr_setup_inp_err_stop( ctrl, 'SOLVER=MATRIXSTORAGE is not supported' )
    else if( solver /= ' ' ) then
      call fstr_setup_inp_warn( ctrl, P, 'SOLVER='//trim(solver)// &
        ' is skipped. The linear solver is given by *SOLVER, DIRECTmf by default' )
    endif

    storage = 2
    rcode = fstr_ctrl_get_param_ex( ctrl, 'STORAGE ', 'YES,NO ', 0, 'P', storage )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( storage == 1 ) call fstr_setup_inp_warn( ctrl, P, &
      'STORAGE=YES is skipped. The eigenmodes are written to the result file of FrontISTR' )

    fmin = 0.0d0
    fmax = -1.0d0
    rcode = fstr_ctrl_get_data_ex( ctrl, 1, 'Irr ', P%EIGEN%nget, fmin, fmax )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( P%EIGEN%nget < 1 ) call fstr_setup_inp_err_stop( ctrl, 'the number of eigenfrequencies must be positive' )
    if( fmin > 0.0d0 .or. fmax >= 0.0d0 ) call fstr_setup_inp_warn( ctrl, P, &
      'the frequency range is skipped. All the eigenmodes computed are written to the result file' )

    P%PARAM%solution_type = kstEIGEN
  end subroutine fstr_setup_inp_FREQUENCY

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
  !> Read in *DLOAD                                                             !
  !-----------------------------------------------------------------------------!
  !> A data line is an element or an element set, the load label and its values: Pn and the pressure on the face n
  !> (P on a shell), GRAV and the size and the direction of the acceleration, or CENTRIF, the square of the rotational
  !> speed, a point on the axis and the direction of the axis. They are the loads 10*n, 4 and 5 of !DLOAD, which
  !> takes the rotational speed and normalizes the direction of GRAV only for the solids.
  subroutine fstr_setup_inp_DLOAD( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n, i, op
    character(len=HECMW_NAME_LEN) :: amp, fmt_label
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:), labels(:)
    integer(kind=kint), allocatable :: lids(:)
    real(kind=kreal), allocatable :: params(:,:)
    real(kind=kreal) :: norm

    amp = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AMPLITUDE ', '# ', 0, 'S', amp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( amp /= ' ' .and. P%STATE%step > 1 ) call fstr_setup_inp_amp_err_stop( ctrl )
    op = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'OP ', 'MOD,NEW ', 0, 'P', op )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_names(n), labels(n), lids(n), params(0:6,n) )
    params = 0.0d0
    write( fmt_label, '(a,i0,a,i0,a)' ) 'S', HECMW_NAME_LEN, 'S', HECMW_NAME_LEN, ' '
    do i = 1, n
      rcode = fstr_ctrl_get_data_ex( ctrl, i, fmt_label, grp_names(i), labels(i) )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
      select case( labels(i) )
      case( 'P', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6' )
        lids(i) = 10
        if( labels(i) /= 'P' ) lids(i) = 10*( iachar(labels(i)(2:2)) - iachar('0') )
        rcode = fstr_ctrl_get_data_ex( ctrl, i, trim(fmt_label)//'R ', grp_names(i), labels(i), params(0,i) )
      case( 'GRAV' )
        lids(i) = 4
        rcode = fstr_ctrl_get_data_ex( ctrl, i, trim(fmt_label)//'RRRR ', grp_names(i), labels(i), &
          params(0,i), params(1,i), params(2,i), params(3,i) )
      case( 'CENTRIF' )
        lids(i) = 5
        rcode = fstr_ctrl_get_data_ex( ctrl, i, trim(fmt_label)//'RRRRRRR ', grp_names(i), labels(i), &
          params(0,i), params(1,i), params(2,i), params(3,i), params(4,i), params(5,i), params(6,i) )
      case default
        call fstr_setup_inp_err_stop( ctrl, 'the load label '//trim(labels(i))//' is not supported yet' )
      end select
      if( rcode /= 0 ) call fstr_ctrl_err_stop
      if( lids(i) == 4 ) then
        norm = sqrt( sum( params(1:3,i)**2 ) )
        if( norm <= 0.0d0 ) call fstr_setup_inp_err_stop( ctrl, 'the direction of GRAV must not be zero' )
        params(1:3,i) = params(1:3,i) / norm
      else if( lids(i) == 5 ) then
        if( params(0,i) < 0.0d0 ) call fstr_setup_inp_err_stop( ctrl, &
          'the square of the rotational speed of CENTRIF must not be negative' )
        if( sum( params(4:6,i)**2 ) <= 0.0d0 ) call fstr_setup_inp_err_stop( ctrl, &
          'the direction of the axis of CENTRIF must not be zero' )
        params(0,i) = sqrt( params(0,i) )
      endif
    enddo

    call fstr_setup_inp_dload_append( ctrl, P, '*DLOAD', op, amp, n, grp_names, labels, lids, params )
  end subroutine fstr_setup_inp_DLOAD

  !-----------------------------------------------------------------------------!
  !> Read in *DSLOAD                                                            !
  !-----------------------------------------------------------------------------!
  !> A data line is a surface, the load label P and the pressure, the load S of !DLOAD on the surface group
  subroutine fstr_setup_inp_DSLOAD( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n, i
    character(len=HECMW_NAME_LEN) :: amp, data_fmt
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:), labels(:)
    integer(kind=kint), allocatable :: lids(:)
    real(kind=kreal), allocatable :: params(:,:)

    amp = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AMPLITUDE ', '# ', 0, 'S', amp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( amp /= ' ' .and. P%STATE%step > 1 ) call fstr_setup_inp_amp_err_stop( ctrl )

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_names(n), labels(n), lids(n), params(0:6,n) )
    params = 0.0d0
    lids = 100
    write( data_fmt, '(a,i0,a,i0,a)' ) 'S', HECMW_NAME_LEN, 'S', HECMW_NAME_LEN, 'R '
    do i = 1, n
      rcode = fstr_ctrl_get_data_ex( ctrl, i, data_fmt, grp_names(i), labels(i), params(0,i) )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
      if( labels(i) /= 'P' ) call fstr_setup_inp_err_stop( ctrl, &
        'the load label '//trim(labels(i))//' is not supported yet' )
    enddo

    call fstr_setup_inp_dload_append( ctrl, P, '*DSLOAD', 1, amp, n, grp_names, labels, lids, params )
  end subroutine fstr_setup_inp_DSLOAD

  !> Appends the n loads of the *DLOAD or the *DSLOAD at the current header of ctrl to the loads of !DLOAD. The group
  !> (GRPID) of a load is the step it starts in, and a load is active to the end, as in fstr_setup_inp_CLOAD. A load
  !> replacing those of the earlier steps (fstr_setup_inp_same_dload), and OP=NEW of the first *DLOAD in the step,
  !> add the loads that cancel them by fstr_setup_inp_cancel_dload. The pressures follow the deformation in a
  !> geometrically nonlinear analysis, as the default FOLLOW=YES of !DLOAD.
  subroutine fstr_setup_inp_dload_append( ctrl, P, header_name, op_given, amp, n, grp_names, labels, lids, params )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    character(len=*), intent(in) :: header_name
    integer(kind=kint), intent(in) :: op_given, n
    character(len=HECMW_NAME_LEN), intent(in) :: amp, grp_names(:), labels(:)
    integer(kind=kint), intent(in) :: lids(:)
    real(kind=kreal), intent(in) :: params(0:,:)

    integer(kind=kint) :: op, step, follow, old_size, e, f
    integer(kind=kint), allocatable :: mark(:)

    op = op_given
    step = max( P%STATE%step, 1 )
    follow = 0
    if( P%PARAM%nlgeom ) follow = 1
    old_size = P%SOLID%DLOAD_ngrp_tot
    if( P%STATE%step > 1 ) then
      allocate( mark(P%MESH%n_elem) )
      mark = 0
    endif
    if( op == 2 .and. P%STATE%step > 1 ) &
      call fstr_setup_inp_first_op( ctrl, P, header_name, P%SOLID%DLOAD_ngrp_GRPID, old_size, step, op )
    if( op == 2 .and. P%STATE%step > 1 ) then
      do f = 1, old_size
        if( P%SOLID%DLOAD_ngrp_GRPID(f) >= step ) cycle
        call fstr_setup_inp_cancel_dload( ctrl, P, f, step, P%SOLID%DLOAD_ngrp_tot + 1, mark )
      enddo
    endif
    old_size = P%SOLID%DLOAD_ngrp_tot
    call fstr_solid_append_dload( P%MESH, P%SOLID, header_name, step, follow, amp, n, grp_names, lids, params )
    call fstr_setup_inp_check_dload( ctrl, P, old_size, labels )

    if( P%STATE%step <= 1 ) return
    do e = old_size + 1, old_size + n
      call fstr_setup_inp_cancel_dload( ctrl, P, e, step, e, mark )
    enddo
  end subroutine fstr_setup_inp_dload_append

  !> Stops at a load of the card at the current header of ctrl, the loads old_size+1 to the last of !DLOAD given by
  !> labels, on an element that FrontISTR cannot load so: process_distributed_loads skips the elements other than
  !> those below, DL_Shell takes any load of a face of a shell as P, and DL_C2 integrates over the thickness 1 and
  !> takes only the plane stress, the plane strain and the axisymmetric options, in the x-y plane.
  !> An amplitude of VALUE=ABSOLUTE would replace the rotational speed of CENTRIF, not its square.
  subroutine fstr_setup_inp_check_dload( ctrl, P, old_size, labels )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: old_size
    character(len=HECMW_NAME_LEN), intent(in) :: labels(:)

    integer(kind=kint) :: i, j, is, ie, ig, lid, icel, face, nface, isect, amp
    real(kind=kreal) :: params(0:6), thick
    character(len=HECMW_MSG_LEN) :: msg, reason
    character(len=HECMW_NAME_LEN) :: label

    do i = old_size + 1, P%SOLID%DLOAD_ngrp_tot
      label = labels(i - old_size)
      lid = P%SOLID%DLOAD_ngrp_LID(i)
      ig = P%SOLID%DLOAD_ngrp_ID(i)
      params(:) = P%SOLID%DLOAD_ngrp_params(:,i)
      amp = P%SOLID%DLOAD_ngrp_amp(i)
      if( lid == 5 .and. amp > 0 ) then
        if( P%MESH%amp%amp_type_value(amp) == HECMW_AMP_TYPEVAL_ABSOLUTE ) &
          call fstr_setup_inp_err_stop( ctrl, 'an amplitude of VALUE=ABSOLUTE with CENTRIF is not supported yet' )
      endif
      if( lid == 100 ) then
        is = P%MESH%surf_group%grp_index(ig-1) + 1
        ie = P%MESH%surf_group%grp_index(ig)
      else
        is = P%MESH%elem_group%grp_index(ig-1) + 1
        ie = P%MESH%elem_group%grp_index(ig)
      endif
      do j = is, ie
        if( lid == 100 ) then
          icel = P%MESH%surf_group%grp_item(2*j-1)
          face = P%MESH%surf_group%grp_item(2*j)
        else
          icel = P%MESH%elem_group%grp_item(j)
          face = lid / 10
        endif
        select case( P%MESH%elem_type(icel) )
        case( 341, 342 )
          nface = 4
        case( 351, 352 )
          nface = 5
        case( 361, 362 )
          nface = 6
        case( 231, 232 )
          nface = 3
        case( 241, 242 )
          nface = 4
        case( 731, 741 )
          nface = 0
        case default
          nface = -1
        end select
        reason = ' '
        if( nface < 0 ) then
          reason = 'the element type is not supported yet'
        else if( nface == 0 ) then
          if( lid == 100 .or. ( lid >= 10 .and. label /= 'P' ) ) reason = 'the pressure on a shell is given only by P'
        else if( lid /= 100 .and. label == 'P' ) then
          reason = 'P is the pressure on a shell; give the face of the element by Pn'
        else if( face > nface ) then
          reason = 'the element has no such face'
        else if( P%MESH%elem_type(icel) < 300 ) then
          isect = P%MESH%section_ID(icel)
          thick = P%MESH%section%sect_R_item( P%MESH%section%sect_R_index(isect-1) + 1 )
          select case( P%MESH%section%sect_opt(isect) )
          case( HECMW_SECT_OPT_PSTRESS, HECMW_SECT_OPT_PSTRAIN )
            if( abs( thick - 1.0d0 ) > 0.0d0 ) &
              reason = 'a two-dimensional element of a thickness other than 1 is not supported yet'
          case( HECMW_SECT_OPT_ASYMMETRY )
          case default
            reason = 'a two-dimensional element of the reduced integration is not supported yet'
          end select
          if( lid == 4 .and. abs( params(3) ) > 0.0d0 ) &
            reason = 'the direction of GRAV on a two-dimensional element must be in the x-y plane'
          if( lid == 5 .and. abs( params(6) ) > 0.0d0 ) &
            reason = 'the axis of CENTRIF on a two-dimensional element must be in the x-y plane'
        endif
        if( reason /= ' ' ) then
          write( msg, '(a,i0,a,i0,a,a)' ) 'the load '//trim(label)//' on element ', P%MESH%global_elem_ID(icel), &
            ' (type ', P%MESH%elem_type(icel), ') : ', trim(reason)
          call fstr_setup_inp_err_stop( ctrl, trim(msg) )
        endif
      enddo
    enddo
  end subroutine fstr_setup_inp_check_dload

  !> Applies op to the faces of the load i of !DLOAD, the bits of mark at the elements set at the face numbers (0 for
  !> GRAV and CENTRIF): op 1 sets them and op 2 clears them, and the number of the faces is returned. op 3 returns the
  !> number of the faces of i whose bits are set.
  integer(kind=kint) function fstr_setup_inp_dload_faces( P, i, op, mark )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: i, op
    integer(kind=kint) :: mark(:)

    integer(kind=kint) :: ig, lid, j, is, ie, icel, face

    ig = P%SOLID%DLOAD_ngrp_ID(i)
    lid = P%SOLID%DLOAD_ngrp_LID(i)
    if( lid == 100 ) then
      is = P%MESH%surf_group%grp_index(ig-1) + 1
      ie = P%MESH%surf_group%grp_index(ig)
    else
      is = P%MESH%elem_group%grp_index(ig-1) + 1
      ie = P%MESH%elem_group%grp_index(ig)
    endif
    fstr_setup_inp_dload_faces = 0
    do j = is, ie
      if( lid == 100 ) then
        icel = P%MESH%surf_group%grp_item(2*j-1)
        face = P%MESH%surf_group%grp_item(2*j)
      else
        icel = P%MESH%elem_group%grp_item(j)
        face = lid / 10
      endif
      select case( op )
      case( 1 )
        mark(icel) = ibset( mark(icel), face )
        fstr_setup_inp_dload_faces = fstr_setup_inp_dload_faces + 1
      case( 2 )
        mark(icel) = ibclr( mark(icel), face )
        fstr_setup_inp_dload_faces = fstr_setup_inp_dload_faces + 1
      case( 3 )
        if( btest( mark(icel), face ) ) fstr_setup_inp_dload_faces = fstr_setup_inp_dload_faces + 1
      end select
    enddo
  end function fstr_setup_inp_dload_faces

  !> Whether a load j of !DLOAD of a later step replaces the load i, as in the inp format: both are pressures on the
  !> same faces, GRAV in the same direction or CENTRIF about the same axis on the same elements. mark has the bits of
  !> the faces of i (fstr_setup_inp_dload_faces), and nface is their number.
  logical function fstr_setup_inp_same_dload( P, i, j, nface, mark )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: i, j, nface
    integer(kind=kint) :: mark(:)

    integer(kind=kint) :: li, lj
    real(kind=kreal) :: pi(0:6), pj(0:6), axis(3), cross(3)

    fstr_setup_inp_same_dload = .false.
    li = P%SOLID%DLOAD_ngrp_LID(i)
    lj = P%SOLID%DLOAD_ngrp_LID(j)
    if( min( li, 10 ) /= min( lj, 10 ) ) return
    pi(:) = P%SOLID%DLOAD_ngrp_params(:,i)
    pj(:) = P%SOLID%DLOAD_ngrp_params(:,j)
    if( li == 4 ) then
      if( any( abs( pi(1:3) - pj(1:3) ) > 1.0d-12 ) ) return
    else if( li == 5 ) then
      axis = pi(4:6) / sqrt( sum( pi(4:6)**2 ) )
      cross = (/ axis(2)*pj(6) - axis(3)*pj(5), axis(3)*pj(4) - axis(1)*pj(6), axis(1)*pj(5) - axis(2)*pj(4) /)
      if( any( abs( cross ) > 1.0d-12*sqrt( sum( pj(4:6)**2 ) ) ) ) return
      pj(1:3) = pj(1:3) - pi(1:3)
      cross = (/ axis(2)*pj(3) - axis(3)*pj(2), axis(3)*pj(1) - axis(1)*pj(3), axis(1)*pj(2) - axis(2)*pj(1) /)
      if( any( abs( cross ) > 1.0d-12*max( sqrt( sum( pi(1:3)**2 ) ), sqrt( sum( pj(1:3)**2 ) ), 1.0d0 ) ) ) return
    endif
    fstr_setup_inp_same_dload =fstr_setup_inp_dload_faces( P, j, 3, mark ) == nface .and. &
      fstr_setup_inp_dload_size( P, j ) == nface
  end function fstr_setup_inp_same_dload

  !> The number of the faces or the elements of the load i of !DLOAD
  integer(kind=kint) function fstr_setup_inp_dload_size( P, i )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: i

    integer(kind=kint) :: ig

    ig = P%SOLID%DLOAD_ngrp_ID(i)
    if( P%SOLID%DLOAD_ngrp_LID(i) == 100 ) then
      fstr_setup_inp_dload_size = P%MESH%surf_group%grp_index(ig) - P%MESH%surf_group%grp_index(ig-1)
    else
      fstr_setup_inp_dload_size = P%MESH%elem_group%grp_index(ig) - P%MESH%elem_group%grp_index(ig-1)
    endif
  end function fstr_setup_inp_dload_size

  !> The value of the load i of !DLOAD that the loads replacing each other sum: the pressure, the size of GRAV, or the
  !> square of the rotational speed of CENTRIF
  real(kind=kreal) function fstr_setup_inp_dload_value( P, i )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: i

    fstr_setup_inp_dload_value = P%SOLID%DLOAD_ngrp_params(0,i)
    if( P%SOLID%DLOAD_ngrp_LID(i) == 5 ) fstr_setup_inp_dload_value = fstr_setup_inp_dload_value*abs( fstr_setup_inp_dload_value )
  end function fstr_setup_inp_dload_value

  !> Adds the load that cancels, from step on, the loads of the earlier steps that the load src replaces
  !> (fstr_setup_inp_same_dload), as a load of the group step on the group of src whose value is minus their sum
  !> (fstr_setup_inp_dload_value). Nothing is added if a load of step before the load upto replaces them, since they
  !> are cancelled already. mark is all 0 on entry and on return.
  !> DL_C2 and DL_Shell take the square of the rotational speed of CENT, which a negative speed cannot cancel.
  subroutine fstr_setup_inp_cancel_dload( ctrl, P, src, step, upto, mark )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: src, step, upto
    integer(kind=kint) :: mark(:)

    integer(kind=kint) :: f, n, nface, ig, j, etype
    real(kind=kreal) :: val, scale
    logical :: replaced

    nface = fstr_setup_inp_dload_faces( P, src, 1, mark )
    val = 0.0d0
    scale = 0.0d0
    replaced = .false.
    do f = 1, P%SOLID%DLOAD_ngrp_tot
      if( P%SOLID%DLOAD_ngrp_GRPID(f) > step ) cycle
      if( P%SOLID%DLOAD_ngrp_GRPID(f) == step .and. f >= upto ) cycle
      if( .not. fstr_setup_inp_same_dload( P, src, f, nface, mark ) ) cycle
      if( P%SOLID%DLOAD_ngrp_GRPID(f) == step ) then
        replaced = .true.
        exit
      endif
      val = val + fstr_setup_inp_dload_value( P, f )
      scale = scale + abs( fstr_setup_inp_dload_value( P, f ) )
    enddo
    nface = fstr_setup_inp_dload_faces( P, src, 2, mark )
    if( replaced .or. fstr_setup_inp_cancelled( val, scale ) ) return
    if( P%SOLID%DLOAD_ngrp_LID(src) == 5 ) then
      ig = P%SOLID%DLOAD_ngrp_ID(src)
      do j = P%MESH%elem_group%grp_index(ig-1) + 1, P%MESH%elem_group%grp_index(ig)
        etype = P%MESH%elem_type( P%MESH%elem_group%grp_item(j) )
        if( etype < 300 .or. etype > 700 ) call fstr_setup_inp_err_stop( ctrl, &
          'changing or removing CENTRIF on a two-dimensional element or a shell in a later step is not supported yet' )
      enddo
    endif

    n = P%SOLID%DLOAD_ngrp_tot + 1
    call hecmw_expand_integer_array( P%SOLID%DLOAD_ngrp_GRPID, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%DLOAD_ngrp_ID, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%DLOAD_ngrp_LID, n-1, n )
    call hecmw_expand_integer_array( P%SOLID%DLOAD_ngrp_amp, n-1, n )
    call fstr_expand_dload_array( P%SOLID%DLOAD_ngrp_params, n-1, n )
    P%SOLID%DLOAD_ngrp_tot = n
    P%SOLID%DLOAD_ngrp_GRPID(n) = step
    P%SOLID%DLOAD_ngrp_ID(n) = P%SOLID%DLOAD_ngrp_ID(src)
    P%SOLID%DLOAD_ngrp_LID(n) = P%SOLID%DLOAD_ngrp_LID(src)
    P%SOLID%DLOAD_ngrp_amp(n) = P%SOLID%DLOAD_ngrp_amp(src)
    P%SOLID%DLOAD_ngrp_params(:,n) = P%SOLID%DLOAD_ngrp_params(:,src)
    if( P%SOLID%DLOAD_ngrp_LID(src) == 5 ) then
      P%SOLID%DLOAD_ngrp_params(0,n) = -sign( sqrt( abs(val) ), val )
    else
      P%SOLID%DLOAD_ngrp_params(0,n) = -val
    endif
  end subroutine fstr_setup_inp_cancel_dload

  !-----------------------------------------------------------------------------!
  !> Read in *TEMPERATURE                                                       !
  !-----------------------------------------------------------------------------!
  !> A data line is a node or a node set and its temperature at the end of the step. The temperature gradients of
  !> shells and beams after it are not read. The nodes not given keep their initial temperature.
  subroutine fstr_setup_inp_TEMPERATURE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n, i, is_given
    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:)
    real(kind=kreal), allocatable :: vals(:)
    character(len=9), parameter :: params(4) = [ 'AMPLITUDE', 'USER     ', 'SUBMODEL ', 'FILE     ' ]

    do i = 1, size(params)
      is_given = 0
      rcode = fstr_ctrl_get_param_ex( ctrl, params(i)//' ', '# ', 0, 'E', is_given )
      if( is_given /= 0 ) call fstr_setup_inp_err_stop( ctrl, trim(params(i))//' is not supported yet' )
    end do
    if( P%STATE%step > 1 ) call fstr_setup_inp_err_stop( ctrl, 'given in a step after the first step is not supported yet' )
    call fstr_setup_inp_check_constant_materials( ctrl, P )

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_names(n), vals(n) )
    if( n > 0 ) then
      write( data_fmt, '(a,i0,a)' ) 'S', HECMW_NAME_LEN, 'R '
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_names, vals )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif

    call fstr_solid_append_temperature( P%MESH, P%SOLID, '*TEMPERATURE', 1, n, grp_names, vals )
  end subroutine fstr_setup_inp_TEMPERATURE

  !-----------------------------------------------------------------------------!
  !> Read in *INITIAL CONDITIONS                                                !
  !-----------------------------------------------------------------------------!
  !> hecmw1 reads the initial temperatures of TYPE=TEMPERATURE into the mesh, and FrontISTR uses them as given to every node.
  subroutine fstr_setup_inp_INITIAL_CONDITIONS( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n
    character(len=HECMW_NAME_LEN) :: type_name

    type_name = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', '# ', 0, 'S', type_name )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( type_name /= 'TEMPERATURE' ) &
      call fstr_setup_inp_err_stop( ctrl, 'TYPE='//trim(type_name)//' is not supported yet' )

    n = 0
    if( P%MESH%hecmw_flag_initcon == 1 ) n = P%MESH%node_init_val_index(P%MESH%n_node)
    if( n /= P%MESH%n_node ) call fstr_setup_inp_err_stop( ctrl, 'the initial temperature must be given to all the nodes' )
    call fstr_setup_inp_check_constant_materials( ctrl, P )
  end subroutine fstr_setup_inp_INITIAL_CONDITIONS

  !> Stops at the card of ctrl if a material constant read by hecmw1 has a table of temperatures, since FrontISTR
  !> takes one value independent of the temperature from the materials of the mesh
  subroutine fstr_setup_inp_check_constant_materials( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: i

    do i = 1, P%MESH%material%n_mat_subitem
      if( P%MESH%material%mat_table_index(i) - P%MESH%material%mat_table_index(i-1) > 1 ) then
        call fstr_setup_inp_err_stop( ctrl, &
          'not supported yet with a material constant depending on the temperature (*ELASTIC, *DENSITY)' )
      endif
    end do
  end subroutine fstr_setup_inp_check_constant_materials

  !-----------------------------------------------------------------------------!
  !> Read in *EXPANSION                                                         !
  !-----------------------------------------------------------------------------!
  !> Stored by the reader of !EXPANSION_COEF in the material of the last *MATERIAL. Only a coefficient independent of
  !> the temperature is read. Its thermal strain alpha (T - T0) does not depend on ZERO, the reference temperature of
  !> the coefficient, which is not used.
  subroutine fstr_setup_inp_EXPANSION( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, is_type

    if( P%STATE%cid == 0 ) call fstr_setup_inp_err_stop( ctrl, 'given before *MATERIAL' )
    is_type = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', '# ', 0, 'E', is_type )
    if( is_type /= 0 ) call fstr_setup_inp_err_stop( ctrl, &
      'TYPE is not supported yet; omit it for an isotropic coefficient' )
    if( fstr_ctrl_get_data_line_n( ctrl ) /= 1 ) call fstr_setup_inp_err_stop( ctrl, &
      'only a coefficient independent of the temperature, in one data line, is supported yet' )

    call fstr_setup_EXPANSION_COEFF( ctrl, P )
  end subroutine fstr_setup_inp_EXPANSION

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
    call fstr_setup_inp_dload_overlap( P )

    do e = 1, P%SOLID%BOUNDARY_ngrp_tot
      s = P%SOLID%BOUNDARY_ngrp_GRPID(e)
      last = P%STATE%boundary_end(e)
      if( last == 0 ) last = nstep + 1
      P%SOLID%BOUNDARY_ngrp_GRPID(e) = fstr_setup_inp_group_id( s, last, nstep )
    enddo
    if( associated(P%STATE%boundary_end) ) deallocate( P%STATE%boundary_end )
  end subroutine fstr_setup_inp_groups

  !> Stops at a pressure of a later step on some faces of a pressure of an earlier step on other faces, unless the
  !> earlier one is cancelled in the step, as for the loads of *CLOAD in fstr_setup_inp_groups: the inp format
  !> replaces the pressure on each face, while a load of !DLOAD has a value for all the faces of its group.
  subroutine fstr_setup_inp_dload_overlap( P )
    implicit none
    type(fstr_param_pack) :: P

    integer(kind=kint) :: e, f, g, k, ne, nf
    integer(kind=kint), allocatable :: mark_e(:), mark_f(:)
    real(kind=kreal) :: val, scale
    character(len=HECMW_MSG_LEN) :: msg

    if( P%SOLID%DLOAD_ngrp_tot == 0 ) return
    allocate( mark_e(P%MESH%n_elem), mark_f(P%MESH%n_elem) )
    mark_e = 0
    mark_f = 0
    do e = 1, P%SOLID%DLOAD_ngrp_tot
      k = P%SOLID%DLOAD_ngrp_GRPID(e)
      if( k == 1 .or. P%SOLID%DLOAD_ngrp_LID(e) < 10 ) cycle
      ne = fstr_setup_inp_dload_faces( P, e, 1, mark_e )
      do f = 1, P%SOLID%DLOAD_ngrp_tot
        if( P%SOLID%DLOAD_ngrp_GRPID(f) >= k .or. P%SOLID%DLOAD_ngrp_LID(f) < 10 ) cycle
        if( fstr_setup_inp_same_dload( P, e, f, ne, mark_e ) ) cycle
        if( fstr_setup_inp_dload_faces( P, f, 3, mark_e ) == 0 ) cycle
        ! the pressures on the faces of f, unless a load of step k cancels them
        nf = fstr_setup_inp_dload_faces( P, f, 1, mark_f )
        val = 0.0d0
        scale = 0.0d0
        do g = 1, P%SOLID%DLOAD_ngrp_tot
          if( P%SOLID%DLOAD_ngrp_GRPID(g) > k ) cycle
          if( .not. fstr_setup_inp_same_dload( P, f, g, nf, mark_f ) ) cycle
          if( P%SOLID%DLOAD_ngrp_GRPID(g) == k ) then
            val = 0.0d0
            exit
          endif
          val = val + P%SOLID%DLOAD_ngrp_params(0,g)
          scale = scale + abs( P%SOLID%DLOAD_ngrp_params(0,g) )
        enddo
        nf = fstr_setup_inp_dload_faces( P, f, 2, mark_f )
        if( fstr_setup_inp_cancelled( val, scale ) ) cycle
        write( msg, '(a,i0,a)' ) '### Error: *DLOAD : a pressure of step ', k, ' is given to some faces of a '// &
          'pressure of an earlier step given to other faces, which is not supported yet'
        call fstr_setup_util_err_stop( trim(msg) )
      enddo
      ne = fstr_setup_inp_dload_faces( P, e, 2, mark_e )
    enddo
  end subroutine fstr_setup_inp_dload_overlap

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
