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
  use fstr_setup_util
  use m_fstr_solid_setup, only &
    : fstr_setup_SECTION &
    , fstr_setup_WRITE_flags &
    , fstr_setup_WRITE_store &
    , fstr_setup_MATERIAL &
    , fstr_setup_EXPANSION_COEFF &
    , fstr_solid_append_boundary &
    , fstr_solid_append_cload &
    , fstr_solid_append_temperature
  use m_fstr_param_setup, only &
    : fstr_setup_SOLVER
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

    call fstr_setup_inp_read_cards( ctrl, T, 1, total, P )

    ! ----- counts stay 0: the materials are those of the mesh,
    ! and no other card of the table is stored in the arrays that fstr_setup_prepare allocates
    call fstr_setup_prepare( P, counts )

    ! -----
    rcode = fstr_ctrl_rewind( ctrl )
    call fstr_setup_inp_read_cards( ctrl, T, 2, total, P )

    call fstr_setup_finish( P, counts )
    rcode = fstr_ctrl_close( ctrl )

  end subroutine fstr_setup_inp

  !> Adds the cards of an inp file to the table T, by the names without the '*' of the header.
  !> The cards of the inp format are read in the first reading, so that fstr_setup_prepare knows the analysis type
  !> and all the boundary conditions. The material constants are read in the second reading into the materials
  !> that fstr_setup_prepare makes from the mesh.
  subroutine fstr_setup_inp_cards( T )
    implicit none
    type(fstr_card_table) :: T

    ! The cards of the inp format. The conditions go into the default step, as in a cnt file without !STEP.
    call fstr_card_add( T, 'STEP',     first=fstr_setup_inp_STEP )
    call fstr_card_add( T, 'END STEP' )
    call fstr_card_add( T, 'STATIC',   first=fstr_setup_inp_STATIC )
    call fstr_card_add( T, 'FREQUENCY', first=fstr_setup_inp_FREQUENCY )
    call fstr_card_add( T, 'BOUNDARY', first=fstr_setup_inp_BOUNDARY )
    call fstr_card_add( T, 'CLOAD',    first=fstr_setup_inp_CLOAD )
    call fstr_card_add( T, 'TEMPERATURE', first=fstr_setup_inp_TEMPERATURE )
    call fstr_card_add( T, 'INITIAL CONDITIONS', first=fstr_setup_inp_INITIAL_CONDITIONS )

    ! The material constants that hecmw1 does not read, stored in the material of the last *MATERIAL
    call fstr_card_add( T, 'MATERIAL',  second=fstr_setup_MATERIAL )
    call fstr_card_add( T, 'EXPANSION', second=fstr_setup_inp_EXPANSION )

    ! FrontISTR's own cards, read as the cnt cards
    call fstr_card_add( T, 'SOLVER',   first=fstr_setup_SOLVER )
    call fstr_card_add( T, 'SECTION',  second=fstr_setup_SECTION )
    call fstr_card_add( T, 'WRITE',    first=fstr_setup_inp_WRITE, second=fstr_setup_WRITE_store )

    ! Read by hecmw1
    call fstr_card_add( T, 'HEADING' )
    call fstr_card_add( T, 'NODE' )
    call fstr_card_add( T, 'ELEMENT' )
    call fstr_card_add( T, 'NSET' )
    call fstr_card_add( T, 'ELSET' )
    call fstr_card_add( T, 'ELASTIC' )
    call fstr_card_add( T, 'DENSITY' )
    call fstr_card_add( T, 'CONDUCTIVITY' )
    call fstr_card_add( T, 'SPECIFIC HEAT' )
    call fstr_card_add( T, 'SOLID SECTION' )
    call fstr_card_add( T, 'SHELL SECTION' )
    call fstr_card_add( T, 'EQUATION' )
    call fstr_card_add( T, 'AMPLITUDE' )

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
        case( '*STEP', '*END STEP', '*STATIC', '*FREQUENCY', '*BOUNDARY', '*CLOAD', &
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
              i = scan( param, '= ' )
              if( i > 0 ) param = param(:i-1)
              do i = 1, len_trim(param)
                c = iachar( param(i:i) )
                if( c >= iachar('a') .and. c <= iachar('z') ) param(i:i) = achar( c - iachar('a') + iachar('A') )
              enddo
              supported = .false.
              select case( trim(header_name) )
              case( '*STEP' )
                supported = param == 'NLGEOM'
              case( '*FREQUENCY' )
                supported = param == 'SOLVER' .or. param == 'STORAGE'
              case( '*BOUNDARY', '*CLOAD' )
                supported = param == 'AMPLITUDE'
              case( '*TEMPERATURE' )
                supported = param == 'AMPLITUDE' .or. param == 'USER' .or. param == 'SUBMODEL' .or. param == 'FILE'
              case( '*INITIAL CONDITIONS' )
                supported = param == 'TYPE'
              case( '*EXPANSION' )
                supported = param == 'TYPE' .or. param == 'ZERO'
              end select
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

  !-----------------------------------------------------------------------------!
  !> Read in *STEP                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_inp_STEP( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, is_nlgeom, nlgeom

    if( P%STATE%slot%first > 1 ) call fstr_setup_inp_err_stop( ctrl, 'more than one step is not supported yet' )

    is_nlgeom = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'NLGEOM ', '# ', 0, 'E', is_nlgeom )
    if( is_nlgeom == 0 ) return
    nlgeom = 1
    rcode = fstr_ctrl_get_param_ex( ctrl, 'NLGEOM ', 'YES,NO ', 0, 'P', nlgeom )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    if( nlgeom == 1 ) call fstr_setup_inp_err_stop( ctrl, 'NLGEOM is not supported yet' )
  end subroutine fstr_setup_inp_STEP

  !-----------------------------------------------------------------------------!
  !> Read in *STATIC                                                            !
  !-----------------------------------------------------------------------------!
  !> The data line of the time increments is not used in a linear analysis.
  subroutine fstr_setup_inp_STATIC( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    P%PARAM%solution_type = kstSTATIC
  end subroutine fstr_setup_inp_STATIC

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
  !> and the value (0 if omitted).
  subroutine fstr_setup_inp_BOUNDARY( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n
    character(len=HECMW_NAME_LEN) :: amp, rotc_name, data_fmt
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:)
    integer(kind=kint), allocatable :: dof_ids(:), dof_ide(:)
    real(kind=kreal), allocatable :: vals(:)

    amp = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AMPLITUDE ', '# ', 0, 'S', amp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
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

    call fstr_solid_append_boundary( P%MESH, P%SOLID, '*BOUNDARY', 1, 0, rotc_name, &
      amp, n, grp_names, dof_ids, dof_ide, vals )
  end subroutine fstr_setup_inp_BOUNDARY

  !-----------------------------------------------------------------------------!
  !> Read in *CLOAD                                                             !
  !-----------------------------------------------------------------------------!
  !> A data line is a node or a node set, the dof and the value.
  subroutine fstr_setup_inp_CLOAD( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode, n
    character(len=HECMW_NAME_LEN) :: amp, rotc_name, data_fmt
    character(len=HECMW_NAME_LEN), allocatable :: grp_names(:)
    integer(kind=kint), allocatable :: dofs(:)
    real(kind=kreal), allocatable :: vals(:)

    amp = ' '
    rcode = fstr_ctrl_get_param_ex( ctrl, 'AMPLITUDE ', '# ', 0, 'S', amp )
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    rotc_name = ' '

    n = fstr_ctrl_get_data_line_n( ctrl )
    allocate( grp_names(n), dofs(n), vals(n) )
    if( n > 0 ) then
      write( data_fmt, '(a,i0,a)' ) 'S', HECMW_NAME_LEN, 'IR '
      rcode = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, grp_names, dofs, vals )
      if( rcode /= 0 ) call fstr_ctrl_err_stop
    endif

    call fstr_solid_append_cload( P%MESH, P%SOLID, '*CLOAD', 1, rotc_name, &
      amp, n, grp_names, dofs, vals )
  end subroutine fstr_setup_inp_CLOAD

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

end module m_fstr_setup_inp
