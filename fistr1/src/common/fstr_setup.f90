!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to read in data from control file
!! and do necessary preparation for following calculation
module m_fstr_setup
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex &
    , elem_grp_name_to_id_ex &
    , surf_grp_name_to_id_ex &
    , hecmw_streqr
  use fstr_setup_util
  use m_fstr_solid_setup
  use m_fstr_heat_setup
  use m_fstr_freq_setup
  use m_fstr_eigen_setup
  use m_fstr_dynamic_setup
  use m_fstr_param_setup
  use m_fstr_setup_finish
  use fstr_ctrl_common
  use fstr_ctrl_static
  use fstr_ctrl_heat
  use fstr_ctrl_eigen
  use fstr_ctrl_dynamic
  use fstr_ctrl_material
  use mContactParam
  use m_static_get_prop
  use m_out
  use m_step
  use m_utilities
  use fstr_ctrl_util_f
  implicit none

contains

  !=============================================================================!
  !> Read in and initialize control data                                        !
  !=============================================================================!
  subroutine fstr_setup( cntl_filename, hecMESH, hecMAT, fstrPARAM,  &
      fstrSOLID, fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL, fstrFREQ )
    character(len=HECMW_FILENAME_LEN) :: cntl_filename, input_filename
    type(hecmwST_local_mesh),target :: hecMESH
    type(hecmwST_matrix),target :: hecMAT
    type(fstr_param),target   :: fstrPARAM
    type(fstr_solid),target   :: fstrSOLID
    type(fstr_eigen),target     :: fstrEIG
    type(fstr_heat),target    :: fstrHEAT
    type(fstr_dynamic),target :: fstrDYNAMIC
    type(fstr_couple),target  :: fstrCPL
    type(fstr_freqanalysis), target :: fstrFREQ

    integer(kind=kint) :: ctrl, ctrl_list(20), ictrl
    type(fstr_param_pack) :: P

    integer, parameter :: MAXOUTFILE = 10
    double precision, parameter :: dpi = 3.14159265358979323846D0

    integer(kind=kint) :: version, result, visual, femap
    integer(kind=kint) :: rcode, n, cid, nout, ierror
    character(len=HECMW_NAME_LEN) :: header_name
    type(t_output_ctrl) :: outctrl

    ! counters
    integer(kind=kint) :: c_solution, c_solver, c_nlsolver, c_step, c_echo, c_amplitude
    integer(kind=kint) :: c_static, c_boundary, c_cload, c_dload, c_temperature, c_reftemp, c_spring, c_elemact
    integer(kind=kint) :: c_heat, c_fixtemp, c_cflux, c_dflux, c_sflux, c_film, c_sfilm, c_radiate, c_sradiate
    integer(kind=kint) :: c_eigen, c_contact, c_contactparam, c_embed, c_contact_if
    integer(kind=kint) :: c_dynamic, c_velocity, c_acceleration
    integer(kind=kint) :: c_fload, c_facc, c_eigenread
    integer(kind=kint) :: c_couple, c_material
    integer(kind=kint) :: c_mpc, c_weldline, c_initial
    integer(kind=kint) :: c_istep, c_localcoord, c_section
    integer(kind=kint) :: c_elemopt, c_aincparam, c_timepoints, c_cnvparam
    integer(kind=kint) :: c_output, islog
    integer(kind=kint) :: cache = 1

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

    c_solution = 0; c_solver   = 0; c_nlsolver = 0; c_step   = 0; c_output = 0; c_echo = 0; c_amplitude = 0
    c_static   = 0; c_boundary = 0; c_cload  = 0; c_dload = 0; c_temperature = 0; c_reftemp = 0; c_spring = 0;
    c_elemact = 0;
    c_heat     = 0; c_fixtemp  = 0; c_cflux  = 0; c_dflux = 0; c_sflux = 0
    c_film     = 0; c_sfilm    = 0; c_radiate= 0; c_sradiate = 0
    c_eigen    = 0; c_contact  = 0; c_contactparam = 0; c_embed  = 0; c_contact_if = 0
    c_dynamic  = 0; c_velocity = 0; c_acceleration = 0
    c_couple   = 0; c_material = 0; c_section =0
    c_mpc      = 0; c_weldline = 0; c_initial = 0
    c_istep    = 0; c_localcoord = 0
    c_fload    = 0; c_facc = 0; c_eigenread = 0
    c_elemopt  = 0;
    c_aincparam= 0; c_timepoints = 0; c_cnvparam = 0

    ctrl_list = 0
    ictrl = 1
    ctrl  = fstr_ctrl_open( cntl_filename )
    if( ctrl < 0 ) then
      write(*,*) '### Error: Cannot open FSTR control file : ', cntl_filename
      write(ILOG,*) '### Error: Cannot open FSTR control file : ', cntl_filename
      stop HECMW_EXIT_INPUT
    end if

    version =0
    do
      rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )
      if(     header_name == '!VERSION' ) then
        rcode = fstr_ctrl_get_data_ex( ctrl, 1, 'i ', version )
      else if(     header_name == '!SOLUTION' ) then
        c_solution = c_solution + 1
        call fstr_setup_SOLUTION( ctrl, c_solution, P )
      else if(     header_name == '!NONLINEAR_SOLVER' ) then
        c_nlsolver = c_nlsolver + 1
        call fstr_setup_NONLINEAR_SOLVER( ctrl, c_nlsolver, P )
      else if( header_name == '!SOLVER' ) then
        c_solver = c_solver + 1
        call fstr_setup_SOLVER( ctrl, c_solver, P )
      else if( header_name == '!ISTEP' ) then
        c_istep = c_istep + 1
      else if( header_name == '!STEP' ) then
        if( version==0 ) then
          c_step = c_step + 1
          call fstr_setup_STEP( ctrl, c_step, P )
        else
          c_istep = c_istep + 1
        endif
      else if( header_name == '!WRITE' ) then
        call fstr_ctrl_get_output( ctrl, outctrl, islog, result, visual, femap )
        if( visual==1 ) P%PARAM%fg_visual= 1
        if( result==1 ) P%PARAM%fg_result = 1
        c_output = c_output+1
      else if( header_name == '!ECHO' ) then
        c_echo = c_echo + 1
        call fstr_setup_ECHO( ctrl, c_echo, P )
      else if( header_name == '!RESTART' ) then
        call fstr_setup_RESTART( ctrl, nout, P%PARAM%restart_version )
        fstrSOLID%restart_nout= nout
        fstrDYNAMIC%restart_nout= nout
        fstrHEAT%restart_nout= nout
      else if( header_name == '!ORIENTATION' ) then
        c_localcoord = c_localcoord + 1
      else if( header_name == '!AUTOINC_PARAM' ) then
        c_aincparam = c_aincparam + 1
      else if( header_name == '!CONVERG_PARAM' ) then
        c_cnvparam = c_cnvparam + 1
      else if( header_name == '!TIME_POINTS' ) then
        c_timepoints = c_timepoints + 1
      else if( header_name == '!OUTPUT_SSTYPE' ) then
        call fstr_setup_OUTPUT_SSTYPE( ctrl, P )
      else if( header_name == '!INITIAL_CONDITION' ) then
        c_initial = c_initial + 1
      else if( header_name == '!AMPLITUDE' ) then
        c_amplitude = c_amplitude + 1
        call fstr_setup_AMPLITUDE( ctrl, P )
      else if( header_name == '!ELEMENT_ACTIVATION' ) then
        c_elemact = c_elemact + 1
        call fstr_setup_ELEMENT_ACTIVATION( ctrl, c_elemact, P )

        !--------------- for static -------------------------

      else if( header_name == '!STATIC' ) then
        c_static = c_static + 1
        call fstr_setup_STATIC( ctrl, c_static, P )
      else if( header_name == '!BOUNDARY' ) then
        c_boundary = c_boundary + 1
        call fstr_setup_BOUNDARY( ctrl, c_boundary, P )
      else if( header_name == '!CLOAD' ) then
        c_cload = c_cload + 1
        call fstr_setup_CLOAD( ctrl, c_cload, P )
        n = fstr_ctrl_get_data_line_n( ctrl )
      else if( header_name == '!DLOAD' ) then
        c_dload = c_dload + 1
        call fstr_setup_DLOAD( ctrl, c_dload, P )
      else if( header_name == '!CONTACT_ALGO' ) then
        call fstr_setup_CONTACTALGO( ctrl, P )
      else if( header_name == '!CONTACT' ) then
        n = fstr_ctrl_get_data_line_n( ctrl )
        c_contact = c_contact + n
      else if( header_name == '!EMBED' ) then
        n = fstr_ctrl_get_data_line_n( ctrl )
        c_embed = c_embed + n
      else if( header_name == '!CONTACT_PARAM' ) then
        c_contactparam = c_contactparam + 1
      else if( header_name == '!CONTACT_INTERFERENCE' ) then
        n = fstr_ctrl_get_data_line_n( ctrl )
        c_contact_if = c_contact_if + n
      else if( header_name == '!MATERIAL' ) then
        c_material = c_material + 1
      else if( header_name == '!TEMPERATURE' ) then
        c_temperature = c_temperature + 1
        call fstr_setup_TEMPERATURE( ctrl, c_temperature, P )
      else if( header_name == '!SPRING' ) then
        c_spring = c_spring + 1
        call fstr_setup_SPRING( ctrl, c_spring, P )
      else if( header_name == '!REFTEMP' ) then
        c_reftemp = c_reftemp + 1
        call fstr_setup_REFTEMP( ctrl, c_reftemp, P )

        !--------------- for heat -------------------------

      else if( header_name == '!HEAT' ) then
        c_heat = c_heat + 1
      else if( header_name == '!FIXTEMP' ) then
        c_fixtemp = c_fixtemp + 1
        call fstr_setup_FIXTEMP( ctrl, c_fixtemp, P )
      else if( header_name == '!CFLUX' ) then
        c_cflux = c_cflux + 1
        call fstr_setup_CFLUX( ctrl, c_cflux, P )
      else if( header_name == '!DFLUX' ) then
        c_dflux = c_dflux + 1
        call fstr_setup_DFLUX( ctrl, c_dflux, P )
      else if( header_name == '!SFLUX' ) then
        c_sflux = c_sflux + 1
        call fstr_setup_SFLUX( ctrl, c_sflux, P )
      else if( header_name == '!FILM' ) then
        c_film = c_film + 1
        call fstr_setup_FILM( ctrl, c_film, P )
      else if( header_name == '!SFILM' ) then
        c_sfilm = c_sfilm + 1
        call fstr_setup_SFILM( ctrl, c_sfilm, P )
      else if( header_name == '!RADIATE' ) then
        c_radiate = c_radiate + 1
        call fstr_setup_RADIATE( ctrl, c_radiate, P )
      else if( header_name == '!SRADIATE' ) then
        c_sradiate = c_sradiate + 1
        call fstr_setup_SRADIATE( ctrl, c_sradiate, P )
      else if( header_name == '!WELD_LINE' ) then
        c_weldline = c_weldline + 1

        !--------------- for eigen -------------------------

      else if( header_name == '!EIGEN' ) then
        c_eigen = c_eigen + 1
        call fstr_setup_EIGEN( ctrl, c_eigen, P )

        !--------------- for dynamic -------------------------

      else if( header_name == '!DYNAMIC' ) then
        c_dynamic = c_dynamic + 1
        call fstr_setup_DYNAMIC( ctrl, c_eigen, P )
      else if( header_name == '!VELOCITY' ) then
        c_velocity = c_velocity + 1
        call fstr_setup_VELOCITY( ctrl, c_eigen, P )
      else if( header_name == '!ACCELERATION' ) then
        c_acceleration = c_acceleration + 1
        call fstr_setup_ACCELERATION( ctrl, c_eigen, P )
      else if( header_name == '!FLOAD' ) then
        c_fload = c_fload + 1
        call fstr_setup_FLOAD( ctrl , c_fload, P )
      else if( header_name == '!FACCELERATION' ) then
        ! acceleration input of the frequency response analysis
        c_facc = c_facc + 1
        call fstr_setup_FACCELERATION( ctrl, c_facc, P )
      else if( header_name == '!EIGENREAD' ) then
        c_eigenread = c_eigenread + 1
        call fstr_setup_eigenread( ctrl, c_eigenread, P )

        !--------------- for couple -------------------------

      else if( header_name == '!COUPLE' ) then
        c_couple = c_couple + 1
        call fstr_setup_COUPLE( ctrl, c_couple, P )

        !--------------- for mpc -------------------------

      else if( header_name == '!MPC' ) then
        c_mpc = c_mpc + 1
        call fstr_setup_MPC( ctrl, c_mpc, P )

        !--------------------- for input -------------------------

      else if( header_name == '!INCLUDE' ) then
        ctrl_list(ictrl) = ctrl
        input_filename   = ""
        ierror = fstr_ctrl_get_param_ex( ctrl, 'INPUT ', '# ', 0, 'F', input_filename )
        ctrl   = fstr_ctrl_open( input_filename )
        if( ctrl < 0 ) then
          write(*,*) '### Error: Cannot open FSTR control file : ', input_filename
          write(ILOG,*) '### Error: Cannot open FSTR control file : ', input_filename
          stop HECMW_EXIT_INPUT
        end if
        ictrl = ictrl + 1
        cycle

        !--------------------- END -------------------------

      else if( header_name == '!END' ) then
        exit
      end if

      ! next
      if( fstr_ctrl_seek_next_header(ctrl) == 0 )then
        if( ictrl == 1 )then
          exit
        else
          ierror= fstr_ctrl_close( ctrl )
          ictrl = ictrl - 1
          ctrl  = ctrl_list(ictrl)
          if( fstr_ctrl_seek_next_header(ctrl) == 0 ) exit
        endif
      endif
    end do

    ! -----
    call fstr_setup_alloc( P, c_contact, c_embed, c_weldline, c_initial, c_istep, c_localcoord, &
      c_aincparam, c_cnvparam, c_timepoints, c_contactparam, c_contact_if, c_material )

    ! -----
    rcode = fstr_ctrl_rewind( ctrl )

    c_istep    = 0
    c_heat     = 0
    c_material = 0
    cid = 0
    c_output = 0
    c_contact  = 0
    c_contactparam  = 0
    c_contact_if  = 0
    c_embed = 0
    c_initial = 0
    c_localcoord = 0
    c_section = 0
    fstrHEAT%WL_tot = 0
    c_elemopt = 0
    c_aincparam = 0
    c_timepoints = 0
    c_cnvparam = 0
    ictrl = 1
    do
      rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )

      if( header_name == '!ORIENTATION' ) then
        c_localcoord = c_localcoord + 1
        if( fstr_setup_ORIENTATION( ctrl, hecMESH, g_LocalCoordSys(c_localcoord) )/=0 ) then
          write(*,*) '### Error: Fail in read in ORIENTATION definition : ', c_localcoord
          write(ILOG,*) '### Error: Fail in read in ORIENTATION definition : ', c_localcoord
          stop HECMW_EXIT_INPUT
        endif

        ! ----- CONTACT condition setting
      elseif( header_name == '!CONTACT' ) then
        call fstr_setup_CONTACT( ctrl, c_contact, P )

        ! ----- EMBED condition setting
      elseif( header_name == '!EMBED' ) then
        call fstr_setup_EMBED( ctrl, c_embed, P )

      else if( header_name == '!ISTEP'  ) then
        c_istep = c_istep+1
        call fstr_setup_ISTEP( ctrl, c_istep, P )
      else if( header_name == '!STEP' .and. version>=1 ) then
        c_istep = c_istep+1
        call fstr_setup_STEP_store( ctrl, c_istep, P )

      else if( header_name == '!HEAT'  ) then
        c_heat = c_heat + 1
        call fstr_setup_HEAT( ctrl, c_heat, P )

      else if( header_name == '!WELD_LINE'  ) then
        call fstr_setup_WELD_LINE( ctrl, P )

      else if( header_name == '!INITIAL_CONDITION' .or. header_name == '!INITIAL CONDITION' ) then
        c_initial = c_initial+1
        if( fstr_setup_INITIAL( ctrl, g_InitialCnd(c_initial), P%MESH )/=0 ) then
           write(*,*) '### Error: Fail in read in INITIAL CONDITION definition : ' ,c_initial
           write(ILOG,*) '### Error: Fail in read in INITIAL CONDITION definition : ', c_initial
           stop HECMW_EXIT_INPUT
        endif

      else if( header_name == '!SECTION'  ) then
        c_section = c_section+1
        call fstr_setup_SECTION( ctrl, c_section, P )

      else if( header_name == '!ELEMOPT'  ) then
        c_elemopt = c_elemopt+1
        call fstr_setup_ELEMOPT( ctrl, c_elemopt, P )

        !== following material properties ==
      else if( header_name == '!MATERIAL' ) then
        c_material = c_material+1
        call fstr_setup_MATERIAL( ctrl, c_material, cache, cid, P )

      else if( header_name == '!ELASTIC' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_ELASTIC( ctrl, cid, P )
      else if( header_name == '!PLASTIC' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_PLASTIC( ctrl, cid, P )
      else if( header_name == '!HYPERELASTIC' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_HYPERELASTIC( ctrl, cid, P )
      else if( header_name == '!VISCOELASTIC' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_VISCOELASTIC( ctrl, cid, P )
      else if( header_name == '!TRS' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_TRS( ctrl, cid, P )
      else if( header_name == '!CREEP' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_CREEP( ctrl, cid, P )
      else if( header_name == '!DENSITY' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_DENSITY( ctrl, cid, P )
      else if( header_name == '!EXPANSION_COEF' .or. header_name == '!EXPANSION_COEFF' .or. &
          header_name == '!EXPANSION') then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_EXPANSION_COEFF( ctrl, cid, P )
      else if( header_name == '!DAMPING') then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_DAMPING( ctrl, cid, P )
      else if( header_name == '!FLUID' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_FLUID( ctrl, cid, P )
      else if( header_name == '!SPRING_D' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_SPRING_D( ctrl, cid, P )
      else if( header_name == '!SPRING_A' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_SPRING_A( ctrl, cid, P )
      else if( header_name == '!DASHPOT_D' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_DASHPOT_D( ctrl, cid, P )
      else if( header_name == '!DASHPOT_A' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_DASHPOT_A( ctrl, cid, P )
      else if( header_name == '!USER_MATERIAL' ) then
        if( cid == 0 ) call material_card_err_stop( header_name )
        call fstr_setup_USER_MATERIAL( ctrl, cid, P )


        ! == Following output control ==
      else if( header_name == '!WRITE' ) then
        call fstr_setup_WRITE_store( ctrl, P )

      else if( header_name == '!OUTPUT_RES' ) then
        call fstr_setup_OUTPUT_RES( ctrl, P )
      else if( header_name == '!OUTPUT_VIS' ) then
        call fstr_setup_OUTPUT_VIS( ctrl, P )
      else if( header_name == '!AUTOINC_PARAM' ) then
        c_aincparam = c_aincparam + 1
        call fstr_setup_AUTOINC_PARAM( ctrl, c_aincparam, P )
      else if( header_name == '!CONVERG_PARAM' ) then
        c_cnvparam = c_cnvparam + 1
        call fstr_setup_CONVERG_PARAM( ctrl, c_cnvparam, P )
      else if( header_name == '!TIME_POINTS'  ) then
        c_timepoints = c_timepoints + 1
        call fstr_setup_TIME_POINTS( ctrl, c_timepoints, P )
      else if( header_name == '!CONTACT_PARAM' ) then
        c_contactparam = c_contactparam + 1
        call fstr_setup_CONTACT_PARAM( ctrl, c_contactparam, P )
      else if( header_name == '!CONTACT_INTERFERENCE' ) then
        call fstr_setup_CONTACT_INTERFERENCE( ctrl, c_contact_if, P )
      else if( header_name == '!ULOAD' ) then
        call fstr_setup_ULOAD( ctrl )

      else if( header_name == '!INCLUDE' ) then
        ctrl_list(ictrl) = ctrl
        input_filename   = ""
        ierror = fstr_ctrl_get_param_ex( ctrl, 'INPUT ', '# ', 0, 'F', input_filename )
        ctrl   = fstr_ctrl_open( input_filename )
        if( ctrl < 0 ) then
          write(*,*) '### Error: Cannot open FSTR control file : ', input_filename
          write(ILOG,*) '### Error: Cannot open FSTR control file : ', input_filename
          stop HECMW_EXIT_INPUT
        end if
        ictrl = ictrl + 1
        cycle

      else if( header_name == '!END' ) then
        exit
      endif

      ! next
      if( fstr_ctrl_seek_next_header(ctrl) == 0 )then
        if( ictrl == 1 )then
          exit
        else
          ierror= fstr_ctrl_close( ctrl )
          ictrl = ictrl - 1
          ctrl  = ctrl_list(ictrl)
          if( fstr_ctrl_seek_next_header(ctrl) == 0 ) exit
        endif
      endif

    end do

    call fstr_setup_finish( P, c_material )
    if( P%PARAM%fg_visual == kON ) call fstr_setup_visualize( ctrl, P%MESH )
    rcode = fstr_ctrl_close( ctrl )

  end subroutine fstr_setup

  subroutine fstr_couple_finalize( C )
    implicit none
    type( fstr_couple ) :: C
    if(associated(C%coupled_node)) deallocate(C%coupled_node)
    if(associated(C%trac)) deallocate(C%trac)
    if(associated(C%disp)) deallocate(C%disp)
    if(associated(C%velo)) deallocate(C%velo)
    if(associated(C%accel)) deallocate(C%accel)
    if(associated(C%index)) deallocate(C%index)
  end subroutine fstr_couple_finalize

  subroutine material_card_err_stop( header_name )
    implicit none
    character(len=*) :: header_name

    call fstr_setup_util_err_stop( '### Error: '//trim(header_name)// &
      ' : given before the first !MATERIAL; put it after !MATERIAL, NAME=<material name>' )
  end subroutine material_card_err_stop

end module m_fstr_setup
