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
    integer(kind=kint) :: rcode, n, i, cid, nout, ierror
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

    ! !DENSITY, !PLASTIC and the others checked by cid before the first !MATERIAL go to the material of the last section
    cid = 0
    do i=1,hecMESH%section%n_sect
      if( hecMESH%section%sect_type(i) /= 4 ) cid = hecMESH%section%sect_mat_ID_item(i)
    enddo

    ! -----
    rcode = fstr_ctrl_rewind( ctrl )

    c_istep    = 0
    c_heat     = 0
    c_material = 0
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
        if( fstr_setup_ORIENTATION( ctrl, hecMESH, c_localcoord, g_LocalCoordSys(c_localcoord) )/=0 ) then
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
        if( c_material >0 ) then
          call fstr_setup_ELASTIC( ctrl, cid, P )
        endif
      else if( header_name == '!PLASTIC' ) then
        if( cid >0 ) then
          call fstr_setup_PLASTIC( ctrl, cid, P )
        endif
      else if( header_name == '!HYPERELASTIC' ) then
        if( cid >0 ) then
          call fstr_setup_HYPERELASTIC( ctrl, cid, P )
        endif
      else if( header_name == '!VISCOELASTIC' ) then
        if( cid >0 ) then
          call fstr_setup_VISCOELASTIC( ctrl, cid, P )
        endif
      else if( header_name == '!TRS' ) then
        if( cid >0 ) then
          call fstr_setup_TRS( ctrl, cid, P )
        endif
      else if( header_name == '!CREEP' ) then
        if( cid >0 ) then
          call fstr_setup_CREEP( ctrl, cid, P )
        endif
      else if( header_name == '!DENSITY' ) then
        if( cid >0 ) then
          call fstr_setup_DENSITY( ctrl, cid, P )
        endif
      else if( header_name == '!EXPANSION_COEF' .or. header_name == '!EXPANSION_COEFF' .or. &
          header_name == '!EXPANSION') then
        if( cid >0 ) then
          call fstr_setup_EXPANSION_COEFF( ctrl, cid, P )
        endif
      else if( header_name == '!DAMPING') then
        if( cid >0 ) then
          call fstr_setup_DAMPING( ctrl, cid, P )
        endif
      else if( header_name == '!FLUID' ) then
        if( c_material >0 ) then
          call fstr_setup_FLUID( ctrl, cid, P )
        endif
      else if( header_name == '!SPRING_D' ) then
        if( c_material >0 ) then
          call fstr_setup_SPRING_D( ctrl, cid, P )
        endif
      else if( header_name == '!SPRING_A' ) then
        if( c_material >0 ) then
          call fstr_setup_SPRING_A( ctrl, cid, P )
        endif
      else if( header_name == '!DASHPOT_D' ) then
        if( c_material >0 ) then
          call fstr_setup_DASHPOT_D( ctrl, cid, P )
        endif
      else if( header_name == '!DASHPOT_A' ) then
        if( c_material >0 ) then
          call fstr_setup_DASHPOT_A( ctrl, cid, P )
        endif
      else if( header_name == '!USER_MATERIAL' ) then
        if( cid >0 ) then
          call fstr_setup_USER_MATERIAL( ctrl, cid, P )
        endif


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

    call fstr_setup_finish( ctrl, P, c_material )
    rcode = fstr_ctrl_close( ctrl )

  end subroutine fstr_setup

  !> Allocation after the cards are counted: the arrays sized by the counts, the materials
  !> and sections from the mesh, and the output controls. The section defaults depend on
  !> the analysis type, and fstr_convert_rotation_dof rewrites the boundary conditions,
  !> so this is called after the analysis type and all the boundary conditions are read.
  subroutine fstr_setup_alloc( P, c_contact, c_embed, c_weldline, c_initial, c_istep, c_localcoord, &
      c_aincparam, c_cnvparam, c_timepoints, c_contactparam, c_contact_if, c_material )
    use mMaterial
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: c_contact, c_embed, c_weldline, c_initial, c_istep, c_localcoord
    integer(kind=kint), intent(in) :: c_aincparam, c_cnvparam, c_timepoints, c_contactparam, c_contact_if
    integer(kind=kint), intent(in) :: c_material

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    type(fstr_heat), pointer :: fstrHEAT
    type(fstr_dynamic), pointer :: fstrDYNAMIC
    integer(kind=kint) :: n, i, cid, n_totlyr
    real(kind=kreal) :: ee, pp, rho, alpha, thick, alpha_over_mu
    real(kind=kreal) :: beam_radius,                          &
      beam_angle1, beam_angle2, beam_angle3,&
      beam_angle4, beam_angle5, beam_angle6
    type(tshellmat),pointer :: shmat(:)
    character(len=HECMW_FILENAME_LEN) :: logfileNAME

    hecMESH => P%MESH
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID
    fstrHEAT => P%HEAT
    fstrDYNAMIC => P%DYN
    write( logfileNAME, '(i5,''.log'')' ) myrank

    fstrSOLID%n_contacts = c_contact
    if( c_contact>0 ) then
      allocate( fstrSOLID%contacts( c_contact ) )
      ! convert SURF_SURF contact to NODE_SURF contact
      call fstr_convert_contact_type( P%MESH )
    endif
    call fstr_convert_rotation_dof( P )
    fstrSOLID%n_embeds = c_embed
    if( c_embed>0 )  allocate( fstrSOLID%embeds( c_embed ) )
    if( c_weldline>0 ) allocate( fstrHEAT%weldline( c_weldline ) )
    if( c_initial>0 ) allocate( g_InitialCnd( c_initial ) )
    if( c_istep>0 ) then
      allocate( fstrSOLID%step_ctrl( c_istep ) )
      do i=1, c_istep
        call init_stepInfo( fstrSOLID%step_ctrl(i) )
        if( p%PARAM%solution_type==kstDYNAMIC ) then
          fstrSOLID%step_ctrl(i)%num_substep = fstrDYNAMIC%n_step
          fstrSOLID%step_ctrl(i)%initdt      = fstrDYNAMIC%t_delta
          fstrSOLID%step_ctrl(i)%elapsetime  = dble(fstrDYNAMIC%n_step) * fstrDYNAMIC%t_delta
          fstrSOLID%step_ctrl(i)%mindt       = fstrDYNAMIC%t_delta
          fstrSOLID%step_ctrl(i)%maxdt       = fstrDYNAMIC%t_delta
        endif
      end do
    endif
    if( c_localcoord>0 ) allocate( g_LocalCoordSys(c_localcoord) )
    allocate( fstrPARAM%ainc(0:c_aincparam) )
    do i=0,c_aincparam
      call init_AincParam( fstrPARAM%ainc(i) )
    end do
    allocate( fstrPARAM%cnvparam(0:c_cnvparam) )
    do i=0,c_cnvparam
      call init_ConvergParam( fstrPARAM%cnvparam(i) )
    end do
    if( c_timepoints>0 ) allocate( fstrPARAM%timepoints(c_timepoints) )
    allocate( fstrPARAM%contactparam(0:c_contactparam) )
    do i=0,c_contactparam
      call init_ContactParam( fstrPARAM%contactparam(i) )
    end do
    if( c_contact_if>0 )then
      allocate( fstrPARAM%contact_if( c_contact_if ) )
      do i=1,c_contact_if
        call init_Contact_IF( fstrPARAM%contact_if(i) )
      end do
    end if

    P%SOLID%is_33shell = 0
    P%SOLID%is_33beam  = 0

    do i=1,hecMESH%n_elem_type
      n =  hecMESH%elem_type_item(i)
      if (n == 781 .or. n == 761)then
        P%SOLID%is_33shell = 1
      elseif (n == 641)then
        P%SOLID%is_33beam  = 1
      endif
    enddo

    n = c_material
    if( hecMESH%material%n_mat>n ) n= hecMESH%material%n_mat
    if( n==0 ) stop "material property not defined!"
    allocate( fstrSOLID%materials( n ) )
    do i = 1, n
      call initMaterial(fstrSOLID%materials(i))
    enddo

    ! for section control
    allocate( fstrSOLID%sections(hecMESH%section%n_sect) )
    do i=1,hecMESH%section%n_sect
      ! set default 361 element formulation
      if( p%PARAM%solution_type==kstSTATIC .or. p%PARAM%solution_type==kstDYNAMIC ) then
        if( p%PARAM%nlgeom ) then
          fstrSOLID%sections(i)%elemopt361 = kel361FBAR
        else
          fstrSOLID%sections(i)%elemopt361 = kel361IC
        end if
      else if( p%PARAM%solution_type==kstEIGEN ) then
        fstrSOLID%sections(i)%elemopt361 = kel361IC
      else if( p%PARAM%solution_type==kstSTATICEIGEN ) then
        fstrSOLID%sections(i)%elemopt361 = kel361FBAR
      else
        fstrSOLID%sections(i)%elemopt361 = kel361FI
      end if
      fstrSOLID%sections(i)%elemopt341 = kel341FI
      fstrSOLID%sections(i)%elemopt611 = kel611EULER
      ! INTERFACE sections are skipped in the loop below and keep this thickness
      fstrSOLID%sections(i)%thickness = 1.0d0
    enddo

    if( hecMESH%section%n_sect >0 ) then
      do i=1,hecMESH%section%n_sect
        if( hecMESH%section%sect_type(i) == 4 ) cycle
        cid = hecMESH%section%sect_mat_ID_item(i)
        if( cid>n ) stop "Error in material property definition!"
        if( fstrPARAM%nlgeom .or. fstrPARAM%solution_type==kstSTATICEIGEN ) &
          fstrSOLID%materials(cid)%nlgeom_flag = 1
        nullify(shmat)
        call fstr_get_prop(hecMESH,shmat,i,ee,pp,rho,alpha,thick,&
          n_totlyr,alpha_over_mu, &
          beam_radius,beam_angle1,beam_angle2,beam_angle3, &
          beam_angle4,beam_angle5,beam_angle6)
        fstrSOLID%sections(i)%thickness = thick
        fstrSOLID%materials(cid)%name = hecMESH%material%mat_name(cid)
        fstrSOLID%materials(cid)%variables(M_YOUNGS)=ee
        fstrSOLID%materials(cid)%variables(M_POISSON)=pp
        fstrSOLID%materials(cid)%variables(M_DENSITY)=rho
        fstrSOLID%materials(cid)%variables(M_EXPANSION)=alpha
        fstrSOLID%materials(cid)%variables(M_ALPHA_OVER_MU)= alpha_over_mu
        fstrSOLID%materials(cid)%variables(M_BEAM_RADIUS)=beam_radius
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE1)=beam_angle1
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE2)=beam_angle2
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE3)=beam_angle3
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE4)=beam_angle4
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE5)=beam_angle5
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE6)=beam_angle6
        fstrSOLID%materials(cid)%mtype = ELASTIC
        if( hecMESH%section%sect_type(i) == 2 ) then
          fstrSOLID%materials(cid)%totallyr = n_totlyr
          fstrSOLID%materials(cid)%shell_var => shmat
        endif
      enddo
    endif

    allocate( fstrSOLID%output_ctrl( 4 ) )
    call fstr_init_outctrl(fstrSOLID%output_ctrl(1))
    fstrSOLID%output_ctrl( 1 )%filename = trim(logfileNAME)
    fstrSOLID%output_ctrl( 1 )%filenum = ILOG
    call fstr_init_outctrl(fstrSOLID%output_ctrl(2))
    call fstr_init_outctrl(fstrSOLID%output_ctrl(3))
    call fstr_init_outctrl(fstrSOLID%output_ctrl(4))
  end subroutine fstr_setup_alloc

  !> Processing after all the cards are read: the temperature arrays, the default step,
  !> the element initialization and the checks, and the result arrays (fstr_setup_post).
  !> ctrl is used only for !VISUAL.
  subroutine fstr_setup_finish( ctrl, P, c_material )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: c_material

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    type(fstr_dynamic), pointer :: fstrDYNAMIC
    integer(kind=kint) :: n, i, cid, ierror, n_totlyr

    hecMESH => P%MESH
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID
    fstrDYNAMIC => P%DYN

    ! ----- material type judgement. in case of infinitive analysis, nlgeom_flag=0
    if( .not. P%PARAM%nlgeom ) then
      do i=1, c_material
        fstrSOLID%materials(i)%nlgeom_flag = 0
      enddo
    endif

    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      allocate ( fstrSOLID%temperature( hecMESH%n_node )      ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TEMPERATURE>'
        write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrSOLID%temperature = fstrPARAM%ref_temp
      allocate ( fstrSOLID%last_temp( hecMESH%n_node )      ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, LAST_TEMP>'
        write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrSOLID%last_temp = 0.d0
      allocate ( fstrSOLID%temp_bak( hecMESH%n_node )      ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TEMP_BAK>'
        write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrSOLID%temp_bak = 0.d0
    endif

    if( associated(fstrSOLID%step_ctrl) )  then
      fstrSOLID%nstep_tot = size(fstrSOLID%step_ctrl)
      call setup_stepInfo_starttime( fstrSOLID%step_ctrl )
    else
      if( p%PARAM%solution_type==kstSTATIC .and. P%PARAM%nlgeom ) then
        write( *,* ) " ERROR: STEP not defined!"
        write( idbg,* ) "ERROR: STEP not defined!"
        call flush(idbg)
        call fstr_abort( HECMW_EXIT_MODEL )
      endif

      if( myrank==0 ) write(*,*)"Step control not defined! Using default step=1"
      fstrSOLID%nstep_tot = 1
      allocate( fstrSOLID%step_ctrl(1) )
      call init_stepInfo( fstrSOLID%step_ctrl(1) )
      if( p%PARAM%solution_type==kstDYNAMIC ) then
        fstrSOLID%step_ctrl(1)%num_substep = fstrDYNAMIC%n_step
        fstrSOLID%step_ctrl(1)%initdt      = fstrDYNAMIC%t_delta
        fstrSOLID%step_ctrl(1)%elapsetime  = dble(fstrDYNAMIC%n_step) * fstrDYNAMIC%t_delta
        fstrSOLID%step_ctrl(1)%mindt       = fstrDYNAMIC%t_delta
        fstrSOLID%step_ctrl(1)%maxdt       = fstrDYNAMIC%t_delta
      endif
      n =  fstrSOLID%BOUNDARY_ngrp_tot
      if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%Boundary(n) )
      do i = 1, n
        fstrSOLID%step_ctrl(1)%Boundary(i) = fstrSOLID%BOUNDARY_ngrp_GRPID(i)
      enddo
      n = fstrSOLID%CLOAD_ngrp_tot + fstrSOLID%DLOAD_ngrp_tot + fstrSOLID%TEMP_ngrp_tot + fstrSOLID%SPRING_ngrp_tot
      if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%Load(n) )
      n = 0
      do i = 1, fstrSOLID%CLOAD_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%CLOAD_ngrp_GRPID(i)
      enddo
      do i = 1, fstrSOLID%DLOAD_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%DLOAD_ngrp_GRPID(i)
      enddo
      do i = 1, fstrSOLID%TEMP_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%TEMP_ngrp_GRPID(i)
      enddo
      do i = 1, fstrSOLID%SPRING_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%SPRING_ngrp_GRPID(i)
      enddo
      n = fstrSOLID%elemact%ELEMACT_egrp_tot
      if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%ElemActivation(n) )
      do i = 1, n
        fstrSOLID%step_ctrl(1)%ElemActivation(i) = fstrSOLID%elemact%ELEMACT_egrp_GRPID(i)
      enddo
    endif

    call fstr_element_init( hecMESH, fstrSOLID, p%PARAM%solution_type )
    call fstr_check_shell_j2_inputs(hecMESH, fstrSOLID, fstrPARAM)
    call fstr_check_up_inputs(hecMESH, fstrSOLID)
    call fstr_check_shell_gauss_output(hecMESH, fstrSOLID)
    if( p%PARAM%solution_type==kstSTATIC .or. p%PARAM%solution_type==kstDYNAMIC .or.   &
      p%PARAM%solution_type==kstEIGEN  .or. p%PARAM%solution_type==kstSTATICEIGEN )  &
      call fstr_solid_alloc( hecMESH, fstrSOLID )

    if( p%PARAM%solution_type == kstHEAT) then
      p%PARAM%fg_irres = fstrSOLID%output_ctrl(3)%frequency
      p%PARAM%fg_iwres = fstrSOLID%output_ctrl(4)%frequency
      p%HEAT%elemact = p%SOLID%elemact
    endif

    n_totlyr = 1
    do i=1,hecMESH%section%n_sect
      cid = hecMESH%section%sect_mat_ID_item(i)
      n  = fstrSOLID%materials(cid)%totallyr
      if (n > n_totlyr)then
        n_totlyr = n
      endif
    enddo
    P%SOLID%max_lyr = n_totlyr

    call fstr_setup_post( ctrl, P )
  end subroutine fstr_setup_finish

  subroutine fstr_check_up_inputs(hecMESH, fstrSOLID)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint) :: i, isect, mtype

    do i = 1, hecMESH%n_elem
      if( hecMESH%elem_type(i) /= 361 ) cycle
      isect = hecMESH%section_ID(i)
      if( fstrSOLID%sections(isect)%elemopt361 /= kel361UP ) cycle
      if( .not.associated(fstrSOLID%elements(i)%gausses) ) cycle
      mtype = fstrSOLID%elements(i)%gausses(1)%pMaterial%mtype
      if( .not.isElastoplastic(mtype) .and. .not.isCreep(mtype) ) cycle
      write(*,*) '###ERROR### : Elastoplastic and creep materials are not supported by FORM361=UP'
      write(ILOG,*) '###ERROR### : Elastoplastic and creep materials are not supported by FORM361=UP'
      call hecmw_abort(hecmw_comm_get_comm())
    enddo
  end subroutine fstr_check_up_inputs

  subroutine fstr_check_shell_gauss_output(hecMESH, fstrSOLID)
    use mMechGauss, only: fstr_shell_output_point_count
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint) :: i

    if( hecMESH%n_dof /= 6 ) return
    if( .not.any(fstrSOLID%output_ctrl(3)%outinfo%on(9:10)) ) return
    do i=1,hecMESH%n_elem
      if( fstr_shell_output_point_count(fstrSOLID%elements(i))>0 ) cycle
      write(*,*) '###ERROR### : ISTRAIN/ISTRESS require shell thickness-point output support; element type = ', &
        hecMESH%elem_type(i)
      call hecmw_abort(hecmw_comm_get_comm())
    enddo
  end subroutine fstr_check_shell_gauss_output

  subroutine fstr_check_shell_j2_inputs(hecMESH, fstrSOLID, fstrPARAM)
    use m_ElastoPlastic, only: checkPlaneStressJ2Material, PlaneStressJ2ErrorMessage
    use elementInfo, only: NumOfShellThicknessQuadPoints
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    type(fstr_param), intent(in) :: fstrPARAM
    type(tMaterial), pointer :: material
    integer(kind=kint) :: i, nn, ierr, isect, index, ilayer
    real(kind=kreal) :: elastic(2)
    logical :: missing
    character(len=160) :: message

    do i = 1, hecMESH%n_elem
      if( hecMESH%elem_type(i) /= 741 ) cycle
      if( .not.associated(fstrSOLID%elements(i)%gausses) ) cycle
      material => fstrSOLID%elements(i)%gausses(1)%pMaterial
      if( .not.isElastoplastic(material%mtype) ) cycle
      nn = hecMESH%elem_node_index(i)-hecMESH%elem_node_index(i-1)
      message = ''
      call checkPlaneStressJ2Material(material, ierr)
      if( ierr /= 0 ) then
        message = PlaneStressJ2ErrorMessage(ierr)
      else if( .not.fstr_uses_plane_stress_j2_shell(741, nn, material) ) then
        message = 'Elastoplastic MITC4 shells require INFINITESIMAL material kinematics'
      else if( fstrPARAM%solution_type /= kstSTATIC .or. .not.fstrPARAM%nlgeom ) then
        message = 'Elastoplastic MITC4 shells require !SOLUTION, TYPE=STATIC, NONLINEAR'
      else if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
        message = 'Temperature loading is not supported in models with elastoplastic MITC4 shells'
      endif
      if( len_trim(message) == 0 .and. associated(material%shell_var) ) then
        ! The J2 return mapping uses the common elastic table, not layer-specific constants.
        call fetch_TableData(MC_ISOELASTIC, material%dict, elastic, missing)
        do ilayer = 1, material%totallyr
          if( .not.ieee_is_finite(material%shell_var(ilayer)%ee) .or. &
              .not.ieee_is_finite(material%shell_var(ilayer)%pp) ) then
            message = 'Invalid elastic properties in J2 shell layer'
            exit
          endif
          if( abs(material%shell_var(ilayer)%ee-elastic(1)) > 1.0d-12*elastic(1) .or. &
              abs(material%shell_var(ilayer)%pp-elastic(2)) > 1.0d-12 ) then
            message = 'J2 shell layers must share the !ELASTIC Young modulus and Poisson ratio'
            exit
          endif
        enddo
      endif
      if( len_trim(message) == 0 ) then
        isect = hecMESH%section_ID(i)
        index = hecMESH%section%sect_I_index(isect-1)
        if( hecMESH%section%sect_I_index(isect)-index /= 1 ) then
          message = 'Elastoplastic MITC4 shells require a SHELL section with 2 thickness integration points'
        else if( hecMESH%section%sect_I_item(index+1) /= NumOfShellThicknessQuadPoints(741) ) then
          message = 'Elastoplastic MITC4 shells require 2 thickness integration points in !SECTION'
        endif
      endif
      if( len_trim(message) == 0 ) cycle
      write(*,*) '###ERROR### : ', trim(message)
      write(ILOG,*) '###ERROR### : ', trim(message)
      call hecmw_abort(hecmw_comm_get_comm())
    enddo
  end subroutine fstr_check_shell_j2_inputs


  !-----------------------------------------------------------------------------!
  !> Initial setting of postprecessor

  subroutine fstr_setup_post_phys_alloc(phys, NDOF, n_node, n_elem)
    implicit none
    type(fstr_solid_physic_val), pointer :: phys
    integer(kind=kint) :: NDOF, n_node, n_elem, mdof, istat
    mdof = (NDOF*NDOF+NDOF)/2;
    allocate ( phys%STRAIN  (mdof*n_node), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%STRAIN"
    allocate ( phys%STRESS  (mdof*n_node), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%STRESS"
    allocate ( phys%MISES   (     n_node), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%MISES"
    allocate ( phys%ESTRAIN (mdof*n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%ESTRAIN"
    allocate ( phys%ESTRESS (mdof*n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%ESTRESS"
    allocate ( phys%EMISES  (     n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%EMISES"
    allocate ( phys%EPLSTRAIN (   n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%EPLSTRAIN"
    allocate ( phys%ENQM    (12*n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%ENQM"
  end subroutine fstr_setup_post_phys_alloc

  subroutine fstr_setup_post( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl, i
    type(fstr_param_pack) :: P
    type(fstr_solid_physic_val), pointer :: phys => null()

    if( P%PARAM%solution_type == kstSTATIC &
        .or. P%PARAM%solution_type == kstEIGEN  &
        .or. P%PARAM%solution_type == kstDYNAMIC &
        .or. P%PARAM%solution_type == kstSTATICEIGEN ) then
      ! Memory Allocation for Result Vectors ------------
      if( P%MESH%n_dof == 6 .or. P%SOLID%is_33shell == 1 ) then
        allocate ( P%SOLID%SHELL )
        call fstr_setup_post_phys_alloc(P%SOLID%SHELL,3, P%MESH%n_node,P%MESH%n_elem)
        allocate ( P%SOLID%SHELL%LAYER(P%SOLID%max_lyr) )
        do i=1,P%SOLID%max_lyr
          allocate ( P%SOLID%SHELL%LAYER(i)%PLUS )
          allocate ( P%SOLID%SHELL%LAYER(i)%MINUS )
          call fstr_setup_post_phys_alloc(P%SOLID%SHELL%LAYER(i)%PLUS , 3, P%MESH%n_node, P%MESH%n_elem)
          call fstr_setup_post_phys_alloc(P%SOLID%SHELL%LAYER(i)%MINUS, 3, P%MESH%n_node, P%MESH%n_elem)
        enddo
        phys => P%SOLID%SHELL
      else
        allocate ( P%SOLID%SOLID )
        phys => P%SOLID%SOLID
        call fstr_setup_post_phys_alloc(phys, P%MESH%n_dof, P%MESH%n_node,  P%MESH%n_elem)
      end if
      P%SOLID%STRAIN => phys%STRAIN
      P%SOLID%STRESS => phys%STRESS
      P%SOLID%MISES  => phys%MISES
      P%SOLID%ESTRAIN => phys%ESTRAIN
      P%SOLID%ESTRESS => phys%ESTRESS
      P%SOLID%EMISES  => phys%EMISES
      P%SOLID%EPLSTRAIN  => phys%EPLSTRAIN
      P%SOLID%ENQM    => phys%ENQM
      allocate( P%SOLID%REACTION( P%MESH%n_dof*P%MESH%n_node ), stat=i )
      if( i /= 0 ) stop "Allocation error: REACTION"
      P%SOLID%REACTION = 0.d0
    end if

    if( P%PARAM%fg_visual == kON )then
      call fstr_setup_visualize( ctrl, P%MESH )
    end if

    call hecmw_barrier( P%MESH ) ! JP-7

    if( P%HEAT%STEPtot == 0 ) then ! No !HEAT Input
      if( P%PARAM%analysis_n == 0 ) then  ! No !STATIC Input
        call reallocate_real( P%PARAM%dtime, 1)
        call reallocate_real( P%PARAM%etime, 1)
        call reallocate_real( P%PARAM%dtmin, 1)
        call reallocate_real( P%PARAM%delmax,1)
        call reallocate_integer( P%PARAM%itmax, 1)
        call reallocate_real( P%PARAM%eps,   1)
        P%PARAM%analysis_n = 1
        P%PARAM%dtime = 0
        P%PARAM%etime = 0
        P%PARAM%dtmin = 0
        P%PARAM%delmax = 0
        P%PARAM%itmax = 20
        P%PARAM%eps = 1.0e-6
      end if
      P%HEAT%STEPtot = 1
      call reallocate_real( P%HEAT%STEP_DLTIME, 1)
      call reallocate_real( P%HEAT%STEP_EETIME, 1)
      call reallocate_real( P%HEAT%STEP_DELMIN, 1)
      call reallocate_real( P%HEAT%STEP_DELMAX, 1)
      P%HEAT%STEP_DLTIME = 0
      P%HEAT%STEP_EETIME = 0
      P%HEAT%STEP_DELMIN = 0
      P%HEAT%STEP_DELMAX = 0
    end if
  end subroutine fstr_setup_post

  !-----------------------------------------------------------------------------!
  !> Convert SURF-SURF contact to NODE-SURF contact                             !
  !-----------------------------------------------------------------------------!

  subroutine fstr_convert_contact_type( hecMESH )
    implicit none
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    integer(kind=kint) :: n, i, sgrp_id, ngrp_id
    ! convert SURF_SURF to NODE_SURF
    n = hecMESH%contact_pair%n_pair
    do i = 1,n
      if( hecMESH%contact_pair%type(i) /= HECMW_CONTACT_TYPE_SURF_SURF ) cycle
      sgrp_id = hecMESH%contact_pair%slave_grp_id(i)
      call append_node_grp_from_surf_grp( hecMESH, sgrp_id, ngrp_id )
      ! change type of contact and slave group ID
      hecMESH%contact_pair%type(i) = HECMW_CONTACT_TYPE_NODE_SURF
      hecMESH%contact_pair%slave_grp_id(i) = ngrp_id
    enddo
  end subroutine fstr_convert_contact_type

  !-----------------------------------------------------------------------------!
  !> Move dof 4-6 of 641/761/781 nodes to their rotational dummy nodes          !
  !-----------------------------------------------------------------------------!

  !> A translational node of 641/761/781 has dof 1-3 only, and its rotation is dof 1-3
  !> of the paired dummy node. Dof 4-6 of !BOUNDARY and !CLOAD given to the translational
  !> nodes are moved to the dummy nodes generated from 611/731/741; the other conditions
  !> cannot give dof 4-6 to them. With explicit dummy nodes in the mesh, the rotation is
  !> given only as dof 1-3 of the dummy nodes.
  subroutine fstr_convert_rotation_dof( P )
    use hecmw_setup_util, only : append_new_groups, get_grp_id
    implicit none
    type(fstr_param_pack) :: P
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(kind=kint), allocatable :: rot_node(:), n_rot(:), n_rot_all(:), n_explicit(:), rot_grp(:), src(:)
    integer(kind=kint), allocatable :: new_index(:), new_list(:)
    character(len=HECMW_NAME_LEN), allocatable :: new_name(:)
    logical, allocatable :: is_rot(:), in_all(:), is_explicit(:)
    integer(kind=kint) :: icel, is, ie, nn, i, j, k, m, n, ig, n_new, first_id, ids, ide

    hecMESH => P%MESH
    if( hecMESH%n_dof /= 3 ) return

    ! explicit dummy nodes are in ALL, and the generated ones are not
    allocate( in_all(hecMESH%n_node) )
    in_all = .false.
    ig = get_grp_id( hecMESH, 'node_grp', 'ALL' )
    if( ig > 0 ) then
      is = hecMESH%node_group%grp_index(ig-1)
      in_all(hecMESH%node_group%grp_item(is+1:hecMESH%node_group%grp_index(ig))) = .true.
    endif

    allocate( rot_node(hecMESH%n_node), is_explicit(hecMESH%n_node) )
    rot_node = 0
    is_explicit = .false.
    do icel = 1, hecMESH%n_elem
      if( .not. hecmw_is_etype_33struct(hecMESH%elem_type(icel)) ) cycle
      nn = hecmw_get_max_node(hecMESH%elem_type(icel))/2
      is = hecMESH%elem_node_index(icel-1)
      do j = 1, nn
        rot_node(hecMESH%elem_node_item(is+j)) = hecMESH%elem_node_item(is+nn+j)
        if( in_all(hecMESH%elem_node_item(is+nn+j)) ) is_explicit(hecMESH%elem_node_item(is+j)) = .true.
      enddo
    enddo
    ! a rank may hold an external translational node without its 641/761/781 element,
    ! so the conditions are converted on every rank by the counts over all the ranks
    n = count( rot_node > 0 )
    call hecmw_allreduce_I1( hecMESH, n, hecmw_max )
    if( n == 0 ) return

    allocate( n_rot(hecMESH%node_group%n_grp), n_rot_all(hecMESH%node_group%n_grp) )
    allocate( n_explicit(hecMESH%node_group%n_grp) )
    do ig = 1, hecMESH%node_group%n_grp
      is = hecMESH%node_group%grp_index(ig-1)
      n_rot(ig) = count( rot_node(hecMESH%node_group%grp_item(is+1:hecMESH%node_group%grp_index(ig))) > 0 )
      n_explicit(ig) = count( is_explicit(hecMESH%node_group%grp_item(is+1:hecMESH%node_group%grp_index(ig))) )
    enddo
    n_rot_all = n_rot
    call hecmw_allreduce_I( hecMESH, n_rot_all, hecMESH%node_group%n_grp, hecmw_max )

    do i = 1, P%SOLID%SPRING_ngrp_tot
      if( P%SOLID%SPRING_ngrp_DOF(i) > 3 .and. n_rot_all(P%SOLID%SPRING_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!SPRING' )
    enddo
    do i = 1, P%SOLID%VELOCITY_ngrp_tot
      if( P%SOLID%VELOCITY_ngrp_rotID(i) > 0 ) cycle
      if( mod(P%SOLID%VELOCITY_ngrp_type(i), 10) > 3 .and. n_rot_all(P%SOLID%VELOCITY_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!VELOCITY' )
    enddo
    do i = 1, P%SOLID%ACCELERATION_ngrp_tot
      if( mod(P%SOLID%ACCELERATION_ngrp_type(i), 10) > 3 .and. n_rot_all(P%SOLID%ACCELERATION_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!ACCELERATION' )
    enddo
    do i = 1, P%FREQ%FLOAD_ngrp_tot
      if( P%FREQ%FLOAD_ngrp_TYPE(i) /= kFLOADTYPE_NODE ) cycle
      if( P%FREQ%FLOAD_ngrp_DOF(i) > 3 .and. n_rot_all(P%FREQ%FLOAD_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!FLOAD' )
    enddo
    !the rotational DOFs of !FACCELERATION are also not allowed (the end DOF is mod(type,10))
    do i = 1, P%FREQ%FACC_ngrp_tot
      if( mod(P%FREQ%FACC_ngrp_type(i), 10) > 3 .and. n_rot_all(P%FREQ%FACC_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!FACCELERATION' )
    enddo
    do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if( hecMESH%mpc%mpc_dof(j) > 3 .and. rot_node(hecMESH%mpc%mpc_item(j)) > 0 ) &
          call rotation_dof_err_stop( '!EQUATION' )
      enddo
    enddo

    ! node groups whose rotation is moved: rot_grp(ig) = -1, then the id of the dummy node group
    allocate( rot_grp(hecMESH%node_group%n_grp) )
    rot_grp = 0
    do i = 1, P%SOLID%BOUNDARY_ngrp_tot
      if( P%SOLID%BOUNDARY_ngrp_rotID(i) > 0 ) cycle
      ig = P%SOLID%BOUNDARY_ngrp_ID(i)
      if( mod(P%SOLID%BOUNDARY_ngrp_type(i), 10) <= 3 .or. n_rot_all(ig) == 0 ) cycle
      if( n_explicit(ig) > 0 ) call explicit_dummy_node_err_stop( '!BOUNDARY' )
      rot_grp(ig) = -1
    enddo
    do i = 1, P%SOLID%CLOAD_ngrp_tot
      if( P%SOLID%CLOAD_ngrp_rotID(i) > 0 .or. P%SOLID%CLOAD_ngrp_DOF(i) <= 3 ) cycle
      ig = P%SOLID%CLOAD_ngrp_ID(i)
      is = hecMESH%node_group%grp_index(ig-1)
      ie = hecMESH%node_group%grp_index(ig)
      if( any( rot_node(hecMESH%node_group%grp_item(is+1:ie)) == 0 .and. &
        hecMESH%node_group%grp_item(is+1:ie) <= hecMESH%nn_internal ) ) call fstr_setup_util_err_stop( &
        '### Error: !CLOAD : dof 4-6 is given to nodes without rotational dof' )
      if( n_explicit(ig) > 0 ) call explicit_dummy_node_err_stop( '!CLOAD' )
      rot_grp(ig) = -1
    enddo

    n_new = count( rot_grp == -1 )
    if( n_new == 0 ) return
    allocate( new_name(n_new), new_index(0:n_new), new_list(sum(n_rot, mask=rot_grp == -1)) )
    new_index(0) = 0
    k = 0
    m = 0
    do ig = 1, hecMESH%node_group%n_grp
      if( rot_grp(ig) /= -1 ) cycle
      k = k + 1
      write( new_name(k), '(a,i0)' ) 'FSTR_ROT_', ig
      do j = hecMESH%node_group%grp_index(ig-1)+1, hecMESH%node_group%grp_index(ig)
        if( rot_node(hecMESH%node_group%grp_item(j)) == 0 ) cycle
        m = m + 1
        new_list(m) = rot_node(hecMESH%node_group%grp_item(j))
      enddo
      new_index(k) = m
      rot_grp(ig) = k
    enddo
    call append_new_groups( hecMESH, 'node_grp', n_new, new_name, new_index, new_list, first_id )
    where( rot_grp > 0 ) rot_grp = rot_grp + first_id - 1

    ! split each !BOUNDARY entry into dof 1-3 of its nodes and dof 1-3 of their dummy nodes
    n = P%SOLID%BOUNDARY_ngrp_tot
    if( n > 0 ) then
      allocate( src(2*n), is_rot(2*n) )
      m = 0
      do i = 1, n
        ide = mod(P%SOLID%BOUNDARY_ngrp_type(i), 10)
        if( P%SOLID%BOUNDARY_ngrp_rotID(i) <= 0 .and. ide > 3 .and. rot_grp(P%SOLID%BOUNDARY_ngrp_ID(i)) > 0 ) then
          if( P%SOLID%BOUNDARY_ngrp_type(i)/10 <= 3 ) then
            m = m + 1
            src(m) = i
            is_rot(m) = .false.
          endif
          m = m + 1
          src(m) = i
          is_rot(m) = .true.
        else
          m = m + 1
          src(m) = i
          is_rot(m) = .false.
        endif
      enddo
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_GRPID, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_ID, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_type, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_amp, src, m )
      call remap_real_array( P%SOLID%BOUNDARY_ngrp_val, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_istot, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_rotID, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_centerID, src, m )
      P%SOLID%BOUNDARY_ngrp_tot = m
      do k = 1, m
        if( P%SOLID%BOUNDARY_ngrp_rotID(k) > 0 ) cycle
        ig = P%SOLID%BOUNDARY_ngrp_ID(k)
        ids = P%SOLID%BOUNDARY_ngrp_type(k)/10
        ide = mod(P%SOLID%BOUNDARY_ngrp_type(k), 10)
        if( ide <= 3 .or. rot_grp(ig) <= 0 ) cycle
        if( is_rot(k) ) then
          P%SOLID%BOUNDARY_ngrp_ID(k) = rot_grp(ig)
          P%SOLID%BOUNDARY_ngrp_type(k) = 10 * (max(ids, 4) - 3) + (ide - 3)
        else
          P%SOLID%BOUNDARY_ngrp_type(k) = 10 * ids + 3
        endif
      enddo
    endif

    do i = 1, P%SOLID%CLOAD_ngrp_tot
      if( P%SOLID%CLOAD_ngrp_rotID(i) > 0 .or. P%SOLID%CLOAD_ngrp_DOF(i) <= 3 ) cycle
      P%SOLID%CLOAD_ngrp_ID(i) = rot_grp(P%SOLID%CLOAD_ngrp_ID(i))
      P%SOLID%CLOAD_ngrp_DOF(i) = P%SOLID%CLOAD_ngrp_DOF(i) - 3
    enddo
  end subroutine fstr_convert_rotation_dof

  subroutine rotation_dof_err_stop( header_name )
    implicit none
    character(len=*) :: header_name

    call fstr_setup_util_err_stop( '### Error: '//header_name// &
      ' : dof 4-6 of beam/shell nodes in a mesh with solid elements is not supported' )
  end subroutine rotation_dof_err_stop

  subroutine explicit_dummy_node_err_stop( header_name )
    implicit none
    character(len=*) :: header_name

    call fstr_setup_util_err_stop( '### Error: '//header_name// &
      ' : dof 4-6 is given to beam/shell nodes with explicit dummy nodes; give it as dof 1-3 of the dummy nodes' )
  end subroutine explicit_dummy_node_err_stop

  subroutine remap_integer_array( array, src, n )
    implicit none
    integer(kind=kint), pointer :: array(:)
    integer(kind=kint), intent(in) :: src(:), n
    integer(kind=kint), pointer :: new_array(:)

    allocate( new_array(n) )
    new_array(1:n) = array(src(1:n))
    deallocate( array )
    array => new_array
  end subroutine remap_integer_array

  subroutine remap_real_array( array, src, n )
    implicit none
    real(kind=kreal), pointer :: array(:)
    integer(kind=kint), intent(in) :: src(:), n
    real(kind=kreal), pointer :: new_array(:)

    allocate( new_array(n) )
    new_array(1:n) = array(src(1:n))
    deallocate( array )
    array => new_array
  end subroutine remap_real_array

end module m_fstr_setup
