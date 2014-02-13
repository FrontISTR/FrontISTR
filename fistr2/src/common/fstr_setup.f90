!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : I/O and Utility                                   !
!                                                                      !
!            Written by Noboru Imai (Univ. of Tokyo)                   !
!                       Xi YUAN ( Advancesoft )                        !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module provides functions to read in data from control file
!! and do neccessary preparation for following calculation
module m_fstr_setup
use m_fstr
use fstr_setup_util
use fstr_ctrl_common
use fstr_ctrl_static
use fstr_ctrl_heat
use fstr_ctrl_eigen
use fstr_ctrl_dynamic
use fstr_ctrl_material
use m_static_get_prop
use m_out
use m_step
use m_utilities

implicit none


include 'fstr_ctrl_util_f.inc'

        !> Package of all data needs to initilize
        type fstr_param_pack
                type(hecmwST_local_mesh), pointer :: MESH     
                type(fstr_param), pointer         :: PARAM   
                type(fstr_solid), pointer         :: SOLID
                type(fstr_heat), pointer          :: HEAT
                type(lczparam), pointer           :: EIGEN
                type(fstr_dynamic), pointer       :: DYN
                type(fstr_couple), pointer        :: CPL
        end type fstr_param_pack

contains

!=============================================================================!
!> Read in and initialize control data                                                             !
!=============================================================================!
subroutine fstr_setup( cntl_filename, hecMESH, fstrPARAM,  &
        fstrSOLID, fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL )
        use mMaterial
        character(len=*) :: cntl_filename
        type(hecmwST_local_mesh),target :: hecMESH
        type(fstr_param),target   :: fstrPARAM
        type(fstr_solid),target   :: fstrSOLID
        type(lczparam),target     :: fstrEIG
        type(fstr_heat),target    :: fstrHEAT
        type(fstr_dynamic),target :: fstrDYNAMIC
        type(fstr_couple),target  :: fstrCPL
		
        integer(kind=kint) :: ctrl
        type(fstr_param_pack) :: P

        integer, parameter :: MAXOUTFILE = 10

        external fstr_ctrl_get_c_h_name
        integer(kind=kint) :: fstr_ctrl_get_c_h_name

        integer(kind=kint) :: version, resul, visual, femap
        integer(kind=kint) :: rcode, n, i, j, ierror
        character(len=HECMW_NAME_LEN) :: header_name, fname(MAXOUTFILE)
        real(kind=kreal) :: ee, pp, rho, alpha, thick
        logical          :: isOK
        type(t_output_ctrl) :: outctrl
        character(len=HECMW_FILENAME_LEN) :: logfileNAME, mName

        ! counters
        integer(kind=kint) :: c_solution, c_solver, c_step, c_write
        integer(kind=kint) :: c_static, c_boundary, c_cload, c_dload, c_temperature, c_reftemp
        integer(kind=kint) :: c_heat, c_fixtemp, c_cflux, c_dflux, c_sflux, c_film, c_sfilm, c_radiate, c_sradiate
        integer(kind=kint) :: c_eigen, c_contact
        integer(kind=kint) :: c_dynamic, c_velocity, c_acceleration
        integer(kind=kint) :: c_couple, c_material
        integer(kind=kint) :: c_mpc, c_amp
        integer(kind=kint) :: c_istep, c_section
        integer(kind=kint) :: c_output, cr_output, islog

        write( logfileNAME, '(i5,''.log'')') myrank

        ! packaging
        P%MESH   => hecMESH
        P%PARAM  => fstrPARAM
        P%SOLID  => fstrSOLID
        P%EIGEN  => fstrEIG
        P%HEAT   => fstrHEAT
        P%DYN    => fstrDYNAMIC
        P%CPL    => fstrCPL

        c_solution = 0; c_solver   = 0; c_step   = 0; c_write = 0
        c_static   = 0; c_boundary = 0; c_cload  = 0; c_dload = 0; c_temperature = 0; c_reftemp = 0
        c_heat     = 0; c_fixtemp  = 0; c_cflux  = 0; c_dflux = 0; c_sflux = 0
        c_film     = 0; c_sfilm    = 0; c_radiate= 0; c_sradiate = 0
        c_eigen    = 0; c_contact  = 0
        c_dynamic  = 0; c_velocity = 0; c_acceleration = 0
        c_couple   = 0; c_material = 0; c_section = 0
        c_mpc      = 0; c_amp      = 0;
        c_istep    = 0
        c_output   = 0
		
        call fstr_setup_init( P )

        ctrl = fstr_ctrl_open( cntl_filename//char(0) )
        if( ctrl < 0 ) then
                write(*,*) '### Error: Cannot open FSTR control file : ', cntl_filename
                write(ILOG,*) '### Error: Cannot open FSTR control file : ', cntl_filename
                STOP
        end if

        version =0
        do
                rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )
                if(     header_name == '!VERSION' ) then
                        rcode = fstr_ctrl_get_data_array_ex( ctrl, 'i ', version )
                else if(     header_name == '!SOLUTION' ) then
                        c_solution = c_solution + 1
                        call fstr_setup_SOLUTION( ctrl, c_solution, P )
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
        !        else if( header_name == '!WRITE' ) then
        !                c_write = c_write + 1
        !                call fstr_setup_WRITE( ctrl, c_write, P )
                else if( header_name == '!ECHO' ) then
                        P%PARAM%fg_echo = 1
                else if( header_name == '!RESTART' ) then
                        call fstr_setup_RESTART( ctrl, n, restartfilNAME )
                        fstrSOLID%restart_nout= n
!
                !--------------- for static -------------------------
!
                else if( header_name == '!STATIC' ) then
                        c_static = c_static + 1
                        call fstr_setup_STATIC( ctrl, c_static, P )
                else if( header_name == '!BOUNDARY' ) then
                  n = fstr_ctrl_get_data_line_n( ctrl )
                  c_boundary = c_boundary + n
                else if( header_name == '!CLOAD' ) then
                  n = fstr_ctrl_get_data_line_n( ctrl )
                  c_cload = c_cload + n
                else if( header_name == '!DLOAD' ) then
                        n = fstr_ctrl_get_data_line_n( ctrl )
                        c_dload = c_dload + n
                else if( header_name == '!SECTION' ) then
                        c_section = c_section + 1
                else if( header_name == '!AMPLITUDE' ) then
                        c_amp = c_amp + 1
                else if( header_name == '!MATERIAL' ) then
                        c_material = c_material + 1
                else if( header_name == '!TEMPERATURE' ) then
                        n = fstr_ctrl_get_data_line_n( ctrl )
                        c_temperature = c_temperature + n
                else if( header_name == '!REFTEMP' ) then
                        c_reftemp = c_reftemp + 1
                        call fstr_setup_REFTEMP( ctrl, c_reftemp, P )
                else if( header_name == '!WRITE' ) then
                        call fstr_ctrl_get_output( ctrl, outctrl, islog, resul, visual, femap )
                        if( visual==1 ) P%PARAM%fg_visual=1
                        if( femap==1 ) P%PARAM%fg_neutral=1
                        if( resul==1 ) then
                          P%PARAM%fg_result = 1
                          if( islog == 1 ) outctrl%filename = trim(logfileNAME)
                          isOK= .true.
                          do i=1,c_output
                            if( trim( fname(i) ) == trim( outctrl%filename ) ) then
                              isOK=.false.; exit
                            endif
                          enddo
                          if( isOK ) then
                            c_output = c_output+1
                            if( c_output>MAXOUTFILE ) stop "Too many output files(>10)"
                            fname( c_output ) = trim( outctrl%filename )
                          endif
                        endif

                !--------------- for heat -------------------------

                else if( header_name == '!HEAT' ) then
                        c_heat = c_heat + 1
                        call fstr_setup_HEAT( ctrl, c_heat, P )
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

                !--------------- for couple -------------------------

                else if( header_name == '!COUPLE' ) then
                        c_couple = c_couple + 1
                        call fstr_setup_COUPLE( ctrl, c_couple, P )
                !--------------- for mpc -------------------------

                else if( header_name == '!MPC' ) then
                        c_mpc = c_mpc + 1
                        call fstr_setup_MPC( ctrl, c_mpc, P )

                !--------------------- END -------------------------

                else if( header_name == '!END' ) then
                        exit
                end if

                ! next
                if( fstr_ctrl_seek_next_header(ctrl) == 0) exit
        end do
!
! ----- 
        if( c_section==0 ) stop "SECTION not defined!"
        allocate( MWSections( c_section ) )
        do i=1,c_section
          MWSections(i)%sect_R_item = 0.d0
        enddo
        if( c_amp>0 ) allocate( MWAmplitudes( c_amp ) )
        if( c_istep>0 ) allocate( fstrSOLID%step_ctrl( c_istep ) )
        if( c_output > 0 ) allocate( fstrSOLID%output_ctrl( c_output ) )
        if( c_material==0 ) stop "material property not defined!"
        allocate( fstrSOLID%materials( c_material ) )
        if( c_cload>0 ) then
          allocate(fstrSOLID%cload_grp(c_cload) )
          do i=1,c_cload
            call init_boundary_grp( fstrSOLID%cload_grp(i) )
          enddo
        endif
        if( c_dload>0 ) then
          allocate(fstrSOLID%dload_grp(c_dload) )
          do i=1,c_dload
            call init_dload_grp( fstrSOLID%dload_grp(i) )
          enddo
        endif
        if( c_boundary>0 ) then
          allocate(fstrSOLID%boundary_grp(c_boundary) )
          do i=1,c_boundary
            call init_boundary_grp( fstrSOLID%boundary_grp(i) )
          enddo
        endif
        if( c_temperature>0 ) then
          allocate(fstrSOLID%temp_grp(c_temperature) )
          do i=1,c_temperature
            call init_ndscalar_grp( fstrSOLID%temp_grp(i) )
          enddo
        endif
!
! ----- 
        rcode = fstr_ctrl_rewind( ctrl )
!
        c_istep    = 0
        c_material = 0
        c_output   = 0
        c_section = 0
        c_cload = 0
        c_dload = 0
        c_boundary = 0
        c_amp = 0
        c_temperature = 0
        do
          rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )

           if( header_name == '!ISTEP'  ) then
            c_istep = c_istep+1
            if( .not. fstr_ctrl_get_ISTEP( ctrl, hecMESH, fstrSOLID%step_ctrl(c_istep) ) ) then
                write(*,*) '### Error: Fail in read in step definition : ' , c_istep
                write(ILOG,*) '### Error: Fail in read in step definition : ', c_istep
                stop
            endif
          else if( header_name == '!STEP' .and. version>=1 ) then
            c_istep = c_istep+1
            if( .not. fstr_ctrl_get_ISTEP( ctrl, hecMESH, fstrSOLID%step_ctrl(c_istep) ) ) then
                write(*,*) '### Error: Fail in read in step definition : ' , c_istep
                write(ILOG,*) '### Error: Fail in read in step definition : ', c_istep
                stop
            endif
			
          else if( header_name == '!CLOAD' ) then
            if( fstr_ctrl_get_CLOAD( ctrl, c_cload, fstrSOLID%cload_grp )/=0 ) then
                write(*,*) '### Error: Fail in read in CLOAD definition : ' , c_cload
                write(ILOG,*) '### Error: Fail in read in CLOAD definition : ', c_cload
                stop
            endif
			
          else if( header_name == '!DLOAD' ) then
            if( fstr_ctrl_get_DLOAD( ctrl, c_dload, fstrSOLID%dload_grp )/=0 ) then
                write(*,*) '### Error: Fail in read in DLOAD definition : ' , c_dload
                write(ILOG,*) '### Error: Fail in read in DLOAD definition : ', c_dload
                stop
            endif
			
          else if( header_name == '!BOUNDARY' ) then
            if( fstr_ctrl_get_BOUNDARY( ctrl, c_boundary, fstrSOLID%boundary_grp )/=0 ) then
                write(*,*) '### Error: Fail in read in BOUNDARY definition : ' , c_boundary
                write(ILOG,*) '### Error: Fail in read in BOUNDARY definition : ', c_boundary
                stop
            endif
			
          else if( header_name == '!TEMPERATURE' ) then
            if( fstr_ctrl_get_TEMPERATURE( ctrl, c_temperature, fstrSOLID%temp_grp,  &
                      fstrSOLID%TEMP_irres, fstrSOLID%TEMP_tstep)/=0 ) then
                write(*,*) '### Error: Fail in read in TEMPERATURE definition : ' , c_temperature
                write(ILOG,*) '### Error: Fail in read in TEMPERATURE definition : ', c_temperature
                stop
            endif
			
          else if( header_name == '!AMPLITUDE' ) then
            c_amp = c_amp + 1
            if( fstr_ctrl_get_AMPLITUDE( ctrl, MWAmplitudes(c_amp) )/=0 ) then
                write(*,*) '### Error: Fail in read in AMPLITUDE definition : ' , c_amp
                write(ILOG,*) '### Error: Fail in read in AMPLITUDE definition : ', c_amp
                stop
            endif
			
          else if( header_name == '!SECTION' ) then
            c_section = c_section + 1
            if( fstr_ctrl_get_SECTION( ctrl, MWSections(c_section) )/=0 ) then
                write(*,*) '### Error: Fail in read in section definition : ' , c_section
                write(ILOG,*) '### Error: Fail in read in section definition : ', c_section
                stop
            endif

          else if( header_name == '!MATERIAL' ) then
            c_material = c_material+1
            if( fstr_ctrl_get_MATERIAL( ctrl, mName )/=0 ) then
                write(*,*) '### Error: Fail in read in material definition : ' , c_material
                write(ILOG,*) '### Error: Fail in read in material definition : ', c_material
                stop
            endif
            fstrSOLID%materials(c_material)%name = mName
            call initMaterial( fstrSOLID%materials(c_material) )
          else if( header_name == '!ELASTIC' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_ELASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(c_material)%mtype,       &
                                             fstrSOLID%materials(c_material)%nlgeom_flag, &
                                             fstrSOLID%materials(c_material)%variables,   &
                                             fstrSOLID%materials(c_material)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in elasticity definition : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in elasticity definition : ', c_material
                 stop
              endif
            endif
          else if( header_name == '!PLASTIC' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_PLASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(c_material)%mtype,       &
                                             fstrSOLID%materials(c_material)%nlgeom_flag, &
                                             fstrSOLID%materials(c_material)%variables,   &
                                             fstrSOLID%materials(c_material)%table,   &
                                             fstrSOLID%materials(c_material)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in plasticity definition : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in plasticity definition : ', c_material
                 stop
              endif
            endif
          else if( header_name == '!HYPERELASTIC' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_HYPERELASTIC( ctrl,                                        &
                                             fstrSOLID%materials(c_material)%mtype,       &
                                             fstrSOLID%materials(c_material)%nlgeom_flag, &
                                             fstrSOLID%materials(c_material)%variables )/=0 ) then
                 write(*,*) '### Error: Fail in read in elasticity definition : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in elasticity definition : ', c_material
                 stop
              endif
            endif
          else if( header_name == '!VISCOELASTIC' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_VISCOELASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(c_material)%mtype,       &
                                             fstrSOLID%materials(c_material)%nlgeom_flag, &
                                             fstrSOLID%materials(c_material)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in plasticity definition : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in plasticity definition : ', c_material
                 stop
              endif
            endif
          else if( header_name == '!DENSITY' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_DENSITY( ctrl, fstrSOLID%materials(c_material)%variables )/=0 ) then
                 write(*,*) '### Error: Fail in read in density definition : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in density definition : ', c_material
                 stop
               endif
            endif
          else if( header_name == '!EXPANSION_COEFF' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_EXPANSION_COEFF( ctrl, fstrSOLID%materials(c_material)%variables, &
                                             fstrSOLID%materials(c_material)%dict)/=0 )  then
                 write(*,*) '### Error: Fail in read in expansion coefficient definition : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in expansion coefficient definition : ', c_material
                 stop
               endif
            endif
          else if( header_name == '!USER_MATERIAL' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_USERMATERIAL( ctrl, fstrSOLID%materials(c_material)%mtype,   &
                 fstrSOLID%materials(c_material)%nlgeom_flag, fstrSOLID%materials(c_material)%nfstatus, &
                 fstrSOLID%materials(c_material)%variables(101:) )/=0 ) then
                 write(*,*) '### Error: Fail in read in user defined material : ' , c_material
                 write(ILOG,*) '### Error: Fail in read in user defined material : ', c_material
                 stop
               endif
            endif


! == Following output control ==
          else if( header_name == '!WRITE' ) then
            call fstr_ctrl_get_output( ctrl, outctrl, islog, resul, visual, femap)
            if( resul==1 ) then
              if( islog == 1 ) then
                outctrl%filename = trim(logfileNAME)
                outctrl%filenum = ILOG
              endif
              isOK = .true.
              do i=1,c_output-1
               if( trim( fstrSOLID%output_ctrl(i)%filename ) == trim( outctrl%filename ) ) then
                    isOK=.false.;  cr_output=i;  exit
               endif
              enddo
              if( isOK ) then 
                c_output = c_output + 1
                cr_output=c_output
                if( islog/=1 ) then 
                  outctrl%filenum = IRESOUT+c_output
                  call append_int2name( myrank, outctrl%filename )
                endif
                call fstr_init_outctrl(fstrSOLID%output_ctrl(cr_output))
                call fstr_copy_outctrl(fstrSOLID%output_ctrl(cr_output), outctrl)
                if( islog/=1 ) open( unit=outctrl%filenum, file=outctrl%filename, status='REPLACE' )
!              call print_output_ctrl( 6, fstrSOLID%output_ctrl(cr_output) )
              endif
            endif
			
          else if( header_name == '!OUTPUT_TYPE=VTK' ) then
		    gVisType=1
			
          else if( header_name == '!OUTPUT_TYPE=MICROAVS' ) then
		    gVisType=2
			
          else if( header_name == '!OUTPUT_TYPE=FVUNS' ) then
		    gVisType=3
            
          else if( header_name == '!NODE_OUTPUT' ) then
            if( c_output >0 ) then
               if( .not. fstr_ctrl_get_outnode( ctrl, hecMESH, fstrSOLID%output_ctrl(cr_output)%outinfo ) ) then
                 write(*,*) '### Error: Fail in read in node output definition : ' , cr_output
                 write(ILOG,*) '### Error: Fail in read in node output definition : ', cr_output
                 stop
               endif
            endif
          else if( header_name == '!ELEMENT_OUTPUT' ) then
            if( c_output >0 ) then
               if( .not. fstr_ctrl_get_outelem( ctrl, hecMESH, fstrSOLID%output_ctrl(cr_output)%outinfo ) ) then
                 write(*,*) '### Error: Fail in read in element output definition : ' , cr_output
                 write(ILOG,*) '### Error: Fail in read in element output definition : ', cr_output
                 stop
               endif
            endif
          else if( header_name == '!ULOAD' ) then
            if( fstr_ctrl_get_USERLOAD( ctrl )/=0 ) then
                 write(*,*) '### Error: Fail in read in ULOAD definition : ' 
                 write(ILOG,*) '### Error: Fail in read in ULOAD definition : '
                 stop
            endif

          else if( header_name == '!END' ) then
                        exit
          endif
!
          if( fstr_ctrl_seek_next_header(ctrl) == 0) exit
        end do
!
! ----- material type judgement. in case of infinitive analysis, nlgeom_flag=0
        if( P%PARAM%solution_type == kstSTATIC ) then
          do i=1, c_material
            fstrSOLID%materials(c_material)%nlgeom_flag = 0
     !     call printMaterial( 6, fstrSOLID%materials(i) ); pause
          enddo
        endif

        if( c_temperature>0 .or. fstrSOLID%TEMP_irres == 1 ) then 
          allocate ( fstrSOLID%temperature( total_node )      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TEMPERATURE>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          fstrSOLID%temperature = REF_TEMP
          allocate ( fstrSOLID%reftemp( total_node )      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, REFTEMP>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          fstrSOLID%reftemp = REF_TEMP
        endif

        fstrSOLID%nstep_tot = 1
        if( associated(fstrSOLID%step_ctrl) )  then
           fstrSOLID%nstep_tot = size(fstrSOLID%step_ctrl)
           !call fstr_print_steps( 6, fstrSOLID%step_ctrl )
        else
           print *, "Step control not defined! Using defualt step=1"
           allocate( fstrSOLID%step_ctrl(1) )
           fstrSOLID%nstep_tot = 1
           allocate( fstrSOLID%step_ctrl(1) )
           call init_stepInfo( fstrSOLID%step_ctrl(1) )
           if( c_boundary>0 ) allocate( fstrSOLID%step_ctrl(1)%Boundary(c_boundary) )
           do i = 1, c_boundary
             fstrSOLID%step_ctrl(1)%Boundary(i) = fstrSOLID%BOUNDARY_grp(i)%gid
           enddo
           n = c_cload + c_dload + c_temperature
           if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%Load(n) )
           do i = 1, c_cload
             fstrSOLID%step_ctrl(1)%Load(i) = fstrSOLID%CLOAD_grp(i)%gid
           enddo
           do i = 1, c_dload
             fstrSOLID%step_ctrl(1)%Load(i+c_cload) = fstrSOLID%DLOAD_grp(i)%gid
           enddo
            do i = 1, c_temperature
             fstrSOLID%step_ctrl(1)%Load(i+c_cload+c_dload) = fstrSOLID%temp_grp(i)%gid
           enddo
        endif

        if( p%PARAM%solution_type /= kstHEAT) call fstr_element_init( fstrSOLID )
        call fstr_setup_post( ctrl, P )
        rcode = fstr_ctrl_close( ctrl )
!
end subroutine fstr_setup

!-----------------------------------------------------------------------------!
!> General initializer
subroutine fstr_setup_init( P )
        implicit none
        type(fstr_param_pack), target :: P

        ! ----  default setting of global params ---

        DT = 1
        ETIME = 1
        ITMAX = 20
        EPS = 1.0e-6

        ! -------  grobal pointer setting ----------
!        IPROC    => P%PARAM%iproc
        INCMAX   => P%PARAM%incmax
        REF_TEMP => P%PARAM%ref_temp
        IECHO    => P%PARAM%fg_echo
        IRESULT  => P%PARAM%fg_result
        IVISUAL  => P%PARAM%fg_visual

        ! for heat ...
        INEUTRAL => P%PARAM%fg_neutral
        IRRES    => P%PARAM%fg_irres
        IWRES    => P%PARAM%fg_iwres
        NRRES    => P%PARAM%nrres
        NPRINT   => P%PARAM%nprint
 
        call fstr_param_init( P%PARAM, P%MESH )
        call fstr_solid_init( P%SOLID )
        call fstr_eigen_init( P%EIGEN )
        call fstr_heat_init ( P%HEAT  )
        call fstr_dynamic_init ( P%DYN  )

end subroutine fstr_setup_init

!> Initializer of structure fstr_solid
subroutine fstr_solid_init( fstrSOLID )
        use m_fstr
        type(fstr_solid)                :: fstrSOLID

        fstrSOLID%file_type  = kbcfFSTR

        fstrSOLID%BOUNDARY_ngrp_tot = 0
        fstrSOLID%CLOAD_ngrp_tot    = 0
        fstrSOLID%DLOAD_ngrp_tot    = 0
        fstrSOLID%TEMP_ngrp_tot     = 0
        fstrSOLID%TEMP_irres        = 0
        fstrSOLID%TEMP_tstep        = 1
        fstrSOLID%VELOCITY_ngrp_tot = 0
        fstrSOLID%ACCELERATION_ngrp_tot = 0
        fstrSOLID%COUPLE_ngrp_tot   = 0

        fstrSOLID%restart_nout= 0
end subroutine fstr_solid_init

!> Allocator of structure fstr_solid	
subroutine fstr_solid_alloc( fstrSOLID )
        use m_fstr
        include "HEC_MW3_For.h"
        type(fstr_solid)                :: fstrSOLID
		
		integer :: ntotal, ierror

        ntotal=assdof(1)*total_node
        allocate ( fstrSOLID%GL( ntotal )          ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, GL>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%unode( ntotal )  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, unode>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%dunode( ntotal )  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, dunode>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%ddunode( ntotal )  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, ddunode>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%QFORCE( ntotal )      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, QFORCE>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

        fstrSOLID%GL(:)=0.d0
        fstrSOLID%unode(:)      = 0.d0
        fstrSOLID%dunode(:)     = 0.d0
        fstrSOLID%ddunode(:)    = 0.d0
        fstrSOLID%QFORCE(:)     = 0.d0
        fstrSOLID%FACTOR( 1:2 ) = 0.d0


        ! initialize for linear static problems 
        fstrSOLID%FACTOR(2)=1.d0
        fstrSOLID%FACTOR(1)=0.d0
end subroutine fstr_solid_alloc

!> Initialize elements info in static calculation
subroutine fstr_element_init( fstrSOLID )
        use elementInfo
        use mMechGauss
        use m_fstr
        type(fstr_solid)                :: fstrSOLID
		
        include "HEC_MW3_For.h"
        integer :: ii,i, j, ng, ndof, tcnt, iAss, iPart, iElem, ic_type, cid
        integer :: igrp, ngrp, isect, csect, scid
        character(len=HECMW_NAME_LEN) :: header_name
        logical :: found
		
        call hecMW_init_section( fstrSOLID%materials )

        tcnt = 0 		
        do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            tcnt = tcnt+ mw_get_num_of_element()
         enddo
        enddo

        if( tcnt <=0 ) then
             stop "no element defined!"
        endif
        allocate( fstrSOLID%elements(tcnt) )
        tcnt = 0 		
        do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            ngrp = mw_get_num_of_elementgroup()
            do iElem = 0, mw_get_num_of_element()-1
               call mw_select_element( iElem )
               ic_type = mw_get_element_type()
               ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
               tcnt=tcnt+1
               fstrSOLID%elements(tcnt)%iAss = iAss
               fstrSOLID%elements(tcnt)%iPart = iPart
             !  fstrSOLID%elements(tcnt)%iID = mw_get_element_id()
			 
               fstrSOLID%elements(tcnt)%etype = ic_type
!          if (hecmw_is_etype_link(fstrSOLID%elements(i)%etype)) cycle
               
               found = .false.
               do igrp = 0, ngrp-1
                 ng = mw_get_elementgroup_name_length(igrp)
                 header_name(:) =''
                 call mw_get_elementgroup_name(igrp, header_name(1:ng), ng)

                 csect = -1
                 do isect=1, size(MWSections)
                     if( MWSections(isect)%egroup_name == header_name ) then
                       csect = isect;  exit
                     endif
                 enddo
                 if( csect==-1 ) then
                     write(IMSG, *) "Property of element group",iAss,iPart,igrp, "not defined"
                     write(*, *) "Property of element group",iAss,iPart,igrp, "not defined"
                     stop
                 endif
                 do i=0, mw_get_num_of_element_id(igrp)-1			   
                   cid = mw_get_element_id_with_elementgroup( igrp, i )
                   if( cid == mw_get_element_id(iElem) ) then
                     fstrSOLID%elements(tcnt)%iSect = csect
					 found=.true.; exit
                   endif
                 enddo
                 if( found ) exit
               enddo
			   if( .not. found ) then
                  write(IMSG, *) "Property of element",iAss,iPart,iElem, "not defined"
                  write(*, *) "Property of element",iAss,iPart,iElem, "not defined"
                  stop
               endif
			   
			   ndof = getSpaceDimension( fstrSOLID%elements(tcnt)%etype )
               if (ndof == 2) then              ! why do this???
                 assDOF(1) = 2
                 scid=MWSections(csect)%sect_opt
                 if( scid==0 ) then
                   fstrSOLID%elements(tcnt)%iset=1
                 else if( scid==1) then
                   fstrSOLID%elements(tcnt)%iset=0
                 else if( scid==2) then
                   fstrSOLID%elements(tcnt)%iset=2
                 endif
               else
                 assDOF(1) = 3
                 if( fstrSOLID%elements(tcnt)%etype==731 .or.            &
                     fstrSOLID%elements(tcnt)%etype==741 ) assDOF(1)=6
               endif
			   
			   cid = MWSections(csect)%sect_mat_ID
               ng = NumOfQuadPoints( fstrSOLID%elements(tcnt)%etype )
               if(ng>0) allocate( fstrSOLID%elements(tcnt)%gausses( ng ) )     
               do j=1,ng
                 fstrSOLID%elements(tcnt)%gausses(j)%pMaterial => fstrSOLID%materials(cid)
                 call fstr_init_gauss( fstrSOLID%elements(tcnt)%gausses( j )  )
               enddo
            enddo
         enddo
        enddo

end subroutine

!> Finalizer of fstr_solid
subroutine fstr_solid_finalize( fstrSOLID )
        type(fstr_solid) :: fstrSOLID
        integer :: i, j, ierror
		
        if( associated(fstrSOLID%materials) )    &
          deallocate( fstrSOLID%materials )
        if( .not. associated(fstrSOLID%elements ) ) return
        do i=1,size(fstrSOLID%elements)
          if( .not. associated(fstrSOLID%elements(i)%gausses) ) cycle
          do j=1,size(fstrSOLID%elements(i)%gausses)
            call fstr_finalize_gauss(fstrSOLID%elements(i)%gausses(j))
          enddo
          deallocate( fstrSOLID%elements(i)%gausses )
        enddo
        deallocate( fstrSOLID%elements )

        if( associated( fstrSOLID%mpc_const ) ) then
          deallocate( fstrSOLID%mpc_const )
        endif
        if( associated(fstrSOLID%step_ctrl) ) then
          do i=1,size(fstrSOLID%step_ctrl)
            call free_stepInfo( fstrSOLID%step_ctrl(i) )
          enddo
          deallocate( fstrSOLID%step_ctrl )
        endif
        if(associated(fstrSOLID%output_ctrl) ) then
           do i=1,size(fstrSOLID%output_ctrl)
              if( fstrSOLID%output_ctrl(i)%filenum/=ILOG )  &
                close(fstrSOLID%output_ctrl(i)%filenum)
           enddo
           deallocate(fstrSOLID%output_ctrl)
        endif
        
        if( associated(fstrSOLID%GL) ) then
            deallocate(fstrSOLID%GL               ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, GL>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%unode) ) then
            deallocate(fstrSOLID%unode       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, unode>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%dunode) ) then
            deallocate(fstrSOLID%dunode       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, dunode>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%ddunode) ) then
            deallocate(fstrSOLID%ddunode       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, ddunode>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%QFORCE) ) then
            deallocate(fstrSOLID%QFORCE           ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, QFORCE>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%temperature) ) then
            deallocate(fstrSOLID%temperature       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, temperature>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%reftemp) ) then
            deallocate(fstrSOLID%reftemp       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, reftemp>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif

end subroutine

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

end subroutine fstr_heat_init

!> Initial setting of eigen calculation
subroutine fstr_eigen_init( fstrEIG )
        implicit none
        type(lczparam) :: fstrEIG
        ! type(fstr_eigen) :: fstrEIG

        fstrEIG%eqset       = 0
        fstrEIG%nget        = 5
        fstrEIG%lczsgm      = 0.0
        fstrEIG%lczmax      = 60
        fstrEIG%lcztol      = 1.0e-8
        fstrEIG%lczrod      = 1.0
        fstrEIG%lczrot      = 0.0
        fstrEIG%iluetol     = 0
end subroutine fstr_eigen_init

!> Initial setting of dynamic calculation
subroutine fstr_dynamic_init( fstrDYNAMIC )
        use m_fstr
        type(fstr_dynamic) :: fstrDYNAMIC
		
        fstrDYNAMIC%idx_eqa  = 1
        fstrDYNAMIC%idx_resp = 1
        fstrDYNAMIC%n_step   = 1
        fstrDYNAMIC%t_start  = 0.0
        fstrDYNAMIC%t_end    = 1.0
        fstrDYNAMIC%t_delta  = 1.0
        fstrDYNAMIC%ganma    = 0.5
        fstrDYNAMIC%beta     = 0.25
        fstrDYNAMIC%idx_mas  = 1
        fstrDYNAMIC%idx_dmp  = 1
        fstrDYNAMIC%ray_m    = 0.0
        fstrDYNAMIC%ray_k    = 0.0
        fstrDYNAMIC%restart_nout = 0
        fstrDYNAMIC%nout         = 100
        fstrDYNAMIC%node_monit_1 = 1
        fstrDYNAMIC%nout_monit   = 1
        fstrDYNAMIC%iout_list(1) = 0
        fstrDYNAMIC%iout_list(2) = 0
        fstrDYNAMIC%iout_list(3) = 0
        fstrDYNAMIC%iout_list(4) = 0
        fstrDYNAMIC%iout_list(5) = 0
        fstrDYNAMIC%iout_list(6) = 0
		
end subroutine fstr_dynamic_init


!> Initial setting of dynamic calculation
subroutine fstr_dynamic_alloc( fstrDYNAMIC )
        use m_fstr
        include "HEC_MW3_For.h"
        type(fstr_dynamic) :: fstrDYNAMIC
		
        integer :: ierror, ndof,nnod, iAss, iPart
		
		nnod = total_node 
        ndof=assDOF(1)
        if(fstrDYNAMIC%idx_eqa == 11) then
          allocate( fstrDYNAMIC%DISP(ndof*nnod,3)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%VEL (ndof*nnod,1)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%ACC (ndof*nnod,1)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        else
          allocate( fstrDYNAMIC%DISP(ndof*nnod,2)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%VEL (ndof*nnod,2)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%ACC (ndof*nnod,2)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        
        
        allocate( fstrDYNAMIC%VEC1(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC1>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate( fstrDYNAMIC%VEC2(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC2>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate( fstrDYNAMIC%VEC3(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC3>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

end subroutine fstr_dynamic_alloc

!> Finalizer of fstr_solid
 subroutine fstr_dynamic_finalize( fstrDYNAMIC )
        type(fstr_dynamic) :: fstrDYNAMIC
		
        integer :: ierror
        if( associated(fstrDYNAMIC%DISP) )	    &
        deallocate( fstrDYNAMIC%DISP    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        if( associated(fstrDYNAMIC%VEL) )	    &
        deallocate( fstrDYNAMIC%VEL     ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        if( associated(fstrDYNAMIC%ACC) )	    &
        deallocate( fstrDYNAMIC%ACC     ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        if( associated(fstrDYNAMIC%VEC1) )	    &
        deallocate( fstrDYNAMIC%VEC1    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEC1>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        if( associated(fstrDYNAMIC%VEC2) )	    &
        deallocate( fstrDYNAMIC%VEC2    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEC2>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        if( associated(fstrDYNAMIC%VEC3) )	    &
        deallocate( fstrDYNAMIC%VEC3    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, VEC3>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        
end subroutine		


!-----------------------------------------------------------------------------!
!> Initial setting of postprecessor

 subroutine fstr_setup_post( ctrl, P )
        implicit none                      
        integer(kind=kint) :: ctrl
        type(fstr_param_pack) :: P

                                           ! JP-6
        if( P%PARAM%solution_type == kstSTATIC &
         .or. P%PARAM%solution_type == kstNLSTATIC  &
         .or. P%PARAM%solution_type == kstEIGEN  &
         .or. P%PARAM%solution_type == kstDYNAMIC ) then
                ! Memory Allocation for Result Vectors ------------
                if( P%MESH%n_dof == 6 ) then
                        allocate ( P%SOLID%STRAIN  (10*total_node))
                        allocate ( P%SOLID%STRESS  (12*total_node))
                        allocate ( P%SOLID%ESTRAIN  (10*total_elem))
                        allocate ( P%SOLID%ESTRESS  (12*total_elem))
                else
                        allocate ( P%SOLID%STRAIN  (6*total_node))
                        allocate ( P%SOLID%STRESS  (7*total_node))
                        allocate ( P%SOLID%ESTRAIN  (6*total_elem))
                        allocate ( P%SOLID%ESTRESS  (7*total_elem))
                end if
                P%SOLID%STRAIN = 0.d0
                P%SOLID%STRESS = 0.d0
                P%SOLID%ESTRAIN = 0.d0
                P%SOLID%ESTRESS = 0.d0
        end if

        P%EIGEN%iluetol = svRarray(1) ! solver tolerance
        if( P%PARAM%solution_type == kstEIGEN ) then
            P%EIGEN%eqset = 1
        else
            P%EIGEN%eqset = 0
        endif

        if( P%PARAM%fg_visual == kON .and. P%MESH%my_rank == 0) then
     !           call fstr_setup_visualize( ctrl )
        end if

    !    call hecmw_barrier( P%MESH ) ! JP-7

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

!*****************************************************************************!
!* GENERAL HEADERS ***********************************************************!
!*****************************************************************************!

!-----------------------------------------------------------------------------!
!> Read in !SOLUTION                                                         !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_SOLUTION( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode

        rcode = fstr_ctrl_get_SOLUTION( ctrl, P%PARAM%solution_type )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

end subroutine fstr_setup_SOLUTION

!-----------------------------------------------------------------------------!
!> Read in !SOLVER                                                           !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_SOLVER( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack),target :: P

        integer(kind=kint) :: rcode

        if( counter >= 2 ) then
                write(ILOG,*) '### Error : !SOLVER exists twice in FSTR control file.'
                STOP
        endif

     !   nier       => svIarray(1)
     !   method     => svIarray(2)
     !   precond    => svIarray(3)
     !   nset       => svIarray(4)
     !   iterpremax => svIarray(5)
     !   nrest      => svIarray(6)
     !   iterlog    => svIarray(21)
     !   timelog    => svIarray(22)

     !   resid      => svRarray(1)
     !   sigma_diag => svRarray(2)
     !   sigma      => svRarray(3)
     !   thresh     => svRarray(4)
     !   filter     => svRarray(5)

        rcode = fstr_ctrl_get_SOLVER( ctrl,                      &
                        svIarray(2), svIarray(3), svIarray(4), svIarray(21), svIarray(22), &
                        svIarray(1), svIarray(5), svIarray(6),                 &
                        svRarray(1), svRarray(2), svRarray(3),                &
                        svRarray(4), svRarray(5) )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        if( svIarray(2) < 100 ) then
                svIarray(99) = 1  ! indirect method
        else
                svIarray(99) = 2  ! direct method
        end if

end subroutine fstr_setup_SOLVER


!-----------------------------------------------------------------------------!
!> Read in !STEP                                                             !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_STEP( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id

        integer(kind=kint) :: rcode, iproc

        amp = ' '
        rcode = fstr_ctrl_get_STEP( ctrl, amp, iproc, P%PARAM%incmax )
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        call amp_name_to_id( P%MESH, '!STEP', amp, amp_id )
    !    P%SOLID%NLSTATIC_ngrp_amp = amp_id;

end subroutine fstr_setup_STEP


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


!> Read in !RESTART
subroutine fstr_setup_RESTART( ctrl, n, fname )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: n
        character(len=HECMW_FILENAME_LEN) :: fname

        integer(kind=kint) :: rcode
        n = 0
        rcode = fstr_ctrl_get_param_ex( ctrl, 'FREQENCY ',  '# ',    0,   'I',   n  )
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        fname=""
        rcode = fstr_ctrl_get_param_ex( ctrl, 'NAME ',     '# ',  0, 'S', fname )
        if( fname=="" .or. rcode /= 0 ) stop "STOP:You must input the restart file name!"
end subroutine fstr_setup_RESTART


!-----------------------------------------------------------------------------!
!> Read in !COUPLE                                                          !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_COUPLE( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P
        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        integer(kind=kint) :: i, n, old_size, new_size

        if( P%SOLID%file_type /= kbcfFSTR ) return

        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return
        old_size = P%SOLID%COUPLE_ngrp_tot
        new_size = old_size + n
        P%SOLID%COUPLE_ngrp_tot = new_size

        call fstr_expand_integer_array ( P%SOLID%COUPLE_ngrp_ID,  old_size, new_size )

        allocate( grp_id_name(n))
        rcode = fstr_ctrl_get_COUPLE( ctrl,           &
                        P%PARAM%fg_couple_first,      &
                        grp_id_name, HECMW_NAME_LEN  )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        call surf_grp_name_to_id_ex( P%MESH, '!COUPLE', &
                        n, grp_id_name, P%SOLID%COUPLE_ngrp_ID(old_size+1:))

        deallocate( grp_id_name )
        P%PARAM%fg_couple = 1

end subroutine fstr_setup_COUPLE


!*****************************************************************************!
!* HEADERS FOR STATIC ANALYSIS ***********************************************!
!*****************************************************************************!

!-----------------------------------------------------------------------------!
!> Read in !STATIC(old)                                                           !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_STATIC( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P
        integer(kind=kint) :: rcode

        integer :: nout, nout_monit,node_monit_1 ,elem_monit_1 ,intg_monit_1  
        integer :: ipt, idx_elpl, iout_list(6) 
        real(kind=kreal) :: sig_y0, h_dash

        if( counter > 1 ) then
                write(*,*) 
        endif
		
        ipt = 0
        if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ', 'INFINITE,NLGEOM ', 0, 'P', ipt  )/=0 )  &
           return
        if( ipt == 2 ) P%PARAM%solution_type = kstNLSTATIC

        rcode = fstr_ctrl_get_STATIC( ctrl, &
                        DT, ETIME, ITMAX, EPS, P%SOLID%restart_nout, &
                        idx_elpl, &
                        iout_list, &
                        sig_y0, h_dash, &
                        nout, nout_monit, node_monit_1, &
                        elem_monit_1, intg_monit_1 )

        if( rcode /= 0 ) call fstr_ctrl_err_stop
        
end subroutine fstr_setup_STATIC

!-----------------------------------------------------------------------------!
!> Reset !DLOAD                                                        !
!-----------------------------------------------------------------------------!

subroutine fstr_expand_dload_array( array, old_size, new_size )
        implicit none
        real(kind=kreal), pointer :: array(:,:)
        integer(kind=kint) :: old_size, new_size, i, j
        real(kind=kreal), pointer :: temp(:,:)

        if( old_size >= new_size ) then
                return
        end if

        if( associated( array ) ) then
                allocate(temp(0:6, old_size))
                temp = array
                deallocate(array)
                allocate(array(0:6, new_size))
                array = 0
                do i=1,old_size
                        do j=0,6
                                array(j,i) = temp(j,i)
                        end do
                end do
                deallocate(temp)
        else
                allocate(array(0:6, new_size))
                array = 0
        end if
end subroutine fstr_expand_dload_array

!-----------------------------------------------------------------------------!
!> Read in !REFTEMP                                                          !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_REFTEMP( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode

        rcode = fstr_ctrl_get_REFTEMP( ctrl, P%PARAM%ref_temp )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

end subroutine fstr_setup_REFTEMP


!*****************************************************************************!
!* HEADERS FOR HEAT ANALYSIS *************************************************!
!*****************************************************************************!

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

        rcode = fstr_ctrl_get_HEAT(   ctrl,        &
                                P%PARAM%fg_irres,  &
                                P%PARAM%fg_iwres,  &
                                P%PARAM%dtime,     &
                                P%PARAM%etime,     &
                                P%PARAM%dtmin,     &
                                P%PARAM%delmax,    &
                                P%PARAM%itmax,     &
                                P%PARAM%eps )
        if( rcode /= 0 ) then
                call fstr_ctrl_err_stop
        end if

        call reallocate_real( P%HEAT%STEP_DLTIME, n)
        call reallocate_real( P%HEAT%STEP_EETIME, n)
        call reallocate_real( P%HEAT%STEP_DELMIN, n)
        call reallocate_real( P%HEAT%STEP_DELMAX, n)
        P%HEAT%STEPtot = n

        P%HEAT%STEP_DLTIME = P%PARAM%dtime
        P%HEAT%STEP_EETIME = P%PARAM%etime
        P%HEAT%STEP_DELMIN = P%PARAM%dtmin
        P%HEAT%STEP_DELMAX = P%PARAM%delmax

end subroutine fstr_setup_HEAT

!-----------------------------------------------------------------------------!
!> Read in !FIXTEMP                                                          !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_FIXTEMP( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack),target :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member(:)
        integer(kind=kint) :: local_id, rtc
        ! ------------------------------------------------

        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return

        allocate( grp_id_name(n))
        allocate( value(n))

        amp = ' '
        rcode = fstr_ctrl_get_FIXTEMP( ctrl, amp, &
                        grp_id_name, HECMW_NAME_LEN, value )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        call amp_name_to_id( P%MESH, '!FIXTEMP', amp, amp_id ) 

        m = 0
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'node', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        m = m + 1
                else if( rtc < 0 ) then
                        m = m + get_grp_member_n( P%MESH, 'node_grp', grp_id_name(i) )
                end if
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( value )
                return
        endif

        ! JP-8
        old_size = P%HEAT%T_FIX_tot
        new_size = old_size + m
        call fstr_expand_integer_array( P%HEAT%T_FIX_node, old_size, new_size )
        call fstr_expand_integer_array( P%HEAT%T_FIX_ampl, old_size, new_size )
        call fstr_expand_real_array(    P%HEAT%T_FIX_val,  old_size, new_size )
        P%HEAT%T_FIX_tot = new_size

        head = old_size + 1
        member => P%HEAT%T_FIX_node(head:)
        id = head
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'node', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        member(1) = local_id
                        member_n = 1
                else if( rtc < 0 ) then
                        member_n = get_grp_member( P%MESH, 'node_grp', grp_id_name(i), member )
                else
                        cycle
                end if
                if( i<n ) then
                        member => member( member_n+1 : )
                endif
                do j = 1, member_n
                        P%HEAT%T_FIX_val  (id) = value(i)
                        P%HEAT%T_FIX_ampl (id) = amp_id
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( value )
end subroutine fstr_setup_FIXTEMP


!-----------------------------------------------------------------------------!
!> Read in !CFLUX                                                            !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_CFLUX( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member(:)
        integer(kind=kint) :: local_id, rtc
        ! ------------------------------------------------

        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return

        allocate( grp_id_name(n))
        allocate( value(n))

        amp = ' '
        rcode = fstr_ctrl_get_CFLUX( ctrl, amp, &
                        grp_id_name, HECMW_NAME_LEN, value )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        call amp_name_to_id( P%MESH, '!CFLUX', amp, amp_id ) 

        m = 0

        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'node', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        m = m + 1
                else if( rtc < 0 ) then
                        m = m + get_grp_member_n( P%MESH, 'node_grp', grp_id_name(i) )
                end if
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( value )
                return
        endif

        ! JP-9
        old_size = P%HEAT%Q_NOD_tot
        new_size = old_size + m
        call fstr_expand_integer_array( P%HEAT%Q_NOD_node, old_size, new_size )
        call fstr_expand_integer_array( P%HEAT%Q_NOD_ampl, old_size, new_size )
        call fstr_expand_real_array(    P%HEAT%Q_NOD_val,  old_size, new_size )
        P%HEAT%Q_NOD_tot = new_size

        head = old_size + 1
        member => P%HEAT%Q_NOD_node(head:)
        id = head
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'node', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        member(1) = local_id
                        member_n = 1
                else if( rtc < 0 ) then
                        member_n = get_grp_member( P%MESH, 'node_grp', grp_id_name(i), member )
                else
                        cycle
                end if
                if( i<n ) member => member( member_n+1 : )
                do j = 1, member_n
                        P%HEAT%Q_NOD_val  (id) = value(i)
                        P%HEAT%Q_NOD_ampl (id) = amp_id
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( value )
end subroutine fstr_setup_CFLUX


!-----------------------------------------------------------------------------!
!> Read in !DFLUX                                                            !
!-----------------------------------------------------------------------------!


subroutine fstr_setup_DFLUX( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        integer(kind=kint),pointer :: load_type(:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member(:)
        integer(kind=kint) :: local_id, rtc
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

        call amp_name_to_id( P%MESH, '!DFLUX', amp, amp_id ) 

        m = 0
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'element', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        m = m + 1
                else if( rtc < 0 ) then
                        m = m + get_grp_member_n( P%MESH, 'elem_grp', grp_id_name(i) )
                end if
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( load_type )
                deallocate( value )
                return
        endif

        ! JP-10
        old_size = P%HEAT%Q_SUF_tot
        new_size = old_size + m
        call fstr_expand_integer_array( P%HEAT%Q_SUF_elem, old_size, new_size )
        call fstr_expand_integer_array( P%HEAT%Q_SUF_ampl, old_size, new_size )
        call fstr_expand_integer_array( P%HEAT%Q_SUF_surf, old_size, new_size )
        call fstr_expand_real_array(    P%HEAT%Q_SUF_val,  old_size, new_size )
        P%HEAT%Q_SUF_tot = new_size

        head = old_size + 1
        member => P%HEAT%Q_SUF_elem(head:)
        id = head
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'element', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        member(1) = local_id
                        member_n = 1
                else if( rtc < 0 ) then
                        member_n = get_grp_member( P%MESH, 'elem_grp', grp_id_name(i), member )
                else
                        cycle
                end if
                if( i<n ) member => member( member_n+1 : )
                do j = 1, member_n
                        P%HEAT%Q_SUF_surf (id) = load_type(i)
                        P%HEAT%Q_SUF_val  (id) = value(i)
                        P%HEAT%Q_SUF_ampl (id) = amp_id
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( load_type )
        deallocate( value )
end subroutine fstr_setup_DFLUX


!-----------------------------------------------------------------------------!
!> Read in !SFLUX                                                            !
!-----------------------------------------------------------------------------!


subroutine fstr_setup_SFLUX( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: value(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member1(:), member2(:)
        ! ------------------------------------------------

        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return

        allocate( grp_id_name(n))
        allocate( value(n))

        amp = ' '
        rcode = fstr_ctrl_get_SFLUX( ctrl, amp, &
                        grp_id_name, HECMW_NAME_LEN, value )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        call amp_name_to_id( P%MESH, '!SFLUX', amp, amp_id ) 

        m = 0
        do i = 1, n
                m = m + get_grp_member_n( P%MESH, 'surf_grp', grp_id_name(i) )
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( value )
                return
        endif

        ! JP-11
        old_size = P%HEAT%Q_SUF_tot
        new_size = old_size + m
        call fstr_expand_integer_array( P%HEAT%Q_SUF_elem, old_size, new_size )
        call fstr_expand_integer_array( P%HEAT%Q_SUF_ampl, old_size, new_size )
        call fstr_expand_integer_array( P%HEAT%Q_SUF_surf, old_size, new_size )
        call fstr_expand_real_array(    P%HEAT%Q_SUF_val,  old_size, new_size )
        P%HEAT%Q_SUF_tot = new_size

        head = old_size + 1
        member1 => P%HEAT%Q_SUF_elem(head:)
        member2 => P%HEAT%Q_SUF_surf(head:)
        id = head
        do i = 1, n
                member_n = get_grp_member( P%MESH, 'surf_grp', grp_id_name(i), member1, member2 )
                if( i<n ) then
                        member1 => member1( member_n+1 : )
                        member2 => member2( member_n+1 : )
                end if
                do j = 1, member_n
                        P%HEAT%Q_SUF_val  (id) = value(i)
                        P%HEAT%Q_SUF_ampl (id) = amp_id
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( value )
end subroutine fstr_setup_SFLUX


!-----------------------------------------------------------------------------!
!> Read in !FILM                                                             !
!-----------------------------------------------------------------------------!


subroutine fstr_setup_FILM( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp1, amp2
        integer(kind=kint) :: amp_id1, amp_id2
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        integer(kind=kint),pointer :: load_type(:)
        real(kind=kreal),pointer :: value(:)
        real(kind=kreal),pointer :: shink(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member(:)
        integer(kind=kint) :: local_id, rtc
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

        call amp_name_to_id( P%MESH, '!FILM', amp1, amp_id1 ) 
        call amp_name_to_id( P%MESH, '!FILM', amp2, amp_id2 ) 

        m = 0
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'element', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        m = m + 1
                else if( rtc < 0 ) then
                        m = m + get_grp_member_n( P%MESH, 'elem_grp', grp_id_name(i) )
                end if
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( load_type )
                deallocate( value )
                deallocate( shink )
                return
        endif

        ! JP-12
        old_size = P%HEAT%H_SUF_tot
        new_size = old_size + m
        call fstr_expand_integer_array(  P%HEAT%H_SUF_elem,    old_size, new_size )
        call fstr_expand_integer_array2( P%HEAT%H_SUF_ampl, 2, old_size, new_size )
        call fstr_expand_integer_array(  P%HEAT%H_SUF_surf,    old_size, new_size )
        call fstr_expand_real_array2(    P%HEAT%H_SUF_val,  2, old_size, new_size )
        P%HEAT%H_SUF_tot = new_size

        head = old_size + 1
        member => P%HEAT%H_SUF_elem(head:)
        id = head
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'element', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        member(1) = local_id
                        member_n = 1
                else if( rtc < 0 ) then
                        member_n = get_grp_member( P%MESH, 'elem_grp', grp_id_name(i), member )
                else
                        cycle
                end if
                if( i<n ) member => member( member_n+1 : )
                do j = 1, member_n
                        P%HEAT%H_SUF_surf (id)   = load_type(i)
                        P%HEAT%H_SUF_val  (id,1) = value(i)
                        P%HEAT%H_SUF_val  (id,2) = shink(i)
                        P%HEAT%H_SUF_ampl (id,1) = amp_id1
                        P%HEAT%H_SUF_ampl (id,2) = amp_id2
                        id= id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( load_type )
        deallocate( value )
        deallocate( shink )
end subroutine fstr_setup_FILM


!-----------------------------------------------------------------------------!
!> Read in !SFILM                                                            !
!-----------------------------------------------------------------------------!


subroutine fstr_setup_SFILM( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp1, amp2
        integer(kind=kint) :: amp_id1, amp_id2
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: value(:)
        real(kind=kreal),pointer :: shink(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member1(:), member2(:)
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

        call amp_name_to_id( P%MESH, '!SFILM', amp1, amp_id1 ) 
        call amp_name_to_id( P%MESH, '!SFILM', amp2, amp_id2 ) 

        m = 0
        do i = 1, n
                m = m + get_grp_member_n( P%MESH, 'surf_grp', grp_id_name(i) )
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( value )
                deallocate( shink )
                return
        endif

        ! JP-13
        old_size = P%HEAT%H_SUF_tot
        new_size = old_size + m
        call fstr_expand_integer_array(  P%HEAT%H_SUF_elem,    old_size, new_size )
        call fstr_expand_integer_array2( P%HEAT%H_SUF_ampl, 2, old_size, new_size )
        call fstr_expand_integer_array(  P%HEAT%H_SUF_surf,    old_size, new_size )
        call fstr_expand_real_array2(    P%HEAT%H_SUF_val,  2, old_size, new_size )
        P%HEAT%H_SUF_tot = new_size

        head = old_size + 1
        member1 => P%HEAT%H_SUF_elem(head:)
        member2 => P%HEAT%H_SUF_surf(head:)
        id = head
        do i = 1, n
                member_n = get_grp_member( P%MESH, 'surf_grp', grp_id_name(i), member1, member2 )
                if( i<n ) then
                        member1 => member1( member_n+1 : )
                        member2 => member2( member_n+1 : )
                end if
                do j = 1, member_n
                        P%HEAT%H_SUF_val  (id,1) = value(i)
                        P%HEAT%H_SUF_val  (id,2) = shink(i)
                        P%HEAT%H_SUF_ampl (id,1) = amp_id1
                        P%HEAT%H_SUF_ampl (id,2) = amp_id2
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( value )
        deallocate( shink )
end subroutine fstr_setup_SFILM


!-----------------------------------------------------------------------------!
!> Read in !RADIATE                                                          !
!-----------------------------------------------------------------------------!


subroutine fstr_setup_RADIATE( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp1, amp2
        integer(kind=kint) :: amp_id1, amp_id2
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        integer(kind=kint),pointer :: load_type(:)
        real(kind=kreal),pointer :: value(:)
        real(kind=kreal),pointer :: shink(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member(:)
        integer(kind=kint) :: local_id, rtc
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

        call amp_name_to_id( P%MESH, '!RADIATE', amp1, amp_id1 ) 
        call amp_name_to_id( P%MESH, '!RADIATE', amp2, amp_id2 ) 

        m = 0
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'element', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        m = m + 1
                else if( rtc < 0 ) then
                        m = m + get_grp_member_n( P%MESH, 'elem_grp', grp_id_name(i) )
                end if
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( load_type )
                deallocate( value )
                deallocate( shink )
                return
        endif

        ! JP-14
        old_size = P%HEAT%R_SUF_tot
        new_size = old_size + m
        call fstr_expand_integer_array(  P%HEAT%R_SUF_elem,    old_size, new_size )
        call fstr_expand_integer_array2( P%HEAT%R_SUF_ampl, 2, old_size, new_size )
        call fstr_expand_integer_array(  P%HEAT%R_SUF_surf,    old_size, new_size )
        call fstr_expand_real_array2(    P%HEAT%R_SUF_val,  2, old_size, new_size )
        P%HEAT%R_SUF_tot = new_size

        head = old_size + 1
        member => P%HEAT%R_SUF_elem(head:)
        id = head
        do i = 1, n
                rtc = get_local_member_index( P%MESH, 'element', grp_id_name(i), local_id )
                if( rtc > 0 ) then
                        member(1) = local_id
                        member_n = 1
                else if( rtc < 0 ) then
                        member_n = get_grp_member( P%MESH, 'elem_grp', grp_id_name(i), member )
                else
                        cycle
                end if
                if( i<n ) member => member( member_n+1 : )
                do j = 1, member_n
                        P%HEAT%R_SUF_surf (id)   = load_type(i)
                        P%HEAT%R_SUF_val  (id,1) = value(i)
                        P%HEAT%R_SUF_val  (id,2) = shink(i)
                        P%HEAT%R_SUF_ampl (id,1) = amp_id1
                        P%HEAT%R_SUF_ampl (id,2) = amp_id2
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( load_type )
        deallocate( value )
        deallocate( shink )
end subroutine fstr_setup_RADIATE


!-----------------------------------------------------------------------------!
!> Read in !SRADIATE                                                         !
!-----------------------------------------------------------------------------!


subroutine fstr_setup_SRADIATE( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp1, amp2
        integer(kind=kint) :: amp_id1, amp_id2
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: value(:)
        real(kind=kreal),pointer :: shink(:)
        integer(kind=kint) :: i, j, n, m, head, id, member_n, old_size, new_size
        integer(kind=kint),pointer :: member1(:), member2(:)
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

        call amp_name_to_id( P%MESH, '!SRADIATE', amp1, amp_id1 ) 
        call amp_name_to_id( P%MESH, '!SRADIATE', amp2, amp_id2 ) 

        m = 0
        do i = 1, n
                m = m + get_grp_member_n( P%MESH, 'surf_grp', grp_id_name(i) )
        end do

        if (m == 0) then
                deallocate( grp_id_name )
                deallocate( value )
                deallocate( shink )
                return
        endif

        ! JP-15
        old_size = P%HEAT%R_SUF_tot
        new_size = old_size + m
        call fstr_expand_integer_array(  P%HEAT%R_SUF_elem,    old_size, new_size )
        call fstr_expand_integer_array2( P%HEAT%R_SUF_ampl, 2, old_size, new_size )
        call fstr_expand_integer_array(  P%HEAT%R_SUF_surf,    old_size, new_size )
        call fstr_expand_real_array2(    P%HEAT%R_SUF_val,  2, old_size, new_size )
        P%HEAT%R_SUF_tot = new_size

        head = old_size + 1
        member1 => P%HEAT%R_SUF_elem(head:)
        member2 => P%HEAT%R_SUF_surf(head:)
        id = head
        do i = 1, n
                member_n = get_grp_member( P%MESH, 'surf_grp', grp_id_name(i), member1, member2 )
                if( i<n ) then
                        member1 => member1( member_n+1 : )
                        member2 => member2( member_n+1 : )
                end if
                do j = 1, member_n
                        P%HEAT%R_SUF_val  (id,1) = value(i)
                        P%HEAT%R_SUF_val  (id,2) = shink(i)
                        P%HEAT%R_SUF_ampl (id,1) = amp_id1
                        P%HEAT%R_SUF_ampl (id,2) = amp_id2
                        id = id + 1
                end do
        end do

        deallocate( grp_id_name )
        deallocate( value )
        deallocate( shink )
end subroutine fstr_setup_SRADIATE


!*****************************************************************************!
!* HEADERS FOR EIGEN ANALYSIS ************************************************!
!*****************************************************************************!

!-----------------------------------------------------------------------------!
!> Read in !EIGEN                                                            !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_EIGEN( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode

        rcode = fstr_ctrl_get_EIGEN( ctrl, P%EIGEN%nget, P%EIGEN%lcztol, P%EIGEN%lczmax)
        if( rcode /= 0) call fstr_ctrl_err_stop

end subroutine fstr_setup_EIGEN


!*****************************************************************************!
!* HEADERS FOR DYNAMIC ANALYSIS **********************************************!
!*****************************************************************************!

!-----------------------------------------------------------------------------!
!> Read in !DYNAMIC                                                          !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_DYNAMIC( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P
        integer(kind=kint) :: rcode

        rcode = fstr_ctrl_get_DYNAMIC( ctrl, &
                P%DYN%idx_eqa, &
                P%DYN%idx_resp,&
                P%DYN%n_step,  &
                P%DYN%t_start, &
                P%DYN%t_end,   &
                P%DYN%t_delta, &
                P%DYN%restart_nout, &
                P%DYN%ganma,   &
                P%DYN%beta,    &
                P%DYN%idx_mas, &
                P%DYN%idx_dmp, &
                P%DYN%ray_m,   &
                P%DYN%ray_k,   &
                P%DYN%nout,    &
                P%DYN%node_monit_1,&
                P%DYN%nout_monit,  &
                P%DYN%iout_list )

        if( rcode /= 0) call fstr_ctrl_err_stop

end subroutine fstr_setup_DYNAMIC


!-----------------------------------------------------------------------------!
!> Read in !VELOCITY                                                         !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_VELOCITY( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        integer(kind=kint),pointer :: dof_ids (:)
        integer(kind=kint),pointer :: dof_ide (:)
        real(kind=kreal),pointer :: val_ptr(:)
        integer(kind=kint) :: i, j, n, old_size, new_size
                
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return
        old_size = P%SOLID%VELOCITY_ngrp_tot
        new_size = old_size + n
        P%SOLID%VELOCITY_ngrp_tot = new_size

        call fstr_expand_integer_array (P%SOLID%VELOCITY_ngrp_ID  , old_size, new_size )
        call fstr_expand_integer_array (P%SOLID%VELOCITY_ngrp_type, old_size, new_size )
        call fstr_expand_real_array    (P%SOLID%VELOCITY_ngrp_val , old_size, new_size )
        call fstr_expand_integer_array (P%SOLID%VELOCITY_ngrp_amp , old_size, new_size )

        allocate( grp_id_name(n))
        allocate( dof_ids (n))
        allocate( dof_ide (n))

        amp = ''
        val_ptr => P%SOLID%VELOCITY_ngrp_val(old_size+1:)
        val_ptr = 0
        rcode = fstr_ctrl_get_VELOCITY( ctrl,  amp,   &
                grp_id_name, HECMW_NAME_LEN,  &
                dof_ids, dof_ide, val_ptr )
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        call amp_name_to_id( P%MESH, '!VELOCITY', amp, amp_id ) 
        call node_grp_name_to_id_ex( P%MESH, '!VELOCITY', &
                n, grp_id_name, P%SOLID%VELOCITY_ngrp_ID(old_size+1:))

        j = old_size+1
        do i = 1, n
                if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < 1).or.(6 < dof_ide(i)) ) then
                         write(ILOG,*) 'fstr contol file error : !VELOCITY : range of dof_ids and dof_ide is from 1 to 6'
                              STOP
                end if
                P%SOLID%VELOCITY_ngrp_type(j) = 10 * dof_ids(i) + dof_ide(i)
                P%SOLID%VELOCITY_ngrp_amp(j) = amp_id
                j = j+1
        end do

        deallocate( grp_id_name )
        deallocate( dof_ids )
        deallocate( dof_ide )

end subroutine fstr_setup_VELOCITY


!-----------------------------------------------------------------------------!
!> Read in !ACCELERATION                                                     !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_ACCELERATION( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        integer(kind=kint),pointer :: dof_ids (:)
        integer(kind=kint),pointer :: dof_ide (:)
        real(kind=kreal),pointer :: val_ptr(:)
        integer(kind=kint) :: i, j, n, old_size, new_size

                
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return
        old_size = P%SOLID%ACCELERATION_ngrp_tot
        new_size = old_size + n
        P%SOLID%ACCELERATION_ngrp_tot = new_size

        call fstr_expand_integer_array (P%SOLID%ACCELERATION_ngrp_ID  , old_size, new_size )
        call fstr_expand_integer_array (P%SOLID%ACCELERATION_ngrp_type, old_size, new_size )
        call fstr_expand_real_array    (P%SOLID%ACCELERATION_ngrp_val , old_size, new_size )
        call fstr_expand_integer_array (P%SOLID%ACCELERATION_ngrp_amp , old_size, new_size )

        allocate( grp_id_name(n))
        allocate( dof_ids (n))
        allocate( dof_ide (n))

        amp = ' '
        val_ptr => P%SOLID%ACCELERATION_ngrp_val(old_size+1:)
        val_ptr = 0
        rcode = fstr_ctrl_get_ACCELERATION( ctrl,  amp,   &
                grp_id_name, HECMW_NAME_LEN,  &
                dof_ids, dof_ide,  val_ptr)
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        call amp_name_to_id( P%MESH, '!ACCELERATION', amp, amp_id ) 
        call node_grp_name_to_id_ex( P%MESH, '!ACCELERATION', &
                n, grp_id_name, P%SOLID%ACCELERATION_ngrp_ID(old_size+1:))

        j = old_size+1
        do i = 1, n
                if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < 1).or.(6 < dof_ide(i)) ) then
                         write(ILOG,*) 'fstr contol file error : !ACCELERATION : range of dof_ids and dof_ide is from 1 to 6'
                              STOP
                end if
                P%SOLID%ACCELERATION_ngrp_type(j) = 10 * dof_ids(i) + dof_ide(i)
                P%SOLID%ACCELERATION_ngrp_amp(j) = amp_id
                j = j+1
        end do

        deallocate( grp_id_name )
        deallocate( dof_ids )
        deallocate( dof_ide )
end subroutine fstr_setup_ACCELERATION


!-----------------------------------------------------------------------------!
!> Read in !MPC                                                              !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_MPC( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack), target :: P

        integer(kind=kint) :: rcode

        rcode = fstr_ctrl_get_MPC( ctrl, svRarray(11))
        if( rcode /= 0) call fstr_ctrl_err_stop
end subroutine fstr_setup_MPC


subroutine hecMW_init_section( matls )
        type(tMaterial)    :: matls(:)
        integer :: i, j, nsect
        logical :: isfind
        if( .not. associated(MWSections) ) stop "SECTION not defined!"
        nsect = size(MWSections)
        do i=1, nsect
          isfind = .false.
          MWSections(i)%sect_R_item = 0.d0
          do j=1, size(matls)
            if( MWSections(i)%mat_name == matls(j)%name ) then
               MWSections(i)%sect_mat_ID = j
               isfind = .true.
            endif
          enddo
          if( .not. isfind ) stop "Error in SECTION definition: Material not found!"
        enddo
end subroutine hecMW_init_section


end module m_fstr_setup


