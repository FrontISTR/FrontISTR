!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : I/O and Utility                                   !
!                                                                      !
!            Written by Noboru Imai (Univ. of Tokyo)                   !
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
use mContact
use m_static_get_prop
use m_out
use m_step
use m_utilities
use m_fstr_freqdata, only : fstr_freqanalysis, kFLOADCASE_RE, kFLOADCASE_IM
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
                type(fstr_freqanalysis), pointer  :: FREQ
        end type fstr_param_pack

contains

!=============================================================================!
!> Read in and initialize control data                                        !
!=============================================================================!
subroutine fstr_setup( cntl_filename, hecMESH, fstrPARAM,  &
        fstrSOLID, fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL, fstrFREQ )
        use mMaterial
        character(len=HECMW_FILENAME_LEN) :: cntl_filename
        type(hecmwST_local_mesh),target :: hecMESH
        type(fstr_param),target   :: fstrPARAM
        type(fstr_solid),target   :: fstrSOLID
        type(lczparam),target     :: fstrEIG
        type(fstr_heat),target    :: fstrHEAT
        type(fstr_dynamic),target :: fstrDYNAMIC
        type(fstr_couple),target  :: fstrCPL
        type(fstr_freqanalysis), target :: fstrFREQ

        integer(kind=kint) :: ctrl
        type(fstr_param_pack) :: P

        integer, parameter :: MAXOUTFILE = 10

        external fstr_ctrl_get_c_h_name
        integer(kind=kint) :: fstr_ctrl_get_c_h_name

        integer(kind=kint) :: version, resul, visual, femap, n_totallayer_ls, shell_matltype, aa, bb, mixedflag
        integer(kind=kint) :: rcode, n, i, j, it, cid, nout, nin, ierror                     
        character(len=HECMW_NAME_LEN) :: header_name, fname(MAXOUTFILE)
        real(kind=kreal) :: ee, pp, rho, alpha, thick, alpha_over_mu, shell_variables(200)
        real(kind=kreal) :: beam_radius,                          &
                            beam_angle1, beam_angle2, beam_angle3,&
                            beam_angle4, beam_angle5, beam_angle6 
        logical          :: isOK
        type(t_output_ctrl) :: outctrl
        character(len=HECMW_FILENAME_LEN) :: logfileNAME, mName

        ! counters
        integer(kind=kint) :: c_solution, c_solver, c_step, c_write, c_echo
        integer(kind=kint) :: c_static, c_boundary, c_cload, c_dload, c_temperature, c_reftemp, c_spring
        integer(kind=kint) :: c_heat, c_fixtemp, c_cflux, c_dflux, c_sflux, c_film, c_sfilm, c_radiate, c_sradiate
        integer(kind=kint) :: c_eigen, c_contact
        integer(kind=kint) :: c_dynamic, c_velocity, c_acceleration
        integer(kind=kint) :: c_fload, c_eigenread
        integer(kind=kint) :: c_couple, c_material
        integer(kind=kint) :: c_mpc, c_weldline
        integer(kind=kint) :: c_istep, c_localcoord, c_section
        integer(kind=kint) :: c_output, islog
        integer(kind=kint) :: k

        write( logfileNAME, '(i5,''.log'')' ) myrank

        ! packaging
        P%MESH   => hecMESH
        P%PARAM  => fstrPARAM
        P%SOLID  => fstrSOLID
        P%EIGEN  => fstrEIG
        P%HEAT   => fstrHEAT
        P%DYN    => fstrDYNAMIC
        P%CPL    => fstrCPL
        P%FREQ   => fstrFREQ

        fstrPARAM%contact_algo = kcaALagrange

        c_solution = 0; c_solver   = 0; c_step   = 0; c_output = 0; c_echo = 0;
        c_static   = 0; c_boundary = 0; c_cload  = 0; c_dload = 0; c_temperature = 0; c_reftemp = 0; c_spring = 0;
        c_heat     = 0; c_fixtemp  = 0; c_cflux  = 0; c_dflux = 0; c_sflux = 0
        c_film     = 0; c_sfilm    = 0; c_radiate= 0; c_sradiate = 0
        c_eigen    = 0; c_contact  = 0
        c_dynamic  = 0; c_velocity = 0; c_acceleration = 0
        c_couple   = 0; c_material = 0; c_section =0
        c_mpc      = 0; c_weldline = 0
        c_istep    = 0; c_localcoord = 0
        c_fload    = 0; c_eigenread = 0

        ctrl = fstr_ctrl_open( cntl_filename)
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
                else if( header_name == '!WRITE' ) then
                        call fstr_ctrl_get_output( ctrl, outctrl, islog, resul, visual, femap )
                        if( visual==1 ) P%PARAM%fg_visual= 1
                        if( resul==1 ) P%PARAM%fg_result = 1
                        c_output = c_output+1
                else if( header_name == '!ECHO' ) then
                        c_echo = c_echo + 1
                        call fstr_setup_ECHO( ctrl, c_echo, P )
                else if( header_name == '!RESTART' ) then
                        call fstr_setup_RESTART( ctrl, nout )
                        fstrSOLID%restart_nout= nout
                        fstrDYNAMIC%restart_nout= nout
                        fstrHEAT%restart_nout= nout
                else if( header_name == '!ORIENTATION' ) then
                        c_localcoord = c_localcoord + 1

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

                !--------------- for dynamic -------------------------
                else if( header_name == '!FLOAD' ) then
                        c_fload = c_fload + 1
                        call fstr_setup_FLOAD( ctrl , c_fload, P )
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

                !--------------------- END -------------------------

                else if( header_name == '!END' ) then
                        exit
                end if

                ! next
                if( fstr_ctrl_seek_next_header(ctrl) == 0) exit
        end do

! ----- 
        if( c_contact>0 ) allocate( fstrSOLID%contacts( c_contact ) )
        if( c_weldline>0 ) allocate( fstrHEAT%weldline( c_weldline ) )	
        if( c_istep>0 ) allocate( fstrSOLID%step_ctrl( c_istep ) )
        if( c_localcoord>0 ) allocate( g_LocalCoordSys(c_localcoord) )

        n = c_material
        if( hecMESH%material%n_mat>n ) n= hecMESH%material%n_mat
        if( n==0 ) stop "material property not defined!"
        allocate( fstrSOLID%materials( n ) )
        if( hecMESH%section%n_sect >0 ) then
           do i=1,hecMESH%section%n_sect
              if( hecMESH%section%sect_type(i) == 4 ) cycle
              cid = hecMESH%section%sect_mat_ID_item(i)
              if( cid>n ) stop "Error in material property definition!"
              call initMaterial(  fstrSOLID%materials(cid) )
              if( fstrPARAM%solution_type == kstNLSTATIC &
                   .or. fstrPARAM%solution_type==6 ) &
                      fstrSOLID%materials(cid)%nlgeom_flag = 1
              call fstr_get_prop(hecMESH,i,ee,pp,rho,alpha,thick,&
                            alpha_over_mu,n_totallayer_ls,&
                            shell_matltype,shell_variables, &
                            beam_radius,beam_angle1,beam_angle2,beam_angle3,   &
                            beam_angle4,beam_angle5,beam_angle6)
              fstrSOLID%materials(cid)%name = hecMESH%material%mat_name(cid)
              fstrSOLID%materials(cid)%variables(M_YOUNGS)=ee
              fstrSOLID%materials(cid)%variables(M_POISSON)=pp
              fstrSOLID%materials(cid)%variables(M_DENSITY)=rho
              fstrSOLID%materials(cid)%variables(M_EXAPNSION)=alpha
              fstrSOLID%materials(cid)%variables(M_THICK)=thick
              fstrSOLID%materials(cid)%variables(M_ALPHA_OVER_MU)= alpha_over_mu
              fstrSOLID%materials(cid)%variables(M_BEAM_RADIUS)=beam_radius
              fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE1)=beam_angle1
              fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE2)=beam_angle2
              fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE3)=beam_angle3
              fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE4)=beam_angle4
              fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE5)=beam_angle5
              fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE6)=beam_angle6
              fstrSOLID%materials(cid)%mtype = ELASTIC
              !<**********  Laminated shell  **********
              fstrSOLID%materials(cid)%variables(M_TOTAL_LAYER)=n_totallayer_ls*1.0D0
              fstrSOLID%materials(cid)%variables(M_SHELL_MATLTYPE)=shell_matltype*1.0D0
              if (shell_matltype == 0)then
                do k=1,n_totallayer_ls
                  fstrSOLID%materials(cid)%variables(100+3*k-2)=shell_variables(3*k-2)		!<ee
                  fstrSOLID%materials(cid)%variables(100+3*k-1)=shell_variables(3*k-1)		!<pp
                  fstrSOLID%materials(cid)%variables(100+3*k  )=shell_variables(3*k  )		!<thickness
                end do
              elseif (shell_matltype == 1)then
                do k=1,n_totallayer_ls
                  fstrSOLID%materials(cid)%variables(100+8*k-7)=shell_variables(8*k-7)		!<ee
                  fstrSOLID%materials(cid)%variables(100+8*k-6)=shell_variables(8*k-6)		!<pp
                  fstrSOLID%materials(cid)%variables(100+8*k-5)=shell_variables(8*k-5)		!<thickness
                  fstrSOLID%materials(cid)%variables(100+8*k-4)=shell_variables(8*k-4)		!<ee2
                  fstrSOLID%materials(cid)%variables(100+8*k-3)=shell_variables(8*k-3)		!<g12
                  fstrSOLID%materials(cid)%variables(100+8*k-2)=shell_variables(8*k-2)		!<g23
                  fstrSOLID%materials(cid)%variables(100+8*k-1)=shell_variables(8*k-1)		!<g31
                  fstrSOLID%materials(cid)%variables(100+8*k  )=(shell_variables(8*k )/180.0d0)*datan(1.0d0)*4d0		!<theta
                enddo
              endif
              !>**********  Laminated shell  **********
           enddo
         endif

         allocate( fstrSOLID%output_ctrl( 4 ) )
         call fstr_init_outctrl(fstrSOLID%output_ctrl(1))
         fstrSOLID%output_ctrl( 1 )%filename = trim(logfileNAME)
         fstrSOLID%output_ctrl( 1 )%filenum = ILOG
         call fstr_init_outctrl(fstrSOLID%output_ctrl(2))
         call fstr_init_outctrl(fstrSOLID%output_ctrl(3))
         call fstr_init_outctrl(fstrSOLID%output_ctrl(4))

! ----- 
        rcode = fstr_ctrl_rewind( ctrl )

        c_istep    = 0
        c_material = 0
        c_output = 0
        c_contact  = 0
        c_localcoord = 0
        c_section = 0
        fstrHEAT%WL_tot = 0
        do
          rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )

          if( header_name == '!ORIENTATION' ) then
                        c_localcoord = c_localcoord + 1
                        if( fstr_setup_ORIENTATION( ctrl, hecMESH, c_localcoord, g_LocalCoordSys(c_localcoord) )/=0 ) then
                          write(*,*) '### Error: Fail in read in ORIENTATION definition : ', c_localcoord
                          write(ILOG,*) '### Error: Fail in read in ORIENTATION definition : ', c_localcoord
                          stop
                        endif

          ! ----- CONTACT condtion setting
          elseif( header_name == '!CONTACT' ) then                          
            n = fstr_ctrl_get_data_line_n( ctrl )
            if( .not. fstr_ctrl_get_CONTACT( ctrl, n, fstrSOLID%contacts(c_contact+1:c_contact+n)   &    
                ,ee, pp, rho, alpha, P%PARAM%contact_algo ) ) then                                    
                write(*,*) '### Error: Fail in read in contact condition : ', c_contact
                write(ILOG,*) '### Error: Fail in read in contact condition : ', c_contact
                STOP
            endif
            ! initialize contact condition
            if( ee>0.d0 ) cdotp = ee
            if( pp>0.d0 ) mut = pp
            if( rho>0.d0 ) cgn = rho
            if( alpha>0.d0 ) cgt = alpha
            do i=1,n
              if( .not. fstr_contact_check( fstrSOLID%contacts(c_contact+i), P%MESH ) ) then
                write(*,*) '### Error: Inconsistence in contact and surface definition : ' , i+c_contact
                write(ILOG,*) '### Error: Inconsistence in contact and surface definition : ', i+c_contact
                stop
              else
                if(paraContactFlag) then
                  isOK = fstr_contact_init( fstrSOLID%contacts(c_contact+i), P%MESH   &
!#ifdef PARA_CONTACT
                                                                                  ,myrank &
!#endif
                                                                                    )
                else
                  isOK = fstr_contact_init( fstrSOLID%contacts(c_contact+i), P%MESH)
                endif
         !       call fstr_write_contact( 6, fstrSOLID%contacts(c_contact+i) )
              endif
            enddo
            c_contact = c_contact+n

          else if( header_name == '!ISTEP'  ) then
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
			
          else if( header_name == '!WELD_LINE'  ) then
            fstrHEAT%WL_tot = fstrHEAT%WL_tot+1
            if( fstr_ctrl_get_WELDLINE( ctrl, hecMESH, HECMW_NAME_LEN, fstrHEAT%weldline(fstrHEAT%WL_tot) )/=0 ) then
                write(*,*) '### Error: Fail in read in Weld Line definition : ' , fstrHEAT%WL_tot
                write(ILOG,*) '### Error: Fail in read in Weld Line definition : ', fstrHEAT%WL_tot
                stop
            endif
			
          else if( header_name == '!SECTION'  ) then
            c_section = c_section+1
            if( fstr_ctrl_get_SECTION( ctrl, hecMESH )/=0 ) then
                write(*,*) '### Error: Fail in read in SECTION definition : ' , c_section
                write(ILOG,*) '### Error: Fail in read in SECTION definition : ', c_section
                stop
            endif

!== following material proerties ==
          else if( header_name == '!MATERIAL' ) then
            c_material = c_material+1
            if( fstr_ctrl_get_MATERIAL( ctrl, mName )/=0 ) then
                write(*,*) '### Error: Fail in read in material definition : ' , c_material
                write(ILOG,*) '### Error: Fail in read in material definition : ', c_material
                stop
            endif
            cid = c_material
            do i=1,hecMESH%material%n_mat
              if( trim(hecMESH%material%mat_name(i)) == trim(mName) ) then
                cid = i; exit
              endif
            enddo
            fstrSOLID%materials(cid)%name = mName
            if(c_material>hecMESH%material%n_mat) call initMaterial( fstrSOLID%materials(cid) )
          else if( header_name == '!ELASTIC' ) then
            if( c_material >0 ) then
               if( fstr_ctrl_get_ELASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(cid)%mtype,       &
                                             fstrSOLID%materials(cid)%nlgeom_flag, &
                                             fstrSOLID%materials(cid)%variables,   &
                                             fstrSOLID%materials(cid)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in elasticity definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in elasticity definition : ', cid
                 stop
              endif
            endif
          else if( header_name == '!PLASTIC' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_PLASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(cid)%mtype,       &
                                             fstrSOLID%materials(cid)%nlgeom_flag, &
                                             fstrSOLID%materials(cid)%variables,   &
                                             fstrSOLID%materials(cid)%table,   &
                                             fstrSOLID%materials(cid)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in plasticity definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in plasticity definition : ', cid
                 stop
              endif
            endif
          else if( header_name == '!HYPERELASTIC' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_HYPERELASTIC( ctrl,                                        &
                                             fstrSOLID%materials(cid)%mtype,       &
                                             fstrSOLID%materials(cid)%nlgeom_flag, &
                                             fstrSOLID%materials(cid)%variables )/=0 ) then
                 write(*,*) '### Error: Fail in read in elasticity definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in elasticity definition : ', cid
                 stop
              endif
            endif
          else if( header_name == '!VISCOELASTIC' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_VISCOELASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(cid)%mtype,       &
                                             fstrSOLID%materials(cid)%nlgeom_flag, &
                                             fstrSOLID%materials(cid)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in plasticity definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in plasticity definition : ', cid
                 stop
              endif
            endif
          else if( header_name == '!TRS' ) then
            if( cid >0 ) then
              if( fstrSOLID%materials(cid)%mtype/=VISCOELASTIC ) then
                 write(*,*) '### WARNING: TRS can only be defined for viscoelastic material! It is ignored! ' 
                 write(ILOG,*) '### WARNING: TRS can only be defined for viscoelastic material! It is ignored! ' 
              else
                 if( fstr_ctrl_get_TRS( ctrl, fstrSOLID%materials(cid)%mtype, fstrSOLID%materials(cid)%variables)/=0 ) then
                   write(*,*) '### Error: Fail in read in TRS definition : ' , cid
                   write(ILOG,*) '### Error: Fail in read in TRS definition : ', cid
                   stop
                 endif
              endif
            endif
          else if( header_name == '!CREEP' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_VISCOPLASTICITY( ctrl,                                        &
                                             fstrSOLID%materials(cid)%mtype,       &
                                             fstrSOLID%materials(cid)%nlgeom_flag, &
                                             fstrSOLID%materials(cid)%dict)/=0 ) then
                 write(*,*) '### Error: Fail in read in plasticity definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in plasticity definition : ', cid
                 stop
              endif
            endif
          else if( header_name == '!DENSITY' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_DENSITY( ctrl, fstrSOLID%materials(cid)%variables )/=0 ) then
                 write(*,*) '### Error: Fail in read in density definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in density definition : ', cid
                 stop
               endif
            endif
          else if( header_name == '!EXPANSION_COEFF' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_EXPANSION_COEFF( ctrl, fstrSOLID%materials(cid)%variables, &
                                             fstrSOLID%materials(cid)%dict)/=0 )  then
                 write(*,*) '### Error: Fail in read in expansion coefficient definition : ' , cid
                 write(ILOG,*) '### Error: Fail in read in expansion coefficient definition : ', cid
                 stop
               endif
            endif
          else if( header_name == '!USER_MATERIAL' ) then
            if( cid >0 ) then
               if( fstr_ctrl_get_USERMATERIAL( ctrl, fstrSOLID%materials(cid)%mtype,   &
                 fstrSOLID%materials(cid)%nlgeom_flag, fstrSOLID%materials(cid)%nfstatus, &
                 fstrSOLID%materials(cid)%variables(101:) )/=0 ) then
                 write(*,*) '### Error: Fail in read in user defined material : ' , cid
                 write(ILOG,*) '### Error: Fail in read in user defined material : ', cid
                 stop
               endif
            endif


! == Following output control ==
          else if( header_name == '!WRITE' ) then
            call fstr_ctrl_get_output( ctrl, outctrl, islog, resul, visual, femap )
            if( islog == 1 ) then
                c_output=1
                outctrl%filename = trim(logfileNAME)
                outctrl%filenum = ILOG
                call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
            endif
            if( femap == 1 ) then
                c_output=2
                write( outctrl%filename, *) 'utable.',myrank,".dat"
                outctrl%filenum = IUTB
                call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
                open( unit=outctrl%filenum, file=outctrl%filename, status='REPLACE' )
            endif
            if( resul == 1 ) then
                c_output=3
                call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
            endif
            if( visual == 1 ) then
                c_output=4
                call fstr_copy_outctrl(fstrSOLID%output_ctrl(c_output), outctrl)
            endif

          else if( header_name == '!OUTPUT_RES' ) then
            c_output=3
            if( .not. fstr_ctrl_get_outitem( ctrl, hecMESH, fstrSOLID%output_ctrl(c_output)%outinfo ) ) then
               write(*,*) '### Error: Fail in read in node output definition : ' , c_output
               write(ILOG,*) '### Error: Fail in read in node output definition : ', c_output
               stop
            endif
			if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name /= 'ALL' ) then
               c_output=2
               do i=1,hecMESH%node_group%n_grp
                 if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name == hecMESH%node_group%grp_name(i) ) then
                   fstrSOLID%output_ctrl(c_output)%outinfo%grp_id = i; exit
                 endif
               enddo
            endif
          else if( header_name == '!OUTPUT_VIS' ) then
            c_output=4
            if( .not. fstr_ctrl_get_outitem( ctrl, hecMESH, fstrSOLID%output_ctrl(c_output)%outinfo ) ) then
              write(*,*) '### Error: Fail in read in element output definition : ' , c_output
              write(ILOG,*) '### Error: Fail in read in element output definition : ', c_output
              stop
            endif
			if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name /= 'ALL' ) then
               c_output=2
               do i=1,hecMESH%node_group%n_grp
                 if( fstrSOLID%output_ctrl(c_output)%outinfo%grp_id_name == hecMESH%node_group%grp_name(i) ) then
                   fstrSOLID%output_ctrl(c_output)%outinfo%grp_id = i; exit
                 endif
               enddo
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

          if( fstr_ctrl_seek_next_header(ctrl) == 0) exit
        end do

! ----- material type judgement. in case of infinitive analysis, nlgeom_flag=0
        if( P%PARAM%solution_type == kstSTATIC ) then
          do i=1, c_material
            fstrSOLID%materials(i)%nlgeom_flag = 0
          enddo
        endif

        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then 
          allocate ( fstrSOLID%temperature( hecMESH%n_node )      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TEMPERATURE>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          fstrSOLID%temperature = REF_TEMP
          allocate ( fstrSOLID%last_temp( hecMESH%n_node )      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, LAST_TEMP>'
              write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          fstrSOLID%last_temp = 0.d0
          allocate ( fstrSOLID%temp_bak( hecMESH%n_node )      ,STAT=ierror )
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
           !call fstr_print_steps( 6, fstrSOLID%step_ctrl )
        else
           if( P%PARAM%solution_type==kstNLSTATIC .or. P%DYN%nlflag/=0 ) then
              write( *,* ) " ERROR: STEP not defined!"
              write( idbg,* ) "ERROR: STEP not defined!"
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
           endif
			  
           print *, "Step control not defined! Using default step=1"
           fstrSOLID%nstep_tot = 1
           allocate( fstrSOLID%step_ctrl(1) )
           call init_stepInfo( fstrSOLID%step_ctrl(1) )
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
        endif

        if( p%PARAM%solution_type /= kstHEAT) call fstr_element_init( hecMESH, fstrSOLID )
        if( p%PARAM%solution_type==kstSTATIC .or. p%PARAM%solution_type==kstNLSTATIC .or. &
           p%PARAM%solution_type==kstDYNAMIC .or. p%PARAM%solution_type==kstEIGEN .or.   &
           p%PARAM%solution_type==6 )            &
          call fstr_solid_alloc( hecMESH, fstrSOLID )

        if( p%PARAM%solution_type == kstHEAT) then
          p%PARAM%fg_irres = fstrSOLID%output_ctrl(3)%freqency
          p%PARAM%fg_iwres = fstrSOLID%output_ctrl(4)%freqency
        endif
		
		n_totallayer_ls = 1
		mixedflag = 0
		
      do it=1,hecMESH%section%n_sect
         cid = hecMESH%section%sect_mat_ID_item(it)
         aa =  int(fstrSOLID%materials(cid)%variables(M_TOTAL_LAYER))
	     if (aa > n_totallayer_ls)then
		    n_totallayer_ls = aa
         endif
      enddo
      do it=1,hecMESH%n_elem_type
         bb =  hecMESH%elem_type_item(it)
         if (bb == 781 .or. bb == 761)then
		    mixedflag = 1
         endif
      enddo
	  
        call fstr_setup_post( ctrl, P, n_totallayer_ls, mixedflag)
        rcode = fstr_ctrl_close( ctrl )

end subroutine fstr_setup


!> Initializer of structure fstr_solid
subroutine fstr_solid_init( hecMESH, fstrSOLID )
        use m_fstr
        type(hecmwST_local_mesh),target :: hecMESH
        type(fstr_solid)                :: fstrSOLID

        integer :: ndof, ntotal, ierror, ic_type

        fstrSOLID%file_type  = kbcfFSTR

        fstrSOLID%BOUNDARY_ngrp_tot = 0
        fstrSOLID%CLOAD_ngrp_tot    = 0
        fstrSOLID%DLOAD_ngrp_tot    = 0
        fstrSOLID%DLOAD_follow      = 1
        fstrSOLID%TEMP_ngrp_tot     = 0
        fstrSOLID%SPRING_ngrp_tot   = 0
        fstrSOLID%TEMP_irres        = 0
        fstrSOLID%TEMP_tstep        = 1
        fstrSOLID%TEMP_interval     = 1
        fstrSOLID%VELOCITY_ngrp_tot = 0
        fstrSOLID%ACCELERATION_ngrp_tot = 0
        fstrSOLID%COUPLE_ngrp_tot   = 0

        fstrSOLID%restart_nout= 0

end subroutine fstr_solid_init

!> Initializer of structure fstr_solid
subroutine fstr_solid_alloc( hecMESH, fstrSOLID )
        use m_fstr
        type(hecmwST_local_mesh),target :: hecMESH
        type(fstr_solid)                :: fstrSOLID

        integer :: ndof, ntotal, ierror, ic_type

        ndof=hecMESH%n_dof
        ic_type= hecMESH%elem_type_item(1)
        if( ic_type==741 .or. ic_type==743 .or. ic_type==731 ) ndof=6

        ntotal=ndof*hecMESH%n_node
        allocate ( fstrSOLID%GL( ntotal )          ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, GL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
!        allocate ( fstrSOLID%TOTAL_DISP( ntotal )  ,STAT=ierror )
!            if( ierror /= 0 ) then
!              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TOTAL_DISP>'
!              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
!              call flush(idbg)
!              call hecmw_abort( hecmw_comm_get_comm())
!            end if
        allocate ( fstrSOLID%unode( ntotal )  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, unode>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%dunode( ntotal )  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, dunode>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%ddunode( ntotal )  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, ddunode>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate ( fstrSOLID%QFORCE( ntotal )      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <FSTR_SOLID, QFORCE>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

        fstrSOLID%GL(:)=0.d0
!        fstrSOLID%TOTAL_DISP(:)=0.d0
        fstrSOLID%unode(:)      = 0.d0
        fstrSOLID%dunode(:)     = 0.d0
        fstrSOLID%ddunode(:)    = 0.d0
        fstrSOLID%QFORCE(:)     = 0.d0
        fstrSOLID%FACTOR( 1:2 ) = 0.d0

        ! for MPC
        fstrSOLID%n_fix_mpc = hecMESH%mpc%n_mpc
        if( fstrSOLID%n_fix_mpc>0 ) then
          allocate( fstrSOLID%mpc_const( fstrSOLID%n_fix_mpc ) )
          fstrSOLID%mpc_const(:) = hecMESH%mpc%mpc_const(:)
        endif

        ! initialize for linear static problems 
        fstrSOLID%FACTOR(2)=1.d0
        fstrSOLID%FACTOR(1)=0.d0
end subroutine fstr_solid_alloc

!> Initialize elements info in static calculation
subroutine fstr_element_init( hecMESH, fstrSOLID )
        use elementInfo
        use mMechGauss
        use m_fstr
        type(hecmwST_local_mesh),target :: hecMESH
        type(fstr_solid)                :: fstrSOLID

        integer :: i, j, ng, isect, ndof, id, nn            

        if( hecMESH%n_elem <=0 ) then
             stop "no element defined!"
        endif
        allocate( fstrSOLID%elements(hecMESH%n_elem) )
        do i=1,hecMESH%n_elem
          fstrSOLID%elements(i)%etype = hecMESH%elem_type(i)
          if( hecMESH%elem_type(i)==301 ) fstrSOLID%elements(i)%etype=111
          if (hecmw_is_etype_link(fstrSOLID%elements(i)%etype)) cycle
          ng = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          if(ng>0) allocate( fstrSOLID%elements(i)%gausses( ng ) )      
          
          isect= hecMESH%section_ID(i)
          ndof = getSpaceDimension( fstrSOLID%elements(i)%etype )
          if (ndof == 2) then              ! why do this???
            id=hecMESH%section%sect_opt(isect)
            if( id==0 ) then
              fstrSOLID%elements(i)%iset=1
            else if( id==1) then
              fstrSOLID%elements(i)%iset=0
            else if( id==2) then
              fstrSOLID%elements(i)%iset=2
            endif
          endif

          if( isect<0 .or. isect>hecMESH%section%n_sect )   &
             stop "Error in element's section definition"
          id = hecMESH%section%sect_mat_ID_item(isect)
          fstrSOLID%materials(id)%cdsys_ID = hecMESH%section%sect_orien_ID(isect)
          do j=1,ng
             fstrSOLID%elements(i)%gausses(j)%pMaterial => fstrSOLID%materials(id)
             call fstr_init_gauss( fstrSOLID%elements(i)%gausses( j )  )
          enddo
          
          nn = hecmw_get_max_node(hecMESH%elem_type(i))           
          allocate(fstrSOLID%elements(i)%equiForces(nn*ndof))   
          fstrSOLID%elements(i)%equiForces = 0.0d0                 
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
          if( associated(fstrSOLID%elements(i)%gausses) ) then           
            do j=1,size(fstrSOLID%elements(i)%gausses)
              call fstr_finalize_gauss(fstrSOLID%elements(i)%gausses(j))
            enddo
            deallocate( fstrSOLID%elements(i)%gausses )
          endif                                                              
          if(associated(fstrSOLID%elements(i)%equiForces) ) then            
            deallocate(fstrSOLID%elements(i)%equiForces)    
          endif                                                               
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
              if( fstrSOLID%output_ctrl(i)%filenum==IUTB )  &
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
        if( associated(fstrSOLID%last_temp) ) then
            deallocate(fstrSOLID%last_temp       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <FSTR_SOLID, reftemp>'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(fstrSOLID%temp_bak) ) then
            deallocate(fstrSOLID%temp_bak       ,STAT=ierror)
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
        fstrHEAT%WL_tot      = 0
end subroutine fstr_heat_init

!> Initial setting of eigen ca;culation
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
        fstrDYNAMIC%nlflag   = 0
        fstrDYNAMIC%idx_eqa  = 1
        fstrDYNAMIC%idx_resp = 1
        fstrDYNAMIC%n_step   = 1
        fstrDYNAMIC%t_start  = 0.0
        fstrDYNAMIC%t_curr = 0.0d0             
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
        fstrDYNAMIC%ngrp_monit   = 0
        fstrDYNAMIC%nout_monit   = 1
        fstrDYNAMIC%iout_list(1) = 0
        fstrDYNAMIC%iout_list(2) = 0
        fstrDYNAMIC%iout_list(3) = 0
        fstrDYNAMIC%iout_list(4) = 0
        fstrDYNAMIC%iout_list(5) = 0
        fstrDYNAMIC%iout_list(6) = 0 
		
end subroutine fstr_dynamic_init


!> Initial setting of dynamic calculation
subroutine fstr_dynamic_alloc( hecMESH, fstrDYNAMIC )
        use m_fstr
        type(hecmwST_local_mesh),target :: hecMESH
        type(fstr_dynamic) :: fstrDYNAMIC
		
        integer :: ierror, ndof,nnod
        
        ndof=hecMESH%n_dof
        nnod=hecMESH%n_node
        if(fstrDYNAMIC%idx_eqa == 11) then
          allocate( fstrDYNAMIC%DISP(ndof*nnod,3)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%VEL (ndof*nnod,1)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%ACC (ndof*nnod,1)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        else
          allocate( fstrDYNAMIC%DISP(ndof*nnod,2)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, DISP>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%VEL (ndof*nnod,2)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
          allocate( fstrDYNAMIC%ACC (ndof*nnod,2)  ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, ACC>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        
        
        allocate( fstrDYNAMIC%VEC1(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC1>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate( fstrDYNAMIC%VEC2(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC2>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        allocate( fstrDYNAMIC%VEC3(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, VEC3>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
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

 subroutine fstr_setup_post( ctrl, P, n_totallayer_ls, mixedflag)
        implicit none                      
        integer(kind=kint) :: ctrl, n_totallayer_ls, mixedflag
        type(fstr_param_pack) :: P

                                           ! JP-6
        if( P%PARAM%solution_type == kstSTATIC &
         .or. P%PARAM%solution_type == kstNLSTATIC  &
         .or. P%PARAM%solution_type == kstEIGEN  &
         .or. P%PARAM%solution_type == kstDYNAMIC &
         .or. P%PARAM%solution_type == 6 ) then
                ! Memory Allocation for Result Vectors ------------
                if( P%MESH%n_dof == 6 ) then
                        allocate ( P%SOLID%STRAIN  (12*n_totallayer_ls*P%MESH%n_node))
                        allocate ( P%SOLID%STRESS  (14*n_totallayer_ls*P%MESH%n_node))
                        allocate ( P%SOLID%ESTRAIN  (12*n_totallayer_ls*P%MESH%n_elem))
                        allocate ( P%SOLID%ESTRESS  ((14*n_totallayer_ls+6)*P%MESH%n_elem))
			    elseif( mixedflag == 1) then
                        allocate ( P%SOLID%STRAIN  (12*n_totallayer_ls*P%MESH%n_node))
                        allocate ( P%SOLID%STRESS  (14*n_totallayer_ls*P%MESH%n_node))
                        allocate ( P%SOLID%ESTRAIN  (12*n_totallayer_ls*P%MESH%n_elem))
                        allocate ( P%SOLID%ESTRESS  ((14*n_totallayer_ls+6)*P%MESH%n_elem))
                else
                        allocate ( P%SOLID%STRAIN  (6*P%MESH%n_node))
                        allocate ( P%SOLID%STRESS  (7*P%MESH%n_node))
                        allocate ( P%SOLID%ESTRAIN  (6*P%MESH%n_elem))
                        allocate ( P%SOLID%ESTRESS  (7*P%MESH%n_elem))
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
                call fstr_setup_visualize( ctrl )
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
     !   scaling    => svIarray(7)
     !   iterlog    => svIarray(21)
     !   timelog    => svIarray(22)
     !   dumptype   => svIarray(31)
     !   dumpexit   => svIarray(32)
     !   usejad     => svIarray(33)
     !   ncolor_in  => svIarray(34)

     !   resid      => svRarray(1)
     !   sigma_diag => svRarray(2)
     !   sigma      => svRarray(3)
     !   thresh     => svRarray(4)
     !   filter     => svRarray(5)

        rcode = fstr_ctrl_get_SOLVER( ctrl,                      &
                        svIarray(2), svIarray(3), svIarray(4), svIarray(21), svIarray(22), &
                        svIarray(1), svIarray(5), svIarray(6), svIarray(7), &
                        svIarray(31), svIarray(32), svIarray(33), svIarray(34), &
                        svRarray(1), svRarray(2), svRarray(3),                &
                        svRarray(4), svRarray(5) )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        if( svIarray(2) <= 100 ) then                      
                svIarray(99) = 1  ! indirect method
        else
                svIarray(99) = svIarray(2)-99 !2  ! direct method          
        end if

end subroutine fstr_setup_SOLVER

!* ----------------------------------------------------------------------------------------------- *!
!> Read in !ORIENTATION       
!* ----------------------------------------------------------------------------------------------- *!

integer function fstr_setup_ORIENTATION( ctrl, hecMESH, cnt, coordsys )
        implicit none
        integer(kind=kint)         :: ctrl
        type( hecmwST_local_mesh ) :: hecMESH
        integer                    :: cnt
        type( tLocalCoordSys )     :: coordsys

        integer                   :: j, iS, iE, grp_id(1)
        character(len=HECMW_NAME_LEN) :: grp_id_name(1)

        integer :: nid, dtype
        character(len=HECMW_NAME_LEN) :: data_fmt
        real(kind=kreal) :: fdum, xyza(3), xyzb(3), xyzc(3), ff1(3), ff2(3), ff3(3)

        fstr_setup_ORIENTATION = -1

		nid = 1
        coordsys%sys_type = 10

        nid = 1	
        data_fmt = 'COORDINATES,NODES '
        if( fstr_ctrl_get_param_ex( ctrl, 'DEFINITION ', data_fmt, 0, 'P', nid )/=0 ) return
        dtype = nid-1
        coordsys%sys_type = coordsys%sys_type + dtype

        if( fstr_ctrl_get_param_ex( ctrl, 'NAME ',  '# ',  1, 'S', grp_id_name(1) )/= 0) return
        coordsys%sys_name = grp_id_name(1)
        call fstr_strupr( coordsys%sys_name )

        if( dtype==0 ) then		
          data_fmt = "RRRRRRrrr "
          xyzc(:) = 0.d0
          if( fstr_ctrl_get_data_array_ex( ctrl, data_fmt, xyza(1), xyza(2),  &
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
          if( fstr_ctrl_get_data_array_ex( ctrl, data_fmt, coordsys%node_ID(1),  &
                coordsys%node_ID(2), coordsys%node_ID(3) )/=0 ) return
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

        fstr_setup_ORIENTATION = 0
end function fstr_setup_ORIENTATION


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
        rcode = fstr_ctrl_get_STEP( ctrl, amp, iproc )
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


!-----------------------------------------------------------------------------!
!> Read in !ECHO                                                            !
!-----------------------------------------------------------------------------!
subroutine fstr_setup_ECHO( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode

        rcode = fstr_ctrl_get_ECHO( ctrl,        &
                                P%PARAM%fg_echo )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

end subroutine fstr_setup_ECHO


!-----------------------------------------------------------------------------!
!> Read in !RESTART                                                         !
!-----------------------------------------------------------------------------!
subroutine fstr_setup_RESTART( ctrl, nout )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: nout

        integer(kind=kint) :: rcode
        nout = 0
        rcode = fstr_ctrl_get_param_ex( ctrl, 'FREQUENCY ', '# ', 0, 'I', nout )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

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
                        P%PARAM%fg_couple_type,       &
                        P%PARAM%fg_couple_first,      &
                        P%PARAM%fg_couple_window,     &
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
!> Read in !BOUNDARY                                                    !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_BOUNDARY( ctrl, counter, P )
  implicit none
  integer(kind=kint) :: ctrl
  integer(kind=kint) :: counter
  type(fstr_param_pack) :: P

  integer(kind=kint) :: rcode
  integer(kind=kint) :: type = 0
  character(HECMW_NAME_LEN) :: amp
  integer(kind=kint) :: amp_id
  character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
  integer(kind=kint),pointer :: dof_ids (:)
  integer(kind=kint),pointer :: dof_ide (:)
  real(kind=kreal),pointer :: val_ptr(:)
  integer(kind=kint) :: i, n, old_size, new_size

  integer(kind=kint) :: gid

  gid = 1
  rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )
!  rcode = fstr_ctrl_get_param_ex( ctrl, 'TYPE ', 'FSTR,NASTRAN ', 0, 'P', type )
!  if( rcode < 0 ) call fstr_ctrl_err_stop
!  if( rcode == 1 ) type = 0 ! PARAM_NOTHING

!  if( type == 0 ) then

    ! ENTIRE -----------------------------------------------
    P%SOLID%file_type = kbcfFSTR

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return
    old_size = P%SOLID%BOUNDARY_ngrp_tot
    new_size = old_size + n
    P%SOLID%BOUNDARY_ngrp_tot = new_size
    call fstr_expand_integer_array (P%SOLID%BOUNDARY_ngrp_GRPID, old_size, new_size )
    call fstr_expand_integer_array (P%SOLID%BOUNDARY_ngrp_ID,    old_size, new_size )
    call fstr_expand_integer_array (P%SOLID%BOUNDARY_ngrp_type,  old_size, new_size )
    call fstr_expand_real_array    (P%SOLID%BOUNDARY_ngrp_val,   old_size, new_size )
    call fstr_expand_integer_array (P%SOLID%BOUNDARY_ngrp_amp,   old_size, new_size )

    allocate( grp_id_name(n) )
    allocate( dof_ids (n) )
    allocate( dof_ide (n) )

    amp = ' '
    val_ptr => P%SOLID%BOUNDARY_ngrp_val(old_size+1:)
    val_ptr = 0
    rcode = fstr_ctrl_get_BOUNDARY( ctrl, amp, grp_id_name, HECMW_NAME_LEN, dof_ids, dof_ide, val_ptr) 
    if( rcode /= 0 ) call fstr_ctrl_err_stop
    call amp_name_to_id( P%MESH, '!BOUNDARY', amp, amp_id )
    P%SOLID%BOUNDARY_ngrp_GRPID(old_size+1:new_size) = gid
    call node_grp_name_to_id_ex( P%MESH, '!BOUNDARY', n, grp_id_name, P%SOLID%BOUNDARY_ngrp_ID(old_size+1:))

    do i = 1, n
      if( (dof_ids(i) < 1).or.(6 < dof_ids(i)).or.(dof_ide(i) < 1).or.(6 < dof_ide(i)) ) then
        write(*,*) 'fstr contol file error : !BOUNDRAY : range of dof_ids and dof_ide is from 1 to 6'
        write(ILOG,*) 'fstr contol file error : !BOUNDRAY : range of dof_ids and dof_ide is from 1 to 6'
        call fstr_ctrl_err_stop
      end if
      P%SOLID%BOUNDARY_ngrp_type(old_size+i) = 10 * dof_ids(i) + dof_ide(i)
      P%SOLID%BOUNDARY_ngrp_amp(old_size+i) = amp_id
    end do

    deallocate( grp_id_name )
    deallocate( dof_ids )
    deallocate( dof_ide )
!  else
!   ! NASTRAN ---------------------------------------------
!
!     P%SOLID%file_type = kbcfNASTRAN
!     call fstr_setup_solid_nastran( ctrl, P%MESH, P%SOLID )
!  end if

end subroutine fstr_setup_BOUNDARY


!-----------------------------------------------------------------------------!
!> Read in !CLOAD                                                           !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_CLOAD( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: val_ptr(:)
        integer(kind=kint),pointer :: id_ptr(:)
        integer(kind=kint) :: i, n, old_size, new_size
        integer(kind=kint) :: gid

        if( P%SOLID%file_type /= kbcfFSTR ) return
        gid = 1
        rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )
!        P%SOLID%CLOAD_ngrp_tot = 0
!
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return
        old_size = P%SOLID%CLOAD_ngrp_tot
        new_size = old_size + n
        P%SOLID%CLOAD_ngrp_tot = new_size
        call fstr_expand_integer_array ( P%SOLID%CLOAD_ngrp_GRPID, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%CLOAD_ngrp_ID,  old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%CLOAD_ngrp_DOF, old_size, new_size )
        call fstr_expand_real_array    ( P%SOLID%CLOAD_ngrp_val, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%CLOAD_ngrp_amp, old_size, new_size )

        allocate( grp_id_name(n))
        amp = ' '
        val_ptr => P%SOLID%CLOAD_ngrp_val(old_size+1:)
        id_ptr =>P%SOLID%CLOAD_ngrp_DOF(old_size+1:)
        val_ptr = 0
        rcode = fstr_ctrl_get_CLOAD( ctrl, amp, grp_id_name, HECMW_NAME_LEN, id_ptr, val_ptr )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        call amp_name_to_id( P%MESH, '!CLOAD', amp, amp_id ) 
        do i=1,n
                P%SOLID%CLOAD_ngrp_amp(old_size+i) = amp_id
        end do
        P%SOLID%CLOAD_ngrp_GRPID(old_size+1:new_size) = gid
        call node_grp_name_to_id_ex( P%MESH, '!CLOAD', n, grp_id_name, P%SOLID%CLOAD_ngrp_ID(old_size+1:))

        deallocate( grp_id_name )

end subroutine fstr_setup_CLOAD

!-----------------------------------------------------------------------------!
!> Read !FLOAD                                                        !
!-----------------------------------------------------------------------------!
include 'fstr_ctrl_freq.f90'

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

!> Read in !DLOAD
subroutine fstr_setup_DLOAD( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        integer(kind=kint) :: follow
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: new_params(:,:)
        logical,pointer :: fg_surface(:)
        integer(kind=kint),pointer :: lid_ptr(:)
        integer(kind=kint) :: i, j, n, old_size, new_size
        integer(kind=kint) :: gid

        if( P%SOLID%file_type /= kbcfFSTR ) return

        gid = 1
        rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )

        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return
        old_size = P%SOLID%DLOAD_ngrp_tot
        new_size = old_size + n
        P%SOLID%DLOAD_ngrp_tot = new_size
        call fstr_expand_integer_array ( P%SOLID%DLOAD_ngrp_GRPID, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%DLOAD_ngrp_ID,  old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%DLOAD_ngrp_LID, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%DLOAD_ngrp_amp, old_size, new_size )
        call fstr_expand_dload_array ( P%SOLID%DLOAD_ngrp_params, old_size, new_size )

        allocate( grp_id_name(n))
        allocate( new_params(0:6,n))
        allocate( fg_surface(n))
        new_params = 0
        amp = ' '
        follow = P%SOLID%DLOAD_follow
        lid_ptr => P%SOLID%DLOAD_ngrp_LID(old_size+1:)
        rcode = fstr_ctrl_get_DLOAD( ctrl, amp, follow, &
                        grp_id_name, HECMW_NAME_LEN,    &
                        lid_ptr, new_params )
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        call amp_name_to_id( P%MESH, '!DLOAD', amp, amp_id ) 
        P%SOLID%DLOAD_follow = follow
        do i=1,n
                P%SOLID%DLOAD_ngrp_amp(old_size+i) = amp_id
                do j=0, 6
                        P%SOLID%DLOAD_ngrp_params(j,old_size+i) = new_params(j,i)
                end do
                fg_surface(i) =  ( lid_ptr(i) == 100 )
        end do
        P%SOLID%DLOAD_ngrp_GRPID(old_size+1:new_size) = gid
        call dload_grp_name_to_id_ex( P%MESH, n, grp_id_name, fg_surface, P%SOLID%DLOAD_ngrp_ID(old_size+1:))
        deallocate( grp_id_name )
        deallocate( new_params )
        deallocate( fg_surface )
end subroutine fstr_setup_DLOAD


!-----------------------------------------------------------------------------!
!> Read in !TEMPERATURE                                                  !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_TEMPERATURE( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode, gid
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: val_ptr(:)
        integer(kind=kint) :: n, old_size, new_size

        if( P%SOLID%file_type /= kbcfFSTR ) return

        gid = 1
        rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )

        n = fstr_ctrl_get_data_line_n( ctrl )
        old_size = P%SOLID%TEMP_ngrp_tot
        if( n > 0 ) then
          new_size = old_size + n
        else
          new_size = old_size + 1
        endif
        call fstr_expand_integer_array ( P%SOLID%TEMP_ngrp_GRPID, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%TEMP_ngrp_ID, old_size, new_size )
        call fstr_expand_real_array    ( P%SOLID%TEMP_ngrp_val,old_size, new_size )

        allocate( grp_id_name(n))
        val_ptr => P%SOLID%TEMP_ngrp_val( old_size+1: )

        rcode = fstr_ctrl_get_TEMPERATURE( ctrl,      &
                        P%SOLID%TEMP_irres,           &
                        P%SOLID%TEMP_tstep,           &
                        P%SOLID%TEMP_interval,        &
                        grp_id_name, HECMW_NAME_LEN,  &
                        val_ptr )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        P%SOLID%TEMP_ngrp_GRPID(old_size+1:new_size) = gid
        if( n > 0 ) then
          if( P%SOLID%TEMP_irres == 0 ) then
            P%SOLID%TEMP_ngrp_tot = new_size
            call node_grp_name_to_id_ex( P%MESH, '!TEMPERATURE', &
                    n, grp_id_name, P%SOLID%TEMP_ngrp_ID(old_size+1:))
          endif
          deallocate( grp_id_name )
        endif

end subroutine fstr_setup_TEMPERATURE


!-----------------------------------------------------------------------------!
!> Read in !SPRING                                                            !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_SPRING( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
        character(HECMW_NAME_LEN) :: amp
        integer(kind=kint) :: amp_id
        character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
        real(kind=kreal),pointer :: val_ptr(:)
        integer(kind=kint),pointer :: id_ptr(:)
        integer(kind=kint) :: i, n, old_size, new_size
        integer(kind=kint) :: gid

        if( P%SOLID%file_type /= kbcfFSTR ) return
        gid = 1
        rcode = fstr_ctrl_get_param_ex( ctrl, 'GRPID ',  '# ',            0, 'I', gid  )
        n = fstr_ctrl_get_data_line_n( ctrl )
        if( n == 0 ) return
        old_size = P%SOLID%SPRING_ngrp_tot
        new_size = old_size + n
        P%SOLID%SPRING_ngrp_tot = new_size
        call fstr_expand_integer_array ( P%SOLID%SPRING_ngrp_GRPID, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%SPRING_ngrp_ID,  old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%SPRING_ngrp_DOF, old_size, new_size )
        call fstr_expand_real_array    ( P%SOLID%SPRING_ngrp_val, old_size, new_size )
        call fstr_expand_integer_array ( P%SOLID%SPRING_ngrp_amp, old_size, new_size )

        allocate( grp_id_name(n))
        amp = ' '
        val_ptr => P%SOLID%SPRING_ngrp_val(old_size+1:)
        id_ptr =>P%SOLID%SPRING_ngrp_DOF(old_size+1:)
        val_ptr = 0
        rcode = fstr_ctrl_get_SPRING( ctrl, amp, grp_id_name, HECMW_NAME_LEN, id_ptr, val_ptr )
        if( rcode /= 0 ) call fstr_ctrl_err_stop

        call amp_name_to_id( P%MESH, '!SPRING', amp, amp_id ) 
        do i=1,n
                P%SOLID%SPRING_ngrp_amp(old_size+i) = amp_id
        end do
        P%SOLID%SPRING_ngrp_GRPID(old_size+1:new_size) = gid
        call node_grp_name_to_id_ex( P%MESH, '!SPRING', n, grp_id_name, P%SOLID%SPRING_ngrp_ID(old_size+1:))

        deallocate( grp_id_name )

end subroutine fstr_setup_SPRING


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
        character(HECMW_NAME_LEN) :: grp_id_name(1)
        integer(kind=kint) :: grp_id(1)

        rcode = fstr_ctrl_get_DYNAMIC( ctrl, &
                P%DYN%nlflag,  &
                P%DYN%idx_eqa, &
                P%DYN%idx_resp,&
                P%DYN%n_step,  &
                P%DYN%t_start, &
                P%DYN%t_end,   &
                P%DYN%t_delta, &
                P%DYN%ganma,   &
                P%DYN%beta,    &
                P%DYN%idx_mas, &
                P%DYN%idx_dmp, &
                P%DYN%ray_m,   &
                P%DYN%ray_k,   &
                P%DYN%nout,    &
                grp_id_name(1), HECMW_NAME_LEN,  &
                P%DYN%nout_monit,  &
                P%DYN%iout_list )

        if( rcode /= 0) call fstr_ctrl_err_stop

        if (P%DYN%idx_resp == 1) then
          call node_grp_name_to_id_ex( P%MESH, '!DYNAMIC', 1, grp_id_name, grp_id)
          P%DYN%ngrp_monit = grp_id(1)
        else
          read(grp_id_name,*) P%DYN%ngrp_monit
        endif

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
        integer(kind=kint) :: vType                     
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
        rcode = fstr_ctrl_get_VELOCITY( ctrl,  vType, amp,   &
                grp_id_name, HECMW_NAME_LEN,  &
                dof_ids, dof_ide, val_ptr )
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        P%SOLID%VELOCITY_type = vType 
        if( vType == kbcInitial ) P%DYN%VarInitialize = .true.                      
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
        integer(kind=kint) :: aType                             
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
        rcode = fstr_ctrl_get_ACCELERATION( ctrl,  aType, amp,   &
                grp_id_name, HECMW_NAME_LEN,  &
                dof_ids, dof_ide,  val_ptr)
        if( rcode /= 0 ) call fstr_ctrl_err_stop
        P%SOLID%ACCELERATION_type = aType
        if( aType == kbcInitial )P%DYN%VarInitialize = .true.                    
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


!*****************************************************************************!
!* MPC ***********************************************************************!
!*****************************************************************************!

!-----------------------------------------------------------------------------!
!> Read in !MPC                                                              !
!-----------------------------------------------------------------------------!

subroutine fstr_setup_MPC( ctrl, counter, P )
        implicit none
        integer(kind=kint) :: ctrl
        integer(kind=kint) :: counter
        type(fstr_param_pack), target :: P

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
!        call fstr_expand_integer_array ( P%MPC_RD%node1,  old_size, new_size )
!        call fstr_expand_integer_array ( P%MPC_RD%node2,  old_size, new_size )
!        call fstr_expand_integer_array ( P%MPC_RD%dof,    old_size, new_size )
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
        rcode = fstr_ctrl_get_MPC( ctrl, svRarray(11))
        if( rcode /= 0) call fstr_ctrl_err_stop
end subroutine fstr_setup_MPC


!*****************************************************************************!
!* IMPORTING NASTRAN BOUNDARY CONDITIONS *************************************!
!*****************************************************************************!

subroutine fstr_setup_solid_nastran( ctrl, hecMESH, fstrSOLID )
        implicit none
        integer(kind=kint) :: ctrl
        type (hecmwST_local_mesh) :: hecMESH
        type (fstr_solid        ) :: fstrSOLID
        write(ILOG,*) '### Error : In !BOUNNDARY, TYPE=NASTRAN is not supported.'
        call hecmw_abort( hecmw_comm_get_comm())
end subroutine fstr_setup_solid_nastran

!-----------------------------------------------------------------------------!
!> Read in !CONTACT                                                           !    
!-----------------------------------------------------------------------------!

subroutine fstr_setup_CONTACTALGO( ctrl, P )                                                                    
        implicit none
        integer(kind=kint) :: ctrl
!        integer(kind=kint) :: counter
        type(fstr_param_pack) :: P

        integer(kind=kint) :: rcode
       

        rcode = fstr_ctrl_get_CONTACTALGO( ctrl, P%PARAM%contact_algo )       
        if( rcode /= 0 ) call fstr_ctrl_err_stop

end subroutine fstr_setup_CONTACTALGO                                       


end module m_fstr_setup


