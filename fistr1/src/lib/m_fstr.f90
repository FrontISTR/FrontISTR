!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
! If new header is supported, change according to following method.    !
!  1) Increase FSTR_CTRL_HEADER_NUMBER                                 !
!  2) Add new header name to fstr_ctrl_header_names                    !
!  3) Describe new function to get parameters from control file        !
!     in fstr_ctrl.f90                                                 !
!  4) Describe new subroutine to set values of the parameter           !
!     in fstr_setup.f90                                                !
!  5) If initial values are necessary, set the value                   !
!     in subroutine fstr_setup_init in fstr_setup.f90                  !
!> This module defined coomon data and basic structures for analysis
module m_fstr
use hecmw
use lczparm
use m_common_struct
use m_step
use m_out
use m_timepoint
use mMechGauss
use mContactDef

implicit none

public

!C
!C-- CONSTANTS
!C
        ! general
        integer(kind=kint),parameter :: kYES       =   1
        integer(kind=kint),parameter :: kNO        =   0
        integer(kind=kint),parameter :: kON        =   1
        integer(kind=kint),parameter :: kOFF       =   0

        ! solution type (st)
        integer(kind=kint),parameter :: kstPRECHECK  =   0
        integer(kind=kint),parameter :: kstSTATIC    =   1
        integer(kind=kint),parameter :: kstEIGEN     =   2
        integer(kind=kint),parameter :: kstHEAT      =   3
        integer(kind=kint),parameter :: kstDYNAMIC   =   4
        !integer(kind=kint),parameter :: kstNLSTATIC  =   5
        integer(kind=kint),parameter :: kstSTATICEIGEN = 6
        integer(kind=kint),parameter :: kstNZPROF    =   7

        ! solver method (sm)    !CAUTION : (<=100):indirect, (>100):direct
        integer(kind=kint),parameter :: ksmCG       =   1
        integer(kind=kint),parameter :: ksmBiCGSTAB =   2
        integer(kind=kint),parameter :: ksmGMRES    =   3
        integer(kind=kint),parameter :: ksmGPBiCG   =   4
        integer(kind=kint),parameter :: ksmDIRECT   = 101

        ! contact analysis algorithm
        integer(kind=kint),parameter :: kcaSLagrange = 1
        integer(kind=kint),parameter :: kcaALagrange = 2

        ! boundary condition file type (bcf)
        integer(kind=kint),parameter :: kbcfFSTR    =   0  ! BC described in fstr control file (default)
        integer(kind=kint),parameter :: kbcfNASTRAN =   1  ! nastran file

        !
        integer(kind=kint),parameter :: kbcInitial    =  1
        integer(kind=kint),parameter :: kbcTransit    =  2

        ! restart type
        integer(kind=kint),parameter :: restart_outLast = 1
        integer(kind=kint),parameter :: restart_outAll  = 2

        ! statistics of newton iteration
        integer(kind=kint),parameter :: knstMAXIT  = 1 ! maximum number of newton iteration
        integer(kind=kint),parameter :: knstSUMIT  = 2 ! total number of newton iteration
        integer(kind=kint),parameter :: knstCITER  = 3 ! number of contact iteration

        ! section control
        integer(kind=kint),parameter :: kel361FI     =  1
        integer(kind=kint),parameter :: kel361BBAR   =  2
        integer(kind=kint),parameter :: kel361IC     =  3

!C
!C-- PARALLEL EXECUTION
!C
        integer(kind = kint) :: myrank
        integer(kind = kint) :: nprocs
!C
!C-- PARALLEL CONTACT FLAG
!C
        logical   ::  paraContactFlag = .false.
!C
!C-- FILE NAME
!C
        character(len=HECMW_FILENAME_LEN):: cntfilNAME
        character(len=HECMW_FILENAME_LEN):: restartfilNAME
!C
!C-- FILE HANDLER
!C
        integer(kind=kint),parameter :: ILOG=16 ! log
        integer(kind=kint),parameter :: ISTA=17 ! status
        integer                      :: IUTB=18 ! utable
        integer(kind=kint),parameter :: IMSG=51 ! message (myrank == 0 only)
        integer(kind=kint),parameter :: IDBG=52 ! debug
        integer(kind=kint),parameter :: IFVS=53 ! visual.ini file
        integer(kind=kint),parameter :: INEU=54 ! neutral file (heat)

        integer(kind=kint),parameter :: IRESOUT=100 ! ~110, keeping for result output file

!C
!C-- SOLVER CONTROL
!C
        integer(kind=kint) :: svIarray(100)
        real(kind=kreal)   :: svRarray(100)

!C
!C-- FLAG for ECHO/RESULT/POST
!C
        integer(kind=kint),pointer :: IECHO
        integer(kind=kint),pointer :: IRESULT
        integer(kind=kint),pointer :: IVISUAL
        ! for heat ...
        integer(kind=kint),pointer :: INEUTRAL  ! flag for femap neutral file
        integer(kind=kint),pointer :: IRRES     ! flag for restart, read
        integer(kind=kint),pointer :: IWRES     ! flag for restart, write
        integer(kind=kint),pointer :: NRRES     ! position of restart read
        integer(kind=kint),pointer :: NPRINT    ! interval of write
!C
!C-- REFTEMP
!C-
        real(kind=kreal),pointer :: REF_TEMP
!C
!C-- ANALYSIS CONTROL for NLGEOM and HEAT
!C
        !for static
        !CAUTION)
        !   DT,ETIME and EPS will be integrated in fstr_param
        !   when !STATIC supports mult-lines.
        real(kind=kreal)           :: DT    ! /=fstr_param%dtime
        real(kind=kreal)           :: ETIME ! /=fstr_param%etime
        integer(kind=kint)         :: ITMAX
        real(kind=kreal)           :: EPS   ! /=fstr_param%eps

        !for heat
        real(kind=kreal)   :: TEMTOL
!C
!C ----------------------------------------------------------------------------
!C
!>  FSTR INNER CONTROL PARAMETERS  (fstrPARAM)
!C   Caution: global pointer parameter IECHO, IRESULT and IVISUAL, etc.
!C            points to corresponded parameter in the structure
!C
        type fstr_param
                integer(kind=kint) :: solution_type !< solution type number
                logical            :: nlgeom        !< is geometrical nonlinearity considered
                integer(kind=kint) :: solver_method !< solver method number

                !!STATIC !HEAT
                integer(kind=kint) :: analysis_n      !< Number of analysis
                real(kind=kreal),pointer  :: dtime(:) !< (=DT)    STEP_DLTIME
                real(kind=kreal),pointer  :: etime(:) !< (/=ETIME) STEP_EETIME
                real(kind=kreal),pointer  :: dtmin(:) !< (=DTMIN) STEP_DELMIN
                real(kind=kreal),pointer  :: delmax(:)!< (=DTMAX) STEP_DELMAX
                integer(kind=kint),pointer:: itmax(:) !< (/=ITMAX)
                real(kind=kreal),pointer  :: eps(:)   !< (/=ESP)
                real(kind=kreal)          :: ref_temp !< (=REF_TEMP)

                ! output control
                integer(kind=kint) :: fg_echo       !< output echo   (kYES/kNO) (=IECHO)
                integer(kind=kint) :: fg_result     !< output result (kYES/kNO) (=IRESULT)
                integer(kind=kint) :: fg_visual     !< visualization (kYES/kNO) (=IVISUAL)

                ! for heat ...
                integer(kind=kint) :: fg_neutral    !< write by neutral (=INEUTRAL)
                integer(kind=kint) :: fg_irres      !< restart read     (=IRRES)
                integer(kind=kint) :: fg_iwres      !< restart write    (=IWRES)
                integer(kind=kint) :: nrres         !< NRRES
                integer(kind=kint) :: nprint        !< NPRINT

                ! index table for global node ID sorting
                integer(kind=kint) :: n_node        !< hecMESH%n_node
                integer(kind=kint) :: nn_internal   !< hecMESH%nn_internal
                integer(kind=kint),pointer :: global_local_ID(:,:)  !> (2:nn_internal) (1,:):global, (2,:):local

                ! for couple analysis
                integer( kind=kint ) :: fg_couple          !< (default:0)
                integer( kind=kint ) :: fg_couple_type     !< (default:0)
                integer( kind=kint ) :: fg_couple_first    !< (default:0)
                integer( kind=kint ) :: fg_couple_window   !< (default:0)

                ! for restart control
                integer( kind=kint ) :: restart_out_type   !< output type of restart file
                integer( kind=kint ) :: restart_version    !< version of restart file

                ! for contact analysis
                integer( kind=kint ) :: contact_algo       !< contact analysis algorithm number(SLagrange or Alagrange)

                ! for auto increment and cutback
                real(kind=kreal)     :: ainc_Rs            !< time increment decreasing ratio
                real(kind=kreal)     :: ainc_Rl            !< time increment increasing ratio
                integer( kind=kint ) :: NRbound_s(10)      !< # of NR iteration bound to decrease time increment
                integer( kind=kint ) :: NRbound_l(10)      !< # of NR iteration bound to increase time increment
                integer( kind=kint ) :: NRtimes_s          !< # of times that decreasing condition is satisfied
                integer( kind=kint ) :: NRtimes_l          !< # of times that increasing condition is satisfied
                real(kind=kreal)     :: ainc_Rc            !< time increment decreasing ratio for cutback
                integer( kind=kint ) :: CBbound            !< maximum # of successive cutback
                type(time_points), pointer :: timepoints(:)  !< time points data
!
        end type fstr_param
!
!
!
        ! **** GLOBAL VARIABLE INITIALIZED IN FSTR_SETUP *****
        type( fstr_param ),target :: fstrPR
!
!C ----------------------------------------------------------------------------
!C
!> Data for STATIC ANSLYSIS  (fstrSOLID)
!C
        type fstr_solid_physic_val
                ! If using shell, substitute average value for this structure. and Substitute Plus and Minus for SHELL_PLUS and SHELL_MINUS
                ! If using laminated shell, substitute average whole laminated elements value for this structure. and Substitute whole average plus and minus for SHELL_PLUS and SHELL_MINUS
                ! And substitute each laminated shell's value for layer(1,2,3,4,)

                real(kind=kreal), pointer :: STRESS(:)    !< nodal stress
                real(kind=kreal), pointer :: STRAIN(:)    !< nodal strain
                real(kind=kreal), pointer :: MISES(:)    !< nodal MISES

                real(kind=kreal), pointer :: PSTRESS(:)   !< nodal principal stress
                real(kind=kreal), pointer :: PSTRAIN(:)   !< nodal principal strain
                real(kind=kreal), pointer :: PSTRESS_VECT(:,:)   !< nodal principal stress vector
                real(kind=kreal), pointer :: PSTRAIN_VECT(:,:)   !< nodal principal strain vector

                real(kind=kreal), pointer :: ESTRESS(:)   !< elemental stress
                real(kind=kreal), pointer :: ESTRAIN(:)   !< elemental strain
                real(kind=kreal), pointer :: EMISES(:)    !< elemental MISES

                real(kind=kreal), pointer :: EPSTRESS(:)   !< elemental principal stress
                real(kind=kreal), pointer :: EPSTRAIN(:)   !< elemental principal strain
                real(kind=kreal), pointer :: EPSTRESS_VECT(:,:)   !< elemental principal stress vector
                real(kind=kreal), pointer :: EPSTRAIN_VECT(:,:)   !< elemental principal strain vector

                type(fstr_solid_physic_val), pointer :: LAYER(:)    !< Laminated Shell's layer (1,2,3,4,5,...)
                type(fstr_solid_physic_val), pointer :: PLUS    !< for SHELL PLUS
                type(fstr_solid_physic_val), pointer :: MINUS    !< for SHELL MINUS

        end type fstr_solid_physic_val

        type fstr_solid

                integer(kind=kint) :: file_type  ! kbcfFSTR or kbcfNASTRAN
                integer(kind=kint) :: StaticType ! 1:Total, 2:Updated, 3:Infinite
!
                integer(kind=kint) :: nstep_tot
                type(step_info), pointer       :: step_ctrl(:)  =>null()   !< step information
                type(t_output_ctrl),  pointer  :: output_ctrl(:)=>null()   !< output  information

                !!BOUNDARY
                integer(kind=kint) :: BOUNDARY_ngrp_tot                    !< Following boundary conditions
                integer(kind=kint), pointer :: BOUNDARY_ngrp_GRPID  (:)  =>null()
                integer(kind=kint), pointer :: BOUNDARY_ngrp_ID     (:)  =>null()
                integer(kind=kint), pointer :: BOUNDARY_ngrp_type   (:)  =>null()
!               integer(kind=kint), pointer :: BOUNDARY_ngrp_iftype (:)  =>null()
                integer(kind=kint), pointer :: BOUNDARY_ngrp_amp    (:)  =>null()
                real(kind=kreal), pointer   :: BOUNDARY_ngrp_val    (:)  =>null()
                integer(kind=kint) :: BOUNDARY_ngrp_rot                   !< number of rotational boundary conditions
                integer(kind=kint), pointer :: BOUNDARY_ngrp_rotID     (:) =>null()
                integer(kind=kint), pointer :: BOUNDARY_ngrp_centerID  (:) =>null()

                !!VELOCITY
                integer(kind=kint) :: VELOCITY_type
                integer(kind=kint) :: VELOCITY_ngrp_tot                    !< Following velocity boundary condition
                integer(kind=kint), pointer :: VELOCITY_ngrp_GRPID  (:)  =>null()
                integer(kind=kint), pointer :: VELOCITY_ngrp_ID     (:)  =>null()
                integer(kind=kint), pointer :: VELOCITY_ngrp_type   (:)  =>null()
                integer(kind=kint), pointer :: VELOCITY_ngrp_amp    (:)  =>null()
                real(kind=kreal), pointer   :: VELOCITY_ngrp_val    (:)  =>null()

                !!ACCELERATION
                integer(kind=kint) :: ACCELERATION_type
                integer(kind=kint) :: ACCELERATION_ngrp_tot                !< Following accelerate boundary condition
                integer(kind=kint), pointer :: ACCELERATION_ngrp_GRPID  (:)  =>null()
                integer(kind=kint), pointer :: ACCELERATION_ngrp_ID     (:)  =>null()
                integer(kind=kint), pointer :: ACCELERATION_ngrp_type   (:)  =>null()
                integer(kind=kint), pointer :: ACCELERATION_ngrp_amp    (:)  =>null()
                real(kind=kreal), pointer   :: ACCELERATION_ngrp_val    (:)  =>null()

                !!CLOAD
                integer(kind=kint) :: CLOAD_ngrp_tot                       !< Following concetrated external load
                integer(kind=kint), pointer :: CLOAD_ngrp_GRPID     (:) =>null()
                integer(kind=kint), pointer :: CLOAD_ngrp_ID        (:)
                integer(kind=kint), pointer :: CLOAD_ngrp_DOF       (:)
                integer(kind=kint), pointer :: CLOAD_ngrp_amp       (:)
                real(kind=kreal), pointer   :: CLOAD_ngrp_val       (:)
                integer(kind=kint) :: CLOAD_ngrp_rot                   !< number of torque load conditions
                integer(kind=kint), pointer :: CLOAD_ngrp_rotID     (:) =>null()
                integer(kind=kint), pointer :: CLOAD_ngrp_centerID  (:) =>null()

                !!DLOAD
                integer(kind=kint) :: DLOAD_ngrp_tot                       !< Following distrubuted external load
                integer(kind=kint) :: DLOAD_follow
                integer(kind=kint), pointer :: DLOAD_ngrp_GRPID     (:) =>null()
                integer(kind=kint), pointer :: DLOAD_ngrp_ID        (:)
                integer(kind=kint), pointer :: DLOAD_ngrp_LID       (:)
                integer(kind=kint), pointer :: DLOAD_ngrp_amp       (:)
                real(kind=kreal), pointer   :: DLOAD_ngrp_params    (:,:)

                !!TEMPERATURE
                integer(kind=kint) :: TEMP_ngrp_tot                        !< Following tempearture conditions
                integer(kind=kint) :: TEMP_irres
                integer(kind=kint) :: TEMP_tstep
                integer(kind=kint) :: TEMP_interval
                real(kind=kreal)   :: TEMP_FACTOR
                integer(kind=kint), pointer :: TEMP_ngrp_GRPID     (:) =>null()
                integer(kind=kint), pointer :: TEMP_ngrp_ID        (:)
                real(kind=kreal), pointer   :: TEMP_ngrp_val       (:)

                !!SPRING
                integer(kind=kint) :: SPRING_ngrp_tot                      !< Following spring boundary conditions
                integer(kind=kint), pointer :: SPRING_ngrp_GRPID    (:) =>null()
                integer(kind=kint), pointer :: SPRING_ngrp_ID       (:)
                integer(kind=kint), pointer :: SPRING_ngrp_DOF      (:)
                integer(kind=kint), pointer :: SPRING_ngrp_amp      (:)
                real(kind=kreal), pointer   :: SPRING_ngrp_val      (:)

                ! for couple analysis
                integer( kind=kint ) :: COUPLE_ngrp_tot                   !< Following for coupling analysis
                integer( kind=kint ),pointer :: COUPLE_ngrp_ID(:)

                ! VALUE
                real(kind=kreal), pointer :: STRESS(:)    !< nodal stress
                real(kind=kreal), pointer :: STRAIN(:)    !< nodal strain
                real(kind=kreal), pointer :: MISES(:)     !< nodal MISES

                real(kind=kreal), pointer :: PSTRESS(:)   !< nodal principal stress
                real(kind=kreal), pointer :: PSTRAIN(:)   !< nodal principal strain
                real(kind=kreal), pointer :: PSTRESS_VECT(:,:)   !< nodal principal stress vector
                real(kind=kreal), pointer :: PSTRAIN_VECT(:,:)   !< nodal principal strain vector

                real(kind=kreal), pointer :: ESTRESS(:)   !< elemental stress
                real(kind=kreal), pointer :: ESTRAIN(:)   !< elemental strain
                real(kind=kreal), pointer :: EMISES(:)    !< elemental MISES

                real(kind=kreal), pointer :: EPSTRESS(:)   !< elemental principal stress
                real(kind=kreal), pointer :: EPSTRAIN(:)   !< elemental principal strain
                real(kind=kreal), pointer :: EPSTRESS_VECT(:,:)   !< elemental principal stress vector
                real(kind=kreal), pointer :: EPSTRAIN_VECT(:,:)   !< elemental principal strain vector

                real(kind=kreal), pointer :: TNSTRAIN(:)   !< thermal nodal strain
                real(kind=kreal), pointer :: TESTRAIN(:)   !< thermal elemental strain

                type(fstr_solid_physic_val), pointer       :: SOLID=>null()     !< for solid physical value stracture
                type(fstr_solid_physic_val), pointer       :: SHELL=>null()     !< for shell physical value stracture
                type(fstr_solid_physic_val), pointer       :: BEAM =>null()     !< for beam physical value stracture

                ! ANALYSIS CONTROL for NLGEOM
                integer(kind=kint) :: restart_nout   !< output interval of restart file
                                                     !< (if  .gt.0) restart file write
                                                     !< (if  .lt.0) restart file read and write
                integer(kind=kint) :: restart_nin    !< input number of restart

                integer(kind=kint) :: max_lyr        !< maximun num of layer
                integer(kind=kint) :: is_33shell
                integer(kind=kint) :: is_33beam
                integer(kind=kint) :: is_heat
                integer(kind=kint), pointer :: is_rot(:)
                integer(kind=kint) :: elemopt361

                real(kind=kreal)          :: FACTOR     (2)      !< factor of incrementation
                                                                 !< 1:time t  2: time t+dt
                ! for increment control
                integer(kind=kint) :: NRstat_i(10)     !< statistics of newton iteration (integer)
                real(kind=kreal)   :: NRstat_r(10)     !< statistics of newton iteration (real)
                integer(kind=kint) :: AutoINC_stat     !< status of auto-increment control
                integer(kind=kint) :: CutBack_stat     !< status of cutback control

                real(kind=kreal), pointer :: GL          (:)           !< exnternal force
                real(kind=kreal), pointer :: EFORCE      (:)           !< exnternal force
                real(kind=kreal), pointer :: QFORCE      (:)           !< equivalent nodal force
                real(kind=kreal), pointer :: unode(:)      => null()   !< disp at the beginning of curr step
                real(kind=kreal), pointer :: dunode(:)     => null()   !< curr total disp
                real(kind=kreal), pointer :: ddunode(:)    => null()   !< =hecMESH%X, disp icrement
                real(kind=kreal), pointer :: temperature(:)=> null()   !< =temperature
                real(kind=kreal), pointer :: temp_bak(:)   => null()
                real(kind=kreal), pointer :: last_temp(:)  => null()

                type( tElement ), pointer :: elements(:)   =>null()  !< elements information
                type( tMaterial ),pointer :: materials(:)  =>null()  !< material properties
                type( tContact ), pointer :: contacts(:)   =>null()  !< contact information
                integer                   :: n_fix_mpc               !< number mpc conditions user defined
                real(kind=kreal), pointer :: mpc_const(:)  =>null()  !< bakeup of hecmwST_mpc%mpc_const
                type(tSection), pointer   :: sections(:)   =>null()  !< definition of setion referred by elements(i)%sectionID

                ! for cutback
                ! ####################### Notice #######################
                ! # If you add new variables to store analysis status, #
                ! # - backup variables with postfix "_bkup" here       #
                ! # - backup process to module m_fstr_Cutback          #
                ! # must be added if necessary.                        #
                ! ######################################################
                real(kind=kreal), pointer :: unode_bkup(:)     => null() !< disp at the beginning of curr step (backup)
                real(kind=kreal), pointer :: QFORCE_bkup(:)    => null() !< equivalent nodal force (backup)
                real(kind=kreal), pointer :: last_temp_bkup(:) => null()
                type( tElement ), pointer :: elements_bkup(:)  =>null()  !< elements information (backup)
                type( tContact ), pointer :: contacts_bkup(:)  =>null()  !< contact information (backup)
!
        end type fstr_solid
!
!C
!C ----------------------------------------------------------------------------
!C
!> Data for HEAT ANSLYSIS  (fstrHEAT)
!C
        type fstr_heat

                ! TIME CONTROL
                real(kind=kreal) :: TEMTOL
                integer :: STEPtot
                integer(kind=kint) :: restart_nout
                real(kind=kreal), pointer :: STEP_DLTIME(:),STEP_EETIME(:)
                real(kind=kreal), pointer :: STEP_DELMIN(:),STEP_DELMAX(:)

                ! MATERIAL
                integer :: MATERIALtot
                real(kind=kreal), pointer :: RHO(:,:), RHOtemp (:,:)
                real(kind=kreal), pointer :: CP(:,:),  CPtemp (:,:)
                real(kind=kreal), pointer :: COND(:,:),CONDtemp (:,:)

                integer, pointer :: RHOtab(:), CPtab(:), CONDtab(:)

                real(kind=kreal), pointer :: RHOfuncA(:,:),  RHOfuncB(:,:)
                real(kind=kreal), pointer :: CPfuncA (:,:),  CPfuncB(:,:)
                real(kind=kreal), pointer :: CONDfuncA (:,:),CONDfuncB(:,:)

                ! AMPLITUDE
                integer :: AMPLITUDEtot
                real(kind=kreal), pointer :: AMPL(:,:), AMPLtime (:,:)
                integer, pointer :: AMPLtab(:)
                real(kind=kreal), pointer :: AMPLfuncA(:,:),  AMPLfuncB(:,:)


                ! VALUE
                real(kind=kreal), pointer :: TEMP0(:)
                real(kind=kreal), pointer :: TEMPC(:)
                real(kind=kreal), pointer :: TEMP (:)
                real(kind=kreal), pointer :: TEMPW(:)

                ! Residual
                real(kind=kreal), pointer :: RE(:)
                real(kind=kreal), pointer :: QV(:)
                real(kind=kreal), pointer :: RR(:)
                real(kind=kreal), pointer :: RL(:)
                real(kind=kreal), pointer :: RU(:)
                real(kind=kreal), pointer :: RD(:)
                real(kind=kreal), pointer :: IWKX(:,:)

                ! BOUNDARY CONDTIONS -------

                !!FIXTEMP
                integer :: T_FIX_tot
                integer, pointer          :: T_FIX_node(:)
                integer, pointer          :: T_FIX_ampl(:)
                real(kind=kreal), pointer :: T_FIX_val(:)

                !!CFLUX
                integer :: Q_NOD_tot
                integer, pointer          :: Q_NOD_node(:)
                integer, pointer          :: Q_NOD_ampl(:)
                real(kind=kreal), pointer :: Q_NOD_val(:)

                !!DFLUX (not used)
                integer :: Q_VOL_tot
                integer, pointer          :: Q_VOL_elem(:)
                integer, pointer          :: Q_VOL_ampl(:)
                real(kind=kreal), pointer :: Q_VOL_val(:)

                !!DFLUX, !SFLUX
                integer :: Q_SUF_tot
                integer, pointer          :: Q_SUF_elem(:)
                integer, pointer          :: Q_SUF_ampl(:)
                integer, pointer          :: Q_SUF_surf(:)
                real(kind=kreal), pointer :: Q_SUF_val(:)

                !!RADIATE, !SRADIATE
                integer :: R_SUF_tot
                integer, pointer          :: R_SUF_elem(:)
                integer, pointer          :: R_SUF_ampl(:,:)
                integer, pointer          :: R_SUF_surf(:)
                real(kind=kreal), pointer :: R_SUF_val(:,:)

                !!FILM, SFILM
                integer :: H_SUF_tot
                integer, pointer          :: H_SUF_elem(:)
                integer, pointer          :: H_SUF_ampl(:,:)
                integer, pointer          :: H_SUF_surf(:)
                real(kind=kreal), pointer :: H_SUF_val(:,:)

                integer  :: WL_tot
                type(tWeldLine), pointer :: weldline(:) => null()

        end type fstr_heat


!C ----------------------------------------------------------------------------
!C
!> Data for DYNAMIC ANSLYSIS  (fstrDYNAMIC)
!C
        type fstr_dynamic

                ! ANALYSIS TYPE CONTROL
                integer(kind=kint) :: idx_eqa       ! implicit or explicit
                integer(kind=kint) :: idx_resp      ! time history or steady-state harmonic response analysis

                ! TIME CONTROL
                integer(kind=kint) :: n_step        ! total step number of analysis
                real(kind=kreal)   :: t_start       ! start time of analysis
                real(kind=kreal)   :: t_curr        ! current time of analysis
                real(kind=kreal)   :: t_end         ! end time of analysis
                real(kind=kreal)   :: t_delta       ! time increment
                integer(kind=kint) :: restart_nout  ! output interval of restart file
                                                    ! (if  .gt.0) restart file write
                                                    ! (if  .lt.0) restart file read and write
                integer(kind=kint) :: restart_nin  !input number of restart file

                ! Newmark-beta parameter
                real(kind=kreal)   :: ganma         ! Newmark-beta parameter ganma
                real(kind=kreal)   :: beta          ! Newmark-beta parameter beta

                ! mass matrix control
                integer(kind=kint) :: idx_mas       ! mass matrix type

                ! damping control
                integer(kind=kint) :: idx_dmp      ! damping type
                real(kind=kreal)   :: ray_m         ! Rayleigh damping parameter Rm
                real(kind=kreal)   :: ray_k         ! Rayleigh damping parameter Rk

                ! initialization control
                logical           :: VarInitialize  ! initialization flag

                ! OUTPUT CONTROL
                integer(kind=kint) :: nout           ! output interval of result
                integer(kind=kint) :: ngrp_monit     ! node of monitoring result
                integer(kind=kint) :: nout_monit     ! output interval of result monitoring
                integer(kind=kint) :: i_step         ! step number
                integer(kind=kint) :: iout_list(6)   ! 0:not output  1:output
                                                     ! iout_list(1): displacement
                                                     ! iout_list(2): velocity
                                                     ! iout_list(3): acceleration
                                                     ! iout_list(4): reaction force
                                                     ! iout_list(5): strain
                                                     ! iout_list(6): stress

                ! VALUE
                real(kind=kreal), pointer :: DISP  (:,:)     !> Displacement, U(t+dt), U(t), U(t-dt)
                real(kind=kreal), pointer :: VEL   (:,:)     !> Velocity
                real(kind=kreal), pointer :: ACC   (:,:)     !> Acceleration

                real(kind=kreal) :: kineticEnergy
                real(kind=kreal) :: strainEnergy
                real(kind=kreal) :: totalEnergy

                ! temporary quantity
                real(kind=kreal), pointer :: VEC1  (:)
                real(kind=kreal), pointer :: VEC2  (:)
                real(kind=kreal), pointer :: VEC3  (:)

                integer(kind=kint) :: dynamic_IW4        =   204
                integer(kind=kint) :: dynamic_IW5        =   205
                integer(kind=kint) :: dynamic_IW6        =   206
                integer(kind=kint) :: dynamic_IW7        =   207
                integer(kind=kint) :: dynamic_IW8        =   208
                integer(kind=kint) :: dynamic_IW9        =   209

        end type fstr_dynamic
!C ----------------------------------------------------------------------------
!C
!> Data for coupling analysis
!C
    type fstr_couple
        ! for revocap
        integer( kind=kint )   :: dof              ! == 3
        integer( kind=kint )   :: ndof             ! total dof (coupled_node_n*dof)
        integer( kind=kint )   :: coupled_node_n
        ! note) following types depend on revocap
        integer, pointer       :: coupled_node(:)  ! local node id sent from revocap
        real( kind=8 ),pointer :: trac(:)          ! input  (x,y,z,x,y,z ... )
        real( kind=8 ),pointer :: disp(:)          ! output (x,y,z,x,y,z ... )
        real( kind=8 ),pointer :: velo(:)          ! output (x,y,z,x,y,z ... )
        real( kind=8 ),pointer :: accel(:)          ! output (x,y,z,x,y,z ... )
        ! for inner use
        integer( kind=kint ),pointer :: index(:)   ! size:total node num.
        !  -1:not relation, >1:index of coupled_node
    end type fstr_couple

!C ----------------------------------------------------------------------------
!C
!> Data for weld line
!C
    type tWeldLine
        integer              :: egrpid
        real( kind=kreal )   :: I
        real( kind=kreal )   :: U
        real( kind=kreal )   :: coe
        real( kind=kreal )   :: v
        integer              :: xyz
        real(kind=kreal)     :: n1, n2
        real(kind=kreal)     :: distol
        real(kind=kreal)     :: tstart
    end type tWeldLine

!C ----------------------------------------------------------------------------
!> Data for section control
!C
    type tSection
        !integer              :: mat_ID
        !integer              :: iset
        !integer              :: orien_ID
        !real(kind=kreal)     :: thickness
        !integer              :: elemopt341
        !integer              :: elemopt342
        !integer              :: elemopt351
        !integer              :: elemopt352
        integer              :: elemopt361
        !integer              :: elemopt362
    end type tSection

!C ----------------------------------------------------------------------------
!C-- GROBAL VARIABLES for MEMORY MANAGEMENT or TIME RECORDING
!C

        ! for memory management
        integer :: masbr,mlczi,mlczr,mreorr
        logical :: minit

        ! for time recording (for subroutine time_log )
        logical :: tinit,teachiter, tenditer
        real(kind=kint) :: tmaints,tmainte,tmainus,tmainue,tmainss,tmainse
        real(kind=kint) :: tasbts,tasbte,tasbus,tasbue,tasbss,tasbse
        real(kind=kint) :: tlczts,tlczte,tlczus,tlczue,tlczss,tlczse
        real(kind=kint) :: tlczavt,tlczavu,tlczavs
        real(kind=kint) :: treorts,treorte,treorus,treorue,treorss,treorse
        real(kind=kint) :: treavt,treavu,treavs
        real(kind=kint) :: tsolts,tsolte,tsolus,tsolue,tsolss,tsolse
        real(kind=kint) :: tsolavt,tsolavu,tsolavs

        integer :: idummy
        integer(kind=kint) :: iexit  ! for eigen ?

contains

!C ----------------------------------------------------------------------------
!C NULL POINTER SETTING TO AVOID RUNTIME ERROR
!C ----------------------------------------------------------------------------

        subroutine fstr_nullify_fstr_param( P )
        implicit none
        type( fstr_param ) :: P
        nullify( P%dtime )
        nullify( P%etime )
        nullify( P%dtmin )
        nullify( P%delmax )
        nullify( P%itmax )
        nullify( P%eps )
        nullify( P%global_local_ID)
        nullify( P%timepoints )
        end subroutine fstr_nullify_fstr_param
!C ----------------------------------------------------------------------------
        subroutine fstr_nullify_fstr_solid( S )
        implicit none
        type( fstr_solid ) :: S
        nullify( S%BOUNDARY_ngrp_ID )
        nullify( S%BOUNDARY_ngrp_type )
!       nullify( S%BOUNDARY_ngrp_iftype )
        nullify( S%BOUNDARY_ngrp_amp )
        nullify( S%BOUNDARY_ngrp_val)
        nullify( S%BOUNDARY_ngrp_rotID )
        nullify( S%BOUNDARY_ngrp_centerID )
        nullify( S%CLOAD_ngrp_ID )
        nullify( S%CLOAD_ngrp_DOF )
        nullify( S%CLOAD_ngrp_amp )
        nullify( S%CLOAD_ngrp_rotID )
        nullify( S%CLOAD_ngrp_centerID )
        nullify( S%CLOAD_ngrp_val )
        nullify( S%DLOAD_ngrp_ID )
        nullify( S%DLOAD_ngrp_LID )
        nullify( S%DLOAD_ngrp_amp )
        nullify( S%DLOAD_ngrp_params )
        nullify( S%TEMP_ngrp_ID )
        nullify( S%TEMP_ngrp_val )
        nullify( S%SPRING_ngrp_ID )
        nullify( S%SPRING_ngrp_DOF )
        nullify( S%SPRING_ngrp_amp )
        nullify( S%SPRING_ngrp_val )
        nullify( S%STRESS )
        nullify( S%STRAIN )
        nullify( S%MISES )
        nullify( S%ESTRESS )
        nullify( S%ESTRAIN )
        nullify( S%EMISES )

        !        nullify( S%ESECTSTRESS )
        nullify( S%GL          )
!        nullify( S%TOTAL_DISP  )
        nullify( S%QFORCE      )

        nullify( S%VELOCITY_ngrp_ID )
        nullify( S%VELOCITY_ngrp_type )
        nullify( S%VELOCITY_ngrp_amp )
        nullify( S%VELOCITY_ngrp_val )
        nullify( S%ACCELERATION_ngrp_ID )
        nullify( S%ACCELERATION_ngrp_type )
        nullify( S%ACCELERATION_ngrp_amp )
        nullify( S%ACCELERATION_ngrp_val )

    ! for couple analysis
        nullify( S%COUPLE_ngrp_ID )

        end subroutine fstr_nullify_fstr_solid

!C ----------------------------------------------------------------------------
        subroutine fstr_nullify_fstr_heat( H )
        implicit none
        type( fstr_heat ) :: H

        nullify( H%STEP_DLTIME )
        nullify( H%STEP_EETIME )
        nullify( H%STEP_DELMIN )
        nullify( H%STEP_DELMAX )
        nullify( H%RHO )
        nullify( H%RHOtemp  )
        nullify( H%CP )
        nullify( H%CPtemp  )
        nullify( H%COND )
        nullify( H%CONDtemp  )
        nullify( H%RHOtab )
        nullify( H%CPtab )
        nullify( H%CONDtab )
        nullify( H%RHOfuncA )
        nullify( H%RHOfuncB )
        nullify( H%CPfuncA  )
        nullify( H%CPfuncB )
        nullify( H%CONDfuncA  )
        nullify( H%CONDfuncB )
        nullify( H%AMPL )
        nullify( H%AMPLtime  )
        nullify( H%AMPLtab )
        nullify( H%AMPLfuncA )
        nullify( H%AMPLfuncB )
        nullify( H%TEMP0 )
        nullify( H%TEMPC )
        nullify( H%TEMP  )
        nullify( H%TEMPW )
        nullify( H%RE )
        nullify( H%QV )
        nullify( H%RR )
        nullify( H%RL )
        nullify( H%RU )
        nullify( H%RD )
        nullify( H%IWKX )
        nullify( H%T_FIX_node )
        nullify( H%T_FIX_ampl )
        nullify( H%T_FIX_val )
        nullify( H%Q_NOD_node )
        nullify( H%Q_NOD_ampl )
        nullify( H%Q_NOD_val )
        nullify( H%Q_VOL_elem )
        nullify( H%Q_VOL_ampl )
        nullify( H%Q_VOL_val )
        nullify( H%Q_SUF_elem )
        nullify( H%Q_SUF_ampl )
        nullify( H%Q_SUF_surf )
        nullify( H%Q_SUF_val )
        nullify( H%R_SUF_elem )
        nullify( H%R_SUF_ampl )
        nullify( H%R_SUF_surf )
        nullify( H%R_SUF_val )
        nullify( H%H_SUF_elem )
        nullify( H%H_SUF_ampl )
        nullify( H%H_SUF_surf )
        nullify( H%H_SUF_val )

        end subroutine fstr_nullify_fstr_heat
!C ----------------------------------------------------------------------------
        subroutine fstr_nullify_fstr_dynamic( DY )
        implicit none
        type( fstr_dynamic ) :: DY

!!         nullify( DY%iout_list )
         nullify( DY%DISP )
         nullify( DY%VEL  )
         nullify( DY%ACC  )
         nullify( DY%VEC1 )
         nullify( DY%VEC2 )
         nullify( DY%VEC3 )

        end subroutine fstr_nullify_fstr_dynamic
!C ----------------------------------------------------------------------------
        subroutine fstr_nullify_fstr_couple( C )
        type( fstr_couple ) :: C
        nullify( C%coupled_node )
        nullify( C%trac )
        nullify( C%velo )
        nullify( C%index )
        end subroutine fstr_nullify_fstr_couple

!> Initializer of structure hecmwST_matrix
subroutine fstr_mat_init( hecMAT )
        implicit none
        type(hecmwST_matrix) :: hecMAT

        hecMAT%Iarray(1) =  100    ! = nier
        hecMAT%Iarray(2) =    1    ! = method
        hecMAT%Iarray(3) =    1    ! = precond
        hecMAT%Iarray(4) =    0    ! = nset
        hecMAT%Iarray(5) =    1    ! = iterpremax
        hecMAT%Iarray(6) =   10    ! = nrest
        hecMAT%Iarray(7) =    0    ! = scaling
        hecMAT%Iarray(21)=  kNO    ! = iterlog
        hecMAT%Iarray(22)=  kNO    ! = timelog
        hecMAT%Iarray(31)=    0    ! = dumptype
        hecMAT%Iarray(32)=    0    ! = dumpexit
        hecMAT%Iarray(33)=    0    ! = usejad
        hecMAT%Iarray(34)=   10    ! = ncolor_in
        hecMAT%Iarray(13)=    3    ! = mpc_method
        hecMAT%Iarray(14)=    0    ! = estcond
        hecMAT%Iarray(35)=    3    ! = maxrecycle_precond

        hecMAT%Rarray(1) =  1.0e-8 ! = resid
        hecMAT%Rarray(2) =  1.0    ! = sigma_diag
        hecMAT%Rarray(3) =  0.0    ! = sigma
        hecMAT%Rarray(4) =  0.1    ! = thresh
        hecMAT%Rarray(5) =  0.1    ! = filter
        hecMAT%Rarray(11)=  1.0e+4 ! = penalty

        hecMAT%Iarray(96) =   0    ! nrecycle_precond
        hecMAT%Iarray(97) = kYES   ! flag_numfact
        hecMAT%Iarray(98) = kYES   ! flag_symbfact
        hecMAT%Iarray(99) = kYES   ! indirect method

end subroutine fstr_mat_init

!C ----------------------------------------------------------------------------
        subroutine hecMAT_init( hecMAT )
! ----------------------------------------------------------------------------
!  Purpose: Memeory allocation
!           July.6, 2009  YUAN Xi
!  Notes: 1. Should be a function of middleware !
!         2. Not available for mixed-dof problems
!  This subroutine should completely rewritten next year
! ----------------------------------------------------------------------------
        type( hecmwST_matrix ) :: hecMAT
        integer ::  ndof, nn, ierror
        ndof = hecMAT%NDOF
        nn = ndof*ndof
        allocate (hecMAT%AL(nn*hecMAT%NPL)        ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%AU(nn*hecMAT%NPU)        ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%B(ndof*hecMAT%NP)          ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        hecMAT%B(:)=0.d0
        allocate (hecMAT%D(nn*hecMAT%NP)          ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%X(ndof*hecMAT%NP)          ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        end if
        allocate (hecMAT%ALU(nn*hecMAT%N)         ,STAT=ierror )
        if( ierror /= 0 ) then
            write(*,*) "##ERROR : not enough memory"
            write(idbg,*) 'stop due to allocation error'
            call flush(idbg)
            call hecmw_abort( hecmw_comm_get_comm() )
        endif
        call hecmw_cmat_init( hecMAT%cmat )
        end subroutine hecMAT_init
! ----------------------------------------------------------------------------
        subroutine hecMAT_finalize( hecMAT )
! ----------------------------------------------------------------------------
!  Purpose: Memeory allocation
!           July.6, 2009  YUAN Xi
!  Notes: 1. Should be a function of middleware !
!         2. Not available for mixed-dof problems
!  This subroutine should completely rewritten next year
! ----------------------------------------------------------------------------
        type( hecmwST_matrix ) :: hecMAT
        integer ::  ndof, nn, ierror
        ndof = hecMAT%NDOF
        nn = ndof*ndof
        if( associated(hecMAT%AL) ) then
            deallocate(hecMAT%AL                  ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%AU) ) then
            deallocate(hecMAT%AU                  ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%B) ) then
            deallocate(hecMAT%B                   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
             call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%D) ) then
            deallocate(hecMAT%D                   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(HECMAT%X) ) then
            deallocate(hecMAT%X                   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        if( associated(hecMAT%ALU) ) then
            deallocate(hecMAT%ALU                 ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error'
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        endif
        call hecmw_cmat_finalize(hecmAT%cmat)
        end subroutine hecMAT_finalize

!> Initializer of structure fstr_param
subroutine fstr_param_init( fstrPARAM, hecMESH )
        implicit none
        type(fstr_param) :: fstrPARAM
        type(hecmwST_local_mesh) :: hecMESH
        integer(kind=kint) :: i
        external fstr_sort_index

        fstrPARAM%solution_type = kstSTATIC
        fstrPARAM%nlgeom        = .false.
        fstrPARAM%solver_method = ksmCG

        !!STATIC !HEAT
        fstrPARAM%analysis_n = 0
        fstrPARAM%ref_temp   = 0

        ! output control
        fstrPARAM%fg_echo       = kOFF
        fstrPARAM%fg_result     = kOFF
        fstrPARAM%fg_visual     = kOFF

        ! for heat ...
        fstrPARAM%fg_neutral    = kOFF
        fstrPARAM%fg_irres      = kNO
        fstrPARAM%fg_iwres      = kNO
        fstrPARAM%nrres         = 1
        fstrPARAM%nprint        = 1

        ! for couple
        fstrPARAM%fg_couple      = 0
        fstrPARAM%fg_couple_type = 0
        fstrPARAM%fg_couple_first= 0
        fstrPARAM%fg_couple_window= 0

        ! for restart control
        fstrPARAM%restart_version = 5

        ! for auto increment and cutback
        fstrPARAM%ainc_Rs   = 0.25d0
        fstrPARAM%ainc_Rl   = 1.25d0
        fstrPARAM%NRbound_s = 0
        fstrPARAM%NRbound_s(knstMAXIT) = 10
        fstrPARAM%NRbound_s(knstSUMIT) = 50
        fstrPARAM%NRbound_s(knstCITER) = 10
        fstrPARAM%NRbound_l = 0
        fstrPARAM%NRbound_l(knstMAXIT) = 1
        fstrPARAM%NRbound_l(knstSUMIT) = 1
        fstrPARAM%NRbound_l(knstCITER) = 1
        fstrPARAM%NRtimes_s = 1
        fstrPARAM%NRtimes_l = 2
        fstrPARAM%ainc_Rc   = 0.25d0
        fstrPARAM%CBbound   = 5

        ! index table for global node ID sorting

        fstrPARAM%n_node = hecMESH%n_node;
        fstrPARAM%nn_internal = hecMESH%nn_internal;
        allocate( fstrPARAM%global_local_ID(2,hecMESH%nn_internal))
        do i = 1, hecMESH%nn_internal
                fstrPARAM%global_local_ID(1,i) = hecMESH%global_node_ID(i)
                fstrPARAM%global_local_ID(2,i) = i
        end do
        call fstr_sort_index( fstrPARAM%global_local_ID, hecMESH%nn_internal )

end subroutine fstr_param_init


        logical function fstr_isBoundaryActive( fstrSOLID, nbc, cstep )
          type(fstr_solid)    :: fstrSOLID
          integer, intent(in) :: nbc    !< group id of boundary condition
          integer, intent(in) :: cstep  !< current step number
          fstr_isBoundaryActive = .true.
          if( .not. associated(fstrSOLID%step_ctrl) ) return
          if( cstep>fstrSOLID%nstep_tot ) return
          fstr_isBoundaryActive = isBoundaryActive( nbc, fstrSOLID%step_ctrl(cstep) )
        end function
!
        logical function fstr_isLoadActive( fstrSOLID, nbc, cstep )
          type(fstr_solid)    :: fstrSOLID
          integer, intent(in) :: nbc    !< group id of boundary condition
          integer, intent(in) :: cstep  !< current step number
          fstr_isLoadActive = .true.
          if( .not. associated(fstrSOLID%step_ctrl) ) return
          if( cstep>fstrSOLID%nstep_tot ) return
          fstr_isLoadActive = isLoadActive( nbc, fstrSOLID%step_ctrl(cstep) )
        end function

        logical function fstr_isContactActive( fstrSOLID, nbc, cstep )
          type(fstr_solid)     :: fstrSOLID !< type fstr_solid
          integer, intent(in) :: nbc       !< group id of boundary condition
          integer, intent(in) :: cstep     !< current step number
          fstr_isContactActive = .true.
          if( .not. associated(fstrSOLID%step_ctrl) ) return
          if( cstep>fstrSOLID%nstep_tot ) return
          fstr_isContactActive = isContactActive( nbc, fstrSOLID%step_ctrl(cstep) )
        end function

    !< This subroutine fetch coords defined by local coordinate system
        subroutine get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )
           integer, intent(in)             :: cdsys_ID      !< id of local coordinate
           type (hecmwST_local_mesh)       :: hecMESH       !< mesh information
           type (fstr_solid),intent(inout) :: fstrSOLID     !< fstr_solid
           real(kind=kreal), intent(out)   :: coords(3,3)

           integer :: ik

           coords =0.d0
           if( cdsys_ID>0 ) then
                 if( isCoordNeeds(g_LocalCoordSys(cdsys_ID)) ) then
                    coords=g_LocalCoordSys(cdsys_ID)%CoordSys
                 else
                    ik=g_LocalCoordSys(cdsys_ID)%node_ID(1)
                    coords(1,:)= hecMESH%node(3*ik-2:3*ik)+fstrSOLID%unode(3*ik-2:3*ik)  &
                                + fstrSOLID%dunode(3*ik-2:3*ik)
                    ik=g_LocalCoordSys(cdsys_ID)%node_ID(2)
                    coords(2,:)= hecMESH%node(3*ik-2:3*ik)+fstrSOLID%unode(3*ik-2:3*ik)  &
                                + fstrSOLID%dunode(3*ik-2:3*ik)
                    ik=g_LocalCoordSys(cdsys_ID)%node_ID(3)
                    if(ik>0) coords(3,:)= hecMESH%node(3*ik-2:3*ik)+fstrSOLID%unode(3*ik-2:3*ik)  &
                                + fstrSOLID%dunode(3*ik-2:3*ik)
                 endif
           endif
        end subroutine

        subroutine fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
          implicit none
          type(hecmwST_local_mesh), intent(inout) :: hecMESH
          type (fstr_solid), intent(in) :: fstrSOLID
          real(kind=kreal), pointer :: coord(:)
          integer(kind=kint) :: i
          if(hecMESH%n_dof == 4) return
          do i = 1, hecMESH%n_node*min(hecMESH%n_dof,3)
            coord(i) = hecMESH%node(i)
            hecMESH%node(i) = coord(i)+fstrSOLID%unode(i)+fstrSOLID%dunode(i)
          enddo
        end subroutine fstr_set_current_config_to_mesh


        subroutine fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)
          implicit none
          type(hecmwST_local_mesh), intent(inout) :: hecMESH
          type (fstr_solid), intent(in) :: fstrSOLID
          real(kind=kreal), pointer :: coord(:)
          integer(kind=kint) :: i
          if(hecMESH%n_dof == 4) return
          do i = 1, hecMESH%n_node*min(hecMESH%n_dof,3)
            hecMESH%node(i) = coord(i)
          enddo
        end subroutine fstr_recover_initial_config_to_mesh


subroutine fstr_solid_phys_zero(phys)
  implicit none
  type(fstr_solid_physic_val), pointer :: phys
  phys%STRAIN  = 0.0d0
  phys%STRESS  = 0.0d0
  phys%MISES   = 0.0d0
  phys%ESTRAIN = 0.0d0
  phys%ESTRESS = 0.0d0
  phys%EMISES  = 0.0d0
end subroutine fstr_solid_phys_zero

subroutine fstr_solid_phys_clear(fstrSOLID)

  implicit none
  type (fstr_solid)         :: fstrSOLID
  type(fstr_solid_physic_val), pointer :: phys
  integer(kind=kint) :: i
  
  if (associated(fstrSOLID%SOLID)) then
    call fstr_solid_phys_zero(fstrSOLID%SOLID)
  end if 
  if (associated(fstrSOLID%SHELL)) then
    call fstr_solid_phys_zero(fstrSOLID%SHELL)
    do i=1,fstrSOLID%max_lyr
      call fstr_solid_phys_zero(fstrSOLID%SHELL%LAYER(i)%PLUS)
      call fstr_solid_phys_zero(fstrSOLID%SHELL%LAYER(i)%MINUS)
    end do
  end if 
  if (associated(fstrSOLID%BEAM)) then
    call fstr_solid_phys_zero(fstrSOLID%BEAM)
  end if 
  
end subroutine fstr_solid_phys_clear
end module m_fstr
