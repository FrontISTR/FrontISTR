!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
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
!> \brief  This module defines common data and basic structures for analysis
module m_fstr
  use hecmw
  use m_common_struct
  use m_step
  use m_out
  use m_timepoint
  use mMechGauss
  use mContactDef

  implicit none

  public

  !> CONSTANTS
  !> general
  integer(kind=kint), parameter :: kYES = 1
  integer(kind=kint), parameter :: kNO  = 0
  integer(kind=kint), parameter :: kON  = 1
  integer(kind=kint), parameter :: kOFF = 0

  !> solution type (st)
  integer(kind=kint), parameter :: kstPRECHECK = 0
  integer(kind=kint), parameter :: kstSTATIC   = 1
  integer(kind=kint), parameter :: kstEIGEN    = 2
  integer(kind=kint), parameter :: kstHEAT     = 3
  integer(kind=kint), parameter :: kstDYNAMIC  = 4
  !integer(kind=kint), parameter :: kstNLSTATIC  =   5
  integer(kind=kint), parameter :: kstSTATICEIGEN = 6
  integer(kind=kint), parameter :: kstNZPROF      = 7

  !> solver method (sm)    !CAUTION : (<=100):indirect, (>100):direct
  integer(kind=kint), parameter :: ksmCG       = 1
  integer(kind=kint), parameter :: ksmBiCGSTAB = 2
  integer(kind=kint), parameter :: ksmGMRES    = 3
  integer(kind=kint), parameter :: ksmGPBiCG   = 4
  integer(kind=kint), parameter :: ksmGMRESR   = 5
  integer(kind=kint), parameter :: ksmGMRESREN = 6
  integer(kind=kint), parameter :: ksmDIRECT   = 101

  !> contact analysis algorithm
  integer(kind=kint), parameter :: kcaSLagrange = 1
  integer(kind=kint), parameter :: kcaALagrange = 2

  !> boundary condition file type (bcf)
  integer(kind=kint), parameter :: kbcfFSTR     =   0  ! BC described in fstr control file (default)
  integer(kind=kint), parameter :: kbcfNASTRAN  =   1  ! nastran file

  integer(kind=kint), parameter :: kbcInitial   =  1
  integer(kind=kint), parameter :: kbcTransit   =  2

  !> restart type
  integer(kind=kint), parameter :: restart_outLast = 1
  integer(kind=kint), parameter :: restart_outAll  = 2

  !> section control
  integer(kind=kint), parameter :: kel341FI     =  1
  integer(kind=kint), parameter :: kel341SESNS  =  2

  integer(kind=kint), parameter :: kel361FI     =  1
  integer(kind=kint), parameter :: kel361BBAR   =  2
  integer(kind=kint), parameter :: kel361IC     =  3
  integer(kind=kint), parameter :: kel361FBAR   =  4

  integer(kind=kint), parameter :: kFLOADTYPE_NODE = 1
  integer(kind=kint), parameter :: kFLOADTYPE_SURF = 2

  integer(kind=kint), parameter :: kFLOADCASE_RE = 1
  integer(kind=kint), parameter :: kFLOADCASE_IM = 2

  !> PARALLEL EXECUTION
  integer(kind = kint) :: myrank
  integer(kind = kint) :: nprocs

  !> PARALLEL CONTACT FLAG
  logical :: paraContactFlag = .false.

  !> FILE NAME
  character(len=HECMW_FILENAME_LEN) :: cntfilNAME
  character(len=HECMW_FILENAME_LEN) :: restartfilNAME

  !> FILE HANDLER
  integer(kind=kint), parameter :: ILOG = 16 ! log
  integer(kind=kint), parameter :: ISTA = 17 ! status
  integer(kind=kint), parameter :: IUTB = 18 ! utable
  integer(kind=kint), parameter :: IMSG = 51 ! message (myrank == 0 only)
  integer(kind=kint), parameter :: IDBG = 52 ! debug
  integer(kind=kint), parameter :: IFVS = 53 ! visual.ini file
  integer(kind=kint), parameter :: INEU = 54 ! neutral file (heat)
  integer(kind=kint), parameter :: IRESOUT = 100 ! ~110, keeping for result output file

  !> SOLVER CONTROL
  integer(kind=kint) :: svIarray(100)
  real(kind=kreal)   :: svRarray(100)

  !> FLAG for ECHO/RESULT/POST
  integer(kind=kint), pointer :: IECHO
  integer(kind=kint), pointer :: IRESULT
  integer(kind=kint), pointer :: IVISUAL
  integer(kind=kint), pointer :: INEUTRAL  ! flag for femap neutral file
  integer(kind=kint), pointer :: IRRES     ! flag for restart, read
  integer(kind=kint), pointer :: IWRES     ! flag for restart, write
  integer(kind=kint), pointer :: NRRES     ! position of restart read
  integer(kind=kint), pointer :: NPRINT    ! interval of write

  integer(kind=kint), parameter :: kOPSS_SOLUTION = 1
  integer(kind=kint), parameter :: kOPSS_MATERIAL = 2
  integer(kind=kint)            :: OPSSTYPE = kOPSS_SOLUTION ! output stress/strain type


  !> REFTEMP
  real(kind=kreal), pointer :: REF_TEMP

  !> ANALYSIS CONTROL for NLGEOM and HEAT
  real(kind=kreal)   :: DT    ! /=fstr_param%dtime
  real(kind=kreal)   :: ETIME ! /=fstr_param%etime
  integer(kind=kint) :: ITMAX
  real(kind=kreal)   :: EPS   ! /=fstr_param%eps

  type tInitialCondition
     character(len=HECMW_FILENAME_LEN)          :: cond_name
     integer                    :: node_elem                    !< node =0;  element =1
     integer                    :: grpid
     integer, pointer           :: intval(:)     => null()      !< if -1, not initialized, otherwise dof number
     real(kind=kreal), pointer  :: realval(:)    => null()      !< initial value
  end type
  type( tInitialCondition ), pointer, save :: g_InitialCnd(:) => null()

  !>  FSTR INNER CONTROL PARAMETERS  (fstrPARAM)
  type fstr_param
    integer(kind=kint) :: solution_type !< solution type number
    integer(kind=kint) :: solver_method !< solver method number
    logical            :: nlgeom        !< is geometrical nonlinearity considered

    !> STATIC !HEAT
    integer(kind=kint) :: analysis_n      !< Number of analysis
    real(kind=kreal), pointer  :: dtime(:) !< (=DT)    STEP_DLTIME
    real(kind=kreal), pointer  :: etime(:) !< (/=ETIME) STEP_EETIME
    real(kind=kreal), pointer  :: dtmin(:) !< (=DTMIN) STEP_DELMIN
    real(kind=kreal), pointer  :: delmax(:)!< (=DTMAX) STEP_DELMAX
    integer(kind=kint), pointer:: itmax(:) !< (/=ITMAX)
    real(kind=kreal), pointer  :: eps(:)   !< (/=ESP)
    real(kind=kreal)           :: ref_temp !< (=REF_TEMP)
    integer(kind=kint)         :: timepoint_id !< time point ID for heat analysis

    !> output control
    integer(kind=kint) :: fg_echo       !< output echo   (kYES/kNO) (=IECHO)
    integer(kind=kint) :: fg_result     !< output result (kYES/kNO) (=IRESULT)
    integer(kind=kint) :: fg_visual     !< visualization (kYES/kNO) (=IVISUAL)

    !> for heat ...
    integer(kind=kint) :: fg_neutral    !< write by neutral (=INEUTRAL)
    integer(kind=kint) :: fg_irres      !< restart read     (=IRRES)
    integer(kind=kint) :: fg_iwres      !< restart write    (=IWRES)
    integer(kind=kint) :: nrres         !< NRRES
    integer(kind=kint) :: nprint        !< NPRINT

    !> index table for global node ID sorting
    integer(kind=kint) :: n_node        !< hecMESH%n_node
    integer(kind=kint) :: nn_internal   !< hecMESH%nn_internal
    integer(kind=kint), pointer :: global_local_ID(:,:)  !> (2:nn_internal) (1,:):global, (2,:):local

    !> for couple analysis
    integer( kind=kint ) :: fg_couple          !< (default:0)
    integer( kind=kint ) :: fg_couple_type     !< (default:0)
    integer( kind=kint ) :: fg_couple_first    !< (default:0)
    integer( kind=kint ) :: fg_couple_window   !< (default:0)

    !> for restart control
    integer( kind=kint ) :: restart_out_type   !< output type of restart file
    integer( kind=kint ) :: restart_version    !< version of restart file

    !> for contact analysis
    integer( kind=kint ) :: contact_algo       !< contact analysis algorithm number(SLagrange or Alagrange)
    type(tContactParam), pointer :: contactparam(:)  !< parameter sets for contact scan

    !> for auto increment and cutback
    type(tParamAutoInc), pointer :: ainc(:)        !< auto increment control
    type(time_points), pointer :: timepoints(:)  !< time points data
  end type fstr_param

  !> GLOBAL VARIABLE INITIALIZED IN FSTR_SETUP
  type( fstr_param ),target :: fstrPR

  !> Data for STATIC ANSLYSIS  (fstrSOLID)
  type fstr_solid_physic_val
    real(kind=kreal), pointer :: STRESS(:) => null()   !< nodal stress
    real(kind=kreal), pointer :: STRAIN(:) => null()   !< nodal strain
    real(kind=kreal), pointer :: MISES(:) => null()    !< nodal MISES

    real(kind=kreal), pointer :: PSTRESS(:) => null()   !< nodal principal stress
    real(kind=kreal), pointer :: PSTRAIN(:) => null()   !< nodal principal strain
    real(kind=kreal), pointer :: PSTRESS_VECT(:,:) => null()  !< nodal principal stress vector
    real(kind=kreal), pointer :: PSTRAIN_VECT(:,:) => null()  !< nodal principal strain vector

    real(kind=kreal), pointer :: ESTRESS(:) => null()  !< elemental stress
    real(kind=kreal), pointer :: ESTRAIN(:) => null()  !< elemental strain
    real(kind=kreal), pointer :: EMISES(:) => null()   !< elemental MISES
    real(kind=kreal), pointer :: EPLSTRAIN(:) => null()   !< elemental plastic strain

    real(kind=kreal), pointer :: EPSTRESS(:) => null()  !< elemental principal stress
    real(kind=kreal), pointer :: EPSTRAIN(:) => null()  !< elemental principal strain
    real(kind=kreal), pointer :: EPSTRESS_VECT(:,:) => null()  !< elemental principal stress vector
    real(kind=kreal), pointer :: EPSTRAIN_VECT(:,:) => null()  !< elemental principal strain vector
    real(kind=kreal), pointer :: ENQM(:) => null()     !< elemental NQM


    type(fstr_solid_physic_val), pointer :: LAYER(:) => null()   !< Laminated Shell's layer (1,2,3,4,5,...)
    type(fstr_solid_physic_val), pointer :: PLUS => null()   !< for SHELL PLUS
    type(fstr_solid_physic_val), pointer :: MINUS => null()  !< for SHELL MINUS
  end type fstr_solid_physic_val

  type fstr_solid
    integer(kind=kint) :: file_type  ! kbcfFSTR or kbcfNASTRAN
    integer(kind=kint) :: StaticType ! 1:Total, 2:Updated, 3:Infinitesimal
    integer(kind=kint) :: nstep_tot

    type(step_info), pointer       :: step_ctrl(:)  =>null()   !< step information
    type(t_output_ctrl),  pointer  :: output_ctrl(:)=>null()   !< output  information

    !> BOUNDARY
    integer(kind=kint) :: BOUNDARY_ngrp_tot                    !< Following boundary conditions
    integer(kind=kint), pointer :: BOUNDARY_ngrp_GRPID  (:)  =>null()
    integer(kind=kint), pointer :: BOUNDARY_ngrp_ID     (:)  =>null()
    integer(kind=kint), pointer :: BOUNDARY_ngrp_type   (:)  =>null()
    integer(kind=kint), pointer :: BOUNDARY_ngrp_amp    (:)  =>null()
    real(kind=kreal), pointer   :: BOUNDARY_ngrp_val    (:)  =>null()
    integer(kind=kint), pointer :: BOUNDARY_ngrp_istot  (:)  =>null()
    integer(kind=kint) :: BOUNDARY_ngrp_rot                   !< number of rotational boundary conditions
    integer(kind=kint), pointer :: BOUNDARY_ngrp_rotID     (:) =>null()
    integer(kind=kint), pointer :: BOUNDARY_ngrp_centerID  (:) =>null()

    !> VELOCITY
    integer(kind=kint) :: VELOCITY_type
    integer(kind=kint) :: VELOCITY_ngrp_tot                    !< Following velocity boundary condition
    integer(kind=kint), pointer :: VELOCITY_ngrp_GRPID  (:)  =>null()
    integer(kind=kint), pointer :: VELOCITY_ngrp_ID     (:)  =>null()
    integer(kind=kint), pointer :: VELOCITY_ngrp_type   (:)  =>null()
    integer(kind=kint), pointer :: VELOCITY_ngrp_amp    (:)  =>null()
    real(kind=kreal), pointer   :: VELOCITY_ngrp_val    (:)  =>null()

    !> ACCELERATION
    integer(kind=kint) :: ACCELERATION_type
    integer(kind=kint) :: ACCELERATION_ngrp_tot                !< Following accelerate boundary condition
    integer(kind=kint), pointer :: ACCELERATION_ngrp_GRPID  (:)  =>null()
    integer(kind=kint), pointer :: ACCELERATION_ngrp_ID     (:)  =>null()
    integer(kind=kint), pointer :: ACCELERATION_ngrp_type   (:)  =>null()
    integer(kind=kint), pointer :: ACCELERATION_ngrp_amp    (:)  =>null()
    real(kind=kreal), pointer   :: ACCELERATION_ngrp_val    (:)  =>null()

    !> CLOAD
    integer(kind=kint) :: CLOAD_ngrp_tot                       !< Following concetrated external load
    integer(kind=kint), pointer :: CLOAD_ngrp_GRPID     (:) =>null()
    integer(kind=kint), pointer :: CLOAD_ngrp_ID        (:)
    integer(kind=kint), pointer :: CLOAD_ngrp_DOF       (:)
    integer(kind=kint), pointer :: CLOAD_ngrp_amp       (:)
    real(kind=kreal), pointer   :: CLOAD_ngrp_val       (:)
    integer(kind=kint) :: CLOAD_ngrp_rot                   !< number of torque load conditions
    integer(kind=kint), pointer :: CLOAD_ngrp_rotID     (:) =>null()
    integer(kind=kint), pointer :: CLOAD_ngrp_centerID  (:) =>null()

    !> DLOAD
    integer(kind=kint) :: DLOAD_ngrp_tot                       !< Following distributed external load
    integer(kind=kint) :: DLOAD_follow
    integer(kind=kint), pointer :: DLOAD_ngrp_GRPID     (:) =>null()
    integer(kind=kint), pointer :: DLOAD_ngrp_ID        (:)
    integer(kind=kint), pointer :: DLOAD_ngrp_LID       (:)
    integer(kind=kint), pointer :: DLOAD_ngrp_amp       (:)
    real(kind=kreal), pointer   :: DLOAD_ngrp_params    (:,:)

    !> TEMPERATURE
    integer(kind=kint) :: TEMP_ngrp_tot                        !< Following temperature conditions
    integer(kind=kint) :: TEMP_irres
    integer(kind=kint) :: TEMP_tstep
    integer(kind=kint) :: TEMP_interval
    integer(kind=kint) :: TEMP_rtype      ! type of reading result; 1: step-based; 2: time-based
    real(kind=kreal)   :: TEMP_FACTOR
    integer(kind=kint), pointer :: TEMP_ngrp_GRPID     (:) =>null()
    integer(kind=kint), pointer :: TEMP_ngrp_ID        (:)
    real(kind=kreal), pointer   :: TEMP_ngrp_val       (:)

    !> SPRING
    integer(kind=kint) :: SPRING_ngrp_tot                      !< Following spring boundary conditions
    integer(kind=kint), pointer :: SPRING_ngrp_GRPID    (:) =>null()
    integer(kind=kint), pointer :: SPRING_ngrp_ID       (:)
    integer(kind=kint), pointer :: SPRING_ngrp_DOF      (:)
    integer(kind=kint), pointer :: SPRING_ngrp_amp      (:)
    real(kind=kreal), pointer   :: SPRING_ngrp_val      (:)

    !> for couple analysis
    integer( kind=kint ) :: COUPLE_ngrp_tot                   !< Following for coupling analysis
    integer( kind=kint ),pointer :: COUPLE_ngrp_ID(:)

    !> VALUE
    integer(kind=kint) :: maxn_gauss

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

    real(kind=kreal), pointer :: YIELD_RATIO(:)    !< yield ratio

    real(kind=kreal), pointer :: ENQM(:)      !< elemental NQM
    real(kind=kreal), pointer :: REACTION(:)    !< reaction_force

    real(kind=kreal), pointer :: CONT_NFORCE(:)  !< contact normal force for output
    real(kind=kreal), pointer :: CONT_FRIC(:)    !< contact friction force for output
    real(kind=kreal), pointer :: CONT_RELVEL(:)  !< contact ralative velocity for output
    real(kind=kreal), pointer :: CONT_STATE(:)   !< contact state for output
    integer(kind=kint), pointer :: CONT_SGRP_ID(:) !< contact element surf ids for output
    real(kind=kreal), pointer :: CONT_AREA(:)    !< contact area
    real(kind=kreal), pointer :: CONT_NTRAC(:)   !< contact normal traction force for output
    real(kind=kreal), pointer :: CONT_FTRAC(:)   !< contact friction traction force for output
    real(kind=kreal), pointer :: EMBED_NFORCE(:)  !< embed force for output

    type(fstr_solid_physic_val), pointer :: SOLID=>null()     !< for solid physical value stracture
    type(fstr_solid_physic_val), pointer :: SHELL=>null()     !< for shell physical value stracture
    type(fstr_solid_physic_val), pointer :: BEAM =>null()     !< for beam physical value stracture

    !> ANALYSIS CONTROL for NLGEOM
    integer(kind=kint) :: restart_nout   !< output interval of restart file
    !< (if  .gt.0) restart file write
    !< (if  .lt.0) restart file read and write
    integer(kind=kint) :: restart_nin    !< input number of restart
    type(step_info)    :: step_ctrl_restart  !< step information for restart

    integer(kind=kint) :: max_lyr          !< maximum num of layer
    integer(kind=kint) :: is_33shell
    integer(kind=kint) :: is_33beam
    integer(kind=kint) :: is_heat
    integer(kind=kint) :: max_ncon_stf     !< maximum num of stiffness matrix size
    integer(kind=kint) :: max_ncon         !< maximum num of element connectivity
    integer(kind=kint), pointer :: is_rot(:) => null()
    integer(kind=kint) :: elemopt361
    logical            :: is_smoothing_active
    real(kind=kreal)   :: FACTOR     (2)   !< factor of incrementation
    !< 1:time t  2: time t+dt
    !> for increment control
    integer(kind=kint) :: NRstat_i(10)     !< statistics of newton iteration (integer)
    real(kind=kreal)   :: NRstat_r(10)     !< statistics of newton iteration (real)
    integer(kind=kint) :: AutoINC_stat     !< status of auto-increment control
    integer(kind=kint) :: CutBack_stat     !< status of cutback control

    real(kind=kreal), pointer :: GL          (:)           !< external force
    real(kind=kreal), pointer :: EFORCE      (:)           !< external force
    real(kind=kreal), pointer :: QFORCE      (:)           !< equivalent nodal force
    real(kind=kreal), pointer :: unode(:)      => null()   !< disp at the beginning of curr step
    real(kind=kreal), pointer :: unode_bak(:)  => null()   !< disp at the beginning of curr step
    real(kind=kreal), pointer :: dunode(:)     => null()   !< curr total disp
    real(kind=kreal), pointer :: ddunode(:)    => null()   !< =hecMESH%X, disp increment
    real(kind=kreal), pointer :: temperature(:)=> null()   !< =temperature
    real(kind=kreal), pointer :: temp_bak(:)   => null()
    real(kind=kreal), pointer :: last_temp(:)  => null()

    type( tElement ), pointer :: elements(:)   =>null()  !< elements information
    type( tMaterial ),pointer :: materials(:)  =>null()  !< material properties
    integer                   :: n_contacts              !< number of contact conditions
    type( tContact ), pointer :: contacts(:)   =>null()  !< contact information
    integer                   :: n_embeds               !< number of embed conditions
    type( tContact ), pointer :: embeds(:)   =>null()  !< contact information
    integer                   :: n_fix_mpc               !< number mpc conditions user defined
    real(kind=kreal), pointer :: mpc_const(:)  =>null()  !< bakeup of hecmwST_mpc%mpc_const
    type(tSection), pointer   :: sections(:)   =>null()  !< definition of section referred by elements(i)%sectionID

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
    type( tContact ), pointer :: embeds_bkup(:)  =>null()  !< contact information (backup)
  end type fstr_solid

  !> Data for HEAT ANSLYSIS  (fstrHEAT)
  type fstr_heat
    !> Crank-Nicolson parameter
    integer(kind=kint) :: is_steady
    real(kind=kreal)   :: beta
    logical :: is_iter_max_limit

    !> TIME CONTROL
    integer(kind=kint) :: STEPtot
    integer(kind=kint) :: restart_nout
    real(kind=kreal), pointer :: STEP_DLTIME(:), STEP_EETIME(:)
    real(kind=kreal), pointer :: STEP_DELMIN(:), STEP_DELMAX(:)
    integer(kind=kint) :: timepoint_id

    !> MATERIAL
    integer(kind=kint) :: MATERIALtot
    integer(kind=kint), pointer :: RHOtab(:), CPtab(:), CONDtab(:)
    real(kind=kreal), pointer :: RHO(:,:), RHOtemp (:,:)
    real(kind=kreal), pointer :: CP(:,:),  CPtemp (:,:)
    real(kind=kreal), pointer :: COND(:,:),CONDtemp (:,:)

    real(kind=kreal), pointer :: RHOfuncA(:,:),  RHOfuncB(:,:)
    real(kind=kreal), pointer :: CPfuncA (:,:),  CPfuncB(:,:)
    real(kind=kreal), pointer :: CONDfuncA (:,:),CONDfuncB(:,:)

    !> AMPLITUDE
    integer(kind=kint) :: AMPLITUDEtot
    integer(kind=kint), pointer :: AMPLtab(:)
    real(kind=kreal), pointer :: AMPL(:,:), AMPLtime (:,:)
    real(kind=kreal), pointer :: AMPLfuncA(:,:),  AMPLfuncB(:,:)

    !> VALUE
    real(kind=kreal), pointer :: TEMP0(:)
    real(kind=kreal), pointer :: TEMPC(:)
    real(kind=kreal), pointer :: TEMP (:)

    !> FIXTEMP
    integer(kind=kint) :: T_FIX_tot
    integer(kind=kint), pointer :: T_FIX_node(:)
    integer(kind=kint), pointer :: T_FIX_ampl(:)
    real(kind=kreal), pointer :: T_FIX_val(:)

    !> CFLUX
    integer(kind=kint) :: Q_NOD_tot
    integer(kind=kint), pointer :: Q_NOD_node(:)
    integer(kind=kint), pointer :: Q_NOD_ampl(:)
    real(kind=kreal), pointer :: Q_NOD_val(:)

    !> DFLUX (not used)
    integer(kind=kint) :: Q_VOL_tot
    integer(kind=kint), pointer :: Q_VOL_elem(:)
    integer(kind=kint), pointer :: Q_VOL_ampl(:)
    real(kind=kreal), pointer :: Q_VOL_val(:)

    !> DFLUX, !SFLUX
    integer(kind=kint) :: Q_SUF_tot
    integer(kind=kint), pointer :: Q_SUF_elem(:)
    integer(kind=kint), pointer :: Q_SUF_ampl(:)
    integer(kind=kint), pointer :: Q_SUF_surf(:)
    real(kind=kreal), pointer :: Q_SUF_val(:)

    !> RADIATE, !SRADIATE
    integer(kind=kint) :: R_SUF_tot
    integer(kind=kint), pointer :: R_SUF_elem(:)
    integer(kind=kint), pointer :: R_SUF_ampl(:,:)
    integer(kind=kint), pointer :: R_SUF_surf(:)
    real(kind=kreal), pointer :: R_SUF_val(:,:)

    !> FILM, SFILM
    integer(kind=kint) :: H_SUF_tot
    integer(kind=kint), pointer :: H_SUF_elem(:)
    integer(kind=kint), pointer :: H_SUF_ampl(:,:)
    integer(kind=kint), pointer :: H_SUF_surf(:)
    real(kind=kreal), pointer :: H_SUF_val(:,:)

    integer(kind=kint)  :: WL_tot
    type(tWeldLine), pointer :: weldline(:) => null()
  end type fstr_heat

  !> Data for DYNAMIC ANSLYSIS  (fstrDYNAMIC)
  type fstr_dynamic
    !> ANALYSIS TYPE CONTROL
    integer(kind=kint) :: idx_eqa       ! implicit or explicit
    integer(kind=kint) :: idx_resp      ! time history or steady-state harmonic response analysis

    !> TIME CONTROL
    integer(kind=kint) :: n_step        ! total step number of analysis
    real(kind=kreal)   :: t_start       ! start time of analysis
    real(kind=kreal)   :: t_curr        ! current time of analysis
    real(kind=kreal)   :: t_end         ! end time of analysis
    real(kind=kreal)   :: t_delta       ! time increment
    integer(kind=kint) :: restart_nout  ! output interval of restart file
    ! (if  .gt.0) restart file write
    ! (if  .lt.0) restart file read and write
    integer(kind=kint) :: restart_nin  !input number of restart file

    !> Newmark-beta parameter
    real(kind=kreal)   :: ganma         ! Newmark-beta parameter ganma
    real(kind=kreal)   :: beta          ! Newmark-beta parameter beta

    !> mass matrix control
    integer(kind=kint) :: idx_mas       ! mass matrix type

    !> damping control
    integer(kind=kint) :: idx_dmp      ! damping type
    real(kind=kreal)   :: ray_m         ! Rayleigh damping parameter Rm
    real(kind=kreal)   :: ray_k         ! Rayleigh damping parameter Rk

    !> initialization control
    logical           :: VarInitialize  ! initialization flag

    !> OUTPUT CONTROL
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

    !> VALUE
    real(kind=kreal), pointer :: DISP  (:,:)     !> Displacement, U(t+dt), U(t), U(t-dt)
    real(kind=kreal), pointer :: VEL   (:,:)     !> Velocity
    real(kind=kreal), pointer :: ACC   (:,:)     !> Acceleration

    real(kind=kreal) :: kineticEnergy
    real(kind=kreal) :: strainEnergy
    real(kind=kreal) :: totalEnergy

    !> temporary quantity
    real(kind=kreal), pointer :: VEC1  (:)
    real(kind=kreal), pointer :: VEC2  (:)
    real(kind=kreal), pointer :: VEC3  (:)

    integer(kind=kint) :: dynamic_IW4        =   204
    integer(kind=kint) :: dynamic_IW5        =   205
    integer(kind=kint) :: dynamic_IW6        =   206
    integer(kind=kint) :: dynamic_IW7        =   207
    integer(kind=kint) :: dynamic_IW8        =   208
    integer(kind=kint) :: dynamic_IW9        =   209
    integer(kind=kint) :: dynamic_IW10       =   210
  end type fstr_dynamic

  type fstr_freqanalysis
    integer(kind=kint)                :: FLOAD_ngrp_tot
    integer(kind=kint), pointer       :: FLOAD_ngrp_GRPID(:) => null()
    integer(kind=kint), pointer       :: FLOAD_ngrp_ID(:)    => null()
    integer(kind=kint), pointer       :: FLOAD_ngrp_TYPE(:)  => null()
    integer(kind=kint), pointer       :: FLOAD_ngrp_DOF(:)   => null()
    real(kind=kreal), pointer         :: FLOAD_ngrp_valre(:) => null()
    real(kind=kreal), pointer         :: FLOAD_ngrp_valim(:) => null()
    character(len=HECMW_FILENAME_LEN) :: eigenlog_filename
    integer(kind=kint)                :: start_mode
    integer(kind=kint)                :: end_mode
  end type fstr_freqanalysis

  type fstr_freqanalysis_data
    integer(kind=kint)        :: numMode
    integer(kind=kint)        :: numNodeDOF
    real(kind=kreal), pointer :: eigOmega(:)    => null()
    real(kind=kreal), pointer :: eigVector(:,:) => null()
    real(kind=kreal)           :: rayAlpha, rayBeta
  end type fstr_freqanalysis_data

  !> Package of data used by Lanczos eigenvalue solver
  type fstr_eigen
    !> Allocatable array, used or Lanczos eigenvalue analysis
    integer(kind=kint)  :: nget      ! Solved eigen value number (default:5)
    integer(kind=kint)  :: maxiter   ! Max. Lcz iterations (default:60)
    integer(kind=kint)  :: iter      ! Max. Lcz iterations (default:60)
    real   (kind=kreal) :: sigma     ! 0.0
    real   (kind=kreal) :: tolerance ! Lcz tolerance (default:1.0e-8)
    real   (kind=kreal) :: totalmass
    real   (kind=kreal), pointer :: eigval(:)
    real   (kind=kreal), pointer :: eigvec(:,:)
    real   (kind=kreal), pointer :: filter(:)
    real   (kind=kreal), pointer :: mass(:)
    real   (kind=kreal), pointer :: effmass(:)
    real   (kind=kreal), pointer :: partfactor(:)
    logical :: is_free = .false.
  end type fstr_eigen

  !> Data for coupling analysis
  type fstr_couple
    !> for revocap
    integer( kind=kint )   :: dof              ! == 3
    integer( kind=kint )   :: ndof             ! total dof (coupled_node_n*dof)
    integer( kind=kint )   :: coupled_node_n
    !> note) following types depend on revocap
    integer, pointer       :: coupled_node(:)  ! local node id sent from revocap
    real( kind=8 ),pointer :: trac(:)          ! input  (x,y,z,x,y,z ... )
    real( kind=8 ),pointer :: disp(:)          ! output (x,y,z,x,y,z ... )
    real( kind=8 ),pointer :: velo(:)          ! output (x,y,z,x,y,z ... )
    real( kind=8 ),pointer :: accel(:)         ! output (x,y,z,x,y,z ... )
    !> for inner use
    integer( kind=kint ),pointer :: index(:)   ! size:total node num.
    !>  -1:not relation, >1:index of coupled_node
  end type fstr_couple

  !> Data for weld line
  type tWeldLine
    integer(kind=kint) :: egrpid
    real( kind=kreal ) :: I
    real( kind=kreal ) :: U
    real( kind=kreal ) :: coe
    real( kind=kreal ) :: v
    integer(kind=kint) :: xyz
    real(kind=kreal)   :: n1, n2
    real(kind=kreal)   :: distol
    real(kind=kreal)   :: tstart
  end type tWeldLine

  !> Data for section control
  type tSection
    !integer              :: mat_ID
    !integer              :: iset
    !integer              :: orien_ID
    !real(kind=kreal)     :: thickness
    integer              :: elemopt341
    !integer              :: elemopt342
    !integer              :: elemopt351
    !integer              :: elemopt352
    integer              :: elemopt361
    !integer              :: elemopt362
  end type tSection

contains

  !> NULL POINTER SETTING TO AVOID RUNTIME ERROR
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

  subroutine fstr_nullify_fstr_solid( S )
    implicit none
    type( fstr_solid ) :: S

    nullify( S%BOUNDARY_ngrp_ID )
    nullify( S%BOUNDARY_ngrp_type )
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
    nullify( S%PSTRESS )
    nullify( S%PSTRAIN )
    nullify( S%PSTRESS_VECT )
    nullify( S%PSTRAIN_VECT )
    nullify( S%REACTION )
    nullify( S%ESTRESS )
    nullify( S%ESTRAIN )
    nullify( S%EMISES )
    nullify( S%EPSTRESS )
    nullify( S%EPSTRAIN )
    nullify( S%EPSTRESS_VECT )
    nullify( S%EPSTRAIN_VECT )
    nullify( S%ENQM )
    nullify( S%GL          )
    nullify( S%QFORCE      )
    nullify( S%VELOCITY_ngrp_ID )
    nullify( S%VELOCITY_ngrp_type )
    nullify( S%VELOCITY_ngrp_amp )
    nullify( S%VELOCITY_ngrp_val )
    nullify( S%ACCELERATION_ngrp_ID )
    nullify( S%ACCELERATION_ngrp_type )
    nullify( S%ACCELERATION_ngrp_amp )
    nullify( S%ACCELERATION_ngrp_val )
    nullify( S%COUPLE_ngrp_ID )
  end subroutine fstr_nullify_fstr_solid

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

  subroutine fstr_nullify_fstr_dynamic( DY )
    implicit none
    type( fstr_dynamic ) :: DY

    nullify( DY%DISP )
    nullify( DY%VEL  )
    nullify( DY%ACC  )
    nullify( DY%VEC1 )
    nullify( DY%VEC2 )
    nullify( DY%VEC3 )
  end subroutine fstr_nullify_fstr_dynamic

  subroutine fstr_nullify_fstr_freqanalysis( f )
    implicit none
    type( fstr_freqanalysis ), intent(inout) :: f

    f%FLOAD_ngrp_tot = 0
    nullify( f%FLOAD_ngrp_GRPID )
    nullify( f%FLOAD_ngrp_ID )
    nullify( f%FLOAD_ngrp_TYPE )
    nullify( f%FLOAD_ngrp_DOF )
    nullify( f%FLOAD_ngrp_valre )
    nullify( f%FLOAD_ngrp_valim )
  end subroutine fstr_nullify_fstr_freqanalysis

  subroutine fstr_nullify_fstr_eigen( E )
    implicit none
    type( fstr_eigen ) :: E

    nullify( E%mass )
  end subroutine fstr_nullify_fstr_eigen

  subroutine fstr_nullify_fstr_couple( C )
    implicit none
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
    hecMAT%Iarray(13)=    0    ! = mpc_method
    hecMAT%Iarray(14)=    0    ! = estcond
    hecMAT%Iarray(35)=    3    ! = maxrecycle_precond
    hecMAT%Iarray(41)=    0    ! = solver_opt1
    hecMAT%Iarray(42)=    0    ! = solver_opt2
    hecMAT%Iarray(43)=    0    ! = solver_opt3
    hecMAT%Iarray(44)=    0    ! = solver_opt4
    hecMAT%Iarray(45)=    0    ! = solver_opt5
    hecMAT%Iarray(46)=    0    ! = solver_opt6

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

  subroutine hecMAT_init( hecMAT )
    implicit none
    type( hecmwST_matrix ) :: hecMAT
    integer ::  ndof, nn, ierror
    ndof = hecMAT%NDOF
    nn = ndof*ndof
    allocate (hecMAT%AL(nn*hecMAT%NPL)        ,stat=ierror )
    if( ierror /= 0 ) then
      write(*,*) "##ERROR : not enough memory"
      write(idbg,*) 'stop due to allocation error'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm() )
    end if
    allocate (hecMAT%AU(nn*hecMAT%NPU)        ,stat=ierror )
    if( ierror /= 0 ) then
      write(*,*) "##ERROR : not enough memory"
      write(idbg,*) 'stop due to allocation error'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm() )
    end if
    allocate (hecMAT%B(ndof*hecMAT%NP)          ,stat=ierror )
    if( ierror /= 0 ) then
      write(*,*) "##ERROR : not enough memory"
      write(idbg,*) 'stop due to allocation error'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm() )
    end if
    hecMAT%B(:)=0.d0
    allocate (hecMAT%D(nn*hecMAT%NP)          ,stat=ierror )
    if( ierror /= 0 ) then
      write(*,*) "##ERROR : not enough memory"
      write(idbg,*) 'stop due to allocation error'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm() )
    end if
    allocate (hecMAT%X(ndof*hecMAT%NP)          ,stat=ierror )
    if( ierror /= 0 ) then
      write(*,*) "##ERROR : not enough memory"
      write(idbg,*) 'stop due to allocation error'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm() )
    end if
    allocate (hecMAT%ALU(nn*hecMAT%N)         ,stat=ierror )
    if( ierror /= 0 ) then
      write(*,*) "##ERROR : not enough memory"
      write(idbg,*) 'stop due to allocation error'
      call flush(idbg)
      call hecmw_abort( hecmw_comm_get_comm() )
    endif
    call hecmw_cmat_init( hecMAT%cmat )
    hecMAT%D  = 0.0d0
    hecMAT%AL = 0.0d0
    hecMAT%AU = 0.0d0
    hecMAT%B  = 0.0d0
    hecMAT%X  = 0.0d0
    hecMAT%ALU = 0.0d0
  end subroutine hecMAT_init

  subroutine hecMAT_finalize( hecMAT )
    implicit none
    type( hecmwST_matrix ) :: hecMAT
    integer ::  ndof, nn, ierror
    ndof = hecMAT%NDOF
    nn = ndof*ndof
    if( associated(hecMAT%AL) ) then
      deallocate(hecMAT%AL                  ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(hecMAT%AU) ) then
      deallocate(hecMAT%AU                  ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(hecMAT%B) ) then
      deallocate(hecMAT%B                   ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(hecMAT%D) ) then
      deallocate(hecMAT%D                   ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(HECMAT%X) ) then
      deallocate(hecMAT%X                   ,stat=ierror)
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to deallocation error'
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
    endif
    if( associated(hecMAT%ALU) ) then
      deallocate(hecMAT%ALU                 ,stat=ierror)
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

  logical function fstr_isLoadActive( fstrSOLID, nbc, cstep )
    type(fstr_solid)    :: fstrSOLID
    integer, intent(in) :: nbc    !< group id of boundary condition
    integer, intent(in) :: cstep  !< current step number
    fstr_isLoadActive = .true.
    if( cstep > 0 ) then
      if( .not. associated(fstrSOLID%step_ctrl) ) return
      if( cstep>fstrSOLID%nstep_tot ) return
      fstr_isLoadActive = isLoadActive( nbc, fstrSOLID%step_ctrl(cstep) )
    else
      fstr_isLoadActive = isLoadActive( nbc, fstrSOLID%step_ctrl_restart )
    endif
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

  logical function fstr_isEmbedActive( fstrSOLID, nbc, cstep )
    type(fstr_solid)     :: fstrSOLID !< type fstr_solid
    integer, intent(in) :: nbc       !< group id of boundary condition
    integer, intent(in) :: cstep     !< current step number
    fstr_isEmbedActive = .true.
    if( .not. associated(fstrSOLID%step_ctrl) ) return
    if( cstep>fstrSOLID%nstep_tot ) return
    fstr_isEmbedActive = isContactActive( nbc, fstrSOLID%step_ctrl(cstep) )
  end function

  !> This subroutine fetch coords defined by local coordinate system
  subroutine get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )
    integer, intent(in)             :: cdsys_ID      !< id of local coordinate
    type(hecmwST_local_mesh)       :: hecMESH       !< mesh information
    type(fstr_solid), intent(inout) :: fstrSOLID     !< fstr_solid
    real(kind=kreal), intent(out)   :: coords(3,3)
    integer :: ik

    coords = 0.d0
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
  end subroutine get_coordsys

  subroutine fstr_set_current_config_to_mesh(hecMESH,fstrSOLID,coord)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type (fstr_solid), intent(in) :: fstrSOLID
    real(kind=kreal), intent(out) :: coord(:)
    integer(kind=kint) :: i
    if(hecMESH%n_dof == 4) return
    do i = 1, hecMESH%nn_internal*min(hecMESH%n_dof,3)
      coord(i) = hecMESH%node(i)
      hecMESH%node(i) = coord(i)+fstrSOLID%unode(i)+fstrSOLID%dunode(i)
    enddo
  end subroutine fstr_set_current_config_to_mesh

  subroutine fstr_recover_initial_config_to_mesh(hecMESH,fstrSOLID,coord)
    implicit none
    type(hecmwST_local_mesh), intent(inout) :: hecMESH
    type (fstr_solid), intent(in) :: fstrSOLID
    real(kind=kreal), intent(in) :: coord(:)
    integer(kind=kint) :: i
    if(hecMESH%n_dof == 4) return
    do i = 1, hecMESH%nn_internal*min(hecMESH%n_dof,3)
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
    phys%ENQM    = 0.0d0
  end subroutine fstr_solid_phys_zero

  subroutine fstr_solid_phys_clear(fstrSOLID)
    implicit none
    type (fstr_solid)         :: fstrSOLID
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
