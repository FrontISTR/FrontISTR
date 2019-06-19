!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module m_fstr_main

  use hecmw
  use m_fstr
  use m_hecmw2fstr_mesh_conv
  use m_fstr_setup
  use m_fstr_solve_heat
  use m_fstr_solve_nlgeom
  use m_fstr_solve_eigen
  use m_static_echo
  use m_heat_init
  use m_heat_echo
  use m_fstr_precheck
  use m_fstr_rcap_io
  use fstr_solver_dynamic
  use fstr_debug_dump
  use fstr_matrix_con_contact

  type(hecmwST_local_mesh), save             :: hecMESH
  type(hecmwST_matrix), save                 :: hecMAT
  type(hecmwST_matrix), save                 :: conMAT
  type(fstr_solid), save                     :: fstrSOLID
  type(fstrST_matrix_contact_lagrange), save :: fstrMAT
  type(fstr_heat), save                      :: fstrHEAT
  type(fstr_eigen), save                     :: fstrEIG
  type(fstr_dynamic), save                   :: fstrDYNAMIC
  type(hecmwST_result_data), save            :: fstrRESULT
  type(fstr_couple), save                    :: fstrCPL
  type(fstr_freqanalysis), save              :: fstrFREQ
  character(len=HECMW_FILENAME_LEN)          :: name_ID

contains

  subroutine fstr_main() bind(C,NAME='fstr_main')
    implicit none
    real(kind=kreal) :: T1, T2, T3

    T1=0.0d0; T2=0.0d0; T3=0.0d0

    ! =============== INITIALIZE ===================

    call hecmw_init
    myrank = hecmw_comm_get_rank()
    nprocs = hecmw_comm_get_size()

    T1 = hecmw_Wtime()

    name_ID = 'fstrMSH'
    call hecmw_get_mesh( name_ID , hecMESH )

    if( hecMESH%contact_pair%n_pair > 0 ) then
      if( nprocs > 1 .and. &
          hecMESH%hecmw_flag_partcontact /= HECMW_FLAG_PARTCONTACT_AGGREGATE ) then
        paraContactFlag = .true.
      endif
      if( myrank == 0 ) then
        print *,'paraContactFlag',paraContactFlag
      endif
    endif

    call hecmw2fstr_mesh_conv( hecMESH )

    call fstr_init

    call fstr_rcap_initialize( hecMESH, fstrPR, fstrCPL )

    T2 = hecmw_Wtime()

    ! =============== ANALYSIS =====================

    select case( fstrPR%solution_type )
      case( kstSTATIC )
        call fstr_static_analysis
      case( kstDYNAMIC )
        call fstr_dynamic_analysis
      case( kstEIGEN )
        call fstr_eigen_analysis
      case( kstHEAT )
        call fstr_heat_analysis
      case( kstSTATICEIGEN )
        call fstr_static_eigen_analysis
      case( kstPRECHECK, kstNZPROF )
        call fstr_precheck( hecMESH, hecMAT, fstrPR%solution_type )
    end select

    T3 = hecmw_Wtime()

    if(hecMESH%my_rank==0) then
      write(*,*)
      write(*,*)           '===================================='
      write(*,'(a,f10.2)') '    TOTAL TIME (sec) :', T3 - T1
      write(*,'(a,f10.2)') '           pre (sec) :', T2 - T1
      write(*,'(a,f10.2)') '         solve (sec) :', T3 - T2
      write(*,*)           '===================================='

      write(IMSG,*)           '===================================='
      write(IMSG,'(a,f10.2)') '    TOTAL TIME (sec) :', T3 - T1
      write(IMSG,'(a,f10.2)') '           pre (sec) :', T2 - T1
      write(IMSG,'(a,f10.2)') '         solve (sec) :', T3 - T2
      write(IMSG,*)           '===================================='
    endif

    ! =============== FINALIZE =====================

    call fstr_rcap_finalize( fstrPR, fstrCPL )
    call fstr_finalize()
    call hecmw_dist_free(hecMESH)
    call hecmw_finalize
    if(hecMESH%my_rank==0) write(*,*) 'FrontISTR Completed !!'

  end subroutine fstr_main

  !=============================================================================!
  !> Initializer                                                                !
  !=============================================================================!

  subroutine fstr_init
    implicit none

    ! set pointer to null
    call hecmw_nullify_matrix     ( hecMAT      )
    call hecmw_nullify_matrix     ( conMAT      )
    call hecmw_nullify_result_data( fstrRESULT  )
    call fstr_nullify_fstr_param  ( fstrPR      )
    call fstr_nullify_fstr_solid  ( fstrSOLID   )
    call fstr_nullify_fstr_heat   ( fstrHEAT    )
    call fstr_nullify_fstr_eigen  ( fstrEIG     )
    call fstr_nullify_fstr_dynamic( fstrDYNAMIC )
    call fstr_nullify_fstr_couple ( fstrCPL     )
    call fstr_init_file

    ! ----  default setting of global params ---
    DT    = 1
    ETIME = 1
    ITMAX = 20
    EPS   = 1.0d-6

    ! -------  grobal pointer setting ----------
    REF_TEMP => fstrPR%ref_temp
    IECHO    => fstrPR%fg_echo
    IRESULT  => fstrPR%fg_result
    IVISUAL  => fstrPR%fg_visual

    ! for heat ...
    INEUTRAL => fstrPR%fg_neutral
    IRRES    => fstrPR%fg_irres
    IWRES    => fstrPR%fg_iwres
    NRRES    => fstrPR%nrres
    NPRINT   => fstrPR%nprint

    call hecmw_mat_con(hecMESH, hecMAT)

    ! ------- initial value setting -------------
    call fstr_mat_init  ( hecMAT   )
    call fstr_param_init( fstrPR, hecMESH )

    call fstr_solid_init( hecMESH, fstrSOLID )
    call fstr_eigen_init( fstrEIG )
    call fstr_heat_init ( fstrHEAT  )
    call fstr_dynamic_init( fstrDYNAMIC  )

    call fstr_init_condition
    hecMAT%NDOF = hecMESH%n_dof
    if( kstHEAT == fstrPR%solution_type ) then
      call heat_init_material (hecMESH,fstrHEAT)
      call heat_init_amplitude(hecMESH,fstrHEAT)
      hecMAT%NDOF = 1
    endif
    call hecMAT_init( hecMAT )

  end subroutine fstr_init

  !------------------------------------------------------------------------------
  !> Open all files preparing calculation
  subroutine fstr_init_file
    implicit none
    character(len=HECMW_FILENAME_LEN) :: s, r
    character(len=HECMW_FILENAME_LEN) :: stafileNAME
    character(len=HECMW_FILENAME_LEN) :: logfileNAME
    character(len=HECMW_FILENAME_LEN) :: msgfileNAME
    character(len=HECMW_FILENAME_LEN) :: dbgfileNAME
    integer :: stat, flag, limit, irank

    ! set file name --------------------------------
    call hecmw_ctrl_is_subdir( flag, limit )
    write(s,*) myrank
    if( flag == 0 ) then
      write( logfileNAME, '(a,a)') trim(adjustl(s)), '.log'
      logfileNAME = adjustl(logfileNAME)
      write( dbgfileNAME, '(a,a)') 'FSTR.dbg.', trim(adjustl(s))
      dbgfileNAME = adjustl(dbgfileNAME)
    else
      if( nprocs > limit ) then
        irank = myrank / limit
        write(r,*) irank
        write( logfileNAME, '(a,a,a,a,a)') 'LOG/TRUNK', trim(adjustl(r)), '/', trim(adjustl(s)), '.log'
        logfileNAME = adjustl(logfileNAME)
        call hecmw_ctrl_make_subdir( logfileNAME, stat )
        if( stat /= 0 ) call fstr_setup_util_err_stop( '### Cannot create directory' )
        write( dbgfileNAME, '(a,a,a,a,a)') 'DBG/TRUNK', trim(adjustl(r)), '/', 'FSTR.dbg.', trim(adjustl(s))
        dbgfileNAME = adjustl(dbgfileNAME)
        call hecmw_ctrl_make_subdir( dbgfileNAME, stat )
        if( stat /= 0 ) call fstr_setup_util_err_stop( '### Cannot create directory' )
      else
        write( logfileNAME, '(a,a,a)') 'LOG/', trim(adjustl(s)), '.log'
        logfileNAME = adjustl(logfileNAME)
        call hecmw_ctrl_make_subdir( logfileNAME, stat )
        if( stat /= 0 ) call fstr_setup_util_err_stop( '### Cannot create directory' )
        write( dbgfileNAME, '(a,a,a)') 'DBG/', 'FSTR.dbg.', trim(adjustl(s))
        dbgfileNAME = adjustl(dbgfileNAME)
        call hecmw_ctrl_make_subdir( dbgfileNAME, stat )
        if( stat /= 0 ) call fstr_setup_util_err_stop( '### Cannot create directory' )
      endif
    endif
    stafileNAME = 'FSTR.sta'
    msgfileNAME = 'FSTR.msg'

    ! open & opening message out -------------------
    ! MSGFILE
    if( myrank == 0) then
      open(IMSG, file=msgfileNAME, status='replace', iostat=stat)
      if( stat /= 0 ) then
        call fstr_setup_util_err_stop( '### Cannot open message file :'//msgfileNAME )
      endif
      write(IMSG,*) ':========================================:'
      write(IMSG,*) ':**   BEGIN FSTR Structural Analysis   **:'
      write(IMSG,*) ':========================================:'
      write(IMSG,*) '        Total no. of processors: ',nprocs
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*) ' *    STAGE Initialization and input   **'
    endif

    ! LOGFILE & STAFILE
    open (ILOG, file = logfileNAME, status = 'replace', iostat=stat )
    if( stat /= 0 ) then
      call fstr_setup_util_err_stop( '### Cannot open log file :'//logfileNAME )
    endif

    if( myrank == 0 ) then
      open (ISTA,file = stafileNAME, status = 'replace', iostat=stat )
      write(ISTA,'(''####''a80)') stafileNAME
      if( stat /= 0 ) then
        call fstr_setup_util_err_stop( '### Cannot open status file :'//stafileNAME )
      endif
    endif

    open (IDBG,file = dbgfileNAME, status = 'replace')
    write(IDBG,'(''####''a80)') dbgfileNAME
    if( stat /= 0 ) then
      call fstr_setup_util_err_stop( '### Cannot open debug file :'//dbgfileNAME )
    endif
  end subroutine fstr_init_file

  !------------------------------------------------------------------------------
  !> Read in control file and do all preparation
  subroutine fstr_init_condition
    implicit none
    character(len=HECMW_FILENAME_LEN) :: cntfileNAME

    name_ID='fstrCNT'
    call hecmw_ctrl_get_control_file( name_ID, cntfileNAME )

    ! loading boundary conditions etc. from fstr control file or nastran mesh file
    ! and setup parameters ...
    svRarray(:) = hecMAT%Rarray(:)
    svIarray(:) = hecMAT%Iarray(:)

    call fstr_setup( cntfileNAME, hecMESH, fstrPR, fstrSOLID, fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL, fstrFREQ )

    hecMAT%Rarray(:) = svRarray(:)
    hecMAT%Iarray(:) = svIarray(:)

    if( myrank == 0) write(*,*) 'fstr_setup: OK'
    write(ILOG,*) 'fstr_setup: OK'
    call flush(6)

  end subroutine fstr_init_condition

  !=============================================================================!
  !> Master subroutine of linear/nonlinear static analysis                      !
  !=============================================================================!

  subroutine fstr_static_analysis
    implicit none

    if( IECHO.eq.1 ) call fstr_echo(hecMESH)

    if(myrank .EQ. 0) then
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
    endif

    if( fstrPR%nlgeom ) then
      if( myrank == 0)  write(IMSG,*) ' ***   STAGE Non Linear static analysis   **'
    else
      if( myrank == 0 ) write(IMSG,*) ' ***   STAGE Linear static analysis   **'
    endif

    if( paraContactFlag ) then
      call fstr_solve_NLGEOM( hecMESH, hecMAT, fstrSOLID, fstrMAT, fstrPR, conMAT )
    else
      call fstr_solve_NLGEOM( hecMESH, hecMAT, fstrSOLID, fstrMAT, fstrPR )
    endif

    call fstr_solid_finalize( fstrSOLID )

  end subroutine fstr_static_analysis

  !=============================================================================!
  !> Master subroutine of eigen analysis                                        !
  !=============================================================================!

  subroutine fstr_eigen_analysis
    use hecmw
    use m_fstr
    implicit none

    if( IECHO.eq.1 ) call fstr_echo(hecMESH)
    if(myrank .EQ. 0) then
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*) ' ***   STAGE Eigenvalue analysis     **'
    endif

    call fstr_solve_EIGEN( hecMESH, hecMAT, fstrEIG, fstrSOLID, fstrRESULT, fstrPR, fstrMAT )

  end subroutine fstr_eigen_analysis

  !=============================================================================!
  !> Master subroutine of heat analysis                                         !
  !=============================================================================!

  subroutine fstr_heat_analysis
    implicit none

    if( IECHO.eq.1 ) call heat_echo(fstrPR,hecMESH,fstrHEAT)
    if(myrank .EQ. 0) then
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*) ' ***   STAGE Heat analysis    **'
    endif

    call fstr_solve_HEAT( hecMESH, hecMAT, fstrRESULT, fstrPR, fstrHEAT )

  end subroutine fstr_heat_analysis

  !=============================================================================!
  !> Master subroutine of dynamic analysis                                      !
  !=============================================================================!

  subroutine fstr_dynamic_analysis
    implicit none

    if( IECHO.eq.1 ) call fstr_echo(hecMESH)

    if(myrank == 0) then
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
      if( fstrPR%nlgeom ) then
        write(IMSG,*) ' ***   STAGE Nonlinear dynamic analysis   **'
      else
        write(IMSG,*) ' ***   STAGE Linear dynamic analysis   **'
      endif
    endif

    if( paraContactFlag ) then
      call fstr_solve_dynamic( hecMESH, hecMAT, fstrSOLID, fstrEIG, &
        fstrDYNAMIC, fstrRESULT, fstrPR, fstrCPL, fstrFREQ, fstrMAT, &
        conMAT )
    else
      call fstr_solve_dynamic( hecMESH, hecMAT, fstrSOLID, fstrEIG, &
        fstrDYNAMIC, fstrRESULT, fstrPR, fstrCPL, fstrFREQ, fstrMAT)
    endif

  end subroutine fstr_dynamic_analysis

  !=============================================================================!
  !> Master subroutine of static -> eigen anaylsis                              !
  !=============================================================================!

  subroutine fstr_static_eigen_analysis
    implicit none

    if( IECHO==1 ) call fstr_echo(hecMESH)

    if(myrank == 0) then
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*) ' ***   STAGE Static -> Eigen analysis   **'
      write(*,*) ' ***   STAGE Static -> Eigen analysis   **'
      write(IMSG,*)
      write(IMSG,*) ' ***   Stage 1: Nonlinear dynamic analysis  **'
      write(*,*) ' ***   Stage 1: Nonlinear dynamic analysis   **'
    endif

    call fstr_solve_NLGEOM( hecMESH, hecMAT, fstrSOLID, fstrMAT, fstrPR )

    if(myrank == 0) then
      write(IMSG,*)
      write(IMSG,*) ' ***   Stage 2: Eigenvalue analysis  **'
      write(*,*)
      write(*,*) ' ***   Stage 2: Eigenvalue analysis   **'
    endif

    call fstr_solve_EIGEN( hecMESH, hecMAT, fstrEIG, fstrSOLID, fstrRESULT, fstrPR, fstrMAT )

    call fstr_solid_finalize( fstrSOLID )

  end subroutine fstr_static_eigen_analysis

  !=============================================================================!
  !> Finalizer                                                                  !
  !=============================================================================!

  subroutine fstr_finalize
    implicit none

    if( myrank == 0 ) then
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*)
      write(IMSG,*) ':========================================:'
      write(IMSG,*) ':**            END of FSTR             **:'
      write(IMSG,*) ':========================================:'
      close(IMSG)
      close(ISTA)
    endif

    call fstr_solid_finalize( fstrSOLID )
    call hecMAT_finalize( hecMAT )

    close(ILOG)
    close(IDBG)
  end subroutine fstr_finalize

end module m_fstr_main
