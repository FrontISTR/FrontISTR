!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Main subroutine                                   !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

program fstr_main
     use hecmw
     use m_fstr
     use m_fstr_setup
     use m_fstr_precheck
     use m_static_echo
     use m_fstr_precheck
     use m_fstr_rcap_io
     use fstr_solver_dynamic
     use fstr_debug_dump
     implicit none

     type (hecmwST_local_mesh) :: hecMESH
     type (hecmwST_matrix )    :: hecMAT
	 type (fstr_solid )        :: fstrSOLID
     type (fstr_heat )         :: fstrHEAT
     type (lczparam)           :: fstrEIG
     type (fstr_dynamic )      :: fstrDYNAMIC
     type (fstr_couple )       :: fstrCPL
	 
	 include "HEC_MW3_For.h"
	 
	 real :: T1,T2,T3
	 integer :: ierr, iAss, iPart, nAss, nPart, nRefine, refinetype
 
	 integer :: i, narg, argv_len, path_len, cntn, cnte
	 character(len=100) :: fname, argv(10)
	 character(len=100) :: cntpath = "./hecmw_ctrl.dat"// CHAR(0)
	
  !   argv(1) = "./fstr"// CHAR(0)
     narg = command_argument_count()
     if( narg>10 ) narg=10
    ! if( narg<1 ) stop "You should input control file name"

     do i=1, narg
         call get_command_argument( i, argv(i) )
         argv(i) = argv(i)// CHAR(0)
     enddo
	 
	! ierr = mw_initialize( narg, argv, path )
	 ierr = mw_initialize_fstr( narg, argv, cntpath )
	 
     print *, "end of mw init", ierr
	 
     myrank = mw_get_rank()
     nprocs = mw_get_num_of_process()
     hecmw_comm = mw_mpi_comm()
	 
     ierr = mw_file_read_fstr()
     print *, "end of mw read", ierr
	 
     nRefine = mw_get_fstr_refine_num()
     refinetype = mw_get_fstr_refine_type()
	 
     if( refinetype/=0 ) then
         argv_len = mw_get_fstr_filename_length_cadfit()
         call mw_get_fstr_filename_cadfit(fname, argv_len)
         fname = fname// CHAR(0)
         ierr =  mw_revocap_refine(fname, nRefine)
     else
         ierr = mw_refine(nRefine)
     endif
      print *, "end of mw refine with refined level", nRefine, "and refine type", refinetype
    ! ierr = mw_file_write()
	 
     nAss = mw_get_num_of_assemble_model()
     call mw_select_assemble_model( 0 )
	 nPart = mw_get_num_of_mesh_part()
	 
	 print *, "Number of Parts, Assmble:", nPart, nAss
     allocate( part_nodes(nAss, nPart+1) )	 
     allocate( part_elems(nAss, nPart+1) )
	 total_node = 0	; total_elem=0 
     print *, "Number of elements & nodes:"
     do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         part_nodes( iAss+1, 1 ) = total_node
         part_elems( iAss+1, 1 ) = total_elem
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            total_node = total_node+ mw_get_num_of_node()
            total_elem = total_elem+ mw_get_num_of_element()
            print *, "-Assmble, Part:",iAss, iPart, mw_get_num_of_element(), mw_get_num_of_node()
            part_nodes(iAss+1, iPart+2) = total_node
            part_elems(iAss+1, iPart+2) = total_elem
         enddo
     enddo
     print *, "Total number of elements & nodes:", total_elem, total_node
     allocate( global_node_ID(total_node) )
     allocate( global_elem_ID(total_elem) )
     cntn = 0; cnte = 0
     do iAss = 0, mw_get_num_of_assemble_model()-1
       call mw_select_assemble_model( iAss )
       do iPart = 0, mw_get_num_of_mesh_part()-1
         call mw_select_mesh_part( iPart )
         do i=0,mw_get_num_of_node()-1
           cntn = cntn+1
           global_node_ID(cntn) = mw_get_node_id(i)
         enddo
         do i=0,mw_get_num_of_element()-1
           cnte = cnte+1
           global_elem_ID(cnte) = mw_get_element_id(i)
         enddo
       enddo
     enddo
	 
     call cpu_time(T1)
	 
	 call fstr_init()
     call fstr_rcap_initialize( hecMESH, fstrPR, fstrCPL )

     nAss = mw_get_num_of_assemble_model()
   
	 call mw_gene_linear_algebra(1, 1, assDOF)
    ! do iAss = 0, mw_get_num_of_assemble_model()-1
    !     call mw_select_assemble_model( iAss )
    !     call mw_select_algebra(0)
    ! enddo
	 
	 CALL CPU_TIME(T2)
	 
	 ! =============== ANALYSIS =====================
        select case( fstrPR%solution_type )
        case ( kstPRECHECK )
                call fstr_precheck( )
        case ( kstSTATIC )
                call fstr_linear_static_analysis				
        case ( kstNLSTATIC )
                call fstr_nonlinear_static_analysis
        case ( kstEIGEN )
      !          call fstr_eigen_analysis
        case ( kstHEAT )
      !          call fstr_heat_analysis
        case ( kstDYNAMIC )
      !          call fstr_linear_dynamic_analysis
        end select

        CALL CPU_TIME(T3)

        write(ILOG,*)
        write(ILOG,*)           '===================================='
        write(ILOG,'(a,f10.2)') '    TOTAL TIME (sec) :', T3 - T1
        write(ILOG,'(a,f10.2)') '           pre (sec) :', T2 - T1
        write(ILOG,'(a,f10.2)') '         solve (sec) :', T3 - T2
        write(ILOG,*)           '===================================='


        ! =============== FINALIZE =====================
        call fstr_rcap_finalize( fstrPR, fstrCPL )
        call fstr_finalize()
        if(myrank==0) write(*,*) 'FrontISTR completed successfully!'
	 
	 
	 ierr = mw_finalize()
	 
contains

  subroutine fstr_init( )
        use m_heat_init
        implicit none
        character(len=HECMW_NAME_LEN) :: cntfileNAME 
        include "HEC_MW3_For.h"	
        integer :: nlen
		
        nlen = mw_get_fstr_filename_length_control()
        cntfileNAME(1:nlen)=' '
        call mw_get_fstr_filename_control(cntfileNAME(1:nlen), nlen)

        call fstr_nullify_fstr_param ( fstrPR     )
        call fstr_nullify_fstr_solid ( fstrSOLID  )
        call fstr_nullify_fstr_heat  ( fstrHEAT   )
        call fstr_nullify_lczparam   ( fstrEIG    )
        call fstr_nullify_fstr_dynamic ( fstrDYNAMIC  )
        call fstr_nullify_fstr_couple ( fstrCPL    )

        call fstr_init_file
        call fstr_init_condition( cntFileName(1:nlen) )
   
  end subroutine fstr_init


!> Open all files preparing calculation
  subroutine fstr_init_file
        implicit none
        character(len=HECMW_FILENAME_LEN) :: s
        character(len=HECMW_FILENAME_LEN) :: stafileNAME
        character(len=HECMW_FILENAME_LEN) :: logfileNAME
        character(len=HECMW_FILENAME_LEN) :: msgfileNAME
        character(len=HECMW_FILENAME_LEN) :: dbgfileNAME
        integer :: stat

        ! set file name --------------------------------

        write( logfileNAME, '(i5,''.log'')') myrank
        logfileNAME = adjustl(logfileNAME)
        stafileNAME = 'FSTR.sta'
        msgfileNAME = 'FSTR.msg'
        write(s,*) myrank
        write( dbgfileNAME, '(a,a)') 'FSTR.dbg.', trim(adjustl(s))
        dbgfileNAME = adjustl(dbgfileNAME)

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
        end if

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

!> Read in control file and do all preparation
  subroutine fstr_init_condition( cntfileName )
        implicit none
        character(len=*) :: cntfileNAME

        ! loading boundary conditions etc. from fstr control file or nastran mesh file
        ! and setup parameters ...
        call fstr_setup( cntfileNAME, hecMESH, fstrPR, fstrSOLID,  &
                           fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL ) 
						   
        write(*,*) 'fstr_setup: OK'; call flush(6)

  end subroutine fstr_init_condition
  

!=============================================================================!
!> Master subroutine of linear static analysis
!=============================================================================!

  subroutine fstr_linear_static_analysis
        use m_fstr_solve_LINEAR

        if( IECHO==1 ) call fstr_echo(fstrSOLID)

        if(myrank == 0) THEN
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*) ' ***   STAGE Linear static analysis   **'
        end if
        call fstr_solid_alloc( fstrSOLID )
        call fstr_solve_LINEAR( fstrEIG, fstrSOLID )
        call fstr_solid_finalize( fstrSOLID )

  end subroutine fstr_linear_static_analysis

!=============================================================================!
!> Master subroutine of nonlinear static analysis                                           !
!=============================================================================!

  subroutine fstr_nonlinear_static_analysis
        use m_fstr_solve_NLGEOM

        if( IECHO==1 ) call fstr_echo(fstrSOLID)
        if(myrank == 0) THEN
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*) ' ***   STAGE Non linear static analysis   **'
        end if
        call fstr_solid_alloc( fstrSOLID )
        call fstr_solve_NLGEOM( fstrSOLID )
        call fstr_solid_finalize( fstrSOLID )

  end subroutine fstr_nonlinear_static_analysis

!=============================================================================!
!> Master subroutine of eigen analysis                                                    !
!=============================================================================!

  subroutine fstr_eigen_analysis
        use m_fstr_solve_eigen
        if( IECHO==1 ) call fstr_echo(fstrSOLID)
        if(myrank == 0) THEN
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*) ' ***   STAGE Eigenvalue analysis     **'
        end if

        call fstr_solve_EIGEN( hecMESH, hecMAT, fstrEIG, fstrSOLID, fstrPR )

  end subroutine fstr_eigen_analysis

!=============================================================================!
!> Master subroutine of heat analysis                                                      !
!=============================================================================!

  subroutine fstr_heat_analysis
        use m_heat_echo
        use m_fstr_solve_heat
		
     !   call heat_init_material (hecMESH,fstrHEAT)
      !  call heat_init_amplitude(hecMESH,fstrHEAT)
		
        if( IECHO==1 ) call heat_echo(fstrPR,hecMESH,fstrHEAT)
        if(myrank ==0) THEN
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*) ' ***   STAGE Heat analysis    **'
        end if

        call fstr_solve_HEAT( hecMESH, hecMAT, fstrPR, fstrHEAT )

  end subroutine fstr_heat_analysis


!=============================================================================!
!> Master subroutine of dynamic analysis                                              !
!=============================================================================!

  subroutine fstr_linear_dynamic_analysis
        use fstr_solver_dynamic

        if( IECHO==1 ) call fstr_echo(fstrSOLID)

        if(myrank == 0) THEN
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*) ' ***   STAGE Linear dynamic analysis   **'
        end if

        call fstr_solve_dynamic( hecMESH, hecMAT,fstrSOLID,fstrEIG    &
                                        ,fstrDYNAMIC,fstrPR,fstrCPL)

  end subroutine fstr_linear_dynamic_analysis

!=============================================================================!
!> Finalizer                                                             !
!=============================================================================!
  subroutine fstr_finalize

        close(ILOG)

        if( myrank==0 ) then
                close(ISTA)
        end if
        if( myrank == 0 ) then
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*)
                write(IMSG,*) ':========================================:'
                write(IMSG,*) ':**            END of FSTR             **:'
                write(IMSG,*) ':========================================:'
                close(IMSG)

        end if

        call fstr_solid_finalize( fstrSOLID )
        if( associated(MWSections) ) deallocate(MWSections)
        if( allocated(part_nodes) ) deallocate(part_nodes)

        close(IDBG)
        close(IMSG)

  end subroutine fstr_finalize
	 
end program
