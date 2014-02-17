!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                       Zhigang Sun(ASTOM)                             !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines controlling dynamic calculation

module fstr_solver_dynamic

use m_fstr
use fstr_dynamic_implicit
use fstr_dynamic_explicit
use fstr_dynamic_nlexplicit
use fstr_dynamic_nlimplicit
use fstr_matrix_con_contact
use fstr_frequency_analysis  !Frequency analysis module

contains

!C================================================================C
!> Master subroutine for dynamic analysis
!C================================================================C
  subroutine fstr_solve_DYNAMIC(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &      
                                      ,fstrCPL,fstrFREQ, fstrMAT  &
                                      ,conMAT )                  
      use m_fstr_setup
      implicit none
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( lczparam            ) :: myEIG
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM
      type ( fstr_dynamic        ) :: fstrDYNAMIC
      type ( fstr_couple         ) :: fstrCPL         !for COUPLE
      type ( fstr_freqanalysis   ) :: fstrFREQ
      type (fstrST_matrix_contact_lagrange)  :: fstrMAT      !< type fstrST_matrix_contact_lagrange                                 
      type (fstr_info_contactChange)         :: infoCTChange !< type fstr_info_contactChange
      type ( hecmwST_matrix      ),optional :: conMAT
      integer(kind=kint) :: i,j,num_monit,ig,iS,iE,ik,in,ing,iunitS,iunit,ierror
      character(len=HECMW_FILENAME_LEN) :: fname
      integer(kind=kint) :: restrt_step_num
      integer(kind=kint) :: restrt_step(1)


      if(dabs(fstrDYNAMIC%t_delta) < 1.0e-20) then
              if( hecMESH%my_rank == 0 ) then
                 write(imsg,*) 'stop due to fstrDYNAMIC%t_delta = 0'
              end if
              call hecmw_abort( hecmw_comm_get_comm())
      end if
      call fstr_dynamic_alloc( hecMESH, fstrDYNAMIC )
      allocate( myEIG%mass(hecMAT%NDOF*hecMESH%n_node)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, mass>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

!C
!C-- file open for local use
!C
      if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
      num_monit = 0
      ig = fstrDYNAMIC%ngrp_monit
      iS = hecMESH%node_group%grp_index(ig-1)+1
      iE = hecMESH%node_group%grp_index(ig)
      do ik=iS,iE
        in = hecMESH%node_group%grp_item(ik)
        if (in > hecMESH%nn_internal) cycle
        num_monit = num_monit+1
        ing = hecMESH%global_node_id(in)
        iunitS = 10*(num_monit-1)

        iunit = iunitS + fstrDYNAMIC%dynamic_IW4
        write(fname,'(a,i0,a)') 'dyna_disp_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          OPEN(iunit,FILE=fname, access = 'append', iostat=ierror)
        else
          OPEN(iunit,FILE=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW5
        write(fname,'(a,i0,a)') 'dyna_velo_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          OPEN(iunit,FILE=fname, access = 'append', iostat=ierror)
        else
          OPEN(iunit,FILE=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW6
        write(fname,'(a,i0,a)') 'dyna_acce_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          OPEN(iunit,FILE=fname, access = 'append', iostat=ierror)
        else
          OPEN(iunit,FILE=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW8
        write(fname,'(a,i0,a)') 'dyna_strain_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          OPEN(iunit,FILE=fname, access = 'append', iostat=ierror)
        else
          OPEN(iunit,FILE=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW9
        write(fname,'(a,i0,a)') 'dyna_stress_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          OPEN(iunit,FILE=fname, access = 'append', iostat=ierror)
        else
          OPEN(iunit,FILE=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if
      enddo
      endif

!C
!C-- initial condition
!C
	  fstrDYNAMIC%DISP = 0.d0
      fstrDYNAMIC%VEL  = 0.d0
      fstrDYNAMIC%ACC  = 0.d0
      fstrSOLID%unode(:)  =0.d0
	  fstrSOLID%QFORCE(:) =0.d0
	  fstrDYNAMIC%kineticEnergy=0.d0
	  fstrDYNAMIC%strainEnergy=0.d0
	  fstrDYNAMIC%totalEnergy=0.d0
      call fstr_UpdateState( hecMESH, fstrSOLID, 0.d0 )

! ---- restart 

      restrt_step_num = 1
      fstrDYNAMIC%i_step = 0
      infoCTChange%contactNode_previous = 0

      if(fstrDYNAMIC%restart_nout >= 0 ) then
        call dynamic_bc_init   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        call dynamic_bc_init_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        call dynamic_bc_init_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
      endif

      if(fstrDYNAMIC%nlflag==0) then
        if(fstrDYNAMIC%restart_nout < 0 ) then
          call fstr_read_restart_dyna_linear(restrt_step_num,fstrDYNAMIC)
          restrt_step_num = restrt_step_num + 1
          fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
          hecMAT%Iarray(98) = 1
        endif
        if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                 call fstr_solve_dynamic_implicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, restrt_step_num )
             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    !write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                    call fstr_solve_frequency_analysis(hecMESH, hecMAT, fstrSOLID, myEIG, fstrDYNAMIC, &
                                                       fstrRESULT, fstrPARAM, fstrCPL, fstrFREQ, fstrMAT, &
                                                       restrt_step_num)
                 end if
             end if

        else if(fstrDYNAMIC%idx_eqa == 11) then  ! explicit dynamic analysis
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                call fstr_solve_dynamic_explicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, restrt_step_num )
             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    !write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                    call fstr_solve_frequency_analysis(hecMESH, hecMAT, fstrSOLID, myEIG, fstrDYNAMIC, &
                                                       fstrRESULT, fstrPARAM, fstrCPL, fstrFREQ, fstrMAT, &
                                                       restrt_step_num)
                 end if
             end if
        end if

      else
        if(fstrDYNAMIC%restart_nout < 0 ) then
          if( .not. associated( fstrSOLID%contacts ) ) then
            call fstr_read_restart_dyna_nl(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM)
          else                                                
            call fstr_read_restart_dyna_nl(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
                                        infoCTChange%contactNode_previous)
          endif
          restrt_step_num = restrt_step_num + 1
          fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
          hecMAT%Iarray(98) = 1
        end if
        if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
               if(.not. associated( fstrSOLID%contacts ) ) then
                 call fstr_solve_dynamic_nlimplicit(1, hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, restrt_step_num )
               elseif( fstrPARAM%contact_algo == kcaSLagrange ) then
!                ----  For Parallel Contact with Multi-Partition Domains
                 if(paraContactFlag.and.present(conMAT)) then
                   call fstr_solve_dynamic_nlimplicit_contactSLag(1, hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL,fstrMAT,restrt_step_num,infoCTChange   &
                                      ,conMAT )
                 else
                   call fstr_solve_dynamic_nlimplicit_contactSLag(1, hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL,fstrMAT,restrt_step_num,infoCTChange )
                 endif
               endif                                                                                                              
             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if

        else if(fstrDYNAMIC%idx_eqa == 11) then  ! explicit dynamic analysis
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                call fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, restrt_step_num )

             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if
        end if
      endif

!C-- file close for local use
      if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
      do i=1,num_monit
        iunitS = 10*(i-1)
        CLOSE(iunitS + fstrDYNAMIC%dynamic_IW4)
        CLOSE(iunitS + fstrDYNAMIC%dynamic_IW5)
        CLOSE(iunitS + fstrDYNAMIC%dynamic_IW6)
        CLOSE(iunitS + fstrDYNAMIC%dynamic_IW8)
        CLOSE(iunitS + fstrDYNAMIC%dynamic_IW9)
      enddo
      endif
!C-- end of finalization

      call fstr_dynamic_finalize( fstrDYNAMIC )
      call hecMAT_finalize( hecMAT )

	      deallocate( myEIG%mass    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, mass>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

  end subroutine fstr_solve_DYNAMIC

end module fstr_solver_dynamic
