!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.4                                   !
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
      integer(kind=kint) :: i,j,i_flg_1, my_rank_monit_1, ierror
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

      my_rank_monit_1 = 999999999
      i_flg_1 = 0

      if( i_flg_1 .ne. 1 ) then
        do i= 1, hecMESH%nn_internal
        if( fstrPARAM%global_local_id(1,i) .eq. fstrDYNAMIC%node_monit_1 ) then
          fstrDYNAMIC%node_monit_1 = i
          i_flg_1 = 1
          do j= 1, hecMESH%nn_internal
            if( j == i ) then
              my_rank_monit_1 = hecMESH%node_ID(2*j)
              exit
            end if
          end do
          exit
        end if
        enddo
      end if

!C
!C-- file open for local use
!C
      if( hecMESH%my_rank == my_rank_monit_1 ) then
        OPEN(fstrDYNAMIC%dynamic_IW4,FILE='dyna_disp_p1.txt', status = 'replace', iostat=ierror)
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error <dyna_disp_p1.txt>'
          call hecmw_abort( hecmw_comm_get_comm())
        end if
        OPEN(fstrDYNAMIC%dynamic_IW5,FILE='dyna_velo_p1.txt', status = 'replace', iostat=ierror)
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error <dyna_velo_p1.txt>'
          call hecmw_abort( hecmw_comm_get_comm())
        end if
        OPEN(fstrDYNAMIC%dynamic_IW6,FILE='dyna_acce_p1.txt', status = 'replace', iostat=ierror)
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error <dyna_acce_p1.txt>'
          call hecmw_abort( hecmw_comm_get_comm())
        end if
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
        if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
             if(fstrDYNAMIC%restart_nout < 0 ) then
               call hecmw_restart_open()
               call hecmw_restart_read_int(restrt_step)
               restrt_step_num = restrt_step(1) + 1
               call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
               call hecmw_restart_read_real(fstrDYNAMIC%VEL (:,1))
               call hecmw_restart_read_real(fstrDYNAMIC%ACC (:,1))
               call hecmw_restart_close()
               fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
               hecMAT%Iarray(98) = 1
             endif
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                 call fstr_solve_dynamic_implicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, my_rank_monit_1, restrt_step_num )
             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    !write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                    call fstr_solve_frequency_analysis(hecMESH, hecMAT, fstrSOLID, myEIG, fstrDYNAMIC, &
                                                       fstrRESULT, fstrPARAM, fstrCPL, fstrFREQ, fstrMAT, my_rank_monit_1, &
                                                       restrt_step_num)
                 end if
             end if

        else if(fstrDYNAMIC%idx_eqa == 11) then  ! explicit dynamic analysis
             if(fstrDYNAMIC%restart_nout < 0 ) then
               call hecmw_restart_open()
               call hecmw_restart_read_int(restrt_step)
               restrt_step_num = restrt_step(1) + 1
               call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
               call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,3))
               call hecmw_restart_close()
               fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
               hecMAT%Iarray(98) = 1
             endif
             if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
                call fstr_solve_dynamic_explicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, my_rank_monit_1, restrt_step_num )
             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    !write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                    call fstr_solve_frequency_analysis(hecMESH, hecMAT, fstrSOLID, myEIG, fstrDYNAMIC, &
                                                       fstrRESULT, fstrPARAM, fstrCPL, fstrFREQ, fstrMAT, my_rank_monit_1, &
                                                       restrt_step_num)
                 end if
             end if
        end if

      else
        if(fstrDYNAMIC%restart_nout < 0 ) then
          if( .not. associated( fstrSOLID%contacts ) ) then
            call fstr_read_restart_dyna(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM)
          else                                                
            call fstr_read_restart_dyna(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
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
                                      ,fstrCPL, my_rank_monit_1, restrt_step_num )
               elseif( fstrPARAM%contact_algo == kcaSLagrange ) then
!                ----  For Parallel Contact with Multi-Partition Domains
                 if(paraContactFlag.and.present(conMAT)) then
                   call fstr_solve_dynamic_nlimplicit_contactSLag(1, hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL,fstrMAT,my_rank_monit_1,restrt_step_num,infoCTChange   &
                                      ,conMAT )
                 else
                   call fstr_solve_dynamic_nlimplicit_contactSLag(1, hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL,fstrMAT,my_rank_monit_1,restrt_step_num,infoCTChange )
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
                                      ,fstrCPL, my_rank_monit_1, restrt_step_num )

             else if(fstrDYNAMIC%idx_resp == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if
        end if
      endif

!C-- file close for local use
      if( hecMESH%my_rank == my_rank_monit_1 ) then
        CLOSE(fstrDYNAMIC%dynamic_IW4)
        CLOSE(fstrDYNAMIC%dynamic_IW5)
        CLOSE(fstrDYNAMIC%dynamic_IW6)
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
