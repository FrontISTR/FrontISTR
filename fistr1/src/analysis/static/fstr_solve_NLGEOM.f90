!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Sato(Advancesoft), X. YUAN(AdavanceSoft)    !
!                       Z. Sun(ASTOM)                                  !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief This module provides main suboruitne for nonliear calculation.
!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \date       2009/01/30
!>  \author     Z. Sun(ASTOM)
!>  \date       2010/11
!>  \version    0.00
!
!======================================================================!
module m_fstr_solve_NLGEOM

use m_fstr
use m_static_lib
use m_static_output

use m_fstr_NonLinearMethod
use m_fstr_Restart

use fstr_matrix_con_contact

    implicit none

contains

!======================================================================!
!> \brief This module provides main suborutine for nonlinear calculation.
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>              Z. Sun(ASTOM)/2010/11
!>  \date       2009/08/31
!>  \version    0.00

      subroutine FSTR_SOLVE_NLGEOM(hecMESH,hecMAT,fstrSOLID,fstrMAT,fstrPARAM,conMAT)    
!======================================================================!
      type (hecmwST_local_mesh)              :: hecMESH      !< mesh information
      type (hecmwST_matrix    )              :: hecMAT       !< linear equation, its right side modified here
      type (fstr_param       )               :: fstrPARAM    !< analysis control parameters
      type (fstr_solid       )               :: fstrSOLID    !< we need boundary conditions of curr step
      type (fstrST_matrix_contact_lagrange)  :: fstrMAT      !< type fstrST_matrix_contact_lagrange                                
      type (fstr_info_contactChange)         :: infoCTChange !< type fstr_info_contactChange                                            
      type (hecmwST_matrix    ),optional     :: conMAT

      integer(kind=kint) :: ndof, nn

      integer(kind=kint) :: j, tot_step, step_count
      integer(kind=kint) :: sub_step
      real(kind=kreal) :: tt, factor
      real(kind=kreal) :: time_1, time_2
      logical          :: ctchanged                             
      integer(kind=kint) :: restart_step_num, restart_substep_num

      hecMAT%NDOF = hecMESH%n_dof

      ndof = hecMAT%NDOF
      nn = ndof*ndof
	  
      if( fstrSOLID%TEMP_ngrp_tot>0 .and. hecMESH%hecmw_flag_initcon==1 ) then
        do j=1, hecMESH%n_node
          fstrSOLID%last_temp(j) = hecMESH%node_init_val_item(j)
          fstrSOLID%temperature(j) = hecMESH%node_init_val_item(j)
        end do
      endif

      if( fstrSOLID%restart_nout == 0 ) then
        fstrSOLID%restart_nout = 999999999
      end if

      restart_step_num = 1
      restart_substep_num = 1
      fstrSOLID%unode = 0.0
      step_count = 0 !**
      infoCTChange%contactNode_previous = 0
      if(fstrSOLID%restart_nout <0 ) then
        if( .not. associated( fstrSOLID%contacts ) ) then   
          call fstr_read_restart(restart_step_num,restart_substep_num,step_count,hecMESH,fstrSOLID,fstrPARAM)            
        else                                                
          call fstr_read_restart(restart_step_num,restart_substep_num,step_count,hecMESH,fstrSOLID,fstrPARAM, & 
                                 infoCTChange%contactNode_previous)     
        endif
        hecMAT%Iarray(98) = 1
      endif

      fstrSOLID%FACTOR    =0.0
                 
      do tot_step=restart_step_num, fstrSOLID%nstep_tot
        if(hecMESH%my_rank==0) write(*,*) ''
        if(hecMESH%my_rank==0) write(*,'(a,i5)') ' loading step=',tot_step
		
        if( fstrSOLID%TEMP_ngrp_tot>0 ) then
           do j=1, hecMESH%n_node
             fstrSOLID%temp_bak(j) = fstrSOLID%temperature(j)
           end do
        endif
        
! -------------------------------------------------------------------------
!      STEP LOOP
! -------------------------------------------------------------------------
        do sub_step = restart_substep_num, fstrSOLID%step_ctrl(tot_step)%num_substep

! ----- time history of factor
          tt = sub_step*fstrSOLID%step_ctrl(tot_step)%initdt          !**
          call table_nlsta(hecMESH,fstrSOLID,tot_step,sub_step-1,factor)
            fstrSOLID%FACTOR(1) = factor
          call table_nlsta(hecMESH,fstrSOLID,tot_step,sub_step, factor)
            fstrSOLID%FACTOR(2) = factor

          if( dabs(fstrSOLID%FACTOR(2)-fstrSOLID%FACTOR(1))<1.0e-20 ) then
            if( hecMESH%my_rank==0) then
               write(imsg,*) 'loading increment finished due to f_2-f_1 = 0.0'
            end if
            exit
          end if
		  	  
          if(hecMESH%my_rank==0) write(*,*) ' substep=',sub_step,fstrSOLID%FACTOR
          step_count = step_count+1
          call cpu_time(time_1)
		  
          call fstr_UpdateState( hecMESH, fstrSOLID, 0.d0 )

!       analysis algorithm ( Newton-Rapshon Method )
          if( .not. associated( fstrSOLID%contacts ) ) then        
            call fstr_Newton( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM,   &           
                              restart_step_num, sub_step  )   
          else
            if( fstrPARAM%contact_algo == kcaSLagrange ) then
!     ----  For Parallel Contact with Multi-Partition Domains
              if(paraContactFlag.and.present(conMAT)) then
                call fstr_Newton_contactSLag( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT,  &            
                                    restart_step_num, restart_substep_num, sub_step, infoCTChange, conMAT )               
              else
                call fstr_Newton_contactSLag( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT,  &            
                                    restart_step_num, restart_substep_num, sub_step, infoCTChange )                   
              endif
            else if( fstrPARAM%contact_algo == kcaALagrange ) then                              
              call fstr_Newton_contactALag( tot_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM,            &     
                                  restart_step_num, restart_substep_num, sub_step, infoCTChange )                                  
            endif                                                                                
          endif   

          if( fstrSOLID%TEMP_irres==0 ) then
            if(fstrSOLID%restart_nout<0) then
              fstrSOLID%restart_nout = - fstrSOLID%restart_nout
            end if
            if( mod(step_count,fstrSOLID%restart_nout) == 0 ) then
              if( .not. associated( fstrSOLID%contacts ) ) then                              
                call fstr_write_restart(tot_step,sub_step,step_count,hecMESH,fstrSOLID,fstrPARAM)      
              else                                                                                                                          
                call fstr_write_restart(tot_step,sub_step,step_count,hecMESH,fstrSOLID,fstrPARAM, &    
                                      infoCTChange%contactNode_current)                  
              endif    
            end if  

! ----- Result output (include visualize output)
            call fstr_static_Output( tot_step, step_count, hecMESH, fstrSOLID, fstrPR%solution_type )
          endif

          call cpu_time(time_2)
          if( hecMESH%my_rank==0) then
            write(ISTA,'(a,f10.2)') '         solve (sec) :', time_2 - time_1
          end if

        enddo    !--- end of substep  loop 

	    restart_substep_num = 1

      enddo      !--- end of tot_step loop
!
!  message
! 
      IF(myrank == 0)THEN
        WRITE(IMSG,'("### FSTR_SOLVE_NLGEOM FINISHED!")')
        WRITE(*,'("### FSTR_SOLVE_NLGEOM FINISHED!")')
      ENDIF

      end subroutine FSTR_SOLVE_NLGEOM
	  

!C================================================================C
!> \brief This subroutine decide the loading increment considering
!>        the amplitude definition
!C================================================================C
    subroutine table_nlsta(hecMESH, fstrSOLID, cstep, substep, f_t)
    type ( hecmwST_local_mesh ), intent(in) :: hecMESH    !< hecmw mesh
    type ( fstr_solid         ), intent(in) :: fstrSOLID  !< fstr_solid
    integer(kind=kint), intent(in)          :: cstep      !< curr loading step
    integer(kind=kint), intent(in)          :: substep    !< curr substep number
    real(kind=kreal), intent(out)           :: f_t        !< loading factor

    integer(kind=kint) :: i
    integer(kind=kint) :: jj_n_amp, jj1, jj2
    integer(kind=kint) :: s1, s2, flag
    real(kind=kreal) :: t_1, t_2, t_t, f_1, f_2, tincre


      jj_n_amp = fstrSOLID%step_ctrl( cstep )%amp_id

      if( jj_n_amp <= 0 ) then  ! Amplitude not defined
          tincre = 1.d0 / fstrSOLID%step_ctrl( cstep )%num_substep
          f_t = tincre*substep
          if( f_t>1.d0 ) f_t=1.d0

      else
          tincre = fstrSOLID%step_ctrl( cstep )%initdt
          jj1 = hecMESH%amp%amp_index(jj_n_amp - 1)
          jj2 = hecMESH%amp%amp_index(jj_n_amp)

          jj1 = jj1 + 2
          t_t = tincre*substep

!      if(jj2 .eq. 0) then
!         f_t = 1.0
          if(t_t .gt. hecMESH%amp%amp_table(jj2)) then
             f_t = hecMESH%amp%amp_val(jj2)
          else if(t_t .le. hecMESH%amp%amp_table(jj2)) then
             flag=0
             do i = jj1, jj2
                if(t_t .le. hecMESH%amp%amp_table(i)) then
                  s2 = i
                  s1 = i - 1
                  flag = 1
                end if
                if( flag == 1 ) exit
             end do

             t_2 = hecMESH%amp%amp_table(s2)
             t_1 = hecMESH%amp%amp_table(s1)
             f_2 = hecMESH%amp%amp_val(s2)
             f_1 = hecMESH%amp%amp_val(s1)
             if( t_2-t_1 .lt. 1.0e-20) then
              if(myrank == 0) then
                 write(imsg,*) 'stop due to t_2-t_1 <= 0'
              end if
              call hecmw_abort( hecmw_comm_get_comm())
            end if
            f_t = ((t_2*f_1 - t_1*f_2) + (f_2 - f_1)*t_t) / (t_2 - t_1)
          end if

      end if

  end subroutine table_nlsta

!
end module m_fstr_solve_NLGEOM
