!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                       Zhigang Sun(ASTOM)                                  !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines for nonlinear implicit dynamic analysis

module fstr_dynamic_nlimplicit

use m_fstr
use lczparm
use m_dynamic_output
use m_fstr_EIG_setMASS
use m_dynamic_mat_ass_bc_ac
use m_dynamic_mat_ass_bc
use m_dynamic_mat_ass_bc_vl
use m_dynamic_mat_ass_load
use m_fstr_StiffMatrix
use m_dynamic_post
use m_fstr_Update
use m_fstr_Restart
use fstr_matrix_con_contact                   
use m_fstr_Residual                           

!-------- for couple -------
use m_dynamic_mat_ass_couple
use m_fstr_rcap_io


contains

!> \brief This subroutine provides function of nonlinear implicit dynamic analysis using the Newmark method.

  subroutine fstr_solve_dynamic_nlimplicit(cstep, hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, my_rank_monit_1, restrt_step_num )

    implicit none
!C
!C-- global variable
!C
    integer, intent(in)          :: cstep     !< current step  
    type ( hecmwST_local_mesh  ) :: hecMESH
    type ( hecmwST_matrix      ) :: hecMAT
    type ( lczparam            ) :: myEIG
    type ( fstr_solid          ) :: fstrSOLID
    type ( hecmwST_result_data ) :: fstrRESULT
    type ( fstr_param          ) :: fstrPARAM
    type ( fstr_dynamic        ) :: fstrDYNAMIC
    type (fstrST_matrix_contact_lagrange)  :: fstrMAT      !< type fstrST_matrix_contact_lagrange  
    type ( fstr_couple         ) :: fstrCPL         !for COUPLE        

!C
!C-- local variable
!C
    

    integer(kind=kint) :: nnod, ndof, numnp, nn
    integer(kind=kint) :: i, j, ids, ide, ims, ime, kk, idm, imm
    integer(kind=kint) :: iter
    integer(kind=kint) :: my_rank_monit_1
    

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize, res
    real :: time_1, time_2

    integer(kind=kint) :: restrt_step_num


    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof
!!
!!-- initial value
!!
      call cpu_time(time_1)

!C
!C-- check parameters
!C
         if(dabs(fstrDYNAMIC%beta) < 1.0e-20) then
              if( hecMESH%my_rank == 0 ) then
                 write(imsg,*) 'stop due to Newmark-beta = 0'
              end if
              call hecmw_abort( hecmw_comm_get_comm())
         end if


!C-- matrix [M]
!C-- lumped mass matrix
             if(fstrDYNAMIC%idx_mas == 1) then

                call setMASS(IDBG,hecMESH,hecMAT,myEIG)

!C-- consistent mass matrix
               else if(fstrDYNAMIC%idx_mas == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: consistent mass matrix is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if
!C--
      hecMAT%Iarray(98) = 1   !Assmebly complete
      hecMAT%Iarray(97) = 1   !Need numerical factorization
	  
!C
!C
!C-- time step loop
!C
    a1 = .5d0/fstrDYNAMIC%beta - 1.d0
    a2 = 1.d0/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)
    a3 = 1.d0/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta*fstrDYNAMIC%t_delta)
    b1 = ( .5d0*fstrDYNAMIC%ganma/fstrDYNAMIC%beta - 1.d0 )*fstrDYNAMIC%t_delta
    b2 = fstrDYNAMIC%ganma/fstrDYNAMIC%beta - 1.d0
    b3 = fstrDYNAMIC%ganma/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)
    c1 = 1.d0 + fstrDYNAMIC%ray_k*b3
    c2 = a3 + fstrDYNAMIC%ray_m*b3


!C-- output of initial state
    if( restrt_step_num == 1 ) then
      call dynamic_nloutput(0,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, my_rank_monit_1)   
    end if
	
	fstrDYNAMIC%VEC3(:) =0.d0
	hecMAT%X(:) =0.d0
	
!!
!!    step = 1,2,....,fstrDYNAMIC%n_step
!!
    do i= restrt_step_num, fstrDYNAMIC%n_step

       fstrDYNAMIC%i_step = i                                                  
       fstrDYNAMIC%t_curr = fstrDYNAMIC%t_curr + fstrDYNAMIC%t_delta                      
	
       if(hecMESH%my_rank==0) then
         write(ISTA,'('' time step='',i10,'' time='',1pe13.4e3)') i,fstrDYNAMIC%t_curr   
         if( mod(i,int(fstrDYNAMIC%nout/10)) == 0 )   &                        
            write(*,'('' time step='',i10,'' time='',1pe13.4e3)') i,fstrDYNAMIC%t_curr    
       endif

       fstrSOLID%dunode(:) =0.d0
      ! call fstr_UpdateEPState( hecMESH, fstrSOLID )  
	   
       do j = 1 ,ndof*nnod
           fstrDYNAMIC%VEC1(j) = a1*fstrDYNAMIC%ACC(j,1) + a2*fstrDYNAMIC%VEL(j,1)
           fstrDYNAMIC%VEC2(j) = b1*fstrDYNAMIC%ACC(j,1) + b2*fstrDYNAMIC%VEL(j,1)
       end do
	   
       do iter = 1, fstrSOLID%step_ctrl(cstep)%max_iter
         call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC%t_delta )

         if( fstrDYNAMIC%ray_k/=0.d0 .or. fstrDYNAMIC%ray_m/=0.d0 ) then 
           do j = 1 ,ndof*nnod
              hecMAT%X(j) = fstrDYNAMIC%VEC2(j) - b3*fstrSOLID%dunode(j)
           end do		 
         endif
         if( fstrDYNAMIC%ray_k/=0.d0 ) then
           if( hecMESH%n_dof == 3 ) then
             call hecmw_matvec_33 (hecMESH, hecMAT, hecMAT%X, fstrDYNAMIC%VEC3)
           else if( hecMESH%n_dof == 2 ) then
             call hecmw_matvec_22 (hecMESH, hecMAT, hecMAT%X, fstrDYNAMIC%VEC3, nnod)
           else if( hecMESH%n_dof == 6 ) then
             call matvec(fstrDYNAMIC%VEC3, hecMAT%X, hecMAT, ndof, hecMAT%D, hecMAT%AU, hecMAT%AL)
           end if	   
         endif
!C
!C-- mechanical boundary condition
         call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
         do j=1, hecMESH%n_node*  hecMESH%n_dof
           hecMAT%B(j)=hecMAT%B(j)- fstrSOLID%QFORCE(j) + myEIG%mass(j)*( fstrDYNAMIC%VEC1(j)-a3*fstrSOLID%dunode(j)   &
		     + fstrDYNAMIC%ray_m* hecMAT%X(j) ) + fstrDYNAMIC%ray_k*fstrDYNAMIC%VEC3(j)
         end do
		 
         do j = 1 ,nn*hecMAT%NP
           hecMAT%D(j)  = c1* hecMAT%D(j)
         end do
         do j = 1 ,nn*hecMAT%NPU
           hecMAT%AU(j) = c1* hecMAT%AU(j)
         end do
         do j = 1 ,nn*hecMAT%NPL
           hecMAT%AL(j) = c1*hecMAT%AL(j)
         end do
         do j=1,nnod
         do kk=1,ndof
           idm = nn*(j-1)+1 + (ndof+1)*(kk-1)
           imm = ndof*(j-1) + kk
           hecMAT%D(idm) = hecMAT%D(idm) + c2*myEIG%mass(imm)  
         end do
         end do


!C ********************************************************************************
!C for couple analysis
        if( fstrPARAM%fg_couple == 1) then
            if( fstrDYNAMIC%i_step > 1 .or. &
                (fstrDYNAMIC%i_step==1 .and. fstrPARAM%fg_couple_first==1 )) then
                  call fstr_rcap_get( fstrCPL )
                  call dynamic_mat_ass_couple( hecMESH, hecMAT, fstrSOLID, fstrCPL )
            endif
        endif
!C ********************************************************************************

!C-- geometrical boundary condition

        call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, iter)
        call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, iter)
        call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, iter)

! ----- check convergence
		res = dot_product( hecMAT%B, hecMAT%B )
        res = sqrt(res)/hecMESH%n_node
        if( hecMESH%my_rank==0 ) then
          if( mod(i,int(fstrDYNAMIC%nout/10)) == 0 )   &
             write(*,'(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res
          write(ISTA,'(''iter='',I5,''- Residual'',2E15.7)')iter,res
        endif
        if( res<fstrSOLID%step_ctrl(cstep)%converg ) exit
		
        CALL solve_LINEQ(hecMESH,hecMAT,imsg)
	
! ----- update the strain, stress, and internal force
        call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,fstrDYNAMIC%t_delta,1,fstrDYNAMIC%strainEnergy ) 
        do j=1,hecMESH%n_node*ndof
          fstrSOLID%dunode(j)  = fstrSOLID%dunode(j)+hecMAT%X(j)
        enddo
		
      enddo
	  
! -----  not convergence
      if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
         if( hecMESH%my_rank==0) then
           write(ILOG,*) '### Fail to Converge  : at step=', i
           write(ISTA,*) '### Fail to Converge  : at step=', i
           write(   *,*) '     ### Fail to Converge  : at step=', i
         end if
         stop
      end if
	  
!
!C-- new displacement, velocity and accelaration
!C
      fstrDYNAMIC%kineticEnergy = 0.0d0
      do j = 1 ,ndof*nnod
          fstrDYNAMIC%ACC (j,2) = -a1*fstrDYNAMIC%ACC(j,1) - a2*fstrDYNAMIC%VEL(j,1) + &
                               a3*fstrSOLID%dunode(j)
          fstrDYNAMIC%VEL (j,2) = -b1*fstrDYNAMIC%ACC(j,1) - b2*fstrDYNAMIC%VEL(j,1) + &
                               b3*fstrSOLID%dunode(j)
          fstrDYNAMIC%ACC (j,1) = fstrDYNAMIC%ACC (j,2)
          fstrDYNAMIC%VEL (j,1) = fstrDYNAMIC%VEL (j,2)

          fstrSOLID%unode(j)  = fstrSOLID%unode(j)+fstrSOLID%dunode(j)
          fstrDYNAMIC%DISP(j,2) = fstrSOLID%unode(j)
          
          fstrDYNAMIC%kineticEnergy = fstrDYNAMIC%kineticEnergy + &                        
                                      0.5d0*myEIG%mass(j)*fstrDYNAMIC%VEL(j,2)*fstrDYNAMIC%VEL(j,2) 
      enddo

!---  Restart info
      if( mod(i,fstrDYNAMIC%restart_nout) == 0 ) then
        call fstr_write_restart_dyna(i,0,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM)    
      end if
!
!C-- output new displacement, velocity and accelaration
      call dynamic_nloutput(i,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
!C
!C-- output result of monitoring node
!C
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, my_rank_monit_1)   
	  
      call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYNAMIC%t_delta )

    enddo
!C
!C-- end of time step loop

    call cpu_time(time_2)
    if( hecMESH%my_rank == 0 ) then
        write(ISTA,'(a,f10.2,a)') '         solve (sec) :', time_2 - time_1, 's'
    end if

  end subroutine fstr_solve_dynamic_nlimplicit
  
!> \brief This subroutine provides function of nonlinear implicit dynamic analysis using the Newmark method. 
!> Standard Lagrange multiplier algorithm for contact analysis is included in this subroutine.  
 subroutine fstr_solve_dynamic_nlimplicit_contactSLag(cstep, hecMESH,hecMAT,fstrSOLID,myEIG   &  
                                                       ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                                       ,fstrCPL,fstrMAT,my_rank_monit_1,restrt_step_num,infoCTChange )
  
    use mContact                                                                                                                   
    use m_addContactStiffness                              
    use m_solve_LINEQ_contact 
    use m_dynamic_init_variables                                                
      
    implicit none
!C
!C-- global variable
!C
    integer, intent(in)          :: cstep     !< current step     
    type ( hecmwST_local_mesh  ) :: hecMESH
    type ( hecmwST_matrix      ) :: hecMAT
    type ( lczparam            ) :: myEIG
    type ( fstr_solid          ) :: fstrSOLID
    type ( hecmwST_result_data ) :: fstrRESULT
    type ( fstr_param          ) :: fstrPARAM
    type ( fstr_dynamic        ) :: fstrDYNAMIC
    type ( fstr_couple         ) :: fstrCPL         !for COUPLE
    type (fstrST_matrix_contact_lagrange)  :: fstrMAT      !< type fstrST_matrix_contact_lagrange  
    type (fstr_info_contactChange)         :: infoCTChange !< fstr_info_contactChange                              

!C
!C-- local variable
!C
    
    integer(kind=kint) :: nnod, ndof, numnp, nn
    integer(kind=kint) :: i, j, ids, ide, ims, ime, kk, idm, imm
    integer(kind=kint) :: iter
    integer(kind=kint) :: my_rank_monit_1    
    

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize, res, res1, rf                                
    real :: time_1, time_2

    integer(kind=kint) :: restrt_step_num
    
    integer(kind=kint) :: ctAlgo                                      
    integer(kind=kint) :: max_iter_contact, count_step                 
    integer(kind=kint) :: stepcnt                                      
    real(kind=kreal)    :: maxDLag                                             

    logical :: is_mat_symmetric
    integer(kind=kint) :: n_node_global
    integer(kind=kint) :: contact_changed_global

! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

    ctAlgo = fstrPARAM%contact_algo                     
       
    if( hecMAT%Iarray(99)==4 .and. .not.fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH) ) then
      write(*,*) ' This type of direct solver is not yet available in such case ! '
      write(*,*) ' Please use intel MKL direct solver !'
      call  hecmw_abort(hecmw_comm_get_comm())
    endif   

    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof      
!!
!!-- initial value
!!
      call cpu_time(time_1)
!C
!C-- check parameters
!C
         if(dabs(fstrDYNAMIC%beta) < 1.0e-20) then
              if( hecMESH%my_rank == 0 ) then
                 write(imsg,*) 'stop due to Newmark-beta = 0'
              end if
              call hecmw_abort( hecmw_comm_get_comm())
         end if


!C-- matrix [M]
!C-- lumped mass matrix
             if(fstrDYNAMIC%idx_mas == 1) then

                call setMASS(IDBG,hecMESH,hecMAT,myEIG)

!C-- consistent mass matrix
               else if(fstrDYNAMIC%idx_mas == 2) then
                 if( hecMESH%my_rank .eq. 0 ) then
                    write(imsg,*) 'stop: consistent mass matrix is not yet available !'
                 end if
                 call hecmw_abort( hecmw_comm_get_comm())
             end if
!C--
      hecMAT%Iarray(98) = 1   !Assmebly complete
      hecMAT%Iarray(97) = 1   !Need numerical factorization
!C
!C-- initialize variables 
!C     
      if( restrt_step_num == 1 .and. fstrDYNAMIC%VarInitialize .and. fstrDYNAMIC%ray_m /= 0.0d0 ) & 
      call dynamic_init_varibles( hecMESH, hecMAT, fstrSOLID, myEIG, fstrDYNAMIC )                     
!C
!C
!C-- time step loop
!C
    a1 = .5d0/fstrDYNAMIC%beta - 1.d0
    a2 = 1.d0/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)
    a3 = 1.d0/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta*fstrDYNAMIC%t_delta)
    b1 = ( .5d0*fstrDYNAMIC%ganma/fstrDYNAMIC%beta - 1.d0 )*fstrDYNAMIC%t_delta
    b2 = fstrDYNAMIC%ganma/fstrDYNAMIC%beta - 1.d0
    b3 = fstrDYNAMIC%ganma/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)
    c1 = 1.d0 + fstrDYNAMIC%ray_k*b3
    c2 = a3 + fstrDYNAMIC%ray_m*b3
      

!C-- output of initial state
    if( restrt_step_num == 1 ) then
      call dynamic_nloutput(0,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, my_rank_monit_1)    
    end if
	
	fstrDYNAMIC%VEC3(:) =0.d0  
	hecMAT%X(:) =0.d0

    call fstr_save_originalMatrixStructure(hecMAT)   
    call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )   
    if ( fstr_is_contact_active() ) then                                                 
       call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange)    
    elseif( hecMAT%Iarray(99)==4 ) then                                                   
       write(*,*) ' This type of direct solver is not yet available in such case ! '
       write(*,*) ' Please change solver type to intel MKL direct solver !'
       call  hecmw_abort(hecmw_comm_get_comm()) 
    endif  
    is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)
    call solve_LINEQ_contact_init(hecMESH,hecMAT,fstrMAT,is_mat_symmetric)
	
!!
!!    step = 1,2,....,fstrDYNAMIC%n_step
!!
    do i= restrt_step_num, fstrDYNAMIC%n_step             

       fstrDYNAMIC%i_step = i                                                  
       fstrDYNAMIC%t_curr = fstrDYNAMIC%t_curr + fstrDYNAMIC%t_delta                      
	
       if(hecMESH%my_rank==0) then
         write(ISTA,'('' time step='',i10,'' time='',1pe13.4e3)') i,fstrDYNAMIC%t_curr   
         if( mod(i,int(fstrDYNAMIC%nout/10)) == 0 )   &                        
            write(*,'('' time step='',i10,'' time='',1pe13.4e3)') i,fstrDYNAMIC%t_curr      
       endif
       
       fstrSOLID%dunode(:) =0.d0
      ! call fstr_UpdateEPState( hecMESH, fstrSOLID )  
 	   
       do j = 1 ,ndof*nnod                  
           fstrDYNAMIC%VEC1(j) = a1*fstrDYNAMIC%ACC(j,1) + a2*fstrDYNAMIC%VEL(j,1)
           fstrDYNAMIC%VEC2(j) = b1*fstrDYNAMIC%ACC(j,1) + b2*fstrDYNAMIC%VEL(j,1)  
       end do                              

       max_iter_contact = 6 !1              
       count_step = 0                                              
       stepcnt = 0                                                 
       loopFORcontactAnalysis: DO WHILE( .TRUE. )                                                                                 
         count_step = count_step + 1                               
         do iter = 1, fstrSOLID%step_ctrl(cstep)%max_iter
           stepcnt=stepcnt+1                                      
           call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC%t_delta )
           if( fstrDYNAMIC%ray_k/=0.d0 .or. fstrDYNAMIC%ray_m/=0.d0 ) then     
             do j = 1 ,ndof*nnod                                   
                hecMAT%X(j) = fstrDYNAMIC%VEC2(j) - b3*fstrSOLID%dunode(j)    
             end do		 
           endif
           if( fstrDYNAMIC%ray_k/=0.d0 ) then    
             if( hecMESH%n_dof == 3 ) then
               call hecmw_matvec_33 (hecMESH, hecMAT, hecMAT%X, fstrDYNAMIC%VEC3)
             else if( hecMESH%n_dof == 2 ) then
               call hecmw_matvec_22 (hecMESH, hecMAT, hecMAT%X, fstrDYNAMIC%VEC3, nnod)
             else if( hecMESH%n_dof == 6 ) then
               call matvec(fstrDYNAMIC%VEC3, hecMAT%X, hecMAT, ndof, hecMAT%D, hecMAT%AU, hecMAT%AL)
             end if	   
           endif
!C
!C-- mechanical boundary condition
           call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)  
           do j=1, hecMESH%n_node*  hecMESH%n_dof 
             hecMAT%B(j)=hecMAT%B(j)- fstrSOLID%QFORCE(j) + myEIG%mass(j)*( fstrDYNAMIC%VEC1(j)-a3*fstrSOLID%dunode(j)   & 
		       + fstrDYNAMIC%ray_m* hecMAT%X(j) ) + fstrDYNAMIC%ray_k*fstrDYNAMIC%VEC3(j)  
           end do
           do j = 1 ,nn*hecMAT%NP             
             hecMAT%D(j)  = c1* hecMAT%D(j)
           end do
           do j = 1 ,nn*hecMAT%NPU
             hecMAT%AU(j) = c1* hecMAT%AU(j)
           end do
           do j = 1 ,nn*hecMAT%NPL
             hecMAT%AL(j) = c1*hecMAT%AL(j)
           end do                                            
           do j=1,nnod               
           do kk=1,ndof
             idm = nn*(j-1)+1 + (ndof+1)*(kk-1)
             imm = ndof*(j-1) + kk
             hecMAT%D(idm) = hecMAT%D(idm) + c2*myEIG%mass(imm)
           end do
           end do                   

           if( fstr_is_contact_active() ) then        
             call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID)    
             call fstr_AddContactStiffness(cstep,iter,hecMAT,fstrMAT,fstrSOLID)   
           endif
!              
!C ********************************************************************************
!C for couple analysis
          if( fstrPARAM%fg_couple == 1) then
              if( fstrDYNAMIC%i_step > 1 .or. &
                  (fstrDYNAMIC%i_step==1 .and. fstrPARAM%fg_couple_first==1 )) then
                    call fstr_rcap_get( fstrCPL )
                    call dynamic_mat_ass_couple( hecMESH, hecMAT, fstrSOLID, fstrCPL )
              endif
          endif
!C ********************************************************************************

!C-- geometrical boundary condition

          call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt)     
          call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt) 
          call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt) 

! ----- check convergence
          call hecmw_innerProduct_R(hecMESH,ndof,hecMAT%B,hecMAT%B,res)
          res = sqrt(res)/n_node_global
          if( hecMESH%my_rank==0 ) then
            if( mod(i,int(fstrDYNAMIC%nout/10)) == 0 )   &
            write(*,'(a,i3,a,2e15.7)') ' - Resiual(',iter,') =',res
            write(ISTA,'(''iter='',I5,''- Residual'',2E15.7)')iter,res
          endif

! ----- check convergence
          if( .not.fstr_is_contact_active() ) maxDLag= 0.0d0                    
          call hecmw_allreduce_R1(hecMESH, maxDlag, HECMW_MAX)
          if( res<fstrSOLID%step_ctrl(cstep)%converg .and. maxDLag<1.0d-5 .and. iter>1 ) exit   
          rf=1.0d0                                                   
          if( iter>1 .and. res>res1 )rf=0.5d0*rf                         
          res1=res                                                                          

          call solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,rf)

! ----- update the strain, stress, and internal force
          call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,fstrDYNAMIC%t_delta,1,fstrDYNAMIC%strainEnergy ) 
          do j=1,hecMESH%n_node*ndof
            fstrSOLID%dunode(j)  = fstrSOLID%dunode(j)+hecMAT%X(j)
          enddo
    
! ----- update the Lagrange multipliers                                         
          if( fstr_is_contact_active() ) then  
            maxDLag = 0.0d0                            
            do j=1,fstrMAT%num_lagrange
              fstrMAT%lagrange(j) = fstrMAT%lagrange(j) + hecMAT%X(hecMESH%n_node*ndof+j) 
              if(dabs(hecMAT%X(hecMESH%n_node*ndof+j))>maxDLag) maxDLag=dabs(hecMAT%X(hecMESH%n_node*ndof+j)) 
!              write(*,*)'Lagrange:', j,fstrMAT%lagrange(j),hecMAT%X(hecMESH%n_node*ndof+j)      
            enddo                                     
          endif 
        enddo                             
	  
! -----  not convergence
        if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
           if( hecMESH%my_rank==0) then
             write(ILOG,*) '### Fail to Converge  : at step=', i
             write(ISTA,*) '### Fail to Converge  : at step=', i
             write(   *,*) '     ### Fail to Converge  : at step=', i
           end if
           stop       
        end if  

        call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B ) 
              
        if( hecMAT%Iarray(99)==4 .and. .not. fstr_is_contact_active() ) then                      
          write(*,*) ' This type of direct solver is not yet available in such case ! '
          write(*,*) ' Please use intel MKL direct solver !'
          call  hecmw_abort(hecmw_comm_get_comm())
        endif

        is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)
        contact_changed_global=0
        if( fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH) ) then
          exit loopFORcontactAnalysis
        elseif( fstr_is_matrixStructure_changed(infoCTChange) ) then  
          call fstr_mat_con_contact(cstep,hecMAT,fstrSOLID,fstrMAT,infoCTChange)      
          contact_changed_global=1
        endif
        call hecmw_allreduce_I1(hecMESH,contact_changed_global,HECMW_MAX)
        if (contact_changed_global > 0) then
          call solve_LINEQ_contact_init(hecMESH,hecMAT,fstrMAT,is_mat_symmetric)
        endif

        if( count_step > max_iter_contact ) exit loopFORcontactAnalysis                            
        

      ENDDO loopFORcontactAnalysis                   
!
!C-- new displacement, velocity and accelaration
!C
      fstrDYNAMIC%kineticEnergy = 0.0d0                                               
      do j = 1 ,ndof*nnod
          fstrDYNAMIC%ACC (j,2) = -a1*fstrDYNAMIC%ACC(j,1) - a2*fstrDYNAMIC%VEL(j,1) + &
                               a3*fstrSOLID%dunode(j)
          fstrDYNAMIC%VEL (j,2) = -b1*fstrDYNAMIC%ACC(j,1) - b2*fstrDYNAMIC%VEL(j,1) + &
                               b3*fstrSOLID%dunode(j)
          fstrDYNAMIC%ACC (j,1) = fstrDYNAMIC%ACC (j,2)
          fstrDYNAMIC%VEL (j,1) = fstrDYNAMIC%VEL (j,2)

          fstrSOLID%unode(j)  = fstrSOLID%unode(j)+fstrSOLID%dunode(j)
          fstrDYNAMIC%DISP(j,2) = fstrSOLID%unode(j)
          
          fstrDYNAMIC%kineticEnergy = fstrDYNAMIC%kineticEnergy + &                        
                                      0.5d0*myEIG%mass(j)*fstrDYNAMIC%VEL(j,2)*fstrDYNAMIC%VEL(j,2) 
      enddo

!---  Restart info
      if( mod(i,fstrDYNAMIC%restart_nout) == 0 ) then
        call fstr_write_restart_dyna(i,0,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
                                     infoCTChange%contactNode_current)               
      end if
!
!C-- output new displacement, velocity and accelaration
      call dynamic_nloutput(i,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
!C
!C-- output result of monitoring node
!C
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, my_rank_monit_1)    

      call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYNAMIC%t_delta )

    enddo
!C
!C-- end of time step loop

    call cpu_time(time_2)
    if( hecMESH%my_rank == 0 ) then
        write(ISTA,'(a,f10.2,a)') '         solve (sec) :', time_2 - time_1, 's'
    end if

  end subroutine fstr_solve_dynamic_nlimplicit_contactSLag 
  
  


end module fstr_dynamic_nlimplicit
