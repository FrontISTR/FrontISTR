!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by  X. YUAN(AdavanceSoft)                         !
!                        Z. Sun(ASTOM)                                 !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions on nonlinear analysis
!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \date       2010/02/01
!>  \author     Z. Sun(ASTOM)
!>  \author     2010/11
!>  \version    0.00
!!
!======================================================================!
module m_fstr_NonLinearMethod

use m_fstr
use m_fstr_para_contact
use m_static_lib
use m_static_output

use m_fstr_StiffMatrix
use m_fstr_Update
use m_fstr_ass_load
use m_fstr_AddBC
use m_fstr_Residual
use m_fstr_Restart
use fstr_matrix_con_contact

implicit none

contains


!> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
!> method
subroutine fstr_Newton( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM,      &
                         restrt_step_num, sub_step  )

  integer, intent(in)                   :: cstep     !< current loading step
  type (hecmwST_local_mesh)              :: hecMESH   !< hecmw mesh
  type (hecmwST_matrix    )              :: hecMAT    !< hecmw matrix
  type (fstr_solid        )              :: fstrSOLID !< fstr_solid
  integer, intent(in)                   :: sub_step  !< substep number of current loading step
  type (fstr_param)                      :: fstrPARAM !< type fstr_param
  type (fstrST_matrix_contact_lagrange)  :: fstrMAT   !< type fstrST_matrix_contact_lagrange

  integer(kind=kint) :: ndof
  integer(kind=kint) :: i, iter, itemp, ttemp, tintl
  integer(kind=kint) :: al_step, stepcnt
  real(kind=kreal)   :: tt0,tt, res, res0, res1, maxv, relres, tincr
  integer(kind=kint) :: restrt_step_num
  logical            :: convg
  integer(kind=kint) :: n_node_global

! sum of n_node among all subdomains (to be used to calc res)
  n_node_global = hecMESH%nn_internal
  call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

  hecMAT%NDOF = hecMESH%n_dof
  ndof = hecMAT%NDOF
  tincr = fstrSOLID%step_ctrl(cstep)%initdt
  if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0
  call cpu_time(tt0)

  hecMAT%X = 0.d0

  stepcnt = 0
  ttemp = 1
  tintl = 1
  if( fstrSOLID%TEMP_irres>1 ) then
    ttemp = fstrSOLID%TEMP_irres
    tintl = fstrSOLID%TEMP_interval
  endif
  do itemp = fstrSOLID%TEMP_tstep, ttemp, tintl
    if( fstrSOLID%TEMP_irres>0 ) then
      if( hecMESH%my_rank==0 ) then
        write(*,*) " - Read in temperature in time step", fstrSOLID%TEMP_tstep
        write(ISTA,*) " - Read in temperature in time step", fstrSOLID%TEMP_tstep
      endif
    endif

    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM )
    fstrSOLID%dunode(:)  = 0.d0
! ----- Inner Iteration, lagrange multiplier constant
    res1=0.d0
    relres = 1.d0
    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt=stepcnt+1
      call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr )
      call fstr_AddSPRING(cstep, sub_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM) 

! ----- Set Boundary condition
      call fstr_AddBC(cstep, sub_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT, stepcnt)
!----- SOLVE [Kt]{du}={R}
      if( sub_step == restrt_step_num .and. iter == 1 ) hecMAT%Iarray(98) = 1
      if( iter >= 1 ) then
        hecMAT%Iarray(97) = 1   !Need numerical factorization
      end if
      CALL solve_LINEQ(hecMESH,hecMAT,imsg)

          if( hecMESH%n_dof==3 ) then
             call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
          else if( hecMESH%n_dof==2 ) then
             call hecmw_update_2_R (hecMESH, hecMAT%X, hecMAT%NP)
          endif
		  
!   ----- update the small displacement and the displacement for 1step 
!            \delta u^k => solver's solution
!            \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i=1,hecMESH%n_node*ndof
        fstrSOLID%dunode(i)  = fstrSOLID%dunode(i) + hecMAT%X(i)
      enddo

! ----- update the strain, stress, and internal force
      call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,tincr,iter )  

! ----- Set residual
      if (fstrSOLID%DLOAD_follow /= 0) &
           call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM )
      call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,sub_step ) 
      res = fstr_get_residual( hecMAT%B, hecMESH )

! ----- Gather global residual
      res = sqrt(res)/n_node_global
      if( iter==1 ) res0=res
      if( res0==0.d0 ) then
        res0 =1.d0
      else
        relres = dabs(res1-res)/res0
      endif
      if( hecMESH%my_rank==0 ) then
        write(*,'(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
        write(ISTA,'(''iter='',I5,''res/res0='',2E15.7)')iter,res,relres
      endif

! ----- check convergence
      if( res<fstrSOLID%step_ctrl(cstep)%converg  .or.    &
          relres<fstrSOLID%step_ctrl(cstep)%converg ) exit
      res1 = res

    enddo
! ----- end of inner loop

! -----  not convergence
    if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
       if( hecMESH%my_rank==0) then
         write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
         write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
         write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
       end if
       stop
    end if

!   ----- update the total displacement
!             u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i)  = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo

    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )

    call cpu_time(tt)
    if( hecMESH%my_rank==0) then
      write(ISTA,'("### Converged in NR ietration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter
    endif

    if( fstrSOLID%TEMP_irres>0 ) then
          if(fstrSOLID%restart_nout<0) then
            fstrSOLID%restart_nout = - fstrSOLID%restart_nout
          end if
          if( mod(itemp,fstrSOLID%restart_nout) == 0 .or. (fstrSOLID%TEMP_irres>1 .and. itemp==ttemp) ) then
            call fstr_write_restart(cstep,1,itemp,hecMESH,fstrSOLID,fstrPARAM) 
          end if  

! ----- Result output (include visualize output)
          call fstr_static_Output( cstep, itemp, hecMESH, fstrSOLID, fstrPR%solution_type )
    endif

    if( fstrSOLID%TEMP_irres>1 ) fstrSOLID%TEMP_tstep = fstrSOLID%TEMP_tstep + tintl
  enddo

end subroutine fstr_Newton

!> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
!> method combined with Nested iteration of augmentation calculation as suggested
!> by Simo & Laursen (Compu & Struct, Vol42, pp97-116, 1992 )
subroutine fstr_Newton_contactALag( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM,     &      
                                restart_step_num, restart_substep_num, sub_step, infoCTChange  )  
  use mContact                                

  integer, intent(in)                   :: cstep     !< current loading step
  type (hecmwST_local_mesh)              :: hecMESH   !< hecmw mesh
  type (hecmwST_matrix    )              :: hecMAT    !< hecmw matrix
  type (fstr_solid        )              :: fstrSOLID !< fstr_solid
  integer, intent(in)                   :: sub_step  !< substep number of current loading step
  type (fstr_param)                      :: fstrPARAM !< type fstr_param
  type (fstr_info_contactChange)         :: infoCTChange  !< fstr_info_contactChange         
  type (fstrST_matrix_contact_lagrange)  :: fstrMAT !< type fstrST_matrix_contact_lagrange     

  integer(kind=kint) :: ndof
  integer(kind=kint) :: ctAlgo                                                           
  integer(kind=kint) :: i, iter, itemp, ttemp
  integer(kind=kint) :: al_step, n_al_step, stepcnt
  real(kind=kreal)   :: tt0,tt, res, res0, res1, maxv, relres, tincr
  integer(kind=kint) :: restart_step_num,  restart_substep_num
  logical            :: convg, ctchange
  integer(kind=kint) :: n_node_global

! sum of n_node among all subdomains (to be used to calc res)
  n_node_global = hecMESH%nn_internal
  call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

  ctAlgo = fstrPARAM%contact_algo                     

  hecMAT%NDOF = hecMESH%n_dof
  ndof = hecMAT%NDOF
  tincr = fstrSOLID%step_ctrl(cstep)%initdt                                  
  if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0       
  call cpu_time(tt0)

  if( cstep == restart_step_num .and. sub_step == restart_substep_num ) then       
    if(hecMESH%my_rank==0) write(*,*) "---Scanning initial contact state---"
    call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange )
    if(hecMESH%my_rank==0) write(*,*)
  endif

  hecMAT%X = 0.d0

  stepcnt = 0                                                               
  ttemp = 1
  if( fstrSOLID%TEMP_irres>1 ) ttemp = fstrSOLID%TEMP_irres
  do itemp = fstrSOLID%TEMP_tstep, ttemp                                   
    if( fstrSOLID%TEMP_irres>0 ) then
      if( hecMESH%my_rank==0 ) then
        write(*,*) " - Read in temperature in time step", fstrSOLID%TEMP_tstep
        write(ISTA,*) " - Read in temperature in time step", fstrSOLID%TEMP_tstep
      endif
    endif                                                                  
	
    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM )   
    
    ! ----- Augmentation loop. In case of no contact, it is inactive
    n_al_step = 10
    if( .not. fstr_is_contact_active() ) n_al_step = 1
    do al_step=1,n_al_step
      if( hecMESH%my_rank==0) then
         if( n_al_step>1 ) then
           print *, "Augmentation step:",al_step
           write(IMSG,*) "Augmentation step:",al_step
         endif
      end if
    
      fstrSOLID%dunode(:)  = 0.d0
      
! ----- Inner Iteration, lagrange multiplier constant
      res1=0.d0
      relres = 1.d0
      do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
        stepcnt=stepcnt+1
        call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr )           
        call fstr_AddSPRING(cstep, sub_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM)
! ----- Contact
        if( fstr_is_contact_active() .and. stepcnt==1 )  then
          maxv = hecmw_mat_diag_max( hecMAT, hecMESH )
          call fstr_set_contact_penalty( maxv )
        endif
        if( fstr_is_contact_active() ) then
          call fstr_contactBC( iter, hecMESH, hecMAT, fstrSOLID )
        endif

! ----- Set Boundary condition
        call fstr_AddBC(cstep, sub_step, hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrMAT,stepcnt)  
!----- SOLVE [Kt]{du}={R}
        if( sub_step == restart_step_num .and. iter == 1 ) hecMAT%Iarray(98) = 1   
        if( iter >= 1 ) then
          hecMAT%Iarray(97) = 1   !Need numerical factorization
        end if
        CALL solve_LINEQ(hecMESH,hecMAT,imsg)
	  
          if( hecMESH%n_dof==3 ) then
             call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
             if( hecMESH%my_rank==0) then
                write(IMSG,*) 'hecmw_update_3_R: OK'
             end if
          else if( hecMESH%n_dof==2 ) then
             call hecmw_update_2_R (hecMESH, hecMAT%X, hecMAT%NP)
             if( hecMESH%my_rank==0) then
                write(IMSG,*) 'hecmw_update_2_R: OK'
             end if
          endif

!   ----- update the small displacement and the displacement for 1step 
!            \delta u^k => solver's solution
!            \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i=1,hecMESH%n_node*ndof
        fstrSOLID%dunode(i)  = fstrSOLID%dunode(i) + hecMAT%X(i)
      enddo

! ----- update the strain, stress, and internal force
      call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,tincr,iter )

! ----- Set residual
        if (fstrSOLID%DLOAD_follow /= 0) &
             call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM )
        call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,sub_step )       
        if( fstr_is_contact_active() ) then                                     
          call fstr_update_contact0( hecMESH, fstrSOLID, hecMAT%B ) 
        endif
        res = fstr_get_residual( hecMAT%B, hecMESH )

! ----- Gather global residual
        res = sqrt(res)/n_node_global
        if( iter==1 ) res0=res
        if( res0==0.d0 ) then
          res0 =1.d0
        else
          relres = dabs(res1-res)/res0
        endif
        if( hecMESH%my_rank==0 ) then
          write(*,'(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
          write(ISTA,'(''iter='',I5,''res/res0='',2E15.7)')iter,res,relres
        endif

! ----- check convergence
        if( res<fstrSOLID%step_ctrl(cstep)%converg  .or.    &
            relres<fstrSOLID%step_ctrl(cstep)%converg ) exit
        res1 = res

      enddo  
! ----- end of inner loop

! -----  not convergence
      if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
         if( hecMESH%my_rank==0) then
           write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
           write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
           write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
         end if
         stop
      end if
      
! ----- deal with contact boundary
      convg = .true.
      ctchange = .false.
      if( associated(fstrSOLID%contacts) ) then
        call fstr_update_contact_multiplier( hecMESH, fstrSOLID, ctchange )
        call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )   
        if( infoCTChange%contact2free+infoCTChange%contact2neighbor+infoCTChange%free2contact > 0 ) & 
        ctchange = .true.
      endif
      if( fstr_is_contact_active() ) then    
        gnt(2)=gnt(2)/iter
        convg = fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH)
      endif      
	
!   ----- update the total displacement
!             u_{n+1} = u_{n} + \Delta u_{n+1}
      do i=1,hecMESH%n_node*ndof
        fstrSOLID%unode(i)  = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
      enddo
      if( convg .and. (.not.ctchange) ) exit  
    enddo
! ----- end of augmentation loop    
	
    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )     

    call cpu_time(tt)
    if( hecMESH%my_rank==0) then
      write(ISTA,'("### Converged in NR ietration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter
    endif
	
    if( fstrSOLID%TEMP_irres>0 ) then                                  
          if(fstrSOLID%restart_nout<0) then
            fstrSOLID%restart_nout = - fstrSOLID%restart_nout
          end if
          if( mod(itemp,fstrSOLID%restart_nout) == 0 ) then
            call fstr_write_restart(cstep,1,itemp,hecMESH,fstrSOLID,fstrPARAM,infoCTChange%contactNode_current) 
          end if  

! ----- Result output (include visualize output)
          call fstr_static_Output( cstep, itemp, hecMESH, fstrSOLID, fstrPR%solution_type )
    endif                                                             
  
    if( fstrSOLID%TEMP_irres>1 ) fstrSOLID%TEMP_tstep = fstrSOLID%TEMP_tstep+1  
  enddo                                                                        

end subroutine fstr_Newton_contactALag


!> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson method.        
!> Standard Lagrange multiplier algorithm for contact analysis is incoluded in this subroutine. 
subroutine fstr_Newton_contactSLag( cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT,     &      
                                 restart_step_num, restart_substep_num, sub_step, infoCTChange, conMAT )    

  use mContact                                                                            
  use m_addContactStiffness  
  use m_solve_LINEQ_contact                                                                                                                            

  integer, intent(in)                   :: cstep        !< current loading step
  type (hecmwST_local_mesh)              :: hecMESH      !< hecmw mesh
  type (hecmwST_matrix    )              :: hecMAT       !< hecmw matrix
  type (fstr_solid        )              :: fstrSOLID    !< fstr_solid
  integer, intent(in)                   :: sub_step     !< substep number of current loading step
  type (fstr_param)                      :: fstrPARAM    !< type fstr_param
  type (fstr_info_contactChange)         :: infoCTChange !< fstr_info_contactChange               
  type (fstrST_matrix_contact_lagrange)  :: fstrMAT      !< type fstrST_matrix_contact_lagrange 
  type (hecmwST_matrix    ),optional     :: conMAT

  integer(kind=kint) :: ndof
  integer(kind=kint) :: ctAlgo                                                           
  integer(kind=kint) :: i, iter, itemp, ttemp, max_iter_contact                          
  integer(kind=kint) :: stepcnt, count_step                                             
  real(kind=kreal)    :: tt0,tt, res, res0, res1, maxv, relres, tincr,resX,ierr
  integer(kind=kint) :: restart_step_num,  restart_substep_num   
  logical :: is_mat_symmetric
  integer(kind=kint) :: n_node_global
  integer(kind=kint) :: contact_changed_global
  integer(kint)   ::  nndof,npdof
  real(kreal),allocatable :: tmp_conB(:)
  integer   ::  istat,i0
  real(kreal) ::  q_residual,x_residual

! sum of n_node among all subdomains (to be used to calc res)
  n_node_global = hecMESH%nn_internal
  call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)
 
  if( hecMAT%Iarray(99)==4 .and. .not.fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH) ) then
    write(*,*) ' This type of direct solver is not yet available in such case ! '
    write(*,*) ' Please use intel MKL direct solver !'
    call  hecmw_abort(hecmw_comm_get_comm())
  endif        
  ctAlgo = fstrPARAM%contact_algo                                                        

  hecMAT%NDOF = hecMESH%n_dof
  ndof = hecMAT%NDOF
  tincr = fstrSOLID%step_ctrl(cstep)%initdt                                  
  if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0       
  call cpu_time(tt0)

  fstrSOLID%dunode(:)  = 0.d0                                                       
                                            
  if( cstep==restart_step_num.and.sub_step==restart_substep_num  ) then       
    call fstr_save_originalMatrixStructure(hecMAT)   
    call fstr_scan_contact_state( cstep, ctAlgo, hecMESH, fstrSOLID, infoCTChange, hecMAT%B )   
    if(paraContactFlag.and.present(conMAT)) then
      call copyClearMatrix(hecMAT,conMAT)
    endif
    if ( fstr_is_contact_active() ) then                                                 
!     ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange, conMAT) 
      else
        call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange)    
      endif
    elseif( hecMAT%Iarray(99)==4 ) then                                                   
      write(*,*) ' This type of direct solver is not yet available in such case ! '
      write(*,*) ' Please change the solver type to intel MKL direct solver !'
      call  hecmw_abort(hecmw_comm_get_comm()) 
    endif  
    is_mat_symmetric = fstr_is_matrixStruct_symmetric(fstrSOLID,hecMESH)
    call solve_LINEQ_contact_init(hecMESH,hecMAT,fstrMAT,is_mat_symmetric)
  endif                                                                                
  
  stepcnt = 0                                                               
  ttemp = 1
  if( fstrSOLID%TEMP_irres>1 ) ttemp = fstrSOLID%TEMP_irres
  do itemp = fstrSOLID%TEMP_tstep, ttemp                                   
    if( fstrSOLID%TEMP_irres>0 ) then
      if( hecMESH%my_rank==0 ) then
        write(*,*) " - Read in temperature in time step", fstrSOLID%TEMP_tstep
        write(ISTA,*) " - Read in temperature in time step", fstrSOLID%TEMP_tstep
      endif
    endif                                                                  
	  
    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM )   
      
    if(paraContactFlag.and.present(conMAT)) then
      conMAT%B(:) = 0.0D0
    endif
    if( fstr_is_contact_active() )  then
!    ----  For Parallel Contact with Multi-Partition Domains
      if(paraContactFlag.and.present(conMAT)) then
        call fstr_ass_load_contact(cstep, hecMESH, conMAT, fstrSOLID, fstrMAT)
      else
        call fstr_ass_load_contact(cstep, hecMESH, hecMAT, fstrSOLID, fstrMAT)
      endif
      
    endif 
  
    fstrSOLID%dunode(:)  = 0.d0
    
    
    max_iter_contact = 10                       
    count_step = 0
    loopFORcontactAnalysis: DO WHILE( .TRUE. )                                        
 
      count_step = count_step + 1  
! ----- Inner Iteration
      res1=0.d0
      relres = 1.d0
      do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      call hecmw_BARRIER(hecMESH)
      if(myrank == 0)print *,'-------------------------------------------------'
      call hecmw_BARRIER(hecMESH)
        stepcnt=stepcnt+1
        call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr )           
        call fstr_AddSPRING(cstep, sub_step, hecMESH, hecMAT, fstrSOLID, fstrPARAM)
        
        if(paraContactFlag.and.present(conMAT)) then
          call hecMAT_clear( conMAT )
        endif
        if( fstr_is_contact_active() )  then
!    ----  For Parallel Contact with Multi-Partition Domains
          if(paraContactFlag.and.present(conMAT)) then
            call fstr_AddContactStiffness(cstep,iter,conMAT,fstrMAT,fstrSOLID)
          else
            call fstr_AddContactStiffness(cstep,iter,hecMAT,fstrMAT,fstrSOLID)
          endif
        endif

! ----- Set Boundary condition

        if(paraContactFlag.and.present(conMAT)) then
          call fstr_AddBC(cstep, sub_step, hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrMAT,stepcnt,conMAT)
        else
          call fstr_AddBC(cstep, sub_step, hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrMAT,stepcnt)
        endif
        
        nndof = hecMAT%N*hecMAT%ndof
        
!----- SOLVE [Kt]{du}={R}
!   ----  For Parallel Contact with Multi-Partition Domains
        if(paraContactFlag.and.present(conMAT)) then
          q_residual = fstr_get_norm_para_contact(hecMAT,fstrMAT,conMAT,hecMESH)
          call solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,1.0D0,conMAT)
        else
          q_residual = fstr_get_norm_contact('residualForce',hecMESH,hecMAT,fstrSOLID,fstrMAT)
          call solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT)
        endif
        x_residual = fstr_get_x_norm_contact(hecMAT,fstrMAT,hecMESH)
! ----- update external nodal displacement increments
        if(paraContactFlag.and.present(conMAT)) then
          call paraContact_update_3_R(hecMESH,hecMAT%X)
        else
          call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
        endif
        call hecmw_innerProduct_R(hecMESH,ndof,hecMAT%X,hecMAT%X,resX)
        resX = sqrt(resX)/n_node_global
        if( hecMESH%my_rank==0 ) then
          write(*,'(a,i3,a,e15.7)') ' - ResidualX    (',iter,') =',resX
          write(*,'(a,i3,a,e15.7)') ' - ResidualX+LAG(',iter,') =',sqrt(x_residual)/n_node_global
          write(*,'(a,i3,a,e15.7)') ' - ResidualQ    (',iter,') =',sqrt(q_residual)/n_node_global
        endif
!   ----- update the small displacement and the displacement for 1step 
        do i=1,hecMESH%n_node*ndof
          fstrSOLID%dunode(i)  = fstrSOLID%dunode(i) + hecMAT%X(i)
        enddo
 
 ! ----- update the Lagrange multipliers   
        if( fstr_is_contact_active() ) then                             
          do i=1,fstrMAT%num_lagrange                 
            fstrMAT%lagrange(i) = fstrMAT%lagrange(i) + hecMAT%X(hecMESH%n_node*ndof+i)  
          enddo                                     
        endif         
! ----- update the strain, stress, and internal force (only QFORCE)
        call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,tincr,iter )

! ----- Set residual
        if(paraContactFlag.and.present(conMAT)) then
          call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,sub_step,conMAT )
        else
          call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID,sub_step )      
        endif 
        
        if( fstr_is_contact_active() )  then
          if(paraContactFlag.and.present(conMAT)) then
            conMAT%B(1:3*hecMAT%NP) = 0.0D0
            call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID,conMAT)
          else
            call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID)
          endif
        endif 
                    
        if(paraContactFlag.and.present(conMAT)) then
          res = fstr_get_norm_para_contact(hecMAT,fstrMAT,conMAT,hecMESH) 
        else
          res = fstr_get_norm_contact('residualForce',hecMESH,hecMAT,fstrSOLID,fstrMAT)
        endif

        res = sqrt(res)/n_node_global
        if( iter==1 ) res0=res
        if( res0==0.d0 ) then
          res0 =1.d0
        else
          relres = dabs(res1-res)/res0
        endif
        if( hecMESH%my_rank==0 ) then
          write(*,'(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
          write(ISTA,'(''iter='',I5,''res/res0='',2E15.7)')iter,res,relres
        endif

! ----- check convergence
        if( res<fstrSOLID%step_ctrl(cstep)%converg  .or.    &
            relres<fstrSOLID%step_ctrl(cstep)%converg ) exit
        res1 = res

      enddo  
! ----- end of inner loop

! -----  not convergence
      if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
         if( hecMESH%my_rank==0) then
           write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
           write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
           write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
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
!     ----  For Parallel Contact with Multi-Partition Domains
        
        if(paraContactFlag.and.present(conMAT)) then
          call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange, conMAT) 
        else
          call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange)    
        endif
        contact_changed_global=1
      else

      endif
      call hecmw_allreduce_I1(hecMESH,contact_changed_global,HECMW_MAX)
      if (contact_changed_global > 0) then
        hecMAT%B(:) = 0.0D0
        if(paraContactFlag.and.present(conMAT)) then
          conMAT%B(:) = 0.0D0
        endif
        call solve_LINEQ_contact_init(hecMESH,hecMAT,fstrMAT,is_mat_symmetric)
      endif 
      if( count_step > max_iter_contact ) exit loopFORcontactAnalysis                   

    ENDDO loopFORcontactAnalysis                                             
    
	
!   ----- update the total displacement
!             u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,hecMESH%n_node*ndof
      fstrSOLID%unode(i)  = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo
	
    call fstr_UpdateState( hecMESH, fstrSOLID, tincr )     

    call cpu_time(tt)
    if( hecMESH%my_rank==0) then
      write(ISTA,'("### Converged in contact iteration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter  
    endif
	
    if( fstrSOLID%TEMP_irres>0 ) then                                  
          if(fstrSOLID%restart_nout<0) then
            fstrSOLID%restart_nout = - fstrSOLID%restart_nout
          end if
          if( mod(itemp,fstrSOLID%restart_nout) == 0 ) then
            call fstr_write_restart(cstep,1,itemp,hecMESH,fstrSOLID,fstrPARAM,infoCTChange%contactNode_current) 
          end if  

! ----- Result output (include visualize output)
          call fstr_static_Output( cstep, itemp, hecMESH, fstrSOLID, fstrPR%solution_type )
    endif                                                             
  
    if( fstrSOLID%TEMP_irres>1 ) fstrSOLID%TEMP_tstep = fstrSOLID%TEMP_tstep+1 
  enddo                                                                        

end subroutine fstr_Newton_contactSLag



end module m_fstr_NonLinearMethod
