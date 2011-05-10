!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.1                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by  X. YUAN(AdavanceSoft)                         !
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
!>  \version    0.00
!!
!======================================================================!
module m_fstr_NonLinearMethod

use m_static_lib

use m_fstr_StiffMatrix
use m_fstr_Update
use m_fstr_AddBC
use m_fstr_Residual
use m_fstr_Result
use m_fstr_Restart

implicit none

contains


!> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
!> method 
subroutine fstr_Newton( cstep, hecMESH, hecMAT, fstrSOLID,     &
                        restrt_step_num, sub_step, fstrPARAM   )
  use m_fstr
  use hecmw
  use m_fstr_ass_load

  integer, intent(in)           :: cstep     !< current loading step
  type (hecmwST_local_mesh)     :: hecMESH   !< hecmw mesh
  type (hecmwST_matrix    )     :: hecMAT    !< hecmw matrix
  type (fstr_solid       )      :: fstrSOLID !< fstr_solid
  integer, intent(in)           :: sub_step  !< substep number of current loading step
  type (fstr_param),intent(in)  :: fstrPARAM

  integer(kind=kint) :: ndof
  integer(kind=kint) :: i, iter, itemp, ttemp
  integer(kind=kint) :: al_step, stepcnt
  real(kind=kreal)   :: tt0,tt, res, res0, res1, maxv, relres, tincr
  integer(kind=kint) :: restrt_step_num
  logical            :: ctchange, convg

  hecMAT%NDOF = hecMESH%n_dof
  ndof = hecMAT%NDOF
  tincr = fstrSOLID%step_ctrl(cstep)%initdt
  if( fstrSOLID%step_ctrl(cstep)%solution == stepStatic ) tincr = 0.d0
  call cpu_time(tt0)

  
  
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
	
    call fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID )
    fstrSOLID%dunode(:)  = 0.d0
! ----- Inner Iteration, lagrange multiplier constant
    res1=0.d0
    relres = 1.d0
    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt=stepcnt+1
      call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr )

! ----- Set Boundary condition
      call fstr_AddBC(cstep, sub_step, hecMESH,hecMAT,fstrSOLID,stepcnt)
!----- SOLVE [Kt]{du}={R}
      if( sub_step == restrt_step_num .and. iter == 1 ) hecMAT%Iarray(98) = 1   !Assmebly complete
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

! ----- update the strain, stress, and internal force
      call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,tincr,iter )

!   ----- update the small displacement and the displacement for 1step 
!            \delta u^k => solver's solution
!            \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i=1,hecMESH%n_node*ndof
        fstrSOLID%dunode(i)  = fstrSOLID%dunode(i) + hecMAT%X(i)
      enddo

! ----- Set residual
      call fstr_Update_NDForce(cstep,hecMESH,hecMAT,fstrSOLID )
      res = fstr_get_residual( hecMAT%B, hecMESH )

! ----- Gather global residual
      call hecmw_allREDUCE_R1(hecMESH,res,hecmw_sum)
      res = sqrt(res)/hecMESH%n_node
      if( iter==1 ) res0=res
      if( res0==0.d0 ) then
        res0 =1.d0
      else
        relres = dabs(res1-res)/res0
      endif
      if( hecMESH%my_rank==0 ) then
        write(*,'(a,i3,a,2e15.7)') ' - Resiual(',iter,') =',res,relres
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
      write(ISTA,'("### Converged in NR ietration : CPU time=",E10.4,"   iter=",I6)') tt-tt0,iter
    endif
	
    if( fstrSOLID%TEMP_irres>0 ) then
          if(fstrSOLID%restart_nout<0) then
            fstrSOLID%restart_nout = - fstrSOLID%restart_nout
          end if
          if( mod(itemp,fstrSOLID%restart_nout) == 0 ) then
            call fstr_write_restart(cstep,1,itemp,hecMESH,fstrSOLID)
          end if  

! ----- Result output (include visualize output)
          call fstr_OutputResult( cstep, itemp, hecMESH, hecMAT, fstrSOLID, fstrPARAM, tt )
    endif
  
    if( fstrSOLID%TEMP_irres>1 ) fstrSOLID%TEMP_tstep = fstrSOLID%TEMP_tstep+1
  enddo

end subroutine fstr_Newton



end module m_fstr_NonLinearMethod
