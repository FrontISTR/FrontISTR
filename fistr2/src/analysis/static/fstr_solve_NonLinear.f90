!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2010/01/16
!>  \version    0.00
!!
!======================================================================!
module m_fstr_NonLinearMethod

use m_static_lib

use m_fstr_StiffMatrix
use m_fstr_Update
use m_fstr_AddBC
use m_fstr_Residual

implicit none

contains


!> \brief This subroutine solve nonlinear solid mechanics problems by Newton-Raphson
!> method 
subroutine fstr_Newton( cstep, fstrSOLID, restrt_step_num, sub_step, factor   )
  use m_fstr
  use hecmw
  use m_fstr_ass_load
  
  include "HEC_MW3_For.h"

  integer, intent(in)           :: cstep     !< current loading step
  type (fstr_solid       )      :: fstrSOLID !< fstr_solid
  integer, intent(in)           :: sub_step  !< substep number of current loading step
  real(kind=kreal), intent(in)  :: factor    !< loading factor of current substep

  integer(kind=kint) :: ndof
  integer(kind=kint) :: i, iter, iErr
  integer(kind=kint) :: al_step, n_al_step, stepcnt
  real(kind=kreal)   :: tt0,tt, res, res0, res1, maxv, relres, tincr
  integer(kind=kint) :: restrt_step_num
  integer :: iAss, iPart, iElem, iNode, iGrp, ndID, ik, snode, enode
  logical            :: ctchange, convg

  ndof = assDOF(1)
  tincr = fstrSOLID%step_ctrl(cstep)%initdt
  call cpu_time(tt0)
  
  do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         call mw_select_algebra(0)
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            call mw_vector_clear( iPart )
         enddo
  enddo

  call fstr_UpdateEPState( fstrSOLID )
! ----- Set Load Vector
  call fstr_ass_load(cstep, fstrSOLID,factor)
  
  stepcnt = 0
    
    fstrSOLID%dunode(:)  = 0.d0
    res1=0.d0
    relres = 1.d0
    do iter=1,fstrSOLID%step_ctrl(cstep)%max_iter
      stepcnt=stepcnt+1
      call fstr_StiffMatrix( fstrSOLID, cstep, tincr, factor )

! ----- Set Boundary condition
      call fstr_AddBC(cstep,fstrSOLID,stepcnt)
!----- SOLVE [Kt]{du}={R}
      iErr=mw_solve(svIarray(1), svRarray(1), svIarray(2), svIarray(3))
      if( iErr==0 ) then
        write(*,*) ' Fatal error: Fails in solver! '
        write(IMSG,*) ' Fatal error: Fails in solver! '
        call hecmw_abort(hecmw_comm_get_comm())
      endif

     ! call mw_send_recv()
      write(IMSG,*) 'hecmw_update: OK'
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part_with_id( iPart )
            snode = part_nodes(iAss+1,iPart+1)
            enode = part_nodes(iAss+1,iPart+2)
            call mw_get_solution_vector(fstrSOLID%ddunode(snode*ndof+1:  &
               enode*ndof), iPart)
         enddo
      enddo

! ----- update the strain, stress, and internal force
      call fstr_UpdateNewton( fstrSOLID,sub_step,tincr,factor,iter )

!   ----- update the small displacement and the displacement for 1step 
!            \delta u^k => solver's solution
!            \Delta u_{n+1}^{k} = \Delta u_{n+1}^{k-1} + \delta u^k
      do i=1,total_node*ndof
        fstrSOLID%dunode(i)  = fstrSOLID%dunode(i) + fstrSOLID%ddunode(i)
      enddo

! ----- Set residual
      call fstr_Update_NDForce(cstep,fstrSOLID,factor )

      res = fstr_get_residual( fstrSOLID )

! ----- Gather global residual
      call mw_allreduce_r( res, 1, mw_mpi_sum() )

      res = sqrt(res)/total_node
      if( iter==1 ) res0=res
      if( res0==0.d0 ) then
        res0 =1.d0
      else
        relres = dabs(res1-res)/res0
      endif
      if( myrank==0 ) then
        write(*,'(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
        write(ISTA,'(''iter='',I5,''res/res0='',2E15.7)')iter,res,relres
      endif

! ----- check convergence
      if( res<fstrSOLID%step_ctrl(cstep)%converg  .or.    &
          relres<fstrSOLID%step_ctrl(cstep)%converg ) exit
      res1 = res

    enddo  

! -----  not convergence
    if( iter>fstrSOLID%step_ctrl(cstep)%max_iter ) then
       if( myrank==0) then
         write(ILOG,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
         write(ISTA,'(a,i5,a,i5)') '### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
         write(   *,'(a,i5,a,i5)') '     ### Fail to Converge  : at total_step=', cstep, '  sub_step=', sub_step
       end if
       stop
    end if
	
!   ----- update the total displacement
!             u_{n+1} = u_{n} + \Delta u_{n+1}
    do i=1,total_node*ndof
      fstrSOLID%unode(i)  = fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo
	

  call cpu_time(tt)
  if( myrank==0) then
   write(ISTA,'("### Converged in NR ietration : CPU time=",E11.4,"   iter=",I6)') tt-tt0,iter
  endif

end subroutine fstr_Newton



end module m_fstr_NonLinearMethod
