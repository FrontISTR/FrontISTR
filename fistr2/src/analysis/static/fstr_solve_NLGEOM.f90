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
!> \brief This module provides main suboruitne for nonliear calculation.
!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2011/02/01
!>  \version    0.00
!
!======================================================================!
module m_fstr_solve_NLGEOM

use m_fstr
use hecmw
use m_static_lib
use m_static_post

use m_fstr_NonLinearMethod
use m_fstr_Result
use m_fstr_Restart

    implicit none

contains

!======================================================================!
!> \brief This module provides main suborutine for nonlinear calculation.
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2010/01/31
!>  \version    0.00

      subroutine FSTR_SOLVE_NLGEOM(fstrSOLID)
!======================================================================!
      type (fstr_solid       )  :: fstrSOLID   !< we need boundary conditions of curr step

      integer(kind=kint) :: ndof

      integer(kind=kint) :: j, tot_step, step_count
      integer(kind=kint) :: sub_step
      real(kind=kreal) :: tt, factor
      real(kind=kreal) :: time_1, time_2
      logical          :: ctchange
      integer(kind=kint) :: restart_step_num, restart_substep_num

      ndof = assDOF(1)

      if( fstrSOLID%restart_nout == 0 ) then
        fstrSOLID%restart_nout = 999999999
      end if

      restart_step_num = 1
      restart_substep_num = 1
      fstrSOLID%unode = 0.0
      if(fstrSOLID%restart_nout <0 ) then
        call fstr_read_restart(restart_step_num,restart_substep_num,fstrSOLID)
      endif

      fstrSOLID%FACTOR    =0.0

      step_count = 0
      do tot_step=restart_step_num, fstrSOLID%nstep_tot
        if(myrank==0) write(*,*) ''
        if(myrank==0) write(*,'(a,i5)') ' loading step=',tot_step
        
! -------------------------------------------------------------------------
!      STEP LOOP
! -------------------------------------------------------------------------
        do sub_step = restart_substep_num, fstrSOLID%step_ctrl(tot_step)%num_substep

! ----- time history of factor
          tt = sub_step*fstrSOLID%step_ctrl(tot_step)%initdt
          call table_nlsta(fstrSOLID,tot_step,sub_step-1,factor)
            fstrSOLID%FACTOR(1) = factor
          call table_nlsta(fstrSOLID,tot_step,sub_step, factor)
            fstrSOLID%FACTOR(2) = factor

          if( dabs(fstrSOLID%FACTOR(2)-fstrSOLID%FACTOR(1))<1.0e-20 ) then
            if( myrank==0) then
               write(imsg,*) 'loading increment finished due to f_2-f_1 = 0.0'
            end if
            exit
          end if
		  
          if( fstrSOLID%step_ctrl(tot_step)%solution==2 ) then
		    if( sub_step==1 ) then
             fstrSOLID%FACTOR(1)=0.d0; fstrSOLID%FACTOR(2)=1.d0
            else
			 fstrSOLID%FACTOR(:)=0.d0
            endif
          endif
		  
          if(myrank==0) write(*,*) ' substep=',sub_step,fstrSOLID%FACTOR
          step_count = step_count+1
          call cpu_time(time_1)

!       analysis algorithm ( Newton-Rapshon Method )
          call fstr_Newton( tot_step, fstrSOLID, restart_step_num, sub_step, factor   )

          if(fstrSOLID%restart_nout<0) then
            fstrSOLID%restart_nout = - fstrSOLID%restart_nout
          end if
          if( mod(step_count,fstrSOLID%restart_nout) == 0 ) then
            call fstr_write_restart(tot_step,sub_step,fstrSOLID)
          end if  

! ----- Result output (include visualize output)
          call fstr_OutputResult( tot_step, step_count, fstrSOLID, tt )

          call cpu_time(time_2)
          if( myrank==0) then
            write(ISTA,'(a,f10.2)') '         solve (sec) :', time_2 - time_1
          end if

        enddo    !--- end of substep  loop 
		  
      enddo      !--- end of tot_step loop
!
!  message
! 
      IF(myrank == 0)THEN
        WRITE(IMSG,'("### FSTR_SOLVE_NLGEOM FINISHED!")')
      ENDIF

      end subroutine FSTR_SOLVE_NLGEOM
	  

!C================================================================C
!> \brief This subroutine decide the loading increment considering
!>        the amplitude definition
!C================================================================C
    subroutine table_nlsta(fstrSOLID, cstep, substep, f_t)
    type ( fstr_solid         ), intent(in) :: fstrSOLID  !< fstr_solid
    integer(kind=kint), intent(in)          :: cstep      !< curr loading step
    integer(kind=kint), intent(in)          :: substep    !< curr substep number
    real(kind=kreal), intent(out)           :: f_t        !< loading factor

    integer(kind=kint) :: i, nitem
    integer(kind=kint) :: jj_n_amp, jj1, jj2
    integer(kind=kint) :: s1, s2, flag
    real(kind=kreal) :: t_1, t_2, t_t, f_1, f_2, tincre


      jj_n_amp = fstrSOLID%step_ctrl( cstep )%amp_id
      tincre = fstrSOLID%step_ctrl( cstep )%initdt

      if( jj_n_amp <= 0 ) then  ! Amplitude not defined
          f_t = tincre*substep

      else
          nitem = size( MWAmplitudes( jj_n_amp )%table, 1 )
          t_t = tincre*substep

          if(t_t > MWAmplitudes( jj_n_amp )%table(nitem,2)) then
             f_t = MWAmplitudes( jj_n_amp )%table(nitem,1)
          else if(t_t < MWAmplitudes( jj_n_amp )%table(1,2)) then
             f_t = MWAmplitudes( jj_n_amp )%table(1,1)
          else 
             do i = 2, nitem
                if(t_t <= MWAmplitudes( jj_n_amp )%table(i,2)) then
                  t_2 = MWAmplitudes( jj_n_amp )%table(i,2)
                  t_1 = MWAmplitudes( jj_n_amp )%table(i-1,2)
                  f_2 = MWAmplitudes( jj_n_amp )%table(i,1)
                  f_1 = MWAmplitudes( jj_n_amp )%table(i-1,1)
                  f_t = ((t_2*f_1 - t_1*f_2) + (f_2 - f_1)*t_t) / (t_2 - t_1)
                  exit
                endif
             end do
          end if

      end if

  end subroutine table_nlsta

!
end module m_fstr_solve_NLGEOM
