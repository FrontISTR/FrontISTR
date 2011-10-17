!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                                                                      !
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

!-------- for couple -------
use m_dynamic_mat_ass_couple
use m_fstr_rcap_io


contains

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
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, my_rank_monit_1)
    end if
	
	fstrDYNAMIC%VEC3(:) =0.d0
	hecMAT%X(:) =0.d0
	
!!
!!    step = 1,2,....,fstrDYNAMIC%n_step
!!
    do i= restrt_step_num, fstrDYNAMIC%n_step
	
       if(hecMESH%my_rank==0) then
         write(ISTA,*) ' time step=',i, ' time=',real(i*fstrDYNAMIC%t_delta)
         if( mod(i,int(fstrDYNAMIC%nout/10)) == 0 )   &
            write(*,*) ' time step=',i, ' time=',real(i*fstrDYNAMIC%t_delta)
       endif

       fstrDYNAMIC%i_step = i
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

        call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, iter)
        call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, iter)
        call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, iter)

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
        call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,fstrDYNAMIC%t_delta,1 )
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
      do j = 1 ,ndof*nnod
          fstrDYNAMIC%ACC (j,2) = -a1*fstrDYNAMIC%ACC(j,1) - a2*fstrDYNAMIC%VEL(j,1) + &
                               a3*fstrSOLID%dunode(j)
          fstrDYNAMIC%VEL (j,2) = -b1*fstrDYNAMIC%ACC(j,1) - b2*fstrDYNAMIC%VEL(j,1) + &
                               b3*fstrSOLID%dunode(j)
          fstrDYNAMIC%ACC (j,1) = fstrDYNAMIC%ACC (j,2)
          fstrDYNAMIC%VEL (j,1) = fstrDYNAMIC%VEL (j,2)

          fstrSOLID%unode(j)  = fstrSOLID%unode(j)+fstrSOLID%dunode(j)
          fstrDYNAMIC%DISP(j,2) = fstrSOLID%unode(j)
      enddo

!---  Restart info
      if( mod(i,fstrDYNAMIC%restart_nout) == 0 ) then
        call fstr_write_restart_dyna(i,0,hecMESH,fstrSOLID,fstrDYNAMIC)
      end if
!
!C-- output new displacement, velocity and accelaration
      call dynamic_nloutput(i,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
!C
!C-- output result of monitoring node
!C
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, my_rank_monit_1)
	  
      call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYNAMIC%t_delta )

    enddo
!C
!C-- end of time step loop

    call cpu_time(time_2)
    if( hecMESH%my_rank == 0 ) then
        write(ISTA,'(a,f10.2,a)') '         solve (sec) :', time_2 - time_1, 's'
    end if

  end subroutine fstr_solve_dynamic_nlimplicit


end module fstr_dynamic_nlimplicit
