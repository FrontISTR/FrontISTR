!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                       Zhigang Sun(ASTOM)                             !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines for nonlinear implicit dynamic analysis

module fstr_dynamic_nlimplicit

use m_fstr
use m_fstr_para_contact
use lczparm
use m_dynamic_output
use m_fstr_EIG_setMASS
use m_dynamic_mat_ass_bc_ac
use m_dynamic_mat_ass_bc
use m_dynamic_mat_ass_bc_vl
use m_dynamic_mat_ass_load
use m_fstr_StiffMatrix
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
                                      ,fstrCPL, restrt_step_num )

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
    integer(kind=kint) :: iiii5, iexit
    integer(kind=kint) :: revocap_flag
    integer(kind=kint) :: kkk0, kkk1

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize, res
    real :: time_1, time_2

    integer(kind=kint) :: restrt_step_num
    integer(kind=kint) :: n_node_global

    real(kind=kreal), parameter :: PI = 3.14159265358979323846D0


!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if

! sum of n_node among all subdomains (to be used to calc res)
    n_node_global = hecMESH%nn_internal
    call hecmw_allreduce_I1(hecMESH,n_node_global,HECMW_SUM)

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
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)
    end if
	
	fstrDYNAMIC%VEC3(:) =0.d0
	hecMAT%X(:) =0.d0
	
!!
!!    step = 1,2,....,fstrDYNAMIC%n_step
!!
    do i= restrt_step_num, fstrDYNAMIC%n_step

       fstrDYNAMIC%i_step = i                                                  
       fstrDYNAMIC%t_curr = fstrDYNAMIC%t_delta * i
	
       if(hecMESH%my_rank==0) then
         write(ISTA,'('' time step='',i10,'' time='',1pe13.4e3)') i,fstrDYNAMIC%t_curr   
       endif

       do j = 1 ,ndof*nnod
           fstrDYNAMIC%VEC1(j) = a1*fstrDYNAMIC%ACC(j,1) + a2*fstrDYNAMIC%VEL(j,1)
           fstrDYNAMIC%VEC2(j) = b1*fstrDYNAMIC%ACC(j,1) + b2*fstrDYNAMIC%VEL(j,1)
       end do

!C ********************************************************************************
!C for couple analysis
      do
       fstrSOLID%dunode(:) =0.d0
      ! call fstr_UpdateEPState( hecMESH, fstrSOLID )  

       if( fstrPARAM%fg_couple == 1) then
         if( fstrPARAM%fg_couple_type==1 .or. &
             fstrPARAM%fg_couple_type==3 .or. &
             fstrPARAM%fg_couple_type==5 ) call fstr_rcap_get( fstrCPL )
       endif
!C ********************************************************************************

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

!C ********************************************************************************
!C for couple analysis
         if( fstrPARAM%fg_couple == 1) then
           if( fstrPARAM%fg_couple_first /= 0 ) then
             bsize = DFLOAT( i ) / DFLOAT( fstrPARAM%fg_couple_first )
             if( bsize > 1.0 ) bsize = 1.0
             do kkk0 = 1, fstrCPL%coupled_node_n
               kkk1 = 3 * kkk0
               fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
               fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
               fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
             enddo
           endif
           if( fstrPARAM%fg_couple_window > 0 ) then
             j = i - restrt_step_num + 1
             kk = fstrDYNAMIC%n_step - restrt_step_num + 1
             bsize = 0.5*(1.0-cos(2.0*PI*DFLOAT(j)/DFLOAT(kk)))
             do kkk0 = 1, fstrCPL%coupled_node_n
               kkk1 = 3 * kkk0
               fstrCPL%trac(kkk1-2) = bsize * fstrCPL%trac(kkk1-2)
               fstrCPL%trac(kkk1-1) = bsize * fstrCPL%trac(kkk1-1)
               fstrCPL%trac(kkk1  ) = bsize * fstrCPL%trac(kkk1  )
             enddo
           endif
           call dynamic_mat_ass_couple( hecMESH, hecMAT, fstrSOLID, fstrCPL )
         endif
!C ********************************************************************************

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

!C-- geometrical boundary condition

        call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, iter)
        call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, iter)
        call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, iter)

!C
!C-- RHS LOAD VECTOR CHECK
!C
        numnp=hecMAT%NP
        bsize=0.0
        do iiii5 = 1,numnp*ndof
          bsize=bsize+hecMAT%B(iiii5)**2
        enddo
!C-- Gather RHS vector
        if (.not. ds) then !In case of Direct Solver prevent MPI
          call hecmw_allREDUCE_R1( hecMESH,bsize,hecmw_sum )
        end if

        iexit = 0
        if( bsize < 1.0e-31 ) then
          iexit = 1
          if( hecMESH%my_rank .eq. 0 ) then
            WRITE(IMSG,*) '###Load Vector Error!'
          endif
        endif

        res = dsqrt(bsize)/n_node_global
        if( hecMESH%my_rank==0 ) then
          write(ISTA,'(''iter='',I5,''- Residual'',E15.7)')iter,res
        endif
        if( res<fstrSOLID%step_ctrl(cstep)%converg ) exit
!C
!C-- linear solver [A]{X} = {B}
!C
        if( iexit .eq. 1 ) then
          hecMAT%X = 0.0
        else
          CALL solve_LINEQ(hecMESH,hecMAT,imsg)
        end if

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

!C *****************************************************
!C for couple analysis
      if( fstrPARAM%fg_couple == 1 ) then
        if( fstrPARAM%fg_couple_type>1 ) then
          do j=1, fstrCPL%coupled_node_n
            if( fstrCPL%dof == 3 ) then
              kkk0 = j*3
              kkk1 = fstrCPL%coupled_node(j)*3

              fstrCPL%disp (kkk0-2) = fstrSOLID%unode(kkk1-2) + fstrSOLID%dunode(kkk1-2)
              fstrCPL%disp (kkk0-1) = fstrSOLID%unode(kkk1-1) + fstrSOLID%dunode(kkk1-1)
              fstrCPL%disp (kkk0  ) = fstrSOLID%unode(kkk1  ) + fstrSOLID%dunode(kkk1  )

              fstrCPL%velo (kkk0-2) = -b1*fstrDYNAMIC%ACC(kkk1-2,1) - b2*fstrDYNAMIC%VEL(kkk1-2,1) + &
                                       b3*fstrSOLID%dunode(kkk1-2)
              fstrCPL%velo (kkk0-1) = -b1*fstrDYNAMIC%ACC(kkk1-1,1) - b2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       b3*fstrSOLID%dunode(kkk1-1)
              fstrCPL%velo (kkk0  ) = -b1*fstrDYNAMIC%ACC(kkk1,1) - b2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       b3*fstrSOLID%dunode(kkk1)
              fstrCPL%accel(kkk0-2) = -a1*fstrDYNAMIC%ACC(kkk1-2,1) - a2*fstrDYNAMIC%VEL(kkk1-2,1) + &
                                       a3*fstrSOLID%dunode(kkk1-2)
              fstrCPL%accel(kkk0-1) = -a1*fstrDYNAMIC%ACC(kkk1-1,1) - a2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       a3*fstrSOLID%dunode(kkk1-1)
              fstrCPL%accel(kkk0  ) = -a1*fstrDYNAMIC%ACC(kkk1,1) - a2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       a3*fstrSOLID%dunode(kkk1)
            else
              kkk0 = j*2
              kkk1 = fstrCPL%coupled_node(j)*2

              fstrCPL%disp (kkk0-1) = fstrSOLID%unode(kkk1-1) + fstrSOLID%dunode(kkk1-1)
              fstrCPL%disp (kkk0  ) = fstrSOLID%unode(kkk1  ) + fstrSOLID%dunode(kkk1  )

              fstrCPL%velo (kkk0-1) = -b1*fstrDYNAMIC%ACC(kkk1-1,1) - b2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       b3*fstrSOLID%dunode(kkk1-1)
              fstrCPL%velo (kkk0  ) = -b1*fstrDYNAMIC%ACC(kkk1,1) - b2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       b3*fstrSOLID%dunode(kkk1)
              fstrCPL%accel(kkk0-1) = -a1*fstrDYNAMIC%ACC(kkk1-1,1) - a2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       a3*fstrSOLID%dunode(kkk1-1)
              fstrCPL%accel(kkk0  ) = -a1*fstrDYNAMIC%ACC(kkk1,1) - a2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       a3*fstrSOLID%dunode(kkk1)
            endif
          end do
          call fstr_rcap_send( fstrCPL )
        endif

        select case ( fstrPARAM%fg_couple_type )
          case (4)
            call fstr_rcap_get( fstrCPL )
          case (5)
            call fstr_get_convergence( revocap_flag )
            if( revocap_flag==0 ) cycle
          case (6)
            call fstr_get_convergence( revocap_flag )
            if( revocap_flag==0 ) then
              call fstr_rcap_get( fstrCPL )
              cycle
            else
              if( i /= fstrDYNAMIC%n_step ) call fstr_rcap_get( fstrCPL )
            endif
        end select
      endif
      exit
     enddo
!C *****************************************************

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
      if( fstrDYNAMIC%restart_nout > 0 .and. &
          (mod(i,fstrDYNAMIC%restart_nout).eq.0 .or. i.eq.fstrDYNAMIC%n_step) ) then
        call fstr_write_restart_dyna_nl(i,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM)
      end if

!C-- output new displacement, velocity and accelaration
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)

!C-- output result of monitoring node
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)

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
                                                       ,fstrCPL,fstrMAT,restrt_step_num,infoCTChange  &
                                                       ,conMAT )
  
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
    type ( hecmwST_matrix      ),optional :: conMAT

!C
!C-- local variable
!C
    
    integer(kind=kint) :: nnod, ndof, numnp, nn
    integer(kind=kint) :: i, j, ids, ide, ims, ime, kk, idm, imm
    integer(kind=kint) :: iter
    

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize, res, res1, rf                                
    real(kind=kreal) :: res0, relres
    real :: time_1, time_2

    integer(kind=kint) :: restrt_step_num
    
    integer(kind=kint) :: ctAlgo                                      
    integer(kind=kint) :: max_iter_contact, count_step                 
    integer(kind=kint) :: stepcnt                                      
    real(kind=kreal)    :: maxDLag                                             

    logical :: is_mat_symmetric
    integer(kind=kint) :: n_node_global
    integer(kind=kint) :: contact_changed_global
    

    integer(kint)   ::  nndof,npdof
    real(kreal),allocatable :: tmp_conB(:)
    integer   ::  istat

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
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)
    end if
	
	fstrDYNAMIC%VEC3(:) =0.d0  
	hecMAT%X(:) =0.d0

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
       fstrDYNAMIC%t_curr = fstrDYNAMIC%t_delta * i
	
       if(hecMESH%my_rank==0) then
         write(ISTA,'('' time step='',i10,'' time='',1pe13.4e3)') i,fstrDYNAMIC%t_curr   
         write(*,'(A)')'-------------------------------------------------'                      
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

       
           if(paraContactFlag.and.present(conMAT)) then
             call hecMAT_clear( conMAT )
             conMAT%B(:) = 0.0D0
           endif
           if( fstr_is_contact_active() ) then
!    ----  For Parallel Contact with Multi-Partition Domains
             if(paraContactFlag.and.present(conMAT)) then
               call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID,conMAT)
               call fstr_AddContactStiffness(cstep,iter,conMAT,fstrMAT,fstrSOLID) 
             else       
               call fstr_Update_NDForce_contact(cstep,hecMESH,hecMAT,fstrMAT,fstrSOLID)    
               call fstr_AddContactStiffness(cstep,iter,hecMAT,fstrMAT,fstrSOLID)   
             endif
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
          if(paraContactFlag.and.present(conMAT)) then
            call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt, conMAT=conMAT)     
            call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt, conMAT=conMAT) 
            call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt, conMAT=conMAT) 
          else
            call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt)     
            call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt) 
            call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT, stepcnt) 
          endif

! ----- check convergence
!   ----  For Parallel Contact with Multi-Partition Domains
          if(paraContactFlag.and.present(conMAT)) then
            res = fstr_get_norm_para_contact(hecMAT,fstrMAT,conMAT,hecMESH) 
          else
            call hecmw_innerProduct_R(hecMESH,ndof,hecMAT%B,hecMAT%B,res)
          endif
          res = sqrt(res)/n_node_global
          
          if( iter==1 ) res0=res
          if( res0==0.d0 ) then
            res0 =1.d0
          else
            relres = dabs(res1-res)/res0
          endif
          
! ----- check convergence
          if( .not.fstr_is_contact_active() ) then
            maxDLag = 0.0d0
          elseif( maxDLag  == 0.0D0) then
            maxDLag = 1.0D0
          endif
          call hecmw_allreduce_R1(hecMESH, maxDlag, HECMW_MAX)
          
          if( hecMESH%my_rank==0 ) then
!            if( mod(i,max(int(fstrDYNAMIC%nout/10),1)) == 0 )   &
            write(*,'(a,i3,a,2e15.7)') ' - Residual(',iter,') =',res,relres
            write(*,'(a,1e15.7)') ' - MaxDLag =',maxDLag
            write(ISTA,'(''iter='',I5,''res/res0='',2E15.7)')iter,res,relres
            write(ISTA,'(a,1e15.7)') ' - MaxDLag =',maxDLag
          endif

! ----- check convergence
!          if( .not.fstr_is_contact_active() ) maxDLag= 0.0d0                    
!          call hecmw_allreduce_R1(hecMESH, maxDlag, HECMW_MAX)
!          if( res<fstrSOLID%step_ctrl(cstep)%converg .and. maxDLag<1.0d-5 .and. iter>1 ) exit   
          if( (res<fstrSOLID%step_ctrl(cstep)%converg  .or.    &
            relres<fstrSOLID%step_ctrl(cstep)%converg) .and. maxDLag<1.0d-4 ) exit
          res1 = res
          rf=1.0d0                                                   
          if( iter>1 .and. res>res1 )rf=0.5d0*rf                         
          res1=res                                                                          

!          call solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,rf)

!   ----  For Parallel Contact with Multi-Partition Domains
          if(paraContactFlag.and.present(conMAT)) then
            call solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,1.0D0,conMAT)
          else
            call solve_LINEQ_contact(hecMESH,hecMAT,fstrMAT,rf)
          endif
                   
! ----- update external nodal displacement increments
          call hecmw_update_3_R (hecMESH, hecMAT%X, hecMAT%NP)
          
! ----- update the strain, stress, and internal force
!          call fstr_UpdateNewton( hecMESH, hecMAT, !fstrSOLID,fstrDYNAMIC%t_delta,1,fstrDYNAMIC%strainEnergy ) 
          do j=1,hecMESH%n_node*ndof
            fstrSOLID%dunode(j)  = fstrSOLID%dunode(j)+hecMAT%X(j)
          enddo
          call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,fstrDYNAMIC%t_delta,1,fstrDYNAMIC%strainEnergy ) 

    
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
!     ----  For Parallel Contact with Multi-Partition Domains
          if(paraContactFlag.and.present(conMAT)) then
            call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange, conMAT) 
          else
            call fstr_mat_con_contact( cstep, hecMAT, fstrSOLID, fstrMAT, infoCTChange)    
          endif  
          contact_changed_global=1
        endif
        call hecmw_allreduce_I1(hecMESH,contact_changed_global,HECMW_MAX)
        if (contact_changed_global > 0) then
          hecMAT%B(:) = 0.0D0
          if(paraContactFlag.and.present(conMAT)) conMAT%B(:) = 0.0D0
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
      if( fstrDYNAMIC%restart_nout > 0 .and. &
          (mod(i,fstrDYNAMIC%restart_nout).eq.0 .or. i.eq.fstrDYNAMIC%n_step) ) then
        call fstr_write_restart_dyna_nl(i,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
                                     infoCTChange%contactNode_current)
      end if

!C-- output new displacement, velocity and accelaration
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)

!C-- output result of monitoring node
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)

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
