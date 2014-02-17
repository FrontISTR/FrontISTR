!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!> \brief This module contains subroutines for nonlinear explicit dynamic analysis

module fstr_dynamic_nlexplicit

use m_fstr
use lczparm
use m_static_lib
use m_dynamic_output
use m_fstr_EIG_setMASS
use m_dynamic_mat_ass_bc_ac
use m_dynamic_mat_ass_bc
use m_dynamic_mat_ass_bc_vl
use m_dynamic_mat_ass_load
use m_fstr_Update
use m_fstr_Restart
use fstr_matrix_con_contact

!-------- for couple -------
use m_dynamic_mat_ass_couple
use m_fstr_rcap_io


contains

!C================================================================C
!C-- subroutine  fstr_solve_LINEAR_DYNAMIC
!C================================================================C
  subroutine fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, restrt_step_num )

    implicit none
!C
!C-- global variable
!C
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( lczparam            ) :: myEIG
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM
      type ( fstr_dynamic        ) :: fstrDYNAMIC
      type (fstrST_matrix_contact_lagrange)  :: fstrMAT  !< type fstrST_matrix_contact_lagrange  
      type ( fstr_couple         ) :: fstrCPL         !for COUPLE      

!C
!C-- local variable
!C
    integer(kind=kint) :: nnod, ndof, nn, numnp
    integer(kind=kint) :: i, j, ids, ide, kk
    integer(kind=kint) :: kkk0, kkk1
    integer(kind=kint) :: ierror
    integer(kind=kint) :: iiii5, iexit
    integer(kind=kint) :: revocap_flag
    real(kind=kreal),pointer :: prevB(:)

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize, res
    real(kind=kreal) :: time_1, time_2

    integer(kind=kint) :: restrt_step_num

    real(kind=kreal), parameter :: PI = 3.14159265358979323846D0


    call cpu_time( time_1 )

!--
    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof

!C--
    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        allocate( prevB(hecMAT%NP*ndof)      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            endif
      endif
    endif

	fstrSOLID%dunode(:) =0.d0

    a1 = 1.d0/fstrDYNAMIC%t_delta**2
    a2 = 1.d0/(2.d0*fstrDYNAMIC%t_delta)
	
    call setMASS(IDBG,hecMESH,hecMAT,myEIG)

    do j = 1 ,ndof*nnod
        fstrDYNAMIC%VEC1(j) = (a1 + a2 *fstrDYNAMIC%ray_m) * myEIG%mass(j)
        if(dabs(fstrDYNAMIC%VEC1(j)) < 1.0e-20) then
          if( hecMESH%my_rank == 0 ) then
            write(*,*) 'stop due to fstrDYNAMIC%VEC(j) = 0 ,  j = ', j
            write(imsg,*) 'stop due to fstrDYNAMIC%VEC(j) = 0 ,  j = ', j
          end if
          call hecmw_abort( hecmw_comm_get_comm())
        endif
    end do
	

!C-- output of initial state
    if( restrt_step_num == 1 ) then
      do j = 1 ,ndof*nnod
        fstrDYNAMIC%DISP(j,3) = fstrDYNAMIC%DISP(j,1) - fstrDYNAMIC%VEL (j,1)/(2.d0*a2) &
                            + fstrDYNAMIC%ACC (j,1)/ (2.d0*a1)
        fstrDYNAMIC%DISP(j,2) = fstrDYNAMIC%DISP(j,1) - fstrDYNAMIC%VEL (j,1)/ a2 &
                            + fstrDYNAMIC%ACC (j,1)/ (2.d0*a1) * 4.d0
      end do

      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)
    end if
	

    do i= restrt_step_num, fstrDYNAMIC%n_step

       fstrDYNAMIC%i_step = i
       fstrDYNAMIC%t_curr = fstrDYNAMIC%t_delta * i
!C
!C-- mechanical boundary condition

        call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        do j=1, hecMESH%n_node*  hecMESH%n_dof 
          hecMAT%B(j)=hecMAT%B(j)-fstrSOLID%QFORCE(j)
        end do

        do j = 1 ,ndof*nnod
          hecMAT%B(j) = hecMAT%B(j) + 2.d0*a1* myEIG%mass(j) * fstrDYNAMIC%DISP(j,1)  &
                 + (- a1 + a2 * fstrDYNAMIC%ray_m) * myEIG%mass(j) * fstrDYNAMIC%DISP(j,3)
        end do

!C ********************************************************************************
!C for couple analysis
        if( fstrPARAM%fg_couple == 1 ) then
          if( fstrPARAM%fg_couple_type==5 .or. &
              fstrPARAM%fg_couple_type==6 ) then
            do j = 1, hecMAT%NP * ndof
              prevB(j) = hecMAT%B(j)
            enddo
          endif
        endif
    do
        if( fstrPARAM%fg_couple == 1 ) then
          if( fstrPARAM%fg_couple_type==1 .or. &
              fstrPARAM%fg_couple_type==3 .or. &
              fstrPARAM%fg_couple_type==5 ) call fstr_rcap_get( fstrCPL )
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

!C
!C-- geometrical boundary condition

        call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT) 
        call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT) 
        call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, fstrMAT) 

! Finish the calculation 
        do j = 1 ,ndof*nnod
          hecMAT%X(j) = hecMAT%B(j) / fstrDYNAMIC%VEC1(j)
          if(dabs(hecMAT%X(j)) > 1.0d+5) then
              if( hecMESH%my_rank == 0 ) then
                 print *, 'Displacement increment too large, please adjust your step size!',i
                 write(imsg,*) 'Displacement increment too large, please adjust your step size!',i,hecMAT%B(j),fstrDYNAMIC%VEC1(j)
              end if
              call hecmw_abort( hecmw_comm_get_comm())
          end if
        end do

!C *****************************************************
!C for couple analysis
      if( fstrPARAM%fg_couple == 1 ) then
        if( fstrPARAM%fg_couple_type>1 ) then
          do j=1, fstrCPL%coupled_node_n
            if( fstrCPL%dof == 3 ) then
              kkk0 = j*3
              kkk1 = fstrCPL%coupled_node(j)*3

              fstrCPL%disp (kkk0-2) = hecMAT%X(kkk1-2)
              fstrCPL%disp (kkk0-1) = hecMAT%X(kkk1-1)
              fstrCPL%disp (kkk0  ) = hecMAT%X(kkk1  )

              fstrCPL%velo (kkk0-2) = -b1*fstrDYNAMIC%ACC(kkk1-2,1) - b2*fstrDYNAMIC%VEL(kkk1-2,1) + &
                                       b3*( hecMAT%X(kkk1-2) - fstrDYNAMIC%DISP(kkk1-2,1) )
              fstrCPL%velo (kkk0-1) = -b1*fstrDYNAMIC%ACC(kkk1-1,1) - b2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       b3*( hecMAT%X(kkk1-1) - fstrDYNAMIC%DISP(kkk1-1,1) )
              fstrCPL%velo (kkk0  ) = -b1*fstrDYNAMIC%ACC(kkk1,1) - b2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       b3*( hecMAT%X(kkk1) - fstrDYNAMIC%DISP(kkk1,1) )
              fstrCPL%accel(kkk0-2) = -a1*fstrDYNAMIC%ACC(kkk1-2,1) - a2*fstrDYNAMIC%VEL(kkk1-2,1) + &
                                       a3*( hecMAT%X(kkk1-2) - fstrDYNAMIC%DISP(kkk1-2,1) )
              fstrCPL%accel(kkk0-1) = -a1*fstrDYNAMIC%ACC(kkk1-1,1) - a2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       a3*( hecMAT%X(kkk1-1) - fstrDYNAMIC%DISP(kkk1-1,1) )
              fstrCPL%accel(kkk0  ) = -a1*fstrDYNAMIC%ACC(kkk1,1) - a2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       a3*( hecMAT%X(kkk1) - fstrDYNAMIC%DISP(kkk1,1) )
            else
              kkk0 = j*2
              kkk1 = fstrCPL%coupled_node(j)*2

              fstrCPL%disp (kkk0-1) = hecMAT%X(kkk1-1)
              fstrCPL%disp (kkk0  ) = hecMAT%X(kkk1  )

              fstrCPL%velo (kkk0-1) = -b1*fstrDYNAMIC%ACC(kkk1-1,1) - b2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       b3*( hecMAT%X(kkk1-1) - fstrDYNAMIC%DISP(kkk1-1,1) )
              fstrCPL%velo (kkk0  ) = -b1*fstrDYNAMIC%ACC(kkk1,1) - b2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       b3*( hecMAT%X(kkk1) - fstrDYNAMIC%DISP(kkk1,1) )
              fstrCPL%accel(kkk0-1) = -a1*fstrDYNAMIC%ACC(kkk1-1,1) - a2*fstrDYNAMIC%VEL(kkk1-1,1) + &
                                       a3*( hecMAT%X(kkk1-1) - fstrDYNAMIC%DISP(kkk1-1,1) )
              fstrCPL%accel(kkk0  ) = -a1*fstrDYNAMIC%ACC(kkk1,1) - a2*fstrDYNAMIC%VEL(kkk1,1) + &
                                       a3*( hecMAT%X(kkk1) - fstrDYNAMIC%DISP(kkk1,1) )
            endif
          end do
          call fstr_rcap_send( fstrCPL )
        endif

        select case ( fstrPARAM%fg_couple_type )
        case (4)
          call fstr_rcap_get( fstrCPL )
        case (5)
          call fstr_get_convergence( revocap_flag )
          if( revocap_flag==0 ) then
            do j = 1, hecMAT%NP * ndof
              hecMAT%B(j) = prevB(j)
            enddo
            cycle
          endif
        case (6)
          call fstr_get_convergence( revocap_flag )
          if( revocap_flag==0 ) then
            do j = 1, hecMAT%NP * ndof
              hecMAT%B(j) = prevB(j)
            enddo
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

!C
!C-- new displacement, velocity and accelaration
!C
      do j = 1 ,ndof*nnod
        fstrDYNAMIC%ACC (j,1) = a1*(hecMAT%X(j) - 2.d0*fstrDYNAMIC%DISP(j,1) &
                                  + fstrDYNAMIC%DISP(j,3))
        fstrDYNAMIC%VEL (j,1) = a2*(hecMAT%X(j) - fstrDYNAMIC%DISP(j,3))

        fstrDYNAMIC%DISP(j,3) = fstrDYNAMIC%DISP(j,1)
        fstrDYNAMIC%DISP(j,1) = hecMAT%X(j)
		
        fstrSOLID%unode(j)  = fstrDYNAMIC%DISP(j,3)
        hecMAT%X(j) = fstrDYNAMIC%DISP(j,1)-fstrDYNAMIC%DISP(j,3)
      end do
	
! ----- update strain, stress, and internal force
      call fstr_UpdateNewton( hecMESH, hecMAT, fstrSOLID,fstrDYNAMIC%t_delta,1 )

      if( fstrDYNAMIC%restart_nout > 0 .and. &
          (mod(i,fstrDYNAMIC%restart_nout).eq.0 .or. i.eq.fstrDYNAMIC%n_step) ) then
        call fstr_write_restart_dyna_nl(i,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM)
      end if
!
!C-- output new displacement, velocity and accelaration
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)
	  
      call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYNAMIC%t_delta )
	  
    enddo

    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        deallocate( prevB      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_NONLINEAR_DYNAMIC, prevB>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            endif
      endif
    endif

	call cpu_time(time_2)
    if( hecMESH%my_rank == 0 ) then
        write(ISTA,'(a,f10.2)') '         solve (sec) :', time_2 - time_1
    end if

  end subroutine fstr_solve_dynamic_nlexplicit

end module fstr_dynamic_nlexplicit
