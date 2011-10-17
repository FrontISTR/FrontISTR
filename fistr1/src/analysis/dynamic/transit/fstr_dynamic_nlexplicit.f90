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
!> \brief This module contains subroutines for nonlinear explicit dynamic analysis

module fstr_dynamic_nlexplicit

use m_fstr
use lczparm
use m_static_lib
use m_static_make_result
use m_dynamic_output
use m_fstr_EIG_setMASS
use m_dynamic_mat_ass_bc_ac
use m_dynamic_mat_ass_bc
use m_dynamic_mat_ass_bc_vl
use m_dynamic_mat_ass_load
use m_dynamic_post
use m_fstr_Update
use m_fstr_Restart

!-------- for couple -------
use m_dynamic_mat_ass_couple
use m_fstr_rcap_io


contains

!C================================================================C
!C-- subroutine  fstr_solve_LINEAR_DYNAMIC
!C================================================================C
  subroutine fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL, my_rank_monit_1, restrt_step_num )

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
      type ( fstr_couple         ) :: fstrCPL         !for COUPLE

!C
!C-- local variable
!C
    integer(kind=kint) :: nnod, ndof, nn, numnp
    integer(kind=kint) :: i, j, ids, ide, kk
    integer(kind=kint) :: kkk0, kkk1
    integer(kind=kint) :: ierror
    integer(kind=kint) :: iiii5, iexit
    integer(kind=kint) :: my_rank_monit_1

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize, res
    real(kind=kreal) :: time_1, time_2

    integer(kind=kint) :: restrt_step_num

	  
    call cpu_time( time_1 )
	  
!--
    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof
	
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

      call dynamic_nloutput(0,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, my_rank_monit_1)
    end if
	

    do i= restrt_step_num, fstrDYNAMIC%n_step

       fstrDYNAMIC%i_step = i
!C
!C-- mechanical boundary condition

        call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        do j=1, hecMESH%n_node*  hecMESH%n_dof 
          hecMAT%B(j)=hecMAT%B(j)-fstrSOLID%QFORCE(j)
        end do

!C ********************************************************************************
!C for couple analysis
        if( fstrPARAM%fg_couple ==1 ) then
            if( fstrDYNAMIC%i_step > 1 .or. &
               (fstrDYNAMIC%i_step==1 .and. fstrPARAM%fg_couple_first==1 )) then
                  call fstr_rcap_get( fstrCPL )
                  call dynamic_mat_ass_couple( hecMESH, hecMAT, fstrSOLID, fstrCPL )
            endif
        endif
!C ********************************************************************************

        do j = 1 ,ndof*nnod
          hecMAT%B(j) = hecMAT%B(j) + 2.d0*a1* myEIG%mass(j) * fstrDYNAMIC%DISP(j,1)  &
                 + (- a1 + a2 * fstrDYNAMIC%ray_m) * myEIG%mass(j) * fstrDYNAMIC%DISP(j,3)
        end do
!C
!C-- geometrical boundary condition

        call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)

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


      if( mod(i,fstrDYNAMIC%restart_nout) == 0 ) then
        call fstr_write_restart_dyna(i,0,hecMESH,fstrSOLID,fstrDYNAMIC)
      end if
!
!C-- output new displacement, velocity and accelaration
      call dynamic_nloutput(i,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, my_rank_monit_1)
	  
      call fstr_UpdateState( hecMESH, fstrSOLID, fstrDYNAMIC%t_delta )
	  
      if( hecMESH%my_rank==0 ) then
          if( mod(i,int(fstrDYNAMIC%nout/100)) == 0 ) then
             res = maxval( dabs(fstrDYNAMIC%DISP(:,1)) )
             write(*,*) "Time:",real(i*fstrDYNAMIC%t_delta),"Max disp.=",real(res)
             write(ISTA,*) "Time:",real(i*fstrDYNAMIC%t_delta),"Max disp.=",real(res)
          endif
      endif
    enddo
	
	call cpu_time(time_2)
    if( hecMESH%my_rank == 0 ) then
        write(ISTA,'(a,f10.2)') '         solve (sec) :', time_2 - time_1
    end if

  end subroutine fstr_solve_dynamic_nlexplicit


end module fstr_dynamic_nlexplicit
