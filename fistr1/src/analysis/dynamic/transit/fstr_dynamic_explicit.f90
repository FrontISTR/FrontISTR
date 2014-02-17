!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Tomotaka Ogasawara (Univ. of Tokyo)            !
!                       Xi YUAN( AdvanceSoft )                         !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains subroutines for explicit dynamic analysis

module fstr_dynamic_explicit

use m_fstr
use lczparm
use m_static_lib
use m_dynamic_output
use m_fstr_EIG_setMASS
use m_dynamic_mat_ass_bc_ac
use m_dynamic_mat_ass_bc
use m_dynamic_mat_ass_bc_vl
use m_dynamic_mat_ass_load
use m_static_mat_ass_main
use fstr_matrix_con_contact
use m_fstr_restart

!-------- for couple -------
use m_dynamic_mat_ass_couple
use m_fstr_rcap_io


contains

!C================================================================C
!C-- subroutine  fstr_solve_LINEAR_DYNAMIC
!C================================================================C
  subroutine fstr_solve_dynamic_explicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
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
    integer(kind=kint) :: nnod, ndof, nn, numnp, n_nod_dof
    integer(kind=kint) :: i, j, ids, ide, ims, ime, kk, idm, imm
    integer(kind=kint) :: kkk0, kkk1
    integer(kind=kint) :: ierror, idummy
    integer(kind=kint) :: iiii5, iexit
    integer(kind=kint) :: revocap_flag
    real(kind=kreal),pointer :: prevB(:)

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize
    real(kind=kreal) :: time_1, time_2

    integer(kind=kint) :: restrt_step_num
    integer(kind=kint) :: restrt_step(1)

    real(kind=kreal), parameter :: PI = 3.14159265358979323846D0

!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if
	  
      call cpu_time( time_1 )

!--
    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof
    n_nod_dof = nnod*ndof

!C
!C-- matrix [KL] & [M]
!C-- matrix [KL]

      call fstr_mat_ass_main (hecMESH, hecMAT, fstrSOLID)

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
    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        allocate( prevB(hecMAT%NP*ndof)      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, prevB>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            endif
      endif
    endif

!C
!C-- time step loop
!C
    a1 = 1./fstrDYNAMIC%t_delta**2
    a2 = 1./(2.0*fstrDYNAMIC%t_delta)

!C-- output initial condition
    if( restrt_step_num .eq. 1 ) then

      do j = 1 ,ndof*nnod
        fstrDYNAMIC%DISP(j,3) = fstrDYNAMIC%DISP(j,1) - fstrDYNAMIC%VEL (j,1)/(2.0*a2) &
                            + fstrDYNAMIC%ACC (j,1)/ (2.0*a1)
        fstrDYNAMIC%DISP(j,2) = fstrDYNAMIC%DISP(j,1) - fstrDYNAMIC%VEL (j,1)/ a2 &
                            + fstrDYNAMIC%ACC (j,1)/ (2.0*a1) * 4.0
      end do

!C-- output new displacement, velocity and accelaration
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)

!C-- output result of monitoring node
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)

    end if
!!
!!    step = 1,2,....,fstrDYNAMIC%n_step
!!

    do i= restrt_step_num, fstrDYNAMIC%n_step

       fstrDYNAMIC%i_step = i
       fstrDYNAMIC%t_curr = fstrDYNAMIC%t_delta * i

!C-- {vec3}=[KL]{U(t)}, second part of right hand of Eq.(1.1.22)

      do j = 1 ,ndof*nnod
         fstrDYNAMIC%VEC1(j) = fstrDYNAMIC%DISP(j,1)
      end do
      if( hecMESH%n_dof .eq. 3 ) then
        call hecmw_matvec_33 (hecMESH, hecMAT, fstrDYNAMIC%VEC1, fstrDYNAMIC%VEC3)
      else if( hecMESH%n_dof .eq. 2 ) then
        call hecmw_matvec_22 (hecMESH, hecMAT, fstrDYNAMIC%VEC1, fstrDYNAMIC%VEC3, nnod)
      else if( hecMESH%n_dof .eq. 6 ) then
        call matvec(fstrDYNAMIC%VEC3, fstrDYNAMIC%VEC1, hecMAT, ndof, hecMAT%D, hecMAT%AU, hecMAT%AL)
      end if
!
!C-- matrix [A] = {VEC1}, Eq.(1.1.23)

    do j = 1 ,ndof*nnod
        fstrDYNAMIC%VEC1(j) = (a1 + a2 *fstrDYNAMIC%ray_m) * myEIG%mass(j)
    end do

!C
!C-- mechanical boundary condition

        call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)

!C
!C Eq. (1.1.22) but without Aj
    do j = 1 ,ndof*nnod
      hecMAT%B(j) = hecMAT%B(j) - fstrDYNAMIC%VEC3(j) + 2.0*a1* myEIG%mass(j) * fstrDYNAMIC%DISP(j,1)  &
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

!C
!C-- RHS LOAD VECTOR CHECK
!C
      numnp=hecMAT%NP
      bsize=0.0
      do iiii5 = 1,numnp*ndof
        bsize=bsize+hecMAT%B(iiii5)**2
      enddo

!C-- Gather RHS vector
!C
        call hecmw_allREDUCE_R1( hecMESH,bsize,hecmw_sum )

        if( hecMESH%my_rank .eq. 0 ) then
          write(IMSG,*) 'Total RHS size=',bsize
        endif

        iexit = 0

        if( bsize < 1.0e-31 ) then
          iexit = 1
!         if( hecMESH%my_rank .eq. 0 ) then
!           WRITE(IMSG,*) '###Load Vector Error!'
!         endif
!!!       call hecmw_abort( hecmw_comm_get_comm())
        endif
!C
!C-- check parameters
!C
!!
      do j = 1 ,ndof*nnod
         if(dabs(fstrDYNAMIC%VEC1(j)) .lt. 1.0e-20) then
              if( hecMESH%my_rank .eq. 0 ) then
                 write(imsg,*) 'stop due to fstrDYNAMIC%VEC(j) = 0 ,  j = ', j
              end if
              call hecmw_abort( hecmw_comm_get_comm())
         end if
      end do
!
!!        if( iexit .eq. 1 ) then
!!          hecMAT%X = 0.0
!!        else

! Finish the calculation indicited by Eq. (1.1.22)
      do j = 1 ,ndof*nnod
          hecMAT%X(j) = hecMAT%B(j) / fstrDYNAMIC%VEC1(j)
          if(dabs(hecMAT%X(j)) > 1.0d+5) then
              if( hecMESH%my_rank == 0 ) then
                 print *, 'Displacement increment too large, please adjust your step size!'
                 write(imsg,*) 'Displacement increment too large, please adjust your step size!'
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
	  ! Eq. (1.1.18)
      fstrDYNAMIC%ACC (j,1) = a1*(hecMAT%X(j) - 2.d0*fstrDYNAMIC%DISP(j,1) &
                                  + fstrDYNAMIC%DISP(j,3))
      ! Eq.(1.1.17)
      fstrDYNAMIC%VEL (j,1) = a2*(hecMAT%X(j) - fstrDYNAMIC%DISP(j,3))

      fstrDYNAMIC%DISP(j,3) = fstrDYNAMIC%DISP(j,1)
      fstrDYNAMIC%DISP(j,1) = hecMAT%X(j)
    end do

!!! restart  !!!
    if( fstrDYNAMIC%restart_nout > 0 .and. &
        (mod(i,fstrDYNAMIC%restart_nout).eq.0 .or. i.eq.fstrDYNAMIC%n_step) ) then
      call fstr_write_restart_dyna_linear(i,fstrDYNAMIC)
    end if

!C-- output new displacement, velocity and accelaration
      call fstr_dynamic_Output(hecMESH, fstrSOLID, fstrDYNAMIC)

!C-- output result of monitoring node
      call dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)

    enddo

    if( fstrPARAM%fg_couple == 1) then
      if( fstrPARAM%fg_couple_type==5 .or. &
          fstrPARAM%fg_couple_type==6 ) then
        deallocate( prevB      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, prevB>'
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
!C
!C-- end of time step loop

  end subroutine fstr_solve_dynamic_explicit


end module fstr_dynamic_explicit
