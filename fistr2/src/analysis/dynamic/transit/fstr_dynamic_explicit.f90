!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Tomotaka Ogasawara (Univ. of Tokyo)            !
!                                                                      !
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
use m_static_make_result
use m_dynamic_output
use m_fstr_EIG_setMASS
use m_dynamic_mat_ass_bc_ac
use m_dynamic_mat_ass_bc
use m_dynamic_mat_ass_bc_vl
use m_dynamic_mat_ass_load
use m_static_mat_ass_main
use m_dynamic_post

!-------- for couple -------
use m_dynamic_mat_ass_couple
use m_fstr_rcap_io


contains

!C================================================================C
!C-- subroutine  fstr_solve_LINEAR_DYNAMIC
!C================================================================C
  subroutine fstr_solve_dynamic_explicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
                                      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
                                      ,fstrCPL )

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
    integer(kind=kint),parameter :: dynamic_IR1        =   70
    integer(kind=kint),parameter :: dynamic_IW1        =   71
    integer(kind=kint),parameter :: dynamic_IW2        =   72
    integer(kind=kint),parameter :: dynamic_IW3        =   73
    integer(kind=kint),parameter :: dynamic_IW4        =   74
    integer(kind=kint),parameter :: dynamic_IW5        =   75
    integer(kind=kint),parameter :: dynamic_IW6        =   76

    integer(kind=kint) :: nnod, ndof, nn, numnp, n_nod_dof
    integer(kind=kint) :: i, j, ids, ide, ims, ime, kk, idm, imm
    integer(kind=kint) :: kkk0, kkk1
    integer(kind=kint) :: ierror, istep, idummy
    integer(kind=kint) :: iiii5, iexit, stat
    integer(kind=kint) :: i_flg_1,i_flg_2
    integer(kind=kint) :: my_rank_monit_1, my_rank_monit_2

    real(kind=kreal) :: a1, a2, a3, b1, b2, b3, c1, c2
    real(kind=kreal) :: bsize
    real :: time_1, time_2

    integer(kind=kint) :: restrt_step_num
    integer(kind=kint) :: restrt_step(1)

!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if

!--
    hecMAT%NDOF=hecMESH%n_dof

    nnod=hecMESH%n_node
    ndof=hecMAT%NDOF
    nn=ndof*ndof
    n_nod_dof = nnod*ndof
!!
!!-- initial value
!!
      my_rank_monit_1 = 999999999
      my_rank_monit_2 = 999999999
      i_flg_1 = 0
      i_flg_2 = 0

      if( i_flg_1 .ne. 1 ) then
        do i= 1, hecMESH%nn_internal
        if( fstrPARAM%global_local_id(1,i) .eq. fstrDYNAMIC%node_monit_1 ) then
          fstrDYNAMIC%node_monit_1 = i
          i_flg_1 = 1
          do j= 1, hecMESH%nn_internal
!!          if( hecMESH%node_ID(2*j-1) .eq. i ) then
            if( j .eq. i ) then
              my_rank_monit_1 = hecMESH%node_ID(2*j)
              exit
            end if
          end do
          exit
        end if
        enddo
      end if
!
      if( fstrDYNAMIC%restart_nout .eq. 0 ) then
        fstrDYNAMIC%restart_nout = 999999999
      end if
!C
!C-- file open for local use
!C
      if( hecMESH%my_rank .eq. my_rank_monit_1 ) then
        OPEN(dynamic_IW4,FILE='dyna_disp_p1.out', status = 'replace', iostat=stat)
        if( stat /= 0 ) then
          write(*,*) 'stop due to file opening error <dyna_disp_p1.out>'
          call hecmw_abort( hecmw_comm_get_comm())
        end if
        OPEN(dynamic_IW5,FILE='dyna_velo_p1.out', status = 'replace', iostat=stat)
        if( stat /= 0 ) then
          write(*,*) 'stop due to file opening error <dyna_velo_p1.out>'
          call hecmw_abort( hecmw_comm_get_comm())
        end if
        OPEN(dynamic_IW6,FILE='dyna_acce_p1.out', status = 'replace', iostat=stat)
        if( stat /= 0 ) then
          write(*,*) 'stop due to file opening error <dyna_acce_p1.out>'
          call hecmw_abort( hecmw_comm_get_comm())
        end if
      endif


!C
!C-- ALLOCATE
!C
    allocate( myEIG%mass(ndof*nnod)    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, mass>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if


!C
!C-- initial condition
!C
!!! restart  !!!
      restrt_step(1) = 1
      restrt_step_num = 1
      fstrDYNAMIC%i_step = 0

      fstrDYNAMIC%DISP = 0.0
      fstrDYNAMIC%VEL  = 0.0
      fstrDYNAMIC%ACC  = 0.0

    if(fstrDYNAMIC%restart_nout .ge. 0 ) then
      call dynamic_bc_init   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
      call dynamic_bc_init_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
      call dynamic_bc_init_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)

    else if(fstrDYNAMIC%restart_nout .lt. 0 ) then
  !    call hecmw_restart_open()
  !    call hecmw_restart_read_int(restrt_step)
  !    call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
  !    call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,3))
  !    call hecmw_restart_close()
      restrt_step_num = restrt_step(1) + 1
    end if
!C
!C-- matrix [KL] & [M]
!C-- matrix [KL]

      call fstr_mat_ass_main(fstrSOLID)

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

!C
!C-- time step loop
!C
    a1 = 1./fstrDYNAMIC%t_delta**2
    a2 = 1./(2.0*fstrDYNAMIC%t_delta)
	
! Restart output
    if( restrt_step_num .eq. 1 ) then

      do j = 1 ,ndof*nnod
        fstrDYNAMIC%DISP(j,3) = fstrDYNAMIC%DISP(j,1) - fstrDYNAMIC%VEL (j,1)/(2.0*a2) &
                            + fstrDYNAMIC%ACC (j,1)/ (2.0*a1)
        fstrDYNAMIC%DISP(j,2) = fstrDYNAMIC%DISP(j,1) - fstrDYNAMIC%VEL (j,1)/ a2 &
                            + fstrDYNAMIC%ACC (j,1)/ (2.0*a1) * 4.0
      end do

!C-- output initial condition

      istep = 0
      call dynamic_output_and_post(istep,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
!
!C-- output result of monitoring node
!
      call dynamic_output_monit   (hecMESH, fstrPARAM, fstrDYNAMIC, &
                                   dynamic_IW4, dynamic_IW5, dynamic_IW6, my_rank_monit_1)
    end if
!!
!!    step = 1,2,....,fstrDYNAMIC%n_step
!!

    do i= restrt_step_num, fstrDYNAMIC%n_step
!
       call cpu_time(time_1)
!
       fstrDYNAMIC%i_step = i

!C-- {vec3}=[KL]{U(t)}, second part of right hand of Eq.(1.1.22)

      do j = 1 ,ndof*nnod
         fstrDYNAMIC%VEC1(j) = fstrDYNAMIC%DISP(j,1)
      end do
      if( hecMESH%n_dof .eq. 3 ) then
    !    call hecmw_matvec_33 (hecMESH, hecMAT, fstrDYNAMIC%VEC1, fstrDYNAMIC%VEC3)
      else if( hecMESH%n_dof .eq. 2 ) then
    !    call hecmw_matvec_22 (hecMESH, hecMAT, fstrDYNAMIC%VEC1, fstrDYNAMIC%VEC3, nnod)
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
!C Eq. (1.1.22) but without Aj
    do j = 1 ,ndof*nnod
      hecMAT%B(j) = hecMAT%B(j) - fstrDYNAMIC%VEC3(j) + 2.0*a1* myEIG%mass(j) * fstrDYNAMIC%DISP(j,1)  &
                 + (- a1 + a2 * fstrDYNAMIC%ray_m) * myEIG%mass(j) * fstrDYNAMIC%DISP(j,3)
    end do
!C
!C-- geometrical boundary condition

        call dynamic_mat_ass_bc   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        call dynamic_mat_ass_bc_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
        call dynamic_mat_ass_bc_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)

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
       ! call hecmw_allREDUCE_R1( hecMESH,bsize,hecmw_sum )

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
!C
!C-- new displacement, velocity and accelaration
!C
    do j = 1 ,ndof*nnod
	  ! Eq. (1.1.18)
      fstrDYNAMIC%ACC (j,1) = a1*(hecMAT%X(j) - 2.0*fstrDYNAMIC%DISP(j,1) &
                                  + fstrDYNAMIC%DISP(j,3))
      ! Eq.(1.1.17)
      fstrDYNAMIC%VEL (j,1) = a2*(hecMAT%X(j) - fstrDYNAMIC%DISP(j,3))

      fstrDYNAMIC%DISP(j,3) = fstrDYNAMIC%DISP(j,1)
      fstrDYNAMIC%DISP(j,1) = hecMAT%X(j)
    end do

!!! restart  !!!
    if(fstrDYNAMIC%restart_nout .lt. 0) then
      fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
    end if
    if( mod(i,fstrDYNAMIC%restart_nout) == 0 ) then
      restrt_step(1) = i
   !   call hecmw_restart_add_int(restrt_step,size(restrt_step))
   !   call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,1),size(fstrDYNAMIC%DISP(:,1)))
   !   call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,3),size(fstrDYNAMIC%DISP(:,3)))
   !   call hecmw_restart_write()
    end if
!
!C-- output new displacement, velocity and accelaration
      istep = i
      call dynamic_output_and_post(istep,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
!C
!C-- output result of monitoring node
!C
      call dynamic_output_monit   (hecMESH, fstrPARAM, fstrDYNAMIC, &
                                   dynamic_IW4, dynamic_IW5, dynamic_IW6, my_rank_monit_1)
!!
      call cpu_time(time_2)
      if( hecMESH%my_rank .eq. 0 ) then
        write(ISTA,'(a,f10.2)') '         solve (sec) :', time_2 - time_1
      end if

    enddo
!C
!C-- end of time step loop

!C-- deallocate
    deallocate( myEIG%mass    ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, mass>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

!C-- file close for local use
      if( hecMESH%my_rank .eq. my_rank_monit_1 ) then
        CLOSE(dynamic_IW4)
        CLOSE(dynamic_IW5)
        CLOSE(dynamic_IW6)
      endif
!C-- end of finalization
!C
  end subroutine fstr_solve_dynamic_explicit


end module fstr_dynamic_explicit
