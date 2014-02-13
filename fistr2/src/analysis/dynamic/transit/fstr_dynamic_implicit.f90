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
!> \brief This module contains subroutines for implicit dynamic analysis

module fstr_dynamic_implicit

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

  subroutine fstr_solve_dynamic_implicit(hecMESH,hecMAT,fstrSOLID,myEIG   &
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

      type ( hecmwST_matrix      ) :: MAT0
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
!    real(kind=kreal) :: time_1, time_2

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


!C=============================C
!C-- implicit dynamic analysis
!C=============================C
!
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
	MAT0%NDOF=hecMESH%n_dof
	MAT0%N=hecMAT%N

    allocate( MAT0%D (nn*hecMAT%NP )        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, D>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    allocate( MAT0%AU(nn*hecMAT%NPU)        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, AU>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    allocate( MAT0%AL(nn*hecMAT%NPL)        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, AL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    allocate( MAT0%indexU(0:hecMAT%NP)      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, indexU>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    allocate( MAT0%indexL(0:hecMAT%NP)      ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, indexL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    allocate( MAT0%itemU(hecMAT%NPU)        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, itemU>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    allocate( MAT0%itemL(hecMAT%NPL)        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <fstr_solve_LINEAR_DYNAMIC, itemL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
 !   call hecmw_cmat_init( MAT0%cmat )

!C
!C-- check parameters
!C
         if(dabs(fstrDYNAMIC%beta) .lt. 1.0e-20) then
              if( hecMESH%my_rank .eq. 0 ) then
                 write(imsg,*) 'stop due to fstrDYNAMIC%beta = 0'
              end if
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
    !  call hecmw_restart_open()
    !  call hecmw_restart_read_int(restrt_step)
    !  call hecmw_restart_read_real(fstrDYNAMIC%DISP(:,1))
    !  call hecmw_restart_read_real(fstrDYNAMIC%VEL (:,1))
    !  call hecmw_restart_read_real(fstrDYNAMIC%ACC (:,1))
    !  call hecmw_restart_close()
      restrt_step_num = restrt_step(1) + 1
    end if
!C
!C-- matrix [KL] & [M]
!C-- matrix [KL]

      call fstr_mat_ass_main(fstrSOLID)

    do j = 1 ,nn*hecMAT%NP
      MAT0%D(j)   = hecMAT%D(j)
    end do
    do j = 1 ,nn*hecMAT%NPU
      MAT0%AU(j)  = hecMAT%AU(j)
    end do
    do j = 1 ,nn*hecMAT%NPL
      MAT0%AL(j)  = hecMAT%AL(j)
    end do
    do j = 0 ,hecMAT%NP
      MAT0%indexU(j) = hecMAT%indexU(j)
    end do
    do j = 0 ,hecMAT%NP
      MAT0%indexL(j) = hecMAT%indexL(j)
    end do
    do j = 1 ,hecMAT%NPU
      MAT0%itemU(j)  = hecMAT%itemU(j)
    end do
    do j = 1 ,hecMAT%NPL
      MAT0%itemL(j)  = hecMAT%itemL(j)
    end do
!C-- matrix [M]
!C-- lumped mass matrix
             if(fstrDYNAMIC%idx_mas .eq. 1) then

                call setMASS(IDBG,hecMESH,hecMAT,myEIG)

!C-- consistent mass matrix
               else if(fstrDYNAMIC%idx_mas .eq. 2) then
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
    a1 = .5/fstrDYNAMIC%beta - 1.
    a2 = 1./(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)
    a3 = 1./(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta**2)
    b1 = ( .5*fstrDYNAMIC%ganma/fstrDYNAMIC%beta - 1. )*fstrDYNAMIC%t_delta
    b2 = fstrDYNAMIC%ganma/fstrDYNAMIC%beta - 1.
    b3 = fstrDYNAMIC%ganma/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)

    c1 = 1. + fstrDYNAMIC%ray_k*fstrDYNAMIC%ganma/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)
    c2 = 1./(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta**2) + &
        fstrDYNAMIC%ray_m*fstrDYNAMIC%ganma/(fstrDYNAMIC%beta*fstrDYNAMIC%t_delta)

!!
!!    step = 0
!!
!C-- output initial condition
    if( restrt_step_num .eq. 1 ) then

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
!        time_1 =  hecmw_Wtime()
!
       fstrDYNAMIC%i_step = i

!C-- matrix [A]
    do j = 1 ,nn*hecMAT%NP
      hecMAT%D(j)  = c1*MAT0%D(j)
    end do
    do j = 1 ,nn*hecMAT%NPU
      hecMAT%AU(j) = c1*MAT0%AU(j)
    end do
    do j = 1 ,nn*hecMAT%NPL
      hecMAT%AL(j) = c1*MAT0%AL(j)
    end do
!
    do j=1,nnod
      do kk=1,ndof
        idm = nn*(j-1)+1 + (ndof+1)*(kk-1)
        imm = ndof*(j-1) + kk
        hecMAT%D(idm) = hecMAT%D(idm) + c2*myEIG%mass(imm)
      end do
    end do
!C
!C-- {vec1}={vec m1}, {vec2}={vec m2}
!C--
    do j = 1 ,ndof*nnod
      fstrDYNAMIC%VEC1(j) = a1*fstrDYNAMIC%ACC(j,1) + a2*fstrDYNAMIC%VEL(j,1) + &
                            a3*fstrDYNAMIC%DISP(j,1)
      fstrDYNAMIC%VEC2(j) = b1*fstrDYNAMIC%ACC(j,1) + b2*fstrDYNAMIC%VEL(j,1) + &
                            b3*fstrDYNAMIC%DISP(j,1)
    end do
!C-- {vec3}={vec k}=[KL]{vec2}

!       call matvec(fstrDYNAMIC%VEC3,fstrDYNAMIC%VEC2,hecMAT,ndof,D,AU,AL)
      if( hecMESH%n_dof .eq. 3 ) then
    !    call hecmw_matvec_33 (hecMESH, MAT0, fstrDYNAMIC%VEC2, fstrDYNAMIC%VEC3)
      else if( hecMESH%n_dof .eq. 2 ) then
    !    call hecmw_matvec_22 (hecMESH, MAT0, fstrDYNAMIC%VEC2, fstrDYNAMIC%VEC3, nnod)
      else if( hecMESH%n_dof .eq. 6 ) then
!!      call hecmw_matvec_11 (hecMESH, MAT0, fstrDYNAMIC%VEC2, fstrDYNAMIC%VEC3, n_nod_dof)
        call matvec(fstrDYNAMIC%VEC3, fstrDYNAMIC%VEC2, hecMAT, ndof, MAT0%D, MAT0%AU, MAT0%AL)
      end if
!C
!C-- mechanical boundary condition

        call dynamic_mat_ass_load (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)


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
!C
    do j = 1 ,ndof*nnod
      hecMAT%B(j) = hecMAT%B(j) + myEIG%mass(j)*( fstrDYNAMIC%VEC1(j) + fstrDYNAMIC%ray_m*  &
                    fstrDYNAMIC%VEC2(j) ) + fstrDYNAMIC%ray_k*fstrDYNAMIC%VEC3(j)
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
        if (.not. ds) then !In case of Direct Solver prevent MPI
      !    call hecmw_allREDUCE_R1( hecMESH,bsize,hecmw_sum )
        end if

        iexit = 0
        if( bsize < 1.0e-31 ) then
          iexit = 1
          if( hecMESH%my_rank .eq. 0 ) then
            WRITE(IMSG,*) '###Load Vector Error!'
          endif
!!!       call hecmw_abort( hecmw_comm_get_comm())
        endif

!C
!C-- linear solver [A]{X} = {B}
!C
        if( iexit .eq. 1 ) then
          hecMAT%X = 0.0
        else
!          CALL solve_LINEQ(hecMESH,hecMAT,imsg)
        end if
!C
!C-- new displacement, velocity and accelaration
!C
    do j = 1 ,ndof*nnod
      fstrDYNAMIC%DISP(j,2) = hecMAT%X(j)
      fstrDYNAMIC%ACC (j,2) = -a1*fstrDYNAMIC%ACC(j,1) - a2*fstrDYNAMIC%VEL(j,1) + &
                               a3*( fstrDYNAMIC%DISP(j,2) - fstrDYNAMIC%DISP(j,1) )
      fstrDYNAMIC%VEL (j,2) = -b1*fstrDYNAMIC%ACC(j,1) - b2*fstrDYNAMIC%VEL(j,1) + &
                               b3*( fstrDYNAMIC%DISP(j,2) - fstrDYNAMIC%DISP(j,1) )

      fstrDYNAMIC%DISP(j,1) = fstrDYNAMIC%DISP(j,2)
      fstrDYNAMIC%ACC (j,1) = fstrDYNAMIC%ACC (j,2)
      fstrDYNAMIC%VEL (j,1) = fstrDYNAMIC%VEL (j,2)
    end do

!!! restart  !!!
    if(fstrDYNAMIC%restart_nout .lt. 0) then
      fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
    end if
    if( mod(i,fstrDYNAMIC%restart_nout) == 0 ) then
      restrt_step(1) = i
  !    call hecmw_restart_add_int(restrt_step,size(restrt_step))
  !    call hecmw_restart_add_real(fstrDYNAMIC%DISP(:,1),size(fstrDYNAMIC%DISP(:,1)))
  !    call hecmw_restart_add_real(fstrDYNAMIC%VEL (:,1),size(fstrDYNAMIC%VEL (:,1)))
  !    call hecmw_restart_add_real(fstrDYNAMIC%ACC (:,1),size(fstrDYNAMIC%ACC (:,1)))
  !    call hecmw_restart_write()
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
!      time_2 =  hecmw_Wtime()
!     if( hecMESH%my_rank .eq. 0 ) then
!       write(ISTA,'(a,f10.2)') '         solve (sec) :', time_2 - time_1
!     end if
!

!C *****************************************************
!C for couple analysis
!      if( fstrPARAM%fg_couple ==1 ) then
!          if( fstrDYNAMIC%i_step < fstrDYNAMIC%n_step .or. &
!                (fstrDYNAMIC%i_step==fstrDYNAMIC%n_step .and. fstrPARAM%fg_couple_first==2)) then
!                  do j=1, fstrCPL%coupled_node_n
!                        if( fstrCPL%coupled_node(j) > hecMESH%n_node ) then
!                           write(IDBG,*) "fstr_solve_linear_dynamic: warning: ", &
!                                & "data in overlap region not set"
!                           cycle
!                        endif
!                        if( fstrCPL%dof == 3 ) then
!                              kkk0 = j*3;
!                              kkk1 = fstrCPL%coupled_node(j)*3
!
!                              fstrCPL%disp (kkk0-2) = fstrDYNAMIC%DISP(kkk1-2,1)
!                              fstrCPL%disp (kkk0-1) = fstrDYNAMIC%DISP(kkk1-1,1)
!                              fstrCPL%disp (kkk0  ) = fstrDYNAMIC%DISP(kkk1  ,1)
!
!                              fstrCPL%velo (kkk0-2) = fstrDYNAMIC%VEL (kkk1-2,1)
!                              fstrCPL%velo (kkk0-1) = fstrDYNAMIC%VEL (kkk1-1,1)
!                              fstrCPL%velo (kkk0  ) = fstrDYNAMIC%VEL (kkk1  ,1)
!
!                              fstrCPL%accel(kkk0-2) = fstrDYNAMIC%ACC (kkk1-2,1)
!                              fstrCPL%accel(kkk0-1) = fstrDYNAMIC%ACC (kkk1-1,1)
!                              fstrCPL%accel(kkk0  ) = fstrDYNAMIC%ACC (kkk1  ,1)
!                        else
!                              kkk0 = j*2;
!                              kkk1 = fstrCPL%coupled_node(j)*2
!
!                              fstrCPL%disp (kkk0-1) = fstrDYNAMIC%DISP(kkk1-1,1)
!                              fstrCPL%disp (kkk0  ) = fstrDYNAMIC%DISP(kkk1  ,1)
!
!                              fstrCPL%velo (kkk0-1) = fstrDYNAMIC%VEL (kkk1-1,1)
!                              fstrCPL%velo (kkk0  ) = fstrDYNAMIC%VEL (kkk1  ,1)
!
!                              fstrCPL%accel(kkk0-1) = fstrDYNAMIC%ACC (kkk1-1,1)
!                              fstrCPL%accel(kkk0  ) = fstrDYNAMIC%ACC (kkk1  ,1)
!                        endif
!                  end do
!                  call fstr_rcap_send( fstrCPL )
!            endif
!      endif
!C *****************************************************


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

    deallocate( MAT0%D             ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, D>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    deallocate( MAT0%AU            ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, AU>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    deallocate( MAT0%AL            ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, AL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    deallocate( MAT0%indexU        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, indexU>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    deallocate( MAT0%indexL        ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, indexL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    deallocate( MAT0%itemU         ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, itemU>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
    deallocate( MAT0%itemL         ,STAT=ierror )
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <fstr_solve_LINEAR_DYNAMIC, itemL>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if


!C-- end of implicit dynamic analysis
!C

!C-- file close for local use
      if( hecMESH%my_rank .eq. my_rank_monit_1 ) then
        CLOSE(dynamic_IW4)
        CLOSE(dynamic_IW5)
        CLOSE(dynamic_IW6)
      endif
!C-- end of finalization
!C
  end subroutine fstr_solve_dynamic_implicit


end module fstr_dynamic_implicit
