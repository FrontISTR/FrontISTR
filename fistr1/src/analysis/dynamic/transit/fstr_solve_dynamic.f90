!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains subroutines controlling dynamic calculation

module fstr_solver_dynamic

  use m_fstr
  use fstr_dynamic_nlexplicit
  use fstr_dynamic_nlimplicit
  use fstr_frequency_analysis  !Frequency analysis module

contains

  !C================================================================C
  !> Master subroutine for dynamic analysis
  !C================================================================C
  subroutine fstr_solve_DYNAMIC(hecMESH,hecMAT,fstrSOLID,fstrEIG   &
      ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
      ,fstrCPL,fstrFREQ, hecLagMAT  &
      ,conMAT )
    use m_fstr_setup
    implicit none
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_eigen)                     :: fstrEIG
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_result_data)            :: fstrRESULT
    type(fstr_param)                     :: fstrPARAM
    type(fstr_dynamic)                   :: fstrDYNAMIC
    type(fstr_couple)                    :: fstrCPL !for COUPLE
    type(fstr_freqanalysis)              :: fstrFREQ
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    type(fstr_info_contactChange)        :: infoCTChange !< type fstr_info_contactChange
    type(hecmwST_matrix)                 :: conMAT
    integer(kind=kint) :: i, j, num_monit, ig, is, iE, ik, in, ing, iunitS, iunit, ierror, flag, limit
    character(len=HECMW_FILENAME_LEN) :: fname, header
    integer(kind=kint) :: restrt_step_num, ndof
    integer(kind=kint) :: restrt_step(1)

    num_monit = 0

    if(dabs(fstrDYNAMIC%t_delta) < 1.0e-20) then
      if( hecMESH%my_rank == 0 ) then
        write(imsg,*) 'stop due to fstrDYNAMIC%t_delta = 0'
      end if
      call hecmw_abort( hecmw_comm_get_comm())
    end if
    call fstr_dynamic_alloc( hecMESH, fstrDYNAMIC )

    !C-- file open for local use
    !C
    if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
      call hecmw_ctrl_is_subdir( flag, limit )
      if( flag == 0)then
        header = ""
      else
        header = adjustl("MONITOR/")
        call hecmw_ctrl_make_subdir( adjustl("MONITOR/test.txt"), limit )
      endif
      ig = fstrDYNAMIC%ngrp_monit
      is = hecMESH%node_group%grp_index(ig-1)+1
      iE = hecMESH%node_group%grp_index(ig)
      do ik=is,iE
        in = hecMESH%node_group%grp_item(ik)
        if (in > hecMESH%nn_internal) cycle
        num_monit = num_monit+1
        ing = hecMESH%global_node_id(in)
        iunitS = 10*(num_monit-1)

        iunit = iunitS + fstrDYNAMIC%dynamic_IW4
        write(fname,'(a,i0,a)') trim(header)//'dyna_disp_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          open(iunit,file=fname, position = 'append', iostat=ierror)
        else
          open(iunit,file=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error:',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW5
        write(fname,'(a,i0,a)') trim(header)//'dyna_velo_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          open(iunit,file=fname, position = 'append', iostat=ierror)
        else
          open(iunit,file=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW6
        write(fname,'(a,i0,a)') trim(header)//'dyna_acce_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          open(iunit,file=fname, position = 'append', iostat=ierror)
        else
          open(iunit,file=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW7
        write(fname,'(a,i0,a)') trim(header)//'dyna_force_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          open(iunit,file=fname, position = 'append', iostat=ierror)
        else
          open(iunit,file=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if
        iunit = iunitS + fstrDYNAMIC%dynamic_IW8
        write(fname,'(a,i0,a)') trim(header)//'dyna_strain_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          open(iunit,file=fname, position = 'append', iostat=ierror)
        else
          open(iunit,file=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if

        iunit = iunitS + fstrDYNAMIC%dynamic_IW9
        write(fname,'(a,i0,a)') trim(header)//'dyna_stress_',ing,'.txt'
        if(fstrDYNAMIC%restart_nout < 0 ) then
          open(iunit,file=fname, position = 'append', iostat=ierror)
        else
          open(iunit,file=fname, status = 'replace', iostat=ierror)
        endif
        if( ierror /= 0 ) then
          write(*,*) 'stop due to file opening error',trim(fname)
          call hecmw_abort( hecmw_comm_get_comm())
        end if
      enddo
    endif

    !C
    !C-- initial condition
    !C
    fstrDYNAMIC%DISP = 0.d0
    fstrDYNAMIC%VEL  = 0.d0
    fstrDYNAMIC%ACC  = 0.d0
    fstrSOLID%unode(:)  =0.d0
    fstrSOLID%QFORCE(:) =0.d0
    fstrDYNAMIC%kineticEnergy=0.d0
    fstrDYNAMIC%strainEnergy=0.d0
    fstrDYNAMIC%totalEnergy=0.d0
    call fstr_UpdateState( hecMESH, fstrSOLID, 0.d0 )

    ! ---- restart

    restrt_step_num = 1
    fstrDYNAMIC%i_step = 0
    infoCTChange%contactNode_previous = 0

    if(associated(g_InitialCnd))then
      ndof = HECMAT%NDOF
      do j = 1, size(g_InitialCnd)
        if(g_InitialCnd(j)%cond_name == "velocity")then
          do i= 1, hecMESH%n_node
            ing = g_InitialCnd(j)%intval(i)
            if(ing <= 0) cycle
            fstrDYNAMIC%VEL(ndof*i-(ndof-ing),1) = g_InitialCnd(j)%realval(i)
          end do
        elseif(g_InitialCnd(j)%cond_name == "acceleration")then
          do i = 1, hecMESH%n_node
            ing = g_InitialCnd(j)%intval(i)
            if(ing <= 0) cycle
            fstrDYNAMIC%ACC(ndof*i-(ndof-ing),1) = g_InitialCnd(j)%realval(i)
          enddo
        endif
      enddo
    endif

    if(fstrDYNAMIC%restart_nout >= 0 ) then
      call dynamic_bc_init   (hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
      call dynamic_bc_init_vl(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
      call dynamic_bc_init_ac(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC)
    endif

    !restart
    if(fstrDYNAMIC%restart_nout < 0 ) then
      if( fstrDYNAMIC%idx_eqa == 1 ) then
        call fstr_read_restart_dyna_nl(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
          infoCTChange%contactNode_previous)
      elseif(fstrDYNAMIC%idx_eqa == 11) then
        if( .not. associated( fstrSOLID%contacts ) ) then
          call fstr_read_restart_dyna_nl(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM)
        else
          call fstr_read_restart_dyna_nl(restrt_step_num,hecMESH,fstrSOLID,fstrDYNAMIC,fstrPARAM,&
            infoCTChange%contactNode_previous)
        endif
      endif
      restrt_step_num = restrt_step_num + 1
      fstrDYNAMIC%restart_nout = - fstrDYNAMIC%restart_nout
      hecMAT%Iarray(98) = 1
    end if

    if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis

      if(fstrDYNAMIC%idx_eqa == 1) then     ! implicit dynamic analysis
        call fstr_solve_dynamic_nlimplicit_contactSLag(1, hecMESH,hecMAT,fstrSOLID,fstrEIG   &
          ,fstrDYNAMIC,fstrRESULT,fstrPARAM &
          ,fstrCPL,hecLagMAT,restrt_step_num,infoCTChange   &
          ,conMAT )

      else if(fstrDYNAMIC%idx_eqa == 11) then  ! explicit dynamic analysis
        call fstr_solve_dynamic_nlexplicit(hecMESH,hecMAT,fstrSOLID,fstrEIG   &
          ,fstrDYNAMIC,fstrRESULT,fstrPARAM,infoCTChange &
          ,fstrCPL, restrt_step_num )
      endif

    else if(fstrDYNAMIC%idx_resp == 2) then  ! frequency response analysis

      if( fstrPARAM%nlgeom ) then
        if( hecMESH%my_rank .eq. 0 ) then
          write(imsg,*) 'stop: steady-state harmonic response analysis is not yet available !'
        end if
        call hecmw_abort( hecmw_comm_get_comm())
      end if

      if( hecMESH%my_rank .eq. 0 ) then
        call fstr_solve_frequency_analysis(hecMESH, hecMAT, fstrSOLID, fstrEIG, fstrDYNAMIC, &
          fstrRESULT, fstrPARAM, fstrCPL, fstrFREQ, hecLagMAT, &
          restrt_step_num)
      end if
    end if

    !C-- file close for local use
    if(fstrDYNAMIC%idx_resp == 1) then   ! time history analysis
      do i=1,num_monit
        iunitS = 10*(i-1)
        close(iunitS + fstrDYNAMIC%dynamic_IW4)
        close(iunitS + fstrDYNAMIC%dynamic_IW5)
        close(iunitS + fstrDYNAMIC%dynamic_IW6)
        close(iunitS + fstrDYNAMIC%dynamic_IW7)
        close(iunitS + fstrDYNAMIC%dynamic_IW8)
        close(iunitS + fstrDYNAMIC%dynamic_IW9)
      enddo
    endif
    !C-- end of finalization

    call fstr_dynamic_finalize( fstrDYNAMIC )
    call hecMAT_finalize( hecMAT )

    deallocate( fstrEIG%mass )

  end subroutine fstr_solve_DYNAMIC

end module fstr_solver_dynamic
