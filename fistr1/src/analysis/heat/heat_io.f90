!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides a function to control heat analysis
module m_heat_io

contains

  subroutine heat_input_restart(fstrHEAT, hecMESH, tstep, tt)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_heat)           :: fstrHEAT
    integer(kind=kint) :: restart_step(1)
    real(kind=kreal)   :: restart_time(1)
    integer(kind=kint) :: i, tstep
    real(kind=kreal)   :: tt

    if(fstrHEAT%restart_nout < 0)then
      fstrHEAT%restart_nout = -fstrHEAT%restart_nout
      call hecmw_restart_open()
      call hecmw_restart_read_int(restart_step)
      call hecmw_restart_read_real(restart_time)
      call hecmw_restart_read_real(fstrHEAT%TEMP0)
      call hecmw_restart_close()
      tstep = restart_step(1)
      tt = restart_time(1)

      do i = 1, hecMESH%n_node
        fstrHEAT%TEMPC(i)= fstrHEAT%TEMP0(i)
        fstrHEAT%TEMP (i)= fstrHEAT%TEMP0(i)
      enddo
      write(ILOG,*) ' Restart read of temperatures: OK'
    endif
  end subroutine heat_input_restart

  subroutine heat_output_log(hecMESH, fstrPARAM, fstrHEAT, tstep, ctime)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_heat)           :: fstrHEAT
    type(fstr_param)          :: fstrPARAM
    integer(kind=kint) :: i, in, inod, tstep, nmax, nmin
    real(kind=kreal)   :: temp, ctime, tmax, tmin

    tmax = -1.0d10
    tmin =  1.0d10
    nmax = -1
    nmin = -1

    write(ILOG,*)
    write(ILOG,'(a,i6)')    ' ISTEP =', tstep
    write(ILOG,'(a,f10.3)') ' Time  =', ctime

    do i = 1, hecMESH%nn_internal
      inod = fstrPARAM%global_local_id(1,i)
      in = fstrPARAM%global_local_id(2,i)
      temp = fstrHEAT%TEMP(in)
      if(tmax < temp)then
        tmax = temp
        nmax = inod
      endif
      if(temp < tmin)then
        tmin = temp
        nmin = inod
      endif
    enddo

    write(ILOG,'(a,f10.3,i10)') ' Maximum Temperature :', tmax
    write(ILOG,'(a,i10)')       ' Maximum Node No.    :', nmax
    write(ILOG,'(a,f10.3,i10)') ' Minimum Temperature :', tmin
    write(ILOG,'(a,i10)')       ' Minimum Node No.    :', nmin

    !global temperature
    call hecmw_allreduce_R1 (hecMESH, tmax, hecmw_max)
    call hecmw_allreduce_R1 (hecMESH, tmin, hecmw_min)
    if( myrank == 0 ) then
      write(ILOG,'(a,f10.3,i10)') ' Maximum Temperature(global) :', tmax
      write(ILOG,'(a,f10.3,i10)') ' Minimum Temperature(global) :', tmin
    end if


  end subroutine heat_output_log

  subroutine heat_output_result(hecMESH, fstrHEAT, tstep, ctime, outflag)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_heat)           :: fstrHEAT
    integer(kind=kint) :: restart_step(1)
    real(kind=kreal)   :: restart_time(1)
    integer(kind=kint) :: i, tstep
    real(kind=kreal)   :: ctime, work(1)
    logical, intent(in)       :: outflag     !< if true, result will be output regardless of istep
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_MSG_LEN)    :: comment
    character(len=HECMW_NAME_LEN)   :: label
    character(len=HECMW_NAME_LEN)   :: nameID

    if(IRESULT == 1 .and. (mod(tstep, IRRES) == 0 .or. outflag))then
      header = '*fstrresult'
      comment = 'nonsteady_heat_result'
      call hecmw_result_init(hecMESH, tstep, header, comment)
      work(1) = ctime
      label = 'TOTALTIME'
      call hecmw_result_add(3, 1, label, work)
      label = 'TEMPERATURE'
      call hecmw_result_add(1, 1, label, fstrHEAT%TEMP)
      nameID = 'fstrRES'
      call hecmw_result_write_by_name(nameID)
      call hecmw_result_finalize
    endif
  end subroutine heat_output_result

  subroutine heat_output_visual(hecMESH, fstrRESULT, fstrHEAT, tstep, ctime, outflag)
    use m_fstr
    use m_hecmw2fstr_mesh_conv
    implicit none
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_heat)           :: fstrHEAT
    type(hecmwST_result_data) :: fstrRESULT
    integer(kind=kint) :: i, tstep
    real(kind=kreal)   :: ctime
    logical, intent(in)       :: outflag     !< if true, result will be output regardless of istep

    if(IVISUAL == 1 .and. (mod(tstep, IWRES) == 0 .or. outflag))then
      call hecmw_nullify_result_data(fstrRESULT)
      fstrRESULT%ng_component = 1
      fstrRESULT%nn_component = 1
      fstrRESULT%ne_component = 0
      allocate(fstrRESULT%ng_dof(1))
      allocate(fstrRESULT%global_label(1))
      allocate(fstrRESULT%global_val_item(1))
      fstrRESULT%ng_dof(1) = 1
      fstrRESULT%global_label(1) = 'TOTALTIME'
      fstrRESULT%global_val_item(1) = ctime
      allocate(fstrRESULT%nn_dof(1))
      allocate(fstrRESULT%node_label(1))
      allocate(fstrRESULT%node_val_item(hecMESH%n_node))
      fstrRESULT%nn_dof(1) = 1
      fstrRESULT%node_label(1) = 'TEMPERATURE'
      fstrRESULT%node_val_item = fstrHEAT%TEMP
      call fstr2hecmw_mesh_conv(hecMESH)
      call hecmw_visualize_init
      call hecmw_visualize( hecMESH, fstrRESULT, tstep )
      call hecmw_visualize_finalize
      call hecmw2fstr_mesh_conv(hecMESH)
      call hecmw_result_free(fstrRESULT)
    endif
  end subroutine heat_output_visual

  subroutine heat_output_restart(hecMESH, fstrHEAT, tstep, outflag, current_time)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh)  :: hecMESH
    type(fstr_heat)           :: fstrHEAT
    integer(kind=kint) :: restart_step(1)
    real(kind=kreal)   :: restart_time(1)
    integer(kind=kint) :: restrt_data_size
    integer(kind=kint) :: tstep
    logical, intent(in)       :: outflag     !< if true, result will be output regardless of istep
    real(kind=kreal)   :: current_time

    if(0 < fstrHEAT%restart_nout .and. (mod(tstep, fstrHEAT%restart_nout) == 0 .or. outflag))then
      restart_step(1) = tstep
      restart_time(1) = current_time
      restrt_data_size = size(restart_step)
      call hecmw_restart_add_int(restart_step, restrt_data_size)
      restrt_data_size = size(restart_time)
      call hecmw_restart_add_real(restart_time, restrt_data_size)
      restrt_data_size = size(fstrHEAT%TEMP)
      call hecmw_restart_add_real(fstrHEAT%TEMP, restrt_data_size)
      call hecmw_restart_write()
      if( hecMESH%my_rank.eq.0 ) then
        write(IMSG,*) '### FSTR output Restart_File.'
        call flush(IMSG)
      endif
    endif
  end subroutine heat_output_restart
end module m_heat_io
