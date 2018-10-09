!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provide a function to control nonsteady heat analysis
module m_heat_solve_TRAN
contains

  subroutine heat_solve_TRAN ( hecMESH,hecMAT,fstrRESULT,fstrPARAM,fstrHEAT,ISTEP )
    use m_fstr
    use m_heat_mat_ass_conductivity
    use m_heat_mat_ass_capacity
    use m_heat_mat_ass_boundary
    use m_heat_init
    use m_hecmw2fstr_mesh_conv
    use m_solve_lineq
    implicit none
    integer(kind=kint) :: ISTEP, ITM, incr, iend, iterALL, i, inod, mnod, max_step, tstep, interval, bup_n_dof
    real(kind=kreal)   :: CTIME, ST, DTIME, EETIME, DELMAX, DELMIN, TT, BETA, val, CHK, tmpmax, dltmp, tmpmax_myrank
    type(hecmwST_local_mesh)  :: hecMESH
    type(hecmwST_matrix)      :: hecMAT
    type(hecmwST_result_data) :: fstrRESULT
    type(fstr_param)          :: fstrPARAM
    type(fstr_heat)           :: fstrHEAT
    type(hecmwST_matrix), pointer :: hecMATmpc
    integer(kind=kint) :: restart_step(1)
    real(kind=kreal)   :: restart_time(1)
    integer(kind=kint) :: restrt_data_size
    integer, parameter :: miniter = 4
    character(len=HECMW_HEADER_LEN) :: header
    character(len=HECMW_NAME_LEN)   :: label
    character(len=HECMW_NAME_LEN)   :: nameID

    call hecmw_mpc_mat_init(hecMESH, hecMAT, hecMATmpc)

    if( ISTEP .eq. 1 ) then
      ST = 0.0d0
    else
      ST = 0.0d0
      do i = 1, ISTEP - 1
        ST = ST + fstrHEAT%STEP_EETIME(i)
      enddo
    endif

    DTIME  = fstrHEAT%STEP_DLTIME(ISTEP)
    EETIME = fstrHEAT%STEP_EETIME(ISTEP)
    DELMIN = fstrHEAT%STEP_DELMIN(ISTEP)
    DELMAX = fstrHEAT%STEP_DELMAX(ISTEP)
    ITM    = fstrPARAM%itmax(ISTEP)
    EPS    = fstrPARAM%eps(ISTEP)
    TT     = 0.0d0
    write(*,*) ' DTIME=',DTIME
    write(*,*) 'EETIME=',EETIME
    write(*,*) 'DELMIN=',DELMIN
    write(*,*) 'DELMAX=',DELMAX
    write(*,*) '   ITM=',ITM
    write(*,*) '   EPS=',EPS

    max_step = ( EETIME + 0.1d0 * DTIME ) / DTIME
    BETA = 0.5d0
    incr = 0
    iend = 0
    tstep = 0
    hecMAT%NDOF = 1
    hecMAT%Iarray(98) = 1   !Assmebly complete
    hecMAT%X = 0.0d0

    !C Restart read
    if( fstrHEAT%restart_nout < 0 ) then
      fstrHEAT%restart_nout = -fstrHEAT%restart_nout
      call hecmw_restart_open()
      call hecmw_restart_read_int(restart_step)
      call hecmw_restart_read_real(restart_time)
      call hecmw_restart_read_real(fstrHEAT%TEMP0)
      call hecmw_restart_close()
      tstep = restart_step(1)
      TT = restart_time(1)

      do i= 1, hecMESH%n_node
        fstrHEAT%TEMPC(i)= fstrHEAT%TEMP0(i)
        fstrHEAT%TEMP (i)= fstrHEAT%TEMP0(i)
      enddo
      write(ILOG,*) ' Restart read of temperatures: OK'
    endif

    !C--------------------   START TRANSIET LOOP   ------------------------
    tr_loop: do
      incr = incr + 1
      if( TT+DTIME*1.0000001d0 >= EETIME ) then
        DTIME = EETIME - TT
        TT = EETIME
        iend = 1
      else
        TT = TT + DTIME
        iend = 0
      endif
      CTIME = TT

      if( hecMESH%my_rank.eq.0 ) then
        write(IMSG,*) '// INCREMENT NO. =', incr, TT, DTIME
        write(*,*) '// INCREMENT NO. =', incr, TT, DTIME
        write(IDBG,*) '// INCREMENT NO. =',incr, TT, DTIME
      endif

      if( DTIME .lt. DELMIN ) then
        if( hecMESH%my_rank.eq.0 ) then
          write(IMSG,*) ' !!! DELTA TIME EXCEEDED TOLERANCE OF TIME INCREMENT'
        endif
        call hecmw_abort( hecmw_comm_get_comm() )
      endif

      iterALL = 0
      !C==============  START OF ITERATION LOOP  ===========
      do
        iterALL= iterALL + 1

        !C-- MATRIX ASSEMBLING -----
        call heat_mat_ass_conductivity ( hecMESH,hecMAT,fstrHEAT,BETA )
        call heat_mat_ass_capacity ( hecMESH,hecMAT,fstrHEAT,DTIME )
        call heat_mat_ass_boundary ( hecMESH,hecMAT,hecMATmpc,fstrHEAT,TT,ST, DTIME )

        !C-- SOLVER
        hecMATmpc%Iarray(97) = 1   !Need numerical factorization
        bup_n_dof = hecMESH%n_dof
        hecMESH%n_dof = 1
        call solve_LINEQ(hecMESH,hecMATmpc)
        hecMESH%n_dof=bup_n_dof
        call hecmw_mpc_tback_sol(hecMESH, hecMAT, hecMATmpc)

        !C-- UPDATE -----
        do i= 1, hecMESH%n_node
          fstrHEAT%TEMPC(i) = fstrHEAT%TEMP(i)
          fstrHEAT%TEMP (i) = hecMAT%X(i)
        enddo

        !C-- GLOBAL RESIDUAL -----
        val= 0.0d0
        do i= 1, hecMESH%nn_internal
          val= val + (fstrHEAT%TEMP(i)-fstrHEAT%TEMPC(i))**2
        enddo
        call hecmw_allREDUCE_R1 ( hecMESH,val,hecmw_sum )
        !C
        CHK= dsqrt(val)
        if ( hecMESH%my_rank.eq.0 ) then
          write(*,'(i8,1pe16.6)') iterALL, CHK
          write(IMSG,'(i8,1pe16.6)') iterALL, CHK
          call flush(IMSG)
        endif
        !C
        if ( CHK.lt.EPS ) then
          if ( hecMESH%my_rank==0 ) then
            write(*,*) ' !!! CONVERGENCE ACHIEVED '
            write(IMSG,*) ' !!! CONVERGENCE ACHIEVED '
            call flush(IMSG)
          endif
          exit
        endif
        !C
        if ( iterALL>ITM ) then
          if ( hecMESH%my_rank==0 ) then
            write(*,*) ' ** ITERATION COUNT OVER : MAX = ', ITM
            write(IMSG,*) ' *** ITERATION COUNT OVER : MAX = ', ITM
            call flush(IMSG)
          endif

          if( DELMIN .gt. 0.d0) then
            TT = TT - DTIME
            DTIME = 0.5*DTIME
            cycle tr_loop
          else
            call hecmw_abort( hecmw_comm_get_comm() )
          endif
        endif
      enddo
      !C==============  END OF ITERATION LOOP  ===========

      if( DELMIN .gt. 0.d0 ) then
        tmpmax = 0.d0
        do i= 1, hecMESH%nn_internal
          inod=fstrPARAM%global_local_id(1,i)
          dltmp  = fstrHEAT%TEMP0(i) - fstrHEAT%TEMP(i)
          if( dabs(dltmp).gt.tmpmax ) then
            mnod = inod
            tmpmax = dabs(dltmp)
          endif
        enddo

        tmpmax_myrank = tmpmax
        call hecmw_allREDUCE_R1(hecMESH, tmpmax, hecmw_max)

        if (tmpmax .gt. tmpmax_myrank) mnod = -1
        call hecmw_allREDUCE_I1(hecMESH, mnod, hecmw_max)

        if ( tmpmax .gt. DELMAX ) then
          if( hecMESH%my_rank.eq.0 ) then
            write(IMSG,*) ' *** EXCEEDED TOLERANCE OF VARIATION IN TEMPERATUTE.'
            write(IMSG,*) ' : NODE NUMBER  = ', mnod, ' : DELTA TEMP   = ', tmpmax
            write(*,*) ' *** EXCEEDED TOLERANCE OF VARIATION IN TEMPERATUTE.'
            write(*,*) ' : NODE NUMBER  = ', mnod, ' : DELTA TEMP   = ', tmpmax
            call flush(IMSG)
          endif
          TT    = TT - DTIME
          DTIME = 0.5*DTIME
          cycle
        endif

        if ( iterALL <=  miniter ) DTIME = DTIME * 1.5d0
      endif

      do i= 1, hecMESH%n_node
        fstrHEAT%TEMP0(i) = fstrHEAT%TEMP(i)
      enddo

      !C
      !C=== OUTPUT
      !C
      tstep = tstep+1
      write(ILOG,*)
      write(ILOG,'(a,i6, a,f10.3)') ' STEP =', tstep, ' Time  =', CTIME

      if( IRESULT.eq.1 .and. (mod(tstep,IRRES).eq.0 .or. tstep.eq.max_step) ) then
        header = '*fstrresult'
        call hecmw_result_init(hecMESH,max_step,tstep,header)
        label = 'TEMPERATURE'
        call hecmw_result_add(1,1,label,fstrHEAT%TEMP)
        nameID = 'fstrRES'
        call hecmw_result_write_by_name(nameID)
        call hecmw_result_finalize
      endif

      if( IVISUAL.eq.1 .and. (mod(tstep,IWRES).eq.0 .or. tstep.eq.max_step) ) then
        interval = IWRES
        call hecmw_nullify_result_data(fstrRESULT)
        fstrRESULT%nn_component = 1
        fstrRESULT%ne_component = 0
        allocate(fstrRESULT%nn_dof(1))
        allocate(fstrRESULT%node_label(1))
        allocate(fstrRESULT%node_val_item(hecMESH%n_node))
        fstrRESULT%nn_dof(1) = 1
        fstrRESULT%node_label(1) = 'TEMPERATURE'
        fstrRESULT%node_val_item = fstrHEAT%TEMP
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize (hecMESH,fstrRESULT,tstep,max_step,interval)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif

      !C Restart write
      if( fstrHEAT%restart_nout > 0 .and. (mod(tstep,fstrHEAT%restart_nout).eq.0 .or. tstep.eq.max_step) ) then
        restart_step(1) = tstep
        restart_time(1) = CTIME
        restrt_data_size = size(restart_step)
        call hecmw_restart_add_int(restart_step,restrt_data_size)
        restrt_data_size = size(restart_time)
        call hecmw_restart_add_real(restart_time,restrt_data_size)
        restrt_data_size = size(fstrHEAT%TEMP)
        call hecmw_restart_add_real(fstrHEAT%TEMP,restrt_data_size)
        call hecmw_restart_write()
        if( hecMESH%my_rank.eq.0 ) then
          write(IMSG,*) '### FSTR output Restart_File.'
          call flush(IMSG)
        endif
      endif

      if( iend.ne.0 ) exit
    enddo tr_loop
    !C--------------------   END TRANSIET LOOP   ------------------------

    call hecmw_mpc_mat_finalize(hecMESH, hecMAT, hecMATmpc)

  end subroutine heat_solve_TRAN
end module m_heat_solve_TRAN
