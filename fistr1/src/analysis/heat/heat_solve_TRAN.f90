!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : Heat Analysis                                     !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This module provide a function to control nonsteady heat analysis
module m_heat_solve_TRAN
   contains
!C
!C** TRANSIENT
!C
   subroutine heat_solve_TRAN ( hecMESH,hecMAT,fstrRESULT,fstrPARAM,fstrHEAT,ISTEP,CTIME, maxincr )

      use m_fstr
      use m_heat_mat_ass_conductivity
      use m_heat_mat_ass_capacity
      use m_heat_mat_ass_boundary
      use m_heat_init
      use m_heat_make_result
      use m_hecmw2fstr_mesh_conv

      implicit none
      integer(kind=kint) ISTEP,ITM,incr,iend,iterALL,i,inod,mnod,nd,id,ndof,maxincr,tstep
      real(kind=kreal)   CTIME,ST,DTIME,EETIME,DELMAX,DELMIN,TT,BETA,VAL,CHK,tmpmax,dltmp,tmpmax_myrank
!C file name
      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN) :: label
      character(len=HECMW_NAME_LEN) :: nameID

      type (hecmwST_local_mesh ) :: hecMESH
      type (hecmwST_matrix     ) :: hecMAT
      type (hecmwST_result_data) :: fstrRESULT
      type (fstr_param         ) :: fstrPARAM
      type (fstr_heat          ) :: fstrHEAT

      integer(kind=kint) :: restart_step(1)
      real(kind=kreal)   :: restart_time(1)
      integer(kind=kint) :: restrt_data_size

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
      TT     = CTIME
      if( NPRINT.eq.0 ) NPRINT = 1
      write(*,*) ' DTIME=',DTIME
      write(*,*) 'EETIME=',EETIME
      write(*,*) 'DELMIN=',DELMIN
      write(*,*) 'DELMAX=',DELMAX
      write(*,*) '   ITM=',ITM
      write(*,*) '   EPS=',EPS
      write(*,*) '    TT=',TT
      write(*,*) 'NPRINT=',NPRINT

      BETA = 0.5
      incr = 0
      iend = 0
      tstep = 0
      hecMAT%NDOF = 1
      hecMAT%Iarray(98) = 1   !Assmebly complete

!C Restart read
      if( fstrPARAM%fg_irres .eq. kYES ) then
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
!C--------------------
        incr = incr + 1
        if( TT+DTIME >= EETIME ) then
          DTIME = EETIME - TT
          TT = EETIME
          iend = 1
        else
          TT = TT + DTIME
          iend = 0
        endif
        CTIME = TT

        if( hecMESH%my_rank.eq.0 ) then
            write(IMSG,*)
            write(IMSG,*) '// INCREMENT NO. =',incr, TT, DTIME
            write(*,*)
            write(*,*) '// INCREMENT NO. =',incr, TT, DTIME
        endif

        if( DTIME .lt.DELMIN ) then
          if( hecMESH%my_rank.eq.0 ) then
            write(IMSG,*) ' !!! DELTA TIME EXCEEDED TOLERANCE OF TIME INCREMENT'
            call flush(IMSG)
          endif
          call hecmw_abort( hecmw_comm_get_comm() )
        endif

        if( incr > maxincr ) then
          if( hecMESH%my_rank==0 ) then
            write(IMSG,*) ' !!! NUMBER OF INCREMENTS EXCEEDED MAXIMUM INCREMENTS'
            call flush(IMSG)
          endif
          call hecmw_abort( hecmw_comm_get_comm() )
        endif

        iterALL= 0
!C==============  START OF ITERATION LOOP  ===========
        do
!C==============
          iterALL= iterALL + 1
!C
!C-- MATRIX ASSEMBLING -----

          call heat_mat_ass_conductivity ( hecMESH,hecMAT,fstrHEAT,BETA )
              write(IDBG,*) 'mat_ass_conductivity: OK'
              call flush(IDBG)

          call heat_mat_ass_capacity ( hecMESH,hecMAT,fstrHEAT,DTIME )
              write(IDBG,*) 'mat_ass_capacity : OK'
              call flush(IDBG)

          call heat_mat_ass_boundary ( hecMESH,hecMAT,fstrHEAT,TT,ST, DTIME )
              write(IDBG,*) 'mat_ass_boundary: OK'
              call flush(IDBG)
			  
!C
!C-- SOLVER -----

          if( hecMAT%Iarray(99).eq. 1 ) then
            call hecmw_solve_11 ( hecMESH,hecMAT )
            write(IDBG,*) 'hecmw_solve_11: OK'
            call flush(IDBG)
          else
            hecMAT%Iarray(97) = 1   !Need numerical factorization
            call hecmw_solve_direct ( hecMESH,hecMAT,IMSG )
             do i=1,hecMAT%NP
                 hecMAT%X(i) = hecMAT%B(i)
             end do
            write(IDBG,*) 'hecmw_solve_direct: OK'
            call flush(IDBG)
          endif
!C
!C-- UPDATE -----

          do i= 1, hecMESH%n_node
            fstrHEAT%TEMPC(i)= fstrHEAT%TEMP(i)
            fstrHEAT%TEMP (i)= hecMAT%X(i)
          enddo
		  
!C
!C-- GLOBAL RESIDUAL -----

          VAL= 0.d0
          do i= 1, hecMESH%nn_internal
            VAL= VAL + (fstrHEAT%TEMP(i)-fstrHEAT%TEMPC(i))**2
          enddo
          call hecmw_allREDUCE_R1 ( hecMESH,VAL,hecmw_sum )
!C
          CHK= dsqrt(VAL)
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
!C==============
        enddo
!C==============  END OF ITERATION LOOP  ===========
!C

        if( DELMIN .gt. 0.d0 ) then
          tmpmax = 0.d0
          do i= 1, hecMESH%nn_internal
            inod=fstrPARAM%global_local_id(1,i)
            dltmp  = fstrHEAT%TEMP0(i) - fstrHEAT%TEMP(i)
            if( DABS(dltmp).gt.tmpmax ) then
              mnod = inod
              tmpmax = DABS(dltmp)
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

          if ( iterALL .le. 2 ) DTIME = DTIME * 2.0
        endif

        do i= 1, hecMESH%n_node
          fstrHEAT%TEMP0(i) = fstrHEAT%TEMP(i)
        enddo

!C
!C=== OUTPUT
!C
        tstep = tstep+1
        if( MOD(tstep,NPRINT)==0 .or. iend==1 ) then
          write(ILOG,*)
          write(ILOG,'(a,i6, a,f10.3)')    ' STEP =',tstep, ' Time  =',CTIME
          write(ILOG,*) '     Node   Temperature   '
          write(ILOG,*) '--------------------------'
          do i= 1, hecMESH%nn_internal
            inod=fstrPARAM%global_local_id(1,i)
            write (ILOG,'(i8,f12.4)') inod, fstrHEAT%TEMP(i)
          enddo
          call flush(ILOG)

          if( IRESULT==1 ) then
            nd = tstep
            header = '*fstrresult'
            call hecmw_result_init( hecMESH,nd,header )
            id    = 1
            ndof  = 1
            label = 'TEMPERATURE'
            call hecmw_result_add(id,ndof,label,fstrHEAT%TEMP)
            nameID = 'fstrRES'
            call hecmw_result_write_by_name(nameID)
            call hecmw_result_finalize
            if( hecMESH%my_rank.eq.0 ) then
              write(IMSG,*) '### FSTR output Result_File.'
              call flush(IMSG)
            endif
          endif

          if( IVISUAL.eq.1 ) then
            call heat_init_result ( hecMESH, fstrRESULT )
            call heat_make_result ( hecMESH, fstrHEAT, fstrRESULT )
            call fstr2hecmw_mesh_conv(hecMESH)
            call hecmw_visualize_init
            call hecmw_visualize ( hecMESH,fstrRESULT,tstep,tstep,0 )
            call hecmw_visualize_finalize
            call hecmw2fstr_mesh_conv(hecMESH)
            call hecmw_result_free(fstrRESULT)
            if( hecMESH%my_rank.eq.0 ) then
              write(IMSG,*) '### FSTR output Visual_File.'
              call flush(IMSG)
            endif
          endif

!C Restart write
          if( fstrPARAM%fg_iwres .eq. kYES ) then
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

        endif

        if( iend.ne.0 ) exit
      enddo tr_loop
!C--------------------   START TRANSIET LOOP   ------------------------

   end subroutine heat_solve_TRAN
end module m_heat_solve_TRAN
