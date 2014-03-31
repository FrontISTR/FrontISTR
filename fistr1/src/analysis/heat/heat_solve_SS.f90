!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!> This module provides a function for stationary heat analysis
module m_heat_solve_SS
   contains

   subroutine heat_solve_SS ( hecMESH,hecMAT,fstrRESULT,fstrPARAM,fstrHEAT,ISTEP,CTIME )

      use m_fstr
      use m_heat_mat_ass_conductivity
      use m_heat_mat_ass_boundary
      use m_heat_init
      use m_hecmw2fstr_mesh_conv
      use m_solve_lineq

      implicit none
      integer(kind=kint) ISTEP,iterALL,ITM,i,INCR,LMAX,LMIN,inod,ii,bup_n_dof
      real(kind=kreal)   CTIME,BETA,STIME,VAL,CHK,TMAX,TMIN,temp

      type (hecmwST_local_mesh  ) :: hecMESH
      type (hecmwST_matrix      ) :: hecMAT
      type (hecmwST_result_data ) :: fstrRESULT
      type (fstr_param          ) :: fstrPARAM
      type (fstr_heat           ) :: fstrHEAT

      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN)   :: label
      character(len=HECMW_NAME_LEN)   :: nameID

      BETA    = 1.0d0
      STIME   = 0.0d0
      ETIME   = 0.0d0
      iterALL = 0

      ITM = fstrPARAM%itmax(ISTEP)
      EPS = fstrPARAM%eps(ISTEP)

      hecMAT%NDOF = 1
      hecMAT%Iarray(98) = 1   !Assmebly complete

!C--------------------  START OF STEADY STATE  ----------------------------
      do
!C--------------------
        iterALL= iterALL + 1
        hecMAT%X = 0.0d0
!C
!C-- MATRIX ASSEMBLING

        call hecmw_barrier(hecMESH)
        call heat_mat_ass_conductivity( hecMESH, hecMAT, fstrHEAT, BETA )  
             write(IDBG,*) 'mat_ass_conductivity: OK'
             ! do i = 1, hecMESH%nn_internal
             !   write(IDBG,*) i, hecMAT%D(i)
             ! enddo
        call flush(IDBG)

        call heat_mat_ass_boundary( hecMESH, hecMAT, fstrHEAT, STIME, ETIME, 0.d0 )
             write(IDBG,*) 'mat_ass_boundary: OK'
             ! do i = 1, hecMESH%nn_internal
             !   write(IDBG,*) i, hecMAT%D(i)
             ! enddo
             call flush(IDBG)

        call hecmw_barrier(hecMESH)

!C
!C-- SOLVER
        hecMAT%Iarray(97) = 1   !Need numerical factorization
        bup_n_dof = hecMESH%n_dof
        hecMESH%n_dof = 1
        CALL solve_LINEQ(hecMESH,hecMAT,imsg)
        hecMESH%n_dof=bup_n_dof
        write(IDBG,*) 'solve_LINEQ: OK'
        call flush(IDBG)
!C
!C-- UPDATE

        do i= 1, hecMESH%n_node
          fstrHEAT%TEMPC(i)= fstrHEAT%TEMP(i)
        enddo

        do i= 1, hecMESH%n_node
          fstrHEAT%TEMP (i)= hecMAT%X(i)
        enddo

        VAL= 0.d0 
        do i= 1, hecMESH%nn_internal
          VAL= VAL + (fstrHEAT%TEMP(i) - fstrHEAT%TEMPC(i))**2
        enddo

        call hecmw_allREDUCE_R1 ( hecMESH, VAL, hecmw_sum )

        CHK = dsqrt(VAL)
        if( hecMESH%my_rank.eq.0 ) then 
          !write(*,'(i8,1p2e16.6,i10)')    iterALL,CHK,hecMAT%RESIDactual,hecMAT%ITERactual 
          !write(IMSG,'(i8,1p2e16.6,i10)') iterALL,CHK,hecMAT%RESIDactual,hecMAT%ITERactual 
          write(*,'(i8,1p1e16.6)')    iterALL,CHK
          write(IMSG,'(i8,1p1e16.6)') iterALL,CHK
          call flush(IMSG)
        endif

        if( CHK.lt.EPS ) then
          if( hecMESH%my_rank.eq.0 ) then
            write(*,*) 
            write(*,*) ' !!! CONVERGENCE ACHIEVED '
            write(IMSG,*) 
            write(IMSG,*) ' !!! CONVERGENCE ACHIEVED '
            INCR = 0
            !write(ISTA,'(3i8,1pE15.7,i8)') ISTEP,INCR,iterALL-1,CHK,hecMAT%ITERactual
            write(ISTA,'(3i8,1pE15.7)') ISTEP,INCR,iterALL-1,CHK
          endif
          exit
        endif

        if ( iterALL.ge.ITM ) then
          if( hecMESH%my_rank.eq.0 ) then
            write(*,*) 
            write(*,*) ' !!! ITERATION COUNT OVER : MAX = ', ITM
            write(IMSG,*) 
            write(IMSG,*) ' !!! ITERATION COUNT OVER : MAX = ', ITM
          endif
          call hecmw_abort( hecmw_comm_get_comm() )
        endif

!C--------------------
      enddo
!C--------------------  END OF STEADY STATE  ------------------

      TMAX = -1.0d10
      TMIN =  1.0d10
      LMAX = -1
      LMIN = -1

      write(ILOG,*)
      write(ILOG,'(a,i6)')    ' ISTEP =',ISTEP
      write(ILOG,'(a,f10.3)') ' Time  =',CTIME

      do i = 1, hecMESH%nn_internal
        inod = fstrPARAM%global_local_id(1,i)
        ii = fstrPARAM%global_local_id(2,i)
        temp = fstrHEAT%TEMP(ii)
        if( temp .gt. TMAX ) THEN
          TMAX = temp
          LMAX = inod
        endif
        if( temp .lt. TMIN ) THEN
          TMIN = temp
          LMIN = inod
        endif
      enddo

      write(ILOG,'(a,f10.3,i10)') ' Maximum Temperature :',TMAX
      write(ILOG,'(a,i10)')       ' Maximum Node No.    :',LMAX
      write(ILOG,'(a,f10.3,i10)') ' Minimum Temperature :',TMIN
      write(ILOG,'(a,i10)')       ' Minimum Node No.    :',LMIN

      if( IRESULT.eq.1 ) then
        header = '*fstrresult'
        call hecmw_result_init(hecMESH,1,1,header)
        label = 'TEMPERATURE'
        call hecmw_result_add(1,1,label,fstrHEAT%TEMP)
        nameID = 'fstrRES'
        call hecmw_result_write_by_name(nameID)
        call hecmw_result_finalize
      endif

      if( IVISUAL.eq.1 ) then
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
        call hecmw_visualize (hecMESH,fstrRESULT,1,1,1)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif

   end subroutine heat_solve_SS
end module m_heat_solve_SS
