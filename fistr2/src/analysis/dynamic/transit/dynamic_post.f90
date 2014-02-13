!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                       Tomotaka Ogasawara (Univ. of Tokyo)            !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!C================================================================C
!> This module contains functions to ouput max and min values of calculated result
!C================================================================C

module m_dynamic_post

contains



!C***
!C***  FSTR POST
!C***
      subroutine DYNAMIC_POST(hecMESH,hecMAT,fstrSOLID,istep,fstrDYNAMIC)

      use m_fstr
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix     ) :: hecMAT
      type (fstr_solid)         :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC

!C file name
      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN) :: label
      character(len=HECMW_NAME_LEN) :: nameID 
!C Max Min
      real(kind=kreal)    Umax(6), Umin(6), Emax(10), Emin(10),Smax(12), Smin(12)
      integer(kind=kint) IUmax(6),IUmin(6),IEmax(10),IEmin(10),ISmax(12),ISmin(12)

      real(kind=kreal)    EEmax(10), EEmin(10) ,ESmax(12), ESmin(12)
      integer(kind=kint) IEEmax(10),IEEmin(10),IESmax(12),IESmin(12)
!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if
!C Clear
      Umax=0.0d0
      Umin=0.0d0
      Emax=0.0d0
      Emin=0.0d0
      Smax=0.0d0
      Smin=0.0d0
      EEmax=0.0d0
      EEmin=0.0d0
      ESmax=0.0d0
      ESmin=0.0d0
      IUmax=0
      IUmin=0
      IEmax=0
      IEmin=0
      ISmax=0
      ISmin=0
      IEEmax=0
      IEEmin=0
      IESmax=0
      IESmin=0
!C*** Update Strain and Stress
      NDOF = hecMESH%n_dof
!      if( NDOF .eq. 2 ) NDOF = 3
!!      if( NDOF .eq. 6 ) then
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%STRAIN,hecMESH%n_node,10)
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%STRESS,hecMESH%n_node,12)
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%ESTRAIN,hecMESH%n_elem,10)
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%ESTRESS,hecMESH%n_elem,12)
!!      else 
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%STRAIN,hecMESH%n_node,6)
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%STRESS,hecMESH%n_node,7)
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%ESTRAIN,hecMESH%n_elem,6)
!!        call hecmw_update_m_R(hecMESH,fstrSOLID%ESTRESS,hecMESH%n_elem,7)
!!      endif
!C*** Show Displacement 
      write(IDBG,*) '#### DISPLACEMENT@POST'
      do i= 1, hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        if( i.eq.1) then
          do k=1,NDOF 
            Umax(k)=fstrSOLID%unode(NDOF*(i-1)+k)
            Umin(k)=fstrSOLID%unode(NDOF*(i-1)+k)
            IUmax(k)=j
            IUmin(k)=j
          enddo
        else
          do k=1,NDOF 
            if( fstrSOLID%unode(NDOF*(i-1)+k) .gt. Umax(k) ) then
              Umax(k)=fstrSOLID%unode(NDOF*(i-1)+k)
              IUmax(k)=j
            endif
            if( fstrSOLID%unode(NDOF*(i-1)+k) .lt. Umin(k) ) then
              Umin(k)=fstrSOLID%unode(NDOF*(i-1)+k)
              IUmin(k)=j
            endif
          enddo
        endif
      enddo
!C*** Show Strain
      if( NDOF .eq. 3 ) MDOF = 6
      if( NDOF .eq. 6 ) MDOF = 10
      write(IDBG,*) '#### STRAIN@POST'
!C @node
      do i= 1, hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        if( i.eq.1) then
          do k=1,MDOF 
            Emax(k)=fstrSOLID%STRAIN(MDOF*(i-1)+k)
            Emin(k)=fstrSOLID%STRAIN(MDOF*(i-1)+k)
            IEmax(k)=j
            IEmin(k)=j
          enddo
        else
          do k=1,MDOF
            if( fstrSOLID%STRAIN(MDOF*(i-1)+k) .gt. Emax(k) ) then
              Emax(k)=fstrSOLID%STRAIN(MDOF*(i-1)+k)
              IEmax(k)=j
            endif
            if( fstrSOLID%STRAIN(MDOF*(i-1)+k) .lt. Emin(k) ) then
              Emin(k)=fstrSOLID%STRAIN(MDOF*(i-1)+k)
              IEmin(k)=j
            endif
          enddo
        endif
      enddo
!C @element
      do i= 1, hecMESH%n_elem
!!      ii = hecMESH%elem_ID(i*2-1)
        ii = i
        ID_area = hecMESH%elem_ID(i*2)
       if( ID_area.eq.hecMESH%my_rank ) then
          j=hecMESH%global_elem_ID(ii)
        if( i.eq.1) then
          do k=1,MDOF 
            EEmax(k)=fstrSOLID%ESTRAIN(MDOF*(i-1)+k)
            EEmin(k)=fstrSOLID%ESTRAIN(MDOF*(i-1)+k)
            IEEmax(k)=j
            IEEmin(k)=j
          enddo
        else
          do k=1,MDOF
            if( fstrSOLID%ESTRAIN(MDOF*(i-1)+k) .gt. EEmax(k) ) then
              EEmax(k)=fstrSOLID%ESTRAIN(MDOF*(i-1)+k)
              IEEmax(k)=j
            endif
            if( fstrSOLID%ESTRAIN(MDOF*(i-1)+k) .lt. EEmin(k) ) then
              EEmin(k)=fstrSOLID%ESTRAIN(MDOF*(i-1)+k)
              IEEmin(k)=j
            endif
          enddo
        endif
       endif
      enddo
!C*** Show Stress
      if( NDOF .eq. 3 ) MDOF = 7
      if( NDOF .eq. 6 ) MDOF = 12
      write(IDBG,*) '#### STRESS@POST'
!C @node
      do i= 1, hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        if( i.eq.1) then
          do k=1,MDOF
            Smax(k)=fstrSOLID%STRESS(MDOF*(i-1)+k)
            Smin(k)=fstrSOLID%STRESS(MDOF*(i-1)+k)
            ISmax(k)=j
            ISmin(k)=j
          enddo
        else
          do k=1,MDOF
            if( fstrSOLID%STRESS(MDOF*(i-1)+k) .gt. Smax(k) ) then
              Smax(k)=fstrSOLID%STRESS(MDOF*(i-1)+k)
              ISmax(k)=j
            endif
            if( fstrSOLID%STRESS(MDOF*(i-1)+k) .lt. Smin(k) ) then
              Smin(k)=fstrSOLID%STRESS(MDOF*(i-1)+k)
              ISmin(k)=j
            endif
          enddo
        endif
      enddo
!C @element
      do i= 1, hecMESH%n_elem
!!      ii = hecMESH%elem_ID(i*2-1)
        ii = i
        ID_area = hecMESH%elem_ID(i*2)
       if( ID_area.eq.hecMESH%my_rank ) then
          j=hecMESH%global_elem_ID(ii)
        if( i.eq.1) then
          do k=1,MDOF
            ESmax(k)=fstrSOLID%ESTRESS(MDOF*(i-1)+k)
            ESmin(k)=fstrSOLID%ESTRESS(MDOF*(i-1)+k)
            IESmax(k)=j
            IESmin(k)=j
          enddo
        else
          do k=1,MDOF
            if( fstrSOLID%ESTRESS(MDOF*(i-1)+k) .gt. ESmax(k) ) then
              ESmax(k)=fstrSOLID%ESTRESS(MDOF*(i-1)+k)
              IESmax(k)=j
            endif
            if( fstrSOLID%ESTRESS(MDOF*(i-1)+k) .lt. ESmin(k) ) then
              ESmin(k)=fstrSOLID%ESTRESS(MDOF*(i-1)+k)
              IESmin(k)=j
            endif
          enddo
        endif
       endif
      enddo
!C
      if( NDOF .eq. 3 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//U1 ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009) '//U2 ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009) '//U3 ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,1009) '//E11',Emax(1),IEmax(1),Emin(1),IEmin(1)
        write(ILOG,1009) '//E22',Emax(2),IEmax(2),Emin(2),IEmin(2)
        write(ILOG,1009) '//E33',Emax(3),IEmax(3),Emin(3),IEmin(3)
        write(ILOG,1009) '//E12',Emax(4),IEmax(4),Emin(4),IEmin(4)
        write(ILOG,1009) '//E23',Emax(5),IEmax(5),Emin(5),IEmin(5)
        write(ILOG,1009) '//E13',Emax(6),IEmax(6),Emin(6),IEmin(6)
        write(ILOG,1009) '//S11',Smax(1),ISmax(1),Smin(1),ISmin(1)
        write(ILOG,1009) '//S22',Smax(2),ISmax(2),Smin(2),ISmin(2)
        write(ILOG,1009) '//S33',Smax(3),ISmax(3),Smin(3),ISmin(3)
        write(ILOG,1009) '//S12',Smax(4),ISmax(4),Smin(4),ISmin(4)
        write(ILOG,1009) '//S23',Smax(5),ISmax(5),Smin(5),ISmin(5)
        write(ILOG,1009) '//S13',Smax(6),ISmax(6),Smin(6),ISmin(6)
        write(ILOG,1009) '//SMS',Smax(7),ISmax(7),Smin(7),ISmin(7)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//E11',EEmax(1),IEEmax(1),EEmin(1),IEEmin(1)
        write(ILOG,1009) '//E22',EEmax(2),IEEmax(2),EEmin(2),IEEmin(2)
        write(ILOG,1009) '//E33',EEmax(3),IEEmax(3),EEmin(3),IEEmin(3)
        write(ILOG,1009) '//E12',EEmax(4),IEEmax(4),EEmin(4),IEEmin(4)
        write(ILOG,1009) '//E23',EEmax(5),IEEmax(5),EEmin(5),IEEmin(5)
        write(ILOG,1009) '//E13',EEmax(6),IEEmax(6),EEmin(6),IEEmin(6)
        write(ILOG,1009) '//S11',ESmax(1),IESmax(1),ESmin(1),IESmin(1)
        write(ILOG,1009) '//S22',ESmax(2),IESmax(2),ESmin(2),IESmin(2)
        write(ILOG,1009) '//S33',ESmax(3),IESmax(3),ESmin(3),IESmin(3)
        write(ILOG,1009) '//S12',ESmax(4),IESmax(4),ESmin(4),IESmin(4)
        write(ILOG,1009) '//S23',ESmax(5),IESmax(5),ESmin(5),IESmin(5)
        write(ILOG,1009) '//S13',ESmax(6),IESmax(6),ESmin(6),IESmin(6)
        write(ILOG,1009) '//SMS',ESmax(7),IESmax(7),ESmin(7),IESmin(7)
!C*** Show Summary
        if (.not. ds) then !In case of Direct Solver prevent MPI
      !    call hecmw_allREDUCE_R(hecMESH,Umax,3,hecmw_max)
     !     call hecmw_allREDUCE_R(hecMESH,Umin,3,hecmw_min)
     !     call hecmw_allREDUCE_R(hecMESH,Emax,6,hecmw_max)
     !     call hecmw_allREDUCE_R(hecMESH,Emin,6,hecmw_min)
      !    call hecmw_allREDUCE_R(hecMESH,Smax,7,hecmw_max)
      !    call hecmw_allREDUCE_R(hecMESH,Smin,7,hecmw_min)
      !    call hecmw_allREDUCE_R(hecMESH,EEmax,6,hecmw_max)
      !    call hecmw_allREDUCE_R(hecMESH,EEmin,6,hecmw_min)
       !   call hecmw_allREDUCE_R(hecMESH,ESmax,7,hecmw_max)
       !   call hecmw_allREDUCE_R(hecMESH,ESmin,7,hecmw_min)
        end if
        if( hecMESH%my_rank .eq. 0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1 ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2 ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3 ',Umax(3),Umin(3)
          write(ILOG,1019) '//E11',Emax(1),Emin(1)
          write(ILOG,1019) '//E22',Emax(2),Emin(2)
          write(ILOG,1019) '//E33',Emax(3),Emin(3)
          write(ILOG,1019) '//E12',Emax(4),Emin(4)
          write(ILOG,1019) '//E23',Emax(5),Emin(5)
          write(ILOG,1019) '//E13',Emax(6),Emin(6)
          write(ILOG,1019) '//S11',Smax(1),Smin(1)
          write(ILOG,1019) '//S22',Smax(2),Smin(2)
          write(ILOG,1019) '//S33',Smax(3),Smin(3)
          write(ILOG,1019) '//S12',Smax(4),Smin(4)
          write(ILOG,1019) '//S23',Smax(5),Smin(5)
          write(ILOG,1019) '//S13',Smax(6),Smin(6)
          write(ILOG,1019) '//SMS',Smax(7),Smin(7)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11',EEmax(1),EEmin(1)
          write(ILOG,1019) '//E22',EEmax(2),EEmin(2)
          write(ILOG,1019) '//E33',EEmax(3),EEmin(3)
          write(ILOG,1019) '//E12',EEmax(4),EEmin(4)
          write(ILOG,1019) '//E23',EEmax(5),EEmin(5)
          write(ILOG,1019) '//E13',EEmax(6),EEmin(6)
          write(ILOG,1019) '//S11',ESmax(1),ESmin(1)
          write(ILOG,1019) '//S22',ESmax(2),ESmin(2)
          write(ILOG,1019) '//S33',ESmax(3),ESmin(3)
          write(ILOG,1019) '//S12',ESmax(4),ESmin(4)
          write(ILOG,1019) '//S23',ESmax(5),ESmin(5)
          write(ILOG,1019) '//S13',ESmax(6),ESmin(6)
          write(ILOG,1019) '//SMS',ESmax(7),ESmin(7)
        endif
        if( IRESULT.eq.0 ) return
!C*** Write Result File
!C*** INITIALIZE
        header='*fstrresult'
        nd = istep
     !   call hecmw_result_init(hecMESH,nd,header)
!C*** ADD 
!C-- implicit dynamic analysis
        if(fstrDYNAMIC%idx_eqa .eq. 1) then
!C*** DISPLACEMENT
          if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
            id = 1
            ndof=3
            label='DISPLACEMENT'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,2))
            end if
          end if
!C*** VELOCITY
          if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
            id = 1
            ndof=3
            label='VELOCITY'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,2))
            end if
          end if
!C*** ACCELERATION
          if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
            id = 1
            ndof=3
            label='ACCELERATION'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,2))
            end if
          end if
!C-- explicit dynamic analysis
        else if(fstrDYNAMIC%idx_eqa .eq. 11) then
!C*** DISPLACEMENT
          if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
            id = 1
            ndof=3
            label='DISPLACEMENT'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,1))
            end if
          end if
!C*** VELOCITY
          if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
            id = 1
            ndof=3
            label='VELOCITY'
            if(fstrDYNAMIC%i_step .eq. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,1))
            end if
          end if
!C*** ACCELERATION
          if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
            id = 1
            ndof=3
            label='ACCELERATION'
            if(fstrDYNAMIC%i_step .eq. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,1))
            end if
          end if
        end if
!
        if(fstrDYNAMIC%i_step .gt. 0) then
          if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
!C*** STRAIN @node
            id = 1
            ndof=6
            label='STRAIN'
         !   call hecmw_result_add(id,ndof,label,fstrSOLID%STRAIN)
!C*** STRAIN @element 
            id = 2
            ndof=6
            label='ESTRAIN'
        !    call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRAIN)
          end if
!
          if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
!C*** STRESS @node
            id = 1
            ndof=7
            label='STRESS'
         !   call hecmw_result_add(id,ndof,label,fstrSOLID%STRESS)
!C*** STRESS @element
            id = 2
            ndof=7
            label='ESTRESS'
         !   call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRESS)
          end if
        end if
!C*** WRITE NOW
        nameID='fstrRES'
      !  call hecmw_result_write_by_name(nameID)
!C*** FINALIZE 
     !   call hecmw_result_finalize
!C
!C*for 6D
      elseif( NDOF .eq. 6 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//U1    ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009)'//U2    ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009)'//U3    ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,1009)'//R1    ',Umax(4),IUmax(4),Umin(4),IUmin(4)
        write(ILOG,1009)'//R2    ',Umax(5),IUmax(5),Umin(5),IUmin(5)
        write(ILOG,1009)'//R3    ',Umax(6),IUmax(6),Umin(6),IUmin(6)
        write(ILOG,1009)'//E11(+)',Emax( 1),IEmax( 1),Emin( 1),IEmin( 1)
        write(ILOG,1009)'//E22(+)',Emax( 2),IEmax( 2),Emin( 2),IEmin( 2)
        write(ILOG,1009)'//E12(+)',Emax( 3),IEmax( 3),Emin( 3),IEmin( 3)
        write(ILOG,1009)'//E23(+)',Emax( 4),IEmax( 4),Emin( 4),IEmin( 4)
        write(ILOG,1009)'//E13(+)',Emax( 5),IEmax( 5),Emin( 5),IEmin( 5)
        write(ILOG,1009)'//E11(-)',Emax( 6),IEmax( 6),Emin( 6),IEmin( 6)
        write(ILOG,1009)'//E22(-)',Emax( 7),IEmax( 7),Emin( 7),IEmin( 7)
        write(ILOG,1009)'//E12(-)',Emax( 8),IEmax( 8),Emin( 8),IEmin( 8)
        write(ILOG,1009)'//E23(-)',Emax( 9),IEmax( 9),Emin( 9),IEmin( 9)
        write(ILOG,1009)'//E13(-)',Emax(10),IEmax(10),Emin(10),IEmin(10)
        write(ILOG,1009)'//S11(+)',Smax( 1),ISmax( 1),Smin( 1),ISmin( 1)
        write(ILOG,1009)'//S22(+)',Smax( 2),ISmax( 2),Smin( 2),ISmin( 2)
        write(ILOG,1009)'//S12(+)',Smax( 3),ISmax( 3),Smin( 3),ISmin( 3)
        write(ILOG,1009)'//S23(+)',Smax( 4),ISmax( 4),Smin( 4),ISmin( 4)
        write(ILOG,1009)'//S13(+)',Smax( 5),ISmax( 5),Smin( 5),ISmin( 5)
        write(ILOG,1009)'//SMS(+)',Smax( 6),ISmax( 6),Smin( 6),ISmin( 6)
        write(ILOG,1009)'//S11(-)',Smax( 7),ISmax( 7),Smin( 7),ISmin( 7)
        write(ILOG,1009)'//S22(-)',Smax( 8),ISmax( 8),Smin( 8),ISmin( 8)
        write(ILOG,1009)'//S12(-)',Smax( 9),ISmax( 9),Smin( 9),ISmin( 9)
        write(ILOG,1009)'//S23(-)',Smax(10),ISmax(10),Smin(10),ISmin(10)
        write(ILOG,1009)'//S13(-)',Smax(11),ISmax(11),Smin(11),ISmin(11)
        write(ILOG,1009)'//SMS(-)',Smax(12),ISmax(12),Smin(12),ISmin(12)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//E11(+)',EEmax( 1),IEEmax( 1),EEmin( 1),IEEmin( 1)
        write(ILOG,1009)'//E22(+)',EEmax( 2),IEEmax( 2),EEmin( 2),IEEmin( 2)
        write(ILOG,1009)'//E12(+)',EEmax( 3),IEEmax( 3),EEmin( 3),IEEmin( 3)
        write(ILOG,1009)'//E23(+)',EEmax( 4),IEEmax( 4),EEmin( 4),IEEmin( 4)
        write(ILOG,1009)'//E13(+)',EEmax( 5),IEEmax( 5),EEmin( 5),IEEmin( 5)
        write(ILOG,1009)'//E11(-)',EEmax( 6),IEEmax( 6),EEmin( 6),IEEmin( 6)
        write(ILOG,1009)'//E22(-)',EEmax( 7),IEEmax( 7),EEmin( 7),IEEmin( 7)
        write(ILOG,1009)'//E12(-)',EEmax( 8),IEEmax( 8),EEmin( 8),IEEmin( 8)
        write(ILOG,1009)'//E23(-)',EEmax( 9),IEEmax( 9),EEmin( 9),IEEmin( 9)
        write(ILOG,1009)'//E13(-)',EEmax(10),IEEmax(10),EEmin(10),IEEmin(10)
        write(ILOG,1009)'//S11(+)',ESmax( 1),IESmax( 1),ESmin( 1),IESmin( 1)
        write(ILOG,1009)'//S22(+)',ESmax( 2),IESmax( 2),ESmin( 2),IESmin( 2)
        write(ILOG,1009)'//S12(+)',ESmax( 3),IESmax( 3),ESmin( 3),IESmin( 3)
        write(ILOG,1009)'//S23(+)',ESmax( 4),IESmax( 4),ESmin( 4),IESmin( 4)
        write(ILOG,1009)'//S13(+)',ESmax( 5),IESmax( 5),ESmin( 5),IESmin( 5)
        write(ILOG,1009)'//SMS(+)',ESmax( 6),IESmax( 6),ESmin( 6),IESmin( 6)
        write(ILOG,1009)'//S11(-)',ESmax( 7),IESmax( 7),ESmin( 7),IESmin( 7)
        write(ILOG,1009)'//S22(-)',ESmax( 8),IESmax( 8),ESmin( 8),IESmin( 8)
        write(ILOG,1009)'//S12(-)',ESmax( 9),IESmax( 9),ESmin( 9),IESmin( 9)
        write(ILOG,1009)'//S23(-)',ESmax(10),IESmax(10),ESmin(10),IESmin(10)
        write(ILOG,1009)'//S13(-)',ESmax(11),IESmax(11),ESmin(11),IESmin(11)
        write(ILOG,1009)'//SMS(-)',ESmax(12),IESmax(12),ESmin(12),IESmin(12)
!C*** Show Summary
        if (.not. ds) then !In case of Direct Solver prevent MPI
      !    call hecmw_allREDUCE_R(hecMESH,Umax, 6,hecmw_max)
      !    call hecmw_allREDUCE_R(hecMESH,Umin, 6,hecmw_min)
      !    call hecmw_allREDUCE_R(hecMESH,Emax,10,hecmw_max)
      !    call hecmw_allREDUCE_R(hecMESH,Emin,10,hecmw_min)
      !    call hecmw_allREDUCE_R(hecMESH,Smax,12,hecmw_max)
       !   call hecmw_allREDUCE_R(hecMESH,Smin,12,hecmw_min)
      !    call hecmw_allREDUCE_R(hecMESH,EEmax,10,hecmw_max)
      !    call hecmw_allREDUCE_R(hecMESH,EEmin,10,hecmw_min)
     !     call hecmw_allREDUCE_R(hecMESH,ESmax,12,hecmw_max)
      !    call hecmw_allREDUCE_R(hecMESH,ESmin,12,hecmw_min)
        end if
        if( hecMESH%my_rank .eq. 0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1    ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2    ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3    ',Umax(3),Umin(3)
          write(ILOG,1019) '//R1    ',Umax(4),Umin(4)
          write(ILOG,1019) '//R2    ',Umax(5),Umin(5)
          write(ILOG,1019) '//R3    ',Umax(6),Umin(6)
          write(ILOG,1019) '//E11(+)',Emax( 1),Emin( 1)
          write(ILOG,1019) '//E22(+)',Emax( 2),Emin( 2)
          write(ILOG,1019) '//E12(+)',Emax( 3),Emin( 3)
          write(ILOG,1019) '//E23(+)',Emax( 4),Emin( 4)
          write(ILOG,1019) '//E13(+)',Emax( 5),Emin( 5)
          write(ILOG,1019) '//E11(-)',Emax( 6),Emin( 6)
          write(ILOG,1019) '//E22(-)',Emax( 7),Emin( 7)
          write(ILOG,1019) '//E12(-)',Emax( 8),Emin( 8)
          write(ILOG,1019) '//E23(-)',Emax( 9),Emin( 9)
          write(ILOG,1019) '//E13(-)',Emax(10),Emin(10)
          write(ILOG,1019) '//S11(+)',Smax( 1),Smin( 1)
          write(ILOG,1019) '//S22(+)',Smax( 2),Smin( 2)
          write(ILOG,1019) '//S12(+)',Smax( 3),Smin( 3)
          write(ILOG,1019) '//S23(+)',Smax( 4),Smin( 4)
          write(ILOG,1019) '//S13(+)',Smax( 5),Smin( 5)
          write(ILOG,1019) '//SMS(+)',Smax( 6),Smin( 6)
          write(ILOG,1019) '//S11(-)',Smax( 7),Smin( 7)
          write(ILOG,1019) '//S22(-)',Smax( 8),Smin( 8)
          write(ILOG,1019) '//S12(-)',Smax( 9),Smin( 9)
          write(ILOG,1019) '//S23(-)',Smax(10),Smin(10)
          write(ILOG,1019) '//S13(-)',Smax(11),Smin(11)
          write(ILOG,1019) '//SMS(-)',Smax(12),Smin(12)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11(+)',EEmax( 1),EEmin( 1)
          write(ILOG,1019) '//E22(+)',EEmax( 2),EEmin( 2)
          write(ILOG,1019) '//E12(+)',EEmax( 3),EEmin( 3)
          write(ILOG,1019) '//E23(+)',EEmax( 4),EEmin( 4)
          write(ILOG,1019) '//E13(+)',EEmax( 5),EEmin( 5)
          write(ILOG,1019) '//E11(-)',EEmax( 6),EEmin( 6)
          write(ILOG,1019) '//E22(-)',EEmax( 7),EEmin( 7)
          write(ILOG,1019) '//E12(-)',EEmax( 8),EEmin( 8)
          write(ILOG,1019) '//E23(-)',EEmax( 9),EEmin( 9)
          write(ILOG,1019) '//E13(-)',EEmax(10),EEmin(10)
          write(ILOG,1019) '//S11(+)',ESmax( 1),ESmin( 1)
          write(ILOG,1019) '//S22(+)',ESmax( 2),ESmin( 2)
          write(ILOG,1019) '//S12(+)',ESmax( 3),ESmin( 3)
          write(ILOG,1019) '//S23(+)',ESmax( 4),ESmin( 4)
          write(ILOG,1019) '//S13(+)',ESmax( 5),ESmin( 5)
          write(ILOG,1019) '//SMS(+)',ESmax( 6),ESmin( 6)
          write(ILOG,1019) '//S11(-)',ESmax( 7),ESmin( 7)
          write(ILOG,1019) '//S22(-)',ESmax( 8),ESmin( 8)
          write(ILOG,1019) '//S12(-)',ESmax( 9),ESmin( 9)
          write(ILOG,1019) '//S23(-)',ESmax(10),ESmin(10)
          write(ILOG,1019) '//S13(-)',ESmax(11),ESmin(11)
          write(ILOG,1019) '//SMS(-)',ESmax(12),ESmin(12)
        endif
        if( IRESULT.eq.0 ) return
!C*** Write Result File
!C*** INITIALIZE
        header='*fstrresult'
        nd = istep
     !   call hecmw_result_init(hecMESH,nd,header)
!C*** ADD 

!C--- shell DISPLACEMENT
!        id = 1
!        ndof=6
!        label='DISPLACEMENT'
!        call hecmw_result_add(id,ndof,label,fstrSOLID%unode)
!C--- shell DISPLACEMENT

!C-- implicit dynamic analysis
        if(fstrDYNAMIC%idx_eqa .eq. 1) then
!C*** DISPLACEMENT
          if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
            id = 1
            ndof=6
            label='DISPLACEMENT'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,2))
            end if
          end if
!C*** VELOCITY
          if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
            id = 1
            ndof=6
            label='VELOCITY'
            if(fstrDYNAMIC%i_step .eq. 0) then
        !       call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
        !       call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,2))
            end if
          end if
!C*** ACCELERATION
          if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
            id = 1
            ndof=6
            label='ACCELERATION'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,2))
            end if
          end if
!C-- explicit dynamic analysis
        else if(fstrDYNAMIC%idx_eqa .eq. 11) then
!C*** DISPLACEMENT
          if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
            id = 1
            ndof=6
            label='DISPLACEMENT'
            if(fstrDYNAMIC%i_step .eq. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
         !      call hecmw_result_add(id,ndof,label,fstrDYNAMIC%DISP(:,1))
            end if
          end if
!C*** VELOCITY
          if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
            id = 1
            ndof=6
            label='VELOCITY'
            if(fstrDYNAMIC%i_step .eq. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%VEL(:,1))
            end if
          end if
!C*** ACCELERATION
          if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
            id = 1
            ndof=6
            label='ACCELERATION'
            if(fstrDYNAMIC%i_step .eq. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,1))
            else if(fstrDYNAMIC%i_step .gt. 0) then
          !     call hecmw_result_add(id,ndof,label,fstrDYNAMIC%ACC(:,1))
            end if
          end if
        end if
!
        if(fstrDYNAMIC%i_step .gt. 0) then
          if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
!C*** STRAIN @node
            id = 1
            ndof=10
            label='STRAIN'
        !    call hecmw_result_add(id,ndof,label,fstrSOLID%STRAIN)
!C*** STRAIN @element 
            id = 2
            ndof=10
            label='ESTRAIN'
        !    call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRAIN)
          end if
!
          if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
!C*** STRESS @node
            id = 1
            ndof=12
            label='STRESS'
         !   call hecmw_result_add(id,ndof,label,fstrSOLID%STRESS)
!C*** STRESS @element
            id = 2
            ndof=12
            label='ESTRESS'
        !    call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRESS)
          end if
        end if
!C*** WRITE NOW
        nameID='fstrRES'
      !  call hecmw_result_write_by_name(nameID)
!C*** FINALIZE 
      !  call hecmw_result_finalize
!C
      endif
 1009 format(a8,1pe12.4,i10,1pe12.4,i10)
 1019 format(a8,1pe12.4,1pe12.4)
      end subroutine DYNAMIC_POST

end module m_dynamic_post
