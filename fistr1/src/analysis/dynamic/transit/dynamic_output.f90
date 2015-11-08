!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.6                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Suemits(Advancesoft)                        !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions to output result.
!!
!>  \author     K. Suemitsu(AdavanceSoft)
!>  \date       2012/01/16
!>  \version    0.00
!!
!======================================================================!
module m_dynamic_output
  implicit none

  contains

!> Output result
!----------------------------------------------------------------------*
  subroutine fstr_dynamic_Output( hecMESH, fstrSOLID, fstrDYNAMIC )
!----------------------------------------------------------------------*
    use m_fstr
    use m_fstr_Update
    use m_fstr_NodalStress
    use m_dynamic_make_result
    use m_hecmw2fstr_mesh_conv
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(inout)      :: fstrSOLID
    type (fstr_dynamic), intent(in)       :: fstrDYNAMIC

    type ( hecmwST_result_data ) :: fstrRESULT
    integer(kind=kint) :: i, j, ndof, maxstep, interval, fnum, iS, iE, gid, istep, idx

    ndof = hecMESH%n_dof

!C-- SET DISPLACEMENT etc.
    istep = fstrDYNAMIC%i_step
    if( fstrDYNAMIC%idx_eqa==1 .and. istep>0 ) then
      idx = 2
    else
      idx = 1
    endif
    if( fstrDYNAMIC%nlflag==0 ) then
      do i = 1, hecMESH%n_node
        do j = 1, ndof
          fstrSOLID%unode(ndof*(i-1)+j) = fstrDYNAMIC%DISP(ndof*(i-1)+j,idx)
        enddo
      enddo
      if( ndof==3 ) then
        call fstr_Update3D( hecMESH, fstrSOLID )
      else if( ndof==2 ) then
        call fstr_Update2D( hecMESH, fstrSOLID )
      else if( ndof==6) then
        call fstr_Update6D( hecMESH, fstrSOLID )
      endif
    endif

    if( fstrSOLID%TEMP_ngrp_tot>0 .or. fstrSOLID%TEMP_irres>0 ) then
       if( ndof==3 ) then
          allocate( fstrSOLID%tnstrain(hecMESH%n_node*6) )
          allocate( fstrSOLID%testrain(hecMESH%n_elem*6) )
       else if( ndof==2 ) then
          allocate( fstrSOLID%tnstrain(hecMESH%n_node*3) )
          allocate( fstrSOLID%testrain(hecMESH%n_elem*3) )
       endif
    endif

    if( ndof==3 ) then
          call fstr_NodalStress3D( hecMESH, fstrSOLID )
    else if( ndof==2 ) then
          call fstr_NodalStress2D( hecMESH, fstrSOLID )
    else if( ndof==6 ) then
          call fstr_NodalStress6D( hecMESH, fstrSOLID )
    endif

    maxstep = fstrDYNAMIC%n_step

    if( (mod(istep,fstrSOLID%output_ctrl(1)%freqency)==0 .or. istep==maxstep) ) then
          fnum = fstrSOLID%output_ctrl(1)%filenum
          call fstr_dynamic_post( fnum, hecMESH, fstrSOLID, fstrDYNAMIC )
    endif

    if( fstrSOLID%output_ctrl(2)%outinfo%grp_id>0 .and. &
        (mod(istep,fstrSOLID%output_ctrl(2)%freqency)==0 .or. istep==maxstep) ) then
      iS = fstrSOLID%output_ctrl(2)%outinfo%grp_id
      fnum = fstrSOLID%output_ctrl(2)%filenum
      do i = hecMESH%node_group%grp_index(iS-1)+1, hecMESH%node_group%grp_index(iS)
        iE = hecMESH%node_group%grp_item(i)
        gid = hecMESH%global_node_ID(iE)
        write(fnum,'(2i10,1p6e15.7)') istep,gid,(fstrSOLID%unode(ndof*(iE-1)+j),j=1,ndof)
      enddo
    endif

    if( IRESULT==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(3)%freqency)==0 .or. istep==maxstep) ) then
          call fstr_write_dynamic_result( hecMESH, fstrSOLID, fstrDYNAMIC, maxstep, istep )
    endif

    if( IVISUAL==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(4)%freqency)==0 .or. istep==maxstep) ) then
          interval = fstrSOLID%output_ctrl(4)%freqency
          call fstr_make_dynamic_result( hecMESH, fstrSOLID, fstrDYNAMIC, fstrRESULT )
          call fstr2hecmw_mesh_conv( hecMESH )
          call hecmw_visualize_init
          call hecmw_visualize( hecMESH, fstrRESULT, istep, maxstep, interval )
          call hecmw_visualize_finalize
          call hecmw2fstr_mesh_conv( hecMESH )
          call hecmw_result_free( fstrRESULT )
    endif

  end subroutine fstr_dynamic_Output

!> Summarizer of output data which prints out max and min output values
!----------------------------------------------------------------------*
  subroutine fstr_dynamic_post( fnum, hecMESH, fstrSOLID, fstrDYNAMIC )
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(in)                   :: fnum
      type (hecmwST_local_mesh), intent(in) :: hecMESH
      type (fstr_solid), intent(in)         :: fstrSOLID
      type (fstr_dynamic), intent(in)       :: fstrDYNAMIC

      real(kind=kreal)   :: Umax(6), Umin(6), Emax(14), Emin(14), Smax(14), Smin(14)
      integer(kind=kint) :: IUmax(6), IUmin(6), IEmax(14), IEmin(14), ISmax(14), ISmin(14)
      real(kind=kreal)   :: EEmax(14), EEmin(14), ESmax(14), ESmin(14)
      integer(kind=kint) :: IEEmax(14), IEEmin(14), IESmax(14), IESmin(14)
      real(kind=kreal)   :: Vmax(6), Vmin(6), Amax(6), Amin(6)
      integer(kind=kint) :: IVmax(6), IVmin(6), IAmax(6), IAmin(6)
      real(kind=kreal)   :: Mmax(1), Mmin(1), EMmax(1), EMmin(1)
      integer(kind=kint) :: IMmax(1), IMmin(1), IEMmax(1), IEMmin(1)
      integer(kind=kint) :: i, j, k, ndof, mdof, ID_area, idx

      if( fstrDYNAMIC%i_step==0 ) return

      if( fstrDYNAMIC%idx_eqa==1 .and. fstrDYNAMIC%i_step>0 ) then
        idx = 2
      else
        idx = 1
      endif

      ndof = hecMESH%n_dof

      write( fnum, '(''#### Result step='',I6)') fstrDYNAMIC%i_step

!C*** Show Displacement
      do i = 1, hecMESH%nn_internal
        if(fstrSOLID%is_rot(i)==1)cycle
        j = hecMESH%global_node_ID(i)
        if( i==1 ) then
          do k = 1, ndof
            Umax(k) = fstrDYNAMIC%DISP(ndof*(i-1)+k,idx)
            Umin(k) = fstrDYNAMIC%DISP(ndof*(i-1)+k,idx)
            IUmax(k)= j
            IUmin(k)= j
            Vmax(k) = fstrDYNAMIC%VEL(ndof*(i-1)+k,idx)
            Vmin(k) = fstrDYNAMIC%VEL(ndof*(i-1)+k,idx)
            IVmax(k)= j
            IVmin(k)= j
            Amax(k) = fstrDYNAMIC%ACC(ndof*(i-1)+k,idx)
            Amin(k) = fstrDYNAMIC%ACC(ndof*(i-1)+k,idx)
            IAmax(k)= j
            IAmin(k)= j
          enddo
        else
          do k = 1, ndof
            if( fstrDYNAMIC%DISP(ndof*(i-1)+k,idx) > Umax(k) ) then
              Umax(k) = fstrDYNAMIC%DISP(ndof*(i-1)+k,idx)
              IUmax(k)= j
            endif
            if( fstrDYNAMIC%DISP(ndof*(i-1)+k,idx) < Umin(k) ) then
              Umin(k) = fstrDYNAMIC%DISP(ndof*(i-1)+k,idx)
              IUmin(k)= j
            endif
            if( fstrDYNAMIC%VEL(ndof*(i-1)+k,idx) > Vmax(k) ) then
              Vmax(k) = fstrDYNAMIC%VEL(ndof*(i-1)+k,idx)
              IVmax(k)= j
            endif
            if( fstrDYNAMIC%VEL(ndof*(i-1)+k,idx) < Vmin(k) ) then
              Vmin(k) = fstrDYNAMIC%VEL(ndof*(i-1)+k,idx)
              IVmin(k)= j
            endif
            if( fstrDYNAMIC%ACC(ndof*(i-1)+k,idx) > Amax(k) ) then
              Amax(k) = fstrDYNAMIC%ACC(ndof*(i-1)+k,idx)
              IAmax(k)= j
            endif
            if( fstrDYNAMIC%ACC(ndof*(i-1)+k,idx) < Amin(k) ) then
              Amin(k) = fstrDYNAMIC%ACC(ndof*(i-1)+k,idx)
              IAmin(k)= j
            endif
          enddo
        endif
      enddo
!C*** Show Strain
      if( ndof==2 ) mdof = 3
      if( ndof==3 ) mdof = 6
      if( ndof==6 ) mdof = 6
!C @node
      do i = 1, hecMESH%nn_internal
        if(fstrSOLID%is_rot(i)==1)cycle
        j = hecMESH%global_node_ID(i)
        if( i==1 ) then
          do k = 1, mdof
            Emax(k) = fstrSOLID%STRAIN(mdof*(i-1)+k)
            Emin(k) = fstrSOLID%STRAIN(mdof*(i-1)+k)
            IEmax(k)= j
            IEmin(k)= j
          enddo
        else
          do k = 1, mdof
            if( fstrSOLID%STRAIN(mdof*(i-1)+k) > Emax(k) ) then
              Emax(k) = fstrSOLID%STRAIN(mdof*(i-1)+k)
              IEmax(k)= j
            endif
            if( fstrSOLID%STRAIN(mdof*(i-1)+k) < Emin(k) ) then
              Emin(k) = fstrSOLID%STRAIN(mdof*(i-1)+k)
              IEmin(k)= j
            endif
          enddo
        endif
      enddo
!C @element
      do i = 1, hecMESH%n_elem
        ID_area = hecMESH%elem_ID(i*2)
        if( ID_area==hecMESH%my_rank ) then
          j = hecMESH%global_elem_ID(i)
          if( hecMESH%elem_ID(i*2-1)==1 ) then
            do k = 1, mdof
              EEmax(k) = fstrSOLID%ESTRAIN(mdof*(i-1)+k)
              EEmin(k) = fstrSOLID%ESTRAIN(mdof*(i-1)+k)
              IEEmax(k)= j
              IEEmin(k)= j
            enddo
          else
            do k = 1, mdof
              if( fstrSOLID%ESTRAIN(mdof*(i-1)+k) > EEmax(k) ) then
                EEmax(k) = fstrSOLID%ESTRAIN(mdof*(i-1)+k)
                IEEmax(k)= j
              endif
              if( fstrSOLID%ESTRAIN(mdof*(i-1)+k) < EEmin(k) ) then
                EEmin(k) = fstrSOLID%ESTRAIN(mdof*(i-1)+k)
                IEEmin(k)= j
              endif
            enddo
          endif
        endif
      enddo
!C*** Show Stress
      if( ndof==2 ) mdof = 3
      if( ndof==3 ) mdof = 6
      if( ndof==6 ) mdof = 6
!C @node
      do i = 1, hecMESH%nn_internal
        if(fstrSOLID%is_rot(i)==1)cycle
        j = hecMESH%global_node_ID(i)
        if( i==1 ) then
          do k = 1, mdof
            Smax(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
            Smin(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
            ISmax(k)= j
            ISmin(k)= j
          enddo
          Mmax(1) = fstrSOLID%MISES(1)
          Mmin(1) = fstrSOLID%MISES(1)
          IMmax(1)= j
          IMmin(1)= j
        else
          do k = 1, mdof
            if( fstrSOLID%STRESS(mdof*(i-1)+k) > Smax(k) ) then
              Smax(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
              ISmax(k)= j
            endif
            if( fstrSOLID%STRESS(mdof*(i-1)+k) < Smin(k) ) then
              Smin(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
              ISmin(k)= j
            endif
          enddo
          if( fstrSOLID%MISES(i) > Mmax(1) ) then
            Mmax(1) = fstrSOLID%MISES(i)
            IMmax(1)= j
          endif
          if( fstrSOLID%MISES(i) < Mmin(1) ) then
            Mmin(1) = fstrSOLID%MISES(i)
            IMmin(1)= j
          endif
        endif
      enddo
!C @element
      do i = 1, hecMESH%n_elem
        ID_area = hecMESH%elem_ID(i*2)
        if( ID_area==hecMESH%my_rank ) then
          j = hecMESH%global_elem_ID(i)
          if( hecMESH%elem_ID(i*2-1)==1 ) then
            do k = 1, mdof
              ESmax(k) = fstrSOLID%ESTRESS(mdof*(i-1)+k)
              ESmin(k) = fstrSOLID%ESTRESS(mdof*(i-1)+k)
              IESmax(k)= j
              IESmin(k)= j
            enddo
            EMmax(1) = fstrSOLID%EMISES(1)
            EMmin(1) = fstrSOLID%EMISES(1)
            IEMmax(1)= j
            IEMmin(1)= j
          else
            do k = 1, mdof
              if( fstrSOLID%ESTRESS(mdof*(i-1)+k) > ESmax(k) ) then
                ESmax(k) = fstrSOLID%ESTRESS(mdof*(i-1)+k)
                IESmax(k)= j
              endif
              if( fstrSOLID%ESTRESS(mdof*(i-1)+k) < ESmin(k) ) then
                ESmin(k) = fstrSOLID%ESTRESS(mdof*(i-1)+k)
                IESmin(k)= j
              endif
            enddo
            if( fstrSOLID%EMISES(i) > EMmax(1) ) then
              EMmax(1) = fstrSOLID%EMISES(i)
              IEMmax(1)= j
            endif
            if( fstrSOLID%EMISES(i) < EMmin(1) ) then
              EMmin(1) = fstrSOLID%EMISES(i)
              IEMmin(1)= j
            endif
          endif
        endif
      enddo

      if( ndof==2 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//U1 ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009) '//U2 ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009) '//V1 ',Vmax(1),IVmax(1),Vmin(1),IVmin(1)
        write(ILOG,1009) '//V2 ',Vmax(2),IVmax(2),Vmin(2),IVmin(2)
        write(ILOG,1009) '//A1 ',Amax(1),IAmax(1),Amin(1),IAmin(1)
        write(ILOG,1009) '//A2 ',Amax(2),IAmax(2),Amin(2),IAmin(2)
        write(ILOG,1009) '//E11',Emax(1),IEmax(1),Emin(1),IEmin(1)
        write(ILOG,1009) '//E22',Emax(2),IEmax(2),Emin(2),IEmin(2)
        write(ILOG,1009) '//E12',Emax(3),IEmax(3),Emin(3),IEmin(3)
        write(ILOG,1009) '//S11',Smax(1),ISmax(1),Smin(1),ISmin(1)
        write(ILOG,1009) '//S22',Smax(2),ISmax(2),Smin(2),ISmin(2)
        write(ILOG,1009) '//S12',Smax(3),ISmax(3),Smin(3),ISmin(3)
        write(ILOG,1009) '//SMS',Mmax(1),IMmax(1),Mmin(1),IMmin(1)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//E11',EEmax(1),IEEmax(1),EEmin(1),IEEmin(1)
        write(ILOG,1009) '//E22',EEmax(2),IEEmax(2),EEmin(2),IEEmin(2)
        write(ILOG,1009) '//E12',EEmax(3),IEEmax(3),EEmin(3),IEEmin(3)
        write(ILOG,1009) '//S11',ESmax(1),IESmax(1),ESmin(1),IESmin(1)
        write(ILOG,1009) '//S22',ESmax(2),IESmax(2),ESmin(2),IESmin(2)
        write(ILOG,1009) '//S12',ESmax(3),IESmax(3),ESmin(3),IESmin(3)
        write(ILOG,1009) '//SMS',EMmax(1),IEMmax(1),EMmin(1),IEMmin(1)
!C*** Show Summary
        call hecmw_allREDUCE_R(hecMESH,Umax,2,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin,2,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Vmax,2,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Vmin,2,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Amax,2,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Amin,2,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,EEmax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,EEmin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,ESmax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,ESmin,3,hecmw_min)
        if( hecMESH%my_rank==0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1 ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2 ',Umax(2),Umin(2)
          write(ILOG,1019) '//V1 ',Vmax(1),Vmin(1)
          write(ILOG,1019) '//V2 ',Vmax(2),Vmin(2)
          write(ILOG,1019) '//A1 ',Amax(1),Amin(1)
          write(ILOG,1019) '//A2 ',Amax(2),Amin(2)
          write(ILOG,1019) '//E11',Emax(1),Emin(1)
          write(ILOG,1019) '//E22',Emax(2),Emin(2)
          write(ILOG,1019) '//E12',Emax(3),Emin(3)
          write(ILOG,1019) '//S11',Smax(1),Smin(1)
          write(ILOG,1019) '//S22',Smax(2),Smin(2)
          write(ILOG,1019) '//S12',Smax(3),Smin(3)
          write(ILOG,1019) '//SMS',Mmax(1),Mmin(1)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11',EEmax(1),EEmin(1)
          write(ILOG,1019) '//E22',EEmax(2),EEmin(2)
          write(ILOG,1019) '//E12',EEmax(3),EEmin(3)
          write(ILOG,1019) '//S11',ESmax(1),ESmin(1)
          write(ILOG,1019) '//S22',ESmax(2),ESmin(2)
          write(ILOG,1019) '//S12',ESmax(3),ESmin(3)
          write(ILOG,1019) '//SMS',EMmax(1),EMmin(1)
        endif
      else if( ndof==3 .or. ndof==6 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//U1 ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009) '//U2 ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009) '//U3 ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,1009) '//V1 ',Vmax(1),IVmax(1),Vmin(1),IVmin(1)
        write(ILOG,1009) '//V2 ',Vmax(2),IVmax(2),Vmin(2),IVmin(2)
        write(ILOG,1009) '//V3 ',Vmax(3),IVmax(3),Vmin(3),IVmin(3)
        write(ILOG,1009) '//A1 ',Amax(1),IAmax(1),Amin(1),IAmin(1)
        write(ILOG,1009) '//A2 ',Amax(2),IAmax(2),Amin(2),IAmin(2)
        write(ILOG,1009) '//A3 ',Amax(3),IAmax(3),Amin(3),IAmin(3)
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
        write(ILOG,1009) '//SMS',Mmax(1),IMmax(1),Mmin(1),IMmin(1)
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
        write(ILOG,1009) '//SMS',EMmax(1),IEMmax(1),EMmin(1),IEMmin(1)
!C*** Show Summary
        call hecmw_allREDUCE_R(hecMESH,Umax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Vmax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Vmin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Amax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Amin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,EEmax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,EEmin,6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,ESmax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,ESmin,6,hecmw_min)
        if( hecMESH%my_rank==0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1 ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2 ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3 ',Umax(3),Umin(3)
          write(ILOG,1019) '//V1 ',Vmax(1),Vmin(1)
          write(ILOG,1019) '//V2 ',Vmax(2),Vmin(2)
          write(ILOG,1019) '//V3 ',Vmax(3),Vmin(3)
          write(ILOG,1019) '//A1 ',Amax(1),Amin(1)
          write(ILOG,1019) '//A2 ',Amax(2),Amin(2)
          write(ILOG,1019) '//A3 ',Amax(3),Amin(3)
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
          write(ILOG,1019) '//SMS',Mmax(1),Mmin(1)
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
          write(ILOG,1019) '//SMS',EMmax(1),EMmin(1)
        endif
      else if( 1 == 0 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//U1    ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009)'//U2    ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009)'//U3    ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,1009)'//V1    ',Vmax(1),IVmax(1),Vmin(1),IVmin(1)
        write(ILOG,1009)'//V2    ',Vmax(2),IVmax(2),Vmin(2),IVmin(2)
        write(ILOG,1009)'//V3    ',Vmax(3),IVmax(3),Vmin(3),IVmin(3)
        write(ILOG,1009)'//A1    ',Amax(1),IAmax(1),Amin(1),IAmin(1)
        write(ILOG,1009)'//A2    ',Amax(2),IAmax(2),Amin(2),IAmin(2)
        write(ILOG,1009)'//A3    ',Amax(3),IAmax(3),Amin(3),IAmin(3)
        write(ILOG,1009)'//R1    ',Umax(4),IUmax(4),Umin(4),IUmin(4)
        write(ILOG,1009)'//R2    ',Umax(5),IUmax(5),Umin(5),IUmin(5)
        write(ILOG,1009)'//R3    ',Umax(6),IUmax(6),Umin(6),IUmin(6)
        write(ILOG,1009)'//RV1   ',Vmax(4),IVmax(4),Vmin(4),IVmin(4)
        write(ILOG,1009)'//RV2   ',Vmax(5),IVmax(5),Vmin(5),IVmin(5)
        write(ILOG,1009)'//RV3   ',Vmax(6),IVmax(6),Vmin(6),IVmin(6)
        write(ILOG,1009)'//RA1   ',Amax(4),IAmax(4),Amin(4),IAmin(4)
        write(ILOG,1009)'//RA2   ',Amax(5),IAmax(5),Amin(5),IAmin(5)
        write(ILOG,1009)'//RA3   ',Amax(6),IAmax(6),Amin(6),IAmin(6)
        write(ILOG,1009)'//E11(+)',Emax( 1),IEmax( 1),Emin( 1),IEmin( 1)
        write(ILOG,1009)'//E22(+)',Emax( 2),IEmax( 2),Emin( 2),IEmin( 2)
        write(ILOG,1009)'//E33(+)',Emax( 3),IEmax( 3),Emin( 3),IEmin( 3)
        write(ILOG,1009)'//E12(+)',Emax( 4),IEmax( 4),Emin( 4),IEmin( 4)
        write(ILOG,1009)'//E23(+)',Emax( 5),IEmax( 5),Emin( 5),IEmin( 5)
        write(ILOG,1009)'//E31(+)',Emax( 6),IEmax( 6),Emin( 6),IEmin( 6)
        write(ILOG,1009)'//E11(-)',Emax( 7),IEmax( 7),Emin( 7),IEmin( 7)
        write(ILOG,1009)'//E22(-)',Emax( 8),IEmax( 8),Emin( 8),IEmin( 8)
        write(ILOG,1009)'//E33(-)',Emax( 9),IEmax( 9),Emin( 9),IEmin( 9)
        write(ILOG,1009)'//E12(-)',Emax(10),IEmax(10),Emin(10),IEmin(10)
        write(ILOG,1009)'//E23(-)',Emax(11),IEmax(11),Emin(11),IEmin(11)
        write(ILOG,1009)'//E31(-)',Emax(12),IEmax(12),Emin(12),IEmin(12)
        write(ILOG,1009)'//S11(+)',Smax( 1),ISmax( 1),Smin( 1),ISmin( 1)
        write(ILOG,1009)'//S22(+)',Smax( 2),ISmax( 2),Smin( 2),ISmin( 2)
        write(ILOG,1009)'//S33(+)',Smax( 3),ISmax( 3),Smin( 3),ISmin( 3)
        write(ILOG,1009)'//S12(+)',Smax( 4),ISmax( 4),Smin( 4),ISmin( 4)
        write(ILOG,1009)'//S23(+)',Smax( 5),ISmax( 5),Smin( 5),ISmin( 5)
        write(ILOG,1009)'//S31(+)',Smax( 6),ISmax( 6),Smin( 6),ISmin( 6)
        write(ILOG,1009)'//S11(-)',Smax( 7),ISmax( 7),Smin( 7),ISmin( 7)
        write(ILOG,1009)'//S22(-)',Smax( 8),ISmax( 8),Smin( 8),ISmin( 8)
        write(ILOG,1009)'//S33(-)',Smax( 9),ISmax( 9),Smin( 9),ISmin( 9)
        write(ILOG,1009)'//S12(-)',Smax(10),ISmax(10),Smin(10),ISmin(10)
        write(ILOG,1009)'//S23(-)',Smax(11),ISmax(11),Smin(11),ISmin(11)
        write(ILOG,1009)'//S31(-)',Smax(12),ISmax(12),Smin(12),ISmin(12)
        write(ILOG,1009)'//SMS(+)',Mmax( 1),IMmax( 1),Mmin( 1),IMmin( 1)
        write(ILOG,1009)'//SMS(-)',Mmax( 1),IMmax( 1),Mmin( 1),IMmin( 1)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//E11(+)',EEmax( 1),IEEmax( 1),EEmin( 1),IEEmin( 1)
        write(ILOG,1009)'//E22(+)',EEmax( 2),IEEmax( 2),EEmin( 2),IEEmin( 2)
        write(ILOG,1009)'//E33(+)',EEmax( 3),IEEmax( 3),EEmin( 3),IEEmin( 3)
        write(ILOG,1009)'//E12(+)',EEmax( 4),IEEmax( 4),EEmin( 4),IEEmin( 4)
        write(ILOG,1009)'//E23(+)',EEmax( 5),IEEmax( 5),EEmin( 5),IEEmin( 5)
        write(ILOG,1009)'//E31(+)',EEmax( 6),IEEmax( 6),EEmin( 6),IEEmin( 6)
        write(ILOG,1009)'//E11(-)',EEmax( 7),IEEmax( 7),EEmin( 7),IEEmin( 7)
        write(ILOG,1009)'//E22(-)',EEmax( 8),IEEmax( 8),EEmin( 8),IEEmin( 8)
        write(ILOG,1009)'//E33(-)',EEmax( 9),IEEmax( 9),EEmin( 9),IEEmin( 9)
        write(ILOG,1009)'//E12(-)',EEmax(10),IEEmax(10),EEmin(10),IEEmin(10)
        write(ILOG,1009)'//E23(-)',EEmax(11),IEEmax(11),EEmin(11),IEEmin(11)
        write(ILOG,1009)'//E31(-)',EEmax(12),IEEmax(12),EEmin(12),IEEmin(12)
        write(ILOG,1009)'//S11(+)',ESmax( 1),IESmax( 1),ESmin( 1),IESmin( 1)
        write(ILOG,1009)'//S22(+)',ESmax( 2),IESmax( 2),ESmin( 2),IESmin( 2)
        write(ILOG,1009)'//S33(+)',ESmax( 3),IESmax( 3),ESmin( 3),IESmin( 3)
        write(ILOG,1009)'//S12(+)',ESmax( 4),IESmax( 4),ESmin( 4),IESmin( 4)
        write(ILOG,1009)'//S23(+)',ESmax( 5),IESmax( 5),ESmin( 5),IESmin( 5)
        write(ILOG,1009)'//S31(+)',ESmax( 6),IESmax( 6),ESmin( 6),IESmin( 6)
        write(ILOG,1009)'//S11(-)',ESmax( 7),IESmax( 7),ESmin( 7),IESmin( 7)
        write(ILOG,1009)'//S22(-)',ESmax( 8),IESmax( 8),ESmin( 8),IESmin( 8)
        write(ILOG,1009)'//S33(-)',ESmax( 9),IESmax( 9),ESmin( 9),IESmin( 9)
        write(ILOG,1009)'//S12(-)',ESmax(10),IESmax(10),ESmin(10),IESmin(10)
        write(ILOG,1009)'//S23(-)',ESmax(11),IESmax(11),ESmin(11),IESmin(11)
        write(ILOG,1009)'//S31(-)',ESmax(12),IESmax(12),ESmin(12),IESmin(12)
        write(ILOG,1009)'//SMS(+)',EMmax( 1),IEMmax( 1),EMmin( 1),IEMmin( 1)
        write(ILOG,1009)'//SMS(-)',EMmax( 1),IEMmax( 1),EMmin( 1),IEMmin( 1)
!C*** Show Summary
        call hecmw_allREDUCE_R(hecMESH,Umax, 6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin, 6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Vmax, 6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Vmin, 6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Amax, 6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Amin, 6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,12,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,12,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,12,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,12,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,EEmax,12,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,EEmin,12,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,ESmax,12,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,ESmin,12,hecmw_min)
        if( hecMESH%my_rank==0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1    ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2    ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3    ',Umax(3),Umin(3)
          write(ILOG,1019) '//V1    ',Vmax(1),Vmin(1)
          write(ILOG,1019) '//V2    ',Vmax(2),Vmin(2)
          write(ILOG,1019) '//V3    ',Vmax(3),Vmin(3)
          write(ILOG,1019) '//A1    ',Amax(1),Amin(1)
          write(ILOG,1019) '//A2    ',Amax(2),Amin(2)
          write(ILOG,1019) '//A3    ',Amax(3),Amin(3)
          write(ILOG,1019) '//R1    ',Umax(4),Umin(4)
          write(ILOG,1019) '//R2    ',Umax(5),Umin(5)
          write(ILOG,1019) '//R3    ',Umax(6),Umin(6)
          write(ILOG,1019) '//RV1   ',Vmax(4),Vmin(4)
          write(ILOG,1019) '//RV2   ',Vmax(5),Vmin(5)
          write(ILOG,1019) '//RV3   ',Vmax(6),Vmin(6)
          write(ILOG,1019) '//RA1   ',Amax(4),Amin(4)
          write(ILOG,1019) '//RA2   ',Amax(5),Amin(5)
          write(ILOG,1019) '//RA3   ',Amax(6),Amin(6)
          write(ILOG,1019) '//E11(+)',Emax( 1),Emin( 1)
          write(ILOG,1019) '//E22(+)',Emax( 2),Emin( 2)
          write(ILOG,1019) '//E33(+)',Emax( 3),Emin( 3)
          write(ILOG,1019) '//E12(+)',Emax( 4),Emin( 4)
          write(ILOG,1019) '//E23(+)',Emax( 5),Emin( 5)
          write(ILOG,1019) '//E31(+)',Emax( 6),Emin( 6)
          write(ILOG,1019) '//E11(-)',Emax( 7),Emin( 7)
          write(ILOG,1019) '//E22(-)',Emax( 8),Emin( 8)
          write(ILOG,1019) '//E33(-)',Emax( 9),Emin( 9)
          write(ILOG,1019) '//E12(-)',Emax(10),Emin(10)
          write(ILOG,1019) '//E23(-)',Emax(11),Emin(11)
          write(ILOG,1019) '//E31(-)',Emax(12),Emin(12)
          write(ILOG,1019) '//S11(+)',Smax( 1),Smin( 1)
          write(ILOG,1019) '//S22(+)',Smax( 2),Smin( 2)
          write(ILOG,1019) '//S33(+)',Smax( 3),Smin( 3)
          write(ILOG,1019) '//S12(+)',Smax( 4),Smin( 4)
          write(ILOG,1019) '//S23(+)',Smax( 5),Smin( 5)
          write(ILOG,1019) '//S31(+)',Smax( 6),Smin( 6)
          write(ILOG,1019) '//S11(-)',Smax( 7),Smin( 7)
          write(ILOG,1019) '//S22(-)',Smax( 8),Smin( 8)
          write(ILOG,1019) '//S33(-)',Smax( 9),Smin( 9)
          write(ILOG,1019) '//S12(-)',Smax(10),Smin(10)
          write(ILOG,1019) '//S23(-)',Smax(11),Smin(11)
          write(ILOG,1019) '//S31(-)',Smax(12),Smin(12)
          write(ILOG,1019) '//SMS(+)',Mmax( 1),Mmin( 1)
          write(ILOG,1019) '//SMS(-)',Mmax( 1),Mmin( 1)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11(+)',EEmax( 1),EEmin( 1)
          write(ILOG,1019) '//E22(+)',EEmax( 2),EEmin( 2)
          write(ILOG,1019) '//E33(+)',EEmax( 3),EEmin( 3)
          write(ILOG,1019) '//E12(+)',EEmax( 4),EEmin( 4)
          write(ILOG,1019) '//E23(+)',EEmax( 5),EEmin( 5)
          write(ILOG,1019) '//E31(+)',EEmax( 6),EEmin( 6)
          write(ILOG,1019) '//E11(-)',EEmax( 7),EEmin( 7)
          write(ILOG,1019) '//E22(-)',EEmax( 8),EEmin( 8)
          write(ILOG,1019) '//E33(-)',EEmax( 9),EEmin( 9)
          write(ILOG,1019) '//E12(-)',EEmax(10),EEmin(10)
          write(ILOG,1019) '//E23(-)',EEmax(11),EEmin(11)
          write(ILOG,1019) '//E31(-)',EEmax(12),EEmin(12)
          write(ILOG,1019) '//S11(+)',ESmax( 1),ESmin( 1)
          write(ILOG,1019) '//S22(+)',ESmax( 2),ESmin( 2)
          write(ILOG,1019) '//S33(+)',ESmax( 3),ESmin( 3)
          write(ILOG,1019) '//S12(+)',ESmax( 4),ESmin( 4)
          write(ILOG,1019) '//S23(+)',ESmax( 5),ESmin( 5)
          write(ILOG,1019) '//S31(+)',ESmax( 6),ESmin( 6)
          write(ILOG,1019) '//S11(-)',ESmax( 7),ESmin( 7)
          write(ILOG,1019) '//S22(-)',ESmax( 8),ESmin( 8)
          write(ILOG,1019) '//S33(-)',ESmax( 9),ESmin( 9)
          write(ILOG,1019) '//S12(-)',ESmax(10),ESmin(10)
          write(ILOG,1019) '//S23(-)',ESmax(11),ESmin(11)
          write(ILOG,1019) '//S31(-)',ESmax(12),ESmin(12)
          write(ILOG,1019) '//SMS(+)',EMmax( 1),EMmin( 1)
          write(ILOG,1019) '//SMS(-)',EMmax( 1),EMmin( 1)
        endif
      endif

 1009 format(a8,1pe12.4,i10,1pe12.4,i10)
 1019 format(a8,1pe12.4,1pe12.4)
  end subroutine fstr_dynamic_post

!C================================================================C
!C-- subroutine dynamic_output_monit
!C================================================================C
  subroutine dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, fstrSOLID)
    use m_fstr
    type ( hecmwST_local_mesh  ) :: hecMESH
    type ( fstr_param          ) :: fstrPARAM
    type ( fstr_dynamic        ) :: fstrDYNAMIC
    type ( lczparam            ) :: myEIG
    type ( fstr_solid          ) :: fstrSOLID

    integer(kind=kint) :: idx, ii, jj, ierr, ncmp
    integer(kind=kint) :: num_monit, ig, iS, iE, ik, iunitS, iunit
    logical :: yes

    if( mod(fstrDYNAMIC%i_step,fstrDYNAMIC%nout_monit)/=0 ) return

    if( fstrDYNAMIC%idx_eqa==1 .and. fstrDYNAMIC%i_step>0 ) then
      idx = 2
    else
      idx = 1
    endif

    num_monit = 0
    ig = fstrDYNAMIC%ngrp_monit
    iS = hecMESH%node_group%grp_index(ig-1)+1
    iE = hecMESH%node_group%grp_index(ig)
    do ik=iS,iE
      ii = hecMESH%node_group%grp_item(ik)
      if (ii > hecMESH%nn_internal) cycle
      num_monit = num_monit+1
      jj = hecMESH%global_node_id(ii)
      iunitS = 10*(num_monit-1)

!C-- displacement
      if( fstrDYNAMIC%iout_list(1)==1 ) then
        iunit = iunitS + fstrDYNAMIC%dynamic_IW4
        write( iunit, '(i10,1pe13.4e3,i10,1p6e13.4e3)' ) &
               fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, jj, &
               fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , idx )
        end if
!C-- velocity
      if( fstrDYNAMIC%iout_list(2)==1 ) then
        iunit = iunitS + fstrDYNAMIC%dynamic_IW5
        write( iunit, '(i10,1pe13.4e3,i10,1p6e13.4e3)' ) &
               fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, jj, &
               fstrDYNAMIC%VEL( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , idx )
        end if
!C-- acceleration
      if( fstrDYNAMIC%iout_list(3)==1 ) then
        iunit = iunitS + fstrDYNAMIC%dynamic_IW6
        write( iunit, '(i10,1pe13.4e3,i10,1p6e13.4e3)' ) &
               fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, jj, &
               fstrDYNAMIC%ACC( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , idx )
        end if

!C-- strain
      if( fstrDYNAMIC%iout_list(5) > 0 ) then
        if (hecMESH%n_dof == 3 .or. hecMESH%n_dof == 2) then
          ncmp = 6
        else
          ncmp = 12
        endif
        iunit = iunitS + fstrDYNAMIC%dynamic_IW8
        write( iunit, '(i10,1pe13.4e3,i10,1p6e13.4e3)') &
               fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, jj, &
               fstrSOLID%STRAIN( ncmp*(ii-1)+1 : ncmp*ii )
      end if
!C-- stress
      if( fstrDYNAMIC%iout_list(6) > 0 ) then
        if (hecMESH%n_dof == 3 .or. hecMESH%n_dof == 2) then
          ncmp = 6
        else
          ncmp = 12
        endif
        iunit = iunitS + fstrDYNAMIC%dynamic_IW9
        write( iunit, '(i10,1pe13.4e3,i10,1p7e13.4e3)') &
               fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, jj, &
               fstrSOLID%STRESS( ncmp*(ii-1)+1 : ncmp*ii )
      end if
    enddo

    if( hecMESH%my_rank==0 ) then
      if( any(fstrDYNAMIC%iout_list(1:3)==1) ) then
        inquire( file='dyna_energy.txt', opened=yes )
        if( .not. yes ) then
          open( fstrDYNAMIC%dynamic_IW7, file='dyna_energy.txt', status='replace', iostat=ierr )
          if( ierr/=0 ) then
            write(*,*) 'stop due to file opening error <dyna_enrgy.txt>'
            call hecmw_abort( hecmw_comm_get_comm() )
          endif
          write( fstrDYNAMIC%dynamic_IW7, * ) &
                 ' time step', '     time    ', '  kinetic energy', '   strain energy', '   total energy'
        endif
        if(fstrDYNAMIC%i_step==0) then
           fstrDYNAMIC%kineticEnergy = 0.0d0
           do ii = 1, hecMESH%n_node*hecMESH%n_dof
             fstrDYNAMIC%kineticEnergy = fstrDYNAMIC%kineticEnergy &
                        + 0.5d0 * myEIG%mass(ii) * fstrDYNAMIC%VEL(ii,idx) * fstrDYNAMIC%VEL(ii,idx)
           enddo
        endif
        fstrDYNAMIC%totalEnergy = fstrDYNAMIC%kineticEnergy + fstrDYNAMIC%strainEnergy
        write( fstrDYNAMIC%dynamic_IW7, '(i10,1pe13.4e3,1p3e16.4e3)' ) &
               fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, &
               fstrDYNAMIC%kineticEnergy, fstrDYNAMIC%strainEnergy, fstrDYNAMIC%totalEnergy
      endif
      if( fstrDYNAMIC%i_step==fstrDYNAMIC%n_step ) close(fstrDYNAMIC%dynamic_IW7)
    endif
  end subroutine dynamic_output_monit

!C================================================================C
!C-- subroutine  matvec
!C================================================================C
  subroutine matvec(y,x,hecMAT,ndof,D,AU,AL)
    use m_fstr
    type (hecmwST_matrix) :: hecMAT

    integer(kind=kint) :: ndof,i,is,ie,j,icol,ki,kj,ix,iy,ip,nn
    real(kind=kreal) :: D(ndof*ndof*hecMAT%NP)
    real(kind=kreal) :: AU(ndof*ndof*hecMAT%NPU)
    real(kind=kreal) :: AL(ndof*ndof*hecMAT%NPL)
    real(kind=kreal) :: x(ndof*hecMAT%NP)
    real(kind=kreal) :: y(ndof*hecMAT%NP)

    nn=ndof*ndof

    y=0.0d0

    do i=1,hecMAT%NP
      is=hecMAT%indexU(i-1)+1
      ie=hecMAT%indexU(i)
      do j=is,ie
        icol=hecMAT%itemU(j)
        do ki=1,ndof
          iy=ndof*(i-1)+ki
          do kj=1,ndof
            ix=ndof*(icol-1)+kj
            ip=nn*(j-1)+ndof*(ki-1)+kj
            y(iy)=y(iy)+AU(ip)*x(ix)
          enddo
        enddo
      enddo
    enddo

    do i=1,hecMAT%NP
      is=hecMAT%indexL(i-1)+1
      ie=hecMAT%indexL(i)
      do j=is,ie
        icol=hecMAT%itemL(j)
        do ki=1,ndof
          iy=ndof*(i-1)+ki
          do kj=1,ndof
            ix=ndof*(icol-1)+kj
            ip=nn*(j-1)+ndof*(ki-1)+kj
            y(iy)=y(iy)+AL(ip)*x(ix)
          enddo
        enddo
      enddo
    enddo

    do i=1,hecMAT%NP
      do ki=1,ndof
        iy=ndof*(i-1)+ki
        do kj=1,ndof
          ix=ndof*(i-1)+kj
          ip=nn*(i-1)+ndof*(ki-1)+kj
          y(iy)=y(iy)+D(ip)*x(ix)
        enddo
      enddo
    enddo
  end subroutine matvec

end module m_dynamic_output
