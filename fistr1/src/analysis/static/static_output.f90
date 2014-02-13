!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Suemitsu(Advancesoft)                       !
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
module m_static_output
	implicit none

	contains

!> Output result
!----------------------------------------------------------------------*
  subroutine fstr_static_Output( cstep, istep, hecMESH, fstrSOLID, flag )
!----------------------------------------------------------------------*
    use m_fstr
    use m_fstr_NodalStress
    use m_static_make_result
    use m_hecmw2fstr_mesh_conv
    integer, intent(in)                   :: cstep       !< current step number
    integer, intent(in)                   :: istep       !< current substep number
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(inout)      :: fstrSOLID
    integer, intent(in)                   :: flag

    type ( hecmwST_result_data ) :: fstrRESULT
    integer(kind=kint) :: i, j, ndof, maxstep, interval, fnum, iS, iE, gid
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)

    ndof = hecMESH%n_dof

    nullify( tnstrain )
    nullify( testrain )
	
    if( fstrSOLID%TEMP_ngrp_tot>0 .or. fstrSOLID%TEMP_irres>0 ) then
       if( ndof==3 ) then
          allocate( tnstrain(hecMESH%n_node*6) )
          allocate( testrain(hecMESH%n_elem*6) )
       else if( ndof==2 ) then
          allocate( tnstrain(hecMESH%n_node*3) )
          allocate( testrain(hecMESH%n_elem*3) )
       endif
    endif
	
    if( ndof==3 ) then
          call fstr_NodalStress3D( hecMESH, fstrSOLID, tnstrain, testrain )
    else if( ndof==2 ) then
          call fstr_NodalStress2D( hecMESH, fstrSOLID, tnstrain, testrain )
    else if( ndof==6 ) then
          call fstr_NodalStress6D( hecMESH, fstrSOLID )
    endif

    if( fstrSOLID%TEMP_irres>0 ) then
          maxstep = fstrSOLID%TEMP_irres
    else
          maxstep = 0
          do i = 1, cstep
            maxstep = maxstep + fstrSOLID%step_ctrl(i)%num_substep
          end do
    endif

    if( flag==6 ) then
      if( IRESULT==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(3)%freqency)==0 .or. istep==maxstep) ) then
          call fstr_write_static_result( hecMESH, fstrSOLID, maxstep, istep, tnstrain, testrain, 1 )
      endif
      if( associated(tnstrain) ) deallocate( tnstrain )
      if( associated(testrain) ) deallocate( testrain )
      return
    endif

    if( IRESULT==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(3)%freqency)==0 .or. istep==maxstep) ) then
          call fstr_write_static_result( hecMESH, fstrSOLID, maxstep, istep, tnstrain, testrain, 0 )
    endif

    if( IVISUAL==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(4)%freqency)==0 .or. istep==maxstep) ) then
          interval = fstrSOLID%output_ctrl(4)%freqency

          call fstr_make_static_result( hecMESH, fstrSOLID, fstrRESULT, tnstrain, testrain )
          call fstr2hecmw_mesh_conv( hecMESH )
          call hecmw_visualize_init
          call hecmw_visualize( hecMESH, fstrRESULT, istep, maxstep, interval )
          call hecmw_visualize_finalize
          call hecmw2fstr_mesh_conv( hecMESH )
          call hecmw_result_free( fstrRESULT )
    endif

    if( (mod(istep,fstrSOLID%output_ctrl(1)%freqency)==0 .or. istep==maxstep) ) then
          fnum = fstrSOLID%output_ctrl(1)%filenum
          call fstr_static_post( fnum, hecMESH, fstrSOLID, istep )
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

    if( associated(tnstrain) ) deallocate( tnstrain )
    if( associated(testrain) ) deallocate( testrain )

  end subroutine fstr_static_Output

!> Summarizer of output data which prints out max and min output values
!----------------------------------------------------------------------*
  subroutine fstr_static_post( fnum, hecMESH, fstrSOLID, istep )
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(in)                   :: fnum, istep
      type (hecmwST_local_mesh), intent(in) :: hecMESH
      type (fstr_solid), intent(in)         :: fstrSOLID

      real(kind=kreal)   :: Umax(6), Umin(6), Emax(14), Emin(14), Smax(14), Smin(14)
      integer(kind=kint) :: IUmax(6), IUmin(6), IEmax(14), IEmin(14), ISmax(14), ISmin(14)
      real(kind=kreal)   :: EEmax(14), EEmin(14), ESmax(14), ESmin(14)
      integer(kind=kint) :: IEEmax(14), IEEmin(14), IESmax(14), IESmin(14)
      integer(kind=kint) :: i, j, k, ndof, mdof, ID_area

      ndof = hecMESH%n_dof

      write( fnum, '(''#### Result step='',I6)') istep

!C*** Show Displacement
      do i = 1, hecMESH%nn_internal
        j = hecMESH%global_node_ID(i)
        if( i==1 ) then
          do k = 1, ndof 
            Umax(k) = fstrSOLID%unode(ndof*(i-1)+k)
            Umin(k) = fstrSOLID%unode(ndof*(i-1)+k)
            IUmax(k)= j
            IUmin(k)= j
          enddo
        else
          do k = 1, ndof 
            if( fstrSOLID%unode(ndof*(i-1)+k) > Umax(k) ) then
              Umax(k) = fstrSOLID%unode(ndof*(i-1)+k)
              IUmax(k)= j
            endif
            if( fstrSOLID%unode(ndof*(i-1)+k) < Umin(k) ) then
              Umin(k) = fstrSOLID%unode(ndof*(i-1)+k)
              IUmin(k)= j
            endif
          enddo
        endif
      enddo
!C*** Show Strain
      if( ndof==2 ) mdof = 3
      if( ndof==3 ) mdof = 6
      if( ndof==6 ) mdof = 12
!C @node
      do i = 1, hecMESH%nn_internal
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
      if( ndof==2 ) mdof = 4
      if( ndof==3 ) mdof = 7
      if( ndof==6 ) mdof = 14
!C @node
      do i = 1, hecMESH%nn_internal
        j = hecMESH%global_node_ID(i)
        if( i==1 ) then
          do k = 1, mdof
            Smax(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
            Smin(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
            ISmax(k)= j
            ISmin(k)= j
          enddo
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
          endif
        endif
      enddo

      if( ndof==2 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//U1 ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009) '//U2 ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009) '//E11',Emax(1),IEmax(1),Emin(1),IEmin(1)
        write(ILOG,1009) '//E22',Emax(2),IEmax(2),Emin(2),IEmin(2)
        write(ILOG,1009) '//E12',Emax(3),IEmax(3),Emin(3),IEmin(3)
        write(ILOG,1009) '//S11',Smax(1),ISmax(1),Smin(1),ISmin(1)
        write(ILOG,1009) '//S22',Smax(2),ISmax(2),Smin(2),ISmin(2)
        write(ILOG,1009) '//S12',Smax(3),ISmax(3),Smin(3),ISmin(3)
        write(ILOG,1009) '//SMS',Smax(4),ISmax(4),Smin(4),ISmin(4)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//E11',EEmax(1),IEEmax(1),EEmin(1),IEEmin(1)
        write(ILOG,1009) '//E22',EEmax(2),IEEmax(2),EEmin(2),IEEmin(2)
        write(ILOG,1009) '//E12',EEmax(3),IEEmax(3),EEmin(3),IEEmin(3)
        write(ILOG,1009) '//S11',ESmax(1),IESmax(1),ESmin(1),IESmin(1)
        write(ILOG,1009) '//S22',ESmax(2),IESmax(2),ESmin(2),IESmin(2)
        write(ILOG,1009) '//S12',ESmax(3),IESmax(3),ESmin(3),IESmin(3)
        write(ILOG,1009) '//SMS',ESmax(4),IESmax(4),ESmin(4),IESmin(4)
!C*** Show Summary
        call hecmw_allREDUCE_R(hecMESH,Umax,2,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin,2,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,4,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,4,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,EEmax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,EEmin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,ESmax,4,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,ESmin,4,hecmw_min)
        if( hecMESH%my_rank==0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1 ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2 ',Umax(2),Umin(2)
          write(ILOG,1019) '//E11',Emax(1),Emin(1)
          write(ILOG,1019) '//E22',Emax(2),Emin(2)
          write(ILOG,1019) '//E12',Emax(3),Emin(3)
          write(ILOG,1019) '//S11',Smax(1),Smin(1)
          write(ILOG,1019) '//S22',Smax(2),Smin(2)
          write(ILOG,1019) '//S12',Smax(3),Smin(3)
          write(ILOG,1019) '//SMS',Smax(4),Smin(4)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11',EEmax(1),EEmin(1)
          write(ILOG,1019) '//E22',EEmax(2),EEmin(2)
          write(ILOG,1019) '//E12',EEmax(3),EEmin(3)
          write(ILOG,1019) '//S11',ESmax(1),ESmin(1)
          write(ILOG,1019) '//S22',ESmax(2),ESmin(2)
          write(ILOG,1019) '//S12',ESmax(3),ESmin(3)
          write(ILOG,1019) '//SMS',ESmax(4),ESmin(4)
        endif
      else if( ndof==3 ) then
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
        call hecmw_allREDUCE_R(hecMESH,Umax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,7,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,7,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,EEmax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,EEmin,6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,ESmax,7,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,ESmin,7,hecmw_min)
        if( hecMESH%my_rank==0 ) then
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
      else if( ndof==6 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//U1    ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009)'//U2    ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009)'//U3    ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,1009)'//R1    ',Umax(4),IUmax(4),Umin(4),IUmin(4)
        write(ILOG,1009)'//R2    ',Umax(5),IUmax(5),Umin(5),IUmin(5)
        write(ILOG,1009)'//R3    ',Umax(6),IUmax(6),Umin(6),IUmin(6)
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
        write(ILOG,1009)'//SMS(+)',Smax(13),ISmax(13),Smin(13),ISmin(13)
        write(ILOG,1009)'//SMS(-)',Smax(14),ISmax(14),Smin(14),ISmin(14)
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
        write(ILOG,1009)'//SMS(+)',ESmax(13),IESmax(13),ESmin(13),IESmin(13)
        write(ILOG,1009)'//SMS(-)',ESmax(14),IESmax(14),ESmin(14),IESmin(14)
!C*** Show Summary
        call hecmw_allREDUCE_R(hecMESH,Umax, 6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin, 6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,14,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,14,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,14,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,14,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,EEmax,14,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,EEmin,14,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,ESmax,14,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,ESmin,14,hecmw_min)
        if( hecMESH%my_rank==0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1    ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2    ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3    ',Umax(3),Umin(3)
          write(ILOG,1019) '//R1    ',Umax(4),Umin(4)
          write(ILOG,1019) '//R2    ',Umax(5),Umin(5)
          write(ILOG,1019) '//R3    ',Umax(6),Umin(6)
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
          write(ILOG,1019) '//SMS(+)',Smax(13),Smin(13)
          write(ILOG,1019) '//SMS(-)',Smax(14),Smin(14)
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
          write(ILOG,1019) '//SMS(+)',ESmax(13),ESmin(13)
          write(ILOG,1019) '//SMS(-)',ESmax(14),ESmin(14)
        endif
      endif

 1009 format(a8,1pe12.4,i10,1pe12.4,i10)
 1019 format(a8,1pe12.4,1pe12.4)
  end subroutine fstr_static_post

end module m_static_output
