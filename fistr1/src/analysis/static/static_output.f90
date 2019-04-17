!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to output result.
module m_static_output
  implicit none
contains

  !> Output result
  !----------------------------------------------------------------------*
  subroutine fstr_static_Output( cstep, istep, time, hecMESH, fstrSOLID, fstrPARAM, flag, outflag )
    !----------------------------------------------------------------------*
    use m_fstr
    use m_fstr_NodalStress
    use m_static_make_result
    use m_hecmw2fstr_mesh_conv
    integer, intent(in)                   :: cstep       !< current step number
    integer, intent(in)                   :: istep       !< current substep number
    real(kind=kreal), intent(in)          :: time        !< current time
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (fstr_solid), intent(inout)      :: fstrSOLID
    type (fstr_param       )              :: fstrPARAM    !< analysis control parameters
    integer, intent(in)                   :: flag
    logical, intent(in)                   :: outflag     !< if true, result will be output regardless of istep

    type ( hecmwST_result_data ) :: fstrRESULT
    integer(kind=kint) :: i, j, ndof, fnum, is, iE, gid

    ndof = hecMESH%n_dof

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

    if( flag==kstSTATICEIGEN ) then
      if( IRESULT==1 .and. &
          (mod(istep,fstrSOLID%output_ctrl(3)%freqency)==0 .or. outflag) ) then
        call fstr_write_static_result( hecMESH, fstrSOLID, fstrPARAM, istep, time, 1 )
      endif
      return
    endif

    if( IRESULT==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(3)%freqency)==0 .or. outflag) ) then
      call fstr_write_static_result( hecMESH, fstrSOLID, fstrPARAM, istep, time, 0 )
    endif

    if( IVISUAL==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(4)%freqency)==0 .or. outflag) ) then

      call fstr_make_static_result( hecMESH, fstrSOLID, fstrRESULT )
      call fstr2hecmw_mesh_conv( hecMESH )
      call hecmw_visualize_init
      call hecmw_visualize( hecMESH, fstrRESULT, istep )
      call hecmw_visualize_finalize
      call hecmw2fstr_mesh_conv( hecMESH )
      call hecmw_result_free( fstrRESULT )
    endif

    if( (mod(istep,fstrSOLID%output_ctrl(1)%freqency)==0 .or. outflag) ) then
      fnum = fstrSOLID%output_ctrl(1)%filenum
      call fstr_static_post( fnum, hecMESH, fstrSOLID, istep )
    endif

    if( fstrSOLID%output_ctrl(2)%outinfo%grp_id>0 .and. &
        (mod(istep,fstrSOLID%output_ctrl(2)%freqency)==0 .or. outflag) ) then
      is = fstrSOLID%output_ctrl(2)%outinfo%grp_id
      fnum = fstrSOLID%output_ctrl(2)%filenum
      do i = hecMESH%node_group%grp_index(is-1)+1, hecMESH%node_group%grp_index(is)
        iE = hecMESH%node_group%grp_item(i)
        gid = hecMESH%global_node_ID(iE)
        write(fnum,'(2i10,1p6e15.7)') istep,gid,(fstrSOLID%unode(ndof*(iE-1)+j),j=1,ndof)
      enddo
    endif

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
    real(kind=kreal)   :: Mmax(1), Mmin(1), EMmax(1), EMmin(1)
    real(kind=kreal)   :: EEmax(14), EEmin(14), ESmax(14), ESmin(14)

    real(kind=kreal)   :: GUmax(6), GUmin(6), GEmax(14), GEmin(14), GSmax(14), GSmin(14)
    real(kind=kreal)   :: GMmax(1), GMmin(1), GEMmax(1), GEMmin(1)
    real(kind=kreal)   :: GEEmax(14), GEEmin(14), GESmax(14), GESmin(14)

    integer(kind=kint) :: IUmax(6), IUmin(6), IEmax(14), IEmin(14), ISmax(14), ISmin(14)
    integer(kind=kint) :: IMmax(1), IMmin(1), IEMmax(1), IEMmin(1)
    integer(kind=kint) :: IEEmax(14), IEEmin(14), IESmax(14), IESmin(14)
    integer(kind=kint) :: i, j, k, ndof, mdof, ID_area
    integer(kind=kint) :: label(6)

    ndof = hecMESH%n_dof

    write( fnum, '(''#### Result step='',I6)') istep
    select case (ndof)
      case (2)
        mdof = 3
        label(1)=11;label(2)=22
        label(3)=12
      case (3,6)
        mdof = 6
        label(1)=11;label(2)=22;label(3)=33;
        label(4)=12;label(5)=23;label(6)=31;
    end select

    j = hecMESH%global_node_ID(1)
    do k = 1, ndof
      Umax(k) = fstrSOLID%unode(k); Umin(k) = fstrSOLID%unode(k)
      IUmax(k)= j; IUmin(k)= j
    enddo
    do k = 1, mdof
      Emax(k) = fstrSOLID%STRAIN(k);   Emin(k) = fstrSOLID%STRAIN(k)
      Smax(k) = fstrSOLID%STRESS(k);   Smin(k) = fstrSOLID%STRESS(k)
      EEmax(k) = fstrSOLID%ESTRAIN(k); EEmin(k) = fstrSOLID%ESTRAIN(k)
      ESmax(k) = fstrSOLID%ESTRESS(k); ESmin(k) = fstrSOLID%ESTRESS(k)
      IEmax(k)= j;  IEmin(k)= j;  ISmax(k)= j;  ISmin(k)= j
      IEEmax(k)= j; IEEmin(k)= j; IESmax(k)= j; IESmin(k)= j
    enddo
    Mmax(1) = fstrSOLID%MISES(1); Mmin(1) = fstrSOLID%MISES(1)
    EMmax(1) = fstrSOLID%EMISES(1); EMmin(1) = fstrSOLID%EMISES(1)
    IMmax(1)= j; IMmin(1)= j; IEMmax(1)= j; IEMmin(1)= j

    !C*** Show Displacement
    do i = 1, hecMESH%nn_internal
      if(fstrSOLID%is_rot(i)==1)cycle
      j = hecMESH%global_node_ID(i)
      do k = 1, ndof
        if     ( fstrSOLID%unode(ndof*(i-1)+k) > Umax(k) ) then
          Umax(k) = fstrSOLID%unode(ndof*(i-1)+k)
          IUmax(k)= j
        else if( fstrSOLID%unode(ndof*(i-1)+k) < Umin(k) ) then
          Umin(k) = fstrSOLID%unode(ndof*(i-1)+k)
          IUmin(k)= j
        endif
      enddo
    enddo
    !C*** Nodal Strain / Stress / MISES
    !C @node
    do i = 1, hecMESH%nn_internal
      if(fstrSOLID%is_rot(i)==1)cycle
      j = hecMESH%global_node_ID(i)
      do k = 1, mdof
        if     ( fstrSOLID%STRAIN(mdof*(i-1)+k) > Emax(k) ) then
          Emax(k) = fstrSOLID%STRAIN(mdof*(i-1)+k)
          IEmax(k)= j
        else if( fstrSOLID%STRAIN(mdof*(i-1)+k) < Emin(k) ) then
          Emin(k) = fstrSOLID%STRAIN(mdof*(i-1)+k)
          IEmin(k)= j
        endif
        if     ( fstrSOLID%STRESS(mdof*(i-1)+k) > Smax(k) ) then
          Smax(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
          ISmax(k)= j
        else if( fstrSOLID%STRESS(mdof*(i-1)+k) < Smin(k) ) then
          Smin(k) = fstrSOLID%STRESS(mdof*(i-1)+k)
          ISmin(k)= j
        endif
      enddo
      if     ( fstrSOLID%MISES(i) > Mmax(1) ) then
        Mmax(1) = fstrSOLID%MISES(i)
        IMmax(1)= j
      else if( fstrSOLID%MISES(i) < Mmin(1) ) then
        Mmin(1) = fstrSOLID%MISES(i)
        IMmin(1)= j
      endif
    enddo
    !C*** Elemental Strain / STRESS
    do i = 1, hecMESH%n_elem
      ID_area = hecMESH%elem_ID(i*2)
      if( ID_area==hecMESH%my_rank ) then
        j = hecMESH%global_elem_ID(i)
        do k = 1, mdof
          if( fstrSOLID%ESTRAIN(mdof*(i-1)+k) > EEmax(k) ) then
            EEmax(k) = fstrSOLID%ESTRAIN(mdof*(i-1)+k)
            IEEmax(k)= j
          else if( fstrSOLID%ESTRAIN(mdof*(i-1)+k) < EEmin(k) ) then
            EEmin(k) = fstrSOLID%ESTRAIN(mdof*(i-1)+k)
            IEEmin(k)= j
          endif
          if( fstrSOLID%ESTRESS(mdof*(i-1)+k) > ESmax(k) ) then
            ESmax(k) = fstrSOLID%ESTRESS(mdof*(i-1)+k)
            IESmax(k)= j
          else if( fstrSOLID%ESTRESS(mdof*(i-1)+k) < ESmin(k) ) then
            ESmin(k) = fstrSOLID%ESTRESS(mdof*(i-1)+k)
            IESmin(k)= j
          endif
        enddo
        if( fstrSOLID%EMISES(i) > EMmax(1) ) then
          EMmax(1) = fstrSOLID%EMISES(i)
          IEMmax(1)= j
        else if( fstrSOLID%EMISES(i) < EMmin(1) ) then
          EMmin(1) = fstrSOLID%EMISES(i)
          IEMmin(1)= j
        endif
      endif
    enddo


    write(ILOG,*)    '##### Local Summary @Node    :Max/IdMax/Min/IdMin####'
    do i = 1, ndof; write(ILOG,1029) ' //U',i,      '  ',Umax(i),IUmax(i),Umin(i),IUmin(i);     end do
    do i = 1, mdof; write(ILOG,1029) ' //E',label(i),' ',Emax(i),IEmax(i),Emin(i),IEmin(i);     end do
    do i = 1, mdof; write(ILOG,1029) ' //S',label(i),' ',Smax(i),ISmax(i),Smin(i),ISmin(i);     end do
    write(ILOG,1009) '//SMS '           ,Mmax(1),IMmax(1),Mmin(1),IMmin(1)
    write(ILOG,*)    '##### Local Summary @Element :Max/IdMax/Min/IdMin####'
    do i = 1, mdof; write(ILOG,1029) ' //E',label(i),' ',EEmax(i),IEEmax(i),EEmin(i),IEEmin(i); end do
    do i = 1, mdof; write(ILOG,1029) ' //S',label(i),' ',ESmax(i),IESmax(i),ESmin(i),IESmin(i); end do
    write(ILOG,1009) '//SMS '           ,EMmax(1),IEMmax(1),EMmin(1),IEMmin(1)

    !C*** Show Summary
    GUmax  = Umax; GUmin  = Umin;
    GEmax  = Emax; GEmin  = Emin; GEEmax = EEmax; GEEmin = EEmin;
    GSmax  = Smax; GSmin  = Smin; GESmax = ESmax; GESmin = ESmin;
    GMmax  = Mmax; GMmin  = Mmin; GEMmax = EMmax; GEMmin = EMmin;

    call hecmw_allREDUCE_R(hecMESH,GUmax,ndof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GUmin,ndof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GEmax,mdof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GEmin,mdof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GSmax,mdof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GSmin,mdof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GMmax,1,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GMmin,1,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GEEmax,mdof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GEEmin,mdof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GESmax,mdof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GESmin,mdof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GEMmax,1,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GEMmin,1,hecmw_min)

    do i=1,ndof
      if(GUmax (i) > Umax (i)) IUmax (i) = 0
      if(GUmin (i) < Umin (i)) IUmin (i) = 0
    enddo
    do i=1,mdof
      if(GEmax (i) > Emax (i)) IEmax (i) = 0
      if(GSmax (i) > Smax (i)) ISmax (i) = 0
      if(GEEmax(i) > EEmax(i)) IEEmax(i) = 0
      if(GESmax(i) > ESmax(i)) IESmax(i) = 0
      if(GEmin (i) < Emin (i)) IEmin (i) = 0
      if(GSmin (i) < Smin (i)) ISmin (i) = 0
      if(GEEmin(i) < EEmin(i)) IEEmin(i) = 0
      if(GESmin(i) < ESmin(i)) IESmin(i) = 0
    enddo
    do i=1,1
      if(GMmax (i) > Mmax (i)) IMmax (i) = 0
      if(GEMmax(i) > EMmax(i)) IEMmax(i) = 0
      if(GMmin (i) < Mmin (i)) IMmin (i) = 0
      if(GEMmin(i) < EMmin(i)) IEMmin(i) = 0
    enddo

    call hecmw_allREDUCE_I(hecMESH,IUmax,ndof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IUmin,ndof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IEmax,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IEmin,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,ISmax,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,ISmin,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IMmax,1,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IMmin,1,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IEEmax,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IEEmin,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IESmax,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IESmin,mdof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IEMmax,1,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IEMmin,1,hecmw_max)

    if( hecMESH%my_rank==0 ) then
      write(ILOG,*)    '##### Global Summary @Node    :Max/IdMax/Min/IdMin####'
      do i = 1, ndof; write(ILOG,1029) ' //U',i,      '  ',GUmax(i),IUmax(i),GUmin(i),IUmin(i);     end do
      do i = 1, mdof; write(ILOG,1029) ' //E',label(i),' ',GEmax(i),IEmax(i),GEmin(i),IEmin(i);     end do
      do i = 1, mdof; write(ILOG,1029) ' //S',label(i),' ',GSmax(i),ISmax(i),GSmin(i),ISmin(i);     end do
      write(ILOG,1009) '//SMS '           ,GMmax(1),IMmax(1),GMmin(1),IMmin(1)
      write(ILOG,*)    '##### Global Summary @Element :Max/IdMax/Min/IdMin####'
      do i = 1, mdof; write(ILOG,1029) ' //E',label(i),' ',GEEmax(i),IEEmax(i),GEEmin(i),IEEmin(i); end do
      do i = 1, mdof; write(ILOG,1029) ' //S',label(i),' ',GESmax(i),IESmax(i),GESmin(i),IESmin(i); end do
      write(ILOG,1009) '//SMS '           ,GEMmax(1),IEMmax(1),GEMmin(1),IEMmin(1)
    endif

    1009 format(a7,1pe12.4,i10,1pe12.4,i10)
    1029 format(a,i0,a,1pe12.4,i10,1pe12.4,i10)
  end subroutine fstr_static_post

end module m_static_output
