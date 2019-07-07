!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to output result.
module m_dynamic_output
  implicit none
contains

  !> Output result
  !----------------------------------------------------------------------*
  subroutine fstr_dynamic_Output( hecMESH, fstrSOLID, fstrDYNAMIC, fstrPARAM )
    !----------------------------------------------------------------------*
    use m_fstr
    use m_fstr_Update
    use m_fstr_NodalStress
    use m_dynamic_make_result
    use m_hecmw2fstr_mesh_conv
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(inout)      :: fstrSOLID
    type(fstr_dynamic), intent(in)       :: fstrDYNAMIC
    type(fstr_param), intent(in)         :: fstrPARAM

    type(hecmwST_result_data) :: fstrRESULT
    integer(kind=kint) :: i, j, ndof, maxstep, interval, fnum, is, iE, gid, istep, idx

    ndof = hecMESH%n_dof

    !C-- SET DISPLACEMENT etc.
    istep = fstrDYNAMIC%i_step
    if( fstrDYNAMIC%idx_eqa==1 .and. istep>0 ) then
      idx = 2
    else
      idx = 1
    endif

    if( fstrSOLID%TEMP_ngrp_tot>0 .or. fstrSOLID%TEMP_irres>0 ) then
      if( ndof==3 ) then
        allocate( fstrSOLID%tnstrain(hecMESH%n_node*6) )
        allocate( fstrSOLID%testrain(hecMESH%n_elem*6) )
      else if( ndof==2 ) then
        allocate( fstrSOLID%tnstrain(hecMESH%n_node*3) )
        allocate( fstrSOLID%testrain(hecMESH%n_elem*3) )
      else if( ndof== 4) then
        write(*,*)'Error: This routine is not implemented'
        stop
      endif
    endif

    if( ndof==3 .or. ndof == 4 ) then
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
      is = fstrSOLID%output_ctrl(2)%outinfo%grp_id
      fnum = fstrSOLID%output_ctrl(2)%filenum
      do i = hecMESH%node_group%grp_index(is-1)+1, hecMESH%node_group%grp_index(is)
        iE = hecMESH%node_group%grp_item(i)
        gid = hecMESH%global_node_ID(iE)
        write(fnum,'(2i10,1p6e15.7)') istep,gid,(fstrSOLID%unode(ndof*(iE-1)+j),j=1,ndof)
      enddo
    endif

    if( IRESULT==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(3)%freqency)==0 .or. istep==maxstep) ) then
      call fstr_write_dynamic_result( hecMESH, fstrSOLID, fstrDYNAMIC, maxstep, istep, fstrDYNAMIC%t_curr )
    endif

    if( IVISUAL==1 .and. &
        (mod(istep,fstrSOLID%output_ctrl(4)%freqency)==0 .or. istep==maxstep) ) then
      call fstr_make_dynamic_result( hecMESH, fstrSOLID, fstrDYNAMIC, fstrRESULT )
      call fstr2hecmw_mesh_conv( hecMESH )
      call hecmw_visualize_init
      call hecmw_visualize( hecMESH, fstrRESULT, istep )
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
    integer, intent(in)                  :: fnum
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in)         :: fstrSOLID
    type(fstr_dynamic), intent(in)       :: fstrDYNAMIC

    real(kind=kreal)   :: Umax(6), Umin(6), Vmax(6), Vmin(6), Amax(6), Amin(6)
    real(kind=kreal)   :: Emax(14), Emin(14), Smax(14), Smin(14)
    real(kind=kreal)   :: Mmax(1), Mmin(1), EMmax(1), EMmin(1)
    real(kind=kreal)   :: EEmax(14), EEmin(14), ESmax(14), ESmin(14)

    real(kind=kreal)   :: GUmax(6), GUmin(6), GVmax(6), GVmin(6), GAmax(6), GAmin(6)
    real(kind=kreal)   :: GEmax(14), GEmin(14), GSmax(14), GSmin(14)
    real(kind=kreal)   :: GMmax(1), GMmin(1), GEMmax(1), GEMmin(1)
    real(kind=kreal)   :: GEEmax(14), GEEmin(14), GESmax(14), GESmin(14)

    integer(kind=kint) :: IUmax(6), IUmin(6), IVmax(6), IVmin(6), IAmax(6), IAmin(6)
    integer(kind=kint) :: IEmax(14), IEmin(14), ISmax(14), ISmin(14)
    integer(kind=kint) :: IMmax(1), IMmin(1), IEMmax(1), IEMmin(1)
    integer(kind=kint) :: IEEmax(14), IEEmin(14), IESmax(14), IESmin(14)
    integer(kind=kint) :: i, j, k, ndof, mdof, ID_area, idx
    integer(kind=kint) :: label(6)

    if( fstrDYNAMIC%i_step==0 ) return
    if( fstrDYNAMIC%idx_eqa==1 .and. fstrDYNAMIC%i_step>0 ) then
      idx = 2
    else
      idx = 1
    endif
    ndof = hecMESH%n_dof

    write( fnum, '(''#### Result step='',I6)') fstrDYNAMIC%i_step
    select case (ndof)
      case (2)
        mdof = 3
        label(1)=11;label(2)=22
        label(3)=12
      case (3,4,6)
        mdof = 6
        label(1)=11;label(2)=22;label(3)=33;
        label(4)=12;label(5)=23;label(6)=31;
    end select

    j = hecMESH%global_node_ID(1)
    do k = 1, ndof
      Umax(k) = fstrDYNAMIC%DISP(k,idx); Umin(k) = fstrDYNAMIC%DISP(k,idx)
      Vmax(k) = fstrDYNAMIC%VEL(k,idx); Vmin(k) = fstrDYNAMIC%VEL(k,idx)
      Amax(k) = fstrDYNAMIC%ACC(k,idx); Amin(k) = fstrDYNAMIC%ACC(k,idx)
      IUmax(k)= j; IUmin(k)= j; IVmax(k)= j; IVmin(k)= j; IAmax(k)= j; IAmin(k)= j
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

    !C*** Show Displacement / Velosity / Acc
    do i = 1, hecMESH%nn_internal
      if(fstrSOLID%is_rot(i)==1)cycle
      j = hecMESH%global_node_ID(i)
      do k = 1, ndof
        if     ( fstrDYNAMIC%DISP(ndof*(i-1)+k,idx) > Umax(k) ) then
          Umax(k) = fstrDYNAMIC%DISP(ndof*(i-1)+k,idx)
          IUmax(k)= j
        else if( fstrDYNAMIC%DISP(ndof*(i-1)+k,idx) < Umin(k) ) then
          Umin(k) = fstrDYNAMIC%DISP(ndof*(i-1)+k,idx)
          IUmin(k)= j
        endif
        if     ( fstrDYNAMIC%VEL(ndof*(i-1)+k,idx) > Vmax(k) ) then
          Vmax(k) = fstrDYNAMIC%VEL(ndof*(i-1)+k,idx)
          IVmax(k)= j
        else if( fstrDYNAMIC%VEL(ndof*(i-1)+k,idx) < Vmin(k) ) then
          Vmin(k) = fstrDYNAMIC%VEL(ndof*(i-1)+k,idx)
          IVmin(k)= j
        endif
        if     ( fstrDYNAMIC%ACC(ndof*(i-1)+k,idx) > Amax(k) ) then
          Amax(k) = fstrDYNAMIC%ACC(ndof*(i-1)+k,idx)
          IAmax(k)= j
        else if( fstrDYNAMIC%ACC(ndof*(i-1)+k,idx) < Amin(k) ) then
          Amin(k) = fstrDYNAMIC%ACC(ndof*(i-1)+k,idx)
          IAmin(k)= j
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
    if (ndof == 4)  write(ILOG,1009) ' //P  '           , Umax(4),IUmax(4),Umin(4),IUmin(4)
    do i = 1, ndof; write(ILOG,1029) ' //V',i,      '  ',Vmax(i),IVmax(i),Vmin(i),IVmin(i);     end do
    do i = 1, ndof; write(ILOG,1029) ' //A',i,      '  ',Amax(i),IAmax(i),Amin(i),IAmin(i);     end do
    do i = 1, mdof; write(ILOG,1029) ' //E',label(i),' ',Emax(i),IEmax(i),Emin(i),IEmin(i);     end do
    do i = 1, mdof; write(ILOG,1029) ' //S',label(i),' ',Smax(i),ISmax(i),Smin(i),ISmin(i);     end do
    write(ILOG,1009) '//SMS '           ,Mmax(1),IMmax(1),Mmin(1),IMmin(1)
    write(ILOG,*)    '##### Local Summary @Element :Max/IdMax/Min/IdMin####'
    do i = 1, mdof; write(ILOG,1029) ' //E',label(i),' ',EEmax(i),IEEmax(i),EEmin(i),IEEmin(i); end do
    do i = 1, mdof; write(ILOG,1029) ' //S',label(i),' ',ESmax(i),IESmax(i),ESmin(i),IESmin(i); end do
    write(ILOG,1009) '//SMS '           ,EMmax(1),IEMmax(1),EMmin(1),IEMmin(1)

    !C*** Show Summary
    GUmax  = Umax; GUmin  = Umin; GVmax  = Vmax; GVmin  = Vmin; GAmax  = Amax; GAmin  = Amin;
    GEmax  = Emax; GEmin  = Emin; GEEmax = EEmax; GEEmin = EEmin;
    GSmax  = Smax; GSmin  = Smin; GESmax = ESmax; GESmin = ESmin;
    GMmax  = Mmax; GMmin  = Mmin; GEMmax = EMmax; GEMmin = EMmin;

    call hecmw_allREDUCE_R(hecMESH,GUmax,ndof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GUmin,ndof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GVmax,ndof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GVmin,ndof,hecmw_min)
    call hecmw_allREDUCE_R(hecMESH,GAmax,ndof,hecmw_max)
    call hecmw_allREDUCE_R(hecMESH,GAmin,ndof,hecmw_min)
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
      if(GVmax (i) > Vmax (i)) IVmax (i) = 0
      if(GAmax (i) > Amax (i)) IAmax (i) = 0
      if(GUmin (i) < Umin (i)) IUmin (i) = 0
      if(GVmin (i) < Vmin (i)) IVmin (i) = 0
      if(GAmin (i) < Amin (i)) IAmin (i) = 0
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
    call hecmw_allREDUCE_I(hecMESH,IVmax,ndof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IVmin,ndof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IAmax,ndof,hecmw_max)
    call hecmw_allREDUCE_I(hecMESH,IAmin,ndof,hecmw_max)
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
      if (ndof == 4)  write(ILOG,1009) ' //P  '           ,GUmax(4),IUmax(4),GUmin(4),IUmin(4)
      do i = 1, ndof; write(ILOG,1029) ' //V',i,      '  ',GVmax(i),IVmax(i),GVmin(i),IVmin(i);     end do
      do i = 1, ndof; write(ILOG,1029) ' //A',i,      '  ',GAmax(i),IAmax(i),GAmin(i),IAmin(i);     end do
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
  end subroutine fstr_dynamic_post

  !C================================================================C
  !C-- subroutine dynamic_output_monit
  !C================================================================C
  subroutine dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, fstrEIG, fstrSOLID)
    use m_fstr
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_param)         :: fstrPARAM
    type(fstr_dynamic)       :: fstrDYNAMIC
    type(fstr_eigen)         :: fstrEIG
    type(fstr_solid)         :: fstrSOLID

    integer(kind=kint) :: idx, ii, jj, ierr, ncmp
    integer(kind=kint) :: num_monit, ig, is, iE, ik, iunitS, iunit
    logical :: yes

    if( mod(fstrDYNAMIC%i_step,fstrDYNAMIC%nout_monit)/=0 ) return

    if( fstrDYNAMIC%idx_eqa==1 .and. fstrDYNAMIC%i_step>0 ) then
      idx = 2
    else
      idx = 1
    endif

    num_monit = 0
    ig = fstrDYNAMIC%ngrp_monit
    is = hecMESH%node_group%grp_index(ig-1)+1
    iE = hecMESH%node_group%grp_index(ig)
    do ik=is,iE
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
      !C-- nodal force
      if( fstrDYNAMIC%iout_list(4)==1 ) then
        iunit = iunitS + fstrDYNAMIC%dynamic_IW7
        write( iunit, '(i10,1pe13.4e3,i10,1p6e13.4e3)' ) &
          fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, jj, &
          fstrSOLID%QFORCE( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii )
      end if
      !C-- strain
      if( fstrDYNAMIC%iout_list(5) > 0 ) then
        if (hecMESH%n_dof == 2 .or. hecMESH%n_dof == 3 .or. hecMESH%n_dof == 4) then
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
        if (hecMESH%n_dof == 2 .or. hecMESH%n_dof == 3 .or. hecMESH%n_dof == 4) then
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
          open( fstrDYNAMIC%dynamic_IW10, file='dyna_energy.txt', status='replace', iostat=ierr )
          if( ierr/=0 ) then
            write(*,*) 'stop due to file opening error <dyna_enrgy.txt>'
            call hecmw_abort( hecmw_comm_get_comm() )
          endif
          write( fstrDYNAMIC%dynamic_IW10, * ) &
            ' time step', '     time    ', '  kinetic energy', '   strain energy', '   total energy'
        endif
        if(fstrDYNAMIC%i_step==0) then
          fstrDYNAMIC%kineticEnergy = 0.0d0
          do ii = 1, hecMESH%n_node*hecMESH%n_dof
            fstrDYNAMIC%kineticEnergy = fstrDYNAMIC%kineticEnergy &
              + 0.5d0 * fstrEIG%mass(ii) * fstrDYNAMIC%VEL(ii,idx) * fstrDYNAMIC%VEL(ii,idx)
          enddo
        endif
        fstrDYNAMIC%totalEnergy = fstrDYNAMIC%kineticEnergy + fstrDYNAMIC%strainEnergy
        write( fstrDYNAMIC%dynamic_IW10, '(i10,1pe13.4e3,1p3e16.4e3)' ) &
          fstrDYNAMIC%i_step, fstrDYNAMIC%t_curr, &
          fstrDYNAMIC%kineticEnergy, fstrDYNAMIC%strainEnergy, fstrDYNAMIC%totalEnergy
      endif
      if( fstrDYNAMIC%i_step==fstrDYNAMIC%n_step ) close(fstrDYNAMIC%dynamic_IW10)
    endif
  end subroutine dynamic_output_monit

  !C================================================================C
  !C-- subroutine  matvec
  !C================================================================C
  subroutine matvec(y,x,hecMAT,ndof,D,AU,AL)
    use m_fstr
    type(hecmwST_matrix) :: hecMAT

    integer(kind=kint) :: ndof, i, is, ie, j, icol, ki, kj, ix, iy, ip, nn
    real(kind=kreal)   :: D(ndof*ndof*hecMAT%NP)
    real(kind=kreal)   :: AU(ndof*ndof*hecMAT%NPU)
    real(kind=kreal)   :: AL(ndof*ndof*hecMAT%NPL)
    real(kind=kreal)   :: x(ndof*hecMAT%NP)
    real(kind=kreal)   :: y(ndof*hecMAT%NP)

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
