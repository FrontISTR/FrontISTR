!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Sato(Advancesoft), X. YUAN(AdavanceSoft)    !
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
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \date       2009/09/10
!>  \version    0.00
!!
!======================================================================!
module m_fstr_Result
  use m_fstr_NodalStress
  implicit none

  contains

!> Output result
!----------------------------------------------------------------------*
  subroutine fstr_OutputResult( cstep, totstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, ttime )
!----------------------------------------------------------------------*
    use m_fstr
    use m_static_make_result
    use m_static_output
    use fstr_setup_util
    use m_hecmw2fstr_mesh_conv
    integer, intent(in)                   :: cstep       !< current step number
    integer, intent(in)                   :: totstep     !< step counter
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    type (fstr_solid), intent(inout)      :: fstrSOLID
    type (fstr_param),intent(in)          :: fstrPARAM
    real(kind=kreal),intent(in)           :: ttime

    real(kind=kreal)   :: fdum, ndforce(6), ss(6)
    integer(kind=kint) :: i, k, ii, jj, nctrl, ninfo, ndof, ncomp, fnum
    integer(kind=kint) :: itype, iS, iE, ic_type, icel, maxstep, interval
    real(kind=kreal), allocatable :: ndstrain(:), ndstress(:)
    real(kind=kreal)   :: s11,s22,s33,s12,s23,s13,ps,smises, vdum(6)
    integer(kind=kint),pointer :: member(:)
    type ( hecmwST_result_data ) :: fstrRESULT

    allocate( ndstrain(hecMESH%n_node*6), ndstress(hecMESH%n_node*7) )
    call fstr_NodalStress3D(hecMESH, hecMAT, fstrSOLID,     &
                              ndstrain, ndstress  )
    fstrSOLID%STRAIN = ndstrain
    fstrSOLID%STRESS = ndstress

    if( fstrSOLID%TEMP_irres>0 ) then
        maxstep = fstrSOLID%TEMP_irres
    else
        maxstep = 0
        do i = 1, cstep
            maxstep = maxstep + fstrSOLID%step_ctrl(i)%num_substep
        end do
    endif

    !-- POST PROCESSING VIA MEMORY
    if( IVISUAL==1 .and. &
        (mod(totstep,fstrSOLID%output_ctrl(4)%freqency)==0 .or. totstep==maxstep) ) then
          interval = fstrSOLID%output_ctrl(4)%freqency
          call fstr_make_result(hecMESH,fstrSOLID,fstrRESULT)
          call fstr2hecmw_mesh_conv(hecMESH)
          call hecmw_visualize_init
          call hecmw_visualize(hecMESH,fstrRESULT,totstep,maxstep,interval)
          call hecmw_visualize_finalize
          call hecmw2fstr_mesh_conv(hecMESH)
          call hecmw_result_free(fstrRESULT)
    endif
	
    if( .not. associated(fstrSOLID%output_ctrl) ) return

    ndof = hecMESH%n_dof
    nctrl=1
    if( fstr_output_active( totstep, fstrSOLID%output_ctrl(nctrl) ) ) then
       fnum = fstrSOLID%output_ctrl(1)%filenum
	   ! Summary & Result
       call fstrNLGEOM_Post(fnum,hecMESH,fstrSOLID,ttime,maxstep,totstep,ndstrain,ndstress)
    endif

	   
	if( size(fstrSOLID%output_ctrl)==2 ) then
       iS = fstrSOLID%output_ctrl(2)%outinfo%grp_id
       if( iS<=0 ) return
       fnum = fstrSOLID%output_ctrl(2)%filenum
	   
       do iE= hecMESH%node_group%grp_index(iS-1) + 1, hecMESH%node_group%grp_index(iS  )
          i   = hecMESH%node_group%grp_item(iE)
          jj=fstrPARAM%global_local_id(1,i)
          ii=fstrPARAM%global_local_id(2,i)

          write(fnum,'(2i10,1p6e15.7)') totstep,jj,(fstrSOLID%unode(ndof*ii-ndof+k),k=1,ndof)
       enddo
    endif

    deallocate( ndstrain, ndstress)
  end subroutine fstr_OutputResult

!> Summarizer of output data which prints out max and min output values
!----------------------------------------------------------------------*
  subroutine fstrNLGEOM_Post(fnum,hecMESH,fstrSOLID,tt,maxstep,i_step,ndstrain, ndstress)
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(in)                   :: fnum
      type (hecmwST_local_mesh), intent(in) :: hecMESH
      type (fstr_solid), intent(in)         :: fstrSOLID

      integer(kind=kint), intent(in) :: maxstep
      integer(kind=kint), intent(in) :: i_step
      real(kind=kreal), intent(in)   :: tt
      real(kind=kreal), intent(in)   :: ndstrain(:)
      real(kind=kreal), intent(in)   :: ndstress(:)
! --- file name
      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN) :: label
      character(len=HECMW_NAME_LEN) :: nameID
! --- Local variables
! --- Max Min
      real(kind=kreal) Umax(3),Umin(3),Emax(6),Emin(6),Smax(7),Smin(7)
      integer(kind=kint) IUmax(3),IUmin(3),IEmax(6),IEmin(6),ISmax(7),ISmin(7)
!
      integer(kind=kint) j,i,k,id,ndof,mdof

!
! --- Time
      write(fnum,'(''#### Result time='',E11.4)') tt
! --- Clear
      Umax=0.0d0
      Umin=0.0d0
      Emax=0.0d0
      Emin=0.0d0
      Smax=0.0d0
      Smin=0.0d0
      IUmax=0
      IUmin=0
      IEmax=0
      IEmin=0
      ISmax=0
      ISmin=0
!
! --- define the node of free degree
      ndof = hecMESH%n_dof
!
! --- Displacement
!      write(fnum,*) '#### DISPLACEMENT@POST'
      do i= 1, hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
! --- Max & Min
        if( i.eq.1) then
          do k=1,ndof
            Umax(k)=fstrSOLID%unode(ndof*i-ndof+k)
            Umin(k)=fstrSOLID%unode(ndof*i-ndof+k)
            IUmax(k)=j
            IUmin(k)=j
          enddo
        else
          do k=1,ndof
            if( fstrSOLID%unode(ndof*i-ndof+k) .gt. Umax(k) ) then
              Umax(k)=fstrSOLID%unode(ndof*i-ndof+k)
              IUmax(k)=j
            endif
            if( fstrSOLID%unode(ndof*i-ndof+k) .lt. Umin(k) ) then
              Umin(k)=fstrSOLID%unode(ndof*i-ndof+k)
              IUmin(k)=j
            endif
          enddo
        endif
      enddo
!
!
! ---Strain
      if( ndof == 2 ) ndof =  3
      if( ndof == 3 ) Mdof =  6
      if( ndof == 6 ) Mdof = 10
!      write(fnum,*) '#### STRAIN@POST'
      do i= 1, hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
! --- Max & Min
        if( i==1 ) then
          do k=1,Mdof
            Emax(k)=ndSTRAIN((i-1)*Mdof+k)
            Emin(k)=ndSTRAIN((i-1)*Mdof+k)
            IEmax(k)=j
            IEmin(k)=j
          enddo
        else
          do k=1,Mdof
            if( ndSTRAIN((i-1)*Mdof+k) > Emax(k) ) then
              Emax(k)=ndSTRAIN((i-1)*Mdof+k)
              IEmax(k)=j
            endif
            if( ndSTRAIN((i-1)*Mdof+k) < Emin(k) ) then
              Emin(k)=ndSTRAIN((i-1)*Mdof+k)
              IEmin(k)=j
            endif
          enddo
        endif
      enddo
! --- Stress
!      write(fnum,*) '#### STRESS@POST'
      do i= 1, hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
! --- Max & Min
        if( i.eq.1) then
          do k=1,7
            Smax(k)=ndSTRESS((i-1)*7+k)
            Smin(k)=ndSTRESS((i-1)*7+k)
            ISmax(k)=j
            ISmin(k)=j
          enddo
        else
          do k=1,7
            if( ndSTRESS((i-1)*7+k) > Smax(k) ) then
              Smax(k)=ndSTRESS((i-1)*7+k)
              ISmax(k)=j
            endif
            if( ndSTRESS((i-1)*7+k) < Smin(k) ) then
              Smin(k)=ndSTRESS((i-1)*7+k)
              ISmin(k)=j
            endif
          enddo
        endif
      enddo
!
      if( NDOF==2 .or. NDOF==3 ) then
        write(fnum,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(fnum,*) '//U1 ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(fnum,*) '//U2 ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(fnum,*) '//U3 ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(fnum,*) '//E11',Emax(1),IEmax(1),Emin(1),IEmin(1)
        write(fnum,*) '//E22',Emax(2),IEmax(2),Emin(2),IEmin(2)
        write(fnum,*) '//E33',Emax(3),IEmax(3),Emin(3),IEmin(3)
        write(fnum,*) '//E12',Emax(4),IEmax(4),Emin(4),IEmin(4)
        write(fnum,*) '//E23',Emax(5),IEmax(5),Emin(5),IEmin(5)
        write(fnum,*) '//E13',Emax(6),IEmax(6),Emin(6),IEmin(6)
        write(fnum,*) '//S11',Smax(1),ISmax(1),Smin(1),ISmin(1)
        write(fnum,*) '//S22',Smax(2),ISmax(2),Smin(2),ISmin(2)
        write(fnum,*) '//S33',Smax(3),ISmax(3),Smin(3),ISmin(3)
        write(fnum,*) '//S12',Smax(4),ISmax(4),Smin(4),ISmin(4)
        write(fnum,*) '//S23',Smax(5),ISmax(5),Smin(5),ISmin(5)
        write(fnum,*) '//S13',Smax(6),ISmax(6),Smin(6),ISmin(6)
        write(fnum,*) '//SMS',Smax(7),ISmax(7),Smin(7),ISmin(7)
!
! --- Show Summary
        call hecmw_allREDUCE_R(hecMESH,Umax,3,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Umin,3,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Emax,6,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Emin,6,hecmw_min)
        call hecmw_allREDUCE_R(hecMESH,Smax,7,hecmw_max)
        call hecmw_allREDUCE_R(hecMESH,Smin,7,hecmw_min)
        if( hecMESH%my_rank .eq. 0 ) then
          write(fnum,*) '##### Global Summary :Max/Min####'
          write(fnum,*) '//U1 ',Umax(1),Umin(1)
          write(fnum,*) '//U2 ',Umax(2),Umin(2)
          write(fnum,*) '//U3 ',Umax(3),Umin(3)
          write(fnum,*) '//E11',Emax(1),Emin(1)
          write(fnum,*) '//E22',Emax(2),Emin(2)
          write(fnum,*) '//E33',Emax(3),Emin(3)
          write(fnum,*) '//E12',Emax(4),Emin(4)
          write(fnum,*) '//E23',Emax(5),Emin(5)
          write(fnum,*) '//E13',Emax(6),Emin(6)
          write(fnum,*) '//S11',Smax(1),Smin(1)
          write(fnum,*) '//S22',Smax(2),Smin(2)
          write(fnum,*) '//S33',Smax(3),Smin(3)
          write(fnum,*) '//S12',Smax(4),Smin(4)
          write(fnum,*) '//S23',Smax(5),Smin(5)
          write(fnum,*) '//S13',Smax(6),Smin(6)
          write(fnum,*) '//SMS',Smax(7),Smin(7)
        endif

        if( IRESULT.eq.0 .or. &
            (mod(i_step,fstrSOLID%output_ctrl(3)%freqency).ne.0 .and. i_step.ne.maxstep) ) return
! --- Write Result File
! --- INITIALIZE
        header='*fstrresult'
        call hecmw_result_init(hecMESH,maxstep,i_step,header)
! --- ADD
! --- DISPLACEMENT
!        if( fstrSOLID%iout_list(1) .eq. 1 ) then
          id = 1
          ndof=3
          label='DISPLACEMENT'
          call hecmw_result_add(id,ndof,label,fstrSOLID%unode)
!        end if
!
        if( i_step .gt. 0) then
!        if( fstrSOLID%iout_list(5) .eq. 1 ) then
! --- STRAIN @node
          id = 1
          ndof=6
          label='STRAIN'
          call hecmw_result_add(id,ndof,label,fstrSOLID%STRAIN)
! --- STRAIN @element
!          id = 2
!          ndof=6
!          label='ESTRAIN'
!          call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRAIN)
!        end if
!
!        if( fstrSOLID%iout_list(6) .eq. 1 ) then
! --- STRESS @node
          id = 1
          ndof=7
          label='STRESS'
          call hecmw_result_add(id,ndof,label,fstrSOLID%STRESS)
! --- STRESS @element
!          id = 2
!          ndof=7
!          label='ESTRESS'
!          call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRESS)
!        end if
        end if
!
        nameID='fstrRES'
        call hecmw_result_write_by_name(nameID)
! --- FINALIZE
        call hecmw_result_finalize
!
      endif
      return
!
  end subroutine fstrNLGEOM_Post
!
!
!
end module m_fstr_Result
