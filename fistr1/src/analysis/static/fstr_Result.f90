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
    integer, intent(in)                   :: totstep     !< total steps
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    type (hecmwST_matrix), intent(in)     :: hecMAT
    type (fstr_solid), intent(inout)      :: fstrSOLID
    type (fstr_param),intent(in)          :: fstrPARAM
    real(kind=kreal),intent(in)           :: ttime

    real(kind=kreal)   :: fdum, ndforce(6), ss(6)
    integer(kind=kint) :: i, k, ii, jj, nctrl, ninfo, ndof, ncomp, fnum
    integer(kind=kint) :: maxiter, itype, iS, iE, ic_type, icel
    real(kind=kreal), allocatable :: ndstrain(:), ndstress(:)
    real(kind=kreal)   :: s11,s22,s33,s12,s23,s13,ps,smises, vdum(6)
    integer(kind=kint),pointer :: member(:)
    type ( hecmwST_result_data ) :: fstrRESULT

    allocate( ndstrain(hecMESH%n_node*6), ndstress(hecMESH%n_node*7) )
    call fstr_NodalStress3D(hecMESH, hecMAT, fstrSOLID,     &
                              ndstrain, ndstress  )
    fstrSOLID%STRAIN = ndstrain
    fstrSOLID%STRESS = ndstress

    maxiter = fstrSOLID%step_ctrl(cstep)%max_iter
	!    call hecmw_nullify_result_data( fstrRESULT )
    !-- POST PROCESSING VIA MEMORY
    if( IVISUAL==1 ) then
          call fstr_make_result(hecMESH,fstrSOLID,fstrRESULT)
          call fstr2hecmw_mesh_conv(hecMESH)
          call hecmw_visualize_init

          call hecmw_visualize(hecMESH,fstrRESULT,totstep,maxiter,0)
          call hecmw_visualize_finalize
          call hecmw2fstr_mesh_conv(hecMESH)
          call hecmw_result_free(fstrRESULT)
    endif
	
    if( .not. associated(fstrSOLID%output_ctrl) ) return

    ndof = hecMESH%n_dof
    nctrl=1
    if( fstr_output_active( totstep, fstrSOLID%output_ctrl(nctrl) ) ) then
	   

       ndforce(3:6)=0.d0
       fnum = fstrSOLID%output_ctrl(nctrl)%filenum
       write(fnum,*) "Step:",totstep
       do ninfo=1, fstrSOLID%output_ctrl(nctrl)%outinfo%num_items
          if( .not. fstrSOLID%output_ctrl(nctrl)%outinfo%on(ninfo) ) cycle
                 
          ncomp = n_comp_valtype(fstrSOLID%output_ctrl(nctrl)%outinfo%vtype(ninfo), ndof)
          if( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "DISPLACEMENT" ) then

! ----- DISPLACEMENT
            if( ndof == 3 ) then
              write(fnum,'(a)') '#### DISPLACEMENT 3D'
              write(fnum,'(a)') '    NODE      X-DISP         Y-DISP         Z-DISP     '
              write(fnum,'(a)') '---------+--------------+--------------+---------------'
            else if( ndof == 2 ) then
              write(fnum,'(a)') '#### DISPLACEMENT 2D'
              write(fnum,'(a)') '    NODE      X-DISP         Y-DISP     '
              write(fnum,'(a)') '---------+--------------+---------------'
            else if( ndof == 6 ) then
            endif
            do i= 1, hecMESH%nn_internal
             jj=fstrPARAM%global_local_id(1,i)
             ii=fstrPARAM%global_local_id(2,i)
             write(fnum,'(i10,1p6e15.7)') jj,(fstrSOLID%unode(ndof*ii-ndof+k),k=1,ndof)
            enddo
!
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "REACTIONFORCE" ) then

! ----- REACTION FORCE
            if( ndof == 3 ) then
               write(fnum,'(a)') '#### REACTION FORCE 3D'
               write(fnum,'(a)') '    NODE      X-DIREC       Y-DIREC         Z-DIREC    '
               write(fnum,'(a)') '---------+--------------+--------------+---------------'
            else if( ndof == 2 ) then
               write(fnum,'(a)') '#### REACTION FORCE 2D'
               write(fnum,'(a)') '    NODE     X-DIREC        Y-DIREC     '
               write(fnum,'(a)') '---------+--------------+---------------'
            else if( ndof == 6 ) then
            endif
            IF( fstrSOLID%output_ctrl(nctrl)%outinfo%grp_id_name=="ALL" ) THEN
              do i=1,hecMESH%nn_internal
                jj=fstrPARAM%global_local_id(1,i)
                ii=fstrPARAM%global_local_id(2,i)
                ndforce(1:hecMESH%n_dof)                                                &
                  = fstrSOLID%QFORCE( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*(ii-1)+hecMESH%n_dof)
                write(fnum,'(i10,1p6e15.7)') jj, ndforce(1:hecMESH%n_dof)
              enddo
            ELSE
              ii = get_grp_member_n( hecMESH, 'node_grp', fstrSOLID%output_ctrl(nctrl)%outinfo%grp_id_name )
              allocate( member(ii) )
              ii= get_grp_member( hecMESH, 'node_grp', fstrSOLID%output_ctrl(nctrl)%outinfo%grp_id_name , member )
              if( ii>0 .and. associated(member) ) then
                vdum = 0.d0
                do i=1,size(member)
                  ii =member(i)
                !  jj=fstrPARAM%global_local_id(1,member(i))
                !  ii=fstrPARAM%global_local_id(2,member(i))
                  vdum = vdum +fstrSOLID%QFORCE( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*(ii-1)+hecMESH%n_dof)
                enddo
                write(fnum,'(a,1p6e15.7)') "    SUM:", vdum(1:hecMESH%n_dof)
              endif
              if( associated(member) ) deallocate( member )
            ENDIF
			  

! ------- NODAL STRAIN
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRAIN"   &
             .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) == -1) then
            if( ndof == 3 ) then
              call fstr_NodalStress3D( hecMESH, hecMAT, fstrSOLID, ndstrain, ndstress  )
              write(fnum,*) '#### NODAL STRAIN'
              write(fnum,'(a,a)') '     NODE       E11            E22            E33      ',   &
                                    '      E12            E23            E13      '
              write(fnum,'(a,a)') '---------+--------------+--------------+--------------+',   &
                                    '--------------+--------------+--------------+'
              do i=1,hecMESH%nn_internal
                 jj=fstrPARAM%global_local_id(1,i)
                 ii=fstrPARAM%global_local_id(2,i)
                 write(fnum,'(i10,1p6e15.7)') jj, (ndSTRAIN(6*(ii-1)+k),k=1,6)
              enddo
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif

! ------- NODAL STRESS
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRESS"   &
            .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) == -1) then
            if( ndof == 3 ) then
              call fstr_NodalStress3D(hecMESH, hecMAT, fstrSOLID,     &
                              ndstrain, ndstress  )
              write(fnum,*) '#### NODAL STRESS'
              write(fnum,'(a,a)') '     NODE       S11            S22            S33      ',   &
                                    '      S12            S23            S13      '
              write(fnum,'(a,a)') '---------+--------------+--------------+--------------+',   &
                                    '--------------+--------------+--------------+'
              do i=1,hecMESH%nn_internal
                jj=fstrPARAM%global_local_id(1,i)
                ii=fstrPARAM%global_local_id(2,i)
                write(fnum,'(i10,1p7e15.7)') jj, (ndSTRESS(7*(ii-1)+k),k=1,7)
              enddo
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif
			
! ---------ELMENTAL STRAIN
        else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRAIN"   &
             .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) >= 0) then
            if( ndof == 3 ) then
              write(fnum,*) '#### ELEMENTAL STRAIN'
              do itype= 1, hecMESH%n_elem_type
                iS= hecMESH%elem_type_index(itype-1) + 1
                iE= hecMESH%elem_type_index(itype  )
                ic_type= hecMESH%elem_type_item(itype)

                if (hecmw_is_etype_link(ic_type)) cycle
                if( ic_type==301 ) ic_type=111
                ii = NumOfQuadPoints( ic_type )
                do icel= iS, iE
                    jj = hecMESH%global_elem_ID(icel)
                    write(fnum,*) '--- ELEMENT:', jj
                    select case ( fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) )
                    case (1)  ! average
                      ss(:) = 0.d0
                      do jj=1,ii
                         ss(:) = ss(:) + fstrSOLID%elements(icel)%gausses(jj)%strain(:)
                      enddo
                      ss(:) = ss(:)/ii
                      write(fnum,'(1p6e15.7)') (ss(k),k=1,6)
                    case (2)  ! quadrature point
                      do jj=1,ii
                        write(fnum,'(i10,1p6e15.7)') jj, (fstrSOLID%elements(icel)%gausses(jj)%strain(k),k=1,6)
                      enddo
                    end select
                enddo
              enddo
			  
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif

! ------- ELEMENTAL STRESS
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRESS"   &
            .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) >= 0) then
            if( ndof == 3 ) then
              write(fnum,*) '#### ELEMENTAL STRESS'
              do itype= 1, hecMESH%n_elem_type
                iS= hecMESH%elem_type_index(itype-1) + 1
                iE= hecMESH%elem_type_index(itype  )
                ic_type= hecMESH%elem_type_item(itype)

                if (hecmw_is_etype_link(ic_type)) cycle
                if( ic_type==301 ) ic_type=111
                ii = NumOfQuadPoints( ic_type )
                do icel= iS, iE
                    jj = hecMESH%global_elem_ID(icel)
                    write(fnum,*) '--- ELEMENT:', jj
                    select case ( fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) )
                    case (1)  ! average
                      ss(:) = 0.d0
                      do jj=1,ii
                         ss(:) = ss(:) + fstrSOLID%elements(icel)%gausses(jj)%stress(:)
                      enddo
                      ss(:) = ss(:)/ii
					  s11=ss(1)
                      s22=ss(2)
                      s33=ss(3)
                      s12=ss(4)
                      s23=ss(5)
                      s13=ss(6)
                      ps=(s11+s22+s33)/3.0
                      smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                            +s12**2+s23**2+s13**2 
                      smises=dsqrt(3.d0*smises)
                      write(fnum,'(1p7e15.7)') (ss(k),k=1,6),smises
                    case (2)  ! quadrature point
                      do jj=1,ii
					    s11=fstrSOLID%elements(icel)%gausses(jj)%stress(1)
                        s22=fstrSOLID%elements(icel)%gausses(jj)%stress(2)
                        s33=fstrSOLID%elements(icel)%gausses(jj)%stress(3)
                        s12=fstrSOLID%elements(icel)%gausses(jj)%stress(4)
                        s23=fstrSOLID%elements(icel)%gausses(jj)%stress(5)
                        s13=fstrSOLID%elements(icel)%gausses(jj)%stress(6)
                        ps=(s11+s22+s33)/3.0
                        smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                            +s12**2+s23**2+s13**2 
                        smises=dsqrt(3.d0*smises)
                        write(fnum,'(i10,1p7e15.7)') jj, (fstrSOLID%elements(icel)%gausses(jj)%stress(k),k=1,6),smises
                      enddo
                    end select
                enddo
              enddo
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif
			
! ------- ELEMENTAL PLASTIC STRAIN
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "PLASTIC_STRAIN" ) then
              write(fnum,*) '#### PLASTIC STRAIN'
              do itype= 1, hecMESH%n_elem_type
                iS= hecMESH%elem_type_index(itype-1) + 1
                iE= hecMESH%elem_type_index(itype  )
                ic_type= hecMESH%elem_type_item(itype)

                if (hecmw_is_etype_link(ic_type)) cycle

                ii = NumOfQuadPoints( ic_type )
                do icel= iS, iE
                    jj = hecMESH%global_elem_ID(icel)
                    write(fnum,*) '--- ELEMENT:', jj
                    select case ( fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) )
                    case (1)  ! average
                      fdum = 0.d0
                      do jj=1,ii
                        if( isElastoplastic( fstrSOLID%elements(icel)%gausses(jj)%pMaterial%mtype )   &
                          .and. associated(fstrSOLID%elements(icel)%gausses(jj)%fstatus) )           &
                         fdum = fdum + fstrSOLID%elements(icel)%gausses(jj)%fstatus(1)
                      enddo
                      write(fnum,'(1p6e15.7)') fdum/ii
                    case (2)  ! quadrature point
                      do jj=1,ii
                        if( isElastoplastic( fstrSOLID%elements(icel)%gausses(jj)%pMaterial%mtype )   &
                          .and. associated(fstrSOLID%elements(icel)%gausses(jj)%fstatus) )             &
                        write(fnum,'(i10,1p1e15.7)') jj, (fstrSOLID%elements(icel)%gausses(jj)%fstatus(1))
                      enddo
                    end select
                    
                enddo
              enddo

			
         end if

       enddo
	   
	   ! Summary & Result
       if( ndof == 3 ) call fstr_NodalStress3D( hecMESH, hecMAT, fstrSOLID, ndstrain, ndstress  )
       call fstrNLGEOM_Post(fnum,hecMESH,fstrSOLID,ttime,totstep,ndstrain,ndstress)
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
  subroutine fstrNLGEOM_Post(fnum,hecMESH,fstrSOLID,tt,i_step,ndstrain, ndstress)
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(in)                   :: fnum
      type (hecmwST_local_mesh), intent(in) :: hecMESH
      type (fstr_solid), intent(in)         :: fstrSOLID

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
      integer(kind=kint) j,i,k,nd,id,ndof,mdof
	  
!
! --- Time
      write(fnum,'(''#### Result time='',E10.4)') tt
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

        if( IRESULT.eq.0 ) return
! --- Write Result File
! --- INITIALIZE
        header='*fstrresult'
        nd = i_step
        call hecmw_result_init(hecMESH,nd,header)
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
