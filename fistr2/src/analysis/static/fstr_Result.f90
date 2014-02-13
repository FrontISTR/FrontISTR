!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
  subroutine fstr_OutputResult( cstep, totstep, fstrSOLID, ttime )
!----------------------------------------------------------------------*
    use m_fstr
    use m_static_make_result
    use m_static_output
    use fstr_setup_util
    use m_fstr_visualize

    integer, intent(in)                   :: cstep       !< current step number
    integer, intent(in)                   :: totstep     !< total steps
    type (fstr_solid), intent(in)         :: fstrSOLID
    real(kind=kreal),intent(in)           :: ttime

    real(kind=kreal)   :: fdum, ndforce(6), ss(6)
    integer(kind=kint) :: i, k, ii, jj, nctrl, ninfo, ndof, ncomp, fnum
    integer(kind=kint) :: maxiter, itype, iS, iE, ic_type, icel
    real(kind=kreal)   :: ndstrain(total_node*6), ndstress(total_node*7)
    real(kind=kreal)   :: s11,s22,s33,s12,s23,s13,ps,smises, vdum(6)
    integer(kind=kint),pointer :: member(:)
    type ( hecmwST_result_data ) :: fstrRESULT

    maxiter = fstrSOLID%step_ctrl(cstep)%max_iter
	!    call hecmw_nullify_result_data( fstrRESULT )
    !-- POST PROCESSING VIA MEMORY
    if( IVISUAL==1 ) then
          call fstr_NodalStress3D( fstrSOLID, fstrSOLID%STRAIN, fstrSOLID%STRESS)
      !    call fstr_make_result(fstrSOLID,fstrRESULT)
      !    call hecmw_visualize_init

      !   call hecmw_visualize(hecMESH,fstrRESULT,totstep,maxiter,0)
      !   call hecmw_visualize_finalize
      !    call hecmw_result_free(fstrRESULT)
          call fstr_visualize( cstep, totstep, fstrSOLID, ttime ) 
    endif
	
    if( .not. associated(fstrSOLID%output_ctrl) ) return

    ndof = assDOF(1)
    do nctrl=1, size(fstrSOLID%output_ctrl)
       if( .not. fstr_output_active( totstep, fstrSOLID%output_ctrl(nctrl) ) ) cycle
	   

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
            do i= 1, total_node
             write(fnum,'(i10,1p6e15.7)') i,(fstrSOLID%unode(ndof*i-ndof+k),k=1,ndof)
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
              do i=1,total_node
                ndforce(1:ndof)                                                &
                  = fstrSOLID%QFORCE( ndof*(i-1)+1 : ndof*(i-1)+ndof)
                write(fnum,'(i10,1p6e15.7)') i, ndforce(1:ndof)
              enddo

			  

! ------- NODAL STRAIN
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRAIN"   &
             .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) == -1) then
            if( ndof == 3 ) then
              call fstr_NodalStress3D( fstrSOLID, ndstrain, ndstress  )
              write(fnum,*) '#### NODAL STRAIN'
              write(fnum,'(a,a)') '     NODE       E11            E22            E33      ',   &
                                    '      E12            E23            E13      '
              write(fnum,'(a,a)') '---------+--------------+--------------+--------------+',   &
                                    '--------------+--------------+--------------+'
              do i=1,total_node
                 write(fnum,'(i10,1p6e15.7)') i, (ndSTRAIN(6*(i-1)+k),k=1,6)
              enddo
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif

! ------- NODAL STRESS
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRESS"   &
            .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) == -1) then
            if( ndof == 3 ) then
              call fstr_NodalStress3D( fstrSOLID, ndstrain, ndstress  )
              write(fnum,*) '#### NODAL STRESS'
              write(fnum,'(a,a)') '     NODE       S11            S22            S33      ',   &
                                    '      S12            S23            S13      '
              write(fnum,'(a,a)') '---------+--------------+--------------+--------------+',   &
                                    '--------------+--------------+--------------+'
              do i=1,total_node
                write(fnum,'(i10,1p7e15.7)') i, (ndSTRESS(7*(i-1)+k),k=1,7)
              enddo
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif
			
! ---------ELMENTAL STRAIN
        else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRAIN"   &
             .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) >= 0) then
            if( ndof == 3 ) then
              write(fnum,*) '#### ELEMENTAL STRAIN'
              do i= 1, total_elem
                ii = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
                write(fnum,*) '--- ELEMENT:', i
                select case ( fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) )
                case (1)  ! average
                      ss(:) = 0.d0
                      do jj=1,ii
                         ss(:) = ss(:) + fstrSOLID%elements(i)%gausses(jj)%strain(:)
                      enddo
                      ss(:) = ss(:)/ii
                      write(fnum,'(1p6e15.7)') (ss(k),k=1,6)
                case (2)  ! quadrature point
                      do jj=1,ii
                        write(fnum,'(i10,1p6e15.7)') jj, (fstrSOLID%elements(i)%gausses(jj)%strain(k),k=1,6)
                      enddo
                end select
              enddo
			  
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif

! ------- ELEMENTAL STRESS
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "STRESS"   &
            .and. fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) >= 0) then
            if( ndof == 3 ) then
              write(fnum,*) '#### ELEMENTAL STRESS'
              do i= 1, total_elem
                ii = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
                write(fnum,*) '--- ELEMENT:', i
                select case ( fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) )
                case (1)  ! average
                      ss(:) = 0.d0
                      do jj=1,ii
                         ss(:) = ss(:) + fstrSOLID%elements(i)%gausses(jj)%stress(:)
                      enddo
                      ss(:) = ss(:)/ii
					  s11=ss(1)
                      s22=ss(2)
                      s33=ss(3)
                      s12=ss(4)
                      s23=ss(5)
                      s13=ss(6)
                      ps=(s11+s22+s33)/3.d0
                      smises=dsqrt( 0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                          +s12**2+s23**2+s13**2 )
                      write(fnum,'(1p7e15.7)') (ss(k),k=1,6),dsqrt(3.d0*smises)
                case (2)  ! quadrature point
                      do jj=1,ii
					    s11=fstrSOLID%elements(i)%gausses(jj)%stress(1)
                        s22=fstrSOLID%elements(i)%gausses(jj)%stress(2)
                        s33=fstrSOLID%elements(i)%gausses(jj)%stress(3)
                        s12=fstrSOLID%elements(i)%gausses(jj)%stress(4)
                        s23=fstrSOLID%elements(i)%gausses(jj)%stress(5)
                        s13=fstrSOLID%elements(i)%gausses(jj)%stress(6)
                        ps=(s11+s22+s33)/3.d0
                        smises=dsqrt( 0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                            +s12**2+s23**2+s13**2 )
                        write(fnum,'(i10,1p7e15.7)') jj, (fstrSOLID%elements(i)%gausses(jj)%stress(k),k=1,6),dsqrt(3.d0*smises)
                      enddo
                end select
              enddo
            else if( ndof==2 ) then
                       
            else if( ndof==6 ) then
                       
            endif
			
! ------- ELEMENTAL PLASTIC STRAIN
         else if ( fstrSOLID%output_ctrl(nctrl)%outinfo%keyWord(ninfo) == "PLASTIC_STRAIN" ) then
              write(fnum,*) '#### PLASTIC STRAIN'
              write(fnum,*) '#### ELEMENTAL STRESS'
              do i= 1, total_elem
                ii = NumOfQuadPoints( fstrSOLID%elements(i)%etype )
                write(fnum,*) '--- ELEMENT:', i
                select case ( fstrSOLID%output_ctrl(nctrl)%outinfo%location(ninfo) )
                case (1)  ! average
                      fdum = 0.d0
                      do jj=1,ii
                        if( isElastoplastic( fstrSOLID%elements(i)%gausses(jj)%pMaterial%mtype )   &
                          .and. associated(fstrSOLID%elements(i)%gausses(jj)%fstatus) )           &
                         fdum = fdum + fstrSOLID%elements(i)%gausses(jj)%fstatus(1)
                      enddo
                      write(fnum,'(1p6e15.7)') fdum/ii
                case (2)  ! quadrature point
                      do jj=1,ii
                        if( isElastoplastic( fstrSOLID%elements(i)%gausses(jj)%pMaterial%mtype )   &
                          .and. associated(fstrSOLID%elements(i)%gausses(jj)%fstatus) )             &
                        write(fnum,'(i10,1p1e15.7)') jj, (fstrSOLID%elements(i)%gausses(jj)%fstatus(1))
                      enddo
                end select
              enddo

			
         end if

       enddo
	   
	   ! Summary & Result
       if( ndof == 3 ) call fstr_NodalStress3D( fstrSOLID, ndstrain, ndstress  )
       call fstrNLGEOM_Post(fnum,fstrSOLID,ttime,totstep,ndstrain,ndstress)
    enddo
  end subroutine fstr_OutputResult

!> Summarizer of output data which prints out max and min output values
!----------------------------------------------------------------------*
  subroutine fstrNLGEOM_Post(fnum,fstrSOLID,tt,i_step,ndstrain, ndstress)
!----------------------------------------------------------------------*
      use m_fstr
      integer, intent(in)                   :: fnum
      type (fstr_solid), intent(in)         :: fstrSOLID
	  
      include "HEC_MW3_For.h"

      integer(kind=kint), intent(in) :: i_step
      real(kind=kreal), intent(in)   :: tt
      real(kind=kreal), intent(in)   :: ndstrain(total_node*6)
      real(kind=kreal), intent(in)   :: ndstress(total_node*7)
! --- file name
      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN) :: label
      character(len=HECMW_NAME_LEN) :: nameID
! --- Local variables
! --- Max Min
      real(kind=kreal), target :: Umax(3),Umin(3),Emax(6),Emin(6),Smax(7),Smin(7)
      integer(kind=kint), target :: IUmax(3),IUmin(3),IEmax(6),IEmin(6),ISmax(7),ISmin(7)
      integer, target :: val_size, op
!
      integer(kind=kint) j,i,k,nd,id,ndof,mdof
!
! --- Time
      write(fnum,'(''#### Resul time='',E11.4)') tt
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
      ndof = assDOF(1)
!
! --- Displacement
!      write(fnum,*) '#### DISPLACEMENT@POST'
      do i= 1, total_node
        if( i==1) then
          do k=1,NDOF 
            Umax(k)=fstrSOLID%unode(NDOF*(i-1)+k)
            Umin(k)=fstrSOLID%unode(NDOF*(i-1)+k)
            IUmax(k)=i
            IUmin(k)=i
          enddo
        else
          do k=1,NDOF 
            if( fstrSOLID%unode(NDOF*(i-1)+k) > Umax(k) ) then
              Umax(k)=fstrSOLID%unode(NDOF*(i-1)+k)
              IUmax(k)=i
            endif
            if( fstrSOLID%unode(NDOF*(i-1)+k) < Umin(k) ) then
              Umin(k)=fstrSOLID%unode(NDOF*(i-1)+k)
              IUmin(k)=i
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
      do i= 1, total_node
! --- Max & Min
        if( i==1 ) then
          do k=1,Mdof
            Emax(k)=ndSTRAIN((i-1)*Mdof+k)
            Emin(k)=ndSTRAIN((i-1)*Mdof+k)
            IEmax(k)=i
            IEmin(k)=i
          enddo
        else
          do k=1,Mdof
            if( ndSTRAIN((i-1)*Mdof+k) > Emax(k) ) then
              Emax(k)=ndSTRAIN((i-1)*Mdof+k)
              IEmax(k)=i
            endif
            if( ndSTRAIN((i-1)*Mdof+k) < Emin(k) ) then
              Emin(k)=ndSTRAIN((i-1)*Mdof+k)
              IEmin(k)=i
            endif
          enddo
        endif
      enddo
! --- Stress
!      write(fnum,*) '#### STRESS@POST'
      do i= 1, total_node
! --- Max & Min
        if( i.eq.1) then
          do k=1,7
            Smax(k)=ndSTRESS((i-1)*7+k)
            Smin(k)=ndSTRESS((i-1)*7+k)
            ISmax(k)=i
            ISmin(k)=i
          enddo
        else
          do k=1,7
            if( ndSTRESS((i-1)*7+k) > Smax(k) ) then
              Smax(k)=ndSTRESS((i-1)*7+k)
              ISmax(k)=i
            endif
            if( ndSTRESS((i-1)*7+k) < Smin(k) ) then
              Smin(k)=ndSTRESS((i-1)*7+k)
              ISmin(k)=i
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
        call mw_allreduce_r( Umax, 3, mw_mpi_max() )
        call mw_allreduce_r( Umin, 3, mw_mpi_min() )
        call mw_allreduce_r( Emax, 6, mw_mpi_max() )
        call mw_allreduce_r( Emin, 6, mw_mpi_min() )
        call mw_allreduce_r( Smax, 6, mw_mpi_max() )
        call mw_allreduce_r( Smin, 6, mw_mpi_min() )

        if( myrank == 0 ) then
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
    !    call hecmw_result_init(nd,header, total_node, total_elem, global_node_ID, global_elem_ID)
! --- ADD
! --- DISPLACEMENT
!        if( fstrSOLID%iout_list(1) .eq. 1 ) then
          id = 1
          ndof=3
          label='DISPLACEMENT'
    !      call hecmw_result_add(id,ndof,label,fstrSOLID%unode)
!        end if
!
        if( i_step .gt. 0) then
!        if( fstrSOLID%iout_list(5) .eq. 1 ) then
! --- STRAIN @node
          id = 1
          ndof=6
          label='STRAIN'
     !     call hecmw_result_add(id,ndof,label,fstrSOLID%STRAIN)
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
       !   call hecmw_result_add(id,ndof,label,ndSTRESS)
! --- STRESS @element
!          id = 2
!          ndof=7
!          label='ESTRESS'
!          call hecmw_result_add(id,ndof,label,fstrSOLID%ESTRESS)
!        end if
        end if
!
        nameID='fstrRES'
      !  call hecmw_result_write_by_name(nameID)
! --- FINALIZE
     !   call hecmw_result_finalize
!
      endif

  end subroutine fstrNLGEOM_Post
!
!
!
end module m_fstr_Result
