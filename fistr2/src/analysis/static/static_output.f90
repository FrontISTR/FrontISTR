!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by Xi YUAN (Advancesoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!C***
!> OUTPUT for FSTR solver
!C***
module m_static_output

   use m_fstr
   use m_static_LIB_shell
   use m_static_mat_ass_main
   use elementInfo
   implicit none

   contains
   
    subroutine solid_output(fstrSOLID)
      type (fstr_solid)        :: fstrSOLID 
  
      if( assDOF(1)==3 ) then
        call fstr_output_3d(fstrSOLID)
      else if( assDOF(1)==2 ) then
        call fstr_output_2d(fstrSOLID)
      else if( assDOF(1)==6) THEN
        call fstr_output_6d(fstrSOLID)
      endif 
    end subroutine
   
subroutine fstr_Update( fstrSOLID, substep,tincr,factor,iter)
!=====================================================================*
  use m_fstr
  use m_static_lib

  type (fstr_solid)           :: fstrSOLID !< 
  integer, intent(in)         :: substep   !< current substep
  real(kind=kreal),intent(in) :: tincr     !< time increment
  real(kind=kreal),intent(in) :: factor    !< loading factor
  integer, intent(in)         :: iter      !< NR iterations
  
  include "HEC_MW3_For.h"
	  
  real(kind=kreal) :: x,y,z, stiff(20*6, 20*6), local_stf(1830)
  integer(kind=kint) :: nids(0:20), iwk(60)
  integer(kind=kint) :: itype, cid, ndID, i, j, k, icel
  integer :: iAss, iPart, iElem, iNode, iGrp, iErr, iSect
  REAL(kind=kreal) edisp(6,20),force(60),ecoord(3,20)

  real(kind=kreal)   :: ee,pp,thick, fval, pa1, ri(4), si(4)
  integer(kind=kint) :: ndof, iS, iE, ic_type, nn, iiS

  real(kind=kreal)   :: total_disp(6,20), array(11)
  real(kind=kreal)   :: tt(20), tt0(20), qf(20*3)
  integer            :: ig0, grpid, ig, iS0, iE0,ik,  ierror
  
  data ri / -1.0,  1.0, 1.0, -1.0 /
  data si / -1.0, -1.0, 1.0,  1.0 /

  ndof = assDOF(1)
  fstrSOLID%QFORCE=0.0d0
  tt0=0.d0;  tt=0.d0
  total_disp(:,:) = 0.d0
  
  
! --------------------------------------------------------------------
!      updated 
!         1. stress and strain  : ep^(k) = ep^(k-1)+dep^(k)
!                                 sgm^(k) = sgm^(k-1)+dsgm^(k)
!         2. Internal Force     : Q^(k-1) ( u^(k-1) )
! --------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------
!      calculate the Strain and Stress and Internal Force ( Equivalent Nodal Force )
! ----------------------------------------------------------------------------------
!
  icel = 0	  
  do iAss = 0, mw_get_num_of_assemble_model()-1
     call mw_select_assemble_model( iAss )
     call mw_select_algebra(0)
     do iPart = 0, mw_get_num_of_mesh_part()-1
        call mw_select_mesh_part( iPart )
        do iElem = 0, mw_get_num_of_element()-1
          icel = icel+1
          call mw_select_element( iElem )
          call  mw_get_element_vert_node_id( nids )
          nn = mw_get_num_of_element_vert()
          ic_type = mw_get_element_type()
          ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
          do iNode = 1, nn
            cid = nids(iNode-1)
            call mw_get_node_coord( cid, x,y,z )
            ecoord(1,iNode)=x
            ecoord(2,iNode)=y
            ecoord(3,iNode)=z
            cid= part_nodes(iAss+1,iPart+1)+cid+1
            if( associated(fstrSOLID%temp_grp) .or. fstrSOLID%TEMP_irres == 1 ) then
                     tt0(j)=fstrSOLID%reftemp( cid ) 
                     tt(j) = fstrSOLID%temperature( cid ) 
            endif	
            do i=1,ndof    
               iwk(ndof*(iNode-1)+i)=ndof*(cid-1)+i
               edisp( i,iNode ) = fstrSOLID%unode((cid-1)*ndof+i)
               force( (iNode-1)*ndof+i ) = fstrSOLID%unode((cid-1)*ndof+i)
            enddo
          enddo
			  
          if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
          fstrSOLID%elements(icel)%gausses(1)%pMaterial%nlgeom_flag = INFINITE
          if ( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
            call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:),      &
                        thick,fstrSOLID%elements(icel)%iset,                                  &
                        total_disp(1:2,1:nn), edisp(1:2,1:nn), qf(1:nn*ndof) )
						
          else if ( ic_type==361 ) then
            call STF_C3D8IC( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3))
            qf(1:nn*3)= matmul( stiff(1:nn*3,1:nn*3), force(1:nn*3) )
            call UpdateST_C3D8IC(ic_type,nn,ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),    &
                  tt(1:nn), tt0(1:nn), edisp(1:3,1:nn), fstrSOLID%elements(icel)%gausses(:))
        
          else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.                          &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres == 1 ) then
              call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), edisp(1:3,1:nn), qf(1:nn*ndof)       &
                        ,fstrSOLID%elements(icel)%gausses(:), substep, iter, tincr, tt(1:nn), tt0(1:nn)  )
            else
              call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), edisp(1:3,1:nn)       &
                        , qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), substep, iter, tincr  )
            endif

         else if ( ic_type==731) then
            ee = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_YOUNGS)
            pp = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_POISSON)
            isect = fstrSOLID%elements(icel)%isect
            thick = MWSections(isect)%sect_R_item(1)
            call STF_S3( ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),ee,pp,thick,local_stf )
            call fstr_local_stf_restore(local_stf, nn*ndof, stiff)
            qf(1:nn*6)= matmul( stiff(1:nn*6,1:nn*6), force(1:nn*6) )
            call RCV_S3( ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),ee,pp,thick,force,1.d0,array )
            fstrSOLID%elements(icel)%gausses(1)%strain(1:5) = array(1:5)
            fstrSOLID%elements(icel)%gausses(1)%stress = array(6:11)
            call RCV_S3( ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),ee,pp,thick,force,-1.d0,array )
            fstrSOLID%elements(icel)%gausses(2)%strain(1:5) = array(1:5)
            fstrSOLID%elements(icel)%gausses(2)%stress = array(6:11)
         else if ( ic_type==741) then
            ee = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_YOUNGS)
            pp = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_POISSON)
            isect = fstrSOLID%elements(icel)%isect
            thick = MWSections(isect)%sect_R_item(1)
            call STF_S4( ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),ee,pp,thick,local_stf )
            call fstr_local_stf_restore(local_stf, nn*ndof, stiff)
            qf(1:nn*6)= matmul( stiff(1:nn*6,1:nn*6), force(1:nn*6) )
            do j = 1, 4
               call RCV_S4( ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),ee,pp,thick,                  &
                                  force,ri(j),si(j),1.d0,array )
               fstrSOLID%elements(icel)%gausses(j)%strain(1:5) = array(1:5)
               fstrSOLID%elements(icel)%gausses(j)%stress = array(6:11)
            enddo
            do j = 1, 4
               call RCV_S4( ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),ee,pp,thick,                  &
                                  force,ri(j),si(j),-1.d0,array )
               fstrSOLID%elements(icel)%gausses(j+4)%strain(1:5) = array(1:5)
               fstrSOLID%elements(icel)%gausses(j+4)%stress = array(6:11)
            enddo
         else
           write(*,*) '###ERROR### : Element type not supported for linear static analysis'
           write(*,*) ' ic_type = ', ic_type
           call hecmw_abort(hecmw_comm_get_comm())
         endif
!
! ----- calculate the global internal force ( Q(u_{n+1}^{k-1}) )
           do j=1,nn*ndof
              fstrSOLID%QFORCE( iwk(j) )=fstrSOLID%QFORCE( iwk(j) )+qf(j)
           enddo 
        enddo
     enddo
  enddo

!C
!C Update for fstrSOLID%QFORCE
!C 
  do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part_with_id( iPart )
            do iNode = 0, mw_get_num_of_neibpe(iPart)-1
               ik = mw_get_transrank(iPart, iNode)
               ndID = part_nodes(iAss+1,iPart+2) - part_nodes(iAss+1,iPart+1)
               call mw_send_recv_r(fstrSOLID%QFORCE(part_nodes(iAss+1,iPart+1)+1:part_nodes(iAss+1,iPart+2)), ndID, ndof, ik)
            enddo
         enddo
  enddo
  
end subroutine fstr_Update

!C
!< ouput for 3D SOLID
!C
   subroutine fstr_output_3d(fstrSOLID)
      type (fstr_solid)         :: fstrSOLID
	  
      include "HEC_MW3_For.h"
	  
      integer :: i,j,k
      real(kind=kreal) :: s11,s22,s33,s12,s23,s13,smises, ps
	  
      write(ILOG,*) ""
      write(ILOG,*) ""
      write(ILOG,*) ""

!C
!C-- DISPLACEMENT
!C
      write(ILOG,*) '#### DISPLACEMENT 3D'
      write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      Z-DISP'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, total_node
        write(ILOG,'(i10,3e12.4)') i,(fstrSOLID%unode(3*(i-1)+k),k=1,3)
      enddo
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 3D'
      write(ILOG,'(a)')'    NODE     X-REAC      Y-REAC      Z-REAC'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, total_node
        write(ILOG,'(i10,3e12.4)') i,(fstrSOLID%QFORCE(3*(i-1)+k),k=1,3)
      enddo
	  
!C
!C-- ELEMENTAL STRESS
!C
      write(ILOG,*) '#### ELEMENTAL STRAIN 3D'
      write(ILOG,'(a,a)')                                              &
                  '     ELEMENT    QP     S11         S22         S33     '   &
                 ,'    S12         S23         S13     '
      write(ILOG,'(a,a)')                                              &
                    '  -------+-----------+-----------+-----------+'   &
                 ,'-----------+-----------+-----------+-------------'
      do i= 1, total_elem
        do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          write(ILOG,'(2i10,1p6e15.7)') i,j, (fstrSOLID%elements(i)%gausses(j)%strain(k),k=1,6)
        enddo
      enddo

!C
!C-- ELEMENTAL STRESS
!C
     write(ILOG,*) '#### ELEMENTAL STRESS 3D'
      write(ILOG,'(a,a)')                                              &
                  '     ELEMENT    QP     S11         S22         S33     '   &
                 ,'    S12         S23         S13        MISES'
      write(ILOG,'(a,a)')                                              &
                    '  -------+-----------+-----------+-----------+'   &
                 ,'-----------+-----------+-----------+-------------'
      do i= 1, total_elem
        do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          s11=fstrSOLID%elements(i)%gausses(j)%stress(1)
          s22=fstrSOLID%elements(i)%gausses(j)%stress(2)
          s33=fstrSOLID%elements(i)%gausses(j)%stress(3)
          s12=fstrSOLID%elements(i)%gausses(j)%stress(4)
          s23=fstrSOLID%elements(i)%gausses(j)%stress(5)
          s13=fstrSOLID%elements(i)%gausses(j)%stress(6)
          ps=(s11+s22+s33)/3.d0
          smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                  +s12**2+s23**2+s13**2
          smises=dsqrt(3.d0*smises)
          write(ILOG,'(2i10,1p7e15.7)') i,j, (fstrSOLID%elements(i)%gausses(j)%stress(k),k=1,6),smises
        enddo
      enddo

      end subroutine fstr_output_3d

!C
!< output for 2D SOLID
!C
      subroutine FSTR_OUTPUT_2D(fstrSOLID)
      type (fstr_solid)         :: fstrSOLID
	  
      integer :: i,j,k
      real(kind=kreal) :: s11,s22,s33,s12,s23,s13,smises, ps
	  
      write(ILOG,*) ""
      write(ILOG,*) ""
      write(ILOG,*) ""
!C
!C-- DISPLACEMENT
!C
      write(ILOG,*) '#### DISPLACEMENT 2D'
      write(ILOG,'(a)')                                      &
                    '     NODE    X-DISP      Y-DISP   '
      write(ILOG,'(a)')                                      &
                    ' --------+-----------+------------'
      do i= 1, total_node
        write(ILOG,'(i10,2e12.4)') i,(fstrSOLID%unode(2*(i-1)+k),k=1,2)
      enddo
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 2D'
      do i= 1, total_node
        write(ILOG,'(i10,2e12.4)') i,(fstrSOLID%QFORCE(2*(i-1)+k),k=1,2)
      enddo
!C
!C-- ELEMENTAL STRESS
!C
      write(ILOG,*) '#### ELEMENTAL STRAIN 2D'
      write(ILOG,'(a,a)')                                              &
               '   ELEMENT    QP       E11         E22         E33     '        &
                             ,'    E12' 
      write(ILOG,'(a,a)')                                              &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-------------'
      do i= 1, total_elem
        do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          write(ILOG,'(2i10,1p4e15.7)') i,j, (fstrSOLID%elements(i)%gausses(j)%strain(k),k=1,4)
        enddo
      enddo

!C
!C-- ELEMENTAL STRESS
!C
      write(ILOG,*) '#### ELEMENTAL STRESS 2D'
      write(ILOG,'(a,a)')                                              &
               '   ELEMENT    QP       S11         S22         S33     '        &
                             ,'    S12        MISES' 
      write(ILOG,'(a,a)')                                              &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-----------+-------------'
      do i= 1, total_elem
        do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          s11=fstrSOLID%elements(i)%gausses(j)%stress(1)
          s22=fstrSOLID%elements(i)%gausses(j)%stress(2)
          s33=fstrSOLID%elements(i)%gausses(j)%stress(3)
          s12=fstrSOLID%elements(i)%gausses(j)%stress(4)
          s23=0.d0
          s13=0.d0
          ps=(s11+s22+s33)/3.d0
          smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                  +s12**2+s23**2+s13**2
          smises=dsqrt(3.d0*smises)
          write(ILOG,'(2i10,1p5e15.7)') i,j, (fstrSOLID%elements(i)%gausses(j)%stress(k),k=1,4),smises
        enddo
      enddo

      end subroutine FSTR_OUTPUT_2D

!C***
!< OUTPUT for shell problem
!C***
      subroutine fstr_output_6d( fstrSOLID )
      type ( fstr_solid         ) :: fstrSOLID
	  
      integer :: i,j,k

      write(ILOG,*) ""
      write(ILOG,*) ""
      write(ILOG,*) ""

!C
      write(ILOG,*) '#### DISPLACEMENT SHELL'
      write(ILOG,'(a,a)')                                            & 
                    '     NODE    X-DISP      Y-DISP      Z-DISP   ' &
                             ,'   X-ROT       Y-ROT       Z-ROT'
      write(ILOG,'(a,a)')                                            & 
                    ' --------+-----------+-----------+-----------+' &
                             ,'-----------+-----------+-------------'
      do i= 1, total_node
        write(ILOG,'(i10,6e12.4)') i,(fstrSOLID%unode(6*(i-1)+k),k=1,6)
      enddo

      call flush(ILOG)

!C
!C-- REACTION FORCE
  
      write(ILOG,*) '#### REACTION FORCE SHELL'
      write(ILOG,'(a,a)')                                              &
                    '   NODE      X-REAC      Y-REAC      Z-REAC   '   &
                             ,'  RX-REAC     RY-REAC     RZ-REAC'
      write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      do i= 1, total_node
        write(ILOG,'(i10,6e12.4)') i,(fstrSOLID%QFORCE(6*(i-1)+k),k=1,6)
      enddo

      call flush(ILOG)

	  
      write(ILOG,*) '#### ELEMENTAL STRESS SHELL'
      write(ILOG,'(a,a)')                                             & 
               '     ELEMENT    QP          E11         E22         E12     '  &
                             ,'    E23         E13'
      write(ILOG,'(a,a)')                                             & 
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-----------+-------------'
      do i= 1, total_elem
        do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          write(ILOG,'(2i10,1p6e15.7)') i,j, (fstrSOLID%elements(i)%gausses(j)%stress(k),k=1,6)
        enddo
      enddo
	  
	  write(ILOG,*) '#### ELEMENTAL STRAIN SHELL'
      write(ILOG,'(a,a)')                                             & 
               '     ELEMENT    QP          E11         E22         E12     '  &
                             ,'    E23         E13'
      write(ILOG,'(a,a)')                                             & 
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-----------+-------------'
      do i= 1, total_elem
        do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
          write(ILOG,'(2i10,1p6e15.7)') i,j, (fstrSOLID%elements(i)%gausses(j)%strain(k),k=1,5)
        enddo
      enddo 

      call flush(ILOG)

      end subroutine fstr_output_6d

end module m_static_output
