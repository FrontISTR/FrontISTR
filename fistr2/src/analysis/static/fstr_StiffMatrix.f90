!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by  X. YUAN(AdavanceSoft)                         !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides function to calcualte tangent stiffness matrix.
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2010/08/26
!>  \version    0.00
!!
!======================================================================!
module m_fstr_StiffMatrix
   implicit none

   contains

!---------------------------------------------------------------------*
!> \brief 接線剛性マトリックスを作成するサブルーチン
subroutine fstr_StiffMatrix ( fstrSOLID, cstep, tincr, factor)
!---------------------------------------------------------------------*
  use m_fstr
  use m_static_LIB
  use mMechGauss
  
  include "HEC_MW3_For.h"

  type (fstr_solid)          :: fstrSOLID    !< we need boundary conditions of curr step
  integer, intent(in)        :: cstep        !< current step
  real(kind=kreal),intent(in) :: tincr       !< time increment
  real(kind=kreal),intent(in) :: factor      !< loading factor
  
  type( tMaterial ), pointer :: material     !< material information

  real(kind=kreal)   :: stiffness(20*6, 20*6)
  real(kind=kreal)   :: tt(20), ecoord(3,20)
  real(kind=kreal)   :: thick, pa1
  integer(kind=kint) :: cid, ndof, itype, iS, iE, ic_type, nn, icel, iiS, i
  real(kind=kreal)   :: u(3,20), du(3,20)
  integer            :: ig0, grpid, ig, iS0, iE0,ik

  real(kind=kreal) :: x,y,z
  integer(kind=kint) :: nids(0:20)
  integer :: iAss, iPart, iElem, iNode, iGrp, iErr

  ndof = assDOF(1)
	  
  icel = 0	  
  do iAss = 0, mw_get_num_of_assemble_model()-1
     call mw_select_assemble_model( iAss )
     do iPart = 0, mw_get_num_of_mesh_part()-1
        call mw_select_mesh_part( iPart )
        call mw_matrix_clear( iPart )
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
                cid = part_nodes(iAss+1,iPart+1)+mw_get_node_index( cid )+1
                do i=1,ndof
                  du(i,iNode) = fstrSOLID%dunode(ndof*cid+i-ndof)
                  u(i,iNode)  = fstrSOLID%unode(ndof*cid+i-ndof) + du(i,iNode)
                enddo
                if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres == 1 )  &
                    tt(iNode)=fstrSOLID%temperature( cid )  
              enddo
		  
          material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
          thick = material%variables(M_THICK)
          if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
          if ( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
            call STF_C2( ic_type,nn,ecoord(1:2,1:nn),fstrSOLID%elements(icel)%gausses(:),thick,  &
                     stiffness(1:nn*ndof,1:nn*ndof), fstrSOLID%elements(icel)%iset,          &
                     u(1:2,1:nn) )
					 
          else if ( ic_type==361 ) then
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres == 1 ) then
              call STF_C3D8Bbar( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, u(1:3,1:nn), tt(1:nn) )
            else
              call STF_C3D8Bbar( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, u(1:3,1:nn) )
            endif
        

          else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.                    &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres == 1 ) then
              call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, u(1:3,1:nn), tt(1:nn) )
            else
              call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, u(1:3,1:nn) )
            endif

          else
            write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
            write(*,*) ' ic_type = ', ic_type
            call hecmw_abort(hecmw_comm_get_comm())
          endif
          iErr = mw_matrix_add_elem( iPart, iElem, stiffness(1:nn*ndof,1:nn*ndof) )
        enddo
     enddo
  enddo  

end subroutine fstr_StiffMatrix

end module m_fstr_StiffMatrix
