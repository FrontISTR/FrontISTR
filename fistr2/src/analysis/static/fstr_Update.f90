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
!> \brief  This module provides functions to do update.
!!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \date       2009/08/28
!>  \version    0.00
!!
!======================================================================!
module m_fstr_Update
  implicit none

  contains

!=====================================================================*
!  UPDATE_C3
!>  \brief 変位／応力・ひずみ／内力のアップデート
!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \version    0.00
!!
!> \if AS
!>    \par  サブルーチン構成
!>    -# 変位の更新					\f$ u_{n+1}^{(k)} = u_{n+1}^{(k-1)} + \delta u^{(k)} \f$
!>    -# ひずみ・応力の更新			\f$ \varepsilon_{n+1}^{(k)} = \varepsilon_{n+1}^{(k-1)} + \delta \varepsilon^{(k)} \f$, \f$ \sigma_{n+1}^{(k)} = \sigma_{n+1}^{(k-1)} + \delta \sigma^{(k)} \f$
!>    -# 内力（等価節点力）の計算	\f$ Q_{n+1}^{(k-1)} ( u_{n+1}^{(k-1)} ) \f$
!> \endif
!
subroutine fstr_UpdateNewton ( fstrSOLID, substep,tincr,factor,iter)
!=====================================================================*
  use m_fstr
  use m_static_lib
  use m_static_LIB_3d
  use m_static_LIB_C3D8

  type (fstr_solid)           :: fstrSOLID !< we need boundary conditions of curr step
  integer, intent(in)         :: substep   !< current substep
  real(kind=kreal),intent(in) :: tincr     !< time increment
  real(kind=kreal),intent(in) :: factor    !< loading factor
  integer, intent(in)         :: iter      !< NR iterations
  
  include "HEC_MW3_For.h"
	  
  real(kind=kreal) :: x,y,z, stiffness(20*6, 20*6)
  integer(kind=kint) :: nids(0:20), iwk(60)
  integer(kind=kint) :: itype, cid, ndID, i, j, k, icel
  integer :: iAss, iPart, iElem, iNode, iGrp, iErr
  REAL(kind=kreal) force(60),ecoord(3,20)

  real(kind=kreal)   :: thick, fval, pa1
  integer(kind=kint) :: ndof, iS, iE, ic_type, nn, iiS

  real(kind=kreal)   :: total_disp(3,20),ddu(3,20)
  real(kind=kreal)   :: tt(20), tt0(20), qf(20*3)
  integer            :: ig0, grpid, ig, iS0, iE0,ik, in, ierror, snode, enode

  ndof = assDOF(1)
  fstrSOLID%QFORCE=0.0d0
  
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
     do iPart = 0, mw_get_num_of_mesh_part()-1
        snode = part_nodes(iAss+1, iPart+1)
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
            cid = mw_get_node_index( cid )+1+snode
            if( associated(fstrSOLID%temp_grp) .or. fstrSOLID%TEMP_irres == 1 ) then
                     tt0(iNode)=fstrSOLID%reftemp( cid ) 
                     tt(iNode) = fstrSOLID%temperature( cid ) 
            endif	
            do i=1,ndof    
                   iwk(ndof*(iNode-1)+i)=ndof*(cid-1)+i
                   ddu(i,iNode) = fstrSOLID%ddunode(ndof*(cid-1)+i)
                   total_disp(i,iNode) = fstrSOLID%unode(ndof*(cid-1)+i) +  &
                                       fstrSOLID%dunode(ndof*(cid-1)+i)
            enddo
	!		print *, iAss, icel, iNode, total_disp(:,iNode)
          enddo

          if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
          if ( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
            call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:),      &
                        thick,fstrSOLID%elements(icel)%iset,                                  &
                        total_disp(1:2,1:nn), ddu(1:2,1:nn), qf(1:nn*ndof) )
						
          else if ( ic_type==361 ) then
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres == 1 ) then
              call UPDATE_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), ddu(1:3,1:nn), qf(1:nn*ndof)       &
                        ,fstrSOLID%elements(icel)%gausses(:), substep, iter, tincr, tt(1:nn), tt0(1:nn)  )
            else
              call Update_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), ddu(1:3,1:nn)       &
                        , qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), substep, iter, tincr  )
            endif
        
         else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.                          &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres == 1 ) then
              call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), ddu(1:3,1:nn), qf(1:nn*ndof)       &
                        ,fstrSOLID%elements(icel)%gausses(:), substep, iter, tincr, tt(1:nn), tt0(1:nn)  )
            else
              call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), ddu(1:3,1:nn)       &
                        , qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), substep, iter, tincr  )
           endif

         else if ( ic_type==731) then
      !     call UPDATE_S3(xx,yy,zz,ee,pp,thick,local_stf)
!        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
         else if ( ic_type==741) then
     !      call UPDATE_S4(xx,yy,zz,ee,pp,thick,local_stf)
!        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
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
            call mw_select_mesh_part( iPart )
            do iNode = 0, mw_get_num_of_neibpe(iPart)-1
           !    ik = mw_get_transrank(iPart, iNode)
               snode = part_nodes(iAss+1,iPart+1)
               enode = part_nodes(iAss+1,iPart+2) 
           !    ndID = part_nodes(iAss+1,iPart+2) - part_nodes(iAss+1,iPart+1)
               call mw_sumup(iAss, iPart, fstrSOLID%QFORCE(snode*ndof+1:enode*ndof), ndof)
           !    call mw_send_recv_r(fstrSOLID%QFORCE(snode*ndof+1:enode*ndof), ndID, ndof, ik)
            enddo
         enddo
  enddo

end subroutine fstr_UpdateNewton

!> Update elastiplastic status
subroutine fstr_UpdateEPState( fstrSOLID)
  use m_fstr
  use m_static_lib
  use m_ElastoPlastic
  
  include "HEC_MW3_For.h"
  
  type (fstr_solid)          :: fstrSOLID   !< fstr_solid
  
  integer(kind=kint) :: ic_type, icel, ngauss, i
  integer :: iAss, iPart, iElem, iNode, iGrp
 
  icel = 0
  do iAss = 0, mw_get_num_of_assemble_model()-1
    call mw_select_assemble_model( iAss )
    do iPart = 0, mw_get_num_of_mesh_part()-1
      call mw_select_mesh_part( iPart )
      do iElem = 0, mw_get_num_of_element()-1
        icel = icel+1
        call mw_select_element( iElem )
        ic_type = mw_get_element_type()
        ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
        ngauss = NumOfQuadPoints( ic_type )
        if( .not. isElastoplastic( fstrSOLID%elements(icel)%gausses(1)   &
         %pMaterial%mtype ) ) cycle
        do i=1,ngauss
          call updateEPState( fstrSOLID%elements(icel)%gausses(i) )
        enddo
      enddo
    enddo
  enddo

end subroutine

end module m_fstr_Update
