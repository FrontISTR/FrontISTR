!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by X. YUAN(AdavanceSoft)                          !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions to take into acount external load
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2011/01/16
!>  \version    0.00
!!
!======================================================================!
module m_fstr_ass_load
    use m_fstr
    use m_static_lib
    implicit none
	
	private :: fstr_call_dl
	
    contains
!
!======================================================================!
!> This subroutine assmble following external force into fstrSOLID%GL and hecMAT%B afterwards
!>  -#  concentrated nodal force
!>  -#  surface pressure
!>  -#  volume force
!>  -#  thermal force
!
    subroutine fstr_ass_load(cstep, fstrSOLID,factor)
!======================================================================!
      use mMechGauss
      use mULoad
	  
      include "HEC_MW3_For.h"
	  
      integer, intent(in)                  :: cstep     !< current step
      type (fstr_solid),intent(inout)      :: fstrSOLID !< fstr_solid
      real(kind=kreal),intent(in)          :: factor    !< loading factor

      real(kind=kreal) :: x,y,z, xx(20), yy(20), zz(20)
      real(kind=kreal) :: params(0:6)
      real(kind=kreal) :: vect(60)
      integer(kind=kint) :: iwk(60)
      integer(kind=kint) :: nodLocal(20), nids(0:20)
      real(kind=kreal) :: tt(20), tt0(20)
      integer(kind=kint) :: ndof, ig0, ig, ityp, ltype, ik, i, j, iii
      integer(kind=kint) :: icel, ic_type, nn, iS, isect, iset, nsize
      integer(kind=kint) :: itype, iE, iErr, grpid, npart, pid, snode
      real(kind=kreal) :: fval, rho, pa1
      logical :: fg_surf
      integer(kind=kint) :: tstep, name_len
      character(len=HECMW_NAME_LEN) :: header_name
      type( tMaterial ), pointer :: material     !< material information
	  
      integer :: iAss, iPart, iElem, iNode, iGrp, ndID, cid, mwigrp

      ndof = assDOF(1)
! -------------------------------------------------------------------
!  CLOAD
! -------------------------------------------------------------------
      fstrSOLID%GL(:)=0.d0

      if( associated(fstrSOLID%cload_grp) ) then  
        do iAss = 0, mw_get_num_of_assemble_model()-1
           call mw_select_assemble_model( iAss )
           npart = mw_get_num_of_mesh_part()
		   
           do igrp= 1, size(fstrSOLID%cload_grp)
	         pid = fstrSOLID%cload_grp(igrp)%part_id
             if( pid>npart-1 ) cycle
             grpid= fstrSOLID%cload_grp(igrp)%gid
             if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
             read( fstrSOLID%cload_grp(igrp)%grp_name, *, IOSTAT=iErr ) ndID
             if( iErr==0 ) then
             !    ndID=part_nodes(iAss+1,pid+1)+ndID+1
                 if( all(global_node_ID/=ndID) ) cycle
                 ndID = mw_get_node_index( ndID )+1
                 do ik=1,assDOF(1)
                   if( fstrSOLID%cload_grp(igrp)%dof(ik)==1 ) then
                     fstrSOLID%GL(ndof*(ndID-1)+ik)=fstrSOLID%GL(ndof*(ndID-1)+ik)  &
                       +fstrSOLID%cload_grp(igrp)%fval
                   endif
                 enddo
             else
               call mw_select_mesh_part_with_id(pid)
               snode = part_nodes(iAss+1, pid+1)
			   do mwigrp = 0, mw_get_num_of_boundary_bnode_mesh()-1
                 nn = mw_get_bnode_mesh_namelength( mwigrp )
                 header_name = ''
                 call mw_get_bnode_mesh_name(mwigrp, header_name, nn)
                 if( header_name(1:nn)/=trim(fstrSOLID%cload_grp(igrp)%grp_name) ) cycle
                 do  iNode=0, mw_get_num_of_bnode_in_bnode_mesh(mwigrp)-1
			        ndID= mw_get_node_id_in_bnode_mesh(mwigrp, iNode)
                !    ndID= part_nodes(iAss+1,pid+1)+ndID+1
                    if( all(global_node_ID/=ndID) ) cycle
                    ndID = mw_get_node_index( ndID )+1+snode
                    do ik=1,assDOF(1)
                      if( fstrSOLID%cload_grp(igrp)%dof(ik)==1 ) then
                        fstrSOLID%GL(ndof*(ndID-1)+ik)=fstrSOLID%GL(ndof*(ndID-1)+ik)  &
                          +fstrSOLID%cload_grp(igrp)%fval
                      endif
                    enddo
                 enddo
               enddo
             endif		   
           enddo	   
        enddo
      endif
	  

! -------------------------------------------------------------------
!  DLOAD
! -------------------------------------------------------------------

      if( associated(fstrSOLID%dload_grp) ) then  
        do iAss = 0, mw_get_num_of_assemble_model()-1
           call mw_select_assemble_model( iAss )
           npart = mw_get_num_of_mesh_part()
		   
           do igrp= 1, size(fstrSOLID%dload_grp)
	         pid = fstrSOLID%dload_grp(igrp)%part_id
             if( pid>npart-1 ) cycle
             grpid= fstrSOLID%dload_grp(igrp)%gid
             if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
             ltype = fstrSOLID%dload_grp(igrp)%itype
             params(0:6) = fstrSOLID%dload_grp(igrp)%fval
             call mw_select_mesh_part( pid )
             
             read( fstrSOLID%dload_grp(igrp)%grp_name, *, IOSTAT=iErr ) ndID
             if( iErr==0 ) then   
             !    icel=part_elems(iAss+1,pid+1)+ndID+1
                 if( all(global_elem_ID/=ndID) ) cycle
                 icel = mw_get_element_index( ndID )+1
                 call mw_select_element_with_id( ndID )
                 call mw_get_element_vert_node_id( nids )
                 nn = mw_get_num_of_element_vert()
                 ic_type = mw_get_element_type()
                 ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
                 do iNode = 1, nn
                   cid = nids(iNode-1)
                   call mw_get_node_coord( cid, x,y,z )
                   xx(iNode)=x
                   yy(iNode)=y
                   zz(iNode)=z
                   cid= part_nodes(iAss+1,pid+1)+mw_get_node_index(cid)
                   do i=1,ndof                      
                      iwk(ndof*(iNode-1)+i)=ndof*cid+i
                   enddo
                 enddo
                 material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
                 rho = material%variables(M_DENSITY)
                 pa1 = material%variables(M_THICK)
                 if (assDOF(1)==2) pa1=1.d0
                 call fstr_call_dl( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn),   &
                   rho, pa1, ltype, params, vect, nsize, fstrSOLID%elements(icel)%iset ) 
                 do j=1,nsize
                      fstrSOLID%GL( iwk(j) )=fstrSOLID%GL( iwk(j) )+vect(j)
                 enddo 
             else
               fg_surf = (fstrSOLID%dload_grp(igrp)%itype == 100)
               if( fg_surf ) then
			     do mwigrp = 0, mw_get_num_of_boundary_bface_mesh()-1
                   nn = mw_get_bface_mesh_namelength( mwigrp )
                   header_name=''
                   call mw_get_bface_mesh_name(mwigrp, header_name, nn)
                   if( header_name(1:nn)/=trim(fstrSOLID%dload_grp(igrp)%grp_name) ) cycle
                   do  iElem=0, mw_get_num_of_bface( mwigrp )-1
			           ndID= mw_get_elem_id_bface(mwigrp, iElem)
                       call mw_select_element_with_id( ndID )
                      ! icel= part_elems(iAss+1,pid+1)+ndID+1  
                       icel = mw_get_element_index( ndID )+1		  
                       call  mw_get_element_vert_node_id( nids )
                       nn = mw_get_num_of_element_vert()
                       ic_type = mw_get_element_type()
                       ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
                       do iNode = 1, nn
                         cid = nids(iNode-1)
                         call mw_get_node_coord( cid, x,y,z )
                         xx(iNode)=x
                         yy(iNode)=y
                         zz(iNode)=z    
                         cid= part_nodes(iAss+1,pid+1)+mw_get_node_index(cid)
                         do i=1,ndof
                           iwk(ndof*(iNode-1)+i)=ndof*cid+i
                         enddo
                       enddo
                       material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
                       rho = material%variables(M_DENSITY)
                       pa1 = material%variables(M_THICK)
                       if (ndof==2) pa1=1.d0
                       ltype = mw_get_face_id_bface(mwigrp, iElem)*10
                       call fstr_call_dl( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn),   &
                        rho, pa1, ltype, params, vect, nsize, fstrSOLID%elements(icel)%iset ) 
                       do j=1,nsize
                        fstrSOLID%GL( iwk(j) )=fstrSOLID%GL( iwk(j) )+vect(j)
                       enddo 
                   enddo
                 enddo
               else
                 do mwigrp = 0, mw_get_num_of_boundary_bvolume_mesh()-1
                   nn = mw_get_bvolume_mesh_namelength( mwigrp )
                   call mw_get_bvolume_mesh_name(mwigrp, header_name, nn)
                   if( header_name(1:nn)/=trim(fstrSOLID%dload_grp(igrp)%grp_name) ) cycle
                   do  iElem=0, mw_get_num_of_bvolume( mwigrp )-1
			         ndID= mw_get_elem_id_bvolume(mwigrp, iElem)
                     call mw_select_element_with_id( ndID )
                    ! icel= part_elems(iAss+1,pid+1)+ndID+1   
                     icel = mw_get_element_index( ndID )+1		
                     call  mw_get_element_vert_node_id( nids )
                     nn = mw_get_num_of_element_vert()
                     ic_type = mw_get_element_type()
                     ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
                     do iNode = 1, nn
                       cid = nids(iNode-1)
                       call mw_get_node_coord( cid, x,y,z )
                       xx(iNode)=x
                       yy(iNode)=y
                       zz(iNode)=z    
                       cid= part_nodes(iAss+1,pid+1)+mw_get_node_index(cid)
                       do i=1,ndof
                         iwk(ndof*(iNode-1)+i)=ndof*cid+i
                       enddo
                     enddo
                     material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
                     rho = material%variables(M_DENSITY)
                     pa1 = material%variables(M_THICK)
                     if (ndof==2) pa1=1.d0
                     call fstr_call_dl( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn),   &
                      rho, pa1, ltype, params, vect, nsize, fstrSOLID%elements(icel)%iset ) 
                     do j=1,nsize
                       fstrSOLID%GL( iwk(j) )=fstrSOLID%GL( iwk(j) )+vect(j)
                     enddo 
                   enddo
                 enddo
                 do mwigrp = 0, mw_get_num_of_elementgroup()-1
                   nn = mw_get_elementgroup_name_length( mwigrp )
                   header_name = ''
                   call mw_get_elementgroup_name(mwigrp, header_name, nn)
                   if( header_name(1:nn)/=trim(fstrSOLID%dload_grp(igrp)%grp_name) ) cycle
                   do  iElem=0, mw_get_num_of_element_id( mwigrp )-1
			         ndID= mw_get_element_id_with_elementgroup(mwigrp, iElem)
                     call mw_select_element_with_id( ndID ) 
                     icel = mw_get_element_index( ndID )+1		
                     call  mw_get_element_vert_node_id( nids )
                     nn = mw_get_num_of_element_vert()
                     ic_type = mw_get_element_type()
                     ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)
                     do iNode = 1, nn
                       cid = nids(iNode-1)
                       call mw_get_node_coord( cid, x,y,z )
                       xx(iNode)=x
                       yy(iNode)=y
                       zz(iNode)=z
                       cid= part_nodes(iAss+1,pid+1)+mw_get_node_index(cid)
                       do i=1,ndof
                         iwk(ndof*(iNode-1)+i)=ndof*cid+i
                       enddo
                     enddo
                     material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
                     rho = material%variables(M_DENSITY)
                     pa1 = material%variables(M_THICK)
                     if (ndof==2) pa1=1.d0
                     call fstr_call_dl( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn),   &
                      rho, pa1, ltype, params, vect, nsize, fstrSOLID%elements(icel)%iset )  
                     do j=1,nsize
                       fstrSOLID%GL( iwk(j) )=fstrSOLID%GL( iwk(j) )+vect(j)
                     enddo 
                   enddo
                 enddo
               endif
               
             endif		   
           enddo	   
        enddo
      endif
 	  
! -----Uload
      call uloading( cstep, factor, fstrSOLID%GL )
	  
!C
!C Update for fstrSOLID%GL
!C 
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part_with_id( iPart )
            do iNode = 0, mw_get_num_of_neibpe(iPart)-1
               ik = mw_get_transrank(iPart, iNode)
               ndID = part_nodes(iAss+1,iPart+2) - part_nodes(iAss+1,iPart+1)
     !          call mw_send_recv_r(fstrSOLID%GL(part_nodes(iAss+1,iPart+1)+1:part_nodes(iAss+1,iPart+2)), ndID, ndof, ik)
            enddo
         enddo
      enddo

      iii = -1	  
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part_with_id( iPart )
            do iNode = 0, mw_get_num_of_node()-1
            !  ndID = part_nodes(iAss+1,iPart+1)+iNode
              iii = iii+1
              ndID = mw_get_node_id( iNode )
              do ik = 1, ndof
                if( fstrSOLID%GL(ndof*iii+ik)/=0.d0 ) then
                  iErr= mw_nl_rhs_set_bc(iPart, ndID, ik-1, factor*fstrSOLID%GL(ndof*iii+ik))
                endif
              enddo
            enddo
         enddo
      enddo
!
! -------------------------------------------------------------------
!  TLOAD : THERMAL LOAD USING TEMPERATURE
! -------------------------------------------------------------------
!C
!C Set Temperature
!C
      if( fstrSOLID%TEMP_irres == 1 ) then
          open( unit=IFTEMP, file=fstrSOLID%temp_grp(1)%fname, status="old", action='read', iostat=iErr )
          if( iErr>0 ) then
            WRITE(*,*) ' Cannot find the file that provide temprature information'
            WRITE(IMSG,*) ' Cannot find the file that provide temprature information'
            call hecmw_abort( hecmw_comm_get_comm())
          endif
          read( IFTEMP, iostat=iErr ) ik
          if( iErr>0 ) then
            WRITE(*,*) ' Fail in reading temperature file'
            WRITE(IMSG,*) ' Fail in reading temperature file'
            call hecmw_abort( hecmw_comm_get_comm())
          endif
          do i=1,ik
            read( IFTEMP, iostat=iErr ), ig, pa1
            if( iErr>0 ) then
              WRITE(*,*) ' Fail in reading temperature file'
              WRITE(IMSG,*) ' Fail in reading temperature file'
              call hecmw_abort( hecmw_comm_get_comm())
            endif
            if( ig/=fstrSOLID%TEMP_tstep ) cycle
            do j=1,total_node
              read( IFTEMP, iostat=iErr ) ndID, fstrSOLID%temperature(j)
              if( iErr>0 ) then
                WRITE(*,*) ' Fail in reading temperature file'
                WRITE(IMSG,*) ' Fail in reading temperature file'
                call hecmw_abort( hecmw_comm_get_comm())
              endif
            enddo
          enddo
          close( IFTEMP )
		  
      else if( associated(fstrSOLID%temp_grp) ) then  
        do iAss = 0, mw_get_num_of_assemble_model()-1
           call mw_select_assemble_model( iAss )
           npart = mw_get_num_of_mesh_part()
		   
           do igrp= 1, size(fstrSOLID%temp_grp)
	         pid = fstrSOLID%temp_grp(igrp)%part_id
             if( pid>npart-1 ) cycle
             grpid= fstrSOLID%temp_grp(igrp)%gid
             if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
             read( fstrSOLID%temp_grp(igrp)%grp_name, *, IOSTAT=iErr ) ndID
             if( iErr==0 ) then
                 if( all(global_node_ID/=ndID) ) cycle
                 ndID = part_nodes(iAss+1,pid+1)+mw_get_node_index( ndID )+1
                 fstrSOLID%temperature(ndID)= fstrSOLID%temp_grp(igrp)%fval*factor
             else
               call mw_select_mesh_part_with_id(pid)
			   do mwigrp = 0, mw_get_num_of_boundary_bnode_mesh()-1              
                 nn = mw_get_bnode_mesh_namelength( mwigrp )
                 call mw_get_bnode_mesh_name(mwigrp, header_name, nn)
                 if( header_name(1:nn)/=trim(fstrSOLID%temp_grp(igrp)%grp_name) ) cycle
                 do  iNode=0, mw_get_num_of_bnode_in_bnode_mesh(mwigrp)-1
			        ndID= mw_get_node_id_in_bnode_mesh(mwigrp, iNode)
                    ndID = part_nodes(iAss+1,pid+1)+mw_get_node_index( ndID )+1
                    fstrSOLID%temperature(ndID)= fstrSOLID%temp_grp(igrp)%fval*factor
                 enddo
               enddo
             endif		   
           enddo	   
        enddo
      endif

      if( fstrSOLID%TEMP_irres == 1 .or. associated(fstrSOLID%temp_grp) ) then
    !    fstrSOLID%reftemp(:) = fstrSOLID%temperature(:)
		
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
                xx(iNode)=x
                yy(iNode)=y
                zz(iNode)=z
                cid= part_nodes(iAss+1,iPart+1)+mw_get_node_index(cid)+1
                tt0(iNode)=fstrSOLID%reftemp( cid ) 
                tt(iNode) = fstrSOLID%temperature(cid ) 
              enddo

              pa1 = material%variables(M_THICK)
              if (ndof==2) pa1=1.d0
              if (ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
                call TLOAD_C2(ic_type,nn,xx(1:nn),yy(1:nn),tt(1:nn),tt0(1:nn),  &
			         fstrSOLID%elements(icel)%gausses,pa1,fstrSOLID%elements(icel)%iset,vect(1:nn*2) )

              else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.  &
                     ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
                call TLOAD_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn), &
                            fstrSOLID%elements(icel)%gausses,vect(1:nn*ndof) )

              else if ( ic_type==741 .or. ic_type==731 ) then
                if( myrank == 0 ) then
                  WRITE(IMSG,*) '*------------------------', &
                                 '-------------------*'
                  WRITE(IMSG,*) ' Thermal loading option for shell elements', &
                                 'not yet available.'
                  WRITE(IMSG,*) '*------------------------', &
                               '-------------------*'
                  call hecmw_abort( hecmw_comm_get_comm())
                endif
              endif

              ndID = 0			  
              do iNode = 1, nn
                cid = nids(iNode-1)
                do i=1,ndof
                  ndID  = ndID+1
                  iErr= mw_nl_rhs_add_bc(iPart, cid, i-1, vect(ndID) )
                enddo
              enddo
            enddo
          enddo
        enddo
      endif


    end subroutine fstr_ass_load
	
	subroutine fstr_call_dl( ic_type, nn, xx, yy, zz, rho, pa1, ltype, params, vect, nsize, iset )
       integer, intent(in)          :: ic_type   !< element type
       integer, intent(in)          :: nn        !< number of elemental nodes
       real(kind=kreal), intent(in) :: xx(nn), yy(nn), zz(nn) !< elemental nodes coord
       real(kind=kreal), intent(in) :: rho, pa1  !< material properties
       integer, intent(in)          :: ltype     !< load type
       real(kind=kreal), intent(in) :: params(0:6)
       real(kind=kreal), intent(out):: vect(:)
       integer, intent(out)         :: nsize
       integer, optional            :: iset
       
	      if (ic_type==241 .or.ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
            call DL_C2(ic_type,nn,xx(:),yy(:),rho,pa1,ltype,params,vect(:),nsize,iset)
          else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.     &
                   ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
            call DL_C3(ic_type,nn,xx(:),yy(:),zz(:),rho,ltype,params,vect(:),nsize)
          else if ( ic_type.eq.741 ) then
            call DL_S4( xx,yy,zz,rho,pa1,ltype,params,vect,nsize )
          else if ( ic_type.EQ.731 ) then
            call DL_S3( xx,yy,zz,rho,pa1,ltype,params,vect,nsize )
          endif
    end subroutine

end module m_fstr_ass_load
