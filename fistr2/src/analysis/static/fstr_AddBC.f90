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
!> \brief  This module provides a function to deal with prescribed displacement.
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2009/08/31
!>  \version    0.00
!!
!======================================================================!
module m_fstr_AddBC

   implicit none

   contains
   
!>  Add Essential Boundary Conditions
!----------------------------------------------------------------------*
      subroutine fstr_AddBC(cstep,fstrSOLID,iter)
!----------------------------------------------------------------------*
      use m_fstr
      include "HEC_MW3_For.h"
      integer, intent(in)       :: cstep     !< current step
      type (fstr_solid        ) :: fstrSOLID !< fstr_solid
      integer(kind=kint)        :: iter      !< NR iterations

      integer(kind=kint) :: iAss, igrp, nPart, iErr, ik
      real(kind=kreal) :: RHS,factor 
      integer(kind=kint) :: grpid, ndID, pid, nn, iNode, mwigrp
      character(len=HECMW_NAME_LEN) :: header_name

      factor = fstrSOLID%FACTOR(2)-fstrSOLID%FACTOR(1)

      if( associated(fstrSOLID%boundary_grp) ) then  
        do iAss = 0, mw_get_num_of_assemble_model()-1
           call mw_select_assemble_model( iAss )
           npart = mw_get_num_of_mesh_part()
		   
           do igrp= 1, size(fstrSOLID%boundary_grp)
	         pid = fstrSOLID%boundary_grp(igrp)%part_id
             if( pid>npart-1 ) cycle
             grpid= fstrSOLID%boundary_grp(igrp)%gid
             if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
             RHS= fstrSOLID%boundary_grp(igrp)%fval*factor
             if( iter>1 ) RHS=0.d0
             read( fstrSOLID%boundary_grp(igrp)%grp_name, * , IOSTAT=iErr ) ndID
             if( iErr==0 ) then
                 if( all(global_node_ID/=ndID) ) cycle
                 do ik=1,assDOF(1)
                   if( fstrSOLID%boundary_grp(igrp)%dof(ik)==1 ) then
                     iErr= mw_matrix_rhs_set_bc2(pid, ndID, ik-1, 1.D0, RHS)
                   endif
                 enddo
             else
               call mw_select_mesh_part_with_id(pid)
			   do mwigrp = 0, mw_get_num_of_boundary_bnode_mesh()-1
                 nn = mw_get_bnode_mesh_namelength( mwigrp )
                 header_name=''
                 call mw_get_bnode_mesh_name(mwigrp, header_name(1:nn), nn)
                 if( header_name(1:nn)/=trim(fstrSOLID%boundary_grp(igrp)%grp_name) ) cycle
                 do  iNode=0, mw_get_num_of_bnode_in_bnode_mesh(mwigrp)-1
			        ndID= mw_get_node_id_in_bnode_mesh(mwigrp, iNode)
                    if( all(global_node_ID/=ndID) ) cycle
                    do ik=1,6
                      if( fstrSOLID%boundary_grp(igrp)%dof(ik)==1 ) then
                        iErr= mw_matrix_rhs_set_bc2(pid, ndID, ik-1, 1.D0, RHS)
                      endif
                    enddo
                 enddo
               enddo
             endif		   
           enddo
        enddo
      endif

      if( myrank==0 ) then
         write(IMSG,*) '####fstr_AddBC finished'
      end if

      end subroutine fstr_AddBC

end module m_fstr_AddBC
