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
!> \brief  This module provides function to calcualte residual of nodal force.
!!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \date       2009/09/14
!>  \version    0.00
!!
!======================================================================!
module m_fstr_Residual
  use hecmw
  implicit none
  
  contains

!C---------------------------------------------------------------------*
      subroutine fstr_Update_NDForce(cstep,fstrSOLID,factor )
!C---------------------------------------------------------------------*
!> In this subroutine, nodal force arose from prescribed displacement constarints 
!> are cleared and nodal force residual is calculated. 
!> Those constraints considered here includes:
!!-#  nodal displacement
      use m_fstr
      use mULoad
	  
      include "HEC_MW3_For.h"
	  
      integer(kind=kint), intent(in)       :: cstep      !< current step
      type (fstr_solid), intent(inout)     :: fstrSOLID  !< we need boundary conditions of curr step
      real(kind=kreal), intent(in)         :: factor     !< loading factor
!    Local variables
      integer(kind=kint) ndof,ityp,ik,idof, iErr, nn, npart
      integer(kind=kint) :: grpid  
      real(kind=kreal) :: rhs, lambda
      integer(kind=kint) :: iAss, igrp, iPart, iNode, ndID, pid
      character(len=HECMW_NAME_LEN) :: header_name
	  
      ndof = assDOF(1)

!    Set residual load
      idof = 0
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            do iNode = 0, mw_get_num_of_node()-1
              ndID = mw_get_node_id(iNode)
              do ik = 1, ndof
                idof = idof+1
                rhs = factor*fstrSOLID%GL(idof)-fstrSOLID%QFORCE(idof)
                iErr= mw_nl_rhs_set_bc(iPart, ndID, ik-1, rhs)
              enddo
            enddo
         enddo
      enddo
	  
!    Consider Uload
      call uResidual( cstep, factor )


!    Consider SPC condition
      if( associated(fstrSOLID%boundary_grp) ) then  
        do iAss = 0, mw_get_num_of_assemble_model()-1
           call mw_select_assemble_model( iAss )
           npart = mw_get_num_of_mesh_part()
		   
           do igrp= 1, size(fstrSOLID%boundary_grp)
	         pid = fstrSOLID%boundary_grp(igrp)%part_id
             grpid= fstrSOLID%boundary_grp(igrp)%gid
             if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
             read( fstrSOLID%boundary_grp(igrp)%grp_name, * , IOSTAT=iErr ) ndID
             if( iErr==0 ) then 
                 do ik=1,assDOF(1)
                   if( fstrSOLID%boundary_grp(igrp)%dof(ik)==1 ) then
                     iErr= mw_nl_rhs_set_bc(pid, ndID, ik-1, 0.d0)
                   endif
                 enddo
             else
               if( pid>npart-1 ) cycle
               nn = mw_get_bnode_mesh_namelength( pid )
               call mw_get_bnode_mesh_name(pid, header_name, nn)
               if( trim(header_name)/=trim(fstrSOLID%boundary_grp(igrp)%grp_name) ) cycle
               do  iNode=0, mw_get_num_of_bnode_in_bnode_mesh(pid)-1
			        ndID= mw_get_node_id_in_bnode_mesh(pid, iNode)
                    do ik=1,6
                      if( fstrSOLID%boundary_grp(igrp)%dof(ik)==1 ) then
                        iErr= mw_nl_rhs_set_bc(pid, ndID, ik-1, 0.d0)
                      endif
                    enddo
               enddo
             endif		   
           enddo
        enddo
      endif
	   	  
    !  call mw_send_recv()

      end subroutine fstr_Update_NDForce
	  
!> Calculate magnitude of a real vector
      real(kind=kreal) function fstr_get_residual(fstrSOLID)
      use m_fstr
      type (fstr_solid), intent(inout)     :: fstrSOLID
      include "HEC_MW3_For.h"
      integer :: i, iAss, igrp, iPart, iNode, ndID, pid
      integer :: ndof, snode, enode
	  
	  ndof = assDOF(1)

      fstr_get_residual =0.d0
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part_with_id( iPart )
            snode = part_nodes(iAss+1,iPart+1)
            enode = part_nodes(iAss+1,iPart+2)
            call mw_get_rhs_vector(fstrSOLID%ddunode(snode*ndof+1:  &
               enode*ndof), iPart)
            do i=snode*ndof+1, enode*ndof
	          fstr_get_residual = fstr_get_residual + fstrSOLID%ddunode(i)*fstrSOLID%ddunode(i)
            enddo			  
         enddo
      enddo
      if( myrank==0) then
         write(IMSG,*) '####fstrNLGEOM_SetResidual finished'
      end if
      end function

end module m_fstr_Residual
