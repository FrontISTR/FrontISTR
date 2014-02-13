!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by  X. YUAN(AdavanceSoft)    !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions to visualize result.
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2013/02/10
!>  \version    0.00
!!
!======================================================================!
module m_fstr_visualize
  use m_fstr_NodalStress
  implicit none

  contains

!> Output result
!----------------------------------------------------------------------*
  subroutine fstr_visualize( cstep, totstep, fstrSOLID, ttime )
!----------------------------------------------------------------------*
    use m_fstr
	
	include "HEC_MW3_For.h"

    integer, intent(in)                   :: cstep       !< current step number
    integer, intent(in)                   :: totstep     !< total steps
    type (fstr_solid), intent(in)         :: fstrSOLID
    real(kind=kreal),intent(in)           :: ttime

    integer(kind=kint) :: iAss, iPart, iDof, nDof
    integer(kind=kint) :: snode, enode, nnode
    real(kind=kreal)   :: ndvalue(total_node*6)
    character(len=10)  :: label, cUnit

    iAss = mw_get_num_of_assemble_model()-1
    call mw_select_assemble_model( iAss )
    if( gVisType==1 ) then
       do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
			! displacement
            label = "disp" //char(0)
            cUnit = "unkonwn"//char(0)
            ndof = 3
            call mw_rec_vtk_label(iPart, label, cUnit, ndof)
            snode = part_nodes(iAss+1,iPart+1)+1
            enode = part_nodes(iAss+1,iPart+2)
            nnode = enode-snode+1
            call mw_rec_vtk_variable(iPart, nnode, label, fstrSOLID%unode((snode-1)*ndof+1:enode*ndof))
            ! strain
            label = "strain" //char(0)
            cUnit = "none"//char(0)
            ndof = 6
            call mw_rec_vtk_label(iPart, label, cUnit, ndof)
            call mw_rec_vtk_variable(iPart, nnode, label, fstrSOLID%STRAIN)
			! stress
            label = "stress" //char(0)
            cUnit = "unkonwn"//char(0)
            ndof = 7
            call mw_rec_vtk_label(iPart, label, cUnit, ndof)
            call mw_rec_vtk_variable(iPart, nnode, label, fstrSOLID%STRESS)
       enddo
       call mw_print_vtk_fem()
    else if( gVisType==2 ) then
       do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            label = "disp" //char(0)
            cUnit = "unkonwn"//char(0)
            ndof = 3
            call mw_rec_avs_label(iPart, label, cUnit, ndof)
            snode = part_nodes(iAss+1,iPart+1)+1
            enode = part_nodes(iAss+1,iPart+2)
            nnode = enode-snode+1
            call mw_rec_avs_variable(iPart, nnode, label, fstrSOLID%unode((snode-1)*ndof+1:enode*ndof))
            ! strain
            label = "strain" //char(0)
            cUnit = "none"//char(0)
            ndof = 6
            call mw_rec_avs_label(iPart, label, cUnit, ndof)
            call mw_rec_avs_variable(iPart, nnode, label, fstrSOLID%STRAIN)
            ! stress
            label = "stress" //char(0)
            cUnit = "unkonwn"//char(0)
            ndof = 7
            call mw_rec_vtk_label(iPart, label, cUnit, ndof)
            call mw_rec_vtk_variable(iPart, nnode, label, fstrSOLID%STRESS)
       enddo
       call mw_print_avs_fem()
    else if( gVisType==3 ) then
       do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            label = "disp" //char(0)
            cUnit = "unkonwn"//char(0)
            ndof = 3
            call mw_rec_uns_label(iPart, label, cUnit, ndof)
            snode = part_nodes(iAss+1,iPart+1)+1
            enode = part_nodes(iAss+1,iPart+2)
            nnode = enode-snode+1
            call mw_rec_uns_variable(iPart, nnode, label, fstrSOLID%unode((snode-1)*ndof+1:enode*ndof))
            ! strain
            label = "strain" //char(0)
            cUnit = "none"//char(0)
            ndof = 6
            call mw_rec_uns_label(iPart, label, cUnit, ndof)
            call mw_rec_uns_variable(iPart, nnode, label, fstrSOLID%STRAIN)
            ! stress
            label = "stress" //char(0)
            cUnit = "unkonwn"//char(0)
            ndof = 7
            call mw_rec_vtk_label(iPart, label, cUnit, ndof)
            call mw_rec_vtk_variable(iPart, nnode, label, fstrSOLID%STRESS)
       enddo
       call mw_print_uns_fem()
    endif
      
  end subroutine fstr_visualize

end module m_fstr_visualize
