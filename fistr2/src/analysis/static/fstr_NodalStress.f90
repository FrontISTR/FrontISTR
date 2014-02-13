!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by X. YUAN(AdavanceSoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions to caluclation nodal stress
!!
!>  \author     X. YUAN(AdavanceSoft)
!>  \date       2010/02/01
!>  \version    0.00
!!
!======================================================================!
module m_fstr_NodalStress
  implicit none

  contains

!> Calculate NODAL STRESS of solid elements
!----------------------------------------------------------------------*
  subroutine fstr_NodalStress3D(fstrSOLID, tdstrain, tdstress  )
!----------------------------------------------------------------------*
    use m_fstr
    use m_static_lib
    type (fstr_solid)      :: fstrSOLID                    !< fstr_solid
    real(kind=kreal)       :: tdstrain(total_node*6)   !< nodal strain
    real(kind=kreal)       :: tdstress(total_node*7)   !< nodal stress
	
    include "HEC_MW3_For.h"
	
!** Local variables
    integer(kind=kint) :: i, ic_type, icel, ID_area, ie, ielem, cid
    integer(kind=kint) :: ii, ik , iS, iS0, iE0, itype, j, jj, jS, k, nn
    REAL(kind=kreal) xx(20),yy(20),zz(20)
    REAL(kind=kreal) edisp(60),estrain(6),estress(6)
    REAL(kind=kreal) tt(20),tt0(20)
    integer(kind=kint) nodLOCAL(20)
    REAL(kind=kreal) edstrain(20,6),edstress(20,6)
    integer(kind=kint) :: arrayTotal
    real(kind=kreal) s11, s22, s33, s12, s23, s13, ps, smises
    real(kind=kreal) ddunode(3,20)
	
	integer(kind=kint) :: nids(0:20)
    integer :: iAss, iPart, iNode, iGrp, iErr, npart, snode, enode

    integer(kind=kint) :: nnumber(total_node)
    real(kind=kreal)   :: ndstrain(total_node,6), ndstress(total_node,7)
    real(kind=kreal)   :: temp(total_node)


!*ZERO CLEAR
    ndstrain=0.d0; ndstress=0.d0
    arrayTotal=0
    nnumber=0

	iAss = mw_get_num_of_assemble_model()-1
    nPart = mw_get_num_of_mesh_part()-1
    icel = part_elems(iAss+1, nPart+1)
    snode = part_nodes(iAss+1, nPart+1)
    enode = part_nodes(iAss+1, nPart+2)

     call mw_select_assemble_model( iAss )
     do iPart = 0, nPart
        call mw_select_mesh_part( iPart )
        call mw_matrix_clear( iPart )
        do iElem = 0, mw_get_num_of_element()-1
          icel = icel+1
          call mw_select_element( iElem )
          call  mw_get_element_vert_node_index( nids )
          nn = mw_get_num_of_element_vert()
          ic_type = mw_get_element_type()
          ic_type = mw_mw3_elemtype_to_fistr_elemtype(ic_type)

!  --- calculate the nodal stress
          call NodalStress_C3( ic_type, nn,                                       &
                             fstrSOLID%elements(icel)%gausses,                  &
                             edstrain(1:nn,:), edstress(1:nn,:)      )


          do j=1,nn
            cid= nids(j-1)+1+snode
            ndstrain( cid,: ) = ndstrain(cid,:) + edstrain(j,:)
            ndstress( cid,1:6 ) = ndstress(cid,1:6) + edstress(j,:)
            nnumber( cid )=nnumber( cid )+1
          enddo

          call ElementStress_C3( ic_type,                                       &
                                 fstrSOLID%elements(icel)%gausses,              &
                                 estrain, estress     )
          fstrSOLID%ESTRAIN(6*icel-5:6*icel) = estrain
          fstrSOLID%ESTRESS(7*icel-6:7*icel-1) = estress
          s11=fstrSOLID%ESTRESS(7*icel-6)
          s22=fstrSOLID%ESTRESS(7*icel-5)
          s33=fstrSOLID%ESTRESS(7*icel-4)
          s12=fstrSOLID%ESTRESS(7*icel-3)
          s23=fstrSOLID%ESTRESS(7*icel-2)
          s13=fstrSOLID%ESTRESS(7*icel-1)
          ps=(s11+s22+s33)/3.0
          smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                +s12**2+s23**2+s13**2
          fstrSOLID%ESTRESS(7*icel)=sqrt(3.0*smises)
        enddo

      enddo       ! icel roop.



!** Average over nodes
    do i=snode+1,enode
      if( nnumber(i)<=0 ) cycle
      ndstrain(i,:)=ndstrain(i,:)/nnumber(i)
      ndstress(i,1:6)=ndstress(i,1:6)/nnumber(i)
    enddo
!** CALCULATE von MISES stress
    do i=1,snode+1,enode
      s11=ndstress(i,1)
      s22=ndstress(i,2)
      s33=ndstress(i,3)
      s12=ndstress(i,4)
      s23=ndstress(i,5)
      s13=ndstress(i,6)
      ps=(s11+s22+s33)/3.0
      smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
            +s12**2+s23**2+s13**2
      ndstress(i,7)=sqrt(3.d0*smises)
    enddo
!** Set Array 
    do i=snode+1,enode
      do j=1,6
  !      fstrSOLID%STRAIN(6*(i-1)+j)=ndstrain(i,j)
		tdSTRAIN(6*(i-1)+j)=ndstrain(i,j)
      enddo
      do j=1,7
   !     fstrSOLID%STRESS(7*(i-1)+j)=ndstress(i,j)
		tdSTRESS(7*(i-1)+j)=ndstress(i,j)
      enddo
    enddo

  end subroutine fstr_NodalStress3D


end module m_fstr_NodalStress
