!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
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
!>  \date       2009/09/09
!>  \version    0.00
!!
!======================================================================!
module m_fstr_NodalStress
  implicit none

  contains

!> Calculate NODAL STRESS of solid elements
!----------------------------------------------------------------------*
  subroutine fstr_NodalStress3D( hecMESH, hecMAT, fstrSOLID,         &
                                 tdstrain, tdstress  )
!----------------------------------------------------------------------*
    use m_fstr
    use m_static_lib
    use m_static_LIB_1d

    type (hecmwST_local_mesh) :: hecMESH                   !< hecmw mesh
    type (hecmwST_matrix    ) :: hecMAT                    !< hemcw mesh
    type (fstr_solid)      :: fstrSOLID                    !< fstr_solid
    real(kind=kreal)       :: tdstrain(:)                  !< nodal strain
    real(kind=kreal)       :: tdstress(:)                  !< nodal stress
!** Local variables
    integer(kind=kint) :: i, ic_type, icel, ID_area, ie, ielem, ig, ig0
    integer(kind=kint) :: ii, ik , in, iS, iS0, iE0, itype, j, jj, jS, k, nn
    REAL(kind=kreal) xx(20),yy(20),zz(20)
    REAL(kind=kreal) edisp(60),estrain(6),estress(6)
    REAL(kind=kreal) tt(20),tt0(20)
    integer(kind=kint) nodLOCAL(20)
    REAL(kind=kreal) edstrain(20,6),edstress(20,6)
    integer(kind=kint) :: arrayTotal
    real(kind=kreal) val
    real(kind=kreal) s11, s22, s33, s12, s23, s13, ps, smises
    real(kind=kreal) ddunode(3,20)

    real(kind=kreal), allocatable :: ndstrain(:,:), ndstress(:,:)
    integer(kind=kint), allocatable :: nnumber(:)
    real(kind=kreal), allocatable :: temp(:)

    allocate( ndstrain(hecMESH%n_node,6), ndstress(hecMESH%n_node,7) )
    allocate( nnumber(hecMESH%n_node) )
    allocate( temp(hecMESH%n_node) )

!*ZERO CLEAR
    ndstrain=0.d0; ndstress=0.d0
    arrayTotal=0
    nnumber=0
!
!!  Set Temperature
!    temp=0
!    if( fstrSOLID%TEMP_ngrp_tot > 0 ) then
!      do ig0= 1, fstrSOLID%TEMP_ngrp_tot
!        ig= fstrSOLID%TEMP_ngrp_ID(ig0)
!        val=fstrSOLID%TEMP_ngrp_val(ig0)
!!C START & END
!        iS0= hecMESH%node_group%grp_index(ig-1) + 1
!        iE0= hecMESH%node_group%grp_index(ig  )
!         do ik= iS0, iE0
!           in   = hecMESH%node_group%grp_item(ik)
!           temp( in ) = val
!        enddo
!      enddo
!    endif
!C
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
    do itype= 1, hecMESH%n_elem_type
      iS= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
      nn = hecmw_get_max_node(ic_type)
!C element loop
      do icel= iS, iE

!  --- calculate the nodal stress
        if( ic_type==301 ) then
          call NodalStress_C1( ic_type, nn,                                       &
                             fstrSOLID%elements(icel)%gausses,                  &
                             edstrain(1:nn,:), edstress(1:nn,:)      )
        else
          call NodalStress_C3( ic_type, nn,                                       &
                             fstrSOLID%elements(icel)%gausses,                  &
                             edstrain(1:nn,:), edstress(1:nn,:)      )
        endif

!        write(*,'(a,i5,a)') '   edstrain  ( icel', icel, ')'
!        do j=1,nn
!          write(*,'(i5,1p12e15.7)') j,edstrain(j,1:6),edstress(j,1:6)
!        enddo
!
        jS= hecMESH%elem_node_index(icel-1)
        do j=1,nn
          nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
          ndstrain( nodLocal(j),: ) = ndstrain(nodLOCAL(j),:) + edstrain(j,:)
          ndstress( nodLocal(j),1:6 ) = ndstress(nodLOCAL(j),1:6) + edstress(j,:)
          nnumber( nodLOCAL(j) )=nnumber( nodLOCAL(j) )+1
        enddo
!
!
!C** elem ID
!!!          ielem = hecMESH%elem_ID(icel*2-1)
        ielem = icel
        ID_area = hecMESH%elem_ID(icel*2)
        if( ID_area==hecMESH%my_rank ) then
          if( ic_type==301 ) then
            call ElementStress_C1( ic_type,                                       &
                                 fstrSOLID%elements(icel)%gausses,              &
                                 estrain, estress     )
          else
            call ElementStress_C3( ic_type,                                       &
                                 fstrSOLID%elements(icel)%gausses,              &
                                 estrain, estress     )
          endif
          fstrSOLID%ESTRAIN(6*ielem-5:6*ielem) = estrain
          fstrSOLID%ESTRESS(7*ielem-6:7*ielem-1) = estress
          s11=fstrSOLID%ESTRESS(7*ielem-6)
          s22=fstrSOLID%ESTRESS(7*ielem-5)
          s33=fstrSOLID%ESTRESS(7*ielem-4)
          s12=fstrSOLID%ESTRESS(7*ielem-3)
          s23=fstrSOLID%ESTRESS(7*ielem-2)
          s13=fstrSOLID%ESTRESS(7*ielem-1)
          ps=(s11+s22+s33)/3.0
          smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                +s12**2+s23**2+s13**2
          fstrSOLID%ESTRESS(7*ielem)=sqrt(3.0*smises)
        endif
!
      enddo       ! icel roop.
    enddo         ! ityp roop.

!    stop

!** Average over nodes
    do i=1,hecMESH%n_node
      ndstrain(i,:)=ndstrain(i,:)/nnumber(i)
      ndstress(i,1:6)=ndstress(i,1:6)/nnumber(i)
    enddo
!** CALCULATE von MISES stress
    do i=1,hecMESH%n_node
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
    do i=1,hecMESH%n_node
      do j=1,6
  !      fstrSOLID%STRAIN(6*(i-1)+j)=ndstrain(i,j)
		tdSTRAIN(6*(i-1)+j)=ndstrain(i,j)
      enddo
      do j=1,7
   !     fstrSOLID%STRESS(7*(i-1)+j)=ndstress(i,j)
		tdSTRESS(7*(i-1)+j)=ndstress(i,j)
      enddo
    enddo

    deallocate( ndstrain, ndstress)
    deallocate( nnumber)
    deallocate( temp )
  end subroutine fstr_NodalStress3D


end module m_fstr_NodalStress
