!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!> \brief  This module provides function to calcualte tangent stiffness matrix.
!!
!>  \author                date                  version 
!>  X.Yuan(Advancesoft)    2010/08/26        original
!>  X.Yuan                 2013/03/18        consider anisotropic materials
!
!======================================================================!
module m_fstr_StiffMatrix
   implicit none

   contains

!---------------------------------------------------------------------*
!> \brief 接線剛性マトリックスを作成するサブルーチン
subroutine fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, tincr)
!---------------------------------------------------------------------*
  use m_fstr
  use m_static_LIB
  use mMechGauss

  type (hecmwST_local_mesh)  :: hecMESH      !< mesh information
  type (hecmwST_matrix)      :: hecMAT       !< linear equation, its right side modified here
  type (fstr_solid)          :: fstrSOLID    !< we need boundary conditions of curr step
  real(kind=kreal),intent(in) :: tincr       !< time increment

  type( tMaterial ), pointer :: material     !< material information

  real(kind=kreal)   :: stiffness(20*6, 20*6)
  integer(kind=kint) :: ierror, nodLOCAL(20)
  real(kind=kreal)   :: tt(20), ecoord(3,20)
  real(kind=kreal)   :: thick, val, pa1
  integer(kind=kint) :: ndof, itype, iS, iE, ic_type, nn, icel, iiS, i, j
  real(kind=kreal)   :: u(3,20), du(3,20), coords(3,3)
  integer            :: ig0, grpid, ig, iS0, iE0,ik, in, isect, ihead, cdsys_ID

! ----- initialize
  call hecmw_mat_clear( hecMAT )

  ndof = hecMAT%NDOF
  do itype= 1, hecMESH%n_elem_type
    iS= hecMESH%elem_type_index(itype-1) + 1
    iE= hecMESH%elem_type_index(itype  )
    ic_type= hecMESH%elem_type_item(itype)
! ----- Ignore link elements
    if (hecmw_is_etype_link(ic_type)) cycle
! ----- Set number of nodes
    nn = hecmw_get_max_node(ic_type)

! ----- element loop
    do icel= iS, iE

! ----- nodal coordinate & displacement
      iiS= hecMESH%elem_node_index(icel-1)
      do j=1,nn
        nodLOCAL(j)= hecMESH%elem_node_item (iiS+j)
        do i=1, 3
          ecoord(i,j) = hecMESH%node(3*nodLOCAL(j)+i-3)
        enddo
        do i=1,ndof
          du(i,j) = fstrSOLID%dunode(ndof*nodLOCAL(j)+i-ndof)
          u(i,j)  = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof) + du(i,j)
        enddo
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 )  &
           tt(j)=fstrSOLID%temperature( nodLOCAL(j) )
      enddo
	  
      cdsys_ID = fstrSOLID%elements(icel)%gausses(1)%pMaterial%cdsys_ID
      if( cdsys_ID>0 ) call get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )

      material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
      thick = material%variables(M_THICK)
      if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
      if ( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 .or. ic_type==2322) then
        call STF_C2( ic_type,nn,ecoord(1:2,1:nn),fstrSOLID%elements(icel)%gausses(:),thick,  &
                     stiffness(1:nn*ndof,1:nn*ndof), fstrSOLID%elements(icel)%iset,          &
                     u(1:2,1:nn) )
					 
      else if ( ic_type==301 ) then
        isect= hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        thick = hecMESH%section%sect_R_item(ihead+1)
        call STF_C1( ic_type,nn,ecoord(:,1:nn),thick,fstrSOLID%elements(icel)%gausses(:),   &
            stiffness(1:nn*ndof,1:nn*ndof), u(1:3,1:nn) )

      else if ( ic_type==361 ) then
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          call STF_C3D8Bbar( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, coords, u(1:3,1:nn), tt(1:nn) )
        else
          call STF_C3D8Bbar( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, coords, u(1:3,1:nn) )
        endif


      else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.                    &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, coords, u(1:3,1:nn), tt(1:nn) )
        else
          call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses(:),   &
                     stiffness(1:nn*ndof,1:nn*ndof), tincr, coords, u(1:3,1:nn) )
        endif


!
!      else if ( ic_type==731) then
!        call STF_S3(xx,yy,zz,ee,pp,thick,local_stf)
!        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
!      else if ( ic_type==741) then
!        call STF_S4(xx,yy,zz,ee,pp,thick,local_stf)
!        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
      else
        write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
        write(*,*) ' ic_type = ', ic_type
        call hecmw_abort(hecmw_comm_get_comm())
      endif
!
! ----- CONSTRUCT the GLOBAL MATRIX STARTED
      call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiffness)

    enddo      ! icel
  enddo        ! itype

end subroutine fstr_StiffMatrix

end module m_fstr_StiffMatrix
