!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides function to calculate tangent stiffness matrix.

module m_fstr_StiffMatrix
  use m_fstr
  use m_static_LIB
  use m_fstr_Update
  implicit none

  private
  public :: fstr_StiffMatrix

contains

  !---------------------------------------------------------------------*
  !> \brief This subroutine creates tangential stiffness matrix
  subroutine fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, fstrPARAM, time, tincr)
    !---------------------------------------------------------------------*
    use mMechGauss

    type (hecmwST_local_mesh)  :: hecMESH      !< mesh information
    type (hecmwST_matrix)      :: hecMAT       !< linear equation, its right side modified here
    type (fstr_solid)          :: fstrSOLID    !< we need boundary conditions of curr step
    type (fstr_param)          :: fstrPARAM !< type fstr_param
    real(kind=kreal),intent(in) :: time        !< current time
    real(kind=kreal),intent(in) :: tincr       !< time increment

    real(kind=kreal)   :: stiffness(fstrSOLID%max_ncon_stf*6, fstrSOLID%max_ncon_stf*6)
    integer(kind=kint) :: nodLOCAL(fstrSOLID%max_ncon)
    real(kind=kreal)   :: tt(fstrSOLID%max_ncon), tt0(fstrSOLID%max_ncon), ttn(fstrSOLID%max_ncon), ecoord(3,fstrSOLID%max_ncon)
    real(kind=kreal)   :: thick
    integer(kind=kint) :: ndof, itype, is, iE, ic_type, nn, icel, iiS, i, j
    real(kind=kreal)   :: u(6,fstrSOLID%max_ncon), du(6,fstrSOLID%max_ncon), coords(3,3), u_prev(6,fstrSOLID%max_ncon)
    integer            :: isect, ihead, cdsys_ID, initt
    integer            :: stftype

    ! ----- initialize
    call hecmw_mat_clear( hecMAT )
    stftype = fstrPARAM%stiffness_type
    ndof = hecMAT%NDOF

    if( stftype == kstfNUMERICAL ) initt = is_initial_temperature_exists( g_InitialCnd )

    do itype= 1, hecMESH%n_elem_type
      is= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      ! ----- Ignore link and patch elements
      if (hecmw_is_etype_link(ic_type)) cycle
      if (hecmw_is_etype_patch(ic_type)) cycle
      ! ----- Set number of nodes

      ! ----- element loop
      !$omp parallel default(none), &
        !$omp&  private(icel,iiS,nn,j,nodLOCAL,i,ecoord,du,u,u_prev,tt,tt0,ttn,cdsys_ID,coords, &
        !$omp&          thick,stiffness,isect,ihead), &
        !$omp&  shared(iS,iE,hecMESH,ndof,fstrSOLID,ic_type,hecMAT,time,tincr,initt,g_InitialCnd,stftype)
      !$omp do
      do icel= is, iE

        ! ----- nodal coordinate & displacement
        iiS= hecMESH%elem_node_index(icel-1)
        nn = hecMESH%elem_node_index(icel)-iiS
     !   if( nn>150 ) stop "elemental nodes > 150!"

        do j=1,nn
          nodLOCAL(j)= hecMESH%elem_node_item (iiS+j)
          do i=1, 3
            ecoord(i,j) = hecMESH%node(3*nodLOCAL(j)+i-3)
          enddo
          do i=1,ndof
            du(i,j) = fstrSOLID%dunode(ndof*nodLOCAL(j)+i-ndof)
            u(i,j)  = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof) + du(i,j)
            u_prev(i,j) = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof)
          enddo
          tt(j) = 0.d0
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 )  &
            tt(j)=fstrSOLID%temperature( nodLOCAL(j) )
        enddo

        isect = hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        cdsys_ID = hecMESH%section%sect_orien_ID(isect)
        if( cdsys_ID > 0 ) call get_coordsys(cdsys_ID, hecMESH, fstrSOLID, coords)
        thick = hecMESH%section%sect_R_item(ihead+1)
        if( getSpaceDimension( ic_type )==2 ) thick =1.d0

        if( stftype == kstfANALYTICAL ) then
          call fstr_StiffMatrix_elem( ic_type, nn, ndof, nodLOCAL, ecoord, u, du, u_prev, fstrSOLID%sections(isect), & 
            & hecMESH%section%sect_R_item(ihead+1:), cdsys_ID, coords, stiffness, fstrSOLID%elements(icel), &
            & time, tincr, tt, thick )
        else if( stftype == kstfNUMERICAL ) then
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
            if( isElastoplastic(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype) .or. &
                fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype == NORTON ) then
              tt0(j)=fstrSOLID%last_temp( nodLOCAL(j) )
            else
              tt0(j) = 0.d0
              if( hecMESH%hecmw_flag_initcon == 1 ) tt0(j) = hecMESH%node_init_val_item(nodLOCAL(j))
              if( initt>0 ) tt0(j) = g_InitialCnd(initt)%realval(nodLOCAL(j))
            endif
            ttn(j) = fstrSOLID%last_temp( nodLOCAL(j) )
          else
            tt0 = 0.d0
            ttn = 0.d0
          endif
          call fstr_StiffMatrix_elem_numerical( ic_type, nn, ndof, nodLOCAL, ecoord, u, du, u_prev, fstrSOLID%sections(isect), & 
            & hecMESH%section%sect_R_item(ihead+1:), cdsys_ID, coords, stiffness, fstrSOLID%elements(icel), &
            & time, tincr, tt, tt0, ttn, thick )

        else
          stop
        endif

        if( ic_type==341 .and. fstrSOLID%sections(isect)%elemopt341 == kel341SESNS ) cycle ! skip smoothed fem
        !
        ! ----- CONSTRUCT the GLOBAL MATRIX STARTED
        call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiffness)

      enddo      ! icel
      !$omp end do
      !$omp end parallel
    enddo        ! itype

  end subroutine fstr_StiffMatrix

  subroutine fstr_StiffMatrix_elem( ic_type, nn, ndof, nodLOCAL, ecoord, u, du, u_prev, section, sect_R_item, &
    & cdsys_ID, coords, stiffness, element, time, tincr, tt, thick )
    integer, intent(in)                :: ic_type      !< element type
    integer, intent(inout)             :: nn
    integer, intent(in)                :: ndof
    integer(kind=kint), intent(inout)  :: nodLOCAL(:)
    real(kind=kreal),intent(in)        :: ecoord(:,:) 
    real(kind=kreal),intent(in)        :: u(:,:) 
    real(kind=kreal),intent(in)        :: du(:,:) 
    real(kind=kreal),intent(in)        :: u_prev(:,:)
    type(tSection),intent(in)          :: section
    real(kind=kreal), intent(in)       :: sect_R_item(:)      
    integer, intent(in)                :: cdsys_ID     
    real(kind=kreal), intent(inout)    :: coords(3,3) 
    real(kind=kreal), intent(out)      :: stiffness(:,:) 
    type(tElement), intent(inout)      :: element    
    real(kind=kreal),intent(in)        :: time      !< current time
    real(kind=kreal),intent(in)        :: tincr     !< time increment
    real(kind=kreal),intent(in)        :: tt(:) 
    real(kind=kreal),intent(in)        :: thick

    type( tMaterial ), pointer :: material     !< material information

    material => element%gausses(1)%pMaterial

    if( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 .or. ic_type==2322) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call StiffMat_abort( ic_type, 2 )
      call STF_C2( ic_type,nn,ecoord(1:2,1:nn),element%gausses(:),thick,  &
        stiffness(1:nn*ndof,1:nn*ndof), element%iset,          &
        u(1:2,1:nn) )

    elseif ( ic_type==301 ) then
      call STF_C1( ic_type,nn,ecoord(:,1:nn),thick,element%gausses(:),   &
        stiffness(1:nn*ndof,1:nn*ndof), u(1:3,1:nn) )

    elseif ( ic_type==361 ) then
      if( section%elemopt361 == kel361FI ) then ! full integration element
        call STF_C3                                                                              &
          ( ic_type, nn, ecoord(:, 1:nn), element%gausses(:),                &
          stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
      else if( section%elemopt361 == kel361BBAR ) then ! B-bar element
        call STF_C3D8Bbar                                                                        &
          ( ic_type, nn, ecoord(:, 1:nn), element%gausses(:),                &
          stiffness(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn), tt(1:nn) )
      else if( section%elemopt361 == kel361IC ) then ! incompatible element
        call STF_C3D8IC                                                              &
          ( ic_type, nn, ecoord(:,1:nn), element%gausses(:), &
          stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), &
          element%aux, tt(1:nn) )
      else if( section%elemopt361 == kel361FBAR ) then ! F-bar element
        call STF_C3D8Fbar                                                                        &
          ( ic_type, nn, ecoord(:, 1:nn), element%gausses(:),                &
          stiffness(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn), tt(1:nn) )
      endif

    elseif (ic_type==341 .or. ic_type==351 .or. ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
      if( ic_type==341 .and. section%elemopt341 == kel341SESNS ) return ! skip smoothed fem
      call STF_C3                                                                              &
        ( ic_type, nn, ecoord(:, 1:nn), element%gausses(:),                &
        stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )

    else if( ic_type == 511) then
      call STF_CONNECTOR( ic_type,nn,ecoord(:,1:nn),element%gausses(:),   &
        stiffness(1:nn*ndof,1:nn*ndof), u(1:3,1:nn), tincr, tt(1:nn) )

    else if( ic_type == 611) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call StiffMat_abort( ic_type, 2 )
      call STF_Beam(ic_type, nn, ecoord, sect_R_item, &
        &   material%variables(M_YOUNGS), material%variables(M_POISSON), stiffness(1:nn*ndof,1:nn*ndof))

    else if( ic_type == 641 ) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call StiffMat_abort( ic_type, 2 )
      call STF_Beam_641(ic_type, nn, ecoord, element%gausses(:), &
        &            sect_R_item, stiffness(1:nn*ndof,1:nn*ndof))

    else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call StiffMat_abort( ic_type, 2 )
      call STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:nn), element%gausses(:), &
        &              stiffness(1:nn*ndof, 1:nn*ndof), thick, 0)

    else if( ic_type == 761 ) then   !for shell-solid mixed analysis
      if( material%nlgeom_flag /= INFINITESIMAL ) call StiffMat_abort( ic_type, 2 )
      call STF_Shell_MITC(731, 3, 6, ecoord(1:3, 1:3), element%gausses(:), &
        &              stiffness(1:nn*ndof, 1:nn*ndof), thick, 2)

    else if( ic_type == 781 ) then   !for shell-solid mixed analysis
      if( material%nlgeom_flag /= INFINITESIMAL ) call StiffMat_abort( ic_type, 2 )
      call STF_Shell_MITC(741, 4, 6, ecoord(1:3, 1:4), element%gausses(:), &
        &              stiffness(1:nn*ndof, 1:nn*ndof), thick, 1)

    elseif ( ic_type==3414 ) then
      if( material%mtype /= INCOMP_NEWTONIAN) call StiffMat_abort( ic_type, 3, material%mtype )
      call STF_C3_vp                                                           &
        ( ic_type, nn, ecoord(:, 1:nn),element%gausses(:), &
        stiffness(1:nn*ndof, 1:nn*ndof), tincr, u_prev(1:4, 1:nn) )
    else if ( ic_type == 881 .or. ic_type == 891 ) then  !for selective es/ns smoothed fem
      call STF_C3D4_SESNS                                                                   &
        ( ic_type,nn,nodLOCAL,ecoord(:, 1:nn), element%gausses(:),         &
        stiffness, cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
    else
      call StiffMat_abort( ic_type, 1 )
    endif    
  end subroutine

  subroutine fstr_StiffMatrix_elem_numerical( ic_type, nn, ndof, nodLOCAL, ecoord, u, du, u_prev, section, sect_R_item, &
    & cdsys_ID, coords, stiffness, element, time, tincr, tt, tt0, ttn, thick )
    integer, intent(in)                :: ic_type      !< element type
    integer, intent(inout)             :: nn
    integer, intent(in)                :: ndof
    integer(kind=kint), intent(inout)  :: nodLOCAL(:)
    real(kind=kreal),intent(in)        :: ecoord(:,:) 
    real(kind=kreal),intent(in)        :: u(:,:) 
    real(kind=kreal),intent(in)        :: du(:,:) 
    real(kind=kreal),intent(in)        :: u_prev(:,:)
    type(tSection),intent(in)          :: section
    real(kind=kreal), intent(in)       :: sect_R_item(:)      
    integer, intent(in)                :: cdsys_ID     
    real(kind=kreal), intent(inout)    :: coords(3,3) 
    real(kind=kreal), intent(out)      :: stiffness(:,:) 
    type(tElement), intent(inout)      :: element    
    real(kind=kreal),intent(in)        :: time      !< current time
    real(kind=kreal),intent(in)        :: tincr     !< time increment
    real(kind=kreal),intent(in)        :: tt(:) 
    real(kind=kreal),intent(in)        :: tt0(:) 
    real(kind=kreal),intent(in)        :: ttn(:) 
    real(kind=kreal),intent(in)        :: thick

    integer(kind=kint) :: i, idof, idx, j, jdof, jdx
    real(kind=kreal) :: u_ori(ndof,nn), du_tmp(ndof,nn), ddu(ndof,nn), pdu(1:ndof)
    real(kind=kreal) :: qfp(ndof*nn),qfm(ndof*nn),qfp2(ndof*nn),qfm2(ndof*nn),tmpval,delta_u

    if( ic_type==341 .and. section%elemopt341 == kel341SESNS ) stop "SELECTIVE ESNS does not support numerical stiffness"
    if( ic_type==361 .and. section%elemopt361 == kel361IC ) stop " incompatible element does not support numerical stiffness"

    ddu(1:ndof,1:nn) = 0.d0

    call get_bounding_box(nn,ndof,du,pdu(1:ndof))
    delta_u = max(1.d-6*max(pdu(1),pdu(2),pdu(3)),1.d-8)

    do i=1,nn
      do idof=1,ndof
        idx = (i-1)*ndof+idof

        ! calc f(u+ddu)
        du_tmp(1:ndof,1:nn) = du(1:ndof,1:nn)
        du_tmp(idof,i) = du_tmp(idof,i) + delta_u
        call fstr_UpdateNewton_elem( ic_type, nn, ndof, nodLOCAL, ecoord, u_prev, du_tmp, ddu, section, & 
        & sect_R_item, cdsys_ID, coords, qfp, element, &
        & 1, time, tincr, tt, tt0, ttn, thick )
    
        ! calc f(u-ddu)
        du_tmp(idof,i) = du_tmp(idof,i) - 2.d0*delta_u
        call fstr_UpdateNewton_elem( ic_type, nn, ndof, nodLOCAL, ecoord, u_prev, du_tmp, ddu, section, & 
        & sect_R_item, cdsys_ID, coords, qfm, element, &
        & 1, time, tincr, tt, tt0, ttn, thick )

        ! calc df/du as (f(u+ddu)-f(u-ddu))/2du
        do j=1,nn
          do jdof=1,ndof
            jdx = (j-1)*ndof+jdof
            stiffness(idx,jdx) = (qfp(jdx)-qfm(jdx))/(2.d0*delta_u)
          enddo
        enddo
  
      enddo
    enddo

  contains

    subroutine get_bounding_box(nn,ndof,ecoord,esize)
      integer, intent(in)             :: nn
      integer, intent(in)             :: ndof
        real(kind=kreal),intent(in)   :: ecoord(:,:) 
      real(kind=kreal),intent(out)    :: esize(:)

      integer(kind=kint) :: i, idof
      real(kind=kreal) :: minec(3), maxec(3)

      if( nn == 1 ) then
        esize(1:ndof) = 1.d0
        return
      endif

      minec(1:ndof) = ecoord(1:ndof,1)
      maxec(1:ndof) = ecoord(1:ndof,1)

      do i=2,nn
        do idof=1,ndof
          minec(idof) = min(minec(idof),ecoord(idof,i))
          maxec(idof) = max(maxec(idof),ecoord(idof,i))
        enddo
      enddo

      do idof=1,ndof
        esize(idof) = maxec(idof)-minec(idof)
      enddo
    end subroutine

  end subroutine

  subroutine StiffMat_abort( ic_type, flag, mtype )
    integer(kind=kint), intent(in)           :: ic_type
    integer(kind=kint), intent(in)           :: flag
    integer(kind=kint), intent(in), optional :: mtype

    if( flag == 1 ) then
      write(*,*) '###ERROR### : Element type not supported for static analysis'
    else if( flag == 2 ) then
      write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
    else if( flag == 3 ) then
      write(*,*) '###ERROR### : This element is not supported for this material'
    endif
    write(*,*) ' ic_type = ', ic_type
    if( present(mtype) ) write(*,*) ' mtype = ', mtype
    call hecmw_abort(hecmw_comm_get_comm())
  end subroutine

end module m_fstr_StiffMatrix
