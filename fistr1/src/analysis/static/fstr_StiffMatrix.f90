!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides function to calcualte tangent stiffness matrix.

module m_fstr_StiffMatrix
  use m_fstr
  implicit none

  private
  public :: fstr_StiffMatrix

contains

  !---------------------------------------------------------------------*
  !> \brief 接線剛性マトリックスを作成するサブルーチン
  subroutine fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, time, tincr)
    !---------------------------------------------------------------------*
    use m_static_LIB
    use mMechGauss

    type (hecmwST_local_mesh)  :: hecMESH      !< mesh information
    type (hecmwST_matrix)      :: hecMAT       !< linear equation, its right side modified here
    type (fstr_solid)          :: fstrSOLID    !< we need boundary conditions of curr step
    real(kind=kreal),intent(in) :: time        !< current time
    real(kind=kreal),intent(in) :: tincr       !< time increment

    type( tMaterial ), pointer :: material     !< material information

    real(kind=kreal)   :: stiffness(20*6, 20*6)
    integer(kind=kint) :: ierror, nodLOCAL(20)
    real(kind=kreal)   :: tt(20), ecoord(3,20)
    real(kind=kreal)   :: thick, val, pa1
    integer(kind=kint) :: ndof, itype, is, iE, ic_type, nn, icel, iiS, i, j
    real(kind=kreal)   :: u(6,20), du(6,20), coords(3,3), u_prev(6,20)
    integer            :: ig0, grpid, ig, iS0, iE0,ik, in, isect, ihead, cdsys_ID

    ! ----- initialize
    call hecmw_mat_clear( hecMAT )

    ndof = hecMAT%NDOF
    do itype= 1, hecMESH%n_elem_type
      is= hecMESH%elem_type_index(itype-1) + 1
      iE= hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      ! ----- Ignore link elements
      if (hecmw_is_etype_link(ic_type)) cycle
      ! ----- Set number of nodes
      nn = hecmw_get_max_node(ic_type)

      ! ----- element loop
      !$omp parallel default(none), &
        !$omp&  private(icel,iiS,j,nodLOCAL,i,ecoord,du,u,u_prev,tt,cdsys_ID,coords, &
        !$omp&          material,thick,stiffness,isect,ihead), &
        !$omp&  shared(iS,iE,hecMESH,nn,ndof,fstrSOLID,ic_type,hecMAT,time,tincr)
      !$omp do
      do icel= is, iE

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
            u_prev(i,j) = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof)
          enddo
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 )  &
            tt(j)=fstrSOLID%temperature( nodLOCAL(j) )
        enddo

        isect = hecMESH%section_ID(icel)
        cdsys_ID = hecMESH%section%sect_orien_ID(isect)
        if( cdsys_ID > 0 ) call get_coordsys(cdsys_ID, hecMESH, fstrSOLID, coords)

        material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
        thick = material%variables(M_THICK)
        if( getSpaceDimension( ic_type )==2 ) thick =1.d0
        if( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 .or. ic_type==2322) then
          if( material%nlgeom_flag /= INFINITE ) call StiffMat_abort( ic_type, 2 )
          call STF_C2( ic_type,nn,ecoord(1:2,1:nn),fstrSOLID%elements(icel)%gausses(:),thick,  &
            stiffness(1:nn*ndof,1:nn*ndof), fstrSOLID%elements(icel)%iset,          &
            u(1:2,1:nn) )

        elseif ( ic_type==301 ) then
          isect= hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call STF_C1( ic_type,nn,ecoord(:,1:nn),thick,fstrSOLID%elements(icel)%gausses(:),   &
            stiffness(1:nn*ndof,1:nn*ndof), u(1:3,1:nn) )

        elseif ( ic_type==361 ) then

          if( fstrSOLID%sections(isect)%elemopt361 == kel361FI ) then ! full integration element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
              call STF_C3                                                                              &
                ( ic_type, nn, ecoord(:, 1:nn), fstrSOLID%elements(icel)%gausses(:),                &
                stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
            else
              call STF_C3                                                                     &
                ( ic_type,nn,ecoord(:, 1:nn),fstrSOLID%elements(icel)%gausses(:),          &
                stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn) )
            endif
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361BBAR ) then ! B-bar element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
              call STF_C3D8Bbar                                                                        &
                ( ic_type, nn, ecoord(:, 1:nn), fstrSOLID%elements(icel)%gausses(:),                &
                stiffness(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn), tt(1:nn) )
            else
              call STF_C3D8Bbar                                                               &
                ( ic_type, nn, ecoord(:, 1:nn),fstrSOLID%elements(icel)%gausses(:),        &
                stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn) )
            endif
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361IC ) then ! incompatible element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
              CALL STF_C3D8IC                                                              &
                ( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
                stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), &
                fstrSOLID%elements(icel)%aux, tt(1:nn) )
            else
              CALL STF_C3D8IC                                                              &
                ( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
                stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), &
                fstrSOLID%elements(icel)%aux )
            endif
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361FBAR ) then ! F-bar element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
              call STF_C3D8Fbar                                                                        &
                ( ic_type, nn, ecoord(:, 1:nn), fstrSOLID%elements(icel)%gausses(:),                &
                stiffness(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn), tt(1:nn) )
            else
              call STF_C3D8Fbar                                                               &
                ( ic_type, nn, ecoord(:, 1:nn),fstrSOLID%elements(icel)%gausses(:),        &
                stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn) )
            endif
          endif

        elseif (ic_type==341 .or. ic_type==351 .or.                     &
            ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
            call STF_C3                                                                              &
              ( ic_type, nn, ecoord(:, 1:nn), fstrSOLID%elements(icel)%gausses(:),                &
              stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
          else
            call STF_C3                                                                     &
              ( ic_type,nn,ecoord(:, 1:nn),fstrSOLID%elements(icel)%gausses(:),          &
              stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3, 1:nn) )
          endif

        else if( ic_type == 611) then
          if( material%nlgeom_flag /= INFINITE ) call StiffMat_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          call STF_Beam(ic_type, nn, ecoord, hecMESH%section%sect_R_item(ihead+1:), &
            &   material%variables(M_YOUNGS), material%variables(M_POISSON), stiffness(1:nn*ndof,1:nn*ndof))

        else if( ic_type == 641 ) then
          if( material%nlgeom_flag /= INFINITE ) call StiffMat_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          call STF_Beam_641(ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses(:), &
            &            hecMESH%section%sect_R_item(ihead+1:), stiffness(1:nn*ndof,1:nn*ndof))

        else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
          if( material%nlgeom_flag /= INFINITE ) call StiffMat_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:nn), fstrSOLID%elements(icel)%gausses(:), &
            &              stiffness(1:nn*ndof, 1:nn*ndof), thick, 0)

        else if( ic_type == 761 ) then   !for shell-solid mixed analysis
          if( material%nlgeom_flag /= INFINITE ) call StiffMat_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call STF_Shell_MITC(731, 3, 6, ecoord(1:3, 1:3), fstrSOLID%elements(icel)%gausses(:), &
            &              stiffness(1:nn*ndof, 1:nn*ndof), thick, 2)

        else if( ic_type == 781 ) then   !for shell-solid mixed analysis
          if( material%nlgeom_flag /= INFINITE ) call StiffMat_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call STF_Shell_MITC(741, 4, 6, ecoord(1:3, 1:4), fstrSOLID%elements(icel)%gausses(:), &
            &              stiffness(1:nn*ndof, 1:nn*ndof), thick, 1)

        elseif ( ic_type==3414 ) then
          if(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype /= INCOMP_NEWTONIAN) then
            write(*, *) '###ERROR### : This element is not supported for this material'
            write(*, *) 'ic_type = ', ic_type, ', mtype = ', fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype
            call hecmw_abort(hecmw_comm_get_comm())
          endif
          call STF_C3_vp                                                           &
            ( ic_type, nn, ecoord(:, 1:nn),fstrSOLID%elements(icel)%gausses(:), &
            stiffness(1:nn*ndof, 1:nn*ndof), tincr, u_prev(1:4, 1:nn) )
          !
          !      elseif ( ic_type==731) then
          !        call STF_S3(xx,yy,zz,ee,pp,thick,local_stf)
          !        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
          !      elseif ( ic_type==741) then
          !        call STF_S4(xx,yy,zz,ee,pp,thick,local_stf)
          !        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
        else
          call StiffMat_abort( ic_type, 1 )
        endif
        !
        ! ----- CONSTRUCT the GLOBAL MATRIX STARTED
        call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiffness)

      enddo      ! icel
      !$omp end do
      !$omp end parallel
    enddo        ! itype

  end subroutine fstr_StiffMatrix

  subroutine StiffMat_abort( ic_type, flag )
    integer(kind=kint), intent(in) :: ic_type
    integer(kind=kint), intent(in) :: flag

    if( flag == 1 ) then
      write(*,*) '###ERROR### : Element type not supported for static analysis'
    else if( flag == 2 ) then
      write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
    endif
    write(*,*) ' ic_type = ', ic_type
    call hecmw_abort(hecmw_comm_get_comm())
  end subroutine

end module m_fstr_StiffMatrix
