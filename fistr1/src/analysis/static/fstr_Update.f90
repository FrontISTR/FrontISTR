!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides function to calcualte to do updates
module m_fstr_Update
  use m_fstr
  implicit none

  private :: Update_abort

contains

  !=====================================================================*
  !  UPDATE_C3
  !>  \brief 変位／応力・ひずみ／内力のアップデート
  !>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
  !>  \version    0.00
  !> \if AS
  !>    \par  サブルーチン構成
  !>    -# 変位の更新              \f$ u_{n+1}^{(k)} = u_{n+1}^{(k-1)} + \delta u^{(k)} \f$
  !>    -# ひずみ・応力の更新       \f$ \varepsilon_{n+1}^{(k)} = \varepsilon_{n+1}^{(k-1)} + \delta \varepsilon^{(k)} \f$, \f$ \sigma_{n+1}^{(k)} = \sigma_{n+1}^{(k-1)} + \delta \sigma^{(k)} \f$
  !>    -# 内力（等価節点力）の計算 \f$ Q_{n+1}^{(k-1)} ( u_{n+1}^{(k-1)} ) \f$
  !> \endif
  subroutine fstr_UpdateNewton ( hecMESH, hecMAT, fstrSOLID, time, tincr,iter, strainEnergy)
    !=====================================================================*
    use m_static_lib

    type (hecmwST_matrix)       :: hecMAT    !< linear equation, its right side modified here
    type (hecmwST_local_mesh)   :: hecMESH   !< mesh information
    type (fstr_solid)           :: fstrSOLID !< we need boundary conditions of curr step
    real(kind=kreal),intent(in) :: time      !< current time
    real(kind=kreal),intent(in) :: tincr     !< time increment
    integer, intent(in)         :: iter      !< NR iterations

    integer(kind=kint) :: nodLOCAL(20)
    real(kind=kreal)   :: ecoord(3, 20), stiff(60, 60)
    real(kind=kreal)   :: thick
    integer(kind=kint) :: ndof, itype, is, iE, ic_type, nn, icel, iiS, i, j

    real(kind=kreal)   :: total_disp(6, 20), du(6, 20), ddu(6, 20)
    real(kind=kreal)   :: tt(20), tt0(20), ttn(20), qf(20*6), coords(3, 3)
    integer            :: ig0, ig, ik, in, ierror, isect, ihead, cdsys_ID
    integer            :: ndim

    real(kind=kreal), optional :: strainEnergy
    real(kind=kreal) :: tmp
    real(kind=kreal)   :: ddaux(3,3)

    ndof = hecMAT%NDOF
    fstrSOLID%QFORCE=0.0d0

    tt0 = 0.d0
    ttn = 0.d0
    tt = 0.d0

    ! --------------------------------------------------------------------
    !      updated
    !   1. stress and strain  : ep^(k) = ep^(k-1)+dep^(k)
    !                           sgm^(k) = sgm^(k-1)+dsgm^(k)
    !   2. Internal Force     : Q^(k-1) ( u^(k-1) )
    ! --------------------------------------------------------------------
    ! ----------------------------------------------------------------------------------
    !      calculate the Strain and Stress and Internal Force ( Equivalent Nodal Force )
    ! ----------------------------------------------------------------------------------

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1)+1
      iE = hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      if (hecmw_is_etype_link(ic_type)) cycle
      nn = hecmw_get_max_node(ic_type)
      if( nn>20 ) stop "Elemental nodes>20!"

      !element loop
      !$omp parallel default(none), &
        !$omp&  private(icel,iiS,j,nodLOCAL,i,ecoord,ddu,du,total_disp, &
        !$omp&  cdsys_ID,coords,thick,qf,isect,ihead,tmp,ndim,ddaux), &
        !$omp&  shared(iS,iE,hecMESH,nn,fstrSOLID,ndof,hecMAT,ic_type,fstrPR, &
        !$omp&         strainEnergy,iter,time,tincr), &
        !$omp&  firstprivate(tt0,ttn,tt)
      !$omp do
      do icel = is, iE

        ! ----- nodal coordinate
        iiS = hecMESH%elem_node_index(icel-1)

        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (iiS+j)
          do i = 1, 3
            ecoord(i,j) = hecMESH%node(3*nodLOCAL(j)+i-3)
          enddo
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
            if( isElastoplastic(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype) .or. &
                fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype == NORTON ) then
              tt0(j)=fstrSOLID%last_temp( nodLOCAL(j) )
            else
              tt0(j) = 0.d0
              if( hecMESH%hecmw_flag_initcon == 1 ) tt0(j) = hecMESH%node_init_val_item(nodLOCAL(j))
            endif
            ttn(j) = fstrSOLID%last_temp( nodLOCAL(j) )
            tt(j)  = fstrSOLID%temperature( nodLOCAL(j) )
          endif
        enddo

        ! nodal displacement
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (iiS+j)
          do i = 1, ndof
            ddu(i,j) = hecMAT%X(ndof*nodLOCAL(j)+i-ndof)
            du(i,j)  = fstrSOLID%dunode(ndof*nodLOCAL(j)+i-ndof)
            total_disp(i,j) = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof)
          enddo
        enddo

        isect = hecMESH%section_ID(icel)
        cdsys_ID = hecMESH%section%sect_orien_ID(isect)
        if( cdsys_ID > 0 ) call get_coordsys(cdsys_ID, hecMESH, fstrSOLID, coords)

        ! ===== calculate the Internal Force
        if( getSpaceDimension( ic_type ) == 2 ) thick = 1.0d0

        if( ic_type == 241 .or. ic_type == 242 .or. ic_type == 231 .or. ic_type == 232 .or. ic_type == 2322 ) then

          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
            call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:), &
              thick,fstrSOLID%elements(icel)%iset,                             &
              total_disp(1:2,1:nn), ddu(1:2,1:nn), qf(1:nn*ndof),              &
              tt(1:nn), tt0(1:nn), ttn(1:nn)  )
          else
            call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:), &
              thick,fstrSOLID%elements(icel)%iset,           &
              total_disp(1:2,1:nn), ddu(1:2,1:nn), qf(1:nn*ndof) )
          endif

        else if( ic_type == 301 ) then
          isect= hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call UPDATE_C1( ic_type,nn,ecoord(:,1:nn), thick, total_disp(1:3,1:nn), du(1:3,1:nn), &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:) )

        else if( ic_type == 361 ) then

          if( fstrSOLID%sections(isect)%elemopt361 == kel361FI ) then ! full integration element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
              call UPDATE_C3( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
            else
              call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr )
            endif
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361BBAR ) then ! B-bar element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
              call UPDATE_C3D8Bbar( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords,    &
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
            else
              call Update_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr )
            endif
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361IC ) then ! incompatible element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
              call UPDATE_C3D8IC( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), ddu(1:3,1:nn), cdsys_ID, coords,&
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, &
                fstrSOLID%elements(icel)%aux, ddaux(1:3,1:3), tt(1:nn), tt0(1:nn), ttn(1:nn) )
            else
              call UPDATE_C3D8IC( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), ddu(1:3,1:nn), cdsys_ID, coords,&
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, &
                fstrSOLID%elements(icel)%aux, ddaux(1:3,1:3) )
            endif
            fstrSOLID%elements(icel)%aux(1:3,1:3) = fstrSOLID%elements(icel)%aux(1:3,1:3) + ddaux(1:3,1:3)
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361FBAR ) then ! F-bar element
            if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
              call UPDATE_C3D8Fbar( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords,    &
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
            else
              call Update_C3D8Fbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
                qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr )
            endif
          endif

        else if (ic_type == 341 .or. ic_type == 351 .or. ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
            call UPDATE_C3( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
              qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
          else
            call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
              qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr )
          endif

        else if( ic_type == 611) then
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          !CALL STF_Beam(ic_type, nn, ecoord, hecMESH%section%sect_R_item(ihead+1:), &
            !     &   material%variables(M_YOUNGS), material%variables(M_POISSON), stiffness(1:nn*ndof,1:nn*ndof))

        else if( ic_type == 641 ) then
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          call UpdateST_Beam_641(ic_type, nn, ecoord, total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &    fstrSOLID%elements(icel)%gausses(:), hecMESH%section%sect_R_item(ihead+1:), qf(1:nn*ndof))

        else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call UpdateST_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:nn), total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &              fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof), thick, 0)

        else if( ic_type == 761 ) then   !for shell-solid mixed analysis
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call UpdateST_Shell_MITC33(731, 3, 6, ecoord(1:3, 1:3), total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &              fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof), thick, 2)

        else if( ic_type == 781 ) then   !for shell-solid mixed analysis
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          isect = hecMESH%section_ID(icel)
          ihead = hecMESH%section%sect_R_index(isect-1)
          thick = hecMESH%section%sect_R_item(ihead+1)
          call UpdateST_Shell_MITC33(741, 4, 6, ecoord(1:3, 1:4), total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &              fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof), thick, 1)

        else if ( ic_type == 3414 ) then
          if(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype /= INCOMP_NEWTONIAN) then
            write(*, *) '###ERROR### : This element is not supported for this material'
            write(*, *) 'ic_type = ', ic_type, ', mtype = ', fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype
            call hecmw_abort(hecmw_comm_get_comm())
          endif
          call UPDATE_C3_vp                                                       &
            ( ic_type, nn, ecoord(:,1:nn), total_disp(1:4,1:nn), du(1:4,1:nn), &
            fstrSOLID%elements(icel)%gausses(:) )
          qf = 0.0d0

          !      else if ( ic_type==731) then
          !        call UPDATE_S3(xx,yy,zz,ee,pp,thick,local_stf)
          !        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
          !      else if ( ic_type==741) then
          !        call UPDATE_S4(xx,yy,zz,ee,pp,thick,local_stf)
          !        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)

        else
          write(*, *) '###ERROR### : Element type not supported for nonlinear static analysis'
          write(*, *) ' ic_type = ', ic_type
          call hecmw_abort(hecmw_comm_get_comm())

        endif

        ! ----- calculate the global internal force ( Q(u_{n+1}^{k-1}) )
        do j = 1, nn
          do i = 1, ndof
            !$omp atomic
            fstrSOLID%QFORCE(ndof*(nodLOCAL(j)-1)+i) = fstrSOLID%QFORCE(ndof*(nodLOCAL(j)-1)+i)+qf(ndof*(j-1)+i)
          enddo
        enddo

        ! ----- calculate strain energy
        if(present(strainEnergy))then
          ndim = getSpaceDimension( fstrSOLID%elements(icel)%etype )
          do j = 1, nn
            do i = 1, ndim
              tmp = 0.5d0*( fstrSOLID%elements(icel)%equiForces(ndim*(j-1)+i)+qf(ndim*(j-1)+i) )*ddu(i,j)
              !$omp atomic
              strainEnergy = strainEnergy+tmp
              fstrSOLID%elements(icel)%equiForces(ndim*(j-1)+i) = qf(ndim*(j-1)+i)
            enddo
          enddo
        endif

      enddo ! icel
      !$omp end do
      !$omp end parallel
    enddo   ! itype

    !C
    !C Update for fstrSOLID%QFORCE
    !C
    if( ndof==3 ) then
      call hecmw_update_3_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node)
    else if( ndof==2 ) then
      call hecmw_update_2_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node)
    else if( ndof==4 ) then
      call hecmw_update_4_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node)
    else if( ndof==6 ) then
      call hecmw_update_m_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node,6)
    endif

  end subroutine fstr_UpdateNewton


  !> Update elastiplastic status
  subroutine fstr_UpdateState( hecMESH, fstrSOLID, tincr)
    use m_fstr
    use m_static_lib
    use m_ElastoPlastic
    use mCreep
    use mViscoElastic
    type(hecmwST_local_mesh) :: hecMESH     !< hecmw mesh
    type(fstr_solid) :: fstrSOLID   !< fstr_solid
    real(kind=kreal) :: tincr
    integer(kind=kint) :: itype, is, iE, ic_type, icel, ngauss, i

    if( associated( fstrSOLID%temperature ) ) then
      do i = 1, hecMESH%n_node
        fstrSOLID%last_temp(i) = fstrSOLID%temperature(i)
      end do
    endif

    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type= hecMESH%elem_type_item(itype)
      if( ic_type == 301 ) ic_type = 111
      if( hecmw_is_etype_link(ic_type) ) cycle

      ngauss = NumOfQuadPoints( ic_type )
      do icel = is, iE
        if( isElastoplastic( fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype ) ) then
          do i = 1, ngauss
            call updateEPState( fstrSOLID%elements(icel)%gausses(i) )
          enddo
        elseif( fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype == NORTON ) then
          if( tincr>0.d0 ) then
            do i = 1, ngauss
              call updateViscoState( fstrSOLID%elements(icel)%gausses(i) )
            enddo
          endif
        elseif( isViscoelastic( fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype ) ) then
          if( tincr > 0.d0 ) then
            do i = 1, ngauss
              call updateViscoElasticState( fstrSOLID%elements(icel)%gausses(i) )
            enddo
          endif
        endif

        do i = 1, ngauss
          fstrSOLID%elements(icel)%gausses(i)%strain_bak = fstrSOLID%elements(icel)%gausses(i)%strain
          fstrSOLID%elements(icel)%gausses(i)%stress_bak = fstrSOLID%elements(icel)%gausses(i)%stress
        enddo
      enddo
    enddo
  end subroutine fstr_UpdateState

  subroutine Update_abort( ic_type, flag )
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

end module m_fstr_Update
