!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides function to calculate to do updates
module m_fstr_Update
  use m_fstr
  implicit none

  private :: Update_abort

contains

  !=====================================================================*
  !  UPDATE_C3
  !>  \brief Update displacement, stress, strain and internal forces
  !>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
  !>  \version    0.00
  !> \if AS
  !>    \par Subroutine Structure
  !>    -# Update displacement      \f$ u_{n+1}^{(k)} = u_{n+1}^{(k-1)} + \delta u^{(k)} \f$
  !>    -# Update stress and strain \f$ \varepsilon_{n+1}^{(k)} = \varepsilon_{n+1}^{(k-1)} + \delta \varepsilon^{(k)} \f$, \f$ \sigma_{n+1}^{(k)} = \sigma_{n+1}^{(k-1)} + \delta \sigma^{(k)} \f$
  !>    -# Upcate internal (equivalent nodal) force  \f$ Q_{n+1}^{(k-1)} ( u_{n+1}^{(k-1)} ) \f$
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

    integer(kind=kint) :: nodLOCAL(fstrSOLID%max_ncon)
    real(kind=kreal)   :: ecoord(3, fstrSOLID%max_ncon)
    real(kind=kreal)   :: thick
    integer(kind=kint) :: ndof, itype, is, iE, ic_type, nn, icel, iiS, i, j

    real(kind=kreal)   :: total_disp(6, fstrSOLID%max_ncon), du(6, fstrSOLID%max_ncon), ddu(6, fstrSOLID%max_ncon)
    real(kind=kreal)   :: tt(fstrSOLID%max_ncon), tt0(fstrSOLID%max_ncon), ttn(fstrSOLID%max_ncon)
    real(kind=kreal)   :: qf(fstrSOLID%max_ncon*6), coords(3, 3)
    integer            :: isect, ihead, cdsys_ID
    integer            :: ndim, initt

    real(kind=kreal), optional :: strainEnergy
    real(kind=kreal) :: tmp
    real(kind=kreal)   :: ddaux(3,3)

    ndof = hecMAT%NDOF
    fstrSOLID%QFORCE=0.0d0

    tt0 = 0.d0
    ttn = 0.d0
    tt = 0.d0

    ! if initial temperature exists
    initt = 0
    if( associated(g_InitialCnd) ) then
        do j=1,size(g_InitialCnd)
          if( g_InitialCnd(j)%cond_name=="temperature" ) then
            initt=j
            exit
          endif
        end do
    endif

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
      if (hecmw_is_etype_patch(ic_type)) cycle

      !element loop
      !$omp parallel default(none), &
        !$omp&  private(icel,iiS,j,nn,nodLOCAL,i,ecoord,ddu,du,total_disp, &
        !$omp&  cdsys_ID,coords,thick,qf,isect,ihead,tmp,ndim,ddaux), &
        !$omp&  shared(iS,iE,hecMESH,fstrSOLID,ndof,hecMAT,ic_type,fstrPR, &
        !$omp&         strainEnergy,iter,time,tincr,initt,g_InitialCnd), &
        !$omp&  firstprivate(tt0,ttn,tt)
      !$omp do
      do icel = is, iE

        ! ----- nodal coordinate, displacement and temperature
        iiS = hecMESH%elem_node_index(icel-1)
        nn = hecMESH%elem_node_index(icel)-iiS
        !if( nn>150 ) stop "elemental nodes > 150!"

        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (iiS+j)
          do i = 1, 3
            ecoord(i,j) = hecMESH%node(3*nodLOCAL(j)+i-3)
          enddo
          do i = 1, ndof
            ddu(i,j) = hecMAT%X(ndof*nodLOCAL(j)+i-ndof)
            du(i,j)  = fstrSOLID%dunode(ndof*nodLOCAL(j)+i-ndof)
            total_disp(i,j) = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof)
          enddo

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
            tt(j)  = fstrSOLID%temperature( nodLOCAL(j) )
          endif
        enddo

        isect = hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        thick = hecMESH%section%sect_R_item(ihead+1)
        cdsys_ID = hecMESH%section%sect_orien_ID(isect)
        if( cdsys_ID > 0 ) call get_coordsys(cdsys_ID, hecMESH, fstrSOLID, coords)

        ! ===== calculate the Internal Force
        if( getSpaceDimension( ic_type ) == 2 ) thick = 1.0d0

        if( ic_type == 241 .or. ic_type == 242 .or. ic_type == 231 .or. ic_type == 232 .or. ic_type == 2322 ) then
          call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:), &
            thick,fstrSOLID%elements(icel)%iset,                             &
            total_disp(1:2,1:nn), ddu(1:2,1:nn), qf(1:nn*ndof),              &
            tt(1:nn), tt0(1:nn), ttn(1:nn)  )

        else if( ic_type == 301 ) then
          call UPDATE_C1( ic_type,nn,ecoord(:,1:nn), thick, total_disp(1:3,1:nn), du(1:3,1:nn), &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:) )

        else if( ic_type == 361 ) then
          if( fstrSOLID%sections(isect)%elemopt361 == kel361FI ) then ! full integration element
            call UPDATE_C3( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
              qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361BBAR ) then ! B-bar element
            call UPDATE_C3D8Bbar( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords,    &
              qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361IC ) then ! incompatible element
            call UPDATE_C3D8IC( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), ddu(1:3,1:nn), cdsys_ID, coords,&
              qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, &
              fstrSOLID%elements(icel)%aux, ddaux(1:3,1:3), tt(1:nn), tt0(1:nn), ttn(1:nn) )
            fstrSOLID%elements(icel)%aux(1:3,1:3) = fstrSOLID%elements(icel)%aux(1:3,1:3) + ddaux(1:3,1:3)
          else if( fstrSOLID%sections(isect)%elemopt361 == kel361FBAR ) then ! F-bar element
            call UPDATE_C3D8Fbar( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords,    &
              qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
          endif

        else if (ic_type == 341 .or. ic_type == 351 .or. ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then
          if( ic_type==341 .and. fstrSOLID%sections(isect)%elemopt341 == kel341SESNS ) cycle ! skip smoothed fem
          call UPDATE_C3( ic_type, nn, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
            qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), iter, time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )

        else if( ic_type == 511) then
          call UPDATE_CONNECTOR( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), tincr )
  
        else if( ic_type == 611) then
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          CALL UpdateST_Beam(ic_type, nn, ecoord, total_disp(1:6,1:nn), du(1:6,1:nn), &
                 &   hecMESH%section%sect_R_item(ihead+1:), fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof))

        else if( ic_type == 641 ) then
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          call UpdateST_Beam_641(ic_type, nn, ecoord, total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &    fstrSOLID%elements(icel)%gausses(:), hecMESH%section%sect_R_item(ihead+1:), qf(1:nn*ndof))

        else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          call UpdateST_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:nn), total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &              fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof), thick, 0)

        else if( ic_type == 761 ) then   !for shell-solid mixed analysis
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          call UpdateST_Shell_MITC33(731, 3, 6, ecoord(1:3, 1:3), total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &              fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof), thick, 2)

        else if( ic_type == 781 ) then   !for shell-solid mixed analysis
          if( fstrPR%nlgeom ) call Update_abort( ic_type, 2 )
          call UpdateST_Shell_MITC33(741, 4, 6, ecoord(1:3, 1:4), total_disp(1:ndof,1:nn), du(1:ndof,1:nn), &
            &              fstrSOLID%elements(icel)%gausses(:), qf(1:nn*ndof), thick, 1)

        else if ( ic_type == 3414 ) then
          if(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype /= INCOMP_NEWTONIAN) &
            &  call Update_abort( ic_type, 3, fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype )
          call UPDATE_C3_vp                                                       &
            ( ic_type, nn, ecoord(:,1:nn), total_disp(1:4,1:nn), du(1:4,1:nn), &
            fstrSOLID%elements(icel)%gausses(:) )
          qf = 0.0d0

        else if ( ic_type == 881 .or. ic_type == 891 ) then  !for selective es/ns smoothed fem
          call UPDATE_C3_SESNS( ic_type, nn, nodLOCAL, ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), cdsys_ID, coords, &
            qf(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:), time, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )

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
    call hecmw_update_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node, ndof)
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
      if( hecmw_is_etype_patch(ic_type) ) cycle

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

  subroutine Update_abort( ic_type, flag, mtype )
    integer(kind=kint), intent(in) :: ic_type
    integer(kind=kint), intent(in) :: flag
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

end module m_fstr_Update
