!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to take into acount external load

module m_fstr_ass_load
  implicit none
contains
  !
  !======================================================================!
  !> This subroutine assmble following external force into fstrSOLID%GL and hecMAT%B afterwards
  !>  -#  concentrated nodal force
  !>  -#  surface pressure
  !>  -#  volume force
  !>  -#  thermal force

  subroutine fstr_ass_load(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM)
    !======================================================================!
    use m_fstr
    use m_static_lib
    use m_fstr_precheck
    use mMechGauss
    use mReadTemp
    use mULoad
    use m_fstr_spring
    use m_common_struct
    use m_utilities
    integer, intent(in)                  :: cstep !< current step
    type(hecmwST_matrix), intent(inout)  :: hecMAT !< hecmw matrix
    type(hecmwST_local_mesh), intent(in) :: hecMESH !< hecmw mesh
    type(fstr_solid), intent(inout)      :: fstrSOLID !< fstr_solid
    type(fstr_param), intent(inout)      :: fstrPARAM !< analysis control parameters

    real(kind=kreal)   :: xx(20), yy(20), zz(20)
    real(kind=kreal)   :: params(0:6)
    real(kind=kreal)   :: vect(60)
    integer(kind=kint) :: iwk(60)
    integer(kind=kint) :: nodLocal(20)
    real(kind=kreal)   :: tt(20), tt0(20), coords(3, 3), factor
    integer(kind=kint) :: ndof, ig0, ig, ityp, ltype, iS0, iE0, ik, in, i, j
    integer(kind=kint) :: icel, ic_type, nn, is, isect, id, iset, nsize
    integer(kind=kint) :: itype, iE, ierror, grpid, cdsys_ID
    real(kind=kreal)   :: fval, rho, thick, pa1
    logical :: fg_surf
    integer(kind=kint) :: tstep
    type(tMaterial), pointer :: material !< material information
    integer(kind=kint) :: ihead
    real(kind=kreal)   :: a

    !for torque load
    integer(kind=kint) :: n_rot, rid, n_nodes, idof
    type(tRotInfo)     :: rinfo
    real(kind=kreal)   :: tval, normal(3), direc(3), ccoord(3), cdisp(3), cdiff(3)

    ndof = hecMAT%NDOF

    ! -------------------------------------------------------------------
    !  CLOAD
    ! -------------------------------------------------------------------
    n_rot = fstrSOLID%CLOAD_ngrp_rot
    if( n_rot > 0 ) call fstr_RotInfo_init(n_rot, rinfo)

    fstrSOLID%GL(:) = 0.0d0
    fstrSOLID%EFORCE(:) = 0.0d0
    do ig0 = 1, fstrSOLID%CLOAD_ngrp_tot
      grpid = fstrSOLID%CLOAD_ngrp_GRPID(ig0)
      if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
      factor = fstrSOLID%factor(2)
      if( fstr_isLoadActive( fstrSOLID, grpid, cstep-1 ) ) factor = 1.0d0
      ig = fstrSOLID%CLOAD_ngrp_ID(ig0)
      ityp = fstrSOLID%CLOAD_ngrp_DOF(ig0)
      fval = fstrSOLID%CLOAD_ngrp_val(ig0)
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )

      if( fstrSOLID%CLOAD_ngrp_rotID(ig0) > 0 ) then ! setup torque load information
        rid = fstrSOLID%CLOAD_ngrp_rotID(ig0)
        if( .not. rinfo%conds(rid)%active ) then
          rinfo%conds(rid)%active = .true.
          rinfo%conds(rid)%center_ngrp_id = fstrSOLID%CLOAD_ngrp_centerID(ig0)
          rinfo%conds(rid)%torque_ngrp_id = ig
        endif
        if( ityp>ndof ) ityp = ityp-ndof
        rinfo%conds(rid)%vec(ityp) = factor*fval
        cycle
      endif

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        fstrSOLID%GL(ndof*(in-1)+ityp) = fstrSOLID%GL(ndof*(in-1)+ityp)+factor*fval
      enddo
    enddo

    !Add torque load to fstrSOLID%GL
    do rid = 1, n_rot
      if( .not. rinfo%conds(rid)%active ) cycle
      !get number of slave nodes
      n_nodes = hecmw_ngrp_get_number(hecMESH, rinfo%conds(rid)%torque_ngrp_id)

      !get center node
      ig = rinfo%conds(rid)%center_ngrp_id
      do idof = 1, ndof
        ccoord(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, hecMESH%node)
        cdisp(idof) = hecmw_ngrp_get_totalvalue(hecMESH, ig, ndof, idof, fstrSOLID%unode)
      enddo
      ccoord(1:ndof) = ccoord(1:ndof) + cdisp(1:ndof)

      tval = dsqrt(dot_product(rinfo%conds(rid)%vec(1:ndof),rinfo%conds(rid)%vec(1:ndof)))
      if( tval < 1.d-16 ) then
        write(*,*) '###ERROR### : norm of torque vector must be > 0.0'
        call hecmw_abort( hecmw_comm_get_comm() )
      endif
      normal(1:ndof) = rinfo%conds(rid)%vec(1:ndof)/tval
      tval = tval/dble(n_nodes)

      ig = rinfo%conds(rid)%torque_ngrp_id
      iS0 = hecMESH%node_group%grp_index(ig-1) + 1
      iE0 = hecMESH%node_group%grp_index(ig  )
      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        cdiff(1:ndof) = hecMESH%node(ndof*(in-1)+1:ndof*in)+fstrSOLID%unode(ndof*(in-1)+1:ndof*in)-ccoord(1:ndof)
        call cross_product(normal,cdiff,vect(1:ndof))
        fval = dot_product(vect(1:ndof),vect(1:ndof))
        if( fval < 1.d-16 ) then
          write(*,*) '###ERROR### : torque node is at the same position as that of center node in rotational surface.'
          call hecmw_abort( hecmw_comm_get_comm() )
        endif
        vect(1:ndof) = (tval/fval)*vect(1:ndof)
        fstrSOLID%GL(ndof*(in-1)+1:ndof*in) = fstrSOLID%GL(ndof*(in-1)+1:ndof*in)+vect(1:ndof)
      enddo
    enddo
    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)
    !
    ! -------------------------------------------------------------------
    !  DLOAD
    ! -------------------------------------------------------------------
    do ig0 = 1, fstrSOLID%DLOAD_ngrp_tot
      grpid = fstrSOLID%DLOAD_ngrp_GRPID(ig0)
      if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
      factor = fstrSOLID%factor(2)
      if( fstr_isLoadActive( fstrSOLID, grpid, cstep-1 ) ) factor = 1.0d0
      ig = fstrSOLID%DLOAD_ngrp_ID(ig0)
      ltype = fstrSOLID%DLOAD_ngrp_LID(ig0)
      do i = 0, 6
        params(i)= fstrSOLID%DLOAD_ngrp_params(i,ig0)
      enddo
      ! ----- START & END
      fg_surf = (ltype == 100)
      if( fg_surf ) then                  ! surface group
        iS0 = hecMESH%surf_group%grp_index(ig-1) + 1
        iE0 = hecMESH%surf_group%grp_index(ig  )
      else                                ! element group
        iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
        iE0 = hecMESH%elem_group%grp_index(ig  )
      endif
      do ik = iS0, iE0
        if( fg_surf ) then                ! surface group
          ltype   = hecMESH%surf_group%grp_item(2*ik)*10
          icel    = hecMESH%surf_group%grp_item(2*ik-1)
          ic_type = hecMESH%elem_type(icel)
        else                              ! element group
          icel    = hecMESH%elem_group%grp_item(ik)
          ic_type = hecMESH%elem_type(icel)
        endif
        if( hecmw_is_etype_link(ic_type) ) cycle
        ! if( ic_type==3422 ) ic_type=342
        nn = hecmw_get_max_node(ic_type)
        ! ----- node ID
        is = hecMESH%elem_node_index(icel-1)
        if( fstrSOLID%DLOAD_follow == 0 ) then
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item (is+j)
            ! ----- nodal coordinate
            xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
            yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
            zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
            ! ----- create iwk array ***
            do i = 1, ndof
              iwk( ndof*(j-1)+i ) = ndof*( nodLOCAL(j)-1 )+i
            enddo
          enddo
        else
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item (is+j)
            ! ----- nodal coordinate
            if (ndof==2) then
              xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )+fstrSOLID%unode( 2*nodLOCAL(j)-1 )+fstrSOLID%dunode( 2*nodLOCAL(j)-1 )
              yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )+fstrSOLID%unode( 2*nodLOCAL(j)   )+fstrSOLID%dunode( 2*nodLOCAL(j)   )
            else if (ndof==3) then
              xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )+fstrSOLID%unode( 3*nodLOCAL(j)-2 )+fstrSOLID%dunode( 3*nodLOCAL(j)-2 )
              yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )+fstrSOLID%unode( 3*nodLOCAL(j)-1 )+fstrSOLID%dunode( 3*nodLOCAL(j)-1 )
              zz(j) = hecMESH%node( 3*nodLOCAL(j)   )+fstrSOLID%unode( 3*nodLOCAL(j)   )+fstrSOLID%dunode( 3*nodLOCAL(j)   )
            else if (ndof==6) then
              xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )+fstrSOLID%unode( 6*nodLOCAL(j)-5 )+fstrSOLID%dunode( 6*nodLOCAL(j)-5 )
              yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )+fstrSOLID%unode( 6*nodLOCAL(j)-4 )+fstrSOLID%dunode( 6*nodLOCAL(j)-4 )
              zz(j) = hecMESH%node( 3*nodLOCAL(j)   )+fstrSOLID%unode( 6*nodLOCAL(j)-3 )+fstrSOLID%dunode( 6*nodLOCAL(j)-3 )
            endif
            ! ----- create iwk array ***
            do i = 1, ndof
              iwk( ndof*(j-1)+i ) = ndof*( nodLOCAL(j)-1 )+i
            enddo
          enddo
        end if
        ! ----- section  ID
        isect = hecMESH%section_ID(icel)
        ! ----- Get Properties
        material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
        rho = material%variables(M_DENSITY)
        call fstr_get_thickness(hecMESH,isect,thick)
        ! ----- Section Data
        if( ndof == 2 ) then
          id=hecMESH%section%sect_opt(isect)
          if( id == 0 ) then
            iset = 1
          else if( id == 1 ) then
            iset = 0
          else if( id == 2) then
            iset = 2
          endif
          pa1=1.d0
        endif
        ! ----- Create local stiffness
        if (ic_type==301)then
          ihead = hecMESH%section%sect_R_index(isect-1)
          call DL_C1(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),rho,thick,ltype,params,vect(1:nn*ndof),nsize)

        elseif( ic_type == 241 .or. ic_type == 242 .or. ic_type == 231 .or. ic_type == 232 .or. ic_type == 2322 ) then
          call DL_C2(ic_type,nn,xx(1:nn),yy(1:nn),rho,pa1,ltype,params,vect(1:nn*ndof),nsize,iset)

        else if ( ic_type == 341 .or. ic_type == 351 .or. ic_type == 361 .or.   &
            ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then
          call DL_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),rho,ltype,params,vect(1:nn*ndof),nsize)

        else if ( ic_type == 641 ) then
          ihead = hecMESH%section%sect_R_index(isect-1)
          call DL_Beam_641(ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), rho, ltype, params, &
            hecMESH%section%sect_R_item(ihead+1:), vect(1:nn*ndof), nsize)

        else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
          call DL_Shell(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, fstrSOLID%elements(icel)%gausses)

        else if( ( ic_type==761 ) .or. ( ic_type==781 ) ) then
          call DL_Shell_33(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, &
            fstrSOLID%elements(icel)%gausses)

        else
          nsize = 0
          write(*,*)"### WARNING: DLOAD",ic_type

        endif
        ! ----- Add vector
        do j=1,nsize
          fstrSOLID%GL( iwk(j) )=fstrSOLID%GL( iwk(j) )+factor*vect(j)
        enddo
      enddo
    enddo

    ! -----Uload
    call uloading( cstep, factor, fstrSOLID%GL )

    !C
    !C Update for fstrSOLID%GL
    !C
    if( hecMESH%n_dof == 3 ) then
      call hecmw_update_3_R (hecMESH,fstrSOLID%GL,hecMESH%n_node)
    else if( hecMESH%n_dof == 2 ) then
      call hecmw_update_2_R (hecMESH,fstrSOLID%GL,hecMESH%n_node)
    endif

    call hecmw_mat_clear_b( hecMAT )
    do i=1, hecMESH%n_node*  hecMESH%n_dof
      hecMAT%B(i)=fstrSOLID%GL(i)-fstrSOLID%QFORCE(i)
    enddo

    do i=1, hecMAT%NDOF*hecMAT%NP
      !thermal load is not considered
      fstrSOLID%EFORCE(i) = fstrSOLID%GL(i)
    enddo

    ! -------------------------------------------------------------------
    !  TLOAD : THERMAL LOAD USING TEMPERATURE
    ! -------------------------------------------------------------------
    !C
    !C Set Temperature
    !C
    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      do ig0 = 1, fstrSOLID%TEMP_ngrp_tot
        grpid = fstrSOLID%TEMP_ngrp_GRPID(ig0)
        if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
        factor = fstrSOLID%factor(2)
        if( fstr_isLoadActive( fstrSOLID, grpid, cstep-1 ) ) factor = 1.0d0
        ig = fstrSOLID%TEMP_ngrp_ID(ig0)
        fval =fstrSOLID%TEMP_ngrp_val(ig0)
        iS0 = hecMESH%node_group%grp_index(ig-1)+1
        iE0 = hecMESH%node_group%grp_index(ig  )
        do ik = iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          pa1 = fstrSOLID%temp_bak( in )
          fstrSOLID%temperature( in ) = pa1+(fval-pa1)*factor
        enddo
      enddo

      if( fstrSOLID%TEMP_irres > 0 ) then
        call read_temperature_result(hecMESH, fstrSOLID%TEMP_irres, fstrSOLID%TEMP_tstep, &
          &  fstrSOLID%TEMP_interval, fstrSOLID%TEMP_factor, fstrSOLID%temperature, fstrSOLID%temp_bak)
      endif

      ! ----- element TYPE loop.
      do itype = 1, hecMESH%n_elem_type

        is = hecMESH%elem_type_index(itype-1)+1
        iE = hecMESH%elem_type_index(itype  )
        ic_type = hecMESH%elem_type_item(itype)
        if( hecmw_is_etype_link(ic_type) ) cycle
        ! ----- Set number of nodes
        nn = hecmw_get_max_node(ic_type)

        ! ----- element loop
        do icel = is, iE

          ! ----- node ID
          is= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)=hecMESH%elem_node_item(is+j)
            ! ----- nodal coordinate
            if (ndof==2) then
              xx(j)=hecMESH%node(3*nodLOCAL(j)-2)+fstrSOLID%unode(ndof*nodLOCAL(j)-1)
              yy(j)=hecMESH%node(3*nodLOCAL(j)-1)+fstrSOLID%unode(ndof*nodLOCAL(j)  )
            else if (ndof==3) then
              xx(j)=hecMESH%node(3*nodLOCAL(j)-2)+fstrSOLID%unode(ndof*nodLOCAL(j)-2)
              yy(j)=hecMESH%node(3*nodLOCAL(j)-1)+fstrSOLID%unode(ndof*nodLOCAL(j)-1)
              zz(j)=hecMESH%node(3*nodLOCAL(j)  )+fstrSOLID%unode(ndof*nodLOCAL(j))
            endif
            tt0(j)=fstrSOLID%last_temp( nodLOCAL(j) )
            tt(j) = fstrSOLID%temperature( nodLOCAL(j) )
            ! ----- create iwk array ***
            do i=1,ndof
              iwk(ndof*(j-1)+i)=ndof*(nodLOCAL(j)-1)+i
            enddo
          enddo

          ! ----- section  Data
          isect= hecMESH%section_ID(icel)
          cdsys_ID = hecMESH%section%sect_orien_ID(isect)
          call get_coordsys(cdsys_ID, hecMESH, fstrSOLID, coords)

          if( ndof == 2 ) then
            id=hecMESH%section%sect_opt(isect)
            if( id==0 ) then
              iset = 1
            else if( id == 1 ) then
              iset = 0
            else if( id == 2 ) then
              iset = 2
            endif
            pa1 = 1.0d0
          endif

          if( ic_type == 641 ) then

            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)

            call TLOAD_Beam_641( ic_type, nn, ndof, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),    &
              fstrSOLID%elements(icel)%gausses, hecMESH%section%sect_R_item(ihead+1:), &
              vect(1:nn*ndof) )

            do j = 1, ndof*nn
              hecMAT%B( iwk(j) ) = hecMAT%B( iwk(j) )+vect(j)
            end do
            cycle
          end if

          ! ----- Create local stiffness
          if(ic_type == 241 .or. ic_type == 242 .or. ic_type == 231 .or. ic_type == 232 ) then
            call TLOAD_C2( ic_type, nn, xx(1:nn), yy(1:nn), tt(1:nn), tt0(1:nn),      &
              fstrSOLID%elements(icel)%gausses,pa1, iset, vect(1:nn*2) )

          else if( ic_type == 361 ) then
            if( fstrSOLID%sections(isect)%elemopt361 == kel361FI ) then
              call TLOAD_C3                                                          &
                ( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),       &
                fstrSOLID%elements(icel)%gausses, vect(1:nn*ndof), cdsys_ID, coords )
            else if( fstrSOLID%sections(isect)%elemopt361 == kel361BBAR ) then
              call TLOAD_C3D8Bbar                                                          &
                ( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),       &
                fstrSOLID%elements(icel)%gausses, vect(1:nn*ndof), cdsys_ID, coords )
            else if( fstrSOLID%sections(isect)%elemopt361 == kel361IC ) then
              call TLOAD_C3D8IC                                                            &
                ( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),       &
                fstrSOLID%elements(icel)%gausses, vect(1:nn*ndof), cdsys_ID, coords )
            else if( fstrSOLID%sections(isect)%elemopt361 == kel361FBAR ) then
              call TLOAD_C3D8Fbar                                                            &
                ( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),       &
                fstrSOLID%elements(icel)%gausses, vect(1:nn*ndof), cdsys_ID, coords )
            endif

          else if( ic_type == 341 .or. ic_type == 351 .or.                       &
              ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then
            call TLOAD_C3                                                                &
              ( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),       &
              fstrSOLID%elements(icel)%gausses, vect(1:nn*ndof), cdsys_ID, coords )

          else if( ic_type == 741 .or. ic_type == 743 .or. ic_type == 731 ) then
            if( myrank == 0 ) then
              write(IMSG,*) '*------------------------', &
                '-------------------*'
              write(IMSG,*) ' Thermal loading option for shell elements', &
                'not yet available.'
              write(IMSG,*) '*------------------------', &
                '-------------------*'
              call hecmw_abort( hecmw_comm_get_comm())
            endif

          endif

          ! ----- Add vector
          do j = 1, ndof*nn
            ! fstrSOLID%GL( iwk(j) )=fstrSOLID%GL( iwk(j) )+vect(j)
            hecMAT%B( iwk(j) ) = hecMAT%B( iwk(j) )+vect(j)
          enddo

        enddo
      enddo
    endif

    ! ----- Spring force
    call fstr_Update_NDForce_spring( cstep, hecMESH, fstrSOLID, hecMAT%B )

    if( associated( fstrSOLID%contacts ) .and. fstrPARAM%contact_algo == kcaALagrange ) then
      do i = 1, size(fstrSOLID%contacts)
        call ass_contact_force( fstrSOLID%contacts(i), hecMESH%node, fstrSOLID%unode, hecMAT%B )
      enddo
    endif

  end subroutine fstr_ass_load

end module m_fstr_ass_load
