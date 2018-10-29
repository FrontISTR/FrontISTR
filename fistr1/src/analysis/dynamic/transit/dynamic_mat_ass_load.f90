!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains function to set boundary condition of external load in dynamic analysis

module m_dynamic_mat_ass_load
contains

  !C
  !C***
  !> This function sets boundary condition of external load
  !C***
  !C
  subroutine DYNAMIC_MAT_ASS_LOAD(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC, fstrPARAM, iter )

    use m_fstr
    use m_static_lib
    use m_fstr_precheck
    use m_table_dyn
    use m_common_struct
    use m_utilities

    implicit none
    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    type(fstr_dynamic)       :: fstrDYNAMIC
    type(fstr_param)         :: fstrPARAM

    real(kind=kreal)   :: xx(20), yy(20), zz(20)
    real(kind=kreal)   :: params(0:6)
    real(kind=kreal)   :: vect(60)
    integer(kind=kint) :: iwk(60)
    integer(kind=kint) :: nodLocal(20)
    real(kind=kreal)   :: tt(20), tt0(20), coords(3,3)
    real(kind=kreal),pointer:: temp(:)
    integer(kind=kint) :: ndof, ig0, ig, ityp, ltype, iS0, iE0, ik, in, i, j
    integer(kind=kint) :: icel, ic_type, nn, is, isect, id, iset, nsize
    integer(kind=kint) :: itype, iE, cdsys_ID
    real(kind=kreal)   :: val, rho, thick, pa1
    logical :: fg_surf
    logical, save :: isFirst = .true.

    integer(kind=kint) :: flag_u, ierror
    integer(kind=kint), optional :: iter
    real(kind=kreal) :: f_t

    integer(kind=kint) :: iiS, idofS, idofE
    real(kind=kreal)   :: ecoord(3, 20)
    real(kind=kreal)   :: v(6, 20),  dv(6, 20), r(6*20)
    real(kind=kreal)   :: RHS
    real(kind=kreal)   :: unode_tmp(hecMAT%NDOF*hecMESH%n_node)

    !for torque load
    integer(kind=kint) :: n_rot, rid, n_nodes, idof
    type(tRotInfo)   :: rinfo
    real(kind=kreal) :: tval, normal(3), direc(3), ccoord(3), cdisp(3), cdiff(3)

    ndof = hecMAT%NDOF
    call hecmw_mat_clear_b( hecMAT )
    !C
    !C CLOAD
    !C
    n_rot = fstrSOLID%CLOAD_ngrp_rot
    if( n_rot > 0 ) call fstr_RotInfo_init(n_rot, rinfo)

    do ig0 = 1, fstrSOLID%CLOAD_ngrp_tot
      ig = fstrSOLID%CLOAD_ngrp_ID(ig0)
      ityp = fstrSOLID%CLOAD_ngrp_DOF(ig0)
      val = fstrSOLID%CLOAD_ngrp_val(ig0)

      flag_u = 0
      call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
      val = val*f_t

      iS0= hecMESH%node_group%grp_index(ig-1)+1
      iE0= hecMESH%node_group%grp_index(ig  )

      if( fstrSOLID%CLOAD_ngrp_rotID(ig0) > 0 ) then ! setup torque load information
        rid = fstrSOLID%CLOAD_ngrp_rotID(ig0)
        if( .not. rinfo%conds(rid)%active ) then
          rinfo%conds(rid)%active = .true.
          rinfo%conds(rid)%center_ngrp_id = fstrSOLID%CLOAD_ngrp_centerID(ig0)
          rinfo%conds(rid)%torque_ngrp_id = ig
        endif
        if( ityp>ndof ) ityp = ityp-ndof
        rinfo%conds(rid)%vec(ityp) = val
        cycle
      endif

      do ik = iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        hecMAT%B( ndof*(in-1)+ityp ) = hecMAT%B( ndof*(in-1)+ityp )+val
      enddo
    enddo

    !Add torque load to hecMAT%B
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
        val = dot_product(vect(1:ndof),vect(1:ndof))
        if( val < 1.d-16 ) then
          write(*,*) '###ERROR### : torque node is at the same position as that of center node in rotational surface.'
          call hecmw_abort( hecmw_comm_get_comm() )
        endif
        vect(1:ndof) = (tval/val)*vect(1:ndof)
        hecMAT%B(ndof*(in-1)+1:ndof*in) = hecMAT%B(ndof*(in-1)+1:ndof*in)+vect(1:ndof)
      enddo
    enddo
    if( n_rot > 0 ) call fstr_RotInfo_finalize(rinfo)

    !C
    !C DLOAD
    !C
    do ig0 = 1, fstrSOLID%DLOAD_ngrp_tot
      ig = fstrSOLID%DLOAD_ngrp_ID(ig0)
      ltype = fstrSOLID%DLOAD_ngrp_LID(ig0)
      do i = 0, 6
        params(i) = fstrSOLID%DLOAD_ngrp_params(i,ig0)
      enddo
      !C START & END
      fg_surf = (ltype == 100)
      if( fg_surf ) then ! surface group
        iS0 = hecMESH%surf_group%grp_index(ig-1) + 1
        iE0 = hecMESH%surf_group%grp_index(ig  )
      else ! element group
        iS0 = hecMESH%elem_group%grp_index(ig-1) + 1
        iE0 = hecMESH%elem_group%grp_index(ig  )
      endif

      do ik = iS0, iE0
        if( fg_surf ) then ! surface group
          ltype   = hecMESH%surf_group%grp_item(2*ik)*10
          icel    = hecMESH%surf_group%grp_item(2*ik-1)
          ic_type = hecMESH%elem_type(icel)
        else ! element group
          icel    = hecMESH%elem_group%grp_item(ik)
          ic_type = hecMESH%elem_type(icel)
        endif
        !C** Create local stiffness
        nn = hecmw_get_max_node(ic_type)
        !C** node ID
        is = hecMESH%elem_node_index(icel-1)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (is+j)
          !C** nodal coordinate
          xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
          yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
          zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
          !C** create iwk array ***
          do i = 1, ndof
            iwk(ndof*(j-1)+i) = ndof*(nodLOCAL(j)-1)+i
          enddo
        enddo
        !C** section  ID
        isect = hecMESH%section_ID(icel)
        !C** Get Properties
        rho = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_DENSITY)
        call fstr_get_thickness(hecMESH,isect,thick)

        !C** Section Data
        if( ndof == 2 ) then
          id = hecMESH%section%sect_opt(isect)
          if( id == 0 ) then
            iset = 1
          elseif( id == 1 ) then
            iset = 0
          elseif( id == 2 ) then
            iset = 2
          endif
          pa1 = 1.0
        endif

        !C** Create local stiffness
        if( ic_type == 241 .or.ic_type == 242 .or. ic_type == 231 .or. ic_type == 232 ) then
          call DL_C2(ic_type,nn,xx(1:nn),yy(1:nn),rho,pa1,ltype,params,vect(1:nn*ndof),nsize,iset)

        elseif( ic_type == 341 .or. ic_type == 351 .or. ic_type == 361 .or.   &
            ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then
          call DL_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),rho,ltype,params,vect(1:nn*ndof),nsize)

        elseif( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
          call DL_Shell(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, fstrSOLID%elements(icel)%gausses)
        elseif( ( ic_type==761 ) ) then
          call DL_Shell_33(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, &
            fstrSOLID%elements(icel)%gausses)
        elseif( ( ic_type==781 ) ) then
          call DL_Shell_33(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize, &
            fstrSOLID%elements(icel)%gausses)

        endif
        !
        !!!!!!  time history

        flag_u = 10
        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
        do j=1,nsize
          vect(j) = vect(j)*f_t
        enddo
        !
        !C** Add vector
        do j = 1, nsize
          hecMAT%B( iwk(j) )=hecMAT%B( iwk(j) )+vect(j)
        enddo
      enddo
    enddo

    if ( present(iter) ) then
      if( iter == 1 ) then
        do i = 1, ndof*hecMESH%n_node
          unode_tmp(i) = fstrSOLID%unode(i)
        enddo

        do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
          ig    = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
          RHS   = fstrSOLID%BOUNDARY_ngrp_val(ig0)
          ityp  = fstrSOLID%BOUNDARY_ngrp_type(ig0)
          idofS = ityp/10
          idofE = ityp-idofS*10
          iS0 = hecMESH%node_group%grp_index(ig-1) + 1
          iE0 = hecMESH%node_group%grp_index(ig  )

          do ik = iS0, iE0
            in = hecMESH%node_group%grp_item(ik)
            do idof = idofS, idofE
              unode_tmp( ndof*(in-1)+idof ) = RHS
            enddo
          enddo
        enddo

        do itype = 1, hecMESH%n_elem_type
          ic_type = hecMESH%elem_type_item(itype)
          if( ic_type == 3414 ) then
            nn = hecmw_get_max_node(ic_type)
            if( nn > 20 ) stop "The number of elemental nodes > 20"

            is = hecMESH%elem_type_index(itype-1)+1
            iE = hecMESH%elem_type_index(itype  )
            do icel = is, iE
              if(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype /= INCOMP_NEWTONIAN) then
                write(*, *) '###ERROR### : This element is not supported for this material'
                write(*, *) 'ic_type = ', ic_type, ', mtype = ', fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype
                stop
                call hecmw_abort(hecmw_comm_get_comm())
              endif

              v = 0.0d0
              dv= 0.0d0
              iiS = hecMESH%elem_node_index(icel-1)
              do j = 1, nn
                nodLOCAL(j) = hecMESH%elem_node_item(iiS+j)
                do i = 1, 3
                  ! nodal coordinates
                  ecoord(i,j) = hecMESH%node( 3*nodLOCAL(j)+i-3 )
                  ! nodal velocity
                  v(i,j) = unode_tmp( ndof*nodLOCAL(j)+i-ndof )
                  fstrSOLID%unode( ndof*nodLOCAL(j)+i-ndof ) = v(i,j)
                  ! nodal velocity increment
                  dv(i,j) = fstrSOLID%dunode( ndof*nodLOCAL(j)+i-ndof )
                enddo
              enddo

              call LOAD_C3_vp                                                &
                ( ic_type, nn, ecoord(:,1:nn), v(1:4,1:nn), dv(1:4,1:nn), &
                r(1:nn*ndof), fstrSOLID%elements(icel)%gausses(:),      &
                fstrDYNAMIC%t_delta )

              do j = 1, nn
                do i = 1, ndof
                  hecMAT%B(ndof*(nodLOCAL(j)-1)+i) = hecMAT%B(ndof*(nodLOCAL(j)-1)+i)+r(ndof*(j-1)+i)
                enddo
              enddo
            enddo ! icel
          endif
        enddo ! itype

      else
        do itype = 1, hecMESH%n_elem_type
          ic_type = hecMESH%elem_type_item(itype)
          if( ic_type == 3414 ) then
            nn = hecmw_get_max_node(ic_type)
            if( nn > 20 ) stop "The number of elemental nodes > 20"
            is = hecMESH%elem_type_index(itype-1)+1
            iE = hecMESH%elem_type_index(itype  )
            do icel = is, iE
              iiS = hecMESH%elem_node_index(icel-1)
              do j = 1, nn
                nodLOCAL(j) = hecMESH%elem_node_item(iiS+j)
              enddo
              do j = 1, nn
                do i = 1, ndof
                  hecMAT%B(ndof*(nodLOCAL(j)-1)+i) = 0.0D0
                enddo
              enddo
            enddo
          endif
        enddo
      endif
    endif

    !C
    !C THERMAL LOAD USING TEMPERATURE
    !C
    !C Set Temperature
    !C
    if( fstrSOLID%TEMP_ngrp_tot > 0 ) then

      if( hecMESH%my_rank .eq. 0 ) then
        write(imsg,*) 'stop: THERMAL LOAD is not yet available in dynamic analysis!'
      endif
      call hecmw_abort( hecmw_comm_get_comm())

      allocate ( temp(hecMESH%n_node) )
      temp=0
      do ig0= 1, fstrSOLID%TEMP_ngrp_tot
        ig= fstrSOLID%TEMP_ngrp_ID(ig0)
        val=fstrSOLID%TEMP_ngrp_val(ig0)
        !C START & END
        iS0= hecMESH%node_group%grp_index(ig-1) + 1
        iE0= hecMESH%node_group%grp_index(ig  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          temp( in ) = val
        enddo
      enddo
      !C
      !C +-------------------------------+
      !C | ELEMENT-by-ELEMENT ASSEMBLING |
      !C | according to ELEMENT TYPE     |
      !C +-------------------------------+
      !C===
      do itype = 1, hecMESH%n_elem_type

        is = hecMESH%elem_type_index(itype-1)+1
        iE = hecMESH%elem_type_index(itype  )
        ic_type = hecMESH%elem_type_item(itype)
        if( hecmw_is_etype_link(ic_type) ) cycle
        !C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
        !C element loop
        do icel = is, iE
          !C** node ID
          is= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)=hecMESH%elem_node_item(is+j)
            !C** nodal coordinate
            xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
            zz(j)=hecMESH%node(3*nodLOCAL(j)  )
            tt(j)=temp( nodLOCAL(j) )
            tt0(j)=ref_temp
            !C** create iwk array ***
            do i=1,ndof
              iwk(ndof*(j-1)+i)=ndof*(nodLOCAL(j)-1)+i
            enddo
          enddo

          !C** section  ID
          isect= hecMESH%section_ID(icel)
          cdsys_ID = fstrSOLID%elements(icel)%gausses(1)%pMaterial%cdsys_ID
          call get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )

          !C** Section Data
          if( ndof .eq. 2 ) then
            id=hecMESH%section%sect_opt(isect)
            if( id.eq.0 ) then
              iset=1
            elseif( id.eq.1) then
              iset=0
            elseif( id.eq.2) then
              iset=2
            endif
            pa1=1.0
          endif

          !C** Create local stiffness
          if( ic_type == 241 .or. ic_type == 242 .or. ic_type == 231 .or. ic_type == 232 ) then
            call TLOAD_C2( ic_type, nn, xx(1:nn), yy(1:nn), tt(1:nn), tt0(1:nn),     &
              fstrSOLID%elements(icel)%gausses,pa1,iset, vect(1:nn*2) )

          elseif( ic_type == 361 ) then
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

          elseif (ic_type == 341 .or. ic_type == 351 .or.                       &
              ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then
            call TLOAD_C3                                                                &
              ( ic_type, nn, xx(1:nn), yy(1:nn), zz(1:nn), tt(1:nn), tt0(1:nn),       &
              fstrSOLID%elements(icel)%gausses, vect(1:nn*ndof), cdsys_ID, coords )

          elseif ( ic_type == 741 .or. ic_type == 743 .or. ic_type == 731 ) then
            if( myrank == 0 ) then
              write(IMSG,*) '*------------------------', &
                '-------------------*'
              write(IMSG,*) ' Thermal loading option ', &
                'not yet available.'
              write(IMSG,*) '*------------------------', &
                '-------------------*'
              call hecmw_abort( hecmw_comm_get_comm())
            endif
          endif
          !C** Add vector
          do j = 1, ndof*nn
            hecMAT%B( iwk(j) ) = hecMAT%B( iwk(j) )+vect(j)
          enddo
        enddo
      enddo
      deallocate ( temp )
    endif
  end subroutine DYNAMIC_MAT_ASS_LOAD


end module m_dynamic_mat_ass_load
