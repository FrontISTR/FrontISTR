!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to calculation nodal stress
module m_fstr_NodalStress
  use m_fstr

  implicit none
  private :: NodalStress_INV3, NodalStress_INV2, inverse_func
contains

  !> Calculate NODAL STRESS of solid elements
  !----------------------------------------------------------------------*
  subroutine fstr_NodalStress3D( hecMESH, fstrSOLID )
    !----------------------------------------------------------------------*
    use m_static_lib
    type(hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)         :: fstrSOLID
    real(kind=kreal), pointer   :: tnstrain(:), testrain(:), yield_ratio(:)
    integer(kind=kint), pointer :: is_rot(:)
    !C** local variables
    integer(kind=kint) :: itype, icel, ic, is, iE, jS, i, j, k, m, ic_type, nn, ni, ID_area
    integer(kind=kint) :: nodlocal(20), ntemp
    integer(kind=kint), allocatable :: nnumber(:)
    real(kind=kreal)   :: estrain(6), estress(6), naturalCoord(3)
    real(kind=kreal)   :: enqm(12)
    real(kind=kreal)   :: ndstrain(20,6), ndstress(20,6), tdstrain(20,6)
    real(kind=kreal)   :: ecoord(3, 20), edisp(60), tt(20), t0(20)
    real(kind=kreal), allocatable :: func(:,:), inv_func(:,:)

    !C** Shell33 variables
    integer(kind=kint) :: isect, ihead, ntot_lyr, nlyr, flag33, cid, truss
    real(kind=kreal)   :: thick, thick_lyr, dtot_lyr
    call fstr_solid_phys_clear(fstrSOLID)

    allocate( nnumber(hecMESH%n_node) )
    if( .not. associated(fstrSOLID%is_rot) ) allocate( fstrSOLID%is_rot(hecMESH%n_node) )
    !allocate( fstrSOLID%yield_ratio(hecMESH%n_elem) )
    nnumber = 0
    fstrSOLID%is_rot = 0
    !fstrSOLID%yield_ratio = 0.0d0

    tnstrain => fstrSOLID%tnstrain
    testrain => fstrSOLID%testrain
    is_rot   => fstrSOLID%is_rot
    yield_ratio => fstrSOLID%yield_ratio

    if( associated(tnstrain) ) tnstrain = 0.0d0

    !C**  setting
    ntot_lyr = fstrSOLID%max_lyr
    flag33   = fstrSOLID%is_33shell
    truss    = fstrSOLID%is_33beam

    !C +-------------------------------+
    !C | according to ELEMENT TYPE     |
    !C +-------------------------------+
    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type == fe_tet10nc ) ic_type = fe_tet10n
      if( .not. (hecmw_is_etype_solid(ic_type) .or. ic_type == 781 &
        & .or. ic_type == 761 .or. ic_type == fe_beam341 ) ) cycle
      !C** set number of nodes and shape function
      nn = hecmw_get_max_node( ic_type )
      ni = NumOfQuadPoints( ic_type )
      allocate( func(ni,nn), inv_func(nn,ni) )
      if( ic_type == fe_tet10n ) then
        ic = hecmw_get_max_node( fe_tet4n )
        do i = 1, ni
          call getQuadPoint( ic_type, i, naturalCoord )
          call getShapeFunc( fe_tet4n, naturalCoord, func(i,1:ic) )
        enddo
        call inverse_func( ic, func, inv_func )
      else if( ic_type == fe_hex8n ) then
        do i = 1, ni
          call getQuadPoint( ic_type, i, naturalCoord )
          call getShapeFunc( ic_type, naturalCoord, func(i,1:nn) )
        enddo
        call inverse_func( ni, func, inv_func )
      else if( ic_type == fe_prism15n ) then
        ic = 0
        do i = 1, ni
          if( i==1 .or. i==2 .or. i==3 .or. i==7 .or. i==8 .or. i==9 ) then
            ic = ic + 1
            call getQuadPoint( ic_type, i, naturalCoord )
            call getShapeFunc( fe_prism6n, naturalCoord, func(ic,1:6) )
          endif
        enddo
        call inverse_func( ic, func, inv_func )
        ni = ic
      else if( ic_type == fe_hex20n ) then
        ic = 0
        do i = 1, ni
          if( i==1 .or. i==3 .or. i==7 .or. i==9 .or. &
              i==19 .or. i==21 .or. i==25 .or. i==27 ) then
            ic = ic + 1
            call getQuadPoint( ic_type, i, naturalCoord )
            call getShapeFunc( fe_hex8n, naturalCoord, func(ic,1:8) )
          endif
        enddo
        call inverse_func( ic, func, inv_func )
        ni = ic
      endif
      !C** element loop
      do icel = is, iE
        jS = hecMESH%elem_node_index(icel-1)
        ID_area = hecMESH%elem_ID(icel*2)
        isect= hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        thick = hecMESH%section%sect_R_item(ihead+1)
        !initialize
        enqm     = 0.0d0
        estrain  = 0.0d0
        estress  = 0.0d0
        ndstrain = 0.0d0
        ndstress = 0.0d0
        !if( ID_area == hecMESH%my_rank ) then

        !--- calculate nodal and elemental value
        if( ic_type == 641 ) then !<3*3 beam section
          do j = 1, 4
            nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
            ecoord(1:3,j)    = hecMESH%node(3*nodLOCAL(j)-2:3*nodLOCAL(j))
            edisp(3*j-2:3*j) = fstrSOLID%unode(3*nodLOCAL(j)-2:3*nodLOCAL(j))
          end do
          ntemp = 0
          if( associated( fstrSOLID%temperature ) ) then
            ntemp = 1
            do j = 1, 4
              nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
              t0(j) = fstrSOLID%last_temp( nodLOCAL(j) )
              tt(j) = fstrSOLID%temperature( nodLOCAL(j) )
            end do
          end if
          call NodalStress_Beam_641( ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses, &
            &     hecMESH%section%sect_R_item(ihead+1:), edisp,                               &
            &     ndstrain(1:nn,1:6), ndstress(1:nn,1:6), tt(1:nn), t0(1:nn), ntemp )
          call ElementalStress_Beam_641( fstrSOLID%elements(icel)%gausses, estrain, estress, enqm )
          fstrSOLID%ENQM(icel*12-11:icel*12) = enqm(1:12)


        elseif( ic_type == 781) then !<3*3 shell section
          do j = 1, 4
            nodLOCAL(j  ) = hecMESH%elem_node_item(jS+j  )
            nodLOCAL(j+4) = hecMESH%elem_node_item(jS+j+4)
            is_rot(nodLOCAL(j+4)) = 1
            ecoord(1:3,j  ) = hecMESH%node(3*nodLOCAL(j  )-2:3*nodLOCAL(j  ))
            ecoord(1:3,j+4) = hecMESH%node(3*nodLOCAL(j+4)-2:3*nodLOCAL(j+4))
            edisp(6*j-5:6*j-3) = fstrSOLID%unode(3*nodLOCAL(j  )-2:3*nodLOCAL(j  ))
            edisp(6*j-2:6*j  ) = fstrSOLID%unode(3*nodLOCAL(j+4)-2:3*nodLOCAL(j+4))
          enddo
          ntot_lyr = fstrSOLID%elements(icel)%gausses(1)%pMaterial%totallyr
          do nlyr=1,ntot_lyr
            call ElementStress_Shell_MITC( 741, 4, 6, ecoord, fstrSOLID%elements(icel)%gausses, edisp, &
              & ndstrain(1:4,1:6), ndstress(1:4,1:6), thick, 1.0d0, nlyr, ntot_lyr)
            call fstr_Stress_add_shelllyr(4,fstrSOLID,icel,nodLOCAL,nlyr,ndstrain(1:4,1:6),ndstress(1:4,1:6),1)
            !minus section
            call ElementStress_Shell_MITC( 741, 4, 6, ecoord, fstrSOLID%elements(icel)%gausses, edisp, &
              & ndstrain(1:4,1:6), ndstress(1:4,1:6), thick,-1.0d0, nlyr, ntot_lyr)
            call fstr_Stress_add_shelllyr(4,fstrSOLID,icel,nodLOCAL,nlyr,ndstrain(1:4,1:6),ndstress(1:4,1:6),-1)
          enddo
          call fstr_getavg_shell(4,fstrSOLID,icel,nodLOCAL,ndstrain(1:4,1:6),ndstress(1:4,1:6),estrain,estress)

        elseif( ic_type == 761) then !<3*3 shell section
          do j = 1, 3
            nodLOCAL(j  ) = hecMESH%elem_node_item(jS+j  )
            nodLOCAL(j+3) = hecMESH%elem_node_item(jS+j+3)
            is_rot(nodLOCAL(j+3)) = 1
            ecoord(1:3,j  ) = hecMESH%node(3*nodLOCAL(j  )-2:3*nodLOCAL(j  ))
            ecoord(1:3,j+3) = hecMESH%node(3*nodLOCAL(j+3)-2:3*nodLOCAL(j+3))
            edisp(6*j-5:6*j-3) = fstrSOLID%unode(3*nodLOCAL(j  )-2:3*nodLOCAL(j  ))
            edisp(6*j-2:6*j  ) = fstrSOLID%unode(3*nodLOCAL(j+3)-2:3*nodLOCAL(j+3))
          enddo
          ntot_lyr = fstrSOLID%elements(icel)%gausses(1)%pMaterial%totallyr
          do nlyr=1,ntot_lyr
            call ElementStress_Shell_MITC( 731, 3, 6, ecoord, fstrSOLID%elements(icel)%gausses, edisp, &
              & ndstrain(1:3,1:6), ndstress(1:3,1:6), thick, 1.0d0, nlyr, ntot_lyr)
            call fstr_Stress_add_shelllyr(3,fstrSOLID,icel,nodLOCAL,nlyr,ndstrain(1:3,1:6),ndstress(1:3,1:6),1)
            !minus section
            call ElementStress_Shell_MITC( 731, 3, 6, ecoord, fstrSOLID%elements(icel)%gausses, edisp, &
              & ndstrain(1:3,1:6), ndstress(1:3,1:6), thick,-1.0d0, nlyr, ntot_lyr)
            call fstr_Stress_add_shelllyr(3,fstrSOLID,icel,nodLOCAL,nlyr,ndstrain(1:3,1:6),ndstress(1:3,1:6),-1)
          enddo
          call fstr_getavg_shell(3,fstrSOLID,icel,nodLOCAL,ndstrain(1:3,1:6),ndstress(1:3,1:6),estrain,estress)

        else if( ic_type == 301 ) then
          call NodalStress_C1( ic_type, nn, fstrSOLID%elements(icel)%gausses, &
            ndstrain(1:nn,1:6), ndstress(1:nn,1:6) )
          call ElementStress_C1( ic_type, fstrSOLID%elements(icel)%gausses, estrain, estress )

        else if( ic_type == fe_tet10n .or. ic_type == fe_hex8n .or. &
            ic_type == fe_prism15n .or. ic_type == fe_hex20n ) then
          call NodalStress_INV3( ic_type, ni, fstrSOLID%elements(icel)%gausses, &
            inv_func, ndstrain(1:nn,1:6), ndstress(1:nn,1:6), &
            tdstrain(1:nn,1:6) )
          call ElementStress_C3( ic_type, fstrSOLID%elements(icel)%gausses, estrain, estress )

        else if ( ic_type == 881 .or. ic_type == 891 ) then  !for selective es/ns smoothed fem
          cycle
        else
          if( ic_type == 341 .and. fstrSOLID%sections(isect)%elemopt341 == kel341SESNS ) cycle

          call NodalStress_C3( ic_type, nn, fstrSOLID%elements(icel)%gausses, &
            ndstrain(1:nn,1:6), ndstress(1:nn,1:6) )
          !call NodalStress_C3( ic_type, nn, fstrSOLID%elements(icel)%gausses, &
            !                     ndstrain(1:nn,1:6), ndstress(1:nn,1:6), tdstrain(1:nn,1:6) )
          call ElementStress_C3( ic_type, fstrSOLID%elements(icel)%gausses, estrain, estress )

        endif

        !ADD VALUE and Count node
        do j = 1, nn
          ic = hecMESH%elem_node_item(jS+j)
          fstrSOLID%STRAIN(6*(ic-1)+1:6*(ic-1)+6)  = fstrSOLID%STRAIN(6*(ic-1)+1:6*(ic-1)+6)  + ndstrain(j,1:6)
          fstrSOLID%STRESS(6*(ic-1)+1:6*(ic-1)+6)  = fstrSOLID%STRESS(6*(ic-1)+1:6*(ic-1)+6)  + ndstress(j,1:6)
          if( associated(tnstrain) )then
            tnstrain(6*(ic-1)+1:6*(ic-1)+6) = tnstrain(6*(ic-1)+1:6*(ic-1)+6) + tdstrain(j,1:6)
          endif
          nnumber(ic) = nnumber(ic) + 1
        enddo

        fstrSOLID%ESTRAIN(6*(icel-1)+1:6*(icel-1)+6) = fstrSOLID%ESTRAIN(6*(icel-1)+1:6*(icel-1)+6) + estrain(1:6)
        fstrSOLID%ESTRESS(6*(icel-1)+1:6*(icel-1)+6) = fstrSOLID%ESTRESS(6*(icel-1)+1:6*(icel-1)+6) + estress(1:6)

        !endif
      enddo !<element loop
      deallocate( func, inv_func )
    enddo !<element type loop

    !C** calculate nodal stress and strain
    do i = 1, hecMESH%n_node
      if( nnumber(i) == 0 ) cycle
      fstrSOLID%STRAIN(6*(i-1)+1:6*(i-1)+6) = fstrSOLID%STRAIN(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
      fstrSOLID%STRESS(6*(i-1)+1:6*(i-1)+6) = fstrSOLID%STRESS(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
      if( associated(tnstrain) )then
        tnstrain(6*(i-1)+1:6*(i-1)+6) = tnstrain(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
      endif
    enddo

    if( fstrSOLID%is_smoothing_active ) call fstr_NodalStress3D_C3D4_SESNS( &
      &  hecMESH, fstrSOLID, nnumber, fstrSOLID%STRAIN, fstrSOLID%STRESS, fstrSOLID%ESTRAIN, fstrSOLID%ESTRESS )

    if( flag33 == 1 )then
      do nlyr = 1, ntot_lyr
        do i = 1, hecMESH%n_node
          if( nnumber(i) == 0 ) cycle
          fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRAIN(6*(i-1)+1:6*(i-1)+6)  = &
            & fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRAIN(6*(i-1)+1:6*(i-1)+6)  / nnumber(i)
          fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRESS(6*(i-1)+1:6*(i-1)+6)  = &
            & fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRESS(6*(i-1)+1:6*(i-1)+6)  / nnumber(i)
          fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRAIN(6*(i-1)+1:6*(i-1)+6) = &
            & fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRAIN(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
          fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRESS(6*(i-1)+1:6*(i-1)+6) = &
            & fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRESS(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
        enddo
      enddo
    endif

    !C** calculate von MISES stress
    do i = 1, hecMESH%n_node
      fstrSOLID%MISES(i) = get_mises(fstrSOLID%STRESS(6*(i-1)+1:6*(i-1)+6))
    enddo
    do i = 1, hecMESH%n_elem
      fstrSOLID%EMISES(i) = get_mises(fstrSOLID%ESTRESS(6*(i-1)+1:6*(i-1)+6))
    enddo

    if( flag33 == 1 )then
      if( fstrSOLID%output_ctrl(3)%outinfo%on(27) .or. fstrSOLID%output_ctrl(4)%outinfo%on(27) ) then
        do nlyr = 1, ntot_lyr
          call make_principal(fstrSOLID, hecMESH, fstrSOLID%SHELL%LAYER(nlyr)%PLUS)
          call make_principal(fstrSOLID, hecMESH, fstrSOLID%SHELL%LAYER(nlyr)%MINUS)
        enddo
      endif
      call make_principal(fstrSOLID, hecMESH, fstrSOLID%SHELL)
    else
      call make_principal(fstrSOLID, hecMESH, fstrSOLID%SOLID)
    endif

    deallocate( nnumber )

  end subroutine fstr_NodalStress3D

  integer(kind=kint) function search_idx_SENES( irow, asect, nid, sid )
    integer(kind=kint), allocatable, intent(in) :: irow(:)
    integer(kind=kint), allocatable, intent(in) :: asect(:)
    integer(kind=kint), intent(in) :: nid
    integer(kind=kint), intent(in) :: sid

    integer(kind=kint) :: i

    search_idx_SENES = -1
    do i=irow(nid-1)+1,irow(nid)
      if( asect(i) == sid ) then
        search_idx_SENES = i
        return
      end if
    end do

  end function

  subroutine fstr_NodalStress3D_C3D4_SESNS( hecMESH, fstrSOLID, nnumber, &
      Nodal_STRAIN, Nodal_STRESS, Elemental_STRAIN, Elemental_STRESS )
    type(hecmwST_local_mesh),intent(in) :: hecMESH
    type(fstr_solid),intent(inout)         :: fstrSOLID
    integer(kind=kint), allocatable, intent(inout) :: nnumber(:)
    real(kind=kreal), pointer, intent(inout)   :: Nodal_STRAIN(:)
    real(kind=kreal), pointer, intent(inout)   :: Nodal_STRESS(:)
    real(kind=kreal), pointer, intent(inout)   :: Elemental_STRAIN(:)
    real(kind=kreal), pointer, intent(inout)   :: Elemental_STRESS(:)

    integer(kind=kint) :: itype, iS, iE, jS, ic_type, icel, i, j, isect
    integer(kind=kint) :: nsize, nid(2), idx(2), nd
    integer(kind=kint) :: nnode, nlen
    type(hecmwST_varray_int), allocatable :: nodal_sections(:)
    real(kind=kreal)   :: tmpval(6), hydval, nsecdup
    integer(kind=kint), allocatable :: irow(:), jcol(:), asect(:)
    real(kind=kreal), allocatable :: stress_hyd(:), strain_hyd(:)
    real(kind=kreal), allocatable :: stress_dev(:)
    real(kind=kreal) :: stress_hyd_ndave(6), strain_hyd_ndave(6)
    real(kind=kreal) :: stress_dev_ndave(6), strain_dev_ndave(6)
    real(kind=kreal), allocatable :: n_dup_dev(:), n_dup_hyd(:)
    real(kind=kreal)   :: edstrain(6), edstress(6)

    nnode = hecMESH%n_node
    nsize = size(Nodal_STRAIN)

    ! create section info at node
    call HECMW_varray_int_initialize_all( nodal_sections, nnode, 2 )
    do itype = 1, hecMESH%n_elem_type
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type /= 341 ) cycle

      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )

      do icel=iS,iE
        isect= hecMESH%section_ID(icel)
        if( fstrSOLID%sections(isect)%elemopt341 /= kel341SESNS ) cycle
        jS = hecMESH%elem_node_index(icel-1)
        do i=1,4
          nd = hecMESH%elem_node_item(jS+i)
          call HECMW_varray_int_add_if_not_exits( nodal_sections(nd), isect )
        end do
      end do
    enddo

    ! create CRS arrays of nodal stress/strain with different sections
    allocate(irow(0:nnode))
    irow(0) = 0
    do i=1,nnode
      irow(i) = irow(i-1)+HECMW_varray_int_get_nitem(nodal_sections(i))
    end do
    nlen = irow(nnode)

    allocate(asect(nlen))
    do i=1,nnode
      if( irow(i-1) == irow(i) ) cycle
      call HECMW_varray_int_get_item_all( nodal_sections(i), asect(irow(i-1)+1:irow(i)) )
    end do

    ! add stress/strain from smoothed elements
    allocate(stress_hyd(6*nlen), strain_hyd(6*nlen))
    allocate(stress_dev(6*nlen))
    allocate(n_dup_dev(nlen),n_dup_hyd(nlen))

    stress_hyd(:) = 0.d0
    strain_hyd(:) = 0.d0
    stress_dev(:) = 0.d0
    n_dup_hyd(:) = 0.d0
    n_dup_dev(:) = 0.d0
    do itype = 1, hecMESH%n_elem_type
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type /= 881 .and. ic_type /= 891 ) cycle

      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )

      do icel=iS,iE
        jS = hecMESH%elem_node_index(icel-1)
        isect= hecMESH%section_ID(icel)
        if( ic_type == 881 ) then
          nid(1) = hecMESH%elem_node_item(jS+1)
          idx(1) = search_idx_SENES( irow, asect, nid(1), isect )

          !strain
          strain_hyd(6*idx(1)-5:6*idx(1)) = fstrSOLID%elements(icel)%gausses(1)%strain_out(1:6)
          !stress
          stress_hyd(6*idx(1)-5:6*idx(1)) =  fstrSOLID%elements(icel)%gausses(1)%stress_out(1:6)
          !number of duplication
          n_dup_hyd(idx(1)) = n_dup_hyd(idx(1)) + 1.d0
        else if( ic_type == 891 ) then
          nid(1:2) = hecMESH%elem_node_item(jS+1:jS+2)
          idx(1) = search_idx_SENES( irow, asect, nid(1), isect )
          idx(2) = search_idx_SENES( irow, asect, nid(2), isect )

          !stress
          tmpval(1:6) = fstrSOLID%elements(icel)%gausses(1)%stress_out(1:6)
          stress_dev(6*idx(1)-5:6*idx(1)) = stress_dev(6*idx(1)-5:6*idx(1)) + tmpval(1:6)
          stress_dev(6*idx(2)-5:6*idx(2)) = stress_dev(6*idx(2)-5:6*idx(2)) + tmpval(1:6)
          !number of duplication
          n_dup_dev(idx(1)) = n_dup_dev(idx(1)) + 1.d0
          n_dup_dev(idx(2)) = n_dup_dev(idx(2)) + 1.d0
        end if
      end do
    enddo

    do i=1,nnode
      if( irow(i-1) == irow(i) ) cycle
      do j=irow(i-1)+1,irow(i)
        if( n_dup_dev(j) < 1.0d-8 ) cycle
        stress_dev(6*j-5:6*j) = stress_dev(6*j-5:6*j)/n_dup_dev(j)
      end do
    end do

    ! average at node for nodal output
    do i=1,nnode
      if( irow(i-1) == irow(i) ) cycle
      strain_hyd_ndave(:) = 0.d0
      stress_hyd_ndave(:) = 0.d0
      stress_dev_ndave(:) = 0.d0
      do j=irow(i-1)+1,irow(i)
        strain_hyd_ndave(1:6) = strain_hyd_ndave(1:6) + strain_hyd(6*j-5:6*j)
        stress_hyd_ndave(1:6) = stress_hyd_ndave(1:6) + stress_hyd(6*j-5:6*j)
        stress_dev_ndave(1:6) = stress_dev_ndave(1:6) + stress_dev(6*j-5:6*j)
      end do
      nsecdup = dble(irow(i)-irow(i-1))
      strain_hyd_ndave(1:6) = strain_hyd_ndave(1:6)/nsecdup
      stress_hyd_ndave(1:6) = stress_hyd_ndave(1:6)/nsecdup
      stress_dev_ndave(1:6) = stress_dev_ndave(1:6)/nsecdup

      if( nnumber(i) == 0 ) then
        Nodal_STRAIN(6*i-5:6*i) = strain_hyd_ndave(1:6)
        Nodal_STRESS(6*i-5:6*i) = stress_hyd_ndave(1:6)+stress_dev_ndave(1:6)
      else
        Nodal_STRAIN(6*i-5:6*i) = 0.5d0*(Nodal_STRAIN(6*i-5:6*i)+strain_hyd_ndave(1:6))
        Nodal_STRESS(6*i-5:6*i) = 0.5d0*(Nodal_STRESS(6*i-5:6*i)+stress_hyd_ndave(1:6)+stress_dev_ndave(1:6))
      endif
    end do

    ! ELEMENTAL STRAIN and STRESS
    do itype = 1, hecMESH%n_elem_type
      ic_type = hecMESH%elem_type_item(itype)
      if( ic_type /= 341 ) cycle

      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )

      do icel=iS,iE
        isect= hecMESH%section_ID(icel)
        if( fstrSOLID%sections(isect)%elemopt341 /= kel341SESNS ) cycle
        jS = hecMESH%elem_node_index(icel-1)
        edstrain(1:6) = 0.d0
        edstress(1:6) = 0.d0
        do i=1,4
          nd = hecMESH%elem_node_item(jS+i)
          idx(1) = search_idx_SENES( irow, asect, hecMESH%elem_node_item(jS+i), isect )
          edstrain(1:6) = edstrain(1:6) + strain_hyd(6*idx(1)-5:6*idx(1))
          edstress(1:6) = edstress(1:6) + stress_hyd(6*idx(1)-5:6*idx(1)) + stress_dev(6*idx(1)-5:6*idx(1))
        end do
        edstrain(1:6) = 0.25d0*edstrain(1:6)
        edstress(1:6) = 0.25d0*edstress(1:6)

        Elemental_STRAIN(6*(icel-1)+1:6*(icel-1)+6) = Elemental_STRAIN(6*(icel-1)+1:6*(icel-1)+6) + edstrain(1:6)
        Elemental_STRESS(6*(icel-1)+1:6*(icel-1)+6) = Elemental_STRESS(6*(icel-1)+1:6*(icel-1)+6) + edstress(1:6)

        fstrSOLID%elements(icel)%gausses(1)%strain_out(1:6) = Elemental_STRAIN(6*(icel-1)+1:6*(icel-1)+6)
        fstrSOLID%elements(icel)%gausses(1)%stress_out(1:6) = Elemental_STRESS(6*(icel-1)+1:6*(icel-1)+6)
      end do
    enddo

    deallocate(stress_hyd, strain_hyd)
    deallocate(stress_dev)
    deallocate(n_dup_dev, n_dup_hyd)

  end subroutine

  subroutine fstr_Stress_add_shelllyr(nn,fstrSOLID,icel,nodLOCAL,nlyr,strain,stress,flag)
    implicit none
    type(fstr_solid)   :: fstrSOLID
    integer(kind=kint) :: nodLOCAL(20)
    integer(kind=kint) :: nn, i, j, k, m, nlyr, weight, icel, flag
    real(kind=kreal)   :: strain(nn, 6), stress(nn, 6)
    type(fstr_solid_physic_val), pointer :: layer => null()

    do j = 1, nn
      i = nodLOCAL(j)
      m = nodLOCAL(j+nn)
      if(flag == 1)then
        layer => fstrSOLID%SHELL%LAYER(nlyr)%PLUS
      elseif(flag == -1)then
        layer => fstrSOLID%SHELL%LAYER(nlyr)%MINUS
      endif
      do k = 1, 6
        layer%STRAIN(6*(i-1)+k)   = layer%STRAIN(6*(i-1)+k)  + strain(j,k)
        layer%STRAIN(6*(m-1)+k)   = layer%STRAIN(6*(m-1)+k)  + strain(j,k)
        layer%STRESS(6*(i-1)+k)   = layer%STRESS(6*(i-1)+k)  + stress(j,k)
        layer%STRESS(6*(m-1)+k)   = layer%STRESS(6*(m-1)+k)  + stress(j,k)
        layer%ESTRAIN(6*(icel-1)+k)  = layer%ESTRAIN(6*(icel-1)+k) + strain(j,k)/nn
        layer%ESTRESS(6*(icel-1)+k)  = layer%ESTRESS(6*(icel-1)+k) + stress(j,k)/nn
      enddo
    enddo
  end subroutine fstr_Stress_add_shelllyr

  subroutine fstr_getavg_shell(nn,fstrSOLID,icel,nodLOCAL,strain,stress,estrain,estress)
    implicit none
    type (fstr_solid)  :: fstrSOLID
    integer(kind=kint) :: nodLOCAL(20)
    integer(kind=kint) :: nn, i, j, k, m, nlyr, icel, flag, ntot_lyr
    real(kind=kreal)   :: strain(nn,6), stress(nn,6), estrain(6), estress(6), weight
    type(fstr_solid_physic_val), pointer :: layer => null()

    ntot_lyr = fstrSOLID%elements(icel)%gausses(1)%pMaterial%totallyr
    strain  = 0.0d0
    stress  = 0.0d0
    estrain = 0.0d0
    estress = 0.0d0

    do nlyr = 1, ntot_lyr
      layer => fstrSOLID%SHELL%LAYER(nlyr)
      weight = fstrSOLID%elements(icel)%gausses(1)%pMaterial%shell_var(nlyr)%weight
      do j = 1, nn
        i = nodLOCAL(j)
        do k = 1, 6
          strain(j,k) = strain(j,k) &
            & + weight*(0.5d0*layer%PLUS%STRAIN(6*(i-1)+k) + 0.5d0*layer%MINUS%STRAIN(6*(i-1)+k))
          stress(j,k) = stress(j,k) &
            & + weight*(0.5d0*layer%PLUS%STRESS(6*(i-1)+k) + 0.5d0*layer%MINUS%STRESS(6*(i-1)+k))
        enddo
        estrain(j) = estrain(j) &
          & + weight*(0.5d0*layer%PLUS%ESTRAIN(6*(icel-1)+j) + 0.5d0*layer%MINUS%ESTRAIN(6*(icel-1)+j))
        estress(j) = estress(j) &
          & + weight*(0.5d0*layer%PLUS%ESTRESS(6*(icel-1)+j) + 0.5d0*layer%MINUS%ESTRESS(6*(icel-1)+j))
      enddo
    enddo
  end subroutine fstr_getavg_shell

  !----------------------------------------------------------------------*
  subroutine NodalStress_INV3( etype, ni, gausses, func, edstrain, edstress, tdstrain )
    !----------------------------------------------------------------------*
    use mMechGauss
    integer(kind=kint) :: etype, ni
    type(tGaussStatus) :: gausses(:)
    real(kind=kreal)   :: func(:, :), edstrain(:, :), edstress(:, :), tdstrain(:, :)
    integer :: i, j, k, ic

    edstrain = 0.0d0
    edstress = 0.0d0
    tdstrain = 0.0d0

    if( etype == fe_hex8n ) then
      do i = 1, ni
        do j = 1, ni
          do k = 1, 6
            edstrain(i,k) = edstrain(i,k) + func(i,j) * gausses(j)%strain_out(k)
            edstress(i,k) = edstress(i,k) + func(i,j) * gausses(j)%stress_out(k)
            !            tdstrain(i,k) = tdstrain(i,k) + func(i,j) * gausses(j)%tstrain(k)
          enddo
        enddo
      enddo
    else if( etype == fe_tet10n ) then
      do i = 1, ni
        do j = 1, ni
          do k = 1, 6
            edstrain(i,k) = edstrain(i,k) + func(i,j) * gausses(j)%strain_out(k)
            edstress(i,k) = edstress(i,k) + func(i,j) * gausses(j)%stress_out(k)
            !            tdstrain(i,k) = tdstrain(i,k) + func(i,j) * gausses(j)%tstrain(k)
          enddo
        enddo
      enddo
      edstrain(5,1:6) = ( edstrain(1,1:6) + edstrain(2,1:6) ) / 2.0
      edstress(5,1:6) = ( edstress(1,1:6) + edstress(2,1:6) ) / 2.0
      tdstrain(5,1:6) = ( tdstrain(1,1:6) + tdstrain(2,1:6) ) / 2.0
      edstrain(6,1:6) = ( edstrain(2,1:6) + edstrain(3,1:6) ) / 2.0
      edstress(6,1:6) = ( edstress(2,1:6) + edstress(3,1:6) ) / 2.0
      tdstrain(6,1:6) = ( tdstrain(2,1:6) + tdstrain(3,1:6) ) / 2.0
      edstrain(7,1:6) = ( edstrain(3,1:6) + edstrain(1,1:6) ) / 2.0
      edstress(7,1:6) = ( edstress(3,1:6) + edstress(1,1:6) ) / 2.0
      tdstrain(7,1:6) = ( tdstrain(3,1:6) + tdstrain(1,1:6) ) / 2.0
      edstrain(8,1:6) = ( edstrain(1,1:6) + edstrain(4,1:6) ) / 2.0
      edstress(8,1:6) = ( edstress(1,1:6) + edstress(4,1:6) ) / 2.0
      tdstrain(8,1:6) = ( tdstrain(1,1:6) + tdstrain(4,1:6) ) / 2.0
      edstrain(9,1:6) = ( edstrain(2,1:6) + edstrain(4,1:6) ) / 2.0
      edstress(9,1:6) = ( edstress(2,1:6) + edstress(4,1:6) ) / 2.0
      tdstrain(9,1:6) = ( tdstrain(2,1:6) + tdstrain(4,1:6) ) / 2.0
      edstrain(10,1:6) = ( edstrain(3,1:6) + edstrain(4,1:6) ) / 2.0
      edstress(10,1:6) = ( edstress(3,1:6) + edstress(4,1:6) ) / 2.0
      tdstrain(10,1:6) = ( tdstrain(3,1:6) + tdstrain(4,1:6) ) / 2.0
    else if( etype == fe_prism15n ) then
      do i = 1, ni
        ic = 0
        do j = 1, NumOfQuadPoints(etype)
          if( j==1 .or. j==2 .or. j==3 .or. j==7 .or. j==8 .or. j==9 ) then
            ic = ic + 1
            do k = 1, 6
              edstrain(i,k) = edstrain(i,k) + func(i,ic) * gausses(j)%strain_out(k)
              edstress(i,k) = edstress(i,k) + func(i,ic) * gausses(j)%stress_out(k)
              !              tdstrain(i,k) = tdstrain(i,k) + func(i,ic) * gausses(j)%tstrain(k)
            enddo
          endif
        enddo
      enddo
      edstrain(7,1:6) = ( edstrain(1,1:6) + edstrain(2,1:6) ) / 2.0
      edstress(7,1:6) = ( edstress(1,1:6) + edstress(2,1:6) ) / 2.0
      tdstrain(7,1:6) = ( tdstrain(1,1:6) + tdstrain(2,1:6) ) / 2.0
      edstrain(8,1:6) = ( edstrain(2,1:6) + edstrain(3,1:6) ) / 2.0
      edstress(8,1:6) = ( edstress(2,1:6) + edstress(3,1:6) ) / 2.0
      tdstrain(8,1:6) = ( tdstrain(2,1:6) + tdstrain(3,1:6) ) / 2.0
      edstrain(9,1:6) = ( edstrain(3,1:6) + edstrain(1,1:6) ) / 2.0
      edstress(9,1:6) = ( edstress(3,1:6) + edstress(1,1:6) ) / 2.0
      tdstrain(9,1:6) = ( tdstrain(3,1:6) + tdstrain(1,1:6) ) / 2.0
      edstrain(10,1:6) = ( edstrain(4,1:6) + edstrain(5,1:6) ) / 2.0
      edstress(10,1:6) = ( edstress(4,1:6) + edstress(5,1:6) ) / 2.0
      tdstrain(10,1:6) = ( tdstrain(4,1:6) + tdstrain(5,1:6) ) / 2.0
      edstrain(11,1:6) = ( edstrain(5,1:6) + edstrain(6,1:6) ) / 2.0
      edstress(11,1:6) = ( edstress(5,1:6) + edstress(6,1:6) ) / 2.0
      tdstrain(11,1:6) = ( tdstrain(5,1:6) + tdstrain(6,1:6) ) / 2.0
      edstrain(12,1:6) = ( edstrain(6,1:6) + edstrain(4,1:6) ) / 2.0
      edstress(12,1:6) = ( edstress(6,1:6) + edstress(4,1:6) ) / 2.0
      tdstrain(12,1:6) = ( tdstrain(6,1:6) + tdstrain(4,1:6) ) / 2.0
      edstrain(13,1:6) = ( edstrain(1,1:6) + edstrain(4,1:6) ) / 2.0
      edstress(13,1:6) = ( edstress(1,1:6) + edstress(4,1:6) ) / 2.0
      tdstrain(13,1:6) = ( tdstrain(1,1:6) + tdstrain(4,1:6) ) / 2.0
      edstrain(14,1:6) = ( edstrain(2,1:6) + edstrain(5,1:6) ) / 2.0
      edstress(14,1:6) = ( edstress(2,1:6) + edstress(5,1:6) ) / 2.0
      tdstrain(14,1:6) = ( tdstrain(2,1:6) + tdstrain(5,1:6) ) / 2.0
      edstrain(15,1:6) = ( edstrain(3,1:6) + edstrain(6,1:6) ) / 2.0
      edstress(15,1:6) = ( edstress(3,1:6) + edstress(6,1:6) ) / 2.0
      tdstrain(15,1:6) = ( tdstrain(3,1:6) + tdstrain(6,1:6) ) / 2.0
    else if( etype == fe_hex20n ) then
      do i = 1, ni
        ic = 0
        do j = 1, NumOfQuadPoints(etype)
          if( j==1 .or. j==3 .or. j==7 .or. j==9 .or. &
              j==19 .or. j==21 .or. j==25 .or. j==27 ) then
            ic = ic + 1
            do k = 1, 6
              edstrain(i,k) = edstrain(i,k) + func(i,ic) * gausses(j)%strain_out(k)
              edstress(i,k) = edstress(i,k) + func(i,ic) * gausses(j)%stress_out(k)
              !              tdstrain(i,k) = tdstrain(i,k) + func(i,ic) * gausses(j)%tstrain(k)
            enddo
          endif
        enddo
      enddo
      edstrain(9,1:6) = ( edstrain(1,1:6) + edstrain(2,1:6) ) / 2.0
      edstress(9,1:6) = ( edstress(1,1:6) + edstress(2,1:6) ) / 2.0
      tdstrain(9,1:6) = ( tdstrain(1,1:6) + tdstrain(2,1:6) ) / 2.0
      edstrain(10,1:6) = ( edstrain(2,1:6) + edstrain(3,1:6) ) / 2.0
      edstress(10,1:6) = ( edstress(2,1:6) + edstress(3,1:6) ) / 2.0
      tdstrain(10,1:6) = ( tdstrain(2,1:6) + tdstrain(3,1:6) ) / 2.0
      edstrain(11,1:6) = ( edstrain(3,1:6) + edstrain(4,1:6) ) / 2.0
      edstress(11,1:6) = ( edstress(3,1:6) + edstress(4,1:6) ) / 2.0
      tdstrain(11,1:6) = ( tdstrain(3,1:6) + tdstrain(4,1:6) ) / 2.0
      edstrain(12,1:6) = ( edstrain(4,1:6) + edstrain(1,1:6) ) / 2.0
      edstress(12,1:6) = ( edstress(4,1:6) + edstress(1,1:6) ) / 2.0
      tdstrain(12,1:6) = ( tdstrain(4,1:6) + tdstrain(1,1:6) ) / 2.0
      edstrain(13,1:6) = ( edstrain(5,1:6) + edstrain(6,1:6) ) / 2.0
      edstress(13,1:6) = ( edstress(5,1:6) + edstress(6,1:6) ) / 2.0
      tdstrain(13,1:6) = ( tdstrain(5,1:6) + tdstrain(6,1:6) ) / 2.0
      edstrain(14,1:6) = ( edstrain(6,1:6) + edstrain(7,1:6) ) / 2.0
      edstress(14,1:6) = ( edstress(6,1:6) + edstress(7,1:6) ) / 2.0
      tdstrain(14,1:6) = ( tdstrain(6,1:6) + tdstrain(7,1:6) ) / 2.0
      edstrain(15,1:6) = ( edstrain(7,1:6) + edstrain(8,1:6) ) / 2.0
      edstress(15,1:6) = ( edstress(7,1:6) + edstress(8,1:6) ) / 2.0
      tdstrain(15,1:6) = ( tdstrain(7,1:6) + tdstrain(8,1:6) ) / 2.0
      edstrain(16,1:6) = ( edstrain(8,1:6) + edstrain(5,1:6) ) / 2.0
      edstress(16,1:6) = ( edstress(8,1:6) + edstress(5,1:6) ) / 2.0
      tdstrain(16,1:6) = ( tdstrain(8,1:6) + tdstrain(5,1:6) ) / 2.0
      edstrain(17,1:6) = ( edstrain(1,1:6) + edstrain(5,1:6) ) / 2.0
      edstress(17,1:6) = ( edstress(1,1:6) + edstress(5,1:6) ) / 2.0
      tdstrain(17,1:6) = ( tdstrain(1,1:6) + tdstrain(5,1:6) ) / 2.0
      edstrain(18,1:6) = ( edstrain(2,1:6) + edstrain(6,1:6) ) / 2.0
      edstress(18,1:6) = ( edstress(2,1:6) + edstress(6,1:6) ) / 2.0
      tdstrain(18,1:6) = ( tdstrain(2,1:6) + tdstrain(6,1:6) ) / 2.0
      edstrain(19,1:6) = ( edstrain(3,1:6) + edstrain(7,1:6) ) / 2.0
      edstress(19,1:6) = ( edstress(3,1:6) + edstress(7,1:6) ) / 2.0
      tdstrain(19,1:6) = ( tdstrain(3,1:6) + tdstrain(7,1:6) ) / 2.0
      edstrain(20,1:6) = ( edstrain(4,1:6) + edstrain(8,1:6) ) / 2.0
      edstress(20,1:6) = ( edstress(4,1:6) + edstress(8,1:6) ) / 2.0
      tdstrain(20,1:6) = ( tdstrain(4,1:6) + tdstrain(8,1:6) ) / 2.0
    endif
  end subroutine NodalStress_INV3

  function get_mises(s)
    implicit none
    real(kind=kreal) :: get_mises, s(1:6)
    real(kind=kreal) :: s11, s22, s33, s12, s23, s13, ps, smises

    s11 = s(1)
    s22 = s(2)
    s33 = s(3)
    s12 = s(4)
    s23 = s(5)
    s13 = s(6)
    ps = ( s11 + s22 + s33 ) / 3.0d0
    smises = 0.5d0 * ( (s11-ps)**2 + (s22-ps)**2 + (s33-ps)**2 ) + s12**2 + s23**2 + s13**2
    get_mises = dsqrt( 3.0d0 * smises )

  end function get_mises

  !> Calculate NODAL STRESS of plane elements
  !----------------------------------------------------------------------*
  subroutine fstr_NodalStress2D( hecMESH, fstrSOLID )
    !----------------------------------------------------------------------*
    use m_static_lib
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    real(kind=kreal), pointer :: tnstrain(:), testrain(:)
    !C** local variables
    integer(kind=kint) :: itype, icel, ic, is, iE, jS, i, j, ic_type, nn, ni, ID_area
    real(kind=kreal)   :: estrain(4), estress(4), tstrain(4), naturalCoord(4)
    real(kind=kreal)   :: edstrain(8,4), edstress(8,4), tdstrain(8,4)
    real(kind=kreal)   :: s11, s22, s33, s12, s23, s13, ps, smises
    real(kind=kreal), allocatable :: func(:,:), inv_func(:,:)
    integer(kind=kint), allocatable :: nnumber(:)

    tnstrain => fstrSOLID%tnstrain
    testrain => fstrSOLID%testrain
    call fstr_solid_phys_clear(fstrSOLID)

    allocate( nnumber(hecMESH%n_node) )
    if( .not. associated(fstrSOLID%is_rot) ) allocate( fstrSOLID%is_rot(hecMESH%n_node) )
    nnumber = 0
    fstrSOLID%is_rot = 0

    !C +-------------------------------+
    !C | according to ELEMENT TYPE     |
    !C +-------------------------------+
    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if( .not. hecmw_is_etype_surface(ic_type) ) cycle
      !C** set number of nodes and shape function
      nn = hecmw_get_max_node( ic_type )
      ni = NumOfQuadPoints( ic_type )
      allocate( func(ni,nn), inv_func(nn,ni) )
      if( ic_type == fe_tri6n ) then
        ic = hecmw_get_max_node( fe_tri3n )
        do i = 1, ni
          call getQuadPoint( ic_type, i, naturalCoord )
          call getShapeFunc( fe_tri3n, naturalCoord, func(i,1:ic) )
        enddo
        call inverse_func( ic, func, inv_func )
      else if( ic_type == fe_quad4n ) then
        do i = 1, ni
          call getQuadPoint( ic_type, i, naturalCoord )
          call getShapeFunc( ic_type, naturalCoord, func(i,1:nn) )
        enddo
        call inverse_func( ni, func, inv_func )
      else if( ic_type == fe_quad8n ) then
        ic = 0
        do i = 1, ni
          if( i==1 .or. i==3 .or. i==7 .or. i==9 ) then
            ic = ic + 1
            call getQuadPoint( ic_type, i, naturalCoord )
            call getShapeFunc( fe_quad4n, naturalCoord, func(ic,1:4) )
          endif
        enddo
        call inverse_func( ic, func, inv_func )
        ni = ic
      endif
      !C** element loop
      do icel = is, iE
        jS = hecMESH%elem_node_index(icel-1)
        ID_area = hecMESH%elem_ID(icel*2)
        !--- calculate nodal stress and strain
        if( ic_type == fe_tri6n .or. ic_type == fe_quad4n .or. ic_type == fe_quad8n ) then
          call NodalStress_INV2( ic_type, ni, fstrSOLID%elements(icel)%gausses, &
            inv_func, edstrain(1:nn,1:4), edstress(1:nn,1:4), &
            tdstrain(1:nn,1:4) )
        else
          call NodalStress_C2( ic_type, nn, fstrSOLID%elements(icel)%gausses, &
            edstrain(1:nn,1:4), edstress(1:nn,1:4) )
          !          call NodalStress_C2( ic_type, nn, fstrSOLID%elements(icel)%gausses, &
            !                               edstrain(1:nn,1:4), edstress(1:nn,1:4), tdstrain(1:nn,1:4) )
        endif
        do j = 1, nn
          ic = hecMESH%elem_node_item(jS+j)
          fstrSOLID%STRAIN(3*ic-2) = fstrSOLID%STRAIN(3*ic-2) + edstrain(j,1)
          fstrSOLID%STRAIN(3*ic-1) = fstrSOLID%STRAIN(3*ic-1) + edstrain(j,2)
          fstrSOLID%STRAIN(3*ic-0) = fstrSOLID%STRAIN(3*ic-0) + edstrain(j,3)
          fstrSOLID%STRESS(3*ic-2) = fstrSOLID%STRESS(3*ic-2) + edstress(j,1)
          fstrSOLID%STRESS(3*ic-1) = fstrSOLID%STRESS(3*ic-1) + edstress(j,2)
          fstrSOLID%STRESS(3*ic-0) = fstrSOLID%STRESS(3*ic-0) + edstress(j,3)

          if( associated(tnstrain) ) then
            tnstrain(3*ic-2) = tnstrain(3*ic-2) + tdstrain(j,1)
            tnstrain(3*ic-1) = tnstrain(3*ic-1) + tdstrain(j,2)
            tnstrain(3*ic  ) = tnstrain(3*ic  ) + tdstrain(j,3)
          endif
          nnumber(ic) = nnumber(ic) + 1
        enddo
        !--- calculate elemental stress and strain
        !        if( ID_area == hecMESH%my_rank ) then
        call ElementStress_C2( ic_type, fstrSOLID%elements(icel)%gausses, estrain, estress )
        !          call ElementStress_C2( ic_type, fstrSOLID%elements(icel)%gausses, estrain, estress, tstrain )

        fstrSOLID%ESTRAIN(3*icel-2) = estrain(1)
        fstrSOLID%ESTRAIN(3*icel-1) = estrain(2)
        fstrSOLID%ESTRAIN(3*icel-0) = estrain(3)
        fstrSOLID%ESTRESS(3*icel-2) = estress(1)
        fstrSOLID%ESTRESS(3*icel-1) = estress(2)
        fstrSOLID%ESTRESS(3*icel-0) = estress(3)

        !if( associated(testrain) ) then
        !  testrain(3*icel-2) = tstrain(1)
        !  testrain(3*icel-1) = tstrain(2)
        !  testrain(3*icel  ) = tstrain(3)
        !endif
        s11 = estress(1)
        s22 = estress(2)
        s12 = estress(3)
        smises =  0.5d0 * ((s11-s22)**2+(s11)**2+(s22)**2) + 3*s12**2
        fstrSOLID%EMISES(icel) = sqrt( smises )
        !        endif
      enddo
      deallocate( func, inv_func )
    enddo

    !C** average over nodes
    do i = 1, hecMESH%n_node
      if( nnumber(i) == 0 ) cycle
      fstrSOLID%STRAIN(3*i-2:3*i-0) = fstrSOLID%STRAIN(3*i-2:3*i-0) / nnumber(i)
      fstrSOLID%STRESS(3*i-2:3*i-0) = fstrSOLID%STRESS(3*i-2:3*i-0) / nnumber(i)
      if( associated(tnstrain) ) tnstrain(3*i-2:3*i) = tnstrain(3*i-2:3*i) / nnumber(i)
    enddo
    !C** calculate von MISES stress
    do i = 1, hecMESH%n_node
      s11 = fstrSOLID%STRESS(3*i-2)
      s22 = fstrSOLID%STRESS(3*i-1)
      s12 = fstrSOLID%STRESS(3*i-0)
      smises =  0.5d0 * ((s11-s22)**2+(s11)**2+(s22)**2) + 3*s12**2
      fstrSOLID%MISES(i) = sqrt( smises )
    enddo

    deallocate( nnumber )
  end subroutine fstr_NodalStress2D

  !----------------------------------------------------------------------*
  subroutine NodalStress_INV2( etype, ni, gausses, func, edstrain, edstress, tdstrain )
    !----------------------------------------------------------------------*
    use mMechGauss
    integer(kind=kint) :: etype, ni
    type(tGaussStatus) :: gausses(:)
    real(kind=kreal)   :: func(:,:), edstrain(:,:), edstress(:,:), tdstrain(:,:)
    integer :: i, j, k, ic

    edstrain = 0.0d0
    edstress = 0.0d0
    tdstrain = 0.0d0

    if( etype == fe_quad4n ) then
      do i = 1, ni
        do j = 1, ni
          do k = 1, 4
            edstrain(i,k) = edstrain(i,k) + func(i,j) * gausses(j)%strain_out(k)
            edstress(i,k) = edstress(i,k) + func(i,j) * gausses(j)%stress_out(k)
            !            tdstrain(i,k) = tdstrain(i,k) + func(i,j) * gausses(j)%tstrain(k)
          enddo
        enddo
      enddo
    else if( etype == fe_tri6n ) then
      do i = 1, ni
        do j = 1, ni
          do k = 1, 4
            edstrain(i,k) = edstrain(i,k) + func(i,j) * gausses(j)%strain_out(k)
            edstress(i,k) = edstress(i,k) + func(i,j) * gausses(j)%stress_out(k)
            !            tdstrain(i,k) = tdstrain(i,k) + func(i,j) * gausses(j)%tstrain(k)
          enddo
        enddo
      enddo
      edstrain(4,1:4) = ( edstrain(1,1:4) + edstrain(2,1:4) ) / 2.0
      edstress(4,1:4) = ( edstress(1,1:4) + edstress(2,1:4) ) / 2.0
      tdstrain(4,1:4) = ( tdstrain(1,1:4) + tdstrain(2,1:4) ) / 2.0
      edstrain(5,1:4) = ( edstrain(2,1:4) + edstrain(3,1:4) ) / 2.0
      edstress(5,1:4) = ( edstress(2,1:4) + edstress(3,1:4) ) / 2.0
      tdstrain(5,1:4) = ( tdstrain(2,1:4) + tdstrain(3,1:4) ) / 2.0
      edstrain(6,1:4) = ( edstrain(3,1:4) + edstrain(1,1:4) ) / 2.0
      edstress(6,1:4) = ( edstress(3,1:4) + edstress(1,1:4) ) / 2.0
      tdstrain(6,1:4) = ( tdstrain(3,1:4) + tdstrain(1,1:4) ) / 2.0
    else if( etype == fe_quad8n ) then
      do i = 1, ni
        ic = 0
        do j = 1, NumOfQuadPoints(etype)
          if( j==1 .or. j==3 .or. j==7 .or. j==9 ) then
            ic = ic + 1
            do k = 1, 4
              edstrain(i,k) = edstrain(i,k) + func(i,ic) * gausses(j)%strain_out(k)
              edstress(i,k) = edstress(i,k) + func(i,ic) * gausses(j)%stress_out(k)
              !              tdstrain(i,k) = tdstrain(i,k) + func(i,ic) * gausses(j)%tstrain(k)
            enddo
          endif
        enddo
      enddo
      edstrain(5,1:4) = ( edstrain(1,1:4) + edstrain(2,1:4) ) / 2.0
      edstress(5,1:4) = ( edstress(1,1:4) + edstress(2,1:4) ) / 2.0
      tdstrain(5,1:4) = ( tdstrain(1,1:4) + tdstrain(2,1:4) ) / 2.0
      edstrain(6,1:4) = ( edstrain(2,1:4) + edstrain(3,1:4) ) / 2.0
      edstress(6,1:4) = ( edstress(2,1:4) + edstress(3,1:4) ) / 2.0
      tdstrain(6,1:4) = ( tdstrain(2,1:4) + tdstrain(3,1:4) ) / 2.0
      edstrain(7,1:4) = ( edstrain(3,1:4) + edstrain(4,1:4) ) / 2.0
      edstress(7,1:4) = ( edstress(3,1:4) + edstress(4,1:4) ) / 2.0
      tdstrain(7,1:4) = ( tdstrain(3,1:4) + tdstrain(4,1:4) ) / 2.0
      edstrain(8,1:4) = ( edstrain(4,1:4) + edstrain(1,1:4) ) / 2.0
      edstress(8,1:4) = ( edstress(4,1:4) + edstress(1,1:4) ) / 2.0
      tdstrain(8,1:4) = ( tdstrain(4,1:4) + tdstrain(1,1:4) ) / 2.0
    endif
  end subroutine NodalStress_INV2

  !----------------------------------------------------------------------*
  subroutine inverse_func( n, a, inv_a )
    !----------------------------------------------------------------------*
    integer(kind=kint) :: n
    real(kind=kreal)   :: a(:,:), inv_a(:,:)
    integer(kind=kint) :: i, j, k
    real(kind=kreal)   :: buf

    do i = 1, n
      do j = 1, n
        if( i == j ) then
          inv_a(i,j) = 1.0
        else
          inv_a(i,j) = 0.0
        endif
      enddo
    enddo

    do i = 1, n
      buf = 1.0 / a(i,i)
      do j = 1, n
        a(i,j) = a(i,j) * buf
        inv_a(i,j) = inv_a(i,j) *buf
      enddo
      do j = 1, n
        if( i /= j ) then
          buf = a(j,i)
          do k = 1, n
            a(j,k) = a(j,k) - a(i,k) * buf
            inv_a(j,k) = inv_a(j,k) - inv_a(i,k) * buf
          enddo
        endif
      enddo
    enddo
  end subroutine inverse_func

  !> Calculate NODAL STRESS of shell elements
  !----------------------------------------------------------------------*
  subroutine fstr_NodalStress6D( hecMESH, fstrSOLID )
    !----------------------------------------------------------------------*
    use m_static_lib
    type (hecmwST_local_mesh) :: hecMESH
    type (fstr_solid)         :: fstrSOLID
    !C** local variables
    integer(kind=kint) :: itype, icel, is, iE, jS, i, j, k, it, ic, ic_type, nn, isect, ihead, ID_area
    integer(kind=kint) :: nodLOCAL(20), n_layer, ntot_lyr, nlyr, n_totlyr, com_total_layer, shellmatl
    real(kind=kreal)   :: ecoord(3,9), edisp(6,9), estrain(6), estress(6), ndstrain(9,6), ndstress(9,6)
    real(kind=kreal)   :: thick, thick_layer
    real(kind=kreal)   :: s11, s22, s33, s12, s23, s13, t11, t22, t33, t12, t23, t13, ps, smises, tmises
    integer(kind=kint), allocatable :: nnumber(:)
    type(fstr_solid_physic_val), pointer :: layer => null()

    call fstr_solid_phys_clear(fstrSOLID)

    n_totlyr = fstrSOLID%max_lyr

    allocate( nnumber(hecMESH%n_node) )
    if( .not. associated(fstrSOLID%is_rot) ) allocate( fstrSOLID%is_rot(hecMESH%n_node) )
    nnumber = 0
    fstrSOLID%is_rot = 0

    !C +-------------------------------+
    !C | according to ELEMENT TYPE     |
    !C +-------------------------------+
    do itype = 1, hecMESH%n_elem_type
      is = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if( .not. hecmw_is_etype_shell(ic_type) ) then
        ntot_lyr = 0
        cycle
      end if
      nn = hecmw_get_max_node( ic_type )
      !C** element loop
      do icel = is, iE
        jS = hecMESH%elem_node_index(icel-1)
        ID_area = hecMESH%elem_ID(icel*2)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
          ecoord(1,j) = hecMESH%node(3*nodLOCAL(j)-2)
          ecoord(2,j) = hecMESH%node(3*nodLOCAL(j)-1)
          ecoord(3,j) = hecMESH%node(3*nodLOCAL(j)  )
          edisp(1,j) = fstrSOLID%unode(6*nodLOCAL(j)-5)
          edisp(2,j) = fstrSOLID%unode(6*nodLOCAL(j)-4)
          edisp(3,j) = fstrSOLID%unode(6*nodLOCAL(j)-3)
          edisp(4,j) = fstrSOLID%unode(6*nodLOCAL(j)-2)
          edisp(5,j) = fstrSOLID%unode(6*nodLOCAL(j)-1)
          edisp(6,j) = fstrSOLID%unode(6*nodLOCAL(j)  )
        enddo
        isect = hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        thick = hecMESH%section%sect_R_item(ihead+1)
        !--- calculate elemental stress and strain
        if( ic_type == 731 .or. ic_type == 741 .or. ic_type == 743 ) then
          ntot_lyr = fstrSOLID%elements(icel)%gausses(1)%pMaterial%totallyr
          do nlyr=1,ntot_lyr
            call ElementStress_Shell_MITC( ic_type, nn, 6, ecoord, fstrSOLID%elements(icel)%gausses, edisp, &
              & ndstrain(1:nn,1:6), ndstress(1:nn,1:6), thick, 1.0d0, nlyr, ntot_lyr)
            do j = 1, nn
              i = nodLOCAL(j)
              layer => fstrSOLID%SHELL%LAYER(nlyr)%PLUS
              do k = 1, 6
                layer%STRAIN(6*(i-1)+k)   = layer%STRAIN(6*(i-1)+k)  + ndstrain(j,k)
                layer%STRESS(6*(i-1)+k)   = layer%STRESS(6*(i-1)+k)  + ndstress(j,k)
                layer%ESTRAIN(6*(icel-1)+k)  = layer%ESTRAIN(6*(icel-1)+k) + ndstrain(j,k)/nn
                layer%ESTRESS(6*(icel-1)+k)  = layer%ESTRESS(6*(icel-1)+k) + ndstress(j,k)/nn
              enddo
            enddo
            !minus section
            call ElementStress_Shell_MITC( ic_type, nn, 6, ecoord, fstrSOLID%elements(icel)%gausses, edisp, &
              & ndstrain(1:nn,1:6), ndstress(1:nn,1:6), thick,-1.0d0, nlyr, ntot_lyr)
            do j = 1, nn
              i = nodLOCAL(j)
              layer => fstrSOLID%SHELL%LAYER(nlyr)%MINUS
              do k = 1, 6
                layer%STRAIN(6*(i-1)+k)   = layer%STRAIN(6*(i-1)+k)  + ndstrain(j,k)
                layer%STRESS(6*(i-1)+k)   = layer%STRESS(6*(i-1)+k)  + ndstress(j,k)
                layer%ESTRAIN(6*(icel-1)+k)  = layer%ESTRAIN(6*(icel-1)+k) + ndstrain(j,k)/nn
                layer%ESTRESS(6*(icel-1)+k)  = layer%ESTRESS(6*(icel-1)+k) + ndstress(j,k)/nn
              enddo
            enddo
          enddo
          call fstr_getavg_shell(nn,fstrSOLID,icel,nodLOCAL,ndstrain(1:nn,1:6),ndstress(1:nn,1:6),estrain,estress)
        endif

        !if( ID_area == hecMESH%my_rank ) then
        !ADD VALUE and Count node
        do j = 1, nn
          ic = hecMESH%elem_node_item(jS+j)
          fstrSOLID%STRAIN(6*(ic-1)+1:6*(ic-1)+6)  = fstrSOLID%STRAIN(6*(ic-1)+1:6*(ic-1)+6)  + ndstrain(j,1:6)
          fstrSOLID%STRESS(6*(ic-1)+1:6*(ic-1)+6)  = fstrSOLID%STRESS(6*(ic-1)+1:6*(ic-1)+6)  + ndstress(j,1:6)
          !if( associated(tnstrain) )then
          !  tnstrain(6*(ic-1)+1:6*(ic-1)+6) = tnstrain(6*(ic-1)+1:6*(ic-1)+6) + tdstrain(j,1:6)
          !endif
          nnumber(ic) = nnumber(ic) + 1
        enddo

        fstrSOLID%ESTRAIN(6*(icel-1)+1:6*(icel-1)+6) = fstrSOLID%ESTRAIN(6*(icel-1)+1:6*(icel-1)+6) + estrain(1:6)
        fstrSOLID%ESTRESS(6*(icel-1)+1:6*(icel-1)+6) = fstrSOLID%ESTRESS(6*(icel-1)+1:6*(icel-1)+6) + estress(1:6)
        !endif
      enddo
    enddo

    !C** calculate nodal stress and strain
    do i = 1, hecMESH%n_node
      if( nnumber(i) == 0 ) cycle
      fstrSOLID%STRAIN(6*(i-1)+1:6*(i-1)+6) = fstrSOLID%STRAIN(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
      fstrSOLID%STRESS(6*(i-1)+1:6*(i-1)+6) = fstrSOLID%STRESS(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
      !if( associated(tnstrain) )then
      !  tnstrain(6*(i-1)+1:6*(i-1)+6) = tnstrain(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
      !endif
    enddo

    do nlyr = 1, ntot_lyr
      do i = 1, hecMESH%n_node
        fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRAIN(6*(i-1)+1:6*(i-1)+6)  = &
          & fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRAIN(6*(i-1)+1:6*(i-1)+6)  / nnumber(i)
        fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRESS(6*(i-1)+1:6*(i-1)+6)  = &
          & fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRESS(6*(i-1)+1:6*(i-1)+6)  / nnumber(i)
        fstrSOLID%SHELL%LAYER(nlyr)%PLUS%MISES(i) = &
          & get_mises(fstrSOLID%SHELL%LAYER(nlyr)%PLUS%STRESS(6*(i-1)+1:6*(i-1)+6))

        fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRAIN(6*(i-1)+1:6*(i-1)+6) = &
          & fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRAIN(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
        fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRESS(6*(i-1)+1:6*(i-1)+6) = &
          & fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRESS(6*(i-1)+1:6*(i-1)+6) / nnumber(i)
        fstrSOLID%SHELL%LAYER(nlyr)%MINUS%MISES(i) = &
          & get_mises(fstrSOLID%SHELL%LAYER(nlyr)%MINUS%STRESS(6*(i-1)+1:6*(i-1)+6))
      enddo
    enddo

    !C** calculate von MISES stress
    do i = 1, hecMESH%n_node
      fstrSOLID%MISES(i) = get_mises(fstrSOLID%STRESS(6*(i-1)+1:6*(i-1)+6))
    enddo
    do i = 1, hecMESH%n_elem
      fstrSOLID%EMISES(i) = get_mises(fstrSOLID%ESTRESS(6*(i-1)+1:6*(i-1)+6))
    enddo
    deallocate( nnumber )

  end subroutine fstr_NodalStress6D

  subroutine make_principal(fstrSOLID, hecMESH, RES)
    use hecmw_util
    use m_out
    use m_static_lib

    type(fstr_solid)            :: fstrSOLID
    type(hecmwST_local_mesh)    :: hecMESH
    type(fstr_solid_physic_val) :: RES
    integer(kind=kint) :: i, flag
    real(kind=kreal)   :: tmat(3, 3), tvec(3), strain(6)

    flag=ieor(flag,flag)
    if( fstrSOLID%output_ctrl(3)%outinfo%on(19) .or. fstrSOLID%output_ctrl(4)%outinfo%on(19) ) then
      if ( .not. associated(RES%PSTRESS) ) then
        allocate(RES%PSTRESS( 3*hecMESH%n_node ))
      endif
      flag=ior(flag,B'00000001')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(23) .or. fstrSOLID%output_ctrl(4)%outinfo%on(23) ) then
      if ( .not. associated(RES%PSTRESS_VECT) ) then
        allocate(RES%PSTRESS_VECT( 3*hecMESH%n_node ,3))
      endif
      flag=ior(flag,B'00000010')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(21) .or. fstrSOLID%output_ctrl(4)%outinfo%on(21) ) then
      if ( .not. associated(RES%PSTRAIN) ) then
        allocate(RES%PSTRAIN( 3*hecMESH%n_node ))
      endif
      flag=ior(flag,B'00000100')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(25) .or. fstrSOLID%output_ctrl(4)%outinfo%on(25) ) then
      if ( .not. associated(RES%PSTRAIN_VECT) ) then
        allocate(RES%PSTRAIN_VECT( 3*hecMESH%n_node ,3))
      endif
      flag=ior(flag,B'00001000')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(20) .or. fstrSOLID%output_ctrl(4)%outinfo%on(20) ) then
      if ( .not. associated(RES%EPSTRESS) ) then
        allocate(RES%EPSTRESS( 3*hecMESH%n_elem ))
      endif
      flag=ior(flag,B'00010000')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(24) .or. fstrSOLID%output_ctrl(4)%outinfo%on(24) ) then
      if ( .not. associated(RES%EPSTRESS_VECT) ) then
        allocate(RES%EPSTRESS_VECT( 3*hecMESH%n_elem ,3))
      endif
      flag=ior(flag,B'00100000')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(22) .or. fstrSOLID%output_ctrl(4)%outinfo%on(22) ) then
      if ( .not. associated(RES%EPSTRAIN) ) then
        allocate(RES%EPSTRAIN( 3*hecMESH%n_elem ))
      endif
      flag=ior(flag,B'01000000')
    end if
    if( fstrSOLID%output_ctrl(3)%outinfo%on(26) .or. fstrSOLID%output_ctrl(4)%outinfo%on(26) ) then
      if ( .not. associated(RES%EPSTRAIN_VECT) ) then
        allocate(RES%EPSTRAIN_VECT( 3*hecMESH%n_elem ,3))
      endif
      flag=ior(flag,B'10000000')
    end if

    if (iand(flag,B'00000011') /= 0) then
      do i = 1, hecMESH%n_node
        call get_principal(RES%STRESS(6*i-5:6*i), tvec, tmat)
        if (iand(flag,B'00000001') /= 0) RES%PSTRESS(3*(i-1)+1:3*(i-1)+3)=tvec
        if (iand(flag,B'00000010') /= 0) RES%PSTRESS_VECT(3*(i-1)+1:3*(i-1)+3,1:3)=tmat
      end do
    end if
    if (iand(flag,B'00001100') /= 0) then
      do i = 1, hecMESH%n_node
        strain(1:6) = RES%STRAIN(6*i-5:6*i)
        strain(4:6) = 0.5d0*strain(4:6)
        call get_principal(strain, tvec, tmat)
        if (iand(flag,B'00000100') /= 0) RES%PSTRAIN(3*(i-1)+1:3*(i-1)+3)=tvec
        if (iand(flag,B'00001000') /= 0) RES%PSTRAIN_VECT(3*(i-1)+1:3*(i-1)+3,1:3)=tmat
      end do
    end if

    if (iand(flag,B'00110000') /= 0) then
      do i = 1, hecMESH%n_elem
        call get_principal( RES%ESTRESS(6*i-5:6*i), tvec, tmat)
        if (iand(flag,B'00010000') /= 0) RES%EPSTRESS(3*(i-1)+1:3*(i-1)+3)=tvec
        if (iand(flag,B'00100000') /= 0) RES%EPSTRESS_VECT(3*(i-1)+1:3*(i-1)+3,1:3)=tmat
      end do
    end if
    if (iand(flag,B'11000000') /= 0) then
      do i = 1, hecMESH%n_elem
        strain(1:6) = RES%ESTRAIN(6*i-5:6*i)
        strain(4:6) = 0.5d0*strain(4:6)
        call get_principal(strain, tvec, tmat)
        if (iand(flag,B'01000000') /= 0) RES%EPSTRAIN(3*(i-1)+1:3*(i-1)+3)=tvec
        if (iand(flag,B'10000000') /= 0) RES%EPSTRAIN_VECT(3*(i-1)+1:3*(i-1)+3,1:3)=tmat
      end do
    end if
  end subroutine make_principal

end module m_fstr_NodalStress
