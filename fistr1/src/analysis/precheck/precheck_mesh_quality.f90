!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides mesh quality check functions for ELEMCHECK
module m_precheck_mesh_quality

  use hecmw
  use m_fstr
  use m_precheck_LIB_elements

  implicit none

  private
  public :: precheck_mesh_quality

contains

  !> Get section thickness
  subroutine precheck_get_thickness(hecMESH, mid, thick)
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: mid, ihead
    real(kind=kreal)   :: thick

    ihead = hecMESH%section%sect_R_index(mid-1)
    thick = hecMESH%section%sect_R_item(ihead+1)
  end subroutine precheck_get_thickness

  !> Mesh quality check (element volume, aspect ratio)
  !! Returns per-element volume and aspect ratio in elem_vol/elem_asp arrays
  subroutine precheck_mesh_quality(hecMESH, hecMAT, elem_vol, elem_asp)
    implicit none

    type(hecmwST_matrix)     :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    real(kind=kreal), intent(out) :: elem_vol(:)
    real(kind=kreal), intent(out) :: elem_asp(:)

    !** Local variables
    integer(kind=kint) :: nelem, mid, j, isect, nline, tline, icel, iiS
    integer(kind=kint) :: ndof2, nelem_wo_mpc
    integer(kind=kint) :: ie, ia, jelem, ic_type, nn, jS, jE, itype
    integer(kind=kint) :: nodLOCAL(20), NTOTsum(1)
    integer(kind=kint) :: nonzero
    real(kind=kreal)   :: ntdof2
    real(kind=kreal)   :: al, almin, almax, AA, thick, vol, avvol
    real(kind=kreal)   :: tvol, tvmax, tvmin, tlmax, tlmin, asp, aspmax
    real(kind=kreal)   :: xx(20), yy(20), zz(20)
    real(kind=kreal)   :: TOTsum(1), TOTmax(3), TOTmin(2)

    !C INIT
    elem_vol = 0.0d0
    elem_asp = 0.0d0
    nelem  = 0
    tvol   = 0.0
    tvmax  = 0.0
    tvmin  = 1.0e+20
    tlmax  = 0.0
    tlmin  = 1.0e+20
    aspmax = 0.0

    !C Mesh Summary
    write(ILOG,"(a)") '###  Mesh Summary  ###'
    write(ILOG,"(a,i12)") '  Num of node:',hecMESH%n_node
    write(*   ,"(a,i12)") '  Num of node:',hecMESH%n_node
    write(ILOG,"(a,i12)") '  Num of DOF :',hecMESH%n_node*hecMESH%n_dof
    write(*   ,"(a,i12)") '  Num of DOF :',hecMESH%n_node*hecMESH%n_dof
    ndof2  = hecMESH%n_dof**2
    ntdof2 = dble(hecMESH%n_node*hecMESH%n_dof)**2
    write(ILOG,"(a,i12)") '  Num of elem:',hecMESH%n_elem
    write(*   ,"(a,i12)") '  Num of elem:',hecMESH%n_elem
    nelem_wo_mpc = 0
    do itype = 1, hecMESH%n_elem_type
      jS = hecMESH%elem_type_index(itype-1)
      jE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      if(.not. (hecmw_is_etype_link(ic_type) .or. hecmw_is_etype_patch(ic_type))) &
           nelem_wo_mpc = nelem_wo_mpc + jE-jS
    enddo
    write(ILOG,"(a,i12)") '   ** w/o MPC/Patch:',nelem_wo_mpc
    write(*   ,"(a,i12)") '   ** w/o MPC/Patch:',nelem_wo_mpc
    do itype = 1, hecMESH%n_elem_type
      jS = hecMESH%elem_type_index(itype-1)
      jE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)
      write(ILOG,"(a,i4,a,i12)") '  Num of ',ic_type,':',jE-jS
      write(*   ,"(a,i4,a,i12)") '  Num of ',ic_type,':',jE-jS
    enddo
    nonzero = ndof2*(hecMAT%NP + hecMAT%NPU + hecMAT%NPL)
    write(ILOG,"(a,i12)") '  Num of NZ  :',nonzero
    write(*   ,"(a,i12)") '  Num of NZ  :',nonzero
    write(ILOG,"(a,1pe12.5,a)") '  Sparsity   :',100.0d0*dble(nonzero)/dble(ntdof2),"[%]"
    write(*   ,"(a,1pe12.5,a)") '  Sparsity   :',100.0d0*dble(nonzero)/dble(ntdof2),"[%]"

    !C 3D
    if( hecMESH%n_dof .eq. 3 ) then
      do ie = 1, hecMESH%n_elem
        ia = hecMESH%elem_ID(ie*2)
        if( ia.ne.hecMESH%my_rank ) cycle
        jelem = hecMESH%global_elem_ID(ie)
        ic_type = hecMESH%elem_type(ie)
        if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_solid(ic_type) &
            & .or. HECMW_is_etype_beam(ic_type) .or. HECMW_is_etype_shell(ic_type))) then
          write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
          cycle
        endif
        vol = 0.0d0
        almax = 0.0d0
        almin = 0.0d0
        nn = hecmw_get_max_node(ic_type)
        jS = hecMESH%elem_node_index(ie-1)
        jE = hecMESH%elem_node_index(ie)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
          xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
          yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
          zz(j) = hecMESH%node(3*nodLOCAL(j)  )
        enddo

        if    ( ic_type.eq.111 ) then
          isect = hecMESH%section_ID(ie)
          mid = hecMESH%section%sect_mat_ID_item(isect)
          call precheck_get_thickness( hecMESH,mid,AA )
          al = sqrt( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
          nline = 1
          tline = al
          vol = AA*al
          almax = al
          almin = al
        elseif( hecmw_is_etype_solid(ic_type) ) then
          call precheck_calc_vol_asp(ic_type, nn, xx, yy, zz, 0.0d0, vol, almax, almin)
        elseif( ic_type.eq.641 ) then
          vol = 1.0d-12
        elseif( ic_type.eq.761 .or. ic_type.eq.781 ) then
          vol = 1.0d-12
        endif

        if( vol.le.0.0 ) then
          write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
        endif
        nelem = nelem + 1
        elem_vol(ie) = vol
        tvol = tvol + vol
        if( vol.gt.tvmax ) tvmax = vol
        if( vol.lt.tvmin ) tvmin = vol
        if( almax.gt.tlmax ) tlmax = almax
        if( almin.lt.tlmin ) tlmin = almin
        asp = almax/almin
        elem_asp(ie) = asp
        if( asp.gt.aspmax ) aspmax = asp
        if( asp.gt.50 ) then
          write(ILOG,*) '  %%%  WARNING %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
          write(ILOG,*) '      Maximum length =',almax
          write(ILOG,*) '      Minimum length =',almin
        endif
      enddo

    !C 2D
    elseif( hecMESH%n_dof .eq. 2 ) then
      do ie = 1, hecMESH%n_elem
        ia = hecMESH%elem_ID(ie*2)
        if( ia.ne.hecMESH%my_rank ) cycle
        jelem = hecMESH%global_elem_ID(ie)
        ic_type = hecMESH%elem_type(ie)
        if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_surface(ic_type))) then
          write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
          cycle
        endif
        vol = 0.0d0
        almax = 0.0d0
        almin = 0.0d0
        nn = hecmw_get_max_node(ic_type)
        jS = hecMESH%elem_node_index(ie-1)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
          xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
          yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
        enddo

        isect = hecMESH%section_ID(ie)
        mid = hecMESH%section%sect_mat_ID_item(isect)
        call precheck_get_thickness( hecMESH,mid,AA )

        if    ( ic_type.eq.111 ) then
          al = sqrt( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2 )
          vol = AA*al
          if( al.gt.tlmax ) tlmax = al
          if( al.lt.tlmin ) tlmin = al
          aspmax = 1.0
        elseif( hecmw_is_etype_surface(ic_type) ) then
          call precheck_calc_vol_asp(ic_type, nn, xx, yy, zz, AA, vol, almax, almin)
        else
          vol = 0.0
        endif

        if( vol.le.0.0 ) then
          write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
        endif
        nelem = nelem + 1
        elem_vol(ie) = vol
        tvol = tvol + vol
        if( vol.gt.tvmax ) tvmax = vol
        if( vol.lt.tvmin ) tvmin = vol
        if( almax.gt.tlmax ) tlmax = almax
        if( almin.lt.tlmin ) tlmin = almin
        asp = almax/almin
        elem_asp(ie) = asp
        if( asp.gt.aspmax ) aspmax = asp
        if( asp.gt.50 ) then
          write(ILOG,*) '  %%%  WARNING %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
          write(ILOG,*) '      Maximum length =',almax
          write(ILOG,*) '      Minimum length =',almin
        endif
      enddo

    !C SHELL
    elseif( hecMESH%n_dof .eq. 6 ) then
      do ie = 1, hecMESH%n_elem
        ia = hecMESH%elem_ID(ie*2)
        if( ia.ne.hecMESH%my_rank ) cycle
        jelem = hecMESH%global_elem_ID(ie)
        ic_type = hecMESH%elem_type(ie)
        if (.not. (hecmw_is_etype_beam(ic_type) .or. hecmw_is_etype_shell(ic_type))) then
          write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
          cycle
        endif
        vol = 0.0d0
        almax = 0.0d0
        almin = 0.0d0
        nn = hecmw_get_max_node(ic_type)
        jS = hecMESH%elem_node_index(ie-1)
        do j = 1, nn
          nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
          xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
          yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
          zz(j) = hecMESH%node(3*nodLOCAL(j)  )
        enddo

        isect = hecMESH%section_ID(ie)
        mid = hecMESH%section%sect_mat_ID_item(isect)
        call precheck_get_thickness( hecMESH,mid,AA )

        if    ( ic_type.eq.111 ) then
          al = sqrt( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
          nline = nline + 1
          tline = tline + al
          vol = AA*al
          if( al.gt.tlmax ) tlmax = al
          if( al.lt.tlmin ) tlmin = al
          aspmax = 1.0
        elseif( ic_type.eq.731 .or. ic_type.eq.741 ) then
          call precheck_calc_vol_asp(ic_type, nn, xx, yy, zz, AA, vol, almax, almin)
        endif

        if( vol.le.0.0 ) then
          write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
        endif
        nelem = nelem + 1
        elem_vol(ie) = vol
        tvol = tvol + vol
        if( vol.gt.tvmax ) tvmax = vol
        if( vol.lt.tvmin ) tvmin = vol
        if( almax.gt.tlmax ) tlmax = almax
        if( almin.lt.tlmin ) tlmin = almin
        asp = almax/almin
        elem_asp(ie) = asp
        if( asp.gt.aspmax ) aspmax = asp
        if( asp.gt.50 ) then
          write(ILOG,*) '  %%%  WARNING %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
          write(ILOG,*) '      Maximum length =',almax
          write(ILOG,*) '      Minimum length =',almin
        endif
      enddo
    endif

    !C Summary
    avvol = tvol / nelem
    write(ILOG,*) '###  Sub Summary  ###'
    write(ILOG,*) ' Total Volumes in this region        = ',tvol
    write(ILOG,*) ' Average Volume of elements          = ',avvol
    write(ILOG,*) ' Maximum Volume of elements          = ',tvmax
    write(ILOG,*) ' Minimum Volume of elements          = ',tvmin
    write(ILOG,*) ' Maximum length of element edges     = ',tlmax
    write(ILOG,*) ' Minimum length of element edges     = ',tlmin
    write(ILOG,*) ' Maximum aspect ratio in this region = ',aspmax
    TOTsum(1) = tvol
    call hecmw_allREDUCE_R(hecMESH,TOTsum,1,hecmw_sum)
    NTOTsum(1) = nelem
    call hecmw_allREDUCE_I(hecMESH,NTOTsum,1,hecmw_sum)
    TOTmax(1) = tvmax
    TOTmax(2) = tlmax
    TOTmax(3) = aspmax
    call hecmw_allREDUCE_R(hecMESH,TOTmax,3,hecmw_max)
    TOTmin(1) = tvmin
    TOTmin(2) = tlmin
    call hecmw_allREDUCE_R(hecMESH,TOTmin,2,hecmw_min)
    if( hecMESH%my_rank .eq. 0 ) then
      avvol = TOTsum(1) / NTOTsum(1)
      write(ILOG,*) '###  Global Summary  ###'
      write(ILOG,*) ' TOTAL VOLUME = ',TOTsum(1)
      write(*,*)    ' TOTAL VOLUME = ',TOTsum(1)
      write(ILOG,*) ' AVERAGE VOLUME OF ELEMENTS = ',avvol
      write(*,*)    ' AVERAGE VOLUME OF ELEMENTS = ',avvol
      write(ILOG,*) ' MAXIMUM VOLUME OF ELEMENTS = ',TOTmax(1)
      write(*,*)    ' MAXIMUM VOLUME OF ELEMENTS = ',TOTmax(1)
      write(ILOG,*) ' MINIMUM VOLUME OF ELEMENTS = ',TOTmin(1)
      write(*,*)    ' MINIMUM VOLUME OF ELEMENTS = ',TOTmin(1)
      write(ILOG,*) ' MAXIMUM LENGTH OF ELEMENT EDGES = ',TOTmax(2)
      write(*,*)    ' MAXIMUM LENGTH OF ELEMENT EDGES = ',TOTmax(2)
      write(ILOG,*) ' MINIMUM LENGTH OF ELEMENT EDGES = ',TOTmin(2)
      write(*,*)    ' MINIMUM LENGTH OF ELEMENT EDGES = ',TOTmin(2)
      write(ILOG,*) ' MAXIMUM ASPECT RATIO  = ',TOTmax(3)
      write(*,*)    ' MAXIMUM ASPECT RATIO  = ',TOTmax(3)
    endif

  end subroutine precheck_mesh_quality

end module m_precheck_mesh_quality
