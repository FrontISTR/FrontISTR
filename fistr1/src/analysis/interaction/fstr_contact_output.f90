!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Contact output vector processing (initialization, area calculation, parallel)
!>
!> This module provides:
!>   - Initialization of contact output vectors (CONT_NFORCE, CONT_FRIC, CONT_AREA, etc.)
!>   - Contact area calculation
!>   - Parallel processing setup for contact values
module m_fstr_contact_output
  use hecmw
  use m_fstr
  use mContactDef
  use m_fstr_contact_geom
  implicit none

  public :: initialize_contact_output_vectors
  public :: initialize_embed_vectors
  public :: setup_contact_elesurf_for_area
  public :: calc_contact_area
  public :: fstr_setup_parancon_contactvalue
  public :: update_contact_state_vectors

  private :: calc_nodalarea_surfelement

contains

  subroutine initialize_contact_output_vectors(fstrSOLID,hecMAT)
    type(fstr_solid)       :: fstrSOLID      !< type fstr_solid
    type(hecmwST_matrix)   :: hecMAT         !< type hecmwST_matrix

    if( .not. associated(fstrSOLID%CONT_NFORCE) ) then
      allocate( fstrSOLID%CONT_NFORCE(hecMAT%NP*3) )
      fstrSOLID%CONT_NFORCE(:) = 0.d0
    end if

    if( .not. associated(fstrSOLID%CONT_FRIC) ) then
      allocate( fstrSOLID%CONT_FRIC(hecMAT%NP*3) )
      fstrSOLID%CONT_FRIC(:) = 0.d0
    end if

    if( .not. associated(fstrSOLID%CONT_RELVEL) ) then
      allocate( fstrSOLID%CONT_RELVEL(hecMAT%NP*3) )
      fstrSOLID%CONT_RELVEL(:) = 0.d0
    end if

    if( .not. associated(fstrSOLID%CONT_STATE) ) then
      allocate( fstrSOLID%CONT_STATE(hecMAT%NP*1) )
      fstrSOLID%CONT_STATE(:) = 0.d0
    end if

    if( .not. associated(fstrSOLID%CONT_AREA) ) then
      allocate( fstrSOLID%CONT_AREA(hecMAT%NP) )
      fstrSOLID%CONT_AREA(:) = 0.d0
    end if

    if( .not. associated(fstrSOLID%CONT_NTRAC) ) then
      allocate( fstrSOLID%CONT_NTRAC(hecMAT%NP*3) )
      fstrSOLID%CONT_NTRAC(:) = 0.d0
    end if

    if( .not. associated(fstrSOLID%CONT_FTRAC) ) then
      allocate( fstrSOLID%CONT_FTRAC(hecMAT%NP*3) )
      fstrSOLID%CONT_FTRAC(:) = 0.d0
    end if

  end subroutine

  subroutine initialize_embed_vectors(fstrSOLID,hecMAT)
    type(fstr_solid)       :: fstrSOLID      !< type fstr_solid
    type(hecmwST_matrix)   :: hecMAT         !< type hecmwST_matrix

    if( .not. associated(fstrSOLID%EMBED_NFORCE) ) then
      allocate( fstrSOLID%EMBED_NFORCE(hecMAT%NP*3) )
      fstrSOLID%EMBED_NFORCE(:) = 0.d0
    end if
  end subroutine

  subroutine setup_contact_elesurf_for_area( cstep, hecMESH, fstrSOLID )
    integer(kind=kint), intent(in)               :: cstep         !< current step number
    type( hecmwST_local_mesh ), intent(in)       :: hecMESH       !< type mesh
    type(fstr_solid), intent(inout)              :: fstrSOLID     !< type fstr_solid

    integer(kind=kint) :: i, j, k, sgrp_id, iS, iE, eid, sid, n_cdsurfs
    logical, pointer :: cdef_surf(:,:)

    if( associated(fstrSOLID%CONT_SGRP_ID) ) deallocate(fstrSOLID%CONT_SGRP_ID)

    allocate(cdef_surf(l_max_elem_surf,hecMESH%n_elem))
    cdef_surf(:,:) = .false.

    ! label contact defined surfaces
    n_cdsurfs = 0
    do i=1, fstrSOLID%n_contacts
      !grpid = fstrSOLID%contacts(i)%group
      !if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) then
      !  call clear_contact_state(fstrSOLID%contacts(i));  cycle
      !endif

      do k=1,2 !slave,master
        if( k==1 ) then !slave
          sgrp_id = fstrSOLID%contacts(i)%surf_id1_sgrp
        else if( k==2 ) then !master
          sgrp_id = fstrSOLID%contacts(i)%surf_id2
        end if

        if( sgrp_id <= 0 ) cycle

        iS = hecMESH%surf_group%grp_index(sgrp_id-1) + 1
        iE = hecMESH%surf_group%grp_index(sgrp_id  )
        do j=iS,iE
          eid = hecMESH%surf_group%grp_item(2*j-1)
          sid = hecMESH%surf_group%grp_item(2*j)
          ! only internal and boundary element should be added
          if( .not. cdef_surf(sid,eid) ) n_cdsurfs = n_cdsurfs + 1
          cdef_surf(sid,eid) = .true.
        enddo
      end do
    enddo

    !gather info of contact defined surfaces
    allocate(fstrSOLID%CONT_SGRP_ID(2*n_cdsurfs))
    n_cdsurfs = 0
    do i=1,hecMESH%n_elem
      do j=1,l_max_elem_surf
        if( cdef_surf(j,i) ) then
          n_cdsurfs = n_cdsurfs + 1
          fstrSOLID%CONT_SGRP_ID(2*n_cdsurfs-1) = i
          fstrSOLID%CONT_SGRP_ID(2*n_cdsurfs  ) = j
        endif
      end do
    end do
    deallocate(cdef_surf)

  end subroutine

  subroutine calc_contact_area( hecMESH, fstrSOLID, flag )
    type( hecmwST_local_mesh ), intent(in)       :: hecMESH       !< type mesh
    type(fstr_solid), intent(inout)              :: fstrSOLID     !< type fstr_solid
    integer(kind=kint), intent(in)               :: flag          !< set 1 if called in NR iteration

    integer(kind=kint), parameter :: NDOF=3
    integer(kind=kint) :: i, isuf, icel, sid, etype, nn, iS, idx
    integer(kind=kint) :: ndlocal(l_max_elem_node)
    real(kind=kreal), allocatable :: coord(:)
    real(kind=kreal)   :: ecoord(NDOF,l_max_elem_node), vect(l_max_elem_node)

    fstrSOLID%CONT_AREA(:) = 0.d0

    if( .not. associated(fstrSOLID%CONT_SGRP_ID) ) return

    allocate(coord(NDOF*hecMESH%n_node))
    do i=1,NDOF*hecMESH%n_node
      coord(i) = hecMESH%node(i)+fstrSOLID%unode(i)
    end do
    if( flag == 1 ) then
      do i=1,NDOF*hecMESH%n_node
        coord(i) = coord(i)+fstrSOLID%dunode(i)
      end do
    end if

    do isuf=1,size(fstrSOLID%CONT_SGRP_ID)/2
      icel = fstrSOLID%CONT_SGRP_ID(2*isuf-1)
      sid  = fstrSOLID%CONT_SGRP_ID(2*isuf  )

      etype = hecMESH%elem_type(icel)
      nn = hecmw_get_max_node(etype)
      iS = hecMESH%elem_node_index(icel-1)
      ndlocal(1:nn) = hecMESH%elem_node_item (iS+1:iS+nn)

      do i=1,nn
        idx = NDOF*(ndlocal(i)-1)
        ecoord(1:NDOF,i) = coord(idx+1:idx+NDOF)
      end do

      call calc_nodalarea_surfelement( etype, nn, ecoord, sid, vect )

      do i=1,nn
        idx = ndlocal(i)
        fstrSOLID%CONT_AREA(idx) = fstrSOLID%CONT_AREA(idx) + vect(i)
      end do

    end do

    deallocate(coord)
  end subroutine

  subroutine calc_nodalarea_surfelement( etype, nn, ecoord, sid, vect )
    integer(kind=kint), intent(in)   :: etype
    integer(kind=kint), intent(in)   :: nn
    real(kind=kreal), intent(in)     :: ecoord(:,:)
    integer(kind=kint), intent(in)   :: sid
    real(kind=kreal), intent(out)    :: vect(:)

    integer(kind=kint), parameter :: NDOF=3
    integer(kind=kint) :: nod(l_max_surface_node)
    integer(kind=kint) :: nsur, stype, ig0, i
    real(kind=kreal)   :: localcoord(2), normal(3), area, wg
    real(kind=kreal)   :: scoord(NDOF,l_max_surface_node), H(l_max_surface_node)

    vect(:) = 0.d0

    call getSubFace( etype, sid, stype, nod )
    nsur = getNumberOfNodes( stype )
    do i=1,nsur
      scoord(1:NDOF,i) = ecoord(1:NDOF,nod(i))
    end do

    area = 0.d0
    do ig0=1,NumOfQuadPoints( stype )
      call getQuadPoint( stype, ig0, localcoord(1:2) )
      call getShapeFunc( stype, localcoord(1:2), H(1:nsur) )

      wg=getWeight( stype, ig0 )
      ! normal = dx/dr_1 \times dx/dr_2
      normal(1:3) = SurfaceNormal( stype, nsur, localcoord(1:2), scoord(1:NDOF,1:nsur) )
      area = area + dsqrt(dot_product(normal,normal))*wg
    enddo
    area = area/dble(nsur)
    do i=1,nsur
      vect(nod(i)) = area
    end do

  end subroutine

  subroutine fstr_setup_parancon_contactvalue(hecMESH,ndof,vec,vtype)
    use m_fstr
    implicit none
    type(hecmwST_local_mesh), intent(in)      :: hecMESH
    integer(kind=kint), intent(in)            :: ndof
    real(kind=kreal), pointer, intent(inout)  :: vec(:)
    integer(kind=kint), intent(in)            :: vtype !1:value, 2:state
    !
    integer(kind=kint) ::  i,N,i0,N_loc
    integer(kind=kint) :: offset, pid, lid
    integer(kind=kint), allocatable :: displs(:)
    real(kind=kreal), allocatable   :: vec_all(:)

    !
    N_loc = hecMESH%nn_internal
    allocate(displs(0:nprocs))
    displs(:) = 0
    displs(myrank+1) = N_loc
    call hecmw_allreduce_I(hecMESH, displs, nprocs+1, hecmw_sum)
    do i=1,nprocs
      displs(i) = displs(i-1) + displs(i)
    end do
    offset = displs(myrank)
    N = displs(nprocs)

    allocate(vec_all(ndof*N))

    if( vtype == 1 ) then
      vec_all(:) = 0.d0
      do i= 1, hecMESH%n_node
        pid = hecMESH%node_ID(i*2)
        lid = hecMESH%node_ID(i*2-1)
        i0 = (displs(pid) + (lid-1))*ndof
        vec_all(i0+1:i0+ndof) = vec((i-1)*ndof+1:i*ndof)
        vec((i-1)*ndof+1:i*ndof) = 0.d0
      enddo

      call hecmw_allreduce_R(hecMESH, vec_all, N*ndof, hecmw_sum)
      do i= 1, hecMESH%n_node
        pid = hecMESH%node_ID(i*2)
        lid = hecMESH%node_ID(i*2-1)
        i0 = (displs(pid) + (lid-1))*ndof
        if (dot_product(vec_all(i0+1:i0+ndof),vec_all(i0+1:i0+ndof)) == 0.d0) cycle
        vec((i-1)*ndof+1:i*ndof) = vec_all(i0+1:i0+ndof)
      enddo
      do i=1,ndof*N_loc
        vec(i) = vec_all(offset*ndof+i)
      end do
    else if( vtype == 2 ) then
      vec_all(:) = -1000.d0
      do i= 1, hecMESH%n_node
        if( vec(i) == 0.d0 ) cycle
        pid = hecMESH%node_ID(i*2)
        lid = hecMESH%node_ID(i*2-1)
        i0 = displs(pid) + lid
        vec_all(i0) = vec(i)
      enddo

      call hecmw_allreduce_R(hecMESH, vec_all, N, hecmw_max)
      do i= 1, hecMESH%n_node
        pid = hecMESH%node_ID(i*2)
        lid = hecMESH%node_ID(i*2-1)
        i0 = displs(pid) + lid
        if(vec_all(i0) == -1000.d0) cycle
        if(vec(i) < vec_all(i0)) vec(i) = vec_all(i0)
      enddo
    end if

    deallocate(displs,vec_all)
  end subroutine

  !> Update contact state output vectors (CONT_RELVEL, CONT_STATE)
  subroutine update_contact_state_vectors( contact, dt, relvel_vec, state_vec )
    type( tContact ), intent(in)      :: contact        !< contact info
    real(kind=kreal), intent(in)      :: dt             !< time increment
    real(kind=kreal), intent(inout)   :: relvel_vec(:) !< relative velocity vector
    real(kind=kreal), intent(inout)   :: state_vec(:)  !< contact state vector

    integer(kind=kint)  :: i, slave

    do i= 1, size(contact%slave)
      slave = contact%slave(i)
      if( state_vec(slave) < 0.1d0 .or. contact%states(i)%state > 0 ) &
      &  state_vec(slave) = dble(contact%states(i)%state)

      if( is_contact_free(contact%states(i)%state) ) cycle   ! not in contact or near
      if( dt < 1.d-16 ) cycle ! too small delta t
      relvel_vec(3*slave-2:3*slave) = contact%states(i)%reldisp(1:3)/dt
    enddo

  end subroutine update_contact_state_vectors

end module m_fstr_contact_output
