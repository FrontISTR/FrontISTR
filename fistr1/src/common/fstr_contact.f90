!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provides functions to calculate contact stiff matrix
module mContact

  use mContactDef
  use hecmw
  use m_fstr
  implicit none

  private :: l_contact2mpc, l_tied2mpc
  integer(kind=kint), save :: n_contact_mpc
  logical, private :: active

  real(kind=kreal), save :: mu=1.d10  !< penalty, default value
  real(kind=kreal), save :: mut=1.d6  !< penalty along tangent direction
  real(kind=kreal), save :: cdotp=1.d3  !< mu=cdotp*maxval

  real(kind=kreal), save :: cgn=1.d-5 !< convergent condition of penetration
  real(kind=kreal), save :: cgt=1.d-3 !< convergent condition of relative tangent disp

  real(kind=kreal), save :: gnt(2)    !< 1:current average penetration;
  !< 2:current relative tangent displacement
  real(kind=kreal), save :: bakgnt(2) !< 1:current average penetration;
  !< 2:current relative tangent displacement!

contains

  !> Write out the contact definition read from mesh file
  subroutine print_contatct_pair( file, pair )
    integer(kind=kint), intent(in)           :: file
    type( hecmwST_contact_pair ), intent(in) :: pair

    integer(kind=kint) :: i
    write(file,*) "Number of contact pair", pair%n_pair
    do i=1,pair%n_pair
      write(file,*) trim(pair%name(i)), pair%type(i), pair%slave_grp_id(i)  &
        ,pair%master_grp_id(i), pair%slave_orisgrp_id(i)
    enddo
  end subroutine

  subroutine fstr_set_contact_penalty( maxv )
    real(kind=kreal), intent(in) :: maxv
    mu = cdotp * maxv
    if( gnt(1)<1.d-3 ) mu=cdotp*10.d0* maxv
    bakgnt = 0.d0
  end subroutine

  logical function fstr_is_contact_active()
    fstr_is_contact_active = active
  end function

  subroutine fstr_set_contact_active( a )
    logical, intent(in) :: a
    active = a
  end subroutine

  logical function fstr_is_contact_conv(ctAlgo,infoCTChange,hecMESH)
    integer(kind=kint), intent(in)  :: ctAlgo                    !< contact analysis algorithm
    type (fstr_info_contactChange), intent(in) :: infoCTChange   !< fstr_contactChange
    type (hecmwST_local_mesh), intent(in) :: hecMESH
    integer(kind=kint) :: is_conv

    fstr_is_contact_conv = .false.
    if( infoCTChange%contact2free+infoCTChange%contact2neighbor+      &
    infoCTChange%contact2difflpos+infoCTChange%free2contact == 0 ) &
    fstr_is_contact_conv = .true.

    is_conv = 0
    if (fstr_is_contact_conv) is_conv = 1
    call hecmw_allreduce_I1(hecMESH, is_conv, HECMW_MIN)
    if (is_conv == 0) then
      fstr_is_contact_conv = .false.
    else
      fstr_is_contact_conv = .true.
    endif
  end function

  logical function fstr_is_matrixStructure_changed(infoCTChange)
    type (fstr_info_contactChange)   :: infoCTChange  !< fstr_contactChange
    fstr_is_matrixStructure_changed = .false.
    if( infoCTChange%contact2free+infoCTChange%contact2neighbor+infoCTChange%free2contact > 0 ) &
      fstr_is_matrixStructure_changed = .true.
  end function

  !> Contact state to equation conditions
  subroutine l_contact2mpc( contact, mpcs, nmpc )
    use fstr_ctrl_modifier
    type( tContact ), intent(in)         :: contact  !< current contact state
    type( hecmwST_mpc ), intent(inout)   :: mpcs     !< to who mpc be appended
    integer(kind=kint), intent(out)      :: nmpc     !< number of mpc conditions appended
    integer(kind=kint), parameter :: ndof = 3        ! 3D problem only, currently
    real(kind=kreal), parameter   :: tol =1.d-10
    integer(kind=kint) :: i, j, k, nn, csurf, nenode, etype, tdof
    integer(kind=kint) :: nodes(l_max_surface_node*ndof+ndof), dofs(l_max_surface_node*ndof+ndof)
    real(kind=kreal) :: values(l_max_surface_node*ndof+ndof+1),val(l_max_surface_node*ndof+ndof+1)
    nmpc=0
    do i=1,size(contact%states)
      if( contact%states(i)%state == -1 ) cycle    ! in free
      csurf = contact%states(i)%surface
      if( csurf<=0 ) stop "error in contact state"
      etype =  contact%master(csurf)%etype
      nenode = size(contact%master(csurf)%nodes)
      tdof = nenode*ndof+ndof
      call contact2mpcval( contact%states(i), etype, nenode, values(1:tdof+1) )
      tdof = 0
      do j=1,ndof
        if( dabs(values(j))<tol ) cycle
        tdof = tdof+1
        nodes(tdof) = contact%slave(i)
        dofs(tdof) = j
        val(tdof) = values(j)
      enddo
      do j=1,nenode
        nn =  contact%master(csurf)%nodes(j)
        nodes( j*ndof+1:j*ndof+ndof ) =  nn
        do k=1,ndof
          if( dabs(values(j*ndof+k)) < tol ) cycle
          tdof=tdof+1
          nodes(tdof)=nn
          dofs(tdof ) = k
          val(tdof)=values(j*ndof+k)
        enddo
      enddo
      val(tdof+1) = values(nenode*ndof+ndof+1)

      call fstr_append_mpc( tdof, nodes(1:tdof), dofs(1:tdof), val(1:tdof+1), mpcs )
      nmpc=nmpc+1
    enddo
  end subroutine l_contact2mpc

  !> Rigid connect condition to equation conditions
  subroutine l_tied2mpc( contact, mpcs, nmpc )
    use fstr_ctrl_modifier
    type( tContact ), intent(in)         :: contact  !< current contact state
    type( hecmwST_mpc ), intent(inout)   :: mpcs     !< to who mpc be appended
    integer(kind=kint), intent(out)      :: nmpc     !< number of mpc conditions appended
    integer(kind=kint) :: i, j, csurf, nenode, etype, tdof
    integer(kind=kint) :: nodes(l_max_surface_node+1), dofs(l_max_surface_node+1)
    real(kind=kreal) :: values(l_max_surface_node+2)
    nmpc=0
    do i=1,size(contact%slave)
      csurf = contact%states(i)%surface
      if( csurf<=0 ) cycle                         ! contactor not exists
      nenode = size(contact%master(csurf)%nodes)
      tdof = nenode+1
      nodes(1) = contact%slave(i)
      nodes( 2:tdof ) =  contact%master(csurf)%nodes(:)
      values(1) = -1.d0
      values(2:tdof) = 1.d0
      values(tdof+1) = 0.d0
      etype =  contact%master(csurf)%etype
      do j=1,3
        dofs(1:tdof) = j
        call fstr_append_mpc( tdof, nodes(1:tdof), dofs(1:tdof), values(1:tdof+1), mpcs )
        nmpc=nmpc+1
      enddo
    enddo
  end subroutine l_tied2mpc

  !> Contact states to equation conditions
  subroutine fstr_contact2mpc( contacts, mpcs )
    type( tContact ), intent(in)       :: contacts(:)  !< current contact state
    type( hecmwST_mpc ), intent(inout) :: mpcs         !< to who mpc be appended
    integer(kind=kint) :: i, nmpc
    n_contact_mpc = 0
    do i=1,size(contacts)
      if( contacts(i)%algtype == CONTACTUNKNOWN ) cycle     ! not initialized
      if( contacts(i)%algtype == CONTACTFSLID ) then
        print *, "Cannot deal with finit slip problems by MPC!"
        cycle
      endif
      if( contacts(i)%algtype == CONTACTSSLID ) then
        call l_contact2mpc( contacts(i), mpcs, nmpc )
        n_contact_mpc = n_contact_mpc + nmpc
      elseif( contacts(i)%algtype == CONTACTTIED ) then
        call l_tied2mpc( contacts(i), mpcs, nmpc )
        n_contact_mpc = n_contact_mpc + nmpc
      endif
    enddo
  end subroutine

  !> Delete mpcs derived from contact conditions
  subroutine fstr_del_contactmpc( mpcs )
    use fstr_ctrl_modifier
    type( hecmwST_mpc ), intent(inout) :: mpcs       !<  mpcs to be modified
    call fstr_delete_mpc( n_contact_mpc, mpcs )
  end subroutine

  !> Print out mpc conditions
  subroutine fstr_write_mpc( file, mpcs )
    integer(kind=kint), intent(in)  :: file       !<  file number
    type( hecmwST_mpc ), intent(in) :: mpcs       !<  mpcs to be printed

    integer(kind=kint) :: i,j,n0,n1
    write(file, *) "Number of equation", mpcs%n_mpc
    do i=1,mpcs%n_mpc
      write(file,*) "--Equation",i
      n0=mpcs%mpc_index(i-1)+1
      n1=mpcs%mpc_index(i)
      write(file, *) n0,n1
      write(file,'(30i5)') (mpcs%mpc_item(j),j=n0,n1)
      write(file,'(30i5)') (mpcs%mpc_dof(j),j=n0,n1)
      write(file,'(30f7.2)') (mpcs%mpc_val(j),j=n0,n1),mpcs%mpc_const(i)
    enddo
  end subroutine

  !> Scanning contact state
  subroutine fstr_scan_contact_state( cstep, sub_step, cont_step, dt, ctAlgo, hecMESH, fstrSOLID, infoCTChange, B )
    integer(kind=kint), intent(in)         :: cstep      !< current step number
    integer(kind=kint), intent(in)         :: sub_step   !< current sub-step number
    integer(kind=kint), intent(in)         :: cont_step  !< current contact step number
    real(kind=kreal), intent(in)           :: dt
    integer(kind=kint), intent(in)         :: ctAlgo     !< contact analysis algorithm
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
    type(fstr_info_contactChange), intent(inout):: infoCTChange   !<
    !      logical, intent(inout)                 :: changed     !< if contact state changed
    real(kind=kreal), optional             :: B(:)        !< nodal force residual
    character(len=9)                       :: flag_ctAlgo !< contact analysis algorithm flag
    integer(kind=kint) :: i, grpid
    logical :: iactive, is_init

    if( associated( fstrSOLID%CONT_RELVEL ) ) fstrSOLID%CONT_RELVEL(:) = 0.d0
    if( associated( fstrSOLID%CONT_STATE ) ) fstrSOLID%CONT_STATE(:) = 0.d0

    if( ctAlgo == kcaSLAGRANGE ) then
      flag_ctAlgo = 'SLagrange'
    elseif( ctAlgo == kcaALAGRANGE ) then
      flag_ctAlgo = 'ALagrange'
    endif

    ! P.A. We redefine fstrSOLID%ddunode as current coordinate of every nodes
    !  fstrSOLID%ddunode(:) = fstrSOLID%unode(:) + fstrSOLID%dunode(:)
    do i = 1, size(fstrSOLID%unode)
      fstrSOLID%ddunode(i) = hecMESH%node(i) + fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo
    active = .false.

    infoCTChange%contact2free = 0
    infoCTChange%contact2neighbor = 0
    infoCTChange%contact2diffLpos = 0
    infoCTChange%free2contact = 0
    infoCTChange%contactNode_current = 0

    is_init = ( cstep == 1 .and. sub_step == 1 .and. cont_step == 0 )

    do i=1,fstrSOLID%n_contacts
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) then
        call clear_contact_state(fstrSOLID%contacts(i));  cycle
      endif
      if( present(B) ) then
        call scan_contact_state( flag_ctAlgo, fstrSOLID%contacts(i), fstrSOLID%ddunode(:), fstrSOLID%dunode(:), &
           & fstrSOLID%QFORCE(:), infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive, mu, B )
      else
        call scan_contact_state( flag_ctAlgo, fstrSOLID%contacts(i), fstrSOLID%ddunode(:), fstrSOLID%dunode(:), &
           & fstrSOLID%QFORCE(:), infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive, mu )
      endif
      if( .not. active ) active = iactive
    enddo

    do i=1,fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) then
        call clear_contact_state(fstrSOLID%embeds(i));  cycle
      endif
      call scan_embed_state( flag_ctAlgo, fstrSOLID%embeds(i), fstrSOLID%ddunode(:), fstrSOLID%dunode(:), &
          & fstrSOLID%QFORCE(:), infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive, mu )
      if( .not. active ) active = iactive
    enddo

    if( is_init .and. ctAlgo == kcaSLAGRANGE .and. fstrSOLID%n_contacts > 0 ) &
      &  call remove_duplication_tiedcontact( cstep, hecMESH, fstrSOLID, infoCTChange )

    !for output contact state
    do i=1,fstrSOLID%n_contacts
      call set_contact_state_vector( fstrSOLID%contacts(i), dt, fstrSOLID%CONT_RELVEL, fstrSOLID%CONT_STATE )
    enddo

    infoCTChange%contactNode_current = infoCTChange%contactNode_previous+infoCTChange%free2contact-infoCTChange%contact2free
    infoCTChange%contactNode_previous = infoCTChange%contactNode_current

    if( .not. active ) then
      if( associated( fstrSOLID%CONT_NFORCE ) ) fstrSOLID%CONT_NFORCE(:) = 0.d0
      if( associated( fstrSOLID%CONT_FRIC ) ) fstrSOLID%CONT_FRIC(:) = 0.d0
    end if

  end subroutine

  !> Scanning contact state
  subroutine remove_duplication_tiedcontact( cstep, hecMESH, fstrSOLID, infoCTChange )
    integer(kind=kint), intent(in)         :: cstep      !< current step number
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
    type(fstr_info_contactChange), intent(inout):: infoCTChange   !<

    integer(kind=kint) :: i, j, grpid, slave
    integer(kind=kint) :: k, id, iSS
    integer(kind=kint) :: ig0, ig, iS0, iE0
    integer(kind=kint), allocatable :: states(:)

    allocate(states(hecMESH%n_node))
    states(:) = CONTACTFREE

    ! if a boundary condition is given, the slave
    do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
      grpid = fstrSOLID%BOUNDARY_ngrp_GRPID(ig0)
      if( .not. fstr_isBoundaryActive( fstrSOLID, grpid, cstep ) ) cycle
      ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
      iS0= hecMESH%node_group%grp_index(ig-1) + 1
      iE0= hecMESH%node_group%grp_index(ig  )
      do k= iS0, iE0
        iSS = hecMESH%node_group%grp_item(k)
        !states(iSS) = CONTACTSTICK
      enddo
    enddo

    do i=1,fstrSOLID%n_contacts
      if( fstrSOLID%contacts(i)%algtype /= CONTACTTIED ) cycle
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      do j=1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state==CONTACTFREE ) cycle   ! free
        slave = fstrSOLID%contacts(i)%slave(j)
        if( states(slave) == CONTACTFREE ) then
          states(slave) = fstrSOLID%contacts(i)%states(j)%state
          id = fstrSOLID%contacts(i)%states(j)%surface
          do k=1,size( fstrSOLID%contacts(i)%master(id)%nodes )
            iSS = fstrSOLID%contacts(i)%master(id)%nodes(k)
            states(iSS) = fstrSOLID%contacts(i)%states(j)%state
          enddo
        else !found duplicate tied contact slave node
          fstrSOLID%contacts(i)%states(j)%state = CONTACTFREE
          infoCTChange%free2contact = infoCTChange%free2contact - 1
          write(*,'(A,i10,A,i6,A,i6,A)') "Node",hecMESH%global_node_ID(slave), &
            " in rank",hecmw_comm_get_rank()," freed due to duplication"
        endif
      enddo
    enddo

  end subroutine


  !> Scanning contact state
  subroutine fstr_scan_contact_state_exp( cstep, hecMESH, fstrSOLID, infoCTChange )
    integer(kind=kint), intent(in)               :: cstep         !< current step number
    type( hecmwST_local_mesh ), intent(in)       :: hecMESH       !< type mesh
    type(fstr_solid), intent(inout)              :: fstrSOLID     !< type fstr_solid
    type(fstr_info_contactChange), intent(inout) :: infoCTChange  !<

    integer(kind=kint) :: i
    logical :: iactive, is_init


    ! P.A. We redefine fstrSOLID%ddunode as current coordinate of every nodes
    !  fstrSOLID%ddunode(:) = fstrSOLID%unode(:) + fstrSOLID%dunode(:)
    do i = 1, size(fstrSOLID%unode)
      fstrSOLID%ddunode(i) = hecMESH%node(i) + fstrSOLID%unode(i) + fstrSOLID%dunode(i)
    enddo
    infoCTChange%active = .false.

    infoCTChange%contact2free = 0
    infoCTChange%contact2neighbor = 0
    infoCTChange%contact2diffLpos = 0
    infoCTChange%free2contact = 0
    infoCTChange%contactNode_current = 0

    is_init = ( cstep == 1 )

    do i=1,fstrSOLID%n_contacts
   !   grpid = fstrSOLID%contacts(i)%group
   !   if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) then
   !     call clear_contact_state(fstrSOLID%contacts(i));  cycle
   !   endif

      call scan_contact_state_exp( fstrSOLID%contacts(i), fstrSOLID%ddunode(:), fstrSOLID%dunode(:), &
           & infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive )

      infoCTChange%active = infoCTChange%active .or. iactive
    enddo

    infoCTChange%contactNode_current = infoCTChange%contactNode_previous+infoCTChange%free2contact-infoCTChange%contact2free
    infoCTChange%contactNode_previous = infoCTChange%contactNode_current
    fstrSOLID%ddunode = 0.d0
  end subroutine

  !> Update lagrangian multiplier
  subroutine fstr_update_contact0( hecMESH, fstrSOLID, B )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH     !< type mesh
    type(fstr_solid), intent(inout)        :: fstrSOLID   !< type fstr_solid
    real(kind=kreal), intent(inout)        :: B(:)        !< nodal force residual

    integer(kind=kint) :: i, algtype

    do i=1, fstrSOLID%n_contacts
      !   if( contacts(i)%mpced ) cycle
      algtype = fstrSOLID%contacts(i)%algtype
      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        call calcu_contact_force0( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
        , fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, mu, mut, B )
      else if( algtype == CONTACTTIED ) then
        call calcu_tied_force0( fstrSOLID%contacts(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), mu, B )
      endif
    enddo

    do i=1, fstrSOLID%n_embeds
      call calcu_tied_force0( fstrSOLID%embeds(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), mu, B )
    enddo

  end subroutine

  !> Update lagrangian multiplier
  subroutine fstr_update_contact_multiplier( hecMESH, fstrSOLID, ctchanged )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH
    type(fstr_solid), intent(inout)        :: fstrSOLID
    logical, intent(out)                   :: ctchanged

    integer(kind=kint) :: i, nc, algtype

    gnt = 0.d0;  ctchanged = .false.
    nc = fstrSOLID%n_contacts+fstrSOLID%n_embeds
    do i=1, fstrSOLID%n_contacts
      algtype = fstrSOLID%contacts(i)%algtype
      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        call update_contact_multiplier( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
        , fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, mu, mut, gnt, ctchanged )
      else if( algtype == CONTACTTIED ) then
        call update_tied_multiplier( fstrSOLID%contacts(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), &
                                      &  mu, ctchanged )
      endif
    enddo

    do i=1, fstrSOLID%n_embeds
      call update_tied_multiplier( fstrSOLID%embeds(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), &
      &  mu, ctchanged )
    enddo

    if( nc>0 ) gnt = gnt/nc
  end subroutine

  !> Update lagrangian multiplier
  subroutine fstr_ass_load_contactAlag( hecMESH, fstrSOLID, B )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH
    type(fstr_solid), intent(inout)        :: fstrSOLID
    real(kind=kreal), intent(inout)        :: B(:)        !< nodal force residual

    integer(kind=kint) :: i, algtype

    do i = 1, fstrSOLID%n_contacts
      algtype = fstrSOLID%contacts(i)%algtype
      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        call ass_contact_force( fstrSOLID%contacts(i), hecMESH%node, fstrSOLID%unode, B )
      else if( algtype == CONTACTTIED ) then
        call calcu_tied_force0( fstrSOLID%contacts(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), mu, B )
      endif
    enddo

    do i = 1, fstrSOLID%n_embeds
      call calcu_tied_force0( fstrSOLID%embeds(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), mu, B )
    enddo
  end subroutine

  !> Update tangent force
  subroutine fstr_update_contact_TangentForce( fstrSOLID )
    type(fstr_solid), intent(inout)        :: fstrSOLID

    integer(kind=kint) :: i

    do i=1, fstrSOLID%n_contacts
      call update_contact_TangentForce( fstrSOLID%contacts(i) )
    enddo
  end subroutine

  !> Introduce contact stiff into global stiff matrix or mpc conditions into hecMESH
  subroutine fstr_contactBC( cstep, iter, hecMESH, hecMAT, fstrSOLID )
    use fstr_ctrl_modifier
    integer(kind=kint)                       :: cstep     !< current loading step
    integer(kind=kint), intent(in)           :: iter      !< NR iterations
    type (hecmwST_local_mesh), intent(inout) :: hecMESH   !< type mesh
    type (hecmwST_matrix), intent(inout)     :: hecMAT    !< type matrix
    type(fstr_solid), intent(inout)          :: fstrSOLID !< type fstr_solid

    integer(kind=kint), parameter :: NDOF=3

    integer(kind=kint) :: i, j, k, m, nnode, nd, etype, grpid
    integer(kind=kint) :: ctsurf, ndLocal(l_max_surface_node+1)
    integer(kind=kint) :: algtype    
    real(kind=kreal) :: factor, elecoord( 3, l_max_surface_node)
    real(kind=kreal) :: stiff(l_max_surface_node*3+3, l_max_surface_node*3+3)
    real(kind=kreal) :: nrlforce, force(l_max_surface_node*3+3)

    factor = fstrSOLID%FACTOR(2)

    do i=1,fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      algtype = fstrSOLID%contacts(i)%algtype

      do j=1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state==CONTACTFREE ) cycle   ! free
        ctsurf = fstrSOLID%contacts(i)%states(j)%surface          ! contacting surface
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        do k=1,nnode
          ndLocal(k+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(k)
          elecoord(1:3,k)=hecMESH%node(3*ndLocal(k+1)-2:3*ndLocal(k+1))
        enddo
        if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
          call contact2stiff( algtype, fstrSOLID%contacts(i)%states(j),    &
          etype, nnode, elecoord(:,:), mu, mut, fstrSOLID%contacts(i)%fcoeff,    &
          fstrSOLID%contacts(i)%symmetric, stiff(:,:), force(:) )
        else if( algtype == CONTACTTIED ) then
          call tied2stiff( algtype, fstrSOLID%contacts(i)%states(j),    &
          etype, nnode, mu, mut, stiff(:,:), force(:) )
        endif
        ! ----- CONSTRUCT the GLOBAL MATRIX STARTED
        call hecmw_mat_ass_elem(hecMAT, nnode+1, ndLocal, stiff)

        if( iter>1 ) cycle
        !  if( fstrSOLID%contacts(i)%states(j)%multiplier(1)/=0.d0 ) cycle
        ! In case of new contact nodes, add enforced disp constraint
        fstrSOLID%contacts(i)%states(j)%wkdist = fstrSOLID%contacts(i)%states(j)%distance
        nrlforce = -mu*fstrSOLID%contacts(i)%states(j)%distance
        force(1:nnode*NDOF+NDOF) = force(1:nnode*NDOF+NDOF)*nrlforce
        do m=1,nnode+1
          nd = ndLocal(m)
          do k=1,NDOF
            hecMAT%B(NDOF*(nd-1)+k)=hecMAT%B(NDOF*(nd-1)+k)-force((m-1)*NDOF+k)
          enddo
        enddo

      enddo
    enddo

    do i=1,fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
      do j=1, size(fstrSOLID%embeds(i)%slave)
        if( fstrSOLID%embeds(i)%states(j)%state==CONTACTFREE ) cycle   ! free
        ctsurf = fstrSOLID%embeds(i)%states(j)%surface          ! contacting surface
        etype = fstrSOLID%embeds(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%embeds(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%embeds(i)%slave(j)
        do k=1,nnode
          ndLocal(k+1) = fstrSOLID%embeds(i)%master(ctsurf)%nodes(k)
          elecoord(1:3,k)=hecMESH%node(3*ndLocal(k+1)-2:3*ndLocal(k+1))
        enddo
        call tied2stiff( algtype, fstrSOLID%embeds(i)%states(j),    &
        etype, nnode, mu, mut, stiff(:,:), force(:) )
        ! ----- CONSTRUCT the GLOBAL MATRIX STARTED
        call hecmw_mat_ass_elem(hecMAT, nnode+1, ndLocal, stiff)
      enddo
    enddo

  end subroutine

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
    real(kind=kreal), pointer   :: coord(:)

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
    integer(kind=kint) :: i, isuf, icel, sid, etype, nn, iS, stype, idx
    integer(kind=kint) :: ndlocal(l_max_elem_node)
    real(kind=kreal), pointer   :: coord(:)
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
    real(kind=kreal) ::  rhsB
    integer(kind=kint) ::  i,j,N,i0,N_loc,nndof
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
      do i= hecMESH%nn_internal+1,hecMESH%n_node
        pid = hecMESH%node_ID(i*2)
        lid = hecMESH%node_ID(i*2-1)
        i0 = (displs(pid) + (lid-1))*ndof
        vec_all(i0+1:i0+ndof) = vec((i-1)*ndof+1:i*ndof)
        vec((i-1)*ndof+1:i*ndof) = 0.d0
      enddo
  
      call hecmw_allreduce_R(hecMESH, vec_all, N*ndof, hecmw_sum)
  
      do i=1,ndof*N_loc
        vec(i) = vec(i) + vec_all(offset*ndof+i)
      end do
    else if( vtype == 2 ) then
      vec_all(:) = -1000.d0
      do i= hecMESH%nn_internal+1,hecMESH%n_node
        if( vec(i) == 0.d0 ) cycle
        pid = hecMESH%node_ID(i*2)
        lid = hecMESH%node_ID(i*2-1)
        i0 = displs(pid) + lid
        vec_all(i0) = vec(i)
      enddo
  
      call hecmw_allreduce_R(hecMESH, vec_all, N, hecmw_max)
  
      do i=1,N_loc
        if( vec_all(offset+i) == -1000.d0 ) cycle
        if( vec(i) < vec_all(offset+i) ) vec(i) = vec_all(offset+i)
      end do
    end if
  
    deallocate(displs,vec_all)
    end subroutine

end module mContact
