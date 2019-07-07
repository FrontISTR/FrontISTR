!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provides functions to calcualte contact stiff matrix
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

  real(kind=kreal), save :: gnt(2)    !< 1:current avarage penetration;
  !< 2:current relative tangent displacement
  real(kind=kreal), save :: bakgnt(2) !< 1:current avarage penetration;
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
        ,pair%master_grp_id(i)
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
    if( ctAlgo== kcaSLagrange ) then
      if( infoCTChange%contact2free+infoCTChange%contact2neighbor+      &
        infoCTChange%contact2difflpos+infoCTChange%free2contact == 0 ) &
        fstr_is_contact_conv = .true.
    elseif( ctAlgo==kcaALagrange ) then
      if( gnt(1)<cgn .and. gnt(2)<cgt ) fstr_is_contact_conv = .true.
      write(*,'(a,2e15.7)') "--Relative displacement in contact surface",gnt
      !  if( dabs( bakgnt(1)-gnt(1) )<cgn*1.d-2 .and.    &
        !      dabs( bakgnt(2)-gnt(2) )<cgn ) fstr_is_contact_conv = .true.
      bakgnt = gnt
    endif
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

    fstrSOLID%CONT_RELVEL(:) = 0.d0
    fstrSOLID%CONT_STATE(:) = 0.d0

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

    do i=1,size(fstrSOLID%contacts)
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

        !for output contact state
        call set_contact_state_vector( fstrSOLID%contacts(i), dt, fstrSOLID%CONT_RELVEL, fstrSOLID%CONT_STATE )
    enddo

    infoCTChange%contactNode_current = infoCTChange%contactNode_previous+infoCTChange%free2contact-infoCTChange%contact2free
    infoCTChange%contactNode_previous = infoCTChange%contactNode_current

  end subroutine

  !> Scanning contact state
  subroutine fstr_scan_contact_state_exp( cstep, hecMESH, fstrSOLID, infoCTChange )
    integer(kind=kint), intent(in)               :: cstep         !< current step number
    type( hecmwST_local_mesh ), intent(in)       :: hecMESH       !< type mesh
    type(fstr_solid), intent(inout)              :: fstrSOLID     !< type fstr_solid
    type(fstr_info_contactChange), intent(inout) :: infoCTChange  !<

    integer(kind=kint) :: i, grpid
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

    do i=1,size(fstrSOLID%contacts)
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

    integer(kind=kint) :: i
    do i=1, size(fstrSOLID%contacts)
      !   if( contacts(i)%mpced ) cycle
      call calcu_contact_force0( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
        , fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, mu, mut, B )
    enddo
  end subroutine

  !> Update lagrangian multiplier
  subroutine fstr_update_contact_multiplier( hecMESH, fstrSOLID, ctchanged )
    type( hecmwST_local_mesh ), intent(in) :: hecMESH
    type(fstr_solid), intent(inout)        :: fstrSOLID
    logical, intent(out)                   :: ctchanged

    integer(kind=kint) :: i, nc
    gnt = 0.d0;  ctchanged = .false.
    nc = size(fstrSOLID%contacts)
    do i=1, nc
      !   if( contacts(i)%mpced ) cycle
      call update_contact_multiplier( fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:)  &
        , fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, mu, mut, gnt, ctchanged )
    enddo
    if( nc>0 ) gnt = gnt/nc
  end subroutine

  !> Update tangent force
  subroutine fstr_update_contact_TangentForce( fstrSOLID )
    type(fstr_solid), intent(inout)        :: fstrSOLID

    integer(kind=kint) :: i, nc

    nc = size(fstrSOLID%contacts)
    do i=1, nc
      call update_contact_TangentForce( fstrSOLID%contacts(i) )
    enddo
  end subroutine

  !> Introduce contact stiff into global stiff matrix or mpc conditions into hecMESH
  subroutine fstr_contactBC( iter, hecMESH, hecMAT, fstrSOLID )
    use fstr_ctrl_modifier
    integer(kind=kint), intent(in)           :: iter      !< NR iterations
    type (hecmwST_local_mesh), intent(inout) :: hecMESH   !< type mesh
    type (hecmwST_matrix), intent(inout)     :: hecMAT    !< type matrix
    type(fstr_solid), intent(inout)          :: fstrSOLID !< type fstr_solid

    integer(kind=kint), parameter :: NDOF=3

    integer(kind=kint) :: i, j, k, m, nnode, nd, etype
    integer(kind=kint) :: ctsurf, ndLocal(l_max_surface_node+1)
    real(kind=kreal) :: factor, elecoord( 3, l_max_surface_node)
    real(kind=kreal) :: stiff(l_max_surface_node*3+3, l_max_surface_node*3+3)
    real(kind=kreal) :: nrlforce, force(l_max_surface_node*3+3)
    !   if( n_contact_mpc>0 ) call fstr_delete_mpc( n_contact_mpc, hecMESH%mpc )
    !   call fstr_contact2mpc( fstrSOLID%contacts, hecMESH%mpc )
    ! temp. need modification
    !   do i=1,size(fstrSOLID%contacts)
    !     fstrSOLID%contacts(i)%mpced = .true.
    !     enddo
    !    call fstr_write_mpc( 6, hecMESH%mpc )
    !print *,"Contact to mpc ok!",n_contact_mpc,"equations generated"
    factor = fstrSOLID%FACTOR(2)
    call hecmw_cmat_clear( hecMAT%cmat )
    do i=1,size(fstrSOLID%contacts)
      do j=1, size(fstrSOLID%contacts(i)%slave)
        if( fstrSOLID%contacts(i)%states(j)%state==CONTACTFREE ) cycle   ! free
        ctsurf = fstrSOLID%contacts(i)%states(j)%surface          ! contacting surface
        etype = fstrSOLID%contacts(i)%master(ctsurf)%etype
        nnode = size(fstrSOLID%contacts(i)%master(ctsurf)%nodes)
        ndLocal(1) = fstrSOLID%contacts(i)%slave(j)
        do k=1,nnode
          ndLocal(k+1) = fstrSOLID%contacts(i)%master(ctsurf)%nodes(k)
          elecoord(1,k)=hecMESH%node(3*ndLocal(k+1)-2)
          elecoord(2,k)=hecMESH%node(3*ndLocal(k+1)-1)
          elecoord(3,k)=hecMESH%node(3*ndLocal(k+1))
        enddo
        call contact2stiff( fstrSOLID%contacts(i)%algtype, fstrSOLID%contacts(i)%states(j),    &
          etype, nnode, elecoord(:,:), mu, mut, fstrSOLID%contacts(i)%fcoeff,    &
          fstrSOLID%contacts(i)%symmetric, stiff(:,:), force(:) )
        call hecmw_mat_ass_contact( hecMAT,nnode+1,ndLocal,stiff )

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
  end subroutine


end module mContact
