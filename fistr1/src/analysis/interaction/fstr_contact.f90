!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief Top-level contact analysis module (System level)
!>
!> Contact Analysis Module Hierarchy:
!> ===================================
!>
!> Level 1: Element Level (1 slave-master pair)
!>   m_fstr_contact_element:
!>     - Computes stiffness matrix and force vector for a single contact pair
!>     - Mathematical calculations: Tm/Tt matrices, friction state, etc.
!>
!> Level 2: Contact Object Level (all pairs in one tContact)
!>   m_fstr_contact_assembly:
!>     - Processes all slave-master pairs in one tContact object
!>     - Updates contact state (multipliers, tangent forces)
!>     - Assembles element matrices into global system
!>
!>   m_fstr_contact_search:
!>     - Detects contact state for all pairs in one tContact object
!>     - Functions: scan_contact_state, track_contact_position
!>
!> Level 3: System Level (all contact objects)
!>   mContact (this module):
!>     - Manages all contact objects in fstrSOLID%contacts(:)
!>     - Core processing and state management
!>
!> Supporting modules:
!>   m_fstr_contact_mpc:    MPC (Multi-Point Constraint) processing
!>   m_fstr_contact_output: Output vector initialization and processing
module mContact

  use mContactDef
  use hecmw
  use m_fstr
  use m_fstr_contact_assembly
  use m_fstr_contact_element
  use m_fstr_contact_geom
  use m_fstr_contact_search
  use m_fstr_contact_mpc
  use m_fstr_contact_output
  implicit none

  logical, private :: active

  ! Re-export functions from sub-modules
  public :: initialize_contact_output_vectors
  public :: initialize_embed_vectors
  public :: setup_contact_elesurf_for_area
  public :: calc_contact_area
  public :: fstr_setup_parancon_contactvalue
  public :: update_contact_state_vectors
  public :: fstr_update_contact_state_vectors

contains

  subroutine fstr_AddContactStiffness(cstep,ctAlgo, iter,hecMESH,conMAT,hecLagMAT,fstrSOLID)

    integer(kind=kint)                   :: cstep !< current loading step
    integer(kind=kint)                   :: ctAlgo !< contact algorithm type
    integer(kind=kint)                   :: iter
    type(hecmwST_local_mesh)             :: hecMESH !< type hecmwST_local_mesh
    type(hecmwST_matrix)                 :: conMAT !< type hecmwST_matrix
    type(fstr_solid)                     :: fstrSOLID !< type fstr_solid
    type(hecmwST_matrix_lagrange)        :: hecLagMAT !< type hecmwST_matrix_lagrange
    integer(kind=kint)                   :: i, grpid

    if( associated(hecLagMAT%AL_lagrange) ) hecLagMAT%AL_lagrange = 0.0d0
    if( associated(hecLagMAT%AU_lagrange) ) hecLagMAT%AU_lagrange = 0.0d0

    do i = 1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_stiffness_NodeSurf( ctAlgo, fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:), &
        iter, hecLagMAT%Lagrange(:), conMAT, hecLagMAT)

    enddo

    do i = 1, fstrSOLID%n_embeds

      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle

      call calcu_contact_stiffness_NodeSurf( ctAlgo, fstrSOLID%embeds(i), hecMESH%node(:), fstrSOLID%unode(:), &
        iter, hecLagMAT%Lagrange(:), conMAT, hecLagMAT)

    enddo

  end subroutine fstr_AddContactStiffness

  !> \brief Compute contact forces for residual vector (conMAT%B).
  subroutine fstr_Update_NDForce_contact(cstep,ctAlgo,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)
    integer(kind=kint), intent(in)       :: cstep
    integer(kind=kint), intent(in)       :: ctAlgo
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_matrix_lagrange)        :: hecLagMAT
    type(hecmwST_matrix)                 :: conMAT

    call fstr_contact_ndforce_core(kctForResidual,cstep,ctAlgo,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)

  end subroutine fstr_Update_NDForce_contact

  !> \brief Compute contact forces for output (CONT_NFORCE/CONT_FRIC).
  subroutine fstr_calc_contact_output_force(cstep,ctAlgo,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)
    integer(kind=kint), intent(in)       :: cstep
    integer(kind=kint), intent(in)       :: ctAlgo
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_matrix_lagrange)        :: hecLagMAT
    type(hecmwST_matrix)                 :: conMAT

    call fstr_contact_ndforce_core(kctForOutput,cstep,ctAlgo,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)

  end subroutine fstr_calc_contact_output_force

  !> \brief Core routine: compute contact nodal forces for all contact/embed pairs.
  !! purpose == kctForResidual: assemble into conMAT%B
  !! purpose == kctForOutput:   store into CONT_NFORCE/CONT_FRIC (zero-cleared first)
  subroutine fstr_contact_ndforce_core(purpose,cstep,ctAlgo,hecMESH,hecMAT,hecLagMAT,fstrSOLID,conMAT)
    integer(kind=kint), intent(in)       :: purpose
    integer(kind=kint), intent(in)       :: cstep
    integer(kind=kint), intent(in)       :: ctAlgo
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid)                     :: fstrSOLID
    type(hecmwST_matrix_lagrange)        :: hecLagMAT
    type(hecmwST_matrix)                 :: conMAT
    integer(kind=kint) :: i, grpid

    ! Zero-clear output arrays only when computing for output
    if( purpose == kctForOutput ) then
      if( associated(fstrSOLID%CONT_NFORCE) ) fstrSOLID%CONT_NFORCE(:) = 0.d0
      if( associated(fstrSOLID%CONT_FRIC) )   fstrSOLID%CONT_FRIC(:) = 0.d0
      if( associated(fstrSOLID%EMBED_NFORCE) ) fstrSOLID%EMBED_NFORCE(:) = 0.d0
    endif

    do i = 1, fstrSOLID%n_contacts
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle
      call calcu_contact_ndforce_NodeSurf( purpose, ctAlgo, fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:), &
        fstrSOLID%dunode(:), hecLagMAT%Lagrange(:), conMAT, &
        fstrSOLID%CONT_NFORCE, fstrSOLID%CONT_FRIC, hecLagMAT )
    enddo

    do i = 1, fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) cycle
      call calcu_contact_ndforce_NodeSurf( purpose, ctAlgo, fstrSOLID%embeds(i), hecMESH%node(:), fstrSOLID%unode(:), &
        fstrSOLID%dunode(:), hecLagMAT%Lagrange(:), conMAT, &
        fstrSOLID%EMBED_NFORCE, fstrSOLID%EMBED_NFORCE, hecLagMAT )
    enddo

  end subroutine fstr_contact_ndforce_core

  !> Calculate reference stiffness for all contact pairs (System Level)
  subroutine fstr_calc_contact_refStiff(cstep, hecMESH, hecMAT, fstrSOLID)
    use hecmw_matrix_misc
    integer(kind=kint), intent(in)       :: cstep
    type(hecmwST_local_mesh)             :: hecMESH
    type(hecmwST_matrix)                 :: hecMAT
    type(fstr_solid)                     :: fstrSOLID
    
    integer(kind=kint) :: i, grpid, ndof
    real(kind=kreal), pointer :: diag(:)
    
    ndof = hecMAT%NDOF
    
    ! Extract diagonal components once (efficient, called only once)
    diag => hecmw_mat_diag(hecMAT)
    
    ! Loop over contact pairs
    do i = 1, fstrSOLID%n_contacts
      grpid = fstrSOLID%contacts(i)%group
      if(.not. fstr_isContactActive(fstrSOLID, grpid, cstep)) cycle
      
      call calc_contact_pair_refStiff(fstrSOLID%contacts(i), diag, ndof, hecMESH)
    enddo
    
    ! Same for embeds
    do i = 1, fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if(.not. fstr_isEmbedActive(fstrSOLID, grpid, cstep)) cycle
      
      call calc_contact_pair_refStiff(fstrSOLID%embeds(i), diag, ndof, hecMESH)
    enddo
    
    deallocate(diag)
  end subroutine fstr_calc_contact_refStiff

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
        & fstrSOLID%QFORCE(:), infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive, B )
      else
        call scan_contact_state( flag_ctAlgo, fstrSOLID%contacts(i), fstrSOLID%ddunode(:), fstrSOLID%dunode(:), &
        & fstrSOLID%QFORCE(:), infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive )
      endif
      if( .not. active ) active = iactive
    enddo

    do i=1,fstrSOLID%n_embeds
      grpid = fstrSOLID%embeds(i)%group
      if( .not. fstr_isEmbedActive( fstrSOLID, grpid, cstep ) ) then
        call clear_contact_state(fstrSOLID%embeds(i));  cycle
      endif
      call scan_embed_state( flag_ctAlgo, fstrSOLID%embeds(i), fstrSOLID%ddunode(:), fstrSOLID%dunode(:), &
      & fstrSOLID%QFORCE(:), infoCTChange, hecMESH%global_node_ID(:), hecMESH%global_elem_ID(:), is_init, iactive )
      if( .not. active ) active = iactive
    enddo

    if( is_init .and. ctAlgo == kcaSLAGRANGE .and. fstrSOLID%n_contacts > 0 ) &
    &  call remove_duplication_tiedcontact( cstep, hecMESH, fstrSOLID, infoCTChange )

    infoCTChange%contactNode_current = infoCTChange%contactNode_previous+infoCTChange%free2contact-infoCTChange%contact2free
    infoCTChange%contactNode_previous = infoCTChange%contactNode_current

    if( .not. active ) then
      if( associated( fstrSOLID%CONT_NFORCE ) ) fstrSOLID%CONT_NFORCE(:) = 0.d0
      if( associated( fstrSOLID%CONT_FRIC ) ) fstrSOLID%CONT_FRIC(:) = 0.d0
    end if

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

    fstr_is_contact_conv = .false.
    if( infoCTChange%contact2free+infoCTChange%contact2neighbor+      &
      infoCTChange%contact2difflpos+infoCTChange%free2contact == 0 ) &
      fstr_is_contact_conv = .true.

    call hecmw_allreduce_L1(hecMESH, fstr_is_contact_conv, HECMW_LAND)
  end function

  logical function fstr_is_matrixStructure_changed(infoCTChange)
    type (fstr_info_contactChange)   :: infoCTChange  !< fstr_contactChange
    fstr_is_matrixStructure_changed = .false.
    if( infoCTChange%contact2free+infoCTChange%contact2neighbor+infoCTChange%free2contact > 0 ) &
      fstr_is_matrixStructure_changed = .true.
  end function

  subroutine fstr_update_contact_multiplier( cstep, ctAlgo, hecMESH, hecLagMAT, fstrSOLID, ctchanged )
    integer(kind=kint), intent(in)        :: cstep
    integer(kind=kint), intent(in)        :: ctAlgo
    type( hecmwST_local_mesh ), intent(in) :: hecMESH
    type(hecmwST_matrix_lagrange), intent(in) :: hecLagMAT
    type(fstr_solid), intent(inout)        :: fstrSOLID
    logical, intent(out)                   :: ctchanged

    integer(kind=kint) :: i, nc, algtype, grpid

    gnt = 0.d0;  ctchanged = .false.
    nc = fstrSOLID%n_contacts+fstrSOLID%n_embeds
    do i=1, fstrSOLID%n_contacts
      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      algtype = fstrSOLID%contacts(i)%algtype
      if( algtype == CONTACTSSLID .or. algtype == CONTACTFSLID ) then
        call update_contact_multiplier( ctAlgo, fstrSOLID%contacts(i), hecMESH%node(:), fstrSOLID%unode(:), &
          fstrSOLID%dunode(:), fstrSOLID%contacts(i)%fcoeff, &
          hecMESH, hecLagMAT, gnt, ctchanged )
      else if( algtype == CONTACTTIED ) then
        call update_tied_multiplier( fstrSOLID%contacts(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), &
        &  ctchanged )
      endif
    enddo

    do i=1, fstrSOLID%n_embeds
      call update_tied_multiplier( fstrSOLID%embeds(i), fstrSOLID%unode(:), fstrSOLID%dunode(:), &
      &  ctchanged )
    enddo

    call hecmw_allreduce_L1(hecMESH, ctchanged, HECMW_LOR)

    if( nc>0 ) gnt = gnt/nc
  end subroutine

  !> Update tangent force
  subroutine fstr_update_contact_TangentForce( cstep, fstrSOLID )
    integer(kind=kint), intent(in)        :: cstep
    type(fstr_solid), intent(inout)        :: fstrSOLID

    integer(kind=kint) :: i, grpid

    do i=1, fstrSOLID%n_contacts

      grpid = fstrSOLID%contacts(i)%group
      if( .not. fstr_isContactActive( fstrSOLID, grpid, cstep ) ) cycle

      call update_contact_TangentForce( fstrSOLID%contacts(i) )
    enddo
  end subroutine

  !> Update contact state output vectors for all contacts
  subroutine fstr_update_contact_state_vectors( fstrSOLID, dt )
    type(fstr_solid), intent(inout)        :: fstrSOLID
    real(kind=kreal), intent(in)           :: dt  !< time increment

    integer(kind=kint) :: i

    if( associated( fstrSOLID%CONT_RELVEL ) ) fstrSOLID%CONT_RELVEL(:) = 0.d0
    if( associated( fstrSOLID%CONT_STATE ) ) fstrSOLID%CONT_STATE(:) = 0.d0

    do i=1, fstrSOLID%n_contacts
      call update_contact_state_vectors( fstrSOLID%contacts(i), dt, fstrSOLID%CONT_RELVEL, fstrSOLID%CONT_STATE )
    enddo
  end subroutine

end module mContact
