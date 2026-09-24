!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief  This module manages the data structure for contact calculation
!!
!!  Contact calculation takes into act after calling the following three
!!  subrotuines provided in this module
!!-#      Reading contact definition with subroutine: fstr_ctrl_get_CONTACT
!!-#      Check its consistency with mesh definition: fstr_contact_check
!!-#      Initializing the contact calculation      : fstr_contact_init
module mContactDef

  use hecmw
  use elementInfo
  use mSurfElement
  use m_hecmw_contact_comm
  use bucket_search
  use mContactParam

  implicit none

  real(kind=kreal), save :: cgn=1.d-5 !< convergent condition of penetration
  real(kind=kreal), save :: cgt=1.d-3 !< convergent condition of relative tangent disp

  real(kind=kreal), save :: gnt(2)    !< 1:current average penetration;
  !< 2:current relative tangent displacement
  real(kind=kreal), save :: bakgnt(2) !< 1:current average penetration;
  !< 2:current relative tangent displacement!

  integer, parameter :: CONTACTUNKNOWN = -1
  !> contact state definition
  integer, parameter :: CONTACTFREE = -1
  integer, parameter :: CONTACTNEAR = 0    !< near contact: projection info available, no LM constraint
  integer, parameter :: CONTACTSTICK = 1
  integer, parameter :: CONTACTSLIP = 2
  integer, parameter :: CANDIDATE_INTP = 3  !< SURF-SURF: integration point to be re-projected in this scan

  !> contact state category: states sharing a category contribute to the matrix in the same way
  integer, parameter :: kcatFREE = 1  !< no projection
  integer, parameter :: kcatNEAR = 2  !< projection only, no Lagrange multiplier
  integer, parameter :: kcatCONT = 3  !< Lagrange multiplier held (STICK or SLIP)

  !> contact type or algorithm definition
  integer, parameter :: CONTACTTIED = 1
  integer, parameter :: CONTACTGLUED = 2
  integer, parameter :: CONTACTSSLID = 3
  integer, parameter :: CONTACTFSLID = 4

  !> contact smoothing type
  integer, parameter :: kcsNONE   = 0
  integer, parameter :: kcsNAGATA = 1

  !> contact method
  integer, parameter :: CONTACTN2S = 1
  integer, parameter :: CONTACTS2S = 2

  !> sparsity expansion mode (per contact pair, set by !CONTACT EXPANSION=)
  integer, parameter :: SPARSITY_NONE     = 0  !< register current master only
  integer, parameter :: SPARSITY_NEIGHBOR = 1  !< register current master + neighbors (fewer matrix rebuilds)

  !> upper bound of integration points per SURF-SURF slave segment
  integer, parameter :: MAX_N_INTP = 128

  !> contact interference type
  integer, parameter :: C_IF_SLAVE = 1
  integer, parameter :: C_IF_MASTER = 2

  !> purpose flag for contact force calculation
  integer, parameter :: kctForResidual = 1  !< compute contact force for residual (conMAT%B)
  integer, parameter :: kctForOutput   = 2  !< compute contact force for output (CONT_NFORCE/CONT_FRIC)

  !> This structure records contact status
  type tContactState
    integer          :: state !< -1:free, 1:in contact, or other needed
    integer          :: surface !< contacting surface number
    real(kind=kreal) :: distance !< penetration value
    real(kind=kreal) :: wkdist !< copy of penetration value
    real(kind=kreal) :: lpos(3) !< contact position(local coordinate)
    real(kind=kreal) :: gpos(3) !< contact position(global coordinate)
    real(kind=kreal) :: direction(3) !< contact direction
    real(kind=kreal) :: multiplier(3) !< Lagrangian multiplier or contact force
    !< 1: normal 2:tangent component
    real(kind=kreal) :: tangentForce(3) !< friction force
    real(kind=kreal) :: tangentForce1(3) !< friction force rotated by element(for trial friction force)
    real(kind=kreal) :: tangentForce_trial(3) !< trial friction force
    real(kind=kreal) :: tangentForce_final(3) !< final friction force
    real(kind=kreal)    :: reldisp(3)
    !
    real(kind=kreal)    :: shrink_factor
    real(kind=kreal)    :: time_factor
    real(kind=kreal)    :: init_pos
    real(kind=kreal)    :: end_pos
    integer             :: interference_flag
  end type

  !> Structure to define a slave surface segment of a SURF-SURF (mortar) contact pair
  type tContactSurf
    integer(kind=kint)              :: eid                  !< elemental index(global)
    integer(kind=kint)              :: etype                !< type of surface element
    integer(kind=kint), pointer     :: nodes(:)=>null()     !< nodes index(global)
    integer(kind=kint)              :: state = CONTACTFREE       !< segment contact state (CONTACTFREE until first scan)
    integer(kind=kint)              :: state_prev = CONTACTFREE  !< previous scan's segment state (for seg on/off detection)
    integer(kind=kint)              :: n_intp               !< num of surface integral point
    type(tContactState), pointer    :: states(:)=>null()    !< contact states of slave surf
    integer(kind=kint), pointer     :: nslave_index(:)=>null() !< mapping from surface node to slave node index (surf-surf)
    ! --- lambda transaction buffers ---
    ! begin: previous substep's COMMIT product, immutable during a substep (warm-start source)
    ! working: current substep's accumulator, cleared at BEGIN, written by the augmentation update
    integer(kind=kint), pointer     :: lam_begin_id(:)=>null()  !< begin masterID list, ascending, size n_intp
    real(kind=kreal),   pointer     :: lam_begin_val(:,:)=>null() !< begin lambda_n (node a, rank r) keyed by lam_begin_id
    integer(kind=kint)              :: lam_begin_n = 0          !< begin valid count
    integer(kind=kint), pointer     :: lam_work_id(:)=>null()   !< working masterID list, ascending, size n_intp
    real(kind=kreal),   pointer     :: lam_work_val(:,:)=>null()  !< working lambda_n (node a, rank r) keyed by lam_work_id
    integer(kind=kint)              :: lam_work_n = 0           !< working valid count
    ! --- tangential friction parallel arrays ---
    ! per-node covariant tangent multiplier (2 components, node a, rank r) and per-node
    ! stick/slip state. Rank r is keyed by lam_*_id (master), node a is the slave-surf node;
    ! both mirror lam_*_val.
    real(kind=kreal),   pointer     :: lam_begin_t(:,:,:)=>null()   !< begin lambda_t (2, node a, rank r) keyed by lam_begin_id
    real(kind=kreal),   pointer     :: lam_work_t(:,:,:)=>null()    !< working lambda_t (2, node a, rank r) keyed by lam_work_id
    integer(kind=kint), pointer     :: lam_begin_fstate(:,:)=>null()!< begin per-node friction state, (node a, rank r)
    integer(kind=kint), pointer     :: lam_work_fstate(:,:)=>null() !< working per-node friction state, (node a, rank r)
    ! segment state carried in the same transaction (cutback does not restore slave_surf)
    integer(kind=kint)              :: state_begin = CONTACTFREE       !< committed segment state
    integer(kind=kint)              :: state_prev_begin = CONTACTFREE  !< committed previous segment state
  end type tContactSurf

  !> Structure to includes all info needed by contact calculation
  type tContact
    ! following contact definition
    character(len=HECMW_NAME_LEN) :: name                    !< name
    integer                       :: ctype                   !< 1:node-surface 2: surface-surface
    integer                       :: group                   !< group number
    character(len=HECMW_NAME_LEN) :: pair_name               !< name of contact pair
    integer                       :: surf_id1, surf_id2      !< slave surface, master surface
    integer                       :: surf_id1_sgrp           !< surface group id of slave surface
    type(tSurfElement), pointer   :: master(:)=>null()       !< master surface (element )
    integer                       :: n_master_owned = 0      !< owned-only master face count (SURF-SURF visibility guard)
    integer, pointer              :: slave(:)=>null()        !< slave surface (node)
    real(kind=kreal)              :: fcoeff                  !< coeeficient of friction
    real(kind=kreal)              :: nPenalty                !< normal penalty coefficient
    real(kind=kreal)              :: tPenalty                !< tangential penalty coefficient
    real(kind=kreal)              :: refStiff                !< reference stiffness for penalty calculation
    real(kind=kreal)              :: damp_alpha              !< damping coefficient (dimensionless, scaled by refStiff)
    real(kind=kreal)              :: damp_gact               !< damping activation distance [length] (<=0: disabled)

    type(tContactSurf), pointer   :: slave_surf(:)=>null()   !< slave surface segments (MORTAR=YES only)

    ! !CONTACT_INTERFERENCE data; default-initialized because check_apply_Contact_IF only
    ! writes them when that card is present, while if_type is read unconditionally
    ! (if_flag = contact%if_type /= 0) by the contact force / search paths.
    real(kind=kreal)    :: ctime = 0.d0
    integer(kind=kint)  :: if_type = 0
    real(kind=kreal)    :: if_etime = 0.d0
    real(kind=kreal)    :: initial_pos = 0.d0
    real(kind=kreal)    :: end_pos = 0.d0
    ! following algorithm
    ! -1: not initialized
    ! 1: TIED-Just rigidly fixed the two surfaces
    ! 2: GLUED-Distance between the two surfaces to zero and glue them
    ! 3: SSLID-Small sliding contact( no position but with contact state change)
    ! 4: FSLID-Finite sliding contact (both changes in contact state and position possible)
    integer                       :: algtype                 !< algorithm flag
    integer                       :: smoothing               !< kcsNONE or kcsNAGATA
    integer                       :: method = CONTACTN2S     !< CONTACTN2S (NODE-SURF) or CONTACTS2S (SURF-SURF mortar).
                                                             !< !EMBED has no MORTAR option, so it keeps this default.
    integer                       :: sparsity_expansion = SPARSITY_NONE  !< SPARSITY_NONE / SPARSITY_NEIGHBOR (!CONTACT EXPANSION=)

    logical                       :: mpced                   !< if turns into mpc condition
    logical                       :: symmetric               !< true for FRICTION_CONE=FROZEN: cone radius kept at the multiplier
    real(kind=kreal)              :: eps_fric_band = 0.d0    !< hysteresis half-band of the stick/slip switch (0 = no band)

    ! following contact state
    type(tContactState), pointer  :: states(:)=>null()       !< contact states of each slave nodes

    type(hecmwST_contact_comm)    :: comm                    !< contact communication table
    type(bucketDB)                :: master_bktDB            !< bucket DB for master surface

    type(tContactParam), pointer  :: cparam=>null()          !< contact parameter
  end type tContact

  type fstr_info_contactChange
    logical            :: active
    integer(kind=kint) :: n_statechange(3,3)     !< counter: slave nodes moved from a state category to another
    integer(kind=kint) :: contact2neighbor       !< counter: contact to neighbor state change (within 1-hop)
    integer(kind=kint) :: contact2beyond         !< counter: contact moved beyond neighbor elements
    integer(kind=kint) :: contact2diffLpos       !< counter: contact to different local position state change (NODE-SURF only)
    integer(kind=kint) :: free2contact_new       !< counter: free to contact per SURF-SURF segment (NEIGHBOR: new master only)
    integer(kind=kint) :: contactNode_previous   !< previous number of nodes in contact
    integer(kind=kint) :: contactNode_current    !< current number of nodes in contact
  end type fstr_info_contactChange

  private :: is_active_contact

contains

  !> Initializer
  subroutine contact_state_init(cstate)
    type(tContactState), intent(inout) :: cstate !< contact state
    cstate%state = -1
    cstate%surface = -1
    cstate%distance = 0.0d0
    cstate%wkdist = 0.0d0
    cstate%lpos(:) = 0.0d0
    cstate%gpos(:) = 0.0d0
    cstate%direction(:) = 0.0d0
    cstate%multiplier(:) = 0.0d0
    cstate%tangentForce(:) = 0.0d0
    cstate%tangentForce1(:) = 0.0d0
    cstate%tangentForce_trial(:) = 0.0d0
    cstate%tangentForce_final(:) = 0.0d0
    cstate%reldisp(:) = 0.0d0
    cstate%shrink_factor = 0.0d0
    cstate%time_factor = 0.0d0
    cstate%init_pos = 0.0d0
    cstate%end_pos = 0.0d0
    cstate%interference_flag = 0
  end subroutine

  !> Copy
  subroutine contact_state_copy(cstate1, cstate2)
    type(tContactState), intent(in)    :: cstate1 !< contact state
    type(tContactState), intent(inout) :: cstate2 !< contact state
    cstate2 = cstate1
  end subroutine

  !> Whether the contact state has active LM constraint (STICK or SLIP)
  pure logical function is_contact_active(state)
    integer, intent(in) :: state
    is_contact_active = (state >= CONTACTSTICK)
  end function

  !> Whether the contact state is completely free (no projection info)
  pure logical function is_contact_free(state)
    integer, intent(in) :: state
    is_contact_free = (state == CONTACTFREE)
  end function

  !> Which of the three categories the contact state belongs to
  pure integer function contact_state_category(state)
    integer, intent(in) :: state
    if( is_contact_free(state) ) then
      contact_state_category = kcatFREE
    else if( is_contact_active(state) ) then
      contact_state_category = kcatCONT
    else
      contact_state_category = kcatNEAR
    endif
  end function

  !> Finalizer
  subroutine fstr_contact_finalize( contact )
    type(tContact), intent(inout)     :: contact !< contact definition
    integer  :: i
    if( associated( contact%slave ) ) deallocate(contact%slave)
    if( associated( contact%master ) ) then
      do i=1,size( contact%master )
        call  finalize_surf( contact%master(i) )
      enddo
      deallocate(contact%master)
    endif
    if( associated(contact%states) ) deallocate(contact%states)
    call hecmw_contact_comm_finalize(contact%comm)
    call bucketDB_finalize( contact%master_bktDB )
  end subroutine

  !>  Check the consistency with given mesh of contact definition
  logical function fstr_contact_check( contact, hecMESH )
    type(tContact), intent(inout)     :: contact  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition

    integer  :: i
    logical  :: isfind

    fstr_contact_check = .false.

    ! if contact pair exist?
    isfind = .false.
    do i=1,hecMESH%contact_pair%n_pair
      if( hecMESH%contact_pair%name(i) == contact%pair_name ) then
        contact%ctype = hecMESH%contact_pair%type(i)
        contact%surf_id1 = hecMESH%contact_pair%slave_grp_id(i)
        contact%surf_id2 = hecMESH%contact_pair%master_grp_id(i)
        contact%surf_id1_sgrp = hecMESH%contact_pair%slave_orisgrp_id(i)
        isfind = .true.
      endif
    enddo
    if( .not. isfind ) return;
    if( contact%fcoeff<=0.d0 ) contact%fcoeff=0.d0
    if( contact%ctype < 1 .and. contact%ctype > 3 ) return
    if( contact%group<=0 ) return

    fstr_contact_check = .true.
  end function

  !>  Number of slave nodes of this contact owned by the current rank
  integer(kind=kint) function fstr_count_internal_slaves( contact, hecMESH )
    type(tContact), intent(in)        :: contact  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition

    integer  :: i, is, ie, cgrp

    fstr_count_internal_slaves = 0

    cgrp = contact%surf_id1
    if( cgrp<=0 ) return
    is= hecMESH%node_group%grp_index(cgrp-1) + 1
    ie= hecMESH%node_group%grp_index(cgrp  )
    do i=is,ie
      if( hecMESH%node_group%grp_item(i) <= hecMESH%nn_internal ) then
        fstr_count_internal_slaves = fstr_count_internal_slaves + 1
      endif
    enddo
  end function

  !>  Initializer of tContactState
  logical function fstr_contact_init( contact, hecMESH, cparam )
    type(tContact), intent(inout)     :: contact  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    type(tContactParam), target       :: cparam   !< contact parameter

    integer  :: i, j, is, ie, cgrp, nsurf, nslave, ic, ic_type, iss, nn, ii
    integer  :: count, ID_area
    logical  :: slave_owner, take_master
    integer, allocatable  :: slave_index(:)  !< global node -> slave index (SURF-SURF only)

    fstr_contact_init = .false.

    contact%cparam => cparam

    slave_owner = hecmw_partcontact_get_owner( hecMESH%hecmw_flag_partcontact ) == HECMW_FLAG_PARTCONTACT_OWNER_SLAVE

    ! update_surface_normal normalises vertex normals only after the cross-rank assembly, so a uniform factor cancels
    ! but a partial one does not: a rank owning no slave must take no master surface at all, even though the
    ! partitioner does leave master elements here to keep the communication tables symmetric
    take_master = .true.
    if( slave_owner ) take_master = fstr_count_internal_slaves( contact, hecMESH ) > 0

    !  master surface
    cgrp = contact%surf_id2
    if( cgrp<=0 ) return
    is= hecMESH%surf_group%grp_index(cgrp-1) + 1
    ie= hecMESH%surf_group%grp_index(cgrp  )

    ! Owned-only master face count, independent of take_master: the SURF-SURF visibility
    ! guard in fstr_setup sums it over ranks to get the global unique master face count.
    contact%n_master_owned = 0
    do i=is,ie
      ic   = hecMESH%surf_group%grp_item(2*i-1)
      if( hecMESH%elem_ID(ic*2) == hecMESH%my_rank ) contact%n_master_owned = contact%n_master_owned + 1
    enddo

    count = 0
    if( take_master ) then
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        ID_area = hecMESH%elem_ID(ic*2)
        if( .not. slave_owner .and. ID_area /= hecMESH%my_rank ) cycle
        count = count + 1
      enddo
    endif
    allocate( contact%master(count) )
    count = 0
    if( take_master ) then
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        ID_area = hecMESH%elem_ID(ic*2)
        if( .not. slave_owner .and. ID_area /= hecMESH%my_rank ) cycle
        count = count + 1
        nsurf = hecMESH%surf_group%grp_item(2*i)
        ic_type = hecMESH%elem_type(ic)
        call initialize_surf( ic, ic_type, nsurf, contact%master(count) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( contact%master(count)%nodes )
          nn = contact%master(count)%nodes(j)
          contact%master(count)%nodes(j) = hecMESH%elem_node_item( iss+nn )
        enddo
      enddo
    endif

    call update_surface_reflen( contact%master, hecMESH%node )

    cgrp = contact%surf_id1
    if( cgrp<=0 ) return
    is= hecMESH%node_group%grp_index(cgrp-1) + 1
    ie= hecMESH%node_group%grp_index(cgrp  )
    nslave = 0
    do i=is,ie
      if( slave_owner .and. hecMESH%node_group%grp_item(i) > hecMESH%nn_internal ) cycle
      nslave = nslave + 1
    enddo
    allocate( contact%slave(nslave) )
    allocate( slave_index(hecMESH%n_node) )
    slave_index(:) = 0
    ii = 0
    do i=is,ie
      if( slave_owner .and. hecMESH%node_group%grp_item(i) > hecMESH%nn_internal ) cycle
      ii = ii + 1
      contact%slave(ii) = hecMESH%node_group%grp_item(i)
      slave_index(contact%slave(ii)) = ii
    enddo

    ! contact state
    allocate( contact%states(nslave) )
    do i=1,nslave
      call contact_state_init( contact%states(i) )
    enddo

    ! neighborhood of surface group
    call update_surface_box_info( contact%master, hecMESH%node )
    call bucketDB_init( contact%master_bktDB )
    call update_surface_bucket_info( contact%master, contact%master_bktDB )
    call find_surface_neighbor( contact%master, contact%master_bktDB )

    if(contact%method == CONTACTS2S) then
      ! The mortar integral passes the master shape functions through a length-4 array and
      ! sizes its element vectors for a first-order pair, so a second-order master face
      ! overruns them. initialize_csurf rejects the slave side for the same reason.
      do i=1, size( contact%master )
        if( contact%master(i)%etype /= fe_quad4n .and. contact%master(i)%etype /= fe_tri3n ) then
          write(*,*) '### Error: MORTAR=YES supports first-order surfaces only (quad4/tri3) : etype=', &
            contact%master(i)%etype
          stop HECMW_EXIT_MODEL
        endif
      enddo

      !  slave surface
      cgrp = contact%surf_id1_sgrp
      if( cgrp<=0 ) return
      is= hecMESH%surf_group%grp_index(cgrp-1) + 1
      ie= hecMESH%surf_group%grp_index(cgrp  )

      ! Slave segments are taken owned-only (a serial mesh owns every element, so this is
      ! the full slave surface there). The master surface must be visible in full on every
      ! slave-owning rank; that is what !PARTITION, CONTACT_OWNER=SLAVE gives and what the
      ! visibility guard in fstr_setup checks.
      count = 0
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        if( hecMESH%elem_ID(ic*2) /= hecMESH%my_rank ) cycle
        count = count + 1
      enddo
      allocate( contact%slave_surf(count) )
      count = 0
      do i=is,ie
        ic   = hecMESH%surf_group%grp_item(2*i-1)
        if( hecMESH%elem_ID(ic*2) /= hecMESH%my_rank ) cycle
        count = count + 1
        nsurf = hecMESH%surf_group%grp_item(2*i)
        ic_type = hecMESH%elem_type(ic)
        call initialize_csurf( ic, ic_type, nsurf, contact%slave_surf(count) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( contact%slave_surf(count)%nodes )
          nn = contact%slave_surf(count)%nodes(j)
          contact%slave_surf(count)%nodes(j) = hecMESH%elem_node_item( iss+nn )
          contact%slave_surf(count)%nslave_index(j) = slave_index(hecMESH%elem_node_item( iss+nn ))
        enddo
      enddo

      ! state for each integration points
      do i=1, size( contact%slave_surf )
        nn = contact%slave_surf(i)%n_intp
        allocate( contact%slave_surf(i)%states(nn) )
        do j = 1, contact%slave_surf(i)%n_intp
          contact%slave_surf(i)%states(j)%state = -1
          contact%slave_surf(i)%states(j)%multiplier(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce1(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce_trial(:) = 0.d0
          contact%slave_surf(i)%states(j)%tangentForce_final(:) = 0.d0
          contact%slave_surf(i)%states(j)%reldisp(:) = 0.d0
          contact%slave_surf(i)%states(j)%time_factor = 0.d0
          contact%slave_surf(i)%states(j)%interference_flag = 0
        enddo
      enddo
    endif

    ! initialize contact communication table
    call hecmw_contact_comm_init( contact%comm, hecMESH, 1, nslave, contact%slave )

    contact%symmetric = .true.
    fstr_contact_init = .true.
  end function

  !>  Initializer of tContactState for embed case
  logical function fstr_embed_init( embed, hecMESH, cparam )
    type(tContact), intent(inout)     :: embed  !< contact definition
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    type(tContactParam), target       :: cparam   !< contact parameter

    integer  :: i, j, is, ie, cgrp, nslave, ic, ic_type, iss, nn, ii
    integer  :: count, ID_area
    logical  :: slave_owner, take_master

    fstr_embed_init = .false.

    embed%cparam => cparam

    slave_owner = hecmw_partcontact_get_owner( hecMESH%hecmw_flag_partcontact ) == HECMW_FLAG_PARTCONTACT_OWNER_SLAVE

    take_master = .true.
    if( slave_owner ) take_master = fstr_count_internal_slaves( embed, hecMESH ) > 0

    !  master surface
    cgrp = embed%surf_id2
    if( cgrp<=0 ) return
    is= hecMESH%elem_group%grp_index(cgrp-1) + 1
    ie= hecMESH%elem_group%grp_index(cgrp  )

    count = 0
    if( take_master ) then
      do i=is,ie
        ic   = hecMESH%elem_group%grp_item(i)
        ID_area = hecMESH%elem_ID(ic*2)
        if( .not. slave_owner .and. ID_area /= hecMESH%my_rank ) cycle
        count = count + 1
      enddo
    endif
    allocate( embed%master(count) )
    count = 0
    if( take_master ) then
      do i=is,ie
        ic   = hecMESH%elem_group%grp_item(i)
        ID_area = hecMESH%elem_ID(ic*2)
        if( .not. slave_owner .and. ID_area /= hecMESH%my_rank ) cycle
        count = count + 1
        ic_type = hecMESH%elem_type(ic)
        call initialize_surf( ic, ic_type, 0, embed%master(count) )
        iss = hecMESH%elem_node_index(ic-1)
        do j=1, size( embed%master(count)%nodes )
          nn = embed%master(count)%nodes(j)
          embed%master(count)%nodes(j) = hecMESH%elem_node_item( iss+nn )
        enddo
      enddo
    endif

    ! slave surface
    cgrp = embed%surf_id1
    if( cgrp<=0 ) return
    is= hecMESH%node_group%grp_index(cgrp-1) + 1
    ie= hecMESH%node_group%grp_index(cgrp  )
    nslave = 0
    do i=is,ie
      if( slave_owner .and. hecMESH%node_group%grp_item(i) > hecMESH%nn_internal ) cycle
      nslave = nslave + 1
    enddo
    allocate( embed%slave(nslave) )
    ii = 0
    do i=is,ie
      if( slave_owner .and. hecMESH%node_group%grp_item(i) > hecMESH%nn_internal ) cycle
      ii = ii + 1
      embed%slave(ii) = hecMESH%node_group%grp_item(i)
    enddo

    ! embed state
    allocate( embed%states(nslave) )
    do i=1,nslave
      call contact_state_init( embed%states(i) )
    enddo

    ! neighborhood of surface group
    call update_surface_box_info( embed%master, hecMESH%node )
    call bucketDB_init( embed%master_bktDB )
    call update_surface_bucket_info( embed%master, embed%master_bktDB )
    call find_surface_neighbor( embed%master, embed%master_bktDB )

    ! initialize contact communication table
    call hecmw_contact_comm_init( embed%comm, hecMESH, 1, nslave, embed%slave )

    ! initialize penalty coefficients
    embed%nPenalty = 1.0d0     ! default normal penalty coefficient
    embed%tPenalty = 0.1d0     ! default tangential penalty coefficient
    embed%refStiff = 0.0d0     ! will be calculated after first stiffness assembly
    embed%damp_alpha = 0.0d0
    embed%damp_gact = 0.0d0

    embed%symmetric = .true.
    fstr_embed_init = .true.
  end function

  function check_apply_Contact_IF( contact_if, contacts )
    type(tContactInterference), intent(inout)     :: contact_if  !< contact definition
    type(tContact)     :: contacts(:) !< type fstr_solid
    
    integer  :: i, j
    logical  :: isfind
    integer(kind=kint)            :: check_apply_Contact_IF

    check_apply_Contact_IF = -1
    ! if contact pair exist?
    isfind = .false.
    do i = 1, size(contacts)
      if( contacts(i)%pair_name == contact_if%cp_name ) then
        ! !CONTACT_INTERFERENCE is not implemented for the SURF-SURF mortar path
        if( contacts(i)%method == CONTACTS2S ) then
          write(*,*) '### Error: CONTACT_INTERFERENCE is not supported with MORTAR=YES'
          stop HECMW_EXIT_MODEL
        endif
        contacts(i)%if_type     = contact_if%if_type
        contacts(i)%if_etime    = contact_if%etime
        contacts(i)%initial_pos = contact_if%initial_pos
        contacts(i)%end_pos     = contact_if%end_pos
        do j = 1, size(contacts(i)%states)
          contacts(i)%states(j)%interference_flag = contact_if%if_type
          contacts(i)%states(j)%init_pos = contact_if%initial_pos
          contacts(i)%states(j)%end_pos  = contact_if%end_pos
          if( contact_if%if_type /= C_IF_SLAVE )then
            contacts(i)%states(j)%time_factor = (contact_if%end_pos - contact_if%initial_pos) / contact_if%etime
          else
            contacts(i)%states(j)%time_factor = contact_if%etime
          end if
        end do
        isfind = .true.
        check_apply_Contact_IF = 0; return
      endif
    enddo
    if( .not. isfind ) return;
    check_apply_Contact_IF = 0

  end function

  !> Reset contact state all to free
  subroutine clear_contact_state( contact )
    type(tContact), intent(inout) :: contact    !< contact definition
    integer :: i
    if( .not. associated(contact%states) ) return
    do i=1,size( contact%states )
      contact%states(i)%state = -1
    enddo
  end subroutine

  !> if contact is active is curr step
  logical function is_active_contact( acgrp, contact )
    integer, intent(in)        :: acgrp(:)      !< active contact group numbers
    type(tContact), intent(in) :: contact       !< contact definition
    if( any( acgrp==contact%group ) ) then
      is_active_contact = .true.
    else
      is_active_contact = .false.
    endif
  end function

  !> Initializer of a slave surface segment (SURF-SURF mortar)
  subroutine initialize_csurf( eid, etype, nsurf, surf )
    use elementInfo
    integer(kind=kint), intent(in)    :: eid    !< element ID
    integer(kind=kint), intent(in)    :: etype  !< element type
    integer(kind=kint), intent(in)    :: nsurf  !< surface ID
    type(tContactSurf), intent(inout) :: surf   !< surface element
    integer(kind=kint) :: n, outtype, nodes(100)
    surf%eid = eid

    call getSubFace( etype, nsurf, outtype, nodes )
    surf%etype = outtype
    n=getNumberOfNodes( outtype )
    if(surf%etype == fe_quad4n )then
      surf%n_intp = 16
    else if(surf%etype == fe_tri3n )then
      surf%n_intp = 27
    else
      write(*,*) '### Error: MORTAR=YES supports first-order surfaces only (quad4/tri3) : etype=', surf%etype
      stop HECMW_EXIT_MODEL
    end if
    allocate( surf%nodes(n) )
    allocate( surf%nslave_index(n) )
    surf%nodes(1:n)=nodes(1:n)
    surf%nslave_index(:)= 0
    ! lam_*_val: dim1 = slave-surf node a (size n), dim2 = rank r (master, size n_intp)
    allocate( surf%lam_begin_id(surf%n_intp), surf%lam_begin_val(n,surf%n_intp) )
    allocate( surf%lam_work_id (surf%n_intp), surf%lam_work_val (n,surf%n_intp) )
    surf%lam_begin_id(:) = 0; surf%lam_begin_val(:,:) = 0.0d0; surf%lam_begin_n = 0
    surf%lam_work_id (:) = 0; surf%lam_work_val (:,:) = 0.0d0; surf%lam_work_n  = 0
    ! state_begin / state_prev_begin keep their type default (CONTACTFREE)
    ! tangential friction parallel arrays: (2, node a, rank r) / (node a, rank r)
    allocate( surf%lam_begin_t(2,n,surf%n_intp), surf%lam_work_t(2,n,surf%n_intp) )
    allocate( surf%lam_begin_fstate(n,surf%n_intp), surf%lam_work_fstate(n,surf%n_intp) )
    surf%lam_begin_t(:,:,:) = 0.0d0; surf%lam_work_t(:,:,:) = 0.0d0
    surf%lam_begin_fstate(:,:) = CONTACTSTICK; surf%lam_work_fstate(:,:) = CONTACTSTICK
  end subroutine


end module mContactDef
