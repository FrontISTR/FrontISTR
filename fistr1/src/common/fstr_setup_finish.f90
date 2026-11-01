!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to process the data read from the cards,
!! independently of the input file format: the allocation and the preparation after the cards are counted
!! (fstr_setup_prepare) and the processing after all the cards are read (fstr_setup_finish).
module m_fstr_setup_finish
  use m_fstr
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_setup_util_err_stop &
    , append_node_grp_from_surf_grp &
    , reallocate_integer &
    , reallocate_real
  use m_fstr_solid_setup, only &
    : fstr_solid_alloc &
    , fstr_element_init
  use mContactParam
  use m_static_get_prop
  implicit none

  !> Counts of the input for fstr_setup_prepare, independent of the input file format.
  !> A reader of an input file format sets them from what it counts in the file.
  type fstr_setup_counts
    integer(kind=kint) :: contact_pairs = 0
    integer(kind=kint) :: embed_pairs = 0
    integer(kind=kint) :: contact_interferences = 0
    integer(kind=kint) :: contact_params = 0
    integer(kind=kint) :: weld_lines = 0
    integer(kind=kint) :: initial_conditions = 0
    integer(kind=kint) :: steps = 0
    integer(kind=kint) :: local_coords = 0
    integer(kind=kint) :: autoinc_params = 0
    integer(kind=kint) :: converg_params = 0
    integer(kind=kint) :: time_points = 0
    integer(kind=kint) :: materials = 0       !< materials defined in the input; the mesh can define more
  end type fstr_setup_counts

contains

  !> Allocation and preparation after the cards are counted: the arrays sized by counts,
  !> the materials and sections from the mesh, and the output controls.
  !> The section defaults depend on the analysis type, and fstr_convert_rotation_dof rewrites the boundary conditions,
  !> so this is called after the analysis type and all the boundary conditions are read.
  subroutine fstr_setup_prepare( P, counts )
    use mMaterial
    implicit none
    type(fstr_param_pack) :: P
    type(fstr_setup_counts), intent(in) :: counts

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    type(fstr_heat), pointer :: fstrHEAT
    type(fstr_dynamic), pointer :: fstrDYNAMIC
    integer(kind=kint) :: n, i, cid, n_totlyr
    real(kind=kreal) :: ee, pp, rho, alpha, thick, alpha_over_mu
    real(kind=kreal) :: beam_radius,                          &
      beam_angle1, beam_angle2, beam_angle3,&
      beam_angle4, beam_angle5, beam_angle6
    type(tshellmat),pointer :: shmat(:)
    character(len=HECMW_FILENAME_LEN) :: logfileNAME

    hecMESH => P%MESH
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID
    fstrHEAT => P%HEAT
    fstrDYNAMIC => P%DYN
    write( logfileNAME, '(i5,''.log'')' ) myrank

    fstrSOLID%n_contacts = counts%contact_pairs
    if( counts%contact_pairs>0 ) then
      allocate( fstrSOLID%contacts( counts%contact_pairs ) )
      ! convert SURF_SURF contact to NODE_SURF contact
      call fstr_convert_contact_type( P%MESH )
    endif
    call fstr_convert_rotation_dof( P )
    fstrSOLID%n_embeds = counts%embed_pairs
    if( counts%embed_pairs>0 )  allocate( fstrSOLID%embeds( counts%embed_pairs ) )
    if( counts%weld_lines>0 ) allocate( fstrHEAT%weldline( counts%weld_lines ) )
    if( counts%initial_conditions>0 ) allocate( g_InitialCnd( counts%initial_conditions ) )
    if( counts%steps>0 ) then
      allocate( fstrSOLID%step_ctrl( counts%steps ) )
      do i=1, counts%steps
        call init_stepInfo( fstrSOLID%step_ctrl(i) )
        if( p%PARAM%solution_type==kstDYNAMIC ) then
          fstrSOLID%step_ctrl(i)%num_substep = fstrDYNAMIC%n_step
          fstrSOLID%step_ctrl(i)%initdt      = fstrDYNAMIC%t_delta
          fstrSOLID%step_ctrl(i)%elapsetime  = dble(fstrDYNAMIC%n_step) * fstrDYNAMIC%t_delta
          fstrSOLID%step_ctrl(i)%mindt       = fstrDYNAMIC%t_delta
          fstrSOLID%step_ctrl(i)%maxdt       = fstrDYNAMIC%t_delta
        endif
      end do
    endif
    if( counts%local_coords>0 ) allocate( g_LocalCoordSys(counts%local_coords) )
    allocate( fstrPARAM%ainc(0:counts%autoinc_params) )
    do i=0,counts%autoinc_params
      call init_AincParam( fstrPARAM%ainc(i) )
    end do
    allocate( fstrPARAM%cnvparam(0:counts%converg_params) )
    do i=0,counts%converg_params
      call init_ConvergParam( fstrPARAM%cnvparam(i) )
    end do
    if( counts%time_points>0 ) allocate( fstrPARAM%timepoints(counts%time_points) )
    allocate( fstrPARAM%contactparam(0:counts%contact_params) )
    do i=0,counts%contact_params
      call init_ContactParam( fstrPARAM%contactparam(i) )
    end do
    if( counts%contact_interferences>0 )then
      allocate( fstrPARAM%contact_if( counts%contact_interferences ) )
      do i=1,counts%contact_interferences
        call init_Contact_IF( fstrPARAM%contact_if(i) )
      end do
    end if

    P%SOLID%is_33shell = 0
    P%SOLID%is_33beam  = 0

    do i=1,hecMESH%n_elem_type
      n =  hecMESH%elem_type_item(i)
      if (n == 781 .or. n == 761)then
        P%SOLID%is_33shell = 1
      elseif (n == 641)then
        P%SOLID%is_33beam  = 1
      endif
    enddo

    n = counts%materials
    if( hecMESH%material%n_mat>n ) n= hecMESH%material%n_mat
    if( n==0 ) stop "material property not defined!"
    allocate( fstrSOLID%materials( n ) )
    do i = 1, n
      call initMaterial(fstrSOLID%materials(i))
    enddo

    ! for section control
    allocate( fstrSOLID%sections(hecMESH%section%n_sect) )
    do i=1,hecMESH%section%n_sect
      ! set default 361 element formulation
      if( p%PARAM%solution_type==kstSTATIC .or. p%PARAM%solution_type==kstDYNAMIC ) then
        if( p%PARAM%nlgeom ) then
          fstrSOLID%sections(i)%elemopt361 = kel361FBAR
        else
          fstrSOLID%sections(i)%elemopt361 = kel361IC
        end if
      else if( p%PARAM%solution_type==kstEIGEN ) then
        fstrSOLID%sections(i)%elemopt361 = kel361IC
      else if( p%PARAM%solution_type==kstSTATICEIGEN ) then
        fstrSOLID%sections(i)%elemopt361 = kel361FBAR
      else
        fstrSOLID%sections(i)%elemopt361 = kel361FI
      end if
      fstrSOLID%sections(i)%elemopt341 = kel341FI
      fstrSOLID%sections(i)%elemopt611 = kel611EULER
      ! INTERFACE sections are skipped in the loop below and keep this thickness
      fstrSOLID%sections(i)%thickness = 1.0d0
    enddo

    if( hecMESH%section%n_sect >0 ) then
      do i=1,hecMESH%section%n_sect
        if( hecMESH%section%sect_type(i) == 4 ) cycle
        cid = hecMESH%section%sect_mat_ID_item(i)
        if( cid>n ) stop "Error in material property definition!"
        if( fstrPARAM%nlgeom .or. fstrPARAM%solution_type==kstSTATICEIGEN ) &
          fstrSOLID%materials(cid)%nlgeom_flag = 1
        nullify(shmat)
        call fstr_get_prop(hecMESH,shmat,i,ee,pp,rho,alpha,thick,&
          n_totlyr,alpha_over_mu, &
          beam_radius,beam_angle1,beam_angle2,beam_angle3, &
          beam_angle4,beam_angle5,beam_angle6)
        fstrSOLID%sections(i)%thickness = thick
        fstrSOLID%materials(cid)%name = hecMESH%material%mat_name(cid)
        fstrSOLID%materials(cid)%variables(M_YOUNGS)=ee
        fstrSOLID%materials(cid)%variables(M_POISSON)=pp
        fstrSOLID%materials(cid)%variables(M_DENSITY)=rho
        fstrSOLID%materials(cid)%variables(M_EXPANSION)=alpha
        fstrSOLID%materials(cid)%variables(M_ALPHA_OVER_MU)= alpha_over_mu
        fstrSOLID%materials(cid)%variables(M_BEAM_RADIUS)=beam_radius
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE1)=beam_angle1
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE2)=beam_angle2
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE3)=beam_angle3
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE4)=beam_angle4
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE5)=beam_angle5
        fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE6)=beam_angle6
        fstrSOLID%materials(cid)%mtype = ELASTIC
        if( hecMESH%section%sect_type(i) == 2 ) then
          fstrSOLID%materials(cid)%totallyr = n_totlyr
          fstrSOLID%materials(cid)%shell_var => shmat
        endif
      enddo
    endif

    allocate( fstrSOLID%output_ctrl( 4 ) )
    call fstr_init_outctrl(fstrSOLID%output_ctrl(1))
    fstrSOLID%output_ctrl( 1 )%filename = trim(logfileNAME)
    fstrSOLID%output_ctrl( 1 )%filenum = ILOG
    call fstr_init_outctrl(fstrSOLID%output_ctrl(2))
    call fstr_init_outctrl(fstrSOLID%output_ctrl(3))
    call fstr_init_outctrl(fstrSOLID%output_ctrl(4))
  end subroutine fstr_setup_prepare

  !> Processing after all the cards are read: the temperature arrays, the default step,
  !> the element initialization and the checks, and the result arrays (fstr_setup_post).
  subroutine fstr_setup_finish( P, c_material )
    implicit none
    type(fstr_param_pack) :: P
    integer(kind=kint), intent(in) :: c_material

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_param), pointer :: fstrPARAM
    type(fstr_solid), pointer :: fstrSOLID
    type(fstr_dynamic), pointer :: fstrDYNAMIC
    integer(kind=kint) :: n, i, cid, ierror, n_totlyr

    hecMESH => P%MESH
    fstrPARAM => P%PARAM
    fstrSOLID => P%SOLID
    fstrDYNAMIC => P%DYN

    ! ----- material type judgement. in case of infinitive analysis, nlgeom_flag=0
    if( .not. P%PARAM%nlgeom ) then
      do i=1, c_material
        fstrSOLID%materials(i)%nlgeom_flag = 0
      enddo
    endif

    if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
      allocate ( fstrSOLID%temperature( hecMESH%n_node )      ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TEMPERATURE>'
        write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrSOLID%temperature = fstrPARAM%ref_temp
      allocate ( fstrSOLID%last_temp( hecMESH%n_node )      ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, LAST_TEMP>'
        write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrSOLID%last_temp = 0.d0
      allocate ( fstrSOLID%temp_bak( hecMESH%n_node )      ,stat=ierror )
      if( ierror /= 0 ) then
        write(idbg,*) 'stop due to allocation error <FSTR_SOLID, TEMP_BAK>'
        write(idbg,*) '  rank = ', myrank,'  ierror = ',ierror
        call flush(idbg)
        call hecmw_abort( hecmw_comm_get_comm())
      end if
      fstrSOLID%temp_bak = 0.d0
    endif

    if( associated(fstrSOLID%step_ctrl) )  then
      fstrSOLID%nstep_tot = size(fstrSOLID%step_ctrl)
      call setup_stepInfo_starttime( fstrSOLID%step_ctrl )
    else
      if( p%PARAM%solution_type==kstSTATIC .and. P%PARAM%nlgeom ) then
        write( *,* ) " ERROR: STEP not defined!"
        write( idbg,* ) "ERROR: STEP not defined!"
        call flush(idbg)
        call fstr_abort( HECMW_EXIT_MODEL )
      endif

      if( myrank==0 ) write(*,*)"Step control not defined! Using default step=1"
      fstrSOLID%nstep_tot = 1
      allocate( fstrSOLID%step_ctrl(1) )
      call init_stepInfo( fstrSOLID%step_ctrl(1) )
      if( p%PARAM%solution_type==kstDYNAMIC ) then
        fstrSOLID%step_ctrl(1)%num_substep = fstrDYNAMIC%n_step
        fstrSOLID%step_ctrl(1)%initdt      = fstrDYNAMIC%t_delta
        fstrSOLID%step_ctrl(1)%elapsetime  = dble(fstrDYNAMIC%n_step) * fstrDYNAMIC%t_delta
        fstrSOLID%step_ctrl(1)%mindt       = fstrDYNAMIC%t_delta
        fstrSOLID%step_ctrl(1)%maxdt       = fstrDYNAMIC%t_delta
      endif
      n =  fstrSOLID%BOUNDARY_ngrp_tot
      if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%Boundary(n) )
      do i = 1, n
        fstrSOLID%step_ctrl(1)%Boundary(i) = fstrSOLID%BOUNDARY_ngrp_GRPID(i)
      enddo
      n = fstrSOLID%CLOAD_ngrp_tot + fstrSOLID%DLOAD_ngrp_tot + fstrSOLID%TEMP_ngrp_tot + fstrSOLID%SPRING_ngrp_tot
      if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%Load(n) )
      n = 0
      do i = 1, fstrSOLID%CLOAD_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%CLOAD_ngrp_GRPID(i)
      enddo
      do i = 1, fstrSOLID%DLOAD_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%DLOAD_ngrp_GRPID(i)
      enddo
      do i = 1, fstrSOLID%TEMP_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%TEMP_ngrp_GRPID(i)
      enddo
      do i = 1, fstrSOLID%SPRING_ngrp_tot
        n = n + 1
        fstrSOLID%step_ctrl(1)%Load(n) = fstrSOLID%SPRING_ngrp_GRPID(i)
      enddo
      n = fstrSOLID%elemact%ELEMACT_egrp_tot
      if( n>0 ) allocate( fstrSOLID%step_ctrl(1)%ElemActivation(n) )
      do i = 1, n
        fstrSOLID%step_ctrl(1)%ElemActivation(i) = fstrSOLID%elemact%ELEMACT_egrp_GRPID(i)
      enddo
    endif

    call fstr_element_init( hecMESH, fstrSOLID, p%PARAM%solution_type )
    call fstr_check_shell_j2_inputs(hecMESH, fstrSOLID, fstrPARAM)
    call fstr_check_up_inputs(hecMESH, fstrSOLID)
    call fstr_check_shell_gauss_output(hecMESH, fstrSOLID)
    if( p%PARAM%solution_type==kstSTATIC .or. p%PARAM%solution_type==kstDYNAMIC .or.   &
      p%PARAM%solution_type==kstEIGEN  .or. p%PARAM%solution_type==kstSTATICEIGEN )  &
      call fstr_solid_alloc( hecMESH, fstrSOLID )

    if( p%PARAM%solution_type == kstHEAT) then
      p%PARAM%fg_irres = fstrSOLID%output_ctrl(3)%frequency
      p%PARAM%fg_iwres = fstrSOLID%output_ctrl(4)%frequency
      p%HEAT%elemact = p%SOLID%elemact
    endif

    n_totlyr = 1
    do i=1,hecMESH%section%n_sect
      cid = hecMESH%section%sect_mat_ID_item(i)
      n  = fstrSOLID%materials(cid)%totallyr
      if (n > n_totlyr)then
        n_totlyr = n
      endif
    enddo
    P%SOLID%max_lyr = n_totlyr

    call fstr_setup_post( P )
  end subroutine fstr_setup_finish

  subroutine fstr_check_up_inputs(hecMESH, fstrSOLID)
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint) :: i, isect, mtype

    do i = 1, hecMESH%n_elem
      if( hecMESH%elem_type(i) /= 361 ) cycle
      isect = hecMESH%section_ID(i)
      if( fstrSOLID%sections(isect)%elemopt361 /= kel361UP ) cycle
      if( .not.associated(fstrSOLID%elements(i)%gausses) ) cycle
      mtype = fstrSOLID%elements(i)%gausses(1)%pMaterial%mtype
      if( .not.isElastoplastic(mtype) .and. .not.isCreep(mtype) ) cycle
      write(*,*) '###ERROR### : Elastoplastic and creep materials are not supported by FORM361=UP'
      write(ILOG,*) '###ERROR### : Elastoplastic and creep materials are not supported by FORM361=UP'
      call hecmw_abort(hecmw_comm_get_comm())
    enddo
  end subroutine fstr_check_up_inputs

  subroutine fstr_check_shell_gauss_output(hecMESH, fstrSOLID)
    use mMechGauss, only: fstr_shell_output_point_count
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    integer(kind=kint) :: i

    if( hecMESH%n_dof /= 6 ) return
    if( .not.any(fstrSOLID%output_ctrl(3)%outinfo%on(9:10)) ) return
    do i=1,hecMESH%n_elem
      if( fstr_shell_output_point_count(fstrSOLID%elements(i))>0 ) cycle
      write(*,*) '###ERROR### : ISTRAIN/ISTRESS require shell thickness-point output support; element type = ', &
        hecMESH%elem_type(i)
      call hecmw_abort(hecmw_comm_get_comm())
    enddo
  end subroutine fstr_check_shell_gauss_output

  subroutine fstr_check_shell_j2_inputs(hecMESH, fstrSOLID, fstrPARAM)
    use m_ElastoPlastic, only: checkPlaneStressJ2Material, PlaneStressJ2ErrorMessage
    use elementInfo, only: NumOfShellThicknessQuadPoints
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(hecmwST_local_mesh), intent(in) :: hecMESH
    type(fstr_solid), intent(in) :: fstrSOLID
    type(fstr_param), intent(in) :: fstrPARAM
    type(tMaterial), pointer :: material
    integer(kind=kint) :: i, nn, ierr, isect, index, ilayer
    real(kind=kreal) :: elastic(2)
    logical :: missing
    character(len=160) :: message

    do i = 1, hecMESH%n_elem
      if( hecMESH%elem_type(i) /= 741 ) cycle
      if( .not.associated(fstrSOLID%elements(i)%gausses) ) cycle
      material => fstrSOLID%elements(i)%gausses(1)%pMaterial
      if( .not.isElastoplastic(material%mtype) ) cycle
      nn = hecMESH%elem_node_index(i)-hecMESH%elem_node_index(i-1)
      message = ''
      call checkPlaneStressJ2Material(material, ierr)
      if( ierr /= 0 ) then
        message = PlaneStressJ2ErrorMessage(ierr)
      else if( .not.fstr_uses_plane_stress_j2_shell(741, nn, material) ) then
        message = 'Elastoplastic MITC4 shells require INFINITESIMAL material kinematics'
      else if( fstrPARAM%solution_type /= kstSTATIC .or. .not.fstrPARAM%nlgeom ) then
        message = 'Elastoplastic MITC4 shells require !SOLUTION, TYPE=STATIC, NONLINEAR'
      else if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
        message = 'Temperature loading is not supported in models with elastoplastic MITC4 shells'
      endif
      if( len_trim(message) == 0 .and. associated(material%shell_var) ) then
        ! The J2 return mapping uses the common elastic table, not layer-specific constants.
        call fetch_TableData(MC_ISOELASTIC, material%dict, elastic, missing)
        do ilayer = 1, material%totallyr
          if( .not.ieee_is_finite(material%shell_var(ilayer)%ee) .or. &
              .not.ieee_is_finite(material%shell_var(ilayer)%pp) ) then
            message = 'Invalid elastic properties in J2 shell layer'
            exit
          endif
          if( abs(material%shell_var(ilayer)%ee-elastic(1)) > 1.0d-12*elastic(1) .or. &
              abs(material%shell_var(ilayer)%pp-elastic(2)) > 1.0d-12 ) then
            message = 'J2 shell layers must share the !ELASTIC Young modulus and Poisson ratio'
            exit
          endif
        enddo
      endif
      if( len_trim(message) == 0 ) then
        isect = hecMESH%section_ID(i)
        index = hecMESH%section%sect_I_index(isect-1)
        if( hecMESH%section%sect_I_index(isect)-index /= 1 ) then
          message = 'Elastoplastic MITC4 shells require a SHELL section with 2 thickness integration points'
        else if( hecMESH%section%sect_I_item(index+1) /= NumOfShellThicknessQuadPoints(741) ) then
          message = 'Elastoplastic MITC4 shells require 2 thickness integration points in !SECTION'
        endif
      endif
      if( len_trim(message) == 0 ) cycle
      write(*,*) '###ERROR### : ', trim(message)
      write(ILOG,*) '###ERROR### : ', trim(message)
      call hecmw_abort(hecmw_comm_get_comm())
    enddo
  end subroutine fstr_check_shell_j2_inputs

  !-----------------------------------------------------------------------------!
  !> Initial setting of postprecessor

  subroutine fstr_setup_post_phys_alloc(phys, NDOF, n_node, n_elem)
    implicit none
    type(fstr_solid_physic_val), pointer :: phys
    integer(kind=kint) :: NDOF, n_node, n_elem, mdof, istat
    mdof = (NDOF*NDOF+NDOF)/2;
    allocate ( phys%STRAIN  (mdof*n_node), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%STRAIN"
    allocate ( phys%STRESS  (mdof*n_node), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%STRESS"
    allocate ( phys%MISES   (     n_node), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%MISES"
    allocate ( phys%ESTRAIN (mdof*n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%ESTRAIN"
    allocate ( phys%ESTRESS (mdof*n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%ESTRESS"
    allocate ( phys%EMISES  (     n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%EMISES"
    allocate ( phys%EPLSTRAIN (   n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%EPLSTRAIN"
    allocate ( phys%ENQM    (12*n_elem), stat=istat)
    if( istat /= 0 ) stop "Allocation error: phys%ENQM"
  end subroutine fstr_setup_post_phys_alloc

  subroutine fstr_setup_post( P )
    implicit none
    integer(kind=kint) :: i
    type(fstr_param_pack) :: P
    type(fstr_solid_physic_val), pointer :: phys => null()

    if( P%PARAM%solution_type == kstSTATIC &
        .or. P%PARAM%solution_type == kstEIGEN  &
        .or. P%PARAM%solution_type == kstDYNAMIC &
        .or. P%PARAM%solution_type == kstSTATICEIGEN ) then
      ! Memory Allocation for Result Vectors ------------
      if( P%MESH%n_dof == 6 .or. P%SOLID%is_33shell == 1 ) then
        allocate ( P%SOLID%SHELL )
        call fstr_setup_post_phys_alloc(P%SOLID%SHELL,3, P%MESH%n_node,P%MESH%n_elem)
        allocate ( P%SOLID%SHELL%LAYER(P%SOLID%max_lyr) )
        do i=1,P%SOLID%max_lyr
          allocate ( P%SOLID%SHELL%LAYER(i)%PLUS )
          allocate ( P%SOLID%SHELL%LAYER(i)%MINUS )
          call fstr_setup_post_phys_alloc(P%SOLID%SHELL%LAYER(i)%PLUS , 3, P%MESH%n_node, P%MESH%n_elem)
          call fstr_setup_post_phys_alloc(P%SOLID%SHELL%LAYER(i)%MINUS, 3, P%MESH%n_node, P%MESH%n_elem)
        enddo
        phys => P%SOLID%SHELL
      else
        allocate ( P%SOLID%SOLID )
        phys => P%SOLID%SOLID
        call fstr_setup_post_phys_alloc(phys, P%MESH%n_dof, P%MESH%n_node,  P%MESH%n_elem)
      end if
      P%SOLID%STRAIN => phys%STRAIN
      P%SOLID%STRESS => phys%STRESS
      P%SOLID%MISES  => phys%MISES
      P%SOLID%ESTRAIN => phys%ESTRAIN
      P%SOLID%ESTRESS => phys%ESTRESS
      P%SOLID%EMISES  => phys%EMISES
      P%SOLID%EPLSTRAIN  => phys%EPLSTRAIN
      P%SOLID%ENQM    => phys%ENQM
      allocate( P%SOLID%REACTION( P%MESH%n_dof*P%MESH%n_node ), stat=i )
      if( i /= 0 ) stop "Allocation error: REACTION"
      P%SOLID%REACTION = 0.d0
    end if

    call hecmw_barrier( P%MESH ) ! JP-7

    if( P%HEAT%STEPtot == 0 ) then ! No !HEAT Input
      if( P%PARAM%analysis_n == 0 ) then  ! No !STATIC Input
        call reallocate_real( P%PARAM%dtime, 1)
        call reallocate_real( P%PARAM%etime, 1)
        call reallocate_real( P%PARAM%dtmin, 1)
        call reallocate_real( P%PARAM%delmax,1)
        call reallocate_integer( P%PARAM%itmax, 1)
        call reallocate_real( P%PARAM%eps,   1)
        P%PARAM%analysis_n = 1
        P%PARAM%dtime = 0
        P%PARAM%etime = 0
        P%PARAM%dtmin = 0
        P%PARAM%delmax = 0
        P%PARAM%itmax = 20
        P%PARAM%eps = 1.0e-6
      end if
      P%HEAT%STEPtot = 1
      call reallocate_real( P%HEAT%STEP_DLTIME, 1)
      call reallocate_real( P%HEAT%STEP_EETIME, 1)
      call reallocate_real( P%HEAT%STEP_DELMIN, 1)
      call reallocate_real( P%HEAT%STEP_DELMAX, 1)
      P%HEAT%STEP_DLTIME = 0
      P%HEAT%STEP_EETIME = 0
      P%HEAT%STEP_DELMIN = 0
      P%HEAT%STEP_DELMAX = 0
    end if
  end subroutine fstr_setup_post

  !-----------------------------------------------------------------------------!
  !> Convert SURF-SURF contact to NODE-SURF contact                             !
  !-----------------------------------------------------------------------------!

  subroutine fstr_convert_contact_type( hecMESH )
    implicit none
    type(hecmwST_local_mesh), pointer :: hecMESH  !< mesh definition
    integer(kind=kint) :: n, i, sgrp_id, ngrp_id
    ! convert SURF_SURF to NODE_SURF
    n = hecMESH%contact_pair%n_pair
    do i = 1,n
      if( hecMESH%contact_pair%type(i) /= HECMW_CONTACT_TYPE_SURF_SURF ) cycle
      sgrp_id = hecMESH%contact_pair%slave_grp_id(i)
      call append_node_grp_from_surf_grp( hecMESH, sgrp_id, ngrp_id )
      ! change type of contact and slave group ID
      hecMESH%contact_pair%type(i) = HECMW_CONTACT_TYPE_NODE_SURF
      hecMESH%contact_pair%slave_grp_id(i) = ngrp_id
    enddo
  end subroutine fstr_convert_contact_type

  !-----------------------------------------------------------------------------!
  !> Move dof 4-6 of 641/761/781 nodes to their rotational dummy nodes          !
  !-----------------------------------------------------------------------------!

  !> A translational node of 641/761/781 has dof 1-3 only, and its rotation is dof 1-3
  !> of the paired dummy node. Dof 4-6 of !BOUNDARY and !CLOAD given to the translational
  !> nodes are moved to the dummy nodes generated from 611/731/741; the other conditions
  !> cannot give dof 4-6 to them. With explicit dummy nodes in the mesh, the rotation is
  !> given only as dof 1-3 of the dummy nodes.
  subroutine fstr_convert_rotation_dof( P )
    use hecmw_setup_util, only : append_new_groups, get_grp_id
    implicit none
    type(fstr_param_pack) :: P
    type(hecmwST_local_mesh), pointer :: hecMESH
    integer(kind=kint), allocatable :: rot_node(:), n_rot(:), n_rot_all(:), n_explicit(:), rot_grp(:), src(:)
    integer(kind=kint), allocatable :: new_index(:), new_list(:)
    character(len=HECMW_NAME_LEN), allocatable :: new_name(:)
    logical, allocatable :: is_rot(:), in_all(:), is_explicit(:)
    integer(kind=kint) :: icel, is, ie, nn, i, j, k, m, n, ig, n_new, first_id, ids, ide

    hecMESH => P%MESH
    if( hecMESH%n_dof /= 3 ) return

    ! explicit dummy nodes are in ALL, and the generated ones are not
    allocate( in_all(hecMESH%n_node) )
    in_all = .false.
    ig = get_grp_id( hecMESH, 'node_grp', 'ALL' )
    if( ig > 0 ) then
      is = hecMESH%node_group%grp_index(ig-1)
      in_all(hecMESH%node_group%grp_item(is+1:hecMESH%node_group%grp_index(ig))) = .true.
    endif

    allocate( rot_node(hecMESH%n_node), is_explicit(hecMESH%n_node) )
    rot_node = 0
    is_explicit = .false.
    do icel = 1, hecMESH%n_elem
      if( .not. hecmw_is_etype_33struct(hecMESH%elem_type(icel)) ) cycle
      nn = hecmw_get_max_node(hecMESH%elem_type(icel))/2
      is = hecMESH%elem_node_index(icel-1)
      do j = 1, nn
        rot_node(hecMESH%elem_node_item(is+j)) = hecMESH%elem_node_item(is+nn+j)
        if( in_all(hecMESH%elem_node_item(is+nn+j)) ) is_explicit(hecMESH%elem_node_item(is+j)) = .true.
      enddo
    enddo
    ! a rank may hold an external translational node without its 641/761/781 element,
    ! so the conditions are converted on every rank by the counts over all the ranks
    n = count( rot_node > 0 )
    call hecmw_allreduce_I1( hecMESH, n, hecmw_max )
    if( n == 0 ) return

    allocate( n_rot(hecMESH%node_group%n_grp), n_rot_all(hecMESH%node_group%n_grp) )
    allocate( n_explicit(hecMESH%node_group%n_grp) )
    do ig = 1, hecMESH%node_group%n_grp
      is = hecMESH%node_group%grp_index(ig-1)
      n_rot(ig) = count( rot_node(hecMESH%node_group%grp_item(is+1:hecMESH%node_group%grp_index(ig))) > 0 )
      n_explicit(ig) = count( is_explicit(hecMESH%node_group%grp_item(is+1:hecMESH%node_group%grp_index(ig))) )
    enddo
    n_rot_all = n_rot
    call hecmw_allreduce_I( hecMESH, n_rot_all, hecMESH%node_group%n_grp, hecmw_max )

    do i = 1, P%SOLID%SPRING_ngrp_tot
      if( P%SOLID%SPRING_ngrp_DOF(i) > 3 .and. n_rot_all(P%SOLID%SPRING_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!SPRING' )
    enddo
    do i = 1, P%SOLID%VELOCITY_ngrp_tot
      if( P%SOLID%VELOCITY_ngrp_rotID(i) > 0 ) cycle
      if( mod(P%SOLID%VELOCITY_ngrp_type(i), 10) > 3 .and. n_rot_all(P%SOLID%VELOCITY_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!VELOCITY' )
    enddo
    do i = 1, P%SOLID%ACCELERATION_ngrp_tot
      if( mod(P%SOLID%ACCELERATION_ngrp_type(i), 10) > 3 .and. n_rot_all(P%SOLID%ACCELERATION_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!ACCELERATION' )
    enddo
    do i = 1, P%FREQ%FLOAD_ngrp_tot
      if( P%FREQ%FLOAD_ngrp_TYPE(i) /= kFLOADTYPE_NODE ) cycle
      if( P%FREQ%FLOAD_ngrp_DOF(i) > 3 .and. n_rot_all(P%FREQ%FLOAD_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!FLOAD' )
    enddo
    !the rotational DOFs of !FACCELERATION are also not allowed (the end DOF is mod(type,10))
    do i = 1, P%FREQ%FACC_ngrp_tot
      if( mod(P%FREQ%FACC_ngrp_type(i), 10) > 3 .and. n_rot_all(P%FREQ%FACC_ngrp_ID(i)) > 0 ) &
        call rotation_dof_err_stop( '!FACCELERATION' )
    enddo
    do i = 1, hecMESH%mpc%n_mpc
      do j = hecMESH%mpc%mpc_index(i-1)+1, hecMESH%mpc%mpc_index(i)
        if( hecMESH%mpc%mpc_dof(j) > 3 .and. rot_node(hecMESH%mpc%mpc_item(j)) > 0 ) &
          call rotation_dof_err_stop( '!EQUATION' )
      enddo
    enddo

    ! node groups whose rotation is moved: rot_grp(ig) = -1, then the id of the dummy node group
    allocate( rot_grp(hecMESH%node_group%n_grp) )
    rot_grp = 0
    do i = 1, P%SOLID%BOUNDARY_ngrp_tot
      if( P%SOLID%BOUNDARY_ngrp_rotID(i) > 0 ) cycle
      ig = P%SOLID%BOUNDARY_ngrp_ID(i)
      if( mod(P%SOLID%BOUNDARY_ngrp_type(i), 10) <= 3 .or. n_rot_all(ig) == 0 ) cycle
      if( n_explicit(ig) > 0 ) call explicit_dummy_node_err_stop( '!BOUNDARY' )
      rot_grp(ig) = -1
    enddo
    do i = 1, P%SOLID%CLOAD_ngrp_tot
      if( P%SOLID%CLOAD_ngrp_rotID(i) > 0 .or. P%SOLID%CLOAD_ngrp_DOF(i) <= 3 ) cycle
      ig = P%SOLID%CLOAD_ngrp_ID(i)
      is = hecMESH%node_group%grp_index(ig-1)
      ie = hecMESH%node_group%grp_index(ig)
      if( any( rot_node(hecMESH%node_group%grp_item(is+1:ie)) == 0 .and. &
        hecMESH%node_group%grp_item(is+1:ie) <= hecMESH%nn_internal ) ) call fstr_setup_util_err_stop( &
        '### Error: !CLOAD : dof 4-6 is given to nodes without rotational dof' )
      if( n_explicit(ig) > 0 ) call explicit_dummy_node_err_stop( '!CLOAD' )
      rot_grp(ig) = -1
    enddo

    n_new = count( rot_grp == -1 )
    if( n_new == 0 ) return
    allocate( new_name(n_new), new_index(0:n_new), new_list(sum(n_rot, mask=rot_grp == -1)) )
    new_index(0) = 0
    k = 0
    m = 0
    do ig = 1, hecMESH%node_group%n_grp
      if( rot_grp(ig) /= -1 ) cycle
      k = k + 1
      write( new_name(k), '(a,i0)' ) 'FSTR_ROT_', ig
      do j = hecMESH%node_group%grp_index(ig-1)+1, hecMESH%node_group%grp_index(ig)
        if( rot_node(hecMESH%node_group%grp_item(j)) == 0 ) cycle
        m = m + 1
        new_list(m) = rot_node(hecMESH%node_group%grp_item(j))
      enddo
      new_index(k) = m
      rot_grp(ig) = k
    enddo
    call append_new_groups( hecMESH, 'node_grp', n_new, new_name, new_index, new_list, first_id )
    where( rot_grp > 0 ) rot_grp = rot_grp + first_id - 1

    ! split each !BOUNDARY entry into dof 1-3 of its nodes and dof 1-3 of their dummy nodes
    n = P%SOLID%BOUNDARY_ngrp_tot
    if( n > 0 ) then
      allocate( src(2*n), is_rot(2*n) )
      m = 0
      do i = 1, n
        ide = mod(P%SOLID%BOUNDARY_ngrp_type(i), 10)
        if( P%SOLID%BOUNDARY_ngrp_rotID(i) <= 0 .and. ide > 3 .and. rot_grp(P%SOLID%BOUNDARY_ngrp_ID(i)) > 0 ) then
          if( P%SOLID%BOUNDARY_ngrp_type(i)/10 <= 3 ) then
            m = m + 1
            src(m) = i
            is_rot(m) = .false.
          endif
          m = m + 1
          src(m) = i
          is_rot(m) = .true.
        else
          m = m + 1
          src(m) = i
          is_rot(m) = .false.
        endif
      enddo
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_GRPID, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_ID, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_type, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_amp, src, m )
      call remap_real_array( P%SOLID%BOUNDARY_ngrp_val, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_istot, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_rotID, src, m )
      call remap_integer_array( P%SOLID%BOUNDARY_ngrp_centerID, src, m )
      P%SOLID%BOUNDARY_ngrp_tot = m
      do k = 1, m
        if( P%SOLID%BOUNDARY_ngrp_rotID(k) > 0 ) cycle
        ig = P%SOLID%BOUNDARY_ngrp_ID(k)
        ids = P%SOLID%BOUNDARY_ngrp_type(k)/10
        ide = mod(P%SOLID%BOUNDARY_ngrp_type(k), 10)
        if( ide <= 3 .or. rot_grp(ig) <= 0 ) cycle
        if( is_rot(k) ) then
          P%SOLID%BOUNDARY_ngrp_ID(k) = rot_grp(ig)
          P%SOLID%BOUNDARY_ngrp_type(k) = 10 * (max(ids, 4) - 3) + (ide - 3)
        else
          P%SOLID%BOUNDARY_ngrp_type(k) = 10 * ids + 3
        endif
      enddo
    endif

    do i = 1, P%SOLID%CLOAD_ngrp_tot
      if( P%SOLID%CLOAD_ngrp_rotID(i) > 0 .or. P%SOLID%CLOAD_ngrp_DOF(i) <= 3 ) cycle
      P%SOLID%CLOAD_ngrp_ID(i) = rot_grp(P%SOLID%CLOAD_ngrp_ID(i))
      P%SOLID%CLOAD_ngrp_DOF(i) = P%SOLID%CLOAD_ngrp_DOF(i) - 3
    enddo
  end subroutine fstr_convert_rotation_dof

  subroutine rotation_dof_err_stop( header_name )
    implicit none
    character(len=*) :: header_name

    call fstr_setup_util_err_stop( '### Error: '//header_name// &
      ' : dof 4-6 of beam/shell nodes in a mesh with solid elements is not supported' )
  end subroutine rotation_dof_err_stop

  subroutine explicit_dummy_node_err_stop( header_name )
    implicit none
    character(len=*) :: header_name

    call fstr_setup_util_err_stop( '### Error: '//header_name// &
      ' : dof 4-6 is given to beam/shell nodes with explicit dummy nodes; give it as dof 1-3 of the dummy nodes' )
  end subroutine explicit_dummy_node_err_stop

  subroutine remap_integer_array( array, src, n )
    implicit none
    integer(kind=kint), pointer :: array(:)
    integer(kind=kint), intent(in) :: src(:), n
    integer(kind=kint), pointer :: new_array(:)

    allocate( new_array(n) )
    new_array(1:n) = array(src(1:n))
    deallocate( array )
    array => new_array
  end subroutine remap_integer_array

  subroutine remap_real_array( array, src, n )
    implicit none
    real(kind=kreal), pointer :: array(:)
    integer(kind=kint), intent(in) :: src(:), n
    real(kind=kreal), pointer :: new_array(:)

    allocate( new_array(n) )
    new_array(1:n) = array(src(1:n))
    deallocate( array )
    array => new_array
  end subroutine remap_real_array

end module m_fstr_setup_finish
