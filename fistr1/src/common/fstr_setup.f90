!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to read in data from control file
!! and do necessary preparation for following calculation
module m_fstr_setup
  use m_fstr
  use hecmw_setup_util, only &
    : node_grp_name_to_id_ex &
    , elem_grp_name_to_id_ex &
    , surf_grp_name_to_id_ex &
    , hecmw_streqr
  use fstr_setup_util
  use m_fstr_solid_setup
  use m_fstr_heat_setup
  use m_fstr_freq_setup
  use m_fstr_eigen_setup
  use m_fstr_dynamic_setup
  use m_fstr_param_setup
  use m_fstr_setup_finish
  use fstr_ctrl_common
  use fstr_ctrl_static
  use fstr_ctrl_heat
  use fstr_ctrl_eigen
  use fstr_ctrl_dynamic
  use fstr_ctrl_material
  use mContactParam
  use m_static_get_prop
  use m_out
  use m_step
  use m_utilities
  use fstr_ctrl_util_f
  implicit none

contains

  !=============================================================================!
  !> Read in and initialize control data                                        !
  !=============================================================================!
  subroutine fstr_setup( cntl_filename, hecMESH, hecMAT, fstrPARAM,  &
      fstrSOLID, fstrEIG, fstrHEAT, fstrDYNAMIC, fstrCPL, fstrFREQ )
    character(len=HECMW_FILENAME_LEN) :: cntl_filename
    type(hecmwST_local_mesh),target :: hecMESH
    type(hecmwST_matrix),target :: hecMAT
    type(fstr_param),target   :: fstrPARAM
    type(fstr_solid),target   :: fstrSOLID
    type(fstr_eigen),target     :: fstrEIG
    type(fstr_heat),target    :: fstrHEAT
    type(fstr_dynamic),target :: fstrDYNAMIC
    type(fstr_couple),target  :: fstrCPL
    type(fstr_freqanalysis), target :: fstrFREQ

    integer(kind=kint) :: ctrl, rcode
    type(fstr_param_pack) :: P
    type(fstr_card_table) :: T
    integer(kind=kint) :: total(FSTR_CARD_MAX)
    type(fstr_setup_counts) :: counts

    ! packaging
    P%MESH   => hecMESH
    P%MAT    => hecMAT
    P%PARAM  => fstrPARAM
    P%SOLID  => fstrSOLID
    P%EIGEN  => fstrEIG
    P%HEAT   => fstrHEAT
    P%DYN    => fstrDYNAMIC
    P%CPL    => fstrCPL
    P%FREQ   => fstrFREQ

    fstrPARAM%contact_algo = kcaALagrange

    call fstr_setup_cards( T )

    ctrl  = fstr_ctrl_open( cntl_filename )
    if( ctrl < 0 ) then
      write(*,*) '### Error: Cannot open FSTR control file : ', cntl_filename
      write(ILOG,*) '### Error: Cannot open FSTR control file : ', cntl_filename
      stop HECMW_EXIT_INPUT
    end if

    call fstr_setup_read_cards( ctrl, T, 1, total, P )

    ! -----
    counts = fstr_setup_card_counts( T, total )
    call fstr_setup_prepare( P, counts )

    ! -----
    rcode = fstr_ctrl_rewind( ctrl )
    call fstr_setup_read_cards( ctrl, T, 2, total, P )

    call fstr_setup_finish( P, counts )
    if( P%PARAM%fg_visual == kON ) call fstr_setup_visualize( ctrl, P%MESH )
    rcode = fstr_ctrl_close( ctrl )

  end subroutine fstr_setup

  !> Adds the cards of the control file to the table T, by the names without the '!' of the header.
  !> The first reading reads the settings and the cards added to arrays as they are read, and counts the cards for
  !> fstr_setup_prepare. The second reading reads the cards stored in the arrays that fstr_setup_prepare allocates,
  !> and the cards that refer to them.
  subroutine fstr_setup_cards( T )
    implicit none
    type(fstr_card_table) :: T

    ! Settings, read in the first reading. fstr_setup_prepare and the second reading use them wherever they are
    ! in the file, and the cards of the first reading use those read before them.
    ! fstr_setup_prepare uses !SOLUTION, !STATIC and !DYNAMIC, !CONTACT uses !CONTACT_ALGO,
    ! and in the first reading, !STEP uses !VERSION and !DLOAD uses NLGEOM of !SOLUTION, !STATIC and !DYNAMIC.
    call fstr_card_add( T, 'VERSION',          first=fstr_setup_VERSION )
    call fstr_card_add( T, 'SOLUTION',         first=fstr_setup_SOLUTION )
    call fstr_card_add( T, 'STATIC',           first=fstr_setup_STATIC )
    call fstr_card_add( T, 'DYNAMIC',          first=fstr_setup_DYNAMIC )
    call fstr_card_add( T, 'CONTACT_ALGO',     first=fstr_setup_CONTACTALGO )
    call fstr_card_add( T, 'NONLINEAR_SOLVER', first=fstr_setup_NONLINEAR_SOLVER )
    call fstr_card_add( T, 'SOLVER',           first=fstr_setup_SOLVER )
    call fstr_card_add( T, 'EIGEN',            first=fstr_setup_EIGEN )
    call fstr_card_add( T, 'ECHO',             first=fstr_setup_ECHO )
    call fstr_card_add( T, 'RESTART',          first=fstr_setup_RESTART )
    call fstr_card_add( T, 'OUTPUT_SSTYPE',    first=fstr_setup_OUTPUT_SSTYPE )

    ! Cards added to arrays as they are read. fstr_setup_prepare converts or checks the rotational dofs of
    ! !BOUNDARY, !CLOAD, !SPRING, !VELOCITY, !ACCELERATION, !FLOAD and !FACCELERATION,
    ! and the cards of the first reading refer to !AMPLITUDE by name.
    call fstr_card_add( T, 'AMPLITUDE',          first=fstr_setup_AMPLITUDE )
    call fstr_card_add( T, 'BOUNDARY',           first=fstr_setup_BOUNDARY )
    call fstr_card_add( T, 'CLOAD',              first=fstr_setup_CLOAD )
    call fstr_card_add( T, 'SPRING',             first=fstr_setup_SPRING )
    call fstr_card_add( T, 'VELOCITY',           first=fstr_setup_VELOCITY )
    call fstr_card_add( T, 'ACCELERATION',       first=fstr_setup_ACCELERATION )
    call fstr_card_add( T, 'FLOAD',              first=fstr_setup_FLOAD )
    call fstr_card_add( T, 'FACCELERATION',      first=fstr_setup_FACCELERATION )
    call fstr_card_add( T, 'DLOAD',              first=fstr_setup_DLOAD )
    call fstr_card_add( T, 'TEMPERATURE',        first=fstr_setup_TEMPERATURE )
    call fstr_card_add( T, 'REFTEMP',            first=fstr_setup_REFTEMP )
    call fstr_card_add( T, 'ELEMENT_ACTIVATION', first=fstr_setup_ELEMENT_ACTIVATION )
    call fstr_card_add( T, 'FIXTEMP',            first=fstr_setup_FIXTEMP )
    call fstr_card_add( T, 'CFLUX',              first=fstr_setup_CFLUX )
    call fstr_card_add( T, 'DFLUX',              first=fstr_setup_DFLUX )
    call fstr_card_add( T, 'SFLUX',              first=fstr_setup_SFLUX )
    call fstr_card_add( T, 'FILM',               first=fstr_setup_FILM )
    call fstr_card_add( T, 'SFILM',              first=fstr_setup_SFILM )
    call fstr_card_add( T, 'RADIATE',            first=fstr_setup_RADIATE )
    call fstr_card_add( T, 'SRADIATE',           first=fstr_setup_SRADIATE )
    call fstr_card_add( T, 'EIGENREAD',          first=fstr_setup_eigenread )
    call fstr_card_add( T, 'COUPLE',             first=fstr_setup_COUPLE )
    call fstr_card_add( T, 'MPC',                first=fstr_setup_MPC )

    ! Counted in the first reading to allocate the arrays they are stored in, and stored in the second reading.
    ! !CONTACT also uses !CONTACT_ALGO, and !CONTACT_INTERFERENCE refers to the pairs of !CONTACT.
    call fstr_card_add( T, 'ORIENTATION',          second=fstr_setup_ORIENTATION )
    call fstr_card_add( T, 'CONTACT',              second=fstr_setup_CONTACT, by_line=.true. )
    call fstr_card_add( T, 'EMBED',                second=fstr_setup_EMBED, by_line=.true. )
    call fstr_card_add( T, 'CONTACT_INTERFERENCE', second=fstr_setup_CONTACT_INTERFERENCE, by_line=.true. )
    call fstr_card_add( T, 'CONTACT_PARAM',        second=fstr_setup_CONTACT_PARAM )
    call fstr_card_add( T, 'AUTOINC_PARAM',        second=fstr_setup_AUTOINC_PARAM )
    call fstr_card_add( T, 'CONVERG_PARAM',        second=fstr_setup_CONVERG_PARAM )
    call fstr_card_add( T, 'TIME_POINTS',          second=fstr_setup_TIME_POINTS )
    call fstr_card_add( T, 'INITIAL_CONDITION',    second=fstr_setup_INITIAL_CONDITION )
    call fstr_card_add( T, 'WELD_LINE',            second=fstr_setup_WELD_LINE )
    call fstr_card_add( T, 'MATERIAL',             second=fstr_setup_MATERIAL )
    ! !ISTEP also refers by name to !AUTOINC_PARAM, !CONVERG_PARAM and !TIME_POINTS read before it.
    ! !STEP of !VERSION 0 is read in the first reading. !STEP of !VERSION 1 or later is a step like !ISTEP,
    ! stored in the same array. !STEP of a negative !VERSION is not read, but is counted as a step in the first reading.
    call fstr_card_add( T, 'ISTEP',                second=fstr_setup_ISTEP, seq='steps' )
    call fstr_card_add( T, 'STEP',                 first=fstr_setup_STEP, version_min=0, version_max=0 )
    call fstr_card_add( T, 'STEP',                 second=fstr_setup_STEP_store, seq='steps', version_min=1 )
    call fstr_card_add( T, 'STEP',                 seq='steps', version_max=-1 )

    ! Read in the second reading to refer by name to the cards of the second reading read before them:
    ! !SECTION to !ORIENTATION, and !HEAT to !TIME_POINTS
    call fstr_card_add( T, 'SECTION',              second=fstr_setup_SECTION )
    call fstr_card_add( T, 'HEAT',                 second=fstr_setup_HEAT )

    ! Stored in the output controls that fstr_setup_prepare allocates.
    ! !WRITE also sets the output flags in the first reading.
    call fstr_card_add( T, 'WRITE',                first=fstr_setup_WRITE_flags, second=fstr_setup_WRITE_store )
    call fstr_card_add( T, 'OUTPUT_RES',           second=fstr_setup_OUTPUT_RES )
    call fstr_card_add( T, 'OUTPUT_VIS',           second=fstr_setup_OUTPUT_VIS )

    ! Stored in the material of the last !MATERIAL
    call fstr_card_add( T, 'ELASTIC',              second=fstr_setup_ELASTIC )
    call fstr_card_add( T, 'PLASTIC',              second=fstr_setup_PLASTIC )
    call fstr_card_add( T, 'HYPERELASTIC',         second=fstr_setup_HYPERELASTIC )
    call fstr_card_add( T, 'VISCOELASTIC',         second=fstr_setup_VISCOELASTIC )
    call fstr_card_add( T, 'TRS',                  second=fstr_setup_TRS )
    call fstr_card_add( T, 'CREEP',                second=fstr_setup_CREEP )
    call fstr_card_add( T, 'DENSITY',              second=fstr_setup_DENSITY )
    call fstr_card_add( T, 'EXPANSION_COEF',       second=fstr_setup_EXPANSION_COEFF )
    call fstr_card_alias( T, 'EXPANSION_COEFF', 'EXPANSION_COEF' )
    call fstr_card_alias( T, 'EXPANSION',       'EXPANSION_COEF' )
    call fstr_card_add( T, 'DAMPING',              second=fstr_setup_DAMPING )
    call fstr_card_add( T, 'FLUID',                second=fstr_setup_FLUID )
    call fstr_card_add( T, 'SPRING_D',             second=fstr_setup_SPRING_D )
    call fstr_card_add( T, 'SPRING_A',             second=fstr_setup_SPRING_A )
    call fstr_card_add( T, 'DASHPOT_D',            second=fstr_setup_DASHPOT_D )
    call fstr_card_add( T, 'DASHPOT_A',            second=fstr_setup_DASHPOT_A )
    call fstr_card_add( T, 'USER_MATERIAL',        second=fstr_setup_USER_MATERIAL )

    ! Read in the second reading
    call fstr_card_add( T, 'ELEMOPT',              second=fstr_setup_ELEMOPT )
    call fstr_card_add( T, 'ULOAD',                second=fstr_setup_ULOAD )
  end subroutine fstr_setup_cards

  !> Counts for fstr_setup_prepare, from the totals of the sequences of the table T in the first reading
  function fstr_setup_card_counts( T, total ) result( counts )
    implicit none
    type(fstr_card_table), intent(in) :: T
    integer(kind=kint), intent(in) :: total(:)
    type(fstr_setup_counts) :: counts

    counts%contact_pairs         = fstr_card_total( T, total, 'CONTACT' )
    counts%embed_pairs           = fstr_card_total( T, total, 'EMBED' )
    counts%contact_interferences = fstr_card_total( T, total, 'CONTACT_INTERFERENCE' )
    counts%contact_params        = fstr_card_total( T, total, 'CONTACT_PARAM' )
    counts%weld_lines            = fstr_card_total( T, total, 'WELD_LINE' )
    counts%initial_conditions    = fstr_card_total( T, total, 'INITIAL_CONDITION' )
    counts%steps                 = fstr_card_total( T, total, 'steps' )
    counts%local_coords          = fstr_card_total( T, total, 'ORIENTATION' )
    counts%autoinc_params        = fstr_card_total( T, total, 'AUTOINC_PARAM' )
    counts%converg_params        = fstr_card_total( T, total, 'CONVERG_PARAM' )
    counts%time_points           = fstr_card_total( T, total, 'TIME_POINTS' )
    counts%materials             = fstr_card_total( T, total, 'MATERIAL' )
  end function fstr_setup_card_counts

  !> Reads the cards of the control file in order, opening the files of !INCLUDE and stopping at !END,
  !> and calls the reader of each card in the table T for the reading ipass (1 or 2), setting the places of the card
  !> in its sequence to P%STATE%slot. The cards are found in T by the names without the '!' of the header.
  !> Each reading counts the sequences from 0. The first reading counts all the cards in T and returns the totals of
  !> the sequences in total. The second reading counts only the cards it reads,
  !> and stops if a sequence goes beyond its total, which the arrays are allocated by. The cards not in T are skipped.
  !> On return, ctrl is the file where the reading ended.
  subroutine fstr_setup_read_cards( ctrl, T, ipass, total, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_card_table), intent(in) :: T
    integer(kind=kint), intent(in) :: ipass
    integer(kind=kint), intent(inout) :: total(:)
    type(fstr_param_pack) :: P

    integer(kind=kint) :: ctrl_list(20), ictrl, rcode, ierror, k, s
    integer(kind=kint) :: pos(size(total))
    character(len=HECMW_NAME_LEN) :: header_name
    character(len=HECMW_FILENAME_LEN) :: input_filename

    pos(:) = 0
    ctrl_list = 0
    ictrl = 1
    do
      rcode = fstr_ctrl_get_c_h_name( ctrl, header_name, HECMW_NAME_LEN )
      if( header_name == '!INCLUDE' ) then
        ctrl_list(ictrl) = ctrl
        input_filename   = ""
        ierror = fstr_ctrl_get_param_ex( ctrl, 'INPUT ', '# ', 0, 'F', input_filename )
        ctrl   = fstr_ctrl_open( input_filename )
        if( ctrl < 0 ) then
          write(*,*) '### Error: Cannot open FSTR control file : ', input_filename
          write(ILOG,*) '### Error: Cannot open FSTR control file : ', input_filename
          stop HECMW_EXIT_INPUT
        end if
        ictrl = ictrl + 1
        cycle
      else if( header_name == '!END' ) then
        exit
      end if

      k = fstr_card_find( T, header_name(2:), P%STATE%version )
      if( k > 0 ) then
        if( ipass == 1 .or. associated(T%card(k)%second) ) then
          s = T%card(k)%seq
          P%STATE%slot%first = pos(s) + 1
          if( T%card(k)%by_line ) then
            pos(s) = pos(s) + fstr_ctrl_get_data_line_n( ctrl )
          else
            pos(s) = pos(s) + 1
          endif
          P%STATE%slot%last = pos(s)
          if( ipass == 2 .and. pos(s) > total(s) ) call fstr_setup_util_err_stop( '### Error: '// &
            trim(header_name)//' : more than counted in the first reading of the control file' )
          if( ipass == 1 .and. associated(T%card(k)%first) ) call T%card(k)%first( ctrl, P )
          if( ipass == 2 ) call T%card(k)%second( ctrl, P )
        endif
      endif

      ! next
      if( fstr_ctrl_seek_next_header(ctrl) == 0 )then
        if( ictrl == 1 )then
          exit
        else
          ierror= fstr_ctrl_close( ctrl )
          ictrl = ictrl - 1
          ctrl  = ctrl_list(ictrl)
          if( fstr_ctrl_seek_next_header(ctrl) == 0 ) exit
        endif
      endif
    end do
    if( ipass == 1 ) total(:) = pos(:)
  end subroutine fstr_setup_read_cards

  subroutine fstr_couple_finalize( C )
    implicit none
    type( fstr_couple ) :: C
    if(associated(C%coupled_node)) deallocate(C%coupled_node)
    if(associated(C%trac)) deallocate(C%trac)
    if(associated(C%disp)) deallocate(C%disp)
    if(associated(C%velo)) deallocate(C%velo)
    if(associated(C%accel)) deallocate(C%accel)
    if(associated(C%index)) deallocate(C%index)
  end subroutine fstr_couple_finalize

end module m_fstr_setup
