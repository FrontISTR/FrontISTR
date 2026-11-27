!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides subroutines to initialize fstr_heat and to read
!! the cards of heat analysis into it. Boundary conditions are appended
!! independently of the input file format.
!! header_name is the card name shown in error messages. A blank name of
!! amplitude means that it is not given.
module m_fstr_heat_setup
  use m_fstr
  use hecmw_setup_util, only &
    : hecmw_streqr &
    , get_grp_member &
    , get_grp_member_n &
    , hecmw_expand_integer_array &
    , hecmw_expand_real_array &
    , hecmw_expand_integer_array2 &
    , hecmw_expand_real_array2
  use fstr_setup_util, only &
    : fstr_param_pack &
    , fstr_card_slot &
    , fstr_ctrl_err_stop &
    , reallocate_integer &
    , reallocate_real &
    , amp_name_to_id &
    , get_local_member_index &
    , get_sorted_local_member_index
  use fstr_ctrl_heat
  use fstr_ctrl_util_f
  implicit none

contains

  !> Initial setting of heat analysis
  subroutine fstr_heat_init( fstrHEAT )
    implicit none
    type(fstr_heat) :: fstrHEAT

    fstrHEAT%STEPtot     = 0
    fstrHEAT%MATERIALtot = 0
    fstrHEAT%AMPLITUDEtot= 0
    fstrHEAT%T_FIX_tot   = 0
    fstrHEAT%Q_NOD_tot   = 0
    fstrHEAT%Q_VOL_tot   = 0
    fstrHEAT%Q_SUF_tot   = 0
    fstrHEAT%R_SUF_tot   = 0
    fstrHEAT%H_SUF_tot   = 0
    fstrHEAT%WL_tot      = 0
    fstrHEAT%beta        = -1.0d0
  end subroutine fstr_heat_init

  !C***
  !C*** INIT_AMPLITUDE
  !C***
  subroutine heat_init_amplitude (hecMESH, fstrHEAT)

    use m_fstr

    implicit none
    integer(kind=kint) :: namax, i, nn, is, iE, icou, j, k
    real(kind=kreal)   :: x1, y1, x2, y2
    type(fstr_heat)    :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH
    !C
    !C===
    namax = 0
    do i = 1, hecMESH%amp%n_amp
      nn = hecMESH%amp%amp_index(i) - hecMESH%amp%amp_index(i-1)
      namax = max(nn,namax)
    enddo

    fstrHEAT%AMPLITUDEtot= hecMESH%amp%n_amp

    if(associated(fstrHEAT%AMPLtab)) deallocate(fstrHEAT%AMPLtab)
    if(associated(fstrHEAT%AMPL)) deallocate(fstrHEAT%AMPL)
    if(associated(fstrHEAT%AMPLtime)) deallocate(fstrHEAT%AMPLtime)

    allocate (fstrHEAT%AMPLtab (fstrHEAT%AMPLITUDEtot) )
    allocate (fstrHEAT%AMPL    (fstrHEAT%AMPLITUDEtot,namax) )
    allocate (fstrHEAT%AMPLtime(fstrHEAT%AMPLITUDEtot,namax) )

    fstrHEAT%AMPLtab  = 0
    fstrHEAT%AMPL     = 0.d0
    fstrHEAT%AMPLtime = 0.d0

    do i = 1, fstrHEAT%AMPLITUDEtot
      is = hecMESH%amp%amp_index(i-1) + 1
      iE = hecMESH%amp%amp_index(i)
      nn = iE - is + 1

      fstrHEAT%AMPLtab(i) = nn

      icou = 0
      do j = is, iE
        icou = icou + 1
        fstrHEAT%AMPL    (i,icou) = hecMESH%amp%amp_val  (j)
        fstrHEAT%AMPLtime(i,icou) = hecMESH%amp%amp_table(j)
      enddo
    enddo
    !C===

    !C
    !C +-----------+
    !C | AMP-TABLE |
    !C +-----------+
    !C===
    if(associated(fstrHEAT%AMPLfuncA)) deallocate(fstrHEAT%AMPLfuncA)
    if(associated(fstrHEAT%AMPLfuncB)) deallocate(fstrHEAT%AMPLfuncB)

    allocate ( fstrHEAT%AMPLfuncA( fstrHEAT%AMPLITUDEtot,namax+1 ) )
    allocate ( fstrHEAT%AMPLfuncB( fstrHEAT%AMPLITUDEtot,namax+1 ) )

    fstrHEAT%AMPLfuncA  = 0.d0
    fstrHEAT%AMPLfuncB  = 0.d0
    !C
    !C--
    do i = 1, fstrHEAT%AMPLITUDEtot

      fstrHEAT%AMPLfuncA(i,1) = 0.d0
      fstrHEAT%AMPLfuncB(i,1) = fstrHEAT%AMPL(i,1)

      nn = fstrHEAT%AMPLtab(i)

      do k = 2, nn

        x1 = fstrHEAT%AMPLtime(i,k-1)
        y1 = fstrHEAT%AMPL    (i,k-1)
        x2 = fstrHEAT%AMPLtime(i,k)
        y2 = fstrHEAT%AMPL    (i,k)

        fstrHEAT%AMPLfuncA(i,k) =  (y2-y1)/(x2-x1)
        fstrHEAT%AMPLfuncB(i,k) = -(y2-y1)/(x2-x1)*x1 + y1

      enddo

      fstrHEAT%AMPLfuncA(i,nn+1) = 0.d0
      fstrHEAT%AMPLfuncB(i,nn+1) = fstrHEAT%AMPL(i,nn)

    enddo

    !C===
  end subroutine heat_init_amplitude

  !C***
  !C*** INIT_MATERIAL
  !C***
  subroutine heat_init_material (hecMESH, fstrHEAT)

    use m_fstr

    implicit none
    integer(kind=kint) :: m1max, m2max, m3max, icou, im, jm, nn, ic, jS, jE, kc, km, k
    real(kind=kreal)   :: aa, bb
    type(fstr_heat)    :: fstrHEAT
    type(hecmwST_local_mesh) :: hecMESH

    !C
    !C +----------+
    !C | MATERIAL |
    !C +----------+
    !C===
    fstrHEAT%MATERIALtot= hecMESH%material%n_mat

    if(associated(fstrHEAT%mat_name)) deallocate(fstrHEAT%mat_name)
    allocate(fstrHEAT%mat_name(fstrHEAT%MATERIALtot))
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%mat_name(im) = hecMESH%material%mat_name(im)
    enddo

    m1max = 0
    m2max = 0
    m3max = 0
    icou  = 0
    do im = 1, hecMESH%material%n_mat
      do jm = 1, 3
        icou = icou + 1
        nn = hecMESH%material%mat_TABLE_index(icou) - hecMESH%material%mat_TABLE_index(icou-1)
        if( jm.eq.1 ) m1max = max(nn,m1max)
        if( jm.eq.2 ) m2max = max(nn,m2max)
        if( jm.eq.3 ) m3max = max(nn,m3max)
      enddo
    enddo

    if(associated(fstrHEAT%RHOtab)) deallocate(fstrHEAT%RHOtab)
    if(associated(fstrHEAT%CPtab)) deallocate(fstrHEAT%CPtab)
    if(associated(fstrHEAT%CONDtab)) deallocate(fstrHEAT%CONDtab)
    if(associated(fstrHEAT%RHO)) deallocate(fstrHEAT%RHO)
    if(associated(fstrHEAT%RHOtemp)) deallocate(fstrHEAT%RHOtemp)
    if(associated(fstrHEAT%CP)) deallocate(fstrHEAT%CP)
    if(associated(fstrHEAT%CPtemp)) deallocate(fstrHEAT%CPtemp)
    if(associated(fstrHEAT%COND)) deallocate(fstrHEAT%COND)
    if(associated(fstrHEAT%CONDtemp)) deallocate(fstrHEAT%CONDtemp)
    if(associated(fstrHEAT%RHOfuncA)) deallocate(fstrHEAT%RHOfuncA)
    if(associated(fstrHEAT%RHOfuncB)) deallocate(fstrHEAT%RHOfuncB)
    if(associated(fstrHEAT%CPfuncA)) deallocate(fstrHEAT%CPfuncA)
    if(associated(fstrHEAT%CPfuncB)) deallocate(fstrHEAT%CPfuncB)
    if(associated(fstrHEAT%CONDfuncA)) deallocate(fstrHEAT%CONDfuncA)
    if(associated(fstrHEAT%CONDfuncB)) deallocate(fstrHEAT%CONDfuncB)

    allocate (fstrHEAT%RHOtab  (fstrHEAT%MATERIALtot) )
    allocate (fstrHEAT%CPtab   (fstrHEAT%MATERIALtot) )
    allocate (fstrHEAT%CONDtab (fstrHEAT%MATERIALtot) )
    allocate (fstrHEAT%RHO     (fstrHEAT%MATERIALtot,m1max) )
    allocate (fstrHEAT%RHOtemp (fstrHEAT%MATERIALtot,m1max) )
    allocate (fstrHEAT%CP      (fstrHEAT%MATERIALtot,m2max) )
    allocate (fstrHEAT%CPtemp  (fstrHEAT%MATERIALtot,m2max) )
    allocate (fstrHEAT%COND    (fstrHEAT%MATERIALtot,m3max) )
    allocate (fstrHEAT%CONDtemp(fstrHEAT%MATERIALtot,m3max) )

    fstrHEAT%RHO  = 0.d0
    fstrHEAT%CP   = 0.d0
    fstrHEAT%COND = 0.d0

    fstrHEAT%RHOtemp  = 0.d0
    fstrHEAT%CPtemp   = 0.d0
    fstrHEAT%CONDtemp = 0.d0
    fstrHEAT%RHOtab   = 0
    fstrHEAT%CPtab    = 0
    fstrHEAT%CONDtab  = 0

    ic = 0
    do im = 1, fstrHEAT%MATERIALtot
      do jm = 1, 3
        ic = ic + 1
        jS = hecMESH%material%mat_TABLE_index(ic-1) + 1
        jE = hecMESH%material%mat_TABLE_index(ic  )
        nn = jE - jS + 1
        if( jm.eq.1 ) fstrHEAT%RHOtab (im) = nn
        if( jm.eq.2 ) fstrHEAT%CPtab  (im) = nn
        if( jm.eq.3 ) fstrHEAT%CONDtab(im) = nn

        kc = 0
        do km = jS, jE
          kc = kc + 1
          if( jm.eq.1 ) then
            fstrHEAT%RHO     (im,kc) = hecMESH%material%mat_VAL (km)
            fstrHEAT%RHOtemp (im,kc) = hecMESH%material%mat_TEMP(km)
          endif
          if( jm.eq.2 ) then
            fstrHEAT%CP      (im,kc) = hecMESH%material%mat_VAL (km)
            fstrHEAT%CPtemp  (im,kc) = hecMESH%material%mat_TEMP(km)
          endif
          if( jm.eq.3 ) then
            fstrHEAT%COND    (im,kc) = hecMESH%material%mat_VAL (km)
            fstrHEAT%CONDtemp(im,kc) = hecMESH%material%mat_TEMP(km)
          endif
        enddo
      enddo
    enddo
    !C===

    !C
    !C +-----------+
    !C | MAT-TABLE |
    !C +-----------+
    !C===
    allocate (fstrHEAT%RHOfuncA (fstrHEAT%MATERIALtot, m1max+1)              &
      ,fstrHEAT%RHOfuncB (fstrHEAT%MATERIALtot, m1max+1))
    allocate (fstrHEAT%CPfuncA  (fstrHEAT%MATERIALtot, m2max+1)              &
      ,fstrHEAT%CPfuncB  (fstrHEAT%MATERIALtot, m2max+1))
    allocate (fstrHEAT%CONDfuncA(fstrHEAT%MATERIALtot, m3max+1)              &
      ,fstrHEAT%CONDfuncB(fstrHEAT%MATERIALtot, m3max+1))

    fstrHEAT%RHOfuncA  = 0.d0
    fstrHEAT%RHOfuncB  = 0.d0
    fstrHEAT%CPfuncA   = 0.d0
    fstrHEAT%CPfuncB   = 0.d0
    fstrHEAT%CONDfuncA = 0.d0
    fstrHEAT%CONDfuncB = 0.d0
    !C
    !C--RHO
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%RHOfuncB(im,1) = fstrHEAT%RHO(im,1)
      do k = 2, fstrHEAT%RHOtab(im)
        bb= fstrHEAT%RHO    (im,k) - fstrHEAT%RHO    (im,k-1)
        aa= fstrHEAT%RHOtemp(im,k) - fstrHEAT%RHOtemp(im,k-1)
        fstrHEAT%RHOfuncA(im,k) =   bb/aa
        fstrHEAT%RHOfuncB(im,k) = -(bb/aa)*fstrHEAT%RHOtemp(im,k-1) + fstrHEAT%RHO(im,k-1)
      enddo
      fstrHEAT%RHOfuncB(im,fstrHEAT%RHOtab(im)+1) = fstrHEAT%RHO(im,fstrHEAT%RHOtab(im))
    enddo
    !C
    !C-- CP
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%CPfuncB(im,1) = fstrHEAT%CP(im,1)
      do k = 2, fstrHEAT%CPtab(im)
        bb= fstrHEAT%CP    (im,k) - fstrHEAT%CP    (im,k-1)
        aa= fstrHEAT%CPtemp(im,k) - fstrHEAT%CPtemp(im,k-1)
        fstrHEAT%CPfuncA(im,k) =   bb/aa
        fstrHEAT%CPfuncB(im,k) = -(bb/aa)*fstrHEAT%CPtemp(im,k-1) + fstrHEAT%CP(im,k-1)
      enddo
      fstrHEAT%CPfuncB(im,fstrHEAT%CPtab(im)+1) =  fstrHEAT%CP(im,fstrHEAT%CPtab(im))
    enddo
    !C
    !C-- COND.
    do im = 1, fstrHEAT%MATERIALtot
      fstrHEAT%CONDfuncB(im,1)= fstrHEAT%COND(im,1)
      do k = 2, fstrHEAT%CONDtab(im)
        bb = fstrHEAT%COND    (im,k) - fstrHEAT%COND    (im,k-1)
        aa = fstrHEAT%CONDtemp(im,k) - fstrHEAT%CONDtemp(im,k-1)
        fstrHEAT%CONDfuncA(im,k) =   bb/aa
        fstrHEAT%CONDfuncB(im,k) = -(bb/aa)*fstrHEAT%CONDtemp(im,k-1) + fstrHEAT%COND(im,k-1)
      enddo
      fstrHEAT%CONDfuncB(im,fstrHEAT%CONDtab(im)+1) = fstrHEAT%COND(im,fstrHEAT%CONDtab(im))
    enddo
    !C===
  end subroutine heat_init_material

  !-----------------------------------------------------------------------------!
  !> Read in !HEAT                                                             !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_HEAT( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    integer(kind=kint) :: n
    character(len=HECMW_NAME_LEN) :: mName
    integer(kind=kint) :: i

    n = fstr_ctrl_get_data_line_n( ctrl )

    if( n == 0 ) return

    call reallocate_real( P%PARAM%dtime, n)
    call reallocate_real( P%PARAM%etime, n)
    call reallocate_real( P%PARAM%dtmin, n)
    call reallocate_real( P%PARAM%delmax,n)
    call reallocate_integer( P%PARAM%itmax, n)
    call reallocate_real( P%PARAM%eps,   n)
    P%PARAM%analysis_n = n

    P%PARAM%dtime = 0
    P%PARAM%etime = 0
    P%PARAM%dtmin = 0
    P%PARAM%delmax = 0
    P%PARAM%itmax = 20
    P%PARAM%eps = 1.0e-6
    P%PARAM%timepoint_id = 0

    rcode = fstr_ctrl_get_HEAT(   ctrl,        &
      P%PARAM%dtime,     &
      P%PARAM%etime,     &
      P%PARAM%dtmin,     &
      P%PARAM%delmax,    &
      P%PARAM%itmax,     &
      P%PARAM%eps,       &
      mName,             &
      P%HEAT%beta)
    if( rcode /= 0 ) then
      call fstr_ctrl_err_stop
    end if

    if( associated(P%PARAM%timepoints) ) then
      do i=1,size(P%PARAM%timepoints)
        if( hecmw_streqr( P%PARAM%timepoints(i)%name, mName ) ) then
          P%PARAM%timepoint_id = i; exit
        endif
      enddo
    endif

    call reallocate_real( P%HEAT%STEP_DLTIME, n)
    call reallocate_real( P%HEAT%STEP_EETIME, n)
    call reallocate_real( P%HEAT%STEP_DELMIN, n)
    call reallocate_real( P%HEAT%STEP_DELMAX, n)
    P%HEAT%STEPtot = n

    P%HEAT%STEP_DLTIME = P%PARAM%dtime
    P%HEAT%STEP_EETIME = P%PARAM%etime
    P%HEAT%STEP_DELMIN = P%PARAM%dtmin
    P%HEAT%STEP_DELMAX = P%PARAM%delmax
    P%HEAT%timepoint_id = P%PARAM%timepoint_id

  end subroutine fstr_setup_HEAT

  !-----------------------------------------------------------------------------!
  !> Read in !FIXTEMP                                                          !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_FIXTEMP( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_FIXTEMP( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_fixtemp( P%MESH, P%PARAM, P%HEAT, '!FIXTEMP', amp, n, grp_id_name, value )

    deallocate( grp_id_name )
    deallocate( value )
  end subroutine fstr_setup_FIXTEMP

  !-----------------------------------------------------------------------------!
  !> Read in !CFLUX                                                            !
  !-----------------------------------------------------------------------------!

  subroutine fstr_setup_CFLUX( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_CFLUX( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_cflux( P%MESH, P%HEAT, '!CFLUX', amp, n, grp_id_name, value )

    deallocate( grp_id_name )
    deallocate( value )
  end subroutine fstr_setup_CFLUX

  !-----------------------------------------------------------------------------!
  !> Read in !DFLUX                                                            !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_DFLUX( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( load_type(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_DFLUX( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, load_type, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_dflux( P%MESH, P%HEAT, '!DFLUX', amp, n, grp_id_name, load_type, value )

    deallocate( grp_id_name )
    deallocate( load_type )
    deallocate( value )
  end subroutine fstr_setup_DFLUX

  !-----------------------------------------------------------------------------!
  !> Read in !SFLUX                                                            !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_SFLUX( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))

    amp = ' '
    rcode = fstr_ctrl_get_SFLUX( ctrl, amp, &
      grp_id_name, HECMW_NAME_LEN, value )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_sflux( P%MESH, P%HEAT, '!SFLUX', amp, n, grp_id_name, value )

    deallocate( grp_id_name )
    deallocate( value )
  end subroutine fstr_setup_SFLUX

  !-----------------------------------------------------------------------------!
  !> Read in !FILM                                                             !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_FILM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( load_type(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '

    rcode = fstr_ctrl_get_FILM( ctrl, amp1, amp2, &
      grp_id_name, HECMW_NAME_LEN, load_type, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_film( P%MESH, P%HEAT, '!FILM', amp1, amp2, &
      n, grp_id_name, load_type, value, shink )

    deallocate( grp_id_name )
    deallocate( load_type )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_FILM

  !-----------------------------------------------------------------------------!
  !> Read in !SFILM                                                            !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_SFILM( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '
    rcode = fstr_ctrl_get_SFILM( ctrl, amp1, amp2, &
      grp_id_name, HECMW_NAME_LEN, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_sfilm( P%MESH, P%HEAT, '!SFILM', amp1, amp2, &
      n, grp_id_name, value, shink )

    deallocate( grp_id_name )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_SFILM

  !-----------------------------------------------------------------------------!
  !> Read in !RADIATE                                                          !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_RADIATE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    integer(kind=kint),pointer :: load_type(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( load_type(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '
    rcode = fstr_ctrl_get_RADIATE( ctrl, amp1, amp2, &
      grp_id_name, HECMW_NAME_LEN, load_type, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_radiate( P%MESH, P%HEAT, '!RADIATE', amp1, amp2, &
      n, grp_id_name, load_type, value, shink )

    deallocate( grp_id_name )
    deallocate( load_type )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_RADIATE

  !-----------------------------------------------------------------------------!
  !> Read in !SRADIATE                                                         !
  !-----------------------------------------------------------------------------!


  subroutine fstr_setup_SRADIATE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    integer(kind=kint) :: rcode
    character(HECMW_NAME_LEN) :: amp1, amp2
    character(HECMW_NAME_LEN), pointer :: grp_id_name(:)
    real(kind=kreal),pointer :: value(:)
    real(kind=kreal),pointer :: shink(:)
    integer(kind=kint) :: n
    ! ------------------------------------------------

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return

    allocate( grp_id_name(n))
    allocate( value(n))
    allocate( shink(n))

    amp1 = ' '
    amp2 = ' '
    rcode = fstr_ctrl_get_SRADIATE( ctrl, amp1, amp2, grp_id_name, HECMW_NAME_LEN, value, shink )
    if( rcode /= 0 ) call fstr_ctrl_err_stop

    call fstr_heat_append_sradiate( P%MESH, P%HEAT, '!SRADIATE', amp1, amp2, &
      n, grp_id_name, value, shink )

    deallocate( grp_id_name )
    deallocate( value )
    deallocate( shink )
  end subroutine fstr_setup_SRADIATE

  !-----------------------------------------------------------------------------!
  !> Read in !WELD_LINE                                                         !
  !-----------------------------------------------------------------------------!
  subroutine fstr_setup_WELD_LINE( ctrl, P )
    implicit none
    integer(kind=kint) :: ctrl
    type(fstr_param_pack) :: P

    type(hecmwST_local_mesh), pointer :: hecMESH
    type(fstr_heat), pointer :: fstrHEAT
    type(fstr_card_slot) :: slot

    slot = P%STATE%slot
    hecMESH => P%MESH
    fstrHEAT => P%HEAT

    fstrHEAT%WL_tot = slot%last
    if( fstr_ctrl_get_WELDLINE( ctrl, hecMESH, HECMW_NAME_LEN, fstrHEAT%weldline(fstrHEAT%WL_tot) )/=0 ) then
      write(*,*) '### Error: Fail in read in Weld Line definition : ' , fstrHEAT%WL_tot
      write(ILOG,*) '### Error: Fail in read in Weld Line definition : ', fstrHEAT%WL_tot
      stop HECMW_EXIT_INPUT
    endif
  end subroutine fstr_setup_WELD_LINE

  !-----------------------------------------------------------------------------!
  !> Append !FIXTEMP                                                            !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_fixtemp( hecMESH, fstrPARAM, fstrHEAT, header_name, &
      amp_name, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_param), target                  :: fstrPARAM
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0
    do i = 1, n
      !rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      rtc = get_sorted_local_member_index( hecMESH, fstrPARAM, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'node_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-8
    old_size = fstrHEAT%T_FIX_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%T_FIX_node, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%T_FIX_ampl, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%T_FIX_val,  old_size, new_size )
    fstrHEAT%T_FIX_tot = new_size

    head = old_size + 1
    member => fstrHEAT%T_FIX_node(head:)
    id = head
    do i = 1, n
      !rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      rtc = get_sorted_local_member_index( hecMESH, fstrPARAM, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'node_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) then
        member => member( member_n+1 : )
      endif
      do j = 1, member_n
        fstrHEAT%T_FIX_val  (id) = vals(i)
        fstrHEAT%T_FIX_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_fixtemp

  !-----------------------------------------------------------------------------!
  !> Append !CFLUX                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_cflux( hecMESH, fstrHEAT, header_name, amp_name, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0

    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'node_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-9
    old_size = fstrHEAT%Q_NOD_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%Q_NOD_node, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_NOD_ampl, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%Q_NOD_val,  old_size, new_size )
    fstrHEAT%Q_NOD_tot = new_size

    head = old_size + 1
    member => fstrHEAT%Q_NOD_node(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'node', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'node_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%Q_NOD_val  (id) = vals(i)
        fstrHEAT%Q_NOD_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_cflux

  !-----------------------------------------------------------------------------!
  !> Append !DFLUX                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_dflux( hecMESH, fstrHEAT, header_name, amp_name, n, grp_names, &
      load_types, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: load_types(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'elem_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-10
    old_size = fstrHEAT%Q_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_elem, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_ampl, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_surf, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%Q_SUF_val,  old_size, new_size )
    fstrHEAT%Q_SUF_tot = new_size

    head = old_size + 1
    member => fstrHEAT%Q_SUF_elem(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'elem_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%Q_SUF_surf (id) = load_types(i)
        fstrHEAT%Q_SUF_val  (id) = vals(i)
        fstrHEAT%Q_SUF_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_dflux

  !-----------------------------------------------------------------------------!
  !> Append !SFLUX                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_sflux( hecMESH, fstrHEAT, header_name, amp_name, n, grp_names, vals )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)

    integer(kind=kint) :: amp_id
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member1(:), member2(:)

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name, amp_id )

    m = 0
    do i = 1, n
      m = m + get_grp_member_n( hecMESH, 'surf_grp', grp_names(i) )
    end do

    if (m == 0) return

    ! JP-11
    old_size = fstrHEAT%Q_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_elem, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_ampl, old_size, new_size )
    call hecmw_expand_integer_array( fstrHEAT%Q_SUF_surf, old_size, new_size )
    call hecmw_expand_real_array(    fstrHEAT%Q_SUF_val,  old_size, new_size )
    fstrHEAT%Q_SUF_tot = new_size

    head = old_size + 1
    member1 => fstrHEAT%Q_SUF_elem(head:)
    member2 => fstrHEAT%Q_SUF_surf(head:)
    id = head
    do i = 1, n
      member_n = get_grp_member( hecMESH, 'surf_grp', grp_names(i), member1, member2 )
      if( i<n ) then
        member1 => member1( member_n+1 : )
        member2 => member2( member_n+1 : )
      end if
      do j = 1, member_n
        fstrHEAT%Q_SUF_val  (id) = vals(i)
        fstrHEAT%Q_SUF_ampl (id) = amp_id
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_sflux

  !-----------------------------------------------------------------------------!
  !> Append !FILM                                                               !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_film( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, load_types, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: load_types(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'elem_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-12
    old_size = fstrHEAT%H_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%H_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%H_SUF_val,  2, old_size, new_size )
    fstrHEAT%H_SUF_tot = new_size

    head = old_size + 1
    member => fstrHEAT%H_SUF_elem(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'elem_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%H_SUF_surf (id)   = load_types(i)
        fstrHEAT%H_SUF_val  (id,1) = vals(i)
        fstrHEAT%H_SUF_val  (id,2) = sinks(i)
        fstrHEAT%H_SUF_ampl (id,1) = amp_id1
        fstrHEAT%H_SUF_ampl (id,2) = amp_id2
        id= id + 1
      end do
    end do

  end subroutine fstr_heat_append_film

  !-----------------------------------------------------------------------------!
  !> Append !SFILM                                                              !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_sfilm( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member1(:), member2(:)

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      m = m + get_grp_member_n( hecMESH, 'surf_grp', grp_names(i) )
    end do

    if (m == 0) return

    ! JP-13
    old_size = fstrHEAT%H_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%H_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%H_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%H_SUF_val,  2, old_size, new_size )
    fstrHEAT%H_SUF_tot = new_size

    head = old_size + 1
    member1 => fstrHEAT%H_SUF_elem(head:)
    member2 => fstrHEAT%H_SUF_surf(head:)
    id = head
    do i = 1, n
      member_n = get_grp_member( hecMESH, 'surf_grp', grp_names(i), member1, member2 )
      if( i<n ) then
        member1 => member1( member_n+1 : )
        member2 => member2( member_n+1 : )
      end if
      do j = 1, member_n
        fstrHEAT%H_SUF_val  (id,1) = vals(i)
        fstrHEAT%H_SUF_val  (id,2) = sinks(i)
        fstrHEAT%H_SUF_ampl (id,1) = amp_id1
        fstrHEAT%H_SUF_ampl (id,2) = amp_id2
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_sfilm

  !-----------------------------------------------------------------------------!
  !> Append !RADIATE                                                            !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_radiate( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, load_types, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    integer(kind=kint), intent(in)            :: load_types(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member(:)
    integer(kind=kint) :: local_id, rtc

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        m = m + 1
      else if( rtc < 0 ) then
        m = m + get_grp_member_n( hecMESH, 'elem_grp', grp_names(i) )
      end if
    end do

    if (m == 0) return

    ! JP-14
    old_size = fstrHEAT%R_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%R_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%R_SUF_val,  2, old_size, new_size )
    fstrHEAT%R_SUF_tot = new_size

    head = old_size + 1
    member => fstrHEAT%R_SUF_elem(head:)
    id = head
    do i = 1, n
      rtc = get_local_member_index( hecMESH, 'element', grp_names(i), local_id )
      if( rtc > 0 ) then
        member(1) = local_id
        member_n = 1
      else if( rtc < 0 ) then
        member_n = get_grp_member( hecMESH, 'elem_grp', grp_names(i), member )
      else
        cycle
      end if
      if( i<n ) member => member( member_n+1 : )
      do j = 1, member_n
        fstrHEAT%R_SUF_surf (id)   = load_types(i)
        fstrHEAT%R_SUF_val  (id,1) = vals(i)
        fstrHEAT%R_SUF_val  (id,2) = sinks(i)
        fstrHEAT%R_SUF_ampl (id,1) = amp_id1
        fstrHEAT%R_SUF_ampl (id,2) = amp_id2
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_radiate

  !-----------------------------------------------------------------------------!
  !> Append !SRADIATE                                                           !
  !-----------------------------------------------------------------------------!
  subroutine fstr_heat_append_sradiate( hecMESH, fstrHEAT, header_name, amp_name1, amp_name2, &
      n, grp_names, vals, sinks )
    implicit none
    type(hecmwST_local_mesh), target          :: hecMESH
    type(fstr_heat)                           :: fstrHEAT
    character(len=*), intent(in)              :: header_name
    character(len=HECMW_NAME_LEN), intent(in) :: amp_name1, amp_name2
    integer(kind=kint), intent(in)            :: n
    character(len=HECMW_NAME_LEN), intent(in) :: grp_names(:)
    real(kind=kreal), intent(in)              :: vals(:)
    real(kind=kreal), intent(in)              :: sinks(:)

    integer(kind=kint) :: amp_id1, amp_id2
    integer(kind=kint) :: i, j, m, head, id, member_n, old_size, new_size
    integer(kind=kint),pointer :: member1(:), member2(:)

    if( n == 0 ) return
    call amp_name_to_id( hecMESH, header_name, amp_name1, amp_id1 )
    call amp_name_to_id( hecMESH, header_name, amp_name2, amp_id2 )

    m = 0
    do i = 1, n
      m = m + get_grp_member_n( hecMESH, 'surf_grp', grp_names(i) )
    end do

    if (m == 0) return

    ! JP-15
    old_size = fstrHEAT%R_SUF_tot
    new_size = old_size + m
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_elem,    old_size, new_size )
    call hecmw_expand_integer_array2( fstrHEAT%R_SUF_ampl, 2, old_size, new_size )
    call hecmw_expand_integer_array(  fstrHEAT%R_SUF_surf,    old_size, new_size )
    call hecmw_expand_real_array2(    fstrHEAT%R_SUF_val,  2, old_size, new_size )
    fstrHEAT%R_SUF_tot = new_size

    head = old_size + 1
    member1 => fstrHEAT%R_SUF_elem(head:)
    member2 => fstrHEAT%R_SUF_surf(head:)
    id = head
    do i = 1, n
      member_n = get_grp_member( hecMESH, 'surf_grp', grp_names(i), member1, member2 )
      if( i<n ) then
        member1 => member1( member_n+1 : )
        member2 => member2( member_n+1 : )
      end if
      do j = 1, member_n
        fstrHEAT%R_SUF_val  (id,1) = vals(i)
        fstrHEAT%R_SUF_val  (id,2) = sinks(i)
        fstrHEAT%R_SUF_ampl (id,1) = amp_id1
        fstrHEAT%R_SUF_ampl (id,2) = amp_id2
        id = id + 1
      end do
    end do

  end subroutine fstr_heat_append_sradiate

end module m_fstr_heat_setup
