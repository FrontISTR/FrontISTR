!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module manages read in of various material properties
module fstr_ctrl_material
  use hecmw
  use mMaterial
  use m_table
  implicit none

  private :: read_user_matl

  include 'fstr_ctrl_util_f.inc'

contains


  !----------------------------------------------------------------------
  !> Read in !MATERIAL
  integer function fstr_ctrl_get_MATERIAL( ctrl, matname )
    integer(kind=kint), intent(in) :: ctrl
    character(len=*), intent(out) :: matname

    matname=""
    fstr_ctrl_get_MATERIAL = fstr_ctrl_get_param_ex( ctrl, 'NAME ',  '# ',  0, 'S', matname )
  end function fstr_ctrl_get_MATERIAL

  !----------------------------------------------------------------------
  !> Read in !USER_MATERIAL
  integer function fstr_ctrl_get_USERMATERIAL( ctrl, mattype, nlgeom, nstatus, matval )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    integer(kind=kint), intent(out)   :: nstatus
    real(kind=kreal),intent(out)      :: matval(:)

    integer(kind=kint) :: ipt
    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=256) :: s, fname

    fstr_ctrl_get_USERMATERIAL = -1
    mattype = USERMATERIAL
    nlgeom = UPDATELAG   !default value
    nstatus = 1
    if( fstr_ctrl_get_param_ex( ctrl, 'NSTATUS ',  '# ',    0,   'I',   nstatus )/= 0) return
    if( fstr_ctrl_get_param_ex( ctrl, 'KIRCHHOFF ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = TOTALLAG

    fstr_ctrl_get_USERMATERIAL = read_user_matl( ctrl, matval )
  end function fstr_ctrl_get_USERMATERIAL

  !----------------------------------------------------------------------
  !> Read in !ELASTIC
  integer function fstr_ctrl_get_ELASTICITY( ctrl, mattype, nlgeom, matval, dict )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    real(kind=kreal),intent(out)      :: matval(:)
    type(DICT_STRUCT), pointer        :: dict

    integer(kind=kint) :: i,j, rcode, depends, ipt, n
    real(kind=kreal),pointer :: fval(:,:)
    character(len=HECMW_NAME_LEN) :: data_fmt
    type( tTable )        :: mattable
    logical            :: isok
    character(len=256) :: s

    fstr_ctrl_get_ELASTICITY = -1
    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends=1   ! temperature depends only currently
    if( depends > 3 ) stop "We cannot read dependencies>3 right now"
    nlgeom = TOTALLAG   !default value

    if( fstr_ctrl_get_param_ex( ctrl, 'CAUCHY ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = UPDATELAG
    if( fstr_ctrl_get_param_ex( ctrl, 'INFINITESIMAL ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = INFINITESIMAL

    ! for backward compatibility
    if( fstr_ctrl_get_param_ex( ctrl, 'INFINITE ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) then
      write(*,*) "Warning : !ELASTIC : parameter 'INFINITE' is deprecated." &
           & //  " Please use the replacement parameter 'INFINITESIMAL'"
      nlgeom = INFINITESIMAL
    endif

    ipt=1
    s = 'ISOTROPIC,ORTHOTROPIC,USER '
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',  s, 0, 'P',   ipt    ) /= 0 ) return

    n = fstr_ctrl_get_data_line_n( ctrl )
    ! ISOTROPIC
    if( ipt==1 ) then
      allocate( fval(2+depends,n) )
      if( depends==0 ) then
        data_fmt = "RR "
        fstr_ctrl_get_ELASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
      endif
      if( depends==1 ) then
        data_fmt = "RRR "
        fstr_ctrl_get_ELASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
      endif
      if( fstr_ctrl_get_ELASTICITY ==0 ) then
        matval(M_YOUNGS) = fval(1,1)
        matval(M_POISSON) = fval(2,1)
        call init_table( mattable, depends, 2+depends,n, fval )
        call dict_add_key( dict, MC_ISOELASTIC, mattable )
        !         call print_table( mattable, 6 ); pause
      endif
      mattype = ELASTIC

      ! ORTHOTROPIC
    else if( ipt==2 ) then
      allocate( fval(9+depends,n) )
      if( depends==0 ) then
        data_fmt = "RRRRRRRRR "
        fstr_ctrl_get_ELASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt,    &
          fval(1,:), fval(2,:), fval(3,:), fval(4,:), fval(5,:), fval(6,:),           &
          fval(7,:), fval(8,:), fval(9,:) )
      else if( depends==1 ) then
        data_fmt = "RRRRRRRRRR "
        fstr_ctrl_get_ELASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt,    &
          fval(1,:), fval(2,:), fval(3,:), fval(4,:), fval(5,:), fval(6,:),           &
          fval(7,:), fval(8,:), fval(9,:), fval(10,:) )
      endif
      if( fstr_ctrl_get_ELASTICITY ==0 ) then
        isok = .true.
        do i=1,n
          if( fval(1,i)<=0.d0 .or. fval(2,i)<=0.d0 .or. fval(3,i)<=0.d0  .or.  &
              fval(7,i)<=0.d0 .or. fval(8,i)<=0.d0 .or. fval(9,i)<=0.d0 ) then
            isok = .false.;  fstr_ctrl_get_ELASTICITY=1; exit
          endif
        enddo
        if( isok ) then
          call init_table( mattable, depends, 9+depends,n, fval )
          call dict_add_key( dict, MC_ORTHOELASTIC, mattable )
          mattype = MN_ORTHOELASTIC
        endif
      endif

    else if( ipt==3 ) then
      fstr_ctrl_get_ELASTICITY = read_user_matl( ctrl, matval(101:200))
      mattype = USERELASTIC
      nlgeom = INFINITESIMAL

    else
      stop "ERROR: Material type not supported"

    endif

    call finalize_table( mattable )
    if( associated(fval) ) deallocate(fval)
  end function fstr_ctrl_get_ELASTICITY


  !----------------------------------------------------------------------
  !> Read in !HYPERELASTIC
  integer function fstr_ctrl_get_HYPERELASTIC( ctrl, mattype, nlgeom, matval )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    real(kind=kreal),intent(out)      :: matval(:)

    integer(kind=kint) :: i,j, rcode, depends, ipt
    real(kind=kreal),pointer :: fval(:,:)
    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=256) :: s

    fstr_ctrl_get_HYPERELASTIC = -1
    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends > 3 ) stop "We cannot read dependencies>3 right now"
    nlgeom = TOTALLAG   !default value
    if( fstr_ctrl_get_param_ex( ctrl, 'CAUCHY ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = UPDATELAG

    ipt=1
    s = 'NEOHOOKE,MOONEY-RIVLIN,ARRUDA-BOYCE,USER,MOONEY-RIVLIN-ANISO '
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',  s, 0, 'P',   ipt    ) /= 0 ) return

    ! NEOHOOKE
    if( ipt==1 ) then
      allocate( fval(2,depends+1) )
      fval =0.0d0
      if( depends==0 ) then
        data_fmt = "RR "
        fstr_ctrl_get_HYPERELASTIC =                           &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
      endif
      if( fstr_ctrl_get_HYPERELASTIC ==0 ) then
        if( fval(2,1)==0.d0 ) stop "We cannot deal with incompressible material currently"
        matval(M_PLCONST1) = fval(1,1)
        matval(M_PLCONST2) = 0.d0
        matval(M_PLCONST3) = fval(2,1)
        ! matval(M_YOUNGS) = fval(1,1)
        ! matval(M_POISSON) = fval(2,1)
      endif
      mattype = NEOHOOKE

      ! MOONEY
    else if( ipt==2 ) then
      allocate( fval(3,depends+1) )
      fval =0.0d0
      if( depends==0 ) then
        data_fmt = "RRR "
        fstr_ctrl_get_HYPERELASTIC =                       &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
      endif
      if( fstr_ctrl_get_HYPERELASTIC ==0 ) then
        matval(M_PLCONST1) = fval(1,1)
        matval(M_PLCONST2) = fval(2,1)
        matval(M_PLCONST3) = fval(3,1)
      endif
      mattype = MOONEYRIVLIN

      ! ARRUDA
    else if( ipt==3 ) then
      allocate( fval(3,depends+1) )
      fval =0.0d0
      if( depends==0 ) then
        data_fmt = "RRR "
        fstr_ctrl_get_HYPERELASTIC = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
      endif
      if( fstr_ctrl_get_HYPERELASTIC ==0 ) then
        matval(M_PLCONST1) = fval(1,1)
        matval(M_PLCONST2) = fval(2,1)
        matval(M_PLCONST3) = fval(3,1)
      endif
      mattype = ARRUDABOYCE

    else if( ipt==4 ) then    !User
      fstr_ctrl_get_HYPERELASTIC = read_user_matl( ctrl, matval(101:200))
      mattype = USERHYPERELASTIC

      ! MOONEY-ORTHO
    else if( ipt==5 ) then
      allocate( fval(10,depends+1) )
      fval =0.0d0
      if( depends==0 ) then
        data_fmt = "RRRRRrrrrr "
        fstr_ctrl_get_HYPERELASTIC = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, &
          fval(1,:), fval(2,:), fval(3,:), fval(4,:), fval(5,:), &
          fval(6,:), fval(7,:), fval(8,:), fval(9,:), fval(10,:) )
      endif
      if( fstr_ctrl_get_HYPERELASTIC ==0 ) then
        matval(M_PLCONST1) = fval(1,1)
        matval(M_PLCONST2) = fval(2,1)
        matval(M_PLCONST3) = fval(3,1)
        matval(M_PLCONST4) = fval(4,1)
        matval(M_PLCONST5) = fval(5,1)
        matval(M_PLCONST6) = fval(6,1)
        matval(M_PLCONST7) = fval(7,1)
        matval(M_PLCONST8) = fval(8,1)
        matval(M_PLCONST9) = fval(9,1)
        matval(M_PLCONST10) = fval(10,1)
      endif
      mattype = MOONEYRIVLIN_ANISO

    endif

    if( associated(fval) ) deallocate(fval)
  end function fstr_ctrl_get_HYPERELASTIC


  !----------------------------------------------------------------------
  !> Read in !VISCOELASTIC
  integer function fstr_ctrl_get_VISCOELASTICITY( ctrl, mattype, nlgeom, dict )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    type(DICT_STRUCT), pointer        :: dict

    integer(kind=kint) :: i,j, rcode, depends, ipt, n
    real(kind=kreal),pointer :: fval(:,:)
    character(len=HECMW_NAME_LEN) :: data_fmt
    type( tTable )        :: mattable
    character(len=256) :: s

    fstr_ctrl_get_VISCOELASTICITY = -1
    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends=1   ! temperature depends only currently
    !depends = 0
    nlgeom = TOTALLAG   !default value
    if( fstr_ctrl_get_param_ex( ctrl, 'INFINITESIMAL ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = INFINITESIMAL

    ! for backward compatibility
    if( fstr_ctrl_get_param_ex( ctrl, 'INFINITE ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) then
      write(*,*) "Warning : !VISCOELASTIC : parameter 'INFINITE' is deprecated." &
           & //  " Please use the replacement parameter 'INFINITESIMAL'"
      nlgeom = INFINITESIMAL
    endif

    ipt=1
    s = 'ISOTROPIC,USER '
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',  s, 0, 'P',   ipt    ) /= 0 ) return
    ipt = 1

    ! ISOTROPIC
    if( ipt==1 ) then
      n = fstr_ctrl_get_data_line_n( ctrl )
      allocate( fval(2+depends,n) )
      if( depends==0 ) then
        data_fmt = "RR "
        fstr_ctrl_get_VISCOELASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
        if( fval(2,1)==0.d0 ) stop "Error in defining viscoelasticity: Relaxation time cannot be zero!"
      endif
      if( depends==1 ) then
        data_fmt = "RRR "
        fstr_ctrl_get_VISCOELASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
      endif
      if( fstr_ctrl_get_VISCOELASTICITY ==0 ) then
        call init_table( mattable, 1, 2+depends,n, fval )
        call dict_add_key( dict, MC_VISCOELASTIC, mattable )
        !         call print_table( mattable, 6 ); pause
      endif
      mattype = VISCOELASTIC

    else
      stop "ERROR: Material type not supported"

    endif

    call finalize_table( mattable )
    if( associated(fval) ) deallocate(fval)
  end function fstr_ctrl_get_VISCOELASTICITY


  !> Read in !TRS
  integer function fstr_ctrl_get_TRS( ctrl, mattype, matval )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    real(kind=kreal),intent(out)      :: matval(:)

    integer :: ipt
    character(len=256) :: s

    ipt=1
    s = 'WLF,ARRHENIUS '
    if( fstr_ctrl_get_param_ex( ctrl, 'DEFINITION ',  s, 0, 'P',   ipt    ) /= 0 ) return

    fstr_ctrl_get_TRS                                                           &
      = fstr_ctrl_get_data_ex( ctrl, 1, "RRR ", matval(1), matval(2), matval(3) )
    if( fstr_ctrl_get_TRS/=0 ) return
    mattype = mattype+ipt

  end function fstr_ctrl_get_TRS


  !----------------------------------------------------------------------
  !> Read in !PLASTIC
  integer function fstr_ctrl_get_PLASTICITY( ctrl, mattype, nlgeom, matval, mattable, dict )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    real(kind=kreal),intent(out)      :: matval(:)
    real(kind=kreal), pointer         :: mattable(:)
    type(DICT_STRUCT), pointer        :: dict

    integer(kind=kint) :: i, n, rcode, depends, ipt, hipt
    real(kind=kreal),pointer :: fval(:,:) => null()
    real(kind=kreal) :: phi, psi
    character(len=HECMW_NAME_LEN) :: data_fmt
    character(len=256)    :: s
    type( tTable )        :: mttable
    real(kind=kreal), parameter :: PI=3.14159265358979d0

    fstr_ctrl_get_PLASTICITY = -1
    ipt = 0; hipt = 0

    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends = 1 ! we consider temperature dependence only currently
    if( depends > 3 ) stop "We cannot read dependencies>3 right now"
    nlgeom = UPDATELAG   !default value
    if( fstr_ctrl_get_param_ex( ctrl, 'KIRCHHOFF ',  '# ',    0,   'E',   ipt )/= 0) return
    ! rcode = fstr_ctrl_get_param_ex( ctrl, 'FILE  ', '# ',           0,   'S',   fname )
    if( ipt/=0 ) nlgeom = TOTALLAG
    if( fstr_ctrl_get_param_ex( ctrl, 'INFINITESIMAL ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = INFINITESIMAL

    ! for backward compatibility
    if( fstr_ctrl_get_param_ex( ctrl, 'INFINITE ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) then
      write(*,*) "Warning : !PLASTIC : parameter 'INFINITE' is deprecated." &
           & //  " Please use the replacement parameter 'INFINITESIMAL'"
      nlgeom = INFINITESIMAL
    endif

    call setDigit( 1, 1, mattype )
    call setDigit( 2, 2, mattype )

    ! hardening
    s = 'BILINEAR,MULTILINEAR,SWIFT,RAMBERG-OSGOOD,KINEMATIC,COMBINED '
    if( fstr_ctrl_get_param_ex( ctrl, 'HARDEN ',  s , 0, 'P',   hipt    ) /= 0 ) return
    if( hipt==0 ) hipt=1  ! default: linear hardening
    call setDigit( 5, hipt-1, mattype )

    ! yield function
    s = 'MISES,MOHR-COULOMB,DRUCKER-PRAGER,USER '
    call setDigit( 2, 2, mattype )
    if( fstr_ctrl_get_param_ex( ctrl, 'YIELD ',  s , 0, 'P',   ipt    ) /= 0 ) return
    if( ipt==0 ) ipt=1  ! default: mises yield function
    call setDigit( 4, ipt-1, mattype )

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return               ! fail in reading plastic
    if( hipt==2 .and. n<2 ) return    ! not enough data
    if( ( ipt==2 .or. ipt==3 ) .and. hipt>2 ) hipt = 1

    select case (ipt)
      case (1)  !Mises
        select case (hipt)
          case (1,5)  ! linear hardening, kinematic hardening
            allocate( fval(2,n) )
            data_fmt = "RR "
            fstr_ctrl_get_PLASTICITY = &
              fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
            if( fstr_ctrl_get_PLASTICITY ==0 ) then
              matval(M_PLCONST1) = fval(1,1)
              if(hipt==1) then
                matval(M_PLCONST2) = fval(2,1)
              else
                matval(M_PLCONST2) = 0.d0
                matval(M_PLCONST3) = fval(2,1)
              endif
            endif
          case (2)  ! multilinear approximation
            allocate( fval(depends+2,n) )
            if( depends==0 ) then
              data_fmt = "RR "
              fstr_ctrl_get_PLASTICITY = &
                fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
              if( fstr_ctrl_get_PLASTICITY ==0 ) then
                if( fval(2,1)/=0.d0 ) then
                  print *, "Multilinear hardening: First plastic strain must be zero"
                  stop
                endif
                do i=1,n
                  if( fval(2,i)<0.0 ) &
                    stop "Multilinear hardening: Error in plastic strain definition"
                enddo
                call init_table( mttable,1, 2+depends, n, fval )
                call dict_add_key( dict, MC_YIELD, mttable )

              endif
            else  ! depends==1
              data_fmt = "RRR "
              fstr_ctrl_get_PLASTICITY = &
                fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
              if( fstr_ctrl_get_PLASTICITY ==0 ) then
                call init_table( mttable,2, 2+depends,n, fval )
                call dict_add_key( dict, MC_YIELD, mttable )
              endif
            endif
          case (3, 4, 6)  ! swift, Ramberg-Osgood, Combined
            allocate( fval(3,1) )
            data_fmt = "RRR "
            fstr_ctrl_get_PLASTICITY = &
              fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
            if( fstr_ctrl_get_PLASTICITY ==0 ) then
              matval(M_PLCONST1) = fval(1,1)
              matval(M_PLCONST2) = fval(2,1)
              matval(M_PLCONST3) = fval(3,1)
            endif
          case default
            print *, "Error in hardening definition!"
            stop
        end select
      case (2, 3)  ! Mohr-Coulomb, Drucker-Prager
        select case (hipt)
        case (1)  ! linear hardening
          allocate( fval(4,n) )
          data_fmt = "RRrr "
          fval(4,:) = -1.d0
          fstr_ctrl_get_PLASTICITY                                                         &
              = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:), fval(4,:) )
          if( fstr_ctrl_get_PLASTICITY ==0 ) then
            matval(M_PLCONST1) = fval(1,1)    ! c
            matval(M_PLCONST2) = fval(3,1)    ! H
            phi = fval(2,1)*PI/180.d0
            if( fval(4,1) >= 0.d0 ) then
              psi = fval(4,1)*PI/180.d0
            else
              psi = phi
            endif
          endif
        case (2)  ! multilinear hardening
          if( depends>0 ) then
            stop "Mohr-Coulomb and Drucker-Prager do not support temperature dependency"
          endif
          allocate( fval(2,n) )
          data_fmt = "Rr "
          fval(2,:) = -1.d0
          fstr_ctrl_get_PLASTICITY = &
              fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
          if( fstr_ctrl_get_PLASTICITY ==0 ) then
              phi =fval(1,1)*PI/180.d0
              if( fval(2,1) >= 0.d0 ) then
                psi = fval(2,1)*PI/180.d0
              else
                psi = phi
              endif
              if( fval(2,2)/=0.d0 ) then
                print *, "Multilinear hardening: First plastic strain must be zero"
                stop
              endif
              do i=2,n
                if( fval(2,i)<0.0 ) &
                    stop "Multilinear hardening: Error in plastic strain definition"
              enddo
            call init_table( mttable,1, 2, n-1, fval(1:2,2:n) )
            call dict_add_key( dict, MC_YIELD, mttable )
          endif
        end select
        if( ipt==3 ) then     ! Drucker-Prager
          matval(M_PLCONST3) = 2.d0*sin(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! eta
          matval(M_PLCONST4) = 6.d0*cos(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! xi
          matval(M_PLCONST5) = 2.d0*sin(psi)/ ( sqrt(3.d0)*(3.d0+sin(psi)) ) ! etabar
        else                  ! Mohr-Coulomb
          matval(M_PLCONST3) = phi
          matval(M_PLCONST4) = psi
        endif

      case(4)
        fstr_ctrl_get_PLASTICITY = read_user_matl( ctrl, matval(101:200) )

      case default
        stop "Yield function not supported"
    end select

    if( associated(fval) ) deallocate(fval)
    call finalize_table( mttable )
  end function fstr_ctrl_get_PLASTICITY


  !----------------------------------------------------------------------
  !> Read in !CREEP
  integer function fstr_ctrl_get_VISCOPLASTICITY( ctrl, mattype, nlgeom, dict )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    type(DICT_STRUCT), pointer        :: dict

    integer(kind=kint) :: i,j, rcode, depends, ipt, n
    real(kind=kreal),pointer :: fval(:,:)
    character(len=HECMW_NAME_LEN) :: data_fmt
    type( tTable )        :: mattable
    character(len=256) :: s

    fstr_ctrl_get_VISCOPLASTICITY = -1
    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends=1   ! temperature depends only currently
    nlgeom = UPDATELAG   !default value
    if( fstr_ctrl_get_param_ex( ctrl, 'KIRCHHOFF ',  '# ',    0,   'E',   ipt )/= 0) return
    if( ipt/=0 ) nlgeom = TOTALLAG

    ipt=1
    s = 'NORTON,USER '
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',  s, 0, 'P',   ipt    ) /= 0 ) return
    ipt = 1

    ! NORTON
    if( ipt==1 ) then
      n = fstr_ctrl_get_data_line_n( ctrl )
      allocate( fval(3+depends,n) )
      if( depends==0 ) then
        data_fmt = "RRR "
        fstr_ctrl_get_VISCOPLASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
      endif
      if( depends==1 ) then
        data_fmt = "RRRR "
        fstr_ctrl_get_VISCOPLASTICITY = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:), fval(4,:) )
      endif
      if( fstr_ctrl_get_VISCOPLASTICITY ==0 ) then
        call init_table( mattable, depends, 3+depends,n, fval )
        call dict_add_key( dict, MC_NORTON, mattable )
        !  call print_table( mattable, 6 ); pause
      endif
      mattype = NORTON

    else
      stop "ERROR: Material type not supported"

    endif

    call finalize_table( mattable )
    if( associated(fval) ) deallocate(fval)
  end function fstr_ctrl_get_VISCOPLASTICITY

  !----------------------------------------------------------------------
  !> Read in !DENSITY
  integer function fstr_ctrl_get_DENSITY( ctrl, matval )
    integer(kind=kint), intent(in) :: ctrl
    real(kind=kreal),intent(out)   :: matval(:)

    integer(kind=kint) :: i, rcode, depends
    real(kind=kreal),pointer :: fval(:,:)
    character(len=HECMW_NAME_LEN) :: data_fmt

    data_fmt = "R "

    fstr_ctrl_get_DENSITY = -1

    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends = 1 ! we consider temperature dependence only currently

    allocate( fval(1,depends+1) )
    do i=2,1+depends
      data_fmt = data_fmt //"R "
    enddo
    fstr_ctrl_get_DENSITY                                      &
      = fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:) )
    if( fstr_ctrl_get_DENSITY==0 ) matval(M_DENSITY) = fval(1,1)

    if( associated(fval) ) deallocate(fval)

  end function fstr_ctrl_get_DENSITY


  !----------------------------------------------------------------------
  !> Read in !EXPANSION_COEFF
  integer function fstr_ctrl_get_EXPANSION_COEFF( ctrl, matval, dict )
    integer(kind=kint), intent(in) :: ctrl
    real(kind=kreal),intent(out)   :: matval(:)
    type(DICT_STRUCT), pointer     :: dict

    integer(kind=kint) :: i, n, rcode, depends, ipt
    real(kind=kreal),pointer :: fval(:,:)
    type( tTable )           :: mttable
    character(len=HECMW_NAME_LEN) :: data_fmt, ss

    data_fmt = "R "

    fstr_ctrl_get_EXPANSION_COEFF = -1
    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n == 0 ) return               ! fail in reading plastic

    ss = 'ISOTROPIC,ORTHOTROPIC '
    ipt = 1  !default
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',  ss, 0, 'P',   ipt    ) /= 0 ) return

    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends = 1 ! we consider temperature dependence only currently

    if( ipt==1 ) then
      allocate( fval(depends+1, n) )
      do i=2,1+depends
        data_fmt = data_fmt //"R "
      enddo
      if( depends==0 ) then
        fstr_ctrl_get_EXPANSION_COEFF = &
          fstr_ctrl_get_data_array_ex( ctrl, "R ", fval(1,:) )
      else
        fstr_ctrl_get_EXPANSION_COEFF = &
          fstr_ctrl_get_data_array_ex( ctrl, "RR ", fval(1,:), fval(2,:) )
      endif
      if( fstr_ctrl_get_EXPANSION_COEFF==0 ) then
        matval(M_EXAPNSION) = fval(1,1)
        call init_table( mttable,depends, 1+depends, n, fval )
        call dict_add_key( dict, MC_THEMOEXP, mttable )
      endif
    else
      allocate( fval(3+depends,n) )
      do i=2,3+depends
        data_fmt = trim(data_fmt) //"R "
      enddo
      if( depends==0 ) then
        fstr_ctrl_get_EXPANSION_COEFF = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:) )
      elseif( depends==1 ) then
        fstr_ctrl_get_EXPANSION_COEFF = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:), fval(3,:), fval(4,:) )
      endif
      if( fstr_ctrl_get_EXPANSION_COEFF==0 ) then
        call init_table( mttable, depends, 3+depends,n, fval )
        if( fstr_ctrl_get_EXPANSION_COEFF==0 ) call dict_add_key( dict, MC_ORTHOEXP, mttable )
      endif
    endif

    call finalize_table( mttable )
    if( associated(fval) ) deallocate(fval)
  end function fstr_ctrl_get_EXPANSION_COEFF


  integer function read_user_matl( ctrl, matval )
    integer(kind=kint), intent(in)    :: ctrl
    real(kind=kreal),intent(out)      :: matval(:)

    integer(kind=kint) :: n, i, j
    real(kind=kreal)   :: fval(10,10)

    read_user_matl = -1

    n = fstr_ctrl_get_data_line_n( ctrl )
    if( n > 10 ) stop "Num of data lines for user-defined material exceeds 10"
    fval =0.d0
    if( fstr_ctrl_get_data_array_ex( ctrl, 'rrrrrrrrrr ', fval(1,:), fval(2,:), fval(3,:),  &
      fval(4,:), fval(5,:), fval(6,:), fval(7,:), fval(8,:), fval(9,:), fval(10,:) ) /= 0 ) return
    do i=1,10
      do j=1,10
        matval((i-1)*10+j)=fval(j,i)
      enddo
    enddo

    read_user_matl = 0
  end function read_user_matl

  !----------------------------------------------------------------------
  !> Read in !FLUID
  integer function fstr_ctrl_get_FLUID( ctrl, mattype, nlgeom, matval, dict )
    integer(kind=kint), intent(in)    :: ctrl
    integer(kind=kint), intent(inout) :: mattype
    integer(kind=kint), intent(out)   :: nlgeom
    real(kind=kreal),intent(out)      :: matval(:)
    type(DICT_STRUCT), pointer        :: dict

    integer(kind=kint) :: i,j, rcode, depends, ipt, n
    real(kind=kreal),pointer :: fval(:,:)
    character(len=HECMW_NAME_LEN) :: data_fmt
    type( tTable )        :: mattable
    logical            :: isok
    character(len=256) :: s

    fstr_ctrl_get_FLUID = -1
    depends = 0
    rcode = fstr_ctrl_get_param_ex( ctrl, 'DEPENDENCIES  ', '# ',           0,   'I',   depends )
    if( depends>1 ) depends=1   ! temperature depends only currently
    if( depends > 3 ) stop "We cannot read dependencies>3 right now"
    nlgeom = TOTALLAG   !default value

    ipt=1
    s = 'INCOMP_NEWTONIAN '
    if( fstr_ctrl_get_param_ex( ctrl, 'TYPE ',  s, 0, 'P',   ipt    ) /= 0 ) return

    n = fstr_ctrl_get_data_line_n( ctrl )
    ! ISOTROPIC
    if( ipt==1 ) then
      allocate( fval(1+depends,n) )
      if( depends==0 ) then
        data_fmt = "R "
        fstr_ctrl_get_FLUID = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:) )
      endif
      if( depends==1 ) then
        data_fmt = "RR "
        fstr_ctrl_get_FLUID = &
          fstr_ctrl_get_data_array_ex( ctrl, data_fmt, fval(1,:), fval(2,:) )
      endif
      if( fstr_ctrl_get_FLUID ==0 ) then
        matval(M_VISCOCITY) = fval(1,1)
        call init_table( mattable, depends, 1+depends,n, fval )
        call dict_add_key( dict, MC_INCOMP_NEWTONIAN, mattable )
        !         call print_table( mattable, 6 ); pause
      endif
      mattype = INCOMP_NEWTONIAN

    else
      stop "ERROR: Material type not supported"

    endif

    call finalize_table( mattable )
    if( associated(fval) ) deallocate(fval)
  end function fstr_ctrl_get_FLUID

end module fstr_ctrl_material
