!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module provide functions for elastoplastic calculation
module m_ElastoPlastic
  use hecmw_util
  use mMaterial
  use m_ElasticLinear
  use mUYield

  implicit none

  private
  public :: calElastoPlasticMatrix
  public :: BackwardEuler
  public :: updateEPState

  real(kind=kreal), private, parameter :: Id(6,6) = reshape( &
    & (/  2.d0/3.d0, -1.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &    -1.d0/3.d0,  2.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &    -1.d0/3.d0, -1.d0/3.d0,  2.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0, 0.5d0,  0.d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0,  0.d0, 0.5d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0,  0.d0,  0.d0, 0.5d0/), &
    & (/6, 6/))

contains

  !> This subroutine calculates elastoplastic constitutive relation
  subroutine calElastoPlasticMatrix( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag )
    type( tMaterial ), intent(in) :: matl      !< material properties
    integer, intent(in)           :: sectType  !< not used currently
    real(kind=kreal), intent(in)  :: stress(6) !< stress
    real(kind=kreal), intent(in)  :: extval(:) !< plastic strain, back stress
    real(kind=kreal), intent(in)  :: plstrain  !< plastic strain
    integer, intent(in)           :: istat     !< plastic state
    real(kind=kreal), intent(out) :: D(:,:)    !< constitutive relation
    real(kind=kreal), intent(in)  :: temperature   !> temperature
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: ytype,hdflag_in

    hdflag_in = 0
    if( present(hdflag) ) hdflag_in = hdflag

    ytype = getYieldFunction( matl%mtype )
    select case (ytype)
    case (0)
      call calElastoPlasticMatrix_VM( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag_in )
    case (1)
      call calElastoPlasticMatrix_MC( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag_in )
    case (2)
      call calElastoPlasticMatrix_DP( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag_in )
    case (3)
      call uElastoPlasticMatrix( matl%variables, stress, istat, extval, plstrain, D, temperature, hdflag_in )
    end select
  end subroutine calElastoPlasticMatrix

  !> This subroutine calculates elastoplastic constitutive relation
  subroutine calElastoPlasticMatrix_VM( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag )
    type( tMaterial ), intent(in) :: matl      !< material properties
    integer, intent(in)           :: sectType  !< not used currently
    real(kind=kreal), intent(in)  :: stress(6) !< stress
    real(kind=kreal), intent(in)  :: extval(:) !< plastic strain, back stress
    real(kind=kreal), intent(in)  :: plstrain  !< plastic strain
    integer, intent(in)           :: istat     !< plastic state
    real(kind=kreal), intent(out) :: D(:,:)    !< constitutive relation
    real(kind=kreal), intent(in)  :: temperature   !> temperature
    integer(kind=kint), intent(in) :: hdflag  !> return only hyd and dev term if specified

    integer :: i,j
    logical :: kinematic
    real(kind=kreal) :: dum, a(6), G, dlambda
    real(kind=kreal) :: C1,C2,C3, back(6)
    real(kind=kreal) :: J1,J2, harden, khard, devia(6)

    if( sectType /=D3 ) stop "Elastoplastic calculation support only Solid element currently"

    call calElasticMatrix( matl, sectTYPE, D, temperature, hdflag=hdflag )
    if( istat == 0 ) return   ! elastic state
    if( hdflag == 2 ) return

    harden = calHardenCoeff( matl, extval(1), temperature )

    kinematic = isKinematicHarden( matl%mtype )
    khard = 0.d0
    if( kinematic ) then
      back(1:6) = extval(2:7)
      khard = calKinematicHarden( matl, extval(1) )
    endif

    J1 = (stress(1)+stress(2)+stress(3))
    devia(1:3) = stress(1:3)-J1/3.d0
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    a(1:6) = devia(1:6)/sqrt(2.d0*J2)
    G = D(4,4)
    dlambda = extval(1)-plstrain
    C3 = sqrt(3.d0*J2)+3.d0*G*dlambda !trial mises stress
    C1 = 6.d0*dlambda*G*G/C3
    dum = 3.d0*G+khard+harden
    C2 = 6.d0*G*G*(dlambda/C3-1.d0/dum)

    do i=1,6
      do j=1,6
        D(i,j) = D(i,j) - C1*Id(i,j) + C2*a(i)*a(j)
      enddo
    enddo
  end subroutine calElastoPlasticMatrix_VM

  !> This subroutine calculates elastoplastic constitutive relation
  subroutine calElastoPlasticMatrix_MC( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag )
    type( tMaterial ), intent(in) :: matl      !< material properties
    integer, intent(in)           :: sectType  !< not used currently
    real(kind=kreal), intent(in)  :: stress(6) !< stress
    real(kind=kreal), intent(in)  :: extval(:) !< plastic strain, back stress
    real(kind=kreal), intent(in)  :: plstrain  !< plastic strain
    integer, intent(in)           :: istat     !< plastic state
    real(kind=kreal), intent(out) :: D(:,:)    !< constitutive relation
    real(kind=kreal), intent(in)  :: temperature   !> temperature
    integer(kind=kint), intent(in) :: hdflag  !> return only hyd and dev term if specified

    integer :: i,j
    logical :: kinematic
    real(kind=kreal) :: dum, dj1(6), dj2(6), dj3(6), a(6)
    real(kind=kreal) :: C1,C2,C3, back(6)
    real(kind=kreal) :: J1,J2,J3, phi, sin3theta, theta, harden, khard, da(6), devia(6)

    if( sectType /=D3 ) stop "Elastoplastic calculation support only Solid element currently"

    call calElasticMatrix( matl, sectTYPE, D, temperature, hdflag=hdflag )
    if( istat == 0 ) return   ! elastic state
    if( hdflag == 2 ) return

    harden = calHardenCoeff( matl, extval(1), temperature )

    kinematic = isKinematicHarden( matl%mtype )
    khard = 0.d0
    if( kinematic ) then
      back(1:6) = extval(2:7)
      khard = calKinematicHarden( matl, extval(1) )
    endif

    J1 = (stress(1)+stress(2)+stress(3))
    devia(1:3) = stress(1:3)-J1/3.d0
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    !derivative of J2
    dj2(1:3) = devia(1:3)
    dj2(4:6) = 2.d0*devia(4:6)
    dj2 = dj2/( 2.d0*sqrt(j2) )

    phi = matl%variables(M_PLCONST3)
    J3 = devia(1)*devia(2)*devia(3)                    &
        +2.d0* devia(4)*devia(5)*devia(6)                     &
        -devia(6)*devia(2)*devia(6)                           &
        -devia(4)*devia(4)*devia(3)                           &
        -devia(1)*devia(5)*devia(5)
    sin3theta = -3.d0*sqrt(3.d0)*J3/( 2.d0*(J2**1.5d0) )
    if( abs( abs(sin3theta)-1.d0 ) <1.d-8 ) then
      C1 = 0.d0
      C2 = sqrt(3.d0)
      C3 = 0.d0
    else
      if( abs(sin3theta) >1.d0 ) stop "Math Error in Mohr-Coulomb calculation"
      theta = asin( sin3theta )/3.d0
      C2 = cos(theta)*( 1.d0*tan(theta)*tan(3.d0*theta) + sin(phi)* &
          ( tan(3.d0*theta)-tan(theta )/sqrt(3.d0) ) )
      C1 = sin(phi)/3.d0
      C3 = sqrt(3.d0)*sin(theta)+cos(theta)*sin(phi)/(2.d0*J2*cos(3.d0*theta))
    endif
    ! deirivative of j1
    dj1(1:3) = 1.d0
    dj1(4:6) = 0.d0
    ! deirivative of j3
    dj3(1) = devia(2)*devia(3)-devia(5)*devia(5)+J2/3.d0
    dj3(2) = devia(1)*devia(3)-devia(6)*devia(6)+J2/3.d0
    dj3(3) = devia(1)*devia(2)-devia(4)*devia(4)+J2/3.d0
    dj3(4) = 2.d0*(devia(5)*devia(6)-devia(3)*devia(4))
    dj3(5) = 2.d0*(devia(4)*devia(6)-devia(1)*devia(5))
    dj3(6) = 2.d0*(devia(4)*devia(5)-devia(2)*devia(6))
    a(:) = C1*dj1 + C2*dj2 + C3*dj3

    da = matmul( D, a )
    dum = harden + khard+ dot_product( da, a )
    do i=1,6
      do j=1,6
        D(i,j) = D(i,j) - da(i)*da(j)/dum
      enddo
    enddo
  end subroutine calElastoPlasticMatrix_MC

  !> This subroutine calculates elastoplastic constitutive relation
  subroutine calElastoPlasticMatrix_DP( matl, sectType, stress, istat, extval, plstrain, D, temperature, hdflag )
    type( tMaterial ), intent(in) :: matl      !< material properties
    integer, intent(in)           :: sectType  !< not used currently
    real(kind=kreal), intent(in)  :: stress(6) !< stress
    real(kind=kreal), intent(in)  :: extval(:) !< plastic strain, back stress
    real(kind=kreal), intent(in)  :: plstrain  !< plastic strain
    integer, intent(in)           :: istat     !< plastic state
    real(kind=kreal), intent(out) :: D(:,:)    !< constitutive relation
    real(kind=kreal), intent(in)  :: temperature   !> temperature
    integer(kind=kint), intent(in) :: hdflag  !> return only hyd and dev term if specified

    integer :: i,j
    logical :: kinematic
    real(kind=kreal) :: dum, dj1(6), dj2(6), a(6)
    real(kind=kreal) :: back(6)
    real(kind=kreal) :: J1,J2, eta, harden, khard, da(6), devia(6)

    if( sectType /=D3 ) stop "Elastoplastic calculation support only Solid element currently"

    call calElasticMatrix( matl, sectTYPE, D, temperature, hdflag=hdflag )
    if( istat == 0 ) return   ! elastic state
    if( hdflag == 2 ) return

    harden = calHardenCoeff( matl, extval(1), temperature )

    kinematic = isKinematicHarden( matl%mtype )
    khard = 0.d0
    if( kinematic ) then
      back(1:6) = extval(2:7)
      khard = calKinematicHarden( matl, extval(1) )
    endif

    J1 = (stress(1)+stress(2)+stress(3))
    devia(1:3) = stress(1:3)-J1/3.d0
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    !derivative of J2
    dj2(1:3) = devia(1:3)
    dj2(4:6) = 2.d0*devia(4:6)
    dj2 = dj2/( 2.d0*sqrt(j2) )

    eta = matl%variables(M_PLCONST3)
    ! deirivative of j1
    dj1(1:3) = 1.d0
    dj1(4:6) = 0.d0
    a(:) = eta*dj1(:) + dj2(:)

    da = matmul( D, a )
    dum = harden + khard+ dot_product( da, a )
    do i=1,6
      do j=1,6
        D(i,j) = D(i,j) - da(i)*da(j)/dum
      enddo
    enddo
  end subroutine calElastoPlasticMatrix_DP

  !> This function calculates hardening coefficient
  real(kind=kreal) function calHardenCoeff( matl, pstrain, temp )
    type( tMaterial ), intent(in)          :: matl    !< material property
    real(kind=kreal), intent(in)           :: pstrain !< plastic strain
    real(kind=kreal), intent(in)           :: temp !< temperature

    integer :: htype
    logical :: ierr
    real(kind=kreal) :: s0, s1,s2, ef, ina(2)

    calHardenCoeff = -1.d0
    htype = getHardenType( matl%mtype )
    select case (htype)
      case (0)  ! Linear hardening
        calHardenCoeff = matl%variables(M_PLCONST2)
      case (1)  ! Multilinear approximation
        ina(1) = temp;  ina(2)=pstrain
        call fetch_TableGrad( MC_YIELD, ina, matl%dict, calHardenCoeff, ierr )
      case (2)  ! Swift
        s0= matl%variables(M_PLCONST1)
        s1= matl%variables(M_PLCONST2)
        s2= matl%variables(M_PLCONST3)
        calHardenCoeff = s1*s2*( s0+pstrain )**(s2-1)
      case (3)  ! Ramberg-Osgood
        s0= matl%variables(M_PLCONST1)
        s1= matl%variables(M_PLCONST2)
        s2= matl%variables(M_PLCONST3)
        ef = calCurrYield( matl, pstrain, temp )
        calHardenCoeff = s1*(ef/s1)**(1.d0-s2) /(s0*s2)
      case(4)   ! Prager
        calHardenCoeff = 0.d0
      case(5)   ! Prager+linear
        calHardenCoeff = matl%variables(M_PLCONST2)
    end select
  end function

  !> This function calculates kinematic hardening coefficient
  real(kind=kreal) function calKinematicHarden( matl, pstrain )
    type( tMaterial ), intent(in) :: matl    !< material property
    real(kind=kreal), intent(in)  :: pstrain !< plastic strain

    integer :: htype
    htype = getHardenType( matl%mtype )
    select case (htype)
      case(4, 5)   ! Prager
        calKinematicHarden = matl%variables(M_PLCONST3)
      case default
        calKinematicHarden = 0.d0
    end select
  end function

  !> This function calculates state of kinematic hardening
  real(kind=kreal) function calCurrKinematic( matl, pstrain )
    type( tMaterial ), intent(in) :: matl    !< material property
    real(kind=kreal), intent(in)  :: pstrain !< plastic strain

    integer :: htype
    htype = getHardenType( matl%mtype )
    select case (htype)
      case(4, 5)   ! Prager
        calCurrKinematic = matl%variables(M_PLCONST3)*pstrain
      case default
        calCurrKinematic = 0.d0
    end select
  end function

  !> This function calculates current yield stress
  real(kind=kreal) function calCurrYield( matl, pstrain, temp )
    type( tMaterial ), intent(in) :: matl    !< material property
    real(kind=kreal), intent(in)  :: pstrain !< plastic strain
    real(kind=kreal), intent(in)  :: temp  !< temperature

    integer :: htype
    real(kind=kreal) :: s0, s1,s2, ina(2), outa(1)
    logical :: ierr
    calCurrYield = -1.d0
    htype = getHardenType( matl%mtype )

    select case (htype)
      case (0, 5)  ! Linear hardening, Linear+Parger hardening
        calCurrYield = matl%variables(M_PLCONST1)+matl%variables(M_PLCONST2)*pstrain
      case (1)  ! Multilinear approximation
        ina(1) = temp;  ina(2)=pstrain
        call fetch_TableData(MC_YIELD, matl%dict, outa, ierr, ina)
        if( ierr ) stop "Fail to get yield stress!"
        calCurrYield = outa(1)
      case (2)  ! Swift
        s0= matl%variables(M_PLCONST1)
        s1= matl%variables(M_PLCONST2)
        s2= matl%variables(M_PLCONST3)
        calCurrYield = s1*( s0+pstrain )**s2
      case (3)  ! Ramberg-Osgood
        s0= matl%variables(M_PLCONST1)
        s1= matl%variables(M_PLCONST2)
        s2= matl%variables(M_PLCONST3)
        if( pstrain<=s0 ) then
          calCurrYield = s1
        else
          calCurrYield = s1*( pstrain/s0 )**(1.d0/s2)
        endif
      case (4)  ! Parger hardening
        calCurrYield = matl%variables(M_PLCONST1)
    end select
  end function

  !> This subroutine does backward-Euler return calculation
  subroutine BackwardEuler( matl, stress, plstrain, istat, fstat, temp, hdflag )
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(in)     :: temp  !< temperature
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: ytype, hdflag_in

    hdflag_in = 0
    if( present(hdflag) ) hdflag_in = hdflag

    ytype = getYieldFunction( matl%mtype )
    select case (ytype)
    case (0)
      call BackwardEuler_VM( matl, stress, plstrain, istat, fstat, temp, hdflag_in )
    case (1)
      call BackwardEuler_MC( matl, stress, plstrain, istat, fstat, temp, hdflag_in )
    case (2)
      call BackwardEuler_DP( matl, stress, plstrain, istat, fstat, temp, hdflag_in )
    case (3)
      call uBackwardEuler( matl%variables, stress, plstrain, istat, fstat, temp, hdflag_in )
    end select
  end subroutine BackwardEuler

  !> This subroutine does backward-Euler return calculation for von Mises
  subroutine BackwardEuler_VM( matl, stress, plstrain, istat, fstat, temp, hdflag )
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(in)     :: temp  !< temperature
    integer(kind=kint), intent(in)   :: hdflag  !> return only hyd and dev term if specified

    real(kind=kreal), parameter :: tol =1.d-3
    integer, parameter          :: MAXITER = 5
    real(kind=kreal) :: dlambda, f
    integer :: i
    real(kind=kreal) :: youngs, poisson, pstrain, ina(1), ee(2)
    real(kind=kreal) :: J1, J2, H, KH, KK, dd, eqvs, yd, G, K, devia(6)
    logical          :: kinematic, ierr
    real(kind=kreal) :: betan, back(6)

    kinematic = isKinematicHarden( matl%mtype )
    if( kinematic ) back(1:6) = fstat(8:13)

    J1 = (stress(1)+stress(2)+stress(3))
    devia(1:3) = stress(1:3)-J1/3.d0
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )
    eqvs = sqrt( 3.d0*J2 )
    yd = calCurrYield( matl, plstrain, temp )
    f = eqvs - yd

    if( abs(f)<tol ) then  ! yielded
      istat = 1
      return
    elseif( f<0.d0 ) then   ! not yielded or unloading
      istat =0
      return
    endif
    if( hdflag == 2 ) return

    istat = 1           ! yielded
    KH = 0.d0; KK=0.d0; betan=0.d0; back(:)=0.d0
    if( kinematic ) betan = calCurrKinematic( matl, plstrain )

    ina(1) = temp
    call fetch_TableData(MC_ISOELASTIC, matl%dict, ee, ierr, ina)
    if( ierr ) then
      stop " fail to fetch young's modulus in elastoplastic calculation"
    else
      youngs = ee(1)
      poisson = ee(2)
    endif
    if( youngs==0.d0 ) stop "YOUNG's ratio==0"
    G = youngs/ ( 2.d0*(1.d0+poisson) )
    K = youngs/ ( 3.d0*(1.d0-2.d0*poisson) )
    dlambda = 0.d0
    pstrain = plstrain

    do i=1,MAXITER
      H= calHardenCoeff( matl, pstrain, temp )
      if( kinematic ) then
        KH = calKinematicHarden( matl, pstrain )
      endif
      dd= 3.d0*G+H+KH
      dlambda = dlambda+f/dd
      if( dlambda<0.d0 ) then
        dlambda = 0.d0
        pstrain = plstrain
        istat=0; exit
      endif
      pstrain = plstrain+dlambda
      yd = calCurrYield( matl, pstrain, temp )
      if( kinematic ) then
        KK = calCurrKinematic( matl, pstrain )
      endif
      f = eqvs-3.d0*G*dlambda-yd -(KK-betan)
      if( abs(f)<tol*tol ) exit
    enddo
    if( kinematic ) then
      KK = calCurrKinematic( matl, pstrain )
      fstat(2:7) = back(:)+(KK-betan)*devia(:)/eqvs
    endif
    devia(:) = (1.d0-3.d0*dlambda*G/eqvs)*devia(:)
    stress(1:3) = devia(1:3)+J1/3.d0
    stress(4:6) = devia(4:6)
    stress(:)= stress(:)+back(:)

    fstat(1) = pstrain
  end subroutine BackwardEuler_VM

  !> This subroutine does backward-Euler return calculation for Mohr-Coulomb
  subroutine BackwardEuler_MC( matl, stress, plstrain, istat, fstat, temp, hdflag )
    use m_utilities, only : eigen3
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(in)     :: temp  !< temperature
    integer(kind=kint), intent(in)   :: hdflag  !> return only hyd and dev term if specified

    real(kind=kreal), parameter :: tol =1.d-3
    integer, parameter          :: MAXITER = 5
    real(kind=kreal) :: dlambda, f, mat(3,3)
    integer :: i, maxp(1), minp(1), mm
    real(kind=kreal) :: youngs, poisson, pstrain, ina(1), ee(2)
    real(kind=kreal) :: J1,J2,J3, H, KH, KK, dd, eqvs, cohe, G, K, devia(6)
    real(kind=kreal) :: prnstre(3), prnprj(3,3), tstre(3,3)
    real(kind=kreal) :: sin3theta, theta, phi, trialprn(3)
    logical          :: kinematic, ierr
    real(kind=kreal) :: betan, back(6)

    phi = matl%variables(M_PLCONST3)

    kinematic = isKinematicHarden( matl%mtype )
    if( kinematic ) back(1:6) = fstat(8:13)

    J1 = (stress(1)+stress(2)+stress(3))
    devia(1:3) = stress(1:3)-J1/3.d0
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    J3 = devia(1)*devia(2)*devia(3)                    &
        +2.d0* devia(4)*devia(5)*devia(6)                     &
        -devia(2)*devia(6)*devia(6)                           &
        -devia(3)*devia(4)*devia(4)                           &
        -devia(1)*devia(5)*devia(5)
    sin3theta = -3.d0*sqrt(3.d0)*J3/( 2.d0*(J2**1.5d0) )
    if( abs( abs(sin3theta)-1.d0 ) <1.d-8 ) sin3theta=sign(1.d0, sin3theta)
    if( abs(sin3theta) >1.d0 ) stop "Math Error in Mohr-Coulomb calculation"
    theta = asin( sin3theta )/3.d0
    eqvs = (cos(theta)-sin(theta)*sin(phi)/sqrt(3.d0))*sqrt(J2)  &
        +J1*sin(phi)/3.d0
    cohe = calCurrYield( matl, plstrain, temp )
    f = eqvs - cohe*cos(phi)

    if( abs(f)<tol ) then  ! yielded
      istat = 1
      return
    elseif( f<0.d0 ) then   ! not yielded or unloading
      istat =0
      return
    endif
    if( hdflag == 2 ) return

    istat = 1           ! yielded
    KH = 0.d0; KK=0.d0; betan=0.d0; back(:)=0.d0

    if( kinematic ) betan = calCurrKinematic( matl, plstrain )

    ina(1) = temp
    call fetch_TableData(MC_ISOELASTIC, matl%dict, ee, ierr, ina)
    if( ierr ) then
      stop " fail to fetch young's modulus in elastoplastic calculation"
    else
      youngs = ee(1)
      poisson = ee(2)
    endif
    if( youngs==0.d0 ) stop "YOUNG's ratio==0"
    G = youngs/ ( 2.d0*(1.d0+poisson) )
    K = youngs/ ( 3.d0*(1.d0-2.d0*poisson) )
    dlambda = 0.d0
    pstrain = plstrain

    do mm=1,6
      if( abs(stress(mm))<1.d-10 ) stress(mm)=0.d0
    enddo
    call eigen3( stress, prnstre, prnprj )
    trialprn = prnstre
    maxp = maxloc( prnstre )
    minp = minloc( prnstre )
    mm = 1
    if( maxp(1)==1 .or. minp(1)==1 ) mm =2
    if( maxp(1)==2 .or. minp(1)==2 ) mm =3
    do i=1,MAXITER
      H= calHardenCoeff( matl, pstrain, temp )
      dd= 4.d0*G*( 1.d0+sin(phi)*sin(theta)/3.d0 )+4.d0*K         &
          *sin(phi)*sin(theta)+4.d0*H*cos(phi)*cos(phi)
      dlambda = dlambda+f/dd
      if( 2.d0*dlambda*cos(phi)<0.d0 ) then
        if( cos(phi)==0.d0 ) stop "Math error in return mapping"
        dlambda = 0.d0
        pstrain = plstrain
        istat=0; exit
      endif
      pstrain = plstrain + 2.d0*dlambda*cos(phi)
      cohe = calCurrYield( matl, pstrain, temp )
      f = prnstre(maxp(1))-prnstre(minp(1))+                     &
          (prnstre(maxp(1))+prnstre(minp(1)))*sin(phi)-            &
          (4.d0*G*(1.d0+sin(phi)*sin(theta)/3.d0)+4.d0*K*sin(phi)   &
          *sin(theta))*dlambda-2.d0*cohe*cos(phi)
      if( abs(f)<tol ) exit
    enddo
    prnstre(maxp(1)) = prnstre(maxp(1))-(2.d0*G*(1.d0+sin(phi)/3.d0)  &
        + 2.d0*K*sin(phi) )*dlambda
    prnstre(minp(1)) = prnstre(minp(1))+(2.d0*G*(1.d0-sin(phi)/3.d0)  &
        - 2.d0*K*sin(phi) )*dlambda
    prnstre(mm) = prnstre(mm)+(4.d0*G/3.d0-2.d0*K)*sin(phi)*dlambda

    tstre(:,:) = 0.d0
    tstre(1,1)= prnstre(1); tstre(2,2)=prnstre(2); tstre(3,3)=prnstre(3)
    mat= matmul( prnprj, tstre )
    mat= matmul( mat, transpose(prnprj) )
    stress(1) = mat(1,1)
    stress(2) = mat(2,2)
    stress(3) = mat(3,3)
    stress(4) = mat(1,2)
    stress(5) = mat(2,3)
    stress(6) = mat(3,1)

    fstat(1) = pstrain
  end subroutine BackwardEuler_MC

  !> This subroutine does backward-Euler return calculation for Drucker-Prager
  subroutine BackwardEuler_DP( matl, stress, plstrain, istat, fstat, temp, hdflag )
    use m_utilities, only : eigen3
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(in)     :: temp  !< temperature
    integer(kind=kint), intent(in)   :: hdflag  !> return only hyd and dev term if specified

    real(kind=kreal), parameter :: tol =1.d-3
    integer, parameter          :: MAXITER = 5
    real(kind=kreal) :: dlambda, f
    integer :: i
    real(kind=kreal) :: youngs, poisson, pstrain, xi, ina(1), ee(2)
    real(kind=kreal) :: J1,J2,H, KH, KK, dd, eqvs, cohe, G, K, devia(6), eta, p
    logical          :: kinematic, ierr
    real(kind=kreal) :: betan, back(6)

    eta = matl%variables(M_PLCONST3)
    xi = matl%variables(M_PLCONST4)

    kinematic = isKinematicHarden( matl%mtype )
    if( kinematic ) back(1:6) = fstat(8:13)

    J1 = (stress(1)+stress(2)+stress(3))
    p = J1/3.d0
    devia(1:3) = stress(1:3)-p
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    eqvs = sqrt(J2)
    cohe = calCurrYield( matl, plstrain, temp )
    f = eqvs + eta*J1 - cohe*xi

    if( abs(f)<tol ) then  ! yielded
      istat = 1
      return
    elseif( f<0.d0 ) then   ! not yielded or unloading
      istat =0
      return
    endif
    if( hdflag == 2 ) return

    istat = 1           ! yielded
    KH = 0.d0; KK=0.d0; betan=0.d0; back(:)=0.d0

    if( kinematic ) betan = calCurrKinematic( matl, plstrain )

    ina(1) = temp
    call fetch_TableData(MC_ISOELASTIC, matl%dict, ee, ierr, ina)
    if( ierr ) then
      stop " fail to fetch young's modulus in elastoplastic calculation"
    else
      youngs = ee(1)
      poisson = ee(2)
    endif
    if( youngs==0.d0 ) stop "YOUNG's ratio==0"
    G = youngs/ ( 2.d0*(1.d0+poisson) )
    K = youngs/ ( 3.d0*(1.d0-2.d0*poisson) )
    dlambda = 0.d0
    pstrain = plstrain

    do i=1,MAXITER
      H= calHardenCoeff( matl, pstrain, temp )
      dd= G+K*eta*eta+H*xi*xi
      dlambda = dlambda+f/dd
      if( xi*dlambda<0.d0 ) then
        if( xi==0.d0 ) stop "Math error in return mapping"
        dlambda = 0.d0
        pstrain = plstrain
        istat=0; exit
      endif
      pstrain = plstrain+xi*dlambda
      cohe = calCurrYield( matl, pstrain, temp  )
      f = eqvs-G*dlambda+eta*(p-K*eta*dlambda)- xi*cohe
      if( abs(f)<tol*tol ) exit
    enddo
    devia(:) = (1.d0-G*dlambda/eqvs)*devia(:)
    p = p-K*eta*dlambda
    stress(1:3) = devia(1:3)+p
    stress(4:6) = devia(4:6)

    fstat(1) = pstrain
  end subroutine BackwardEuler_DP

  !> Clear elatoplastic state
  subroutine updateEPState( gauss )
    use mMechGauss
    type(tGaussStatus), intent(inout) :: gauss  ! status of curr gauss point
    gauss%plstrain= gauss%fstatus(1)
    if(isKinematicHarden(gauss%pMaterial%mtype)) then
      gauss%fstatus(8:13) =gauss%fstatus(2:7)
    endif
  end subroutine

end module m_ElastoPlastic
