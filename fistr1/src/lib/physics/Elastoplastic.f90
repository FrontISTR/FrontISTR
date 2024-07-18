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

  real(kind=kreal), parameter :: Id(6,6) = reshape( &
    & (/  2.d0/3.d0, -1.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &    -1.d0/3.d0,  2.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &    -1.d0/3.d0, -1.d0/3.d0,  2.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0, 0.5d0,  0.d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0,  0.d0, 0.5d0,  0.d0,   &
    &          0.d0,       0.d0,       0.d0,  0.d0,  0.d0, 0.5d0/), &
    & (/6, 6/))
  real(kind=kreal), parameter :: I2(6) = (/ 1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0 /)

  integer, parameter :: VM_ELASTIC = 0
  integer, parameter :: VM_PLASTIC = 1

  integer, parameter :: MC_ELASTIC       = 0
  integer, parameter :: MC_PLASTIC_SURF  = 1
  integer, parameter :: MC_PLASTIC_RIGHT = 2
  integer, parameter :: MC_PLASTIC_LEFT  = 3
  integer, parameter :: MC_PLASTIC_APEX  = 4

  integer, parameter :: DP_ELASTIC      = 0
  integer, parameter :: DP_PLASTIC_SURF = 1
  integer, parameter :: DP_PLASTIC_APEX = 2

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
    if( istat == VM_ELASTIC ) return
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
    use m_utilities, only : eigen3,deriv_general_iso_tensor_func_3d
    type( tMaterial ), intent(in) :: matl      !< material properties
    integer, intent(in)           :: sectType  !< not used currently
    real(kind=kreal), intent(in)  :: stress(6) !< stress
    real(kind=kreal), intent(in)  :: extval(:) !< plastic strain, back stress
    real(kind=kreal), intent(in)  :: plstrain  !< plastic strain
    integer, intent(in)           :: istat     !< plastic state
    real(kind=kreal), intent(out) :: D(:,:)    !< constitutive relation
    real(kind=kreal), intent(in)  :: temperature   !> temperature
    integer(kind=kint), intent(in) :: hdflag  !> return only hyd and dev term if specified

    real(kind=kreal) :: G, K, harden, r2G, r2Gd3, r4Gd3, r2K, youngs, poisson
    real(kind=kreal) :: phi, psi, cosphi, sinphi, cotphi, sinpsi, sphsps, r2cosphi, r4cos2phi
    real(kind=kreal) :: prnstre(3), prnprj(3,3), tstre(3,3), prnstra(3)
    integer(kind=kint) :: m1, m2, m3
    real(kind=kreal) :: C1,C2,C3, CA1, CA2, CA3, CAm, CAp, CD1, CD2, CD3, Cdiag, Coffd
    real(kind=kreal) :: CK1, CK2, CK3
    real(kind=kreal) :: dum, da, db, dc, dd, detinv
    real(kind=kreal) :: dpsdpe(3,3)

    if( sectType /=D3 ) stop "Elastoplastic calculation support only Solid element currently"

    call calElasticMatrix( matl, sectTYPE, D, temperature, hdflag=hdflag )
    if( istat == MC_ELASTIC ) return
    if( hdflag == 2 ) return

    harden = calHardenCoeff( matl, extval(1), temperature )
    G = D(4,4)
    K = D(1,1)-(4.d0/3.d0)*G
    r2G = 2.d0*G
    r2K = 2.d0*K
    r2Gd3 = r2G/3.d0
    r4Gd3 = 2.d0*r2Gd3
    youngs = 9.d0*K*G/(3.d0*K+G)
    poisson = (3.d0*K-r2G)/(6.d0*K+r2G)

    phi = matl%variables(M_PLCONST3)
    psi = matl%variables(M_PLCONST4)
    sinphi = sin(phi)
    cosphi = cos(phi)
    sinpsi = sin(psi)
    sphsps = sinphi*sinpsi
    r2cosphi = 2.d0*cosphi
    r4cos2phi = r2cosphi*r2cosphi

    call eigen3( stress, prnstre, prnprj )
    m1 = maxloc( prnstre, 1 )
    m3 = minloc( prnstre, 1 )
    if( m1 == m3 ) then
      m1 = 1; m2 = 2; m3 = 3
    else
      m2 = 6 - (m1 + m3)
    endif

    C1 = 4.d0*(G*(1.d0+sphsps/3.d0)+K*sphsps)
    if( istat==MC_PLASTIC_SURF ) then
      dd= C1 + r4cos2phi*harden
      CD1 = (r2G*(1.d0+sinpsi/3.d0) + r2K*sinpsi)/dd
      CD2 = (r4Gd3-r2K)*sinpsi/dd
      CD3 = (r2G*(1.d0-sinpsi/3.d0) - r2K*sinpsi)/dd
      CAp = 1.d0+sinphi/3.d0
      CAm = 1.d0-sinphi/3.d0
      CK1 = 1.d0-2.d0*CD1*sinphi
      CK2 = 1.d0+2.d0*CD2*sinphi
      CK3 = 1.d0+2.d0*CD3*sinphi
      dpsdpe(m1,m1) = r2G*( 2.d0/3.d0-CD1*CAp)+K*CK1
      dpsdpe(m1,m2) = (K-r2Gd3)*CK1
      dpsdpe(m1,m3) = r2G*(-1.d0/3.d0+CD1*CAm)+K*CK1
      dpsdpe(m2,m1) = r2G*(-1.d0/3.d0+CD2*CAp)+K*CK2
      dpsdpe(m2,m2) = r4Gd3*( 1.d0-CD2*sinphi)+K*CK2
      dpsdpe(m2,m3) = r2G*(-1.d0/3.d0-CD2*CAm)+K*CK2
      dpsdpe(m3,m1) = r2G*(-1.d0/3.d0+CD3*CAp)+K*CK3
      dpsdpe(m3,m2) = (K-r2Gd3)*CK3
      dpsdpe(m3,m3) = r2G*( 2.d0/3.d0-CD3*CAm)+K*CK3
    else if( istat==MC_PLASTIC_APEX ) then
      cotphi = cosphi/sinphi
      dpsdpe(:,:) = K*(1.d0-(K/(K+harden*cotphi*cosphi/sinpsi)))
    else ! EDGE
      if( istat==MC_PLASTIC_RIGHT ) then
        C2 = r2G*(1.d0+sinphi+sinpsi-sphsps/3.d0) + 4.d0*K*sphsps
      else if( istat==MC_PLASTIC_LEFT ) then
        C2 = r2G*(1.d0-sinphi-sinpsi-sphsps/3.d0) + 4.d0*K*sphsps
      endif
      dum = r4cos2phi*harden
      da = C1 + dum
      db = C2 + dum
      dc = db
      dd = da
      detinv = 1.d0/(da*dd-db*dc)
      CA1 = r2G*(1.d0+sinphi/3.d0)+r2K*sinpsi
      CA2 = (r4Gd3-r2K)*sinpsi
      CA3 = r2G*(1.d0-sinpsi/3.d0)-r2K*sinpsi
      Cdiag = K+r4Gd3
      Coffd = K-r2Gd3
      if( istat==MC_PLASTIC_RIGHT ) then
        dpsdpe(m1,m1) = Cdiag+CA1*(db-dd-da+dc)*(r2G+(r2K+r2Gd3)*sinphi)*detinv
        dpsdpe(m1,m2) = Coffd+CA1*(r2G*(da-db)+((db-dd-da+dc)*(r2K+r2Gd3)+(dd-dc)*r2G)*sinphi)*detinv
        dpsdpe(m1,m3) = Coffd+CA1*(r2G*(dd-dc)+((db-dd-da+dc)*(r2K+r2Gd3)+(da-db)*r2G)*sinphi)*detinv
        dpsdpe(m2,m1) = Coffd+(CA2*(dd-db)+CA3*(da-dc))*(r2G+(r2K+r2Gd3)*sinphi)*detinv
        dpsdpe(m2,m2) = Cdiag+(CA2*((r2K*(dd-db)-(db*r2Gd3+dd*r4Gd3))*sinphi+db*r2G) &
            &                 +CA3*((r2K*(da-dc)+(da*r2Gd3+dc*r4Gd3))*sinphi-da*r2G))*detinv
        dpsdpe(m2,m3) = Coffd+(CA2*((r2K*(dd-db)+(db*r4Gd3+dd*r2Gd3))*sinphi-dd*r2G) &
            &                 +CA3*((r2K*(da-dc)-(da*r4Gd3+dc*r2Gd3))*sinphi+dc*r2G))*detinv
        dpsdpe(m3,m1) = Coffd+(CA2*(da-dc)+CA3*(dd-db))*(r2G+(r2K+r2Gd3)*sinphi)*detinv
        dpsdpe(m3,m2) = Coffd+(CA2*((r2K*(da-dc)+(da*r2Gd3+dc*r4Gd3))*sinphi-da*r2G) &
            &                 +CA3*((r2K*(dd-db)-(db*r2Gd3+dd*r4Gd3))*sinphi+db*r2G))*detinv
        dpsdpe(m3,m3) = Cdiag+(CA2*((r2K*(da-dc)-(da*r4Gd3+dc*r2Gd3))*sinphi+dc*r2G) &
            &                 +CA3*((r2K*(dd-db)+(db*r4Gd3+dd*r2Gd3))*sinphi-dd*r2G))*detinv
      else if( istat==MC_PLASTIC_LEFT ) then
        dpsdpe(m1,m1) = Cdiag+(CA1*((r2K*(db-dd)-(db*r4Gd3+dd*r2Gd3))*sinphi-dd*r2G) &
            &                 +CA2*((r2K*(da-dc)-(da*r4Gd3+dc*r2Gd3))*sinphi-dc*r2G))*detinv
        dpsdpe(m1,m2) = Coffd+(CA1*((r2K*(db-dd)+(db*r2Gd3+dd*r4Gd3))*sinphi+db*r2G) &
            &                 +CA2*((r2K*(da-dc)+(da*r2Gd3+dc*r4Gd3))*sinphi+da*r2G))*detinv
        dpsdpe(m1,m3) = Coffd+(CA1*(db-dd)+CA2*(da-dc))*(-r2G+(r2K+r2Gd3)*sinphi)*detinv
        dpsdpe(m2,m1) = Coffd+(CA1*((r2K*(dc-da)+(da*r4Gd3+dc*r2Gd3))*sinphi+dc*r2G) &
            &                 +CA2*((r2K*(dd-db)+(db*r4Gd3+dd*r2Gd3))*sinphi+dd*r2G))*detinv
        dpsdpe(m2,m2) = Cdiag+(CA1*((r2K*(dc-da)-(da*r2Gd3+dc*r4Gd3))*sinphi-da*r2G) &
            &                 +CA2*((r2K*(dd-db)-(db*r2Gd3+dd*r4Gd3))*sinphi-db*r2G))*detinv
        dpsdpe(m2,m3) = Coffd+(CA1*(dc-da)+CA2*(dd-db))*(-r2G+(r2K+r2Gd3)*sinphi)*detinv
        dpsdpe(m3,m1) = Coffd+CA3*((r2K*(-db+dd+da-dc)+(db-da)*r4Gd3+(dd-dc)*r2Gd3)*sinphi+(dd-dc)*r2G)*detinv
        dpsdpe(m3,m2) = Coffd+CA3*((r2K*(-db+dd+da-dc)+(da-db)*r2Gd3+(dd-dc)*r4Gd3)*sinphi+(da-db)*r2G)*detinv
        dpsdpe(m3,m3) = Cdiag+CA3*(-db+dd+da-dc)*(-r2G+(r2K+r2Gd3)*sinphi)*detinv
      endif
    endif
    ! compute principal elastic strain from principal stress
    prnstra(1) = (prnstre(1)-poisson*(prnstre(2)+prnstre(3)))/youngs
    prnstra(2) = (prnstre(2)-poisson*(prnstre(1)+prnstre(3)))/youngs
    prnstra(3) = (prnstre(3)-poisson*(prnstre(1)+prnstre(2)))/youngs
    call deriv_general_iso_tensor_func_3d(dpsdpe, D, prnprj, prnstra, prnstre)
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
    real(kind=kreal) :: dum, a(6), dlambda, G, K
    real(kind=kreal) :: J1,J2, eta, xi, etabar, harden, devia(6)
    real(kind=kreal) :: alpha, beta, C1, C2, C3, C4, CA, devia_norm

    if( sectType /=D3 ) stop "Elastoplastic calculation support only Solid element currently"

    call calElasticMatrix( matl, sectTYPE, D, temperature, hdflag=hdflag )
    if( istat == DP_ELASTIC ) return   ! elastic state
    if( hdflag == 2 ) return

    harden = calHardenCoeff( matl, extval(1), temperature )

    G = D(4,4)
    K = D(1,1)-(4.d0/3.d0)*G

    eta = matl%variables(M_PLCONST3)
    xi = matl%variables(M_PLCONST4)
    etabar = matl%variables(M_PLCONST5)

    if( istat==DP_PLASTIC_SURF ) then
      J1 = (stress(1)+stress(2)+stress(3))
      devia(1:3) = stress(1:3)-J1/3.d0
      devia(4:6) = stress(4:6)
      J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
          dot_product( devia(4:6), devia(4:6) )

      devia_norm = sqrt(2.d0*J2)
      a(1:6) = devia(1:6)/devia_norm
      dlambda = extval(1)-plstrain
      CA = 1.d0 / (G + K*eta*etabar + xi*xi*harden)
      dum = sqrt(2.d0)*devia_norm
      C1 = 4.d0*G*G*dlambda/dum
      C2 = 2.d0*G*(2.d0*G*dlambda/dum - G*CA)
      C3 = sqrt(2.d0)*G*CA*K
      C4 = K*K*eta*etabar*CA
      do j=1,6
        do i=1,6
          D(i,j) = D(i,j) - C1*Id(i,j) + C2*a(i)*a(j) &
              - C3*(eta*a(i)*I2(j) + etabar*I2(i)*a(j)) &
              - C4*I2(i)*I2(j)
        enddo
      enddo
    else ! istat==DP_PLASTIC_APEX
      alpha = xi/etabar
      beta = xi/eta
      C1 = K*(1.d0 - K/(K + alpha*beta*harden))
      do j=1,6
        do i=1,6
          D(i,j) = C1*I2(i)*I2(j)
        enddo
      enddo
    endif
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

    real(kind=kreal), parameter :: tol =1.d-8
    integer, parameter          :: MAXITER = 10
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

    if( abs(f/yd)<tol ) then  ! yielded
      istat = VM_PLASTIC
      return
    elseif( f<0.d0 ) then   ! not yielded or unloading
      istat = VM_ELASTIC
      return
    endif
    if( hdflag == 2 ) return

    istat = VM_PLASTIC      ! yielded
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
        istat=VM_ELASTIC; exit
      endif
      pstrain = plstrain+dlambda
      yd = calCurrYield( matl, pstrain, temp )
      if( kinematic ) then
        KK = calCurrKinematic( matl, pstrain )
      endif
      f = eqvs-3.d0*G*dlambda-yd -(KK-betan)
      if( abs(f/yd)<tol ) exit
      ! if( i==MAXITER ) then
      !   stop 'ERROR: BackwardEuler_VM: convergence failure'
      ! endif
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

    real(kind=kreal), parameter :: tol =1.d-8
    integer, parameter          :: MAXITER = 10
    real(kind=kreal) :: dlambda, f, mat(3,3)
    integer :: i, m1, m2, m3
    real(kind=kreal) :: youngs, poisson, pstrain, ina(1), ee(2)
    real(kind=kreal) :: H, dd, eqvs, cohe, G, K
    real(kind=kreal) :: prnstre(3), prnprj(3,3), tstre(3,3)
    real(kind=kreal) :: phi, psi, trialprn(3)
    logical          :: ierr
    real(kind=kreal) :: C1, C2, CS1, CS2, CS3
    real(kind=kreal) :: sinphi, cosphi, sinpsi, sphsps, r2cosphi, r4cos2phi, cotphi
    real(kind=kreal) :: da, db, dc, depv, detinv, dlambdb, dum, eps, eqvsb, fb
    real(kind=kreal) :: pt, p, resid

    phi = matl%variables(M_PLCONST3)
    psi = matl%variables(M_PLCONST4)
    sinphi = sin(phi)
    cosphi = cos(phi)
    r2cosphi = 2.d0*cosphi

    call eigen3( stress, prnstre, prnprj )
    trialprn = prnstre
    m1 = maxloc( prnstre, 1 )
    m3 = minloc( prnstre, 1 )
    if( m1 == m3 ) then
      m1 = 1; m2 = 2; m3 = 3
    else
      m2 = 6 - (m1 + m3)
    endif

    eqvs = prnstre(m1)-prnstre(m3) + (prnstre(m1)+prnstre(m3))*sinphi
    cohe = calCurrYield( matl, plstrain, temp )
    f = eqvs - r2cosphi*cohe

    if( abs(f/cohe)<tol ) then  ! yielded
      istat = MC_PLASTIC_SURF
      return
    elseif( f<0.d0 ) then   ! not yielded or unloading
      istat = MC_ELASTIC
      return
    endif
    if( hdflag == 2 ) return

    istat = MC_PLASTIC_SURF   ! yielded

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

    sinpsi = sin(psi)
    sphsps = sinphi*sinpsi
    r4cos2phi = r2cosphi*r2cosphi
    C1 = 4.d0*(G*(1.d0+sphsps/3.d0)+K*sphsps)
    do i=1,MAXITER
      H= calHardenCoeff( matl, pstrain, temp )
      dd= C1 + r4cos2phi*H
      dlambda = dlambda+f/dd
      if( r2cosphi*dlambda<0.d0 ) then
        if( cosphi==0.d0 ) stop "Math error in return mapping"
        dlambda = 0.d0
        pstrain = plstrain
        istat = MC_ELASTIC; exit
      endif
      pstrain = plstrain + r2cosphi*dlambda
      cohe = calCurrYield( matl, pstrain, temp )
      f = eqvs - C1*dlambda - r2cosphi*cohe
      if( abs(f/cohe)<tol ) exit
      ! if( i==MAXITER ) then
      !   stop 'ERROR: BackwardEuler_MC: convergence failure'
      ! endif
    enddo
    CS1 =2.d0*G*(1.d0+sinpsi/3.d0) + 2.d0*K*sinpsi
    CS2 =(4.d0*G/3.d0-2.d0*K)*sinpsi
    CS3 =2.d0*G*(1.d0-sinpsi/3.d0) - 2.d0*K*sinpsi
    prnstre(m1) = prnstre(m1)-CS1*dlambda
    prnstre(m2) = prnstre(m2)+CS2*dlambda
    prnstre(m3) = prnstre(m3)+CS3*dlambda
    eps = (abs(prnstre(m1))+abs(prnstre(m2))+abs(prnstre(m3)))*tol
    if( prnstre(m1) < prnstre(m2)-eps .or. prnstre(m2) < prnstre(m3)-eps ) then
      ! return mapping to EDGE
      prnstre = trialprn
      dlambda = 0.d0
      dlambdb = 0.d0
      if( (1.d0-sinpsi)*prnstre(m1) - 2*prnstre(m2) + (1.d0+sinpsi)*prnstre(m3) > 0) then
        istat = MC_PLASTIC_RIGHT
        eqvsb = prnstre(m1)-prnstre(m2) + (prnstre(m1)+prnstre(m2))*sinphi
        C2 = 2.d0*G*(1.d0+sinphi+sinpsi-sphsps/3.d0) + 4.d0*K*sphsps
      else
        istat = MC_PLASTIC_LEFT
        eqvsb = prnstre(m2)-prnstre(m3) + (prnstre(m2)+prnstre(m3))*sinphi
        C2 = 2.d0*G*(1.d0-sinphi-sinpsi-sphsps/3.d0) + 4.d0*K*sphsps
      endif
      cohe = calCurrYield( matl, plstrain, temp )
      f = eqvs - r2cosphi*cohe
      fb = eqvsb - r2cosphi*cohe
      pstrain = plstrain
      do i=1,MAXITER
        H= calHardenCoeff( matl, pstrain, temp )
        dum = r4cos2phi*H
        da = C1 + dum
        db = C2 + dum
        dc = db
        dd = da
        detinv = 1.d0/(da*dd-db*dc)
        dlambda = dlambda + detinv*( dd*f - db*fb)
        dlambdb = dlambdb + detinv*(-dc*f + da*fb)
        pstrain = plstrain + r2cosphi*(dlambda+dlambdb)
        cohe = calCurrYield( matl, pstrain, temp )
        f = eqvs - C1*dlambda - C2*dlambdb - r2cosphi*cohe
        fb = eqvsb - C2*dlambda - C1*dlambdb - r2cosphi*cohe
        if( (abs(f)+abs(fb))/(abs(eqvs)+abs(eqvsb)) < tol ) exit
        ! if( i==MAXITER ) then
        !   stop 'ERROR: BackwardEuler_MC: convergence failure(2)'
        ! endif
      enddo
      if( istat==MC_PLASTIC_RIGHT ) then
        prnstre(m1) = prnstre(m1)-CS1*(dlambda+dlambdb)
        prnstre(m2) = prnstre(m2)+CS2*dlambda+CS3*dlambdb
        prnstre(m3) = prnstre(m3)+CS3*dlambda+CS2*dlambdb
      else
        prnstre(m1) = prnstre(m1)-CS1*dlambda+CS2*dlambdb
        prnstre(m2) = prnstre(m2)+CS2*dlambda-CS1*dlambdb
        prnstre(m3) = prnstre(m3)+CS3*(dlambda+dlambdb)
      endif
      eps = (abs(prnstre(m1))+abs(prnstre(m2))+abs(prnstre(m3)))*tol
      if( prnstre(m1) < prnstre(m2)-eps .or. prnstre(m2) < prnstre(m3)-eps ) then
        ! return mapping to APEX
        prnstre = trialprn
        istat = MC_PLASTIC_APEX
        if( sinphi==0.d0 ) stop 'ERROR: BackwardEuler_MC: phi==0.0'
        if( sinpsi==0.d0 ) stop 'ERROR: BackwardEuler_MC: psi==0.0'
        depv = 0.d0
        cohe = calCurrYield( matl, plstrain, temp )
        cotphi = cosphi/sinphi
        pt = (stress(1)+stress(2)+stress(3))/3.d0
        resid = cotphi*cohe - pt
        pstrain = plstrain
        do i=1,MAXITER
          H= calHardenCoeff( matl, pstrain, temp )
          dd= cosphi*cotphi*H/sinpsi + K
          depv = depv - resid/dd
          pstrain = plstrain + cosphi*depv/sinpsi
          cohe = calCurrYield( matl, pstrain,temp )
          p = pt-K*depv
          resid = cotphi*cohe-p
          if( abs(resid/cohe)<tol ) exit
          ! if( i==MAXITER ) then
          !   stop 'ERROR: BackwardEuler_MC: convergence failure(3)'
          ! endif
        enddo
        prnstre(m1) = p
        prnstre(m2) = p
        prnstre(m3) = p
      endif
    endif
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
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(in)     :: temp  !< temperature
    integer(kind=kint), intent(in)   :: hdflag  !> return only hyd and dev term if specified

    real(kind=kreal), parameter :: tol =1.d-8
    integer, parameter          :: MAXITER = 10
    real(kind=kreal) :: dlambda, f
    integer :: i
    real(kind=kreal) :: youngs, poisson, pstrain, xi, ina(1), ee(2)
    real(kind=kreal) :: J1,J2,H, dd, eqvst, eqvs, cohe, G, K, devia(6), eta, etabar, pt, p
    logical          :: ierr
    real(kind=kreal) :: alpha, beta, depv, factor, resid

    eta = matl%variables(M_PLCONST3)
    xi = matl%variables(M_PLCONST4)
    etabar = matl%variables(M_PLCONST5)

    J1 = (stress(1)+stress(2)+stress(3))
    pt = J1/3.d0
    devia(1:3) = stress(1:3)-pt
    devia(4:6) = stress(4:6)
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    eqvst = sqrt(J2)
    cohe = calCurrYield( matl, plstrain, temp )
    f = eqvst + eta*pt - xi*cohe

    if( abs(f/cohe)<tol ) then  ! yielded
      istat = DP_PLASTIC_SURF
      return
    elseif( f<0.d0 ) then   ! not yielded or unloading
      istat = DP_ELASTIC
      return
    endif
    if( hdflag == 2 ) return

    istat = DP_PLASTIC_SURF

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
      dd= G+K*etabar*eta+H*xi*xi
      dlambda = dlambda+f/dd
      if( xi*dlambda<0.d0 ) then
        if( xi==0.d0 ) stop "Math error in return mapping"
        dlambda = 0.d0
        pstrain = plstrain
        istat=0; exit
      endif
      pstrain = plstrain+xi*dlambda
      cohe = calCurrYield( matl, pstrain, temp  )
      eqvs = eqvst-G*dlambda
      p = pt-K*etabar*dlambda
      f = eqvs + eta*p- xi*cohe
      if( abs(f/cohe)<tol ) exit
      ! if( i==MAXITER ) then
      !   stop 'ERROR: BackwardEuler_DP: convergence failure'
      ! endif
    enddo
    if( eqvs>=0.d0 ) then ! converged
      factor = 1.d0-G*dlambda/eqvst
    else                  ! return mapping to APEX
      istat = DP_PLASTIC_APEX
      if( eta==0.d0 ) stop 'ERROR: BackwardEuler_DP: eta==0.0'
      if( etabar==0.d0 ) stop 'ERROR: BackwardEuler_DP: etabar==0.0'
      alpha = xi/etabar
      beta = xi/eta
      depv=0.d0
      pstrain = plstrain
      cohe = calCurrYield( matl, pstrain, temp )
      resid = beta*cohe - pt
      do i=1,MAXITER
        H= calHardenCoeff( matl, pstrain, temp )
        dd= alpha*beta*H + K
        depv = depv - resid/dd
        pstrain = plstrain+alpha*depv
        cohe = calCurrYield( matl, pstrain, temp )
        p = pt-K*depv
        resid = beta*cohe - p
        if( abs(resid/cohe)<tol ) then
          dlambda=depv/etabar
          factor=0.d0
          exit
        endif
        ! if( i==MAXITER ) then
        !   stop 'ERROR: BackwardEuler_DP: convergence failure(2)'
        ! endif
      enddo
    endif
    devia(:) = factor*devia(:)
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
