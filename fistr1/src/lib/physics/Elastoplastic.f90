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
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private
  public :: calElastoPlasticMatrix
  public :: BackwardEuler
  public :: updateEPState
  public :: Update_PlaneStressJ2
  public :: PlaneStressJ2Tangent
  public :: checkPlaneStressJ2Material, PlaneStressJ2ErrorMessage
  integer, parameter, public :: J2_UNSUPPORTED=1, J2_BAD_STATE=2, J2_BAD_PROPERTY=3
  integer, parameter, public :: J2_SINGULAR=4, J2_NO_CONVERGENCE=5, J2_NONFINITE=6, J2_TEMPERATURE=7

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

  real(kind=kreal), parameter :: SHELL_SHEAR_CORRECTION = 5.0d0/6.0d0
  integer, parameter :: PLANE_STRESS_COMPONENTS(3) = (/ 1, 2, 4 /)

contains

  logical function isPlaneStressJ2PerfectPlasticity( matl )
    type(tMaterial), intent(in) :: matl

    isPlaneStressJ2PerfectPlasticity = .false.
    if( .not. isElastoplastic(matl%mtype) ) return
    if( getElasticType(matl%mtype) /= 0 ) return
    if( getYieldFunction(matl%mtype) /= 0 ) return
    if( getHardenType(matl%mtype) /= 0 ) return
    if( abs(matl%variables(M_PLCONST2)) > tiny(1.0d0) ) return
    isPlaneStressJ2PerfectPlasticity = .true.
  end function isPlaneStressJ2PerfectPlasticity

  subroutine checkPlaneStressJ2Material(matl, ierr)
    type(tMaterial), intent(in) :: matl
    integer(kind=kint), intent(out) :: ierr
    type(tTable), pointer :: table
    real(kind=kreal) :: constants(4), elastic(2)
    logical :: missing

    ierr = J2_UNSUPPORTED
    if( .not. isPlaneStressJ2PerfectPlasticity(matl) ) return
    if( associated(matl%shell_var) ) then
      if( any(matl%shell_var%ortho /= 0) ) return
    endif
    call fetch_Table(MC_ISOELASTIC, matl%dict, table, missing)
    if( .not. missing ) then
      if( table%ndepends /= 0 ) then
        ierr = J2_TEMPERATURE
        return
      endif
    endif
    call fetch_TableData(MC_ISOELASTIC, matl%dict, elastic, missing)
    ierr = J2_BAD_PROPERTY
    if( missing ) return
    constants = (/ elastic, matl%variables(M_PLCONST1), matl%variables(M_PLCONST2) /)
    if( .not. all(ieee_is_finite(constants)) ) return
    if( constants(1) <= 0.0d0 .or. constants(2) <= -1.0d0 .or. constants(2) >= 0.5d0 ) return
    if( constants(3) <= 0.0d0 ) return
    ierr = 0
  end subroutine checkPlaneStressJ2Material

  function PlaneStressJ2ErrorMessage(ierr) result(message)
    integer(kind=kint), intent(in) :: ierr
    character(len=128) :: message

    select case(ierr)
    case(J2_UNSUPPORTED)
      message = 'Plane-stress J2 requires isotropic elasticity and perfect plasticity'
    case(J2_BAD_STATE)
      message = 'Invalid plane-stress J2 history state'
    case(J2_BAD_PROPERTY)
      message = 'Invalid or missing plane-stress J2 properties: require E > 0, -1 < nu < 0.5, yield stress > 0'
    case(J2_SINGULAR)
      message = 'Singular thickness tangent in plane-stress J2 condensation'
    case(J2_NO_CONVERGENCE)
      message = 'Plane-stress J2 local iteration did not converge'
    case(J2_NONFINITE)
      message = 'Non-finite value in plane-stress J2 material update'
    case(J2_TEMPERATURE)
      message = 'Temperature-dependent elasticity is not supported by plane-stress J2 shells'
    case default
      message = 'Unknown plane-stress J2 error'
    end select
  end function PlaneStressJ2ErrorMessage

  subroutine condensePlaneStressJ2Tangent( tangent3d, shear_stiffness, tangent, ierr )
    real(kind=kreal), intent(in) :: tangent3d(6,6)
    real(kind=kreal), intent(in) :: shear_stiffness
    real(kind=kreal), intent(out) :: tangent(5,5)
    integer(kind=kint), intent(out) :: ierr

    integer :: i, j, ii, jj

    ierr = 0
    tangent = 0.0d0
    if( .not. all(ieee_is_finite(tangent3d)) ) then
      ierr = J2_NONFINITE
      return
    endif
    if( abs(tangent3d(3,3)) <= 100.0d0*epsilon(1.0d0)*maxval(abs(tangent3d)) ) then
      ierr = J2_SINGULAR
      return
    endif

    do i = 1, 3
      ii = PLANE_STRESS_COMPONENTS(i)
      do j = 1, 3
        jj = PLANE_STRESS_COMPONENTS(j)
        tangent(i,j) = tangent3d(ii,jj)-tangent3d(ii,3)*tangent3d(3,jj)/tangent3d(3,3)
      enddo
    enddo
    tangent(4,4) = shear_stiffness
    tangent(5,5) = shear_stiffness
  end subroutine condensePlaneStressJ2Tangent

  !> Return the plane-stress response of the existing three-dimensional J2
  !> perfect-plastic material update. The strain increment and stress use the
  !> shell-local engineering ordering (11, 22, 33, 12, 23, 31).
  !> Transverse shear remains elastic and is excluded from the yield function.
  subroutine Update_PlaneStressJ2( matl, strain_increment, stress_bak, plstrain, &
      fstatus, stress, strain_increment_out, tangent, istat, fstatus_out, plpotential, energy_increment, ierr )
    type(tMaterial), intent(in) :: matl
    real(kind=kreal), intent(in) :: strain_increment(6), stress_bak(6)
    real(kind=kreal), intent(in) :: plstrain, fstatus(:)
    real(kind=kreal), intent(out) :: stress(6), strain_increment_out(6), tangent(5,5)
    integer, intent(out) :: istat
    real(kind=kreal), intent(out) :: fstatus_out(:), plpotential, energy_increment
    integer(kind=kint), intent(out) :: ierr

    integer, parameter :: maxiter = 25
    integer :: iter
    real(kind=kreal) :: elastic(6,6), tangent3d(6,6)
    real(kind=kreal) :: strain_increment_work(6)
    real(kind=kreal) :: stress_work(6), stress_bak_work(6)
    real(kind=kreal) :: shear_stiffness, residual, tolerance

    ierr = 0
    stress = 0.0d0
    strain_increment_out = strain_increment
    tangent = 0.0d0
    istat = 0
    fstatus_out = 0.0d0
    plpotential = 0.0d0
    energy_increment = 0.0d0

    call checkPlaneStressJ2Material(matl, ierr)
    if( ierr /= 0 ) return
    if( size(fstatus) /= size(fstatus_out) .or. size(fstatus) < 1 ) then
      ierr = J2_BAD_STATE
      return
    endif
    if( .not. all(ieee_is_finite(strain_increment)) .or. .not. all(ieee_is_finite(stress_bak)) .or. &
        .not. all(ieee_is_finite(fstatus)) .or. .not. ieee_is_finite(plstrain) ) then
      ierr = J2_NONFINITE
      return
    endif
    fstatus_out = fstatus

    ! Only the in-plane response enters the three-dimensional J2 update.
    strain_increment_work = strain_increment
    strain_increment_work(5:6) = 0.0d0
    stress_bak_work = stress_bak
    stress_bak_work(5:6) = 0.0d0
    strain_increment_work(3) = 0.0d0
    call calElasticMatrix( matl, D3, elastic, 0.0d0 )
    shear_stiffness = SHELL_SHEAR_CORRECTION*elastic(5,5)

    do iter = 1, maxiter
      ! Each plane-stress iteration starts from the committed J2 state.
      fstatus_out = fstatus
      fstatus_out(1) = plstrain
      stress_work = stress_bak_work+matmul(elastic, strain_increment_work)
      call BackwardEuler( matl, stress_work, plstrain, istat, fstatus_out, plpotential, 0.0d0 )
      call calElastoPlasticMatrix( matl, D3, stress_work, istat, fstatus_out, plstrain, tangent3d, 0.0d0 )

      if( .not. all(ieee_is_finite(stress_work)) .or. .not. all(ieee_is_finite(tangent3d)) .or. &
          .not. all(ieee_is_finite(fstatus_out)) ) then
        ierr = J2_NONFINITE
        return
      endif

      residual = stress_work(3)
      tolerance = 1.0d-10*max(matl%variables(M_PLCONST1), maxval(abs(stress_work)))
      if( abs(residual) <= tolerance ) exit
      if( abs(tangent3d(3,3)) <= 100.0d0*epsilon(1.0d0)*maxval(abs(tangent3d)) ) then
        ierr = J2_SINGULAR
        return
      endif
      strain_increment_work(3) = strain_increment_work(3)-residual/tangent3d(3,3)
    enddo
    if( iter > maxiter ) then
      ierr = J2_NO_CONVERGENCE
      return
    endif

    call condensePlaneStressJ2Tangent(tangent3d, shear_stiffness, tangent, ierr)
    if( ierr /= 0 ) return

    stress = stress_work
    stress(3) = 0.0d0
    stress(5:6) = stress_bak(5:6)+shear_stiffness*strain_increment(5:6)
    strain_increment_out = strain_increment_work
    strain_increment_out(5:6) = strain_increment(5:6)
    ! Use the elastic predictor work before applying the plastic correction.
    energy_increment = dot_product(stress_bak_work, strain_increment_work) &
      +0.5d0*dot_product(matmul(elastic, strain_increment_work), strain_increment_work)+plpotential &
      +dot_product(stress_bak(5:6)+0.5d0*shear_stiffness*strain_increment(5:6), strain_increment(5:6))
    if( .not. ieee_is_finite(energy_increment) ) ierr = J2_NONFINITE
  end subroutine Update_PlaneStressJ2

  !> Condense the in-plane J2 tangent and retain the elastic transverse-shear stiffness.
  subroutine PlaneStressJ2Tangent( matl, stress, istat, fstatus, plstrain, tangent, ierr )
    type(tMaterial), intent(in) :: matl
    real(kind=kreal), intent(in) :: stress(6), fstatus(:), plstrain
    integer, intent(in) :: istat
    real(kind=kreal), intent(out) :: tangent(5,5)
    integer(kind=kint), intent(out) :: ierr

    real(kind=kreal) :: tangent3d(6,6), stress_work(6), elastic(6,6)

    ierr = 0
    tangent = 0.0d0
    call checkPlaneStressJ2Material(matl, ierr)
    if( ierr /= 0 ) return
    if( size(fstatus) < 1 ) then
      ierr = J2_BAD_STATE
      return
    endif
    if( .not. all(ieee_is_finite(stress)) .or. .not. all(ieee_is_finite(fstatus)) .or. &
        .not. ieee_is_finite(plstrain) ) then
      ierr = J2_NONFINITE
      return
    endif
    stress_work = stress
    stress_work(5:6) = 0.0d0
    call calElasticMatrix( matl, D3, elastic, 0.0d0 )
    call calElastoPlasticMatrix( matl, D3, stress_work, istat, fstatus, plstrain, tangent3d, 0.0d0 )
    call condensePlaneStressJ2Tangent(tangent3d, SHELL_SHEAR_CORRECTION*elastic(5,5), tangent, ierr)
  end subroutine PlaneStressJ2Tangent

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
    real(kind=kreal) :: prnstre(3), prnprj(3,3), prnstra(3)
    integer(kind=kint) :: m1, m2, m3
    real(kind=kreal) :: C1,C2, CA1, CA2, CA3, CAm, CAp, CD1, CD2, CD3, Cdiag, Coffd
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
  subroutine BackwardEuler( matl, stress, plstrain, istat, fstat, plpotential, temp, hdflag )
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(inout)  :: plpotential    !< plastic potential
    real(kind=kreal), intent(in)     :: temp  !< temperature
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: ytype, hdflag_in

    hdflag_in = 0
    if( present(hdflag) ) hdflag_in = hdflag

    ytype = getYieldFunction( matl%mtype )
    select case (ytype)
    case (0)
      call BackwardEuler_VM( matl, stress, plstrain, istat, fstat, plpotential, temp, hdflag_in )
    case (1)
      call BackwardEuler_MC( matl, stress, plstrain, istat, fstat, temp, hdflag_in )
    case (2)
      call BackwardEuler_DP( matl, stress, plstrain, istat, fstat, temp, hdflag_in )
    case (3)
      call uBackwardEuler( matl%variables, stress, plstrain, istat, fstat, temp, hdflag_in )
    end select
  end subroutine BackwardEuler

  !> This subroutine does backward-Euler return calculation for von Mises
  subroutine BackwardEuler_VM( matl, stress, plstrain, istat, fstat, plpotential, temp, hdflag )
    type( tMaterial ), intent(in)    :: matl        !< material properties
    real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
    real(kind=kreal), intent(in)     :: plstrain    !< plastic strain till current substep
    integer, intent(inout)           :: istat       !< plastic state
    real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
    real(kind=kreal), intent(inout)  :: plpotential    !< plastic potential
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

    plpotential = 0.d0
    kinematic = isKinematicHarden( matl%mtype )
    if( kinematic ) back(1:6) = fstat(8:13)

    J1 = (stress(1)+stress(2)+stress(3))
    devia(1:3) = stress(1:3)-J1/3.d0
    devia(4:6) = stress(4:6)
    if( kinematic ) devia = devia-back
    J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
      dot_product( devia(4:6), devia(4:6) )

    eqvs = dsqrt( 3.d0*J2 )
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
    KH = 0.d0; KK=0.d0; betan=0.d0
    if( kinematic ) then
      betan = calCurrKinematic( matl, plstrain )  ! keep back = alpha_n (loaded above) so it accumulates at fstat update and is restored into stress
    else
      back(:)=0.d0
    endif

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

    H= calHardenCoeff( matl, pstrain, temp ) !a
    yd = calCurrYield( matl, plstrain, temp ) !b
    plpotential = -0.5d0*(eqvs-yd)*(eqvs-yd)/(H+3.d0*G)

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
