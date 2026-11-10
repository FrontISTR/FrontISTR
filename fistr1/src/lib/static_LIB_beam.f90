!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provide common functions of beam elements
module m_static_LIB_beam

  use hecmw
  use mMechGauss
  use m_utilities
  use elementinfo
  use m_fstr, only: kel611TIMOSHENKO

  implicit none

contains

  !< This subroutine build up coordinate frame for beam elements
  subroutine framtr(refx, xl, le, t)
    real(kind=kreal), intent(in)  :: refx(3)      !< Reference Vector
    real(kind=kreal), intent(in)  :: xl(3,2)      !< Element coordinates
    real(kind=kreal), intent(out) :: le           !< length of the element
    real(kind=kreal), intent(out) :: t(3,3)       !< Transformation array

    real(kind=kreal) ::       dl
    real(kind=kreal), parameter ::  tol = 1.d-08

    t(1,1) = xl(1,2) - xl(1,1)
    t(1,2) = xl(2,2) - xl(2,1)
    t(1,3) = xl(3,2) - xl(3,1)
    le  = sqrt(t(1,1)*t(1,1)+t(1,2)*t(1,2)+t(1,3)*t(1,3))
    dl  = 1.0d0/le
    t(1,1) = t(1,1)*dl
    t(1,2) = t(1,2)*dl
    t(1,3) = t(1,3)*dl

    t(3,1) = refx(1)
    t(3,2) = refx(2)
    t(3,3) = refx(3)


    t(2,1) = (t(3,2)*t(1,3) - t(3,3)*t(1,2))
    t(2,2) = (t(3,3)*t(1,1) - t(3,1)*t(1,3))
    t(2,3) = (t(3,1)*t(1,2) - t(3,2)*t(1,1))
    dl  = sqrt(t(2,1)*t(2,1)+t(2,2)*t(2,2)+t(2,3)*t(2,3))
    if(dl<tol*le) then
      stop  "Bad reference for beam element!"
    else
      t(2,1) = t(2,1)/dl
      t(2,2) = t(2,2)/dl
      t(2,3) = t(2,3)/dl
      t(3,1) = t(1,2)*t(2,3) - t(1,3)*t(2,2)
      t(3,2) = t(1,3)*t(2,1) - t(1,1)*t(2,3)
      t(3,3) = t(1,1)*t(2,2) - t(1,2)*t(2,1)
    endif

  end subroutine framtr

  !> Calculate the local stiffness matrix of 611 beam elements
  subroutine STF_Beam_local(le, section, E, P, formulation, stiff)
    real(kind=kreal), intent(in)  :: le
    real(kind=kreal), intent(in)  :: section(:)
    real(kind=kreal), intent(in)  :: E, P
    integer(kind=kint), intent(in), optional :: formulation
    real(kind=kreal), intent(out) :: stiff(12,12)

    real(kind=kreal) :: G, L2, L3, A, Iy, Iz, Jx, EA
    real(kind=kreal) :: phi_y, phi_z
    real(kind=kreal) :: kyy, kyz, kzz, kryz
    real(kind=kreal) :: kww, kwy, krr, krz

    L2 = le*le
    L3 = L2*le
    G = E/(2.d0*(1.d0 + P))

    A = section(4); Iy = section(5); Iz = section(6); Jx = section(7)

    phi_y = 0.d0
    phi_z = 0.d0
    if( present(formulation) ) then
      if( formulation == kel611TIMOSHENKO ) then
        ! The BEAM section has no independent shear areas; use A in both planes.
        phi_y = 12.d0*E*Iz/(G*A*L2)
        phi_z = 12.d0*E*Iy/(G*A*L2)
      endif
    endif

    EA = E*A/le
    kyy = 12.d0*E*Iz/(L3*(1.d0 + phi_y))
    kyz = 6.d0*E*Iz/(L2*(1.d0 + phi_y))
    kzz = (4.d0 + phi_y)*E*Iz/(le*(1.d0 + phi_y))
    kryz = (2.d0 - phi_y)*E*Iz/(le*(1.d0 + phi_y))
    kww = 12.d0*E*Iy/(L3*(1.d0 + phi_z))
    kwy = 6.d0*E*Iy/(L2*(1.d0 + phi_z))
    krr = (4.d0 + phi_z)*E*Iy/(le*(1.d0 + phi_z))
    krz = (2.d0 - phi_z)*E*Iy/(le*(1.d0 + phi_z))

    stiff = 0.d0
    stiff(1,1) = EA
    stiff(7,1) = -EA

    stiff(2,2) = kyy
    stiff(6,2) = kyz
    stiff(8,2) = -kyy
    stiff(12,2) = kyz

    stiff(3,3) = kww
    stiff(5,3) = -kwy
    stiff(9,3) = -kww
    stiff(11,3) = -kwy

    stiff(4,4) = G*Jx/le
    stiff(10,4) = -G*Jx/le

    stiff(3,5) = -kwy
    stiff(5,5) = krr
    stiff(9,5) = kwy
    stiff(11,5) = krz

    stiff(2,6) = kyz
    stiff(6,6) = kzz
    stiff(8,6) = -kyz
    stiff(12,6) = kryz

    stiff(1,7) = -EA
    stiff(7,7) = EA

    stiff(2,8) = -kyy
    stiff(6,8) = -kyz
    stiff(8,8) = kyy
    stiff(12,8) = -kyz

    stiff(3,9) = -kww
    stiff(5,9) = kwy
    stiff(9,9) = kww
    stiff(11,9) = kwy

    stiff(4,10) = -G*Jx/le
    stiff(10,10) = G*Jx/le

    stiff(3,11) = -kwy
    stiff(5,11) = krz
    stiff(9,11) = kwy
    stiff(11,11) = krr

    stiff(2,12) = kyz
    stiff(6,12) = kryz
    stiff(8,12) = -kyz
    stiff(12,12) = kzz
  end subroutine STF_Beam_local

  !> Calculate stiff matrix of BEAM elements
  subroutine STF_Beam(nn,ecoord,section,E,P,STIFF,formulation)
    integer, intent(in)            :: nn           !< number of element's nodes
    real(kind=kreal), intent(in)   :: ecoord(3,nn) !< coordinates of elemental nodes
    real(kind=kreal), intent(in)   :: section(:)   !< section parameters
    real(kind=kreal), intent(in)   :: E,P   !< status of qudrature points
    real(kind=kreal), intent(out)  :: STIFF(nn*6,nn*6)   !< elemental stiff matrix
    integer(kind=kint), intent(in), optional :: formulation

    real(kind=kreal) :: le, trans(3,3), refv(3), transt(3,3)

    refv = section(1:3)
    call framtr(refv, ecoord, le, trans)
    transT= transpose(trans)

    call STF_Beam_local(le, section, E, P, formulation, stiff)

    stiff(1:3,:) = matmul( transT, stiff(1:3,:) )
    stiff(4:6,:) = matmul( transT, stiff(4:6,:) )
    stiff(7:9,:) = matmul( transT, stiff(7:9,:) )
    stiff(10:12,:) = matmul( transT, stiff(10:12,:) )

    stiff(:,1:3) = matmul( stiff(:,1:3), trans )
    stiff(:,4:6) = matmul( stiff(:,4:6), trans )
    stiff(:,7:9) = matmul( stiff(:,7:9), trans )
    stiff(:,10:12) = matmul( stiff(:,10:12), trans )

  end subroutine STF_Beam

  !####################################################################
  subroutine UpdateST_Beam(nn,ecoord,u,du,section,gausses,QF,formulation)
    integer, intent(in)            :: nn           !< number of element's nodes
    real(kind=kreal), intent(in)   :: ecoord(3,nn) !< coordinates of elemental nodes
    real(kind=kreal), intent(in)   :: u(6,nn)      !< displacement of elemental nodes
    real(kind=kreal), intent(in)   :: du(6,nn)     !< coordinates of elemental nodes
    real(kind=kreal), intent(in)   :: section(:)   !< section parameters
    type(tGaussStatus), intent(in) :: gausses(:)         !< status of Gaussian qudrature points
    real(kind=kreal), intent(out)  :: QF(nn*6)     !< elemental force matrix
    integer(kind=kint), intent(in), optional :: formulation

    real(kind=kreal)   :: stiff(nn*6, nn*6), totaldisp(nn*6)
    integer(kind=kint) :: i, j
    real(kind=kreal)   :: E,P   !< status of qudrature points

    E = gausses(1)%pMaterial%variables(M_YOUNGS)
    P = gausses(1)%pMaterial%variables(M_POISSON)

    call STF_Beam(nn,ecoord,section,E,P,STIFF,formulation)

    do i=1,nn
      do j=1,6
        totaldisp(6*(i-1)+j) = u(j,i) + du(j,i)
      end do
    end do

    qf = matmul(stiff,totaldisp)

  end subroutine UpdateST_Beam

  subroutine STF_Beam_641_from_611(ecoord, gausses, section, stiff, formulation)
    use mMechGauss
    implicit none

    real(kind=kreal), intent(in) :: ecoord(3, 4), section(:)
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: stiff(12, 12)
    integer(kind=kint), intent(in) :: formulation

    integer(kind=kint) :: i, j
    integer(kind=kint), parameter :: mixed_to_natural(12) = &
      (/ 1, 2, 3, 7, 8, 9, 4, 5, 6, 10, 11, 12 /)
    real(kind=kreal) :: natural_stiff(12, 12), ee, pp

    ee = gausses(1)%pMaterial%variables(M_YOUNGS)
    pp = gausses(1)%pMaterial%variables(M_POISSON)
    call STF_Beam(2, ecoord(1:3,1:2), section, ee, pp, natural_stiff, formulation)

    do j = 1, 12
      do i = 1, 12
        stiff(i,j) = natural_stiff(mixed_to_natural(i), mixed_to_natural(j))
      enddo
    enddo
  end subroutine STF_Beam_641_from_611

  subroutine UpdateST_Beam_641_from_611(ecoord, u, du, gausses, section, qf, formulation)
    use mMechGauss
    implicit none

    real(kind=kreal), intent(in) :: ecoord(3, 4), u(3, 4), du(3, 4), section(:)
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out) :: qf(12)
    integer(kind=kint), intent(in) :: formulation

    integer(kind=kint) :: i
    integer(kind=kint), parameter :: mixed_to_natural(12) = &
      (/ 1, 2, 3, 7, 8, 9, 4, 5, 6, 10, 11, 12 /)
    real(kind=kreal) :: beam_u(6, 2), beam_du(6, 2), natural_qf(12)

    beam_u(1:3,1:2) = u(1:3,1:2)
    beam_u(4:6,1:2) = u(1:3,3:4)
    beam_du(1:3,1:2) = du(1:3,1:2)
    beam_du(4:6,1:2) = du(1:3,3:4)
    call UpdateST_Beam(2, ecoord(1:3,1:2), beam_u, beam_du, section, &
      gausses, natural_qf, formulation)

    do i = 1, 12
      qf(i) = natural_qf(mixed_to_natural(i))
    enddo
  end subroutine UpdateST_Beam_641_from_611

  !> Calculate elemental section force (N, Q, M) of 611 beam elements
  subroutine NQM_Beam(nn, ecoord, gausses, section, ul, rnqm, formulation)
    use mMechGauss
    integer(kind=kint), intent(in)    :: nn              !< number of element's nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn)   !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)    :: gausses(:)      !< status of Gaussian quadrature points
    real(kind=kreal), intent(in)      :: section(:)      !< section parameters
    real(kind=kreal), intent(in)      :: ul(nn*6)        !< nodal displacement in the local frame
    real(kind=kreal), intent(out)     :: rnqm(nn*6)      !< elemental NQM (local frame)
    integer(kind=kint), intent(in), optional :: formulation

    real(kind=kreal) :: ee, pp, le, refv(3), trans(3,3), ec(3,2)
    real(kind=kreal) :: stiff(nn*6, nn*6)

    ee = gausses(1)%pMaterial%variables(M_YOUNGS)
    pp = gausses(1)%pMaterial%variables(M_POISSON)

    refv(1:3)  = section(1:3)
    ec(1:3, 1) = ecoord(1:3, 1)
    ec(1:3, 2) = ecoord(1:3, 2)
    call framtr(refv, ec, le, trans)

    call STF_Beam_local(le, section, ee, pp, formulation, stiff)

    rnqm = matmul(stiff, ul)

  end subroutine NQM_Beam

  !> Calculate NODAL STRESS and STRAIN of 611 beam elements
  subroutine NodalStress_Beam(nn, ecoord, gausses, section, edisp, ndstrain, ndstress, formulation)
    use mMechGauss
    integer(kind=kint), intent(in)    :: nn              !< number of element's nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn)   !< coordinates of elemental nodes
    type(tGaussStatus), intent(inout) :: gausses(:)      !< status of Gaussian quadrature points
    real(kind=kreal), intent(in)      :: section(:)      !< section parameters
    real(kind=kreal), intent(in)      :: edisp(6, nn)    !< nodal displacement (global frame)
    real(kind=kreal), intent(out)     :: ndstrain(nn, 6) !< nodal strain
    real(kind=kreal), intent(out)     :: ndstress(nn, 6) !< nodal stress
    integer(kind=kint), intent(in), optional :: formulation

    integer(kind=kint) :: k
    real(kind=kreal)   :: ee, pi
    real(kind=kreal)   :: radius, angle(6)
    real(kind=kreal)   :: refv(3), ec(3,2), trans(3,3), le
    real(kind=kreal)   :: ul(nn*6), rnqm(nn*6)
    real(kind=kreal)   :: x2_hat, x3_hat, eps
    real(kind=kreal)   :: stress_i, stress_j

    pi = 4.0D0*datan(1.0D0)

    ee = gausses(1)%pMaterial%variables(M_YOUNGS)

    refv(1:3)  = section(1:3)
    ec(1:3, 1) = ecoord(1:3, 1)
    ec(1:3, 2) = ecoord(1:3, 2)
    call framtr(refv, ec, le, trans)

    radius   = gausses(1)%pMaterial%variables(M_BEAM_RADIUS)
    angle(1) = gausses(1)%pMaterial%variables(M_BEAM_ANGLE1)
    angle(2) = gausses(1)%pMaterial%variables(M_BEAM_ANGLE2)
    angle(3) = gausses(1)%pMaterial%variables(M_BEAM_ANGLE3)
    angle(4) = gausses(1)%pMaterial%variables(M_BEAM_ANGLE4)
    angle(5) = gausses(1)%pMaterial%variables(M_BEAM_ANGLE5)
    angle(6) = gausses(1)%pMaterial%variables(M_BEAM_ANGLE6)

    !< nodal displacements/rotations in the local frame
    ul(1:3)   = matmul(trans, edisp(1:3, 1))
    ul(4:6)   = matmul(trans, edisp(4:6, 1))
    ul(7:9)   = matmul(trans, edisp(1:3, 2))
    ul(10:12) = matmul(trans, edisp(4:6, 2))

    eps = (ul(7)-ul(1))/le
    call NQM_Beam(nn, ecoord, gausses, section, ul, rnqm, formulation)

    ndstrain = 0.0D0
    ndstress = 0.0D0

    do k = 1, 6

      angle(k) = angle(k)/180.0D0*pi
      x2_hat = radius*dcos(angle(k))
      x3_hat = radius*dsin(angle(k))

      stress_i = ee*eps + x2_hat*rnqm(6)/section(6) - x3_hat*rnqm(5)/section(5)
      stress_j = ee*eps - x2_hat*rnqm(12)/section(6) + x3_hat*rnqm(11)/section(5)

      !< elemental value (evaluated at the element center)
      gausses(1)%strain(k) = eps
      gausses(1)%stress(k) = 0.5D0*(stress_i + stress_j)
      gausses(1)%strain_out(k) = gausses(1)%strain(k)
      gausses(1)%stress_out(k) = gausses(1)%stress(k)

      !< nodal value at node i
      ndstrain(1, k) = eps
      ndstress(1, k) = stress_i

      !< nodal value at node j
      ndstrain(2, k) = eps
      ndstress(2, k) = stress_j

    end do

    !< section force (N, Q, M)
    gausses(1)%nqm(1:nn*6) = rnqm(1:nn*6)

  end subroutine NodalStress_Beam

  !> Copy elemental stress/strain/NQM of 611 beam elements from Gauss status
  subroutine ElementalStress_Beam(gausses, estrain, estress, enqm)
    use mMechGauss
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: estrain(6)
    real(kind=kreal), intent(out)  :: estress(6)
    real(kind=kreal), intent(out)  :: enqm(12)

    estrain(1:6) = gausses(1)%strain_out(1:6)
    estress(1:6) = gausses(1)%stress_out(1:6)
    enqm(1:12)   = gausses(1)%nqm(1:12)

  end subroutine ElementalStress_Beam

  ! (Gaku Hashimoto, The University of Tokyo, 2014/02/06) <
  !####################################################################
  subroutine DL_Beam_641(etype, nn, xx, yy, zz, rho, ltype, params, &
      section, vect, nsize)
    !####################################################################
    !**
    !**  SET DLOAD
    !**
    !   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
    !   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
    !   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
    !   GRAV LTYPE=4  :GRAVITY FORCE
    !   CENT LTYPE=5  :CENTRIFUGAL LOAD
    !   P1   LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR FACE-1
    !   P2   LTYPE=20 :TRACTION IN NORMAL-DIRECTION FOR FACE-2
    !   P3   LTYPE=30 :TRACTION IN NORMAL-DIRECTION FOR FACE-3
    !   P4   LTYPE=40 :TRACTION IN NORMAL-DIRECTION FOR FACE-4
    !   P5   LTYPE=50 :TRACTION IN NORMAL-DIRECTION FOR FACE-5
    !   P6   LTYPE=60 :TRACTION IN NORMAL-DIRECTION FOR FACE-6
    ! I/F VARIABLES
    integer(kind = kint), intent(in)  :: etype, nn
    real(kind = kreal), intent(in)    :: xx(:), yy(:), zz(:)
    real(kind = kreal), intent(in)    :: params(0:6)
    real(kind = kreal), intent(in)    :: section(:)
    real(kind = kreal), intent(inout) :: vect(:)
    real(kind = kreal) :: rho
    integer(kind = kint) :: ltype, nsize
    ! LOCAL VARIABLES
    integer(kind = kint) :: ndof
    parameter(ndof = 3)
    integer(kind = kint) :: ivol, isuf, nod(nn)
    integer(kind = kint) :: i ,surtype, nsur
    real(kind = kreal) :: vx, vy, vz, val, a, aa

    !--------------------------------------------------------------------

    val = params(0)

    !--------------------------------------------------------------

    ivol = 0
    isuf = 0

    if( ltype .LT. 10 ) then

      ivol = 1

    else if( ltype .GE. 10 ) then

      isuf = 1

      call getSubFace(etype, ltype/10, surtype, nod)

      nsur = getNumberOfNodes(surtype)

    end if

    !--------------------------------------------------------------------

    nsize = nn*ndof

    !--------------------------------------------------------------------

    vect(1:nsize) = 0.0D0

    !--------------------------------------------------------------

    ! Volume force

    if( ivol .EQ. 1 ) then

      if( ltype .EQ. 4 ) then

        AA = dsqrt( ( xx(2)-xx(1) )*( xx(2)-xx(1) ) &
          +( yy(2)-yy(1) )*( yy(2)-yy(1) ) &
          +( zz(2)-zz(1) )*( zz(2)-zz(1) ) )

        a = section(4)

        vx = params(1)
        vy = params(2)
        vz = params(3)
        vx = vx/dsqrt( params(1)**2+params(2)**2+params(3)**2 )
        vy = vy/dsqrt( params(1)**2+params(2)**2+params(3)**2 )
        vz = vz/dsqrt( params(1)**2+params(2)**2+params(3)**2 )

        do i = 1, 2

          vect(3*i-2) = val*rho*a*0.5D0*AA*vx
          vect(3*i-1) = val*rho*a*0.5D0*AA*vy
          vect(3*i  ) = val*rho*a*0.5D0*AA*vz

        end do

        do i = 3, 4

          vect(3*i-2) = 0.0D0
          vect(3*i-1) = 0.0D0
          vect(3*i  ) = 0.0D0

        end do

      end if

    end if

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine DL_Beam_641
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2014/02/06)


  ! (Gaku Hashimoto, The University of Tokyo, 2014/02/06) <
  !####################################################################
  subroutine TLOAD_Beam_641 &
      (nn, ndof, xx, yy, zz, tt, t0, &
      gausses, section, reftemp, vect)
    !####################################################################

    use hecmw
    use m_fstr
    use m_utilities
    use mMechGauss

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: nn
    integer(kind = kint), intent(in) :: ndof
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(in)   :: section(:)
    real(kind = kreal), intent(in)   :: xx(nn), yy(nn), zz(nn)
    real(kind = kreal), intent(in)   :: tt(nn), t0(nn)
    real(kind = kreal), intent(in)   :: reftemp
    real(kind = kreal), intent(out)  :: vect(nn*ndof)

    !--------------------------------------------------------------------

    real(kind = kreal) :: tempc, temp0
    real(kind = kreal) :: ecoord(3, nn)
    real(kind = kreal) :: ec(3, 2)
    real(kind = kreal) :: ina1(1), outa1(2)
    real(kind = kreal) :: ina2(1), outa2(1)
    real(kind = kreal) :: alp, alp0
    real(kind = kreal) :: ee, pp
    real(kind = kreal) :: a
    real(kind = kreal) :: refv(3)
    real(kind = kreal) :: g
    real(kind = kreal) :: le
    real(kind = kreal) :: trans(3, 3), transt(3, 3)

    logical :: ierr

    !--------------------------------------------------------------------

    ecoord(1, 1:nn) = xx(1:nn)
    ecoord(2, 1:nn) = yy(1:nn)
    ecoord(3, 1:nn) = zz(1:nn)

    !--------------------------------------------------------------------

    tempc = 0.5D0*( tt(1)+tt(2) )
    temp0 = 0.5D0*( t0(1)+t0(2) )

    !--------------------------------------------------------------

    ina1(1) = tempc

    call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa1, ierr, ina1 )

    if( ierr ) then

      ee = gausses(1)%pMaterial%variables(M_YOUNGS)
      pp = gausses(1)%pMaterial%variables(M_POISSON)

    else

      ee = outa1(1)
      pp = outa1(2)

    end if

    !--------------------------------------------------------------

    ina2(1) = tempc

    call fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa2(:), ierr, ina2 )

    if( ierr ) stop "Fails in fetching expansion coefficient!"

    alp = outa2(1)

    !--------------------------------------------------------------

    ina2(1) = temp0

    call fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa2(:), ierr, ina2 )

    if( ierr ) stop "Fails in fetching expansion coefficient!"

    alp0 = outa2(1)

    !--------------------------------------------------------------------

    refv(1) = section(1)
    refv(2) = section(2)
    refv(3) = section(3)

    ec(1, 1) = ecoord(1, 1)
    ec(2, 1) = ecoord(2, 1)
    ec(3, 1) = ecoord(3, 1)
    ec(1, 2) = ecoord(1, 2)
    ec(2, 2) = ecoord(2, 2)
    ec(3, 2) = ecoord(3, 2)

    call framtr(refv, ec, le, trans)

    transt= transpose( trans )

    !--------------------------------------------------------------------

    a = section(4)

    g = ee/( 2.0D0*( 1.0D0+pp ))

    !--------------------------------------------------------------------

    vect( 1) = -a*ee*( alp*( tempc-reftemp )-alp0*( temp0-reftemp ) )
    vect( 2) =  0.0D0
    vect( 3) =  0.0D0

    vect( 4) =  a*ee*( alp*( tempc-reftemp )-alp0*( temp0-reftemp ) )
    vect( 5) =  0.0D0
    vect( 6) =  0.0D0

    vect( 7) =  0.0D0
    vect( 8) =  0.0D0
    vect( 9) =  0.0D0

    vect(10) =  0.0D0
    vect(11) =  0.0D0
    vect(12) =  0.0D0

    !--------------------------------------------------------------------

    vect(  1:3) = matmul( transt,   vect(1:3) )
    vect(  4:6) = matmul( transt,   vect(4:6) )
    vect(  7:9) = matmul( transt,   vect(7:9) )
    vect(10:12) = matmul( transt, vect(10:12) )

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine TLOAD_Beam_641
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2013/09/13)

end module
