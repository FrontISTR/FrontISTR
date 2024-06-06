!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module m_static_LIB_3d_vp

  use hecmw, only : kint, kreal
  use elementInfo

  implicit none

contains
  !--------------------------------------------------------------------
  subroutine STF_C3_vp                                  &
      (etype, nn, ecoord, gausses, stiff, tincr, &
      v, temperature)
    !--------------------------------------------------------------------

    use mMechGauss

    !--------------------------------------------------------------------

    integer(kind=kint), intent(in) :: etype                     !< element type
    integer(kind=kint), intent(in) :: nn                        !< the number of elemental nodes
    real(kind=kreal),   intent(in) :: ecoord(3, nn)             !< coordinates of elemental nodes
    type(tGaussStatus), intent(in) :: gausses(:)                !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:,:)               !< stiff matrix
    real(kind=kreal),   intent(in) :: tincr                     !< time increment
    real(kind=kreal),   intent(in), optional :: v(:, :)         !< nodal velocity
    real(kind=kreal),   intent(in), optional :: temperature(nn) !< temperature

    !--------------------------------------------------------------------

    integer(kind=kint) :: i, j
    integer(kind=kint) :: na, nb
    integer(kind=kint) :: isize, jsize
    integer(kind=kint) :: LX

    real(kind=kreal) :: MM(nn, nn), AA(nn, nn), DD(nn, nn, 3, 3), &
      trD(nn, nn), BB(nn, nn), CC(nn, nn, 3),   &
      MS(nn, nn), AS(nn, nn), CS(nn, nn, 3),    &
      MP(nn, nn, 3), AP(nn, nn, 3), CP(nn, nn)
    real(kind=kreal) :: spfunc(nn), gderiv(nn, 3)
    real(kind=kreal) :: elem(3, nn)
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: dndx(nn, 3)
    real(kind=kreal) :: tincr_inv
    real(kind=kreal) :: volume, volume_inv
    real(kind=kreal) :: mu
    real(kind=kreal) :: rho, rho_inv
    real(kind=kreal) :: vx, vy, vz
    real(kind=kreal) :: t1, t2, t3
    real(kind=kreal) :: v_dot_v
    real(kind=kreal) :: d
    real(kind=kreal) :: det, wg
    real(kind=kreal) :: tau
    real(kind=kreal), parameter :: gamma = 0.5D0

    !--------------------------------------------------------------------

    tincr_inv = 1.0D0/tincr

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)

    !--------------------------------------------------------------------

    t1 = 2.0D0*tincr_inv

    !---------------------------------------------------------------

    volume = 0.0D0

    loopVolume: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      volume = volume+wg

      !----------------------------------------------------------

    end do loopVolume

    volume_inv = 1.0D0/volume

    !---------------------------------------------------------------

    naturalCoord(1) = 0.25D0
    naturalCoord(2) = 0.25D0
    naturalCoord(3) = 0.25D0

    call getShapeFunc(etype, naturalCoord, spfunc)

    vx = 0.0D0
    vy = 0.0D0
    vz = 0.0D0

    do na = 1, nn

      vx = vx+spfunc(na)*v(1, na)
      vy = vy+spfunc(na)*v(2, na)
      vz = vz+spfunc(na)*v(3, na)

    end do

    v_dot_v = vx*vx+vy*vy+vz*vz

    !---------------------------------------------------------------

    mu  = 0.0D0
    rho = 0.0D0

    do na = 1, nn

      dndx(na, 1) = 0.0D0
      dndx(na, 2) = 0.0D0
      dndx(na, 3) = 0.0D0

    end do

    loopGlobalDeriv: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(1:3) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      mu = mu+wg*gausses(LX)%pMaterial%variables(M_VISCOCITY)

      rho = rho+wg*gausses(LX)%pMaterial%variables(M_DENSITY)

      !----------------------------------------------------------

      do na = 1, nn

        dndx(na, 1) = dndx(na, 1)+wg*gderiv(na, 1)
        dndx(na, 2) = dndx(na, 2)+wg*gderiv(na, 2)
        dndx(na, 3) = dndx(na, 3)+wg*gderiv(na, 3)

      end do

      !----------------------------------------------------------

    end do loopGlobalDeriv

    mu  = volume_inv*mu
    rho = volume_inv*rho

    do na = 1, nn

      dndx(na, 1) = volume_inv*dndx(na, 1)
      dndx(na, 2) = volume_inv*dndx(na, 2)
      dndx(na, 3) = volume_inv*dndx(na, 3)

    end do

    !---------------------------------------------------------------

    d = 0.0D0

    do na = 1, nn

      d = d+dabs( vx*dndx(na, 1)+vy*dndx(na, 2)+vz*dndx(na, 3) )

    end do

    ! h_es3d = 2.0D0/( d/DSQRT( v_dot_v ) )

    !---------------------------------------------------------------

    ! t2 = 2.0D0*DSQRT( v_dot_v )/h_es3d
    t2 = d

    !----------------------------------------------------------------

    if( v_dot_v .LT. 1.0D-15 ) then

      t3 = 4.0D0*mu/( rho*volume**(2.0D0/3.0D0) )

    else

      t3 = mu*d*d/( rho*v_dot_v )

    end if

    !----------------------------------------------------------------

    tau = 1.0D0/dsqrt( t1*t1+t2*t2+t3*t3 )

    !--------------------------------------------------------------------

    stiff(:, :) = 0.0D0

    loopGauss: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      mu = gausses(LX)%pMaterial%variables(M_VISCOCITY)

      rho = gausses(LX)%pMaterial%variables(M_DENSITY)
      rho_inv = 1.0D0/rho

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(1:3) )
      call getShapeFunc( etype, naturalCoord(:), spfunc(:) )
      call getGlobalDeriv( etype, nn, naturalCoord(:), elem(:,:), &
        det, gderiv(:,:) )

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      vx = 0.0D0
      vy = 0.0D0
      vz = 0.0D0

      do na = 1, nn

        vx = vx+spfunc(na)*v(1, na)
        vy = vy+spfunc(na)*v(2, na)
        vz = vz+spfunc(na)*v(3, na)

      end do

      !----------------------------------------------------------

      do nb = 1,nn
        do na = 1,nn
          MM(na, nb) = spfunc(na)*spfunc(nb)
          AA(na, nb) = vx*( spfunc(na)*gderiv(nb, 1) ) &
            +vy*( spfunc(na)*gderiv(nb, 2) ) &
            +vz*( spfunc(na)*gderiv(nb, 3) )
          DD(na, nb, 1, 1) = gderiv(na, 1)*gderiv(nb, 1)
          DD(na, nb, 1, 2) = gderiv(na, 1)*gderiv(nb, 2)
          DD(na, nb, 1, 3) = gderiv(na, 1)*gderiv(nb, 3)
          DD(na, nb, 2, 1) = gderiv(na, 2)*gderiv(nb, 1)
          DD(na, nb, 2, 2) = gderiv(na, 2)*gderiv(nb, 2)
          DD(na, nb, 2, 3) = gderiv(na, 2)*gderiv(nb, 3)
          DD(na, nb, 3, 1) = gderiv(na, 3)*gderiv(nb, 1)
          DD(na, nb, 3, 2) = gderiv(na, 3)*gderiv(nb, 2)
          DD(na, nb, 3, 3) = gderiv(na, 3)*gderiv(nb, 3)
          trD(na, nb) = DD(na, nb, 1, 1) &
            +DD(na, nb, 2, 2) &
            +DD(na, nb, 3, 3)
          BB(na, nb) = ( vx*vx )*DD(na, nb, 1, 1) &
            +( vx*vy )*DD(na, nb, 1, 2) &
            +( vx*vz )*DD(na, nb, 1, 3) &
            +( vy*vx )*DD(na, nb, 2, 1) &
            +( vy*vy )*DD(na, nb, 2, 2) &
            +( vy*vz )*DD(na, nb, 2, 3) &
            +( vz*vx )*DD(na, nb, 3, 1) &
            +( vz*vy )*DD(na, nb, 3, 2) &
            +( vz*vz )*DD(na, nb, 3, 3)
          CC(na, nb, 1) = gderiv(na, 1)*spfunc(nb)
          CC(na, nb, 2) = gderiv(na, 2)*spfunc(nb)
          CC(na, nb, 3) = gderiv(na, 3)*spfunc(nb)

          MS(nb, na) = AA(na, nb)
          AS(na, nb) = BB(na, nb)
          CS(na, nb, 1) = vx*DD(na, nb, 1, 1) &
            +vy*DD(na, nb, 2, 1) &
            +vz*DD(na, nb, 3, 1)
          CS(na, nb, 2) = vx*DD(na, nb, 1, 2) &
            +vy*DD(na, nb, 2, 2) &
            +vz*DD(na, nb, 3, 2)
          CS(na, nb, 3) = vx*DD(na, nb, 1, 3) &
            +vy*DD(na, nb, 2, 3) &
            +vz*DD(na, nb, 3, 3)
          MP(na, nb, 1) = spfunc(nb)*gderiv(na, 1)
          MP(na, nb, 2) = spfunc(nb)*gderiv(na, 2)
          MP(na, nb, 3) = spfunc(nb)*gderiv(na, 3)
          AP(nb, na, 1) = CS(na, nb, 1)
          AP(nb, na, 2) = CS(na, nb, 2)
          AP(nb, na, 3) = CS(na, nb, 3)
          CP(na, nb) = trD(na, nb)
        enddo
      enddo

      !----------------------------------------------------------

      do nb = 1, nn

        do na = 1, nn

          i = 1
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                              &
            = stiff(isize, jsize)                            &
            +wg                                             &
            *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
            +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
            +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )

          i = 1
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 2, 1) )

          i = 1
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 3, 1) )

          i = 1
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                     &
            = stiff(isize, jsize)                   &
            +wg                                    &
            *( -CC(na, nb, 1)+tau*CS(na, nb, 1) )

          i = 2
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 1, 2) )

          i = 2
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                              &
            = stiff(isize, jsize)                            &
            +wg                                             &
            *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
            +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
            +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )

          i = 2
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 3, 2) )

          i = 2
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                     &
            = stiff(isize, jsize)                   &
            +wg                                    &
            *( -CC(na, nb, 2)+tau*CS(na, nb, 2) )

          i = 3
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 1, 3) )

          i = 3
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 2, 3) )

          i = 3
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                               &
            = stiff(isize, jsize)                             &
            +wg                                              &
            *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) )  &
            +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )      &
            +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )

          i = 3
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                     &
            = stiff(isize, jsize)                   &
            +wg                                    &
            *( -CC(na, nb, 3)+tau*CS(na, nb, 3) )

          i = 4
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( CC(nb, na, j)               &
            +tincr_inv*tau*MP(na, nb, j) &
            +gamma*tau*AP(na, nb, j) )

          i = 4
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( CC(nb, na, j)               &
            +tincr_inv*tau*MP(na, nb, j) &
            +gamma*tau*AP(na, nb, j) )

          i = 4
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( CC(nb, na, j)               &
            +tincr_inv*tau*MP(na, nb, j) &
            +gamma*tau*AP(na, nb, j) )

          i = 4
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)            &
            = stiff(isize, jsize)          &
            +wg                           &
            *( rho_inv*tau*trD(na, nb) )

        end do

      end do

      !----------------------------------------------------------

    end do loopGauss

    !--------------------------------------------------------------------
  end subroutine STF_C3_vp
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  subroutine UPDATE_C3_vp                        &
      (etype, nn, ecoord, v, dv, gausses)
    !--------------------------------------------------------------------

    use mMechGauss

    !--------------------------------------------------------------------

    integer(kind=kint), intent(in) :: etype            !< element type
    integer(kind=kint), intent(in) :: nn               !< the number of elemental nodes
    real(kind=kreal), intent(in) :: ecoord(3, nn)      !< coordinates of elemental nodes
    real(kind=kreal), intent(in) :: v(4, nn)           !< nodal velcoity
    real(kind=kreal), intent(in) :: dv(4, nn)          !< nodal velocity increment
    type(tGaussStatus), intent(inout) :: gausses(:)    !< status of qudrature points

    !--------------------------------------------------------------------

    integer(kind=kint) :: LX

    real(kind=kreal) :: elem(3, nn)
    real(kind=kreal) :: totalvelo(4, nn)
    real(kind=kreal) :: spfunc(nn), gderiv(nn, 3)
    real(kind=kreal) :: gveloderiv(3, 3)
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: det
    real(kind=kreal) :: mu
    real(kind=kreal) :: p

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)

    totalvelo(:, :) = v(:, :)+dv(:, :)

    !--------------------------------------------------------------------

    loopMatrix: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      mu = gausses(LX)%pMaterial%variables(M_VISCOCITY)

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      ! Deformation rate tensor
      gveloderiv(1:3, 1:3) = matmul( totalvelo(1:3, 1:nn), gderiv(1:nn, 1:3) )
      gausses(LX)%strain(1) = gveloderiv(1, 1)
      gausses(LX)%strain(2) = gveloderiv(2, 2)
      gausses(LX)%strain(3) = gveloderiv(3, 3)
      gausses(LX)%strain(4) = 0.5D0*( gveloderiv(1, 2)+gveloderiv(2, 1) )
      gausses(LX)%strain(5) = 0.5D0*( gveloderiv(2, 3)+gveloderiv(3, 2) )
      gausses(LX)%strain(6) = 0.5D0*( gveloderiv(3, 1)+gveloderiv(1, 3) )

      !----------------------------------------------------------

      ! Pressure
      p = dot_product(totalvelo(4, 1:nn), spfunc(1:nn))

      ! Cauchy stress tensor
      gausses(LX)%stress(1) = -p+2.0D0*mu*gausses(LX)%strain(1)
      gausses(LX)%stress(2) = -p+2.0D0*mu*gausses(LX)%strain(2)
      gausses(LX)%stress(3) = -p+2.0D0*mu*gausses(LX)%strain(3)
      gausses(LX)%stress(4) = 2.0D0*mu*gausses(LX)%strain(4)
      gausses(LX)%stress(5) = 2.0D0*mu*gausses(LX)%strain(5)
      gausses(LX)%stress(6) = 2.0D0*mu*gausses(LX)%strain(6)

      !----------------------------------------------------------

      !set stress and strain for output
      gausses(LX)%strain_out(1:6) = gausses(LX)%strain(1:6)
      gausses(LX)%stress_out(1:6) = gausses(LX)%stress(1:6)

    end do loopMatrix

    !----------------------------------------------------------------

    !--------------------------------------------------------------------
  end subroutine UPDATE_C3_vp
  !--------------------------------------------------------------------


  !--------------------------------------------------------------------
  subroutine LOAD_C3_vp                                    &
      (etype, nn, ecoord, v, dv, r, gausses, tincr)
    !--------------------------------------------------------------------

    use mMechGauss

    !--------------------------------------------------------------------

    integer(kind=kint), intent(in)    :: etype         !< element type
    integer(kind=kint), intent(in)    :: nn            !< the number of elemental nodes
    real(kind=kreal), intent(in)      :: ecoord(3, nn) !< coordinates of elemental nodes
    real(kind=kreal), intent(in)      :: v(4, nn)      !< nodal dislplacements
    real(kind=kreal), intent(in)      :: dv(4, nn)     !< nodal velocity increment
    real(kind=kreal), intent(out)     :: r(4*nn)       !< elemental residual
    type(tGaussStatus), intent(inout) :: gausses(:)    !< status of qudrature points
    real(kind=kreal), intent(in)      :: tincr         !< time increment

    !--------------------------------------------------------------------

    integer(kind=kint) :: i, j, k
    integer(kind=kint) :: na, nb
    integer(kind=kint) :: isize, jsize
    integer(kind=kint) :: LX

    real(kind=kreal) :: elem(3, nn)
    real(kind=kreal) :: velo_new(4*nn)
    real(kind=kreal) :: stiff(4*nn, 4*nn)
    real(kind=kreal) :: b(4*nn)
    real(kind=kreal) :: MM(nn, nn), AA(nn, nn), DD(nn, nn, 3, 3), &
      trD(nn, nn), BB(nn, nn), CC(nn, nn, 3),   &
      MS(nn, nn), AS(nn, nn), CS(nn, nn, 3),    &
      MP(nn, nn, 3), AP(nn, nn, 3), CP(nn, nn)
    real(kind=kreal) :: spfunc(nn), gderiv(nn, 3)
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: dndx(nn, 3)
    real(kind=kreal) :: tincr_inv
    real(kind=kreal) :: volume, volume_inv
    real(kind=kreal) :: mu
    real(kind=kreal) :: rho, rho_inv
    real(kind=kreal) :: vx, vy, vz
    real(kind=kreal) :: t1, t2, t3
    real(kind=kreal) :: v_a_dot_v_a
    real(kind=kreal) :: d
    real(kind=kreal) :: det, wg
    real(kind=kreal) :: tau
    real(kind=kreal) :: m_v(3), a_v(3), d_v(3, 3, 3),        &
      ms_v(3), as_v(3), mp_dot_v, ap_dot_v
    real(kind=kreal) :: stiff_velo
    real(kind=kreal), parameter :: gamma = 0.5D0

    !--------------------------------------------------------------------

    tincr_inv = 1.0D0/tincr

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)

    do i = 1,4
      do na = 1,nn
        velo_new( 4*(na-1)+i ) = v(i, na)+dv(i, na)
      enddo
    enddo

    !--------------------------------------------------------------------

    t1 = 2.0D0*tincr_inv

    !---------------------------------------------------------------

    volume = 0.0D0

    loopVolume: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      volume = volume+wg

      !----------------------------------------------------------

    end do loopVolume

    volume_inv = 1.0D0/volume

    !---------------------------------------------------------------

    naturalCoord(1) = 0.25D0
    naturalCoord(2) = 0.25D0
    naturalCoord(3) = 0.25D0

    call getShapeFunc(etype, naturalCoord, spfunc)

    vx = 0.0D0
    vy = 0.0D0
    vz = 0.0D0

    do na = 1, nn

      vx = vx+spfunc(na)*v(1, na)
      vy = vy+spfunc(na)*v(2, na)
      vz = vz+spfunc(na)*v(3, na)

    end do

    v_a_dot_v_a = vx*vx+vy*vy+vz*vz

    !---------------------------------------------------------------

    do na = 1, nn

      dndx(na, 1) = 0.0D0
      dndx(na, 2) = 0.0D0
      dndx(na, 3) = 0.0D0

    end do

    mu = 0.0d0
    rho = 0.0d0

    loopGlobalDeriv: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(1:3) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      mu  = mu +wg*gausses(LX)%pMaterial%variables(M_VISCOCITY)
      rho = rho+wg*gausses(LX)%pMaterial%variables(M_DENSITY)

      !----------------------------------------------------------

      do na = 1, nn

        dndx(na, 1) = dndx(na, 1)+wg*gderiv(na, 1)
        dndx(na, 2) = dndx(na, 2)+wg*gderiv(na, 2)
        dndx(na, 3) = dndx(na, 3)+wg*gderiv(na, 3)

      end do

      !----------------------------------------------------------

    end do loopGlobalDeriv

    mu  = volume_inv*mu
    rho = volume_inv*rho

    do na = 1, nn

      dndx(na, 1) = volume_inv*dndx(na, 1)
      dndx(na, 2) = volume_inv*dndx(na, 2)
      dndx(na, 3) = volume_inv*dndx(na, 3)

    end do

    !---------------------------------------------------------------

    d = 0.0D0

    do na = 1, nn

      d = d+dabs( vx*dndx(na, 1)+vy*dndx(na, 2)+vz*dndx(na, 3) )

    end do

    ! h_es3d = 2.0D0/( d/DSQRT( v_dot_v ) )

    !---------------------------------------------------------------

    ! t2 = 2.0D0*DSQRT( v_dot_v )/h_es3d
    t2 = d

    !----------------------------------------------------------------

    if( v_a_dot_v_a .LT. 1.0D-15 ) then

      t3 = 4.0D0*mu/( rho*volume**(2.0D0/3.0D0) )

    else

      t3 = mu*d*d/( rho*v_a_dot_v_a )

    end if

    !----------------------------------------------------------------

    tau = 1.0D0/dsqrt( t1*t1+t2*t2+t3*t3 )

    !--------------------------------------------------------------------

    stiff(:, :) = 0.0D0

    loopMatrix: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      mu = gausses(LX)%pMaterial%variables(M_VISCOCITY)

      rho = gausses(LX)%pMaterial%variables(M_DENSITY)
      rho_inv = 1.0D0/rho

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      vx = 0.0D0
      vy = 0.0D0
      vz = 0.0D0

      do na = 1, nn

        vx = vx+spfunc(na)*v(1, na)
        vy = vy+spfunc(na)*v(2, na)
        vz = vz+spfunc(na)*v(3, na)

      end do

      !----------------------------------------------------------

      do nb = 1,nn
        do na = 1,nn

          MM(na, nb) = spfunc(na)*spfunc(nb)
          AA(na, nb) = vx*( spfunc(na)*gderiv(nb, 1) ) &
            +vy*( spfunc(na)*gderiv(nb, 2) ) &
            +vz*( spfunc(na)*gderiv(nb, 3) )
          DD(na, nb, 1, 1) = gderiv(na, 1)*gderiv(nb, 1)
          DD(na, nb, 1, 2) = gderiv(na, 1)*gderiv(nb, 2)
          DD(na, nb, 1, 3) = gderiv(na, 1)*gderiv(nb, 3)
          DD(na, nb, 2, 1) = gderiv(na, 2)*gderiv(nb, 1)
          DD(na, nb, 2, 2) = gderiv(na, 2)*gderiv(nb, 2)
          DD(na, nb, 2, 3) = gderiv(na, 2)*gderiv(nb, 3)
          DD(na, nb, 3, 1) = gderiv(na, 3)*gderiv(nb, 1)
          DD(na, nb, 3, 2) = gderiv(na, 3)*gderiv(nb, 2)
          DD(na, nb, 3, 3) = gderiv(na, 3)*gderiv(nb, 3)
          trD(na, nb) = DD(na, nb, 1, 1) &
            +DD(na, nb, 2, 2) &
            +DD(na, nb, 3, 3)
          BB(na, nb) = ( vx*vx )*DD(na, nb, 1, 1) &
            +( vx*vy )*DD(na, nb, 1, 2) &
            +( vx*vz )*DD(na, nb, 1, 3) &
            +( vy*vx )*DD(na, nb, 2, 1) &
            +( vy*vy )*DD(na, nb, 2, 2) &
            +( vy*vz )*DD(na, nb, 2, 3) &
            +( vz*vx )*DD(na, nb, 3, 1) &
            +( vz*vy )*DD(na, nb, 3, 2) &
            +( vz*vz )*DD(na, nb, 3, 3)
          CC(na, nb, 1) = gderiv(na, 1)*spfunc(nb)
          CC(na, nb, 2) = gderiv(na, 2)*spfunc(nb)
          CC(na, nb, 3) = gderiv(na, 3)*spfunc(nb)

          MS(nb, na) = AA(na, nb)
          AS(na, nb) = BB(na, nb)
          CS(na, nb, 1) = vx*DD(na, nb, 1, 1) &
            +vy*DD(na, nb, 2, 1) &
            +vz*DD(na, nb, 3, 1)
          CS(na, nb, 2) = vx*DD(na, nb, 1, 2) &
            +vy*DD(na, nb, 2, 2) &
            +vz*DD(na, nb, 3, 2)
          CS(na, nb, 3) = vx*DD(na, nb, 1, 3) &
            +vy*DD(na, nb, 2, 3) &
            +vz*DD(na, nb, 3, 3)
          MP(na, nb, 1) = spfunc(nb)*gderiv(na, 1)
          MP(na, nb, 2) = spfunc(nb)*gderiv(na, 2)
          MP(na, nb, 3) = spfunc(nb)*gderiv(na, 3)
          AP(nb, na, 1) = CS(na, nb, 1)
          AP(nb, na, 2) = CS(na, nb, 2)
          AP(nb, na, 3) = CS(na, nb, 3)
          CP(na, nb) = trD(na, nb)
        enddo
      enddo

      !----------------------------------------------------------

      do nb = 1, nn

        do na = 1, nn

          i = 1
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                              &
            = stiff(isize, jsize)                            &
            +wg                                             &
            *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
            +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
            +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )

          i = 1
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 2, 1) )

          i = 1
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 3, 1) )

          i = 1
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                     &
            = stiff(isize, jsize)                   &
            +wg                                    &
            *( -CC(na, nb, 1)+tau*CS(na, nb, 1) )

          i = 2
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 1, 2) )

          i = 2
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                              &
            = stiff(isize, jsize)                            &
            +wg                                             &
            *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
            +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
            +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )

          i = 2
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 3, 2) )

          i = 2
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                     &
            = stiff(isize, jsize)                   &
            +wg                                    &
            *( -CC(na, nb, 2)+tau*CS(na, nb, 2) )

          i = 3
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 1, 3) )

          i = 3
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( gamma*mu*DD(na, nb, 2, 3) )

          i = 3
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                              &
            = stiff(isize, jsize)                            &
            +wg                                             &
            *( tincr_inv*rho*( MM(na, nb)+tau*MS(na, nb) ) &
            +gamma*rho*( AA(na, nb)+tau*AS(na, nb) )     &
            +gamma*mu*( DD(na, nb, i, j)+trD(na, nb) ) )

          i = 3
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)                     &
            = stiff(isize, jsize)                   &
            +wg                                    &
            *( -CC(na, nb, 3)+tau*CS(na, nb, 3) )

          i = 4
          j = 1
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( CC(nb, na, j)               &
            +tincr_inv*tau*MP(na, nb, j) &
            +gamma*tau*AP(na, nb, j) )

          i = 4
          j = 2
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( CC(nb, na, j)               &
            +tincr_inv*tau*MP(na, nb, j) &
            +gamma*tau*AP(na, nb, j) )

          i = 4
          j = 3
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)              &
            = stiff(isize, jsize)            &
            +wg                             &
            *( CC(nb, na, j)               &
            +tincr_inv*tau*MP(na, nb, j) &
            +gamma*tau*AP(na, nb, j) )

          i = 4
          j = 4
          isize = 4*(na-1)+i
          jsize = 4*(nb-1)+j

          stiff(isize, jsize)            &
            = stiff(isize, jsize)          &
            +wg                           &
            *( rho_inv*tau*trD(na, nb) )

        end do

      end do

      !----------------------------------------------------------

    end do loopMatrix

    !--------------------------------------------------------------------

    b(:) = 0.0D0

    loopVector: do LX = 1, NumOfQuadPoints(etype)

      !----------------------------------------------------------

      mu = gausses(LX)%pMaterial%variables(M_VISCOCITY)

      rho = gausses(LX)%pMaterial%variables(M_DENSITY)
      rho_inv = 1.0D0/rho

      !----------------------------------------------------------

      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getShapeFunc(etype, naturalCoord, spfunc)
      call getGlobalDeriv(etype, nn, naturalCoord, elem, det, gderiv)

      !----------------------------------------------------------

      wg = getWeight(etype, LX)*det

      !----------------------------------------------------------

      vx = 0.0D0
      vy = 0.0D0
      vz = 0.0D0

      do na = 1, nn

        vx = vx+spfunc(na)*v(1, na)
        vy = vy+spfunc(na)*v(2, na)
        vz = vz+spfunc(na)*v(3, na)

      end do

      !----------------------------------------------------------

      do nb = 1,nn
        do na = 1,nn
          MM(na, nb) = spfunc(na)*spfunc(nb)
          AA(na, nb) = vx*( spfunc(na)*gderiv(nb, 1) ) &
            +vy*( spfunc(na)*gderiv(nb, 2) ) &
            +vz*( spfunc(na)*gderiv(nb, 3) )
          DD(na, nb, 1, 1) = gderiv(na, 1)*gderiv(nb, 1)
          DD(na, nb, 1, 2) = gderiv(na, 1)*gderiv(nb, 2)
          DD(na, nb, 1, 3) = gderiv(na, 1)*gderiv(nb, 3)
          DD(na, nb, 2, 1) = gderiv(na, 2)*gderiv(nb, 1)
          DD(na, nb, 2, 2) = gderiv(na, 2)*gderiv(nb, 2)
          DD(na, nb, 2, 3) = gderiv(na, 2)*gderiv(nb, 3)
          DD(na, nb, 3, 1) = gderiv(na, 3)*gderiv(nb, 1)
          DD(na, nb, 3, 2) = gderiv(na, 3)*gderiv(nb, 2)
          DD(na, nb, 3, 3) = gderiv(na, 3)*gderiv(nb, 3)
          BB(na, nb) = ( vx*vx )*DD(na, nb, 1, 1) &
            +( vx*vy )*DD(na, nb, 1, 2) &
            +( vx*vz )*DD(na, nb, 1, 3) &
            +( vy*vx )*DD(na, nb, 2, 1) &
            +( vy*vy )*DD(na, nb, 2, 2) &
            +( vy*vz )*DD(na, nb, 2, 3) &
            +( vz*vx )*DD(na, nb, 3, 1) &
            +( vz*vy )*DD(na, nb, 3, 2) &
            +( vz*vz )*DD(na, nb, 3, 3)

          MS(nb, na) = AA(na, nb)
          AS(na, nb) = BB(na, nb)
          CS(na, nb, 1) = vx*DD(na, nb, 1, 1) &
            +vy*DD(na, nb, 2, 1) &
            +vz*DD(na, nb, 3, 1)
          CS(na, nb, 2) = vx*DD(na, nb, 1, 2) &
            +vy*DD(na, nb, 2, 2) &
            +vz*DD(na, nb, 3, 2)
          CS(na, nb, 3) = vx*DD(na, nb, 1, 3) &
            +vy*DD(na, nb, 2, 3) &
            +vz*DD(na, nb, 3, 3)
          MP(na, nb, 1) = spfunc(nb)*gderiv(na, 1)
          MP(na, nb, 2) = spfunc(nb)*gderiv(na, 2)
          MP(na, nb, 3) = spfunc(nb)*gderiv(na, 3)
          AP(nb, na, 1) = CS(na, nb, 1)
          AP(nb, na, 2) = CS(na, nb, 2)
          AP(nb, na, 3) = CS(na, nb, 3)
        enddo
      enddo

      !----------------------------------------------------------

      do na = 1, nn

        !----------------------------------------------------

        do i = 1, 3

          m_v(i) = 0.0D0
          a_v(i) = 0.0D0
          do j = 1, 3
            do k = 1, 3
              d_v(j, k, i) = 0.0D0
            end do
          end do
          ms_v(i) = 0.0D0
          as_v(i) = 0.0D0
          mp_dot_v = 0.0D0
          ap_dot_v = 0.0D0

          do nb = 1, nn

            ! Unsteady term
            m_v(i) = m_v(i)+MM(na, nb)*v(i, nb)
            ! Advection term
            a_v(i) = a_v(i)+AA(na, nb)*v(i, nb)
            ! Diffusion term
            do j = 1, 3
              do k = 1, 3
                d_v(j, k, i) = d_v(j, k, i)+DD(na, nb, j, k)*v(i, nb)
              end do
            end do
            ! Unsteady term (SUPG)
            ms_v(i) = ms_v(i)+MS(na, nb)*v(i, nb)
            ! Advection term (SUPG)
            as_v(i) = as_v(i)+AS(na, nb)*v(i, nb)
            ! Unsteady term (PSPG)
            mp_dot_v = mp_dot_v+( MP(na, nb, 1)*v(1, nb)   &
              +MP(na, nb, 2)*v(2, nb)   &
              +MP(na, nb, 3)*v(3, nb) )
            ! Advection term (PSPG)
            ap_dot_v = ap_dot_v+( AP(na, nb, 1)*v(1, nb)   &
              +AP(na, nb, 2)*v(2, nb)   &
              +AP(na, nb, 3)*v(3, nb) )
          end do

        end do

        !----------------------------------------------------

        do i = 1, 3

          isize = 4*(na-1)+i

          b(isize)                                                  &
            = b(isize)                                                &
            +wg                                                      &
            *( tincr_inv*rho*( m_v(i)+tau*ms_v(i) )                 &
            -( 1.0D0-gamma )*rho*( a_v(i)+tau*as_v(i) )           &
            -( 1.0D0-gamma )                                      &
            *mu*( ( d_v(1, 1, i)+d_v(2, 2, i)+d_v(3, 3, i) )     &
            +( d_v(1, i, 1)+d_v(2, i, 2)+d_v(3, i, 3) ) ) )

        end do

        i = 4
        isize = 4*(na-1)+i

        b(isize)                                &
          = b(isize)                              &
          +wg                                    &
          *( tincr_inv*tau*( mp_dot_v )         &
          -( 1.0D0-gamma )*tau*( ap_dot_v ) )

        !----------------------------------------------------

      end do

      !----------------------------------------------------------

    end do loopVector

    !----------------------------------------------------------------

    do isize = 1, 4*nn

      stiff_velo = 0.0D0

      do jsize = 1, 4*nn

        stiff_velo = stiff_velo+stiff(isize, jsize)*velo_new(jsize)

      end do

      r(isize) = b(isize)-stiff_velo

    end do

    !--------------------------------------------------------------------
  end subroutine LOAD_C3_vp
  !--------------------------------------------------------------------


end module m_static_LIB_3d_vp
