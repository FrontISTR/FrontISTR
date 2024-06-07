!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module ...
module m_static_LIB_shell
  use hecmw, only : kint, kreal
  use elementInfo

  !--------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------

  !--------------------------------------------------------------
  !
  ! (Programmer)
  ! Gaku Hashimoto
  ! Department of Human and Engineered Environmental Studies
  ! Graduate School of Frontier Sciences, The University of Tokyo
  ! 5-1-5 Kashiwanoha, Kashiwa, Chiba 277-8563 JAPAN
  !
  ! (Ref.)
  ! [1] Noguchi, H. and Hisada, T.,
  !     "Sensitivity analysis in post-buckling problems of shell
  !     structures,"
  !     Computers & Structures, Vol.47, No.4, pp.699-710, (1993).
  ! [2] Dvorkin, E.N. and Bathe, K.J.,
  !     "A Continuum Mechanics Based Four-node Shell Element for
  !     General Non-linear Analysis,"
  !     Engineering Computations, Vol.1, pp.77-88, (1984).
  ! [3] Bucalem, M.L. and Bathe, K.J.,
  !     "Higher-order MITC general shell element,"
  !     International Journal for Numerical Methods in
  !     Engineering, Vol.36, pp.3729-3754, (1993).
  ! [4] Lee, P.S. and Bathe, K.J.,
  !     "Development of MITC Isotropic Triangular Shell Finite
  !     Elements,"
  !     Computers & Structures, Vol.82, pp.945-962, (2004).
  !
  ! Xi YUAN
  !   Apr. 13, 2019: Introduce mass matrix calculation
  ! (Ref.)
  ! [5] E.Hinton, T.A.Rock, O.C.Zienkiewicz(1976): A Note on Mass Lumping
  ! and Related Process in FEM. International Journal on Earthquake Eng
  ! and structural dynamics, 4, pp245-249
  !
  !--------------------------------------------------------------

  !--------------------------------------------------------------------

contains


  !####################################################################
  subroutine STF_Shell_MITC                                           &
      (etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, nddisp)
    !####################################################################

    use mMechGauss
    use gauss_integration
    use m_MatMatrix

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn, mixflag
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(out)  :: stiff(:, :)
    real(kind = kreal), intent(in)   :: thick

    real(kind = kreal), intent(in), optional :: nddisp(3, nn)

    !--------------------------------------------------------------------

    integer :: flag, flag_dof
    integer :: i, j, m
    integer :: lx, ly
    integer :: fetype
    integer :: ny
    integer :: ntying
    integer :: npoints_tying(3)
    integer :: it, ip
    integer :: na, nb
    integer :: isize, jsize
    integer :: jsize1, jsize2, jsize3, &
      jsize4, jsize5, jsize6
    integer :: n_layer,n_totlyr, sstable(24)

    real(kind = kreal) :: D(5, 5), B(5, ndof*nn), DB(5, ndof*nn)
    real(kind = kreal) :: tmpstiff(ndof*nn, ndof*nn)
    real(kind = kreal) :: elem(3, nn)
    real(kind = kreal) :: unode(3, nn)
    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: w_w_lx, w_ly
    real(kind = kreal) :: B_di(5, ndof*nn, 6, 3, 7)
    real(kind = kreal) :: B1(3, ndof*nn), B2(3, ndof*nn), &
      B3(3, ndof*nn)
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: tpcoord(6, 2, 3)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: aa1(3), aa2(3), aa3(3)
    real(kind = kreal) :: bb1(3), bb2(3), bb3(3)
    real(kind = kreal) :: cc1(3), cc2(3)
    real(kind = kreal) :: alpha
    real(kind = kreal) :: xxi_lx, eeta_lx
    real(kind = kreal) :: xxi_di(6, 3), eeta_di(6, 3)
    real(kind = kreal) :: h(nn, 3)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_i(3), v2_i(3), v3_i(3)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
      dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: g3_abs
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: cg1(3), cg2(3), cg3(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_inv
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: w_w_w_det
    real(kind = kreal) :: e1_hat(3), e2_hat(3), e3_hat(3)
    real(kind = kreal) :: e1_hat_abs, e2_hat_abs
    real(kind = kreal) :: Cv12(ndof*nn), Cv13(ndof*nn), &
      Cv21(ndof*nn), Cv23(ndof*nn), &
      Cv31(ndof*nn), Cv32(ndof*nn)
    real(kind = kreal) :: Cv_theta(ndof*nn), Cv_w(ndof*nn)
    real(kind = kreal) :: Cv(ndof*nn)
    real(kind = kreal) :: sumlyr

    sstable = 0
    flag_dof = 0
    ny = 0

    !--------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      fetype = fe_mitc4_shell

      ny = 2

      ntying = 1
      npoints_tying(1)= 4

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      fetype = fe_mitc9_shell

      ny = 3

      ntying = 3
      npoints_tying(1)= 6
      npoints_tying(2)= 6
      npoints_tying(3)= 4

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      fetype = fe_mitc3_shell

      ny = 2

      ntying = 1
      npoints_tying(1)= 3

    end if

    !--------------------------------------------------------------------

    if( present( nddisp ) ) then

      unode(:, 1:nn) = nddisp(:, :)

    end if

    !--------------------------------------------------------------------

    flag = gausses(1)%pMaterial%nlgeom_flag

    if( .not. present( nddisp ) ) flag = INFINITESIMAL

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)

    !--------------------------------------------------------------------

    tmpstiff(:, :) = 0.0D0

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    !-------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) =  0.0D0
      tpcoord(2, 1, 1) =  1.0D0
      tpcoord(3, 1, 1) =  0.0D0
      tpcoord(4, 1, 1) = -1.0D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0
      tpcoord(2, 2, 1) =  0.0D0
      tpcoord(3, 2, 1) =  1.0D0
      tpcoord(4, 2, 1) =  0.0D0

      !--------------------------------------------------------

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      !--------------------------------------------------------

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 1)  = -1.0D0
      xxi_di(2, 1)  =  1.0D0
      xxi_di(3, 1)  =  1.0D0
      xxi_di(4, 1)  = -1.0D0
      xxi_di(5, 1)  =  1.0D0
      xxi_di(6, 1)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 1) = -1.0D0
      eeta_di(2, 1) = -1.0D0
      eeta_di(3, 1) =  1.0D0
      eeta_di(4, 1) =  1.0D0
      eeta_di(5, 1) =  0.0D0
      eeta_di(6, 1) =  0.0D0

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 2)  = -1.0D0
      xxi_di(2, 2)  =  0.0D0
      xxi_di(3, 2)  =  1.0D0
      xxi_di(4, 2)  =  1.0D0
      xxi_di(5, 2)  =  0.0D0
      xxi_di(6, 2)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 2) = -1.0D0
      eeta_di(2, 2) = -1.0D0
      eeta_di(3, 2) = -1.0D0
      eeta_di(4, 2) =  1.0D0
      eeta_di(5, 2) =  1.0D0
      eeta_di(6, 2) =  1.0D0

      !--------------------------------------------------------

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1)  = 0.5D0
      tpcoord(2, 1, 1)  = 0.0D0
      tpcoord(3, 1, 1)  = 0.5D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = 0.0D0
      tpcoord(2, 2, 1) = 0.5D0
      tpcoord(3, 2, 1) = 0.5D0

      !--------------------------------------------------------

    end if

    !--------------------------------------------------------------------

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    do i = 1, 3

      g1(i) = 0.0D0

      do na = 1, nn

        g1(i) = g1(i)+shapederiv(na, 1) &
          *elem(i, na)

      end do

    end do

    e_0(1) = g1(1)
    e_0(2) = g1(2)
    e_0(3) = g1(3)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1) &
            *elem(i, na)
          g2(i) = g2(i)+shapederiv(na, 2) &
            *elem(i, na)

        end do

      end do

      !------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(1, nb) = det_cg3(1)/det_cg3_abs
      v3(2, nb) = det_cg3(2)/det_cg3_abs
      v3(3, nb) = det_cg3(3)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs .GT. 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else    ! YX: impossible

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !---------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------
    !     MODIFIED to LAMINATED SHELL ANALYSIS
    !--------------------------------------------------------------------
    zeta_ly = 0.0D0

    n_totlyr =  gausses(1)%pMaterial%totallyr
    do n_layer=1,n_totlyr
      do ly = 1, ny

        !--------------------------------------------------------

        ! MITC4
        if( etype .EQ. fe_mitc4_shell ) then

          zeta_ly = 0.0D0

          ! MITC9
        else if( etype .EQ. fe_mitc9_shell ) then

          zeta_ly = xg(ny, ly)

          ! MITC3
        else if( etype .EQ. fe_mitc3_shell )then

          zeta_ly = 0.0D0

        end if

        !--------------------------------------------------------

        do it = 1, ntying

          do ip = 1, npoints_tying(it)

            !-------------------------------------------------

            naturalcoord(1) = tpcoord(ip, 1, it)
            naturalcoord(2) = tpcoord(ip, 2, it)

            call getShapeFunc(fetype, naturalcoord, shapefunc)

            call getShapeDeriv(fetype, naturalcoord, shapederiv)

            !-------------------------------------------------

            do na = 1, nn

              do i = 1, 3

                u_rot(i, na)                      &
                  = shapefunc(na)                   &
                  *( zeta_ly*a_over_2_v3(i, na) )

                dudxi_rot(i, na)                  &
                  = shapederiv(na, 1)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudeta_rot(i, na)                 &
                  = shapederiv(na, 2)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudzeta_rot(i, na)                &
                  = shapefunc(na)                   &
                  *( a_over_2_v3(i, na) )

              end do

            end do

            !-------------------------------------------------

            ! Covariant basis vector
            do i = 1, 3

              g1(i) = 0.0D0
              g2(i) = 0.0D0
              g3(i) = 0.0D0

              do na = 1, nn

                g1(i) = g1(i)+shapederiv(na, 1)  &
                  *elem(i, na)       &
                  +dudxi_rot(i, na)
                g2(i) = g2(i)+shapederiv(na, 2)  &
                  *elem(i, na)       &
                  +dudeta_rot(i, na)
                g3(i) = g3(i)+dudzeta_rot(i, na)

              end do

            end do

            !-------------------------------------------------

            ! [ B L ] matrix
            do nb = 1, nn

              jsize1 = ndof*(nb-1)+1
              jsize2 = ndof*(nb-1)+2
              jsize3 = ndof*(nb-1)+3
              jsize4 = ndof*(nb-1)+4
              jsize5 = ndof*(nb-1)+5
              jsize6 = ndof*(nb-1)+6

              aa1(1) = dudxi_rot(2, nb)  *g1(3)-dudxi_rot(3, nb)  *g1(2)
              aa1(2) = dudxi_rot(3, nb)  *g1(1)-dudxi_rot(1, nb)  *g1(3)
              aa1(3) = dudxi_rot(1, nb)  *g1(2)-dudxi_rot(2, nb)  *g1(1)

              aa2(1) = dudxi_rot(2, nb)  *g2(3)-dudxi_rot(3, nb)  *g2(2)
              aa2(2) = dudxi_rot(3, nb)  *g2(1)-dudxi_rot(1, nb)  *g2(3)
              aa2(3) = dudxi_rot(1, nb)  *g2(2)-dudxi_rot(2, nb)  *g2(1)

              aa3(1) = dudxi_rot(2, nb)  *g3(3)-dudxi_rot(3, nb)  *g3(2)
              aa3(2) = dudxi_rot(3, nb)  *g3(1)-dudxi_rot(1, nb)  *g3(3)
              aa3(3) = dudxi_rot(1, nb)  *g3(2)-dudxi_rot(2, nb)  *g3(1)

              bb1(1) = dudeta_rot(2, nb) *g1(3)-dudeta_rot(3, nb) *g1(2)
              bb1(2) = dudeta_rot(3, nb) *g1(1)-dudeta_rot(1, nb) *g1(3)
              bb1(3) = dudeta_rot(1, nb) *g1(2)-dudeta_rot(2, nb) *g1(1)

              bb2(1) = dudeta_rot(2, nb) *g2(3)-dudeta_rot(3, nb) *g2(2)
              bb2(2) = dudeta_rot(3, nb) *g2(1)-dudeta_rot(1, nb) *g2(3)
              bb2(3) = dudeta_rot(1, nb) *g2(2)-dudeta_rot(2, nb) *g2(1)

              bb3(1) = dudeta_rot(2, nb) *g3(3)-dudeta_rot(3, nb) *g3(2)
              bb3(2) = dudeta_rot(3, nb) *g3(1)-dudeta_rot(1, nb) *g3(3)
              bb3(3) = dudeta_rot(1, nb) *g3(2)-dudeta_rot(2, nb) *g3(1)

              cc1(1) = dudzeta_rot(2, nb)*g1(3)-dudzeta_rot(3, nb)*g1(2)
              cc1(2) = dudzeta_rot(3, nb)*g1(1)-dudzeta_rot(1, nb)*g1(3)
              cc1(3) = dudzeta_rot(1, nb)*g1(2)-dudzeta_rot(2, nb)*g1(1)

              cc2(1) = dudzeta_rot(2, nb)*g2(3)-dudzeta_rot(3, nb)*g2(2)
              cc2(2) = dudzeta_rot(3, nb)*g2(1)-dudzeta_rot(1, nb)*g2(3)
              cc2(3) = dudzeta_rot(1, nb)*g2(2)-dudzeta_rot(2, nb)*g2(1)

              B_di(1, jsize1, ip, it, ly) = shapederiv(nb, 1)*g1(1)
              B_di(2, jsize1, ip, it, ly) = shapederiv(nb, 2)*g2(1)
              B_di(3, jsize1, ip, it, ly) = shapederiv(nb, 1)*g2(1) &
                +shapederiv(nb, 2)*g1(1)
              B_di(4, jsize1, ip, it, ly) = shapederiv(nb, 2)*g3(1)
              B_di(5, jsize1, ip, it, ly) = shapederiv(nb, 1)*g3(1)

              B_di(1, jsize2, ip, it, ly) = shapederiv(nb, 1)*g1(2)
              B_di(2, jsize2, ip, it, ly) = shapederiv(nb, 2)*g2(2)
              B_di(3, jsize2, ip, it, ly) = shapederiv(nb, 1)*g2(2) &
                +shapederiv(nb, 2)*g1(2)
              B_di(4, jsize2, ip, it, ly) = shapederiv(nb, 2)*g3(2)
              B_di(5, jsize2, ip, it, ly) = shapederiv(nb, 1)*g3(2)

              B_di(1, jsize3, ip, it, ly) = shapederiv(nb, 1)*g1(3)
              B_di(2, jsize3, ip, it, ly) = shapederiv(nb, 2)*g2(3)
              B_di(3, jsize3, ip, it, ly) = shapederiv(nb, 1)*g2(3) &
                +shapederiv(nb, 2)*g1(3)
              B_di(4, jsize3, ip, it, ly) = shapederiv(nb, 2)*g3(3)
              B_di(5, jsize3, ip, it, ly) = shapederiv(nb, 1)*g3(3)

              B_di(1, jsize4, ip, it, ly) = aa1(1)
              B_di(2, jsize4, ip, it, ly) = bb2(1)
              B_di(3, jsize4, ip, it, ly) = aa2(1)+bb1(1)
              B_di(4, jsize4, ip, it, ly) = bb3(1)+cc2(1)
              B_di(5, jsize4, ip, it, ly) = aa3(1)+cc1(1)

              B_di(1, jsize5, ip, it, ly) = aa1(2)
              B_di(2, jsize5, ip, it, ly) = bb2(2)
              B_di(3, jsize5, ip, it, ly) = aa2(2)+bb1(2)
              B_di(4, jsize5, ip, it, ly) = bb3(2)+cc2(2)
              B_di(5, jsize5, ip, it, ly) = aa3(2)+cc1(2)

              B_di(1, jsize6, ip, it, ly) = aa1(3)
              B_di(2, jsize6, ip, it, ly) = bb2(3)
              B_di(3, jsize6, ip, it, ly) = aa2(3)+bb1(3)
              B_di(4, jsize6, ip, it, ly) = bb3(3)+cc2(3)
              B_di(5, jsize6, ip, it, ly) = aa3(3)+cc1(3)

            end do

            !-------------------------------------------------

          end do

        end do

        !--------------------------------------------------------

        sumlyr = 0.0d0
        do m = 1, n_layer
          sumlyr = sumlyr +2 *gausses(1)%pMaterial%shell_var(m)%weight
        enddo
        zeta_ly = -1 +sumlyr -gausses(1)%pMaterial%shell_var(n_layer)%weight*(1-xg(ny, ly))

        w_ly    = wgt(ny, ly)

        !--------------------------------------------------------

        do lx = 1, NumOfQuadPoints(fetype)

          !--------------------------------------------------

          call getQuadPoint(fetype, lx, naturalcoord)

          xi_lx  = naturalcoord(1)
          eta_lx = naturalcoord(2)

          w_w_lx = getWeight(fetype, lx)

          call getShapeFunc(fetype, naturalcoord, shapefunc)

          call getShapeDeriv(fetype, naturalcoord, shapederiv)

          !--------------------------------------------------

          do i = 1, 3

            v1_i(i) = 0.0D0
            v2_i(i) = 0.0D0
            v3_i(i) = 0.0D0

            do na = 1, nn

              v1_i(i) = v1_i(i)+shapefunc(na)*v1(i, na)
              v2_i(i) = v2_i(i)+shapefunc(na)*v2(i, na)
              v3_i(i) = v3_i(i)+shapefunc(na)*v3(i, na)

            end do

          end do

          !--------------------------------------------------

          do na = 1, nn

            do i = 1, 3

              u_rot(i, na)                      &
                = shapefunc(na)                   &
                *( zeta_ly*a_over_2_v3(i, na) )

              dudxi_rot(i, na)                  &
                = shapederiv(na, 1)               &
                *( zeta_ly*a_over_2_v3(i, na) )
              dudeta_rot(i, na)                 &
                = shapederiv(na, 2)               &
                *( zeta_ly*a_over_2_v3(i, na) )
              dudzeta_rot(i, na)                &
                = shapefunc(na)                   &
                *( a_over_2_v3(i, na) )

            end do

          end do

          !--------------------------------------------------

          ! Covariant basis vector
          do i = 1, 3

            g1(i) = 0.0D0
            g2(i) = 0.0D0
            g3(i) = 0.0D0

            do na = 1, nn

              g1(i) = g1(i)+shapederiv(na, 1)  &
                *elem(i, na)       &
                +dudxi_rot(i, na)
              g2(i) = g2(i)+shapederiv(na, 2)  &
                *elem(i, na)       &
                +dudeta_rot(i, na)
              g3(i) = g3(i)+dudzeta_rot(i, na)

            end do

          end do

          !--------------------------------------------------

          ! Jacobian
          det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
            +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
            +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

          det_inv = 1.0D0/det

          !--------------------------------------------------

          ! Contravariant basis vector
          cg1(1) = det_inv                      &
            *( g2(2)*g3(3)-g2(3)*g3(2) )
          cg1(2) = det_inv                      &
            *( g2(3)*g3(1)-g2(1)*g3(3) )
          cg1(3) = det_inv                      &
            *( g2(1)*g3(2)-g2(2)*g3(1) )
          cg2(1) = det_inv                      &
            *( g3(2)*g1(3)-g3(3)*g1(2) )
          cg2(2) = det_inv                      &
            *( g3(3)*g1(1)-g3(1)*g1(3) )
          cg2(3) = det_inv                      &
            *( g3(1)*g1(2)-g3(2)*g1(1) )
          cg3(1) = det_inv                      &
            *( g1(2)*g2(3)-g1(3)*g2(2) )
          cg3(2) = det_inv                      &
            *( g1(3)*g2(1)-g1(1)*g2(3) )
          cg3(3) = det_inv                      &
            *( g1(1)*g2(2)-g1(2)*g2(1) )

          !--------------------------------------------------

          g3_abs = dsqrt( g3(1)*g3(1)   &
            +g3(2)*g3(2)   &
            +g3(3)*g3(3) )

          !--------------------------------------------------

          ! Orthonormal vectors

          e3_hat(1) = g3(1)/g3_abs
          e3_hat(2) = g3(2)/g3_abs
          e3_hat(3) = g3(3)/g3_abs

          e1_hat(1) = g2(2)*e3_hat(3) &
            -g2(3)*e3_hat(2)
          e1_hat(2) = g2(3)*e3_hat(1) &
            -g2(1)*e3_hat(3)
          e1_hat(3) = g2(1)*e3_hat(2) &
            -g2(2)*e3_hat(1)
          e1_hat_abs = dsqrt( e1_hat(1)*e1_hat(1)   &
            +e1_hat(2)*e1_hat(2)   &
            +e1_hat(3)*e1_hat(3) )
          e1_hat(1) = e1_hat(1)/e1_hat_abs
          e1_hat(2) = e1_hat(2)/e1_hat_abs
          e1_hat(3) = e1_hat(3)/e1_hat_abs

          e2_hat(1) = e3_hat(2)*e1_hat(3) &
            -e3_hat(3)*e1_hat(2)
          e2_hat(2) = e3_hat(3)*e1_hat(1) &
            -e3_hat(1)*e1_hat(3)
          e2_hat(3) = e3_hat(1)*e1_hat(2) &
            -e3_hat(2)*e1_hat(1)
          e2_hat_abs = dsqrt( e2_hat(1)*e2_hat(1)   &
            +e2_hat(2)*e2_hat(2)   &
            +e2_hat(3)*e2_hat(3) )
          e2_hat(1) = e2_hat(1)/e2_hat_abs
          e2_hat(2) = e2_hat(2)/e2_hat_abs
          e2_hat(3) = e2_hat(3)/e2_hat_abs

          !--------------------------------------------------

          call MatlMatrix_Shell                        &
            (gausses(lx), Shell, D,                 &
            e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
            alpha, n_layer)

          !--------------------------------------------------

          ! [ B L ] matrix
          do nb = 1, nn

            jsize1 = ndof*(nb-1)+1
            jsize2 = ndof*(nb-1)+2
            jsize3 = ndof*(nb-1)+3
            jsize4 = ndof*(nb-1)+4
            jsize5 = ndof*(nb-1)+5
            jsize6 = ndof*(nb-1)+6

            aa1(1) = dudxi_rot(2, nb)  *g1(3)-dudxi_rot(3, nb)  *g1(2)
            aa1(2) = dudxi_rot(3, nb)  *g1(1)-dudxi_rot(1, nb)  *g1(3)
            aa1(3) = dudxi_rot(1, nb)  *g1(2)-dudxi_rot(2, nb)  *g1(1)

            aa2(1) = dudxi_rot(2, nb)  *g2(3)-dudxi_rot(3, nb)  *g2(2)
            aa2(2) = dudxi_rot(3, nb)  *g2(1)-dudxi_rot(1, nb)  *g2(3)
            aa2(3) = dudxi_rot(1, nb)  *g2(2)-dudxi_rot(2, nb)  *g2(1)

            aa3(1) = dudxi_rot(2, nb)  *g3(3)-dudxi_rot(3, nb)  *g3(2)
            aa3(2) = dudxi_rot(3, nb)  *g3(1)-dudxi_rot(1, nb)  *g3(3)
            aa3(3) = dudxi_rot(1, nb)  *g3(2)-dudxi_rot(2, nb)  *g3(1)

            bb1(1) = dudeta_rot(2, nb) *g1(3)-dudeta_rot(3, nb) *g1(2)
            bb1(2) = dudeta_rot(3, nb) *g1(1)-dudeta_rot(1, nb) *g1(3)
            bb1(3) = dudeta_rot(1, nb) *g1(2)-dudeta_rot(2, nb) *g1(1)

            bb2(1) = dudeta_rot(2, nb) *g2(3)-dudeta_rot(3, nb) *g2(2)
            bb2(2) = dudeta_rot(3, nb) *g2(1)-dudeta_rot(1, nb) *g2(3)
            bb2(3) = dudeta_rot(1, nb) *g2(2)-dudeta_rot(2, nb) *g2(1)

            bb3(1) = dudeta_rot(2, nb) *g3(3)-dudeta_rot(3, nb) *g3(2)
            bb3(2) = dudeta_rot(3, nb) *g3(1)-dudeta_rot(1, nb) *g3(3)
            bb3(3) = dudeta_rot(1, nb) *g3(2)-dudeta_rot(2, nb) *g3(1)

            cc1(1) = dudzeta_rot(2, nb)*g1(3)-dudzeta_rot(3, nb)*g1(2)
            cc1(2) = dudzeta_rot(3, nb)*g1(1)-dudzeta_rot(1, nb)*g1(3)
            cc1(3) = dudzeta_rot(1, nb)*g1(2)-dudzeta_rot(2, nb)*g1(1)

            cc2(1) = dudzeta_rot(2, nb)*g2(3)-dudzeta_rot(3, nb)*g2(2)
            cc2(2) = dudzeta_rot(3, nb)*g2(1)-dudzeta_rot(1, nb)*g2(3)
            cc2(3) = dudzeta_rot(1, nb)*g2(2)-dudzeta_rot(2, nb)*g2(1)

            B(1, jsize1) = shapederiv(nb, 1)*g1(1)
            B(2, jsize1) = shapederiv(nb, 2)*g2(1)
            B(3, jsize1) = shapederiv(nb, 1)*g2(1) &
              +shapederiv(nb, 2)*g1(1)
            B(4, jsize1) = shapederiv(nb, 2)*g3(1)
            B(5, jsize1) = shapederiv(nb, 1)*g3(1)

            B(1, jsize2) = shapederiv(nb, 1)*g1(2)
            B(2, jsize2) = shapederiv(nb, 2)*g2(2)
            B(3, jsize2) = shapederiv(nb, 1)*g2(2) &
              +shapederiv(nb, 2)*g1(2)
            B(4, jsize2) = shapederiv(nb, 2)*g3(2)
            B(5, jsize2) = shapederiv(nb, 1)*g3(2)

            B(1, jsize3) = shapederiv(nb, 1)*g1(3)
            B(2, jsize3) = shapederiv(nb, 2)*g2(3)
            B(3, jsize3) = shapederiv(nb, 1)*g2(3) &
              +shapederiv(nb, 2)*g1(3)
            B(4, jsize3) = shapederiv(nb, 2)*g3(3)
            B(5, jsize3) = shapederiv(nb, 1)*g3(3)

            B(1, jsize4) = aa1(1)
            B(2, jsize4) = bb2(1)
            B(3, jsize4) = aa2(1)+bb1(1)
            B(4, jsize4) = bb3(1)+cc2(1)
            B(5, jsize4) = aa3(1)+cc1(1)

            B(1, jsize5) = aa1(2)
            B(2, jsize5) = bb2(2)
            B(3, jsize5) = aa2(2)+bb1(2)
            B(4, jsize5) = bb3(2)+cc2(2)
            B(5, jsize5) = aa3(2)+cc1(2)

            B(1, jsize6) = aa1(3)
            B(2, jsize6) = bb2(3)
            B(3, jsize6) = aa2(3)+bb1(3)
            B(4, jsize6) = bb3(3)+cc2(3)
            B(5, jsize6) = aa3(3)+cc1(3)

          end do

          !--------------------------------------------------

          ! MITC4
          if( etype .EQ. fe_mitc4_shell ) then

            do jsize = 1, ndof*nn

              B(4, jsize) = 0.0D0
              B(5, jsize) = 0.0D0

              ! B_as(4, jsize)
              B(4, jsize)                                       &
                = 0.5D0*( 1.0D0-xi_lx  )*B_di(4, jsize, 4, 1, ly) &
                +0.5D0*( 1.0D0+xi_lx  )*B_di(4, jsize, 2, 1, ly)
              ! B_as(5, jsize)
              B(5, jsize)                                       &
                = 0.5D0*( 1.0D0-eta_lx )*B_di(5, jsize, 1, 1, ly) &
                +0.5D0*( 1.0D0+eta_lx )*B_di(5, jsize, 3, 1, ly)

            end do

            ! MITC9
          else if( etype .EQ. fe_mitc9_shell ) then

            xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
            eeta_lx = eta_lx/dsqrt( 3.0D0/5.0D0 )

            do ip = 1, npoints_tying(1)

              h(ip, 1)                                     &
                = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )   &
                *( ( 0.5D0*eeta_di(ip, 1)*eeta_lx )        &
                *( 1.0D0+eeta_di(ip, 1)*eeta_lx )       &
                +( 1.0D0-eeta_di(ip, 1)*eeta_di(ip, 1) ) &
                *( 1.0D0-eeta_lx*eeta_lx ) )

            end do

            xxi_lx  = xi_lx /dsqrt( 3.0D0/5.0D0 )
            eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

            do ip = 1, npoints_tying(2)

              h(ip, 2)                                      &
                = ( ( 0.5D0*xxi_di(ip, 2) *xxi_lx  )          &
                *( 1.0D0+xxi_di(ip, 2) *xxi_lx  )         &
                +( 1.0D0-xxi_di(ip, 2) *xxi_di(ip, 2)  )   &
                *( 1.0D0-xxi_lx*xxi_lx ) )                &
                *( 0.5D0*( 1.0D0+eeta_di(ip, 2)*eeta_lx ) )

            end do

            xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
            eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

            do ip = 1, npoints_tying(3)

              h(ip, 3)                                      &
                = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )    &
                *( 0.5D0*( 1.0D0+eeta_di(ip, 1)*eeta_lx ) )

            end do

            do jsize = 1, ndof*nn

              B(1, jsize) = 0.0D0
              B(2, jsize) = 0.0D0
              B(3, jsize) = 0.0D0
              B(4, jsize) = 0.0D0
              B(5, jsize) = 0.0D0

              do ip = 1, npoints_tying(1)

                ! B_as(1, jsize)
                B(1, jsize)                                      &
                  = B(1, jsize)+h(ip, 1)*B_di(1, jsize, ip, 1, ly)
                ! B_as(5, jsize)
                B(5, jsize)                                      &
                  = B(5, jsize)+h(ip, 1)*B_di(5, jsize, ip, 1, ly)

              end do

              do ip = 1, npoints_tying(2)

                ! B_as(2, jsize)
                B(2, jsize)                                      &
                  = B(2, jsize)+h(ip, 2)*B_di(2, jsize, ip, 2, ly)
                ! B_as(4, jsize)
                B(4, jsize)                                      &
                  = B(4, jsize)+h(ip, 2)*B_di(4, jsize, ip, 2, ly)

              end do

              do ip = 1, npoints_tying(3)

                ! B_as(3, jsize)
                B(3, jsize)                                      &
                  = B(3, jsize)+h(ip, 3)*B_di(3, jsize, ip, 3, ly)

              end do

            end do

            ! MITC3
          else if( etype .EQ. fe_mitc3_shell ) then

            do jsize = 1, ndof*nn

              B(4, jsize) = 0.0D0
              B(5, jsize) = 0.0D0

              ! B_as(4, jsize)
              B(4, jsize)                                 &
                = ( 1.0D0-xi_lx  )*B_di(4, jsize, 2, 1, ly) &
                +xi_lx *B_di(5, jsize, 1, 1, ly)           &
                +xi_lx *( B_di(4, jsize, 3, 1, ly)         &
                -B_di(5, jsize, 3, 1, ly)  )

              ! B_as(5, jsize)
              B(5, jsize)                                 &
                = eta_lx*B_di(4, jsize, 2, 1, ly)           &
                +( 1.0D0-eta_lx )*B_di(5, jsize, 1, 1, ly) &
                -eta_lx*( B_di(4, jsize, 3, 1, ly)         &
                -B_di(5, jsize, 3, 1, ly)  )

            end do

          end if

          !--------------------------------------------------

          w_w_w_det = w_w_lx*w_ly*det

          !--------------------------------------------------

          DB(1:5, 1:ndof*nn) = matmul( D, B(1:5, 1:ndof*nn ) )

          !--------------------------------------------------

          do jsize=1,ndof*nn 
            do isize=1,ndof*nn
              tmpstiff(isize, jsize)                               &
                = tmpstiff(isize, jsize)                             &
                +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight*dot_product( B(:, isize), DB(:, jsize) )
            end do
          end do

          !--------------------------------------------------

          ! [ B_{i} ] matrix
          do nb = 1, nn

            jsize1 = ndof*(nb-1)+1
            jsize2 = ndof*(nb-1)+2
            jsize3 = ndof*(nb-1)+3
            jsize4 = ndof*(nb-1)+4
            jsize5 = ndof*(nb-1)+5
            jsize6 = ndof*(nb-1)+6

            B1(1, jsize1) =  shapederiv(nb, 1)
            B1(2, jsize1) =  0.0D0
            B1(3, jsize1) =  0.0D0
            B1(1, jsize2) =  0.0D0
            B1(2, jsize2) =  shapederiv(nb, 1)
            B1(3, jsize2) =  0.0D0
            B1(1, jsize3) =  0.0D0
            B1(2, jsize3) =  0.0D0
            B1(3, jsize3) =  shapederiv(nb, 1)
            B1(1, jsize4) =  0.0D0
            B1(2, jsize4) = -dudxi_rot(3, nb)
            B1(3, jsize4) =  dudxi_rot(2, nb)
            B1(1, jsize5) =  dudxi_rot(3, nb)
            B1(2, jsize5) =  0.0D0
            B1(3, jsize5) = -dudxi_rot(1, nb)
            B1(1, jsize6) = -dudxi_rot(2, nb)
            B1(2, jsize6) =  dudxi_rot(1, nb)
            B1(3, jsize6) =  0.0D0

            B2(1, jsize1) =  shapederiv(nb, 2)
            B2(2, jsize1) =  0.0D0
            B2(3, jsize1) =  0.0D0
            B2(1, jsize2) =  0.0D0
            B2(2, jsize2) =  shapederiv(nb, 2)
            B2(3, jsize2) =  0.0D0
            B2(1, jsize3) =  0.0D0
            B2(2, jsize3) =  0.0D0
            B2(3, jsize3) =  shapederiv(nb, 2)
            B2(1, jsize4) =  0.0D0
            B2(2, jsize4) = -dudeta_rot(3, nb)
            B2(3, jsize4) =  dudeta_rot(2, nb)
            B2(1, jsize5) =  dudeta_rot(3, nb)
            B2(2, jsize5) =  0.0D0
            B2(3, jsize5) = -dudeta_rot(1, nb)
            B2(1, jsize6) = -dudeta_rot(2, nb)
            B2(2, jsize6) =  dudeta_rot(1, nb)
            B2(3, jsize6) =  0.0D0

            B3(1, jsize1) =  0.0D0
            B3(2, jsize1) =  0.0D0
            B3(3, jsize1) =  0.0D0
            B3(1, jsize2) =  0.0D0
            B3(2, jsize2) =  0.0D0
            B3(3, jsize2) =  0.0D0
            B3(1, jsize3) =  0.0D0
            B3(2, jsize3) =  0.0D0
            B3(3, jsize3) =  0.0D0
            B3(1, jsize4) =  0.0D0
            B3(2, jsize4) = -dudzeta_rot(3, nb)
            B3(3, jsize4) =  dudzeta_rot(2, nb)
            B3(1, jsize5) =  dudzeta_rot(3, nb)
            B3(2, jsize5) =  0.0D0
            B3(3, jsize5) = -dudzeta_rot(1, nb)
            B3(1, jsize6) = -dudzeta_rot(2, nb)
            B3(2, jsize6) =  dudzeta_rot(1, nb)
            B3(3, jsize6) =  0.0D0

          end do

          !--------------------------------------------------

          ! { C_{ij} } vector
          do jsize = 1, ndof*nn

            Cv12(jsize) = ( cg1(1)*B1(2, jsize)   &
              +cg2(1)*B2(2, jsize)   &
              +cg3(1)*B3(2, jsize) ) &
              -( cg1(2)*B1(1, jsize)   &
              +cg2(2)*B2(1, jsize)   &
              +cg3(2)*B3(1, jsize) )
            Cv13(jsize) = ( cg1(1)*B1(3, jsize)   &
              +cg2(1)*B2(3, jsize)   &
              +cg3(1)*B3(3, jsize) ) &
              -( cg1(3)*B1(1, jsize)   &
              +cg2(3)*B2(1, jsize)   &
              +cg3(3)*B3(1, jsize) )
            Cv21(jsize) = ( cg1(2)*B1(1, jsize)   &
              +cg2(2)*B2(1, jsize)   &
              +cg3(2)*B3(1, jsize) ) &
              -( cg1(1)*B1(2, jsize)   &
              +cg2(1)*B2(2, jsize)   &
              +cg3(1)*B3(2, jsize) )
            Cv23(jsize) = ( cg1(2)*B1(3, jsize)   &
              +cg2(2)*B2(3, jsize)   &
              +cg3(2)*B3(3, jsize) ) &
              -( cg1(3)*B1(2, jsize)   &
              +cg2(3)*B2(2, jsize)   &
              +cg3(3)*B3(2, jsize) )
            Cv31(jsize) = ( cg1(3)*B1(1, jsize)   &
              +cg2(3)*B2(1, jsize)   &
              +cg3(3)*B3(1, jsize) ) &
              -( cg1(1)*B1(3, jsize)   &
              +cg2(1)*B2(3, jsize)   &
              +cg3(1)*B3(3, jsize) )
            Cv32(jsize) = ( cg1(3)*B1(2, jsize)   &
              +cg2(3)*B2(2, jsize)   &
              +cg3(3)*B3(2, jsize) ) &
              -( cg1(2)*B1(3, jsize)   &
              +cg2(2)*B2(3, jsize)   &
              +cg3(2)*B3(3, jsize) )

          end do

          !--------------------------------------------------

          ! { Cw } vector
          do nb = 1, nn

            do j = 1, ndof

              jsize = ndof*(nb-1)+j

              Cv_w(jsize)                       &
                = v1_i(1)*Cv12(jsize)*v2_i(2) &
                +v1_i(1)*Cv13(jsize)*v2_i(3) &
                +v1_i(2)*Cv21(jsize)*v2_i(1) &
                +v1_i(2)*Cv23(jsize)*v2_i(3) &
                +v1_i(3)*Cv31(jsize)*v2_i(1) &
                +v1_i(3)*Cv32(jsize)*v2_i(2)

            end do

          end do

          ! { Ctheta } vector
          do nb = 1, nn

            jsize1 = ndof*(nb-1)+1
            jsize2 = ndof*(nb-1)+2
            jsize3 = ndof*(nb-1)+3
            jsize4 = ndof*(nb-1)+4
            jsize5 = ndof*(nb-1)+5
            jsize6 = ndof*(nb-1)+6

            Cv_theta(jsize1) = 0.0D0
            Cv_theta(jsize2) = 0.0D0
            Cv_theta(jsize3) = 0.0D0
            Cv_theta(jsize4) = v3_i(1)*shapefunc(nb)
            Cv_theta(jsize5) = v3_i(2)*shapefunc(nb)
            Cv_theta(jsize6) = v3_i(3)*shapefunc(nb)

          end do

          ! { C } vector
          do jsize = 1, ndof*nn

            Cv(jsize) = Cv_theta(jsize)-0.5D0*Cv_w(jsize)

          end do

          !--------------------------------------------------

          ! [ K L ] matrix
          do jsize = 1, ndof*nn
            do isize = 1, ndof*nn

              tmpstiff(isize, jsize)                &
                = tmpstiff(isize, jsize)              &
                +w_w_w_det*gausses(1)%pMaterial%shell_var(n_layer)%weight*alpha*Cv(isize)*Cv(jsize)

            end do
          end do


          !--------------------------------------------------

        end do

        !--------------------------------------------------------

      end do

      !--------------------------------------------------------------

      stiff(1:nn*ndof, 1:nn*ndof) = tmpstiff(1:nn*ndof, 1:nn*ndof)

      !--------------------------------------------------------------------
    end do     !< LAMINATED SHELL ANALYSIS

    !   write(*,"(24E25.16)") stiff

    !******************** Shell-Solid mixed analysis ********************
    !		mixglaf = 0 < natural shell (6 dof)
    !		mixglaf = 1 < mixed 361 (3*2 dof)*8 nod
    !		mixglaf = 2 < mixed 351 (3*2 dof)*6 nod
    if( mixflag == 1 )then

      !       write(*,*) 'convert for shell-solid mixed analysis'
      sstable(1) = 1
      sstable(2) = 2
      sstable(3) = 3
      sstable(4) = 7
      sstable(5) = 8
      sstable(6) = 9
      sstable(7) = 13
      sstable(8) = 14
      sstable(9) = 15
      sstable(10)= 19
      sstable(11)= 20
      sstable(12)= 21
      sstable(13)= 4
      sstable(14)= 5
      sstable(15)= 6
      sstable(16)= 10
      sstable(17)= 11
      sstable(18)= 12
      sstable(19)= 16
      sstable(20)= 17
      sstable(21)= 18
      sstable(22)= 22
      sstable(23)= 23
      sstable(24)= 24

      tmpstiff(1:nn*ndof, 1:nn*ndof) = stiff(1:nn*ndof, 1:nn*ndof)

      do i = 1, nn*ndof
        do j = 1, nn*ndof
          stiff(i,j) = tmpstiff(sstable(i),sstable(j))
        enddo
      enddo

    elseif( mixflag == 2 )then
      !		write(*,*) 'convert for shell-solid mixed analysis 351'
      sstable(1) = 1
      sstable(2) = 2
      sstable(3) = 3
      sstable(4) = 7
      sstable(5) = 8
      sstable(6) = 9
      sstable(7) = 13
      sstable(8) = 14
      sstable(9) = 15
      sstable(10)= 4
      sstable(11)= 5
      sstable(12)= 6
      sstable(13)= 10
      sstable(14)= 11
      sstable(15)= 12
      sstable(16)= 16
      sstable(17)= 17
      sstable(18)= 18

      tmpstiff(1:nn*ndof, 1:nn*ndof) = stiff(1:nn*ndof, 1:nn*ndof)

      do i = 1, nn*ndof
        do j = 1, nn*ndof
          stiff(i,j) = tmpstiff(sstable(i),sstable(j))
        enddo
      enddo

    endif

    return

    !####################################################################
  end subroutine STF_Shell_MITC
  !####################################################################


  !####################################################################
  subroutine ElementStress_Shell_MITC                  &
      (etype, nn, ndof, ecoord, gausses, edisp, &
      strain, stress, thick, zeta, n_layer, n_totlyr)
    !####################################################################

    use mMechGauss
    use gauss_integration
    use m_MatMatrix

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(in)   :: edisp(6, nn)
    real(kind = kreal), intent(out)  :: strain(:,:)
    real(kind = kreal), intent(out)  :: stress(:,:)
    real(kind = kreal), intent(in)   :: thick
    real(kind = kreal), intent(in)   :: zeta

    !--------------------------------------------------------------------

    integer :: i, m
    integer :: lx
    integer :: fetype
    integer :: ntying
    integer :: npoints_tying(3)
    integer :: it, ip
    integer :: na, nb
    integer :: n_layer, n_totlyr

    real(kind = kreal) :: D(5, 5)
    real(kind = kreal) :: elem(3, nn)
    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: tpcoord(6, 2, 3)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: alpha
    real(kind = kreal) :: xxi_lx, eeta_lx
    real(kind = kreal) :: xxi_di(6, 3), eeta_di(6, 3)
    real(kind = kreal) :: h(nn, 3)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: a_over_2_theta_cross_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: theta(3, nn)
    real(kind = kreal) :: dudxi(3), dudeta(3), dudzeta(3)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
      dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: g3_abs
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: cg1(3), cg2(3), cg3(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_inv
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: e1_hat(3), e2_hat(3), e3_hat(3)
    real(kind = kreal) :: e1_hat_abs, e2_hat_abs
    real(kind = kreal) :: e11, e22, e12_2, e23_2, e31_2
    real(kind = kreal) :: e11_di(6, 3), e22_di(6, 3),     &
      e12_di_2(6, 3), e23_di_2(6, 3), &
      e31_di_2(6, 3)
    real(kind = kreal) :: E(3, 3), Ev(5)
    real(kind = kreal) :: S(3, 3), Sv(5)
    real(kind = kreal) :: sumlyr


    zeta_ly = 0.0d0

    !--------------------------------------------------------------------

    ! for lamina stress

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      fetype = fe_mitc4_shell

      ntying = 1
      npoints_tying(1)= 4

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      fetype = fe_mitc9_shell

      ntying = 3
      npoints_tying(1)= 6
      npoints_tying(2)= 6
      npoints_tying(3)= 4

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      fetype = fe_mitc3_shell

      ntying = 1
      npoints_tying(1)= 3

    end if

    !--------------------------------------------------------------------

    elem(:, :) = ecoord(:, :)

    !--------------------------------------------------------------------

    do na = 1, nn

      theta(1, na) = edisp(4, na)
      theta(2, na) = edisp(5, na)
      theta(3, na) = edisp(6, na)

    end do

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    !-------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) =  0.0D0
      tpcoord(2, 1, 1) =  1.0D0
      tpcoord(3, 1, 1) =  0.0D0
      tpcoord(4, 1, 1) = -1.0D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0
      tpcoord(2, 2, 1) =  0.0D0
      tpcoord(3, 2, 1) =  1.0D0
      tpcoord(4, 2, 1) =  0.0D0

      !--------------------------------------------------------

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 1, 1) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 1, 1) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 2, 1) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 2, 1) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 2, 1) =  0.0D0*dsqrt( 3.0D0/5.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(2, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(3, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(4, 1, 2) =  1.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(5, 1, 2) =  0.0D0*dsqrt( 3.0D0/5.0D0 )
      tpcoord(6, 1, 2) = -1.0D0*dsqrt( 3.0D0/5.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 2) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(5, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(6, 2, 2) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 1, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 1, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(2, 2, 3) = -1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(3, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )
      tpcoord(4, 2, 3) =  1.0D0*dsqrt( 1.0D0/3.0D0 )

      !--------------------------------------------------------

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 1)  = -1.0D0
      xxi_di(2, 1)  =  1.0D0
      xxi_di(3, 1)  =  1.0D0
      xxi_di(4, 1)  = -1.0D0
      xxi_di(5, 1)  =  1.0D0
      xxi_di(6, 1)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 1) = -1.0D0
      eeta_di(2, 1) = -1.0D0
      eeta_di(3, 1) =  1.0D0
      eeta_di(4, 1) =  1.0D0
      eeta_di(5, 1) =  0.0D0
      eeta_di(6, 1) =  0.0D0

      ! Xi-coordinate at a tying point in a local element
      xxi_di(1, 2)  = -1.0D0
      xxi_di(2, 2)  =  0.0D0
      xxi_di(3, 2)  =  1.0D0
      xxi_di(4, 2)  =  1.0D0
      xxi_di(5, 2)  =  0.0D0
      xxi_di(6, 2)  = -1.0D0
      ! Eta-coordinate at a tying point in a local element
      eeta_di(1, 2) = -1.0D0
      eeta_di(2, 2) = -1.0D0
      eeta_di(3, 2) = -1.0D0
      eeta_di(4, 2) =  1.0D0
      eeta_di(5, 2) =  1.0D0
      eeta_di(6, 2) =  1.0D0

      !--------------------------------------------------------

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      !--------------------------------------------------------

      ! xi-coordinate at a tying point in a local element
      tpcoord(1, 1, 1)  = 0.5D0
      tpcoord(2, 1, 1)  = 0.0D0
      tpcoord(3, 1, 1)  = 0.5D0
      ! eta-coordinate at a tying point in a local element
      tpcoord(1, 2, 1) = 0.0D0
      tpcoord(2, 2, 1) = 0.5D0
      tpcoord(3, 2, 1) = 0.5D0

      !--------------------------------------------------------

    end if

    !--------------------------------------------------------------------

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    do i = 1, 3

      g1(i) = 0.0D0

      do na = 1, nn

        g1(i) = g1(i)+shapederiv(na, 1) &
          *elem(i, na)

      end do

    end do

    e_0(1) = g1(1)
    e_0(2) = g1(2)
    e_0(3) = g1(3)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1) &
            *elem(i, na)
          g2(i) = g2(i)+shapederiv(na, 2) &
            *elem(i, na)

        end do

      end do

      !--------------------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(1, nb) = det_cg3(1)/det_cg3_abs
      v3(2, nb) = det_cg3(2)/det_cg3_abs
      v3(3, nb) = det_cg3(3)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs .GT. 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !--------------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)

      !--------------------------------------------------------

      a_over_2_theta_cross_v3(1, nb)    &
        = theta(2, nb)*a_over_2_v3(3, nb) &
        -theta(3, nb)*a_over_2_v3(2, nb)
      a_over_2_theta_cross_v3(2, nb)    &
        = theta(3, nb)*a_over_2_v3(1, nb) &
        -theta(1, nb)*a_over_2_v3(3, nb)
      a_over_2_theta_cross_v3(3, nb)    &
        = theta(1, nb)*a_over_2_v3(2, nb) &
        -theta(2, nb)*a_over_2_v3(1, nb)

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------
    !     Modified stress in laminated shell
    !--------------------------------------------------------------------
    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      zeta_ly = 0.0D0

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      zeta_ly = zeta

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      zeta_ly = 0.0D0

    end if

    !---------------------------------------------------------

    do it = 1, ntying

      do ip = 1, npoints_tying(it)

        !-------------------------------------------------

        naturalcoord(1) = tpcoord(ip, 1, it)
        naturalcoord(2) = tpcoord(ip, 2, it)

        call getShapeFunc(fetype, naturalcoord, shapefunc)

        call getShapeDeriv(fetype, naturalcoord, shapederiv)

        !-------------------------------------------------

        do na = 1, nn

          do i = 1, 3

            u_rot(i, na)                      &
              = shapefunc(na)                   &
              *( zeta_ly*a_over_2_v3(i, na) )

            dudxi_rot(i, na)                  &
              = shapederiv(na, 1)               &
              *( zeta_ly*a_over_2_v3(i, na) )
            dudeta_rot(i, na)                 &
              = shapederiv(na, 2)               &
              *( zeta_ly*a_over_2_v3(i, na) )
            dudzeta_rot(i, na)                &
              = shapefunc(na)                   &
              *( a_over_2_v3(i, na) )

          end do

        end do

        !-------------------------------------------------

        ! Covariant basis vector
        do i = 1, 3

          g1(i) = 0.0D0
          g2(i) = 0.0D0
          g3(i) = 0.0D0

          do na = 1, nn

            g1(i) = g1(i)+shapederiv(na, 1)  &
              *elem(i, na)       &
              +dudxi_rot(i, na)
            g2(i) = g2(i)+shapederiv(na, 2)  &
              *elem(i, na)       &
              +dudeta_rot(i, na)
            g3(i) = g3(i)+dudzeta_rot(i, na)

          end do

        end do

        !---------------------------------------------

        do i = 1, 3

          dudxi(i)   = 0.0D0
          dudeta(i)  = 0.0D0
          dudzeta(i) = 0.0D0

          do na = 1, nn

            dudxi(i)                                      &
              = dudxi(i)                                    &
              +shapederiv(na, 1)                           &
              *( edisp(i, na)                             &
              +zeta_ly*a_over_2_theta_cross_v3(i, na) )
            dudeta(i)                                     &
              = dudeta(i)                                   &
              +shapederiv(na, 2)                           &
              *( edisp(i, na)                             &
              +zeta_ly*a_over_2_theta_cross_v3(i, na) )
            dudzeta(i)                                    &
              = dudzeta(i)                                  &
              +shapefunc(na)                               &
              *( a_over_2_theta_cross_v3(i, na) )

          end do

        end do

        !---------------------------------------------

        ! Infinitesimal strain tensor
        e11_di(ip, it)                             &
          = 0.5D0                                    &
          *( ( g1(1)*dudxi(1) +dudxi(1) *g1(1) )   &
          +( g1(2)*dudxi(2) +dudxi(2) *g1(2) )   &
          +( g1(3)*dudxi(3) +dudxi(3) *g1(3) ) )
        e22_di(ip, it)                             &
          = 0.5D0                                    &
          *( ( g2(1)*dudeta(1)+dudeta(1)*g2(1) )   &
          +( g2(2)*dudeta(2)+dudeta(2)*g2(2) )   &
          +( g2(3)*dudeta(3)+dudeta(3)*g2(3) ) )
        e12_2                                   &
          = ( g1(1)*dudeta(1) +dudxi(1)  *g2(1) ) &
          +( g1(2)*dudeta(2) +dudxi(2)  *g2(2) ) &
          +( g1(3)*dudeta(3) +dudxi(3)  *g2(3) )
        e23_di_2(ip, it)                        &
          = ( g2(1)*dudzeta(1)+dudeta(1) *g3(1) ) &
          +( g2(2)*dudzeta(2)+dudeta(2) *g3(2) ) &
          +( g2(3)*dudzeta(3)+dudeta(3) *g3(3) )
        e31_di_2(ip, it)                        &
          = ( g3(1)*dudxi(1)  +dudzeta(1)*g1(1) ) &
          +( g3(2)*dudxi(2)  +dudzeta(2)*g1(2) ) &
          +( g3(3)*dudxi(3)  +dudzeta(3)*g1(3) )

        !-------------------------------------------------

      end do

    end do

    !--------------------------------------------------------

    sumlyr = 0.0D0
    do m = 1, n_layer
      sumlyr = sumlyr +2*gausses(1)%pMaterial%shell_var(m)%weight
    end do
    zeta_ly = -1 +sumlyr -gausses(1)%pMaterial%shell_var(n_layer)%weight*(1-zeta)

    !--------------------------------------------------------

    do lx = 1, nn

      !--------------------------------------------------

      naturalcoord(1) = nncoord(lx, 1)
      naturalcoord(2) = nncoord(lx, 2)

      xi_lx  = naturalcoord(1)
      eta_lx = naturalcoord(2)

      call getShapeFunc(fetype, naturalcoord, shapefunc)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------

      do na = 1, nn

        do i = 1, 3

          u_rot(i, na)                      &
            = shapefunc(na)                   &
            *( zeta_ly*a_over_2_v3(i, na) )

          dudxi_rot(i, na)                  &
            = shapederiv(na, 1)               &
            *( zeta_ly*a_over_2_v3(i, na) )
          dudeta_rot(i, na)                 &
            = shapederiv(na, 2)                &
            *( zeta_ly*a_over_2_v3(i, na) )
          dudzeta_rot(i, na)                &
            = shapefunc(na)                   &
            *( a_over_2_v3(i, na) )

        end do

      end do

      !--------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0
        g3(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1)  &
            *elem(i, na)       &
            +dudxi_rot(i, na)
          g2(i) = g2(i)+shapederiv(na, 2)  &
            *elem(i, na)       &
            +dudeta_rot(i, na)
          g3(i) = g3(i)+dudzeta_rot(i, na)

        end do
      end do

      !--------------------------------------------------

      ! Jacobian
      det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
        +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
        +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

      if(det == 0.0d0) then
        write(*,*)"ERROR:LIB Shell in l2009 Not Jacobian"
        stop
      endif

      det_inv = 1.0D0/det

      !--------------------------------------------------

      ! Contravariant basis vector
      cg1(1) = det_inv                      &
        *( g2(2)*g3(3)-g2(3)*g3(2) )
      cg1(2) = det_inv                      &
        *( g2(3)*g3(1)-g2(1)*g3(3) )
      cg1(3) = det_inv                      &
        *( g2(1)*g3(2)-g2(2)*g3(1) )
      cg2(1) = det_inv                      &
        *( g3(2)*g1(3)-g3(3)*g1(2) )
      cg2(2) = det_inv                      &
        *( g3(3)*g1(1)-g3(1)*g1(3) )
      cg2(3) = det_inv                      &
        *( g3(1)*g1(2)-g3(2)*g1(1) )
      cg3(1) = det_inv                      &
        *( g1(2)*g2(3)-g1(3)*g2(2) )
      cg3(2) = det_inv                      &
        *( g1(3)*g2(1)-g1(1)*g2(3) )
      cg3(3) = det_inv                      &
        *( g1(1)*g2(2)-g1(2)*g2(1) )

      !--------------------------------------------------

      g3_abs = dsqrt( g3(1)*g3(1)   &
        +g3(2)*g3(2)   &
        +g3(3)*g3(3) )

      !--------------------------------------------------

      ! Orthonormal vectors

      e3_hat(1) = g3(1)/g3_abs
      e3_hat(2) = g3(2)/g3_abs
      e3_hat(3) = g3(3)/g3_abs

      e1_hat(1) = g2(2)*e3_hat(3) &
        -g2(3)*e3_hat(2)
      e1_hat(2) = g2(3)*e3_hat(1) &
        -g2(1)*e3_hat(3)
      e1_hat(3) = g2(1)*e3_hat(2) &
        -g2(2)*e3_hat(1)
      e1_hat_abs = dsqrt( e1_hat(1)*e1_hat(1)   &
        +e1_hat(2)*e1_hat(2)   &
        +e1_hat(3)*e1_hat(3) )
      e1_hat(1) = e1_hat(1)/e1_hat_abs
      e1_hat(2) = e1_hat(2)/e1_hat_abs
      e1_hat(3) = e1_hat(3)/e1_hat_abs

      e2_hat(1) = e3_hat(2)*e1_hat(3) &
        -e3_hat(3)*e1_hat(2)
      e2_hat(2) = e3_hat(3)*e1_hat(1) &
        -e3_hat(1)*e1_hat(3)
      e2_hat(3) = e3_hat(1)*e1_hat(2) &
        -e3_hat(2)*e1_hat(1)
      e2_hat_abs = dsqrt( e2_hat(1)*e2_hat(1)   &
        +e2_hat(2)*e2_hat(2)   &
        +e2_hat(3)*e2_hat(3) )
      e2_hat(1) = e2_hat(1)/e2_hat_abs
      e2_hat(2) = e2_hat(2)/e2_hat_abs
      e2_hat(3) = e2_hat(3)/e2_hat_abs

      !--------------------------------------------------

      do i = 1, 3

        dudxi(i)   = 0.0D0
        dudeta(i)  = 0.0D0
        dudzeta(i) = 0.0D0

        do na = 1, nn

          dudxi(i)                                      &
            = dudxi(i)                                    &
            +shapederiv(na, 1)                           &
            *( edisp(i, na)                             &
            +zeta_ly*a_over_2_theta_cross_v3(i, na) )
          dudeta(i)                                     &
            = dudeta(i)                                   &
            +shapederiv(na, 2)                           &
            *( edisp(i, na)                             &
            +zeta_ly*a_over_2_theta_cross_v3(i, na) )
          dudzeta(i)                                    &
            = dudzeta(i)                                  &
            +shapefunc(na)                               &
            *( a_over_2_theta_cross_v3(i, na) )

        end do

      end do

      !--------------------------------------------------

      ! Infinitesimal strain tensor
      e11                                        &
        = 0.5D0                                    &
        *( ( g1(1)*dudxi(1) +dudxi(1) *g1(1) )   &
        +( g1(2)*dudxi(2) +dudxi(2) *g1(2) )   &
        +( g1(3)*dudxi(3) +dudxi(3) *g1(3) ) )
      e22                                        &
        = 0.5D0                                    &
        *( ( g2(1)*dudeta(1)+dudeta(1)*g2(1) )   &
        +( g2(2)*dudeta(2)+dudeta(2)*g2(2) )   &
        +( g2(3)*dudeta(3)+dudeta(3)*g2(3) ) )
      e12_2                                   &
        = ( g1(1)*dudeta(1) +dudxi(1)  *g2(1) ) &
        +( g1(2)*dudeta(2) +dudxi(2)  *g2(2) ) &
        +( g1(3)*dudeta(3) +dudxi(3)  *g2(3) )
      e23_2                                   &
        = ( g2(1)*dudzeta(1)+dudeta(1) *g3(1) ) &
        +( g2(2)*dudzeta(2)+dudeta(2) *g3(2) ) &
        +( g2(3)*dudzeta(3)+dudeta(3) *g3(3) )
      e31_2                                   &
        = ( g3(1)*dudxi(1)  +dudzeta(1)*g1(1) ) &
        +( g3(2)*dudxi(2)  +dudzeta(2)*g1(2) ) &
        +( g3(3)*dudxi(3)  +dudzeta(3)*g1(3) )

      ! MITC4
      if( etype .EQ. fe_mitc4_shell ) then

        e23_2 = 0.0D0
        e31_2 = 0.0D0

        ! e23_as_2
        e23_2                                   &
          = 0.5D0*( 1.0D0-xi_lx  )*e23_di_2(4, 1) &
          +0.5D0*( 1.0D0+xi_lx  )*e23_di_2(2, 1)
        ! e31_as_2
        e31_2                                   &
          = 0.5D0*( 1.0D0-eta_lx )*e31_di_2(1, 1) &
          +0.5D0*( 1.0D0+eta_lx )*e31_di_2(3, 1)

        ! MITC9
      else if( etype .EQ. fe_mitc9_shell ) then

        xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
        eeta_lx = eta_lx/dsqrt( 3.0D0/5.0D0 )

        do ip = 1, npoints_tying(1)

          h(ip, 1)                                     &
            = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )   &
            *( ( 0.5D0*eeta_di(ip, 1)*eeta_lx )        &
            *( 1.0D0+eeta_di(ip, 1)*eeta_lx )       &
            +( 1.0D0-eeta_di(ip, 1)*eeta_di(ip, 1) ) &
            *( 1.0D0-eeta_lx*eeta_lx ) )

        end do

        xxi_lx  = xi_lx /dsqrt( 3.0D0/5.0D0 )
        eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

        do ip = 1, npoints_tying(2)

          h(ip, 2)                                      &
            = ( ( 0.5D0*xxi_di(ip, 2) *xxi_lx  )          &
            *( 1.0D0+xxi_di(ip, 2) *xxi_lx  )         &
            +( 1.0D0-xxi_di(ip, 2) *xxi_di(ip, 2)  )   &
            *( 1.0D0-xxi_lx*xxi_lx ) )                &
            *( 0.5D0*( 1.0D0+eeta_di(ip, 2)*eeta_lx ) )

        end do

        xxi_lx  = xi_lx /dsqrt( 1.0D0/3.0D0 )
        eeta_lx = eta_lx/dsqrt( 1.0D0/3.0D0 )

        do ip = 1, npoints_tying(3)

          h(ip, 3)                                      &
            = ( 0.5D0*( 1.0D0+xxi_di(ip, 1)*xxi_lx ) )    &
            *( 0.5D0*( 1.0D0+eeta_di(ip, 1)*eeta_lx ) )

        end do

        e11   = 0.0D0
        e31_2 = 0.0D0

        ! e11_as, e31_as_2
        do ip = 1, npoints_tying(1)

          e11   = e11  +h(ip, 1)*e11_di(ip, 1)
          e31_2 = e31_2+h(ip, 1)*e31_di_2(ip, 1)

        end do

        e22   = 0.0D0
        e23_2 = 0.0D0

        ! e22_as, e23_as_2
        do ip = 1, npoints_tying(2)

          e22   = e22  +h(ip, 2)*e22_di(ip, 2)
          e23_2 = e23_2+h(ip, 2)*e23_di_2(ip, 2)

        end do

        e12_2 = 0.0D0

        ! e12_as_2
        do ip = 1, npoints_tying(3)

          e12_2 = e12_2+h(ip, 3)*e12_di_2(ip, 3)

        end do

        ! MITC3
      else if( etype .EQ. fe_mitc3_shell ) then

        e23_2 = 0.0D0
        e31_2 = 0.0D0

        ! e23_as_2
        e23_2                                       &
          = ( 1.0D0-xi_lx )*e23_di_2(2, 1)           &
          +xi_lx *e31_di_2(1, 1)                    &
          +xi_lx *( e23_di_2(3, 1)-e31_di_2(3, 1) )
        ! e31_as_2
        e31_2                                       &
          = eta_lx*e23_di_2(2, 1)                    &
          +( 1.0D0-eta_lx )*e31_di_2(1, 1)          &
          -eta_lx*( e23_di_2(3, 1)-e31_di_2(3, 1) )

      end if

      !--------------------------------------------------

      ! { E } vector
      Ev(1) = e11
      Ev(2) = e22
      Ev(3) = e12_2
      Ev(4) = e23_2
      Ev(5) = e31_2

      ! Infinitesimal strain tensor
      ! [ E ] matrix
      E(1, 1) = Ev(1)
      E(2, 2) = Ev(2)
      E(3, 3) = 0.0D0
      E(1, 2) = 0.5D0*Ev(3)
      E(2, 1) = 0.5D0*Ev(3)
      E(2, 3) = 0.5D0*Ev(4)
      E(3, 2) = 0.5D0*Ev(4)
      E(3, 1) = 0.5D0*Ev(5)
      E(1, 3) = 0.5D0*Ev(5)

      !--------------------------------------------------
      !   write(*,*) 'Stress_n_layer', n_layer
      call MatlMatrix_Shell                        &
        (gausses(lx), Shell, D,                 &
        e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
        alpha, n_layer)

      !--------------------------------------------------

      Sv = matmul( D, Ev )

      ! Infinitesimal stress tensor
      ! [ S ] matrix
      S(1, 1) = Sv(1)
      S(2, 2) = Sv(2)
      S(3, 3) = 0.0D0
      S(1, 2) = Sv(3)
      S(2, 1) = Sv(3)
      S(2, 3) = Sv(4)
      S(3, 2) = Sv(4)
      S(3, 1) = Sv(5)
      S(1, 3) = Sv(5)

      !------------------------------------------------

      stress(lx,1)                  &
        = ( S(1, 1)*g1(1)*g1(1)   &
        +S(1, 2)*g1(1)*g2(1)   &
        +S(1, 3)*g1(1)*g3(1)   &
        +S(2, 1)*g2(1)*g1(1)   &
        +S(2, 2)*g2(1)*g2(1)   &
        +S(2, 3)*g2(1)*g3(1)   &
        +S(3, 1)*g3(1)*g1(1)   &
        +S(3, 2)*g3(1)*g2(1) )
      stress(lx,2)                 &
        = ( S(1, 1)*g1(2)*g1(2)   &
        +S(1, 2)*g1(2)*g2(2)   &
        +S(1, 3)*g1(2)*g3(2)   &
        +S(2, 1)*g2(2)*g1(2)   &
        +S(2, 2)*g2(2)*g2(2)   &
        +S(2, 3)*g2(2)*g3(2)   &
        +S(3, 1)*g3(2)*g1(2)   &
        +S(3, 2)*g3(2)*g2(2) )
      stress(lx,3)                 &
        = ( S(1, 1)*g1(3)*g1(3)   &
        +S(1, 2)*g1(3)*g2(3)   &
        +S(1, 3)*g1(3)*g3(3)   &
        +S(2, 1)*g2(3)*g1(3)   &
        +S(2, 2)*g2(3)*g2(3)   &
        +S(2, 3)*g2(3)*g3(3)   &
        +S(3, 1)*g3(3)*g1(3)   &
        +S(3, 2)*g3(3)*g2(3) )
      stress(lx,4)                 &
        = ( S(1, 1)*g1(1)*g1(2)   &
        +S(1, 2)*g1(1)*g2(2)   &
        +S(1, 3)*g1(1)*g3(2)   &
        +S(2, 1)*g2(1)*g1(2)   &
        +S(2, 2)*g2(1)*g2(2)   &
        +S(2, 3)*g2(1)*g3(2)   &
        +S(3, 1)*g3(1)*g1(2)   &
        +S(3, 2)*g3(1)*g2(2) )
      stress(lx,5)                 &
        = ( S(1, 1)*g1(2)*g1(3)   &
        +S(1, 2)*g1(2)*g2(3)   &
        +S(1, 3)*g1(2)*g3(3)   &
        +S(2, 1)*g2(2)*g1(3)   &
        +S(2, 2)*g2(2)*g2(3)   &
        +S(2, 3)*g2(2)*g3(3)   &
        +S(3, 1)*g3(2)*g1(3)   &
        +S(3, 2)*g3(2)*g2(3) )
      stress(lx,6)                 &
        = ( S(1, 1)*g1(3)*g1(1)   &
        +S(1, 2)*g1(3)*g2(1)   &
        +S(1, 3)*g1(3)*g3(1)   &
        +S(2, 1)*g2(3)*g1(1)   &
        +S(2, 2)*g2(3)*g2(1)   &
        +S(2, 3)*g2(3)*g3(1)   &
        +S(3, 1)*g3(3)*g1(1)   &
        +S(3, 2)*g3(3)*g2(1) )

      strain(lx,1)                   &
        = ( E(1, 1)*cg1(1)*cg1(1)   &
        +E(1, 2)*cg1(1)*cg2(1)   &
        +E(1, 3)*cg1(1)*cg3(1)   &
        +E(2, 1)*cg2(1)*cg1(1)   &
        +E(2, 2)*cg2(1)*cg2(1)   &
        +E(2, 3)*cg2(1)*cg3(1)   &
        +E(3, 1)*cg3(1)*cg1(1)   &
        +E(3, 2)*cg3(1)*cg2(1) )
      strain(lx,2)                   &
        = ( E(1, 1)*cg1(2)*cg1(2)   &
        +E(1, 2)*cg1(2)*cg2(2)   &
        +E(1, 3)*cg1(2)*cg3(2)   &
        +E(2, 1)*cg2(2)*cg1(2)   &
        +E(2, 2)*cg2(2)*cg2(2)   &
        +E(2, 3)*cg2(2)*cg3(2)   &
        +E(3, 1)*cg3(2)*cg1(2)   &
        +E(3, 2)*cg3(2)*cg2(2) )
      strain(lx,3)                   &
        = ( E(1, 1)*cg1(3)*cg1(3)   &
        +E(1, 2)*cg1(3)*cg2(3)   &
        +E(1, 3)*cg1(3)*cg3(3)   &
        +E(2, 1)*cg2(3)*cg1(3)   &
        +E(2, 2)*cg2(3)*cg2(3)   &
        +E(2, 3)*cg2(3)*cg3(3)   &
        +E(3, 1)*cg3(3)*cg1(3)   &
        +E(3, 2)*cg3(3)*cg2(3) )
      strain(lx,4)                   &
        = ( E(1, 1)*cg1(1)*cg1(2)   &
        +E(1, 2)*cg1(1)*cg2(2)   &
        +E(1, 3)*cg1(1)*cg3(2)   &
        +E(2, 1)*cg2(1)*cg1(2)   &
        +E(2, 2)*cg2(1)*cg2(2)   &
        +E(2, 3)*cg2(1)*cg3(2)   &
        +E(3, 1)*cg3(1)*cg1(2)   &
        +E(3, 2)*cg3(1)*cg2(2) )
      strain(lx,5)                   &
        = ( E(1, 1)*cg1(2)*cg1(3)   &
        +E(1, 2)*cg1(2)*cg2(3)   &
        +E(1, 3)*cg1(2)*cg3(3)   &
        +E(2, 1)*cg2(2)*cg1(3)   &
        +E(2, 2)*cg2(2)*cg2(3)   &
        +E(2, 3)*cg2(2)*cg3(3)   &
        +E(3, 1)*cg3(2)*cg1(3)   &
        +E(3, 2)*cg3(2)*cg2(3) )
      strain(lx,6)                   &
        = ( E(1, 1)*cg1(3)*cg1(1)   &
        +E(1, 2)*cg1(3)*cg2(1)   &
        +E(1, 3)*cg1(3)*cg3(1)   &
        +E(2, 1)*cg2(3)*cg1(1)   &
        +E(2, 2)*cg2(3)*cg2(1)   &
        +E(2, 3)*cg2(3)*cg3(1)   &
        +E(3, 1)*cg3(3)*cg1(1)   &
        +E(3, 2)*cg3(3)*cg2(1) )

      !--------------------------------------------------

    end do
    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine ElementStress_Shell_MITC
  !####################################################################


  !####################################################################
  subroutine DL_Shell                                  &
      (etype, nn, ndof, xx, yy, zz, rho, thick, &
      ltype, params, vect, nsize, gausses)
    !####################################################################

    use hecmw
    use m_utilities
    use mMechGauss
    use gauss_integration

    type(tGaussStatus), intent(in)   :: gausses(:)
    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in) :: xx(*), yy(*), zz(*)
    real(kind = kreal), intent(in) :: rho
    real(kind = kreal), intent(in) :: thick
    real(kind = kreal), intent(in) :: params(*)
    real(kind = kreal), intent(out) :: vect(*)
    integer(kind = kint), intent(out) :: nsize

    !--------------------------------------------------------------------

    integer :: ivol, isurf
    integer :: lx, ly
    integer :: fetype
    integer :: ny
    integer :: i, m
    integer :: na, nb
    integer :: isize
    integer :: jsize1, jsize2, jsize3, &
      jsize4, jsize5, jsize6
    integer :: ltype
    integer :: n_totlyr, n_layer

    real(kind = kreal) :: elem(3, nn)
    real(kind = kreal) :: val
    real(kind = kreal) :: ax, ay, az
    real(kind = kreal) :: rx, ry, rz
    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: w_w_lx, w_ly
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), &
      dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: g1_cross_g2(3)
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: w_w_w_det
    real(kind = kreal) :: N(3, ndof*nn)
    real(kind = kreal) :: hx, hy, hz
    real(kind = kreal) :: phx, phy, phz
    real(kind = kreal) :: coefx, coefy, coefz
    real(kind = kreal) :: x, y, z
    real(kind = kreal) :: sumlyr

    ny = 0

    !--------------------------------------------------------------------

    !   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
    !   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
    !   BZ   LTYPE=3  :BODY FORCE IN Z-DIRECTION
    !   CRAV LTYPE=4  :GRAVITY FORCE
    !   CENT LTYPE=5  :CENTRIFUGAL LOAD
    !   P    LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR SHELL SURFACE

    !--------------------------------------------------------------------

    nsize = ndof*nn

    !--------------------------------------------------------------------

    val = params(1)
    ax  = params(2)
    ay  = params(3)
    az  = params(4)
    rx  = params(5)
    ry  = params(6)
    rz  = params(7)

    !--------------------------------------------------------------------

    ! MITC4
    if( etype .EQ. fe_mitc4_shell ) then

      fetype = fe_mitc4_shell

      ny = 2

      ! MITC9
    else if( etype .EQ. fe_mitc9_shell ) then

      fetype = fe_mitc9_shell

      ny = 3

      ! MITC3
    else if( etype .EQ. fe_mitc3_shell ) then

      fetype = fe_mitc3_shell

      ny = 2

    end if

    !--------------------------------------------------------------------

    do na = 1, nn

      elem(1, na) = xx(na)
      elem(2, na) = yy(na)
      elem(3, na) = zz(na)

    end do

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    !--------------------------------------------------------------------

    ! Local load vector
    do isize = 1, ndof*nn

      vect(isize) = 0.0D0

    end do

    !--------------------------------------------------------------------

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    do i = 1, 3

      g1(i) = 0.0D0

      do na = 1, nn

        g1(i) = g1(i)+shapederiv(na, 1) &
          *elem(i, na)

      end do

    end do

    e_0(1) = g1(1)
    e_0(2) = g1(2)
    e_0(3) = g1(3)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)

      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      do i = 1, 3

        g1(i) = 0.0D0
        g2(i) = 0.0D0

        do na = 1, nn

          g1(i) = g1(i)+shapederiv(na, 1) &
            *elem(i, na)
          g2(i) = g2(i)+shapederiv(na, 2) &
            *elem(i, na)

        end do

      end do

      !--------------------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(1, nb) = det_cg3(1)/det_cg3_abs
      v3(2, nb) = det_cg3(2)/det_cg3_abs
      v3(3, nb) = det_cg3(3)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs .GT. 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !--------------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------

    ! Selection of load type

    ivol = 0
    isurf = 0

    if( ltype .LT. 10 ) then

      ivol = 1

    else if( ltype .GE. 10 ) then

      isurf = 1

    end if

    !--------------------------------------------------------------------

    !** Surface load
    if( isurf .EQ. 1 ) then

      !--------------------------------------------------------

      do lx = 1, NumOfQuadPoints(fetype)

        !--------------------------------------------------

        call getQuadPoint(fetype, lx, naturalcoord)

        xi_lx  = naturalcoord(1)
        eta_lx = naturalcoord(2)

        w_w_lx = getWeight(fetype, lx)

        call getShapeFunc(fetype, naturalcoord, shapefunc)

        call getShapeDeriv(fetype, naturalcoord, shapederiv)

        !--------------------------------------------------

        do na = 1, nn

          do i = 1, 3

            u_rot(i, na)                    &
              = shapefunc(na)                 &
              *( 0.0D0*a_over_2_v3(i, na) )

            dudxi_rot(i, na)                &
              = shapederiv(na, 1)             &
              *( 0.0D0*a_over_2_v3(i, na) )
            dudeta_rot(i, na)               &
              = shapederiv(na, 2)             &
              *( 0.0D0*a_over_2_v3(i, na) )
            dudzeta_rot(i, na)              &
              = shapefunc(na)                 &
              *( a_over_2_v3(i, na) )

          end do

        end do

        !--------------------------------------------------

        ! Covariant basis vector
        do i = 1, 3

          g1(i) = 0.0D0
          g2(i) = 0.0D0
          !g3(i) = 0.0D0

          do na = 1, nn

            g1(i) = g1(i)+shapederiv(na, 1)  &
              *elem(i, na)       &
              +dudxi_rot(i, na)
            g2(i) = g2(i)+shapederiv(na, 2)  &
              *elem(i, na)       &
              +dudeta_rot(i, na)
            !g3(i) = g3(i)+dudzeta_rot(i, na)

          end do

        end do

        !--------------------------------------------------

        !g3_abs = DSQRT( g3(1)*g3(1)   &
          !               +g3(2)*g3(2)   &
          !               +g3(3)*g3(3) )

        !--------------------------------------------------

        !e3_hat(1) = g3(1)/g3_abs
        !e3_hat(2) = g3(2)/g3_abs
        !e3_hat(3) = g3(3)/g3_abs

        !--------------------------------------------------

        ! Jacobian
        !det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
          !     +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
          !     +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

        !--------------------------------------------------

        g1_cross_g2(1) = g1(2)*g2(3)-g1(3)*g2(2)
        g1_cross_g2(2) = g1(3)*g2(1)-g1(1)*g2(3)
        g1_cross_g2(3) = g1(1)*g2(2)-g1(2)*g2(1)

        !--------------------------------------------------

        do nb = 1, nn

          jsize1 = ndof*(nb-1)+1
          jsize2 = ndof*(nb-1)+2
          jsize3 = ndof*(nb-1)+3
          jsize4 = ndof*(nb-1)+4
          jsize5 = ndof*(nb-1)+5
          jsize6 = ndof*(nb-1)+6

          N(1, jsize1) = shapefunc(nb)
          N(1, jsize2) = 0.0D0
          N(1, jsize3) = 0.0D0
          N(1, jsize4) = 0.0D0
          N(1, jsize5) = 0.0D0
          N(1, jsize6) = 0.0D0
          N(2, jsize1) = 0.0D0
          N(2, jsize2) = shapefunc(nb)
          N(2, jsize3) = 0.0D0
          N(2, jsize4) = 0.0D0
          N(2, jsize5) = 0.0D0
          N(2, jsize6) = 0.0D0
          N(3, jsize1) = 0.0D0
          N(3, jsize2) = 0.0D0
          N(3, jsize3) = shapefunc(nb)
          N(3, jsize4) = 0.0D0
          N(3, jsize5) = 0.0D0
          N(3, jsize6) = 0.0D0

        end do

        do isize = 1, ndof*nn

          vect(isize)                                    &
            = vect(isize)                                  &
            +w_w_lx*( N(1, isize)*g1_cross_g2(1)       &
            +N(2, isize)*g1_cross_g2(2)       &
            +N(3, isize)*g1_cross_g2(3) )*val

        end do

        !--------------------------------------------------

      end do

      !--------------------------------------------------------

    end if

    !--------------------------------------------------------------------

    !** Volume load
    if( ivol .EQ. 1 ) then

      !--------------------------------------------------------
      n_totlyr =  gausses(1)%pMaterial%totallyr
      do n_layer=1,n_totlyr
        do ly = 1, ny

          !--------------------------------------------------


          sumlyr = 0.0D0
          do m = 1, n_layer
            sumlyr = sumlyr + 2 * gausses(1)%pMaterial%shell_var(m)%weight
          end do
          zeta_ly = -1 + sumlyr - gausses(1)%pMaterial%shell_var(n_layer)%weight*(1-xg(ny, ly))

          !zeta_ly = xg(ny, ly)
          w_ly    = wgt(ny, ly)

          !--------------------------------------------------

          do lx = 1, NumOfQuadPoints(fetype)

            !--------------------------------------------

            call getQuadPoint(fetype, lx, naturalcoord)

            xi_lx  = naturalcoord(1)
            eta_lx = naturalcoord(2)

            w_w_lx = getWeight(fetype, lx)

            call getShapeFunc(fetype, naturalcoord, shapefunc)

            call getShapeDeriv(fetype, naturalcoord, shapederiv)

            !--------------------------------------------

            do na = 1, nn

              do i = 1, 3

                u_rot(i, na)                      &
                  = shapefunc(na)                   &
                  *( zeta_ly*a_over_2_v3(i, na) )

                dudxi_rot(i, na)                  &
                  = shapederiv(na, 1)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudeta_rot(i, na)                 &
                  = shapederiv(na, 2)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudzeta_rot(i, na)                &
                  = shapefunc(na)                   &
                  *( a_over_2_v3(i, na) )

              end do

            end do

            !--------------------------------------------

            ! Covariant basis vector
            do i = 1, 3

              g1(i) = 0.0D0
              g2(i) = 0.0D0
              g3(i) = 0.0D0

              do na = 1, nn

                g1(i) = g1(i)+shapederiv(na, 1)  &
                  *elem(i, na)       &
                  +dudxi_rot(i, na)
                g2(i) = g2(i)+shapederiv(na, 2)  &
                  *elem(i, na)       &
                  +dudeta_rot(i, na)
                g3(i) = g3(i)+dudzeta_rot(i, na)

              end do

            end do

            !--------------------------------------------

            ! Jacobian
            det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
              +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
              +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

            !--------------------------------------------

            ! [ N ] matrix
            do nb = 1, nn

              jsize1 = ndof*(nb-1)+1
              jsize2 = ndof*(nb-1)+2
              jsize3 = ndof*(nb-1)+3
              jsize4 = ndof*(nb-1)+4
              jsize5 = ndof*(nb-1)+5
              jsize6 = ndof*(nb-1)+6

              N(1, jsize1) =  shapefunc(nb)
              N(2, jsize1) =  0.0D0
              N(3, jsize1) =  0.0D0
              N(1, jsize2) =  0.0D0
              N(2, jsize2) =  shapefunc(nb)
              N(3, jsize2) =  0.0D0
              N(1, jsize3) =  0.0D0
              N(2, jsize3) =  0.0D0
              N(3, jsize3) =  shapefunc(nb)
              N(1, jsize4) =  0.0D0
              N(2, jsize4) = -u_rot(3, nb)
              N(3, jsize4) =  u_rot(2, nb)
              N(1, jsize5) =  u_rot(3, nb)
              N(2, jsize5) =  0.0D0
              N(3, jsize5) = -u_rot(1, nb)
              N(1, jsize6) = -u_rot(2, nb)
              N(2, jsize6) =  u_rot(1, nb)
              N(3, jsize6) =  0.0D0

            enddo

            !--------------------------------------------

            w_w_w_det = w_w_lx*w_ly*det

            !--------------------------------------------

            if( ltype .EQ. 1 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*N(1, isize)*val

              end do

            else if( ltype .EQ. 2 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*N(2, isize)*val

              end do

            else if( ltype .EQ. 3 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*N(3, isize)*val

              end do

            else if( ltype .EQ. 4 ) then

              do isize = 1, ndof*nn

                vect(isize) = vect(isize)+w_w_w_det*rho*ax*N(1, isize)*val
                vect(isize) = vect(isize)+w_w_w_det*rho*ay*N(2, isize)*val
                vect(isize) = vect(isize)+w_w_w_det*rho*az*N(3, isize)*val

              end do

            else if( ltype .EQ. 5 ) then

              x = 0.0D0
              y = 0.0D0
              z = 0.0D0

              do nb = 1, nn

                x = x+shapefunc(nb)*elem(1, nb)
                y = y+shapefunc(nb)*elem(2, nb)
                z = z+shapefunc(nb)*elem(3, nb)

              end do

              hx = ax+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*rx
              hy = ay+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*ry
              hz = az+( (x-ax)*rx+(y-ay)*ry+(z-az)*rz )/( rx**2+ry**2+rz**2 )*rz

              phx = x-hx
              phy = y-hy
              phz = z-hz

              coefx = phx*val*rho*val
              coefy = phy*val*rho*val
              coefz = phz*val*rho*val

              do isize = 1, ndof*nn

                vect(isize)                       &
                  = vect(isize)                     &
                  +w_w_w_det*( N(1, isize)*coefx   &
                  +N(2, isize)*coefy   &
                  +N(3, isize)*coefz )

              end do

            end if

            !--------------------------------------------

          end do

          !----------------------------------------------

        end do

        !----------------------------------------------

      end do

      !--------------------------------------------------------

    end if
    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine DL_Shell
  !####################################################################


  !####################################################################
  subroutine DL_Shell_33                               &
      (ic_type, nn, ndof, xx, yy, zz, rho, thick, &
      ltype, params, vect, nsize, gausses)
    !####################################################################

    use hecmw
    use m_utilities
    use mMechGauss
    use gauss_integration

    type(tGaussStatus) :: gausses(:)
    !--------------------------------------------------------------------

    integer(kind = kint) :: ic_type
    integer(kind = kint) :: nn
    integer(kind = kint) :: ndof
    real(kind = kreal) :: xx(*), yy(*), zz(*)
    real(kind = kreal) :: rho
    real(kind = kreal) :: thick
    real(kind = kreal) :: params(*)
    real(kind = kreal) :: vect(*)
    integer(kind = kint) :: nsize
    integer :: ltype, i
    real(kind = kreal) :: tmp(24)
    !--------------------------------------------------------------------

    if(ic_type == 761)then
      !ic_type = 731
      !nn = 3
      !ndof = 6
      call DL_Shell(731, 3, 6, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
      !ic_type = 761
      !nn = 6
      !ndof = 3

      tmp = 0.0
      do i=1,18
        tmp(i) = vect(i)
      enddo

      vect( 1) = tmp(1)
      vect( 2) = tmp(2)
      vect( 3) = tmp(3)
      vect( 4) = tmp(7)
      vect( 5) = tmp(8)
      vect( 6) = tmp(9)
      vect( 7) = tmp(13)
      vect( 8) = tmp(14)
      vect( 9) = tmp(15)
      vect(10) = tmp(4)
      vect(11) = tmp(5)
      vect(12) = tmp(6)
      vect(13) = tmp(10)
      vect(14) = tmp(11)
      vect(15) = tmp(12)
      vect(16) = tmp(16)
      vect(17) = tmp(17)
      vect(18) = tmp(18)

    elseif(ic_type == 781)then
      !ic_type = 741
      !nn = 4
      !ndof = 6
      call DL_Shell(741, 4, 6, xx, yy, zz, rho, thick, ltype, params, vect, nsize, gausses)
      !ic_type = 781
      !nn = 8
      !ndof = 3

      tmp = 0.0
      do i=1,24
        tmp(i) = vect(i)
      enddo

      vect( 1) = tmp(1)
      vect( 2) = tmp(2)
      vect( 3) = tmp(3)
      vect( 4) = tmp(7)
      vect( 5) = tmp(8)
      vect( 6) = tmp(9)
      vect( 7) = tmp(13)
      vect( 8) = tmp(14)
      vect( 9) = tmp(15)
      vect(10) = tmp(19)
      vect(11) = tmp(20)
      vect(12) = tmp(21)
      vect(13) = tmp(4)
      vect(14) = tmp(5)
      vect(15) = tmp(6)
      vect(16) = tmp(10)
      vect(17) = tmp(11)
      vect(18) = tmp(12)
      vect(19) = tmp(16)
      vect(20) = tmp(17)
      vect(21) = tmp(18)
      vect(22) = tmp(22)
      vect(23) = tmp(23)
      vect(24) = tmp(24)

    endif

  end subroutine DL_Shell_33

  !####################################################################
  ! this subroutine can be used only for linear analysis
  subroutine UpdateST_Shell_MITC                                           &
      (etype, nn, ndof, ecoord, u, du, gausses, qf, thick, mixflag, nddisp)
    !####################################################################

    use mMechGauss
    use gauss_integration
    use m_MatMatrix

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn, mixflag
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    real(kind = kreal), intent(in)   :: u(:, :)
    real(kind = kreal), intent(in)   :: du(:, :)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(out)  :: qf(:)
    real(kind = kreal), intent(in)   :: thick

    real(kind = kreal), intent(in), optional :: nddisp(3, nn)
    !--------------------------------------------------------------------

    real(kind = kreal)   :: stiff(nn*ndof, nn*ndof), totaldisp(nn*ndof)
    integer(kind = kint) :: i

    call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, nddisp)

    totaldisp = 0.d0
    do i=1,nn
      totaldisp(ndof*(i-1)+1:ndof*i) = u(1:ndof,i) + du(1:ndof,i)
    end do

    qf = matmul(stiff,totaldisp)

  end subroutine UpdateST_Shell_MITC

  !####################################################################
  ! this subroutine can be used only for linear analysis
  subroutine UpdateST_Shell_MITC33                                         &
      (etype, nn, ndof, ecoord, u, du, gausses, qf, thick, mixflag, nddisp)
    !####################################################################

    use mMechGauss
    use gauss_integration
    use m_MatMatrix

    !--------------------------------------------------------------------

    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn, mixflag
    integer(kind = kint), intent(in) :: ndof
    real(kind = kreal), intent(in)   :: ecoord(3, nn)
    real(kind = kreal), intent(in)   :: u(3, nn*2)
    real(kind = kreal), intent(in)   :: du(3, nn*2)
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind = kreal), intent(out)  :: qf(:)
    real(kind = kreal), intent(in)   :: thick

    real(kind = kreal), intent(in), optional :: nddisp(3, nn)
    !--------------------------------------------------------------------

    real(kind = kreal)   :: stiff(nn*ndof, nn*ndof), totaldisp(nn*ndof)
    integer(kind = kint) :: i

    call STF_Shell_MITC(etype, nn, ndof, ecoord, gausses, stiff, thick, mixflag, nddisp)

    totaldisp = 0.d0
    do i=1,nn
      totaldisp(ndof*(i-1)+1:ndof*(i-1)+3) = u(1:3,2*i-1) + du(1:3,2*i-1)
      totaldisp(ndof*(i-1)+4:ndof*(i-1)+6) = u(1:3,2*i) + du(1:3,2*i)
    end do

    qf = matmul(stiff,totaldisp)

  end subroutine UpdateST_Shell_MITC33
  
  !####################################################################
  subroutine mass_Shell(etype, nn, elem, rho, thick, gausses, mass, lumped)
  !####################################################################
    use hecmw
    use m_utilities
    use mMechGauss
    use gauss_integration
    integer(kind = kint), intent(in) :: etype
    integer(kind = kint), intent(in) :: nn
    real(kind = kreal), intent(in)   :: elem(3,nn)
    real(kind = kreal), intent(in)   :: rho
    real(kind = kreal), intent(in)   :: thick
    type(tGaussStatus), intent(in)   :: gausses(:)
    real(kind=kreal), intent(out)    :: mass(:,:)      !< mass matrix
    real(kind=kreal), intent(out)    :: lumped(:)      !< lumped mass matrix

    !--------------------------------------------------------------------

    integer :: lx, ly, nsize, ndof
    integer :: fetype
    integer :: ny
    integer :: i, m
    integer :: na, nb
    integer :: jsize1, jsize2, jsize3, jsize4, jsize5, jsize6
    integer :: n_totlyr, n_layer

    real(kind = kreal) :: xi_lx, eta_lx, zeta_ly
    real(kind = kreal) :: w_w_lx, w_ly
    real(kind = kreal) :: naturalcoord(2)
    real(kind = kreal) :: nncoord(nn, 2)
    real(kind = kreal) :: shapefunc(nn)
    real(kind = kreal) :: shapederiv(nn, 2)
    real(kind = kreal) :: v1(3, nn), v2(3, nn), v3(3, nn)
    real(kind = kreal) :: v1_abs, v2_abs, v3_abs
    real(kind = kreal) :: a_over_2_v3(3, nn)
    real(kind = kreal) :: u_rot(3, nn)
    real(kind = kreal) :: dudxi_rot(3, nn), dudeta_rot(3, nn), dudzeta_rot(3, nn)
    real(kind = kreal) :: g1(3), g2(3), g3(3)
    real(kind = kreal) :: e_0(3)
    real(kind = kreal) :: det
    real(kind = kreal) :: det_cg3(3)
    real(kind = kreal) :: det_cg3_abs
    real(kind = kreal) :: w_w_w_det
    real(kind = kreal) :: N(3, 6*nn)
    real(kind = kreal) :: sumlyr, totalmass, totdiag

    ny = 0; ndof=6
    nsize = ndof*nn

    !--------------------------------------------------------------------

    ! MITC4
    if( etype == fe_mitc4_shell ) then
      fetype = fe_mitc4_shell
      ny = 2

      ! MITC9
    else if( etype == fe_mitc9_shell ) then
      fetype = fe_mitc9_shell
      ny = 3

      ! MITC3
    else if( etype == fe_mitc3_shell ) then
      fetype = fe_mitc3_shell
      ny = 2

    end if

    !-------------------------------------------------------------------

    ! xi-coordinate at a node in a local element
    ! eta-coordinate at a node in a local element
    call getNodalNaturalCoord(fetype, nncoord)

    ! xi-coordinate at the center point in a local element
    ! eta-coordinate at the center point in a local element
    naturalcoord(1) = 0.0D0
    naturalcoord(2) = 0.0D0

    call getShapeDeriv(fetype, naturalcoord, shapederiv)

    !--------------------------------------------------------------

    ! Covariant basis vector
    g1(:) = matmul( elem, shapederiv(:,1) )

    e_0(:) = g1(:)

    !--------------------------------------------------------------

    do nb = 1, nn

      !--------------------------------------------------------

      naturalcoord(1) = nncoord(nb, 1)
      naturalcoord(2) = nncoord(nb, 2)
      call getShapeDeriv(fetype, naturalcoord, shapederiv)

      !--------------------------------------------------------

      ! Covariant basis vector
      g1(:) = matmul( elem, shapederiv(:,1) )
      g2(:) = matmul( elem, shapederiv(:,2) )

      !--------------------------------------------------------

      det_cg3(1) = g1(2)*g2(3)-g1(3)*g2(2)
      det_cg3(2) = g1(3)*g2(1)-g1(1)*g2(3)
      det_cg3(3) = g1(1)*g2(2)-g1(2)*g2(1)

      det_cg3_abs = dsqrt( det_cg3(1)*det_cg3(1)   &
        +det_cg3(2)*det_cg3(2)   &
        +det_cg3(3)*det_cg3(3) )

      v3(:, nb) = det_cg3(:)/det_cg3_abs

      !--------------------------------------------------------

      v2(1, nb) = v3(2, nb)*e_0(3)-v3(3, nb)*e_0(2)
      v2(2, nb) = v3(3, nb)*e_0(1)-v3(1, nb)*e_0(3)
      v2(3, nb) = v3(1, nb)*e_0(2)-v3(2, nb)*e_0(1)

      v2_abs = dsqrt( v2(1, nb)*v2(1, nb)   &
        +v2(2, nb)*v2(2, nb)   &
        +v2(3, nb)*v2(3, nb) )

      if( v2_abs > 1.0D-15 ) then

        v2(1, nb) = v2(1, nb)/v2_abs
        v2(2, nb) = v2(2, nb)/v2_abs
        v2(3, nb) = v2(3, nb)/v2_abs

        v1(1, nb) = v2(2, nb)*v3(3, nb) &
          -v2(3, nb)*v3(2, nb)
        v1(2, nb) = v2(3, nb)*v3(1, nb) &
          -v2(1, nb)*v3(3, nb)
        v1(3, nb) = v2(1, nb)*v3(2, nb) &
          -v2(2, nb)*v3(1, nb)

        v1_abs = dsqrt( v1(1, nb)*v1(1, nb)   &
          +v1(2, nb)*v1(2, nb)   &
          +v1(3, nb)*v1(3, nb) )

        v1(1, nb) = v1(1, nb)/v1_abs
        v1(2, nb) = v1(2, nb)/v1_abs
        v1(3, nb) = v1(3, nb)/v1_abs

      else

        v1(1, nb) =  0.0D0
        v1(2, nb) =  0.0D0
        v1(3, nb) = -1.0D0

        v2(1, nb) = 0.0D0
        v2(2, nb) = 1.0D0
        v2(3, nb) = 0.0D0

      end if

      !--------------------------------------------------------

      v3(1, nb) = v1(2, nb)*v2(3, nb) &
        -v1(3, nb)*v2(2, nb)
      v3(2, nb) = v1(3, nb)*v2(1, nb) &
        -v1(1, nb)*v2(3, nb)
      v3(3, nb) = v1(1, nb)*v2(2, nb) &
        -v1(2, nb)*v2(1, nb)

      v3_abs = dsqrt( v3(1, nb)*v3(1, nb)   &
        +v3(2, nb)*v3(2, nb)   &
        +v3(3, nb)*v3(3, nb) )

      v3(1, nb) = v3(1, nb)/v3_abs
      v3(2, nb) = v3(2, nb)/v3_abs
      v3(3, nb) = v3(3, nb)/v3_abs

      !--------------------------------------------------------

      a_over_2_v3(1, nb) = 0.5D0*thick*v3(1, nb)
      a_over_2_v3(2, nb) = 0.5D0*thick*v3(2, nb)
      a_over_2_v3(3, nb) = 0.5D0*thick*v3(3, nb)

      !--------------------------------------------------------

    end do

    !--------------------------------------------------------------------

    mass(:,:) = 0.0D0
    totalmass  = 0.d0
    n_totlyr =  gausses(1)%pMaterial%totallyr
    do n_layer=1,n_totlyr
        do ly = 1, ny

          !--------------------------------------------------

          sumlyr = 0.0D0
          do m = 1, n_layer
            sumlyr = sumlyr + 2 * gausses(1)%pMaterial%shell_var(m)%weight
          end do
          zeta_ly = -1 + sumlyr - gausses(1)%pMaterial%shell_var(n_layer)%weight*(1-xg(ny, ly))

          !zeta_ly = xg(ny, ly)
          w_ly    = wgt(ny, ly)

          !--------------------------------------------------

          do lx = 1, NumOfQuadPoints(fetype)

            !--------------------------------------------

            call getQuadPoint(fetype, lx, naturalcoord)

            xi_lx  = naturalcoord(1)
            eta_lx = naturalcoord(2)

            w_w_lx = getWeight(fetype, lx)

            call getShapeFunc(fetype, naturalcoord, shapefunc)
            call getShapeDeriv(fetype, naturalcoord, shapederiv)

            !--------------------------------------------

            do na = 1, nn

              do i = 1, 3

                u_rot(i, na) = shapefunc(na)*( zeta_ly*a_over_2_v3(i, na) )

                dudxi_rot(i, na) = shapederiv(na, 1)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudeta_rot(i, na) = shapederiv(na, 2)               &
                  *( zeta_ly*a_over_2_v3(i, na) )
                dudzeta_rot(i, na) = shapefunc(na)                   &
                  *( a_over_2_v3(i, na) )

              end do

            end do

            !--------------------------------------------

            ! Covariant basis vector
            do i = 1, 3
              g1(i) = 0.0D0
              g2(i) = 0.0D0
              g3(i) = 0.0D0
              do na = 1, nn
                g1(i) = g1(i)+shapederiv(na, 1) *elem(i, na)       &
                  +dudxi_rot(i, na)
                g2(i) = g2(i)+shapederiv(na, 2) *elem(i, na)       &
                  +dudeta_rot(i, na)
                g3(i) = g3(i)+dudzeta_rot(i, na)
              end do
            end do

            !--------------------------------------------

            ! Jacobian
            det = g1(1)*( g2(2)*g3(3)-g2(3)*g3(2) ) &
              +g1(2)*( g2(3)*g3(1)-g2(1)*g3(3) ) &
              +g1(3)*( g2(1)*g3(2)-g2(2)*g3(1) )

            !--------------------------------------------

            ! [ N ] matrix
            do nb = 1, nn

              jsize1 = ndof*(nb-1)+1
              jsize2 = ndof*(nb-1)+2
              jsize3 = ndof*(nb-1)+3
              jsize4 = ndof*(nb-1)+4
              jsize5 = ndof*(nb-1)+5
              jsize6 = ndof*(nb-1)+6

              N(1, jsize1) =  shapefunc(nb)
              N(2, jsize1) =  0.0D0
              N(3, jsize1) =  0.0D0
              N(1, jsize2) =  0.0D0
              N(2, jsize2) =  shapefunc(nb)
              N(3, jsize2) =  0.0D0
              N(1, jsize3) =  0.0D0
              N(2, jsize3) =  0.0D0
              N(3, jsize3) =  shapefunc(nb)
              N(1, jsize4) =  0.0D0
              N(2, jsize4) = -u_rot(3, nb)
              N(3, jsize4) =  u_rot(2, nb)
              N(1, jsize5) =  u_rot(3, nb)
              N(2, jsize5) =  0.0D0
              N(3, jsize5) = -u_rot(1, nb)
              N(1, jsize6) = -u_rot(2, nb)
              N(2, jsize6) =  u_rot(1, nb)
              N(3, jsize6) =  0.0D0

            enddo

            !--------------------------------------------

            w_w_w_det = w_w_lx*w_ly*det
            mass(1:nsize,1:nsize) = mass(1:nsize,1:nsize)+ matmul( transpose(N), N )*w_w_w_det*rho
            totalmass = totalmass + w_w_w_det*rho
            !--------------------------------------------

          end do

          !----------------------------------------------

        end do

        !----------------------------------------------

    end do
    totalmass = totalmass*3.d0

    totdiag=0.d0
    do nb = 1, nn
        DO i = 1, 3
          lx = (nb-1)*ndof+i
          totdiag = totdiag + mass(lx,lx)
        END DO
    ENDDO
    do nb = 1, nn
        DO i = 1, 6
          lx = (nb-1)*ndof+i
          lumped(lx) = mass(lx,lx)/totdiag* totalmass
        END DO
    ENDDO

  !####################################################################
  end subroutine mass_Shell
  !####################################################################

end module m_static_LIB_shell
