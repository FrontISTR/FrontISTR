!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions for elastic material
module m_ElasticLinear
  use hecmw_util
  use mMaterial

  implicit none

contains

  !> Calculate isotropic elastic matrix
  subroutine calElasticMatrix( matl, sectType, D, temp  )
    type( tMaterial ), intent(in) :: matl       !> material properties
    integer, intent(in)           :: sectType   !> plane strain/stress or 3D
    real(kind=kreal), intent(out) :: D(:,:)     !> elastic matrix
    real(kind=kreal), intent(in)  :: temp       !> temperature
    real(kind=kreal) :: EE, PP, COEF1, COEF2, ina(1), outa(2)
    logical :: ierr

    D(:,:)=0.d0

    ina(1) = temp
    call fetch_TableData( MC_ISOELASTIC, matl%dict, outa, ierr, ina )
    if( ierr ) then
      ee = matl%variables(M_YOUNGS)
      pp = matl%variables(M_POISSON)
    else
      EE = outa(1)
      PP = outa(2)
    endif

    select case (sectType)
      case (D3)
        D(1,1)=EE*(1.d0-PP)/(1.d0-2.d0*PP)/(1.d0+PP)
        D(1,2)=EE*PP/(1.d0-2.d0*PP)/(1.d0+PP)
        D(1,3)=D(1,2)
        D(2,1)=D(1,2)
        D(2,2)=D(1,1)
        D(2,3)=D(1,2)
        D(3,1)=D(1,3)
        D(3,2)=D(2,3)
        D(3,3)=D(1,1)
        D(4,4)=EE/(1.d0+PP)*0.5d0
        D(5,5)=EE/(1.d0+PP)*0.5d0
        D(6,6)=EE/(1.d0+PP)*0.5d0
      case (PlaneStress)
        COEF1=EE/(1.d0-PP*PP)
        COEF2=0.5d0*(1.d0-PP)
        D(1,1)=COEF1
        D(1,2)=COEF1*PP
        D(1,3)=0.d0
        D(2,1)=D(1,2)
        D(2,2)=D(1,1)
        D(2,3)=0.d0
        D(3,1)=0.d0
        D(3,2)=0.d0
        D(3,3)=COEF1*COEF2
      case (Planestrain)
        COEF1=EE/((1.d0+PP)*(1.d0-2.d0*PP))
        COEF2=EE/(2.d0*(1.d0+PP))
        D(1,1)=COEF1*(1.d0-PP)
        D(1,2)=COEF1*PP
        D(1,3)=0.d0
        D(2,1)=D(1,2)
        D(2,2)=D(1,1)
        D(2,3)=0.d0
        D(3,1)=0.d0
        D(3,2)=0.d0
        D(3,3)=COEF2
      case (AxisSymetric)
        COEF1=EE*(1.d0-PP)/((1.d0+PP)*(1.d0-2.d0*PP))
        COEF2=(1.d0-2.d0*PP)/(2.d0*(1.d0-PP))
        D(1,1)=COEF1
        D(1,2)=COEF1*PP/(1.d0-PP)
        D(1,3)=0.d0
        D(1,4)=D(1,2)
        D(2,1)=D(1,2)
        D(2,2)=D(1,1)
        D(2,3)=0.d0
        D(2,4)=D(1,2)
        D(3,1)=0.d0
        D(3,2)=0.d0
        D(3,3)=COEF1*COEF2
        D(3,4)=0.d0
        D(4,1)=D(1,4)
        D(4,2)=D(2,4)
        D(4,3)=0.d0
        D(4,4)=D(1,1)
      case default
        stop "Section type not defined"
    end select

  end subroutine


  !> Calculate orthotropic elastic matrix
  subroutine calElasticMatrix_ortho( matl, sectType, bij, DMAT, temp  )
    use m_utilities
    type( tMaterial ), intent(in) :: matl       !> material properties
    integer, intent(in)           :: sectType   !> plane strain/stress or 3D
    real(kind=kreal), intent(in)  :: bij(3,3)   !> director
    real(kind=kreal), intent(out) :: DMAT(:,:)  !> elastic matrix
    real(kind=kreal), intent(in)  :: temp       !> temperature
    real(kind=kreal) :: E1, E2, E3, G12, G23, G13, nyu12, nyu23,nyu13
    real(kind=kreal) :: nyu21,nyu32,nyu31, delta1, ina(1), outa(9)
    real(kind=kreal) :: tm(6,6)
    logical :: ierr

    ina(1) = temp
    call fetch_TableData( MC_ORTHOELASTIC, matl%dict, outa, ierr, ina )
    if( ierr ) then
      stop "Fails in fetching orthotropic elastic constants!"
    endif

    E1    = outa(1)
    E2    = outa(2)
    E3    = outa(3)
    nyu12 = outa(4)
    nyu13 = outa(5)
    nyu23 = outa(6)
    G12   = outa(7)
    G13   = outa(8)
    G23   = outa(9)
    nyu21 = E2/E1*nyu12
    nyu32 = E3/E2*nyu23
    nyu31 = E3/E1*nyu13
    delta1 = 1.d0/(1.d0 -nyu12*nyu21 -nyu23*nyu32 -nyu31*nyu13 -2.d0*nyu21*nyu32*nyu13)

    DMAT(:,:)=0.d0
    DMAT(1,1) = E1*(1.d0-nyu23*nyu32)*delta1
    DMAT(2,2) = E2*(1.d0-nyu13*nyu31)*delta1
    DMAT(3,3) = E3*(1.d0-nyu12*nyu21)*delta1
    DMAT(1,2) = E1*(nyu21+nyu31*nyu23)*delta1
    DMAT(1,3) = E1*(nyu31+nyu21*nyu32)*delta1
    DMAT(2,3) = E2*(nyu32+nyu12*nyu31)*delta1
    DMAT(4,4) = G12
    DMAT(5,5) = G23
    DMAT(6,6) = G13

    DMAT(2,1) = DMAT(1,2)
    DMAT(3,2) = DMAT(2,3)
    DMAT(3,1) = DMAT(1,3)

    call transformation(bij, tm)

    dmat = matmul( transpose(tm), dmat)
    dmat = matmul( dmat, (tm) )

  end subroutine

  !< calculate contravariant elastic tensor
  !####################################################################
  subroutine LinearElastic_Shell                     &
      (matl, sectType, c,                     &
      e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
      alpha, n_layer)
    !####################################################################

    type( tMaterial ), intent(in)   :: matl
    integer, intent(in)             :: sectType, n_layer
    real(kind = kreal), intent(out) :: c(:, :, :, :)
    real(kind = kreal), intent(in)  :: e1_hat(3), e2_hat(3), e3_hat(3) !< local Orthonormal frame
    real(kind = kreal), intent(in)  :: cg1(3), cg2(3), cg3(3)            !< contravaraint frame
    real(kind = kreal), intent(out) :: alpha

    !--------------------------------------------------------------------

    real(kind = kreal) :: ee, pp, ee2, g12, g23, g31, theta, pp2
    real(kind = kreal) :: outa(2)
    real(kind = kreal) :: lambda1, lambda2, mu, k_correction
    real(kind = kreal) :: c_hat(3, 3, 3, 3), D(5,5), D_hat(5,5), D_temp(5,5), T(5,5)
    real(kind = kreal) :: e_hat_dot_cg(3, 3)
    real(kind = kreal) :: alpha_over_mu

    integer :: index_i(5), index_j(5), index_k(5), index_l(5)
    integer :: is, js, ii, ij, ik, il
    integer :: i, j, k, l, total_layer, matl_type

    logical :: ierr

    !--------------------------------------------------------------------

    call fetch_TableData(MC_ISOELASTIC, matl%dict, outa, ierr)

    !--------------------------------------------------------------------

    matl_type = matl%shell_var(n_layer)%ortho
    if(matl_type == 0)then
      if( ierr ) then

        ee = matl%shell_var(n_layer)%ee
        pp = matl%shell_var(n_layer)%pp
        alpha_over_mu = matl%variables(M_ALPHA_OVER_MU)

      else

        ee = outa(1)
        pp = outa(2)
        alpha_over_mu = matl%variables(M_ALPHA_OVER_MU)

      end if

      !--------------------------------------------------------------------

      ! Elastic constant
      lambda1 = ee/( 1.0D0-pp*pp )
      lambda2 = pp*lambda1
      mu      = 0.5D0*ee/( 1.0D0+pp )

      !--------------------------------------------------------------------

      ! Shear correction factor
      k_correction = 5.0D0/6.0D0
      !k_correction = 1.0D0

      !--------------------------------------------------------------------

      alpha = alpha_over_mu*mu

      !--------------------------------------------------------------------

      ! Constitutive tensor
      c_hat(:, :, :, :) = 0.0D0

      !--------------------------------------------------------

      c_hat(1, 1, 1, 1) = lambda1
      c_hat(1, 1, 2, 2) = lambda2
      c_hat(2, 2, 1, 1) = lambda2
      c_hat(2, 2, 2, 2) = lambda1
      c_hat(1, 2, 1, 2) = mu
      c_hat(1, 2, 2, 1) = mu
      c_hat(2, 1, 1, 2) = mu
      c_hat(2, 1, 2, 1) = mu
      c_hat(1, 3, 1, 3) = k_correction*mu
      c_hat(1, 3, 3, 1) = k_correction*mu
      c_hat(2, 3, 2, 3) = k_correction*mu
      c_hat(2, 3, 3, 2) = k_correction*mu
      c_hat(3, 1, 3, 1) = k_correction*mu
      c_hat(3, 1, 1, 3) = k_correction*mu
      c_hat(3, 2, 3, 2) = k_correction*mu
      c_hat(3, 2, 2, 3) = k_correction*mu

      !--------------------------------------------------------

    elseif(matl_type == 1)then
      total_layer = matl%totallyr

      ee = matl%shell_var(n_layer)%ee
      pp = matl%shell_var(n_layer)%pp
      ee2 = matl%shell_var(n_layer)%ee2
      g12 = matl%shell_var(n_layer)%g12
      g23 = matl%shell_var(n_layer)%g23
      g31 = matl%shell_var(n_layer)%g31
      theta = matl%shell_var(n_layer)%angle

      alpha_over_mu = matl%variables(M_ALPHA_OVER_MU)

      !--------------------------------------------------------------------

      ! Shear correction factor
      k_correction = 5.0D0/6.0D0

      !--------------------------------------------------------------------

      alpha = alpha_over_mu * 0.5D0 * ee / ( 1.0D0+pp )

      !--------------------------------------------------------------------

      !    write(*,*) ee,pp,ee2,g12,g23,g31,theta,alpha_over_mu,n_layer

      D(:,:) = 0.0D0
      D_hat(:,:) = 0.0D0
      D_temp(:,:) = 0.0D0
      T(:,:) = 0.0D0

      pp2 = pp * ee2 / ee

      D(1,1) = ee / (1.0D0 - pp * pp2)
      D(1,2) = pp2 * ee / (1.0D0- pp * pp2)
      D(2,1) = pp2 * ee / (1.0D0- pp * pp2)
      D(2,2) = ee2 / (1.0D0 - pp * pp2)
      D(3,3) = g12
      D(4,4) = g23
      D(5,5) = g31

      T(1,1) = cos(theta) * cos(theta)
      T(1,2) = sin(theta) * sin(theta)
      T(2,1) = sin(theta) * sin(theta)
      T(2,2) = cos(theta) * cos(theta)
      T(3,3) = cos(theta) * cos(theta) - sin(theta) * sin(theta)
      T(1,3) = sin(theta) * cos(theta)
      T(2,3) = -sin(theta) * cos(theta)
      T(3,1) = -2.0D0 * sin(theta) * cos(theta)
      T(3,2) =  2.0D0 * sin(theta) * cos(theta)
      T(4,4) = cos(theta)
      T(4,5) = sin(theta)
      T(5,4) = -sin(theta)
      T(5,5) = cos(theta)

      !--------------------- D_temp = [D]*[T]

      D_temp(1,1) = D(1,1)*T(1,1)+D(1,2)*T(2,1)
      D_temp(1,2) = D(1,1)*T(1,2)+D(1,2)*T(2,2)
      D_temp(2,1) = D(2,1)*T(1,1)+D(2,2)*T(2,1)
      D_temp(2,2) = D(2,1)*T(1,2)+D(2,2)*T(2,2)
      D_temp(3,1) = D(3,3)*T(3,1)
      D_temp(3,2) = D(3,3)*T(3,2)
      D_temp(1,3) = D(1,1)*T(1,3)+D(1,2)*T(2,3)
      D_temp(2,3) = D(2,1)*T(1,3)+D(2,2)*T(2,3)
      D_temp(3,3) = D(3,3)*T(3,3)
      D_temp(4,4) = D(4,4)*T(4,4)
      D_temp(4,5) = D(4,4)*T(4,5)
      D_temp(5,4) = D(5,5)*T(5,4)
      D_temp(5,5) = D(5,5)*T(5,5)

      !--------------------- D_hat = [trans_T]*[D_temp]

      D_hat(1,1) = T(1,1)*D_temp(1,1)+T(1,2)*D_temp(2,1)+T(3,1)*D_temp(3,1)
      D_hat(1,2) = T(1,1)*D_temp(1,2)+T(1,2)*D_temp(2,2)+T(3,1)*D_temp(3,2)
      D_hat(2,1) = T(2,1)*D_temp(1,1)+T(2,2)*D_temp(2,1)+T(3,2)*D_temp(3,1)
      D_hat(2,2) = T(2,1)*D_temp(1,2)+T(2,2)*D_temp(2,2)+T(3,2)*D_temp(3,2)
      D_hat(3,1) = T(1,3)*D_temp(1,1)+T(2,3)*D_temp(2,1)+T(3,3)*D_temp(3,1)
      D_hat(3,2) = T(1,3)*D_temp(1,2)+T(2,3)*D_temp(2,2)+T(3,3)*D_temp(3,2)
      D_hat(1,3) = T(1,1)*D_temp(1,3)+T(1,2)*D_temp(2,3)+T(3,1)*D_temp(3,3)
      D_hat(2,3) = T(2,1)*D_temp(1,3)+T(2,2)*D_temp(2,3)+T(3,2)*D_temp(3,3)
      D_hat(3,3) = T(1,3)*D_temp(1,3)+T(2,3)*D_temp(2,3)+T(3,3)*D_temp(3,3)
      D_hat(4,4) = T(4,4)*D_temp(4,4)+T(5,4)*D_temp(5,4)
      D_hat(4,5) = T(4,4)*D_temp(4,5)+T(5,4)*D_temp(5,5)
      D_hat(5,4) = T(4,5)*D_temp(4,4)+T(5,5)*D_temp(5,4)
      D_hat(5,5) = T(4,5)*D_temp(4,5)+T(5,5)*D_temp(5,5)

      !--------------------------------------------------------------------

      ! Constitutive tensor
      c_hat(:, :, :, :) = 0.0D0
      !write(*,*) 'Elastic linear.f90 make c_hat'

      !-------D_hat to c_hat

      index_i(1) = 1
      index_i(2) = 2
      index_i(3) = 1
      index_i(4) = 2
      index_i(5) = 3

      index_j(1) = 1
      index_j(2) = 2
      index_j(3) = 2
      index_j(4) = 3
      index_j(5) = 1

      index_k(1) = 1
      index_k(2) = 2
      index_k(3) = 1
      index_k(4) = 2
      index_k(5) = 3

      index_l(1) = 1
      index_l(2) = 2
      index_l(3) = 2
      index_l(4) = 3
      index_l(5) = 1

      !--------------------------------------------------------------------

      do js = 1, 5
        do is = 1, 5

          ii = index_i(is)
          ij = index_j(is)
          ik = index_k(js)
          il = index_l(js)

          c_hat(ii, ij, ik, il) = D_hat(is, js)

        end do
      end do

      !--------------------------------------------------------

    else
      write(*,*) 'shell matl type isnot collect'
      stop
    endif

    select case( sectType )
      case( Shell )

        e_hat_dot_cg(1, 1)                  &
          = e1_hat(1)*cg1(1)+e1_hat(2)*cg1(2) &
          +e1_hat(3)*cg1(3)
        e_hat_dot_cg(2, 1)                  &
          = e2_hat(1)*cg1(1)+e2_hat(2)*cg1(2) &
          +e2_hat(3)*cg1(3)
        e_hat_dot_cg(3, 1) = 0.0D0
        e_hat_dot_cg(1, 2)                  &
          = e1_hat(1)*cg2(1)+e1_hat(2)*cg2(2) &
          +e1_hat(3)*cg2(3)
        e_hat_dot_cg(2, 2)                  &
          = e2_hat(1)*cg2(1)+e2_hat(2)*cg2(2) &
          +e2_hat(3)*cg2(3)
        e_hat_dot_cg(3, 2) = 0.0D0
        e_hat_dot_cg(1, 3)                  &
          = e1_hat(1)*cg3(1)+e1_hat(2)*cg3(2) &
          +e1_hat(3)*cg3(3)
        e_hat_dot_cg(2, 3)                  &
          = e2_hat(1)*cg3(1)+e2_hat(2)*cg3(2) &
          +e2_hat(3)*cg3(3)
        e_hat_dot_cg(3, 3)                  &
          = e3_hat(1)*cg3(1)+e3_hat(2)*cg3(2) &
          +e3_hat(3)*cg3(3)

        !--------------------------------------------------------

        ! Constitutive tensor

        c(1, 1, 1, 1) = 0.0D0
        c(2, 2, 1, 1) = 0.0D0
        c(1, 2, 1, 1) = 0.0D0
        c(2, 2, 2, 2) = 0.0D0
        c(1, 2, 2, 2) = 0.0D0
        c(1, 2, 1, 2) = 0.0D0
        c(3, 1, 1, 1) = 0.0D0
        c(3, 1, 2, 2) = 0.0D0
        c(3, 1, 1, 2) = 0.0D0
        c(2, 3, 1, 1) = 0.0D0
        c(2, 3, 2, 2) = 0.0D0
        c(2, 3, 1, 2) = 0.0D0
        c(3, 1, 3, 1) = 0.0D0
        c(3, 1, 2, 3) = 0.0D0
        c(2, 3, 2, 3) = 0.0D0

        do l = 1, 2

          do k = 1, 2

            do j = 1, 2

              do i = 1, 2

                c(1, 1, 1, 1)                            &
                  = c(1, 1, 1, 1)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j ,1) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1)
                c(2, 2, 1, 1)                            &
                  = c(2, 2, 1, 1)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 2) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1)
                c(1, 2, 1, 1)                            &
                  = c(1, 2, 1, 1)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j, 2) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1)
                c(2, 2, 2, 2)                            &
                  = c(2, 2, 2, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 2) &
                  *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2)

                c(1, 2, 2, 2)                            &
                  = c(1, 2, 2, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j, 2) &
                  *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2)
                c(1, 2, 1, 2)                            &
                  = c(1, 2, 1, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j, 2) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 2)

              end do

              do i = 1, 3

                c(3, 1, 1, 1)                            &
                  = c(3, 1, 1, 1)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1)
                c(3, 1, 2, 2)                            &
                  = c(3, 1, 2, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
                  *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2)
                c(3, 1, 1, 2)                            &
                  = c(3, 1, 1, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 2)

              end do

            end do

            do j = 1, 3

              do i = 1, 2

                c(2, 3, 1, 1)                            &
                  = c(2, 3, 1, 1)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1)
                c(2, 3, 2, 2)                            &
                  = c(2, 3, 2, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
                  *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2)
                c(2, 3, 1, 2)                            &
                  = c(2, 3, 1, 2)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
                  *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 2)

              end do

            end do

          end do

          do k = 1, 3

            do j = 1, 2

              do i = 1, 3

                c(3, 1, 3, 1)                            &
                  = c(3, 1, 3, 1)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
                  *e_hat_dot_cg(k, 3)*e_hat_dot_cg(l, 1)

              end do

            end do

          end do

        end do

        do l = 1, 3

          do k = 1, 2

            do j = 1, 2

              do i = 1, 3

                c(3, 1, 2, 3)                            &
                  = c(3, 1, 2, 3)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
                  *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 3)

              end do

            end do

            do j = 1, 3

              do i = 1, 2

                c(2, 3, 2, 3)                            &
                  = c(2, 3, 2, 3)                          &
                  +c_hat(i, j, k, l)                      &
                  *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
                  *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 3)

              end do

            end do

          end do

        end do

        c(1, 1, 2, 2) = c(2, 2, 1, 1)
        c(1, 1, 1, 2) = c(1, 2, 1, 1)
        c(1, 1, 2, 3) = c(2, 3, 1, 1)
        c(1, 1, 3, 1) = c(3, 1, 1, 1)
        c(2, 2, 1, 2) = c(1, 2, 2, 2)
        c(2, 2, 2, 3) = c(2, 3, 2, 2)
        c(2, 2, 3, 1) = c(3, 1, 2, 2)
        c(1, 2, 2, 3) = c(2, 3, 1, 2)
        c(1, 2, 3, 1) = c(3, 1, 1, 2)
        c(2, 3, 3, 1) = c(3, 1, 2, 3)

        !--------------------------------------------------------

        ! DO l = 1, 3
        !
        !  DO k = 1, 3
        !
        !   DO j = 1, 3
        !
        !    DO i = 1, 3
        !
        !     c(i, j, k, l) = 0.0D0
        !
        !     DO ll = 1, 3
        !
        !      DO kk = 1, 3
        !
        !       DO jj = 1, 3
        !
        !        DO ii = 1, 3
        !
        !         c(i, j, k, l)                              &
          !         = c(i, j, k, l)                            &
          !          +c_hat(ii, jj, kk, ll)                    &
          !           *e_hat_dot_cg(ii, i)*e_hat_dot_cg(jj, j) &
          !           *e_hat_dot_cg(kk, k)*e_hat_dot_cg(ll, l)
        !
        !        END DO
        !
        !       END DO
        !
        !      END DO
        !
        !     END DO
        !
        !    END DO
        !
        !   END DO
        !
        !  END DO
        !
        ! END DO

        !--------------------------------------------------------

    end select

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine LinearElastic_Shell
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)



end module m_ElasticLinear
