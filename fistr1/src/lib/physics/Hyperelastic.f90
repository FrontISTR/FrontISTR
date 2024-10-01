!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions for hyperelastic calculation
module mHyperElastic
  use hecmw_util
  use mMaterial
  implicit none

contains

  !> This subroutine calculates derivative of the invariant with respect to Cauchy-Green tensor
  subroutine cderiv( matl, sectType, ctn, itn,               &
      inv1b, inv2b, inv3b, dibdc, d2ibdc2, strain, direction, inv4b, dibdc_ani, d2ibdc2_ani )
    type( tMaterial ), intent(in) :: matl                  !< material rpoperties
    integer, intent(in)           :: sectType              !< not used currently
    real(kind=kreal), intent(out) :: inv1b                 !< invariants
    real(kind=kreal), intent(out) :: inv2b                 !< invariants
    real(kind=kreal), intent(out) :: inv3b                 !< invariants
    real(kind=kreal), intent(out) :: dibdc(3,3,3)          !< derivative of the invariant with respect to c(i,j)
    real(kind=kreal), intent(out) :: d2ibdc2(3,3,3,3,3)    !< derivative of the invariant with respect to c(i,j)
    real(kind=kreal), intent(in)  :: strain(6)             !< Cauchy-Lagrange strain tensor
    real(kind=kreal), intent(out) :: ctn(3,3)              !< right Cauchy-Green deformation tensor
    real(kind=kreal), intent(out) :: itn(3,3)              !< identity tensor
    real(kind=kreal), intent(in), optional :: direction(3)       !< direction vector of anisotropic materials
    real(kind=kreal), intent(out), optional :: inv4b                 !< invariants
    real(kind=kreal), intent(out), optional :: dibdc_ani(3,3)        !< derivative of the 4th anisotropic invariant
    real(kind=kreal), intent(out), optional :: d2ibdc2_ani(3,3,3,3)  !< 2nd derivative of the 4th anisotropic invariant

    integer :: i, j, k, l, m, n

    real(kind=kreal) :: inv1, inv2, inv3, inv33
    real(kind=kreal) :: delta(3,3)
    real(kind=kreal) :: didc(3,3,3), ctninv(3,3)
    real(kind=kreal) :: d2idc2(3,3,3,3,3)

    !   Kronecker's delta
    delta(:,:) = 0.d0
    delta(1,1) = 1.d0
    delta(2,2) = 1.d0
    delta(3,3) = 1.d0

    ! Green-Lagrange strain to right Cauchy-Green deformation tensor
    ctn(1,1)=strain(1)*2.d0+1.d0
    ctn(2,2)=strain(2)*2.d0+1.d0
    ctn(3,3)=strain(3)*2.d0+1.d0
    ctn(1,2)=strain(4);   ctn(2,1)=ctn(1,2)
    ctn(2,3)=strain(5);   ctn(3,2)=ctn(2,3)
    ctn(3,1)=strain(6);   ctn(1,3)=ctn(3,1)

    itn(:,:) = delta(:,:)

    ! ----- calculate the invariant of C
    inv1 = ctn(1,1)+ctn(2,2)+ctn(3,3)
    inv2 = ctn(2,2)*ctn(3,3)+ctn(1,1)*ctn(3,3)+ctn(1,1)*ctn(2,2)                   &
      -ctn(2,3)*ctn(2,3)-ctn(1,3)*ctn(1,3)-ctn(1,2)*ctn(1,2)
    inv3 =  ctn(1,1)*ctn(2,2)*ctn(3,3)                                             &
      +ctn(2,1)*ctn(3,2)*ctn(1,3)                                             &
      +ctn(3,1)*ctn(1,2)*ctn(2,3)                                             &
      -ctn(3,1)*ctn(2,2)*ctn(1,3)                                             &
      -ctn(2,1)*ctn(1,2)*ctn(3,3)                                             &
      -ctn(1,1)*ctn(3,2)*ctn(2,3)
    inv33 = inv3**(-1.d0/3.d0)

    ! ----- calculate the inverse of C
    ctninv(1,1)=(ctn(2,2)*ctn(3,3)-ctn(2,3)*ctn(2,3))/inv3
    ctninv(2,2)=(ctn(1,1)*ctn(3,3)-ctn(1,3)*ctn(1,3))/inv3
    ctninv(3,3)=(ctn(1,1)*ctn(2,2)-ctn(1,2)*ctn(1,2))/inv3
    ctninv(1,2)=(ctn(1,3)*ctn(2,3)-ctn(1,2)*ctn(3,3))/inv3
    ctninv(1,3)=(ctn(1,2)*ctn(2,3)-ctn(2,2)*ctn(1,3))/inv3
    ctninv(2,3)=(ctn(1,2)*ctn(1,3)-ctn(1,1)*ctn(2,3))/inv3
    ctninv(2,1)=ctninv(1,2)
    ctninv(3,1)=ctninv(1,3)
    ctninv(3,2)=ctninv(2,3)

    ! ----- calculate the derivative of the C-invariant with respect to c(i,j)
    do i=1,3
      do j=1,3
        didc(i,j,1) = delta(i,j)
        didc(i,j,2) = inv1*delta(i,j)-ctn(i,j)
        didc(i,j,3) = inv3*ctninv(i,j)
      enddo
    enddo

    ! ----- calculate the second derivative of the C-invariant
    do k=1,3
      do l=1,3
        do m=1,3
          do n=1,3
            d2idc2(k,l,m,n,1)=0.d0
            d2idc2(k,l,m,n,2)=delta(k,l)*delta(m,n)-             &
              (delta(k,m)*delta(l,n)+delta(k,n)*delta(l,m))/2.d0
            d2idc2(k,l,m,n,3)=inv3*(ctninv(m,n)*ctninv(k,l)-   &
              (ctninv(k,m)*ctninv(n,l)+ctninv(k,n)*ctninv(m,l))/2.d0)
          enddo
        enddo
      enddo
    enddo

    ! ----- derivatives for the reduced invariants
    inv1b = inv1*inv33
    inv2b = inv2*inv33*inv33
    inv3b = dsqrt(inv3)

    !   --- first derivative the reduced c-invarians w.r.t. c(i,j)
    do i=1,3
      do j=1,3
        dibdc(i,j,1) = -inv33**4*inv1*didc(i,j,3)/3.d0 + inv33*didc(i,j,1)
        dibdc(i,j,2) = -2.d0*inv33**5*inv2*didc(i,j,3)/3.d0 + inv33**2*didc(i,j,2)
        dibdc(i,j,3) = didc(i,j,3)/(2.d0*dsqrt(inv3))
      enddo
    enddo

    !   --- second derivative of the reduced c-invariants w.r.t. c(i,j)
    do i=1,3
      do j=1,3
        do k=1,3
          do l=1,3
            d2ibdc2(i,j,k,l,1) =  4.d0/9.d0*inv33**7*inv1*didc(i,j,3)*didc(k,l,3)                       &
              -inv33**4/3.d0*(didc(k,l,1)*didc(i,j,3)+didc(i,j,1)*didc(k,l,3))       &
              -inv33**4/3.d0*inv1*d2idc2(i,j,k,l,3)                                  &
              +inv33*d2idc2(i,j,k,l,1)
            d2ibdc2(i,j,k,l,2) = 10.d0/9.d0*inv33**8*inv2*didc(i,j,3)*didc(k,l,3)                       &
              -2.d0/3.d0*inv33**5*(didc(k,l,2)*didc(i,j,3)+didc(i,j,2)*didc(k,l,3))  &
              -2.d0/3.d0*inv33**5*inv2*d2idc2(i,j,k,l,3)                             &
              +inv33**2*d2idc2(i,j,k,l,2)
            d2ibdc2(i,j,k,l,3) = -didc(i,j,3)*didc(k,l,3)/(4.d0*inv3**1.5d0)                            &
              +d2idc2(i,j,k,l,3)/(2.d0*dsqrt(inv3))
          enddo
        enddo
      enddo
    enddo

    if( present(direction) ) call cderiv_aniso( ctn, inv3, didc(:,:,3), d2idc2(:,:,:,:,3), &
      direction, inv4b, dibdc_ani, d2ibdc2_ani )

  end subroutine cderiv

  !> This subroutine calculates derivative of the 4th reduced invariant I4 with respect to Cauchy-Green tensor
  !> where I4 = I_3^{-1/3}*a_i*C_{ij}*a_j and a_i = direc(i).
  subroutine cderiv_aniso( ctn, inv3, didc_3, d2idc2_3, direction, inv4b, dibdc_ani, d2ibdc2_ani )
    real(kind=kreal), intent(in)  :: ctn(3,3)              !< right Cauchy-Green deformation tensor
    real(kind=kreal), intent(in)  :: inv3                  !< 3rd invariants
    real(kind=kreal), intent(in)  :: didc_3(3,3)           !< derivative of the 3rd invariant with respect to c(i,j)
    real(kind=kreal), intent(in)  :: d2idc2_3(3,3,3,3)     !< derivative of the 3rd invariant with respect to c(i,j)
    real(kind=kreal), intent(in)  :: direction(3)          !< direction vector of anisotropic materials
    real(kind=kreal), intent(out) :: inv4b                 !< 4th reduced invariant I4
    real(kind=kreal), intent(out) :: dibdc_ani(3,3)        !< derivative of the 4th anisotropic invariant
    real(kind=kreal), intent(out) :: d2ibdc2_ani(3,3,3,3)  !< 2nd derivative of the 4th anisotropic invariant

    integer :: i, j, k, l, m, n
    real(kind=kreal) :: inv4, inv33, inv3m43, inv4d3
    real(kind=kreal) :: d2ibdc24

    inv33 = inv3**(-1.d0/3.d0) ! = I_3^{-1/3}
    inv3m43 = inv33/inv3       ! = I_3^{-4/3}

    inv4 = 0.d0
    do i=1,3
      do j=1,3
        inv4 = inv4 + direction(j)*ctn(j,i)*direction(i)
      enddo
    enddo
    inv4b = inv33*inv4 ! I4
    inv4d3 = inv4/inv3   ! I4*I_3^{-1}

    !   --- first derivative the reduced 4th c-invarians w.r.t. c(i,j)
    do i=1,3
      do j=1,3
        dibdc_ani(i,j) = inv33*( (-1.d0/3.d0)*didc_3(i,j)*inv4d3+direction(i)*direction(j) )
      enddo
    enddo

    !   --- second derivative of the reduced c-invariants w.r.t. c(i,j)
    do k=1,3
      do l=1,3
        do m=1,3
          do n=1,3
            d2ibdc24 = (2.d0/3.d0)*didc_3(k,l)*didc_3(m,n)*inv4d3
            d2ibdc24 = d2ibdc24 - d2idc2_3(k,l,m,n)*inv4
            d2ibdc24 = d2ibdc24 - direction(k)*direction(l)*didc_3(m,n)
            d2ibdc24 = d2ibdc24 - didc_3(k,l)*direction(m)*direction(n)
            d2ibdc2_ani(k,l,m,n) = (1.d0/3.d0)*inv3m43*d2ibdc24
          enddo
        enddo
      enddo
    enddo

  end subroutine


  !-------------------------------------------------------------------------------
  !> This subroutine provides elastic tangent coefficient for Arruda-Boyce hyperelastic material
  !
  !-------------------------------------------------------------------------------
  subroutine calElasticArrudaBoyce( matl, sectType, cijkl, strain )
    type( tMaterial ), intent(in) :: matl             !< material rpoperties
    integer, intent(in)           :: sectType         !< not used currently
    real(kind=kreal), intent(out) :: cijkl(3,3,3,3)   !< constitutive relation
    real(kind=kreal), intent(in)  :: strain(6)        !< Cauchy-Lagrange strain tensor

    integer :: i, j, k, l
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: inv1b, inv2b, inv3b
    real(kind=kreal) :: dibdc(3,3,3)
    real(kind=kreal) :: d2ibdc2(3,3,3,3,3)
    real(kind=kreal) :: constant(3), coef

    call cderiv(  matl, sectType, ctn, itn, inv1b, inv2b, inv3b,            &
      dibdc, d2ibdc2, strain    )
    constant(1:3)=matl%variables(M_PLCONST1:M_PLCONST3)
    coef = constant(2)

    do l=1,3 
      do k=1,3
        do j=1,3
          do i=1,3
            cijkl(i,j,k,l) =  constant(1)*(1.d0/(10.d0*coef**2)                            &
              +66.d0*inv1b/(1050.d0*coef**4)                 &
              +228.d0*inv1b**2/(7000.d0*coef**6)             &
              +10380.d0*inv1b**3/(673750.d0*coef**8))        &
              *dibdc(i,j,1)*dibdc(k,l,1)                       &
              +constant(1)*(0.5d0+inv1b/(10.d0*coef**2)                    &
              +33.d0*inv1b**2/(1050.d0*coef**4)              &
              +76.d0*inv1b**3/(7000.d0*coef**6)              &
              +2595.d0*inv1b**4/(673750.d0*coef**8))         &
              *d2ibdc2(i,j,k,l,1)                              &
              +(1.d0+1.d0/inv3b**2)*dibdc(i,j,3)*dibdc(k,l,3)/constant(3)  &
              +(inv3b-1.d0/inv3b)*d2ibdc2(i,j,k,l,3)/constant(3)
          enddo
        enddo
      enddo
    enddo
    cijkl(:,:,:,:) = 4.d0*cijkl(:,:,:,:)

  end subroutine calElasticArrudaBoyce

  !-------------------------------------------------------------------------------
  !> This subroutine provides to update stress and strain for Arrude-Royce material
  subroutine calUpdateElasticArrudaBoyce( matl, sectType, dstrain, dstress )
    type( tMaterial ), intent(in) :: matl         !> material properties
    integer, intent(in)           :: sectType     !> not used currently
    real(kind=kreal), intent(out) :: dstress(6)   !> Cauchy-Green strain
    real(kind=kreal), intent(in)  :: dstrain(6)   !> 2nd Piola-Kirchhoff stress

    integer :: i, j
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: inv1b, inv2b, inv3b
    real(kind=kreal) :: dibdc(3,3,3)
    real(kind=kreal) :: d2ibdc2(3,3,3,3,3)
    real(kind=kreal) :: constant(3), coef
    real(kind=kreal) :: PKstress(3,3)


    constant(1:3)=matl%variables(M_PLCONST1:M_PLCONST3)
    coef = constant(2)
    call cderiv(  matl, sectType, ctn, itn,      &
      inv1b, inv2b, inv3b,           &
      dibdc, d2ibdc2, dstrain    )


    ! ----- calculate the stress
    do i=1,3
      do j=1,3
        PKstress(i,j) = constant(1)*( 0.5d0+inv1b/(10.d0*coef**2)                     &
          +33.d0*inv1b*inv1b/(1050.d0*coef**4)                 &
          +76.d0*inv1b**3/(7000.d0*coef**6)                    &
          +2595.d0*inv1b**4/(673750.d0*coef**8))*dibdc(i,j,1)  &
          +(inv3b-1.d0/inv3b)*dibdc(i,j,3)/constant(3)
      enddo
    enddo

    dstress(1) = 2.d0*PKstress(1,1)
    dstress(2) = 2.d0*PKstress(2,2)
    dstress(3) = 2.d0*PKstress(3,3)
    dstress(4) = PKstress(1,2) + PKstress(2,1)
    dstress(5) = PKstress(2,3) + PKstress(3,2)
    dstress(6) = PKstress(1,3) + PKstress(3,1)

  end subroutine calUpdateElasticArrudaBoyce

  !-------------------------------------------------------------------------------
  !> This subroutine provides elastic tangent coefficient for Mooney-Rivlin hyperelastic material
  !-------------------------------------------------------------------------------
  subroutine calElasticMooneyRivlin( matl, sectType, cijkl, strain, hdflag )
    type( tMaterial ), intent(in) :: matl             !< material rpoperties
    integer, intent(in)           :: sectType         !< not used curr
    real(kind=kreal), intent(out) :: cijkl(3,3,3,3)   !< constitutive relation
    real(kind=kreal), intent(in)  :: strain(6)        !< Cauchy-Lagrange strain tensor
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: k, l, m, n
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: inv1b, inv2b, inv3b
    real(kind=kreal) :: dibdc(3,3,3)
    real(kind=kreal) :: d2ibdc2(3,3,3,3,3)
    real(kind=kreal) :: constant(3)


    constant(1:3)=matl%variables(M_PLCONST1:M_PLCONST3)
    if( present(hdflag) ) then
      if( hdflag == 1 ) then
        constant(3)=0.d0
      else if( hdflag == 2 ) then
        constant(1:2)=0.d0
      end if
    end if
    if( dabs(constant(3))>1.d-12 ) constant(3) = 1.d0/constant(3)

    call cderiv( matl, sectType, ctn, itn, inv1b, inv2b, inv3b,            &
      dibdc, d2ibdc2, strain    )

    do n=1,3
      do m=1,3
        do l=1,3
          do k=1,3
            cijkl(k,l,m,n) = d2ibdc2(k,l,m,n,1)*constant(1) +  &
              d2ibdc2(k,l,m,n,2)*constant(2) +  &
              2.d0*(dibdc(k,l,3)*dibdc(m,n,3)+  &
              (inv3b-1.d0)*d2ibdc2(k,l,m,n,3))*constant(3)
          enddo
        enddo
      enddo
    enddo
    cijkl(:,:,:,:)=4.d0*cijkl(:,:,:,:)

  end subroutine calElasticMooneyRivlin

  !-------------------------------------------------------------------------------
  !> This subroutine provides to update stress and strain for Mooney-Rivlin material
  !-------------------------------------------------------------------------------
  subroutine calUpdateElasticMooneyRivlin( matl, sectType, strain, stress, hdflag )
    type( tMaterial ), intent(in) :: matl        !< material properties
    integer, intent(in)           :: sectType    !< not used currently
    real(kind=kreal), intent(out) :: stress(6)   !< 2nd Piola-Kirchhoff stress
    real(kind=kreal), intent(in)  :: strain(6)   !< Green-Lagrangen strain
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: k, l
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: inv1b, inv2b, inv3b
    real(kind=kreal) :: dibdc(3,3,3)
    real(kind=kreal) :: d2ibdc2(3,3,3,3,3)
    real(kind=kreal) :: constant(3)
    real(kind=kreal) :: dudc(3,3)

    constant(1:3)=matl%variables(M_PLCONST1:M_PLCONST3)
    if( present(hdflag) ) then
      if( hdflag == 1 ) then
        constant(3)=0.d0
      else if( hdflag == 2 ) then
        constant(1:2)=0.d0
      end if
    end if
    if( dabs(constant(3))>1.d-12 ) constant(3) = 1.d0/constant(3)

    call cderiv( matl, sectType, ctn, itn, inv1b, inv2b, inv3b,      &
      dibdc, d2ibdc2, strain    )


    ! ----- stress
    do l=1,3
      do k=1,3
        dudc(k,l) = dibdc(k,l,1)*constant(1)+dibdc(k,l,2)*constant(2)  &
          +2.d0*(inv3b-1.d0)*dibdc(k,l,3)*constant(3)
      enddo
    enddo

    stress(1)=2.d0*dudc(1,1)
    stress(2)=2.d0*dudc(2,2)
    stress(3)=2.d0*dudc(3,3)
    stress(4)=2.d0*dudc(1,2)
    stress(5)=2.d0*dudc(2,3)
    stress(6)=2.d0*dudc(1,3)

  end subroutine calUpdateElasticMooneyRivlin

  !-------------------------------------------------------------------------------
  !> This subroutine provides elastic tangent coefficient for anisotropic Mooney-Rivlin hyperelastic material
  !-------------------------------------------------------------------------------
  subroutine calElasticMooneyRivlinAniso( matl, sectType, cijkl, strain, cdsys, hdflag  )
    type( tMaterial ), intent(in) :: matl             !< material rpoperties
    integer, intent(in)           :: sectType         !< not used curr
    real(kind=kreal), intent(out) :: cijkl(3,3,3,3)   !< constitutive relation
    real(kind=kreal), intent(in)  :: strain(6)        !< Cauchy-Lagrange strain tensor
    real(kind=kreal), intent(in)  :: cdsys(3,3) !> material coordinate system
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: i, j, k, l, m, n, jj
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: inv1b, inv2b, inv3b, inv4b
    real(kind=kreal) :: dibdc(3,3,3)
    real(kind=kreal) :: d2ibdc2(3,3,3,3,3)
    real(kind=kreal) :: constant(5)
    real(kind=kreal) :: dibdc_ani(3,3)
    real(kind=kreal) :: d2ibdc2_ani(3,3,3,3)


    constant(1:5)=matl%variables(M_PLCONST1:M_PLCONST5)
    if( present(hdflag) ) then
      if( hdflag == 1 ) then
        constant(3)=0.d0
      else if( hdflag == 2 ) then
        constant(1:2)=0.d0
        constant(4:5)=0.d0
      end if
    end if
    if( dabs(constant(3))>1.d-12 ) constant(3) = 1.d0/constant(3)

    call cderiv( matl, sectType, ctn, itn, inv1b, inv2b, inv3b,      &
      dibdc, d2ibdc2, strain, cdsys(1,1:3), inv4b, dibdc_ani, d2ibdc2_ani )

    do n=1,3 
      do m=1,3
        do l=1,3
          do k=1,3
            cijkl(k,l,m,n) = d2ibdc2(k,l,m,n,1)*constant(1) +  &
              d2ibdc2(k,l,m,n,2)*constant(2) +  &
              2.d0*(dibdc(k,l,3)*dibdc(m,n,3)+  &
              (inv3b-1.d0)*d2ibdc2(k,l,m,n,3))*constant(3)+ &
              (2.d0*constant(4)+6.d0*(inv4b-1.d0)*constant(5))*dibdc_ani(k,l)*dibdc_ani(m,n)+ &
              (inv4b-1.d0)*(2.d0*constant(4)+3.d0*(inv4b-1.d0)*constant(5))*d2ibdc2_ani(k,l,m,n)
          enddo
        enddo
      enddo
    enddo
    cijkl(:,:,:,:)=4.d0*cijkl(:,:,:,:)

  end subroutine calElasticMooneyRivlinAniso

  !-------------------------------------------------------------------------------
  !> This subroutine provides to update stress and strain for anisotropic Mooney-Rivlin material
  !-------------------------------------------------------------------------------
  subroutine calUpdateElasticMooneyRivlinAniso( matl, sectType, strain, stress, cdsys, hdflag  )
    type( tMaterial ), intent(in) :: matl        !< material properties
    integer, intent(in)           :: sectType    !< not used currently
    real(kind=kreal), intent(out) :: stress(6)   !< 2nd Piola-Kirchhoff stress
    real(kind=kreal), intent(in)  :: strain(6)   !< Green-Lagrangen strain
    real(kind=kreal), intent(in)  :: cdsys(3,3) !> material coordinate system
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: k, l
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: inv1b, inv2b, inv3b, inv4b
    real(kind=kreal) :: dibdc(3,3,3)
    real(kind=kreal) :: d2ibdc2(3,3,3,3,3)
    real(kind=kreal) :: constant(5)
    real(kind=kreal) :: dudc(3,3)
    real(kind=kreal) :: dibdc_ani(3,3)
    real(kind=kreal) :: d2ibdc2_ani(3,3,3,3)

    constant(1:5)=matl%variables(M_PLCONST1:M_PLCONST5)
    if( present(hdflag) ) then
      if( hdflag == 1 ) then
        constant(3)=0.d0
      else if( hdflag == 2 ) then
        constant(1:2)=0.d0
        constant(4:5)=0.d0
      end if
    end if
    if( dabs(constant(3))>1.d-12 ) constant(3) = 1.d0/constant(3)

    call cderiv( matl, sectType, ctn, itn, inv1b, inv2b, inv3b,      &
      dibdc, d2ibdc2, strain, cdsys(1,1:3), inv4b, dibdc_ani, d2ibdc2_ani )

    ! ----- stress
    do l=1,3
      do k=1,3
        dudc(k,l) = dibdc(k,l,1)*constant(1)+dibdc(k,l,2)*constant(2)  &
          +2.d0*(inv3b-1.d0)*dibdc(k,l,3)*constant(3) &
          +(inv4b-1.d0)*(2.d0*constant(4)+3.d0*(inv4b-1.d0)*constant(5))*dibdc_ani(k,l)
      enddo
    enddo

    stress(1)=2.d0*dudc(1,1)
    stress(2)=2.d0*dudc(2,2)
    stress(3)=2.d0*dudc(3,3)
    stress(4)=2.d0*dudc(1,2)
    stress(5)=2.d0*dudc(2,3)
    stress(6)=2.d0*dudc(1,3)

  end subroutine calUpdateElasticMooneyRivlinAniso

end module
