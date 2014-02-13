!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by K. Satoh (Advancesoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides function on NEO-HOOKEAN hyperelastic material
!
!>  \author     K. Sato(Advancesoft)
!>  \date       2009/12/16
!>  \version    0.00
!!
!> \par neo-Hookean hyperelastic material
!>  - potential function : \f$ \Psi = \frac{\mu}{2} ( I_{c}-3 ) - \mu \ln J + \frac{\lambda}{2} ( \ln J )^2 \f$
!>  - 2nd Piola-Kirchhoff stress : \f$ \mathbf{S} = \mu ( \mathbf{I} - \mathbf{C}^{-1} ) + \lambda ( \ln J ) \mathbf{C}^{-1} \f$
!>  - elastic tangent coefficient : \f$ \mathbf{\mathsf{C}} = \lambda \mathbf{C}^{-1} \otimes \mathbf{C}^{-1} + 2 ( \mu - \lambda \ln J ) \mathbf{I} \f$
!
!======================================================================!
module m_ElasticNeoHooke
!
  use mMaterial
!
  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)
!
  contains
!
!
!
!-------------------------------------------------------------------------------
!> \brief This subroutine provides elastic tangent coefficient for neo-Hookean material
!-------------------------------------------------------------------------------
  subroutine calElasticNeoHooke( matl, sectType, cijkl, ftn )
!
    type( tMaterial ), intent(in) :: matl            !> material properties
    integer, intent(in)           :: sectType        !> not used currently
    real(kind=kreal), intent(out) :: cijkl(3,3,3,3)  !> \f$ \mathbf{\mathsf{C}} \f$ Elastic tanget coefficient
    real(kind=kreal), INTENT(IN)  :: ftn(:,:)        !> \f$ \mathbf{F}_{n+1}^{(k)} \f$ deformation gradient tensor ( \f$ \mathbf{F} = \mathbf{I} + \mathbf{u}_{n+1}^{(k)} \otimes \mathbf{\nabla}_{X} \f$ )
!
    integer :: i, j, k, l
    real(kind=kreal) :: ctn(3,3)     !< right Cauchy-Green deformation tensor
    real(kind=kreal) :: ctninv(3,3)  !< inversion of right Cauchy-Green deformation tensor
    real(kind=kreal) :: det, dum
    real(kind=kreal) :: cc(3,3,3,3), ii(3,3,3,3), delta(3,3)
    real(kind=kreal) :: ee, pp, lambda, mu
    real(kind=kreal) :: lnJ
!
! ----- Lame coefficients
    ee = matl%variables(M_YOUNGS)
    pp = matl%variables(M_POISSON)
    lambda = pp*ee/(1.d0+pp)/(1.d0-2.d0*pp)
    mu     = ee/2.d0/(1.d0+pp)
!
! ----- calculate the right Cauchy-Green deformation tensor C
    ctn(:,:)    = 0.d0
    do i=1,3
      do j=1,3
        do k=1,3
          ctn(i,j) = ctn(i,j) + ftn(k,i)*ftn(k,j)
        enddo
      enddo
    enddo
!
! ---- inverse of right Cauchy-Green deformation tensor
!    --- determinant of Jacobian
    det=ctn(1,1)*ctn(2,2)*ctn(3,3)                                             &
       +ctn(2,1)*ctn(3,2)*ctn(1,3)                                             &
       +ctn(3,1)*ctn(1,2)*ctn(2,3)                                             &
       -ctn(3,1)*ctn(2,2)*ctn(1,3)                                             &
       -ctn(2,1)*ctn(1,2)*ctn(3,3)                                             &
       -ctn(1,1)*ctn(3,2)*ctn(2,3)
    if( det==0.d0 ) stop "Math error in ctn(right CG deformation tensor) Determinaut==0.0"
! ----- calculate the \ln J
    lnJ = 0.5d0*log(det)
!
!    --- inversion of Jacobian
    ctninv(:,:) = 0.d0
    dum=1.d0/det
    ctninv(1,1)=dum*( ctn(2,2)*ctn(3,3)-ctn(3,2)*ctn(2,3) )
    ctninv(1,2)=dum*(-ctn(1,2)*ctn(3,3)+ctn(3,2)*ctn(1,3) )
    ctninv(1,3)=dum*( ctn(1,2)*ctn(2,3)-ctn(2,2)*ctn(1,3) )
    ctninv(2,1)=dum*(-ctn(2,1)*ctn(3,3)+ctn(3,1)*ctn(2,3) )
    ctninv(2,2)=dum*( ctn(1,1)*ctn(3,3)-ctn(3,1)*ctn(1,3) )
    ctninv(2,3)=dum*(-ctn(1,1)*ctn(2,3)+ctn(2,1)*ctn(1,3) )
    ctninv(3,1)=dum*( ctn(2,1)*ctn(3,2)-ctn(3,1)*ctn(2,2) )
    ctninv(3,2)=dum*(-ctn(1,1)*ctn(3,2)+ctn(3,1)*ctn(1,2) )
    ctninv(3,3)=dum*( ctn(1,1)*ctn(2,2)-ctn(2,1)*ctn(1,2) )
!
! ----- calculation the \bm{C}^{-1} \otimes \bm{C}^{-1}
    cc(:,:,:,:) = 0.d0
    do i=1,3
      do j=1,3
        do k=1,3
          do l=1,3
            cc(i,j,k,l) = cc(i,j,k,l) + ctninv(i,j)*ctninv(k,l)
          enddo
        enddo
      enddo
    enddo
!
! ----- calculation the \bm{\mathsf{I}}
    ii(:,:,:,:) = 0.d0
    do i=1,3
      do j=1,3
        do k=1,3
          do l=1,3
            ii(i,j,k,l) = 0.5d0*( ctninv(i,k)*ctninv(j,l)         &
                                 + ctninv(i,l)*ctninv(j,k) )
          enddo
        enddo
      enddo
    enddo
!
! ----- elastic tanget coefficients for neo-Hooke hyperelastic model
!    \bm{\mathsf{C}} = \lambda \bm{C}^{-1} \otimes \bm{C}^{-1} + 2 ( \mu - \lambda \ln J ) \bm{\mathsf{I}}
    cijkl(:,:,:,:) = 0.d0
    do i=1,3
      do j=1,3
        do k=1,3
          do l=1,3
            cijkl(i,j,k,l) = lambda*cc(i,j,k,l) + 2.d0*(mu-lambda*lnJ)*ii(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
!
  end subroutine calElasticNeoHooke
!
!
!
!-------------------------------------------------------------------------------
!> \brief This subroutine provides to update stress and strain for neo-Hookean material
!-------------------------------------------------------------------------------
  subroutine calUpdateElasticNeoHooke( matl, sectType, dstrain, dstress, ftn )
    type( tMaterial ), intent(in) :: matl           !> material properties
    integer, intent(in)           :: sectType       !> not used currently
    real(kind=kreal), intent(out) :: dstrain(6)     !> \f$ \mathbf{E} \f$ Cauchy-Green strain
    real(kind=kreal), intent(out) :: dstress(6)     !> \f$ \mathbf{S} \f$ 2nd Piola-Kirchhoff stress
    real(kind=kreal), intent(in)  :: ftn(:,:)       !> \f$ \mathbf{F}_{n+1}^{(k)} \f$ deformation gradient tensor ( \f$ \mathbf{F} = \mathbf{I} + \mathbf{u}_{n+1}^{(k)} \otimes \mathbf{\nabla}_{X} \f$ )
!
    integer :: i, j, k
    real(kind=kreal) :: ctn(3,3), itn(3,3)
    real(kind=kreal) :: ctninv(3,3)
    real(kind=kreal) :: det, dum
    real(kind=kreal) :: ee, pp, lambda, mu
    real(kind=kreal) :: lnJ
    real(kind=kreal) :: CGstrain(3,3), PKstress(3,3)
!
    real(kind=kreal) :: volJ
    real(kind=kreal) :: sf(3,3), cauchy(3,3), mises
!
! ----- Lame coefficients
    ee = matl%variables(M_YOUNGS)
    pp = matl%variables(M_POISSON)
    lambda = pp*ee/(1.d0+pp)/(1.d0-2.d0*pp)
    mu     = ee/2.d0/(1.d0+pp)
!
! ----- calculate the right Cauchy-Green deformation tensor C
    ctn(:,:)    = 0.d0
    do i=1,3
      do j=1,3
        do k=1,3
          ctn(i,j) = ctn(i,j) + ftn(k,i)*ftn(k,j)
        enddo
      enddo
    enddo
!
! ----- calculate the right Cauchy-Green deformation tensor C
    itn(:,:)   = 0.d0
    do i=1,3
      itn(i,i) = 1.d0
    enddo
!
! ---- inverse of right Cauchy-Green deformation tensor
!    --- determinant of Jacobian
    det=ctn(1,1)*ctn(2,2)*ctn(3,3)                                             &
       +ctn(2,1)*ctn(3,2)*ctn(1,3)                                             &
       +ctn(3,1)*ctn(1,2)*ctn(2,3)                                             &
       -ctn(3,1)*ctn(2,2)*ctn(1,3)                                             &
       -ctn(2,1)*ctn(1,2)*ctn(3,3)                                             &
       -ctn(1,1)*ctn(3,2)*ctn(2,3)
    if( det==0.d0 ) stop "Math error in ctn(right CG deformation tensor) Determinaut==0.0"
!    --- inversion of Jacobian
    ctninv(:,:)   = 0.d0
    dum=1.d0/det
    ctninv(1,1)=dum*( ctn(2,2)*ctn(3,3)-ctn(3,2)*ctn(2,3) )
    ctninv(1,2)=dum*(-ctn(1,2)*ctn(3,3)+ctn(3,2)*ctn(1,3) )
    ctninv(1,3)=dum*( ctn(1,2)*ctn(2,3)-ctn(2,2)*ctn(1,3) )
    ctninv(2,1)=dum*(-ctn(2,1)*ctn(3,3)+ctn(3,1)*ctn(2,3) )
    ctninv(2,2)=dum*( ctn(1,1)*ctn(3,3)-ctn(3,1)*ctn(1,3) )
    ctninv(2,3)=dum*(-ctn(1,1)*ctn(2,3)+ctn(2,1)*ctn(1,3) )
    ctninv(3,1)=dum*( ctn(2,1)*ctn(3,2)-ctn(3,1)*ctn(2,2) )
    ctninv(3,2)=dum*(-ctn(1,1)*ctn(3,2)+ctn(3,1)*ctn(1,2) )
    ctninv(3,3)=dum*( ctn(1,1)*ctn(2,2)-ctn(2,1)*ctn(1,2) )
!
! ----- calculate the \ln J
    lnJ = 0.5d0*log(det)
!
! ----- calculate the strain
    CGstrain(:,:) = 0.d0
    CGstrain(1:3,1:3) = 0.5d0*( ctn(1:3,1:3) - itn(1:3,1:3) )
!
! ----- calculate the stress
    PKstress(:,:) = 0.d0
    PKstress(1:3,1:3) = mu*( itn(1:3,1:3) - ctninv(1:3,1:3) ) + lambda*lnJ*ctninv(1:3,1:3)
!
    dstrain(:)    = 0.d0
    dstress(:)    = 0.d0
!
    dstrain(1) = CGstrain(1,1)
    dstrain(2) = CGstrain(2,2)
    dstrain(3) = CGstrain(3,3)
    dstrain(4) = 0.5d0*( CGstrain(1,2) + CGstrain(2,1) )
    dstrain(5) = 0.5d0*( CGstrain(2,3) + CGstrain(3,2) )
    dstrain(6) = 0.5d0*( CGstrain(1,3) + CGstrain(3,1) )
!
    dstress(1) = PKstress(1,1)
    dstress(2) = PKstress(2,2)
    dstress(3) = PKstress(3,3)
    dstress(4) = 0.5d0*( PKstress(1,2) + PKstress(2,1) )
    dstress(5) = 0.5d0*( PKstress(2,3) + PKstress(3,2) )
    dstress(6) = 0.5d0*( PKstress(1,3) + PKstress(3,1) )

  end subroutine

end module m_ElasticNeoHooke
