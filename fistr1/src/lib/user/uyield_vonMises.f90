!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This subroutine read in used-defined material properties
!>  tangent
module mUYield
  use hecmw
  implicit none

  private
  public :: uElastoPlasticNumStatus
  public :: uElastoPlasticMatrix
  public :: uBackwardEuler

  ! system-defined material properties are saved in matl(1:100)
  integer, parameter :: M_YOUNGS   = 1
  integer, parameter :: M_POISSON  = 2

  ! user-defined material properties are saved in matl(101:200)
  integer, parameter :: M_VM_SIGMA_Y0 = 101
  integer, parameter :: M_VM_HARDEN   = 102
  integer, parameter :: RM_VM_ELASTIC = 0
  integer, parameter :: RM_VM_PLASTIC = 1

  real(kind=kreal), parameter :: Id(6,6) = reshape( &
      & (/  2.d0/3.d0, -1.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
      &    -1.d0/3.d0,  2.d0/3.d0, -1.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
      &    -1.d0/3.d0, -1.d0/3.d0,  2.d0/3.d0,  0.d0,  0.d0,  0.d0,   &
      &          0.d0,       0.d0,       0.d0, 0.5d0,  0.d0,  0.d0,   &
      &          0.d0,       0.d0,       0.d0,  0.d0, 0.5d0,  0.d0,   &
      &          0.d0,       0.d0,       0.d0,  0.d0,  0.d0, 0.5d0/), &
      & (/6, 6/))

contains
  !> This function returns the number of real state variables
  integer(kind=kint) function uElastoPlasticNumStatus( matl )
    real(kind=kreal),   intent(in)    :: matl(:)   !< material properties

    uElastoPlasticNumStatus = 1
  end function uElastoPlasticNumStatus

  !> This subroutine calculates elastoplastic constitutive relation
  subroutine uElastoPlasticMatrix( matl, stress, istat, fstat, plstrain, D, temp, hdflag )
    real(kind=kreal),   intent(in)  :: matl(:)   !< material properties
    real(kind=kreal),   intent(in)  :: stress(6) !< stress
    integer(kind=kint), intent(in)  :: istat     !< integer state variable
    real(kind=kreal),   intent(in)  :: fstat(:)  !< real state variables
    real(kind=kreal),   intent(in)  :: plstrain  !< plastic strain at the beginning of current substep
    real(kind=kreal),   intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal),   intent(in)  :: temp      !> temperature
    integer(kind=kint), intent(in)  :: hdflag    !> return total(0), dev term only(1) or hyd term only(2)

    real(kind=kreal) :: youngs, poisson, sigma_y0, harden
    real(kind=kreal) :: K, G, p, sdev_norm
    real(kind=kreal) :: dlambda, eqvs, C1, C2, C3, dum
    real(kind=kreal) :: sdev(6), a(6)
    integer :: i, j

    youngs   = matl(M_YOUNGS)
    poisson  = matl(M_POISSON)

    ! isotropic elastic matrix
    K = youngs/(1.d0-2.d0*poisson)/3.d0
    G = youngs/(1.d0+poisson)*0.5d0
    D(:,:)=0.d0
    D(1,1)=K+(4.d0/3.d0)*G
    D(1,2)=K-(2.d0/3.d0)*G
    D(4,4)=G
    D(1,3)=D(1,2)
    D(2,1)=D(1,2)
    D(2,2)=D(1,1)
    D(2,3)=D(1,2)
    D(3,1)=D(1,3)
    D(3,2)=D(2,3)
    D(3,3)=D(1,1)
    D(5,5)=D(4,4)
    D(6,6)=D(4,4)

    if( istat == RM_VM_ELASTIC ) return

    youngs   = matl(M_YOUNGS)
    poisson  = matl(M_POISSON)
    sigma_y0 = matl(M_VM_SIGMA_Y0)
    harden   = matl(M_VM_HARDEN)

    K = youngs/(1.d0-2.d0*poisson)/3.d0
    G = youngs/(1.d0+poisson)*0.5d0

    p = (stress(1)+stress(2)+stress(3))/3.d0
    sdev(1:3) = stress(1:3)-p
    sdev(4:6) = stress(4:6)
    sdev_norm = sqrt(dot_product( sdev(1:3), sdev(1:3) ) +  &
        2.d0 * dot_product( sdev(4:6), sdev(4:6) ))

    ! Mises
    a(1:6) = sdev(1:6)/sdev_norm
    dlambda = fstat(1) - plstrain
    eqvs = sqrt( 3.d0/2.d0 ) * sdev_norm
    C3 = eqvs + 3.d0*G*dlambda !trial mises stress
    C1 = 6.d0*dlambda*G*G/C3
    dum = 3.d0*G+harden
    C2 = 6.d0*G*G*(dlambda/C3-1.d0/dum)

    do i=1,6
      do j=1,6
        D(i,j) = D(i,j) - C1*Id(i,j) + C2*a(i)*a(j)
      enddo
    enddo
  end subroutine uElastoPlasticMatrix

  !> This subroutine does backward-Euler return calculation
  subroutine uBackwardEuler( matl, stress, plstrain, istat, fstat, temp, hdflag )
    real(kind=kreal),   intent(in)    :: matl(:)   !< material properties
    real(kind=kreal),   intent(inout) :: stress(6) !< trial->real stress
    real(kind=kreal),   intent(in)    :: plstrain  !< plastic strain at the beginning of current substep
    integer(kind=kint), intent(inout) :: istat     !< integer state variable
    real(kind=kreal),   intent(inout) :: fstat(:)  !< real state variables
    real(kind=kreal),   intent(in)    :: temp      !< temperature
    integer(kind=kint), intent(in)    :: hdflag    !> return total(0), dev term only(1) or hyd term only(2)

    real(kind=kreal), parameter :: tol =1.d-6
    integer, parameter          :: MAXITER = 5
    real(kind=kreal) :: youngs, poisson, sigma_y0, harden
    real(kind=kreal) :: f, pstrain, p, sigma_y, G, K, dlambda, dd, eqvs, sdev_norm
    real(kind=kreal) :: sdev(6)
    integer :: i

    youngs   = matl(M_YOUNGS)
    poisson  = matl(M_POISSON)
    sigma_y0 = matl(M_VM_SIGMA_Y0)
    harden   = matl(M_VM_HARDEN)

    p = (stress(1)+stress(2)+stress(3))/3.d0
    sdev(1:3) = stress(1:3)-p
    sdev(4:6) = stress(4:6)
    sdev_norm = sqrt( dot_product( sdev(1:3), sdev(1:3) ) +  &
        2.d0 * dot_product( sdev(4:6), sdev(4:6) ))
    eqvs = sqrt( 3.d0/2.d0 ) * sdev_norm

    sigma_y = calCurrYield_mises(matl, plstrain)
    f = eqvs - sigma_y

    if( abs(f/sigma_y)<tol ) then  ! on yield surface
      istat = RM_VM_PLASTIC
      return
    elseif( f<0.d0 ) then  ! not yielded or unloading
      istat = RM_VM_ELASTIC
      return
    endif

    istat = RM_VM_PLASTIC  ! yielded

    if( youngs==0.d0 ) stop "YOUNG's ratio==0"
    G = youngs/ ( 2.d0*(1.d0+poisson) )
    K = youngs/ ( 3.d0*(1.d0-2.d0*poisson) )

    dlambda = 0.d0
    pstrain = plstrain
    do i=1,MAXITER
      dd= 3.d0*G+harden
      dlambda = dlambda+f/dd
      if( dlambda<0.d0 ) then
        dlambda = 0.d0
        exit
      endif
      pstrain = plstrain+dlambda
      sigma_y = calCurrYield_mises(matl, pstrain)
      f = eqvs - 3.d0*G*dlambda - sigma_y
      if( abs(f/sigma_y)<tol ) exit
      if( i==MAXITER ) then
        write(0,*) 'ERROR: uBackwardEuler_mises: convergence failure'
      endif
    enddo
    sdev(:) = (1.d0-3.d0*dlambda*G/eqvs)*sdev(:)
    stress(1:3) = sdev(1:3)+p
    stress(4:6) = sdev(4:6)

    fstat(1) = pstrain
  end subroutine uBackwardEuler

  real(kind=kreal) function calCurrYield_mises(matl, pstrain)
    real(kind=kreal), intent(in) :: matl(:)
    real(kind=kreal), intent(in) :: pstrain

    real(kind=kreal) :: sigma_y0, harden

    sigma_y0 = matl(M_VM_SIGMA_Y0)
    harden   = matl(M_VM_HARDEN)

    calCurrYield_mises = sigma_y0+harden*pstrain
  end function calCurrYield_mises

end module mUYield
