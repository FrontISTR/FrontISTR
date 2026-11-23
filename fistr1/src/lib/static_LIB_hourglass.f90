!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  Hourglass control of the one point elements 241 and 361.
!>
!> The stiffness form of Flanagan and Belytschko (IJNME 17, 1981), with its
!> coefficient fixed on the reference configuration as in Joldes, Wittek and
!> Miller (CNME 24, 2008). With the hourglass vectors
!>   gamma_Ia = h_Ia - g_I . sum_J X_J h_Ja,   q_a = sum_I gamma_Ia u_I,
!> the energy is k/2 sum_a q_a.q_a with k = 0.005 G0 V0 sum_I g_I.g_I. The force
!> vanishes for an affine displacement, and q_a rotates with a rigid rotation
!> because sum_I gamma_Ia X_I = 0. It is kept apart from the material states.
module m_static_LIB_hourglass
  use hecmw, only : kint, kreal
  use elementInfo
  implicit none

  private
  public :: mean_gradient, hourglass_shear_modulus
  public :: add_hourglass_stiffness, add_hourglass_force

  real(kind=kreal), parameter :: hourglass_ratio = 0.005d0
  !> Natural coordinates of the corner nodes: the first 4 are also those of 241.
  real(kind=kreal), parameter :: corner(3,8) = reshape( (/ &
    -1.d0,-1.d0,-1.d0,   1.d0,-1.d0,-1.d0,   1.d0, 1.d0,-1.d0,  -1.d0, 1.d0,-1.d0, &
    -1.d0,-1.d0, 1.d0,   1.d0,-1.d0, 1.d0,   1.d0, 1.d0, 1.d0,  -1.d0, 1.d0, 1.d0 /), (/3,8/) )

contains

  !> Volume (area for ndim=2) and mean gradient g_I of the shape functions over it.
  !> The rule of the element integrates both exactly.
  subroutine mean_gradient( etype, nn, ndim, ecoord, volume, gbar )
    integer(kind=kint), intent(in) :: etype, nn, ndim
    real(kind=kreal), intent(in)   :: ecoord(ndim,nn)
    real(kind=kreal), intent(out)  :: volume, gbar(nn,ndim)
    real(kind=kreal)   :: pos(3), det, wg, gderiv(nn,ndim)
    integer(kind=kint) :: LX

    volume = 0.d0
    gbar = 0.d0
    do LX = 1, NumOfQuadPoints( etype )
      call getQuadPoint( etype, LX, pos )
      call getGlobalDeriv( etype, nn, pos(1:ndim), ecoord, det, gderiv )
      wg = getWeight( etype, LX )*det
      volume = volume + wg
      gbar = gbar + wg*gderiv
    enddo
    gbar = gbar/volume
  end subroutine mean_gradient

  !> Shear modulus of the undeformed, unyielded material at the given temperature.
  real(kind=kreal) function hourglass_shear_modulus( material, temperature )
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    type(tMaterial), pointer     :: material
    real(kind=kreal), intent(in) :: temperature
    type(tGaussStatus) :: gauss
    real(kind=kreal)   :: D(6,6), axes(3,3)
    integer(kind=kint) :: isEp, i

    gauss%pMaterial => material
    call fstr_init_gauss( gauss )
    isEp = 0
    if( isElastoplastic(material%mtype) .or. isCreep(material%mtype) ) isEp = 1
    axes = 0.d0
    do i = 1, 3
      axes(i,i) = 1.d0
    enddo
    call MatlMatrix( gauss, D3, D, 0.d0, 0.d0, axes, temperature, isEp=isEp )
    hourglass_shear_modulus = ( D(4,4) + D(5,5) + D(6,6) )/3.d0
    call fstr_finalize_gauss( gauss )
  end function hourglass_shear_modulus

  !> Hourglass vectors gamma(:,1:nmode) and the coefficient k per unit G0.
  subroutine hourglass_vectors( etype, nn, ndim, ecoord, thick, gamma, nmode, k )
    integer(kind=kint), intent(in)  :: etype, nn, ndim
    real(kind=kreal), intent(in)    :: ecoord(ndim,nn), thick
    real(kind=kreal), intent(out)   :: gamma(nn,4), k
    integer(kind=kint), intent(out) :: nmode
    real(kind=kreal)   :: h(nn,4), gbar(nn,ndim), volume
    integer(kind=kint) :: a

    if( ndim == 2 ) then
      nmode = 1
      h(:,1) = corner(1,1:nn)*corner(2,1:nn)
    else
      nmode = 4
      h(:,1) = corner(2,1:nn)*corner(3,1:nn)
      h(:,2) = corner(3,1:nn)*corner(1,1:nn)
      h(:,3) = corner(1,1:nn)*corner(2,1:nn)
      h(:,4) = corner(1,1:nn)*corner(2,1:nn)*corner(3,1:nn)
    endif
    call mean_gradient( etype, nn, ndim, ecoord, volume, gbar )
    do a = 1, nmode
      gamma(:,a) = h(:,a) - matmul( gbar, matmul( ecoord, h(:,a) ) )
    enddo
    k = hourglass_ratio*thick*volume*sum( gbar*gbar )
  end subroutine hourglass_vectors

  !> K_IJ = k sum_a gamma_Ia gamma_Ja I
  subroutine add_hourglass_stiffness( etype, nn, ndim, ecoord, thick, shear, stiff )
    integer(kind=kint), intent(in) :: etype, nn, ndim
    real(kind=kreal), intent(in)   :: ecoord(ndim,nn), thick, shear
    real(kind=kreal), intent(inout) :: stiff(:,:)
    real(kind=kreal)   :: gamma(nn,4), k, kij
    integer(kind=kint) :: nmode, i, j, d

    call hourglass_vectors( etype, nn, ndim, ecoord, thick, gamma, nmode, k )
    do j = 1, nn
      do i = 1, nn
        kij = k*shear*dot_product( gamma(i,1:nmode), gamma(j,1:nmode) )
        do d = 1, ndim
          stiff(ndim*(i-1)+d, ndim*(j-1)+d) = stiff(ndim*(i-1)+d, ndim*(j-1)+d) + kij
        enddo
      enddo
    enddo
  end subroutine add_hourglass_stiffness

  !> f_I = k sum_a gamma_Ia q_a, with u the total displacement from the reference.
  subroutine add_hourglass_force( etype, nn, ndim, ecoord, thick, shear, u, qf )
    integer(kind=kint), intent(in) :: etype, nn, ndim
    real(kind=kreal), intent(in)   :: ecoord(ndim,nn), thick, shear, u(ndim,nn)
    real(kind=kreal), intent(inout) :: qf(:)
    real(kind=kreal)   :: gamma(nn,4), k, q(ndim)
    integer(kind=kint) :: nmode, a, i

    call hourglass_vectors( etype, nn, ndim, ecoord, thick, gamma, nmode, k )
    do a = 1, nmode
      q = matmul( u, gamma(:,a) )
      do i = 1, nn
        qf(ndim*(i-1)+1:ndim*i) = qf(ndim*(i-1)+1:ndim*i) + k*shear*gamma(i,a)*q
      enddo
    enddo
  end subroutine add_hourglass_force

end module m_static_LIB_hourglass
