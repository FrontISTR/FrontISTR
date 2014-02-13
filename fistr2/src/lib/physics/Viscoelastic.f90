!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides functions for viscoelastic calculation
!
!>  \author     X.Yuan(Advancesoft)
!>  \date       2010/10/06
!>  \version    0.00
!
!======================================================================!
module mViscoElastic
	
  use mMaterial
  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)
  
  private :: hvisc
	
  contains

  !> Compute integration factor for viscoelastic material: H_visc(x) = [ 1 - exp(-x)]/x
  real(kind=kreal) function hvisc(x,expx)
      real(kind=kreal), intent(in) :: x       !< x
      real(kind=kreal), intent(in) :: expx    !< exp(-x)

      if(x < 1.d-04) then

        hvisc = 1.d0 - 0.5d0*x*(1.d0 - x/3.d0*(1.d0     &
                    - 0.25d0*x*(1.d0 - 0.2d0*x)))

      else

        hvisc = (1.d0 - expx)/x

      endif

  end function
  
!-------------------------------------------------------------------------------
!> This subroutine calculates tangent moduli for isotropic viscoelastic material
!-------------------------------------------------------------------------------
  subroutine calViscoelasticMatrix(matl, sectType, dt, D, temp)
      TYPE( tMaterial ), INTENT(IN) :: matl      !< material properties
      INTEGER, INTENT(IN)           :: sectType  !< not used currently
      REAL(KIND=kreal), INTENT(IN)  :: dt        !< time increment
      REAL(KIND=kreal), INTENT(OUT) :: D(:,:)    !< constitutive relation
      REAL(KIND=kreal), OPTIONAL    :: temp      !> temprature

      integer   i,j, n
      real(kind=kreal) :: G,Gg,K,Kg, gfac,exp_n,mu_0,mu_n,dq_n,dtau
      real(kind=kreal) :: ina(1), outa(2), EE, PP

      type(tTable)  :: dicval
      logical :: ierr


      D(:,:)=0.d0
      if( present(temp) ) then
        ina(1) = temp
        call fetch_TableData( MC_ISOELASTIC, matl%dict, outa, ierr, ina )
        if( ierr ) then
	      EE = matl%variables(M_YOUNGS)
          PP = matl%variables(M_POISSON)
        else
          EE = outa(1)
          PP = outa(2)
        endif
      else
        EE = matl%variables(M_YOUNGS)
        PP = matl%variables(M_POISSON)
      endif

!     Set elastic parameters for G (mu) and lambda

      G = EE/(2.d0*(1.d0 + PP))
      K = EE/(3.d0*(1.d0 - 2.d0*PP))


!     Set properties for integrating the q_i terms

      mu_0 = 0.0d0
      gfac = 0.0d0

      call fetch_Table( MC_VISCOELASTIC, matl%dict, dicval, ierr )
      if( ierr ) stop "Viscoelastic properties not defined"
	  
      do n = 1,dicval%tbrow
        mu_n  = dicval%tbval(1,n)
        dtau  = dt/dicval%tbval(2,n)
        exp_n = dexp(-dtau)

        dq_n = mu_n * hvisc(dtau,exp_n)
        gfac = gfac + dq_n
        mu_0 = mu_0 + mu_n

      end do

      mu_0 = 1.d0 - mu_0
      gfac = gfac + mu_0

!    Set tangent parameters

      Gg = G*gfac
      Kg = K - 0.6666666666666666d0*Gg
      do j =1,3
        do i = 1,3
          D(i,j) = Kg
        end do
        D(j,j) = D(j,j) + 2.d0*Gg
      end do 

      do i = 4,6
        D(i,i) = Gg
      end do 
	  
      call finalize_table( dicval )

  end subroutine
  
!-------------------------------------------------------------------------------
!> This subroutine provides to update stress for viscoelastic material
  subroutine UpdateViscoelastic(matl, sectType, eps, epsn, sig, vsig, dt, temp)
      TYPE( tMaterial ), INTENT(IN)   :: matl      !< material properties
      INTEGER, INTENT(IN)             :: sectType  !< not used currently
      REAL(KIND=kreal), INTENT(IN)    :: eps(6)    !< strain after this step
      REAL(KIND=kreal), INTENT(IN)    :: epsn(6)   !< strain before this step
      REAL(KIND=kreal), INTENT(out)   :: sig(6)    !< stress
      REAL(KIND=kreal), INTENT(inout) :: vsig(:)   !< Visco stress components
      REAL(KIND=kreal), INTENT(IN)    :: dt        !< time increment
      REAL(KIND=kreal), OPTIONAL      :: temp      !> temprature

      integer   i,j, n
      real(kind=kreal) :: G,Gg,K,Kg,Kth, gfac,exp_n,mu_0,mu_n,dq_n,dtau, theta, thetan
      real(kind=kreal) :: ina(1), outa(2), EE, PP

      real(kind=kreal) :: devstrain(6), en(6)
      type(tTable)  :: dicval
      logical :: ierr

!print *,vsig

      if( present(temp) ) then
        ina(1) = temp
        call fetch_TableData( MC_ISOELASTIC, matl%dict, outa, ierr, ina )
        if( ierr ) then
	      EE = matl%variables(M_YOUNGS)
          PP = matl%variables(M_POISSON)
        else
          EE = outa(1)
          PP = outa(2)
        endif
      else
        EE = matl%variables(M_YOUNGS)
        PP = matl%variables(M_POISSON)
      endif

!     Set elastic parameters for G (mu) and lambda

      G = EE/(2.d0*(1.d0 + PP))
      K = EE/(3.d0*(1.d0 - 2.d0*PP))

!     Compute volumetric strain and deviatoric components

      theta = (eps(1) + eps(2) + eps(3))/3.d0
      thetan = (epsn(1) + epsn(2) + epsn(3))/3.d0
      do i = 1,3
        devstrain(i) = eps(i) - theta
        en(i) = epsn(i) - thetan
      end do 
      do i = 4,6
        devstrain(i) = eps(i)*0.5d0
        en(i) = epsn(i)*0.5d0
      end do 

!     Set properties for integrating the q_i terms

      sig(:) = 0.0d0
      mu_0 = 0.0d0
      gfac = 0.0d0

      call fetch_Table( MC_VISCOELASTIC, matl%dict, dicval, ierr )
      if( ierr ) stop "Viscoelastic properties not defined"
	  
      do n = 1,dicval%tbrow
        mu_n  = dicval%tbval(1,n)
        dtau  = dt/dicval%tbval(2,n)
        exp_n = dexp(-dtau)

        dq_n = mu_n * hvisc(dtau,exp_n)
        gfac = gfac + dq_n
        mu_0 = mu_0 + mu_n

!       Update history and compute viscoelastic deviatoric stress

        do i = 1,6
          vsig((n-1)*6+i) = exp_n*vsig((n-1)*6+i) + dq_n*(devstrain(i) - en(i))
          sig(i)  = sig(i) + vsig((n-1)*6+i)
        end do 
      end do

!     Finish stress update

      mu_0 = 1.d0 - mu_0
      gfac = gfac + mu_0
      do i = 1,6
        sig(i) = 2.d0*G*(mu_0*devstrain(i) + sig(i))
      end do

!     Add elastic bulk term

      Kth = K*theta*3.0d0
      do i = 1,3
        sig(i) = sig(i) + Kth
      end do

      call finalize_table( dicval )
!print *,vsig; pause
  end subroutine
      
end module
