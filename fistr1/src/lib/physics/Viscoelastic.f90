!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!>  \author                date              version 
!>  X.Yuan(Advancesoft)    2010/10/06        original
!>  X.Yuan                 2012/09/18        add trs function(WLF, ARRHENIUS)       
!>  
!
!======================================================================!
module mViscoElastic
	
  use mMaterial
  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)
  
  private :: hvisc, trs, trsinc
	
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
  
  !> Compute trs time increment
  real(kind=kreal) function trsinc(tn1,tn,mtype, mvar)
      real(kind=kreal), intent(in) :: tn1     !< temperature at time n+1
      real(kind=kreal), intent(in) :: tn      !< temperature at time n
      integer, intent(in)          :: mtype   !< material type
      real(kind=kreal), intent(in) :: mvar(:) !< material properties
	  
      real(kind=kreal)  ::  hsn1, hsn, asn1, asn

      if(mtype==VISCOELASTIC+2) then   ! Arrhenius

        hsn1 = mvar(2)*( 1.d0/(tn1-mvar(3))-1.d0/(mvar(1)-mvar(3)) )/mvar(4)
        hsn = mvar(2)*( 1.d0/(tn-mvar(3))-1.d0/(mvar(1)-mvar(3)) )/mvar(4)
        asn1 = dexp(-hsn1)
        asn = dexp(-hsn)
        
      else   ! WLF

        hsn1 = mvar(2)*( tn1-mvar(1) )/ (mvar(3)+tn1-mvar(1) )*dlog(10.d0)
        hsn = mvar(2)*( tn-mvar(1) )/ (mvar(3)+tn-mvar(1) )*dlog(10.d0)
        asn1 = dexp(hsn1)
        asn = dexp(hsn)

      endif
	  
      trsinc = (asn1-asn)/(hsn1-hsn)

  end function
  
  !> Compute trs time
  real(kind=kreal) function trs(tn,mtype, mvar)
      real(kind=kreal), intent(in) :: tn      !< temperature at current step
      integer, intent(in)          :: mtype   !< material type
      real(kind=kreal), intent(in) :: mvar(:) !< material properties
	  
      real(kind=kreal)  ::  hsn

      if(mtype==VISCOELASTIC+2) then   ! Arrhenius
        hsn = mvar(2)*( 1.d0/(tn-mvar(3))-1.d0/(mvar(1)-mvar(3)) )/mvar(4)
      else   ! WLF
        hsn = mvar(2)*( tn-mvar(1) )/ (mvar(3)+tn-mvar(1) )*dlog(10.d0)
      endif
	  
      trs = dexp(hsn)

  end function
  
!-------------------------------------------------------------------------------
!> This subroutine calculates tangent moduli for isotropic viscoelastic material
!-------------------------------------------------------------------------------
  subroutine calViscoelasticMatrix(matl, sectType, dt, D, temp)
      use m_ElasticLinear
      TYPE( tMaterial ), INTENT(IN) :: matl      !< material properties
      INTEGER, INTENT(IN)           :: sectType  !< not used currently
      REAL(KIND=kreal), INTENT(IN)  :: dt        !< time increment
      REAL(KIND=kreal), INTENT(OUT) :: D(:,:)    !< constitutive relation
      REAL(KIND=kreal), OPTIONAL    :: temp      !> temprature

      integer   i,j, n
      real(kind=kreal) :: G,Gg,K,Kg, gfac,exp_n,mu_0,mu_n,dq_n,dtau
      real(kind=kreal) :: ina(1), outa(2), EE, PP, ddt

      type(tTable), pointer  :: dicval
      logical :: ierr


      D(:,:)=0.d0
      if( dt==0.d0 ) then
        if( present( temp ) ) then
          call calElasticMatrix( matl, D3, D, temp  )
        else
          call calElasticMatrix( matl, D3, D )
        endif
        return
      endif	  

      ddt = dt
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
        if( matl%mtype>VISCOELASTIC ) ddt=trs(temp,matl%mtype, matl%variables)*dt
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
        dtau  = ddt/dicval%tbval(2,n)
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
  

  end subroutine
  
!-------------------------------------------------------------------------------
!> This subroutine provides to update stress for viscoelastic material
  subroutine UpdateViscoelastic(matl, sectType, eps, sig, vsig, dt, temp, tempn)
      TYPE( tMaterial ), INTENT(IN)   :: matl      !< material properties
      INTEGER, INTENT(IN)             :: sectType  !< not used currently
      REAL(KIND=kreal), INTENT(IN)    :: eps(6)    !< strain after this step
      REAL(KIND=kreal), INTENT(out)   :: sig(6)    !< stress
      REAL(KIND=kreal), INTENT(inout) :: vsig(:)   !< Visco stress components
      REAL(KIND=kreal), INTENT(IN)    :: dt        !< time increment
      REAL(KIND=kreal), OPTIONAL      :: temp      !> current temprature
      REAL(KIND=kreal), OPTIONAL      :: tempn     !> temperature at last step

      integer   i,j, n
      real(kind=kreal) :: G,Gg,K,Kg,Kth, exp_n,mu_0,mu_n,dq_n,dtau, theta
      real(kind=kreal) :: ina(1), outa(2), EE, PP, ddt

      real(kind=kreal) :: devstrain(6), en(6)
      type(tTable), pointer  :: dicval
      logical :: ierr

      ddt = dt
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
        if( matl%mtype>VISCOELASTIC ) ddt=trsinc(temp,tempn, matl%mtype, matl%variables)*dt
      else
        EE = matl%variables(M_YOUNGS)
        PP = matl%variables(M_POISSON)
      endif

!     Set elastic parameters for G (mu) and lambda

      G = EE/(2.d0*(1.d0 + PP))
      K = EE/(3.d0*(1.d0 - 2.d0*PP))
      
!     Compute volumetric strain and deviatoric components

      theta = (eps(1) + eps(2) + eps(3))/3.d0
      do i = 1,3
        devstrain(i) = eps(i) - theta
      end do 
      do i = 4,6
        devstrain(i) = eps(i)*0.5d0
      end do 

!     Set properties for integrating the q_i terms

      sig(:) = 0.0d0
      mu_0 = 0.0d0

      call fetch_Table( MC_VISCOELASTIC, matl%dict, dicval, ierr )
      if( ierr ) stop "Viscoelastic properties not defined"
	  
      en(:) = vsig(12*dicval%tbrow+1:12*dicval%tbrow+6)

      do n = 1,dicval%tbrow
        mu_n  = dicval%tbval(1,n)
        dtau  = ddt/dicval%tbval(2,n)
        exp_n = dexp(-dtau)

        dq_n = mu_n * hvisc(dtau,exp_n)
        mu_0 = mu_0 + mu_n

        do i = 1,6
          vsig((n-1)*12+i+6) = exp_n*vsig((n-1)*12+i) + dq_n*(devstrain(i) - en(i))
          sig(i)  = sig(i) + vsig((n-1)*12+i+6)
        end do 
      end do

!     Finish stress update

      mu_0 = 1.d0 - mu_0
      do i = 1,6
        sig(i) = 2.d0*G*(mu_0*devstrain(i) + sig(i))
      end do

!     Add elastic bulk term

      Kth = K*theta*3.0d0
      do i = 1,3
        sig(i) = sig(i) + Kth
      end do

  end subroutine
  
  !> Update viscoplastic state 
   subroutine updateViscoElasticState( gauss )
      use mMechGauss
      type(tGaussStatus), intent(inout) :: gauss  ! status of curr gauss point 
	  
      integer :: i, n, nrow
      real(kind=kreal) :: thetan, vstrain(6)
      nrow = fetch_TableRow( MC_VISCOELASTIC, gauss%pMaterial%dict )
      do n = 1,nrow
        do i = 1,6
          gauss%fstatus((n-1)*12+i) = gauss%fstatus((n-1)*12+i+6) 
        end do 
      end do
      vstrain = gauss%strain
      thetan = (vstrain(1) + vstrain(2) + vstrain(3))/3.d0
      do i = 1,3
        gauss%fstatus(nrow*12+i) = vstrain(i) - thetan
      end do 
      do i = 4,6
        gauss%fstatus(nrow*12+i) = vstrain(i)*0.5d0
      end do 
   end subroutine
      
end module
