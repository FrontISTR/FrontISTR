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
!>  \brief   This module provide functions for elastoplastic calculation
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/01/12
!>  \version    0.00
!======================================================================!
module m_ElastoPlastic
  use mMaterial
  use m_ElasticLinear

  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)

  contains

   !> This subroutine calculates elastoplastic constitutive relation
   SUBROUTINE calElastoPlasticMatrix( matl, sectType, stress, istat, extval, D, temperature )
     TYPE( tMaterial ), INTENT(IN) :: matl      !< material properties
     INTEGER, INTENT(IN)           :: sectType  !< not used currently
     REAL(KIND=kreal), INTENT(IN)  :: stress(6) !< stress
     REAL(KIND=kreal), INTENT(IN)  :: extval(:) !< plastic strain, back stress
     INTEGER, INTENT(IN)           :: istat     !< plastic state
     REAL(KIND=kreal), INTENT(OUT) :: D(:,:)    !< constitutive relation
     REAL(KIND=kreal), OPTIONAL    :: temperature   !> temprature

     INTEGER :: i,j,ytype
     logical :: kinematic
     REAL(KIND=kreal) :: dum, dj1(6), dj2(6), dj3(6), a(6), De(6,6)
     REAL(KIND=kreal) :: C1,C2,C3, back(6)
     REAL(KIND=kreal) :: J1,J2,J3, fai, sita, harden, khard, da(6), devia(6)

	 ytype = getYieldFunction( matl%mtype )
	 if( ytype==3 ) THEN
       call uElastoPlasticMatrix( matl, stress, istat, extval, D  )
       return
     endif
     IF( sectType /=D3 ) stop "Elastoplastic calculation support only Solid element currently"
     kinematic = isKinematicHarden( matl%mtype )
     khard = 0.d0
     if( kinematic ) then
       back(:) = extval(2:7)
       khard = calKinematicHarden( matl, extval(1) )
     endif
     if( present( temperature ) ) then
          call calElasticMatrix( matl, sectTYPE, De, temperature  )
     else
          call calElasticMatrix( matl, sectTYPE, De )
     endif
   !  call calElasticMatrix( matl, sectType, De  )

     J1 = (stress(1)+stress(2)+stress(3))
     devia(1:3) = stress(1:3)-J1/3.d0
     devia(4:6) = stress(4:6)
     if( kinematic ) devia = devia-back
     J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
          dot_product( devia(4:6), devia(4:6) )
 
     D(:,:) = De(:,:) 
     if( istat == 0 ) return   ! elastic state

     !derivative of J2
     dj2(1:3) = devia(1:3)
     dj2(4:6) = 2.d0*devia(4:6)
     dj2 = dj2/( 2.d0*dsqrt(j2) )
     if( present(temperature) ) then
       harden = calHardenCoeff( matl, extval(1), temperature )
     else
       harden = calHardenCoeff( matl, extval(1) )
     endif
	
     SELECT CASE (yType)
     CASE (0)       ! Mises or. Isotropic
       a(1:6) = dsqrt(3.d0) * dj2
     CASE (1)      ! Mohr-Coulomb
       fai = matl%variables(M_PLCONST3)
       J3 = devia(1)*devia(2)*devia(3)                    &
           +2.d0* devia(4)*devia(5)*devia(6)                     &
           -devia(6)*devia(2)*devia(6)                           &
           -devia(4)*devia(4)*devia(3)                           &
           -devia(1)*devia(5)*devia(5)
       sita = -3.d0*dsqrt(3.d0)*J3/( 2.d0*(J2**1.5d0) )
       if( dabs( dabs(sita)-1.d0 ) <1.d-8 ) THEN
          C1 = 0.d0
          C2 = dsqrt(3.d0)
          C3 = 0.d0
       else
         if( dabs(sita) >1.d0 ) stop "Math Error in Mohr-Coulomb calculation"
         sita = asin( sita )/3.d0
         C2 = cos(sita)*( 1.d0*tan(sita)*tan(3.d0*sita) + sin(fai)* &
             ( tan(3.d0*sita)-tan(sita )/dsqrt(3.d0) ) )
         C1 = sin(fai)/3.d0
         C3 = dsqrt(3.d0)*sin(sita)+cos(sita)*sin(fai)/(2.d0*J2*cos(3.d0*sita))
       endif
       ! deirivative of j1
       dj1(1:3) = 1.d0
       dj1(4:6) = 0.d0
       ! deirivative of j3
       dj3(1) = devia(2)*devia(3)-devia(5)*devia(5)+J2/3.d0
       dj3(2) = devia(1)*devia(3)-devia(6)*devia(6)+J2/3.d0
       dj3(3) = devia(1)*devia(2)-devia(4)*devia(4)+J2/3.d0
       dj3(4) = 2.d0*(devia(5)*devia(6)-devia(3)*devia(4))
       dj3(5) = 2.d0*(devia(4)*devia(6)-devia(1)*devia(5))
       dj3(6) = 2.d0*(devia(4)*devia(5)-devia(2)*devia(6))
       a(:) = C1*dj1 + C2*dj2 + C3*dj3
     CASE (2)      ! Drucker-Prager
       fai = matl%variables(M_PLCONST3)
       ! deirivative of j1
       dj1(1:3) = 1.d0
       dj1(4:6) = 0.d0
       a(:) = fai*dj1(:) + dj2(:)
     END SELECT

     da = matmul( de, a )
     dum = harden + khard+ dot_product( da, a )
     DO i=1,6
     DO j=1,6
       D(i,j) = De(i,j) - da(i)*da(j)/dum
     ENDDO
     ENDDO

   END SUBROUTINE

   !> This subrouitne calculate equivalent stress
   REAL(KIND=kreal) FUNCTION cal_equivalent_stress(matl, stress, extval)
     TYPE( tMaterial ), INTENT(IN) :: matl        !< material property
     REAL(kind=kreal), INTENT(IN)  :: stress(6)   !< stress
     REAL(KIND=kreal), INTENT(IN)  :: extval(:)   !< plastic strain, back stress

     INTEGER :: ytype
     logical :: kinematic
     REAL(kind=kreal) :: eqvs, sita, fai, J1,J2,J3, devia(6)
     REAL(KIND=kreal) :: back(6)
     kinematic = isKinematicHarden( matl%mtype )
     if( kinematic ) back(:) = extval(2:7)
	 
     ytype = getYieldFunction( matl%mtype )
     J1 = (stress(1)+stress(2)+stress(3))
     devia(1:3) = stress(1:3)-J1/3.d0
     devia(4:6) = stress(4:6)
     if( kinematic ) devia = devia-back
     J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
          dot_product( devia(4:6), devia(4:6) )

     SELECT CASE (yType)
     CASE (0)       ! Mises or. Isotropic
        eqvs = dsqrt( 3.d0*J2 )
     CASE (1)      ! Mohr-Coulomb
       fai = matl%variables(M_PLCONST1)
       J3 = devia(1)*devia(2)*devia(3)                    &
           +2.d0* devia(4)*devia(5)*devia(6)                     &
           -devia(6)*devia(2)*devia(6)                           &
           -devia(4)*devia(4)*devia(3)                           &
           -devia(1)*devia(5)*devia(5)
       sita = -3.d0*dsqrt(3.d0)*J3/( 2.d0*(J2**1.5d0) )
       if( dabs( dabs(sita)-1.d0 ) <1.d-8 ) sita=SIGN(1.d0, sita)
       if( dabs(sita) >1.d0 ) stop "Math Error in Mohr-Coulomb calculation"
       sita = asin( sita )/3.d0
       eqvs = (cos(sita)-sin(sita)*sin(fai)/dsqrt(3.d0))*dsqrt(J2)  &
             +J1*sin(fai)/3.d0
     CASE (2)      ! Drucker-Prager
       eqvs = dsqrt(J2) 
     CASE DEFAULT
       eqvs = -1.d0
     END SELECT

     cal_equivalent_stress = eqvs
   END FUNCTION
   
   !> This subrouitne calculate equivalent stress
   REAL(KIND=kreal) FUNCTION cal_mises_strain( strain )
     REAL(kind=kreal), INTENT(IN)  :: strain(6)        !< strain
	 cal_mises_strain = 2.d0*dot_product( strain(1:3), strain(1:3) )
	 cal_mises_strain = cal_mises_strain+ dot_product( strain(4:6), strain(4:6) )
	 cal_mises_strain = dsqrt( cal_mises_strain/3.d0 )
   END FUNCTION
      
   !> This function calcualtes hardening coefficient
   REAL(KIND=kreal) FUNCTION calHardenCoeff( matl, pstrain, temp )
     TYPE( tMaterial ), INTENT(IN)          :: matl    !< material property
     REAL(KIND=kreal), INTENT(IN)           :: pstrain !< plastic strain
     REAL(KIND=kreal), INTENT(IN), OPTIONAL :: temp !< temprature

     INTEGER :: i, nc, htype
     logical :: ierr
     REAL(KIND=kreal) :: s0, s1,s2, ef, ina(2)
	 
     calHardenCoeff = -1.d0
     htype = getHardenType( matl%mtype )
     SELECT CASE (htype)
     CASE (0)  ! Linear hardening
       calHardenCoeff = matl%variables(M_PLCONST2)
     CASE (1)  ! Multilinear approximation
       if( present(temp) ) then
          ina(1) = temp;  ina(2)=pstrain
          call fetch_TableGrad( MC_YIELD, ina, matl%dict, calHardenCoeff, ierr )
	!	  print *, ina, calHardenCoeff; pause
       else
         ina(1)=pstrain
         call fetch_TableGrad( MC_YIELD, ina(1:1), matl%dict, calHardenCoeff, ierr )
       endif
	 !  print *,1, calHardenCoeff; pause
     CASE (2)  ! Swift
       s0= matl%variables(M_PLCONST1)
       s1= matl%variables(M_PLCONST2)
       s2= matl%variables(M_PLCONST3)
       calHardenCoeff = s1*s2*( s0+pstrain )**(s2-1)
     CASE (3)  ! Ramberg-Osgood
       s0= matl%variables(M_PLCONST1)
       s1= matl%variables(M_PLCONST2)
       s2= matl%variables(M_PLCONST3)
       if( present(temp) ) then
         ef = calCurrYield( matl, pstrain, temp )
       else
         ef = calCurrYield( matl, pstrain )
       endif
       calHardenCoeff = s1*(ef/s1)**(1.d0-s2) /(s0*s2)
     CASE(4)   ! Prager
       calHardenCoeff = 0.d0
     CASE(5)   ! Prager+linear
       calHardenCoeff = matl%variables(M_PLCONST2)
     END SELECT
   END FUNCTION

   !> This function calcualtes kinematic hardening coefficient
   REAL(KIND=kreal) FUNCTION calKinematicHarden( matl, pstrain )
     TYPE( tMaterial ), INTENT(IN) :: matl    !< material property
     REAL(KIND=kreal), INTENT(IN)  :: pstrain !< plastic strain
	 
     INTEGER :: htype
     htype = getHardenType( matl%mtype )
     SELECT CASE (htype)
     CASE(4, 5)   ! Prager
       calKinematicHarden = matl%variables(M_PLCONST3)
     CASE DEFAULT
       calKinematicHarden = 0.d0
     END SELECT
   END FUNCTION
   
   !> This function calcualtes state of kinematic hardening
   REAL(KIND=kreal) FUNCTION calCurrKinematic( matl, pstrain )
     TYPE( tMaterial ), INTENT(IN) :: matl    !< material property
     REAL(KIND=kreal), INTENT(IN)  :: pstrain !< plastic strain
	 
     INTEGER :: htype
     htype = getHardenType( matl%mtype )
     SELECT CASE (htype)
     CASE(4, 5)   ! Prager
       calCurrKinematic = matl%variables(M_PLCONST3)*pstrain
     CASE DEFAULT
       calCurrKinematic = 0.d0
     END SELECT
   END FUNCTION
   
   !> This function calcualtes current yield stress
   REAL(KIND=kreal) FUNCTION calCurrYield( matl, pstrain, temp )
     TYPE( tMaterial ), INTENT(IN) :: matl    !< material property
     REAL(KIND=kreal), INTENT(IN)  :: pstrain !< plastic strain
     REAL(kind=kreal), INTENT(IN), optional :: temp  !< temperature

     INTEGER :: i, nc, htype
     REAL(KIND=kreal) :: s0, s1,s2, ina(2), outa(1)
	 logical :: ierr
     calCurrYield = -1.d0
     htype = getHardenType( matl%mtype )

     SELECT CASE (htype)
     CASE (0, 5)  ! Linear hardening, Linear+Parger hardening
       calCurrYield = matl%variables(M_PLCONST1)+matl%variables(M_PLCONST2)*pstrain
     CASE (1)  ! Multilinear approximation
       if( present(temp) ) then
         ina(1) = temp;  ina(2)=pstrain
         CALL fetch_TableData(MC_YIELD, matl%dict, outa, ierr, ina)
       else
         ina(1) = pstrain
         CALL fetch_TableData(MC_YIELD, matl%dict, outa, ierr, ina(1:1))
       endif
       if( ierr ) stop "Fail to get yield stress!"
       calCurrYield = outa(1)
     CASE (2)  ! Swift
       s0= matl%variables(M_PLCONST1)
       s1= matl%variables(M_PLCONST2)
       s2= matl%variables(M_PLCONST3)
       calCurrYield = s1*( s0+pstrain )**s2
     CASE (3)  ! Ramberg-Osgood
       s0= matl%variables(M_PLCONST1)
       s1= matl%variables(M_PLCONST2)
       s2= matl%variables(M_PLCONST3)
       if( pstrain<=s0 ) then
         calCurrYield = s1
       else
         calCurrYield = s1*( pstrain/s0 )**(1.d0/s2)
       endif
     CASE (4)  ! Parger hardening
       calCurrYield = matl%variables(M_PLCONST1)
     END SELECT
   END FUNCTION
   
   !> This function calcualtes yield state 
   REAL(KIND=kreal) FUNCTION calYieldFunc( matl, stress, extval, temp )
     TYPE( tMaterial ), INTENT(IN) :: matl        !< material property
     REAL(KIND=kreal), INTENT(IN)  :: stress(6)   !< stress
     REAL(KIND=kreal), INTENT(IN)  :: extval(:)   !< plastic strain, back stress
     REAL(kind=kreal), INTENT(IN), optional :: temp  !< temperature

     INTEGER :: ytype
     logical :: kinematic
     REAL(kind=kreal) :: eqvs, sita, eta, fai, J1,J2,J3, f, devia(6)
     REAL(KIND=kreal) :: pstrain, back(6)
	 
     kinematic = isKinematicHarden( matl%mtype )
     if( kinematic ) back(:) = extval(2:7)
	 
     pstrain = extval(1)
     ytype = getYieldFunction( matl%mtype )
     J1 = (stress(1)+stress(2)+stress(3))
     devia(1:3) = stress(1:3)-J1/3.d0
     devia(4:6) = stress(4:6)
     if( kinematic ) devia = devia-back
	 
     J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
          dot_product( devia(4:6), devia(4:6) )
     if( present(temp) ) then
       eqvs = calCurrYield( matl, pstrain, temp )
     else
       eqvs = calCurrYield( matl, pstrain )
     endif

     SELECT CASE (yType)
     CASE (0)       ! Mises or. Isotropic
        f = dsqrt( 3.d0*J2 ) - eqvs
     CASE (1)      ! Mohr-Coulomb
       fai = matl%variables(M_PLCONST3)
       J3 = devia(1)*devia(2)*devia(3)                    &
           +2.d0* devia(4)*devia(5)*devia(6)                     &
           -devia(2)*devia(6)*devia(6)                           &
           -devia(3)*devia(4)*devia(4)                           &
           -devia(1)*devia(5)*devia(5)
       sita = -3.d0*dsqrt(3.d0)*J3/( 2.d0*(J2**1.5d0) )
       if( dabs( dabs(sita)-1.d0 ) <1.d-8 ) sita=SIGN(1.d0, sita)
       if( dabs(sita) >1.d0 ) stop "Math Error in Mohr-Coulomb calculation"
       sita = asin( sita )/3.d0
       f = (cos(sita)-sin(sita)*sin(fai)/dsqrt(3.d0))*dsqrt(J2)  &
             +J1*sin(fai)/3.d0 - eqvs*cos(fai)
     CASE (2)      ! Drucker-Prager
       eta = matl%variables(M_PLCONST3)
       f = dsqrt(J2) + eta*J1 - eqvs*matl%variables(M_PLCONST4)
     END SELECT
	 
     calYieldFunc = f
   END FUNCTION

   !> This subroutine does backward-Euler return calculation   
   subroutine BackwardEuler( matl, stress, plstrain, istat, fstat, temp )
      use m_utilities, only : eigen3
      type( tMaterial ), intent(in)    :: matl        !< material properties
      real(kind=kreal), intent(inout)  :: stress(6)   !< trial->real stress
      real(kind=kreal), intent(inout)  :: plstrain    !< plastic strain till current substep
      integer, intent(inout)           :: istat       !< plastic state
      real(kind=kreal), intent(inout)  :: fstat(:)    !< plastic strain, back stress
      REAL(kind=kreal), INTENT(IN), optional :: temp  !< temperature

      real(kind=kreal), parameter :: tol =1.d-3
      integer, parameter          :: MAXITER = 5
      real(kind=kreal) :: dlambda, f, mat(3,3)
      integer :: i,j,ytype, maxp(1), minp(1), mm
      real(kind=kreal) :: youngs, poisson, pstrain, dum, ina(1), ee(2)
      real(kind=kreal) :: J1,J2,J3, H, KH, KK, dd, yd, G, K, devia(6)
      real(kind=kreal) :: prnstre(3), prnprj(3,3), tstre(3,3)
      real(kind=kreal) :: sita, fai, dep, trialprn(3)
      real(kind=kreal) :: a,b,siga,sigb,lamab(2),fab(2)
      real(kind=kreal) :: resi(2,2), invd(2,2)
      logical          :: right, kinematic, ierr
      REAL(KIND=kreal) :: ftrial, betan, back(6)

      ytype = getYieldFunction( matl%mtype )
      if( ytype==3 ) THEN
       call uBackwardEuler( matl, stress, istat, fstat )
       return
      endif
	  
      pstrain = fstat(1)
      if( present(temp) ) then
        f = calYieldFunc( matl, stress, fstat, temp )
      else
        f = calYieldFunc( matl, stress, fstat )
      endif
      if( dabs(f)<tol ) then  ! yielded
        istat = 1
        return
      elseif( f<0.d0 ) then   ! not yielded or unloading
        if( istat==0 ) return
        if( plstrain==0.d0 ) then
          istat =0
          return
        endif
      endif
      istat = 1           ! yielded	 
      KH = 0.d0; KK=0.d0; betan=0.d0; back(:)=0.d0
	  
      kinematic = isKinematicHarden( matl%mtype )
      if( kinematic ) then
        back(:) = fstat(2:7)
        betan = calCurrKinematic( matl, pstrain )
      endif
	  
      J1 = (stress(1)+stress(2)+stress(3))/3.d0
      devia(1:3) = stress(1:3)-J1
      devia(4:6) = stress(4:6)
      if( kinematic ) devia = devia-back
      yd = cal_equivalent_stress(matl, stress, fstat)

      if( present(temp) ) then	  
	    ina(1) = temp
        CALL fetch_TableData(MC_ISOELASTIC, matl%dict, ee, ierr, ina)
      else
        call fetch_TableData(MC_ISOELASTIC, matl%dict, ee, ierr )
      endif
      if( ierr ) then
         stop " fail to fetch young's modulus in elastoplastic calculation"
      else
        youngs = ee(1)
        poisson = ee(2)
      endif
      G = youngs/ ( 2.d0*(1.d0+poisson) )
      K = youngs/ ( 3.d0*(1.d0-2.d0*poisson) )
      dlambda = 0.d0
	  
      if( yType==0 ) then    ! Mises or. Isotropic
        do i=1,MAXITER
          if( present(temp) ) then
            H= calHardenCoeff( matl, pstrain+dlambda, temp )
          else
            H= calHardenCoeff( matl, pstrain+dlambda )
          endif
          if( kinematic ) then
            KH = calKinematicHarden( matl, pstrain+dlambda )
          endif
          dd= 3.d0*G+H+KH
          dlambda = dlambda+f/dd
          if( plstrain+dlambda<0.d0 ) then
            dlambda = -plstrain
            istat=0; exit
          endif
          if( present(temp) ) then
            dum = calCurrYield( matl, pstrain+dlambda, temp )
          else
            dum = calCurrYield( matl, pstrain+dlambda )
          endif
          if( kinematic ) then
            KK = calCurrKinematic( matl, pstrain+dlambda )
          endif
          f = yd-3.d0*G*dlambda-dum -(KK-betan)
          if( dabs(f)<tol*tol ) exit
        enddo
        pstrain = pstrain+dlambda
        plstrain = plstrain+dlambda
        if( kinematic ) then
          KK = calCurrKinematic( matl, pstrain )
          fstat(2:7) = back(:)+(KK-betan)*devia(:)/yd
        endif
        devia(:) = (1.d0-3.d0*dlambda*G/yd)*devia(:)
        stress(1:3) = devia(1:3)+J1
        stress(4:6) = devia(4:6)
        stress(:)= stress(:)+back(:)
      elseif(yType==1) then    ! Mohr-Coulomb
        fai = matl%variables(M_PLCONST3)
		
     !   do j=1,MAXITER
        J2 = 0.5d0* dot_product( devia(1:3), devia(1:3) ) +  &
          dot_product( devia(4:6), devia(4:6) )
        J3 = devia(1)*devia(2)*devia(3)                    &
           +2.d0* devia(4)*devia(5)*devia(6)                     &
           -devia(6)*devia(2)*devia(6)                           &
           -devia(4)*devia(4)*devia(3)                           &
           -devia(1)*devia(5)*devia(5)
        sita = -3.d0*dsqrt(3.d0)*J3/( 2.d0*(J2**1.5d0) )
        if( dabs( dabs(sita)-1.d0 ) <1.d-8 ) sita=SIGN(1.d0, sita)
        if( dabs(sita) >1.d0 ) stop "Math Error in Mohr-Coulomb calculation"
        sita = asin( sita )/3.d0
        do mm=1,6
          if( dabs(stress(mm))<1.d-10 ) stress(mm)=0.d0
        enddo
        call eigen3( stress, prnstre, prnprj )  
        trialprn = prnstre
        maxp = MAXLOC( prnstre )
        minp = MINLOC( prnstre )
        mm = 1
        if( maxp(1)==1 .or. minp(1)==1 ) mm =2
        if( maxp(1)==2 .or. minp(1)==2 ) mm =3
        do i=1,MAXITER
          if( present(temp) ) then
            H= calHardenCoeff( matl, pstrain, temp )
          else
            H= calHardenCoeff( matl, pstrain )
          endif
          dd= 4.d0*G*( 1.d0+sin(fai)*sin(sita)/3.d0 )+4.d0*K         &
             *sin(fai)*sin(sita)+4.d0*H*cos(fai)*cos(fai)
          dlambda = dlambda+f/dd
          if( plstrain + 2.d0*dlambda*cos(fai)<0.d0 ) then
            if( cos(fai)==0.d0 ) STOP "Math error in return mapping"
            dlambda = -0.5d0*plstrain/cos(fai)
            istat=0; exit
          endif
          dum = pstrain + 2.d0*dlambda*cos(fai)
          if( present(temp) ) then
            yd = calCurrYield( matl, dum, temp )
          else
            yd = calCurrYield( matl, dum )
          endif
          f = prnstre(maxp(1))-prnstre(minp(1))+                     &
            (prnstre(maxp(1))+prnstre(minp(1)))*sin(fai)-            &
            (4.d0*G*(1.d0+sin(fai)*sin(sita)/3.d0)+4.d0*K*sin(fai)   &
            *sin(sita))*dlambda-2.d0*yd*cos(fai)
          if( dabs(f)<tol ) exit
        enddo
        pstrain = pstrain + 2.d0*dlambda*cos(fai)
        plstrain = plstrain + 2.d0*dlambda*cos(fai)
        prnstre(maxp(1)) = prnstre(maxp(1))-(2.d0*G*(1.d0+sin(fai)/3.d0)  &
                 + 2.d0*K*sin(fai) )*dlambda
        prnstre(minp(1)) = prnstre(minp(1))+(2.d0*G*(1.d0-sin(fai)/3.d0)  &
                 - 2.d0*K*sin(fai) )*dlambda
        prnstre(mm) = prnstre(mm)+(4.d0*G/3.d0-2.d0*K)*sin(fai)*dlambda

        tstre(:,:) = 0.d0
        tstre(1,1)= prnstre(1); tstre(2,2)=prnstre(2); tstre(3,3)=prnstre(3)
        mat= matmul( prnprj, tstre )
		mat= matmul( mat, transpose(prnprj) )
        stress(1) = mat(1,1)
        stress(2) = mat(2,2)
        stress(3) = mat(3,3)
        stress(4) = mat(1,2)
        stress(5) = mat(2,3)
        stress(6) = mat(3,1)
      elseif(yType==2) then    ! Drucker-Prager
        fai = matl%variables(M_PLCONST3)
        dum = matl%variables(M_PLCONST4)
        do i=1,MAXITER
          if( present(temp) ) then
            H= calHardenCoeff( matl, pstrain, temp )
          else
            H= calHardenCoeff( matl, pstrain )
          endif
          dd= G+K*fai*fai+H*dum*dum
          dlambda = dlambda+f/dd
          if( plstrain+dum*dlambda<0.d0 ) then
            if( dum==0.d0 ) STOP "Math error in return mapping"
            dlambda = -plstrain/dum
            istat=0; exit
          endif
          if( present(temp) ) then
            f = calCurrYield( matl, pstrain+dum*dlambda, temp  )
          else
            f = calCurrYield( matl, pstrain+dum*dlambda  )
          endif
          f = yd-G*dlambda+fai*(J1-K*fai*dlambda)- dum*f
          if( dabs(f)<tol*tol ) exit
        enddo
        pstrain = pstrain+dum*dlambda 
        plstrain = plstrain+dum*dlambda 
        devia(:) = (1.d0-G*dlambda/yd)*devia(:)
        J1 = J1-K*fai*dlambda
        stress(1:3) = devia(1:3)+J1
        stress(4:6) = devia(4:6)
      end if

      fstat(1) = pstrain
   end subroutine BackwardEuler

   ! Clear elatoplastic state 
   subroutine updateEPState( gauss )
      use mMechGauss
      type(tGaussStatus), intent(inout) :: gauss  ! status of curr gauss point
      gauss%plstrain= 0.d0
   end subroutine
		 
end module m_ElastoPlastic
