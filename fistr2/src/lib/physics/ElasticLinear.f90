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
!> \brief  This module provides function on elastic material
module m_ElasticLinear
  use mMaterial

  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)

  contains

!> Calculate elastic matrix
   SUBROUTINE calElasticMatrix( matl, sectType, D, temp  )
   TYPE( tMaterial ), INTENT(IN) :: matl       !> material properties
   INTEGER, INTENT(IN)           :: sectType   !> plane strain/stress or 3D 
   REAL(KIND=kreal), INTENT(OUT) :: D(:,:)     !> elastic matrix
   REAL(KIND=kreal), OPTIONAL    :: temp       !> temprature
   REAL(KIND=kreal) :: EE, PP, COEF1, COEF2, ina(1), outa(2)
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
   !  call fetch_TableData( MC_ISOELASTIC, matl%dict, outa, ierr )
	 EE = matl%variables(M_YOUNGS)
     PP = matl%variables(M_POISSON)
   endif
   
     SELECT CASE (sectType)
     CASE (D3)
       D(1,1)=EE*(1.0-PP)/(1.0-2.0*PP)/(1.0+PP)
       D(1,2)=EE*PP/(1.0-2.0*PP)/(1.0+PP)
       D(1,3)=D(1,2)
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=D(1,2)
       D(3,1)=D(1,3)
       D(3,2)=D(2,3)
       D(3,3)=D(1,1)
       D(4,4)=EE/(1.0+PP)*0.5
       D(5,5)=EE/(1.0+PP)*0.5
       D(6,6)=EE/(1.0+PP)*0.5
     CASE (PlaneStress)
       COEF1=EE/(1.0-PP*PP)
       COEF2=0.5*(1.0-PP)
       D(1,1)=COEF1
       D(1,2)=COEF1*PP
       D(1,3)=0.0
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=0.0
       D(3,1)=0.0
       D(3,2)=0.0
       D(3,3)=COEF1*COEF2
     CASE (PlaneStrain)
       COEF1=EE/((1.0+PP)*(1.0-2.0*PP))
       COEF2=EE/(2.0*(1.0+PP))
       D(1,1)=COEF1*(1.0-PP)
       D(1,2)=COEF1*PP
       D(1,3)=0.0
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=0.0
       D(3,1)=0.0
       D(3,2)=0.0
       D(3,3)=COEF2
     CASE (AxisSymetric)
       COEF1=EE*(1.0-PP)/((1.0+PP)*(1.0-2.0*PP))
       COEF2=(1.0-2.0*PP)/(2.0*(1.0-PP))
       D(1,1)=COEF1
       D(1,2)=COEF1*PP/(1.0-PP)
       D(1,3)=0.0
       D(1,4)=D(1,2)
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=0.0
       D(2,4)=D(1,2)
       D(3,1)=0.0
       D(3,2)=0.0
       D(3,3)=COEF1*COEF2
       D(3,4)=0.0
       D(4,1)=D(1,4)
       D(4,2)=D(2,4)
       D(4,3)=0.0
       D(4,4)=D(1,1)
     CASE (Shell)
     END SELECT

   END SUBROUTINE


end module m_ElasticLinear
