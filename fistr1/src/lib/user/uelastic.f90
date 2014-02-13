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
!>  \brief   This module provides interface for elastic or hyperelastic 
!>  calculation
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/02/16
!>  \version    0.00
!======================================================================!
module mUElastic

implicit none

contains

!> This subroutine calculates constitutive relation
subroutine uElasticMatrix( matl, strain, D )
    use hecmw
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
    real(kind=kreal), intent(out) :: D(6,6)    !< constitutive matrix
   
! following examples of linear elasticicty  
    real(kind=kreal) :: EE, PP

    D(:,:)=0.d0

    EE = matl(1)
    PP = matl(2)
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
end subroutine

!> This subroutine calculate updated strain and stress   
subroutine uElasticUpdate( matl, strain, stress )
    use hecmw
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
	real(kind=kreal), intent(in)  :: strain(6) !< strain
    real(kind=kreal), intent(out) :: stress(6) !< stress
	
! following examples of linear elasticicty   
    real(kind=kreal) :: D(6,6)
	call uElasticMatrix( matl(:), strain, D )
	stress = matmul( D, strain )
end subroutine

end module

