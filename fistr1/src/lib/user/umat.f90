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
!>  \brief   This subroutine read in used-defined material properties
!>  tangent 
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/01/12
!>  \version    0.00
!======================================================================!
module mUmat
implicit none

INTEGER, PARAMETER, PRIVATE :: kreal = kind(0.0d0)

contains

!> This subroutine calculates constitutive matrix
subroutine uMatlMatrix( mname, matl, strain, stress, fstat, D  &
    , dtime, ttime, temperature )
     CHARACTER(len=*), INTENT(IN)  :: mname     !< material name
     REAL(KIND=kreal), INTENT(IN)  :: matl(:)   !< material properties
     real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
     REAL(KIND=kreal), INTENT(IN)  :: stress(6) !< 2nd Piola-Kirchhiff stress tensor
     REAL(KIND=kreal), INTENT(IN)  :: fstat(:)  !< state variables 
     REAL(KIND=kreal), INTENT(OUT) :: D(:,:)    !< strain-stress relation
     REAL(KIND=kreal), INTENT(IN)  :: dtime     !< time increment
     REAL(KIND=kreal), INTENT(IN)  :: ttime     !< total time at the start of the current increment
	 REAL(KIND=kreal), optional    :: temperature !< temprature
     
end subroutine

!> This subroutine calculate strain and stress increment  
subroutine uUpdate(  mname, matl, strain, stress, fstat, dtime, ttime, temperature )
      character(len=*), intent(in)    :: mname      !< material name
      real(KIND=kreal), intent(in)    :: matl(:)    !< material properties
      real(kind=kreal), intent(in)    :: strain(6)  !< strain
      real(kind=kreal), intent(inout) :: stress(6)  !< 2nd Piola-Kirchhiff stress tensor
      real(kind=kreal), intent(inout) :: fstat(:)   !< state variables
      REAL(KIND=kreal), INTENT(IN)    :: dtime     !< time increment
      REAL(KIND=kreal), INTENT(IN)    :: ttime     !< total time at the start of the current increment
      real(KIND=kreal), optional      :: temperature !< temperature

end subroutine

end module
