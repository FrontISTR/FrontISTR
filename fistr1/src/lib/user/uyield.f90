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

!> This subroutine calculates elastoplastic constitutive relation
subroutine uElastoPlasticMatrix( matl, stress, istat, fstat, D  )
     use hecmw
     implicit none
     REAL(KIND=kreal), INTENT(IN)  :: matl(:)   !< material properties
     REAL(KIND=kreal), INTENT(IN)  :: stress(6) !< stress
     INTEGER, INTENT(IN)           :: istat     !< plastic state
     REAL(KIND=kreal), INTENT(IN)  :: fstat(:)  !< plastic strain, back stress
     REAL(KIND=kreal), INTENT(OUT) :: D(:,:)    !< strain-stress relation
	
end subroutine

!> This subroutine does backward-Euler return calculation   
subroutine uBackwardEuler( matl, stress, istat, fstat )
      use hecmw
      implicit none
      REAL(KIND=kreal), intent(in)    :: matl       !< material properties
      real(kind=kreal), intent(inout) :: stress(6)  !< stress
      integer, intent(inout)          :: istat      !< plastic state
      real(kind=kreal), intent(inout) :: fstat(:)   !< plastic strain, back stress
end subroutine
