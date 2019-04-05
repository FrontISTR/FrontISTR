!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This subroutine read in used-defined material properties
!>  tangent
!> This subroutine calculates elastoplastic constitutive relation
subroutine uElastoPlasticMatrix( matl, stress, istat, fstat, D  )
  use hecmw
  implicit none
  real(kind=kreal), intent(in)  :: matl(:)   !< material properties
  real(kind=kreal), intent(in)  :: stress(6) !< stress
  integer, intent(in)           :: istat     !< plastic state
  real(kind=kreal), intent(in)  :: fstat(:)  !< plastic strain, back stress
  real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation
end subroutine

!> This subroutine does backward-Euler return calculation
subroutine uBackwardEuler( matl, stress, istat, fstat )
  use hecmw
  implicit none
  real(kind=kreal), intent(in)    :: matl       !< material properties
  real(kind=kreal), intent(inout) :: stress(6)  !< stress
  integer, intent(inout)          :: istat      !< plastic state
  real(kind=kreal), intent(inout) :: fstat(:)   !< plastic strain, back stress
end subroutine
