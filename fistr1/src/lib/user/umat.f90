!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This subroutine read in used-defined material properties
!>  tangent
module mUmat
  implicit none

  integer, parameter, private :: kreal = kind(0.0d0)

contains

  !> This subroutine calculates constitutive matrix
  subroutine uMatlMatrix( mname, matl, strain, stress, fstat, D  &
      , dtime, ttime, temperature )
    character(len=*), intent(in)  :: mname     !< material name
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
    real(kind=kreal), intent(in)  :: stress(6) !< 2nd Piola-Kirchhiff stress tensor
    real(kind=kreal), intent(in)  :: fstat(:)  !< state variables
    real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal), intent(in)  :: dtime     !< time increment
    real(kind=kreal), intent(in)  :: ttime     !< total time at the start of the current increment
    real(kind=kreal), optional    :: temperature !< temprature

  end subroutine

  !> This subroutine calculate strain and stress increment
  subroutine uUpdate(  mname, matl, strain, stress, fstat, dtime, ttime, temperature )
    character(len=*), intent(in)    :: mname      !< material name
    real(kind=kreal), intent(in)    :: matl(:)    !< material properties
    real(kind=kreal), intent(in)    :: strain(6)  !< strain
    real(kind=kreal), intent(inout) :: stress(6)  !< 2nd Piola-Kirchhiff stress tensor
    real(kind=kreal), intent(inout) :: fstat(:)   !< state variables
    real(kind=kreal), intent(in)    :: dtime     !< time increment
    real(kind=kreal), intent(in)    :: ttime     !< total time at the start of the current increment
    real(kind=kreal), optional      :: temperature !< temperature

  end subroutine

end module
