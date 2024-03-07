!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This subroutine read in used-defined material properties
!>  tangent
module mUYield
  use hecmw
  implicit none

  private
  public :: uElastoPlasticMatrix
  public :: uBackwardEuler

  ! system-defined material properties are saved in matl(1:100)
  integer, parameter :: M_YOUNGS   = 1
  integer, parameter :: M_POISSON  = 2

  ! user-defined material properties are saved in matl(101:200)

contains
  !> This subroutine calculates elastoplastic constitutive relation
  subroutine uElastoPlasticMatrix( matl, stress, istat, fstat, plstrain, D  )
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: stress(6) !< stress
    integer, intent(in)           :: istat     !< plastic state
    real(kind=kreal), intent(in)  :: fstat(:)  !< plastic strain, back stress
    real(kind=kreal), intent(in)  :: plstrain  !< plastic strain
    real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation

  end subroutine uElastoPlasticMatrix

  !> This subroutine does backward-Euler return calculation
  subroutine uBackwardEuler( matl, stress, plstrain, istat, fstat )
    real(kind=kreal), intent(in)    :: matl(:)    !< material properties
    real(kind=kreal), intent(inout) :: stress(6)  !< trial->real stress
    real(kind=kreal), intent(in)    :: plstrain   !< plastic strain till current substep
    integer, intent(inout)          :: istat      !< plastic state
    real(kind=kreal), intent(inout) :: fstat(:)   !< plastic strain, back stress

  end subroutine uBackwardEuler

end module mUYield
