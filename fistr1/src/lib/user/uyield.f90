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
  public :: uElastoPlasticNumStatus
  public :: uElastoPlasticMatrix
  public :: uBackwardEuler

  ! system-defined material properties are saved in matl(1:100)
  integer, parameter :: M_YOUNGS   = 1
  integer, parameter :: M_POISSON  = 2

  ! user-defined material properties are saved in matl(101:200)

contains
  !> This function returns the number of real state variables
  integer(kind=kint) function uElastoPlasticNumStatus( matl )
    real(kind=kreal),   intent(in)    :: matl(:)   !< material properties

    uElastoPlasticNumStatus = 0
  end function uElastoPlasticNumStatus

  !> This subroutine calculates elastoplastic constitutive relation
  subroutine uElastoPlasticMatrix( matl, stress, istat, fstat, plstrain, D, temp, hdflag )
    real(kind=kreal),   intent(in)  :: matl(:)   !< material properties
    real(kind=kreal),   intent(in)  :: stress(6) !< stress
    integer(kind=kint), intent(in)  :: istat     !< integer state variable
    real(kind=kreal),   intent(in)  :: fstat(:)  !< real state variables
    real(kind=kreal),   intent(in)  :: plstrain  !< plastic strain at the beginning of current substep
    real(kind=kreal),   intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal),   intent(in)  :: temp      !> temperature
    integer(kind=kint), intent(in)  :: hdflag    !> return total(0), dev term only(1) or hyd term only(2)

  end subroutine uElastoPlasticMatrix

  !> This subroutine does backward-Euler return calculation
  subroutine uBackwardEuler( matl, stress, plstrain, istat, fstat, temp, hdflag )
    real(kind=kreal),   intent(in)    :: matl(:)   !< material properties
    real(kind=kreal),   intent(inout) :: stress(6) !< trial->real stress
    real(kind=kreal),   intent(in)    :: plstrain  !< plastic strain at the beginning of current substep
    integer(kind=kint), intent(inout) :: istat     !< integer state variable
    real(kind=kreal),   intent(inout) :: fstat(:)  !< real state variables
    real(kind=kreal),   intent(in)    :: temp      !< temperature
    integer(kind=kint), intent(in)    :: hdflag    !> return total(0), dev term only(1) or hyd term only(2)

  end subroutine uBackwardEuler

end module mUYield
