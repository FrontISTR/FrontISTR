!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> MITC assumed-strain tying-point rules.
module MITC_Tying
  use elementInfo, only: fe_mitc3_shell, fe_mitc4_shell, fe_mitc9_shell
  implicit none

  private

  integer, parameter :: kreal = kind(0.0d0)
  real(kind=kreal), parameter :: s13 = sqrt(1.0d0/3.0d0)
  real(kind=kreal), parameter :: s35 = sqrt(3.0d0/5.0d0)

  real(kind=kreal), parameter :: mitc4_points(2, 4) = reshape((/ &
     0.0d0, -1.0d0, &
     1.0d0,  0.0d0, &
     0.0d0,  1.0d0, &
    -1.0d0,  0.0d0  &
    /), (/ 2, 4 /))

  real(kind=kreal), parameter :: mitc9_rr_points(2, 6) = reshape((/ &
    -s13, -s35, &
     s13, -s35, &
     s13,  s35, &
    -s13,  s35, &
     s13,  0.0d0, &
    -s13,  0.0d0  &
    /), (/ 2, 6 /))

  real(kind=kreal), parameter :: mitc9_ss_points(2, 6) = reshape((/ &
    -s35, -s13, &
     0.0d0, -s13, &
     s35, -s13, &
     s35,  s13, &
     0.0d0,  s13, &
    -s35,  s13  &
    /), (/ 2, 6 /))

  real(kind=kreal), parameter :: mitc9_rs_points(2, 4) = reshape((/ &
    -s13, -s13, &
     s13, -s13, &
     s13,  s13, &
    -s13,  s13  &
    /), (/ 2, 4 /))

  real(kind=kreal), parameter :: mitc3_points(2, 3) = reshape((/ &
    0.5d0, 0.0d0, &
    0.0d0, 0.5d0, &
    0.5d0, 0.5d0  &
    /), (/ 2, 3 /))

  !> Sign patterns used by the MITC9 interpolation polynomials.
  real(kind=kreal), parameter, public :: mitc9_xi_sign(6, 2) = reshape((/ &
    -1.0d0,  1.0d0,  1.0d0, -1.0d0,  1.0d0, -1.0d0, &
    -1.0d0,  0.0d0,  1.0d0,  1.0d0,  0.0d0, -1.0d0  &
    /), (/ 6, 2 /))

  real(kind=kreal), parameter, public :: mitc9_eta_sign(6, 2) = reshape((/ &
    -1.0d0, -1.0d0,  1.0d0,  1.0d0,  0.0d0,  0.0d0, &
    -1.0d0, -1.0d0, -1.0d0,  1.0d0,  1.0d0,  1.0d0  &
    /), (/ 6, 2 /))

  public :: NumOfTyingSets
  public :: NumOfTyingPoints
  public :: getTyingPoint

contains

  !> Number of tying-point sets used by an MITC shell element.
  integer function NumOfTyingSets(etype)
    integer, intent(in) :: etype

    select case (etype)
    case (fe_mitc3_shell, fe_mitc4_shell)
      NumOfTyingSets = 1
    case (fe_mitc9_shell)
      NumOfTyingSets = 3
    case default
      stop "Unsupported MITC shell element type"
    end select
  end function NumOfTyingSets

  !> Number of tying points in one tying set.
  integer function NumOfTyingPoints(etype, iset)
    integer, intent(in) :: etype, iset

    if (iset < 1 .or. iset > NumOfTyingSets(etype)) stop "Invalid MITC tying set"

    select case (etype)
    case (fe_mitc3_shell)
      NumOfTyingPoints = 3
    case (fe_mitc4_shell)
      NumOfTyingPoints = 4
    case (fe_mitc9_shell)
      select case (iset)
      case (1, 2)
        NumOfTyingPoints = 6
      case (3)
        NumOfTyingPoints = 4
      end select
    end select
  end function NumOfTyingPoints

  !> Natural coordinate of one MITC tying point.
  subroutine getTyingPoint(etype, iset, ip, pos)
    integer, intent(in) :: etype, iset, ip
    real(kind=kreal), intent(out) :: pos(2)

    if (iset < 1 .or. iset > NumOfTyingSets(etype)) stop "Invalid MITC tying set"
    if (ip < 1 .or. ip > NumOfTyingPoints(etype, iset)) stop "Invalid MITC tying point"

    select case (etype)
    case (fe_mitc3_shell)
      pos = mitc3_points(:, ip)
    case (fe_mitc4_shell)
      pos = mitc4_points(:, ip)
    case (fe_mitc9_shell)
      select case (iset)
      case (1)
        pos = mitc9_rr_points(:, ip)
      case (2)
        pos = mitc9_ss_points(:, ip)
      case (3)
        pos = mitc9_rs_points(:, ip)
      end select
    end select
  end subroutine getTyingPoint

end module MITC_Tying
