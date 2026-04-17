!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module provides unified element quality check functions
!! using the elementInfo common API
module m_precheck_LIB_elements

  use hecmw
  use elementInfo

  implicit none

  private
  public :: precheck_calc_vol_asp

contains

  !> Calculate volume and edge lengths (for aspect ratio) for any supported element type
  subroutine precheck_calc_vol_asp(ic_type, nn, xx, yy, zz, thick, vol, almax, almin)
    integer(kind=kint), intent(in) :: ic_type, nn
    real(kind=kreal), intent(in) :: xx(:), yy(:), zz(:), thick
    real(kind=kreal), intent(out) :: vol, almax, almin

    call calc_volume(ic_type, nn, xx, yy, zz, thick, vol)
    call calc_edgelength(ic_type, nn, xx, yy, zz, almax, almin)
  end subroutine

  !> Calculate element volume using numerical integration via elementInfo API
  subroutine calc_volume(ic_type, nn, xx, yy, zz, thick, vol)
    integer(kind=kint), intent(in) :: ic_type, nn
    real(kind=kreal), intent(in) :: xx(:), yy(:), zz(:), thick
    real(kind=kreal), intent(out) :: vol

    integer :: ndim, nqp, ip, i
    real(kind=kreal) :: pos(3), wt, det
    real(kind=kreal) :: elecoord(3, nn)

    if (HECMW_is_etype_shell(ic_type)) then
      call calc_volume_shell(ic_type, nn, xx, yy, zz, thick, vol)
      return
    endif

    ndim = getSpaceDimension(ic_type)

    elecoord = 0.0d0
    do i = 1, nn
      elecoord(1, i) = xx(i)
      elecoord(2, i) = yy(i)
      if (ndim >= 3) elecoord(3, i) = zz(i)
    enddo

    nqp = NumOfQuadPoints(ic_type)
    vol = 0.0d0
    do ip = 1, nqp
      call getQuadPoint(ic_type, ip, pos)
      wt = getWeight(ic_type, ip)
      det = getDeterminant(ic_type, nn, pos, elecoord)
      vol = vol + wt * det
    enddo

    if (ndim == 2) vol = vol * thick
  end subroutine

  !> Calculate shell element volume as area * thickness
  !! Area is computed by integrating the magnitude of the cross product of tangent vectors
  subroutine calc_volume_shell(ic_type, nn, xx, yy, zz, thick, vol)
    integer(kind=kint), intent(in) :: ic_type, nn
    real(kind=kreal), intent(in) :: xx(:), yy(:), zz(:), thick
    real(kind=kreal), intent(out) :: vol

    integer :: nqp, ip, i
    real(kind=kreal) :: pos(3), wt
    real(kind=kreal) :: deriv(nn, 3)
    real(kind=kreal) :: g1(3), g2(3), g3(3), xsum, area

    nqp = NumOfQuadPoints(ic_type)
    area = 0.0d0
    do ip = 1, nqp
      call getQuadPoint(ic_type, ip, pos)
      wt = getWeight(ic_type, ip)
      call getShapeDeriv(ic_type, pos, deriv)

      g1 = 0.0d0
      g2 = 0.0d0
      do i = 1, nn
        g1(1) = g1(1) + deriv(i,1) * xx(i)
        g1(2) = g1(2) + deriv(i,1) * yy(i)
        g1(3) = g1(3) + deriv(i,1) * zz(i)
        g2(1) = g2(1) + deriv(i,2) * xx(i)
        g2(2) = g2(2) + deriv(i,2) * yy(i)
        g2(3) = g2(3) + deriv(i,2) * zz(i)
      enddo

      ! Cross product magnitude = surface area element
      g3(1) = g1(2)*g2(3) - g1(3)*g2(2)
      g3(2) = g1(3)*g2(1) - g1(1)*g2(3)
      g3(3) = g1(1)*g2(2) - g1(2)*g2(1)
      xsum = sqrt(g3(1)**2 + g3(2)**2 + g3(3)**2)

      area = area + wt * xsum
    enddo

    vol = area * thick
  end subroutine

  !> Calculate max/min edge lengths for aspect ratio computation
  !! 1st-order elements: vertex-to-vertex distance
  !! 2nd-order elements: vertex-to-midnode + midnode-to-vertex distance
  subroutine calc_edgelength(ic_type, nn, xx, yy, zz, almax, almin)
    integer(kind=kint), intent(in) :: ic_type, nn
    real(kind=kreal), intent(in) :: xx(:), yy(:), zz(:)
    real(kind=kreal), intent(out) :: almax, almin

    integer :: ndim

    ndim = getSpaceDimension(ic_type)
    if (HECMW_is_etype_shell(ic_type)) ndim = 3

    almax = 0.0d0
    almin = 1.0d20

    select case(ic_type)
    ! 3D 1st order
    case(341)
      call upd1(1,2); call upd1(2,3); call upd1(1,3)
      call upd1(1,4); call upd1(2,4); call upd1(3,4)
    case(351)
      call upd1(1,2); call upd1(2,3); call upd1(1,3)
      call upd1(4,5); call upd1(5,6); call upd1(4,6)
      call upd1(1,4); call upd1(2,5); call upd1(3,6)
    case(361)
      call upd1(1,2); call upd1(2,3); call upd1(3,4); call upd1(1,4)
      call upd1(5,6); call upd1(6,7); call upd1(7,8); call upd1(5,8)
      call upd1(1,5); call upd1(2,6); call upd1(3,7); call upd1(4,8)
    ! 3D 2nd order
    case(342)
      call upd2(1,5,2);  call upd2(2,6,3);  call upd2(3,7,1)
      call upd2(1,8,4);  call upd2(2,9,4);  call upd2(3,10,4)
    case(352)
      call upd2(1,7,2);  call upd2(2,8,3);  call upd2(1,9,3)
      call upd2(4,10,5); call upd2(5,11,6); call upd2(4,12,6)
      call upd2(1,13,4); call upd2(2,14,5); call upd2(3,15,6)
    case(362)
      call upd2(1,9,2);   call upd2(2,10,3);  call upd2(3,11,4);  call upd2(4,12,1)
      call upd2(5,13,6);  call upd2(6,14,7);  call upd2(7,15,8);  call upd2(8,16,5)
      call upd2(1,17,5);  call upd2(2,18,6);  call upd2(3,19,7);  call upd2(4,20,8)
    ! 2D 1st order
    case(231)
      call upd1(1,2); call upd1(2,3); call upd1(1,3)
    case(241)
      call upd1(1,2); call upd1(2,3); call upd1(3,4); call upd1(1,4)
    ! 2D 2nd order
    case(232)
      call upd2(1,4,2); call upd2(2,5,3); call upd2(3,6,1)
    case(242)
      call upd2(1,5,2); call upd2(2,6,3); call upd2(3,7,4); call upd2(4,8,1)
    ! Shell 1st order
    case(731)
      call upd1(1,2); call upd1(2,3); call upd1(1,3)
    case(741)
      call upd1(1,2); call upd1(2,3); call upd1(3,4); call upd1(1,4)
    case default
      almax = 0.0d0
      almin = 0.0d0
    end select

  contains

    subroutine upd1(n1, n2)
      integer, intent(in) :: n1, n2
      real(kind=kreal) :: al
      if (ndim == 2) then
        al = sqrt((xx(n2)-xx(n1))**2 + (yy(n2)-yy(n1))**2)
      else
        al = sqrt((xx(n2)-xx(n1))**2 + (yy(n2)-yy(n1))**2 + (zz(n2)-zz(n1))**2)
      endif
      if (al > almax) almax = al
      if (al < almin) almin = al
    end subroutine

    subroutine upd2(n1, nm, n2)
      integer, intent(in) :: n1, nm, n2
      real(kind=kreal) :: a1, a2, al
      if (ndim == 2) then
        a1 = sqrt((xx(nm)-xx(n1))**2 + (yy(nm)-yy(n1))**2)
        a2 = sqrt((xx(n2)-xx(nm))**2 + (yy(n2)-yy(nm))**2)
      else
        a1 = sqrt((xx(nm)-xx(n1))**2 + (yy(nm)-yy(n1))**2 + (zz(nm)-zz(n1))**2)
        a2 = sqrt((xx(n2)-xx(nm))**2 + (yy(n2)-yy(nm))**2 + (zz(n2)-zz(nm))**2)
      endif
      al = a1 + a2
      if (al > almax) almax = al
      if (al < almin) almin = al
    end subroutine

  end subroutine calc_edgelength

end module m_precheck_LIB_elements
