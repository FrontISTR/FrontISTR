!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_array_util
  use hecmw_util
  implicit none

  private
  public :: hecmw_qsort_int_array
  public :: hecmw_uniq_int_array
  public :: hecmw_bsearch_int_array

contains

  recursive subroutine hecmw_qsort_int_array(array, istart, iend)
    implicit none
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: istart, iend
    integer(kind=kint) :: left, right, center
    integer(kind=kint) :: pivot, tmp
    if (istart >= iend) return
    center = (istart + iend) / 2
    pivot = array(center)
    left = istart
    right = iend
    do
      do while (array(left) < pivot)
        left = left + 1
      end do
      do while (pivot < array(right))
        right = right - 1
      end do
      if (left >= right) exit
      tmp = array(left)
      array(left) = array(right)
      array(right) = tmp
      left = left + 1
      right = right - 1
    end do
    if (istart < left-1) call hecmw_qsort_int_array(array, istart, left-1)
    if (right+1 < iend) call hecmw_qsort_int_array(array, right+1, iend)
  end subroutine hecmw_qsort_int_array

  subroutine hecmw_uniq_int_array(array, istart, iend, ndup)
    implicit none
    integer(kind=kint), intent(inout) :: array(:)
    integer(kind=kint), intent(in) :: istart, iend
    integer(kind=kint), intent(out) :: ndup
    integer(kind=kint) :: i
    ndup = 0
    do i = istart+1, iend
      if (array(i) == array(i - 1 - ndup)) then
        ndup = ndup + 1
      else if (ndup > 0) then
        array(i - ndup) = array(i)
      endif
    end do
  end subroutine hecmw_uniq_int_array

  subroutine hecmw_bsearch_int_array(array, istart, iend, val, idx)
    implicit none
    integer(kind=kint), intent(in) :: array(:)
    integer(kind=kint), intent(in) :: istart, iend
    integer(kind=kint), intent(in) :: val
    integer(kind=kint), intent(out) :: idx
    integer(kind=kint) :: center, left, right, pivot
    left = istart
    right = iend
    do
      if (left > right) then
        idx = -1
        exit
      end if
      center = (left + right) / 2
      pivot = array(center)
      if (val < pivot) then
        right = center - 1
        cycle
      else if (pivot < val) then
        left = center + 1
        cycle
      else ! if (pivot == val) then
        idx = center
        exit
      end if
    end do
  end subroutine hecmw_bsearch_int_array

end module hecmw_array_util
