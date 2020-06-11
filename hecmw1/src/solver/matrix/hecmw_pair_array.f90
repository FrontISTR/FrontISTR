module hecmw_pair_array
  use hecmw_util

  private

  public :: hecmwST_pair_array
  public :: hecmw_pair_array_init
  public :: hecmw_pair_array_finalize
  public :: hecmw_pair_array_append
  public :: hecmw_pair_array_sort
  public :: hecmw_pair_array_find_id

  type hecmwST_pair
    integer(kind=kint) :: id
    integer(kind=kint) :: i1, i2
  end type hecmwST_pair

  type hecmwST_pair_array
    integer(kind=kint) :: num
    integer(kind=kint) :: max_num
    type (hecmwST_pair), pointer :: pairs(:) => null()
  end type hecmwST_pair_array

contains

  subroutine hecmw_pair_array_init(parray, max_num)
    implicit none
    type (hecmwST_pair_array), intent(inout) :: parray
    integer(kind=kint), intent(in) :: max_num
    !if (associated(parray%pairs)) deallocate(parray%pairs)
    allocate(parray%pairs(max_num))
    parray%max_num = max_num
    parray%num = 0
  end subroutine hecmw_pair_array_init

  subroutine hecmw_pair_array_finalize(parray)
    implicit none
    type (hecmwST_pair_array), intent(inout) :: parray
    if (associated(parray%pairs)) deallocate(parray%pairs)
    parray%max_num = 0
    parray%num = 0
  end subroutine hecmw_pair_array_finalize

  subroutine hecmw_pair_array_append(parray, id, i1, i2)
    implicit none
    type (hecmwST_pair_array), intent(inout) :: parray
    integer(kind=kint), intent(in) :: id, i1, i2
    if (parray%num >= parray%max_num) then
      stop 'ERROR: hecmw_pair_array_append: overflow'
    endif
    parray%num = parray%num + 1
    parray%pairs(parray%num)%id = id
    parray%pairs(parray%num)%i1 = i1
    parray%pairs(parray%num)%i2 = i2
  end subroutine hecmw_pair_array_append

  subroutine hecmw_pair_array_sort(parray)
    implicit none
    type (hecmwST_pair_array), intent(inout) :: parray
    call pairs_sort(parray%pairs, 1, parray%num)
  end subroutine hecmw_pair_array_sort

  function hecmw_pair_array_find_id(parray, i1, i2)
    implicit none
    integer(kind=kint) :: hecmw_pair_array_find_id
    type (hecmwST_pair_array), intent(inout) :: parray
    integer(kind=kint), intent(in) :: i1, i2
    type (hecmwST_pair) :: p
    integer(kind=kint) :: id
    p%i1 = i1
    p%i2 = i2
    call pairs_find(parray%pairs, 1, parray%num, p, id)
    hecmw_pair_array_find_id = id
  end function hecmw_pair_array_find_id

  function pairs_comp(p1, p2)
    implicit none
    integer(kind=kint) :: pairs_comp
    type (hecmwST_pair), intent(in) :: p1, p2
    if (p1%i1 < p2%i1) then
      pairs_comp = -1
    else if (p1%i1 > p2%i1) then
      pairs_comp = 1
    else
      if (p1%i2 < p2%i2) then
        pairs_comp = -1
      else if (p1%i2 > p2%i2) then
        pairs_comp = 1
      else
        pairs_comp = 0
      endif
    endif
  end function pairs_comp

  recursive subroutine pairs_sort(pairs, from, to)
    implicit none
    type (hecmwST_pair), pointer :: pairs(:)
    integer(kind=kint), intent(in) :: from, to
    integer(kind=kint) :: center, left, right
    type (hecmwST_pair) :: pivot, tmp
    if (from >= to) return
    center = (from + to) / 2
    pivot = pairs(center)
    left = from
    right = to
    do
      do while (pairs_comp(pairs(left), pivot) < 0)
        left = left + 1
      enddo
      do while (pairs_comp(pivot, pairs(right)) < 0)
        right = right - 1
      enddo
      if (left >= right) exit
      tmp = pairs(left)
      pairs(left) = pairs(right)
      pairs(right) = tmp
      left = left + 1
      right = right - 1
    enddo
    if (from < left-1) call pairs_sort(pairs, from, left-1)
    if (right+1 < to) call pairs_sort(pairs, right+1, to)
    return
  end subroutine pairs_sort

  recursive subroutine pairs_find(pairs, from, to, p, id)
    implicit none
    type (hecmwST_pair), pointer :: pairs(:)
    integer(kind=kint), intent(in) :: from, to
    type (hecmwST_pair), intent(in) :: p
    integer(kind=kint), intent(out) :: id
    integer(kind=kint) :: center, icomp
    if (from > to) then
      id = -1
      return
    endif
    center = (from + to) / 2
    icomp = pairs_comp(p, pairs(center))
    if (icomp < 0) then
      call pairs_find(pairs, from, center-1, p, id)
      return
    else if (icomp > 0) then
      call pairs_find(pairs, center+1, to, p, id)
      return
    else
      id = pairs(center)%id
      return
    endif
  end subroutine pairs_find

end module hecmw_pair_array
