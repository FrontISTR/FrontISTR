!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module m_hecmw_matrix_ordering_CM
  use hecmw_util
  implicit none

  private
  public :: hecmw_matrix_ordering_CM
  public :: hecmw_matrix_ordering_RCM

  integer(kind=kint), parameter :: DEBUG = 0

contains

  subroutine hecmw_matrix_ordering_CM(N, indexL, itemL, indexU, itemU, &
      perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(out) :: perm(:), iperm(:)
    integer(kind=kint), parameter :: NMINMAX = 5
    integer(kind=kint) :: i, nmin
    integer(kind=kint) :: mins(NMINMAX)
    integer(kind=kint), allocatable :: degs(:)
    integer(kind=kint), allocatable :: nlevel(:)
    integer(kind=kint), allocatable :: lv_index(:,:), lv_item(:,:)
    integer(kind=kint) :: nlevel_max, max_id
    allocate(degs(N))
    call count_degrees(N, indexL, indexU, itemU, degs)
    call find_minimum_degrees(N, degs, NMINMAX, nmin, mins)
    allocate(nlevel(nmin), lv_index(0:N,nmin), lv_item(N,nmin))
    ! perform CM ordering starting from each minimum degree node
    !$omp parallel default(none),private(i), &
      !$omp&  shared(nmin,N,indexL,itemL,indexU,itemU,degs,mins,nlevel,lv_index,lv_item)
    !$omp do
    do i=1,nmin
      call ordering_CM_inner(N, indexL, itemL, indexU, itemU, degs, mins(i), &
        nlevel(i), lv_index(:,i), lv_item(:,i))
      if (DEBUG > 0) write(*,*) 'DEBUG:: hecmw_matrix_ordering_CM: i, nstart, nlevel = ', i, mins(i), nlevel(i)
    end do
    !$omp end do
    !$omp end parallel
    deallocate(degs)
    ! choose CM ordering with maximum nlevel
    nlevel_max = nlevel(1)
    max_id = 1
    do i=2,nmin
      if (nlevel(i) > nlevel_max) then
        nlevel_max = nlevel(i)
        max_id = i
      end if
    end do
    if (DEBUG > 0) write(*,*) 'DEBUG:: hecmw_matrix_ordering_CM: chose ordering',max_id
    do i=1,N
      perm(i) = lv_item(i,max_id)
      iperm(perm(i)) = i
    end do
    deallocate(nlevel, lv_index, lv_item)
  end subroutine hecmw_matrix_ordering_CM

  subroutine hecmw_matrix_ordering_RCM(N, indexL, itemL, indexU, itemU, &
      perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(out) :: perm(:), iperm(:)
    call hecmw_matrix_ordering_CM(N, indexL, itemL, indexU, itemU, perm, iperm)
    call reverse_ordering(N, perm, iperm)
    if (DEBUG > 0) then
      call write_nonzero_profile(N, indexL, itemL, indexU, itemU, perm, iperm)
      call write_perm(N, perm, iperm)
    endif
  end subroutine hecmw_matrix_ordering_RCM

  subroutine ordering_CM_inner(N, indexL, itemL, indexU, itemU, degs, nstart, &
      nlevel, lv_index, lv_item)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(in) :: degs(:)
    integer(kind=kint), intent(in) :: nstart
    integer(kind=kint), intent(out) :: nlevel
    integer(kind=kint), intent(out) :: lv_index(0:)
    integer(kind=kint), intent(out) :: lv_item(:)
    integer(kind=kint), allocatable :: iwk(:)
    integer(kind=kint) :: level, cntall, cnt, j, jnode, k, knode
    allocate(iwk(N))
    iwk(:) = 0
    lv_index(0) = 0
    ! first level
    iwk(nstart) = 1
    cntall = 1
    lv_item(1) = nstart
    lv_index(1) = 1
    ! other levels
    do level=2,N
      cnt = 0
      ! all nodes in previous level
      PRLV: do j = lv_index(level-2)+1, lv_index(level-1)
        jnode = lv_item(j)
        ! all connected nodes
        do k = indexL(jnode-1)+1, indexL(jnode)
          knode = itemL(k)
          if (iwk(knode) == 0) then
            iwk(knode) = level
            cnt = cnt + 1
            cntall = cntall + 1
            lv_item(cntall) = knode
            !if (cntall == N) exit PRLV
          end if
        end do
        do k = indexU(jnode-1)+1, indexU(jnode)
          knode = itemU(k)
          if (knode > N) cycle
          if (iwk(knode) == 0) then
            iwk(knode) = level
            cnt = cnt + 1
            cntall = cntall + 1
            lv_item(cntall) = knode
            !if (cntall == N) exit PRLV
          end if
        end do
      end do PRLV
      if (cnt == 0) then
        if (DEBUG > 0) write(*,*) 'DEBUG: choose any uncolored node..'
        do knode = 1, N
          if (iwk(knode) == 0) then
            iwk(knode) = level
            cnt = cnt + 1
            cntall = cntall + 1
            lv_item(cntall) = knode
            exit
          end if
        end do
      endif
      lv_index(level) = cntall
      call sort_nodes_by_degree(lv_item, lv_index(level-1)+1, lv_index(level), degs)
      if (cntall == N) then
        nlevel = level
        exit
      end if
    end do
    if (DEBUG > 0) then
      level = 0
      do j = 1, N
        jnode = lv_item(j)
        if (jnode <= 0 .or. jnode > N) stop 'ERROR: ordering_CM_inner: out of range'
        if (iwk(jnode) < 0) stop 'ERROR: ordering_CM_inner: duplicated node found'
        if (iwk(jnode) < level) stop 'ERROR: ordering_CM_inner: not in level order'
        level = iwk(jnode)
        iwk(jnode) = -1
      enddo
      do j = 1, N
        if (iwk(j) /= -1) stop 'ERROR: ordering_CM_inner: non-numbered node found'
      enddo
    endif
    deallocate(iwk)
  end subroutine ordering_CM_inner

  subroutine count_degrees(N, indexL, indexU, itemU, degs)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemU(:)
    integer(kind=kint), intent(out) :: degs(:)
    integer(kind=kint) :: i, j
    do i=1,N
      degs(i) = indexL(i) - indexL(i-1)
      do j = indexU(i-1)+1, indexU(i)
        if (itemU(j) <= N) degs(i) = degs(i) + 1
      end do
    end do
  end subroutine count_degrees

  subroutine find_minimum_degrees(N, degs, nminmax, nmin, mins)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: degs(:)
    integer(kind=kint), intent(in) :: nminmax
    integer(kind=kint), intent(out) :: nmin
    integer(kind=kint), intent(out) :: mins(nminmax)
    integer(kind=kint) :: degmin, i, j
    degmin = N
    nmin = 0
    do i=1,N
      ! skip unconnected nodes
      if (degs(i) == 0) cycle
      if (degs(i) < degmin) then
        degmin = degs(i)
        nmin = 1
        mins(1) = i
      else if (degs(i) == degmin) then
        nmin = nmin + 1
        if (nmin <= nminmax) mins(nmin) = i
      end if
    end do
    if (DEBUG > 0) write(*,*) 'DEBUG:: find_minimum_degrees: nmin, deg = ', nmin, degmin
    if (nmin > nminmax) nmin = nminmax
  end subroutine find_minimum_degrees

  recursive subroutine sort_nodes_by_degree(lv_item, istart, iend, degs)
    implicit none
    integer(kind=kint), intent(inout) :: lv_item(:)
    integer(kind=kint), intent(in) :: istart, iend
    integer(kind=kint), intent(in) :: degs(:)
    integer(kind=kint) :: pivot, left, right, itmp
    if (istart >= iend) return
    pivot = degs(lv_item((istart + iend) / 2))
    left = istart
    right = iend
    do
      do while (degs(lv_item(left)) < pivot)
        left = left + 1
      end do
      do while (pivot < degs(lv_item(right)))
        right = right - 1
      end do
      if (left >= right) exit
      itmp = lv_item(left)
      lv_item(left) = lv_item(right)
      lv_item(right) = itmp
      left = left + 1
      right = right - 1
    end do
    call sort_nodes_by_degree(lv_item, istart, left-1, degs)
    call sort_nodes_by_degree(lv_item, right+1, iend, degs)
  end subroutine sort_nodes_by_degree

  subroutine reverse_ordering(N, perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(inout) :: perm(:), iperm(:)
    integer(kind=kint) :: i, N1
    N1 = N + 1
    do i=1,N
      iperm(i) = N1 - iperm(i)
      perm(iperm(i)) = i
    end do
  end subroutine reverse_ordering

  subroutine write_nonzero_profile(N, indexL, itemL, indexU, itemU, perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(in) :: perm(:), iperm(:)
    integer(kind=kint), parameter :: F_ORG = 901
    integer(kind=kint), parameter :: F_NEW = 902
    integer(kind=kint) :: i, j, irow, jcol
    open(F_ORG, file='nzprof_org.txt', status='replace')
    do irow = 1, N
      i = irow
      do j = indexL(i-1)+1, indexL(i)
        jcol = itemL(j)
        write(F_ORG,*) irow, jcol
      end do
      do j = indexU(i-1)+1, indexU(i)
        jcol = itemU(j)
        if (jcol > N) cycle
        write(F_ORG,*) irow, jcol
      end do
    end do
    close(F_ORG)
    open(F_NEW, file='nzprof_new.txt', status='replace')
    do irow = 1, N
      i = perm(irow)
      do j = indexL(i-1)+1, indexL(i)
        jcol = itemL(j)
        write(F_NEW,*) irow, iperm(jcol)
      end do
      do j = indexU(i-1)+1, indexU(i)
        jcol = itemU(j)
        if (jcol > N) cycle
        write(F_NEW,*) irow, iperm(jcol)
      end do
    end do
    close(F_NEW)
  end subroutine write_nonzero_profile

  subroutine write_perm(N, perm, iperm)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: perm(:), iperm(:)
    integer(kind=kint), parameter :: F_PERM = 903
    integer(kind=kint) :: i
    open(F_PERM, file='perm_iperm.txt', status='replace')
    do i = 1, N
      write(F_PERM,*) i, perm(i), iperm(i)
    end do
    close(F_PERM)
  end subroutine write_perm

end module m_hecmw_matrix_ordering_CM
