!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2014/01/25                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kazuya Goto (PExProCS LLC)                     !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_matrix_reorder
  use hecmw_util
  implicit none

  private
  public :: hecmw_matrix_reorder_profile
  public :: hecmw_matrix_reorder_values
  public :: hecmw_matrix_reorder_vector
  public :: hecmw_matrix_reorder_back_vector
  public :: hecmw_matrix_reorder_renum_item

contains

  subroutine hecmw_matrix_reorder_profile(N, perm, iperm, &
       indexL, indexU, itemL, itemU, indexLp, indexUp, itemLp, itemUp)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: perm(:), iperm(:)
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    integer(kind=kint), intent(out) :: indexLp(0:), indexUp(0:)
    integer(kind=kint), intent(out) :: itemLp(:), itemUp(:)
    integer(kind=kint) :: cntL, cntU, inew, iold, j, jold, jnew
    cntL = 0
    cntU = 0
    indexLp(0) = 0
    indexUp(0) = 0
    do inew=1,N
      iold = perm(inew)
      ! original L
      do j = indexL(iold-1)+1, indexL(iold)
        jold = itemL(j)
        jnew = iperm(jold)
        if (jnew < inew) then
          cntL = cntL + 1
          itemLp(cntL) = jnew
        else
          cntU = cntU + 1
          itemUp(cntU) = jnew
        end if
      end do
      ! original U
      do j = indexU(iold-1)+1, indexU(iold)
        jold = itemU(j)
        if (jold > N) cycle
        jnew = iperm(jold)
        if (jnew < inew) then
          cntL = cntL + 1
          itemLp(cntL) = jnew
        else
          cntU = cntU + 1
          itemUp(cntU) = jnew
        end if
      end do
      indexLp(inew) = cntL
      indexUp(inew) = cntU
      call sort_int_array(itemLp, indexLp(inew-1)+1, indexLp(inew))
      call sort_int_array(itemUp, indexUp(inew-1)+1, indexUp(inew))
    end do
  end subroutine hecmw_matrix_reorder_profile

  subroutine hecmw_matrix_reorder_values(N, NDOF, perm, iperm, &
       indexL, indexU, itemL, itemU, AL, AU, D, &
       indexLp, indexUp, itemLp, itemUp, ALp, AUp, Dp)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: perm(:), iperm(:)
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    real(kind=kreal), intent(in) :: AL(:), AU(:), D(:)
    integer(kind=kint), intent(in) :: indexLp(0:), indexUp(0:)
    integer(kind=kint), intent(in) :: itemLp(:), itemUp(:)
    real(kind=kreal), intent(out) :: ALp(:), AUp(:), Dp(:)
    Dp = 0.d0
    ALp = 0.d0
    AUp = 0.d0
    ! call reorder_diag(N, NDOF, perm, D, Dp)
    ! call reorder_off_diag(N, NDOF, perm, &
    !    indexL, indexU, itemL, itemU, AL, AU, &
    !    indexLp, itemLp, ALp)
    ! call reorder_off_diag(N, NDOF, perm, &
    !    indexL, indexU, itemL, itemU, AL, AU, &
    !    indexUp, itemUp, AUp)
    call reorder_diag2(N, NDOF, iperm, D, Dp)
    call reorder_off_diag2(N, NDOF, iperm, &
       indexL, itemL, AL, &
       indexLp, indexUp, itemLp, itemUp, ALp, AUp)
    call reorder_off_diag2(N, NDOF, iperm, &
       indexU, itemU, AU, &
       indexLp, indexUp, itemLp, itemUp, ALp, AUp)
  end subroutine hecmw_matrix_reorder_values

  subroutine hecmw_matrix_reorder_vector(N, NDOF, perm, X, Xp)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: perm(:)
    real(kind=kreal), intent(in) :: X(:)
    real(kind=kreal), intent(out) :: Xp(:)
    integer(kind=kint) :: inew, iold, j0new, j0old, j
    do inew=1,N
      iold = perm(inew)
      j0new = (inew-1)*NDOF
      j0old = (iold-1)*NDOF
      do j=1,NDOF
        Xp(j0new + j) = X(j0old + j)
      end do
    end do
  end subroutine hecmw_matrix_reorder_vector

  subroutine hecmw_matrix_reorder_back_vector(N, NDOF, perm, Xp, X)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: perm(:)
    real(kind=kreal), intent(in) :: Xp(:)
    real(kind=kreal), intent(out) :: X(:)
    integer(kind=kint) :: inew, iold, j0new, j0old, j
    do inew=1,N
      iold = perm(inew)
      j0new = (inew-1)*NDOF
      j0old = (iold-1)*NDOF
      do j=1,NDOF
        X(j0old + j) = Xp(j0new + j)
      end do
    end do
  end subroutine hecmw_matrix_reorder_back_vector

  subroutine hecmw_matrix_reorder_renum_item(N, perm, indexXp, itemXp)
    implicit none
    integer(kind=kint), intent(in) :: N
    integer(kind=kint), intent(in) :: perm(:)
    integer(kind=kint), intent(in) :: indexXp(0:)
    integer(kind=kint), intent(inout) :: itemXp(:)
    integer(kind=kint) :: NPX, i
    NPX = indexXp(N)
    do i=1,NPX
      itemXp(i) = perm( itemXp(i) )
    end do
  end subroutine hecmw_matrix_reorder_renum_item

  subroutine reorder_diag(N, NDOF, perm, D, Dp)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: perm(:)
    real(kind=kreal), intent(in) :: D(:)
    real(kind=kreal), intent(out) :: Dp(:)
    integer(kind=kint) :: NDOF2, inew, iold, j0new, j0old, j
    NDOF2 = NDOF*NDOF
    ! diagonal
    do inew=1,N
      iold = perm(inew)
      j0new = (inew-1)*NDOF2
      j0old = (iold-1)*NDOF2
      do j=1,NDOF2
        Dp(j0new + j) = D(j0old + j)
      end do
    end do
  end subroutine reorder_diag

  subroutine reorder_diag2(N, NDOF, iperm, D, Dp)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: iperm(:)
    real(kind=kreal), intent(in) :: D(:)
    real(kind=kreal), intent(out) :: Dp(:)
    integer(kind=kint) :: NDOF2, inew, iold, j0new, j0old, j
    NDOF2 = NDOF*NDOF
    ! diagonal
    do iold=1,N
      inew = iperm(iold)
      j0old = (iold-1)*NDOF2
      j0new = (inew-1)*NDOF2
      do j=1,NDOF2
        Dp(j0new + j) = D(j0old + j)
      end do
    end do
  end subroutine reorder_diag2

  subroutine reorder_off_diag(N, NDOF, perm, &
       indexL, indexU, itemL, itemU, AL, AU, &
       indexXp, itemXp, AXp)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: perm(:)
    integer(kind=kint), intent(in) :: indexL(0:), indexU(0:)
    integer(kind=kint), intent(in) :: itemL(:), itemU(:)
    real(kind=kreal), intent(in) :: AL(:), AU(:)
    integer(kind=kint), intent(in) :: indexXp(0:)
    integer(kind=kint), intent(in) :: itemXp(:)
    real(kind=kreal), intent(out) :: AXp(:)
    integer(kind=kint) :: NDOF2, inew, iold
    integer(kind=kint) :: jsoldL, jeoldL, jsoldU, jeoldU
    integer(kind=kint) :: jnew, knew, kold, jold, l0new, l0old, l
    NDOF2 = NDOF*NDOF
    ! new L
    do inew=1,N
      iold = perm(inew)
      jsoldL = indexL(iold-1)+1
      jeoldL = indexL(iold)
      jsoldU = indexU(iold-1)+1
      jeoldU = indexU(iold)
      do jnew = indexXp(inew-1)+1, indexXp(inew)
        knew = itemXp(jnew)
        kold = perm(knew)
        if (kold < iold) then
          call bsearch_int_array(itemL, jsoldL, jeoldL, kold, jold)
          if (jold < 0) then
            write(0,*) 'DEBUG:: jold < 0 in reorder_off_diag'
            cycle
          end if
          l0new = (jnew-1)*NDOF2
          l0old = (jold-1)*NDOF2
          do l=1,NDOF2
            AXp(l0new + l) = AL(l0old + l)
          end do
        else
          call bsearch_int_array(itemU, jsoldU, jeoldU, kold, jold)
          if (jold < 0) then
            write(0,*) 'DEBUG:: jold < 0 in reorder_off_diag'
            cycle
          end if
          l0new = (jnew-1)*NDOF2
          l0old = (jold-1)*NDOF2
          do l=1,NDOF2
            AXp(l0new + l) = AU(l0old + l)
          end do
        end if
      end do
    end do
  end subroutine reorder_off_diag

  subroutine reorder_off_diag2(N, NDOF, iperm, &
       indexX, itemX, AX, &
       indexLp, indexUp, itemLp, itemUp, ALp, AUp)
    implicit none
    integer(kind=kint), intent(in) :: N, NDOF
    integer(kind=kint), intent(in) :: iperm(:)
    integer(kind=kint), intent(in) :: indexX(0:)
    integer(kind=kint), intent(in) :: itemX(:)
    real(kind=kreal), intent(in) :: AX(:)
    integer(kind=kint), intent(in) :: indexLp(0:), indexUp(0:)
    integer(kind=kint), intent(in) :: itemLp(:), itemUp(:)
    real(kind=kreal), intent(out) :: ALp(:), AUp(:)
    integer(kind=kint) :: NDOF2, iold, inew
    integer(kind=kint) :: jsnewL, jenewL, jsnewU, jenewU
    integer(kind=kint) :: jold, kold, knew, jnew, l0old, l0new, l
    NDOF2 = NDOF*NDOF
    ! new L
    do iold=1,N
      inew = iperm(iold)
      jsnewL = indexLp(inew-1)+1
      jenewL = indexLp(inew)
      jsnewU = indexUp(inew-1)+1
      jenewU = indexUp(inew)
      do jold = indexX(iold-1)+1, indexX(iold)
        kold = itemX(jold)
        if (kold > N) cycle
        knew = iperm(kold)
        if (knew < inew) then
          call bsearch_int_array(itemLp, jsnewL, jenewL, knew, jnew)
          if (jnew < 0) then
            write(0,*) 'ERROR:: jnew < 0 in reorder_off_diag2'
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
          l0old = (jold-1)*NDOF2
          l0new = (jnew-1)*NDOF2
          do l=1,NDOF2
            ALp(l0new + l) = AX(l0old + l)
          end do
        else
          call bsearch_int_array(itemUp, jsnewU, jenewU, knew, jnew)
          if (jnew < 0) then
            write(0,*) 'ERROR:: jnew < 0 in reorder_off_diag2'
            call hecmw_abort( hecmw_comm_get_comm() )
          end if
          l0old = (jold-1)*NDOF2
          l0new = (jnew-1)*NDOF2
          do l=1,NDOF2
            AUp(l0new + l) = AX(l0old + l)
          end do
        end if
      end do
    end do
  end subroutine reorder_off_diag2

  recursive subroutine sort_int_array(array, istart, iend)
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
    if (istart < left-1) call sort_int_array(array, istart, left-1)
    if (right+1 < iend) call sort_int_array(array, right+1, iend)
  end subroutine sort_int_array

  subroutine bsearch_int_array(array, istart, iend, val, idx)
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
  end subroutine bsearch_int_array

end module hecmw_matrix_reorder
