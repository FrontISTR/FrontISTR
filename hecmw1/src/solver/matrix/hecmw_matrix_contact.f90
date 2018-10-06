!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module hecmw_matrix_contact
  use hecmw_util
  implicit none

  private
  public :: hecmw_cmat_init
  public :: hecmw_cmat_finalize
  public :: hecmw_cmat_clear
  public :: hecmw_cmat_add
  public :: hecmw_cmat_pack
  public :: hecmw_cmat_multvec_add
  public :: hecmw_cmat_ass_bc
  public :: hecmw_cmat_n_newpair
  public :: hecmw_cmat_LU
  public :: hecmw_cmat_LU_free
  public :: hecmw_cmat_substitute

  integer(kind=kint), parameter :: CMAT_MAX_VAL_INIT = 128
  integer(kind=kint), parameter :: CMAT_MAX_VAL_GROW = 2

contains

  subroutine hecmw_cmat_init( cmat )
    type(hecmwST_matrix_contact) :: cmat

    nullify( cmat%pair )
    cmat%max_val = 0
    call hecmw_cmat_clear( cmat)
  end subroutine hecmw_cmat_init

  subroutine hecmw_cmat_finalize( cmat )
    type(hecmwST_matrix_contact) :: cmat

    if( cmat%max_val > 0 ) deallocate( cmat%pair )
    call hecmw_cmat_init( cmat )
  end subroutine hecmw_cmat_finalize

  subroutine hecmw_cmat_clear( cmat )
    type(hecmwST_matrix_contact) :: cmat

    cmat%n_val = 0
    cmat%checked = .true.
    cmat%sorted = .true.
    cmat%max_row = 0
    cmat%max_col = 0
  end subroutine hecmw_cmat_clear

  function compair_pair_by_index( p1, p2 )
    integer(kind=kint) :: compair_pair_by_index
    type(hecmwST_index_value_pair) :: p1, p2

    if( p1%i < p2%i ) then
      compair_pair_by_index = -1
    else if( p1%i > p2%i ) then
      compair_pair_by_index = 1
    else if( p1%j < p2%j ) then
      compair_pair_by_index = -1
    else if( p1%j > p2%j ) then
      compair_pair_by_index = 1
    else
      compair_pair_by_index = 0
    endif
  end function compair_pair_by_index

  subroutine cmat_resize( cmat, newlen )
    type(hecmwST_matrix_contact) :: cmat
    integer(kind=kint) :: newlen
    type(hecmwST_index_value_pair), allocatable :: temp(:)

    integer(kind=kint) :: i

    if( newlen <= cmat%max_val ) return

    if( cmat%max_val > 0 ) then
      allocate( temp( cmat%n_val ) )
      do i = 1, cmat%n_val
        temp(i)%i = cmat%pair(i)%i
        temp(i)%j = cmat%pair(i)%j
        temp(i)%val = cmat%pair(i)%val
      enddo
      deallocate( cmat%pair )

      allocate( cmat%pair( newlen ) )

      do i = 1, cmat%n_val
        cmat%pair(i)%i = temp(i)%i
        cmat%pair(i)%j = temp(i)%j
        cmat%pair(i)%val = temp(i)%val
      enddo
      deallocate( temp )
    endif

    cmat%max_val = newlen
  end subroutine cmat_resize

  subroutine cmat_grow( cmat )
    type(hecmwST_matrix_contact) :: cmat
    integer(kind=kint) :: newlen

    if( cmat%max_val == 0 ) then
      newlen = CMAT_MAX_VAL_INIT
    else
      newlen = cmat%max_val * CMAT_MAX_VAL_GROW
    endif
    call cmat_resize( cmat, newlen )
  end subroutine cmat_grow

  subroutine hecmw_cmat_add( cmat, i, j, val )
    type(hecmwST_matrix_contact) :: cmat
    integer(kind=kint) :: i
    integer(kind=kint) :: j
    real(kind=kreal), dimension(:,:) :: val
    integer(kind=kint) :: cmp

    if( cmat%n_val == cmat%max_val ) then
      call cmat_grow( cmat )
    endif

    cmat%pair( cmat%n_val+1 )%i = i
    cmat%pair( cmat%n_val+1 )%j = j
    cmat%pair( cmat%n_val+1 )%val = val

    if( cmat%n_val > 0 .and. cmat%sorted ) then
      cmp = compair_pair_by_index( cmat%pair( cmat%n_val ), cmat%pair( cmat%n_val+1 ) )
      if( cmp > 0 ) then
        cmat%checked = .false.
        cmat%sorted = .false.
      endif

      if( cmat%checked .and. cmp == 0 ) then
        cmat%checked = .false.
      endif
    endif

    cmat%n_val = cmat%n_val + 1

    if( cmat%max_row < i ) cmat%max_row = i
    if( cmat%max_col < j ) cmat%max_col = j
  end subroutine hecmw_cmat_add

  recursive subroutine sort_pair_by_index( pair, first, last )
    type(hecmwST_index_value_pair), pointer :: pair(:)
    integer(kind=kint) :: first
    integer(kind=kint) :: last
    integer(kind=kint) :: left, right
    type(hecmwST_index_value_pair) :: pivot, temp

    pivot = pair( (first + last) / 2 )
    left = first
    right = last
    do
      do while( compair_pair_by_index( pair(left), pivot ) < 0 )
        left = left + 1
      enddo
      do while( compair_pair_by_index( pivot, pair(right) ) < 0 )
        right = right - 1
      enddo
      if ( left >= right ) exit

      temp        = pair(left)
      pair(left)  = pair(right)
      pair(right) = temp

      left  = left  + 1
      right = right - 1
    enddo
    if( first < left - 1 ) call sort_pair_by_index( pair, first, left - 1 )
    if( right + 1 < last ) call sort_pair_by_index( pair, right + 1, last )
  end subroutine sort_pair_by_index

  subroutine hecmw_cmat_pack( cmat )
    type(hecmwST_matrix_contact) :: cmat
    integer(kind=kint) :: i
    integer(kind=kint) :: n_dup, cmp

    if( cmat%checked .or. cmat%n_val <= 1 ) return

    if( .not. cmat%sorted ) then
      call sort_pair_by_index( cmat%pair, 1, cmat%n_val )
      cmat%sorted = .true.
    endif

    n_dup = 0
    do i = 2,cmat%n_val
      cmp = compair_pair_by_index( cmat%pair( i-1-n_dup ), cmat%pair( i ) )
      if( cmp == 0 ) then
        n_dup = n_dup + 1
        cmat%pair( i-n_dup )%val = cmat%pair( i-n_dup )%val + cmat%pair( i )%val
      else
        cmat%pair( i-n_dup )%i   = cmat%pair( i )%i
        cmat%pair( i-n_dup )%j   = cmat%pair( i )%j
        cmat%pair( i-n_dup )%val = cmat%pair( i )%val
      endif
    enddo
    cmat%n_val = cmat%n_val - n_dup
    cmat%checked = .true.
  end subroutine hecmw_cmat_pack

  subroutine hecmw_cmat_multvec_add( cmat, X, Y, len )
    type(hecmwST_matrix_contact) :: cmat
    real(kind=kreal) :: X(:), Y(:)
    integer(kind=kint) :: len
    integer(kind=kint) :: i, ii, jj, k, l
    integer, parameter :: ndof = 3

    if( cmat%max_row > len .or. cmat%max_col > len ) then
      write(*,*) 'ERROR: cmat_multvec_add: vector too short'
      call hecmw_abort( hecmw_comm_get_comm())
    endif

    do i = 1, cmat%n_val
      ii = ndof * (cmat%pair(i)%i - 1)
      jj = ndof * (cmat%pair(i)%j - 1)
      do l = 1, ndof
        do k = 1, ndof
          Y(ii + k) =  &
            & Y(ii + k) + cmat%pair(i)%val(k, l) * X(jj + l)
        enddo
      enddo
    enddo
  end subroutine hecmw_cmat_multvec_add

  !< Assmble prescribed disp boundary condition into cmat
  subroutine hecmw_cmat_ass_bc(hecMAT, inode, idof, RHS )
    type (hecmwST_matrix) :: hecMAT
    integer(kind=kint)    :: inode, idof
    real(kind=kreal)      :: RHS
    integer(kind=kint) :: nval, i, cnode
    integer(kind=kint), parameter :: NDOF=3

    !  NDOF = hecMAT%NDOF
    if( NDOF < idof ) return

    !-- modify rhs
    if( RHS /= 0.d0 ) then
      do nval=1,hecMAT%cmat%n_val
        if( hecMAT%cmat%pair(nval)%j /= inode ) cycle
        cnode = hecMAT%cmat%pair(nval)%i
        if( cnode == inode ) then
          do i=1,NDOF
            if( i==idof ) cycle
            hecMAT%B(NDOF*(cnode-1)+i) = hecMAT%B(NDOF*(cnode-1)+i)      &
              - hecMAT%cmat%pair(nval)%val(i,idof)*RHS
          enddo
        else
          do i=1, NDOF
            hecMAT%B(NDOF*(cnode-1)+i) = hecMAT%B(NDOF*(cnode-1)+i)       &
              - hecMAT%cmat%pair(nval)%val(i,idof)*RHS
          enddo
        endif
      enddo
    endif

    !-- clear cmat
    do nval=1,hecMAT%cmat%n_val
      if( hecMAT%cmat%pair(nval)%i /= inode ) cycle

      hecMAT%cmat%pair(nval)%val(idof,:)=0.d0
    enddo

    do nval=1,hecMAT%cmat%n_val
      if( hecMAT%cmat%pair(nval)%j /= inode ) cycle

      hecMAT%cmat%pair(nval)%val(:,idof)=0.d0
    enddo

  end subroutine hecmw_cmat_ass_bc

  !< number of new pair introduced by cmat ( symm only )
  integer function hecmw_cmat_n_newpair(hecMAT, symm )
    type(hecmwST_matrix), intent(in) :: hecMAT
    integer, intent(in)              :: symm
    integer(kind=kint) :: nval, ind, jnd, k, m, mnd
    logical :: isfind

    hecmw_cmat_n_newpair = 0
    do nval=1,hecMAT%cmat%n_val
      ind = hecMAT%cmat%pair(nval)%i
      jnd = hecMAT%cmat%pair(nval)%j

      if( (symm==1) .and. (ind<jnd) ) cycle
      if( ind==jnd ) cycle
      isfind = .false.
      do k=1,hecMAT%NP
        if( ind/=k .and. jnd/=k ) cycle
        do m= hecMAT%indexL(k-1)+1, hecMAT%indexL(k)
          mnd= hecMAT%itemL(m)
          if( (ind==k .and. jnd==mnd) .or. (ind==mnd .and. jnd==k) ) then
            isfind = .true.
          endif
          if ( isfind ) exit
        enddo
        if( isfind ) exit
      enddo
      if( .not. isfind ) hecmw_cmat_n_newpair = hecmw_cmat_n_newpair+1
    enddo
  end function hecmw_cmat_n_newpair

  !< Assmble LU-part of cmat
  subroutine hecmw_cmat_LU( hecMAT )
    type (hecmwST_matrix) :: hecMAT
    integer(kind=kint)    :: i, j, k, l, countCAL, countCAU

    allocate (hecMAT%indexCL(0:hecMAT%NP), hecMAT%indexCU(0:hecMAT%NP))

    hecMAT%indexCL = 0
    hecMAT%indexCU = 0

    do i= 1, hecMAT%NP
      do j= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(j)%i == i .and. hecMAT%cmat%pair(j)%j < i) then
          hecMAT%indexCL(i) = hecMAT%indexCL(i) + 1
        endif
        if (hecMAT%cmat%pair(j)%i == i .and. hecMAT%cmat%pair(j)%j > i) then
          hecMAT%indexCU(i) = hecMAT%indexCU(i) + 1
        endif
      enddo
      hecMAT%indexCL(i) = hecMAT%indexCL(i) + hecMAT%indexCL(i-1)
      hecMAT%indexCU(i) = hecMAT%indexCU(i) + hecMAT%indexCU(i-1)
    enddo

    hecMAT%NPCL = hecMAT%indexCL(hecMAT%NP)
    hecMAT%NPCU = hecMAT%indexCU(hecMAT%NP)

    allocate (hecMAT%itemCL(hecMAT%NPCL), hecMAT%itemCU(hecMAT%NPCU))
    allocate (hecMAT%CAL(9*hecMAT%NPCL), hecMAT%CAU(9*hecMAT%NPCU))

    countCAL = 0
    countCAU = 0
    do i= 1, hecMAT%NP
      do j= 1, hecMAT%cmat%n_val
        if (hecMAT%cmat%pair(j)%i == i .and. hecMAT%cmat%pair(j)%j < i) then
          countCAL = countCAL + 1
          hecMAT%itemCL(countCAL) = hecMAT%cmat%pair(j)%j
          do k= 1, 3
            do l= 1, 3
              hecMAT%CAL(9*(countCAL-1)+3*(k-1)+l) = hecMAT%cmat%pair(j)%val(k,l)
            enddo
          enddo
        endif
        if (hecMAT%cmat%pair(j)%i == i .and. hecMAT%cmat%pair(j)%j > i) then
          countCAU = countCAU + 1
          hecMAT%itemCU(countCAU) = hecMAT%cmat%pair(j)%j
          do k= 1, 3
            do l= 1, 3
              hecMAT%CAU(9*(countCAU-1)+3*(k-1)+l) = hecMAT%cmat%pair(j)%val(k,l)
            enddo
          enddo
        endif
      enddo
    enddo
  end subroutine hecmw_cmat_LU

  !< free LU-part of cmat
  subroutine hecmw_cmat_LU_free( hecMAT )
    type (hecmwST_matrix) :: hecMAT

    deallocate (hecMAT%indexCL, hecMAT%itemCL, hecMAT%CAL)
    deallocate (hecMAT%indexCU, hecMAT%itemCU, hecMAT%CAU)
  end subroutine hecmw_cmat_LU_free

  subroutine hecmw_cmat_substitute( dest, src )
    implicit none
    type(hecmwST_matrix_contact) :: dest
    type(hecmwST_matrix_contact) :: src
    dest%n_val = src%n_val
    dest%max_val = src%max_val
    if (associated(src%pair)) dest%pair => src%pair
    dest%checked = src%checked
    dest%sorted = src%sorted
    dest%max_row = src%max_row
    dest%max_col = src%max_col
  end subroutine hecmw_cmat_substitute

end module hecmw_matrix_contact
