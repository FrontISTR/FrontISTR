!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to set arrays for direct sparse solver
!> \in the case of using standard Lagrange multiplier algorithm for contact analysis.

module m_set_arrays_directsolver_contact

  use m_fstr
  use fstr_matrix_con_contact

  implicit none

  integer (kind=kint)               :: numNon0
  integer (kind=kint), allocatable :: pointers(:)        !< ia
  integer (kind=kint), allocatable :: indices(:)         !< ja
  real (kind=kreal)   , allocatable :: values(:)          !< a

  logical :: symmetricMatrixStruc

contains


  !> \brief This subroutine sets index arrays for direct sparse solver from those stored
  !> \in the matrix structures defined in MODULE fstr_matrix_con_contact
  subroutine set_pointersANDindices_directsolver(hecMAT,fstrMAT,is_sym)

    type(hecmwST_matrix)                     :: hecMAT    !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT   !< type fstrST_matrix_contact_lagrange
    logical :: is_sym

    integer (kind=kint)     :: np                         !< total number of nodes
    integer (kind=kint)     :: ndof                       !< degree of freedom
    integer (kind=kint)     :: num_lagrange               !< total number of Lagrange multipliers
    integer (kind=kint)     :: nn                         !< size of pointers
    integer (kind=kint)     :: ierr                       !< error indicateor
    integer (kind=kint)     :: i, j, k, l, countNon0

    symmetricMatrixStruc = is_sym

    np = hecMAT%NP ; ndof = hecMAT%NDOF ; num_lagrange = fstrMAT%num_lagrange
    nn = np*ndof + num_lagrange + 1

    if( symmetricMatrixStruc )then
      numNon0 = hecMAT%NPU*ndof**2+hecMAT%NP*ndof*(ndof+1)/2 &
        + (fstrMAT%numU_lagrange)*ndof+fstrMAT%num_lagrange
    else
      numNon0 = (hecMAT%NPL+hecMAT%NPU+hecMAT%NP)*ndof**2 &
        + (fstrMAT%numL_lagrange+fstrMAT%numU_lagrange)*ndof
    endif

    if(allocated(pointers))deallocate(pointers)
    allocate(pointers(nn), stat=ierr)
    if( ierr /= 0 ) stop " Allocation error, mkl%pointers "
    pointers = 0

    if(allocated(indices))deallocate(indices)
    allocate(indices(numNon0), stat=ierr)
    if( ierr /= 0 ) stop " Allocation error, mkl%indices "
    indices = 0

    pointers(1) = 1
    countNon0 = 1

    do i = 1, np
      do j = 1, ndof
        if( .not. symmetricMatrixStruc )then
          do l = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
            do k = 1, ndof
              indices(countNon0) = (hecMAT%itemL(l)-1)*ndof + k
              countNon0 = countNon0 + 1
            enddo
          enddo
          do k = 1, j-1
            indices(countNon0) = (i-1)*ndof + k
            countNon0 = countNon0 + 1
          enddo
        endif
        do k = j, ndof
          indices(countNon0) = (i-1)*ndof + k
          countNon0 = countNon0 + 1
        enddo
        do l = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
          do k = 1, ndof
            indices(countNon0) = (hecMAT%itemU(l)-1)*ndof + k
            countNon0 = countNon0 + 1
          enddo
        enddo
        if( num_lagrange > 0 )then
          do l = fstrMAT%indexU_lagrange(i-1)+1, fstrMAT%indexU_lagrange(i)
            indices(countNon0) = np*ndof + fstrMAT%itemU_lagrange(l)
            countNon0 = countNon0 + 1
          enddo
        endif
        pointers((i-1)*ndof+j+1) = countNon0
      enddo
    enddo

    if( num_lagrange > 0 )then
      do i = 1, num_lagrange
        if( symmetricMatrixStruc )then
          indices(countNon0) = np*ndof + i
          countNon0 = countNon0 + 1
        else
          do l = fstrMAT%indexL_lagrange(i-1)+1, fstrMAT%indexL_lagrange(i)
            do k = 1, ndof
              indices(countNon0) = (fstrMAT%itemL_lagrange(l)-1)*ndof + k
              countNon0 = countNon0 + 1
            enddo
          enddo
        endif
        pointers(np*ndof+i+1) = countNon0
      enddo
    endif

  end subroutine set_pointersANDindices_directsolver


  !> \brief This subroutine sets the array for direct sparse solver that contains
  !> \the non-zero items(elements)of stiffness matrix from those stored
  !> \in the matrix structures defined in MODULE fstr_matrix_con_contact
  subroutine set_values_directsolver(hecMAT,fstrMAT)

    type(hecmwST_matrix)                    :: hecMAT    !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)   :: fstrMAT   !< type fstrST_matrix_contact_lagrange

    integer (kind=kint)     :: np                        !< total number of nodes
    integer (kind=kint)     :: ndof                      !< degree of freedom
    integer (kind=kint)     :: num_lagrange              !< total number of Lagrange multipliers
    integer (kind=kint)     :: ierr                      !< error indicator
    integer (kind=kint)     :: i, j, k, l
    integer (kind=kint)     :: countNon0, locINal, locINd, locINau, locINal_lag, locINau_lag

    np = hecMAT%NP ; ndof = hecMAT%NDOF ; num_lagrange = fstrMAT%num_lagrange

    if(allocated(values))deallocate(values)
    allocate(values(numNon0), stat=ierr)
    if( ierr /= 0 ) stop " Allocation error, mkl%values "
    values = 0.0D0

    countNon0 = 1
    do i = 1, np
      do j = 1, ndof
        if( .not. symmetricMatrixStruc )then
          do l = hecMAT%indexL(i-1)+1, hecMAT%indexL(i)
            do k = 1, ndof
              locINal = ((l-1)*ndof+j-1)*ndof + k
              values(countNon0) = hecMAT%AL(locINal)
              countNon0 = countNon0 + 1
            enddo
          enddo
          do k = 1, j-1
            locINd = ((i-1)*ndof+j-1)*ndof + k
            values(countNon0) = hecMAT%D(locINd)
            countNon0 = countNon0 + 1
          enddo
        endif
        do k = j, ndof
          locINd = ((i-1)*ndof+j-1)*ndof + k
          values(countNon0) = hecMAT%D(locINd)
          countNon0 = countNon0 + 1
        enddo
        do l = hecMAT%indexU(i-1)+1, hecMAT%indexU(i)
          do k = 1, ndof
            locINau = ((l-1)*ndof+j-1)*ndof + k
            values(countNon0) = hecMAT%AU(locINau)
            countNon0 = countNon0 + 1
          enddo
        enddo
        if( num_lagrange > 0 )then
          do l = fstrMAT%indexU_lagrange(i-1)+1, fstrMAT%indexU_lagrange(i)
            locINau_lag = (l-1)*ndof + j
            values(countNon0) = fstrMAT%AU_lagrange(locINau_lag)
            countNon0 = countNon0 + 1
          enddo
        endif
      enddo
    enddo

    if( .not.symmetricMatrixStruc .and. num_lagrange > 0 )then
      do i = 1, num_lagrange
        do l = fstrMAT%indexL_lagrange(i-1)+1, fstrMAT%indexL_lagrange(i)
          do k = 1, ndof
            locINal_lag = (l-1)*ndof + k
            values(countNon0) = fstrMAT%AL_lagrange(locINal_lag)
            countNon0 = countNon0 + 1
          enddo
        enddo
      enddo
    endif

  end subroutine set_values_directsolver

  !> \brief This subroutine gets the residual vector
  subroutine getApproximateB(ntdf,x,y)

    integer(kind=kint)     :: ntdf                       !< total degree of freedom
    integer(kind=kint)     :: i, j, k
    real(kind=kreal)        :: x(ntdf)                    !< solution vector
    real(kind=kreal)        :: y(ntdf)                    !< residual vector

    y = 0.0d0
    do i = 1, ntdf
      do j = pointers(i), pointers(i+1)-1
        k = indices(j)
        y(i) = y(i) + values(j)*x(k)
        if( symmetricMatrixStruc .and. k/=i )&
          y(k)=y(k)+values(j)*x(i)
      enddo
    enddo

  end subroutine getApproximateB


  subroutine checkResidual(hecMAT,fstrMAT)
    type(hecmwST_matrix)                    :: hecMAT    !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)   :: fstrMAT   !< type fstrST_matrix_contact_lagrange
    integer(kind=kint)                       :: ntdf
    real(kind=kreal), allocatable            :: y(:)           !< right-hand side vector
    real(kind=kreal)                          :: residual_Max   !< maximum residual

    ntdf = hecMAT%NP*hecMAT%NDOF + fstrMAT%num_lagrange

    allocate(y(size(hecMAT%B)))
    y = 0.0d0
    call getApproximateB(ntdf,hecMAT%X,y)
    residual_Max=maxval(dabs(y-hecMAT%B))
    write(*,*)' maximum residual = ',residual_Max
    if( dabs(residual_Max) >= 1.0d-8) then
      write(*,*) ' ###Maximum residual exceeded 1.0d-8---Direct Solver### '
      !        stop
    endif
    deallocate(y)

  end subroutine checkResidual


end module m_set_arrays_directsolver_contact
