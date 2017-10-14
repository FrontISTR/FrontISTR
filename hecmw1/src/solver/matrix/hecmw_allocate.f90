!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module hecmw_allocate
contains

  subroutine hecmw_allocate_matrix(hecMAT, mat, nBlock)
    use hecmw_util

    implicit none

    type (hecmwST_matrix):: hecMAT, mat
    integer:: nBlock
    integer:: ierr

    allocate( mat%AL( nBlock*hecMAT%NPL ), stat=ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'Allocation error: Not enough memory for matrix!'
      stop
    end if
    allocate( mat%AU( nBlock*hecMAT%NPU ), stat=ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'Allocation error: Not enough memory for matrix!'
      stop
    end if
    allocate( mat%D ( nBlock*hecMAT%NP  ), stat=ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'Allocation error: Not enough memory for matrix!'
      stop
    end if

    MAT%AL = 0.0d0
    MAT%AU = 0.0d0
    MAT%D  = 0.0d0

  end subroutine hecmw_allocate_matrix

  subroutine hecmw_allocate_vector_I(vector, size)
    use hecmw_util
    implicit none

    integer(kind=kint), pointer:: vector(:)
    integer(kind=kint):: size
    integer(kind=kint):: ierr

    allocate( vector( size ), stat=ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'Allocation error: Not enough memory for integer array!'
      stop
    end if

    vector = 0
  end subroutine hecmw_allocate_vector_I

  subroutine hecmw_allocate_vector_R(vector, size)
    use hecmw_util
    implicit none

    real(kind=kreal), pointer:: vector(:)
    integer(kind=kint):: size
    integer(kind=kint):: ierr

    allocate( vector( size ), stat=ierr )
    if ( ierr /= 0 ) then
      write(0,*) 'Allocation error: Not enough memory for real array!'
      stop
    end if

    vector = 0.0d0
  end subroutine hecmw_allocate_vector_R

end module hecmw_allocate
