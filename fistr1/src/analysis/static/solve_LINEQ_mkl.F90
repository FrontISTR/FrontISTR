!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to solve sparse system of
!> \linear equitions using intel MKL direct sparse solver
module m_solve_LINEQ_mkl

  use m_fstr
  use m_set_arrays_directsolver_contact

  implicit none
  !< ------------------------------------------------ input parameters for MKL solver(see Intel(R) MKL Reference Manual)
  integer(kind=8)               :: pt(64)
  integer (kind=kint)           :: maxfct
  integer (kind=kint)           :: mnum
  integer (kind=kint)           :: mtype
  integer (kind=kint)           :: phase
  integer (kind=kint)           :: nrhs
  integer (kind=kint)           :: msglvl
  integer (kind=kint)           :: iparm(64)
  real     (kind=kreal)          :: dparm(64)
  !< ----------------------------------------------- input parameters for MKL solver

  integer (kind=kint)           :: ntdf_previous

contains


  subroutine solve_LINEQ_mkl_init(hecMAT,fstrMAT,is_sym)

    type (hecmwST_matrix)                    :: hecMAT       !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT      !< type fstrST_matrix_contact_lagrange
    logical                                  :: is_sym       !< symmetry of matrix

    call set_pointersANDindices_directsolver(hecMAT,fstrMAT,is_sym)
    phase = -1

  end subroutine solve_LINEQ_mkl_init

  !> \brief This subroutine executes phase 11 of intel MKL solver
  !> \(see Intel(R) MKL Reference Manual)
  subroutine initialize_solver_mkl(hecMAT,fstrMAT)

    type (hecmwST_matrix)                    :: hecMAT       !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT      !< type fstrST_matrix_contact_lagrange

    integer(kind=kint)                       :: ntdf         !< total degree of freedom
    integer(kind=kint)                       :: idum(1), ierr
    real(kind=kreal)                          :: ddum(1)

    ntdf = hecMAT%NP*hecMAT%NDOF + fstrMAT%num_lagrange
    ntdf_previous = ntdf

    maxfct = 1; mnum = 1; nrhs = 1; msglvl = 0; ierr = 0
    iparm = 0
    iparm(1)  = 1; iparm(2)  = 2; iparm(3)  = 1; iparm(10) = 13
    iparm(11) = 1; iparm(13) = 1; iparm(18) =-1; iparm(19) =-1
    !        iparm(21) = 1

    if(phase==-1)then
      call pardiso(pt, maxfct, mnum, mtype, phase, ntdf_previous, ddum, idum, idum, &
        idum, nrhs, iparm, msglvl, ddum, ddum, ierr)
    endif

    if( symmetricMatrixStruc )then
      mtype = -2
    else
      mtype = 11
    endif

    phase = 11
    call pardiso(pt, maxfct, mnum, mtype, phase, ntdf, values, pointers, indices, &
      idum, nrhs, iparm, msglvl, ddum, ddum, ierr)
    if(ierr /= 0) then
      write(*,'(" initialize_solver_mkl: ERROR was detected in phase", 2I2)') phase,ierr
      stop
    endif
    write(*,*) ' [Pardiso_MKL]: Initialization completed ! '

  end subroutine initialize_solver_mkl


  !> \brief This subroutine executes phase 22 and phase 33 of the MKL solver
  !> \(see Intel(R) MKL Reference Manual)
  subroutine solve_LINEQ_mkl(hecMESH,hecMAT,fstrMAT,ierr)

    type (hecmwST_local_mesh)                :: hecMESH        !< hecmw mesh
    type (hecmwST_matrix)                    :: hecMAT         !< type hecmwST_matrix
    type (fstrST_matrix_contact_lagrange)    :: fstrMAT        !< type fstrST_matrix_contact_lagrange
    integer(kind=kint), intent(out)          :: ierr
    integer(kind=kint)                       :: ntdf           !< total degree of freedom
    integer(kind=kint)                       :: idum(1)
    real(kind=kreal)                          :: ddum(1)
    real(kind=kreal), allocatable            :: x(:)           !< solution vector

    call hecmw_mat_dump(hecMAT, hecMESH)

    call set_values_directsolver(hecMAT,fstrMAT)
    if( phase == -1 ) call initialize_solver_mkl(hecMAT,fstrMAT)

    ntdf = hecMAT%NP*hecMAT%NDOF + fstrMAT%num_lagrange

    phase = 22
    call pardiso(pt, maxfct, mnum, mtype, phase, ntdf, values, pointers, indices,idum,   &
      nrhs, iparm, msglvl, ddum, ddum, ierr)
    if(ierr /= 0) then
      write(*,'(" solve_LINEQ_mkl: [Error] was detected in phase ", 2I2)')phase,ierr
      return
    endif
    write(*,*) ' [Pardiso_MKL]: Factorization completed ! '

    allocate(x(size(hecMAT%X)))
    x = 0.0d0

    iparm(8) = 6
    phase = 33
    call pardiso(pt, maxfct, mnum, mtype, phase, ntdf, values, pointers, indices, idum,  &
      nrhs, iparm, msglvl, hecMAT%B, X, ierr)
    if(ierr /= 0) then
      write(*,'(" solve_LINEQ_mkl: [Error] was detected in phase ", 2I2)')phase,ierr
      deallocate(x)
      return
    endif
    write(*,*) ' [Pardiso_MKL]: Solve completed ... '

    hecMAT%X = x

    deallocate(x)

    call hecmw_mat_dump_solution(hecMAT)

  end subroutine solve_LINEQ_mkl


#ifndef WITH_MKL

  subroutine pardiso(pt, maxfct, mnum, mtype, phase, ntdf, values, pointers, indices, &
      idum, nrhs, iparm, msglvl, ddum1, ddum2, ierr)

    use m_fstr

    !< ------------------------------------------------ input parameters for MKL solver(see Intel(R) MKL Reference Manual)
    integer (kind=8)              :: pt(64)
    integer (kind=kint)           :: maxfct
    integer (kind=kint)           :: mnum
    integer (kind=kint)           :: mtype
    integer (kind=kint)           :: phase
    integer (kind=kint)           :: ntdf
    integer (kind=kint)           :: pointers(:)        !< ia
    integer (kind=kint)           :: indices(:)         !< ja
    real     (kind=kreal)          :: values(:)          !< a
    integer (kind=kint)           :: idum(:)
    integer (kind=kint)           :: nrhs
    integer (kind=kint)           :: iparm(64)
    integer (kind=kint)           :: msglvl
    real(kind=kreal)               :: ddum1(:), ddum2(:)
    integer (kind=kint)           :: ierr

    write(*,*) "pardiso: ERROR was detected. Please install MKL library and try again."
    stop

  end subroutine

#endif

end module m_solve_LINEQ_mkl
