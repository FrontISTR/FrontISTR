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

module hecmw_matrix_dump
  use hecmw_util
  use hecmw_matrix_misc
  use m_hecmw_comm_f

  private

  public :: HECMW_MAT_DUMP_TYPE_NONE
  public :: HECMW_MAT_DUMP_TYPE_MM
  public :: HECMW_MAT_DUMP_TYPE_CSR
  public :: HECMW_MAT_DUMP_TYPE_BSR

  public :: hecmw_mat_dump
  public :: hecmw_mat_dump_rhs
  public :: hecmw_mat_dump_solution

  integer(kind=kint), parameter :: HECMW_MAT_DUMP_TYPE_NONE = 0
  integer(kind=kint), parameter :: HECMW_MAT_DUMP_TYPE_MM   = 1
  integer(kind=kint), parameter :: HECMW_MAT_DUMP_TYPE_CSR  = 2
  integer(kind=kint), parameter :: HECMW_MAT_DUMP_TYPE_BSR  = 3

  integer, save :: NumCall = 0

contains

  subroutine hecmw_mat_dump( hecMAT, hecMESH )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    type(hecmwST_local_mesh) :: hecMESH
    NumCall = NumCall + 1
    select case( hecmw_mat_get_dump(hecMAT) )
    case (HECMW_MAT_DUMP_TYPE_NONE)
      return
    case (HECMW_MAT_DUMP_TYPE_MM)
      call hecmw_mat_dump_mm(hecMAT)
    case (HECMW_MAT_DUMP_TYPE_CSR)
      call hecmw_mat_dump_csr(hecMAT)
    case (HECMW_MAT_DUMP_TYPE_BSR)
      call hecmw_mat_dump_bsr(hecMAT)
    end select
    call hecmw_mat_dump_rhs(hecMAT)
    if (hecmw_mat_get_dump_exit(hecMAT) /= 0) then
      call hecmw_barrier( hecMESH )
      stop "Exiting program after dumping matrix"
    end if
  end subroutine hecmw_mat_dump

  subroutine make_file_name(ext, fname)
    implicit none
    character(*) :: ext
    character(*) :: fname
    write(fname,"('dump_matrix_',I0,'_',I0,A)") &
         NumCall, hecmw_comm_get_rank(), ext
  end subroutine make_file_name

  subroutine hecmw_mat_dump_mm( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer, parameter :: iDump = 201
    character(len=64) :: fname
    integer :: i, j, i0, j0, idof, jdof, ii, jj
    integer :: idxL0, idxL, idxD, idxU0, idxU
    integer :: n, np, ndof, ndof2, nnz
    character(len=64), parameter :: lineFormat = "(I0,' ',I0,' ',e20.12e3)"
    integer :: stat
    !n = hecMAT%N
    n = hecMAT%NP
    np = hecMAT%NP
    ndof = hecMAT%NDOF
    ndof2 = ndof * ndof
    ! make fname
    call make_file_name('.mm', fname)
    ! open file
    open(iDump, file=fname, status='replace', iostat=stat)
    if (stat /= 0) then
      write(*,*) 'WARNING: cannot open file ', fname, ' for matrix dump'
      return
    end if
    ! header
    write(iDump,"(A)") '%%MatrixMarket matrix coordinate real general'
    nnz = ndof2 * (n + hecMAT%indexL(n) + hecMAT%indexU(n))
    write(iDump,"(I0,' ',I0,' ',I0)") n*ndof, np*ndof, nnz
    idxD = 0
    do i = 1, n
      i0 = (i-1)*ndof
      do idof = 1, ndof
        ii = i0 + idof
        ! Lower
        do j = hecMAT%indexL(i-1)+1,hecMAT%indexL(i)
          j0 = (hecMAT%itemL(j)-1)*ndof
          idxL0 = j0*ndof + (idof-1)*ndof
          do jdof = 1, ndof
            jj = j0 + jdof
            idxL = idxL0 + jdof
            write(iDump,lineFormat) ii, jj, hecMAT%AL(idxL)
          end do
        end do
        ! Diagonal
        j0 = i0
        do jdof = 1, ndof
          jj = j0 + jdof
          idxD = idxD + 1
          write(iDump,lineFormat) ii, jj, hecMAT%D(idxD)
        end do
        ! Upper
        do j = hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
          j0 = (hecMAT%itemU(j)-1)*ndof
          idxU0 = j0*ndof + (idof-1)*ndof
          do jdof = 1, ndof
            jj = j0 + jdof
            idxU = idxU0 + jdof
            write(iDump,lineFormat) ii, jj, hecMAT%AU(idxU)
          end do
        end do
      end do
    end do
    ! close file
    close(iDump)
  end subroutine hecmw_mat_dump_mm

  subroutine hecmw_mat_dump_csr( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer, parameter :: iDump = 201
    character(len=64) :: fname
    integer :: i, j, i0, j0, idof, jdof, ii, jj
    integer :: idx, idxD, idxL, idxU, idxL0, idxU0
    integer :: n, np, ndof, ndof2, nnz, nnz1
    character(len=64), parameter :: lineFormat = "(e20.12e3)"
    integer :: stat
    !n = hecMAT%N
    n = hecMAT%NP
    np = hecMAT%NP
    ndof = hecMAT%NDOF
    ndof2 = ndof * ndof
    ! make fname
    call make_file_name('.csr', fname)
    ! open file
    open(iDump, file=fname, status='replace', iostat=stat)
    if (stat /= 0) then
      write(*,*) 'WARNING: cannot open file ', fname, ' for matrix dump'
      return
    end if
    ! header
    write(iDump,"(A)") '%%CSR matrix real general'
    nnz = ndof2 * (n + hecMAT%indexL(n) + hecMAT%indexU(n))
    write(iDump,"(A)") '% nrow ncol nnonzero'
    write(iDump,"(I0,' ',I0,' ',I0)") n*ndof, np*ndof, nnz
    ! index
    write(iDump,"(A)") '% index(0:nrow)'
    idx = 0
    write(iDump, "(I0)") idx
    do i = 1, n
      nnz1 = ndof * ((hecMAT%indexL(i)-hecMAT%indexL(i-1)) + &
           1 + (hecMAT%indexU(i)-hecMAT%indexU(i-1)))
      do idof = 1, ndof
        idx = idx + nnz1
        write(iDump, "(I0)") idx
      end do
    end do
    ! item
    write(iDump,"(A)") '% item(1:nnonzero)'
    do i = 1, n
      i0 = (i-1)*ndof
      do idof = 1, ndof
        ! Lower
        do j = hecMAT%indexL(i-1)+1,hecMAT%indexL(i)
          j0 = (hecMAT%itemL(j)-1)*ndof
          do jdof = 1, ndof
            jj = j0 + jdof
            write(iDump,"(I0)") jj
          end do
        end do
        ! Diagonal
        j0 = i0
        do jdof = 1, ndof
          jj = j0 + jdof
          write(iDump,"(I0)") jj
        end do
        ! Upper
        do j = hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
          j0 = (hecMAT%itemU(j)-1)*ndof
          do jdof = 1, ndof
            jj = j0 + jdof
            write(iDump,"(I0)") jj
          end do
        end do
      end do
    end do
    ! values
    write(iDump,"(A)") '% value(1:nnonzero)'
    idxD = 0
    do i = 1, n
      i0 = (i-1)*ndof
      do idof = 1, ndof
        ii = i0 + idof
        ! Lower
        do j = hecMAT%indexL(i-1)+1,hecMAT%indexL(i)
          j0 = (hecMAT%itemL(j)-1)*ndof
          idxL0 = j0*ndof + (idof-1)*ndof
          do jdof = 1, ndof
            jj = j0 + jdof
            idxL = idxL0+jdof
            write(iDump,lineFormat) hecMAT%AL(idxL)
          end do
        end do
        ! Diagonal
        j0 = i0
        do jdof = 1, ndof
          jj = j0 + jdof
          idxD = idxD + 1
          write(iDump,lineFormat) hecMAT%D(idxD)
        end do
        ! Upper
        do j = hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
          j0 = (hecMAT%itemU(j)-1)*ndof
          idxU0 = j0*ndof + (idof-1)*ndof
          do jdof = 1, ndof
            jj = j0 + jdof
            idxU = idxU0 + jdof
            write(iDump,lineFormat) hecMAT%AU(idxU)
          end do
        end do
      end do
    end do
    ! close file
    close(iDump)
  end subroutine hecmw_mat_dump_csr

  subroutine hecmw_mat_dump_bsr( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer, parameter :: iDump = 201
    character(len=64) :: fname
    integer :: i, j, jj
    integer :: idx, idxL0, idxD0, idxU0
    integer :: n, np, ndof, ndof2, nnz, nnz1
    character(len=64), parameter :: lineFormat = "(e20.12e3)"
    integer :: stat
    !n = hecMAT%N
    n = hecMAT%NP
    np = hecMAT%NP
    ndof = hecMAT%NDOF
    ndof2 = ndof * ndof
    ! make fname
    call make_file_name('.bsr', fname)
    ! open file
    open(iDump, file=fname, status='replace', iostat=stat)
    if (stat /= 0) then
      write(*,*) 'WARNING: cannot open file ', fname, ' for matrix dump'
      return
    end if
    ! header
    write(iDump,"(A)") '%%Block-CSR matrix real general'
    nnz = n + hecMAT%indexL(n) + hecMAT%indexU(n)
    write(iDump,"(A)") '% nrow ncol nnonzero ndof'
    write(iDump,"(I0,' ',I0,' ',I0,' ',I0)") n, np, nnz, ndof
    ! index
    write(iDump,"(A)") '% index(0:nrow)'
    idx = 0
    write(iDump, "(I0)") idx
    do i = 1, n
      nnz1 = (hecMAT%indexL(i)-hecMAT%indexL(i-1)) + &
           1 + (hecMAT%indexU(i)-hecMAT%indexU(i-1))
      idx = idx + nnz1
      write(iDump, "(I0)") idx
    end do
    ! item
    write(iDump,"(A)") '% item(1:nnonzero)'
    do i = 1, n
      ! Lower
      do j = hecMAT%indexL(i-1)+1,hecMAT%indexL(i)
        write(iDump,"(I0)") hecMAT%itemL(j)
      end do
      ! Diagonal
      write(iDump,"(I0)") i
      ! Upper
      do j = hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
        write(iDump,"(I0)") hecMAT%itemU(j)
      end do
    end do
    ! values
    write(iDump,"(A)") '% value(1:nnonzero*ndof*ndof)'
    idxD0 = 0
    do i = 1, n
      ! Lower
      do j = hecMAT%indexL(i-1)+1,hecMAT%indexL(i)
        jj = hecMAT%itemL(j)
        idxL0 = (jj-1)*ndof2
        write(iDump,lineFormat) hecMAT%AL(idxL0+1:idxL0+ndof2)
      end do
      ! Diagonal
      write(iDump,lineFormat) hecMAT%D(idxD0+1:idxD0+ndof2)
      idxD0 = idxD0 + ndof2
      ! Upper
      do j = hecMAT%indexU(i-1)+1,hecMAT%indexU(i)
        jj = hecMAT%itemU(j)
        idxU0 = (jj-1)*ndof2
        write(iDump,lineFormat) hecMAT%AU(idxU0+1:idxU0+ndof2)
      end do
    end do
    ! close file
    close(iDump)
  end subroutine hecmw_mat_dump_bsr

  subroutine hecmw_mat_dump_rhs( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer, parameter :: iDump = 201
    character(len=64) :: fname
    integer :: i
    integer :: n, np, ndof, ndof2
    character(len=64), parameter :: lineFormat = "(e20.12e3)"
    integer :: stat
    if( hecmw_mat_get_dump(hecMAT) == HECMW_MAT_DUMP_TYPE_NONE) return
    !n = hecMAT%N
    n = hecMAT%NP
    np = hecMAT%NP
    ndof = hecMAT%NDOF
    ndof2 = ndof * ndof
    ! make fname
    call make_file_name('.rhs', fname)
    ! open file
    open(iDump, file=fname, status='replace', iostat=stat)
    if (stat /= 0) then
      write(*,*) 'WARNING: cannot open file ', fname, ' for matrix dump'
      return
    end if
    do i = 1, np*ndof
      write(iDump,lineFormat) hecMAT%B(i)
    end do
    ! close file
    close(iDump)
  end subroutine hecmw_mat_dump_rhs

  subroutine hecmw_mat_dump_solution( hecMAT )
    implicit none
    type(hecmwST_matrix) :: hecMAT
    integer, parameter :: iDump = 201
    character(len=64) :: fname
    integer :: i
    integer :: n, np, ndof, ndof2
    character(len=64), parameter :: lineFormat = "(e20.12e3)"
    integer :: stat
    if( hecmw_mat_get_dump(hecMAT) == HECMW_MAT_DUMP_TYPE_NONE) return
    !n = hecMAT%N
    n = hecMAT%NP
    np = hecMAT%NP
    ndof = hecMAT%NDOF
    ndof2 = ndof * ndof
    ! make fname
    call make_file_name('.sol', fname)
    ! open file
    open(iDump, file=fname, status='replace', iostat=stat)
    if (stat /= 0) then
      write(*,*) 'WARNING: cannot open file ', fname, ' for matrix dump'
      return
    end if
    do i = 1, np*ndof
      write(iDump,lineFormat) hecMAT%X(i)
    end do
    ! close file
    close(iDump)
  end subroutine hecmw_mat_dump_solution

end module hecmw_matrix_dump
