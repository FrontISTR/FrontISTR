!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see License.txt
!-------------------------------------------------------------------------------
! SA-AMG unit test : F3b (VTK aggregate visualization)
!
! Aggregates a small 2D grid, writes the legacy-VTK point cloud, then reads the
! agg_id scalars back and checks they match the aggregation (point count + ids).
!-------------------------------------------------------------------------------
program test_step8
  use hecmw_precond_SAAMG_util
  use hecmw_precond_SAAMG_matrix
  use hecmw_precond_SAAMG_aggregate
  implicit none

  integer(kind=kint), parameter :: nx = 8, ny = 8
  type(hecmwST_saamg_bcsr) :: A
  type(hecmwST_saamg_nodegraph) :: g
  integer(kind=kint), allocatable :: aggr(:), aggr_read(:)
  real(kind=kreal), allocatable :: coord(:), dense(:,:)
  integer(kind=kint) :: n, naggr, ix, iy, id, npts_read, npass, nfail, iu, i
  character(len=64) :: line
  character(len=*), parameter :: fname = 'agg_test.vtk'

  npass = 0; nfail = 0
  n = nx*ny

  ! 2D grid operator + coordinates (z=0)
  allocate(dense(n,n), coord(3*n)); dense = 0.0d0
  do iy = 1, ny
    do ix = 1, nx
      id = (iy-1)*nx + ix
      dense(id,id) = 4.0d0
      if (ix > 1)  dense(id, id-1)  = -1.0d0
      if (ix < nx) dense(id, id+1)  = -1.0d0
      if (iy > 1)  dense(id, id-nx) = -1.0d0
      if (iy < ny) dense(id, id+nx) = -1.0d0
      coord(3*id-2) = real(ix,kreal)
      coord(3*id-1) = real(iy,kreal)
      coord(3*id  ) = 0.0d0
    end do
  end do
  call hecmw_saamg_bcsr_from_dense(dense, n, 1, A)
  call hecmw_saamg_build_nodegraph(A, g)
  call hecmw_saamg_aggregate(g, 3, 9, aggr, naggr)

  call hecmw_saamg_write_vtk(fname, coord, n, aggr)

  ! read back: find POINTS count and the agg_id scalar block
  allocate(aggr_read(n)); aggr_read = -999
  npts_read = -1
  open(newunit=iu, file=fname, status='old', action='read')
  do
    read(iu,'(a)',end=100) line
    if (line(1:7) == 'POINTS ') read(line(8:),*) npts_read
    if (line(1:15) == 'LOOKUP_TABLE de') then
      do i = 1, n
        read(iu,*) aggr_read(i)
      end do
      exit
    end if
  end do
100 close(iu)

  call check('VTK point count matches', npts_read == n, npass, nfail)
  call check('VTK agg_id values match aggregation', &
       all(aggr_read(1:n) == aggr(1:n)), npass, nfail)
  call check('agg ids in 1..naggr (no excluded here)', &
       minval(aggr) >= 1 .and. maxval(aggr) <= naggr, npass, nfail)
  write(*,'(a,i0,a,i0,a)') '   wrote '//fname//' : ', n, ' points, ', naggr, ' aggregates'

  deallocate(dense, coord, aggr, aggr_read)
  call hecmw_saamg_nodegraph_free(g); call hecmw_saamg_bcsr_free(A)

  write(*,'(/,a,i0,a,i0,a)') '=== F3b: ', npass, ' passed, ', nfail, ' failed ==='
  if (nfail > 0) stop 1

contains

  subroutine check(name, cond, npass, nfail)
    character(len=*),   intent(in)    :: name
    logical,            intent(in)    :: cond
    integer(kind=kint), intent(inout) :: npass, nfail
    if (cond) then
      npass = npass + 1; write(*,'(a,a)') '  [PASS] ', name
    else
      nfail = nfail + 1; write(*,'(a,a)') '  [FAIL] ', name
    end if
  end subroutine check

end program test_step8
