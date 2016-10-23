!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------

module m_elap
use hecmw_util
! for elaps time
real(kind=kreal)   :: epocht, curt
integer(kind=kint) :: iunit ! output filehandler
logical            :: lout  ! output elaptime or not

private

public initelap
public elapout

contains

subroutine initelap(t,i)
include 'mpif.h'
logical, intent(in) :: t
integer(kind=kint), intent(in) :: i
integer(kind=kint) :: ierr
iunit = i
lout  = t
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
call ptime(epocht)!ELAP
end subroutine initelap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine elapout(mes)
character(*) mes
call ptime(curt)
if (lout) then
  write(iunit,'(a, 1f15.5, 3x, a)') '#elap ',curt - epocht, mes
  ! call flush(iunit)
end if
return
end subroutine elapout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ptime(cputim)
real(kind=kreal) :: cputim
!real(kind=kreal) :: cputim,elaptime
!real x(2)
! machine dependent cpu time by hour
!     cputim=etime(x)
!     cputim=x(1)
cputim=hecmw_Wtime()
return
end subroutine ptime

end module m_elap
