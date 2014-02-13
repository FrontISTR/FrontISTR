!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2008/03/13                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takeshi Kitayama (Univ. of Tokyo)              !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
module m_crs_matrix_lag
use my_hecmw_util_lag

! compact row storage matrix

private

public crs_matrix
public irjctocrs
public symbolicirjctocrs
public crs_matrix_getvec

type crs_matrix
  integer(kind=kint) :: nrows  ! number of rows in whole matrix.
  integer(kind=kint) :: ncols ! number of columns in whole matrix.
  integer(kind=kint) :: nttbr ! number of non zero elements in whole matrix.
  integer(kind=kint) :: ndeg  ! degree of freedom of each element
  integer(kind=kint), pointer :: ia(:)    ! size is ia(nrows+1); ia(k)..ia(k+1)-1 elements in ja(:) show columns of k'th row in whole matrix.
  integer(kind=kint), pointer :: ja(:)    ! size is ja(nttbr); store column indices of non zero elements.
  real(kind=kreal),   pointer :: val(:,:) ! size is val(ndeg*ndeg, nttbr). store matrix elements.
end type crs_matrix

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symbolicirjctocrs(ndeg, nttbr, irow, jcol, ncols, nrows, c)

implicit none

integer(kind=kint), intent(in) :: ndeg     ! degree of freedom of each element
integer(kind=kint), intent(in) :: nttbr    ! number of non zero elements in whole matrix.
integer(kind=kint), intent(in) :: irow(:)  ! irjc matrix style row number of element
integer(kind=kint), intent(in) :: jcol(:)  ! irjc matrix style column number of element
integer(kind=kint), intent(in) :: ncols    ! C column size
integer(kind=kint), intent(in) :: nrows    ! C row size

type (crs_matrix), intent(out) :: c

! internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind=kint), allocatable :: istat(:), jrapt(:), jcolno(:)
integer(kind=kint) :: ntt

integer(kind=kint) :: ipass, i,j,k,l,m,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(istat(nrows))
allocate(jrapt(nttbr))
allocate(jcolno(nttbr))

istat=0
jrapt=0
jcolno=0
ntt=0

do l=1,nttbr
  i=irow(l)
  j=jcol(l)
  call reserv(i,j,istat,jrapt,jcolno,ntt)
end do

! allocation and set c%ia, c%ja
c%nrows=nrows
c%ncols=ncols
c%nttbr=ntt
c%ndeg=ndeg

call crsallocation(c)

call stiajac(c,istat,jrapt,jcolno)

deallocate(istat,jrapt,jcolno)

end subroutine symbolicirjctocrs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine irjctocrs(ndeg, nttbr, irow, jcol, val, ncols, nrows, c)

implicit none

integer(kind=kint), intent(in) :: ndeg     ! degree of freedom of each element
integer(kind=kint), intent(in) :: nttbr    ! number of non zero elements in whole matrix.
integer(kind=kint), intent(in) :: irow(:)  ! irjc matrix style row number of element
integer(kind=kint), intent(in) :: jcol(:)  ! irjc matrix style column number of element
real(kind=kreal),   intent(in) :: val(:,:) ! store matrix element
integer(kind=kint), intent(in) :: ncols    ! C column size
integer(kind=kint), intent(in) :: nrows    ! C row size

type (crs_matrix), intent(out) :: c

! internal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(kind=kint), allocatable :: istat(:), jrapt(:), jcolno(:)
integer(kind=kint) :: ntt

integer(kind=kint) :: ipass, i,j,k,l,m,n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(istat(nrows))
allocate(jrapt(nttbr))
allocate(jcolno(nttbr))

istat=0
jrapt=0
jcolno=0
ntt=0

do ipass=1,2
  do l=1,nttbr
    i=irow(l)
    j=jcol(l)

    if(ipass.eq.1)then
      call reserv(i,j,istat,jrapt,jcolno,ntt)
    endif
    if(ipass.eq.2)then
      call stval(c,i,j,val(:,l),0)
    endif
  end do

  ! allocation and set c%ia, c%ja
  if(ipass.eq.1)then
    c%nrows=nrows
    c%ncols=ncols
    c%nttbr=ntt
    c%ndeg=ndeg

    call crsallocation(c)

    call stiajac(c,istat,jrapt,jcolno)

    deallocate(istat,jrapt,jcolno)
  endif
end do

end subroutine irjctocrs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine reserv(i,j,istat,jrapt,jcolno,k)
      ! count number of non zero elements
      implicit none
      integer(kind=kint) :: jrapt(:)
      integer(kind=kint) :: jcolno(:)
      integer(kind=kint) :: istat(:)
      integer(kind=kint) :: i,j,k,l, locr, loc
      locr=0
      loc=istat(i)
  110 continue
      if(loc.eq.0) goto 120
      if(jcolno(loc).eq.j) then
         goto 100
      elseif(jcolno(loc).gt.j) then
         goto 130
      endif
      locr=loc
      loc=jrapt(loc)
      goto 110
  120 continue
      k=k+1
      if(locr.eq.0) then
        istat(i)=k
      else
        jrapt(locr)=k
      endif
      jcolno(k)=j
      goto 150
  130 continue
      k=k+1
      if(locr.eq.0) then
        istat(i)=k
      else
        jrapt(locr)=k
      endif
      jrapt(k)=loc
      jcolno(k)=j
  150 continue
  100 continue
      return
      end subroutine reserv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine crsallocation(c)
      implicit none
      type (crs_matrix) :: c

      if (0 >= c%nrows) then
        call errtrp('set positive nrows')
      else if (0 >= c%ndeg) then
        call errtrp('set positive ndeg')
      else if (0 >= c%nttbr) then
        call errtrp('set positive nttbr')
      end if

      allocate(c%ia(c%nrows+1))
      allocate(c%ja(c%nttbr))
      allocate(c%val(c%ndeg*c%ndeg,c%nttbr))

      return
      end subroutine crsallocation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine stiajac(c,istat,jrapt,jcolno)
      implicit none
! set ia, ja
      type (crs_matrix) :: c
      integer(kind=kint) :: istat(:)
      integer(kind=kint) :: jrapt(:)
      integer(kind=kint) :: jcolno(:)
      integer(kind=kint) :: i,j,k,l, ii, loc, idbg

      if (0 >= c%ncols) then
        call errtrp('set positive ncols')
      else if (0 >= c%nrows) then
        call errtrp('set positive nrows')
      else if (0 >= c%nttbr) then
        call errtrp('set positive nttbr')
      else if (.not. associated(c%ia)) then
        call errtrp('ia is not allocated')
      else if (.not. associated(c%ja)) then
        call errtrp('ja is not allocated')
      else if (c%nrows+1 /= size(c%ia)) then
        call errtrp('ia size unmatched with nrows')
      else if (c%nttbr /= size(c%ja)) then
        call errtrp('ja size unmatched with nttbr')
      else if ((c%ndeg*c%ndeg /= size(c%val, DIM=1)) .or. (c%nttbr /= size(c%val, DIM=2))) then
        call errtrp('ja size unmatched with ndeg, nttbr')
      end if

      idbg=0
      c%ia(1)=1
      l=0
      do 100 k=1,c%nrows
         loc=istat(k)
  110    continue
         if(loc.eq.0) goto 120
         ii=jcolno(loc)
         l=l+1
         c%ja(l)=ii
         loc=jrapt(loc)
         goto 110
  120    c%ia(k+1)=l+1
  100 continue
      if(idbg.ne.0) then
         write(20,*) 'c%ia '
         write(20,60) (c%ia(i),i=1,c%nrows+1) ! ; call flush(20)
         write(20,*) 'c%ja '
         write(20,60) (c%ja(i),i=1,c%ja(c%nrows+1))
      end if
   60 format(10i7)
      return
      end subroutine stiajac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stval(c,i,j,val,itrans)
! store matrix element
implicit none

type (crs_matrix) :: c
real(kind=kreal), dimension(:) :: val
integer(kind=kint) :: ndeg,itrans
integer(kind=kint) :: i,j,k,m,n

ndeg=c%ndeg

do k=c%ia(i),c%ia(i+1)-1
  if(j.eq.c%ja(k)) then
    if(itrans.eq.0) then
      c%val(:,k)=val(:)
    else
      do m=1,ndeg
        do n=1,ndeg
          c%val(m + (n-1)*ndeg,k)=val((m-1)*ndeg + n)
        end do
      end do
    end if
    return
  end if
end do
write(6,*) "something wrong"
stop
end subroutine stval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine crs_matrix_getvec(c, k, v)
! return k'th row vector in whole matrix
implicit none
type (crs_matrix), intent(in) :: c
integer(kind=kint), intent(in) :: k ! given in "ndeg opened" numbering.
real(kind=kreal), dimension(:), intent(out) :: v ! output as "ndeg opened" dens vector

integer(kind=kint) :: ndeg, nndeg, i, idx, jcol, iofset

ndeg=c%ndeg
nndeg=ndeg*ndeg
v=0
jcol = (k+ndeg-1) / ndeg     ! row number in sparse matrix 
iofset = mod(k+ndeg-1, ndeg) ! offset in val. 0offset

do i=c%ia(jcol),c%ia(jcol+1)-1
  idx=c%ja(i)                  ! row number in sparse matrix
  v(ndeg*(idx-1)+1:ndeg*(idx-1)+ndeg)  &
& =c%val(ndeg*iofset + 1 : ndeg*iofset + ndeg ,i)
end do

return
end subroutine crs_matrix_getvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine errtrp(mes) 
character(*) mes
write(*,*)  'Error in m_crs_matrix: ', mes
stop
end subroutine errtrp

end module m_crs_matrix_lag
