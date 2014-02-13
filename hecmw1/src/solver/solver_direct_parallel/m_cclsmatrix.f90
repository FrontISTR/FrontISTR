!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2007/11/21                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Takeshi Kitayama  (Univ. of Tokyo)             !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module m_cclsmatrix
! compact column storage matrix

type ccls_matrix
  integer :: neqns ! number of rows in whole matrix.
  integer :: ncol  ! number of columns in whole matrix.
  integer :: nttbr ! number of non zero elements in whole matrix.
  integer :: ndeg  ! degree of freedom of each element
  integer, pointer :: ia(:)  ! size is ia(ncol+1). ia(k)..ia(k+1)-1 elements in ja(:) show rows of k'th column in whole matrix.
  integer, pointer :: ja(:)  ! size is ja(nttbr). store row indices of non zero elements.
  real(8), pointer :: val(:,:) ! size is val(ndeg*ndeg, nttbr). store matrix elements.
end type ccls_matrix

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reserv(i,j,jstat,irpt,irowno,k)
      ! count number of non zero elements
      implicit none
!      integer irpt(*),irowno(*),jstat(*)
      integer, dimension(:)   :: irpt
      integer, dimension(:) :: irowno
      integer, dimension(:)  :: jstat
      integer i,j,k,l, locr, loc
      locr=0
      loc=jstat(j)
  110 continue
      if(loc.eq.0) goto 120
      if(irowno(loc).eq.i) then
         goto 100
      elseif(irowno(loc).gt.i) then
         goto 130
      endif
      locr=loc
      loc=irpt(loc)
      goto 110
  120 continue
      k=k+1
      if(locr.eq.0) then
        jstat(j)=k
      else
        irpt(locr)=k
      endif
      irowno(k)=i
      goto 150
  130 continue
      k=k+1
      if(locr.eq.0) then
        jstat(j)=k
      else
        irpt(locr)=k
      endif
      irpt(k)=loc
      irowno(k)=i
  150 continue
  100 continue
      return
      end subroutine reserv

      subroutine cclsallocation(c)
      implicit none
      type (ccls_matrix) :: c

      if (0 >= c%ncol) then
        call m_cclsmatrix_errtrp('set positive ncol')
      else if (0 >= c%ndeg) then
        call m_cclsmatrix_errtrp('set positive ndeg')
      else if (0 >= c%nttbr) then
        call m_cclsmatrix_errtrp('set positive nttbr')
      end if

      allocate(c%ia(c%ncol+1))
      allocate(c%ja(c%nttbr))
      allocate(c%val(c%ndeg*c%ndeg,c%nttbr))

      return
      end subroutine cclsallocation

      subroutine stiajac(c,jstat,irpt,irowno)
      implicit none
! set ia, ja
      type (ccls_matrix) :: c
!      integer jstat(*),irpt(*),irowno(*)
      integer, dimension(:) :: jstat
      integer, dimension(:) :: irpt
      integer, dimension(:) :: irowno
      integer i,j,k,l, ii, loc, idbg

      if (0 >= c%neqns) then
        call m_cclsmatrix_errtrp('set positive neqns')
      else if (0 >= c%ncol) then
        call m_cclsmatrix_errtrp('set positive ncol')
      else if (0 >= c%nttbr) then
        call m_cclsmatrix_errtrp('set positive nttbr')
      else if (.false. == associated(c%ia)) then
        call m_cclsmatrix_errtrp('ia is not allocated')
      else if (.false. == associated(c%ja)) then
        call m_cclsmatrix_errtrp('ja is not allocated')
      else if (c%ncol+1 /= size(c%ia)) then
        call m_cclsmatrix_errtrp('ia size unmatched with ncol')
      else if (c%nttbr /= size(c%ja)) then
        call m_cclsmatrix_errtrp('ja size unmatched with nttbr')
      else if ((c%ndeg*c%ndeg /= size(c%val, DIM=1)) .or. (c%nttbr /= size(c%val, DIM=2))) then
        call m_cclsmatrix_errtrp('ja size unmatched with ndeg, nttbr')
      end if

      idbg=0
      c%ia(1)=1
      l=0
      do 100 k=1,c%ncol
         loc=jstat(k)
  110    continue
         if(loc.eq.0) goto 120
         ii=irowno(loc)
         l=l+1
         c%ja(l)=ii
         loc=irpt(loc)
         goto 110
  120    c%ia(k+1)=l+1
  100 continue
      if(idbg.ne.0) then
         write(20,*) 'c%ia '
         write(20,60) (c%ia(i),i=1,c%ncol+1) ! ; call flush(20)
         write(20,*) 'c%ja '
         write(20,60) (c%ja(i),i=1,c%ja(c%ncol+1))
      end if
   60 format(10i7)
      return
      end subroutine stiajac

subroutine stval(c,i,j,val,itrans)
! store matrix element
implicit none

type (ccls_matrix) :: c
real(8), dimension(c%ndeg*c%ndeg) :: val
integer ndeg,itrans
integer i,j,k,m,n

ndeg=c%ndeg

do k=c%ia(j),c%ia(j+1)-1
  if(i.eq.c%ja(k)) then
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

subroutine m_cclsmatrix_getvec(c, k, v)
! return k'th column vector in whole matrix
implicit none
type (ccls_matrix), intent(in) :: c
integer, intent(in) :: k ! given in "ndeg opened" numbering.
!real(8), dimension(c%ndeg*c%neqns), intent(out) :: v ! output as "ndeg opened" dens vector
real(8), dimension(:), intent(out) :: v ! output as "ndeg opened" dens vector

integer ndeg, nndeg, i, idx, jcol, iofset

ndeg=c%ndeg
nndeg=ndeg*ndeg
v=0
jcol = (k+ndeg-1) / ndeg     ! column number in sparse matrix 
iofset = mod(k+ndeg-1, ndeg) ! offset in val. 0offset

do i=c%ia(jcol),c%ia(jcol+1)-1
  idx=c%ja(i)                  ! row number in sparse matrix
  v(ndeg*(idx-1)+1:ndeg*(idx-1)+ndeg)  &
& =c%val(ndeg*iofset + 1 : ndeg*iofset + ndeg ,i)
end do

return
end subroutine m_cclsmatrix_getvec

subroutine m_cclsmatrix_errtrp(mes) 
character(*) mes
write(*,*)  'Error in m_cclsmatrix: ', mes
stop
end subroutine m_cclsmatrix_errtrp

end module m_cclsmatrix
