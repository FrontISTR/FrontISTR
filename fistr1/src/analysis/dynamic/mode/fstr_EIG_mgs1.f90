!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> Modified Gram-Schmidt orthonormalization
module m_fstr_EIG_mgs1
contains

      subroutine mgs1(b,a,amass,m,irank,hecMESH,Gntotal)
!c
!c     modified Gram-Schimdt
!c
      use m_fstr
      use hecmw_util
      use m_fstr_EIG_lanczos_util
      implicit none
      REAL(kind=kreal), DIMENSION(:), POINTER :: q,p,qw
      REAL(kind=kreal) :: a(m),b(m),amass(m)
      REAL(kind=kreal) :: pp, t, y, z
      INTEGER(kind=kint) :: irank, Gntotal, ierror, m, i
      TYPE (hecmwST_local_mesh) :: hecMESH

!*Allocate local working arrays
      ALLOCATE(q(m),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, MGS1"
      ALLOCATE(p(m),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, MGS1"
      ALLOCATE(qw(m),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, MGS1"
!c
            t=0.
         do 3 i=1,m
            q(i)=a(i)
            p(i) = b(i)
 3    continue

      do i = 1,m
       qw(i) = 0.0D0
      end do
      CALL MATPRO(qw,amass,b,m,1)
      pp = 0.0D0
      do i = 1,Gntotal
       pp = pp + p(i)*p(i)
      end do
      CALL hecmw_allreduce_R1(hecMESH,pp,hecmw_sum)
         y = 0.0D0
         do i = 1,Gntotal
          y = y+p(i)*qw(i)
         end do
      CALL hecmw_allreduce_R1(hecMESH,y,hecmw_sum)
            z=0.
            do 6 i=1,Gntotal
               z=z+q(i)*qw(i)
 6          continue
      CALL hecmw_allreduce_R1(hecMESH,z,hecmw_sum)
            if(y.ne.0.0D0)then
            do  i=1,m
               q(i)=q(i)-(z/y)*p(i)
            end do
            endif

       do 8 i = 1,m
        a(i) = q(i)
 8    continue

!DEBUG
!*EHM DEBUG 31Mar04: Uncomment for debug
!      WRITE(ILOG,*) 'GSO z,y: ',z,y
!*Deallocations
        if( associated(q) )  DEALLOCATE(q)
        if( associated(p) )  DEALLOCATE(p)
        if( associated(qw) ) DEALLOCATE(qw)
      return
      end subroutine mgs1

end module m_fstr_EIG_mgs1

