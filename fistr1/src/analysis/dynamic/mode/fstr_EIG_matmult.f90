!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This module contains subroutine of matrix multply for skyline matrix
module m_fstr_EIG_matmult
contains

      subroutine MATMULT2(hecMESH,hecMAT,BB,ALPA,XX,NUMNP,NDOF)

      use m_fstr
      use hecmw_util
      implicit none
      integer(kind=kint) ::  NUMNP,NDOF,ntotal,i,j,k
      integer(kind=kint) ::  iiS,iiE,kk,kki,kkj,ierror
      real(kind=kreal) :: XX(NUMNP*NDOF)
      REAL(kind=kreal), POINTER, DIMENSION(:) :: yy       !(NUMNP*NDOF)
      real(kind=kreal) :: BB(NUMNP*NDOF),ALPA
      real(kind=kreal) :: X1,X2,WVAL1,WVAL2
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT

!DEBUG
!*For parallel part reduction
!      INTEGER(kind=kint) nglobal,numn,istart,iend, nprocs, myrank, ierror
!      INTEGER(kind=kint) GID
!      INTEGER(kind=kint), POINTER :: istarray(:,:)
!      INTEGER(kind=kint) GROUP,XSAME,hecGROUP

!*Allocate work array
      ALLOCATE(yy(numnp*ndof),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, matmult2"
      ntotal = NUMNP*NDOF
      CALL hecmw_update_2_R(hecMESH,XX,NUMNP)
!C
!C +-------------------------+
!C | {r0}= {b}+alp*[A]{xini} |
!C +-------------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

!C
!C-- BEGIN calculation

      do j= 1, NUMNP
       X1= XX(2*j-1)
       X2= XX(2*j  )
       WVAL1= BB(2*j-1) + ALPA*(hecMAT%D(4*j-3)*X1+hecMAT%D(4*j-2)*X2)
       WVAL2= BB(2*j  ) + ALPA*(hecMAT%D(4*j-1)*X1+hecMAT%D(4*j  )*X2)
       do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
               i= hecMAT%itemL(k)
               X1= XX(2*i-1)
               X2= XX(2*i  )
          WVAL1= WVAL1 + ALPA*(hecMAT%AL(4*k-3)*X1+hecMAT%AL(4*k-2)*X2)
          WVAL2= WVAL2 + ALPA*(hecMAT%AL(4*k-1)*X1+hecMAT%AL(4*k  )*X2)
       enddo

       do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
               i= hecMAT%itemU(k)
               X1= XX(2*i-1)
               X2= XX(2*i  )
         WVAL1= WVAL1 + ALPA*(hecMAT%AU(4*k-3)*X1+hecMAT%AU(4*k-2)*X2)
         WVAL2= WVAL2 + ALPA*(hecMAT%AU(4*k-1)*X1+hecMAT%AU(4*k  )*X2)
       enddo
       YY(2*j-1)= WVAL1
       YY(2*j)  = WVAL2
      enddo

      DO i = 1,NUMNP*NDOF
       XX(i) = YY(i)
      END DO

!*Deallocate work array
      if( associated(yy) ) DEALLOCATE(yy)
      RETURN
      END subroutine MATMULT2

      subroutine MATMULT3(hecMESH,hecMAT,BB,ALPA,XX,NUMNP,NDOF)

      use m_fstr
      use hecmw_util
      implicit none
      integer(kind=kint) ::  NUMNP,NDOF,ntotal,i,j,k
      integer(kind=kint) ::  iiS,iiE,kk,kki,kkj,ierror
      real(kind=kreal) ::  XX(NUMNP*NDOF)
      REAL(kind=kreal), POINTER, DIMENSION(:) :: yy   !(NUMNP*NDOF)
      real(kind=kreal) ::  BB(NUMNP*NDOF),ALPA
      real(kind=kreal) ::  X1,X2,X3,WVAL1,WVAL2,WVAL3
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT

!*Allocate work array
      ALLOCATE(yy(numnp*ndof),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, matmult3"
      ntotal = NUMNP*NDOF
      CALL hecmw_update_3_R(hecMESH,XX,NUMNP)
!C
!C +-------------------------+
!C | {r0}= {b}+alpa*[A]{xini} |
!C +-------------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

!C
!C-- BEGIN calculation

      do j= 1, NUMNP
           X1= XX(3*j-2)
           X2= XX(3*j-1)
           X3= XX(3*j  )
       WVAL1= BB(3*j-2) + ALPA*(hecMAT%D(9*j-8)*X1+ &
     &                       hecMAT%D(9*j-7)*X2+hecMAT%D(9*j-6)*X3)
       WVAL2= BB(3*j-1) + ALPA*(hecMAT%D(9*j-5)*X1+ &
     &                       hecMAT%D(9*j-4)*X2+hecMAT%D(9*j-3)*X3)
       WVAL3= BB(3*j  ) + ALPA*(hecMAT%D(9*j-2)*X1+ &
     &                       hecMAT%D(9*j-1)*X2+hecMAT%D(9*j  )*X3)

        do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
          i= hecMAT%itemL(k)
          X1= XX(3*i-2)
          X2= XX(3*i-1)
          X3= XX(3*i  )
       WVAL1= WVAL1 + ALPA*(hecMAT%AL(9*k-8)*X1+ &
     &                     hecMAT%AL(9*k-7)*X2+hecMAT%AL(9*k-6)*X3)
       WVAL2= WVAL2 + ALPA*(hecMAT%AL(9*k-5)*X1+ &
     &                     hecMAT%AL(9*k-4)*X2+hecMAT%AL(9*k-3)*X3)
       WVAL3= WVAL3 + ALPA*(hecMAT%AL(9*k-2)*X1+ &
     &                     hecMAT%AL(9*k-1)*X2+hecMAT%AL(9*k  )*X3)
        enddo

        do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
          i= hecMAT%itemU(k)
          X1= XX(3*i-2)
          X2= XX(3*i-1)
          X3= XX(3*i  )
       WVAL1= WVAL1 + ALPA*(hecMAT%AU(9*k-8)*X1+ &
     &                     hecMAT%AU(9*k-7)*X2+hecMAT%AU(9*k-6)*X3)
       WVAL2= WVAL2 + ALPA*(hecMAT%AU(9*k-5)*X1+ &
     &                     hecMAT%AU(9*k-4)*X2+hecMAT%AU(9*k-3)*X3)
       WVAL3= WVAL3 + ALPA*(hecMAT%AU(9*k-2)*X1+ &
     &                     hecMAT%AU(9*k-1)*X2+hecMAT%AU(9*k  )*X3)
        enddo

      YY(3*j-2)= WVAL1
      YY(3*j-1)= WVAL2
      YY(3*j  )= WVAL3
      enddo

      DO i = 1,NUMNP*NDOF
       XX(i) = YY(i)
      END DO

!*Deallocate work array
      if( associated(yy) ) DEALLOCATE(yy)
      RETURN
      END subroutine MATMULT3

      subroutine MATMULT6(hecMESH,hecMAT,BB,ALPA,XX,NUMNP,NDOF)

      use m_fstr
      use hecmw_util
      implicit none
      integer(kind=kint) :: NUMNP,NDOF,ntotal,i,j,k
      integer(kind=kint) :: iiS,iiE,kk,kki,kkj,kkl,ierror
      real(kind=kreal) :: XX(NUMNP*NDOF)
      REAL(kind=kreal), POINTER, DIMENSION(:) :: yy  !(NUMNP*NDOF)
      real(kind=kreal) :: BB(NUMNP*NDOF),ALPA
      real(kind=kreal) :: X(NDOF),WVAL(NDOF)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT

!*Allocate work array
      ALLOCATE(yy(numnp*ndof),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, matmult6"
      ntotal = NUMNP*NDOF
      CALL hecmw_update_m_R(hecMESH,XX,NUMNP,NDOF)
!C
!C +-------------------------+
!C | {r0}= {b}+alp*[A]{xini} |
!C +-------------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE

!C
!C-- BEGIN calculation

      do j= 1, NUMNP
       DO i = NDOF-1,0,-1
        X(NDOF-i)= XX(NDOF*j-i)
       END DO

       kki = NDOF*NDOF
       kkj = NDOF*NDOF-1
       DO i = NDOF-1,0,-1
         WVAL(NDOF-i)= BB(NDOF*j-i)
        DO kk = NDOF-1,0,-1
      WVAL(NDOF-i)= WVAL(NDOF-i)+ALPA*hecMAT%D(kki*j-kkj)*X(NDOF-kk)
         kkj = kkj-1
        END DO
       END DO

       do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
               i= hecMAT%itemL(k)
            DO kk = NDOF-1,0,-1
             X(NDOF-kk)= XX(NDOF*i-kk)
            END DO
       kki = NDOF*NDOF
       kkl = NDOF*NDOF-1
       DO kkj = NDOF-1,0,-1
        DO kk = NDOF-1,0,-1
         WVAL(NDOF-kkj)= WVAL(NDOF-kkj)+ &
     &       ALPA*hecMAT%AL(kki*k-kkl)*X(NDOF-kk)
         kkl = kkl - 1
        END DO
       END DO
       enddo

       do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
               i= hecMAT%itemU(k)
          DO kk = NDOF-1,0,-1
           X(NDOF-kk) = xx(NDOF*i-kk)
          END DO
       kki = NDOF*NDOF
       kkl = NDOF*NDOF-1
       DO kkj = NDOF-1,0,-1
        DO kk = NDOF-1,0,-1
         WVAL(NDOF-kkj)= WVAL(NDOF-kkj)+ &
     &       ALPA*hecMAT%AU(kki*k-kkl)*X(NDOF-kk)
         kkl = kkl - 1
        END DO
       END DO
       enddo

      DO i = NDOF-1,0,-1
       YY(NDOF*j-i)= WVAL(NDOF-i)
      END DO
      enddo

      DO i = 1,NUMNP*NDOF
       XX(i) = YY(i)
      END DO

!*Deallocate work array
      if( associated(yy) ) DEALLOCATE(yy)

      RETURN
      END subroutine MATMULT6

end module m_fstr_EIG_matmult

