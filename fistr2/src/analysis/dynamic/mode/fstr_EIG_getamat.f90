!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Eigen Analysis                                    !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module provides functions of matrix calculation
module m_fstr_EIG_getamat
contains

      subroutine GETAMAT2(AMAT,NUMNP,NDOF)

      use m_fstr
      implicit none
      integer(kind=kint) ::  NUMNP,NDOF,ntotal,i,j,k
      integer(kind=kint) :: iiS,iiE,kk,kki,kkj
      real(kind=kreal) :: AMAT(NUMNP*NDOF,NUMNP*NDOF)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT



      ntotal = NUMNP*NDOF
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
                   AMAT(2*j-1,2*j-1) = hecMAT%D(4*j-3)
                   AMAT(2*j-1,2*j) =   hecMAT%D(4*j-2)
                   AMAT(2*j,2*j-1) =   hecMAT%D(4*j-1)
                   AMAT(2*j,2*j) =     hecMAT%D(4*j)
                do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
                        i= hecMAT%itemL(k)
                   AMAT(2*j-1,2*i-1) = hecMAT%AL(4*k-3)
                   AMAT(2*j-1,2*i) = hecMAT%AL(4*k-2)
                   AMAT(2*j,2*i-1) = hecMAT%AL(4*k-1)
                   AMAT(2*j,2*i) = hecMAT%AL(4*k)
                enddo

                do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
                        i= hecMAT%itemU(k)
                   AMAT(2*j-1,2*i-1) = hecMAT%AU(4*k-3)
                   AMAT(2*j-1,2*i) = hecMAT%AU(4*k-2)
                   AMAT(2*j,2*i-1) = hecMAT%AU(4*k-1)
                   AMAT(2*j,2*i) = hecMAT%AU(4*k)
                enddo
      enddo

      RETURN
      END subroutine GETAMAT2



      subroutine GETAMAT3(AMAT,NUMNP,NDOF)
      use m_fstr
      implicit none
      integer(kind=kint) :: NUMNP,NDOF,ntotal,i,j,k 
      integer(kind=kint) :: iiS,iiE,kk,kki,kkj
      real(kind=kreal) :: AMAT(NUMNP*NDOF,NUMNP*NDOF)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT

       

      ntotal = NUMNP*NDOF
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
                   AMAT(3*j-2,3*j-2) = hecMAT%D(9*j-8)
                   AMAT(3*j-2,3*j-1) = hecMAT%D(9*j-7)
                   AMAT(3*j-2,3*j  ) = hecMAT%D(9*j-6)
                   AMAT(3*j-1,3*j-2) = hecMAT%D(9*j-5)
                   AMAT(3*j-1,3*j-1) = hecMAT%D(9*j-4)
                   AMAT(3*j-1,3*j  ) = hecMAT%D(9*j-3)
                   AMAT(3*j  ,3*j-2) = hecMAT%D(9*j-2)
                   AMAT(3*j  ,3*j-1) = hecMAT%D(9*j-1)
                   AMAT(3*j  ,3*j  ) = hecMAT%D(9*j  )

                do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
                        i= hecMAT%itemL(k)
                   AMAT(3*j-2,3*i-2) = hecMAT%AL(9*k-8)
                   AMAT(3*j-2,3*i-1) = hecMAT%AL(9*k-7)
                   AMAT(3*j-2,3*i  ) = hecMAT%AL(9*k-6)
                   AMAT(3*j-1,3*i-2) = hecMAT%AL(9*k-5)
                   AMAT(3*j-1,3*i-1) = hecMAT%AL(9*k-4)
                   AMAT(3*j-1,3*i  ) = hecMAT%AL(9*k-3)
                   AMAT(3*j  ,3*i-2) = hecMAT%AL(9*k-2)
                   AMAT(3*j  ,3*i-1) = hecMAT%AL(9*k-1)
                   AMAT(3*j  ,3*i  ) = hecMAT%AL(9*k  )
                enddo

                do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
                        i= hecMAT%itemU(k)
                   AMAT(3*j-2,3*i-2) = hecMAT%AU(9*k-8)
                   AMAT(3*j-2,3*i-1) = hecMAT%AU(9*k-7)
                   AMAT(3*j-2,3*i  ) = hecMAT%AU(9*k-6)
                   AMAT(3*j-1,3*i-2) = hecMAT%AU(9*k-5)
                   AMAT(3*j-1,3*i-1) = hecMAT%AU(9*k-4)
                   AMAT(3*j-1,3*i  ) = hecMAT%AU(9*k-3)
                   AMAT(3*j  ,3*i-2) = hecMAT%AU(9*k-2)
                   AMAT(3*j  ,3*i-1) = hecMAT%AU(9*k-1)
                   AMAT(3*j  ,3*i  ) = hecMAT%AU(9*k  )
                enddo
      enddo

      RETURN
      END subroutine GETAMAT3

      subroutine GETAMAT6(AMAT,NUMNP,NDOF)

      use m_fstr
      implicit none
      integer(kind=kint) :: NUMNP,NDOF,ntotal,i,j,k
      integer(kind=kint) :: iiS,iiE,kk,kki,kkj
      real(kind=kreal) :: AMAT(NUMNP*NDOF,NUMNP*NDOF)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT



      ntotal = NUMNP*NDOF
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

                KK = NDOF*NDOF-1
                DO K = NDOF-1,0,-1
                 DO I = NDOF-1,0,-1
                  AMAT(NDOF*j-K,NDOF*j-I) = hecMAT%D(NDOF*j - KK)
                  KK = KK - 1
                 END DO
                END DO

                do k= hecMAT%indexL(j-1)+1, hecMAT%indexL(j)
                        i= hecMAT%itemL(k)
                  KK = NDOF*NDOF-1
                  DO KKj = NDOF-1,0,-1
                   DO KKi = NDOF-1,0,-1
                    AMAT(NDOF*j-KKj,NDOF*i-KKi) = hecMAT%AL(NDOF*k-KK)
                    KK = KK - 1
                   END DO
                  END DO
                enddo

                do k= hecMAT%indexU(j-1)+1, hecMAT%indexU(j)
                        i= hecMAT%itemU(k)
                   KK = NDOF*NDOF - 1
                   DO KKj = NDOF-1,0,-1
                    DO KKi = NDOF-1,0,-1
                     AMAT(NDOF*j-KKj,NDOF*i-KKi) = hecMAT%AU(NDOF*k-KK)
                     KK = KK - 1
                    END DO
                   END DO
                enddo
      enddo

      RETURN
      END subroutine GETAMAT6

end module m_fstr_EIG_getamat

