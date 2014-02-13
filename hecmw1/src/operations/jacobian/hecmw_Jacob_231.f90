!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2008/03/13                                         !
!        Category : Jacobian calculation                               !
!                                                                      !
!            Written by Satoshi Ito (Univ. of Tokyo)                   !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_Jacob231
contains

   subroutine hecmw_Jacob_231 ( hecMESH, iElem, DET, W, N, NX, NY )
      use hecmw_util
      implicit none

      type(hecmwST_local_mesh):: hecMESH
      integer(kind=kint)::       iElem
      real(kind=kreal)::         DET
      real(kind=kreal)::         W(3),N(3,3),NX(3,3),NY(3,3)

      integer(kind=kint):: i, ii, j, iLocal
      real(kind=kreal):: DUM
      real(kind=kreal):: XX(3), YY(3)
      real(kind=kreal):: NR(3), NS(3)
      real(kind=kreal):: Jacob(2,2), Jinv(2,2)

!!
!! ******** Set values of coordinates********
!!
      do i = 1, 3
         ii = 3 - i
         iLocal = hecMESH%elem_node_item( 3*iElem -ii )
         XX(i) = hecMESH%node( iLocal*3 -1 )
         YY(i) = hecMESH%node( iLocal*3    )
         W(i)  = 1.0 / 3.0D0
      end do

!!
!! ******** Set values to shape functions ********
!!
      N(1,1) = 0.5D0
      N(2,1) = 0.5D0
      N(3,1) = 0.0D0

      N(1,2) = 0.0D0
      N(2,2) = 0.5D0
      N(3,2) = 0.5D0

      N(1,3) = 0.5D0
      N(2,3) = 0.0D0
      N(3,3) = 0.5D0
!!
!! ******** Derivative of shape functions ********
!!
      !! ----------- For R-Coordinate -------------
      NR(1) = 1.0D0
      NR(2) = 0.0D0
      NR(3) =-1.0D0
      !! ----------- For S-Coordinate -------------
      NS(1) = 0.0D0
      NS(2) = 1.0D0
      NS(3) =-1.0D0
!!
!! ******** Jacobi matrix calculation********
!!
      Jacob(1,1) = NR(1)*XX(1)+NR(2)*XX(2)+NR(3)*XX(3)
      Jacob(2,1) = NS(1)*XX(1)+NS(2)*XX(2)+NS(3)*XX(3)
      Jacob(1,2) = NR(1)*YY(1)+NR(2)*YY(2)+NR(3)*YY(3)
      Jacob(2,2) = NS(1)*YY(1)+NS(2)*YY(2)+NS(3)*YY(3)

      DET = Jacob(1,1)*Jacob(2,2) - Jacob(1,2)*Jacob(2,1)
!!
!! ******** Inverse Jacobi matrix calculation ********
!!
      DUM = 1.0 / DET
      Jinv(1,1) =  DUM * Jacob(2,2)
      Jinv(2,1) = -DUM * Jacob(1,2)
      Jinv(1,2) = -DUM * Jacob(2,1)
      Jinv(2,2) =  DUM * Jacob(1,1)

      do i = 1, 3
         do j = 1, 3

            NX(i,j) = Jinv(1,1)*NR(i) + Jinv(1,2)*NS(i)
            NY(i,j) = Jinv(2,1)*NR(i) + Jinv(2,2)*NS(i)

         end do
      end do

   end subroutine hecmw_Jacob_231
end module hecmw_Jacob231

