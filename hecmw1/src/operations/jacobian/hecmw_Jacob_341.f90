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

module hecmw_Jacob341
contains

   subroutine hecmw_Jacob_341 ( hecMESH, iElem, DET, W, N, NX, NY, NZ)
      use hecmw_util
      implicit none

      type(hecmwST_local_mesh):: hecMESH
      integer(kind=kint)::       iElem
      real(kind=kreal)::         DET
      real(kind=kreal)::         W(4),N(4,4),NX(4,4),NY(4,4),NZ(4,4)

      integer(kind=kint):: i, ii, j, iLocal
      real(kind=kreal):: DUM
      real(kind=kreal):: XX(4), YY(4), ZZ(4)
      real(kind=kreal):: NR(4), NS(4), NT(4)
      real(kind=kreal):: Jacob(3,3), Jinv(3,3)
      real(kind=kreal),parameter:: alpha = 0.58541020, beta = 0.13819660

!!
!! ******** Set values of coordinates********
!!
      do i = 1, 4
         ii = 4 - i
         iLocal = hecMESH%elem_node_item( 4*iElem -ii )
         XX(i) = hecMESH%node( iLocal*3 -2 )
         YY(i) = hecMESH%node( iLocal*3 -1 )
         ZZ(i) = hecMESH%node( iLocal*3    )
         W(i)  = 0.25D0
      end do

!!
!! ******** Set values to shape functions ********
!!
      N(1,1) = alpha
      N(2,1) = beta
      N(3,1) = beta
      N(4,1) = beta

      N(1,2) = beta
      N(2,2) = alpha
      N(3,2) = beta
      N(4,2) = beta

      N(1,3) = beta
      N(2,3) = beta
      N(3,3) = alpha
      N(4,3) = beta

      N(1,4) = beta
      N(2,4) = beta
      N(3,4) = beta
      N(4,4) = alpha
!!
!! ******** Derivative of shape functions ********
!!
      !! ----------- For R-Coordinate -------------
      NR(1) = 1.0
      NR(2) = 0.0
      NR(3) = 0.0
      NR(4) =-1.0
      !! ----------- For S-Coordinate -------------
      NS(1) = 0.0
      NS(2) = 1.0
      NS(3) = 0.0
      NS(4) =-1.0
      !! ----------- For T-Coordinate -------------
      NT(1) = 0.0
      NT(2) = 0.0
      NT(3) = 1.0
      NT(4) =-1.0
!!
!! ******** Jacobi matrix calculation********
!!
      Jacob(1,1) = NR(1)*XX(1)+NR(2)*XX(2)+NR(3)*XX(3)+NR(4)*XX(4)
      Jacob(2,1) = NS(1)*XX(1)+NS(2)*XX(2)+NS(3)*XX(3)+NS(4)*XX(4)
      Jacob(3,1) = NT(1)*XX(1)+NT(2)*XX(2)+NT(3)*XX(3)+NT(4)*XX(4)
      Jacob(1,2) = NR(1)*YY(1)+NR(2)*YY(2)+NR(3)*YY(3)+NR(4)*YY(4)
      Jacob(2,2) = NS(1)*YY(1)+NS(2)*YY(2)+NS(3)*YY(3)+NS(4)*YY(4)
      Jacob(3,2) = NT(1)*YY(1)+NT(2)*YY(2)+NT(3)*YY(3)+NT(4)*YY(4)
      Jacob(1,3) = NR(1)*ZZ(1)+NR(2)*ZZ(2)+NR(3)*ZZ(3)+NR(4)*ZZ(4)
      Jacob(2,3) = NS(1)*ZZ(1)+NS(2)*ZZ(2)+NS(3)*ZZ(3)+NS(4)*ZZ(4)
      Jacob(3,3) = NT(1)*ZZ(1)+NT(2)*ZZ(2)+NT(3)*ZZ(3)+NT(4)*ZZ(4)

      DET = Jacob(1,1) * Jacob(2,2) * Jacob(3,3)   &
&         + Jacob(1,2) * Jacob(2,3) * Jacob(3,1)   &
&         + Jacob(1,3) * Jacob(2,1) * Jacob(3,2)   &
&         - Jacob(1,3) * Jacob(2,2) * Jacob(3,1)   &
&         - Jacob(1,2) * Jacob(2,1) * Jacob(3,3)   &
&         - Jacob(1,1) * Jacob(2,3) * Jacob(3,2)
!!
!! ******** Inverse Jacobi matrix calculation ********
!!
      DUM = 1.0 / DET
      Jinv(1,1) = DUM*(  Jacob(2,2)*Jacob(3,3) -Jacob(2,3)*Jacob(3,2) )
      Jinv(2,1) = DUM*( -Jacob(2,1)*Jacob(3,3) +Jacob(2,3)*Jacob(3,1) )
      Jinv(3,1) = DUM*(  Jacob(2,1)*Jacob(3,2) -Jacob(2,2)*Jacob(3,1) )
      Jinv(1,2) = DUM*( -Jacob(1,2)*Jacob(3,3) +Jacob(1,3)*Jacob(3,2) )
      Jinv(2,2) = DUM*(  Jacob(1,1)*Jacob(3,3) -Jacob(1,3)*Jacob(3,1) )
      Jinv(3,2) = DUM*( -Jacob(1,1)*Jacob(3,2) +Jacob(1,2)*Jacob(3,1) )
      Jinv(1,3) = DUM*(  Jacob(1,2)*Jacob(2,3) -Jacob(1,3)*Jacob(2,2) )
      Jinv(2,3) = DUM*( -Jacob(1,1)*Jacob(2,3) +Jacob(1,3)*Jacob(2,1) )
      Jinv(3,3) = DUM*(  Jacob(1,1)*Jacob(2,2) -Jacob(1,2)*Jacob(2,1) )

      do i = 1, 4
         do j = 1, 4

            NX(i,j) = Jinv(1,1)*NR(i) + Jinv(1,2)*NS(i) +Jinv(1,3)*NT(i)
            NY(i,j) = Jinv(2,1)*NR(i) + Jinv(2,2)*NS(i) +Jinv(2,3)*NT(i)
            NZ(i,j) = Jinv(3,1)*NR(i) + Jinv(3,2)*NS(i) +Jinv(3,3)*NT(i)

         end do
      end do

   end subroutine hecmw_Jacob_341

end module hecmw_Jacob341
