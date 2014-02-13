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

module hecmw_Jacob241
contains

   subroutine hecmw_Jacob_241 ( hecMESH, iElem, DET, W, N, NX, NY )
      use hecmw_util

      implicit none

      type(hecmwST_local_mesh):: hecMESH
      integer(kind=kint)::       iElem
      real(kind=kreal)::         DET
      real(kind=kreal)::         W(4), N(4,4), NX(4,4), NY(4,4)

      integer(kind=kint):: i, j, jj, iLocal
      integer(kind=kint):: LX, LY
      real(kind=kreal)::   DUM
      real(kind=kreal)::   XX(4), YY(4)
      real(kind=kreal)::   XG(2), HR(2,2,4), HS(2,2,4), HT(2,2,4)
      real(kind=kreal)::   H(2,2,4),BX(2,2,4),BY(2,2,4),BZ(2,2,4)
      real(kind=kreal)::   RI, SI, RP, SP, RM, SM
      real(kind=kreal)::   XJ11,XJ12,XJ21,XJ22
      real(kind=kreal)::   XJI11,XJI12,XJI21,XJI22

      DATA XG/-0.5773502691896D0, 0.5773502691896D0/

      do j = 1, 4
         jj = 4 - j
         iLocal = hecMESH%elem_node_item  ( 4*iElem -jj )
         XX(j)  = hecMESH%node( iLocal*2 -1 )
         YY(j)  = hecMESH%node( iLocal*2    )
          W(j)  = 1.0D0
      end do
!C
!C*LOOP OVER ALL INTEGRATION POINTS
      do LX=1,2
         RI=XG(LX)
         do LY=1,2
          SI=XG(LY)
          RP=1.0+RI
          SP=1.0+SI
          RM=1.0-RI
          SM=1.0-SI
!!
!! ******** Shape functions ********
!!
          H(LX,LY,1)=0.25*RM*SM
          H(LX,LY,2)=0.25*RP*SM
          H(LX,LY,3)=0.25*RP*SP
          H(LX,LY,4)=0.25*RM*SP
!!
!! ******** Derivative of shape functions ********
!!
      !! ----------- For R-Coordinate -------------
          HR(LX,LY,1)=-.25*SM
          HR(LX,LY,2)= .25*SM
          HR(LX,LY,3)= .25*SP
          HR(LX,LY,4)=-.25*SP
      !! ----------- For S-Coordinate -------------
          HS(LX,LY,1)=-.25*RM
          HS(LX,LY,2)=-.25*RP
          HS(LX,LY,3)= .25*RP
          HS(LX,LY,4)= .25*RM
!!
!! ******** Jacobi matrix calculation********
!!
          XJ11=0.0
          XJ21=0.0
          XJ12=0.0
          XJ22=0.0
          do I=1,4
            XJ11=XJ11+HR(LX,LY,I)*XX(I)
            XJ21=XJ21+HS(LX,LY,I)*XX(I)
            XJ12=XJ12+HR(LX,LY,I)*YY(I)
            XJ22=XJ22+HS(LX,LY,I)*YY(I)
         end do
          DET=XJ11*XJ22-XJ21*XJ12
!!
!! ******** Inverse Jacobi matrix calculation ********
!!
          DUM=1.0/DET
          XJI11= XJ22*DUM
          XJI12=-XJ12*DUM
          XJI21=-XJ21*DUM
          XJI22= XJ11*DUM
          do J=1, 4
            BX(LX,LY,J)=XJI11*HR(LX,LY,J)+XJI12*HS(LX,LY,J)
            BY(LX,LY,J)=XJI21*HR(LX,LY,J)+XJI22*HS(LX,LY,J)
         end do
!C
      end do
   end do
!C

  J = 1
  do LX = 1, 2
     do LY = 1, 2

        N(J,1)  = H(LX,LY,1)
        N(J,2)  = H(LX,LY,2)
        N(J,3)  = H(LX,LY,3)
        N(J,4)  = H(LX,LY,4)

        NX(J,1) = BX(LX,LY,1)
        NX(J,2) = BX(LX,LY,2)
        NX(J,3) = BX(LX,LY,3)
        NX(J,4) = BX(LX,LY,4)

        NY(J,1) = BY(LX,LY,1)
        NY(J,2) = BY(LX,LY,2)
        NY(J,3) = BY(LX,LY,3)
        NY(J,4) = BY(LX,LY,4)

        J = J + 1
     end do
  end do

!C
!C
end subroutine hecmw_Jacob_241

end module hecmw_Jacob241
