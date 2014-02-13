!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by X. YUAN                                !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!!======================================================================!
!                                                                      !
!> \brief  This module contains functions for interpolation in 15 node
!!    prism element  (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/27
!>  \version    0.00
!======================================================================!

module shape_prism15n
  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_prism15n(ncoord,shp)
      real(kind=kreal), intent(in)  :: ncoord(3)
      real(kind=kreal) :: shp(15)
      real(kind=kreal) :: xi,et,a,ze
      xi = ncoord(1);  et = ncoord(2)
      a=1.d0-xi-et
      ze = ncoord(3)
      shp(1)=0.5* a*(1.0-ze)*(2.0*a -2.0-ze)
      shp(2)=0.5*xi*(1.0-ze)*(2.0*xi-2.0-ze)
      shp(3)=0.5*et*(1.0-ze)*(2.0*et-2.0-ze)
      shp(4)=0.5* a*(1.0+ze)*(2.0*a -2.0+ze)
      shp(5)=0.5*xi*(1.0+ze)*(2.0*xi-2.0+ze)
      shp(6)=0.5*et*(1.0+ze)*(2.0*et-2.0+ze)
      shp(7)=2.0*xi*a*(1.0-ze)
      shp(8)=2.0*xi*et*(1.0-ze)
      shp(9)=2.0*et*a*(1.0-ze)
      shp(10)=2.0*xi*a*(1.0+ze)
      shp(11)=2.0*xi*et*(1.0+ze)
      shp(12)=2.0*et*a*(1.0+ze)
      shp(13)= a*(1.0-ze*ze)
      shp(14)=xi*(1.0-ze*ze)
      shp(15)=et*(1.0-ze*ze)
    end subroutine

    subroutine ShapeDeriv_prism15n(ncoord,func)
      real(kind=kreal), intent(in)  :: ncoord(3)
      real(kind=kreal) :: func(15,3)
      real(kind=kreal) :: xi,et,a,ze
      xi = ncoord(1);  et = ncoord(2)
      a=1.d0-xi-et
      ze = ncoord(3)
!
!     local derivatives of the shape functions: xi-derivative
!
      func(1,1)= -0.5*(1.0-ze)*(4.0*a -ze-2.0)
      func(2,1)=  0.5*(1.0-ze)*(4.0*xi-ze-2.0)
      func(3,1)=  0.d0
      func(4,1)= -0.5*(1.0+ze)*(4.0*a +ze-2.0)
      func(5,1)=  0.5*(1.0+ze)*(4.0*xi+ze-2.0)
      func(6,1)=  0.d0
      func(7,1)=  2.0*(1.0-ze)*(1.0-2.0*xi-et)
      func(8,1)=  2.0*et*(1.0-ze)
      func(9,1)= -2.0*et*(1.0-ze)
      func(10,1)= 2.0*(1.0+ze)*(1.0-2.0*xi-et)
      func(11,1)= 2.0*et*(1.0+ze)
      func(12,1)= -2.0*et*(1.0+ze)
      func(13,1)= -(1.0-ze*ze)
      func(14,1)=  (1.0-ze*ze)
      func(15,1)=  0.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      func(1,2)=-0.5*(1.0-ze)*(4.0*a -ze-2.0)
      func(2,2)= 0.d0
      func(3,2)= 0.5*(1.0-ze)*(4.0*et-ze-2.0)
      func(4,2)=-0.5*(1.0+ze)*(4.0*a +ze-2.0)
      func(5,2)= 0.d0
      func(6,2)= 0.5*(1.0+ze)*(4.0*et+ze-2.0)
      func(7,2)=-2.0*xi*(1.0-ze)
      func(8,2)= 2.0*xi*(1.0-ze)
      func(9,2)= 2.0*(1.0-ze)*(1.0-xi-2.0*et)
      func(10,2)=-2.0*xi*(1.0+ze)
      func(11,2)= 2.0*xi*(1.0+ze)
      func(12,2)= 2.0*(1.0+ze)*(1.0-xi-2.0*et)
      func(13,2)=-(1.0-ze*ze)
      func(14,2)= 0.0d0
      func(15,2)= (1.0-ze*ze)
!
!     local derivatives of the shape functions: zeta-derivative
!
      func(1,3)=  a*(xi+et+ze-0.5)
      func(2,3)= xi*(-xi+ze+0.5)
      func(3,3)= et*(-et+ze+0.5)
      func(4,3)=  a*(-xi-et+ze+0.5)
      func(5,3)= xi*(xi+ze-0.5)
      func(6,3)= et*(et+ze-0.5)
      func(7,3)= -2*xi*a
      func(8,3)= -2*xi*et
      func(9,3)= -2*et*a
      func(10,3)= 2*xi*a
      func(11,3)= 2*xi*et
      func(12,3)= 2*et*a
      func(13,3)=-2*a*ze
      func(14,3)=-2*xi*ze
      func(15,3)=-2*et*ze
    end subroutine

end module
