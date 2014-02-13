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
!======================================================================!
!                                                                      !
!> \brief  This module contains functions for interpolation in 10 node 
!!    tetrahedron element (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/15
!>  \version    0.00
!======================================================================!

module shape_tet10n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_tet10n(volcoord,shp)
      real(kind=kreal), intent(in) :: volcoord(3)
      real(kind=kreal) :: shp(10)
      real(kind=kreal) :: xi, et, ze, a
      xi= volcoord(1);  et=volcoord(2);  ze=volcoord(3)
      a=1.d0-xi-et-ze
      shp(1)=(2.d0*a-1.d0)*a
      shp(2)=xi*(2.d0*xi-1.d0)
      shp(3)=et*(2.d0*et-1.d0)
      shp(4)=ze*(2.d0*ze-1.d0)
      shp(5)=4.d0*xi*a
      shp(6)=4.d0*xi*et
      shp(7)=4.d0*et*a
      shp(8)=4.d0*ze*a
      shp(9)=4.d0*xi*ze
      shp(10)=4.d0*et*ze
    end subroutine
                     
    subroutine ShapeDeriv_tet10n(volcoord,shp)
      real(kind=kreal), intent(in) :: volcoord(3)
      real(kind=kreal) :: shp(10,3)
      real(kind=kreal) :: xi, et, ze, a
 !     real(kind=kreal) :: x1, x2, x3, x4
      xi= volcoord(1);  et=volcoord(2);  ze=volcoord(3)
      a=1.d0-xi-et-ze
!
!     local derivatives of the shape functions: xi-derivative
!
      shp(1,1)=1.d0-4.d0*a
      shp(2,1)=4.d0*xi-1.d0
      shp(3,1)=0.d0
      shp(4,1)=0.d0
      shp(5,1)=4.d0*(1.d0-2.d0*xi-et-ze)
      shp(6,1)=4.d0*et
      shp(7,1)=-4.d0*et
      shp(8,1)=-4.d0*ze
      shp(9,1)=4.d0*ze
      shp(10,1)=0.d0
!
!     local derivatives of the shape functions: eta-derivative
!
      shp(1,2)=1.d0-4.d0*a
      shp(2,2)=0.d0
      shp(3,2)=4.d0*et-1.d0
      shp(4,2)=0.d0
      shp(5,2)=-4.d0*xi
      shp(6,2)=4.d0*xi
      shp(7,2)=4.d0*(1.d0-xi-2.d0*et-ze)
      shp(8,2)=-4.d0*ze
      shp(9,2)=0.d0
      shp(10,2)=4.d0*ze
!
!     local derivatives of the shape functions: zeta-derivative
!
      shp(1,3)=1.d0-4.d0*a
      shp(2,3)=0.d0
      shp(3,3)=0.d0
      shp(4,3)=4.d0*ze-1.d0
      shp(5,3)=-4.d0*xi
      shp(6,3)=0.d0
      shp(7,3)=-4.d0*et
      shp(8,3)=4.d0*(1.d0-xi-et-2.d0*ze)
      shp(9,3)=4.d0*xi
      shp(10,3)=4.d0*et
    end subroutine

end module
