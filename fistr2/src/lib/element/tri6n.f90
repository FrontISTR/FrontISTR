!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
!> \brief  This module contains functions for interpolation in 6 node 
!!  trianglar element (Langrange  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/05/12
!>  \version    0.00
!======================================================================!

module shape_tri6n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_tri6n(areacoord,func)
      real(kind=kreal), intent(in) :: areacoord(2)
      real(kind=kreal) :: func(6)
      real(kind=kreal) :: xi,et,st
      xi=areacoord(1);  et=areacoord(2); st=1.d0-xi-et

      func(1)=st*(2.d0*st-1.d0)
      func(2)=xi*(2.d0*xi-1.d0)
      func(3)=et*(2.d0*et-1.d0)
      func(4)=4.d0*xi*st
      func(5)=4.d0*xi*et
      func(6)=4.d0*et*st
    end subroutine

    subroutine ShapeDeriv_tri6n(areacoord,func)
      real(kind=kreal), intent(in) :: areacoord(2)
      real(kind=kreal) :: func(6,2)
      real(kind=kreal) :: xi,et,st
      xi=areacoord(1);  et=areacoord(2); st=1.d0-xi-et

      func(1,1)=1.d0-4.d0*st
      func(2,1)=4.d0*xi-1.d0
      func(3,1)=0.d0
      func(4,1)=4.d0*(st-xi)
      func(5,1)=4.d0*et
      func(6,1)=-4.d0*et

      func(1,2)=1.d0-4.d0*st
      func(2,2)=0.d0
      func(3,2)=4.d0*et-1.d0
      func(4,2)=-4.d0*xi
      func(5,2)=4.d0*xi
      func(6,2)=4.d0*(st-et)
      
    end subroutine
    
    subroutine Shape2ndDeriv_tri6n(func)
      real(kind=kreal) :: func(6,2,2)
      func(1,1,1) = 4.d0;  func(1,1,2) = 4.d0
      func(2,1,1) = 4.d0;  func(2,1,2) = 0.d0
      func(3,1,1) = 0.d0;  func(3,1,2) = 0.d0
      func(4,1,1) =-8.d0;  func(4,1,2) = -4.d0
      func(5,1,1) = 0.d0;  func(5,1,2) = 4.d0
      func(6,1,1) = 0.d0;  func(6,1,2) = -4.d0

      func(1,2,1) = 4.d0;  func(1,2,2) = 4.d0
      func(2,2,1) = 0.d0;  func(2,2,2) = 0.d0
      func(3,2,1) = 0.d0;  func(3,2,2) = 4.d0
      func(4,2,1) =-4.d0;  func(4,2,2) = 0.d0
      func(5,2,1) = 4.d0;  func(5,2,2) = 0.d0
      func(6,2,1) =-4.d0;  func(6,2,2) =-8.d0
    end subroutine

end module
