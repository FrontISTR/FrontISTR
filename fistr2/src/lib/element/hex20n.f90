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
!> \brief  This module contains functions for interpolation in 20 node
!!      hexahedral element  (Serendipity  interpolation) 
!                                                                      !
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/04/15
!>  \version    0.00
!======================================================================!


module shape_hex20n
  integer, parameter, private :: kreal = kind(0.0d0)

  contains
    subroutine ShapeFunc_hex20n(localcoord,func)
      real(kind=kreal) :: localcoord(3)
      real(kind=kreal) :: func(20)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      RI=localcoord(1);  SI=localcoord(2);  TI=localcoord(3)
      RP=1.0+RI; SP=1.0+SI; TP=1.0+TI
      RM=1.0-RI; SM=1.0-SI; TM=1.0-TI
      func(1)=-0.125*RM*SM*TM*(2.0+RI+SI+TI)
      func(2)=-0.125*RP*SM*TM*(2.0-RI+SI+TI)
      func(3)=-0.125*RP*SP*TM*(2.0-RI-SI+TI)
      func(4)=-0.125*RM*SP*TM*(2.0+RI-SI+TI)
      func(5)=-0.125*RM*SM*TP*(2.0+RI+SI-TI)
      func(6)=-0.125*RP*SM*TP*(2.0-RI+SI-TI)
      func(7)=-0.125*RP*SP*TP*(2.0-RI-SI-TI)
      func(8)=-0.125*RM*SP*TP*(2.0+RI-SI-TI)
      func(9)=0.25*(1.0-RI**2)*SM*TM
      func(10)=0.25*RP*(1.0-SI**2)*TM
      func(11)=0.25*(1.0-RI**2)*SP*TM
      func(12)=0.25*RM*(1.0-SI**2)*TM
      func(13)=0.25*(1.0-RI**2)*SM*TP
      func(14)=0.25*RP*(1.0-SI**2)*TP
      func(15)=0.25*(1.0-RI**2)*SP*TP
      func(16)=0.25*RM*(1.0-SI**2)*TP
      func(17)=0.25*RM*SM*(1.0-TI**2)
      func(18)=0.25*RP*SM*(1.0-TI**2)
      func(19)=0.25*RP*SP*(1.0-TI**2)
      func(20)=0.25*RM*SP*(1.0-TI**2)
    end subroutine

    subroutine ShapeDeriv_hex20n(localcoord,func)
      real(kind=kreal) :: localcoord(3)
      real(kind=kreal) :: func(20,3)
      REAL(kind=kreal) RI,SI,TI,RP,SP,TP,RM,SM,TM
      RI=localcoord(1);  SI=localcoord(2);  TI=localcoord(3)
      RP=1.d0+RI; SP=1.d0+SI; TP=1.d0+TI
      RM=1.d0-RI; SM=1.d0-SI; TM=1.d0-TI
!  FOR R-COORDINATE
     func(1,1)=-0.125*RM*SM*TM+0.125*SM*TM*(2.0+RI+SI+TI)
     func(2,1)=+0.125*RP*SM*TM-0.125*SM*TM*(2.0-RI+SI+TI)
     func(3,1)=+0.125*RP*SP*TM-0.125*SP*TM*(2.0-RI-SI+TI)
     func(4,1)=-0.125*RM*SP*TM+0.125*SP*TM*(2.0+RI-SI+TI)
     func(5,1)=-0.125*RM*SM*TP+0.125*SM*TP*(2.0+RI+SI-TI)
     func(6,1)=+0.125*RP*SM*TP-0.125*SM*TP*(2.0-RI+SI-TI)
     func(7,1)=+0.125*RP*SP*TP-0.125*SP*TP*(2.0-RI-SI-TI)
     func(8,1)=-0.125*RM*SP*TP+0.125*SP*TP*(2.0+RI-SI-TI)
     func(9,1 )=-0.50*RI*SM*TM
     func(10,1)=+0.25*(1.0-SI**2)*TM
     func(11,1)=-0.50*RI*SP*TM
     func(12,1)=-0.25*(1.0-SI**2)*TM
     func(13,1)=-0.50*RI*SM*TP
     func(14,1)=+0.25*(1.0-SI**2)*TP
     func(15,1)=-0.50*RI*SP*TP
     func(16,1)=-0.25*(1.0-SI**2)*TP
     func(17,1)=-0.25*SM*(1.0-TI**2)
     func(18,1)=+0.25*SM*(1.0-TI**2)
     func(19,1)=+0.25*SP*(1.0-TI**2)
     func(20,1)=-0.25*SP*(1.0-TI**2)
!  FOR S-COORDINATE
     func(1,2)=-0.125*RM*SM*TM+0.125*RM*TM*(2.0+RI+SI+TI)
     func(2,2)=-0.125*RP*SM*TM+0.125*RP*TM*(2.0-RI+SI+TI)
     func(3,2)=+0.125*RP*SP*TM-0.125*RP*TM*(2.0-RI-SI+TI)
     func(4,2)=+0.125*RM*SP*TM-0.125*RM*TM*(2.0+RI-SI+TI)
     func(5,2)=-0.125*RM*SM*TP+0.125*RM*TP*(2.0+RI+SI-TI)
     func(6,2)=-0.125*RP*SM*TP+0.125*RP*TP*(2.0-RI+SI-TI)
     func(7,2)=+0.125*RP*SP*TP-0.125*RP*TP*(2.0-RI-SI-TI)
     func(8,2)=+0.125*RM*SP*TP-0.125*RM*TP*(2.0+RI-SI-TI)
     func(9,2)=-0.25*(1.0-RI**2)*TM
     func(10,2)=-0.50*RP*SI*TM
     func(11,2)=+0.25*(1.0-RI**2)*TM
     func(12,2)=-0.50*RM*SI*TM
     func(13,2)=-0.25*(1.0-RI**2)*TP
     func(14,2)=-0.50*RP*SI*TP
     func(15,2)=+0.25*(1.0-RI**2)*TP
     func(16,2)=-0.50*RM*SI*TP
     func(17,2)=-0.25*RM*(1.0-TI**2)
     func(18,2)=-0.25*RP*(1.0-TI**2)
     func(19,2)=+0.25*RP*(1.0-TI**2)
     func(20,2)=+0.25*RM*(1.0-TI**2)
!  FOR T-COORDINATE
     func(1,3)=-0.125*RM*SM*TM+0.125*RM*SM*(2.0+RI+SI+TI)
     func(2,3)=-0.125*RP*SM*TM+0.125*RP*SM*(2.0-RI+SI+TI)
     func(3,3)=-0.125*RP*SP*TM+0.125*RP*SP*(2.0-RI-SI+TI)
     func(4,3)=-0.125*RM*SP*TM+0.125*RM*SP*(2.0+RI-SI+TI)
     func(5,3)=+0.125*RM*SM*TP-0.125*RM*SM*(2.0+RI+SI-TI)
     func(6,3)=+0.125*RP*SM*TP-0.125*RP*SM*(2.0-RI+SI-TI)
     func(7,3)=+0.125*RP*SP*TP-0.125*RP*SP*(2.0-RI-SI-TI)
     func(8,3)=+0.125*RM*SP*TP-0.125*RM*SP*(2.0+RI-SI-TI)
     func(9,3)=-0.25*(1.0-RI**2)*SM
     func(10,3)=-0.25*RP*(1.0-SI**2)
     func(11,3)=-0.25*(1.0-RI**2)*SP
     func(12,3)=-0.25*RM*(1.0-SI**2)
     func(13,3)=0.25*(1.0-RI**2)*SM
     func(14,3)=0.25*RP*(1.0-SI**2)
     func(15,3)=0.25*(1.0-RI**2)*SP
     func(16,3)=0.25*RM*(1.0-SI**2)
     func(17,3)=-0.5*RM*SM*TI
     func(18,3)=-0.5*RP*SM*TI
     func(19,3)=-0.5*RP*SP*TI
     func(20,3)=-0.5*RM*SP*TI
    end subroutine

end module
