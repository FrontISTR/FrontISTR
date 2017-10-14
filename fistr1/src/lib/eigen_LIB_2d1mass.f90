!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  This module contains subroutines used in 2d eigen analysis for
!!  element 231, 241
module m_eigen_LIB_2d1mass
contains

  !C*--------------------------------------------------------------------*
  subroutine MASS_C2D4(XX,YY,EE,PP,RHO,PARAM1,SS,ISET,fstrEIG)
    !C*--------------------------------------------------------------------*
    !C*
    !C* CALCULATION 2D 4 NODE PLANE ELEMENT
    !C*
    use hecmw
    use m_fstr
    use gauss_integration

    implicit none
    !C* I/F VARIABLES
    real(kind=kreal) XX(*),YY(*),SS(*),EE,PP,PARAM1
    integer(kind=kint) ISET
    !C* LOCAL VARIABLES
    integer(kind=kint) NN,NDOF,ISIZE
    integer(kind=kint) NG
    parameter(NN=4,NDOF=2,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=3)
    real(kind=kreal) D(4,4),B(4,NDOF*NN),DB(4,NDOF*NN)
    real(kind=kreal) H(NN),HR(NN),HS(NN)
    real(kind=kreal) RHO,THICK,COEF1,COEF2,PAI
    real(kind=kreal) RI,SI,RP,SP,RM,SM
    real(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,RR,WG,DUM
    real(kind=kreal) XJI11,XJI21,XJI12,XJI22
    integer(kind=kint) I,J,J2,K,LX,LY,NUM
    integer(kind=kint) ind1, ind2
    real(kind=kreal) totdiag, totmass
    type(fstr_eigen) :: fstrEIG

    totdiag = 0.0
    totmass = 0.0
    !C**************************
    !C*  CONSTANT
    !C**************************
    PAI=4.0*atan(1.0)
    do I=1,ISIZE
      SS(I)=0.0
    enddo
    !C*THICKNESS
    THICK=PARAM1
    !C*FOR AX-SYM. ANALYSIS
    if(ISET.EQ.2) then
      THICK=1.0
    end if
    !C*CLEAR [D] MATRIX
    do J=1,4
      do I=1,4
        D(I,J)=0.0
      enddo
    enddo
    !C*LOOP OVER ALL INTEGRATION POINTS
    do LX=1,NG
      RI=XG(NG,LX)
      do LY=1,NG
        SI=XG(NG,LY)
        RP=1.0+RI
        SP=1.0+SI
        RM=1.0-RI
        SM=1.0-SI
        H(1)=0.25*RM*SM
        H(2)=0.25*RP*SM
        H(3)=0.25*RP*SP
        H(4)=0.25*RM*SP
        HR(1)=-.25*SM
        HR(2)= .25*SM
        HR(3)= .25*SP
        HR(4)=-.25*SP
        HS(1)=-.25*RM
        HS(2)=-.25*RP
        HS(3)= .25*RP
        HS(4)= .25*RM
        XJ11=0.0
        XJ21=0.0
        XJ12=0.0
        XJ22=0.0
        do I=1,NN
          XJ11=XJ11+HR(I)*XX(I)
          XJ21=XJ21+HS(I)*XX(I)
          XJ12=XJ12+HR(I)*YY(I)
          XJ22=XJ22+HS(I)*YY(I)
        enddo
        DET=XJ11*XJ22-XJ21*XJ12
        if(ISET.EQ.2) then
          RR=0.0
          do I=1,NN
            RR=RR+H(I)*XX(I)
          enddo
          WG=WGT(NG,LX)*WGT(NG,LY)*DET*RR*2.0*PAI
        else
          RR=THICK
          WG=WGT(NG,LX)*WGT(NG,LY)*DET*RR
        end if
        DUM=1.0/DET
        XJI11= XJ22*DUM
        XJI12=-XJ12*DUM
        XJI21=-XJ21*DUM
        XJI22= XJ11*DUM

        ind1 = 1
        do J=1,NN*NDOF
          ind2 = 1
          do J2=1,J
            NUM=(J-1)*J/2+J2
            if(J.EQ.J2) then
              SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
              totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
            endif
            if(mod(J2,ndof).EQ.mod(J,ndof)) then
              totmass=totmass+H(ind1)*H(ind2)*WG*RHO
            endif
            if(mod(J2,NDOF).EQ.0) ind2 = ind2 + 1
          enddo
          if(mod(J,NDOF).EQ.0) ind1 = ind1 + 1
        enddo
        !C
      enddo
    enddo
    !C
    do j = 1,nn*ndof
      num = j*(j+1)/2
      ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
    end do
    if(2*totmass-totdiag .LT. 0) stop "Mass error 241!"

    return
  end subroutine MASS_C2D4
  !C*--------------------------------------------------------------------*
  subroutine MASS_C2D3(XX,YY,EE,PP,RHO,PARAM1,SS,ISET,fstrEIG)
    !C*--------------------------------------------------------------------*
    !C*
    !C* CALCULATION 2D 3 NODE PLANE ELEMENT
    !C*
    use hecmw
    use m_fstr
    use gauss_integration
    implicit none
    !C*I/F VARIABLES
    real(kind=kreal) XX(*),YY(*),SS(*),EE,PP,PARAM1
    integer(kind=kint) ISET
    !C*LOCAL VARIABLES
    integer(kind=kint) NN,NDOF,ISIZE
    integer(kind=kint) NG
    parameter(NN=3,NDOF=2,ISIZE=(NN*NDOF)*(NN*NDOF+1)/2,NG=2)
    real(kind=kreal) D(4,4),B(4,NDOF*NN),DB(4,NDOF*NN)
    real(kind=kreal) H(NN),HL1(NN),HL2(NN),HL3(NN)
    integer(kind=kint) L1,L2
    real(kind=kreal) X1,X2,X3,XL1,XL2
    real(kind=kreal) DNDX,DNDY
    real(kind=kreal) RHO,THICK,COEF1,COEF2,PAI
    real(kind=kreal) RI,SI,RP,SP,RM,SM
    real(kind=kreal) XJ11,XJ21,XJ12,XJ22,DET,RR,WG,DUM
    real(kind=kreal) XJI11,XJI21,XJI12,XJI22
    integer(kind=kint) I,J,J2,K,LX,LY,NUM
    integer(kind=kint) ind1,ind2
    real(kind=kreal) totdiag, totmass
    type(fstr_eigen) :: fstrEIG
    totdiag = 0.0
    totmass = 0.0
    !C**************************
    !C*  CONSTANT
    !C**************************
    PAI=4.0*atan(1.0)
    do I=1,ISIZE
      SS(I)=0.0
    enddo
    !C*THICKNESS
    THICK=PARAM1
    !C*FOR AX-SYM. ANALYSIS
    if(ISET.EQ.2) then
      THICK=1.0
    end if
    !C*CLEAR [D] MATRIX
    do J=1,4
      do I=1,4
        D(I,J)=0.0
      enddo
    enddo
    !C*LOOP OVER ALL INTEGRATION POINTS
    do L2=1,NG
      XL2=XG(NG,L2)
      X2 =(XL2+1.0)*0.5
      do L1=1,NG
        XL1=XG(NG,L1)
        X1=0.5*(1.0-X2)*(XL1+1.0)
        !C*INTERPOLATION FUNCTION
        X3=1.0-X1-X2
        H(1)=X1
        H(2)=X2
        H(3)=X3
        !C*DERIVATIVE OF INTERPOLATION FUNCTION
        !C*FOR L1-COORDINATE
        HL1(1)=1.0
        HL1(2)=0.0
        HL1(3)=0.0
        !C*FOR L2-COORDINATE
        HL2(1)=0.0
        HL2(2)=1.0
        HL2(3)=0.0
        !C*FOR L3-COORDINATE
        HL3(1)=0.0
        HL3(2)=0.0
        HL3(3)=1.0
        !C*JACOBI MATRIX
        XJ11=0.0
        XJ21=0.0
        XJ12=0.0
        XJ22=0.0
        do I=1,NN
          XJ11=XJ11+(HL1(I)-HL3(I))*XX(I)
          XJ21=XJ21+(HL2(I)-HL3(I))*XX(I)
          XJ12=XJ12+(HL1(I)-HL3(I))*YY(I)
          XJ22=XJ22+(HL2(I)-HL3(I))*YY(I)
        enddo
        DET=XJ11*XJ22-XJ21*XJ12
        if(ISET.EQ.2) then
          RR=0.0
          do I=1,NN
            RR=RR+H(I)*XX(I)
          enddo
          WG=WGT(NG,L1)*WGT(NG,L2)*DET*RR*2.0*PAI*(1.0-X2)*0.25
        else
          RR=THICK
          WG=WGT(NG,L1)*WGT(NG,L2)*DET*RR*(1.0-X2)*0.25
        end if
        !C
        ind1 = 1
        do J=1,NN*NDOF
          ind2 = 1
          do J2=1,J
            NUM=(J-1)*J/2+J2
            !C
            if(J.EQ.J2) then
              SS(NUM)=SS(NUM)+H(ind1)*H(ind2)*WG*RHO
              totdiag=totdiag+H(ind1)*H(ind2)*WG*RHO
            endif
            if(mod(J2,ndof).EQ.mod(J,ndof)) then
              totmass=totmass+H(ind1)*H(ind2)*WG*RHO
            endif
            if(mod(J2,NDOF).EQ.0) ind2 = ind2 + 1
          enddo
          if(mod(J,NDOF).EQ.0) ind1 = ind1 + 1
        enddo
        !C
      enddo
    enddo
    !C
    do j = 1,nn*ndof
      num = j*(j+1)/2
      ss(num) = ss(num)*(2*totmass-totdiag)/(totdiag)
    end do
    if(2*totmass-totdiag .LT. 0) stop "Mass error 231!"
    !C
    return
  end subroutine MASS_C2D3

end module m_eigen_LIB_2d1mass
