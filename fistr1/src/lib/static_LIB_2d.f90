!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides common functions of Plane deformation elements
module m_static_LIB_2d
  use hecmw, only : kint, kreal
  use elementInfo
  implicit none

contains
  !***********************************************************************
  !  2D Element:
  !  STF_C2   :Calculate stiff matrix of 2d elements
  !  DL_C2    :Deal with DLOAD conditions
  !  TLOAD_C2 :Deal with thermal expansion force
  !  UpdateST_C2  : Update Strain and stress
  !***********************************************************************
  !----------------------------------------------------------------------*
  subroutine STF_C2( ETYPE,NN,ecoord,gausses,PARAM1,stiff         &
      ,ISET, u)
    !----------------------------------------------------------------------*
    use mMechGauss
    use m_MatMatrix
    integer(kind=kint), intent(in) :: ETYPE
    integer(kind=kint), intent(in) :: NN
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(in)   :: ecoord(2,NN),PARAM1
    real(kind=kreal), intent(out)  :: stiff(:,:)
    integer(kind=kint),intent(in)  :: ISET
    real(kind=kreal), intent(in), optional :: u(:,:)

    ! LOCAL VARIABLES
    integer(kind=kint) :: flag
    integer(kind=kint) NDOF
    parameter(NDOF=2)
    real(kind=kreal) D(4,4),B(4,NDOF*NN),DB(4,NDOF*NN)
    real(kind=kreal) H(NN),stress(4)
    real(kind=kreal) THICK,PAI
    real(kind=kreal) DET,RR,WG
    real(kind=kreal) localcoord(2)
    real(kind=kreal) gderiv(NN,2), cdsys(3,3)
    integer(kind=kint) J,LX
    real(kind=kreal) gdispderiv(2,2)
    real(kind=kreal) B1(4,nn*ndof)
    real(kind=kreal) Smat(4,4)
    real(kind=kreal) BN(4,nn*ndof), SBN(4,nn*ndof)

    stiff =0.d0
    flag = gausses(1)%pMaterial%nlgeom_flag
    ! THICKNESS
    THICK=PARAM1
    !  FOR AX-SYM. ANALYSIS
    if(ISET==2) then
      THICK=1.d0
      PAI=4.d0*atan(1.d0)
    endif
    !* LOOP OVER ALL INTEGRATION POINTS
    do LX=1,NumOfQuadPoints( ETYPE )
      call MatlMatrix( gausses(LX), ISET, D, 1.d0, 1.d0, cdsys, 0.d0 )
      if( .not. present(u) ) flag=INFINITESIMAL    ! enforce to infinitesimal deformation analysis

      if( flag==1 .and. ISET == 2 ) then
        write(*,'(a)') '    PROGRAM STOP : non-TL element for axisymmetric element'
        stop
      endif

      call getQuadPoint( ETYPE, LX, localcoord )
      call getGlobalDeriv( etype, nn, localcoord, ecoord, det, gderiv )
      !
      if(ISET==2) then
        call getShapeFunc( ETYPE, localcoord, H(:) )
        RR=dot_product( H(1:NN), ecoord(1,1:NN) )
        WG=getWeight( ETYPE, LX )*DET*RR*2.d0*PAI
      else
        RR=THICK
        H(:)=0.d0
        WG=getWeight( ETYPE, LX )*DET*RR
      end if
      do J=1,NN
        B(1,2*J-1)=gderiv(J,1)
        B(2,2*J-1)=0.d0
        B(3,2*J-1)=gderiv(J,2)
        B(1,2*J  )=0.d0
        B(2,2*J  )=gderiv(J,2)
        B(3,2*J  )=gderiv(J,1)
        B(4,2*J-1)=H(J)/RR
        B(4,2*J  )=0.d0
      enddo
      !
      ! ----- calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
      if( flag==1 ) then
        !       ----- dudx(i,j) ==> gdispderiv(i,j)
        gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
        B1(1:4,1:nn*ndof) = 0.d0
        do j=1,nn
          B1(1,2*j-1) = gdispderiv(1,1)*gderiv(j,1)
          B1(2,2*j-1) = gdispderiv(1,2)*gderiv(j,2)
          B1(3,2*j-1) = gdispderiv(1,2)*gderiv(j,1)+gdispderiv(1,1)*gderiv(j,2)
          B1(1,2*j  ) = gdispderiv(2,1)*gderiv(j,1)
          B1(2,2*j  ) = gdispderiv(2,2)*gderiv(j,2)
          B1(3,2*j  ) = gdispderiv(2,2)*gderiv(j,1)+gdispderiv(2,1)*gderiv(j,2)
          B1(4,2*j-1) = 0.d0
          B1(4,2*j  ) = 0.d0
        enddo
        do j=1,nn*ndof
          B(:,j) = B(:,j)+B1(:,j)
        enddo
      endif
      !
      DB(1:4,1:nn*ndof) = 0.d0
      DB(1:4,1:nn*ndof) = matmul( D(1:4,1:4), B(1:4,1:nn*ndof) )
      stiff(1:nn*ndof,1:nn*ndof) = stiff(1:nn*ndof,1:nn*ndof) +             &
        matmul( transpose( B(1:4,1:nn*ndof) ), DB(1:4,1:nn*ndof) )*WG
      !
      !    ----- calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag==1 ) then
        stress(1:4)=gausses(LX)%stress(1:4)
        BN(1:4,1:nn*ndof) = 0.d0
        do j=1,nn
          BN(1,2*j-1) = gderiv(j,1)
          BN(2,2*j  ) = gderiv(j,1)
          BN(3,2*j-1) = gderiv(j,2)
          BN(4,2*j  ) = gderiv(j,2)
        enddo
        Smat(:,:) = 0.d0
        do j=1,2
          Smat(j  ,j  ) = stress(1)
          Smat(j  ,j+2) = stress(3)
          Smat(j+2,j  ) = stress(3)
          Smat(j+2,j+2) = stress(2)
        enddo
        SBN(1:4,1:nn*ndof) = matmul( Smat(1:4,1:4), BN(1:4,1:nn*ndof) )
        stiff(1:nn*ndof,1:nn*ndof) = stiff(1:nn*ndof,1:nn*ndof) +           &
          matmul( transpose( BN(1:4,1:nn*ndof) ), SBN(1:4,1:nn*ndof) )*WG
      endif
      !
    enddo

  end subroutine STF_C2


  !----------------------------------------------------------------------*
  subroutine DL_C2(ETYPE,NN,XX,YY,RHO,PARAM1,LTYPE,PARAMS,VECT,NSIZE,ISET)
    !----------------------------------------------------------------------*
    !**
    !**  Deal with DLOAD conditions
    !**
    !   BX   LTYPE=1  :BODY FORCE IN X-DIRECTION
    !   BY   LTYPE=2  :BODY FORCE IN Y-DIRECTION
    !   GRAV LTYPE=4  :GRAVITY FORCE
    !   CENT LTYPE=5  :CENTRIFUGAL LOAD
    !   P1   LTYPE=10 :TRACTION IN NORMAL-DIRECTION FOR FACE-1
    !   P2   LTYPE=20 :TRACTION IN NORMAL-DIRECTION FOR FACE-2
    !   P3   LTYPE=30 :TRACTION IN NORMAL-DIRECTION FOR FACE-3
    !   P4   LTYPE=40 :TRACTION IN NORMAL-DIRECTION FOR FACE-4
    ! I/F VARIABLES
    integer(kind=kint), intent(in) :: ETYPE,NN
    real(kind=kreal), intent(in)   :: XX(NN),YY(NN),PARAMS(0:6)
    real(kind=kreal), intent(out)  :: VECT(NN*2)
    real(kind=kreal) RHO,PARAM1
    integer(kind=kint) LTYPE,NSIZE,ISET
    ! LOCAL VARIABLES
    integer(kind=kint), parameter :: NDOF=2
    real(kind=kreal) H(NN)
    real(kind=kreal) PLX(NN),PLY(NN)
    real(kind=kreal) XJ(2,2),DET,RR,WG
    real(kind=kreal) PAI,val,THICK
    integer(kind=kint) IVOL,ISUF
    real(kind=kreal) AX,AY,RX,RY
    real(kind=kreal) normal(2)
    real(kind=kreal) COEFX,COEFY,XCOD,YCOD,HX,HY,PHX,PHY
    integer(kind=kint) NOD(NN)
    integer(kind=kint) LX,I,SURTYPE,NSUR
    real(kind=kreal) VX,VY,localcoord(2),deriv(NN,2),elecoord(2,NN)

    rx = 0.0d0; ry = 0.0d0
    ay = 0.0d0; ax = 0.0d0

    !  CONSTANT
    PAI=4.d0*atan(1.d0)
    ! SET VALUE
    val=PARAMS(0)
    ! THICKNESS
    THICK=PARAM1
    ! CLEAR VECT
    NSIZE=NN*NDOF
    VECT(1:NSIZE)=0.d0
    !
    ! SELECTION OF LOAD TYPE
    !
    IVOL=0
    ISUF=0
    if( LTYPE.LT.10 ) then
      IVOL=1
      if( LTYPE.EQ.5 ) then
        AX=PARAMS(1)
        AY=PARAMS(2)
        RX=PARAMS(4)
        RY=PARAMS(5)
      endif
    else if( LTYPE.GE.10 ) then
      ISUF=1
      call getSubFace( ETYPE, LTYPE/10, SURTYPE, NOD )
      NSUR = getNumberOfNodes( SURTYPE )
    endif
    !** SURFACE LOAD
    if( ISUF==1 ) then
      do I=1,NSUR
        elecoord(1,i)=XX(NOD(I))
        elecoord(2,i)=YY(NOD(i))
      enddo
      !** INTEGRATION OVER SURFACE
      do LX=1,NumOfQuadPoints( SURTYPE )
        call getQuadPoint( SURTYPE, LX, localcoord(1:1) )
        call getShapeFunc( SURTYPE, localcoord(1:1), H(1:NSUR) )
        normal=EdgeNormal( SURTYPE, NSUR, localcoord(1:1), elecoord(:,1:NSUR) )
        WG = getWeight( SURTYPE, LX )
        if( ISET==2 ) then
          RR=0.d0
          do I=1,NSUR
            RR=RR+H(I)*XX(NOD(I))
          enddo
          WG=WG*RR*2.d0*PAI
        else
          WG=WG*THICK
        endif
        do I=1,NSUR
          VECT(2*NOD(I)-1)=VECT(2*NOD(I)-1)+val*WG*H(I)*normal(1)
          VECT(2*NOD(I)  )=VECT(2*NOD(I)  )+val*WG*H(I)*normal(2)
        enddo
      enddo
    endif
    !** VOLUME LOAD
    if( IVOL==1 ) then
      PLX(:)=0.d0
      PLY(:)=0.d0
      do LX=1,NumOfQuadPoints( ETYPE )
        call getQuadPoint( ETYPE, LX, localcoord(1:2) )
        call getShapeDeriv( ETYPE, localcoord(1:2), deriv )
        call getShapeFunc( ETYPE, localcoord(1:2), H(1:NN) )
        XJ(1,1:2)=matmul( XX(1:NN), deriv(1:NN,1:2) )
        XJ(2,1:2)=matmul( YY(1:NN), deriv(1:NN,1:2) )

        DET=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)

        WG = getWeight( ETYPE, LX )
        if(ISET==2) then
          RR=dot_product( H(1:NN),XX(1:NN) )
          WG=WG*DET*RR*2.d0*PAI
        else
          RR=THICK
          WG=WG*DET*RR
        end if
        COEFX=1.d0
        COEFY=1.d0
        ! CENTRIFUGAL LOAD
        if( LTYPE==5 ) then
          XCOD=dot_product( H(1:NN),XX(1:NN) )
          YCOD=dot_product( H(1:NN),YY(1:NN) )
          HX=AX + ( (XCOD-AX)*RX+(YCOD-AY)*RY )/(RX**2+RY**2)*RX
          HY=AY + ( (XCOD-AX)*RX+(YCOD-AY)*RY )/(RX**2+RY**2)*RY
          PHX=XCOD-HX
          PHY=YCOD-HY
          COEFX=RHO*val*val*PHX
          COEFY=RHO*val*val*PHY
        end if
        do I=1,NN
          PLX(I)=PLX(I)+H(I)*WG*COEFX
          PLY(I)=PLY(I)+H(I)*WG*COEFY
        enddo
      enddo
      if( LTYPE.EQ.1) then
        do I=1,NN
          VECT(2*I-1)=val*PLX(I)
          VECT(2*I  )=0.d0
        enddo
      else if( LTYPE.EQ.2 ) then
        do I=1,NN
          VECT(2*I-1)=0.d0
          VECT(2*I  )=val*PLY(I)
        enddo
      else if( LTYPE.EQ.4 ) then
        VX=PARAMS(1)
        VY=PARAMS(2)
        VX=VX/sqrt( PARAMS(1)**2 + PARAMS(2)**2 )
        VY=VY/sqrt( PARAMS(1)**2 + PARAMS(2)**2 )
        do I=1,NN
          VECT(2*I-1)=val*PLX(I)*RHO*VX
          VECT(2*I  )=val*PLY(I)*RHO*VY
        enddo
      else if( LTYPE==5 ) then
        do I=1,NN
          VECT(2*I-1)=PLX(I)
          VECT(2*I  )=PLY(I)
        enddo
      end if
    endif
  end subroutine DL_C2

  !----------------------------------------------------------------------*
  subroutine TLOAD_C2(ETYPE,NN,XX,YY,TT,T0,gausses,PARAM1,ISET,VECT)
    !----------------------------------------------------------------------*
    !
    ! THERMAL LOAD OF PLANE ELEMENT
    !
    use mMaterial
    use m_MatMatrix
    use m_ElasticLinear
    ! I/F VARIABLES
    integer(kind=kint), intent(in) :: ETYPE,NN
    real(kind=kreal),intent(in)    :: XX(NN),YY(NN),TT(NN),T0(NN),PARAM1
    type(tGaussStatus), intent(in) :: gausses(:)          !< status of qudrature points
    real(kind=kreal),intent(out)   :: VECT(NN*2)
    integer(kind=kint) ISET
    ! LOCAL VARIABLES
    integer(kind=kint) NDOF
    parameter(NDOF=2)
    real(kind=kreal) ALP,PP,D(4,4),B(4,NDOF*NN)
    real(kind=kreal) H(NN)
    real(kind=kreal) EPS(4),SGM(4),localcoord(2)
    real(kind=kreal) deriv(NN,2),DNDE(2,NN)
    real(kind=kreal) XJ(2,2),DET,RR,DUM
    real(kind=kreal) XJI(2,2)
    real(kind=kreal) PAI,THICK,WG
    integer(kind=kint) J,LX
    real(kind=kreal) TEMPC,TEMP0,THERMAL_EPS
    !*************************
    !  CONSTANT
    !*************************
    PAI=4.d0*atan(1.d0)
    ! CLEAR VECT
    VECT(1:NN*NDOF)=0.d0
    ! THICKNESS
    THICK=PARAM1
    !  FOR AX-SYM. ANALYSIS
    if(ISET==2) THICK=1.d0
    ! We suppose material properties doesn't varies inside element
    call calElasticMatrix( gausses(1)%pMaterial, ISET, D, 0.d0  )
    ALP = gausses(1)%pMaterial%variables(M_EXAPNSION)
    pp = gausses(1)%pMaterial%variables(M_POISSON)
    !* LOOP OVER ALL INTEGRATION POINTS
    do LX=1,NumOfQuadPoints( ETYPE )
      call getQuadPoint( ETYPE, LX, localcoord )
      call getShapeFunc( ETYPE, localcoord, H(1:NN) )
      call getShapeDeriv( ETYPE, localcoord, deriv(:,:) )
      XJ(1,1:2)=matmul( XX(1:NN), deriv(1:NN,1:2) )
      XJ(2,1:2)=matmul( YY(1:NN), deriv(1:NN,1:2) )

      DET=XJ(1,1)*XJ(2,2)-XJ(2,1)*XJ(1,2)

      TEMPC=dot_product(H(1:NN),TT(1:NN))
      TEMP0=dot_product(H(1:NN),T0(1:NN))
      if(ISET==2) then
        RR=dot_product(H(1:NN),XX(1:NN))
        WG=getWeight( ETYPE, LX )*DET*RR*2.d0*PAI
      else
        RR=THICK
        WG=getWeight( ETYPE, LX )*DET*RR
      end if
      DUM=1.d0/DET
      XJI(1,1)= XJ(2,2)*DUM
      XJI(1,2)=-XJ(2,1)*DUM
      XJI(2,1)=-XJ(1,2)*DUM
      XJI(2,2)= XJ(1,1)*DUM

      DNDE=matmul( XJI, transpose(deriv) )
      do J=1,NN
        B(1,2*J-1)=DNDE(1,J)
        B(2,2*J-1)=0.d0
        B(3,2*J-1)=DNDE(2,J)
        B(1,2*J  )=0.d0
        B(2,2*J  )=DNDE(2,J)
        B(3,2*J  )=DNDE(1,J)
        B(4,2*J-1)=0.d0
        B(4,2*J  )=0.d0
        if( ISET==2 ) then
          B(4,2*J-1)=H(J)/RR
        endif
      enddo
      !**
      !** THERMAL EPS
      !**
      THERMAL_EPS=ALP*(TEMPC-TEMP0)
      if( ISET .EQ. 2 ) then
        EPS(1)=THERMAL_EPS
        EPS(2)=THERMAL_EPS
        EPS(3)=0.d0
        EPS(4)=THERMAL_EPS
      else if( ISET.EQ.0 ) then
        EPS(1)=THERMAL_EPS*(1.d0+PP)
        EPS(2)=THERMAL_EPS*(1.d0+PP)
        EPS(3)=0.d0
        EPS(4)=0.d0
      else if( ISET.EQ.1 ) then
        EPS(1)=THERMAL_EPS
        EPS(2)=THERMAL_EPS
        EPS(3)=0.d0
        EPS(4)=0.d0
      endif
      !**
      !** SET SGM  {s}=[D]{e}
      !**
      SGM = matmul( EPS(1:4), D(1:4,1:4) )
      !**
      !** CALCULATE LOAD {F}=[B]T{e}
      !**
      VECT(1:NN*NDOF)=VECT(1:NN*NDOF)+matmul( SGM(1:4), B(1:4,1:NN*NDOF) )*WG
    enddo

  end subroutine TLOAD_C2
  !
  !
  !> Update strain and stress inside element
  !---------------------------------------------------------------------*
  subroutine UPDATE_C2( etype,nn,ecoord,gausses,PARAM1,iset,         &
      u, ddu, qf, TT, T0, TN )
    !---------------------------------------------------------------------*
    use m_fstr
    use mMechGauss
    use m_MatMatrix
    ! I/F VARIABLES
    integer(kind=kint), intent(in)     :: etype, nn
    real(kind=kreal),   intent(in)     :: ecoord(3,nn),PARAM1
    integer(kind=kint), intent(in)     :: ISET
    type(tGaussStatus), intent(inout)  :: gausses(:)
    real(kind=kreal),   intent(in)     :: u(2,nn)
    real(kind=kreal),   intent(in)     :: ddu(2,nn)
    real(kind=kreal),   intent(out)    :: qf(:)
    real(kind=kreal),   intent(in)     :: TT(nn)   !< current temperature
    real(kind=kreal),   intent(in)     :: T0(nn)   !< reference temperature
    real(kind=kreal),   intent(in)     :: TN(nn)   !< reference temperature


    integer(kind=kint), parameter :: ndof=2
    real(kind=kreal)   :: D(4,4), B(4,ndof*nn)
    real(kind=kreal)   :: H(nn)
    real(kind=kreal)   :: thick,pai,rr
    real(kind=kreal)   :: det, WG
    real(kind=kreal)   :: localCoord(2)
    real(kind=kreal)   :: gderiv(nn,2), ttc, tt0, ttn
    real(kind=kreal)   :: cdsys(3,3)
    integer(kind=kint) :: j, LX
    real(kind=kreal)   :: totaldisp(2*nn)
    real(kind=kreal)   :: EPSTH(4), alp, alp0, ina(1), outa(1)
    logical            :: ierr
    !
    qf(:)              = 0.d0
    do j=1,nn
      totaldisp((j-1)*2+1) = u(1,j) + ddu(1,j)
      totaldisp(j*2)       = u(2,j) + ddu(2,j)
    enddo
    !
    ! THICKNESS
    THICK=PARAM1
    !  FOR AX-SYM. ANALYSIS
    if(ISET==2) then
      THICK=1.d0
      PAI=4.d0*atan(1.d0)
    endif
    !* LOOP OVER ALL INTEGRATION POINTS
    do LX=1, NumOfQuadPoints(etype)
      call MatlMatrix( gausses(LX), ISET, D, 1.d0, 1.d0, cdsys, 0.d0 )

      call getQuadPoint( etype, LX, localCoord(:) )
      call getGlobalDeriv( etype, nn, localcoord, ecoord, det, gderiv )
      !
      EPSTH = 0.d0
      call getShapeFunc( ETYPE, localcoord, H(:) )
      ttc = dot_product(TT(:), H(:))
      tt0 = dot_product(T0(:), H(:))
      ttn = dot_product(TN(:), H(:))
      if( dabs(ttc-tt0) > 1.d-14 ) then
        ina(1) = ttc
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp = outa(1)
        ina(1) = tt0
        call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
        if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
        alp0 = outa(1)
        EPSTH(1:2)=alp*(ttc-ref_temp)-alp0*(tt0-ref_temp)
      end if
      !
      WG=getWeight( etype, LX )*DET
      if(ISET==2) then
        call getShapeFunc( ETYPE, localcoord, H(:) )
        RR=dot_product( H(1:NN), ecoord(1,1:NN) )
      else
        RR=THICK
        H(:)=0.d0
      end if
      do J=1,NN
        B(1,2*J-1)=gderiv(J,1)
        B(2,2*J-1)=0.d0
        B(3,2*J-1)=gderiv(J,2)
        B(1,2*J  )=0.d0
        B(2,2*J  )=gderiv(J,2)
        B(3,2*J  )=gderiv(J,1)
        B(4,2*J-1)=H(J)/RR
        B(4,2*J  )=0.d0
      enddo

      gausses(LX)%strain(1:4) = matmul( B(:,:), totaldisp )
      gausses(LX)%stress(1:4) = matmul( D(1:4, 1:4), gausses(LX)%strain(1:4)-EPSTH(1:4) )

      !set stress and strain for output
      gausses(LX)%strain_out(1:4) = gausses(LX)%strain(1:4)
      gausses(LX)%stress_out(1:4) = gausses(LX)%stress(1:4)

      !
      !    ----- calculate the Internal Force
      qf(1:nn*ndof) = qf(1:nn*ndof) +                                     &
        matmul( gausses(LX)%stress(1:4), B(1:4,1:nn*ndof) )*WG
      !
    enddo
    !
  end subroutine UPDATE_C2
  !
  !----------------------------------------------------------------------*
  subroutine NodalStress_C2(ETYPE,NN,gausses,ndstrain,ndstress)
    !----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    use mMechGauss

    integer(kind=kint), intent(in) :: ETYPE,NN
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: ndstrain(NN,4)
    real(kind=kreal), intent(out)  :: ndstress(NN,4)

    integer          :: i,ic
    real(kind=kreal) :: TEMP(8)

    TEMP(:)=0.d0
    IC = NumOfQuadPoints(etype)
    do i=1,IC
      TEMP(1:4) = TEMP(1:4) + gausses(i)%strain_out(1:4)
      TEMP(5:8) = TEMP(5:8) + gausses(i)%stress_out(1:4)
    enddo
    TEMP(1:8) = TEMP(1:8)/IC
    do i=1,NN 
      ndstrain(i,1:4) = TEMP(1:4)
      ndstress(i,1:4) = TEMP(5:8)
    enddo

  end subroutine

  !----------------------------------------------------------------------*
  subroutine ElementStress_C2(ETYPE,gausses,strain,stress)
    !----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    use mMechGauss
    integer(kind=kint), intent(in) :: ETYPE
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: strain(4)
    real(kind=kreal), intent(out)  :: stress(4)

    integer          :: i,ic

    strain(:)=0.d0; stress(:)=0.d0
    IC = NumOfQuadPoints(etype)
    do i=1,IC
      strain(:) = strain(:) + gausses(i)%strain_out(1:4)
      stress(:) = stress(:) + gausses(i)%stress_out(1:4)
    enddo
    strain(:) = strain(:)/IC
    stress(:) = stress(:)/IC

  end subroutine

end module m_static_LIB_2d
