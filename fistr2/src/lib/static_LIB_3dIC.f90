!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by X. YUAN (AdavanceSoft)                         !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  Eight-node hexagonal element with imcompatible mode
!>
!> \see R.L.Taylor,P.J.Bereford, and E.L.Wilson, "A Nonconforming element
!>  for Stress Analysis", Intl. J. Numer. Methods Engng, 10(6), pp1211-1219
!>  ,1976
!!
!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2009/08/12
!>  \version    0.00
!                                                                    !
!======================================================================!
module m_static_LIB_3dIC

   use hecmw, only : kint, kreal
   use elementInfo
   implicit none
   
   private :: INVER9

   contains

!>  CALCULATION STIFF Matrix for C3D8IC ELEMENT
!----------------------------------------------------------------------*
   SUBROUTINE STF_C3D8IC( etype,nn,ecoord,gausses,stiff,nddisp,ehdisp )
!----------------------------------------------------------------------*
     use mMechGauss
     use m_MatMatrix
     use m_static_LIB_3d, only: GEOMAT_C3
     INTEGER(kind=kint), INTENT(IN) :: etype          !< element type, not used here
     INTEGER(kind=kint), INTENT(IN) :: nn             !< number of elements nodes
     REAL(kind=kreal), INTENT(IN)   :: ecoord(3,nn)   !< nodal coord of curr element
     TYPE(tGaussStatus), INTENT(IN) :: gausses(:)     !< Info of qudrature points
     REAL(kind=kreal), INTENT(OUT)  :: stiff(:,:)     !< stiff matrix
     
     REAL(kind=kreal), INTENT(IN), optional :: nddisp(3,nn) !< nodal displacemwent
     REAL(kind=kreal), INTENT(IN), optional :: ehdisp(3,3)  !< enhanced disp of bending mode

     INTEGER(kind=kint)  :: flag           
     INTEGER(kind=kint), PARAMETER :: NDOF=3
     REAL(kind=kreal) D(6,6),B(6,NDOF*(nn+3)),DB(6,NDOF*(nn+3))
     REAL(kind=kreal) gderiv(nn+3,3),stress(6)
     real(kind=kreal) xj(9,9),jacobian(3,3),inverse(3,3)
     REAL(kind=kreal) tmpstiff((nn+3)*3,(nn+3)*3), tmpk(nn*3,9)
     REAL(kind=kreal) det, wg, elem(3,nn), mat(6,6)
     INTEGER(kind=kint) I,J,LX, fetype
     REAL(kind=kreal) naturalCoord(3),unode(3,nn+3)
     REAL(kind=kreal) gdispderiv(3,3)
     REAL(kind=kreal) B1(6,NDOF*(nn+3))
     REAL(kind=kreal) Smat(9,9)
     REAL(kind=kreal) BN(9,NDOF*(nn+3)), SBN(9,NDOF*(NN+3))

     fetype = fe_hex8n
     if( present(nddisp) .and. present(ehdisp) ) then
       unode(:,1:nn)      = nddisp(:,:)
       unode(:,nn+1:nn+3) = ehdisp(:,:)
     endif
	 
	 ! we suppose the same material type in the element
     flag = gausses(1)%pMaterial%nlgeom_flag
     if( .not. present(nddisp) ) flag=INFINITE    ! enforce to infinite deformation analysis
     elem(:,:) = ecoord(:,:)
     if( flag == UPDATELAG ) elem(:,:) = ecoord(:,:) + unode(:,1:nn)
	
     ! --- Inverse of Jacobian at elemental center
     naturalcoord(:) = 0.d0
     CALL getJacobian( fetype, nn, naturalcoord, elem, det, jacobian, inverse )
     inverse(:,:)= inverse(:,:)*det
     ! ---- We now calculate stiff matrix include imcompatible mode
     !    [   Kdd   Kda ]
     !    [   Kad   Kaa ]
     tmpstiff(:,:) = 0.d0
     B1(1:6,1:(nn+3)*NDOF)=0.d0
     BN(1:9,1:(nn+3)*ndof) = 0.d0
     DO LX=1,NumOfQuadPoints(fetype)
       CALL MatlMatrix( gausses(LX), D3, D, 1.d0 )
       IF( flag==UPDATELAG ) then
         call GEOMAT_C3( gausses(LX)%stress, mat )
         D(:,:) = D(:,:)+mat
       ENDIF
       CALL getQuadPoint( fetype, LX, naturalCoord )
       CALL getGlobalDeriv( fetype, nn, naturalcoord, elem, det, gderiv(1:nn,1:3) )
       ! -- Derivative of shape function of imcompatible mode --
       !     [ -2*a   0,   0   ]
       !     [   0,  -2*b, 0   ]
       !     [   0,   0,  -2*c ]
       ! we don't call getShapeDeriv but use above value directly for computation efficiency
       gderiv(nn+1,:) = -2.d0*naturalcoord(1)*inverse(:,1)/det
       gderiv(nn+2,:) = -2.d0*naturalcoord(2)*inverse(:,2)/det
       gderiv(nn+3,:) = -2.d0*naturalcoord(3)*inverse(:,3)/det

       wg=getWeight( fetype, LX )*det
       B(1:6,1:(nn+3)*NDOF)=0.d0
       DO J=1,nn+3
          B(1,3*J-2)=gderiv(J,1)
          B(2,3*J-1)=gderiv(J,2)
          B(3,3*J  )=gderiv(J,3)
          B(4,3*J-2)=gderiv(J,2)
          B(4,3*J-1)=gderiv(J,1)
          B(5,3*J-1)=gderiv(J,3)
          B(5,3*J  )=gderiv(J,2)
          B(6,3*J-2)=gderiv(J,3)
          B(6,3*J  )=gderiv(J,1)
       ENDDO
!
!    calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
       if( flag==TOTALLAG ) then
!      ---dudx(i,j) ==> gdispderiv(i,j)
         gdispderiv(1:ndof,1:ndof) = matmul( unode(1:ndof,1:nn+3), gderiv(1:nn+3,1:ndof) )
         do j=1,nn+3
           B1(1,3*J-2)=gdispderiv(1,1)*gderiv(J,1)
           B1(1,3*J-1)=gdispderiv(2,1)*gderiv(J,1)
           B1(1,3*J  )=gdispderiv(3,1)*gderiv(J,1)
           B1(2,3*J-2)=gdispderiv(1,2)*gderiv(J,2)
           B1(2,3*J-1)=gdispderiv(2,2)*gderiv(J,2)
           B1(2,3*J  )=gdispderiv(3,2)*gderiv(J,2)
           B1(3,3*J-2)=gdispderiv(1,3)*gderiv(J,3)
           B1(3,3*J-1)=gdispderiv(2,3)*gderiv(J,3)
           B1(3,3*J  )=gdispderiv(3,3)*gderiv(J,3)
           B1(4,3*J-2)=gdispderiv(1,2)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,2)
           B1(4,3*J-1)=gdispderiv(2,2)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,2)
           B1(4,3*J  )=gdispderiv(3,2)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,2)
           B1(5,3*J-2)=gdispderiv(1,2)*gderiv(J,3)+gdispderiv(1,3)*gderiv(J,2)
           B1(5,3*J-1)=gdispderiv(2,2)*gderiv(J,3)+gdispderiv(2,3)*gderiv(J,2)
           B1(5,3*J  )=gdispderiv(3,2)*gderiv(J,3)+gdispderiv(3,3)*gderiv(J,2)
           B1(6,3*J-2)=gdispderiv(1,3)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,3)
           B1(6,3*J-1)=gdispderiv(2,3)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,3)
           B1(6,3*J  )=gdispderiv(3,3)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,3)
         enddo
!    ---BL = BL0 + BL1
         do j=1,(nn+3)*ndof
            B(:,j) = B(:,j)+B1(:,j)
         enddo
       endif

        DB(1:6,1:(nn+3)*NDOF) = matmul( D,B(1:6,1:(nn+3)*NDOF) )
        forall( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
         tmpstiff(i,j) = tmpstiff(i,j) + dot_product( B(:,i), DB(:,j) )*wg
        end forall
      ENDDO
      
!    calculate the stress matrix ( TOTAL LAGRANGE METHOD )
      if( flag==TOTALLAG .or. flag==UPDATELAG ) then
        stress(1:6)=gausses(LX)%stress
        do j=1,nn+3
          BN(1,3*J-2) = gderiv(J,1)
          BN(2,3*J-1) = gderiv(J,1)
          BN(3,3*J  ) = gderiv(J,1)
          BN(4,3*J-2) = gderiv(J,2)
          BN(5,3*J-1) = gderiv(J,2)
          BN(6,3*J  ) = gderiv(J,2)
          BN(7,3*J-2) = gderiv(J,3)
          BN(8,3*J-1) = gderiv(J,3)
          BN(9,3*J  ) = gderiv(J,3)
        enddo
        Smat(:,:) = 0.d0
        do j=1,3
          Smat(j  ,j  ) = stress(1)
          Smat(j  ,j+3) = stress(4)
          Smat(j  ,j+6) = stress(6)
          Smat(j+3,j  ) = stress(4)
          Smat(j+3,j+3) = stress(2)
          Smat(j+3,j+6) = stress(5)
          Smat(j+6,j  ) = stress(6)
          Smat(j+6,j+3) = stress(5)
          Smat(j+6,j+6) = stress(3)
        enddo
        SBN(1:9,1:(nn+3)*ndof) = matmul( Smat(1:9,1:9), BN(1:9,1:(nn+3)*ndof) )
        forall( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
          tmpstiff(i,j) = tmpstiff(i,j)+dot_product( BN(:,i), SBN(:,j) )*wg
        end forall
      endif


      ! -----Condense tmpstiff to stiff
      xj(1:9,1:9)= tmpstiff(nn*NDOF+1:(nn+3)*NDOF,nn*NDOF+1:(nn+3)*NDOF)
      call INVER9(xj)
      tmpk = matmul( tmpstiff( 1:nn*NDOF,nn*NDOF+1:(nn+3)*NDOF ), xj )
      stiff = tmpstiff(1:nn*NDOF,1:nn*NDOF)- matmul( tmpk, tmpstiff(nn*NDOF+1:(nn+3)*NDOF,1:nn*NDOF)  )

   end subroutine STF_C3D8IC

!>  Update Strain stress of this element
!----------------------------------------------------------------------*
   SUBROUTINE UpdateST_C3D8IC(ETYPE,NN,XX,YY,ZZ,TT,T0,EDISP,gauss)
!----------------------------------------------------------------------*
      use mMechGauss
      use m_MatMatrix
!
      INTEGER(kind=kint), PARAMETER :: NDOF=3
      INTEGER(kind=kint), INTENT(IN)     :: ETYPE                 !< element type, not used here
      INTEGER(kind=kint), INTENT(IN)     :: NN                    !< number of element nodes
      REAL(kind=kreal), INTENT(IN)       :: XX(NN),YY(NN),ZZ(NN)  !< nodes coordinate of element
      REAL(kind=kreal), INTENT(IN)       :: TT(NN),T0(NN)         !< current and ref temprature
      REAL(kind=kreal), INTENT(IN)       :: EDISP(NN*NDOF)        !< nodal displacement
      TYPE(tGaussStatus), INTENT(INOUT)  :: gauss(:)              !< info about qudrature points

      REAL(kind=kreal) ALP
      REAL(kind=kreal) D(6,6),B(6,NDOF*(NN+3)),DB(6,NDOF*(nn+3))
      REAL(kind=kreal) H(NN)
      REAL(kind=kreal) DET,wg,ecoord(3,NN)
      REAL(kind=kreal) naturalcoord(3),gderiv(NN+3,3)
      INTEGER(kind=kint) J,K,fetype
      REAL(kind=kreal) EPS(6),SGM(6),xj(9,9)
      real(kind=kreal) jacobian(3,3),inverse(3,3)
      REAL(kind=kreal) stiff((nn+3)*3,(nn+3)*3)
      REAL(kind=kreal) tmpforce(9),cdisp((nn+3)*3)
      REAL(kind=kreal) EPSTH(6),THERMAL_EPS,TEMPC,TEMP0
      INTEGER(kind=kint) IC

      fetype = fe_hex8n
      ecoord(1,:)=XX(:)
      ecoord(2,:)=YY(:)
      ecoord(3,:)=ZZ(:)
! ---- Calculate enhanced displacement at first
       ! --- Inverse of Jacobian at elemental center
      naturalcoord(:) = 0.d0
      CALL getJacobian( fetype, nn, naturalcoord, ecoord, det, jacobian, inverse )
      inverse(:,:)= inverse(:,:)*det
      ! ---- We now calculate stiff matrix include imcompatible mode
      !    [   Kdd   Kda ]
      !    [   Kad   Kaa ]
      stiff(:,:) = 0.d0
      B(1:6,1:(nn+3)*NDOF)=0.d0
      DO IC=1,NumOfQuadPoints(fetype)
        CALL MatlMatrix( gauss(IC), D3, D, 1.d0 )
        CALL getQuadPoint( fetype, IC, naturalCoord )
        CALL getGlobalDeriv( fetype, nn, naturalcoord, ecoord, det, gderiv(1:nn,1:3) )
        ! -- Derivative of shape function of imcompatible mode --
        !     [ -2*a   0,   0   ]
        !     [   0,  -2*b, 0   ]
        !     [   0,   0,  -2*c ]
        ! we don't call getShapeDeriv but use above value directly for computation efficiency
        gderiv(nn+1,:) = -2.d0*naturalcoord(1)*inverse(:,1)/det
        gderiv(nn+2,:) = -2.d0*naturalcoord(2)*inverse(:,2)/det
        gderiv(nn+3,:) = -2.d0*naturalcoord(3)*inverse(:,3)/det

        wg=getWeight( fetype, IC )*det
        DO J=1,nn+3
          B(1,3*J-2)=gderiv(J,1)
          B(2,3*J-1)=gderiv(J,2)
          B(3,3*J  )=gderiv(J,3)
          B(4,3*J-2)=gderiv(J,2)
          B(4,3*J-1)=gderiv(J,1)
          B(5,3*J-1)=gderiv(J,3)
          B(5,3*J  )=gderiv(J,2)
          B(6,3*J-2)=gderiv(J,3)
          B(6,3*J  )=gderiv(J,1)
        ENDDO

        DB(1:6,1:(nn+3)*NDOF) = matmul( D,B(1:6,1:(nn+3)*NDOF) )
        stiff(1:(nn+3)*NDOF,1:(nn+3)*NDOF) = stiff(1:(nn+3)*NDOF,1:(nn+3)*NDOF) +        &
               matmul( transpose(B(1:6,1:(nn+3)*NDOF)), DB(1:6,1:(nn+3)*NDOF) )*wg
      ENDDO
      xj(1:9,1:9)= stiff(nn*NDOF+1:(nn+3)*NDOF,nn*NDOF+1:(nn+3)*NDOF)
      call INVER9(xj)
      ! ---  [Kda]*edisp
      tmpforce(:) = matmul( stiff(nn*3+1:(nn+3)*NDOF,1:nn*NDOF), edisp )
      cdisp(1:nn*3)=edisp(:)
      ! ---  -[Kaa]-1 * [Kda] * edisp
      cdisp(nn*3+1:(nn+3)*3)=-matmul( xj(:,:), tmpforce(:) )

! ---- Now strain and stress calculation
      DO IC=1,NumOfQuadPoints(etype)
        ALP = gauss(IC)%pMaterial%variables(4)
        CALL MatlMatrix( gauss(IC), D3, D, 1.d0  )
        CALL getQuadPoint( etype, IC, naturalCoord(:) )
        CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv(1:nn,1:3) )
        CALL getShapeFunc( etype, naturalcoord, H(1:NN) )
        ! -- Derivative of shape function of imcompatible mode --
        !     [ -2*a   0,   0   ]
        !     [   0,  -2*b, 0   ]
        !     [   0,   0,  -2*c ]
        ! we don't call getShapeDeriv but use above value directly for computation efficiency
        gderiv(nn+1,:) = -2.d0*naturalcoord(1)*inverse(:,1)/det
        gderiv(nn+2,:) = -2.d0*naturalcoord(2)*inverse(:,2)/det
        gderiv(nn+3,:) = -2.d0*naturalcoord(3)*inverse(:,3)/det

            B(1:6,1:(NN+3)*NDOF)=0.d0
            DO J=1,NN+3
              B(1,3*J-2)=gderiv(J,1)
              B(2,3*J-1)=gderiv(J,2)
              B(3,3*J  )=gderiv(J,3)
              B(4,3*J-2)=gderiv(J,2)
              B(4,3*J-1)=gderiv(J,1)
              B(5,3*J-1)=gderiv(J,3)
              B(5,3*J  )=gderiv(J,2)
              B(6,3*J-2)=gderiv(J,3)
              B(6,3*J  )=gderiv(J,1)
            ENDDO

            TEMPC=DOT_PRODUCT( H(1:NN),TT(1:NN) )
            TEMP0=DOT_PRODUCT( H(1:NN),T0(1:NN) )
 !**
!** THERMAL EPS
!**
            THERMAL_EPS=ALP*(TEMPC-TEMP0)
            EPSTH(1)=THERMAL_EPS
            EPSTH(2)=THERMAL_EPS
            EPSTH(3)=THERMAL_EPS
            EPSTH(4)=0.d0
            EPSTH(5)=0.d0
            EPSTH(6)=0.d0
!**
!** SET EPS  {e}=[B]{u}
!**
            EPS(1:6) = MATMUL( B(1:6,:), CDISP(:) )
!**
!** SET SGM  {S}=[D]{e}
!**
            DO J=1,6
              SGM(J)=0.d0
              DO K=1,6
                SGM(J)=SGM(J)+D(J,K)*(EPS(K)-EPSTH(K))
              ENDDO
            ENDDO
!**
!** Adding stress in each gauss points
!**
            gauss(IC)%strain(1:6)=EPS(1:6)
            gauss(IC)%stress(1:6)=SGM(1:6)

     ENDDO

   end subroutine UpdateST_C3D8IC
   
!> TO GET INVERSE MATRIX & DETERMINANT ( SIZE=9)
!--------------------------------------------------------------------*
      SUBROUTINE INVER9(A)
!
!  A(NN,NN):ARRAY FOR MATRIX [A]
!  NN      :DIMENSION OF MATRIX [A]
!  IP(N)   :WORK SPACE
!  EPS     :ZERO LIMIT
!  DET     :DETERMINANT OF MATRIX [A]
!*********************************************************************
      INTEGER, PARAMETER :: NN=9
      INTEGER          :: I, J,K,IW,LR,IP(NN)
      REAL(kind=kreal) :: W,WMAX,PIVOT,API,EPS,DET,A(NN,NN)
      DATA EPS/1.0E-35/
      DET=1.d0
      DO 10 I=1,NN
        IP(I)=I
   10 CONTINUE
      DO 70 K=1,NN
        WMAX=0.d0
        DO 20 I=K,NN
          W=DABS(A(I,K))
          IF (W.GT.WMAX) THEN
            WMAX=W
            LR=I
          ENDIF
   20   CONTINUE
        PIVOT=A(LR,K)
        API=ABS(PIVOT)
        IF(API.LE.EPS) THEN
          WRITE(*,'(''PIVOT ERROR AT'',I5)') K
          STOP
        END IF
        DET=DET*PIVOT
        IF (LR.NE.K) THEN
          DET=-DET
          IW=IP(K)
          IP(K)=IP(LR)
          IP(LR)=IW
          DO 30 J=1,NN
            W=A(K,J)
            A(K,J)=A(LR,J)
            A(LR,J)=W
   30     CONTINUE
        ENDIF
        DO 40 I=1,NN
          A(K,I)=A(K,I)/PIVOT
   40   CONTINUE
        DO 60 I=1,NN
          IF (I.NE.K) THEN
            W=A(I,K)
            IF (W.NE.0.) THEN
              DO 50 J=1,NN
                IF (J.NE.K) A(I,J)=A(I,J)-W*A(K,J)
   50         CONTINUE
              A(I,K)=-W/PIVOT
            ENDIF
          ENDIF
   60   CONTINUE
        A(K,K)=1./PIVOT
   70 CONTINUE
      DO 100 I=1,NN
   80   K=IP(I)
        IF (K.NE.I) THEN
          IW=IP(K)
          IP(K)=IP(I)
          IP(I)=IW
          DO 90 J=1,NN
            W=A(J,I)
            A(J,I)=A(J,K)
            A(J,K)=W
   90     CONTINUE
          GOTO 80
        ENDIF
  100 CONTINUE
      RETURN
      
   end subroutine INVER9
   
end module m_static_LIB_3dIC
