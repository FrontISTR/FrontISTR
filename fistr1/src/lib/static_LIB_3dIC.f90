!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
   use m_utilities
   use elementInfo
   implicit none
   
   contains

!>  CALCULATION STIFF Matrix for C3D8IC ELEMENT
!----------------------------------------------------------------------*
   SUBROUTINE STF_C3D8IC( etype,nn,ecoord,gausses,stiff,nddisp,ehdisp )
!----------------------------------------------------------------------*
     use mMechGauss
     use m_MatMatrix
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
     REAL(kind=kreal) gdispderiv(3,3), coordsys(3,3)
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
       CALL MatlMatrix( gausses(LX), D3, D, 1.d0, coordsys)
       CALL getQuadPoint( fetype, LX, naturalCoord )
       CALL getGlobalDeriv( fetype, nn, naturalcoord, elem, det, gderiv(1:nn,1:3) )
       ! -- Derivative of shape function of imcompatible mode --
       !     [ -2*a   0,   0   ]
       !     [   0,  -2*b, 0   ]
       !     [   0,   0,  -2*c ]
       ! we don't call getShapeDeriv but use above value directly for computation efficiency
       gderiv(nn+1,:) = -2.d0*naturalcoord(1)*inverse(1,:)/det
       gderiv(nn+2,:) = -2.d0*naturalcoord(2)*inverse(2,:)/det
       gderiv(nn+3,:) = -2.d0*naturalcoord(3)*inverse(3,:)/det

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
      call calInverse(9, xj)
      tmpk = matmul( tmpstiff( 1:nn*NDOF,nn*NDOF+1:(nn+3)*NDOF ), xj )
      stiff(1:nn*NDOF,1:nn*NDOF) = tmpstiff(1:nn*NDOF,1:nn*NDOF)- matmul( tmpk, tmpstiff(nn*NDOF+1:(nn+3)*NDOF,1:nn*NDOF)  )

   end subroutine STF_C3D8IC

!>  Update Strain stress of this element
!----------------------------------------------------------------------*
   SUBROUTINE UpdateST_C3D8IC(ETYPE,NN,XX,YY,ZZ,TT,T0,EDISP,gausses,coords)
!----------------------------------------------------------------------*

      use m_fstr
      use mMechGauss
      use m_MatMatrix
!
      INTEGER(kind=kint), PARAMETER :: NDOF=3
      INTEGER(kind=kint), INTENT(IN)     :: ETYPE                 !< element type, not used here
      INTEGER(kind=kint), INTENT(IN)     :: NN                    !< number of element nodes
      REAL(kind=kreal), INTENT(IN)       :: XX(NN),YY(NN),ZZ(NN)  !< nodes coordinate of element
      REAL(kind=kreal), INTENT(IN)       :: TT(NN),T0(NN)         !< current and ref temprature
      REAL(kind=kreal), INTENT(IN)       :: EDISP(NN*NDOF)        !< nodal displacement
      TYPE(tGaussStatus), INTENT(INOUT)  :: gausses(:)            !< info about qudrature points
      REAL(kind=kreal), INTENT(INOUT)    :: coords(3,3)

      REAL(kind=kreal) ALP,ALP0
      REAL(kind=kreal) D(6,6),B(6,NDOF*(NN+3)),DB(6,NDOF*(nn+3))
      REAL(kind=kreal) DET,wg,ecoord(3,NN)
      INTEGER(kind=kint) J,K,IC,cdsys_ID,serr,fetype
      REAL(kind=kreal) EPSA(6),EPSTH(6),SGM(6),H(NN),alpo(3),alpo0(3),coordsys(3,3),xj(9,9)
      REAL(kind=kreal) naturalcoord(3),gderiv(NN+3,3)
      REAL(kind=kreal) jacobian(3,3),inverse(3,3)
      REAL(kind=kreal) stiff((nn+3)*3,(nn+3)*3)
      REAL(kind=kreal) tmpforce(9),cdisp((nn+3)*3)
      REAL(kind=kreal) TEMPC,TEMP0,THERMAL_EPS,tm(6,6),outa(1),ina(1)
      logical   :: ierr, matlaniso

        matlaniso = .false.
        cdsys_ID = gausses(1)%pMaterial%cdsys_ID
        if( cdsys_ID>0 ) then   ! cannot define aniso exapansion when no local coord defined
          ina = TT(1)
          call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
          if( .not. ierr ) matlaniso = .true.
        endif

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
        CALL MatlMatrix( gausses(IC), D3, D, 1.d0, coordsys )
        CALL getQuadPoint( fetype, IC, naturalCoord )
        CALL getGlobalDeriv( fetype, nn, naturalcoord, ecoord, det, gderiv(1:nn,1:3) )
        ! -- Derivative of shape function of imcompatible mode --
        !     [ -2*a   0,   0   ]
        !     [   0,  -2*b, 0   ]
        !     [   0,   0,  -2*c ]
        ! we don't call getShapeDeriv but use above value directly for computation efficiency
        gderiv(nn+1,:) = -2.d0*naturalcoord(1)*inverse(1,:)/det
        gderiv(nn+2,:) = -2.d0*naturalcoord(2)*inverse(2,:)/det
        gderiv(nn+3,:) = -2.d0*naturalcoord(3)*inverse(3,:)/det

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
      call calInverse(9, xj)
      ! ---  [Kda]*edisp
      tmpforce(:) = matmul( stiff(nn*3+1:(nn+3)*NDOF,1:nn*NDOF), edisp )
      cdisp(1:nn*3)=edisp(:)
      ! ---  -[Kaa]-1 * [Kda] * edisp
      cdisp(nn*3+1:(nn+3)*3)=-matmul( xj(:,:), tmpforce(:) )

! ---- Now strain and stress calculation
      DO IC=1,NumOfQuadPoints(etype)
        CALL getQuadPoint( etype, IC, naturalCoord(:) )
        CALL getShapeFunc( etype, naturalcoord, H(1:NN) )
        CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv(1:nn,1:3) )
        ! -- Derivative of shape function of imcompatible mode --
        !     [ -2*a   0,   0   ]
        !     [   0,  -2*b, 0   ]
        !     [   0,   0,  -2*c ]
        ! we don't call getShapeDeriv but use above value directly for computation efficiency
        gderiv(nn+1,:) = -2.d0*naturalcoord(1)*inverse(1,:)/det
        gderiv(nn+2,:) = -2.d0*naturalcoord(2)*inverse(2,:)/det
        gderiv(nn+3,:) = -2.d0*naturalcoord(3)*inverse(3,:)/det

        if( matlaniso ) then
          call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys, serr)
          if( serr==-1 ) stop "Fail to setup local coordinate"
          if( serr==-2 ) then
            write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
          endif
        endif

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

        CALL MatlMatrix( gausses(IC), D3, D, 0.d0, coordsys, TEMPC )

            ina(1) = TEMPC
            if( matlaniso ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
               call fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
               alp = outa(1)
            endif
            ina(1) = TEMP0
            if( matlaniso  ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(IC)%pMaterial%dict, alpo0(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
              call fetch_TableData( MC_THEMOEXP, gausses(IC)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(IC)%pMaterial%variables(M_EXAPNSION)
              alp0 = outa(1)
            endif

!**
!** THERMAL strain
!**
            if( matlaniso ) then
              do j=1,3
                EPSTH(j)=ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
              enddo
              EPSTH(4:6)=0.d0
              call transformation(coordsys, tm) 
              EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord  
              EPSTH(4:6)=EPSTH(4:6)*2.d0 
            else
              THERMAL_EPS=ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
              EPSTH(1:3)=THERMAL_EPS
              EPSTH(4:6)=0.d0
            endif
!**
!** SET EPS  {e}=[B]{u}
!**
            EPSA(1:6) = MATMUL( B(1:6,:), CDISP(:) )
!**
!** SET SGM  {S}=[D]{e}
!**
            DO J=1,6
              SGM(J)=0.d0
              DO K=1,6
                SGM(J)=SGM(J)+D(J,K)*(EPSA(K)-EPSTH(K))
              ENDDO
            ENDDO
!**
!** Adding stress in each gauss points
!**
            gausses(IC)%strain(1:6)=EPSA(1:6)
            gausses(IC)%stress(1:6)=SGM(1:6)

     ENDDO

   end subroutine UpdateST_C3D8IC
   
   
end module m_static_LIB_3dIC
