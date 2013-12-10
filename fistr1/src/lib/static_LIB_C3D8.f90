!======================================================================!
!
!> \brief  This module contains several strategy to free locking problem
!> in Eight-node hexagonal element
!>

!>  \author     Xi YUAN (AdavanceSoft)
!>  \date       2010/04/01
!>  \version    0.00
!                                                                      !
!======================================================================!
module m_static_LIB_C3D8

   use hecmw, only : kint, kreal
   use elementInfo
   implicit none
   

   contains

!>  This subroutine calculate stiff matrix using b-bar method
!> 
!> \see Hughes, T. J. "Generalization of Selective Integration Procedures 
!>  to Anisotropic and Nonlinear Media", Intl. J. Numer. Methods Engng, 15, 
!>  pp1413-1418,1980
!----------------------------------------------------------------------*
   SUBROUTINE STF_C3D8Bbar( etype,nn,ecoord,gausses,stiff, tincr, coords,u, temperature )
    USE mMechGauss
    use m_MatMatrix
    use m_common_struct
    INTEGER(kind=kint), INTENT(IN)  :: etype               !< element type
    INTEGER(kind=kint), INTENT(IN)  :: nn                  !< number of elemental nodes
    REAL(kind=kreal),   INTENT(IN)  :: ecoord(3,nn)        !< coordinates of elemental nodes
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)          !< status of qudrature points
    REAL(kind=kreal),   INTENT(OUT) :: stiff(:,:)          !< stiff matrix
    real(kind=kreal),   intent(in)  :: tincr               !< time increment
    REAL(kind=kreal), INTENT(INOUT) :: coords(3,3)     !< variables to define matreial coordinate system
    REAL(kind=kreal),   INTENT(IN), optional :: u(:,:)     !< nodal displacemwent
    REAL(kind=kreal),   INTENT(IN), optional :: temperature(nn)     !< temperature

    INTEGER             :: flag
    INTEGER(kind=kint), PARAMETER :: NDOF=3
    REAL(kind=kreal) D(6,6),B(6,NDOF*NN),DB(6,NDOF*NN)
    REAL(kind=kreal) gderiv(NN,3),stress(6),mat(6,6)
    REAL(kind=kreal) DET,WG, temp, spfunc(NN)
    INTEGER(kind=kint) I,J,LX, cdsys_ID, serr
    REAL(kind=kreal) naturalCoord(3)
    REAL(kind=kreal) gdispderiv(3,3)
    REAL(kind=kreal) B1(6,NDOF*NN), Bbar(NN,3)
    REAL(kind=kreal) Smat(9,9), elem(3,nn)
    REAL(kind=kreal) BN(9,NDOF*NN), SBN(9,NDOF*NN)
	REAL(kind=kreal) B4,B6,B8,Vol, coordsys(3,3)

    stiff(:,:) = 0.d0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(u) ) flag=INFINITE    ! enforce to infinite deformation analysis
    elem(:,:) = ecoord(:,:)
    if( flag == UPDATELAG ) elem(:,:) = ecoord(:,:) + u(:,:)
	
    cdsys_ID = gausses(1)%pMaterial%cdsys_ID
	
	! dilatation component at centroid
    naturalCoord=0.d0
    CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, Bbar )

    DO LX=1,NumOfQuadPoints(etype)
      CALL getQuadPoint( etype, LX, naturalCoord(:) )
      CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )
	  
      if( cdsys_ID>0 ) then
        call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr)
        if( serr==-1 ) stop "Fail to setup local coordinate"
        if( serr==-2 ) then
           write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        endif
      endif

      if( present(temperature) ) then
        CALL getShapeFunc( etype, naturalcoord, spfunc )
        temp = dot_product( temperature, spfunc )
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, temp )
      else
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys )
      endif
	  
      WG=getWeight( etype, LX )*DET
      B(1:6,1:NN*NDOF)=0.d0
      DO J=1,NN
        B4=(Bbar(J,1)-gderiv(J,1))/3.d0
        B6=(Bbar(J,2)-gderiv(J,2))/3.d0
        B8=(Bbar(J,3)-gderiv(J,3))/3.d0
        B(1,3*J-2)=gderiv(J,1)+B4
        B(1,3*J-1)=B6
        B(1,3*J  )=B8
		
		B(2,3*J-2)=B4
        B(2,3*J-1)=gderiv(J,2)+B6
        B(2,3*J  )=B8
		
		B(3,3*J-2)=B4
        B(3,3*J-1)=B6
        B(3,3*J  )=gderiv(J,3)+B8
		
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
        gdispderiv(1:ndof,1:ndof) = matmul( u(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
        B1(1:6,1:NN*NDOF)=0.d0
        do j=1,nn
          B4=(Bbar(J,1)-gderiv(J,1))/3.d0
          B6=(Bbar(J,2)-gderiv(J,2))/3.d0
          B8=(Bbar(J,3)-gderiv(J,3))/3.d0
          B1(1,3*J-2)=gdispderiv(1,1)*(gderiv(J,1)+B4)
          B1(1,3*J-1)=gdispderiv(2,1)*(gderiv(J,1)+B6)
          B1(1,3*J  )=gdispderiv(3,1)*(gderiv(J,1)+B8)
          B1(2,3*J-2)=gdispderiv(1,2)*(gderiv(J,2)+B4)
          B1(2,3*J-1)=gdispderiv(2,2)*(gderiv(J,2)+B6)
          B1(2,3*J  )=gdispderiv(3,2)*(gderiv(J,2)+B8)
          B1(3,3*J-2)=gdispderiv(1,3)*(gderiv(J,3)+B4)
          B1(3,3*J-1)=gdispderiv(2,3)*(gderiv(J,3)+B6)
          B1(3,3*J  )=gdispderiv(3,3)*(gderiv(J,3)+B8)
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
        do j=1,nn*ndof
          B(:,j) = B(:,j)+B1(:,j)
        enddo
      endif

      DB(1:6,1:NN*NDOF) = matmul( D,B(1:6,1:NN*NDOF) )
      forall( i=1:nn*ndof, j=1:nn*ndof )
        stiff(i,j) = stiff(i,j)+dot_product( B(:,i), DB(:,j) )*WG
      end forall

!    calculate the initial stress matrix 
      if( flag==TOTALLAG .or. flag==UPDATELAG ) then
        stress(1:6)=gausses(LX)%stress
        BN(1:9,1:nn*ndof) = 0.d0
        do j=1,nn
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
        SBN(1:9,1:nn*ndof) = matmul( Smat(1:9,1:9), BN(1:9,1:nn*ndof) )
        forall( i=1:nn*ndof, j=1:nn*ndof )
          stiff(i,j) = stiff(i,j)+dot_product( BN(:,i), SBN(:,j) )*WG
        end forall
      endif
	  
    ENDDO      ! gauss roop
		
   end subroutine STF_C3D8Bbar

!>  Update Strain stress of this element
!----------------------------------------------------------------------*
   SUBROUTINE Update_C3D8Bbar(etype,nn,ecoord, u, du, coords, qf ,gausses, iter, tincr, TT,T0, TN  )
!---------------------------------------------------------------------*
    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use mHyperElastic
    use m_utilities
    integer(kind=kint), INTENT(IN)     :: etype           !< \param [in] element type
    integer(kind=kint), INTENT(IN)     :: nn              !< \param [in] number of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: ecoord(3,nn)    !< \param [in] coordinates of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: u(3,nn)         !< \param [in] nodal dislplacements 
    real(kind=kreal),   INTENT(IN)     :: du(3,nn)        !< \param [in] nodal displacement increment
    REAL(kind=kreal), INTENT(INOUT)    :: coords(3,3)     !< variables to define matreial coordinate system
    real(kind=kreal),   INTENT(OUT)    :: qf(nn*3)        !< \param [out] Internal Force    
    type(tGaussStatus), INTENT(INOUT)  :: gausses(:)      !< \param [out] status of qudrature points
    integer, intent(in)                :: iter
	real(kind=kreal),    intent(in)    :: tincr           !< time increment
    REAL(kind=kreal),   INTENT(IN), optional :: TT(nn)    !< current temperature
    REAL(kind=kreal),   INTENT(IN), optional :: T0(nn)    !< reference temperature
    REAL(kind=kreal),   INTENT(IN), optional :: TN(nn)    !< reference temperature

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof=3
    real(kind=kreal)   :: D(6,6), B(6,ndof*nn), B1(6,ndof*nn)
    real(kind=kreal)   :: gderiv(nn,3), gdispderiv(3,3), det, WG
    integer(kind=kint) :: i, j, k, LX, mtype, cdsys_ID, serr
    real(kind=kreal)   :: naturalCoord(3), rot(3,3), R(3,3), spfunc(nn)
    real(kind=kreal)   :: totaldisp(3,nn), elem(3,nn), elem1(3,nn), coordsys(3,3), tm(6,6)
    real(kind=kreal)   :: dstrain(6),dstress(6),dumstress(3,3),dum(3,3)
    real(kind=kreal)   :: dvol, vol0, Bbar(nn,3), derivdum(1:ndof,1:ndof), BBar2(nn,3)
    real(kind=kreal)   :: B4,B6,B8, ttc,tt0,ttn, alp,alp0, alpo(3),alpo0(3),outa(1),ina(1), EPSTH(6)
    logical            :: ierr, matlaniso

    qf(:)    = 0.d0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    elem(:,:) = ecoord(:,:)
    totaldisp(:,:) = u(:,:)+ du(:,:)
    if( flag == UPDATELAG ) then
        elem(:,:) = (0.5D0*du(:,:)+u(:,:) ) +ecoord(:,:)
        elem1(:,:) = (du(:,:)+u(:,:) ) +ecoord(:,:) 
      !  elem = elem1
        totaldisp(:,:) = du(:,:)
    endif
	
	cdsys_ID = gausses(1)%pMaterial%cdsys_ID
    matlaniso = .false.
    if( cdsys_ID>0 .and. present(TT) ) then
       ina = TT(1)
       call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
       if( .not. ierr ) matlaniso = .true.
    endif

	! dilatation at centroid
	naturalCoord=0.d0
    CALL getGlobalDeriv( etype, nn, naturalcoord, elem, det, Bbar )
    derivdum= matmul( totaldisp(1:ndof,1:nn), Bbar(1:nn,1:ndof) )
    vol0 = (derivdum(1,1)+derivdum(2,2)+derivdum(3,3))/3.d0
    if( flag == UPDATELAG ) CALL getGlobalDeriv( etype, nn, naturalcoord, elem1, det, Bbar2 )

    do LX=1, NumOfQuadPoints(etype)
      mtype = gausses(LX)%pMaterial%mtype
      call getQuadPoint( etype, LX, naturalCoord(:) )
      call getGlobalDeriv( etype, nn, naturalcoord, elem, det, gderiv )
	  
      if( cdsys_ID>0 ) then
        call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr)
        if( serr==-1 ) stop "Fail to setup local coordinate"
        if( serr==-2 ) then
           write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
        endif
      endif
 
	  gdispderiv(1:ndof,1:ndof) = matmul( totaldisp(1:ndof,1:nn), gderiv(1:nn,1:ndof) )
	  dvol = vol0-(gdispderiv(1,1)+gdispderiv(2,2)+gdispderiv(3,3))/3.d0
	  gdispderiv(1,1) = gdispderiv(1,1)+dvol
      gdispderiv(2,2) = gdispderiv(2,2)+dvol
      gdispderiv(3,3) = gdispderiv(3,3)+dvol
!
! ========================================================
!     UPDATE STRAIN and STRESS
! ========================================================

        if( isElastoplastic(gausses(LX)%pMaterial%mtype)    &
          .or. gausses(LX)%pMaterial%mtype==NORTON )        &
           gausses(LX)%pMaterial%mtype = ELASTIC
		   
      EPSTH = 0.d0		   
      if( present(tt) .and. present(t0) ) then
        CALL getShapeFunc( etype, naturalcoord, spfunc )
        ttc = dot_product( TT, spfunc )
        tt0 = dot_product( T0, spfunc )
        ttn = dot_product( TN, spfunc )
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys, ttc )
            
          ina(1) = ttc
          if( matlaniso ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
          else
               call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
               alp = outa(1)
          endif
          ina(1) = tt0
          if( matlaniso  ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
          else
              call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
              if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
              alp0 = outa(1)
          endif
          if( matlaniso ) then
              do j=1,3
                EPSTH(j)=ALPO(j)*(ttc-ref_temp)-alpo0(j)*(tt0-ref_temp)
              enddo
              call transformation(coordsys(:,:), tm)  
              EPSTH(:) = matmul( EPSTH(:), tm  )      ! to global coord 
              EPSTH(4:6)=EPSTH(4:6)*2.d0 
          else
              EPSTH(1:3)=ALP*(ttc-ref_temp)-alp0*(tt0-ref_temp)
          endif

      else
        CALL MatlMatrix( gausses(LX), D3, D, tincr, coordsys )
      endif

!       Small strain
        dstrain(1) = gdispderiv(1,1) 
        dstrain(2) = gdispderiv(2,2) 
        dstrain(3) = gdispderiv(3,3)
        dstrain(4) = ( gdispderiv(1,2)+gdispderiv(2,1) )
        dstrain(5) = ( gdispderiv(2,3)+gdispderiv(3,2) )
        dstrain(6) = ( gdispderiv(3,1)+gdispderiv(1,3) )
        dstrain(:) = dstrain(:)-EPSTH(:)   
		
        if( flag==INFINITE ) then
          gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
          gausses(LX)%stress(1:6) = matmul( D(1:6,1:6), dstrain(1:6) )
          if( isViscoelastic(mtype) .and. tincr/=0.d0 ) then
              if( present(TT) .and. present(T0) ) then
                call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, ttn )
              else
                call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
              endif
              gausses(LX)%stress = real(gausses(LX)%stress)
          endif	
        else if( flag==TOTALLAG ) then
!       Green-Lagrange strain
          dstrain(1) = dstrain(1) +0.5d0*dot_product( gdispderiv(:,1), gdispderiv(:,1) )
          dstrain(2) = dstrain(2) +0.5d0*dot_product( gdispderiv(:,2), gdispderiv(:,2) )
          dstrain(3) = dstrain(3) +0.5d0*dot_product( gdispderiv(:,3), gdispderiv(:,3) )
          dstrain(4) = dstrain(4)+ (gdispderiv(1,1)*gdispderiv(1,2)    &
                      +gdispderiv(2,1)*gdispderiv(2,2)+gdispderiv(3,1)*gdispderiv(3,2) )
          dstrain(5) = dstrain(5)+ (gdispderiv(1,2)*gdispderiv(1,3)    &
                      +gdispderiv(2,2)*gdispderiv(2,3)+gdispderiv(3,2)*gdispderiv(3,3) )
          dstrain(6) = dstrain(6)+ (gdispderiv(1,1)*gdispderiv(1,3)    &
                      +gdispderiv(2,1)*gdispderiv(2,3)+gdispderiv(3,1)*gdispderiv(3,3) )

		  
          if( mtype == NEOHOOKE .or. mtype == MOONEYRIVLIN .or.  mtype == ARRUDABOYCE  .or.      &
            mtype==USERELASTIC .or. mtype==USERHYPERELASTIC .or. mtype==USERMATERIAL ) then
             gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
            call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
          else if( (isViscoelastic(mtype) &
              .or. mtype==NORTON) .and. tincr/=0.d0  ) then
            gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
            gausses(LX)%stress(1:6) = matmul( D(1:6,1:6), dstrain(1:6) )
            gausses(LX)%pMaterial%mtype=mtype
            if( present(TT) .and. present(T0) ) then
              call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr, ttc, ttn )
            else
              call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress, tincr )
            endif
          else
            gausses(LX)%strain(1:6) = dstrain(1:6)+EPSTH(:)
            gausses(LX)%stress(1:6) = matmul( D(1:6,1:6), dstrain(1:6) )
          endif
        else if( flag==UPDATELAG ) then
          rot =0.d0
          rot(1,2)= 0.5d0*(gdispderiv(1,2)-gdispderiv(2,1) );  rot(2,1) = -rot(1,2)
          rot(2,3)= 0.5d0*(gdispderiv(2,3)-gdispderiv(3,2) );  rot(3,2) = -rot(2,3)
          rot(1,3)= 0.5d0*(gdispderiv(1,3)-gdispderiv(3,1) );  rot(3,1) = -rot(1,3)

          gausses(LX)%strain(1:6) = gausses(LX)%strain_bak(1:6)+ dstrain(1:6)+EPSTH(:)
		
          dstress = ( matmul( D(1:6,1:6), dstrain(1:6) ) )
          dumstress(1,1) = gausses(LX)%stress_bak(1)
          dumstress(2,2) = gausses(LX)%stress_bak(2)
          dumstress(3,3) = gausses(LX)%stress_bak(3)
          dumstress(1,2) = gausses(LX)%stress_bak(4);  dumstress(2,1)=dumstress(1,2)
          dumstress(2,3) = gausses(LX)%stress_bak(5);  dumstress(3,2)=dumstress(2,3)
          dumstress(3,1) = gausses(LX)%stress_bak(6);  dumstress(1,3)=dumstress(3,1)
    
          dum(:,:) = matmul( rot,dumstress ) -matmul( dumstress, rot )
          gausses(LX)%stress(1) = gausses(LX)%stress_bak(1)+dstress(1)+ dum(1,1)
          gausses(LX)%stress(2) = gausses(LX)%stress_bak(2)+dstress(2)+ dum(2,2)
          gausses(LX)%stress(3) = gausses(LX)%stress_bak(3)+dstress(3)+ dum(3,3)
          gausses(LX)%stress(4) = gausses(LX)%stress_bak(4)+dstress(4)+ dum(1,2)
          gausses(LX)%stress(5) = gausses(LX)%stress_bak(5)+dstress(5)+ dum(2,3)
          gausses(LX)%stress(6) = gausses(LX)%stress_bak(6)+dstress(6)+ dum(3,1)
		  
          if( mtype==USERMATERIAL ) then
            call StressUpdate( gausses(LX), D3, dstrain, gausses(LX)%stress )
          elseif( mtype==NORTON ) then
	        gausses(LX)%pMaterial%mtype=mtype
		    if( tincr/=0.d0 .and. any(gausses(LX)%stress/=0.d0) ) then
              if( present(TT) .and. present(T0) ) then
                call StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, tincr, ttc, ttn )
              else
                call StressUpdate( gausses(LX), D3, gausses(LX)%strain, gausses(LX)%stress, tincr )
              endif    
            endif
          endif

        endif
        if( isElastoplastic(mtype) ) then
          gausses(LX)%pMaterial%mtype=mtype
          if( present(tt) ) then
            call BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain  &
             , gausses(LX)%istatus(1), gausses(LX)%fstatus, ttc )
          else
            call BackwardEuler( gausses(LX)%pMaterial, gausses(LX)%stress, gausses(LX)%plstrain  &
             , gausses(LX)%istatus(1), gausses(LX)%fstatus )
          endif
        endif
!
! ========================================================
!     calculate the internal force ( equivalent nodal force )
! ========================================================
!       Small strain
      B(1:6, 1:nn*ndof) = 0.0d0
      DO J=1,NN
          B4=(Bbar(J,1)-gderiv(J,1))/3.d0
          B6=(Bbar(J,2)-gderiv(J,2))/3.d0
          B8=(Bbar(J,3)-gderiv(J,3))/3.d0
          B(1,3*J-2)=gderiv(J,1)+B4
          B(1,3*J-1)=B6
          B(1,3*J  )=B8
		
		  B(2,3*J-2)=B4
          B(2,3*J-1)=gderiv(J,2)+B6
          B(2,3*J  )=B8
		
		  B(3,3*J-2)=B4
          B(3,3*J-1)=B6
          B(3,3*J  )=gderiv(J,3)+B8
		  
          B(4,3*J-2) = gderiv(J,2)
          B(4,3*J-1) = gderiv(J,1)
          B(5,3*J-1) = gderiv(J,3)
          B(5,3*J  ) = gderiv(J,2)
          B(6,3*J-2) = gderiv(J,3)
          B(6,3*J  ) = gderiv(J,1)
      ENDDO

      if( flag==INFINITE ) then
	  
      else if( flag==TOTALLAG ) then
        B1(1:6,1:NN*NDOF)=0.d0
        do j=1,nn
          B4=(Bbar(J,1)-gderiv(J,1))/3.d0
          B6=(Bbar(J,2)-gderiv(J,2))/3.d0
          B8=(Bbar(J,3)-gderiv(J,3))/3.d0
          B1(1,3*J-2)=gdispderiv(1,1)*(gderiv(J,1)+B4)
          B1(1,3*J-1)=gdispderiv(2,1)*(gderiv(J,1)+B6)
          B1(1,3*J  )=gdispderiv(3,1)*(gderiv(J,1)+B8)
          B1(2,3*J-2)=gdispderiv(1,2)*(gderiv(J,2)+B4)
          B1(2,3*J-1)=gdispderiv(2,2)*(gderiv(J,2)+B6)
          B1(2,3*J  )=gdispderiv(3,2)*(gderiv(J,2)+B8)
          B1(3,3*J-2)=gdispderiv(1,3)*(gderiv(J,3)+B4)
          B1(3,3*J-1)=gdispderiv(2,3)*(gderiv(J,3)+B6)
          B1(3,3*J  )=gdispderiv(3,3)*(gderiv(J,3)+B8)
          B1(4,3*J-2) = gdispderiv(1,2)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,2)
          B1(4,3*J-1) = gdispderiv(2,2)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,2)
          B1(4,3*J  ) = gdispderiv(3,2)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,2)
          B1(5,3*J-2) = gdispderiv(1,2)*gderiv(J,3)+gdispderiv(1,3)*gderiv(J,2)
          B1(5,3*J-1) = gdispderiv(2,2)*gderiv(J,3)+gdispderiv(2,3)*gderiv(J,2)
          B1(5,3*J  ) = gdispderiv(3,2)*gderiv(J,3)+gdispderiv(3,3)*gderiv(J,2)
          B1(6,3*J-2) = gdispderiv(1,3)*gderiv(J,1)+gdispderiv(1,1)*gderiv(J,3)
          B1(6,3*J-1) = gdispderiv(2,3)*gderiv(J,1)+gdispderiv(2,1)*gderiv(J,3)
          B1(6,3*J  ) = gdispderiv(3,3)*gderiv(J,1)+gdispderiv(3,1)*gderiv(J,3)
        enddo
!       -BL = BL0 + BL1
        do j=1,nn*ndof
            B(:,j) = B(:,j)+B1(:,j)
        enddo
!
      else if( flag == UPDATELAG ) then
        call getGlobalDeriv( etype, nn, naturalcoord, elem1, det, gderiv )
        B(1:6, 1:nn*ndof) = 0.0d0
        DO J=1,NN
          B4=(Bbar2(J,1)-gderiv(J,1))/3.d0
          B6=(Bbar2(J,2)-gderiv(J,2))/3.d0
          B8=(Bbar2(J,3)-gderiv(J,3))/3.d0
          B(1,3*J-2)=gderiv(J,1)+B4
          B(1,3*J-1)=B6
          B(1,3*J  )=B8
		
          B(2,3*J-2)=B4
          B(2,3*J-1)=gderiv(J,2)+B6
          B(2,3*J  )=B8
		
          B(3,3*J-2)=B4
          B(3,3*J-1)=B6
          B(3,3*J  )=gderiv(J,3)+B8
		  
          B(4,3*J-2) = gderiv(J,2)
          B(4,3*J-1) = gderiv(J,1)
          B(5,3*J-1) = gderiv(J,3)
          B(5,3*J  ) = gderiv(J,2)
          B(6,3*J-2) = gderiv(J,3)
          B(6,3*J  ) = gderiv(J,1)
        ENDDO
      endif

!!  calculate the Internal Force 
      WG=getWeight( etype, LX )*DET
      qf(1:nn*ndof) = qf(1:nn*ndof) +                                   &
          matmul( gausses(LX)%stress(1:6), B(1:6,1:nn*ndof) )*WG
    enddo

   end subroutine Update_C3D8Bbar
   
   !> This subroutien calculate thermal loading
!----------------------------------------------------------------------*
   SUBROUTINE TLOAD_C3D8Bbar(ETYPE,NN,XX,YY,ZZ,TT,T0,gausses,VECT,coords)
!----------------------------------------------------------------------*
      use m_fstr
      USE mMechGauss
      use m_MatMatrix
      use m_utilities
      INTEGER(kind=kint), PARAMETER   :: NDOF=3
      INTEGER(kind=kint), INTENT(IN)  :: ETYPE,NN
      TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)          !< status of qudrature points
      REAL(kind=kreal), INTENT(IN)    :: XX(NN),YY(NN),ZZ(NN),TT(NN),T0(NN)
      REAL(kind=kreal), INTENT(OUT)   :: VECT(NN*NDOF)
      REAL(kind=kreal), INTENT(INOUT) :: coords(3,3)       !< variables to define matreial coordinate system

      REAL(kind=kreal) ALP,alp0, D(6,6),B(6,NDOF*NN)
      REAL(kind=kreal) B4,B6,B8,DET,ecoord(3,NN)
      INTEGER(kind=kint) J,LX, cdsys_ID, serr
      REAL(kind=kreal) estrain(6),SGM(6),H(NN)
      REAL(kind=kreal) naturalcoord(3),gderiv(NN,3)
      REAL(kind=kreal) WG, outa(1), ina(1), Bbar(nn,3),alpo(3),alpo0(3), coordsys(3,3)
      REAL(kind=kreal) TEMPC,TEMP0,THERMAL_EPS, tm(6,6)
      logical   :: ierr, matlaniso
	  
        matlaniso = .false.
        cdsys_ID = gausses(1)%pMaterial%cdsys_ID
        if( cdsys_ID>0 ) then   ! cannot define aniso exapansion when no local coord defined
          ina = TT(1)
          call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
          if( .not. ierr ) matlaniso = .true.
        endif

      VECT(:)=0.d0
	  
      ecoord(1,:)=XX(:)
      ecoord(2,:)=YY(:)
      ecoord(3,:)=ZZ(:)
	  
	  naturalCoord=0.d0
      CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, Bbar )

! LOOP FOR INTEGRATION POINTS
      DO LX=1,NumOfQuadPoints(etype)  
        CALL getQuadPoint( etype, LX, naturalCoord(:) )
        CALL getShapeFunc( ETYPE, naturalcoord, H(1:NN) )
        CALL getGlobalDeriv( etype, nn, naturalcoord, ecoord, det, gderiv )
		
        if( matlaniso ) then
          call set_localcoordsys(coords, g_LocalCoordSys(cdsys_ID), coordsys, serr)
          if( serr==-1 ) stop "Fail to setup local coordinate"
          if( serr==-2 ) then
            write(*,*) "WARNING! Cannot setup local coordinate, it is modified automatically"
          endif
        endif

!  WEIGHT VALUE AT GAUSSIAN POINT
        WG=getWeight( etype, LX )*DET
        B(1:6, 1:nn*ndof) = 0.0d0
        DO J=1,NN
          B4=(Bbar(J,1)-gderiv(J,1))/3.d0
          B6=(Bbar(J,2)-gderiv(J,2))/3.d0
          B8=(Bbar(J,3)-gderiv(J,3))/3.d0
          B(1,3*J-2)=gderiv(J,1)+B4
          B(1,3*J-1)=B6
          B(1,3*J  )=B8
		
          B(2,3*J-2)=B4
          B(2,3*J-1)=gderiv(J,2)+B6
          B(2,3*J  )=B8
		
          B(3,3*J-2)=B4
          B(3,3*J-1)=B6
          B(3,3*J  )=gderiv(J,3)+B8
		  
          B(4,3*J-2) = gderiv(J,2)
          B(4,3*J-1) = gderiv(J,1)
          B(5,3*J-1) = gderiv(J,3)
          B(5,3*J  ) = gderiv(J,2)
          B(6,3*J-2) = gderiv(J,3)
          B(6,3*J  ) = gderiv(J,1)
        ENDDO

            TEMPC=DOT_PRODUCT( H(1:NN),TT(1:NN) )
            TEMP0=DOT_PRODUCT( H(1:NN),T0(1:NN) )
			
     !   CALL calElasticMatrix( gausses(LX)%pMaterial, D3, D, TEMPC  )
        CALL MatlMatrix( gausses(LX), D3, D, 0.d0, coordsys, TEMPC )
		
            ina(1) = TEMPC
            if( matlaniso ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
               call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
               alp = outa(1)
            endif
            ina(1) = TEMP0
            if( matlaniso  ) then
               call fetch_TableData( MC_ORTHOEXP, gausses(LX)%pMaterial%dict, alpo0(:), ierr, ina )
               if( ierr ) stop "Fails in fetching orthotropic expansion coefficient!"
            else
               call fetch_TableData( MC_THEMOEXP, gausses(LX)%pMaterial%dict, outa(:), ierr, ina )
               if( ierr ) outa(1) = gausses(LX)%pMaterial%variables(M_EXAPNSION)
               alp0 = outa(1)
            endif

!**
!** THERMAL strain
!**
            if( matlaniso ) then
              do j=1,3
                estrain(j)=ALPO(j)*(TEMPC-ref_temp)-alpo0(j)*(TEMP0-ref_temp)
              enddo
              estrain(4:6)=0.d0
              call transformation(coordsys, tm) 
              estrain(:) = matmul( estrain(:), tm  )      ! to global coord  
              estrain(4:6)=estrain(4:6)*2.d0 
            else
              THERMAL_EPS=ALP*(TEMPC-ref_temp)-alp0*(TEMP0-ref_temp)
              estrain(1:3)=THERMAL_EPS
              estrain(4:6)=0.d0
            endif
			
!**
!** SET SGM  {s}=[D]{e}
!**
            SGM(:)=MATMUL( D(:,:),estrain(:) )
!**
!** CALCULATE LOAD {F}=[B]T{e}
!**
            VECT(1:NN*NDOF)=VECT(1:NN*NDOF)+MATMUL( SGM(:),B(:,:))*WG
      ENDDO
   end subroutine TLOAD_C3D8Bbar
    
end module m_static_LIB_C3D8
