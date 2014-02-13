!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!            Written by X. YUAN,                                       !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!>   This module provide common functions of 3D truss elements
!>
!>  \author     X. YUAN (AdavanceSoft)
!>  \date       2012/01/06
!>  \version    0.00
!!
!======================================================================!
module m_static_LIB_1d
   use hecmw, only : kint, kreal
   use elementInfo
   implicit none

   contains
!
!=====================================================================*
!  STF_C1
!=====================================================================*
!>  This subroutine calculate stiff matrix of 2-nodes truss element
  SUBROUTINE STF_C1( etype,nn,ecoord,area,gausses,stiff, u )
    USE mMechGauss
    INTEGER(kind=kint), INTENT(IN)  :: etype               !< element type
    INTEGER(kind=kint), INTENT(IN)  :: nn                  !< number of elemental nodes
    REAL(kind=kreal),   INTENT(IN)  :: ecoord(3,nn)        !< coordinates of elemental nodes
    real(kind=kreal),   INTENT(IN)  :: area                !< section area
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)          !< status of qudrature points
    REAL(kind=kreal),   INTENT(OUT) :: stiff(:,:)          !< stiff matrix
    REAL(kind=kreal),   INTENT(IN), optional :: u(:,:)     !< nodal displacemwent

    REAL(kind=kreal) trans(6,6),stress(6),cdsys(6,6)
    REAL(kind=kreal) DET,WG, llen, elem(3,nn)
    logical :: ierr
    REAL(kind=kreal) outa(1), direc(3), direcxy(3), vn(3)

    stiff(:,:) = 0.d0
    ! we suppose the same material type in the element
    if( present(u) ) then
      elem(:,:) = ecoord(:,:) + u(:,:)
    else
      elem(:,:) = ecoord(:,:)
    endif
	
    direc = elem(:,2)-elem(:,1)
    llen = dsqrt( dot_product(direc, direc) )
    direc = direc/llen
	
    stiff=0.d0
	stiff(1,1)= 1.d0;  stiff(4,1)=-1.d0
    stiff(1,4)=-1.d0;  stiff(4,4)= 1.d0
	
    call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr )
    if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_YOUNGS)
    stiff= outa(1)*area/llen*stiff

    direcxy(:)=0.d0;  direcxy(1)=1.d0
    vn(1) = direcxy(2)*direc(3) - direcxy(3)*direc(2)
    vn(2) = direcxy(3)*direc(1) - direcxy(1)*direc(3)
    vn(3) = direcxy(1)*direc(2) - direcxy(2)*direc(1)  
    if( dsqrt( dot_product(vn,vn) )<1.d-8 )then
      direcxy(:)=0.d0;  direcxy(2)=0
      vn(1) = direc(2)*direcxy(3) - direc(3)*direcxy(2)
      vn(2) = direc(3)*direcxy(1) - direc(1)*direcxy(3)
      vn(3) = direc(1)*direcxy(2) - direc(2)*direcxy(1)  
    endif

    cdsys=0.d0
    cdsys(1,1:3) = direc
    cdsys(2,1:3) = direcxy
    cdsys(3,1:3) = vn
	cdsys(4,4:6) = direc
    cdsys(5,4:6) = direcxy
    cdsys(6,4:6) = vn
    trans = transpose( cdsys )
	
    stiff= matmul( trans, stiff )
    stiff= matmul( stiff, cdsys )

  end subroutine STF_C1

!
!> Update strain and stress inside element
!---------------------------------------------------------------------*
  SUBROUTINE UPDATE_C1( etype,nn,ecoord, area, u, ddu, qf ,gausses )
!---------------------------------------------------------------------*
    USE mMechGauss
! I/F VARIAVLES
    integer(kind=kint), INTENT(IN)     :: etype           !< \param [in] element type
    integer(kind=kint), INTENT(IN)     :: nn              !< \param [in] number of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: ecoord(3,nn)    !< \param [in] coordinates of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: area            !< section area
    real(kind=kreal),   INTENT(IN)     :: u(3,nn)         !< \param [in] nodal dislplacements 
    real(kind=kreal),   INTENT(IN)     :: ddu(3,nn)       !< \param [in] nodal displacement ( solutions of solver )
    real(kind=kreal),   INTENT(OUT)    :: qf(nn*3)        !< \param [out] Internal Force    
    type(tGaussStatus), INTENT(INOUT)  :: gausses(:)      !< \param [out] status of qudrature points


! LCOAL VARIAVLES
    real(kind=kreal)   :: direc(3)
    real(kind=kreal)   :: llen, outa(1)
    integer(kind=kint) :: i
    real(kind=kreal)   :: elem(3,nn)
    real(kind=kreal)   :: dstrain,dstress,edisp(2)
    logical            :: ierr

    qf(:)              = 0.d0
    ! we suppose the same material type in the element
    elem(:,:) = ecoord(:,:) + u(:,:) + ddu(:,:)

    direc = elem(:,2)-elem(:,1)
    llen = dsqrt( dot_product(direc, direc) )
    direc = direc/llen
	
    do i=1,2
      edisp(i) = dot_product( ddu(:,i), direc )
    enddo	  
    dstrain = (edisp(2)-edisp(1))/llen
    gausses(1)%strain(1) = gausses(1)%strain(1) + dstrain

    call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr )
    if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_YOUNGS)
    dstress = outa(1)*dstrain
    gausses(1)%stress(1) = gausses(1)%stress(1) + dstress
	
    qf(1) = gausses(1)%stress(1)*area
    qf(1:3) = -qf(1)*direc
    qf(4:6) = -qf(1:3)

  end subroutine UPDATE_C1
  
  
 !
!> Update strain and stress inside element
!---------------------------------------------------------------------*
  SUBROUTINE UPDATEST_C1( etype,nn,XX,YY,ZZ, area, EDISP, gausses )
!---------------------------------------------------------------------*
    USE mMechGauss
! I/F VARIAVLES
    integer(kind=kint), INTENT(IN)     :: etype           !< \param [in] element type
    integer(kind=kint), INTENT(IN)     :: nn              !< \param [in] number of elemental nodes
    REAL(kind=kreal), INTENT(IN)   :: XX(NN),YY(NN),ZZ(NN)
    real(kind=kreal),   INTENT(IN)     :: area            !< section area
    real(kind=kreal),   INTENT(IN)     :: EDISP(NN*3)         !< \param [in] nodal dislplacements 
    type(tGaussStatus), INTENT(INOUT)  :: gausses(:)      !< \param [out] status of qudrature points


! LCOAL VARIAVLES
    real(kind=kreal)   :: direc(3)
    real(kind=kreal)   :: llen, outa(1)
    integer(kind=kint) :: i
    real(kind=kreal)   :: elem(3,nn), disp(3,nn)
    real(kind=kreal)   :: dstrain,dstress,eedisp(2)
    logical            :: ierr

    elem(1,:)=XX(:)
    elem(2,:)=YY(:)
    elem(3,:)=ZZ(:)
    do i=1,2
	  disp(1,i)=edisp((i-1)*3+1)
      disp(2,i)=edisp((i-1)*3+2)
      disp(3,i)=edisp((i-1)*3+3) 
    enddo

    direc = elem(:,2)-elem(:,1)
    llen = dsqrt( dot_product(direc, direc) )
    direc = direc/llen
	
    do i=1,2
      eedisp(i) = dot_product( disp(:,i), direc )
    enddo	  
    dstrain = (eedisp(2)-eedisp(1))/llen
    gausses(1)%strain(1) = gausses(1)%strain(1) + dstrain

    call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr )
    if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_YOUNGS)
    dstress = outa(1)*dstrain
    gausses(1)%stress(1) = gausses(1)%stress(1) + dstress
	

  end subroutine UPDATEST_C1


!----------------------------------------------------------------------*
   SUBROUTINE NodalStress_C1(ETYPE,NN,gausses,ndstrain,ndstress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: ndstrain(NN,6)
   REAL(kind=kreal), INTENT(OUT)  :: ndstress(NN,6)

     ndstrain(1,1:6) = gausses(1)%strain(1:6)
     ndstress(1,1:6) = gausses(1)%stress(1:6)
     ndstrain(2,1:6) = gausses(1)%strain(1:6)
     ndstress(2,1:6) = gausses(1)%stress(1:6)

   END SUBROUTINE
!
!
!
!----------------------------------------------------------------------*
   SUBROUTINE ElementStress_C1(ETYPE,gausses,strain,stress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: strain(6)
   REAL(kind=kreal), INTENT(OUT)  :: stress(6)


   strain(:) = gausses(1)%strain(1:6)
   stress(:) = gausses(1)%stress(1:6)

   END SUBROUTINE
   
end module m_static_LIB_1d
