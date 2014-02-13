!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Eigen Analysis                                    !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> Set up lumped mass matrix
module m_fstr_EIG_setMASS
use hecmw
use m_eigen_lib
use m_static_get_prop
use elementInfo
contains

!C---------------------------------------------------------------------*
!> Set up lumped mass matrix
      subroutine setMASS(IOUT,hecMESH,hecMAT,myEIG)
!C---------------------------------------------------------------------*
!C*******************************************
!C* SET MASS MATRIX (Lumped mass)
!C*******************************************
      use m_fstr
      !use m_fstr_lib
      !use lczparm
      implicit none
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type(lczparam)            :: myEIG
!C
      integer(kind=kint) nodLOCAL(20),itype,ic_type,icel,isect
      integer(kind=kint) IOUT
      integer(kind=kint) i,j,k,ii,jj,iS,iE,jS,num
      integer(kind=kint) j1,j2,j3,j4
      integer(kind=kint) ix,jx,ll(4),ielm
      integer(kind=kint) pind(20),istart
      integer(kind=kint) npoin,head,nn,nid,numnp,numn,NDOF
      integer(kind=kint) ierror,kk,iax,jk
      real(kind=kreal) xx(20), yy(20), zz(20), ee,pp,rho,rho1,thick,alfa
      real(kind=kreal) xx1(20),yy1(20),zz1(20),Area
      real(kind=kreal) x(20),  y(20),  z(20),  AA,Volume,val
      real(kind=kreal) smax,chkmass
      real(kind=kreal),POINTER :: ss(:)
!C
!*Allocate work array
      ALLOCATE(ss(2000),STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, setMASS"
!C
      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' ****   STAGE Mass   assembly           **'
      ENDIF
!C
      numnp = hecMAT%NP
      numn  = hecMAT%N
      NDOF  = hecMESH%n_dof
!C
      myEIG%mass = 0.
!C
!C**FOR ALL FINITE ELEMENTS
      head   = 0
      smax   = 0.0
!C
      do itype = 1, hecMESH%n_elem_type
        iS = hecMESH%elem_type_index(itype-1) + 1
        iE = hecMESH%elem_type_index(itype  )
        ic_type = hecMESH%elem_type_item(itype)
!C
        nn = getNumberOfNodes(ic_type)
!C
        do icel = iS, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
          enddo
!C
          isect = hecMESH%section_ID(icel)
          iax = hecMESH%section%sect_opt(isect)
          CALL fstr_get_prop(hecMESH,isect,ee,pp,rho,alfa,thick)
!C
          do j = 1,nn
            nid = hecMESH%global_node_ID(nodLOCAL(j))
            xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
            yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
            zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
          enddo
!C
!C Area/volume calculations
!C
          if    ( ic_type.EQ.231 ) then 
            CALL mass_c2d3 ( xx,yy,ee,pp,rho,thick,ss,iax,myEIG )
          elseif( ic_type.EQ.232 ) then 
            CALL mass_c2d6 ( xx,yy,ee,pp,rho,thick,ss,iax,myEIG )
          elseif( ic_type.EQ.241 ) then 
            CALL mass_c2d4 ( xx,yy,ee,pp,rho,thick,ss,iax,myEIG )
          elseif( ic_type.EQ.242 ) then 
            CALL mass_c2d8 ( xx,yy,ee,pp,rho,thick,ss,iax,myEIG )
          elseif( ic_type.EQ.341 ) then 
            CALL mass_c3d4 ( xx,yy,zz,ee,pp,rho,ss,myEIG )
          elseif( ic_type.EQ.342 ) then 
            CALL mass_c3d10( xx,yy,zz,ee,pp,rho,ss,myEIG )
          elseif( ic_type.EQ.351 ) then 
            CALL MASS_C3D6 ( xx,yy,zz,ee,pp,rho,ss,myEIG )
          elseif( ic_type.EQ.352 ) then 
            CALL MASS_C3D15( xx,yy,zz,ee,pp,rho,ss,myEIG )
          elseif( ic_type.EQ.361 ) then 
            CALL mass_c3d8 ( xx,yy,zz,ee,pp,rho,ss,myEIG )
          elseif( ic_type.EQ.362 ) then 
            CALL mass_c3d20( xx,yy,zz,ee,pp,rho,ss,myEIG )
          elseif( ic_type.EQ.731 ) then 
            call face3( xx,yy,zz,AA )
            val = AA/3.0*thick*rho
          elseif( ic_type.EQ.741 ) then 
            call face4( xx,yy,zz,AA )
            val = AA/4.0*thick*rho
          endif
!C
          do j = 1,nn
            nid = nodLOCAL(j)
            js = NDOF*(nid-1)
!C
            if( ic_type.eq.731 .or. ic_type.eq.741 ) then
              myEIG%mass(js+1) = myEIG%mass(js+1) + val
              myEIG%mass(js+2) = myEIG%mass(js+2) + val
              myEIG%mass(js+3) = myEIG%mass(js+3) + val
              myEIG%mass(js+4) = myEIG%mass(js+4) + val*thick**2/12.0
              myEIG%mass(js+5) = myEIG%mass(js+5) + val*thick**2/12.0
              myEIG%mass(js+6) = myEIG%mass(js+6) + val*thick**2/12.0
            else
              jj = NDOF*( j-1 )
              do k = 1, NDOF
                kk = jj + k
                jk = kk*( kk + 1 )/2
                myEIG%mass(js+k) = myEIG%mass(js+k) + ss(jk)
              enddo
            endif
          enddo
!C
        enddo
      enddo
!C
      if    ( NDOF.eq.3 ) then
    !   CALL hecmw_update_3_R( hecMESH,myEIG%mass,numnp )
      elseif( NDOF.eq.2 ) then
    !   CALL hecmw_update_2_R( hecMESH,myEIG%mass,numnp )
      elseif( NDOF.eq.6 ) then
    !   CALL hecmw_update_m_R( hecMESH,myEIG%mass,numnp,NDOF )
      endif
!C
      chkmass = 0.
      do i = 1, numn
        ii = (i-1)*NDOF + 1
        chkmass = chkmass + myEIG%mass(ii)
      end do
!C
    !  CALL hecmw_allreduce_R1(hecMESH,chkmass,hecmw_sum)
      IF(myrank.EQ.0) THEN
        WRITE(IMSG,*) '+===================+'
        WRITE(IMSG,*) 'Total mass: ',chkmass
        WRITE(IMSG,*) '+===================+'
      ENDIF
!C
!C*Deallocate work array
      DEALLOCATE(ss)
!C
      RETURN
      END subroutine setMASS
	  
!C*--------------------------------------------------------------------*
!>  CALCULATION SURFACE AREA FOR 3 POINTS
      SUBROUTINE FACE3(XX,YY,ZZ,AA)
!C*--------------------------------------------------------------------*
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(3),YY(3),ZZ(3)
!C
      A1 = ( XX(2)-XX(1) )**2 + ( YY(2)-YY(1) )**2 + ( ZZ(2)-ZZ(1) )**2
      A2 = ( XX(1)-XX(3) )*( XX(2)-XX(1) ) &
     &   + ( YY(1)-YY(3) )*( YY(2)-YY(1) ) &
     &   + ( ZZ(1)-ZZ(3) )*( ZZ(2)-ZZ(1) )
      A3 = ( XX(3)-XX(1) )**2 + ( YY(3)-YY(1) )**2 + ( ZZ(3)-ZZ(1) )**2
      AA = 0.5*SQRT( A1*A3 - A2**2 )
!C
      RETURN
      END SUBROUTINE FACE3
!C*--------------------------------------------------------------------*
!>  CALCULATION SURFACE AREA FOR 4 POINTS
      SUBROUTINE FACE4(XX,YY,ZZ,AA)
!C*--------------------------------------------------------------------*
      IMPLICIT REAL(kind=kreal)(A-H,O-Z)
      DIMENSION XX(4),YY(4),ZZ(4)
      DIMENSION XG(2),H(4),HR(4),HS(4)
!C
      XG(1) =-0.5773502691896258D0
      XG(2) =-XG(1)
      AA=0.0
      DO 200 LX=1,2
        RI=XG(LX)
        DO 300 LY=1,2
          SI=XG(LY)
          RP=1.0+RI
          SP=1.0+SI
          RM=1.0-RI
          SM=1.0-SI
!C*INTERPOLATION FUNCTION
          H(1)=0.25*RP*SP
          H(2)=0.25*RM*SP
          H(3)=0.25*RM*SM
          H(4)=0.25*RP*SM
!C*DERIVATIVE OF INTERPOLATION FUNCTION
!C*  FOR R-COORDINATE
          HR(1)= .25*SP
          HR(2)=-.25*SP
          HR(3)=-.25*SM
          HR(4)= .25*SM
!C*  FOR S-COORDINATE
          HS(1)= .25*RP
          HS(2)= .25*RM
          HS(3)=-.25*RM
          HS(4)=-.25*RP
!C*JACOBI MATRIX 
          XR=HR(1)*XX(1)+HR(2)*XX(2)+HR(3)*XX(3)+HR(4)*XX(4)
          XS=HS(1)*XX(1)+HS(2)*XX(2)+HS(3)*XX(3)+HS(4)*XX(4)
          YR=HR(1)*YY(1)+HR(2)*YY(2)+HR(3)*YY(3)+HR(4)*YY(4)
          YS=HS(1)*YY(1)+HS(2)*YY(2)+HS(3)*YY(3)+HS(4)*YY(4)
          ZR=HR(1)*ZZ(1)+HR(2)*ZZ(2)+HR(3)*ZZ(3)+HR(4)*ZZ(4)
          ZS=HS(1)*ZZ(1)+HS(2)*ZZ(2)+HS(3)*ZZ(3)+HS(4)*ZZ(4)
          DET=(YR*ZS-ZR*YS)**2+(ZR*XS-XR*ZS)**2+(XR*YS-YR*XS)**2
          DET=SQRT(DET)
          AA=AA+DET
  300   CONTINUE
  200 CONTINUE
!C
      RETURN
      END SUBROUTINE FACE4

end module m_fstr_EIG_setMASS

