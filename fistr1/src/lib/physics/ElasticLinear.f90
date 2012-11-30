!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief  This module provides function on elastic material
module m_ElasticLinear
  use mMaterial

  implicit none
  integer, parameter, private :: kreal = kind(0.0d0)

  contains

!> Calculate elastic matrix
   SUBROUTINE calElasticMatrix( matl, sectType, D, temp  )
   TYPE( tMaterial ), INTENT(IN) :: matl       !> material properties
   INTEGER, INTENT(IN)           :: sectType   !> plane strain/stress or 3D 
   REAL(KIND=kreal), INTENT(OUT) :: D(:,:)     !> elastic matrix
   REAL(KIND=kreal), OPTIONAL    :: temp       !> temprature
   REAL(KIND=kreal) :: EE, PP, COEF1, COEF2, ina(1), outa(2)
   logical :: ierr

   D(:,:)=0.d0
   if( present(temp) ) then
     ina(1) = temp
     call fetch_TableData( MC_ISOELASTIC, matl%dict, outa, ierr, ina )
     if( ierr ) then
	   EE = matl%variables(M_YOUNGS)
       PP = matl%variables(M_POISSON)
     else
       EE = outa(1)
       PP = outa(2)
     endif
   else
     call fetch_TableData( MC_ISOELASTIC, matl%dict, outa, ierr )
     if( ierr ) then
	   EE = matl%variables(M_YOUNGS)
       PP = matl%variables(M_POISSON)
     else
       EE = outa(1)
       PP = outa(2)
     endif
   endif
   
     SELECT CASE (sectType)
     CASE (D3)
       D(1,1)=EE*(1.0-PP)/(1.0-2.0*PP)/(1.0+PP)
       D(1,2)=EE*PP/(1.0-2.0*PP)/(1.0+PP)
       D(1,3)=D(1,2)
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=D(1,2)
       D(3,1)=D(1,3)
       D(3,2)=D(2,3)
       D(3,3)=D(1,1)
       D(4,4)=EE/(1.0+PP)*0.5
       D(5,5)=EE/(1.0+PP)*0.5
       D(6,6)=EE/(1.0+PP)*0.5
     CASE (PlaneStress)
       COEF1=EE/(1.0-PP*PP)
       COEF2=0.5*(1.0-PP)
       D(1,1)=COEF1
       D(1,2)=COEF1*PP
       D(1,3)=0.0
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=0.0
       D(3,1)=0.0
       D(3,2)=0.0
       D(3,3)=COEF1*COEF2
     CASE (PlaneStrain)
       COEF1=EE/((1.0+PP)*(1.0-2.0*PP))
       COEF2=EE/(2.0*(1.0+PP))
       D(1,1)=COEF1*(1.0-PP)
       D(1,2)=COEF1*PP
       D(1,3)=0.0
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=0.0
       D(3,1)=0.0
       D(3,2)=0.0
       D(3,3)=COEF2
     CASE (AxisSymetric)
       COEF1=EE*(1.0-PP)/((1.0+PP)*(1.0-2.0*PP))
       COEF2=(1.0-2.0*PP)/(2.0*(1.0-PP))
       D(1,1)=COEF1
       D(1,2)=COEF1*PP/(1.0-PP)
       D(1,3)=0.0
       D(1,4)=D(1,2)
       D(2,1)=D(1,2)
       D(2,2)=D(1,1)
       D(2,3)=0.0
       D(2,4)=D(1,2)
       D(3,1)=0.0
       D(3,2)=0.0
       D(3,3)=COEF1*COEF2
       D(3,4)=0.0
       D(4,1)=D(1,4)
       D(4,2)=D(2,4)
       D(4,3)=0.0
       D(4,4)=D(1,1)
     CASE (Shell)
     END SELECT

   END SUBROUTINE


      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE LinearElastic_Shell                     &
                 (matl, sectType, c,                     &
                  e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
                  alpha)                                 
!####################################################################
      
      TYPE( tMaterial ), INTENT(IN)   :: matl
      INTEGER, INTENT(IN)             :: sectType
      REAL(KIND = kreal), INTENT(OUT) :: c(:, :, :, :)
      REAL(KIND = kreal), INTENT(IN)  :: e1_hat(3), e2_hat(3), e3_hat(3)
      REAL(KIND = kreal), INTENT(IN)  :: cg1(3), cg2(3), cg3(3)
      REAL(KIND = kreal), INTENT(OUT) :: alpha
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal) :: ee, pp
      REAL(KIND = kreal) :: coef1, coef2
      REAL(KIND = kreal) :: outa(2)
      REAL(KIND = kreal) :: lambda1, lambda2, mu, k_correction
      REAL(KIND = kreal) :: c_hat(3, 3, 3, 3)
      REAL(KIND = kreal) :: e_hat_dot_cg(3, 3)
      REAL(KIND = kreal) :: alpha_over_mu
      
      INTEGER :: i, j, k, l
      
      LOGICAL :: ierr
      
!--------------------------------------------------------------------
      
      CALL fetch_TableData(MC_ISOELASTIC, matl%dict, outa, ierr)
      
!--------------------------------------------------------------------
      
      IF( ierr ) THEN
       
       ee = matl%variables(M_YOUNGS)
       pp = matl%variables(M_POISSON)
       
       alpha_over_mu = matl%variables(M_ALPHA_OVER_MU)
       
      ELSE
       
       ee = outa(1)
       pp = outa(2)
       
       alpha_over_mu = matl%variables(M_ALPHA_OVER_MU)
       
      END IF
      
!--------------------------------------------------------------------
      
      ! Elastic constant
      lambda1 = ee/( 1.0D0-pp*pp )
      lambda2 = pp*lambda1
      mu      = 0.5D0*ee/( 1.0D0+pp )
      
!--------------------------------------------------------------------
      
      ! Shear correction factor
      k_correction = 5.0D0/6.0D0
      
!--------------------------------------------------------------------
      
      alpha = alpha_over_mu*mu
      
!--------------------------------------------------------------------
      
      ! Constitutive tensor
      c_hat(:, :, :, :) = 0.0D0
      
      SELECT CASE( sectType )
      CASE( Shell )
       
       !--------------------------------------------------------
       
       c_hat(1, 1, 1, 1) = lambda1
       c_hat(1, 1, 2, 2) = lambda2
       c_hat(2, 2, 1, 1) = lambda2
       c_hat(2, 2, 2, 2) = lambda1
       c_hat(1, 2, 1, 2) = mu
       c_hat(1, 2, 2, 1) = mu
       c_hat(2, 1, 1, 2) = mu
       c_hat(2, 1, 2, 1) = mu
       c_hat(1, 3, 1, 3) = k_correction*mu
       c_hat(1, 3, 3, 1) = k_correction*mu
       c_hat(2, 3, 2, 3) = k_correction*mu
       c_hat(2, 3, 3, 2) = k_correction*mu
       c_hat(3, 1, 3, 1) = k_correction*mu
       c_hat(3, 1, 1, 3) = k_correction*mu
       c_hat(3, 2, 3, 2) = k_correction*mu
       c_hat(3, 2, 2, 3) = k_correction*mu
       
       !--------------------------------------------------------
       
       e_hat_dot_cg(1, 1)                  &
       = e1_hat(1)*cg1(1)+e1_hat(2)*cg1(2) &
        +e1_hat(3)*cg1(3)                  
       e_hat_dot_cg(2, 1)                  &
       = e2_hat(1)*cg1(1)+e2_hat(2)*cg1(2) &
        +e2_hat(3)*cg1(3)                  
       e_hat_dot_cg(3, 1) = 0.0D0
       e_hat_dot_cg(1, 2)                  &
       = e1_hat(1)*cg2(1)+e1_hat(2)*cg2(2) &
        +e1_hat(3)*cg2(3)                  
       e_hat_dot_cg(2, 2)                  &
       = e2_hat(1)*cg2(1)+e2_hat(2)*cg2(2) &
        +e2_hat(3)*cg2(3)                  
       e_hat_dot_cg(3, 2) = 0.0D0
       e_hat_dot_cg(1, 3)                  &
       = e1_hat(1)*cg3(1)+e1_hat(2)*cg3(2) &
        +e1_hat(3)*cg3(3)                  
       e_hat_dot_cg(2, 3)                  &
       = e2_hat(1)*cg3(1)+e2_hat(2)*cg3(2) &
        +e2_hat(3)*cg3(3)                  
       e_hat_dot_cg(3, 3)                  &
       = e3_hat(1)*cg3(1)+e3_hat(2)*cg3(2) &
        +e3_hat(3)*cg3(3)                  
       
       !--------------------------------------------------------
       
       ! Constitutive tensor
       
       c(1, 1, 1, 1) = 0.0D0
       c(2, 2, 1, 1) = 0.0D0
       c(1, 2, 1, 1) = 0.0D0
       c(2, 2, 2, 2) = 0.0D0
       c(1, 2, 2, 2) = 0.0D0
       c(1, 2, 1, 2) = 0.0D0
       c(3, 1, 1, 1) = 0.0D0
       c(3, 1, 2, 2) = 0.0D0
       c(3, 1, 1, 2) = 0.0D0
       c(2, 3, 1, 1) = 0.0D0
       c(2, 3, 2, 2) = 0.0D0
       c(2, 3, 1, 2) = 0.0D0
       c(3, 1, 3, 1) = 0.0D0
       c(3, 1, 2, 3) = 0.0D0
       c(2, 3, 2, 3) = 0.0D0
       
       DO l = 1, 2
        
        DO k = 1, 2
         
         DO j = 1, 2
          
          DO i = 1, 2
           
           c(1, 1, 1, 1)                            &
           = c(1, 1, 1, 1)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j ,1) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1) 
           c(2, 2, 1, 1)                            &
           = c(2, 2, 1, 1)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 2) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1) 
           c(1, 2, 1, 1)                            &
           = c(1, 2, 1, 1)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j, 2) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1) 
           c(2, 2, 2, 2)                            &
           = c(2, 2, 2, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 2) &
             *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2) 
           
           c(1, 2, 2, 2)                            &
           = c(1, 2, 2, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j, 2) &
             *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2) 
           c(1, 2, 1, 2)                            &
           = c(1, 2, 1, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 1)*e_hat_dot_cg(j, 2) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 2) 
           
          END DO
          
          DO i = 1, 3
           
           c(3, 1, 1, 1)                            &
           = c(3, 1, 1, 1)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1) 
           c(3, 1, 2, 2)                            &
           = c(3, 1, 2, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
             *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2) 
           c(3, 1, 1, 2)                            &
           = c(3, 1, 1, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 2) 
           
          END DO
          
         END DO
         
         DO j = 1, 3
          
          DO i = 1, 2
           
           c(2, 3, 1, 1)                            &
           = c(2, 3, 1, 1)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 1) 
           c(2, 3, 2, 2)                            &
           = c(2, 3, 2, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
             *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 2) 
           c(2, 3, 1, 2)                            &
           = c(2, 3, 1, 2)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
             *e_hat_dot_cg(k, 1)*e_hat_dot_cg(l, 2) 
           
          END DO
          
         END DO
         
        END DO
        
        DO k = 1, 3
         
         DO j = 1, 2
          
          DO i = 1, 3
           
           c(3, 1, 3, 1)                            &
           = c(3, 1, 3, 1)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
             *e_hat_dot_cg(k, 3)*e_hat_dot_cg(l, 1) 
           
          END DO
          
         END DO
         
        END DO
        
       END DO
       
       DO l = 1, 3
        
        DO k = 1, 2
         
         DO j = 1, 2
          
          DO i = 1, 3
           
           c(3, 1, 2, 3)                            &
           = c(3, 1, 2, 3)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 3)*e_hat_dot_cg(j, 1) &
             *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 3) 
           
          END DO
          
         END DO
         
         DO j = 1, 3
          
          DO i = 1, 2
           
           c(2, 3, 2, 3)                            &
           = c(2, 3, 2, 3)                          &
            +c_hat(i, j, k, l)                      &
             *e_hat_dot_cg(i, 2)*e_hat_dot_cg(j, 3) &
             *e_hat_dot_cg(k, 2)*e_hat_dot_cg(l, 3) 
           
          END DO
          
         END DO
         
        END DO
        
       END DO
       
       c(1, 1, 2, 2) = c(2, 2, 1, 1)
       c(1, 1, 1, 2) = c(1, 2, 1, 1)
       c(1, 1, 2, 3) = c(2, 3, 1, 1)
       c(1, 1, 3, 1) = c(3, 1, 1, 1)
       c(2, 2, 1, 2) = c(1, 2, 2, 2)
       c(2, 2, 2, 3) = c(2, 3, 2, 2)
       c(2, 2, 3, 1) = c(3, 1, 2, 2)
       c(1, 2, 2, 3) = c(2, 3, 1, 2)
       c(1, 2, 3, 1) = c(3, 1, 1, 2)
       c(2, 3, 3, 1) = c(3, 1, 2, 3)
       
       !--------------------------------------------------------
       
       ! DO l = 1, 3
       !  
       !  DO k = 1, 3
       !   
       !   DO j = 1, 3
       !    
       !    DO i = 1, 3
       !     
       !     c(i, j, k, l) = 0.0D0
       !     
       !     DO ll = 1, 3
       !      
       !      DO kk = 1, 3
       !       
       !       DO jj = 1, 3
       !        
       !        DO ii = 1, 3
       !         
       !         c(i, j, k, l)                              &
       !         = c(i, j, k, l)                            &
       !          +c_hat(ii, jj, kk, ll)                    &
       !           *e_hat_dot_cg(ii, i)*e_hat_dot_cg(jj, j) &
       !           *e_hat_dot_cg(kk, k)*e_hat_dot_cg(ll, l) 
       !         
       !        END DO
       !        
       !       END DO
       !       
       !      END DO
       !      
       !     END DO
       !     
       !    END DO
       !    
       !   END DO
       !   
       !  END DO
       !  
       ! END DO
       
       !--------------------------------------------------------
       
      END SELECT
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE LinearElastic_Shell
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


end module m_ElasticLinear
