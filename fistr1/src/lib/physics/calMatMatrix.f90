!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
!                                                                      !
!      Module Name : lib                                               !
!                                                                      !
!                    Written by Xi YUAN (AdavanceSoft)                 !
!                               K. Satoh (Advancesoft)                 !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!>  \brief   This module manages calculation relates with materials
!!
!>  \author     Xi YUAN (AdavanceSoft), K. Satoh (Advancesoft)
!>  \date       2010/01/12
!>  \version    0.00
!======================================================================!

module m_MatMatrix

  use mMaterial
  use mMechGauss
  use m_ElasticLinear
  use mHyperElastic
  use m_ElastoPlastic
  use mViscoElastic
  use mCreep
  use mUElastic
  use mUmat

  implicit none
  INTEGER, PARAMETER, PRIVATE :: kreal = kind(0.0d0)

  contains

!> Fetch the nlgeom flag of the material
    integer function getNlgeomFlag( gauss )
       type( tGaussStatus ), intent(in) :: gauss      !> status of qudrature point
       getNlgeomFlag = gauss%pMaterial%nlgeom_flag
    end function

!> Calculate constituive matrix
    subroutine MatlMatrix( gauss, sectType, matrix, dt, temperature )
      type( tGaussStatus ), intent(in) :: gauss          !> status of qudrature point
      INTEGER, INTENT(IN)              :: sectType       !> plane strain/stress or 3D
      REAL(KIND=kreal), INTENT(OUT)    :: matrix(:,:)    !> constitutive matrix
      REAL(KIND=kreal), INTENT(IN)     :: dt             !> time increment
      REAL(KIND=kreal), INTENT(IN), optional  :: temperature   !> temperature
      
      real(kind=kreal)            :: cijkl(3,3,3,3)
      TYPE( tMaterial ), pointer  :: matl
      matl=>gauss%pMaterial

      if( matl%mtype==USERELASTIC ) then
        call uElasticMatrix( matl%variables(101:), gauss%strain, matrix )
      elseif( matl%mtype==VISCOELASTIC ) then
        call calViscoelasticMatrix( matl, sectTYPE, dt, matrix )
      elseif( isElastic(matl%mtype) ) then
        if( present( temperature ) ) then
          call calElasticMatrix( matl, sectTYPE, matrix, temperature  )
        else
          call calElasticMatrix( matl, sectTYPE, matrix )
        endif
      elseif( matl%mtype==NEOHOOKE .or. matl%mtype==MOONEYRIVLIN ) then 
        call calElasticMooneyRivlin( matl, sectType, cijkl, gauss%strain  )
        call mat_c2d( cijkl, matrix, sectType )
      elseif( matl%mtype==ARRUDABOYCE )  then
        call calElasticArrudaBoyce( matl, sectType, cijkl, gauss%strain )
        call mat_c2d( cijkl, matrix, sectType )
      elseif( matl%mtype==USERHYPERELASTIC )  then
        call uElasticMatrix( matl%variables(101:), gauss%strain, matrix )
      elseif( isElastoplastic(matl%mtype) )  then
        if( present( temperature ) ) then
          call calElastoPlasticMatrix( matl, sectType, gauss%stress,  &
	           gauss%istatus(1), gauss%fstatus, matrix, temperature  )
        else
          call calElastoPlasticMatrix( matl, sectType, gauss%stress,  &
	         gauss%istatus(1), gauss%fstatus, matrix  )
        endif
      elseif( matl%mtype==USERMATERIAL ) then
        call uMatlMatrix( matl%name, matl%variables(101:), gauss%strain,  &
            gauss%stress, gauss%fstatus, matrix, dt, gauss%ttime )
      elseif( matl%mtype==NORTON ) then
        if( present( temperature ) ) then
          call iso_creep( matl, sectTYPE, gauss%stress, gauss%strain, gauss%fstatus,  &
		       gauss%plstrain, dt, gauss%ttime, matrix, temperature  )
        else
           call iso_creep( matl, sectTYPE, gauss%stress, gauss%strain, gauss%fstatus,  &
		       gauss%plstrain, dt, gauss%ttime, matrix  )
        endif
      else
        stop "Material type not supported!"
      endif
      
    end subroutine

!
!> Update strain and stress for elastic and hyperelastic materials
    subroutine StressUpdate( gauss, sectType, strain, stress, dt )
      type( tGaussStatus ), intent(inout) :: gauss      !> status of qudrature point
      integer, intent(in)                 :: sectType   !> plane strain/stress or 3D
      real(kind=kreal), intent(in)        :: strain(6)  !> strain
      real(kind=kreal), intent(out)       :: stress(6)  !> stress
      real(kind=kreal), intent(in), optional  :: dt         !> time increment

      select case ( gauss%pMaterial%mtype )
        case ( NEOHOOKE, MOONEYRIVLIN )  ! Mooney-Livlin Hyperelastic material
          call calUpdateElasticMooneyRivlin( gauss%pMaterial, sectType, strain, stress )
        case ( ARRUDABOYCE )  ! Arruda-Boyce Hyperelastic material
          call calUpdateElasticArrudaBoyce( gauss%pMaterial, sectType, strain, stress )
        case ( USERHYPERELASTIC, USERELASTIC ) ! user-defined 
          call uElasticUpdate( gauss%pMaterial%variables(101:), strain, stress )
        case ( VISCOELASTIC ) 
          if( .not. present(dt) ) stop "error in viscoelastic update!"
          call UpdateViscoelastic( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, dt )
        case ( NORTON ) 
          if( .not. present(dt)  ) stop "error in viscoelastic update!"
          call update_iso_creep( gauss%pMaterial, sectType, strain, stress, gauss%fstatus,gauss%plstrain, dt, gauss%ttime )
        case ( USERMATERIAL)  ! user-defined 
          call uUpdate(  gauss%pMaterial%name, gauss%pMaterial%variables(101:),   &
             strain, stress, gauss%fstatus, dt, gauss%ttime )      
      end select

    end subroutine StressUpdate

!> Transfer rank 4 constituive matrix to rank 2 form
subroutine mat_c2d( cijkl, dij, itype )
  real(kind=kreal),   intent(in)  :: cijkl(3,3,3,3)
  real(kind=kreal),   intent(out) :: dij(6,6)
  integer,            intent(in)  :: itype

  dij(:,:) = 0.d0
  SELECT CASE( itype )
  CASE( D3 )
    dij(1,1) = cijkl(1,1,1,1)   ! ---
    dij(1,2) = cijkl(1,1,2,2)
    dij(1,3) = cijkl(1,1,3,3)
    dij(1,4) = cijkl(1,1,1,2)
    dij(1,5) = cijkl(1,1,2,3)
    dij(1,6) = cijkl(1,1,3,1)
    dij(2,1) = cijkl(2,2,1,1)   ! ---
    dij(2,2) = cijkl(2,2,2,2)
    dij(2,3) = cijkl(2,2,3,3)
    dij(2,4) = cijkl(2,2,1,2)
    dij(2,5) = cijkl(2,2,2,3)
    dij(2,6) = cijkl(2,2,3,1)
    dij(3,1) = cijkl(3,3,1,1)   ! ---
    dij(3,2) = cijkl(3,3,2,2)
    dij(3,3) = cijkl(3,3,3,3)
    dij(3,4) = cijkl(3,3,1,2)
    dij(3,5) = cijkl(3,3,2,3)
    dij(3,6) = cijkl(3,3,3,1)
    dij(4,1) = cijkl(1,2,1,1)   ! ---
    dij(4,2) = cijkl(1,2,2,2)
    dij(4,3) = cijkl(1,2,3,3)
    dij(4,4) = cijkl(1,2,1,2)
    dij(4,5) = cijkl(1,2,2,3)
    dij(4,6) = cijkl(1,2,3,1)
    dij(5,1) = cijkl(2,3,1,1)   ! ---
    dij(5,2) = cijkl(2,3,2,2)
    dij(5,3) = cijkl(2,3,3,3)
    dij(5,4) = cijkl(2,3,1,2)
    dij(5,5) = cijkl(2,3,2,3)
    dij(5,6) = cijkl(2,3,3,1)
    dij(6,1) = cijkl(3,1,1,1)   ! ---
    dij(6,2) = cijkl(3,1,2,2)
    dij(6,3) = cijkl(3,1,3,3)
    dij(6,4) = cijkl(3,1,1,2)
    dij(6,5) = cijkl(3,1,2,3)
    dij(6,6) = cijkl(3,1,3,1)
!
  CASE( PlaneStress, PlaneStrain )
    dij(1,1) = cijkl(1,1,1,1)   ! ---
    dij(1,2) = cijkl(1,1,2,2)
    dij(1,3) = cijkl(1,1,1,2)
    dij(2,1) = cijkl(2,2,1,1)   ! ---
    dij(2,2) = cijkl(2,2,2,2)
    dij(2,3) = cijkl(2,2,1,2)
    dij(3,1) = cijkl(1,2,1,1)   ! ---
    dij(3,2) = cijkl(1,2,2,2)
    dij(3,3) = cijkl(1,2,1,2)
  CASE( AxisSymetric )
    dij(1,1) = cijkl(1,1,1,1)
    dij(1,2) = cijkl(1,1,2,2)
    dij(1,3) = cijkl(1,1,1,2)
    dij(1,4) = cijkl(1,1,3,3)
    dij(2,1) = cijkl(2,2,1,1)
    dij(2,2) = cijkl(2,2,2,2)
    dij(2,3) = cijkl(2,2,1,2)
    dij(2,4) = cijkl(2,2,3,3)
    dij(3,1) = cijkl(1,2,1,1)
    dij(3,2) = cijkl(1,2,2,2)
    dij(3,3) = cijkl(1,2,1,2)
    dij(3,4) = cijkl(1,2,3,3)
    dij(4,1) = cijkl(3,3,1,1)
    dij(4,2) = cijkl(3,3,2,2)
    dij(4,3) = cijkl(3,3,1,2)
    dij(4,4) = cijkl(3,3,3,3)
  CASE( Shell )
  END SELECT

end subroutine mat_c2d


      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE MatlMatrix_Shell                        &
                 (gauss, sectType, D,                    &
                  e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
                  alpha)                                 
!####################################################################
      
      TYPE(tGaussStatus), INTENT(IN)  :: gauss
      INTEGER, INTENT(IN)             :: sectType
      REAL(KIND = kreal), INTENT(OUT) :: D(:, :)
      REAL(KIND = kreal), INTENT(IN)  :: e1_hat(3), e2_hat(3), e3_hat(3)
      REAL(KIND = kreal), INTENT(IN)  :: cg1(3), cg2(3), cg3(3)
      REAL(KIND = kreal), INTENT(OUT) :: alpha
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal)       :: c(3, 3, 3, 3)
      TYPE(tMaterial), POINTER :: matl
      
!--------------------------------------------------------------------
      
      matl => gauss%pMaterial
      
!--------------------------------------------------------------------
      
      IF( isElastic(matl%mtype) ) THEN
       
       CALL LinearElastic_Shell                     &
            (matl, sectType, c,                     &
             e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
             alpha)                                 
       
       CALL mat_c2d_Shell(c, D, sectType)
       
      ELSE
       
       STOP "Material type not supported!"
       
      END IF
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE MatlMatrix_Shell
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)
      
      
      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE mat_c2d_Shell(c, D, itype)
!####################################################################
      
      REAL(KIND = kreal), INTENT(IN)  :: c(:, :, :, :)
      REAL(KIND = kreal), INTENT(OUT) :: D(:, :)
      INTEGER, INTENT(IN)             :: itype
      
!--------------------------------------------------------------------
      
      INTEGER :: index_i(5), index_j(5), &
                 index_k(5), index_l(5)  
      INTEGER :: i, j, k, l
      INTEGER :: is, js
      
!--------------------------------------------------------------------
      
      index_i(1) = 1
      index_i(2) = 2
      index_i(3) = 1
      index_i(4) = 2
      index_i(5) = 3
      
      index_j(1) = 1
      index_j(2) = 2
      index_j(3) = 2
      index_j(4) = 3
      index_j(5) = 1
      
      index_k(1) = 1
      index_k(2) = 2
      index_k(3) = 1
      index_k(4) = 2
      index_k(5) = 3
       
      index_l(1) = 1
      index_l(2) = 2
      index_l(3) = 2
      index_l(4) = 3
      index_l(5) = 1
      
!--------------------------------------------------------------------
      
      D(:, :) = 0.0D0
      
!--------------------------------------------------------------------
      
      SELECT CASE( itype )
      CASE( Shell )
       
       DO js = 1, 5
        
        DO is = 1, 5
         
         i = index_i(is)
         j = index_j(is)
         k = index_k(js)
         l = index_l(js)
         
         D(is, js) = c(i, j, k, l)
         
        END DO
        
       END DO
       
      END SELECT
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE mat_c2d_Shell
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)

end module m_MatMatrix
