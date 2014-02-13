!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
   !   elseif( matl%mtype==NEOHOOKE ) then
   !     call calElasticNeoHooke( matl, sectType, cijkl, ftn )
   !     call mat_c2d( cijkl, matrix, sectType )
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
            gauss%stress, gauss%fstatus, matrix )
      else
        stop "Material type not supported!"
      endif
      
    end subroutine

!
!> Update strain and stress for elastic and hyperelastic materials
    subroutine StressUpdate( gauss, sectType, strain, stress, strainn, dt )
      type( tGaussStatus ), intent(inout) :: gauss      !> status of qudrature point
      integer, intent(in)                 :: sectType   !> plane strain/stress or 3D
      real(kind=kreal), intent(in)        :: strain(6)  !> strain
      real(kind=kreal), intent(out)       :: stress(6)  !> stress
      real(kind=kreal), intent(in), optional  :: strainn(6) !> strain of last step
      real(kind=kreal), intent(in), optional  :: dt         !> time increment

      select case ( gauss%pMaterial%mtype )
        case ( NEOHOOKE, MOONEYRIVLIN )  ! Mooney-Livlin Hyperelastic material
          call calUpdateElasticMooneyRivlin( gauss%pMaterial, sectType, strain, stress )
        case ( ARRUDABOYCE )  ! Arruda-Boyce Hyperelastic material
          call calUpdateElasticArrudaBoyce( gauss%pMaterial, sectType, strain, stress )
        case ( USERHYPERELASTIC, USERELASTIC ) ! user-defined 
          call uElasticUpdate( gauss%pMaterial%variables(101:), strain, stress )
        case ( VISCOELASTIC ) 
          if( .not. present(strainn) .or. (.not. present(dt) ) ) stop "error in viscoelastic update!"
          if( dt>0 ) &
            call UpdateViscoelastic( gauss%pMaterial, sectType, strain, strainn, stress, gauss%fstatus, dt )
        case ( USERMATERIAL)  ! user-defined 
          call uUpdate(  gauss%pMaterial%name, gauss%pMaterial%variables(101:),   &
             strain, stress, gauss%fstatus )      
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
  CASE( PlaneStress )
  CASE( PlaneStrain )
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

end module m_MatMatrix
