!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!======================================================================!
!                                                                      !
!   Software Name : FrontSOLID                                         !
!         Version : 0.00                                               !
!                                                                      !
!     Last Update : 2009/06/10                                         !
!        Category : FE of mechanical analysis                          !
!                                                                      !
!            Written by Xi YUAN (AdvanceSoft)                          !
!                                                                      !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

!> This modules defines a structure to record history dependent parameter in static analysis
MODULE mMechGauss
  use mMaterial
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: kreal = kind(0.0d0)

! ----------------------------------------------------------------------------
        !> All data should be recorded in every quadrature points
        type tGaussStatus
            type(tMaterial), pointer  :: pMaterial => null()    !< point to material property definition
            real(kind=kreal)          :: strain(6)              !< strain
            real(kind=kreal)          :: stress(6)              !< stress
            integer, pointer          :: istatus(:) =>null()    !< status variables (integer type)
            real(kind=kreal), pointer :: fstatus(:) => null()   !< status variables (double precision type)
            real(kind=kreal)          :: plstrain               !< plastic strain
            real(kind=kreal)          :: ttime                  !< total time at the start of the current increment
            real(kind=kreal)          :: strain_bak(6)          !< strain
            real(kind=kreal)          :: stress_bak(6)              !< stress
        end type

! ----------------------------------------------------------------------------
        !> All data should be recorded in every elements
        type tElement
            integer                     :: etype                 !< element's type
            integer                     :: iset                  !< plane strain, stress etc
            real(kind=kreal), pointer   :: equiForces(:) => null()  !< equivalent forces    
            type(tGaussStatus), pointer :: gausses(:) => null()  !< info of qudrature points
        end type

  CONTAINS

  !> Initializer
  subroutine fstr_init_gauss( gauss )
    type( tGaussStatus ), intent(inout) :: gauss
    integer :: n
    gauss%strain=0.d0; gauss%stress=0.d0
    gauss%strain_bak=0.d0; gauss%stress_bak=0.d0
    gauss%plstrain =0.d0
    gauss%ttime = 0.d0
    if( gauss%pMaterial%mtype==USERMATERIAL ) then
      if( gauss%pMaterial%nfstatus> 0 ) then
        allocate( gauss%fstatus(gauss%pMaterial%nfstatus) )
         gauss%fstatus(:) = 0.d0
      endif
    else if( isElastoplastic(gauss%pMaterial%mtype) ) then
      allocate( gauss%istatus(1) )    ! 0:elastic 1:plastic
      if( isKinematicHarden( gauss%pMaterial%mtype ) ) then
        allocate( gauss%fstatus(7) )  ! plastic strain, back stress
      else
        allocate( gauss%fstatus(2) )    ! plastic strain
      endif
      gauss%istatus = 0
      gauss%fstatus = 0.d0
    else if( isViscoelastic(gauss%pMaterial%mtype) ) then
      n = fetch_TableRow( MC_VISCOELASTIC, gauss%pMaterial%dict )
      if( n>0 ) then
        allocate( gauss%fstatus(12*n+6) )    ! visco stress components
        gauss%fstatus = 0.d0
      else
        stop "Viscoelastic properties not defined"
      endif
    else if( gauss%pMaterial%mtype==NORTON ) then
      allocate( gauss%fstatus(2) )        ! effective stress, effective viscoplastic strain
      gauss%fstatus = 0.d0
      gauss%plstrain = 0.d0
    endif
  end subroutine fstr_init_gauss
  
  !> Finializer
  subroutine fstr_finalize_gauss( gauss )
     type( tGaussStatus ), intent(inout) :: gauss
     if( associated( gauss%istatus ) ) deallocate( gauss%istatus )
     if( associated( gauss%fstatus ) ) deallocate( gauss%fstatus )
  end subroutine

END MODULE


  
