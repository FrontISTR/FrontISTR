!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This modules defines a structure to record history dependent parameter in static analysis
module mMechGauss
  use hecmw_util
  use mMaterial
  implicit none

  ! ----------------------------------------------------------------------------
  !> All data should be recorded in every quadrature points
  type tGaussStatus
    type(tMaterial), pointer  :: pMaterial => null()    !< point to material property definition
    real(kind=kreal)          :: strain(6)              !< strain
    real(kind=kreal)          :: stress(6)              !< stress
    integer, pointer          :: istatus(:) =>null()    !< status variables (integer type)
    real(kind=kreal), pointer :: fstatus(:) => null()   !< status variables (double precision type)
    real(kind=kreal)          :: plstrain               !< plastic strain
    real(kind=kreal)          :: strain_bak(6)          !< strain
    real(kind=kreal)          :: stress_bak(6)          !< stress
    real(kind=kreal)          :: nqm(12)                !< NQM
    real(kind=kreal)          :: strain_out(6)          !< strain
    real(kind=kreal)          :: stress_out(6)          !< stress
  end type

  ! ----------------------------------------------------------------------------
  !> All data should be recorded in every elements
  type tElement
    integer                     :: etype                 !< element's type
    integer                     :: iset                  !< plane strain, stress etc
    real(kind=kreal), pointer   :: equiForces(:) => null()  !< equivalent forces
    type(tGaussStatus), pointer :: gausses(:) => null()  !< info of qudrature points
    real(kind=kreal), pointer   :: aux(:,:) => null()    !< nodeless dof for incompatible element
  end type

contains

  !> Initializer
  subroutine fstr_init_gauss( gauss )
    use mUYield, only: uElastoPlasticNumStatus
    type( tGaussStatus ), intent(inout) :: gauss
    integer :: n
    gauss%strain=0.d0; gauss%stress=0.d0
    gauss%strain_bak=0.d0; gauss%stress_bak=0.d0
    gauss%strain_out=0.d0; gauss%stress_out=0.d0
    gauss%plstrain =0.d0
    gauss%nqm =0.d0
    if( gauss%pMaterial%mtype==USERMATERIAL ) then
      if( gauss%pMaterial%nfstatus> 0 ) then
        allocate( gauss%fstatus(gauss%pMaterial%nfstatus) )
        gauss%fstatus(:) = 0.d0
      endif
    else if( isElastoplastic(gauss%pMaterial%mtype) ) then
      allocate( gauss%istatus(1) )    ! 0:elastic 1:plastic
      if( getYieldFunction( gauss%pMaterial%mtype )==3 ) then  ! user defined
        n = uElastoPlasticNumStatus( gauss%pMaterial%variables )
        if( n>0 ) allocate( gauss%fstatus(n) )
      elseif( isKinematicHarden( gauss%pMaterial%mtype ) ) then
        allocate( gauss%fstatus(7+6) )  ! plastic strain, back stress
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

  !> Copy
  subroutine fstr_copy_gauss( gauss1, gauss2 )
    type( tGaussStatus ), intent(in)    :: gauss1
    type( tGaussStatus ), intent(inout) :: gauss2

    gauss2%strain     = gauss1%strain
    gauss2%stress     = gauss1%stress
    gauss2%strain_bak = gauss1%strain_bak
    gauss2%stress_bak = gauss1%stress_bak
    gauss2%plstrain   = gauss1%plstrain

    if( associated(gauss1%istatus) .and. associated(gauss2%istatus) ) then
      gauss2%istatus   = gauss1%istatus
    end if
    if( associated(gauss1%fstatus) .and. associated(gauss2%fstatus) ) then
      gauss2%fstatus   = gauss1%fstatus
    end if
  end subroutine fstr_copy_gauss


end module



