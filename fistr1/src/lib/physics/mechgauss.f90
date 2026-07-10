!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This modules defines a structure to record history dependent parameter in static analysis
module mMechGauss
  use hecmw_util
  use mMaterial
  use elementInfo, only: fe_mitc3_shell, fe_mitc4_shell, fe_mitc9_shell, &
    fe_mitc3_shell361, fe_mitc4_shell361
  use Quadrature, only: gauss1d2, gauss1d3, weight1d2, weight1d3
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
    real(kind=kreal)          :: strain_energy          !< strain energy
    real(kind=kreal)          :: strain_energy_bak      !< strain energy
    real(kind=kreal)          :: plpotential            !< plastic potential
  end type

  ! ----------------------------------------------------------------------------
  !> All data should be recorded in every elements
  type tElement
    integer                     :: etype                 !< element's type
    integer                     :: iset                  !< plane strain, stress etc
    real(kind=kreal), pointer   :: equiForces(:) => null()  !< equivalent forces
    type(tGaussStatus), pointer :: gausses(:) => null()  !< info of qudrature points
    type(tGaussStatus), pointer :: shell_layer_gausses(:) => null()  !< shell layer/thickness quadrature histories
    integer(kind=kint)          :: shell_nlayer = 0       !< number of shell material layers with history
    integer(kind=kint)          :: shell_nthick = 0       !< number of thickness quadrature points per layer
    real(kind=kreal), pointer   :: aux(:,:) => null()    !< nodeless dof for incompatible element
    integer                     :: elemact_flag            !< -1:undefined, 0:active, 1:inactive
    real(kind=kreal)            :: elemact_coeff           !< coefficient to internal force of elemact element
    real(kind=kreal)            :: p(1)                  !< pressure
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
    gauss%strain_energy =0.d0
    gauss%strain_energy_bak =0.d0
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

  !> Number of through-thickness quadrature points used by shell stiffness.
  integer(kind=kint) function fstr_shell_num_thickness_points( etype )
    integer(kind=kint), intent(in) :: etype

    select case( etype )
    case( fe_mitc3_shell, fe_mitc4_shell, fe_mitc3_shell361, fe_mitc4_shell361 )
      fstr_shell_num_thickness_points = 2
    case( fe_mitc9_shell )
      fstr_shell_num_thickness_points = 3
    case default
      fstr_shell_num_thickness_points = 0
    end select
  end function fstr_shell_num_thickness_points

  !> Allocate shell history for every surface Gauss point, layer, and thickness point.
  subroutine fstr_init_shell_layer_gausses( element, ng, nlayer, nthick )
    type( tElement ), intent(inout) :: element
    integer(kind=kint), intent(in)  :: ng, nlayer, nthick
    integer(kind=kint) :: i, nstatus

    if( ng <= 0 .or. nlayer <= 0 .or. nthick <= 0 ) return
    if( .not. associated( element%gausses ) ) return
    if( associated( element%shell_layer_gausses ) ) return

    nstatus = ng*nlayer*nthick
    allocate( element%shell_layer_gausses( nstatus ) )
    element%shell_nlayer = nlayer
    element%shell_nthick = nthick

    do i = 1, nstatus
      element%shell_layer_gausses(i)%pMaterial => element%gausses(1)%pMaterial
      call fstr_init_gauss( element%shell_layer_gausses(i) )
    enddo
  end subroutine fstr_init_shell_layer_gausses

  !> Convert surface Gauss/layer/thickness indices to shell_layer_gausses index.
  integer(kind=kint) function fstr_shell_layer_gauss_index( element, ig, ilayer, ithick )
    type( tElement ), intent(in)    :: element
    integer(kind=kint), intent(in)  :: ig, ilayer, ithick
    integer(kind=kint) :: ngauss

    fstr_shell_layer_gauss_index = 0
    if( .not. associated( element%gausses ) ) return
    if( .not. associated( element%shell_layer_gausses ) ) return
    if( element%shell_nlayer <= 0 .or. element%shell_nthick <= 0 ) return
    if( ilayer < 1 .or. ilayer > element%shell_nlayer ) return
    if( ithick < 1 .or. ithick > element%shell_nthick ) return

    ngauss = size( element%gausses )
    if( ig < 1 .or. ig > ngauss ) return
    if( size( element%shell_layer_gausses ) < ngauss*element%shell_nlayer*element%shell_nthick ) return

    fstr_shell_layer_gauss_index = ((ig-1)*element%shell_nlayer + ilayer-1) &
      *element%shell_nthick + ithick
  end function fstr_shell_layer_gauss_index

  !> Through-thickness quadrature point and weight used by shell elements.
  subroutine fstr_shell_thickness_quadrature( etype, ithick, zeta, weight, ierr )
    integer(kind=kint), intent(in)  :: etype, ithick
    real(kind=kreal), intent(out)   :: zeta, weight
    integer(kind=kint), intent(out) :: ierr

    ierr = 0
    zeta = 0.0d0
    weight = 0.0d0

    select case( etype )
    case( fe_mitc3_shell, fe_mitc4_shell, fe_mitc3_shell361, fe_mitc4_shell361 )
      if( ithick < 1 .or. ithick > 2 ) then
        ierr = 1
        return
      endif
      zeta = gauss1d2(1, ithick)
      weight = weight1d2(ithick)
    case( fe_mitc9_shell )
      if( ithick < 1 .or. ithick > 3 ) then
        ierr = 1
        return
      endif
      zeta = gauss1d3(1, ithick)
      weight = weight1d3(ithick)
    case default
      ierr = 1
    end select
  end subroutine fstr_shell_thickness_quadrature

  !> Layer-local shell thickness coordinate and quadrature weight.
  subroutine fstr_shell_layer_quadrature( element, ilayer, ithick, zeta_layer, weight, ierr )
    type( tElement ), intent(in)    :: element
    integer(kind=kint), intent(in)  :: ilayer, ithick
    real(kind=kreal), intent(out)   :: zeta_layer, weight
    integer(kind=kint), intent(out) :: ierr

    ierr = 0
    zeta_layer = 0.0d0
    weight = 0.0d0

    if( .not. associated( element%gausses ) ) then
      ierr = 1
      return
    endif

    call fstr_shell_layer_quadrature_gauss( element%etype, element%gausses(1), ilayer, ithick, &
      zeta_layer, weight, ierr )
  end subroutine fstr_shell_layer_quadrature

  !> Layer-local shell thickness coordinate and quadrature weight from material status.
  subroutine fstr_shell_layer_quadrature_gauss( etype, gauss, ilayer, ithick, zeta_layer, weight, ierr )
    integer(kind=kint), intent(in)  :: etype, ilayer, ithick
    type( tGaussStatus ), intent(in) :: gauss
    real(kind=kreal), intent(out)   :: zeta_layer, weight
    integer(kind=kint), intent(out) :: ierr
    real(kind=kreal) :: zeta

    ierr = 0
    zeta_layer = 0.0d0
    weight = 0.0d0

    call fstr_shell_thickness_quadrature( etype, ithick, zeta, weight, ierr )
    if( ierr /= 0 ) return
    call fstr_shell_layer_zeta( gauss, ilayer, zeta, zeta_layer, ierr )
  end subroutine fstr_shell_layer_quadrature_gauss

  !> Map layer-local zeta to the whole shell thickness coordinate.
  subroutine fstr_shell_layer_zeta( gauss, ilayer, zeta, zeta_layer, ierr )
    type( tGaussStatus ), intent(in) :: gauss
    integer(kind=kint), intent(in)  :: ilayer
    real(kind=kreal), intent(in)    :: zeta
    real(kind=kreal), intent(out)   :: zeta_layer
    integer(kind=kint), intent(out) :: ierr
    integer(kind=kint) :: i
    real(kind=kreal) :: sumlyr

    ierr = 0
    zeta_layer = 0.0d0

    if( .not. associated( gauss%pMaterial ) ) then
      ierr = 1
      return
    endif
    if( ilayer < 1 .or. ilayer > gauss%pMaterial%totallyr ) then
      ierr = 1
      return
    endif

    sumlyr = 0.0d0
    do i = 1, ilayer
      sumlyr = sumlyr + 2.0d0*gauss%pMaterial%shell_var(i)%weight
    enddo
    zeta_layer = -1.0d0 + sumlyr - gauss%pMaterial%shell_var(ilayer)%weight*(1.0d0-zeta)
  end subroutine fstr_shell_layer_zeta

  !> Release shell layer/thickness history.
  subroutine fstr_finalize_shell_layer_gausses( element )
    type( tElement ), intent(inout) :: element
    integer(kind=kint) :: i

    if( associated( element%shell_layer_gausses ) ) then
      do i = 1, size( element%shell_layer_gausses )
        call fstr_finalize_gauss( element%shell_layer_gausses(i) )
      enddo
      deallocate( element%shell_layer_gausses )
    endif
    element%shell_nlayer = 0
    element%shell_nthick = 0
  end subroutine fstr_finalize_shell_layer_gausses

  !> Copy
  subroutine fstr_copy_gauss( gauss1, gauss2 )
    type( tGaussStatus ), intent(in)    :: gauss1
    type( tGaussStatus ), intent(inout) :: gauss2

    gauss2%strain     = gauss1%strain
    gauss2%stress     = gauss1%stress
    gauss2%strain_bak = gauss1%strain_bak
    gauss2%stress_bak = gauss1%stress_bak
    gauss2%nqm        = gauss1%nqm
    gauss2%strain_out = gauss1%strain_out
    gauss2%stress_out = gauss1%stress_out
    gauss2%plstrain   = gauss1%plstrain
    gauss2%strain_energy = gauss1%strain_energy
    gauss2%strain_energy_bak = gauss1%strain_energy_bak
    gauss2%plpotential = gauss1%plpotential

    if( associated(gauss1%istatus) .and. associated(gauss2%istatus) ) then
      gauss2%istatus   = gauss1%istatus
    end if
    if( associated(gauss1%fstatus) .and. associated(gauss2%fstatus) ) then
      gauss2%fstatus   = gauss1%fstatus
    end if
  end subroutine fstr_copy_gauss

  !> Copy shell layer/thickness history.
  subroutine fstr_copy_shell_layer_gausses( element1, element2 )
    type( tElement ), intent(in)    :: element1
    type( tElement ), intent(inout) :: element2
    integer(kind=kint) :: i

    if( .not. associated( element1%shell_layer_gausses ) ) return
    if( .not. associated( element2%shell_layer_gausses ) ) return
    if( size( element1%shell_layer_gausses ) /= size( element2%shell_layer_gausses ) ) return

    do i = 1, size( element1%shell_layer_gausses )
      call fstr_copy_gauss( element1%shell_layer_gausses(i), element2%shell_layer_gausses(i) )
    enddo
  end subroutine fstr_copy_shell_layer_gausses


end module
