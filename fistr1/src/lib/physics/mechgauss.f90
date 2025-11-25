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
    real(kind=kreal), pointer   :: aux(:,:) => null()    !< nodeless dof for incompatible element
    integer                     :: elemact_flag            !< -1:undefined, 0:active, 1:inactive
    real(kind=kreal)            :: elemact_coeff           !< coefficient to internal force of elemact element
  end type

  type integerBufferType
    integer          :: cap
    integer          :: size
    integer, pointer :: buffer(:) => null()
  end type

  type realBufferType
    integer                   :: cap
    integer                   :: size
    real(kind=kreal), pointer :: buffer(:) => null()
  end type

  type tGaussStatusBufferType
    integer                     :: cap
    integer                     :: size
    type(tGaussStatus), pointer :: buffer(:) => null()
  end type

  type(integerBufferType)      :: integerBuffer(50)
  type(realBufferType)         :: realBuffer(50)
  type(tGaussStatusBufferType) :: tGaussStatusBuffer(50)

contains
  subroutine fstr_module_init_gauss
    integer :: i

    do i = 1, size(integerBuffer)
      integerBuffer(i)%cap = 0
      integerBuffer(i)%size = 0
    enddo
    integerBuffer(1)%cap = 64
    allocate(integerBuffer(1)%buffer(integerBuffer(1)%cap))

    do i = 1, size(realBuffer)
      realBuffer(i)%cap = 0
      realBuffer(i)%size = 0
    enddo
    realBuffer(1)%cap = 64
    allocate(realBuffer(1)%buffer(realBuffer(1)%cap))

    do i = 1, size(tGaussStatusBuffer)
      tGaussStatusBuffer(i)%cap = 0
      tGaussStatusBuffer(i)%size = 0
    enddo
    tGaussStatusBuffer(1)%cap = 64
    allocate(tGaussStatusBuffer(1)%buffer(tGaussStatusBuffer(1)%cap))
  end subroutine

  subroutine fstr_module_finalize_gauss
    integer :: i

    do i = 1, size(integerBuffer)
      if(integerBuffer(i)%cap > 0) then
        deallocate(integerBuffer(i)%buffer)
      endif
    enddo

    do i = 1, size(realBuffer)
      if(realBuffer(i)%cap > 0) then
        deallocate(realBuffer(i)%buffer)
      endif
    enddo

    do i = 1, size(tGaussStatusBuffer)
      if(tGaussStatusBuffer(i)%cap > 0) then
        deallocate(tGaussStatusBuffer(i)%buffer)
      endif
    enddo
  end subroutine

  subroutine fstr_gauss_allocate_integer(ptr, alloc_size)
    integer, pointer, intent(inout) :: ptr(:)
    integer, intent(in) :: alloc_size

    integer :: pre
    integer :: i

    do i = 1, size(integerBuffer)
      if(integerBuffer(i)%cap == 0) then
        integerBuffer(i)%cap = pre * 2
        allocate(integerBuffer(i)%buffer(integerBuffer(i)%cap))
      endif

      if(integerBuffer(i)%cap - integerBuffer(i)%size >= alloc_size) then
        ptr => integerBuffer(i)%buffer(integerBuffer(i)%size+1 : integerBuffer(i)%size+alloc_size)
        integerBuffer(i)%size = integerBuffer(i)%size + alloc_size
        exit
      endif

      pre = integerBuffer(i)%cap
    enddo
  end subroutine

  subroutine fstr_gauss_allocate_real(ptr, alloc_size)
    real(kind=kreal), pointer, intent(inout) :: ptr(:)
    integer, intent(in) :: alloc_size

    integer :: pre
    integer :: i

    do i = 1, size(realBuffer)
      if(realBuffer(i)%cap == 0) then
        realBuffer(i)%cap = pre * 2
        allocate(realBuffer(i)%buffer(realBuffer(i)%cap))
      endif

      if(realBuffer(i)%cap - realBuffer(i)%size >= alloc_size) then
        ptr => realBuffer(i)%buffer(realBuffer(i)%size+1 : realBuffer(i)%size+alloc_size)
        realBuffer(i)%size = realBuffer(i)%size + alloc_size
        exit
      endif

      pre = realBuffer(i)%cap
    enddo
  end subroutine

  subroutine fstr_gauss_allocate_GaussStatus(ptr, alloc_size)
    type(tGaussStatus), pointer, intent(inout) :: ptr(:)
    integer, intent(in) :: alloc_size

    integer :: pre
    integer :: i

    do i = 1, size(tGaussStatusBuffer)
      if(tGaussStatusBuffer(i)%cap == 0) then
        tGaussStatusBuffer(i)%cap = pre * 2
        allocate(tGaussStatusBuffer(i)%buffer(tGaussStatusBuffer(i)%cap))
      endif

      if(tGaussStatusBuffer(i)%cap - tGaussStatusBuffer(i)%size >= alloc_size) then
        ptr => tGaussStatusBuffer(i)%buffer(tGaussStatusBuffer(i)%size+1 : tGaussStatusBuffer(i)%size+alloc_size)
        tGaussStatusBuffer(i)%size = tGaussStatusBuffer(i)%size + alloc_size
        exit
      endif

      pre = tGaussStatusBuffer(i)%cap
    enddo
  end subroutine

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
#ifdef _OPENACC
        call fstr_gauss_allocate_real( gauss%fstatus, gauss%pMaterial%nfstatus )
#else
        allocate( gauss%fstatus(gauss%pMaterial%nfstatus) )
#endif
        gauss%fstatus(:) = 0.d0
      endif
    else if( isElastoplastic(gauss%pMaterial%mtype) ) then
#ifdef _OPENACC
      call fstr_gauss_allocate_integer( gauss%istatus, 1 )    ! 0:elastic 1:plastic
#else
      allocate( gauss%istatus(1) )    ! 0:elastic 1:plastic
#endif
      if( getYieldFunction( gauss%pMaterial%mtype )==3 ) then  ! user defined
        n = uElastoPlasticNumStatus( gauss%pMaterial%variables )
#ifdef _OPENACC
        if( n>0 ) call fstr_gauss_allocate_real( gauss%fstatus, n )
#else
        if( n>0 ) allocate( gauss%fstatus(n) )
#endif
      elseif( isKinematicHarden( gauss%pMaterial%mtype ) ) then
#ifdef _OPENACC
        call fstr_gauss_allocate_real( gauss%fstatus, 7+6 )  ! plastic strain, back stress
#else
        allocate( gauss%fstatus(7+6) )  ! plastic strain, back stress
#endif
      else
#ifdef _OPENACC
        call fstr_gauss_allocate_real( gauss%fstatus, 2 )    ! plastic strain
#else
        allocate( gauss%fstatus(2) )    ! plastic strain
#endif
      endif
      gauss%istatus = 0
      gauss%fstatus = 0.d0
    else if( isViscoelastic(gauss%pMaterial%mtype) ) then
      n = fetch_TableRow( MC_VISCOELASTIC, gauss%pMaterial%dict )
      if( n>0 ) then
#ifdef _OPENACC
        call fstr_gauss_allocate_real( gauss%fstatus, 12*n+6 )    ! visco stress components
#else
        allocate( gauss%fstatus(12*n+6) )    ! visco stress components
#endif
        gauss%fstatus = 0.d0
      else
        stop "Viscoelastic properties not defined"
      endif
    else if( gauss%pMaterial%mtype==NORTON ) then
#ifdef _OPENACC
      call fstr_gauss_allocate_real( gauss%fstatus, 2 )        ! effective stress, effective viscoplastic strain
#else
      allocate( gauss%fstatus(2) )        ! effective stress, effective viscoplastic strain
#endif
      gauss%fstatus = 0.d0
      gauss%plstrain = 0.d0
    endif
  end subroutine fstr_init_gauss

  !> Finializer
  subroutine fstr_finalize_gauss( gauss )
    type( tGaussStatus ), intent(inout) :: gauss
#ifdef _OPENACC
    if( associated( gauss%istatus ) ) nullify( gauss%istatus )
    if( associated( gauss%fstatus ) ) nullify( gauss%fstatus )
#else
    if( associated( gauss%istatus ) ) deallocate( gauss%istatus )
    if( associated( gauss%fstatus ) ) deallocate( gauss%fstatus )
#endif
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
    gauss2%strain_energy = gauss1%strain_energy

    if( associated(gauss1%istatus) .and. associated(gauss2%istatus) ) then
      gauss2%istatus   = gauss1%istatus
    end if
    if( associated(gauss1%fstatus) .and. associated(gauss2%fstatus) ) then
      gauss2%fstatus   = gauss1%fstatus
    end if
  end subroutine fstr_copy_gauss


end module



