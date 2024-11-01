!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This module manages calculation relates with materials
module m_MatMatrix
  use hecmw_util
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

contains

  !> Fetch the nlgeom flag of the material
  integer function getNlgeomFlag( gauss )
    type( tGaussStatus ), intent(in) :: gauss      !> status of qudrature point
    getNlgeomFlag = gauss%pMaterial%nlgeom_flag
  end function

  !> Calculate constituive matrix
  subroutine MatlMatrix( gauss, sectType, matrix, time, dtime, cdsys, temperature, isEp, hdflag )
    type( tGaussStatus ), intent(in) :: gauss          !> status of qudrature point
    integer, intent(in)              :: sectType       !> plane strain/stress or 3D
    real(kind=kreal), intent(out)    :: matrix(:,:)    !> constitutive matrix
    real(kind=kreal), intent(in)     :: time           !> current time
    real(kind=kreal), intent(in)     :: dtime          !> time increment
    real(kind=kreal), intent(in)     :: cdsys(3,3)     !> material coordinate system
    real(kind=kreal), intent(in)     :: temperature   !> temperature
    integer(kind=kint), intent(in), optional :: isEp
    integer(kind=kint), intent(in), optional :: hdflag  !> return only hyd and dev term if specified

    integer :: i
    integer :: flag, hdflag_in
    real(kind=kreal)            :: cijkl(3,3,3,3)
    type( tMaterial ), pointer  :: matl
    matl=>gauss%pMaterial

    flag = 0
    if( present(isEp) )then
      if( isEp == 1 )flag = 1
    endif

    hdflag_in = 0
    if( present(hdflag) ) hdflag_in = hdflag

    if( matl%mtype==USERELASTIC ) then
      call uElasticMatrix( matl%variables(101:), gauss%strain, matrix )
    elseif( isViscoelastic(matl%mtype) ) then
      call calViscoelasticMatrix( matl, sectTYPE, dtime, matrix, temperature )
    elseif( isElastic(matl%mtype) .or. flag==1 ) then
      if(flag==1)then
        i = getElasticType(ELASTIC)
      else
        i = getElasticType(gauss%pMaterial%mtype)
      endif

      if( i==0 ) then
        call calElasticMatrix( matl, sectTYPE, matrix, temperature, hdflag=hdflag_in )
      elseif(  i==1 ) then
        call calElasticMatrix_ortho( gauss%pMaterial, sectTYPE, cdsys, matrix, temperature )
      else
        print *, "Elasticity type", matl%mtype, "not supported"
        stop
      endif

    elseif( matl%mtype==NEOHOOKE .or. matl%mtype==MOONEYRIVLIN ) then
      call calElasticMooneyRivlin( matl, sectType, cijkl, gauss%strain, hdflag=hdflag_in )
      call mat_c2d( cijkl, matrix, sectType )
    elseif( matl%mtype==ARRUDABOYCE )  then
      call calElasticArrudaBoyce( matl, sectType, cijkl, gauss%strain )
      call mat_c2d( cijkl, matrix, sectType )
    elseif( matl%mtype==MOONEYRIVLIN_ANISO ) then
      call calElasticMooneyRivlinAniso( matl, sectType, cijkl, gauss%strain, cdsys, hdflag=hdflag_in )
      call mat_c2d( cijkl, matrix, sectType )
    elseif( matl%mtype==USERHYPERELASTIC )  then
      call uElasticMatrix( matl%variables(101:), gauss%strain, matrix )
    elseif( isElastoplastic(matl%mtype) )  then
      call calElastoPlasticMatrix( matl, sectType, gauss%stress,  &
        gauss%istatus(1), gauss%fstatus, gauss%plstrain, matrix, temperature, hdflag=hdflag_in  )
    elseif( matl%mtype==USERMATERIAL ) then
      call uMatlMatrix( matl%name, matl%variables(101:), gauss%strain,  &
        gauss%stress, gauss%fstatus, matrix, dtime, time )
    elseif( matl%mtype==NORTON ) then
      call iso_creep( matl, sectTYPE, gauss%stress, gauss%strain, gauss%fstatus,  &
        gauss%plstrain, dtime, time, matrix, temperature, hdflag=hdflag_in )
    else
      stop "Material type not supported!"
    endif

  end subroutine

  !
  !> Update strain and stress for elastic and hyperelastic materials
  subroutine StressUpdate( gauss, sectType, strain, stress, cdsys, time, dtime, temp, tempn, hdflag  )
    type( tGaussStatus ), intent(inout) :: gauss      !> status of qudrature point
    integer, intent(in)                 :: sectType   !> plane strain/stress or 3D
    real(kind=kreal), intent(in)        :: strain(6)  !> strain
    real(kind=kreal), intent(out)       :: stress(6)  !> stress
    real(kind=kreal), intent(in)        :: cdsys(3,3) !> material coordinate system
    real(kind=kreal), intent(in), optional  :: time   !> current time
    real(kind=kreal), intent(in), optional  :: dtime  !> time increment
    real(kind=kreal), intent(in)        :: temp       !> current temperature
    real(kind=kreal), intent(in)        :: tempn      !> temperature at last step
    integer(kind=kint), intent(in)      :: hdflag  !> return only hyd and dev term if specified

    if( gauss%pMaterial%mtype==NEOHOOKE .or. gauss%pMaterial%mtype==MOONEYRIVLIN ) then
      call calUpdateElasticMooneyRivlin( gauss%pMaterial, sectType, strain, stress, hdflag=hdflag )
    elseif( gauss%pMaterial%mtype==ARRUDABOYCE ) then ! Arruda-Boyce Hyperelastic material
      call calUpdateElasticArrudaBoyce( gauss%pMaterial, sectType, strain, stress )
    elseif( gauss%pMaterial%mtype==MOONEYRIVLIN_ANISO ) then
      call calUpdateElasticMooneyRivlinAniso( gauss%pMaterial, sectType, strain, stress, cdsys, hdflag=hdflag )
    elseif( gauss%pMaterial%mtype==USERHYPERELASTIC .or. gauss%pMaterial%mtype==USERELASTIC ) then ! user-defined
      call uElasticUpdate( gauss%pMaterial%variables(101:), strain, stress )
    elseif( isViscoelastic( gauss%pMaterial%mtype) ) then
      if( .not. present(dtime) ) stop "error in viscoelastic update!"
      call UpdateViscoelastic( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, dtime, temp, tempn )
    elseif ( gauss%pMaterial%mtype==NORTON ) then
      if( .not. present(dtime)  ) stop "error in viscoelastic update!"
      call update_iso_creep( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, &
        &  gauss%plstrain, dtime, time, temp, hdflag=hdflag )
    elseif ( gauss%pMaterial%mtype==USERMATERIAL)  then ! user-defined
      call uUpdate(  gauss%pMaterial%name, gauss%pMaterial%variables(101:),   &
        strain, stress, gauss%fstatus, dtime, time )
    end if

  end subroutine StressUpdate

  !> Transfer rank 4 constituive matrix to rank 2 form
  subroutine mat_c2d( cijkl, dij, itype )
    real(kind=kreal),   intent(in)  :: cijkl(3,3,3,3)
    real(kind=kreal),   intent(out) :: dij(6,6)
    integer,            intent(in)  :: itype

    dij(:,:) = 0.d0
    select case( itype )
      case( D3 )
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
      case( PlaneStress, PlaneStrain )
        dij(1,1) = cijkl(1,1,1,1)   ! ---
        dij(1,2) = cijkl(1,1,2,2)
        dij(1,3) = cijkl(1,1,1,2)
        dij(2,1) = cijkl(2,2,1,1)   ! ---
        dij(2,2) = cijkl(2,2,2,2)
        dij(2,3) = cijkl(2,2,1,2)
        dij(3,1) = cijkl(1,2,1,1)   ! ---
        dij(3,2) = cijkl(1,2,2,2)
        dij(3,3) = cijkl(1,2,1,2)
      case( AxisSymetric )
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
      case( Shell )
    end select

  end subroutine mat_c2d

  ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
  !####################################################################
  subroutine MatlMatrix_Shell                        &
      (gauss, sectType, D,                    &
      e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
      alpha, n_layer)
    !####################################################################

    type(tGaussStatus), intent(in)  :: gauss
    integer, intent(in)             :: sectType, n_layer
    real(kind = kreal), intent(out) :: D(:, :)
    real(kind = kreal), intent(in)  :: e1_hat(3), e2_hat(3), e3_hat(3)
    real(kind = kreal), intent(in)  :: cg1(3), cg2(3), cg3(3)
    real(kind = kreal), intent(out) :: alpha

    !--------------------------------------------------------------------

    real(kind = kreal)       :: c(3, 3, 3, 3)
    type(tMaterial), pointer :: matl

    !--------------------------------------------------------------------

    matl => gauss%pMaterial

    !--------------------------------------------------------------------

    if( isElastic(matl%mtype) ) then
      call LinearElastic_Shell                     &
        (matl, sectType, c,                     &
        e1_hat, e2_hat, e3_hat, cg1, cg2, cg3, &
        alpha, n_layer)

      call mat_c2d_Shell(c, D, sectType)
    else
      stop "Material type not supported!"
    end if

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine MatlMatrix_Shell
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)


  ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
  !####################################################################
  subroutine mat_c2d_Shell(c, D, itype)
    !####################################################################

    real(kind = kreal), intent(in)  :: c(:, :, :, :)
    real(kind = kreal), intent(out) :: D(:, :)
    integer, intent(in)             :: itype

    !--------------------------------------------------------------------

    integer :: index_i(5), index_j(5), &
      index_k(5), index_l(5)
    integer :: i, j, k, l
    integer :: is, js

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

    select case( itype )
      case( Shell )

        do js = 1, 5

          do is = 1, 5

            i = index_i(is)
            j = index_j(is)
            k = index_k(js)
            l = index_l(js)

            D(is, js) = c(i, j, k, l)

          end do

        end do

    end select

    !--------------------------------------------------------------------

    return

    !####################################################################
  end subroutine mat_c2d_Shell
  !####################################################################
  ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)

  subroutine GetConnectorProperty( gauss, ctype, dofid, rparams, iparams, stretch )
    type( tGaussStatus ), intent(in)         :: gauss          !> status of qudrature point
    integer(kind=kint), intent(in)           :: ctype          !> connector type
    integer(kind=kint), intent(in)           :: dofid          !> dofid ( used only in ***_D type connector )
    real(kind=kreal), intent(out)            :: rparams(:)      !> real paramters
    integer(kind=kint), intent(out)          :: iparams(:)      !> integer paramters
    real(kind=kreal), intent(in), optional   :: stretch        !> current stretch

    real(kind=kreal) ina(1), outa(4)
    logical :: ierr
    character(len=DICT_KEY_LENGTH) :: cnkey

    rparams(:) = 0.d0
    iparams(:) = 0

    if( ctype == M_SPRING_DOF ) then ! Dof Spring
      write(cnkey,'(A,I0)') trim(MC_SPRING),dofid
      iparams(1) = gauss%pMaterial%variables_i(M_SPRING_D_NDOFFSET+2*dofid  )
      iparams(2) = gauss%pMaterial%variables_i(M_SPRING_D_NDOFFSET+2*dofid+1)
    else if( ctype == M_SPRING_AXIAL ) then ! Dof Spring
      write(cnkey,'(A,A)') trim(MC_SPRING),'_A'
    else if( ctype == M_DASHPOT_DOF ) then ! Dof Dashpot
      write(cnkey,'(A,I0)') trim(MC_DASHPOT),dofid
      iparams(1) = gauss%pMaterial%variables_i(M_DASHPOT_D_NDOFFSET+2*dofid  )
      iparams(2) = gauss%pMaterial%variables_i(M_DASHPOT_D_NDOFFSET+2*dofid+1)
    else if( ctype == M_DASHPOT_AXIAL ) then ! Dof Dashpot
      write(cnkey,'(A,A)') trim(MC_DASHPOT),'_A'
    else
      stop "CONNECTOR ctype is not defined"
    endif

    if( present(stretch) ) then
      ina(1) = stretch
      call fetch_TableData( cnkey, gauss%pMaterial%dict, outa(1:1), ierr, ina )
    else
      call fetch_TableData( cnkey, gauss%pMaterial%dict, outa(1:1), ierr )
    endif
    rparams(1) = outa(1)

  end subroutine

end module m_MatMatrix
