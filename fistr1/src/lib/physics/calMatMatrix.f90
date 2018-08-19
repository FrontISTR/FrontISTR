!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
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
  subroutine MatlMatrix( gauss, sectType, matrix, time, dtime, cdsys, temperature, isEp )
    type( tGaussStatus ), intent(in) :: gauss          !> status of qudrature point
    integer, intent(in)              :: sectType       !> plane strain/stress or 3D
    real(kind=kreal), intent(out)    :: matrix(:,:)    !> constitutive matrix
    real(kind=kreal), intent(in)     :: time           !> current time
    real(kind=kreal), intent(in)     :: dtime          !> time increment
    real(kind=kreal), intent(in)     :: cdsys(3,3)     !> material coordinate system
    real(kind=kreal), intent(in), optional  :: temperature   !> temperature
    integer(kind=kint), intent(in), optional :: isEp

    integer :: i
    integer :: flag
    real(kind=kreal)            :: cijkl(3,3,3,3)
    type( tMaterial ), pointer  :: matl
    matl=>gauss%pMaterial

    flag = 0
    if( present(isEp) )then
      if( isEp == 1 )flag = 1
    endif

    if( matl%mtype==USERELASTIC ) then
      call uElasticMatrix( matl%variables(101:), gauss%strain, matrix )
    elseif( isViscoelastic(matl%mtype) ) then
      if( present(temperature) ) then
        call calViscoelasticMatrix( matl, sectTYPE, dtime, matrix, temperature )
      else
        call calViscoelasticMatrix( matl, sectTYPE, dtime, matrix )
      endif
    elseif( isElastic(matl%mtype) .or. flag==1 ) then
      if(flag==1)then
        i = getElasticType(ELASTIC)
      else
        i = getElasticType(gauss%pMaterial%mtype)
      endif

      if( i==0 ) then
        if( present(temperature) ) then
          call calElasticMatrix( matl, sectTYPE, matrix, temperature  )
        else
          call calElasticMatrix( matl, sectTYPE, matrix )
        endif
      elseif(  i==1 ) then
        if( present(temperature) ) then
          call calElasticMatrix_ortho( gauss%pMaterial, sectTYPE, cdsys, matrix, temperature )
        else
          call calElasticMatrix_ortho( gauss%pMaterial, sectTYPE, cdsys, matrix )
        endif
      else
        print *, "Elasticity type", matl%mtype, "not supported"
        stop
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
        gauss%stress, gauss%fstatus, matrix, dtime, time )
    elseif( matl%mtype==NORTON ) then
      if( present( temperature ) ) then
        call iso_creep( matl, sectTYPE, gauss%stress, gauss%strain, gauss%fstatus,  &
          gauss%plstrain, dtime, time, matrix, temperature  )
      else
        call iso_creep( matl, sectTYPE, gauss%stress, gauss%strain, gauss%fstatus,  &
          gauss%plstrain, dtime, time, matrix  )
      endif
    else
      stop "Material type not supported!"
    endif

  end subroutine

  !
  !> Update strain and stress for elastic and hyperelastic materials
  subroutine StressUpdate( gauss, sectType, strain, stress, time, dtime, temp, tempn )
    type( tGaussStatus ), intent(inout) :: gauss      !> status of qudrature point
    integer, intent(in)                 :: sectType   !> plane strain/stress or 3D
    real(kind=kreal), intent(in)        :: strain(6)  !> strain
    real(kind=kreal), intent(out)       :: stress(6)  !> stress
    real(kind=kreal), intent(in), optional  :: time   !> current time
    real(kind=kreal), intent(in), optional  :: dtime  !> time increment
    real(kind=kreal), optional          :: temp       !> current temprature
    real(kind=kreal), optional          :: tempn      !> temperature at last step

    if( gauss%pMaterial%mtype==NEOHOOKE .or. gauss%pMaterial%mtype==MOONEYRIVLIN ) then
      call calUpdateElasticMooneyRivlin( gauss%pMaterial, sectType, strain, stress )
    elseif( gauss%pMaterial%mtype==ARRUDABOYCE ) then ! Arruda-Boyce Hyperelastic material
      call calUpdateElasticArrudaBoyce( gauss%pMaterial, sectType, strain, stress )
    elseif( gauss%pMaterial%mtype==USERHYPERELASTIC .or. gauss%pMaterial%mtype==USERELASTIC ) then ! user-defined
      call uElasticUpdate( gauss%pMaterial%variables(101:), strain, stress )
    elseif( isViscoelastic( gauss%pMaterial%mtype) ) then
      if( .not. present(dtime) ) stop "error in viscoelastic update!"
      if( present(temp) .and. present(tempn) ) then
        call UpdateViscoelastic( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, dtime, temp, tempn )
      else
        call UpdateViscoelastic( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, dtime )
      endif
    elseif ( gauss%pMaterial%mtype==NORTON ) then
      if( .not. present(dtime)  ) stop "error in viscoelastic update!"
      if( present(temp) ) then
        call update_iso_creep( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, gauss%plstrain, dtime, time, temp )
      else
        call update_iso_creep( gauss%pMaterial, sectType, strain, stress, gauss%fstatus, gauss%plstrain, dtime, time )
      endif
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

end module m_MatMatrix
