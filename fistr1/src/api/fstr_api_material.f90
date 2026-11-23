!-------------------------------------------------------------------------------
! Copyright (c) 2026 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
module fstr_api_material
  use iso_c_binding
  use mMaterial
  implicit none

  real(kind=kreal), parameter :: PI = 4.0d0 * datan(1.0d0)

contains

  !> @brief Retrive material handler from the solid analysis derived type
  !! @param[in] solid : A solid analysis handler
  !! @param[in] cid : Index of material property
  !! @return Material property handler of type(c_ptr)
  type(c_ptr) function fstr_api_solid_material(solid,cid) &
    bind(C,name='fstr_api_solid_material')
    use m_fstr, only : fstr_solid
    implicit none
    type(c_ptr), value :: solid
    integer(c_int), value, intent(in) :: cid
    type(fstr_solid), pointer :: fstrSOLID
    call c_f_pointer(cptr=solid, fptr=fstrSOLID)
    fstr_api_solid_material = c_loc(fstrSOLID%materials(cid))
  end function

  !> @brief Retrive material type
  !! @param[in] material : A material property handler
  !! @return Material type
  integer(c_int) function fstr_api_material_type(material) &
    bind(C,name='fstr_api_material_type')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_type = mat%mtype
  end function

  !> @brief Retrive material name
  !! @param[in] material : A material property handler
  !! @param[inout] name : A material name
  !! @param[in] namelen : Length of name
  subroutine fstr_api_material_get_name(material,name,namelen) &
    bind(C,name='fstr_api_material_get_name')
    use hecmw_api_common
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(inout) :: name(*)
    integer(c_int), value, intent(in) :: namelen
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    call f_c_str_copy(mat%name, name, namelen)
  end subroutine

  !> @brief Set material name
  !! @param[in] material : A material property handler
  !! @param[in] name : A material name
  subroutine fstr_api_material_set_name(material,mat_name) &
    bind(C,name='fstr_api_material_set_name')
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw_api_common, only : c_f_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(in) :: mat_name(*)
    type(tMaterial), pointer :: mat
    character(len=HECMW_NAME_LEN) :: mat_name_f

    call c_f_pointer(cptr=material, fptr=mat)
    call c_f_str_copy(mat_name, mat_name_f, HECMW_NAME_LEN)
    mat%name = mat_name_f
  end subroutine

  !> @brief Is material property elastic?
  !! @param[in] material : A material property handler
  !! @return A logical value indicating whether it is elastic
  logical(c_bool) function fstr_api_material_isElastic(material) &
    bind(C,name='fstr_api_material_isElastic')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_isElastic = isElastic(mat%mtype)
  end function

  !> @brief Is material property elastoplastic?
  !! @param[in] material : A material property handler
  !! @return A logical value indicating whether it is elastoplastic
  logical(c_bool) function fstr_api_material_isElastoplastic(material) &
    bind(C,name='fstr_api_material_isElastoplastic')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_isElastoplastic = isElastoplastic(mat%mtype)
  end function

  !> @brief Is material property hyperelastic?
  !! @param[in] material : A material property handler
  !! @return A logical value indicating whether it is hyperelastic
  logical(c_bool) function fstr_api_material_isHyperelastic(material) &
    bind(C,name='fstr_api_material_isHyperelastic')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_isHyperelastic = isHyperelastic(mat%mtype)
  end function

  !> @brief Is material property viscoelastic?
  !! @param[in] material : A material property handler
  !! @return A logical value indicating whether it is viscoelastic
  logical(c_bool) function fstr_api_material_isViscoelastic(material) &
    bind(C,name='fstr_api_material_isViscoelastic')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_isViscoelastic = isViscoelastic(mat%mtype)
  end function

  !> @brief Is material property creep?
  !! @param[in] material : A material property handler
  !! @return A logical value indicating whether it is creep
  logical(c_bool) function fstr_api_material_isCreep(material) &
    bind(C,name='fstr_api_material_isCreep')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_isCreep = isCreep(mat%mtype)
  end function

  !> @brief Retirive the type of geometric nonlinearity
  !! @param[in] material : A material property handler
  !! @return the type of geometric nonlinearity
  integer(c_int) function fstr_api_material_nlgeom(material) &
    bind(C,name='fstr_api_material_nlgeom')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_nlgeom = mat%nlgeom_flag
  end function

  !> @brief Set the type of geometric nonlinearity
  !! @param[in] material : A material property handler
  !! @param[in] nlgeom : The type of geometric nonlinearity
  subroutine fstr_api_material_set_nlgeom(material,nlgeom) &
    bind(C,name='fstr_api_material_set_nlgeom')
    implicit none
    type(c_ptr), value :: material
    integer(c_int), value :: nlgeom
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    mat%nlgeom_flag = nlgeom
  end subroutine

  !> @brief Retirive elastic parameters
  !! @param[in] material : A material property handler
  !! @param[out] youngs : Young's modulus
  !! @param[out] poisson : Poisson's ratio
  subroutine fstr_api_material_get_elastic(material,youngs,poisson) &
    bind(C,name='fstr_api_material_get_elastic')
    implicit none
    type(c_ptr), value :: material
    real(c_double), intent(out) :: youngs
    real(c_double), intent(out) :: poisson

    type(tMaterial), pointer :: mat
    real(kind=kreal) :: ina(1), outa(2)
    logical :: ierr

    call c_f_pointer(cptr=material, fptr=mat)
 
    ina(1) = 0.0
    call fetch_TableData( MC_ISOELASTIC, mat%dict, outa, ierr, ina )
    if( ierr ) then
      youngs = mat%variables(M_YOUNGS)
      poisson = mat%variables(M_POISSON)
    else
      youngs = outa(1)
      poisson = outa(2)
    end if
  end subroutine

  !> @brief Set elastic parameters
  !! @param[in] material : A material property handler
  !! @param[in] youngs : Young's modulus
  !! @param[in] poisson : Poisson's ratio
  subroutine fstr_api_material_set_elastic(material,youngs,poisson) &
    bind(C,name='fstr_api_material_set_elastic')
    implicit none
    type(c_ptr), value :: material
    real(c_double), value :: youngs
    real(c_double), value :: poisson
    type(tMaterial), pointer :: mat
    type( tTable ) :: mattable
    real(kind=kreal), allocatable :: fval(:,:)

    call c_f_pointer(cptr=material, fptr=mat)
    mat%mtype = ELASTIC

    mat%variables(M_YOUNGS) = youngs
    mat%variables(M_POISSON) = poisson

    allocate(fval(2,1))
    fval(1,1) = youngs 
    fval(2,1) = poisson
    call init_table( mattable, 0, 2, 1, fval )
    call dict_add_key( mat%dict, MC_ISOELASTIC, mattable )
    call finalize_table( mattable )
    deallocate(fval)
  end subroutine

  !> @brief Retirive density
  !! @param[in] material : A material property handler
  !! @return Density
  real(c_double) function fstr_api_material_density(material) &
    bind(C,name='fstr_api_material_density')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_density = mat%variables(M_DENSITY)
  end function

  !> @brief Set density
  !! @param[in] material : A material property handler
  !! @param[in] density : Density
  subroutine fstr_api_material_set_density(material,density) &
    bind(C,name='fstr_api_material_set_density')
    implicit none
    type(c_ptr), value :: material
    real(c_double), value :: density
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    mat%variables(M_DENSITY) = density
  end subroutine

  !> @brief Retirive coefficient of linear expansion
  !! @param[in] material : A material property handler
  !! @return Coefficient of linear expansion
  real(c_double) function fstr_api_material_expansion(material) &
    bind(C,name='fstr_api_material_expansion')
    implicit none
    type(c_ptr), value :: material
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    fstr_api_material_expansion = mat%variables(M_EXPANSION)
  end function

  !> @brief Set coefficient of linear expansion
  !! @param[in] material : A material property handler
  !! @return Coefficient of linear expansion
  subroutine fstr_api_material_set_expansion(material,expansion) &
    bind(C,name='fstr_api_material_set_expansion')
    implicit none
    type(c_ptr), value :: material
    real(c_double), value :: expansion
    type(tMaterial), pointer :: mat
    call c_f_pointer(cptr=material, fptr=mat)
    mat%variables(M_EXPANSION) = expansion
  end subroutine

  !> @brief Retirive yield criterion, hardening rule and size of plastic parameters
  !! @param[in] material : A material property handler
  !! @param[out] yield : Yield criterion
  !! @param[out] harden : Hardening rule
  !! @param[out] nvar : Size of parameters
  !! @param[out] rows : Number of rows in the table
  !! @param[out] cols : Number of cols in the table
  subroutine fstr_api_material_get_plastic_size(material,yield,harden,nvar,rows,cols) &
    bind(C,name='fstr_api_material_get_plastic_size')
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(inout) :: yield(*)
    character(kind=c_char), intent(inout) :: harden(*)
    integer(c_int), intent(out) :: nvar
    integer(c_int), intent(out) :: rows
    integer(c_int), intent(out) :: cols

    type(tMaterial), pointer :: mat
    type(DICT_DATA), pointer :: dicval

    call c_f_pointer(cptr=material, fptr=mat)
    if(mat%mtype==120000) then
      call f_c_str_copy("MISES",yield,32)
      call f_c_str_copy("BILINEAR",harden,32)
      nvar = 2
      rows = 0
      cols = 0
    else if(mat%mtype==120010) then
      call f_c_str_copy("MISES",yield,32)
      call f_c_str_copy("MULTILINEAR",harden,32)
      dicval => dict_get_key( mat%dict, MC_YIELD )
      nvar = 0
      rows = dicval%tbrow
      cols = dicval%tbcol
    else if(mat%mtype==120020) then
      call f_c_str_copy("MISES",yield,32)
      call f_c_str_copy("SWIFT",harden,32)
      nvar = 3
      rows = 0
      cols = 0
    else if(mat%mtype==120030) then
      call f_c_str_copy("MISES",yield,32)
      call f_c_str_copy("RAMBERG-OSGOOD",harden,32)
      nvar = 3
      rows = 0
      cols = 0
    else if(mat%mtype==120040) then
      call f_c_str_copy("MISES",yield,32)
      call f_c_str_copy("KINEMATIC",harden,32)
      nvar = 2
      rows = 0
      cols = 0
    else if(mat%mtype==120050) then
      call f_c_str_copy("MISES",yield,32)
      call f_c_str_copy("COMBINED",harden,32)
      nvar = 3
      rows = 0
      cols = 0
    else if(mat%mtype==120100) then
      call f_c_str_copy("MOHR-COULOMB",yield,32)
      call f_c_str_copy("BILINEAR",harden,32)
      nvar = 4
      rows = 0
      cols = 0
    else if(mat%mtype==120110) then
      call f_c_str_copy("MOHR-COULOMB",yield,32)
      call f_c_str_copy("MULTILINEAR",harden,32)
      dicval => dict_get_key( mat%dict, MC_YIELD )
      nvar = 2
      rows = dicval%tbrow
      cols = dicval%tbcol
    else if(mat%mtype==120200) then
      call f_c_str_copy("DRUCKER-PRAGER",yield,32)
      call f_c_str_copy("BILINEAR",harden,32)
      nvar = 4
      rows = 0
      cols = 0
    else if(mat%mtype==120210) then
      call f_c_str_copy("DRUCKER-PRAGER",yield,32)
      call f_c_str_copy("MULTILINEAR",harden,32)
      dicval => dict_get_key( mat%dict, MC_YIELD )
      nvar = 2
      rows = dicval%tbrow
      cols = dicval%tbcol
    else if(mat%mtype==120300) then
      call f_c_str_copy("USER",yield,32)
      call f_c_str_copy("BILINEAR",harden,32)
      nvar = 100
      rows = 0
      cols = 0
    else
      call f_c_str_copy("NOT-DEFINED",yield,32)
      call f_c_str_copy("NOT-DEFINED",harden,32)
      nvar = 0
      rows = 0
      cols = 0
    end if
  end subroutine

  !> @brief Retirive parameters of plastic material
  !! @param[in] material : A material property handler
  !! @param[out] variables : Parameters (nvar + rows*cols)
  subroutine fstr_api_material_get_plastic_vars(material,variables) &
    bind(C,name='fstr_api_material_get_plastic_vars')
    implicit none
    type(c_ptr), value :: material
    real(c_double), intent(out) :: variables(*)
    type(tMaterial), pointer :: mat
    type(DICT_DATA), pointer :: dicval
    real(c_double) :: eta
    real(c_double) :: etabar
    integer :: i, j

    call c_f_pointer(cptr=material, fptr=mat)
    if(mat%mtype==120000) then
      ! MISES BILINEAR
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
    else if(mat%mtype==120010) then
      ! MISES MULTILINEAR
      dicval => dict_get_key( mat%dict, MC_YIELD )
      do j=1, dicval%tbrow
        do i=1, dicval%tbcol
          variables(i+(j-1)*dicval%tbcol) = dicval%tbval(i,j)
        end do
      end do
    else if(mat%mtype==120020) then
      ! MISES SWIFT
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
      variables(3) = mat%variables(M_PLCONST3)
    else if(mat%mtype==120030) then
      ! MISES RAMBERG-OSGOOD
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
      variables(3) = mat%variables(M_PLCONST3)
    else if(mat%mtype==120040) then
      ! MISES KINEMATIC
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST3)
    else if(mat%mtype==120050) then
      ! MISES COMBINED
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
      variables(3) = mat%variables(M_PLCONST3)
    else if(mat%mtype==120100) then
      ! MOHR-COULOMB BILINEAR
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST3)*180.0/PI
      variables(3) = mat%variables(M_PLCONST2)
      variables(4) = mat%variables(M_PLCONST4)*180.0/PI
    else if(mat%mtype==120110) then
      ! MOHR-COULOMB MULTILINEAR
      variables(1) = mat%variables(M_PLCONST3)*180.0/PI
      variables(2) = mat%variables(M_PLCONST4)*180.0/PI
      dicval => dict_get_key( mat%dict, MC_YIELD )
      do j=1, dicval%tbrow
        do i=1, dicval%tbcol
          variables(2+i+(j-1)*dicval%tbcol) = dicval%tbval(i,j)
        end do
      end do
    else if(mat%mtype==120200) then
      ! DRUCKER-PRAGER BILINEAR
      eta = mat%variables(M_PLCONST3) ! 2.d0*sin(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) )
      etabar = mat%variables(M_PLCONST5) ! 2.d0*sin(psi)/ ( sqrt(3.d0)*(3.d0+sin(psi)) )
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = asin(3.0*sqrt(3.d0)*eta/(2.0-sqrt(3.d0)*eta))*180.0/PI
      variables(3) = mat%variables(M_PLCONST2)
      variables(4) = asin(3.0*sqrt(3.d0)*etabar/(2.0-sqrt(3.d0)*etabar))*180.0/PI
    else if(mat%mtype==120210) then
      ! DRUCKER-PRAGER MULTILINEAR
      eta = mat%variables(M_PLCONST3) ! 2.d0*sin(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! eta
      etabar = mat%variables(M_PLCONST5) ! 2.d0*sin(psi)/ ( sqrt(3.d0)*(3.d0+sin(psi)) ) ! etabar
      variables(1) = asin(3.0*sqrt(3.d0)*eta/(2.0-sqrt(3.d0)*eta))*180.0/PI
      variables(2) = asin(3.0*sqrt(3.d0)*etabar/(2.0-sqrt(3.d0)*etabar))*180.0/PI
      dicval => dict_get_key( mat%dict, MC_YIELD )
      do j=1, dicval%tbrow
        do i=1, dicval%tbcol
          variables(2+i+(j-1)*dicval%tbcol) = dicval%tbval(i,j)
        end do
      end do
    else if(mat%mtype==120300) then
      ! USER
      do i=1, 100
        variables(i) = mat%variables(100+i)
      end do
    else
      ! NOT-DEFINED
    end if
  end subroutine

  !> @brief Set yield criterion, hardening rule and parameters
  !! @param[in] material : A material property handler
  !! @param[in] yield : Yield criterion
  !! @param[in] harden : Hardening rule
  !! @param[in] variables : Parameters
  !! @param[in] rows : Number of rows in the table
  !! @param[in] cols : Number of cols in the table
  subroutine fstr_api_material_set_plastic(material,yield,harden,variables,rows,cols) &
    bind(C,name='fstr_api_material_set_plastic')
    use hecmw_api_common, only : c_f_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(in) :: yield(*)
    character(kind=c_char), intent(in) :: harden(*)
    real(c_double), intent(in)         :: variables(*)
    integer(c_int), value :: rows
    integer(c_int), value :: cols

    type(tMaterial), pointer :: mat
    character(len=HECMW_NAME_LEN) :: yield_f
    character(len=HECMW_NAME_LEN) :: harden_f
    real(kind=kreal) :: phi, psi
    integer :: i, j
    type(tTable) :: mattable
    real(kind=kreal) :: fval(cols,rows)

    call c_f_pointer(cptr=material, fptr=mat)
    call c_f_str_copy(yield, yield_f, HECMW_NAME_LEN)
    call c_f_str_copy(harden, harden_f, HECMW_NAME_LEN)

    if(yield_f=="MISES") then
      if(harden_f=="BILINEAR") then
        mat%mtype = EPLASTIC
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = variables(2)
      else if(harden_f=="MULTILINEAR") then
        mat%mtype = EPLASTIC + 10
        fval = reshape(variables(1:cols*rows), [cols,rows], order=[1,2])
        call init_table( mattable, cols-2, cols, rows, fval )
        call dict_add_key( mat%dict, MC_YIELD, mattable )
        call finalize_table( mattable )
      else if(harden_f=="SWIFT") then
        mat%mtype = EPLASTIC + 20
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = variables(2)
        mat%variables(M_PLCONST3) = variables(3)
      else if(harden_f=="RAMBERRG-OSWOOD") then
        mat%mtype = EPLASTIC + 30
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = variables(2)
        mat%variables(M_PLCONST3) = variables(3)
      else if(harden_f=="KINEMATIC") then
        mat%mtype = EPLASTIC + 40
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = 0.0
        mat%variables(M_PLCONST3) = variables(2)
      else if(harden_f=="COMBINED") then
        mat%mtype = EPLASTIC + 50
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = variables(2)
        mat%variables(M_PLCONST3) = variables(3)
      endif
    else if(yield_f=="MOHR-COULOMB") then
      if(harden_f=="BILINEAR") then
        mat%mtype = EPLASTIC + 100
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = variables(3)
        phi = variables(2)*PI/180.d0
        if( variables(4) >= 0.d0 ) then
          psi = variables(4)*PI/180.d0
        else
          psi = phi
        endif
        mat%variables(M_PLCONST3) = phi
        mat%variables(M_PLCONST4) = psi
      else if(harden_f=="MULTILINEAR") then
        mat%mtype = EPLASTIC + 110
        phi = variables(1)*PI/180.d0
        if( variables(1) >= 0.d0 ) then
          psi = variables(2)*PI/180.d0
        else
          psi = phi
        endif
        mat%variables(M_PLCONST3) = phi
        mat%variables(M_PLCONST4) = psi
        fval = reshape(variables(3:cols*rows+2), [cols,rows], order=[1,2])
        call init_table( mattable, 1, cols, rows, fval )
        call dict_add_key( mat%dict, MC_YIELD, mattable )
        call finalize_table( mattable )
      endif
    else if(yield_f=="DRUCKER-PRAGER") then
      if(harden_f=="BILINEAR") then
        mat%mtype = EPLASTIC + 200
        mat%variables(M_PLCONST1) = variables(1)
        mat%variables(M_PLCONST2) = variables(3)
        phi = variables(2)*PI/180.d0
        if( variables(4) >= 0.d0 ) then
          psi = variables(4)*PI/180.d0
        else
          psi = phi
        endif
        mat%variables(M_PLCONST3) = 2.d0*sin(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! eta
        mat%variables(M_PLCONST4) = 6.d0*cos(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! xi
        mat%variables(M_PLCONST5) = 2.d0*sin(psi)/ ( sqrt(3.d0)*(3.d0+sin(psi)) ) ! etabar
      else if(harden_f=="MULTILINEAR") then
        mat%mtype = EPLASTIC + 210
        phi = variables(1)*PI/180.d0
        if( variables(2) >= 0.d0 ) then
          psi = variables(2)*PI/180.d0
        else
          psi = phi
        endif
        mat%variables(M_PLCONST3) = 2.d0*sin(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! eta
        mat%variables(M_PLCONST4) = 6.d0*cos(phi)/ ( sqrt(3.d0)*(3.d0+sin(phi)) ) ! xi
        mat%variables(M_PLCONST5) = 2.d0*sin(psi)/ ( sqrt(3.d0)*(3.d0+sin(psi)) ) ! etabar
        fval = reshape(variables(3:cols*rows+2), [cols,rows], order=[1,2])
        call init_table( mattable, 1, cols, rows, fval )
        call dict_add_key( mat%dict, MC_YIELD, mattable )
        call finalize_table( mattable )
      endif
    else if(yield_f=="USER") then
      mat%mtype = EPLASTIC + 300
      do i=1, 100
        mat%variables(100+i) = variables(i)
      end do
    endif

  end subroutine

  !> @brief Retirive type and parameters of hyperelastic material
  !! @param[in] material : A material property handler
  !! @param[out] type : Hyperelastic type
  !! @param[out] variables : Parameters
  subroutine fstr_api_material_get_hyperelastic(material,type,variables) &
    bind(C,name='fstr_api_material_get_hyperelastic')
    use mMaterial
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(inout) :: type(*)
    real(c_double), intent(out) :: variables(*)
    type(tMaterial), pointer :: mat
    integer :: i

    call c_f_pointer(cptr=material, fptr=mat)
    if(mat%mtype==NEOHOOKE) then
      call f_c_str_copy("NEOHOOKE",type,32)
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST3)
    else if(mat%mtype==MOONEYRIVLIN) then
      call f_c_str_copy("MOONEY-RIVLIN",type,32)
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
      variables(3) = mat%variables(M_PLCONST3)
    else if(mat%mtype==ARRUDABOYCE) then
      call f_c_str_copy("ARRUDA-BOYCE",type,32)
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
      variables(3) = mat%variables(M_PLCONST3)
    else if(mat%mtype==USERHYPERELASTIC) then
      call f_c_str_copy("USER",type,32)
      do i=1, 100
        variables(i) = mat%variables(100+i)
      end do
    else if(mat%mtype==MOONEYRIVLIN_ANISO) then
      call f_c_str_copy("MOONEY-RIVLIN-ANISO",type,32)
      variables(1) = mat%variables(M_PLCONST1)
      variables(2) = mat%variables(M_PLCONST2)
      variables(3) = mat%variables(M_PLCONST3)
      variables(4) = mat%variables(M_PLCONST4)
      variables(5) = mat%variables(M_PLCONST5)
    endif
  end subroutine

  !> @brief Set type and parameters of hyperelastic material
  !! @param[in] material : A material property handler
  !! @param[in] type : Hyperelastic type
  !! @param[in] variables : Parameters
  subroutine fstr_api_material_set_hyperelastic(material,type,variables) &
    bind(C,name='fstr_api_material_set_hyperelastic')
    use mMaterial
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw_api_common, only : c_f_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(in) :: type(*)
    real(c_double), intent(in) :: variables(*)
    type(tMaterial), pointer :: mat

    character(len=HECMW_NAME_LEN) :: type_f
    integer :: i

    call c_f_pointer(cptr=material, fptr=mat)
    call c_f_str_copy(type, type_f, HECMW_NAME_LEN)

    if(type_f=="NEOHOOKE") then
      mat%mtype = NEOHOOKE
      mat%nlgeom_flag = TOTALLAG
      mat%nfstatus = 0
      mat%variables(M_PLCONST1) = variables(1)
      mat%variables(M_PLCONST2) = 0.d0
      mat%variables(M_PLCONST3) = variables(2)
    else if(type_f=="MOONEY-RIVLIN") then
      mat%mtype = MOONEYRIVLIN
      mat%nlgeom_flag = TOTALLAG
      mat%nfstatus = 0
      mat%variables(M_PLCONST1) = variables(1)
      mat%variables(M_PLCONST2) = variables(2)
      mat%variables(M_PLCONST3) = variables(3)
    else if(type_f=="ARRUDA-BOYCE") then
      mat%mtype=ARRUDABOYCE
      mat%variables(M_PLCONST1) = variables(1)
      mat%variables(M_PLCONST2) = variables(2)
      mat%variables(M_PLCONST3) = variables(3)
    else if(type_f=="USER") then
      mat%mtype=USERHYPERELASTIC
      do i=1, 100
        mat%variables(100+i) = variables(i)
      end do
    else if(type_f=="MOONEY-RIVLIN-ANISO") then
      mat%mtype=MOONEYRIVLIN_ANISO
      mat%variables(M_PLCONST1) = variables(1)
      mat%variables(M_PLCONST2) = variables(2)
      mat%variables(M_PLCONST3) = variables(3)
      mat%variables(M_PLCONST4) = variables(4)
      mat%variables(M_PLCONST5) = variables(5)
    endif
  end subroutine

  !> @brief Retirive size of viscoelastic parameter table
  !! @param[in] material : A material property handler
  !! @param[out] rows : Number of rows in the table
  !! @param[out] cols : Number of cols in the table
  subroutine fstr_api_material_get_viscoelastic_size(material,rows,cols) &
    bind(C,name='fstr_api_material_get_viscoelastic_size')
    use mMaterial
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: material
    integer(c_int), intent(out) :: rows
    integer(c_int), intent(out) :: cols

    type(tMaterial), pointer :: mat
    type(DICT_DATA), pointer :: dicval

    call c_f_pointer(cptr=material, fptr=mat)
    dicval => dict_get_key( mat%dict, MC_VISCOELASTIC )
    rows = dicval%tbrow
    cols = dicval%tbcol
  end subroutine

  !> @brief Retrive viscoelastic parameters
  !! @param[in] material : A material property handler
  !! @param[out] variables : Parameters (rows*cols)
  subroutine fstr_api_material_get_viscoelastic_vars(material,variables) &
    bind(C,name='fstr_api_material_get_viscoelastic_vars')
    use mMaterial
    implicit none
    type(c_ptr), value :: material
    real(c_double), intent(out) :: variables(*)

    type(tMaterial), pointer :: mat
    type(DICT_DATA), pointer :: dicval
    integer :: i, j

    call c_f_pointer(cptr=material, fptr=mat)

    dicval => dict_get_key( mat%dict, MC_VISCOELASTIC )
    do j=1, dicval%tbrow
      do i=1, dicval%tbcol
        variables(i+(j-1)*dicval%tbcol) = dicval%tbval(i,j)
      end do
    end do
  end subroutine

  !> @brief Set viscoelastic parameter table
  !! @param[in] material : A material property handler
  !! @param[in] variables : Parameters (rows*cols)
  !! @param[in] rows : Number of rows in the table
  !! @param[in] cols : Number of cols in the table
  subroutine fstr_api_material_set_viscoelastic(material,variables,rows,cols) &
    bind(C,name='fstr_api_material_set_viscoelastic')
    use mMaterial
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw_api_common, only : c_f_str_copy
    implicit none
    type(c_ptr), value :: material
    real(c_double), intent(in) :: variables(*)
    integer(c_int), value :: rows
    integer(c_int), value :: cols

    integer(kind=kint) :: depends
    type(tMaterial), pointer :: mat
    type(tTable) :: mattable
    real(kind=kreal) :: fval(cols,rows)

    call c_f_pointer(cptr=material, fptr=mat)
    fval = reshape(variables(1:cols*rows), [cols,rows], order=[1,2])

    mat%nlgeom_flag = TOTALLAG
    depends = 1
    call init_table( mattable, depends, cols, rows, fval )
    call dict_add_key( mat%dict, MC_VISCOELASTIC, mattable )
    mat%mtype = VISCOELASTIC
    call finalize_table( mattable )

  end subroutine

  !> @brief Retirive type and size of creep parameter table
  !! @param[in] material : A material property handler
  !! @param[out] type : Creep type
  !! @param[out] rows : Number of rows in the table
  !! @param[out] cols : Number of cols in the table
  subroutine fstr_api_material_get_creep_size(material,type,rows,cols) &
    bind(C,name='fstr_api_material_get_creep_size')
    use mMaterial
    use hecmw_api_common, only : f_c_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(inout) :: type(*)
    integer(c_int), intent(out) :: rows
    integer(c_int), intent(out) :: cols

    type(tMaterial), pointer :: mat
    type(DICT_DATA), pointer :: dicval

    call c_f_pointer(cptr=material, fptr=mat)
    rows = 0
    cols = 0
    if(mat%mtype==NORTON) then
      call f_c_str_copy("NORTON",type,32)
      dicval => dict_get_key( mat%dict, MC_NORTON )
      rows = dicval%tbrow
      cols = dicval%tbcol
    endif
  end subroutine

  !> @brief Retrive creep parameters
  !! @param[in] material : A material property handler
  !! @param[out] variables : Parameters (rows*cols)
  subroutine fstr_api_material_get_creep_vars(material,variables) &
    bind(C,name='fstr_api_material_get_creep_vars')
    use mMaterial
    implicit none
    type(c_ptr), value :: material
    real(c_double), intent(out) :: variables(*)

    type(tMaterial), pointer :: mat
    type(DICT_DATA), pointer :: dicval
    integer :: i, j

    call c_f_pointer(cptr=material, fptr=mat)

    if(mat%mtype==NORTON) then
      dicval => dict_get_key( mat%dict, MC_NORTON )
      do j=1, dicval%tbrow
        do i=1, dicval%tbcol
          variables(i+(j-1)*dicval%tbcol) = dicval%tbval(i,j)
        end do
      end do
    end if
  end subroutine

  !> @brief Set creep parameter table
  !! @param[in] material : A material property handler
  !! @param[in] type : Creep type
  !! @param[in] variables : Parameters (rows*cols)
  !! @param[in] rows : Number of rows in the table
  !! @param[in] cols : Number of cols in the table
  subroutine fstr_api_material_set_creep(material,type,variables,rows,cols) &
    bind(C,name='fstr_api_material_set_creep')
    use mMaterial
    use hecmw_util, only : HECMW_NAME_LEN
    use hecmw_api_common, only : c_f_str_copy
    implicit none
    type(c_ptr), value :: material
    character(kind=c_char), intent(in) :: type(*)
    real(c_double), intent(in) :: variables(*)
    integer(c_int), value :: rows
    integer(c_int), value :: cols

    integer(kind=kint) :: depends
    type(tMaterial), pointer :: mat
    type(tTable) :: mattable
    real(kind=kreal) :: fval(cols,rows)

    character(len=HECMW_NAME_LEN) :: type_f

    call c_f_pointer(cptr=material, fptr=mat)
    call c_f_str_copy(type, type_f, HECMW_NAME_LEN)
    fval = reshape(variables(1:cols*rows), [cols,rows], order=[1,2])

    if(type_f=="NORTON") then
      mat%nlgeom_flag = UPDATELAG
      depends = cols-3
      call init_table( mattable, depends, cols, rows, fval )
      call dict_add_key( mat%dict, MC_NORTON, mattable )
      mat%mtype = NORTON
      call finalize_table( mattable )
    end if

  end subroutine

end module
