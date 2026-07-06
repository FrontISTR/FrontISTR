!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module assembles the tangent stiffness matrix and, in the
!> implicit dynamic case, also the effective dynamic system matrix
!> A = c1*K + c2*M (+ b3*C for connectors) together with the corresponding
!> RHS contribution stored in fstrSOLID%DFORCE.

module m_fstr_CreateMatrix_and_DampingForce
  use m_fstr
  implicit none

  private
  public :: fstr_CreateMatrix_and_DampingForce

contains

  !---------------------------------------------------------------------*
  !> \brief Assemble the system matrix and, optionally, the dynamic damping force.
  !!
  !! Static / eigen / frequency (fstrDYNAMIC and coef absent):
  !!   Assemble only the tangent stiffness matrix K into hecMAT.
  !!
  !! Implicit dynamic (fstrDYNAMIC and coef present):
  !!   Assemble the effective dynamic system matrix
  !!     A = c1*K + c2*M  (+ b3*C  for connector elements)
  !!   into hecMAT, and accumulate the mass / Rayleigh damping / connector
  !!   damping contribution to the RHS into fstrSOLID%DFORCE.
  subroutine fstr_CreateMatrix_and_DampingForce( hecMESH, hecMAT, fstrSOLID, time, tincr, fstrDYNAMIC, coef )
  !---------------------------------------------------------------------*
    use m_static_LIB
    use mMechGauss
    use m_dynamic_mass
    use m_elemact
    use m_fstr_NodalKinematics, only: fstr_ensure_finite_rotation_state

    type(hecmwST_local_mesh)              :: hecMESH      !< mesh information
    type(hecmwST_matrix)                  :: hecMAT       !< system matrix
    type(fstr_solid)                      :: fstrSOLID    !< solid state
    real(kind=kreal), intent(in)          :: time         !< current time
    real(kind=kreal), intent(in)          :: tincr        !< time increment
    type(fstr_dynamic), intent(in), optional :: fstrDYNAMIC !< dynamic info (omit for static)
    real(kind=kreal), intent(in), optional :: coef(6)     !< Newmark time-integration coefficients

    integer(kind=kint) :: ndof, itype, iS, iE, ic_type, nn, icel, iiS, i, j, in, jn
    integer(kind=kint) :: nodLOCAL(fstrSOLID%max_ncon)
    real(kind=kreal)   :: stiff_mat(fstrSOLID%max_ncon_stf*6, fstrSOLID%max_ncon_stf*6)
    real(kind=kreal)   :: mass_mat(20*6, 20*6), damp_mat(20*6, 20*6)
    real(kind=kreal)   :: mat(20*6, 20*6), lumped(20*6)
    real(kind=kreal)   :: tt(fstrSOLID%max_ncon), ecoord(3,fstrSOLID%max_ncon)
    real(kind=kreal)   :: u(6,fstrSOLID%max_ncon), du(6*20), u_prev(6,fstrSOLID%max_ncon)
    real(kind=kreal)   :: vecA(20*6), vecB(20*6)
    real(kind=kreal)   :: acc(20*6), vec(20*6), df(20*6)
    real(kind=kreal)   :: a1, a2, a3, b1, b2, b3
    type(tMaterial), pointer :: material
    logical            :: is_dynamic

    is_dynamic = present(fstrDYNAMIC) .and. present(coef)

    call hecmw_mat_clear( hecMAT )
    if( is_dynamic ) then
      fstrSOLID%DFORCE = 0.0d0
      a1 = coef(1); a2 = coef(2); a3 = coef(3)
      b1 = coef(4); b2 = coef(5); b3 = coef(6)
    endif

    ndof = hecMAT%NDOF
    call fstr_ensure_finite_rotation_state( hecMESH, fstrSOLID, ndof )

    do itype = 1, hecMESH%n_elem_type
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      ic_type = hecMESH%elem_type_item(itype)

      ! ----- Ignore link and patch elements
      if (hecmw_is_etype_link(ic_type)) cycle
      if (hecmw_is_etype_patch(ic_type)) cycle

      nn = hecmw_get_max_node(ic_type)

      !$omp parallel default(none), &
        !$omp&  private(icel,iiS,nn,j,nodLOCAL,i,in,jn,ecoord,du,u,u_prev,tt, &
        !$omp&          material,stiff_mat,mass_mat,damp_mat, &
        !$omp&          lumped,mat,df,vec,acc,vecA,vecB), &
        !$omp&  shared(iS,iE,hecMESH,ndof,fstrSOLID,ic_type,hecMAT,time,tincr,fstrDYNAMIC, &
        !$omp&         is_dynamic,a1,a2,a3,b1,b2,b3)
      !$omp do
      do icel = iS, iE

        ! ----- nodal coordinate & displacement
        iiS = hecMESH%elem_node_index(icel-1)
        nn  = hecMESH%elem_node_index(icel) - iiS
        tt(:) = 0.d0
        do j = 1, nn
          in = hecMESH%elem_node_item(iiS+j)
          nodLOCAL(j) = in
          do i = 1, 3
            ecoord(i,j) = hecMESH%node(3*(in-1)+i)
          enddo
          do i = 1, ndof
            du(ndof*(j-1)+i) = fstrSOLID%dunode(ndof*(in-1)+i)
            u(i,j) = fstrSOLID%unode(ndof*(in-1)+i) + du(ndof*(j-1)+i)
            u_prev(i,j) = fstrSOLID%unode(ndof*(in-1)+i)
          enddo
          if( is_dynamic ) then
            do i = 1, ndof
              vec(ndof*(j-1)+i) = fstrDYNAMIC%VEL(ndof*(in-1)+i,1)
              acc(ndof*(j-1)+i) = fstrDYNAMIC%ACC(ndof*(in-1)+i,1)
            enddo
          endif
          if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres > 0 ) then
            tt(j) = fstrSOLID%temperature( nodLOCAL(j) )
          endif
        enddo

        ! ----- inactive element : assemble dummy stiffness and skip the rest
        if( fstrSOLID%elements(icel)%elemact_flag == kELACT_INACTIVE ) then
          stiff_mat = 0.0d0
          call STF_DUMMY( ndof, nn, ecoord(:,1:nn), u(1:3,1:nn), &
            &  stiff_mat(1:nn*ndof, 1:nn*ndof), fstrSOLID%elements(icel) )
          call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiff_mat)
          cycle
        endif

        ! ----- element K, M, C (section / material context resolved inside)
        stiff_mat = 0.0d0
        mass_mat  = 0.0d0
        damp_mat  = 0.0d0

        call calc_stiff_and_mass_elem( ic_type, nn, ndof, ecoord, u, u_prev, tt, &
          time, tincr, is_dynamic, fstrSOLID, hecMESH, icel, nodLOCAL, &
          stiff_mat, mass_mat, damp_mat, lumped )

        ! ----- assemble system matrix (and the dynamic damping force)
        if( is_dynamic ) then
          ! Newmark intermediate vectors
          do i = 1, nn*ndof
            vecA(i) = -a3*du(i) + a2*vec(i) + a1*acc(i)
            vecB(i) = -b3*du(i) + b2*vec(i) + b1*acc(i)
          enddo

          ! effective dynamic system matrix A = c1*K + c2*M (+ b3*C) and damping force df
          material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
          call calc_damping_mat_and_force_elem( ic_type, nn, ndof, material, fstrDYNAMIC, &
            a3, b3, vecA, vecB, stiff_mat, mass_mat, damp_mat, mat, df )

          call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, mat)

          do j = 1, nn
            do i = 1, ndof
              !$omp atomic
              fstrSOLID%DFORCE(ndof*(nodLOCAL(j)-1)+i) = fstrSOLID%DFORCE(ndof*(nodLOCAL(j)-1)+i)+df(ndof*(j-1)+i)
            enddo
          enddo
        else
          ! static / eigen / frequency : assemble K only
          call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiff_mat)
        endif

      enddo      ! icel
      !$omp end do
      !$omp end parallel
    enddo        ! itype

  end subroutine fstr_CreateMatrix_and_DampingForce

  !---------------------------------------------------------------------*
  !> Compute the element tangent stiffness (and, in the dynamic case, the
  !> element mass and connector damping matrices) for the given element type.
  !!
  !! Special case: 341 elements with the selective ES/NS smoothed FEM option
  !! (kel341SESNS) leave stiff_mat as zero here because their stiffness is
  !! assembled separately through the 881/891 element-type path. The mass
  !! matrix is still computed normally so that the dynamic system matrix
  !! contains the proper inertia contribution.
  subroutine calc_stiff_and_mass_elem( ic_type, nn, ndof, ecoord, u, u_prev, tt, &
      time, tincr, is_dynamic, fstrSOLID, hecMESH, icel, nodLOCAL, &
      stiff_mat, mass_mat, damp_mat, lumped )
  !---------------------------------------------------------------------*
    use m_static_LIB
    use m_static_LIB_shell, only: ShellComposeNodalDisplacement
    use mMechGauss
    use m_dynamic_mass
    use m_fstr_NodalKinematics, only: fstr_get_shell_trial_directors, fstr_get_shell_reference_directors
    use m_fstr_FiniteRotationKinematics, only: fstr_uses_finite_rotation_kinematics

    integer(kind=kint), intent(in)    :: ic_type, ndof
    integer(kind=kint), intent(inout) :: nn  ! some legacy STF_* routines mutate this
    real(kind=kreal),   intent(in)    :: ecoord(:,:), u(:,:), u_prev(:,:), tt(:)
    real(kind=kreal),   intent(in)    :: time, tincr
    logical,            intent(in)    :: is_dynamic
    type(fstr_solid),               intent(inout) :: fstrSOLID
    type(hecmwST_local_mesh),       intent(in)    :: hecMESH
    integer(kind=kint),             intent(in)    :: icel
    integer(kind=kint),             intent(inout) :: nodLOCAL(:)
    real(kind=kreal),   intent(inout) :: stiff_mat(:,:), mass_mat(:,:), damp_mat(:,:)
    real(kind=kreal),   intent(inout) :: lumped(:)

    type(tMaterial), pointer :: material
    integer(kind=kint) :: isect, ihead, cdsys_ID, sec_opt, i, j
    real(kind=kreal) :: coords(3,3), thick
    real(kind=kreal) :: rho, length, surf
    real(kind=kreal) :: shell_nddisp(6,fstrSOLID%max_ncon)
    real(kind=kreal) :: shell_du(6,fstrSOLID%max_ncon)
    real(kind=kreal) :: shell_director(3,fstrSOLID%max_ncon)
    real(kind=kreal) :: shell_ref_director(3,fstrSOLID%max_ncon)

    ! ----- section / material context
    isect = hecMESH%section_ID(icel)
    ihead = hecMESH%section%sect_R_index(isect-1)
    cdsys_ID = hecMESH%section%sect_orien_ID(isect)
    sec_opt = hecMESH%section%sect_opt(isect)
    coords = 0.0d0
    if( cdsys_ID > 0 ) call get_coordsys(cdsys_ID, hecMESH, fstrSOLID, coords)

    material => fstrSOLID%elements(icel)%gausses(1)%pMaterial
    thick = hecMESH%section%sect_R_item(ihead+1)

    if( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 .or. ic_type==2322) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call CreateMat_abort( ic_type, 2 )
      call STF_C2( ic_type, nn, ecoord(1:2,1:nn), fstrSOLID%elements(icel)%gausses(:), thick, &
        stiff_mat(1:nn*ndof,1:nn*ndof), fstrSOLID%elements(icel)%iset, u(1:2,1:nn) )

      if( is_dynamic ) call mass_C2(ic_type, nn, ecoord(1:2,1:nn), fstrSOLID%elements(icel)%gausses, &
        sec_opt, thick, mass_mat, lumped)

    elseif( ic_type==301 ) then
      call STF_C1( ic_type, nn, ecoord(:,1:nn), thick, fstrSOLID%elements(icel)%gausses(:), &
        stiff_mat(1:nn*ndof,1:nn*ndof), u(1:3,1:nn) )

    elseif( ic_type==361 ) then
      if( fstrSOLID%sections(isect)%elemopt361 == kel361FI ) then ! full integration element
        call STF_C3( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          stiff_mat(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
      else if( fstrSOLID%sections(isect)%elemopt361 == kel361BBAR ) then ! B-bar element
        call STF_C3D8Bbar( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          stiff_mat(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
      else if( fstrSOLID%sections(isect)%elemopt361 == kel361IC ) then ! incompatible element
        call STF_C3D8IC( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          stiff_mat(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), &
          fstrSOLID%elements(icel)%aux, tt(1:nn) )
      else if( fstrSOLID%sections(isect)%elemopt361 == kel361FBAR ) then ! F-bar element
        call STF_C3D8Fbar( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          stiff_mat(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
      else if( fstrSOLID%sections(isect)%elemopt361 == kel361UP ) then ! UP (u-p mixed) element
        call STF_C3_up( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          stiff_mat(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, &
          1, fstrSOLID%elements(icel)%p, u(1:3,1:nn), tt(1:nn) )
      endif

      if( is_dynamic ) call mass_C3(ic_type, nn, ecoord(1:3,1:nn), fstrSOLID%elements(icel)%gausses, mass_mat, lumped)

    elseif( ic_type==341 .or. ic_type==351 .or. ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
      ! SESNS option: stiffness is assembled separately via ic_type==881/891;
      ! return stiff_mat = 0 here and still compute the element mass below.
      if( .not. (ic_type==341 .and. fstrSOLID%sections(isect)%elemopt341 == kel341SESNS) ) then
        call STF_C3( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          stiff_mat(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )
      endif

      if( is_dynamic ) call mass_C3(ic_type, nn, ecoord(1:3,1:nn), fstrSOLID%elements(icel)%gausses, mass_mat, lumped)

    else if( ic_type == 511 ) then
      call STF_CONNECTOR( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
        stiff_mat(1:nn*ndof,1:nn*ndof), u(1:3,1:nn), tt(1:nn))

      if( is_dynamic ) call DMP_CONNECTOR( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
        damp_mat(1:nn*ndof,1:nn*ndof), u(1:3,1:nn), tt(1:nn))

    else if( ic_type == 611 ) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call CreateMat_abort( ic_type, 2 )
      call STF_Beam(ic_type, nn, ecoord, hecMESH%section%sect_R_item(ihead+1:), &
        &   material%variables(M_YOUNGS), material%variables(M_POISSON), stiff_mat(1:nn*ndof,1:nn*ndof))

      if( is_dynamic ) then
        surf = hecMESH%section%sect_R_item(ihead+4)
        length = get_length(ecoord(1:3,1:nn))
        rho = material%variables(M_DENSITY)
        call mass_Beam(surf, length, rho, mass_mat)
      endif

    else if( ic_type == 641 ) then
      if( material%nlgeom_flag /= INFINITESIMAL ) call CreateMat_abort( ic_type, 2 )
      call STF_Beam_641(ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses(:), &
        &            hecMESH%section%sect_R_item(ihead+1:), stiff_mat(1:nn*ndof,1:nn*ndof))

      if( is_dynamic ) then
        surf = hecMESH%section%sect_R_item(ihead+4)
        length = get_length(ecoord(1:3,1:nn))
        rho = material%variables(M_DENSITY)
        call mass_Beam_33(surf, length, rho, mass_mat)
      endif

    else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
      if( fstr_uses_finite_rotation_kinematics( ic_type, nn, material ) ) then
        shell_du(:, :) = 0.0d0
        do j = 1, nn
          do i = 1, min(ndof, 6)
            shell_du(i,j) = fstrSOLID%dunode(ndof*(nodLOCAL(j)-1)+i)
          enddo
        enddo
        call ShellComposeNodalDisplacement( ndof, nn, u_prev(1:6,1:nn), shell_du(1:6,1:nn), &
          shell_nddisp(1:6,1:nn) )
        call fstr_get_shell_trial_directors( fstrSOLID, thick, nn, nodLOCAL(1:nn), shell_director(1:3,1:nn) )
        call fstr_get_shell_reference_directors( fstrSOLID, thick, nn, nodLOCAL(1:nn), &
          shell_ref_director(1:3,1:nn) )
        call STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          &              stiff_mat(1:nn*ndof,1:nn*ndof), thick, 0, nddisp=shell_nddisp(1:6,1:nn), &
          &              element=fstrSOLID%elements(icel), nddirector=shell_director(1:3,1:nn), &
          &              ndrefdirector=shell_ref_director(1:3,1:nn))
      else
        if( material%nlgeom_flag /= INFINITESIMAL ) call CreateMat_abort( ic_type, 2 )
        call STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3,1:nn), fstrSOLID%elements(icel)%gausses(:), &
          &              stiff_mat(1:nn*ndof,1:nn*ndof), thick, 0)
      endif

      if( is_dynamic ) then
        rho = material%variables(M_DENSITY)
        call mass_shell(ic_type, nn, ecoord(1:3,1:nn), rho, thick, fstrSOLID%elements(icel)%gausses, mass_mat, lumped)
      endif

    else if( ic_type == 761 ) then   ! for shell-solid mixed analysis
      if( material%nlgeom_flag /= INFINITESIMAL ) call CreateMat_abort( ic_type, 2 )
      call STF_Shell_MITC(731, 3, 6, ecoord(1:3,1:3), fstrSOLID%elements(icel)%gausses(:), &
        &              stiff_mat(1:nn*ndof,1:nn*ndof), thick, 2)

      if( is_dynamic ) then
        surf = get_face3(ecoord(1:3,1:nn))
        rho = material%variables(M_DENSITY)
        call mass_S3(surf, thick, rho, mass_mat)
      endif

    else if( ic_type == 781 ) then   ! for shell-solid mixed analysis
      if( material%nlgeom_flag /= INFINITESIMAL ) call CreateMat_abort( ic_type, 2 )
      call STF_Shell_MITC(741, 4, 6, ecoord(1:3,1:4), fstrSOLID%elements(icel)%gausses(:), &
        &              stiff_mat(1:nn*ndof,1:nn*ndof), thick, 1)

      if( is_dynamic ) then
        surf = get_face4(ecoord(1:3,1:nn))
        rho = material%variables(M_DENSITY)
        call mass_S4(surf, thick, rho, mass_mat)
      endif

    elseif( ic_type==3414 ) then
      if( material%mtype /= INCOMP_NEWTONIAN ) call CreateMat_abort( ic_type, 3, material%mtype )
      call STF_C3_vp( ic_type, nn, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
        stiff_mat(1:nn*ndof,1:nn*ndof), tincr, u_prev(1:4,1:nn) )

    else if( ic_type == 881 .or. ic_type == 891 ) then  ! for selective es/ns smoothed fem
      call STF_C3D4_SESNS( ic_type, nn, nodLOCAL, ecoord(:,1:nn), fstrSOLID%elements(icel)%gausses(:), &
        stiff_mat, cdsys_ID, coords, time, tincr, u(1:3,1:nn), tt(1:nn) )

    else
      call CreateMat_abort( ic_type, 1 )
    endif

  end subroutine calc_stiff_and_mass_elem

  !---------------------------------------------------------------------*
  !> Combine the element K, M (and connector C) into the effective dynamic
  !> system matrix mat = c1*K + c2*M (+ b3*C) and compute the corresponding
  !> RHS contribution df.
  !!
  !! Per-material Rayleigh damping (is_elem_Rayleigh_damping) overrides the
  !! global fstrDYNAMIC%ray_{m,k} for the current element.
  subroutine calc_damping_mat_and_force_elem( ic_type, nn, ndof, material, fstrDYNAMIC, &
      a3, b3, vecA, vecB, stiff_mat, mass_mat, damp_mat, mat, df )
  !---------------------------------------------------------------------*
    integer(kind=kint), intent(in) :: ic_type, nn, ndof
    type(tMaterial), pointer, intent(in)    :: material
    type(fstr_dynamic),       intent(in)    :: fstrDYNAMIC
    real(kind=kreal),         intent(in)    :: a3, b3
    real(kind=kreal),         intent(in)    :: vecA(:), vecB(:)
    real(kind=kreal),         intent(in)    :: stiff_mat(:,:), mass_mat(:,:), damp_mat(:,:)
    real(kind=kreal),         intent(out)   :: mat(:,:), df(:)

    integer(kind=kint) :: i, j
    real(kind=kreal)   :: ray_m, ray_k, c1, c2
    real(kind=kreal)   :: vecC(20*6), Kb(20*6)

    if( material%is_elem_Rayleigh_damping ) then
      ray_m = material%variables(M_DAMPING_RM)
      ray_k = material%variables(M_DAMPING_RK)
    else
      ray_m = fstrDYNAMIC%ray_m
      ray_k = fstrDYNAMIC%ray_k
    endif

    ! ----- LHS: A = c1*K + c2*M (+ b3*C for connector)
    c1 = 1.d0 + ray_k*b3
    c2 = a3   + ray_m*b3
    do i = 1, nn*ndof
      do j = 1, nn*ndof
        mat(j,i) = c1*stiff_mat(j,i) + c2*mass_mat(j,i)
      enddo
    enddo
    if( ic_type == 511 ) then
      do i = 1, nn*ndof
        do j = 1, nn*ndof
          mat(j,i) = mat(j,i) + b3*damp_mat(j,i)
        enddo
      enddo
    endif

    ! ----- RHS: df = M*(vecA + ray_m*vecB) + ray_k*K*vecB (+ C*vecB for connector)
    do i = 1, nn*ndof
      vecC(i) = vecA(i) + ray_m*vecB(i)
    enddo

    Kb = 0.0d0
    do i = 1, nn*ndof
      do j = 1, nn*ndof
        Kb(i) = Kb(i) + stiff_mat(i,j)*vecB(j)
      enddo
    enddo

    df = 0.0d0
    df(1:nn*ndof) = matmul(mass_mat(1:nn*ndof,1:nn*ndof), vecC(1:nn*ndof))
    do i = 1, nn*ndof
      df(i) = df(i) + ray_k*Kb(i)
    enddo
    if( ic_type == 511 ) then
      do i = 1, nn*ndof
        do j = 1, nn*ndof
          df(i) = df(i) + damp_mat(i,j)*vecB(j)
        enddo
      enddo
    endif

  end subroutine calc_damping_mat_and_force_elem

  subroutine CreateMat_abort( ic_type, flag, mtype )
    integer(kind=kint), intent(in)           :: ic_type
    integer(kind=kint), intent(in)           :: flag
    integer(kind=kint), intent(in), optional :: mtype

    if( flag == 1 ) then
      write(*,*) '###ERROR### : Element type not supported for static analysis'
    else if( flag == 2 ) then
      write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
    else if( flag == 3 ) then
      write(*,*) '###ERROR### : This element is not supported for this material'
    endif
    write(*,*) ' ic_type = ', ic_type
    if( present(mtype) ) write(*,*) ' mtype = ', mtype
    call hecmw_abort(hecmw_comm_get_comm())
  end subroutine CreateMat_abort

end module m_fstr_CreateMatrix_and_DampingForce
