!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module contains several strategy to free locking problem
!> in Eight-node hexagonal element
module m_static_LIB_C3D4SESNS

  use hecmw, only : kint, kreal
  use elementInfo
  use mMechGauss
  use m_MatMatrix
  use m_common_struct
  use m_static_LIB_3d
  use m_utilities

  implicit none
  private
  public :: STF_C3D4_SESNS, UPDATE_C3_SESNS, Return_nn_comp_C3D4_SESNS

  logical, parameter :: DEBUG=.false.

contains

  subroutine create_compressed_vector_by_node_id(nn,nodlocal,nn_comp,local_nid)
    integer(kind=kint), intent(in)     :: nn           !< number of elemental nodes
    integer(kind=kint), intent(inout)  :: nodlocal(:)  !<
    integer(kind=kint), intent(out)    :: nn_comp      !<
    integer(kind=kint), intent(out)    :: local_nid(:) !<

    integer(kind=kint) :: i, j, nid
    logical            :: id_found
    integer(kind=kint) :: cmp_nodlocal(nn)

    nn_comp = 0
    cmp_nodlocal(1:nn) = 0

    do i=1,nn
      nid = nodlocal(i)
      id_found = .false.
      do j=1,nn_comp
        if( cmp_nodlocal(j) /= nid ) cycle
        id_found = .true.
        local_nid(i) = j
        exit
      end do
      if( id_found ) cycle
      nn_comp = nn_comp + 1
      cmp_nodlocal(nn_comp) = nid
      local_nid(i) = nn_comp
    end do

    if(DEBUG) then
      write(*,*) "nn:",nn
      write(*,*) "nodlocal:",nodlocal(1:nn)
    end if

    nodlocal(1:nn_comp) = cmp_nodlocal(1:nn_comp)

    if(DEBUG) then
      write(*,*) "nn_comp:",nn_comp
      write(*,*) "nodlocal(comp):",nodlocal(1:nn_comp)
      write(*,*) "local_nid:",local_nid(1:nn)
    end if

  end subroutine

  integer(kind=kint) function Return_nn_comp_C3D4_SESNS(nn, nodlocal)
    integer(kind=kint), intent(in)    :: nn                !< number of elemental nodes
    integer(kind=kint), intent(inout) :: nodlocal(:)       !< node id ( used to compress data )

    integer(kind=kint) :: nn_comp, local_nid(nn)

    call create_compressed_vector_by_node_id(nn,nodlocal,nn_comp,local_nid)
    Return_nn_comp_C3D4_SESNS = nn_comp
  end function

  subroutine Get_Compressed_gderiv_and_disp_C3D4_SESNS  &
      &  (flag, etype, nn, nodlocal, ecoord, gderiv, det, u, du, tmpu, tmpdu, gderiv1)
    integer(kind=kint), intent(in)    :: flag
    integer(kind=kint), intent(in)    :: etype             !< element type
    integer(kind=kint), intent(inout) :: nn                !< number of elemental nodes
    integer(kind=kint), intent(inout) :: nodlocal(:)       !< node id ( used to compress data )
    real(kind=kreal),   intent(in)    :: ecoord(3,nn)      !< coordinates of elemental nodes
    real(kind=kreal), intent(out)     :: gderiv(nn, 3)
    real(kind=kreal), intent(out)     :: det
    real(kind=kreal), intent(out)     :: u(:,:)
    real(kind=kreal), intent(out)     :: du(:,:)
    real(kind=kreal), intent(in)      :: tmpu(:,:)          !< nodal displacemwent
    real(kind=kreal), intent(in), optional   :: tmpdu(:,:)         !< nodal displacemwent
    real(kind=kreal), intent(out), optional  :: gderiv1(nn, 3)

    integer(kind=kint) :: i, j, LX
    integer(kind=kint) :: n_subelem, nn_comp, local_nid(nn)
    real(kind=kreal) :: elem(3, 4)
    real(kind=kreal) :: tmpgderiv(2*nn,3), ratio(nn)
    real(kind=kreal) :: naturalCoord(3)
    real(kind=kreal) :: coord_comp(3,nn), F(3,3)

    !compress node id
    if( etype == 881 ) then !node-smoothed element
      n_subelem = (nn-1)/3
    else if( etype == 891 ) then !edge-smoothed element
      n_subelem = (nn-2)/2
    end if
    call create_compressed_vector_by_node_id(nn,nodlocal,nn_comp,local_nid)

    !get smoothed global derivative
    call getQuadPoint( fe_tet4n, 1, naturalCoord(:) )
    det = 0.d0
    ratio(:) = 0.d0

    if( etype == 881 ) then !node-smoothed element
      elem(1:3,1) = ecoord(1:3,1)
      do i=1,n_subelem
        elem(1:3,2:4) = ecoord(1:3,3*i-1:3*i+1)
        call getGlobalDeriv(fe_tet4n, 4, naturalcoord, elem, ratio(i), tmpgderiv(4*i-3:4*i,1:3))
        ratio(i) = ratio(i)/24.d0
        det = det + ratio(i)
      end do
    else if( etype == 891 ) then !edge-smoothed element
      elem(1:3,1:2) = ecoord(1:3,1:2)
      do i=1,n_subelem
        elem(1:3,3:4) = ecoord(1:3,2*i+1:2*i+2)
        call getGlobalDeriv(fe_tet4n, 4, naturalcoord, elem, ratio(i), tmpgderiv(4*i-3:4*i,1:3))
        ratio(i) = ratio(i)/36.d0
        det = det + ratio(i)
      end do
    end if

    ratio(1:n_subelem) = ratio(1:n_subelem)/det
    do i=1,n_subelem
      tmpgderiv(4*i-3:4*i,1:3) = ratio(i)*tmpgderiv(4*i-3:4*i,1:3)
    end do

    ! compress
    !! compress gderiv
    gderiv(:,:) = 0.d0
    if( etype == 881 ) then !node-smoothed element
      do i=1,n_subelem
        gderiv(1,1:3) = gderiv(1,1:3) + tmpgderiv(4*i-3,1:3)
        do j=2,4
          LX = local_nid(3*i-3+j)
          gderiv(LX,1:3) = gderiv(LX,1:3) + tmpgderiv(4*i+j-4,1:3)
        end do
      end do
    else if( etype == 891 ) then !edge-smoothed element
      do i=1,n_subelem
        gderiv(1:2,1:3) = gderiv(1:2,1:3) + tmpgderiv(4*i-3:4*i-2,1:3)
        do j=3,4
          LX = local_nid(2*i-2+j)
          gderiv(LX,1:3) = gderiv(LX,1:3) + tmpgderiv(4*i+j-4,1:3)
        end do
      end do
    end if

    if( flag == UPDATELAG ) then
      do i=1,nn
        LX = local_nid(i)
        u(1:3,LX) = tmpu(1:3,i)
        coord_comp(1:3,LX) = ecoord(1:3,i)
      end do
      if( present(tmpdu) ) then
        do i=1,nn
          LX = local_nid(i)
          du(1:3,LX) = tmpdu(1:3,i)
        end do
      else
        du(1:3,1:nn_comp) = 0.d0
      end if

      F(1:3,1:3) = matmul(coord_comp(1:3,1:nn_comp)+u(1:3,1:nn_comp)+du(1:3,1:nn_comp),gderiv(1:nn_comp,1:3))
      det = det * Determinant33( F )
      call calInverse(3,F)
      ! dN/dx_{n+1}
      gderiv1(1:nn_comp,1:3) = matmul(gderiv(1:nn_comp,1:3),F(1:3,1:3))

      F(1:3,1:3) = matmul(coord_comp(1:3,1:nn_comp)+u(1:3,1:nn_comp)+0.5d0*du(1:3,1:nn_comp),gderiv(1:nn_comp,1:3))
      call calInverse(3,F)
      ! dN/dx_{n+1/2}
      gderiv(1:nn_comp,1:3) = matmul(gderiv(1:nn_comp,1:3),F(1:3,1:3))
    else
      do i=1,nn
        LX = local_nid(i)
        u(1:3,LX) = tmpu(1:3,i)
      end do
    end if

    !! compress num of nodes
    nn = nn_comp
  end subroutine

  subroutine STF_C3D4_SESNS                            &
      (etype, nn, nodlocal, ecoord, gausses, stiff, cdsys_ID, coords, &
      time, tincr, tmpu ,temperature)
    !----------------------------------------------------------------------*

    use mMechGauss
    use m_MatMatrix
    use m_common_struct

    !---------------------------------------------------------------------

    integer(kind=kint), intent(in)  :: etype                  !< element type
    integer(kind=kint), intent(inout)  :: nn                  !< number of elemental nodes
    integer(kind=kint), intent(inout)  :: nodlocal(:)         !< node id ( used to compress data )
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)           !< coordinates of elemental nodes
    type(tGaussStatus), intent(in)  :: gausses(:)             !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:,:)             !< stiff matrix
    integer(kind=kint), intent(in)  :: cdsys_ID
    real(kind=kreal), intent(inout) :: coords(3,3)            !< variables to define matreial coordinate system
    real(kind=kreal), intent(in)    :: time                   !< current time
    real(kind=kreal), intent(in)    :: tincr                  !< time increment
    real(kind=kreal), intent(in)    :: temperature(nn) !< temperature
    real(kind=kreal), intent(in), optional :: tmpu(:,:)       !< nodal displacemwent

    !---------------------------------------------------------------------

    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal) :: D(6, 6), B(6, NDOF*nn), DB(6, NDOF*nn)
    real(kind=kreal) :: gderiv(nn, 3), stress(6), mat(6, 6), det
    real(kind=kreal) :: wg
    integer(kind=kint) :: i, j, serr
    real(kind=kreal) :: temp
    real(kind=kreal) :: gdispderiv(3, 3)
    real(kind=kreal) :: B1(6, NDOF*nn), coordsys(3, 3)
    real(kind=kreal) :: Smat(9, 9)
    real(kind=kreal) :: BN(9, NDOF*nn), SBN(9, NDOF*nn)
    real(kind=kreal) :: gderiv1(nn,3)

    !---------------------------------------------------------------------

    real(kind=kreal) :: u(3,nn), du(3,nn)

    !---------------------------------------------------------------------

    stiff(:, :) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag
    if( .not. present(tmpu) ) flag = INFINITESIMAL    ! enforce to infinitesimal deformation analysis

    ! compress id and calc gderiv
    if( present(tmpu) ) then
      if( flag == UPDATELAG ) then
        call Get_Compressed_gderiv_and_disp_C3D4_SESNS &
          &  (flag, etype, nn, nodlocal, ecoord, gderiv, det, u, du, tmpu, gderiv1=gderiv1)
      else
        call Get_Compressed_gderiv_and_disp_C3D4_SESNS &
          &  (flag, etype, nn, nodlocal, ecoord, gderiv, det, u, du, tmpu)
      end if
    end if

    ! calc stiffness matrix
    if( cdsys_ID > 0 ) then
      call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:, :), serr )
      if( serr == -1 ) stop "Fail to setup local coordinate"
      if( serr == -2 ) then
        write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
      end if
    end if

    if( etype == 881 ) then !node-smoothed element
      temp = temperature(1)
    else if( etype == 891 ) then !edge-smoothed element
      temp = 0.5d0*(temperature(1)+temperature(2))
    end if
    if( etype == 881 ) then !node-smoothed element
      call MatlMatrix( gausses(1), D3, D, time, tincr, coordsys, temp, hdflag=2 )
    else if( etype == 891 ) then !edge-smoothed element
      call MatlMatrix( gausses(1), D3, D, time, tincr, coordsys, temp, hdflag=1 )
    end if

    if( flag == UPDATELAG ) then
      call GEOMAT_C3( gausses(1)%stress, mat )
      D(:, :) = D(:, :)-mat
    endif

    wg = det
    B(1:6, 1:nn*ndof) = 0.0D0
    do j = 1, nn
      B(1, 3*j-2)=gderiv(j, 1)
      B(2, 3*j-1)=gderiv(j, 2)
      B(3, 3*j  )=gderiv(j, 3)
      B(4, 3*j-2)=gderiv(j, 2)
      B(4, 3*j-1)=gderiv(j, 1)
      B(5, 3*j-1)=gderiv(j, 3)
      B(5, 3*j  )=gderiv(j, 2)
      B(6, 3*j-2)=gderiv(j, 3)
      B(6, 3*j  )=gderiv(j, 1)
    enddo

    ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
    if( flag == TOTALLAG ) then
      ! ---dudx(i, j) ==> gdispderiv(i, j)
      gdispderiv(1:ndof, 1:ndof) = matmul( u(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
      B1(1:6, 1:nn*NDOF)=0.0D0
      do j=1, nn
        B1(1, 3*j-2) = gdispderiv(1, 1)*gderiv(j, 1)
        B1(1, 3*j-1) = gdispderiv(2, 1)*gderiv(j, 1)
        B1(1, 3*j  ) = gdispderiv(3, 1)*gderiv(j, 1)
        B1(2, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 2)
        B1(2, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 2)
        B1(2, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 2)
        B1(3, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 3)
        B1(3, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 3)
        B1(3, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 3)
        B1(4, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 2)
        B1(4, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 2)
        B1(4, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 2)
        B1(5, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 3)+gdispderiv(1, 3)*gderiv(j, 2)
        B1(5, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 3)+gdispderiv(2, 3)*gderiv(j, 2)
        B1(5, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 3)+gdispderiv(3, 3)*gderiv(j, 2)
        B1(6, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 3)
        B1(6, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 3)
        B1(6, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 3)
      end do
      ! ---BL = BL0 + BL1
      do j=1, nn*ndof
        B(:, j) = B(:, j)+B1(:, j)
      end do
    end if

    DB(1:6, 1:nn*ndof) = matmul( D, B(1:6, 1:nn*ndof) )
    do j=1,nn*ndof 
      do i=1,nn*ndof
        stiff(i, j) = stiff(i, j)+dot_product( B(:, i),  DB(:, j) )*WG
      enddo
    enddo

    ! calculate the stress matrix ( TOTAL LAGRANGE METHOD )
    if( flag == TOTALLAG .OR. flag==UPDATELAG ) then
      stress(1:6) = gausses(1)%stress(1:6)

      BN(1:9, 1:nn*ndof) = 0.0D0
      do j = 1,  nn
        BN(1, 3*j-2) = gderiv(j, 1)
        BN(2, 3*j-1) = gderiv(j, 1)
        BN(3, 3*j  ) = gderiv(j, 1)
        BN(4, 3*j-2) = gderiv(j, 2)
        BN(5, 3*j-1) = gderiv(j, 2)
        BN(6, 3*j  ) = gderiv(j, 2)
        BN(7, 3*j-2) = gderiv(j, 3)
        BN(8, 3*j-1) = gderiv(j, 3)
        BN(9, 3*j  ) = gderiv(j, 3)
      end do
      Smat(:, :) = 0.0D0
      do j = 1, 3
        Smat(j  , j  ) = stress(1)
        Smat(j  , j+3) = stress(4)
        Smat(j  , j+6) = stress(6)
        Smat(j+3, j  ) = stress(4)
        Smat(j+3, j+3) = stress(2)
        Smat(j+3, j+6) = stress(5)
        Smat(j+6, j  ) = stress(6)
        Smat(j+6, j+3) = stress(5)
        Smat(j+6, j+6) = stress(3)
      end do
      SBN(1:9, 1:nn*ndof) = matmul( Smat(1:9, 1:9),  BN(1:9, 1:nn*ndof) )
      do j=1,nn*ndof 
        do i=1,nn*ndof
          stiff(i, j) = stiff(i, j)+dot_product( BN(:, i),  SBN(:, j) )*WG
        end do
      end do

    end if

  end subroutine

  !> Update strain and stress inside element
  !---------------------------------------------------------------------*
  subroutine UPDATE_C3_SESNS                     &
      (etype, nn, nodlocal, ecoord, u, ddu, cdsys_ID, coords, qf, &
      gausses, time, tincr, TT, T0, TN)
    !---------------------------------------------------------------------*

    use m_fstr
    use mMaterial
    use mMechGauss
    use m_MatMatrix
    use m_ElastoPlastic
    use m_utilities

    integer(kind=kint), intent(in)    :: etype         !< \param [in] element type
    integer(kind=kint), intent(inout) :: nn            !< \param [in] number of elemental nodes
    integer(kind=kint), intent(inout) :: nodlocal(:)   !< \param [inout] node id ( used to compress data )
    real(kind=kreal), intent(in)      :: ecoord(3, nn) !< \param [in] coordinates of elemental nodes
    real(kind=kreal), intent(in)      :: u(3, nn)      !< \param [in] nodal dislplacements
    real(kind=kreal), intent(in)      :: ddu(3, nn)    !< \param [in] nodal displacement
    integer(kind=kint), intent(in)    :: cdsys_ID
    real(kind=kreal), intent(inout)   :: coords(3, 3)  !< variables to define matreial coordinate system
    real(kind=kreal), intent(out)     :: qf(nn*3)      !< \param [out] Internal Force
    type(tGaussStatus), intent(inout) :: gausses(:)    !< \param [out] status of qudrature points
    real(kind=kreal), intent(in)      :: time          !< current time
    real(kind=kreal), intent(in)      :: tincr         !< time increment
    real(kind=kreal), intent(in)      :: TT(nn)   !< current temperature
    real(kind=kreal), intent(in)      :: T0(nn)   !< reference temperature
    real(kind=kreal), intent(in)      :: TN(nn)   !< reference temperature

    ! LCOAL VARIAVLES
    integer(kind=kint) :: flag
    integer(kind=kint), parameter :: ndof = 3
    real(kind=kreal)   :: B(6,ndof*nn), B1(6,ndof*nn), ina(1)
    real(kind=kreal)   :: gderiv(nn,3), gderiv1(nn,3), gdispderiv(3,3), F(3,3), det, WG, ttc,tt0, ttn
    integer(kind=kint) :: j, serr
    real(kind=kreal)   :: rot(3,3), EPSTH(6)
    real(kind=kreal)   :: totaldisp(3,nn), totalddisp(3,nn), coordsys(3,3)
    real(kind=kreal)   :: dstrain(6)
    real(kind=kreal)   :: alpo(3)
    logical            :: ierr, matlaniso

    real(kind=kreal)   :: stress(6)

    qf(:) = 0.0D0
    ! we suppose the same material type in the element
    flag = gausses(1)%pMaterial%nlgeom_flag

    ! compress id and calc gderiv
    if( flag == UPDATELAG ) then
      call Get_Compressed_gderiv_and_disp_C3D4_SESNS &
        &  (flag, etype, nn, nodlocal, ecoord, gderiv, det, totaldisp, totalddisp, u, ddu, gderiv1)
    else
      call Get_Compressed_gderiv_and_disp_C3D4_SESNS &
        &  (flag, etype, nn, nodlocal, ecoord, gderiv, det, totaldisp, totalddisp, u+ddu )
    end if

    matlaniso = .FALSE.
    ina = TT(1)
    call fetch_TableData( MC_ORTHOEXP, gausses(1)%pMaterial%dict, alpo(:), ierr, ina )
    if( .not. ierr ) matlaniso = .true.

    if( cdsys_ID > 0 ) then
      call set_localcoordsys( coords, g_LocalCoordSys(cdsys_ID), coordsys(:,:), serr )
      if( serr == -1 ) stop "Fail to setup local coordinate"
      if( serr == -2 ) then
        write(*, *) "WARNING! Cannot setup local coordinate, it is modified automatically"
      end if
    end if

    ! ========================================================
    ! UPDATE STRAIN and STRESS
    ! ========================================================

    ! Thermal Strain
    EPSTH = 0.0D0
    if( etype == 881 ) then !node-smoothed element
      ttc = TT(1)
      tt0 = T0(1)
      ttn = TN(1)
    else if( etype == 891 ) then !edge-smoothed element
      ttc = 0.5d0*(TT(1)+TT(2))
      tt0 = 0.5d0*(T0(1)+T0(2))
      ttn = 0.5d0*(TN(1)+TN(2))
    end if
    call Cal_Thermal_expansion_C3( tt0, ttc, gausses(1)%pMaterial, coordsys, matlaniso, EPSTH )

    ! Update strain
    ! Small strain
    if( flag == UPDATELAG ) then
      gdispderiv(1:ndof, 1:ndof) = matmul( totalddisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
    else
      gdispderiv(1:ndof, 1:ndof) = matmul( totaldisp(1:ndof, 1:nn), gderiv(1:nn, 1:ndof) )
    end if

    dstrain(1) = gdispderiv(1, 1)
    dstrain(2) = gdispderiv(2, 2)
    dstrain(3) = gdispderiv(3, 3)
    dstrain(4) = ( gdispderiv(1, 2)+gdispderiv(2, 1) )
    dstrain(5) = ( gdispderiv(2, 3)+gdispderiv(3, 2) )
    dstrain(6) = ( gdispderiv(3, 1)+gdispderiv(1, 3) )
    dstrain(:) = dstrain(:)-EPSTH(:)   ! allright?

    F(1:3,1:3) = 0.d0; F(1,1)=1.d0; F(2,2)=1.d0; F(3,3)=1.d0; !deformation gradient
    if( flag == INFINITESIMAL ) then
      gausses(1)%strain(1:6) = dstrain(1:6)+EPSTH(1:6)

    else if( flag == TOTALLAG ) then
      ! Green-Lagrange strain
      dstrain(1) = dstrain(1)+0.5d0*dot_product( gdispderiv(:, 1), gdispderiv(:, 1) )
      dstrain(2) = dstrain(2)+0.5d0*dot_product( gdispderiv(:, 2), gdispderiv(:, 2) )
      dstrain(3) = dstrain(3)+0.5d0*dot_product( gdispderiv(:, 3), gdispderiv(:, 3) )
      dstrain(4) = dstrain(4)+( gdispderiv(1, 1)*gdispderiv(1, 2)                                     &
        +gdispderiv(2, 1)*gdispderiv(2, 2)+gdispderiv(3, 1)*gdispderiv(3, 2) )
      dstrain(5) = dstrain(5)+( gdispderiv(1, 2)*gdispderiv(1, 3)                                     &
        +gdispderiv(2, 2)*gdispderiv(2, 3)+gdispderiv(3, 2)*gdispderiv(3, 3) )
      dstrain(6) = dstrain(6)+( gdispderiv(1, 1)*gdispderiv(1, 3)                                     &
        +gdispderiv(2, 1)*gdispderiv(2, 3)+gdispderiv(3, 1)*gdispderiv(3, 3) )

      gausses(1)%strain(1:6) = dstrain(1:6)+EPSTH(:)
      F(1:3,1:3) = F(1:3,1:3) + gdispderiv(1:3,1:3)

    else if( flag == UPDATELAG ) then
      rot = 0.0D0
      rot(1, 2)= 0.5d0*(gdispderiv(1, 2)-gdispderiv(2, 1) );  rot(2, 1) = -rot(1, 2)
      rot(2, 3)= 0.5d0*(gdispderiv(2, 3)-gdispderiv(3, 2) );  rot(3, 2) = -rot(2, 3)
      rot(1, 3)= 0.5d0*(gdispderiv(1, 3)-gdispderiv(3, 1) );  rot(3, 1) = -rot(1, 3)

      gausses(1)%strain(1:6) = gausses(1)%strain_bak(1:6)+dstrain(1:6)+EPSTH(:)

      F(1:3,1:3) = F(1:3,1:3) + matmul( totaldisp(1:ndof, 1:nn)+totalddisp(1:ndof, 1:nn), gderiv1(1:nn, 1:ndof) )
      gderiv(1:nn, 1:ndof) = gderiv1(1:nn, 1:ndof)
    end if

    ! Update stress
    if( etype == 881 ) then !node-smoothed element
      call Update_Stress3D( flag, gausses(1), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn, hdflag=2 )
    else if( etype == 891 ) then !edge-smoothed element
      call Update_Stress3D( flag, gausses(1), rot, dstrain, F, coordsys, time, tincr, ttc, tt0, ttn, hdflag=1 )
    end if

    ! ========================================================
    ! calculate the internal force ( equivalent nodal force )
    ! ========================================================

    ! Small strain
    B(1:6, 1:nn*ndof) = 0.0D0
    do J=1,nn
      B(1,3*j-2) = gderiv(j, 1)
      B(2,3*j-1) = gderiv(j, 2)
      B(3,3*j  ) = gderiv(j, 3)
      B(4,3*j-2) = gderiv(j, 2)
      B(4,3*j-1) = gderiv(j, 1)
      B(5,3*j-1) = gderiv(j, 3)
      B(5,3*j  ) = gderiv(j, 2)
      B(6,3*j-2) = gderiv(j, 3)
      B(6,3*j  ) = gderiv(j, 1)
    end do

    ! calculate the BL1 matrix ( TOTAL LAGRANGE METHOD )
    if( flag == TOTALLAG ) then
      B1(1:6, 1:nn*ndof)=0.0D0
      do j = 1,nn
        B1(1, 3*j-2) = gdispderiv(1, 1)*gderiv(j, 1)
        B1(1, 3*j-1) = gdispderiv(2, 1)*gderiv(j, 1)
        B1(1, 3*j  ) = gdispderiv(3, 1)*gderiv(j, 1)
        B1(2, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 2)
        B1(2, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 2)
        B1(2, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 2)
        B1(3, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 3)
        B1(3, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 3)
        B1(3, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 3)
        B1(4, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 2)
        B1(4, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 2)
        B1(4, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 2)
        B1(5, 3*j-2) = gdispderiv(1, 2)*gderiv(j, 3)+gdispderiv(1, 3)*gderiv(j, 2)
        B1(5, 3*j-1) = gdispderiv(2, 2)*gderiv(j, 3)+gdispderiv(2, 3)*gderiv(j, 2)
        B1(5, 3*j  ) = gdispderiv(3, 2)*gderiv(j, 3)+gdispderiv(3, 3)*gderiv(j, 2)
        B1(6, 3*j-2) = gdispderiv(1, 3)*gderiv(j, 1)+gdispderiv(1, 1)*gderiv(j, 3)
        B1(6, 3*j-1) = gdispderiv(2, 3)*gderiv(j, 1)+gdispderiv(2, 1)*gderiv(j, 3)
        B1(6, 3*j  ) = gdispderiv(3, 3)*gderiv(j, 1)+gdispderiv(3, 1)*gderiv(j, 3)
      end do
      ! BL = BL0 + BL1
      do j=1,nn*ndof
        B(:,j) = B(:,j)+B1(:,j)
      end do

    end if

    ! calculate the Internal Force
    stress(1:6) = gausses(1)%stress(1:6)

    WG=det
    qf(1:nn*ndof) = qf(1:nn*ndof)+matmul( stress(1:6), B(1:6,1:nn*ndof) )*WG

  end subroutine

end module
