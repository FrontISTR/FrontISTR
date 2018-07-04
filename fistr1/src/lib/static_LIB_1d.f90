!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>   This module provide common functions of 3D truss elements
module m_static_LIB_1d
  use hecmw, only : kint, kreal
  use elementInfo
  implicit none

contains
  !
  !=====================================================================*
  !  STF_C1
  !=====================================================================*
  !>  This subroutine calculate stiff matrix of 2-nodes truss element
  subroutine STF_C1( etype,nn,ecoord,area,gausses,stiff, u ,temperature )
    use mMechGauss
    integer(kind=kint), intent(in)  :: etype               !< element type
    integer(kind=kint), intent(in)  :: nn                  !< number of elemental nodes
    real(kind=kreal),   intent(in)  :: ecoord(3,nn)        !< coordinates of elemental nodes
    real(kind=kreal),   intent(in)  :: area                !< section area
    type(tGaussStatus), intent(in)  :: gausses(:)          !< status of qudrature points
    real(kind=kreal),   intent(out) :: stiff(:,:)          !< stiff matrix
    real(kind=kreal),   intent(in), optional :: u(:,:)     !< nodal displacemwent
    real(kind=kreal),   intent(in), optional :: temperature(nn)     !< temperature

    real(kind=kreal) DET,WG, llen, llen0, elem(3,nn)
    logical :: ierr
    real(kind=kreal) ina(1), outa(1), direc(3), direc0(3), coeff, strain
    integer(kind=kint) :: i,j

    ! we suppose the same material type in the element
    if( present(u) ) then
      elem(:,:) = ecoord(:,:) + u(:,:)
    else
      elem(:,:) = ecoord(:,:)
    endif

    direc = elem(:,2)-elem(:,1)
    llen = dsqrt( dot_product(direc, direc) )
    direc = direc/llen
    direc0 = ecoord(:,2)-ecoord(:,1)
    llen0 = dsqrt( dot_product(direc0, direc0) )

    if( present(temperature) ) then
      ina(1) = 0.5d0*(temperature(1)+temperature(2))
      call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr, ina )
    else
      call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr )
    endif
    if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_YOUNGS)
    coeff = outa(1)*area*llen0/(llen*llen)
    strain = gausses(1)%strain(1)

    stiff(:,:) = 0.d0
    do i=1,3
      stiff(i,i) = coeff*strain
      do j=1,3
        stiff(i,j) = stiff(i,j) + coeff*(1.d0-2.d0*strain)*direc(i)*direc(j)
      enddo
    enddo

    stiff(4:6,1:3) = -stiff(1:3,1:3)
    stiff(1:3,4:6) = transpose(stiff(4:6,1:3))
    stiff(4:6,4:6) = stiff(1:3,1:3)

  end subroutine STF_C1

  !
  !> Update strain and stress inside element
  !---------------------------------------------------------------------*
  subroutine UPDATE_C1( etype, nn, ecoord, area, u, du, qf ,gausses, TT, T0 )
    !---------------------------------------------------------------------*
    use m_fstr
    use mMechGauss
    ! I/F VARIAVLES
    integer(kind=kint), intent(in)     :: etype           !< \param [in] element type
    integer(kind=kint), intent(in)     :: nn              !< \param [in] number of elemental nodes
    real(kind=kreal),   intent(in)     :: ecoord(3,nn)    !< \param [in] coordinates of elemental nodes
    real(kind=kreal),   intent(in)     :: area            !< section area
    real(kind=kreal),   intent(in)     :: u(3,nn)         !< \param [in] nodal dislplacements
    real(kind=kreal),   intent(in)     :: du(3,nn)       !< \param [in] nodal displacement ( solutions of solver )
    real(kind=kreal),   intent(out)    :: qf(nn*3)        !< \param [out] Internal Force
    type(tGaussStatus), intent(inout)  :: gausses(:)      !< \param [out] status of qudrature points
    real(kind=kreal),   intent(in), optional :: TT(nn)    !< current temperature
    real(kind=kreal),   intent(in), optional :: T0(nn)    !< reference temperature

    ! LCOAL VARIAVLES
    real(kind=kreal)   :: direc(3), direc0(3)
    real(kind=kreal)   :: llen, llen0, ina(1), outa(1)
    real(kind=kreal)   :: elem(3,nn)
    real(kind=kreal)   :: young
    real(kind=kreal)   :: ttc, tt0, alp, alp0, epsth
    logical            :: ierr

    qf(:)              = 0.d0
    ! we suppose the same material type in the element
    elem(:,:) = ecoord(:,:) + u(:,:) + du(:,:)

    direc = elem(:,2)-elem(:,1)
    llen = dsqrt( dot_product(direc, direc) )
    direc = direc/llen
    direc0 = ecoord(:,2)-ecoord(:,1)
    llen0 = dsqrt( dot_product(direc0, direc0) )

    epsth = 0.d0
    if( present(tt) .and. present(t0) ) then
      ttc = 0.5d0*(TT(1)+TT(2))
      tt0 = 0.5d0*(T0(1)+T0(2))

      ina(1) = ttc
      call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr, ina )
      if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_YOUNGS)
      young = outa(1)

      call fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa(:), ierr, ina )
      if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_EXAPNSION)
      alp = outa(1)

      ina(1) = tt0
      call fetch_TableData( MC_THEMOEXP, gausses(1)%pMaterial%dict, outa(:), ierr, ina )
      if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_EXAPNSION)
      alp0 = outa(1)

      epsth=alp*(ttc-ref_temp)-alp0*(tt0-ref_temp)
    else
      call fetch_TableData( MC_ISOELASTIC, gausses(1)%pMaterial%dict, outa, ierr )
      if( ierr ) outa(1) = gausses(1)%pMaterial%variables(M_YOUNGS)
      young = outa(1)
    endif

    gausses(1)%strain(1) = dlog(llen/llen0)
    gausses(1)%stress(1) = young*(gausses(1)%strain(1)-epsth)

    !set stress and strain for output
    gausses(1)%strain_out(1) = gausses(1)%strain(1)
    gausses(1)%stress_out(1) = gausses(1)%stress(1)

    qf(1) = gausses(1)%stress(1)*area*llen0/llen
    qf(1:3) = -qf(1)*direc
    qf(4:6) = -qf(1:3)

  end subroutine UPDATE_C1

  !----------------------------------------------------------------------*
  subroutine NodalStress_C1(ETYPE,NN,gausses,ndstrain,ndstress)
    !----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    use mMechGauss
    integer(kind=kint), intent(in) :: ETYPE,NN
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: ndstrain(NN,6)
    real(kind=kreal), intent(out)  :: ndstress(NN,6)

    ndstrain(1,1:6) = gausses(1)%strain(1:6)
    ndstress(1,1:6) = gausses(1)%stress(1:6)
    ndstrain(2,1:6) = gausses(1)%strain(1:6)
    ndstress(2,1:6) = gausses(1)%stress(1:6)

  end subroutine
  !
  !
  !
  !----------------------------------------------------------------------*
  subroutine ElementStress_C1(ETYPE,gausses,strain,stress)
    !----------------------------------------------------------------------*
    !
    ! Calculate Strain and Stress increment of solid elements
    !
    use mMechGauss
    integer(kind=kint), intent(in) :: ETYPE
    type(tGaussStatus), intent(in) :: gausses(:)
    real(kind=kreal), intent(out)  :: strain(6)
    real(kind=kreal), intent(out)  :: stress(6)


    strain(:) = gausses(1)%strain(1:6)
    stress(:) = gausses(1)%stress(1:6)

  end subroutine


  !----------------------------------------------------------------------*
  subroutine DL_C1(etype, nn, xx, yy, zz, rho, thick, ltype, params, &
      vect, nsize)
    !----------------------------------------------------------------------*
    !**  SET DLOAD
    !   GRAV LTYPE=4  :GRAVITY FORCE
    ! I/F VARIABLES
    integer(kind = kint), intent(in)  :: etype, nn
    real(kind = kreal), intent(in)    :: xx(:), yy(:), zz(:)
    real(kind = kreal), intent(in)    :: params(0:6)
    real(kind = kreal), intent(inout) :: vect(:)
    real(kind = kreal) :: rho, thick
    integer(kind = kint) :: ltype, nsize, surtype
    ! LOCAL VARIABLES
    integer(kind = kint) :: ndof = 3
    integer(kind = kint) :: ivol, isuf, nsur, i
    integer(kind = kint) :: nod(nn)
    real(kind = kreal) :: vx, vy, vz, val, a, AA
    !--------------------------------------------------------------------
    val = params(0)
    !--------------------------------------------------------------

    ivol = 0
    isuf = 0
    if( ltype .LT. 10 ) then
      ivol = 1
    else if( ltype .GE. 10 ) then
      isuf = 1
    end if

    !--------------------------------------------------------------------
    nsize = nn*ndof
    !--------------------------------------------------------------------
    vect(1:nsize) = 0.0D0

    ! Volume force
    if( ivol .EQ. 1 ) then
      if( ltype .EQ. 4 ) then
        AA = dsqrt( ( xx(2)-xx(1) )*( xx(2)-xx(1) )   &
          +( yy(2)-yy(1) )*( yy(2)-yy(1) )   &
          +( zz(2)-zz(1) )*( zz(2)-zz(1) ) )

        a = thick
        vx = params(1)
        vy = params(2)
        vz = params(3)
        vx = vx/dsqrt( params(1)**2+params(2)**2+params(3)**2 )
        vy = vy/dsqrt( params(1)**2+params(2)**2+params(3)**2 )
        vz = vz/dsqrt( params(1)**2+params(2)**2+params(3)**2 )

        do i = 1, 2
          vect(3*i-2) = val*rho*a*0.5D0*AA*vx
          vect(3*i-1) = val*rho*a*0.5D0*AA*vy
          vect(3*i  ) = val*rho*a*0.5D0*AA*vz
        end do

      end if
    end if
    !--------------------------------------------------------------------

    return
  end subroutine

  !
  !
  !----------------------------------------------------------------------*
  subroutine truss_diag_modify(hecMAT,hecMESH)
    !----------------------------------------------------------------------*
    !
    use hecmw
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: itype, is, iE, ic_type, icel, jS, j, n

    do itype = 1, hecMESH%n_elem_type
      ic_type = hecMESH%elem_type_item(itype)
      if(ic_type == 301)then
        is = hecMESH%elem_type_index(itype-1) + 1
        iE = hecMESH%elem_type_index(itype  )
        do icel = is, iE
          jS = hecMESH%elem_node_index(icel-1)
          do j=1,2
            n = hecMESH%elem_node_item(jS+j)
            if( hecMAT%D(9*n-8) == 0.0d0)then
              hecMAT%D(9*n-8) = 1.0d0
              !call search_diag_modify(n,1,hecMAT,hecMESH)
            endif
            if( hecMAT%D(9*n-4) == 0.0d0)then
              hecMAT%D(9*n-4) = 1.0d0
              !call search_diag_modify(n,2,hecMAT,hecMESH)
            endif
            if( hecMAT%D(9*n  ) == 0.0d0)then
              hecMAT%D(9*n  ) = 1.0d0
              !call search_diag_modify(n,3,hecMAT,hecMESH)
            endif
          enddo
        enddo
      endif
    enddo

  end subroutine

  !
  !
  !----------------------------------------------------------------------*
  subroutine search_diag_modify(n,nn,hecMAT,hecMESH)
    !----------------------------------------------------------------------*
    !
    use hecmw
    type (hecmwST_matrix)     :: hecMAT
    type (hecmwST_local_mesh) :: hecMESH
    integer :: n, nn, is, iE, i, j, in
    integer :: flagl, flagu, flagb, a

    if(nn == 1)then
      a = 0
      is = hecMAT%IndexL(n-1)+1
      iE = hecMAT%IndexL(n  )
      do i=is,iE
        if(hecMAT%AL(9*i-8) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-7) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-6) /= 0.0d0) a = 1
      enddo
      is = hecMAT%IndexU(n-1)+1
      iE = hecMAT%IndexU(n  )
      do i=is,iE
        if(hecMAT%AU(9*i-8) /= 0.0d0) a = 1
        if(hecMAT%AU(9*i-7) /= 0.0d0) a = 1
        if(hecMAT%AU(9*i-6) /= 0.0d0) a = 1
      enddo
      if(hecMAT%D(9*n-7) /= 0.0d0) a = 1
      if(hecMAT%D(9*n-6) /= 0.0d0) a = 1
      if(a == 0)then
        hecMAT%D(9*n-8) = 1.0d0
        !write(*,"(a,i,a,i,a)")"### FIX DIAGONAL n:",n,", ID:",hecMESH%global_node_ID(n),", dof:1"
      endif
    endif
    if(nn == 2)then
      a = 0
      is = hecMAT%IndexL(n-1)+1
      iE = hecMAT%IndexL(n  )
      do i=is,iE
        if(hecMAT%AL(9*i-5) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-4) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-3) /= 0.0d0) a = 1
      enddo
      is = hecMAT%IndexU(n-1)+1
      iE = hecMAT%IndexU(n  )
      do i=is,iE
        if(hecMAT%AU(9*i-5) /= 0.0d0) a = 1
        if(hecMAT%AU(9*i-4) /= 0.0d0) a = 1
        if(hecMAT%AU(9*i-3) /= 0.0d0) a = 1
      enddo
      if(hecMAT%D(9*n-5) /= 0.0d0) a = 1
      if(hecMAT%D(9*n-3) /= 0.0d0) a = 1
      if(a == 0)then
        hecMAT%D(9*n-4) = 1.0d0
        !write(*,"(a,i,a,i,a)")"### FIX DIAGONAL n:",n,", ID:",hecMESH%global_node_ID(n),", dof:2"
      endif
    endif
    if(nn == 3)then
      a = 0
      is = hecMAT%IndexL(n-1)+1
      iE = hecMAT%IndexL(n  )
      do i=is,iE
        if(hecMAT%AL(9*i-2) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-1) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i  ) /= 0.0d0) a = 1
      enddo
      is = hecMAT%IndexU(n-1)+1
      iE = hecMAT%IndexU(n  )
      do i=is,iE
        if(hecMAT%AU(9*i-2) /= 0.0d0) a = 1
        if(hecMAT%AU(9*i-1) /= 0.0d0) a = 1
        if(hecMAT%AU(9*i  ) /= 0.0d0) a = 1
      enddo
      if(hecMAT%D(9*n-2) /= 0.0d0) a = 1
      if(hecMAT%D(9*n-1) /= 0.0d0) a = 1
      if(a == 0)then
        hecMAT%D(9*n  ) = 1.0d0
        !write(*,"(a,i,a,i,a)")"### FIX DIAGONAL n:",n,", ID:",hecMESH%global_node_ID(n),", dof:3"
      endif
    endif

  end subroutine


end module m_static_LIB_1d
