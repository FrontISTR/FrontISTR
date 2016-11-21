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
  SUBROUTINE STF_C1( etype,nn,ecoord,area,gausses,stiff, u ,temperature )
    USE mMechGauss
    INTEGER(kind=kint), INTENT(IN)  :: etype               !< element type
    INTEGER(kind=kint), INTENT(IN)  :: nn                  !< number of elemental nodes
    REAL(kind=kreal),   INTENT(IN)  :: ecoord(3,nn)        !< coordinates of elemental nodes
    real(kind=kreal),   INTENT(IN)  :: area                !< section area
    TYPE(tGaussStatus), INTENT(IN)  :: gausses(:)          !< status of qudrature points
    REAL(kind=kreal),   INTENT(OUT) :: stiff(:,:)          !< stiff matrix
    REAL(kind=kreal),   INTENT(IN), optional :: u(:,:)     !< nodal displacemwent
    REAL(kind=kreal),   INTENT(IN), optional :: temperature(nn)     !< temperature

    REAL(kind=kreal) DET,WG, llen, llen0, elem(3,nn)
    logical :: ierr
    REAL(kind=kreal) ina(1), outa(1), direc(3), direc0(3), coeff, strain
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
  SUBROUTINE UPDATE_C1( etype, nn, ecoord, area, u, du, qf ,gausses, TT, T0 )
!---------------------------------------------------------------------*
    use m_fstr
    use mMechGauss
! I/F VARIAVLES
    integer(kind=kint), INTENT(IN)     :: etype           !< \param [in] element type
    integer(kind=kint), INTENT(IN)     :: nn              !< \param [in] number of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: ecoord(3,nn)    !< \param [in] coordinates of elemental nodes
    real(kind=kreal),   INTENT(IN)     :: area            !< section area
    real(kind=kreal),   INTENT(IN)     :: u(3,nn)         !< \param [in] nodal dislplacements
    real(kind=kreal),   INTENT(IN)     :: du(3,nn)       !< \param [in] nodal displacement ( solutions of solver )
    real(kind=kreal),   INTENT(OUT)    :: qf(nn*3)        !< \param [out] Internal Force
    type(tGaussStatus), INTENT(INOUT)  :: gausses(:)      !< \param [out] status of qudrature points
    REAL(kind=kreal),   INTENT(IN), optional :: TT(nn)    !< current temperature
    REAL(kind=kreal),   INTENT(IN), optional :: T0(nn)    !< reference temperature

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

    qf(1) = gausses(1)%stress(1)*area*llen0/llen
    qf(1:3) = -qf(1)*direc
    qf(4:6) = -qf(1:3)

  end subroutine UPDATE_C1

 !
!> Update strain and stress inside element
!---------------------------------------------------------------------*
  SUBROUTINE UPDATEST_C1( etype, nn, XX, YY, ZZ, TT, T0, area, EDISP, gausses )
!---------------------------------------------------------------------*
    use m_fstr
    USE mMechGauss
! I/F VARIAVLES
    integer(kind=kint), INTENT(IN)     :: etype           !< \param [in] element type
    integer(kind=kint), INTENT(IN)     :: nn              !< \param [in] number of elemental nodes
    REAL(kind=kreal), INTENT(IN)       :: XX(NN),YY(NN),ZZ(NN)
    real(kind=kreal),   INTENT(IN)     :: area            !< section area
    REAL(kind=kreal), INTENT(IN)       :: TT(NN),T0(NN),EDISP(NN*3)
    type(tGaussStatus), INTENT(INOUT)  :: gausses(:)      !< \param [out] status of qudrature points


! LCOAL VARIAVLES
    real(kind=kreal)   :: direc(3), direc0(3)
    real(kind=kreal)   :: llen, llen0, ina(1), outa(1)
    real(kind=kreal)   :: elem(3,nn), ecoord(3,nn)
    logical            :: ierr
    real(kind=kreal)   :: young
    real(kind=kreal)   :: ttc, tt0, alp, alp0, epsth

    ecoord(1,:) = XX(:)
    ecoord(2,:) = YY(:)
    ecoord(3,:) = ZZ(:)
    elem(:,1) = ecoord(:,1) + edisp(1:3)
    elem(:,2) = ecoord(:,2) + edisp(4:6)

    direc = elem(:,2)-elem(:,1)
    llen = dsqrt( dot_product(direc, direc) )
    direc0 = ecoord(:,2)-ecoord(:,1)
    llen0 = dsqrt( dot_product(direc0, direc0) )

    epsth = 0.d0
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

    gausses(1)%strain(1) = (llen-llen0)/llen0
    gausses(1)%stress(1) = young*(gausses(1)%strain(1)-epsth)

  end subroutine UPDATEST_C1


!----------------------------------------------------------------------*
   SUBROUTINE NodalStress_C1(ETYPE,NN,gausses,ndstrain,ndstress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE,NN
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: ndstrain(NN,6)
   REAL(kind=kreal), INTENT(OUT)  :: ndstress(NN,6)

     ndstrain(1,1:6) = gausses(1)%strain(1:6)
     ndstress(1,1:6) = gausses(1)%stress(1:6)
     ndstrain(2,1:6) = gausses(1)%strain(1:6)
     ndstress(2,1:6) = gausses(1)%stress(1:6)

   END SUBROUTINE
!
!
!
!----------------------------------------------------------------------*
   SUBROUTINE ElementStress_C1(ETYPE,gausses,strain,stress)
!----------------------------------------------------------------------*
!
! Calculate Strain and Stress increment of solid elements
!
   USE mMechGauss
   INTEGER(kind=kint), INTENT(IN) :: ETYPE
   TYPE(tGaussStatus), INTENT(IN) :: gausses(:)
   REAL(kind=kreal), INTENT(OUT)  :: strain(6)
   REAL(kind=kreal), INTENT(OUT)  :: stress(6)


   strain(:) = gausses(1)%strain(1:6)
   stress(:) = gausses(1)%stress(1:6)

   END SUBROUTINE


!----------------------------------------------------------------------*
   SUBROUTINE DL_C1(etype, nn, xx, yy, zz, rho, thick, ltype, params, &
                             vect, nsize)
!----------------------------------------------------------------------*
!**  SET DLOAD
!   GRAV LTYPE=4  :GRAVITY FORCE
! I/F VARIABLES
      INTEGER(KIND = kint), INTENT(IN)  :: etype, nn
      REAL(KIND = kreal), INTENT(IN)    :: xx(:), yy(:), zz(:)
      REAL(KIND = kreal), INTENT(IN)    :: params(0:6)
      REAL(KIND = kreal), INTENT(INOUT) :: vect(:)
      REAL(KIND = kreal) :: rho, thick
      INTEGER(KIND = kint) :: ltype, nsize, surtype
! LOCAL VARIABLES
      INTEGER(KIND = kint) :: ndof = 3
      INTEGER(KIND = kint) :: ivol, isuf, nsur, i
      INTEGER(KIND = kint) :: nod(nn)
      REAL(KIND = kreal) :: vx, vy, vz, val, a, AA
!--------------------------------------------------------------------
      val = params(0)
      !--------------------------------------------------------------

      ivol = 0
      isuf = 0
      IF( ltype .LT. 10 ) THEN
       ivol = 1
      ELSE IF( ltype .GE. 10 ) THEN
       isuf = 1
      END IF

!--------------------------------------------------------------------
      nsize = nn*ndof
!--------------------------------------------------------------------
      vect(1:nsize) = 0.0D0

      ! Volume force
      IF( ivol .EQ. 1 ) THEN
       IF( ltype .EQ. 4 ) THEN
        AA = DSQRT( ( xx(2)-xx(1) )*( xx(2)-xx(1) )   &
                   +( yy(2)-yy(1) )*( yy(2)-yy(1) )   &
                   +( zz(2)-zz(1) )*( zz(2)-zz(1) ) )

        a = thick
        vx = params(1)
        vy = params(2)
        vz = params(3)
        vx = vx/DSQRT( params(1)**2+params(2)**2+params(3)**2 )
        vy = vy/DSQRT( params(1)**2+params(2)**2+params(3)**2 )
        vz = vz/DSQRT( params(1)**2+params(2)**2+params(3)**2 )

        DO i = 1, 2
         vect(3*i-2) = val*rho*a*0.5D0*AA*vx
         vect(3*i-1) = val*rho*a*0.5D0*AA*vy
         vect(3*i  ) = val*rho*a*0.5D0*AA*vz
        END DO

       END IF
      END IF
!--------------------------------------------------------------------

      RETURN
   END SUBROUTINE

!
!
!----------------------------------------------------------------------*
   SUBROUTINE truss_diag_modify(hecMAT,hecMESH)
!----------------------------------------------------------------------*
!
   use hecmw_matrix_misc
   type (hecmwST_matrix)     :: hecMAT
   type (hecmwST_local_mesh) :: hecMESH
   integer(kind=kint) :: itype, iS, iE, ic_type, icel, jS, j, n

   do itype = 1, hecMESH%n_elem_type
    ic_type = hecMESH%elem_type_item(itype)
    if(ic_type == 301)then
      iS = hecMESH%elem_type_index(itype-1) + 1
      iE = hecMESH%elem_type_index(itype  )
      do icel = iS, iE
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

   END SUBROUTINE

!
!
!----------------------------------------------------------------------*
   SUBROUTINE search_diag_modify(n,nn,hecMAT,hecMESH)
!----------------------------------------------------------------------*
!
   use hecmw_matrix_misc
   type (hecmwST_matrix)     :: hecMAT
   type (hecmwST_local_mesh) :: hecMESH
   integer :: n, nn, iS, iE, i, j, in
   integer :: flagl, flagu, flagb, a

    if(nn == 1)then
      a = 0
      iS = hecMAT%IndexL(n-1)+1
      iE = hecMAT%IndexL(n  )
      do i=iS,iE
        if(hecMAT%AL(9*i-8) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-7) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-6) /= 0.0d0) a = 1
      enddo
      iS = hecMAT%IndexU(n-1)+1
      iE = hecMAT%IndexU(n  )
      do i=iS,iE
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
      iS = hecMAT%IndexL(n-1)+1
      iE = hecMAT%IndexL(n  )
      do i=iS,iE
        if(hecMAT%AL(9*i-5) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-4) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-3) /= 0.0d0) a = 1
      enddo
      iS = hecMAT%IndexU(n-1)+1
      iE = hecMAT%IndexU(n  )
      do i=iS,iE
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
      iS = hecMAT%IndexL(n-1)+1
      iE = hecMAT%IndexL(n  )
      do i=iS,iE
        if(hecMAT%AL(9*i-2) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i-1) /= 0.0d0) a = 1
        if(hecMAT%AL(9*i  ) /= 0.0d0) a = 1
      enddo
      iS = hecMAT%IndexU(n-1)+1
      iE = hecMAT%IndexU(n  )
      do i=iS,iE
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

   END SUBROUTINE


end module m_static_LIB_1d
