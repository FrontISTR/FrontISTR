!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Sato(Advancesoft), X. YUAN(AdavanceSoft)    !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!======================================================================!
!
!> \brief  This module provides function to calcualte to do updates
!!
!>  \author                date                  version 
!>  X.Yuan(Advancesoft)    2009/08/28        original
!>  X.Yuan                 2013/03/18        consider anisotropic materials
!
!======================================================================!
module m_fstr_Update
  implicit none

  contains

!=====================================================================*
!  UPDATE_C3
!>  \brief 変位／応力・ひずみ／内力のアップデート
!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \version    0.00
!!
!> \if AS
!>    \par  サブルーチン構成
!>    -# 変位の更新					\f$ u_{n+1}^{(k)} = u_{n+1}^{(k-1)} + \delta u^{(k)} \f$
!>    -# ひずみ・応力の更新			\f$ \varepsilon_{n+1}^{(k)} = \varepsilon_{n+1}^{(k-1)} + \delta \varepsilon^{(k)} \f$, \f$ \sigma_{n+1}^{(k)} = \sigma_{n+1}^{(k-1)} + \delta \sigma^{(k)} \f$
!>    -# 内力（等価節点力）の計算	\f$ Q_{n+1}^{(k-1)} ( u_{n+1}^{(k-1)} ) \f$
!> \endif
!
subroutine fstr_UpdateNewton ( hecMESH, hecMAT, fstrSOLID, tincr,iter, strainEnergy)
!=====================================================================*
  use m_fstr
  use m_static_lib
!#ifdef PARA_CONTACT
  use m_fstr_para_contact
!#endif

  type (hecmwST_matrix)       :: hecMAT    !< linear equation, its right side modified here
  type (hecmwST_local_mesh)   :: hecMESH   !< mesh information
  type (fstr_solid)           :: fstrSOLID !< we need boundary conditions of curr step
  real(kind=kreal),intent(in) :: tincr     !< time increment
  integer, intent(in)         :: iter      !< NR iterations

  integer(kind=kint) :: nodLOCAL(20)
  real(kind=kreal)   :: ecoord(3,20)
  real(kind=kreal)   :: thick
  integer(kind=kint) :: ndof, itype, iS, iE, ic_type, nn, icel, iiS, i, j

  real(kind=kreal)   :: total_disp(6,20),du(6,20),ddu(6,20)
  real(kind=kreal)   :: tt(20), tt0(20), ttn(20), qf(20*6), coords(3,3)
  integer            :: ig0,  ig, ik, in, ierror, isect, ihead, cdsys_ID

  real(kind=kreal), optional :: strainEnergy

  ndof = hecMAT%NDOF
  fstrSOLID%QFORCE=0.0d0


! --------------------------------------------------------------------
!      updated
!         1. stress and strain  : ep^(k) = ep^(k-1)+dep^(k)
!                                 sgm^(k) = sgm^(k-1)+dsgm^(k)
!         2. Internal Force     : Q^(k-1) ( u^(k-1) )
! --------------------------------------------------------------------
!
! ----------------------------------------------------------------------------------
!      calculate the Strain and Stress and Internal Force ( Equivalent Nodal Force )
! ----------------------------------------------------------------------------------
!

  do itype= 1, hecMESH%n_elem_type
    iS= hecMESH%elem_type_index(itype-1) + 1
    iE= hecMESH%elem_type_index(itype  )
    ic_type= hecMESH%elem_type_item(itype)
    if (hecmw_is_etype_link(ic_type)) cycle
    nn = hecmw_get_max_node(ic_type)
    if( nn>20 ) stop "Elemental nodes>20!"

! element loop
    do icel= iS, iE
!
! ----- nodal coordinate
      iiS= hecMESH%elem_node_index(icel-1)
      do j=1,nn
        nodLOCAL(j)= hecMESH%elem_node_item (iiS+j)
        do i=1,3
          ecoord(i,j) = hecMESH%node(3*nodLOCAL(j)+i-3)
        enddo
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          if( isElastoplastic(fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype) .or. &
             fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype==NORTON ) then
            tt0(j)=fstrSOLID%last_temp( nodLOCAL(j) ) 
          else
            tt0(j) = 0.d0
            if( hecMESH%hecmw_flag_initcon==1 ) tt0(j) = hecMESH%node_init_val_item(nodLOCAL(j))
          endif
          ttn(j)=fstrSOLID%last_temp( nodLOCAL(j) )
          tt(j) = fstrSOLID%temperature( nodLOCAL(j) ) 
        endif
      enddo

! nodal displacement
        do j=1,nn
          nodLOCAL(j)= hecMESH%elem_node_item (iiS+j)
          do i=1,ndof
            ddu(i,j) = hecMAT%X(ndof*nodLOCAL(j)+i-ndof)
            du(i,j)  = fstrSOLID%dunode(ndof*nodLOCAL(j)+i-ndof)
            total_disp(i,j) = fstrSOLID%unode(ndof*nodLOCAL(j)+i-ndof) 
          enddo
        enddo
		
      cdsys_ID = fstrSOLID%elements(icel)%gausses(1)%pMaterial%cdsys_ID
      if( cdsys_ID>0 ) call get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )
!
! ===== calculate the Internal Force
      if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
      if ( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 .or. ic_type==2322 ) then
        call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:),      &
                        thick,fstrSOLID%elements(icel)%iset,                                  &
                        total_disp(1:2,1:nn), ddu(1:2,1:nn), qf(1:nn*ndof) )
						
      else if ( ic_type==301 ) then
        isect= hecMESH%section_ID(icel)
        ihead = hecMESH%section%sect_R_index(isect-1)
        thick = hecMESH%section%sect_R_item(ihead+1)
        call UPDATE_C1(  ic_type,nn,ecoord(:,1:nn), thick, total_disp(1:3,1:nn), ddu(1:3,1:nn),  &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:) )
						
      else if ( ic_type==361 ) then
        if( fstrPR%solution_type==kstSTATIC ) then
          call UpdateST_C3D8IC(ic_type,nn,ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),       &
                 tt(1:nn), tt0(1:nn),ddu(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:),coords)
        else
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          call UPDATE_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), coords,   &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), iter, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
        else
          call Update_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), coords       &
                        , qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), iter, tincr  )
        endif
        endif
        
      else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.                          &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), coords, qf(1:nn*ndof)       &
                        ,fstrSOLID%elements(icel)%gausses(:), iter, tincr, tt(1:nn), tt0(1:nn), ttn(1:nn)  )
        else
          call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), coords       &
                        , qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), iter, tincr  )
        endif


 !     else if ( ic_type==731) then
!        call UPDATE_S3(xx,yy,zz,ee,pp,thick,local_stf)
!        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
!      else if ( ic_type==741) then
!        call UPDATE_S4(xx,yy,zz,ee,pp,thick,local_stf)
!        call fstr_local_stf_restore_temp(local_stf, nn*ndof, stiffness)
      else
        write(*,*) '###ERROR### : Element type not supported for nonlinear static analysis'
        write(*,*) ' ic_type = ', ic_type
        call hecmw_abort(hecmw_comm_get_comm())
      endif
!
! ----- calculate the global internal force ( Q(u_{n+1}^{k-1}) )
        do j= 1, nn
          do i = 1, ndof
            fstrSOLID%QFORCE(ndof*(nodLOCAL(j)-1)+i)                          &
               = fstrSOLID%QFORCE(ndof*(nodLOCAL(j)-1)+i)+qf(ndof*(j-1)+i)
          enddo
        enddo
!
! ----- calculate strain energy
        if(present(strainEnergy))then
          do j= 1, nn
            do i=1, ndof
              strainEnergy=strainEnergy+0.5d0*(fstrSOLID%elements(icel)%equiForces(ndof*(j-1)+i)&
                                              +qf(ndof*(j-1)+i))*ddu(i,j)
              fstrSOLID%elements(icel)%equiForces(ndof*(j-1)+i)=qf(ndof*(j-1)+i)
            enddo
          enddo
        endif

    enddo      ! icel
  enddo        ! itype


!C
!C Update for fstrSOLID%QFORCE
!C
      if( ndof==3 ) then
        if(paraContactFlag) then
          call paraContact_update_3_R(hecMESH,fstrSOLID%QFORCE)
        else
          call hecmw_update_3_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node)
        endif
      else if( ndof==2 ) then
        call hecmw_update_2_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node)
      else if( ndof==6 ) then
        call hecmw_update_m_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node,6)
      endif

end subroutine fstr_UpdateNewton

!> Update elastiplastic status
subroutine fstr_UpdateState( hecMESH, fstrSOLID, tincr)
	use m_fstr
	use m_static_lib
	use m_ElastoPlastic
	use mCreep
	use mViscoElastic
	type (hecmwST_local_mesh)  :: hecMESH     !< hecmw mesh
	type (fstr_solid)          :: fstrSOLID   !< fstr_solid
	real(kind=kreal)           :: tincr

	integer(kind=kint) :: itype, iS, iE, ic_type, icel, ngauss, i
  
  
	if( associated( fstrSOLID%temperature ) ) then 
		do i=1, hecMESH%n_node
		fstrSOLID%last_temp(i) = fstrSOLID%temperature(i)
		end do
	endif

	do itype= 1, hecMESH%n_elem_type
		iS= hecMESH%elem_type_index(itype-1) + 1
		iE= hecMESH%elem_type_index(itype  )
		ic_type= hecMESH%elem_type_item(itype)
		if( ic_type==301 ) ic_type=111
		if (hecmw_is_etype_link(ic_type)) cycle

		ngauss = NumOfQuadPoints( ic_type )
		do icel= iS, iE
			if( isElastoplastic( fstrSOLID%elements(icel)%gausses(1)   &
				 %pMaterial%mtype ) ) then
				do i=1,ngauss
					call updateEPState( fstrSOLID%elements(icel)%gausses(i) )
				enddo
			elseif( fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype == NORTON ) then
				if( tincr>0.d0 ) then
					do i=1,ngauss
						fstrSOLID%elements(icel)%gausses(i)%ttime = fstrSOLID%elements(icel)%gausses(i)%ttime+tincr
						call updateViscoState( fstrSOLID%elements(icel)%gausses(i) )
					enddo
				endif
			elseif( isViscoelastic( fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype ) ) then
				if( tincr>0.d0 ) then 
					do i=1,ngauss
						fstrSOLID%elements(icel)%gausses(i)%ttime = fstrSOLID%elements(icel)%gausses(i)%ttime+tincr
						call updateViscoElasticState( fstrSOLID%elements(icel)%gausses(i) )
					enddo
				endif
			endif
	  
			do i=1,ngauss
				fstrSOLID%elements(icel)%gausses(i)%strain_bak = fstrSOLID%elements(icel)%gausses(i)%strain
				fstrSOLID%elements(icel)%gausses(i)%stress_bak = fstrSOLID%elements(icel)%gausses(i)%stress
			enddo
		enddo
	enddo
end subroutine fstr_UpdateState

!> Update at linear static analysis ----------------------------------------------------------------------
subroutine fstr_Update3D( hecMESH, fstrSOLID )
	use m_fstr
	use m_static_lib
	type (hecmwST_local_mesh) :: hecMESH
	type (fstr_solid)         :: fstrSOLID
!C** local variables
	integer(kind=kint) :: itype, icel, iS, iE, jS, i, j, ic_type, ig, nn, isect, ihead, iflag, mixflag
	integer(kind=kint) :: nodLOCAL(20), flag_dof
	real(kind=kreal)   :: xx(20), yy(20), zz(20), tt(20), tt0(20), edisp(60), force(60)
	real(kind=kreal)   :: ecoord(3,20), stiff(60,60)
	real(kind=kreal)   :: thick, coords(3,3)
	real(kind=kreal), allocatable :: temp(:)
	integer(kind=kint), allocatable :: id_spc(:)
	real(kind=kreal) :: a

!C
!C set temperature
!C
	allocate( temp(hecMESH%n_node) )
	temp = 0.0d0
	do i = 1, fstrSOLID%TEMP_ngrp_tot
		ig = fstrSOLID%TEMP_ngrp_ID(i)
		iS = hecMESH%node_group%grp_index(ig-1) + 1
		iE = hecMESH%node_group%grp_index(ig  )
		do j = iS, iE
			temp( hecMESH%node_group%grp_item(j) ) = fstrSOLID%TEMP_ngrp_val(i)
		enddo
	enddo
!C
!C set boundary force
!C
	fstrSOLID%QFORCE = 0.0d0
	allocate ( id_spc(hecMESH%n_node) )
	id_spc = 0
	do i = 1, fstrSOLID%BOUNDARY_ngrp_tot
		ig = fstrSOLID%BOUNDARY_ngrp_ID(i)
		iS = hecMESH%node_group%grp_index(ig-1) + 1
		iE = hecMESH%node_group%grp_index(ig  )
		do j = iS, iE
		id_spc( hecMESH%node_group%grp_item(j) ) = 1
		enddo
	enddo

!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
	do itype = 1, hecMESH%n_elem_type
		iS = hecMESH%elem_type_index(itype-1) + 1
		iE = hecMESH%elem_type_index(itype  )
		ic_type = hecMESH%elem_type_item(itype)
		if( ic_type == fe_tet10nc ) ic_type = fe_tet10n
		if(.not. (hecmw_is_etype_solid(ic_type) .or. ic_type == 781 .or. ic_type == 761 .or. ic_type == 641 )) cycle
		nn = hecmw_get_max_node( ic_type )
!C element loop
		do icel = iS, iE
			jS = hecMESH%elem_node_index(icel-1)
			do j = 1, nn
				nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
				xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
				yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
				zz(j) = hecMESH%node(3*nodLOCAL(j)  )
				tt(j) = temp( nodLOCAL(j) )
				tt0(j)= ref_temp
				ecoord(1,j) = hecMESH%node(3*nodLOCAL(j)-2)
				ecoord(2,j) = hecMESH%node(3*nodLOCAL(j)-1)
				ecoord(3,j) = hecMESH%node(3*nodLOCAL(j)  )
				edisp(3*j-2) = fstrSOLID%unode(3*nodLOCAL(j)-2)
				edisp(3*j-1) = fstrSOLID%unode(3*nodLOCAL(j)-1)
				edisp(3*j  ) = fstrSOLID%unode(3*nodLOCAL(j)  )
			enddo
!--- calculate stress and strain of gauss point
			if ( ic_type == 781 ) then   !for shell-solid mixed analysis
				isect= hecMESH%section_ID(icel)
				ihead = hecMESH%section%sect_R_index(isect-1)
				thick = hecMESH%section%sect_R_item(ihead+1)
				nn = 4;	ic_type = 741; mixflag = 1
				call STF_Shell_MITC(ic_type, nn, 6, ecoord(1:3, 1:8), fstrSOLID%elements(icel)%gausses, stiff, thick, mixflag)
				ic_type = 781; nn = 8
				CYCLE
			elseif ( ic_type == 761 ) then   !for shell-solid mixed analysis
				isect= hecMESH%section_ID(icel)
				ihead = hecMESH%section%sect_R_index(isect-1)
				thick = hecMESH%section%sect_R_item(ihead+1)
				nn = 3; ic_type = 731; mixflag = 2
				call STF_Shell_MITC(ic_type, nn, 6, ecoord(1:3, 1:8), fstrSOLID%elements(icel)%gausses, stiff, thick, mixflag)
				ic_type = 761; nn = 6
				CYCLE
			ELSEIF( ic_type == 641 ) THEN
				isect = hecMESH%section_ID(icel)
				ihead = hecMESH%section%sect_R_index(isect-1)
				CALL STF_Beam_641( ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses(:), hecMESH%section%sect_R_item(ihead+1:), stiff )
				iflag = 0
				DO j = 1, nn
					IF( id_spc( nodLOCAL(j) ) == 1 ) iflag = 1
				END DO
				if( iflag == 1 ) THEN
					IF( ic_type == 641 ) THEN
						isect = hecMESH%section_ID(icel)
						ihead = hecMESH%section%sect_R_index(isect-1)
						CALL STF_Beam_641( ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses(:), hecMESH%section%sect_R_item(ihead+1:), stiff )
					END IF
					force(1:nn*3) = MATMUL( stiff(1:nn*3,1:nn*3), edisp(1:nn*3) )
					DO j = 1, nn
						IF( id_spc( nodLOCAL(j) ) == 1 ) THEN
							fstrSOLID%QFORCE( 3*nodLOCAL(j)-2 ) = fstrSOLID%QFORCE( 3*nodLOCAL(j)-2 )+force(3*j-2)
							fstrSOLID%QFORCE( 3*nodLOCAL(j)-1 ) = fstrSOLID%QFORCE( 3*nodLOCAL(j)-1 )+force(3*j-1)
							fstrSOLID%QFORCE( 3*nodLOCAL(j)   ) = fstrSOLID%QFORCE( 3*nodLOCAL(j)   )+force(3*j  )
						END IF
					END DO
				ENDIF
				CYCLE
			elseif( ic_type == 361 ) then
				call UpdateST_C3D8IC( ic_type, nn, xx, yy, zz, tt, tt0, edisp, fstrSOLID%elements(icel)%gausses, coords )
			else if( ic_type == 301 ) then
				isect = hecMESH%section_ID(icel)
				ihead = hecMESH%section%sect_R_index(isect-1)
				thick = hecMESH%section%sect_R_item(ihead+1)
				call UpdateST_C1( ic_type, nn, xx, yy, zz, thick, edisp, fstrSOLID%elements(icel)%gausses )
			else
				call UpdateST_C3( ic_type, nn, xx, yy, zz, tt, tt0, edisp, fstrSOLID%elements(icel)%gausses, coords )
			endif
!--- calculate reaction force
			iflag = 0
			do j = 1, nn
				if( id_spc( nodLOCAL(j) ) == 1 ) iflag = 1
			enddo
			if( iflag == 1 ) then
				if( ic_type == 361 ) then
					call STF_C3D8IC( ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses, stiff )
				else if( ic_type == 301 ) then
					isect = hecMESH%section_ID(icel)
					ihead = hecMESH%section%sect_R_index(isect-1)
					thick = hecMESH%section%sect_R_item(ihead+1)
					call STF_C1( ic_type, nn, ecoord, thick, fstrSOLID%elements(icel)%gausses, stiff )
				else
					call STF_C3( ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses, stiff, 1.0d0, coords )
				endif
				force(1:nn*3) = matmul( stiff(1:nn*3,1:nn*3), edisp(1:nn*3) )
				do j = 1, nn
					if( id_spc( nodLOCAL(j) ) == 1 ) then
					fstrSOLID%QFORCE(3*nodLOCAL(j)-2) = fstrSOLID%QFORCE(3*nodLOCAL(j)-2) + force(3*j-2)
					fstrSOLID%QFORCE(3*nodLOCAL(j)-1) = fstrSOLID%QFORCE(3*nodLOCAL(j)-1) + force(3*j-1)
					fstrSOLID%QFORCE(3*nodLOCAL(j)  ) = fstrSOLID%QFORCE(3*nodLOCAL(j)  ) + force(3*j  )
					endif
				enddo
			endif
		enddo
	enddo

	call hecmw_update_3_R( hecMESH, fstrSOLID%QFORCE, hecMESH%n_node )

	deallocate( temp )
	deallocate( id_spc )
end subroutine fstr_Update3D

!> Update at linear static analysis
subroutine fstr_Update2D( hecMESH, fstrSOLID )
	use m_fstr
	use m_static_lib
	type (hecmwST_local_mesh) :: hecMESH
	type (fstr_solid)         :: fstrSOLID
!C** local variables
	integer(kind=kint) :: itype, icel, iS, iE, jS, i, j, ic_type, ig, nn, iflag
	integer(kind=kint) :: nodLOCAL(8)
	real(kind=kreal)   :: xx(8), yy(8), tt(8), tt0(8), edisp(16), force(16)
	real(kind=kreal)   :: ecoord(2,8), stiff(16,16)
	real(kind=kreal), allocatable :: temp(:)
	integer(kind=kint), allocatable :: id_spc(:)

!C
!C set temperature
!C
	allocate( temp(hecMESH%n_node) )
	temp = 0.0d0
	do i = 1, fstrSOLID%TEMP_ngrp_tot
		ig = fstrSOLID%TEMP_ngrp_ID(i)
		iS = hecMESH%node_group%grp_index(ig-1) + 1
		iE = hecMESH%node_group%grp_index(ig  )
		do j = iS, iE
			temp( hecMESH%node_group%grp_item(j) ) = fstrSOLID%TEMP_ngrp_val(i)
		enddo
	enddo
!C
!C set boundary force
!C
	fstrSOLID%QFORCE = 0.0d0
	allocate ( id_spc(hecMESH%n_node) )
	id_spc = 0
	do i = 1, fstrSOLID%BOUNDARY_ngrp_tot
		ig = fstrSOLID%BOUNDARY_ngrp_ID(i)
		iS = hecMESH%node_group%grp_index(ig-1) + 1
		iE = hecMESH%node_group%grp_index(ig  )
		do j = iS, iE
			id_spc( hecMESH%node_group%grp_item(j) ) = 1
		enddo
	enddo

!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
	do itype = 1, hecMESH%n_elem_type
		iS = hecMESH%elem_type_index(itype-1) + 1
		iE = hecMESH%elem_type_index(itype  )
		ic_type = hecMESH%elem_type_item(itype)
		if( .not. hecmw_is_etype_surface(ic_type) ) cycle
		nn = hecmw_get_max_node( ic_type )
!C element loop
		do icel = iS, iE
			jS = hecMESH%elem_node_index(icel-1)
			do j = 1, nn
				nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
				xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
				yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
				tt(j) = temp( nodLOCAL(j) )
				tt0(j)= ref_temp
				ecoord(1,j) = hecMESH%node(3*nodLOCAL(j)-2)
				ecoord(2,j) = hecMESH%node(3*nodLOCAL(j)-1)
				edisp(2*j-1) = fstrSOLID%unode(2*nodLOCAL(j)-1)
				edisp(2*j  ) = fstrSOLID%unode(2*nodLOCAL(j)  )
			enddo
!--- calculate stress and strain of gauss points
			call UpdateST_C2( ic_type, nn, xx, yy, tt, tt0, 1.0d0, fstrSOLID%elements(icel)%iset, &
                        edisp, fstrSOLID%elements(icel)%gausses )
!--- calculate reaction force
			iflag = 0
			do j = 1, nn
				if( id_spc( nodLOCAL(j) ) == 1 ) iflag = 1
			enddo
			if( iflag == 1 ) then
				call STF_C2( ic_type, nn, ecoord, fstrSOLID%elements(icel)%gausses, 1.0d0, &
						 stiff, fstrSOLID%elements(icel)%iset )
				force(1:nn*2) = matmul( stiff(1:nn*2,1:nn*2), edisp(1:nn*2) )
				do j = 1, nn
					if( id_spc( nodLOCAL(j) ) == 1 ) then
						fstrSOLID%QFORCE(2*nodLOCAL(j)-1) = fstrSOLID%QFORCE(2*nodLOCAL(j)-1) + force(2*j-1)
						fstrSOLID%QFORCE(2*nodLOCAL(j)  ) = fstrSOLID%QFORCE(2*nodLOCAL(j)  ) + force(2*j  )
					endif
				enddo
			endif
		enddo
	enddo

  call hecmw_update_2_R( hecMESH, fstrSOLID%QFORCE, hecMESH%n_node )

  deallocate( temp )
  deallocate( id_spc )
end subroutine fstr_Update2D

!> Update at linear static analysis
subroutine fstr_Update6D( hecMESH, fstrSOLID )
	use m_fstr
	use m_static_lib
	type (hecmwST_local_mesh) :: hecMESH
	type (fstr_solid)         :: fstrSOLID
!C** local variables
	integer(kind=kint) :: itype, icel, iS, iE, jS, i, j, ic_type, ig, nn, isect, ihead, iflag, mixflag
	integer(kind=kint) :: nodLOCAL(9)
	real(kind=kreal)   :: ecoord(3,9), edisp(54), force(54), stiff(54,54)
	real(kind=kreal)   :: thick, logr
	integer(kind=kint), allocatable :: id_spc(:)

!C
!C set boundary force
!C
	fstrSOLID%QFORCE = 0.0d0
	allocate ( id_spc(hecMESH%n_node) )
	id_spc = 0
	do i = 1, fstrSOLID%BOUNDARY_ngrp_tot
		ig = fstrSOLID%BOUNDARY_ngrp_ID(i)
		iS = hecMESH%node_group%grp_index(ig-1) + 1
		iE = hecMESH%node_group%grp_index(ig  )
		do j = iS, iE
			id_spc( hecMESH%node_group%grp_item(j) ) = 1
		enddo
	enddo

!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
	do itype = 1, hecMESH%n_elem_type
		iS = hecMESH%elem_type_index(itype-1) + 1
		iE = hecMESH%elem_type_index(itype  )
		ic_type = hecMESH%elem_type_item(itype)
		if( .not. hecmw_is_etype_shell(ic_type) ) cycle
		nn = hecmw_get_max_node( ic_type )
!C element loop
		do icel = iS, iE
			jS = hecMESH%elem_node_index(icel-1)
			do j = 1, nn
				nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
				ecoord(1,j) = hecMESH%node(3*nodLOCAL(j)-2)
				ecoord(2,j) = hecMESH%node(3*nodLOCAL(j)-1)
				ecoord(3,j) = hecMESH%node(3*nodLOCAL(j)  )
				edisp(6*j-5) = fstrSOLID%unode(6*nodLOCAL(j)-5)
				edisp(6*j-4) = fstrSOLID%unode(6*nodLOCAL(j)-4)
				edisp(6*j-3) = fstrSOLID%unode(6*nodLOCAL(j)-3)
				edisp(6*j-2) = fstrSOLID%unode(6*nodLOCAL(j)-2)
				edisp(6*j-1) = fstrSOLID%unode(6*nodLOCAL(j)-1)
				edisp(6*j  ) = fstrSOLID%unode(6*nodLOCAL(j)  )
			enddo
			isect = hecMESH%section_ID(icel)
			ihead = hecMESH%section%sect_R_index(isect-1)
			thick = hecMESH%section%sect_R_item(ihead+1)
!--- calculate reaction force
			if( ic_type == 731 .or. ic_type == 741 .or. ic_type == 743 ) then
				iflag = 0
				do j = 1, nn
				if( id_spc( nodLOCAL(j) ) == 1 ) iflag = 1
				enddo
				if( iflag == 1 ) then
					mixflag = 0
					call STF_Shell_MITC( ic_type, nn, 6, ecoord(1:3,1:nn), fstrSOLID%elements(icel)%gausses, &
					stiff(1:nn*6,1:nn*6), thick, mixflag )
					force(1:nn*6) = matmul( stiff(1:nn*6,1:nn*6), edisp(1:nn*6) )
					do j = 1, nn
						if( id_spc( nodLOCAL(j) ) == 1 ) then
						fstrSOLID%QFORCE(6*nodLOCAL(j)-5) = fstrSOLID%QFORCE(6*nodLOCAL(j)-5) + force(6*j-5)
						fstrSOLID%QFORCE(6*nodLOCAL(j)-4) = fstrSOLID%QFORCE(6*nodLOCAL(j)-4) + force(6*j-4)
						fstrSOLID%QFORCE(6*nodLOCAL(j)-3) = fstrSOLID%QFORCE(6*nodLOCAL(j)-3) + force(6*j-3)
						fstrSOLID%QFORCE(6*nodLOCAL(j)-2) = fstrSOLID%QFORCE(6*nodLOCAL(j)-2) + force(6*j-2)
						fstrSOLID%QFORCE(6*nodLOCAL(j)-1) = fstrSOLID%QFORCE(6*nodLOCAL(j)-1) + force(6*j-1)
						fstrSOLID%QFORCE(6*nodLOCAL(j)  ) = fstrSOLID%QFORCE(6*nodLOCAL(j)  ) + force(6*j  )
						endif
					enddo
				endif
			endif
		enddo
	enddo

	call hecmw_update_m_R( hecMESH, fstrSOLID%QFORCE, hecMESH%n_node, 6 )

	deallocate( id_spc )
end subroutine fstr_Update6D

end module m_fstr_Update
