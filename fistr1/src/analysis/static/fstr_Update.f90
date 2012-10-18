!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
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
!> \brief  This module provides functions to do update.
!!
!>  \author     K. Sato(Advancesoft), X. YUAN(AdavanceSoft)
!>  \date       2009/08/28
!>  \version    0.00
!!
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
  real(kind=kreal)   :: tt(20), tt0(20), qf(20*6)
  integer            :: ig0,  ig, ik, in, ierror

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
!
! ===== calculate the Internal Force
      if(  getSpaceDimension( ic_type )==2 ) thick =1.d0
      if ( ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 .or. ic_type==2322 ) then
        call UPDATE_C2( ic_type,nn,ecoord(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:),      &
                        thick,fstrSOLID%elements(icel)%iset,                                  &
                        total_disp(1:2,1:nn), ddu(1:2,1:nn), qf(1:nn*ndof) )
						
      else if ( ic_type==301 ) then
        thick = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_THICK)
        call UPDATE_C1(  ic_type,nn,ecoord(:,1:nn), thick, total_disp(1:3,1:nn), ddu(1:3,1:nn),  &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:) )
						
      else if ( ic_type==361 ) then
        if( fstrPR%solution_type==kstSTATIC ) then
          call UpdateST_C3D8IC(ic_type,nn,ecoord(1,1:nn),ecoord(2,1:nn),ecoord(3,1:nn),       &
                 tt(1:nn), tt0(1:nn),ddu(1:3,1:nn),fstrSOLID%elements(icel)%gausses(:))
        else
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          call UPDATE_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn),   &
            qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), iter, tincr, tt(1:nn), tt0(1:nn)  )
        else
          call Update_C3D8Bbar( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn)       &
                        , qf(1:nn*ndof),fstrSOLID%elements(icel)%gausses(:), iter, tincr  )
        endif
        endif
        
      else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.                          &
               ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
        if( fstrSOLID%TEMP_ngrp_tot > 0 .or. fstrSOLID%TEMP_irres >0 ) then
          call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn), qf(1:nn*ndof)       &
                        ,fstrSOLID%elements(icel)%gausses(:), iter, tincr, tt(1:nn), tt0(1:nn)  )
        else
          call UPDATE_C3( ic_type,nn,ecoord(:,1:nn), total_disp(1:3,1:nn), du(1:3,1:nn)       &
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
        call hecmw_update_3_R(hecMESH,fstrSOLID%QFORCE,hecMESH%n_node)
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
      elseif( fstrSOLID%elements(icel)%gausses(1)%pMaterial%mtype == VISCOELASTIC ) then
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
end subroutine

end module m_fstr_Update
