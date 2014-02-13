!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Dynamic Transit Analysis                          !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                       Tomotaka Ogasawara (Univ. of Tokyo)            !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!C================================================================C
!> \brief This module contains function to set boundary condition of external load in dynamic analysis
!C================================================================C

module m_dynamic_mat_ass_load
    contains
!C
!C***
!> This function sets boundary condition of external load
!C***
!C
    subroutine DYNAMIC_MAT_ASS_LOAD(hecMESH, hecMAT, fstrSOLID, fstrDYNAMIC )

      use m_fstr
      use m_static_lib
      use m_fstr_precheck
      use m_table_dyn

      implicit none
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid        ) :: fstrSOLID
      type ( fstr_dynamic     ) :: fstrDYNAMIC

      real(kind=kreal) :: xx(20), yy(20), zz(20)
      real(kind=kreal) :: params(0:6)
      real(kind=kreal) :: vect(60)
      integer(kind=kint) :: iwk(60)
      integer(kind=kint) :: nodLocal(20)
      real(kind=kreal) :: tt(20), tt0(20), coords(3,3)
      real(kind=kreal),pointer:: temp(:)
      integer(kind=kint) :: ndof, ig0, ig, ityp, ltype, iS0, iE0, ik, in, i, j
      integer(kind=kint) :: icel, ic_type, nn, iS, isect, id, iset, nsize
      integer(kind=kint) :: itype, iE, cdsys_ID
      real(kind=kreal) :: val, rho, thick, pa1
      logical :: fg_surf
 
      integer(kind=kint) :: flag_u, ierror 
      real(kind=kreal) :: f_t


      ndof = hecMAT%NDOF
      hecMAT%B(:) = 0.0d0
!C
!C CLOAD
!C 
      do ig0= 1, fstrSOLID%CLOAD_ngrp_tot
        ig= fstrSOLID%CLOAD_ngrp_ID(ig0)
        ityp= fstrSOLID%CLOAD_ngrp_DOF(ig0)
        val= fstrSOLID%CLOAD_ngrp_val(ig0)

        flag_u = 0
        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
        val = val * f_t

        iS0= hecMESH%node_group%grp_index(ig-1) + 1
        iE0= hecMESH%node_group%grp_index(ig  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          hecMAT%B(ndof*(in-1)+ityp)=hecMAT%B(ndof*(in-1)+ityp)+val
        enddo
      enddo
!C
!C DLOAD
!C 
      do ig0= 1, fstrSOLID%DLOAD_ngrp_tot
        ig= fstrSOLID%DLOAD_ngrp_ID(ig0)
        ltype= fstrSOLID%DLOAD_ngrp_LID(ig0)
        do i=0,6
          params(i)= fstrSOLID%DLOAD_ngrp_params(i,ig0)
        enddo
!C START & END
        fg_surf = (ltype == 100)
        if( fg_surf ) then ! surface group
          iS0= hecMESH%surf_group%grp_index(ig-1) + 1
          iE0= hecMESH%surf_group%grp_index(ig  )
        else ! element group
          iS0= hecMESH%elem_group%grp_index(ig-1) + 1
          iE0= hecMESH%elem_group%grp_index(ig  )
        endif
        do ik= iS0, iE0
          if( fg_surf ) then ! surface group
            ltype  = hecMESH%surf_group%grp_item(2*ik)*10
            icel   = hecMESH%surf_group%grp_item(2*ik-1)
            ic_type= hecMESH%elem_type(icel)
          else ! element group
            icel   = hecMESH%elem_group%grp_item(ik)
            ic_type= hecMESH%elem_type(icel)
          endif
!C** Create local stiffness
          nn = hecmw_get_max_node(ic_type)
!C** node ID
          iS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (iS+j)
!C** nodal coordinate
            xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
            zz(j)=hecMESH%node(3*nodLOCAL(j)  )
!C** create iwk array ***
            do i=1,ndof
              iwk(ndof*(j-1)+i)=ndof*(nodLOCAL(j)-1)+i
            enddo
          enddo  
!C** section  ID
          isect= hecMESH%section_ID(icel)
!C** Get Properties
          rho = fstrSOLID%elements(icel)%gausses(1)%pMaterial%variables(M_DENSITY)
          call fstr_get_thickness(hecMESH,isect,thick)
!C** Section Data
          if (ndof .eq. 2) then
            id=hecMESH%section%sect_opt(isect)
            if( id.eq.0 ) then
              iset=1
            else if( id.eq.1) then
              iset=0
            else if( id.eq.2) then
              iset=2
            endif
            pa1=1.0
          endif
!C** Create local stiffness
          if (ic_type==241 .or.ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
            call DL_C2(ic_type,nn,xx(1:nn),yy(1:nn),rho,pa1,ltype,params,vect(1:nn*ndof),nsize,iset)
!
          else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.     &
                   ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
            call DL_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),rho,ltype,params,vect(1:nn*ndof),nsize)
!
          else if( ( ic_type==741 ) .or. ( ic_type==743 ) .or. ( ic_type==731 ) ) then
            call DL_Shell(ic_type, nn, ndof, xx, yy, zz, rho, thick, ltype, params, vect, nsize)
          endif
!
!!!!!!  time history

          flag_u = 10
          call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, f_t, flag_u)
          do j=1,nsize
             vect(j) = vect(j) * f_t
          end do
!
!C** Add vector
          do j=1,nsize
            hecMAT%B( iwk(j) )=hecMAT%B( iwk(j) )+vect(j)
          enddo
        enddo
      enddo
!C
!C THERMAL LOAD USING TEMPERATURE
!C 
!C Set Temperature
!C 
      if( fstrSOLID%TEMP_ngrp_tot > 0 ) then

        if( hecMESH%my_rank .eq. 0 ) then
          write(imsg,*) 'stop: THERMAL LOAD is not yet available in dynamic analysis!'
        end if
        call hecmw_abort( hecmw_comm_get_comm())
!
        allocate ( temp(hecMESH%n_node) )
        temp=0
        do ig0= 1, fstrSOLID%TEMP_ngrp_tot
          ig= fstrSOLID%TEMP_ngrp_ID(ig0)
          val=fstrSOLID%TEMP_ngrp_val(ig0)
!C START & END
          iS0= hecMESH%node_group%grp_index(ig-1) + 1
          iE0= hecMESH%node_group%grp_index(ig  )
          do ik= iS0, iE0
            in   = hecMESH%node_group%grp_item(ik)
            temp( in ) = val
          enddo
        enddo                        
!C
!C +-------------------------------+
!C | ELEMENT-by-ELEMENT ASSEMBLING | 
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
!C===
        do itype= 1, hecMESH%n_elem_type
          iS= hecMESH%elem_type_index(itype-1) + 1
          iE= hecMESH%elem_type_index(itype  )
          ic_type= hecMESH%elem_type_item(itype)
          if (hecmw_is_etype_link(ic_type)) cycle
!C** Set number of nodes		
          nn = hecmw_get_max_node(ic_type)
!C element loop
          do icel= iS, iE
!C** node ID
            iS= hecMESH%elem_node_index(icel-1)
            do j=1,nn
              nodLOCAL(j)=hecMESH%elem_node_item(iS+j)
!C** nodal coordinate
              xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
              yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
              zz(j)=hecMESH%node(3*nodLOCAL(j)  )
              tt(j)=temp( nodLOCAL(j) ) 
              tt0(j)=ref_temp
!C** create iwk array ***
              do i=1,ndof
                iwk(ndof*(j-1)+i)=ndof*(nodLOCAL(j)-1)+i
              enddo
            enddo
!C** section  ID
            isect= hecMESH%section_ID(icel)
            cdsys_ID = fstrSOLID%elements(icel)%gausses(1)%pMaterial%cdsys_ID
            call get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )

!C** Section Data
            if (ndof .eq. 2) then
              id=hecMESH%section%sect_opt(isect)
              if( id.eq.0 ) then
                iset=1
              else if( id.eq.1) then
                iset=0
              else if( id.eq.2) then
                iset=2
              endif
              pa1=1.0
            endif
!C** Create local stiffness
            if (ic_type==241 .or. ic_type==242 .or. ic_type==231 .or. ic_type==232 ) then
              call TLOAD_C2(ic_type,nn,xx(1:nn),yy(1:nn),tt(1:nn),tt0(1:nn),  &
			         fstrSOLID%elements(icel)%gausses,pa1,iset,vect(1:nn*2) )
!
            else if (ic_type==341 .or. ic_type==351 .or. ic_type==361 .or.  &
                     ic_type==342 .or. ic_type==352 .or. ic_type==362 ) then
              call TLOAD_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn), &
                            fstrSOLID%elements(icel)%gausses,vect(1:nn*ndof), coords )
!
            else if ( ic_type.eq.741 .or. ic_type.eq.743 .or. ic_type.eq.731 ) then
              if( myrank .eq. 0 ) then
                WRITE(IMSG,*) '*------------------------', &
                               '-------------------*'
                WRITE(IMSG,*) ' Thermal loading option ', &
                               'not yet available.'
                WRITE(IMSG,*) '*------------------------', &
                               '-------------------*'
                call hecmw_abort( hecmw_comm_get_comm())
              endif
            endif
!C** Add vector
            do j=1,ndof*nn
              hecMAT%B( iwk(j) )=hecMAT%B( iwk(j) )+vect(j)
            enddo
          enddo
        enddo
        deallocate ( temp )
      endif
    end subroutine DYNAMIC_MAT_ASS_LOAD 
end module m_dynamic_mat_ass_load
