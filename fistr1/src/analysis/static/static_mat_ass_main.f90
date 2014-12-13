!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by K. Sato (Advancesoft)                          !
!                       X. Yuan (Advancesoft)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!C
!C***
!> CONSTRUCT the GLOBAL STIFF MATRIX 
!C***
!C
MODULE m_static_mat_ass_main

    IMPLICIT NONE
    
    CONTAINS
    
    
    SUBROUTINE FSTR_MAT_ASS_MAIN (hecMESH, hecMAT, fstrSOLID)
        
        USE m_fstr
        
        TYPE( hecmwST_matrix )     :: hecMAT
        TYPE( hecmwST_local_mesh ) :: hecMESH
        TYPE( fstr_solid )         :: fstrSOLID
        
!** Local variables
        REAL(kind=kreal) :: xx(20), yy(20), zz(20), stiffness(20*6, 20*6)
        INTEGER(kind=kint) :: nodLOCAL(20)
        INTEGER(kind=kint) :: ndof, itype, iS, iE, ic_type, nn, icel, iiS, j
        INTEGER(kind=kint) :: flag_shell, flag_solid, mixflag
        INTEGER(kind=kint) :: isect, cdsys_ID
        REAL(kind=kreal) :: coords(3,3)
!C
!C +-------+
!C | INIT. |
!C +-------+
!C====
        CALL hecmw_mat_clear(hecMAT)
        hecMAT%X = 0.d0
!        CALL HECMW_setloglv(HECMW_LOG_DEBUG)
!C
!C +-------------------------------+
!C | ELEMENT-by-ELEMENT ASSEMBLING |
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
        ndof = hecMAT%NDOF
        
        do itype = 1, hecMESH%n_elem_type
            iS= hecMESH%elem_type_index(itype-1) + 1
            iE= hecMESH%elem_type_index(itype  )
            ic_type= hecMESH%elem_type_item(itype)
!C** Ignore link elements
            if (hecmw_is_etype_link(ic_type)) cycle
!C** Set number of nodes
            nn = hecmw_get_max_node(ic_type)
!C element loop
!$omp parallel default(none),private(icel,iiS,j,nodLOCAL,xx,yy,zz,stiffness), &
!$omp&         shared(iS,iE,hecMESH,nn,ndof,ic_type,fstrSOLID,hecMAT)
!$omp do
            do icel= iS, iE
!C** node ID
                iiS= hecMESH%elem_node_index(icel-1)
                do j=1,nn
                    nodLOCAL(j)= hecMESH%elem_node_item (iiS+j)
!C** nodal coordinate
                    xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
                    yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
                    zz(j)=hecMESH%node(3*nodLOCAL(j))
                enddo
!C** Create local stiffness
                isect = hecMESH%section_ID(icel)
                cdsys_ID = hecMESH%section%sect_orien_ID(isect)
                if( cdsys_ID > 0 ) call get_coordsys( cdsys_ID, hecMESH, fstrSOLID, coords )
                call fstr_local_stf_create                                                                                       &
                     (hecMESH, ndof, ic_type, icel, xx, yy, zz, fstrSOLID%elements(icel)%gausses,                                &
                      fstrSOLID%elements(icel)%iset, stiffness, cdsys_ID, coords, fstrSOLID%TEMP_ngrp_tot, fstrSOLID%TEMP_irres) 
!== CONSTRUCT the GLOBAL MATRIX STARTED
                call hecmw_mat_ass_elem(hecMAT, nn, nodLOCAL, stiffness)
            enddo
!$omp end do
!$omp end parallel
        enddo
        
!* for EQUATION
!        CALL hecmw_mat_ass_equation ( hecMESH, hecMAT )
        
    END SUBROUTINE FSTR_MAT_ASS_MAIN
    
    
!> Calculate stiff matrix of current element
    SUBROUTINE FSTR_LOCAL_STF_CREATE                                          &
               (hecMESH, ndof, ic_type, icel, xx, yy, zz, gausses,            &
                iset, stiffness, cdsys_ID, coords, TEMP_ngrp_tot, TEMP_irres) 
        
        USE m_fstr
        USE m_static_lib
        USE mMechGauss
        
        TYPE( hecmwST_local_mesh ) :: hecMESH
        INTEGER(kind=kint) :: ndof, ic_type, icel, iset
        REAL(kind=kreal) :: xx(:), yy(:), zz(:)
        TYPE( tGaussStatus ), INTENT(IN) :: gausses(:)
        REAL(kind=kreal) :: stiffness(:, :)
        INTEGER(kind=kint) :: cdsys_ID
        REAL(kind=kreal) :: coords(3, 3)
        INTEGER(kind=kint) :: TEMP_ngrp_tot, TEMP_irres
        
        !** Local variables
        REAL(kind=kreal) :: ee, pp, thick, ecoord(3, 20)
        TYPE( tMaterial ), POINTER :: material
        REAL(kind=kreal) :: local_stf(1830)
        INTEGER(kind=kint) :: nn, isect, ihead, mixflag
        
        nn = hecmw_get_max_node(ic_type)
        ecoord(1,1:nn) = xx(1:nn)
        ecoord(2,1:nn) = yy(1:nn)
        ecoord(3,1:nn) = zz(1:nn)
        material => gausses(1)%pMaterial
        ee = material%variables(M_YOUNGS)
        pp = material%variables(M_POISSON)
        
        if ( ic_type == 241 .or. ic_type == 242 .or.                        &
             ic_type == 231 .or. ic_type == 232 .or. ic_type == 2322 ) then 
            thick =1.d0
            CALL STF_C2( ic_type,nn,ecoord(1:2,1:nn),gausses(:),thick,stiffness(1:nn*ndof,1:nn*ndof),iset)
            
        else if( ic_type == 301 ) then
            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)
            thick = hecMESH%section%sect_R_item(ihead+1)
            CALL STF_C1( ic_type,nn,ecoord(:,1:nn),thick,gausses(:),stiffness(1:nn*ndof,1:nn*ndof) )
            
        else if( ic_type == 781 ) then   !for shell-solid mixed analysis
            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)
            thick = hecMESH%section%sect_R_item(ihead+1)
            nn = 4; ndof = 6; ic_type = 741; mixflag = 1
            CALL STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:4), gausses(:), stiffness(1:nn*ndof, 1:nn*ndof), thick, mixflag)
            ic_type = 781; nn = 8; ndof = 3
            return
            
        else if( ic_type == 761 ) then   !for shell-solid mixed analysis
            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)
            thick = hecMESH%section%sect_R_item(ihead+1)
            nn = 3; ndof = 6; ic_type = 731; mixflag = 2
            CALL STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:3), gausses(:), stiffness(1:nn*ndof, 1:nn*ndof), thick, mixflag)
            ic_type = 761; nn = 6; ndof = 3
            return
            
        else if( ic_type == 361 ) then
            if( TEMP_ngrp_tot > 0 .or. TEMP_irres > 0 ) then
             CALL STF_C3D8Bbar( ic_type, nn, ecoord(:, 1:nn), gausses(:), stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords, 1.0D0 )
            else
             CALL STF_C3D8IC( ic_type, nn, ecoord(:, 1:nn), gausses(:), stiffness(1:nn*ndof, 1:nn*ndof), cdsys_ID, coords )
            endif
            
        else if( ic_type == 341 .or. ic_type == 351 .or.                       &
                 ic_type == 342 .or. ic_type == 352 .or. ic_type == 362 ) then 
            CALL STF_C3( ic_type,nn,ecoord(:,1:nn),gausses(:),stiffness(1:nn*ndof,1:nn*ndof), cdsys_ID, coords, 1.0D0)
            
        else if( ( ic_type == 741 ) .or. ( ic_type == 743 ) .or. ( ic_type == 731 ) ) then
            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)
            thick = hecMESH%section%sect_R_item(ihead+1)
            mixflag = 0
            CALL STF_Shell_MITC(ic_type, nn, ndof, ecoord(1:3, 1:nn), gausses(:), stiffness(1:nn*ndof, 1:nn*ndof), thick, mixflag)
            
        else if( ic_type == 611) then
            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)
            CALL STF_Beam(ic_type,nn,ecoord,hecMESH%section%sect_R_item(ihead+1:),ee, pp,stiffness(1:nn*ndof,1:nn*ndof))
            
        else if( ic_type == 641 ) then
            isect= hecMESH%section_ID(icel)
            ihead = hecMESH%section%sect_R_index(isect-1)
            CALL STF_Beam_641(ic_type, nn, ecoord,gausses(:), hecMESH%section%sect_R_item(ihead+1:), stiffness(1:nn*ndof,1:nn*ndof))
            
        else
            write(*,*) '###ERROR### : Element type not supported for linear static analysis'
            write(*,*) ' ic_type = ', ic_type
            CALL hecmw_abort(hecmw_comm_get_comm())
            
        endif
        
    END SUBROUTINE FSTR_LOCAL_STF_CREATE
    
    
END MODULE m_static_mat_ass_main
