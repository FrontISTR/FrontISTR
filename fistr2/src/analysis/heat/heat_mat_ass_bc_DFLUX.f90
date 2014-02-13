!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Heat Analysis                                     !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module provides a subroutine for setting distributed 
!!  heat flux boundary conditions
module m_heat_mat_ass_bc_DFLUX
     use elementInfo
   contains
!C
!C***
!C*** MAT_ASS_DFLUX
!C***
!C
   subroutine heat_mat_ass_bc_DFLUX( hecMESH, hecMAT, fstrHEAT, CTIME )

      use m_fstr
      use m_heat_get_amplitude
      use m_heat_LIB_DFLUX

      implicit none
      integer(kind=kint) k,icel,ic_type,isect,isuf,iamp,nn,iS,j
      real(kind=kreal)   CTIME,QQ,val,asect,thick
      type (fstr_heat         ) :: fstrHEAT
      type (hecmwST_matrix    ) :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      real(kind=kreal) xx(20),yy(20),zz(20)
      real(kind=kreal) vect(20)
      integer(kind=kint) nodLocal(20)
!C
        do k = 1, fstrHEAT%Q_SUF_tot 

          icel    = fstrHEAT%Q_SUF_elem(k)
          ic_type = hecMESH%elem_type(icel)
          isect   = hecMESH%section_ID(icel)
          isuf    = fstrHEAT%Q_SUF_surf(k)
          iamp    = fstrHEAT%Q_SUF_ampl(k)

          call heat_get_amplitude( fstrHEAT,iamp,CTIME,QQ )
          val     = fstrHEAT%Q_SUF_val (k) * QQ
!C**
          nn = getNumberOfNodes(ic_type)
!C** 
          iS = hecMESH%elem_node_index(icel-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item(iS+j)
            xx(j) = hecMESH%node( 3*nodLOCAL(j)-2 )
            yy(j) = hecMESH%node( 3*nodLOCAL(j)-1 )
            zz(j) = hecMESH%node( 3*nodLOCAL(j)   )
          enddo  
!C**		
          if    ( ic_type.eq.111 ) then
            is = hecMesh%section%sect_R_index(isect)
            asect = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_111(nn,xx,yy,zz,asect,isuf,val,vect)

          elseif( ic_type.eq.231 ) then
            is = hecMesh%section%sect_R_index(isect)
            thick = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_231(nn,xx,yy,zz,thick,isuf,val,vect)

          elseif( ic_type.eq.232 ) then
            is = hecMesh%section%sect_R_index(isect)
            thick = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_232(nn,xx,yy,zz,thick,isuf,val,vect)

          elseif( ic_type.eq.241 ) then
            is = hecMesh%section%sect_R_index(isect)
            thick = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_241(nn,xx,yy,zz,thick,isuf,val,vect)

          elseif( ic_type.eq.242 ) then
            is = hecMesh%section%sect_R_index(isect)
            thick = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_242(nn,xx,yy,zz,thick,isuf,val,vect)

          elseif( ic_type.eq.341 ) then
            call heat_DFLUX_341(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.342 ) then
            call heat_DFLUX_342(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.351 ) then
            call heat_DFLUX_351(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.352 ) then
            call heat_DFLUX_352(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.361 ) then
            call heat_DFLUX_361(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.362 ) then
            call heat_DFLUX_362(nn,xx,yy,zz,isuf,val,vect)

          elseif( ic_type.eq.731 ) then
            is    = hecMesh%section%sect_R_index(isect)
            thick = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_731(nn,xx,yy,zz,thick,isuf,val,vect)

          elseif( ic_type.eq.741 ) then
            is    = hecMesh%section%sect_R_index(isect)
            thick = hecMESH%section%sect_R_item(is)
            call heat_DFLUX_741(nn,xx,yy,zz,thick,isuf,val,vect)

          endif

!C 
          do j = 1, nn
            hecMAT%B( nodLOCAL(j) ) = hecMAT%B( nodLOCAL(j) ) - vect(j) 
          enddo
!C

        enddo

      return
   end subroutine heat_mat_ass_bc_DFLUX 
end module m_heat_mat_ass_bc_DFLUX
