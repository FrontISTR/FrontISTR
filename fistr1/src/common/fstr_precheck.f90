!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Library                                            !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module provides function to check input data of IFSTR solver
module m_fstr_precheck
   contains

   subroutine fstr_get_thickness(hecMESH,mid,thick)
      use hecmw
      use m_fstr
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      TYPE (hecmwST_local_mesh) :: hecMESH

      ihead = hecMESH%section%sect_R_index(mid-1)
      thick = hecMESH%section%sect_R_item(ihead+1)
    !  if(thick.LE.0.0) then
    !      write(*,*) "Zero thickness <= 0 is illegal"
    !      call hecmw_abort( hecmw_comm_get_comm())
    !  endif
   end subroutine fstr_get_thickness
!C
!C***
!C*** Pre Check for FSTR solver
!C***
!C
   subroutine fstr_precheck ( hecMESH )

      use m_fstr

      IMPLICIT DOUBLE PRECISION(a-h,o-z)

      type (hecmwST_local_mesh) :: hecMESH

      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' ****   STAGE PreCheck  **'
      ENDIF

      call fstr_precheck_elem ( hecMESH )
      write(IDBG,*) 'fstr_precheck_elem: OK'

   end subroutine fstr_precheck
!C
!C
   subroutine fstr_precheck_elem ( hecMESH )

      use m_fstr
      use m_precheck_LIB_2d
      use m_precheck_LIB_3d
      use m_precheck_LIB_shell

      implicit REAL(kind=kreal) (A-H,O-Z)

      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)         :: fstrSOLID

!** Local variables
      real(kind=kreal) xx(20),yy(20),zz(20)
      integer(kind=kint) nodLOCAL(20),NTOTsum(1)
      real(kind=kreal) TOTsum(1),TOTmax(3),TOTmin(2)
!C
!C INIT
!C
      nelem  = 0
      tvol   = 0.0
      tvmax  = 0.0
      tvmin  = 1.0e+20
      tlmax  = 0.0
      tlmin  = 1.0e+20
      aspmax = 0.0
!C
!C 3D
!C
      if( hecMESH%n_dof .eq. 3 ) then
        do ie = 1, hecMESH%n_elem
          ia = hecMESH%elem_ID(ie*2)
          if( ia.ne.hecMESH%my_rank ) cycle
!          je = hecMESH%elem_ID(ie*2-1)
          jelem = hecMESH%global_elem_ID(ie)
!C
          ic_type = hecMESH%elem_type(ie)
!C
          if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_solid(ic_type))) then
            write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
            cycle
          endif
          nn = hecmw_get_max_node(ic_type)
!C
          jS = hecMESH%elem_node_index(ie-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
            xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
            yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
            zz(j) = hecMESH%node(3*nodLOCAL(j)  )
          enddo  
!C
          if    ( ic_type.eq.111 ) then
            isect = hecMESH%section_ID(ie)
            mid = hecMESH%section%sect_mat_ID_item(isect)
            CALL fstr_get_thickness( hecMESH,mid,AA )
            al = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
            nline = 1
            tline = al
            vol = AA*al
            almax = al
            almin = al
          elseif( ic_type.eq.341 ) then
            call PRE_341 ( xx,yy,zz,vol,almax,almin )
          elseif( ic_type.eq.351 ) then
            call PRE_351 ( xx,yy,zz,vol,almax,almin )
          elseif( ic_type.eq.361 ) then
            call PRE_361 ( xx,yy,zz,vol,almax,almin )
          elseif( ic_type.eq.342 ) then
            call PRE_342 ( xx,yy,zz,vol,almax,almin )
          elseif( ic_type.eq.352 ) then
            call PRE_352 ( xx,yy,zz,vol,almax,almin )
          elseif( ic_type.eq.362 ) then
            call PRE_362 ( xx,yy,zz,vol,almax,almin )
          endif
!C
          if( vol.le.0.0 ) then
            write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
          endif
          nelem = nelem + 1
          tvol = tvol + vol
          if( vol.gt.tvmax ) tvmax = vol 
          if( vol.lt.tvmin ) tvmin = vol
          if( almax.gt.tlmax ) tlmax = almax
          if( almin.lt.tlmin ) tlmin = almin
          asp = almax/almin
          if( asp.gt.aspmax ) aspmax = asp
          if( asp.gt.50 ) then
            write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
            write(ILOG,*) '      Maximum length =',almax                
            write(ILOG,*) '      Minimum length =',almin                
          endif
        enddo
!C
!C 2D
!C
      elseif( hecMESH%n_dof .eq. 2 ) then
        do ie = 1, hecMESH%n_elem
          ia = hecMESH%elem_ID(ie*2)
          if( ia.ne.hecMESH%my_rank ) cycle
!          je = hecMESH%elem_ID(ie*2-1)
          jelem = hecMESH%global_elem_ID(ie)
!C
          ic_type = hecMESH%elem_type(ie)
!C
          if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_surface(ic_type))) then
            write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
            cycle
          endif
          nn = hecmw_get_max_node(ic_type)
!C
          jS = hecMESH%elem_node_index(ie-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
            xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
            yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
          enddo  
!C
          isect = hecMESH%section_ID(ie)
          mid = hecMESH%section%sect_mat_ID_item(isect)
          CALL fstr_get_thickness( hecMESH,mid,AA )
!C
          if    ( ic_type.eq.111 ) then
            al = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2 )
            vol = AA*al
            if( al.gt.tlmax ) tlmax = al
            if( al.lt.tlmin ) tlmin = al
            aspmax = 1.0 
          elseif( ic_type.eq.231 ) then
            call PRE_231 ( xx,yy,AA,vol,almax,almin )
          elseif( ic_type.eq.241 ) then
            call PRE_241 ( xx,yy,AA,vol,almax,almin )
          elseif( ic_type.eq.232 ) then
            call PRE_232 ( xx,yy,AA,vol,almax,almin )
          elseif( ic_type.eq.242 ) then
            call PRE_242 ( xx,yy,AA,vol,almax,almin )
          else
            vol = 0.0
          endif
!C
          if( vol.le.0.0 ) then
            write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
          endif
          nelem = nelem + 1
          tvol = tvol + vol
          if( vol.gt.tvmax ) tvmax = vol 
          if( vol.lt.tvmin ) tvmin = vol
          if( almax.gt.tlmax ) tlmax = almax
          if( almin.lt.tlmin ) tlmin = almin
          asp = almax/almin
          if( asp.gt.aspmax ) aspmax = asp
          if( asp.gt.50 ) then
            write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
            write(ILOG,*) '      Maximum length =',almax                
            write(ILOG,*) '      Minimum length =',almin                
          endif
        enddo
!C
!C SHELL
!C
      elseif( hecMESH%n_dof .eq. 6 ) then
        do ie = 1, hecMESH%n_elem
          ia = hecMESH%elem_ID(ie*2)
          if( ia.ne.hecMESH%my_rank ) cycle
!          je = hecMESH%elem_ID(ie*2-1)
          jelem = hecMESH%global_elem_ID(ie)
!C
          ic_type = hecMESH%elem_type(ie)
!C
          if (.not. (hecmw_is_etype_beam(ic_type) .or. hecmw_is_etype_shell(ic_type))) then
            write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
            cycle
          endif
          nn = hecmw_get_max_node(ic_type)
!C
          jS = hecMESH%elem_node_index(ie-1)
          do j = 1, nn
            nodLOCAL(j) = hecMESH%elem_node_item (jS+j)
            xx(j) = hecMESH%node(3*nodLOCAL(j)-2)
            yy(j) = hecMESH%node(3*nodLOCAL(j)-1)
            zz(j) = hecMESH%node(3*nodLOCAL(j)  )
          enddo  
!C
          isect = hecMESH%section_ID(ie)
          mid = hecMESH%section%sect_mat_ID_item(isect)
          CALL fstr_get_thickness( hecMESH,mid,AA )
!C
          if    ( ic_type.eq.111 ) then
            al = SQRT( (xx(2)-xx(1))**2+(yy(2)-yy(1))**2+(zz(2)-zz(1))**2 )
            nline = nline + 1
            tline = tline + al
            vol = AA*al
            if( al.gt.tlmax ) tlmax = al
            if( al.lt.tlmin ) tlmin = al
            aspmax = 1.0 
          elseif( ic_type.eq.731 ) then
            call PRE_731 ( xx,yy,zz,AA,vol,almax,almin )
          elseif( ic_type.eq.741 ) then
            call PRE_741 ( xx,yy,zz,AA,vol,almax,almin )
          endif
!C
          if( vol.le.0.0 ) then
            write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',jelem,' is zero or negative.'
          endif
          nelem = nelem + 1
          tvol = tvol + vol
          if( vol.gt.tvmax ) tvmax = vol 
          if( vol.lt.tvmin ) tvmin = vol
          if( almax.gt.tlmax ) tlmax = almax
          if( almin.lt.tlmin ) tlmin = almin
          asp = almax/almin
          if( asp.gt.aspmax ) aspmax = asp
          if( asp.gt.50 ) then
            write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',jelem,' exceeds 50.'
            write(ILOG,*) '      Maximum length =',almax                
            write(ILOG,*) '      Minimum length =',almin                
          endif
        enddo
      endif
!C
      avvol = tvol / nelem 
      write(ILOG,*) '###  Sub Summary  ###'                
      write(ILOG,*) ' Total Volumes in this region        = ',tvol
      write(ILOG,*) ' Average Volume of elements          = ',avvol
      write(ILOG,*) ' Maximum Volume of elements          = ',tvmax
      write(ILOG,*) ' Minimum Volume of elements          = ',tvmin
      write(ILOG,*) ' Maximum length of element edges     = ',tlmax
      write(ILOG,*) ' Minimum length of element edges     = ',tlmin
         
      write(ILOG,*) ' Maximum aspect ratio in this region = ',aspmax
      TOTsum(1) = tvol
      call hecmw_allREDUCE_R(hecMESH,TOTsum,1,hecmw_sum)
      NTOTsum(1) = nelem
      call hecmw_allREDUCE_I(hecMESH,NTOTsum,1,hecmw_sum)
      TOTmax(1) = tvmax
      TOTmax(2) = tlmax
      TOTmax(3) = aspmax
      call hecmw_allREDUCE_R(hecMESH,TOTmax,3,hecmw_max)
      TOTmin(1) = tvmin
      TOTmin(2) = tlmin
      call hecmw_allREDUCE_R(hecMESH,TOTmin,2,hecmw_min)
      if( hecMESH%my_rank .eq. 0 ) then
        avvol = TOTsum(1) / NTOTsum(1)
        write(ILOG,*) '###  Global Summary  ###'
        write(ILOG,*) ' TOTAL VOLUME = ',TOTsum(1)
        write(*,*)    ' TOTAL VOLUME = ',TOTsum(1)
        write(ILOG,*) ' AVERAGE VOLUME OF ELEMENTS = ',avvol
        write(*,*)    ' AVERAGE VOLUME OF ELEMENTS = ',avvol
        write(ILOG,*) ' MAXIMUM VOLUME OF ELEMENTS = ',TOTmax(1)
        write(*,*)    ' MAXIMUM VOLUME OF ELEMENTS = ',TOTmax(1)
        write(ILOG,*) ' MINIMUM VOLUME OF ELEMENTS = ',TOTmin(1)
        write(*,*)    ' MINIMUM VOLUME OF ELEMENTS = ',TOTmin(1)
        write(ILOG,*) ' MAXIMUM LENGTH OF ELEMENT EDGES = ',TOTmax(2)
        write(*,*)    ' MAXIMUM LENGTH OF ELEMENT EDGES = ',TOTmax(2)
        write(ILOG,*) ' MINIMUM LENGTH OF ELEMENT EDGES = ',TOTmin(2)
        write(*,*)    ' MINIMUM LENGTH OF ELEMENT EDGES = ',TOTmin(2)
        write(ILOG,*) ' MAXIMUM ASPECT RATIO  = ',TOTmax(3)
        write(*,*)    ' MAXIMUM ASPECT RATIO  = ',TOTmax(3)
      endif
!C
   end subroutine fstr_precheck_elem
end module m_fstr_precheck
