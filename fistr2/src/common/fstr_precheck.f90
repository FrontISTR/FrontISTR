!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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

!> This module provides function to check geometric shape of inputted elements
module m_fstr_precheck

   use hecmw
   use m_fstr
   
   implicit none
   
   contains


   subroutine fstr_precheck ( )
      IF(myrank == 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' ****   STAGE PreCheck  **'
      ENDIF

      call fstr_precheck_elem ( )
      write(IDBG,*) 'fstr_precheck_elem: OK'

   end subroutine fstr_precheck

   subroutine fstr_precheck_elem (  )
      use m_precheck_LIB_2d
      use m_precheck_LIB_3d
      use m_precheck_LIB_shell
	  
      include "HEC_MW3_For.h"

      real(kind=kreal) xx(0:20),yy(0:20),zz(0:20)
      real(kind=kreal) :: tvol, tvmax, tvmin, tlmax, tlmin, aspmax, AA
      real(kind=kreal) :: almax, almin, asp, vol, avvol
      integer(kind=kint) nodLOCAL(0:20),NTOTsum(1)
      real(kind=kreal) :: TOTsum(1),TOTmax(3),TOTmin(2)
      integer :: iAss, iPart, iElem, iNode
      integer :: ic_type, nelem, gelem


      do iAss = 0, mw_get_num_of_assemble_model()-1
	     nelem  = 0
         tvol   = 0.0
         tvmax  = 0.0
         tvmin  = 1.0e+20
         tlmax  = 0.0
         tlmin  = 1.0e+20
         aspmax = 0.0
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            do iElem = 0, mw_get_num_of_element()-1
              call mw_select_element( iElem )
              ic_type = mw_get_element_type()
			  ic_type=361
	!      if (.not. (hecmw_is_etype_rod(ic_type) .or. hecmw_is_etype_solid(ic_type))) then
    !        write(ILOG,*) jelem, ' This Element cannot be checked. Type=',ic_type
    !        cycle
    !      endif
              call  mw_get_element_vert_node_id( nodLocal )
              do iNode = 0, mw_get_num_of_element_vert()-1
                 call mw_get_node_coord( nodLocal(iNode), xx(iNode), yy(iNode), zz(iNode) )                 
              enddo
			  
              select case (ic_type)
              case ( 341 )
                call PRE_341 ( xx,yy,zz,vol,almax,almin )
              case ( 351 )
                call PRE_351 ( xx,yy,zz,vol,almax,almin )
              case ( 361 )
                call PRE_361 ( xx,yy,zz,vol,almax,almin )
              case( 342 )
                call PRE_342 ( xx,yy,zz,vol,almax,almin )
              case( 352 )
                call PRE_352 ( xx,yy,zz,vol,almax,almin )
              case( 362 ) 
                call PRE_362 ( xx,yy,zz,vol,almax,almin )
				
              case( 231 )
                call PRE_231 ( xx,yy,AA,vol,almax,almin )
              case( 241 )
                call PRE_241 ( xx,yy,AA,vol,almax,almin )
              case( 232 )
                call PRE_232 ( xx,yy,AA,vol,almax,almin )
              case( 242 )
                 call PRE_242 ( xx,yy,AA,vol,almax,almin )
				 
              case( 731 )
                 call PRE_731 ( xx,yy,zz,AA,vol,almax,almin )
              case( 741 )
                 call PRE_741 ( xx,yy,zz,AA,vol,almax,almin )

              case default
                 vol = 0.0
              end select
			  
              if( vol<=0.0 ) then
                write(ILOG,*) '  %%%  ERROR %%%  Volume of Element no.=',gelem,' is zero or negative.'
              endif
			  
              nelem = nelem + 1
              tvol = tvol + vol
              if( vol>tvmax ) tvmax = vol 
              if( vol<tvmin ) tvmin = vol
              if( almax>tlmax ) tlmax = almax
              if( almin<tlmin ) tlmin = almin
              asp = almax/almin
              if( asp>aspmax ) aspmax = asp
              if( asp>50 ) then
                write(ILOG,*) '  %%%  WARNIG %%% Aspect ratio of Element no.=',gelem,' exceeds 50.'
                write(ILOG,*) '      Maximum length =',almax                
                write(ILOG,*) '      Minimum length =',almin                
              endif
            enddo
         enddo
		 
		 avvol = tvol / nelem 
      write(ILOG,*) '###  Sub Summary  ### level=' , iAss               
      write(ILOG,*) ' Total Volumes in this region        = ',tvol
      write(ILOG,*) ' Average Volume of elements          = ',avvol
      write(ILOG,*) ' Maximum Volume of elements          = ',tvmax
      write(ILOG,*) ' Minimum Volume of elements          = ',tvmin
      write(ILOG,*) ' Maximum length of element edges     = ',tlmax
      write(ILOG,*) ' Minimum length of element edges     = ',tlmin
         
      write(ILOG,*) ' Maximum aspect ratio in this region = ',aspmax
      TOTsum(1) = tvol
      call mw_allreduce_r(TOTsum,1,hecmw_sum)
      NTOTsum(1) = nelem
      call mw_allreduce_i(NTOTsum,1,hecmw_sum)
      TOTmax(1) = tvmax
      TOTmax(2) = tlmax
      TOTmax(3) = aspmax
      call mw_allreduce_r(ToTmax, 3, hecmw_max )
      TOTmin(1) = tvmin
      TOTmin(2) = tlmin
      call mw_allreduce_r(TOTmin,2,hecmw_min)
      if( myrank == 0 ) then
        avvol = TOTsum(1) / NTOTsum(1)
        write(ILOG,*) '###  Global Summary  ### level=', iAss
        write(*,*)    '###  Global Summary  ### level=', iAss
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
      enddo

   end subroutine fstr_precheck_elem
end module m_fstr_precheck
