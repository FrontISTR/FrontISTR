!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                       Noboru Imai (Univ. of Tokyo)                   !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module provide a function to ECHO for IFSTR solver
module m_static_echo
   contains
!C
!C***
!> ECHO for IFSTR solver
!C***
!C
   subroutine FSTR_ECHO (fstrSOLID)
      use m_fstr
      implicit none
      type (fstr_solid )        :: fstrSOLID


      integer(kind=kint) :: i,j,nn,ic_type, cnt
      integer(kind=kint) :: icel,isect,ig1,iS0,iE0,ik
      integer(kind=kint) :: nid, itype
      real(kind=kreal) :: x,y,z
      integer(kind=kint):: nids(0:20)
	  
      integer :: iAss, iPart, iElem, iNode, iGrp
      character(len=HECMW_NAME_LEN) :: header_name
	  
      include "HEC_MW3_For.h"
	  
      write(ILOG,*) "Number of assemble model=", mw_get_num_of_assemble_model()
      call mw_select_assemble_model( 0 )
      write(ILOG,*) "Number of parts=", mw_get_num_of_mesh_part()

      do iAss = 0, mw_get_num_of_assemble_model()-1
         write( ILOG, * ) "Assmble model:", iAss+1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            write(ILOG, *) "Part:",iPart+1, part_nodes(iAss+1, iPart+1), part_nodes(iAss+1, iPart+2)  &
                 , part_elems(iAss+1, iPart+1), part_elems(iAss+1, iPart+2)
         enddo
      enddo
      write(ILOG, *) ""
	  
!C====
!C +-------------------------------+
!C | NODE                          | 
!C +-------------------------------+
!C===
!C** nodal coordinate
      write(ILOG,*) '### Number of nodes',total_node
      write(ILOG,*) 'IASS IPART ID X Y Z'
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            do iNode = 0, mw_get_num_of_node()-1
              nid = mw_get_node_id(iNode)
              call mw_get_node_coord( nid, x,y,z )
              write(ILOG,'(3i8,3e15.5)') iAss,iPart,nid,x,y,z
            enddo
         enddo
      enddo

!C
!C +-------------------------------+
!C | ELEMENT                       | 
!C +-------------------------------+
!C===
      write(ILOG,*) '### Elements', total_elem
      write(ILOG,*) 'IASS IPART TYPE X Y Z'
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            do iElem = 0, mw_get_num_of_element()-1
              call mw_select_element( iElem )
              call  mw_get_element_vert_node_id( nids )
              nn = mw_get_num_of_element_vert()
              ic_type = mw_get_element_type()
              write(ILOG,*) iAss, iPart, ic_type, (nids(j),j=0,nn-1)
            enddo
         enddo
      enddo
	  
!C +-------------------------------+
!C | NODE GROUP                    | 
!C +-------------------------------+
      write(ILOG,*) '### Ngroup'
      write(ILOG,*) 'IASS IPART ID '
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            do iGrp = 0,mw_get_num_of_boundary_bnode_mesh()-1      ! node group
               nn = mw_get_bnode_mesh_namelength( iGrp )
               header_name=''
               call mw_get_bnode_mesh_name(iGrp, header_name(1:nn), nn)
               write(ILOG,*) "Node group:",header_name(1:nn)
			   do  iNode=0, mw_get_num_of_bnode_in_bnode_mesh(iGrp)-1
			        nn= mw_get_node_id_in_bnode_mesh(iGrp, iNode)
                    write(ILOG,*) iAss, iPart, nn
               enddo
            enddo
         enddo
      enddo
!C +-------------------------------+
!C | ELEMEN GROUP                  | 
!C +-------------------------------+
      write(ILOG,*) '### Egroup'
      write(ILOG,*) 'IASS IPART ID '
      do iAss = 0, mw_get_num_of_assemble_model()-1
         call mw_select_assemble_model( iAss )
         do iPart = 0, mw_get_num_of_mesh_part()-1
            call mw_select_mesh_part( iPart )
            do iGrp = 0, mw_get_num_of_elementgroup()-1
              nn = mw_get_elementgroup_name_length(iGrp)
              header_name=''
              call mw_get_elementgroup_name(iGrp, header_name(1:nn), nn)
              write(ILOG,*) "Element group:",header_name(1:nn)
              do iElem = 0, mw_get_num_of_element_id(iGrp)-1
                nn = mw_get_element_id_with_elementgroup( iGrp, iElem )
                write(ILOG,*) iAss, iPart, nn
              enddo
            enddo
         enddo
      enddo

!C +-------------------------------+
!C | MATERIAL                      | 
!C +-------------------------------+	
      write(ILOG,*) ""
      write(ILOG,*) '### Material'
      do i=1, size(fstrSOLID%materials)
          call printMaterial( ILOG, fstrSOLID%materials(i) )
      enddo

!C +-------------------------------+
!C | SECTION                       | 
!C +-------------------------------+	
      write(ILOG,*) ""
      write(ILOG,*) '### Section'
      do i=1, size(MWSections)
          write(ILOG,*) "Material=",trim(MWSections(i)%mat_name),"  EGRP=", trim(MWSections(i)%egroup_name)
          write(ILOG,*) MWSections(i)%sect_type,MWSections(i)%sect_opt,MWSections(i)%sect_mat_ID
          write(ILOG,*) MWSections(i)%sect_R_item
      enddo
	  
      call flush(ILOG)
      print *, "end of echo"
   end subroutine FSTR_ECHO
end module m_static_echo
