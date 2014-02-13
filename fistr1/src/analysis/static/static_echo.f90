!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
   subroutine FSTR_ECHO (hecMESH)
      use m_fstr
      use m_hecmw2fstr_mesh_conv
      implicit none
      type (hecmwST_local_mesh):: hecMESH

!** Local variables
      integer(kind=kint) :: i,j,iS,iE,nn,ic_type
      integer(kind=kint) :: icel,isect,ig1,iS0,iE0,ik
      integer(kind=kint) :: in,nid, itype
      real(kind=kreal) :: x,y,z
      integer(kind=kint):: nids(20)
!C====
!C +-------------------------------+
!C | NODE                          | 
!C +-------------------------------+
!C===
!C** nodal coordinate
      write(ILOG,*) '### Number of nodes',hecMESH%n_node
      write(ILOG,*) 'ID X Y Z'
      do i=1,hecMESH%n_node
        nid = hecMESH%global_node_ID(i)
        x = hecMESH%node(3*i-2)
        y = hecMESH%node(3*i-1)
        z = hecMESH%node(3*i)
        write(ILOG,'(i8,e15.5,e15.5,e15.5)') nid,x,y,z
      enddo
!C
!C +-------------------------------+
!C | ELEMENT                       | 
!C +-------------------------------+
!C===
      call fstr2hecmw_mesh_conv( hecMESH )
      write(ILOG,*) '### Elements', hecMESH%n_elem
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
!C** Set number of nodes		
        nn = hecmw_get_max_node(ic_type)
!C** element loop
        do icel= iS, iE
!C** node ID
          iS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            if( hecMESH%n_refine > 0 ) then
              nids(j)= hecMESH%elem_node_item (iS+j)
            else
              nids(j)= hecMESH%global_node_ID( hecMESH%elem_node_item (iS+j))
            endif
          enddo
!C** section  ID
          isect= hecMESH%section_ID(icel)
          write(ILOG,*) '### Element ID=',ic_type,isect,hecMESH%global_elem_id(icel)
          write(ILOG,*) (nids(j),j=1,nn)
        enddo
      enddo
      call hecmw2fstr_mesh_conv( hecMESH )
!C +-------------------------------+
!C | NODE GROUP                    | 
!C +-------------------------------+
      write(ILOG,*) '### Ngroup'
      do ig1= 1, hecMESH%node_group%n_grp
        write(ILOG,*)
        write(ILOG,'(a80)') hecMESH%node_group%grp_name(ig1)
        iS0= hecMESH%node_group%grp_index(ig1-1) + 1
        iE0= hecMESH%node_group%grp_index(ig1  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          write(ILOG,*) hecMESH%global_node_ID(in)
        enddo
      enddo
!C +-------------------------------+
!C | ELEMEN GROUP                  | 
!C +-------------------------------+
      write(ILOG,*) '### Egroup'
      do ig1= 1, hecMESH%elem_group%n_grp
        write(ILOG,*)
        write(ILOG,'(a80)') hecMESH%elem_group%grp_name(ig1)
        iS0= hecMESH%elem_group%grp_index(ig1-1) + 1
        iE0= hecMESH%elem_group%grp_index(ig1  )
        do ik= iS0, iE0
          in   = hecMESH%elem_group%grp_item(ik)
          write(ILOG,*) hecMESH%global_elem_ID(in)
        enddo
      enddo
      write(ILOG,*) '### Reftemp',ref_temp
!C====
   end subroutine FSTR_ECHO
end module m_static_echo
