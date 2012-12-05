!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.2                                   !
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
!> This module contains functions for output in dynamic analysis
!C================================================================C

module m_dynamic_output
use m_fstr
use m_static_lib
use m_static_get_prop
use m_dynamic_post
use m_hecmw2fstr_mesh_conv
use m_static_make_result

implicit none


contains

!C***
!C*** OUTPUT for FSTR solver
!C***
!C
!C 3D SOLID
!C
      subroutine dynamic_output_3d(hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC)

      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_param)         :: fstrPARAM
      type (fstr_solid)         :: fstrSOLID
      type (fstr_dynamic)       :: fstrDYNAMIC
      integer(kind=kint) :: idx, i, ii, jj, k
!C
!C-- displacement, velocity and accelaration
!C
      if((fstrDYNAMIC%idx_eqa .eq. 1) .and. (fstrDYNAMIC%i_step .gt. 0)) then
        idx = 2
      else
        idx = 1
      endif
!C
!C-- SET DISPLACEMENT
!C
      if( fstrDYNAMIC%nlflag==0 ) then
      do i= 1, hecMESH%n_node
        fstrSOLID%unode(3*i-2)= fstrDYNAMIC%DISP(3*i-2, idx)
        fstrSOLID%unode(3*i-1)= fstrDYNAMIC%DISP(3*i-1, idx)
        fstrSOLID%unode(3*i  )= fstrDYNAMIC%DISP(3*i  , idx)
      enddo
      endif

!C
!C-- DISPLACEMENT
!C
      if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then

      write(ILOG,*) '#### DISPLACEMENT 3D'
      write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      Z-DISP'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%DISP(3*(ii-1)+k, idx),k=1,3)
      enddo

      end if
!C
!C-- VELOCITY
!C
      if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### VELOCITY 3D'
      write(ILOG,'(a)')'    NODE     X-VELO      Y-VELO      Z-VELO'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%VEL(3*(ii-1)+k, idx),k=1,3)
      enddo

      end if
!C
!C-- ACCELERATION
!C
      if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### ACCELERATION 3D'
      write(ILOG,'(a)')'    NODE     X-ACCE      Y-ACCE      Z-ACCE'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%ACC(3*(ii-1)+k, idx),k=1,3)
      enddo

      end if
!!
!C
!C-- end of displacement, velocity and accelaration
!C
!C
!C-- reaction force, strain and stress
!C
!      if(fstrDYNAMIC%i_step .gt. 0) then
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 3D'
      call dynamic_reaction_force_3d( hecMESH, fstrSOLID, fstrDYNAMIC)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRESS 3D'
      call dynamic_nodal_stress_3d( hecMESH, fstrSOLID,fstrPARAM, fstrDYNAMIC)
!!
!      end if
!C
!C-- end of reaction force, strain and stress
!C
      end subroutine dynamic_output_3d
!C
!C NODAL STRESS 3D
!C
      subroutine dynamic_nodal_stress_3d(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      REAL(kind=kreal) xx(20),yy(20),zz(20)
      REAL(kind=kreal) edisp(60),estrain(6),estress(6)
      REAL(kind=kreal) val, tt(20),tt0(20)                 
      integer(kind=kint) nodLOCAL(20),ierror 
      REAL(kind=kreal) edstrain(20,6),edstress(20,6)
!** Array for nodal recovery
      integer (kind=KINT), dimension(:), allocatable :: nnumber
      real    (kind=KREAL), dimension(:), allocatable :: temp
      real(kind=KREAL), allocatable :: ndstrain(:,:),ndstress(:,:)
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, ielem, ID_area, ii,jj
      real(kind=KREAL) :: s11, s22, s33, s12, s23, s13, ps, smises

!*Allocate array
      allocate ( ndstrain(hecMESH%n_node,6), ndstress(hecMESH%n_node,7) ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_3d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      allocate ( nnumber(hecMESH%n_node)           ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_3d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
!*ZERO CLEAR 
      ndstrain=0.d0; ndstress=0.d0
      nnumber=0
!C
!C Set Temperature
!C
      allocate ( temp(hecMESH%n_node) ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_3d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      temp=0
      if( fstrSOLID%TEMP_ngrp_tot > 0 ) then
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
      endif
!C
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
            zz(j)=hecMESH%node(3*nodLOCAL(j)  )
            tt(j)=temp( nodLOCAL(j) )
            tt0(j)=ref_temp
            edisp(3*j-2)=fstrSOLID%unode(3*nodLOCAL(j)-2)
            edisp(3*j-1)=fstrSOLID%unode(3*nodLOCAL(j)-1)
            edisp(3*j  )=fstrSOLID%unode(3*nodLOCAL(j)  )
          enddo  
!C** Create local stiffness
          if (ic_type==361) then
            call UpdateST_C3D8IC(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn),        &
                           edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          else
            call UpdateST_C3(ic_type,nn,xx(1:nn),yy(1:nn),zz(1:nn),tt(1:nn),tt0(1:nn)   &
                          ,edisp(1:nn*3),fstrSOLID%elements(icel)%gausses )
          endif
          call NodalStress_C3(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:) )
!C** elem ID
!!        ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            call ElementStress_C3(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            fstrSOLID%ESTRAIN(6*ielem-5:6*ielem) = estrain
            fstrSOLID%ESTRESS(7*ielem-6:7*ielem-1) = estress
            s11=fstrSOLID%ESTRESS(7*ielem-6)
            s22=fstrSOLID%ESTRESS(7*ielem-5)
            s33=fstrSOLID%ESTRESS(7*ielem-4)
            s12=fstrSOLID%ESTRESS(7*ielem-3)
            s23=fstrSOLID%ESTRESS(7*ielem-2)
            s13=fstrSOLID%ESTRESS(7*ielem-1)
            ps=(s11+s22+s33)/3.d0
            smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                  +s12**2+s23**2+s13**2
            fstrSOLID%ESTRESS(7*ielem)=sqrt(3.0*smises)
          end if

          do j=1,nn
            ndstrain( nodLocal(j),: ) = ndstrain(nodLOCAL(j),:) + edstrain(j,:)
            ndstress( nodLocal(j),1:6 ) = ndstress(nodLOCAL(j),1:6) + edstress(j,:)
            nnumber( nodLOCAL(j) )=nnumber( nodLOCAL(j) )+1  
          enddo
        enddo
      enddo
!** Average over nodes
      do i=1,hecMESH%n_node
        ndstrain(i,:)=ndstrain(i,:)/nnumber(i)
        ndstress(i,1:6)=ndstress(i,1:6)/nnumber(i)
      enddo
!** RE-CALCULATE von MISES stress
      do i=1,hecMESH%n_node
        s11=ndstress(i,1)
        s22=ndstress(i,2)
        s33=ndstress(i,3)
        s12=ndstress(i,4)
        s23=ndstress(i,5)
        s13=ndstress(i,6)
        ps=(s11+s22+s33)/3.d0
        smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
              +s12**2+s23**2+s13**2
        ndstress(i,7)=sqrt(3.d0*smises)
      enddo
!** Set Array
      do i=1,hecMESH%n_node
        do j=1,6
          fstrSOLID%STRAIN(6*(i-1)+j)=ndstrain(i,j)
        enddo
        do j=1,7
          fstrSOLID%STRESS(7*(i-1)+j)=ndstress(i,j)
        enddo
      enddo
!** Show

      if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then

      write(ILOG,*) '#### STRAIN'
      write(ILOG,'(a,a)')                                                &
                    '     NODE     E11         E22         E33     '     &
                             ,'    E12         E23         E13     '
      write(ILOG,'(a,a)')                                                &
                    '  -------+-----------+-----------+-----------+'     &
                             ,'-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p6e12.4)') jj,(ndSTRAIN(ii,k),k=1,6)
      enddo

      end if
      if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then

      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                                &
                    '     NODE     S11         S22         S33     '     &
                 ,'    S12         S23         S13        MISES'
      write(ILOG,'(a,a)')                                                &
                    '  -------+-----------+-----------+-----------+'     &
                 ,'-----------+-----------+-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p7e12.4)') jj,(ndSTRESS(ii,k),k=1,7)
      enddo

      end if

!*Deallocate array
      deallocate( ndstrain,ndstress   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_3d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( nnumber      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_3d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( temp         ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_3d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      end subroutine dynamic_nodal_stress_3d
!C
!C REACTION FORCE 3D
!C
      subroutine dynamic_reaction_force_3d(hecMESH,fstrSOLID, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      integer, parameter :: NDOF=3
      REAL(kind=kreal) stiff(120,120)
      REAL(kind=kreal) edisp(60),force(60),ecoord(3,20)
      integer(kind=kint) nodLOCAL(20),ierror
      real    (kind=KREAL),dimension(:), allocatable :: spcForce
      integer (kind=KINT), dimension(:), allocatable :: id_spc
	  integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, iflag, icount, nspc, nid


!C Allocate array
      allocate ( id_spc( hecMESH%n_node )     ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_3d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      do i=1,hecMESH%n_node
        id_spc(i)=0
      enddo
!C
!C  SPC Nodal Check Flag
!C 
      icount=0
      do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        iS0= hecMESH%node_group%grp_index(ig-1) + 1
        iE0= hecMESH%node_group%grp_index(ig  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          if( id_spc(in) .eq. 0 ) then
            icount=icount+1
            id_spc(in)=icount
          endif
        enddo
      enddo

      icount=0
      do i=1,hecMESH%n_node
        if( id_spc(i).ge.1 ) icount=icount+1
      enddo
      nspc=icount
!C Allocate array
      if( nspc .gt. 0 ) then
         allocate ( spcForce( nspc * NDOF )   ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_3d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

         do i=1,nspc*NDOF
           spcForce(i)=0.0
         enddo
      else
        deallocate ( id_spc    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_3d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
        return
      endif
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_solid(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            ecoord(1,j)=hecMESH%node(3*nodLOCAL(j)-2)
            ecoord(2,j)=hecMESH%node(3*nodLOCAL(j)-1)
            ecoord(3,j)=hecMESH%node(3*nodLOCAL(j)  )
            edisp(3*j-2)=fstrSOLID%unode(3*nodLOCAL(j)-2)
            edisp(3*j-1)=fstrSOLID%unode(3*nodLOCAL(j)-1)
            edisp(3*j  )=fstrSOLID%unode(3*nodLOCAL(j)  )
          enddo
!C** check node id
          iflag=0
          do j=1,nn
            if( id_spc(  nodLOCAL(j) ) .ge. 1 ) iflag=1
          enddo
          if( iflag == 1 ) then
!C** Create local stiffness
            if (ic_type==361) then
              call STF_C3D8IC( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3))
            else
              call STF_C3( ic_type,nn,ecoord(:,1:nn),fstrSOLID%elements(icel)%gausses,stiff(1:nn*3,1:nn*3),1.d0)
            endif
            force(1:nn*3)= matmul( stiff(1:nn*3,1:nn*3), edisp(1:nn*3) )
!C*** Add SPC FORCE ****
            do j=1,nn
              nid=id_spc(  nodLOCAL(j) )
              if( nid .ge. 1 ) then
                do k=1,NDOF
                  spcForce(NDOF*(nid-1)+k)=spcForce(NDOF*(nid-1)+k) +force(NDOF*(j-1)+k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!C*** Show

      if( fstrDYNAMIC%iout_list(4) .eq. 1 ) then

      write(ILOG,'(a)')'    NODE     X-REAC      Y-REAC      Z-REAC'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i=1,hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        nid=id_spc(i)
        if( nid.ge.1 ) then
          write(ILOG,'(i10,1p3e12.4)') j,  &                          
          spcForce(NDOF*(nid-1)+1),        &                         
          spcForce(NDOF*(nid-1)+2),        &                        
          spcForce(NDOF*(nid-1)+3)
        endif
      enddo

      end if

!C Deallocate array
      deallocate ( spcForce     ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_3d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate ( id_spc       ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_3d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_reaction_force_3d
!C***
!C*** OUTPUT for FSTR solver
!C***
!C
!C 2D SOLID
!C
      subroutine DYNAMIC_OUTPUT_2D(hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC)
      type (hecmwST_matrix)     :: hecMAT
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_param)         :: fstrPARAM
      type (fstr_solid)         :: fstrSOLID

      type (fstr_dynamic      ) :: fstrDYNAMIC
      integer(kind=kint) :: idx
      integer :: i,k, ii,jj
!C
!C-- displacement, velocity and accelaration
!C
      if((fstrDYNAMIC%idx_eqa .eq. 1) .and. (fstrDYNAMIC%i_step .gt. 0)) then
        idx = 2
      else
        idx = 1
      endif
!C
!C-- SET DISPLACEMENT
!C
      do i= 1, hecMESH%n_node
        fstrSOLID%unode(2*i-1)= fstrDYNAMIC%DISP(2*i-1, idx)
        fstrSOLID%unode(2*i  )= fstrDYNAMIC%DISP(2*i  , idx)
      enddo

!C
!C-- DISPLACEMENT
!C
      if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then

      write(ILOG,*) '#### DISPLACEMENT 2D'
      write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      '
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,2e12.4)') jj,(fstrDYNAMIC%DISP(2*(ii-1)+k, idx),k=1,2)
      enddo

      end if
!C
!C-- VELOCITY
!C
      if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### VELOCITY 2D'
      write(ILOG,'(a)')'    NODE     X-VELO      Y-VELO      '
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,2e12.4)') jj,(fstrDYNAMIC%VEL(2*(ii-1)+k, idx),k=1,2)
      enddo

      end if
!C
!C-- ACCELERATION
!C
      if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then

       if(fstrDYNAMIC%idx_eqa .eq. 11) then
        if(fstrDYNAMIC%i_step .gt. 0) then
          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &
                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if
       end if

      write(ILOG,*) '#### ACCELERATION 2D'
      write(ILOG,'(a)')'    NODE     X-ACCE      Y-ACCE      '
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

        write(ILOG,'(i10,2e12.4)') jj,(fstrDYNAMIC%ACC(2*(ii-1)+k, idx),k=1,2)
      enddo

      end if
!C
!C-- end of displacement, velocity and accelaration
!C
!C
!C-- reaction force, strain and stress
!C
!      if(fstrDYNAMIC%i_step .gt. 0) then
!C
!C-- REACTION FORCE
!C
      write(ILOG,*) '#### REACTION FORCE 2D'
      call dynamic_reaction_force_2d( hecMESH, fstrSOLID, fstrDYNAMIC)
!C
!C-- NODAL STRESS
!C
      write(ILOG,*) '#### NODAL STRESS 2D'
      call dynamic_nodal_stress_2d( hecMESH, fstrSOLID, fstrPARAM, fstrDYNAMIC)
!!
!      end if
!C
!C-- end of reaction force, strain and stress
!C
      end subroutine DYNAMIC_OUTPUT_2D

!C
!C NODAL STRESS 2D
!C
      subroutine dynamic_nodal_stress_2d(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      REAL(kind=kreal) xx(8),yy(8)
      REAL(kind=kreal) edisp(16),estrain(4),estress(4)
      REAL(kind=kreal) val, tt(8),tt0(8)
      integer(kind=kint) nodLOCAL(8),ierror
      REAL(kind=kreal) edstrain(9,4),edstress(9,4)
!** Array for nodal recovery
      real    (kind=KREAL), dimension(:,:), allocatable :: arrayTotal
      integer (kind=KINT), dimension(:), allocatable :: nnumber
      real    (kind=KREAL), dimension(:), allocatable :: temp
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel
      integer :: i,j,k, in, itype, ielem, ID_area, ii,jj
      real(kind=KREAL) :: s11, s22, s33, s12, s23, s13, ps, smises
	  
!*Allocate array
      allocate ( arrayTotal( hecMESH%n_node, 9  )    ,STAT=ierror) 
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_2d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      allocate ( nnumber(hecMESH%n_node)             ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_2d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
!*ZERO CLEAR
      arrayTotal=0
      nnumber=0
!C
!C Set Temperature
!C
      allocate ( temp(hecMESH%n_node)      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_nodal_stress_2d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      temp=0
      if( fstrSOLID%TEMP_ngrp_tot > 0 ) then
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
      endif
!C
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_surface(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            xx(j)=hecMESH%node(3*nodLOCAL(j)-2)
            yy(j)=hecMESH%node(3*nodLOCAL(j)-1)
            tt(j)=temp( nodLOCAL(j) )
            tt0(j)=ref_temp
            edisp(2*j-1)=fstrSOLID%unode(2*nodLOCAL(j)-1)
            edisp(2*j  )=fstrSOLID%unode(2*nodLOCAL(j))
          enddo  

!C** Create local stiffness
          call UpdateST_C2(ic_type,nn,xx(1:nn),yy(1:nn),tt(1:nn),tt0(1:nn),1.d0,              &
                     fstrSOLID%elements(icel)%iset,edisp(1:nn*2),fstrSOLID%elements(icel)%gausses )
          call NodalStress_C2(ic_type,nn,fstrSOLID%elements(icel)%gausses,edstrain(1:nn,:),edstress(1:nn,:) )
!C** elem ID
!!        ielem = hecMESH%elem_ID(icel*2-1)
          ielem = icel
          ID_area = hecMESH%elem_ID(icel*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            call ElementStress_C2(ic_type,fstrSOLID%elements(icel)%gausses,estrain,estress)
            fstrSOLID%ESTRAIN(6*ielem-5:6*ielem-2) = estrain
            fstrSOLID%ESTRESS(7*ielem-6:7*ielem-3) = estress
            s11=fstrSOLID%ESTRESS(7*ielem-6)
            s22=fstrSOLID%ESTRESS(7*ielem-5)
            s33=fstrSOLID%ESTRESS(7*ielem-4)
            s12=fstrSOLID%ESTRESS(7*ielem-3)
            s23=0.0D0
            s13=0.0D0
            ps=(s11+s22+s33)/3.0
            smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )         &
                +s12**2+s23**2+s13**2
            fstrSOLID%ESTRESS(7*ielem)=sqrt(3.0*smises)
          end if

          do j=1,nn
            do k=1,4
              arrayTotal(nodLOCAL(j),k)=                                 &
              arrayTotal(nodLOCAL(j),k)+edstrain(j,k)
            enddo
            do k=1,4
              arrayTotal(nodLOCAL(j),k+4)=                                 &
              arrayTotal(nodLOCAL(j),k+4)+edstress(j,k)
            enddo
            nnumber( nodLOCAL(j) )=nnumber( nodLOCAL(j) )+1
          enddo

        enddo
      enddo
!** Average over nodes
      do i=1,hecMESH%n_node
        do j=1,9
          arrayTotal(i,j)=arrayTotal(i,j)/FLOAT(nnumber(i))
        enddo
      enddo
!** RE-CALCULATE von MISES stress
      do i=1,hecMESH%n_node
        s11=arrayTotal(i,5)
        s22=arrayTotal(i,6)
        s33=arrayTotal(i,8)
        s12=arrayTotal(i,7)
        s23=0.0
        s13=0.0
        ps=(s11+s22+s33)/3.0
        smises=0.5*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 ) +s12**2+s23**2+s13**2
        arrayTotal(i,9)=sqrt(3.0*smises)
      enddo
!** Set Array 
      do i=1,hecMESH%n_node
        fstrSOLID%STRAIN(6*(i-1)+1)=arrayTotal(i,1)
        fstrSOLID%STRAIN(6*(i-1)+2)=arrayTotal(i,2)
        fstrSOLID%STRAIN(6*(i-1)+3)=arrayTotal(i,4)
        fstrSOLID%STRAIN(6*(i-1)+4)=arrayTotal(i,3)
        fstrSOLID%STRAIN(6*(i-1)+5)=0.0
        fstrSOLID%STRAIN(6*(i-1)+6)=0.0
        fstrSOLID%STRESS(7*(i-1)+1)=arrayTotal(i,5)
        fstrSOLID%STRESS(7*(i-1)+2)=arrayTotal(i,6)
        fstrSOLID%STRESS(7*(i-1)+3)=arrayTotal(i,8)
        fstrSOLID%STRESS(7*(i-1)+4)=arrayTotal(i,7)
        fstrSOLID%STRESS(7*(i-1)+5)=0.0
        fstrSOLID%STRESS(7*(i-1)+6)=0.0
        fstrSOLID%STRESS(7*(i-1)+7)=arrayTotal(i,9)
      enddo
!** Show
      if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
      write(ILOG,*) '#### STRIN'
      write(ILOG,'(a,a)')                                         & 
               '   NODE       E11         E22         E33     '   &
                             ,'    E12' 
      write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p4e12.4)') jj,(fstrSOLID%STRAIN(6*(ii-1)+k),k=1,4)
      enddo
      end if

      if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
      write(ILOG,*) '#### STRESS'
      write(ILOG,'(a,a)')                                         &
               '   NODE       S11         S22         S33     '   &
                             ,'    S12        MISES' 
      write(ILOG,'(a,a)')                                             &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-----------+-------------'
      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)
        write(ILOG,'(i10,1p5e12.4)') jj,                               &
                                 (fstrSOLID%STRESS(7*(ii-1)+k),k=1,4)  &
                                 ,fstrSOLID%STRESS(7*(ii-1)+7)
      enddo
      end if
!*Deallocate array
      deallocate( arrayTotal      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_2d, arrayTotal>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( nnumber         ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_2d, nnumber>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate( temp            ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_nodal_stress_2d, temp>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_nodal_stress_2d

!C
!C REACTION FORCE 2D
!C
      subroutine dynamic_reaction_force_2d(hecMESH,fstrSOLID, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_dynamic      ) :: fstrDYNAMIC
!** Local variables
      integer, parameter :: NDOF=3
      REAL(kind=kreal) ecoord(2,8),stiff(18,18)
      REAL(kind=kreal) val, edisp(16),force(16)
      integer(kind=kint) nodLOCAL(8),ierror
      real    (kind=KREAL),dimension(:), allocatable :: spcForce
      integer (kind=KINT), dimension(:), allocatable :: id_spc
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, jS, icel, iflag, nn
      integer :: i,j,k, in, ielem, ii,jj, icount, nspc, nid, itype

!C Allocate array
      allocate ( id_spc( hecMESH%n_node )    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_2d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      do i=1,hecMESH%n_node
        id_spc(i)=0                 
      enddo
!C
!C  SPC Nodal Check Flag
!C
      icount=0
      do ig0= 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig= fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        iS0= hecMESH%node_group%grp_index(ig-1) + 1
        iE0= hecMESH%node_group%grp_index(ig  )
        do ik= iS0, iE0
          in   = hecMESH%node_group%grp_item(ik)
          if( id_spc(in) .eq. 0 ) then
            icount=icount+1
            id_spc(in)=icount
          endif
        enddo
      enddo

      icount=0
      do i=1,hecMESH%n_node
        if( id_spc(i).ge.1 ) icount=icount+1
      enddo
      nspc=icount
!C Allocate array
      if( nspc .gt. 0 ) then
         allocate ( spcForce( nspc * NDOF )    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to allocation error <dynamic_reaction_force_2d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

         do i=1,nspc*NDOF
           spcForce(i)=0.0
         enddo
      else
        deallocate ( id_spc     ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_2d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

        return
      endif
!C +-------------------------------+
!C | according to ELEMENT TYPE     |
!C +-------------------------------+
      do itype= 1, hecMESH%n_elem_type
        iS= hecMESH%elem_type_index(itype-1) + 1
        iE= hecMESH%elem_type_index(itype  )
        ic_type= hecMESH%elem_type_item(itype)
        if (.not. hecmw_is_etype_surface(ic_type)) cycle
!C** Set number of nodes
        nn = hecmw_get_max_node(ic_type)
!C element loop
        do icel= iS, iE
!C** node ID
          jS= hecMESH%elem_node_index(icel-1)
          do j=1,nn
            nodLOCAL(j)= hecMESH%elem_node_item (jS+j)
!C** nodal coordinate
            ecoord(1,j)=hecMESH%node(3*nodLOCAL(j)-2)
            ecoord(2,j)=hecMESH%node(3*nodLOCAL(j)-1)
            edisp(2*j-1)=fstrSOLID%unode(2*nodLOCAL(j)-1)
            edisp(2*j  )=fstrSOLID%unode(2*nodLOCAL(j)  )
          enddo
!C** check node id
          iflag=0
          do j=1,nn
            if( id_spc(  nodLOCAL(j) ) .ge. 1 ) iflag=1
          enddo
          if( iflag == 1 ) then
!C** Create local stiffness
            call STF_C2( ic_type,nn,ecoord(1:2,1:nn),fstrSOLID%elements(icel)%gausses,1.d0        &
                        ,stiff(1:nn*2,1:nn*2),fstrSOLID%elements(icel)%iset)
            force = matmul( stiff(1:nn*2,1:nn*2), edisp(1:nn*2) )
!C*** Add SPC FORCE ****
            do j=1,nn
              nid=id_spc( nodLOCAL(j) )
              if( nid .ge. 1 ) then
                do k=1,NDOF
                  spcForce(NDOF*(nid-1)+k)=spcForce(NDOF*(nid-1)+k) +force(NDOF*(j-1)+k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
!C*** Show
      if( fstrDYNAMIC%iout_list(4) .eq. 1 ) then

      write(ILOG,'(a)') '   NODE    X-REACTION  Y-REACTION'
      write(ILOG,'(a)') ' --------+-----------+------------'
      do i=1,hecMESH%nn_internal
        j=hecMESH%global_node_ID(i)
        nid=id_spc(i)
        if( nid.ge.1 ) then
          write(ILOG,'(i10,1p2e12.4)') j,spcForce(NDOF*(nid-1)+1), spcForce(NDOF*(nid-1)+2)
        endif
      enddo

      end if
!C Deallocate array
      deallocate ( spcForce    ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_2d, spcForce>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if
      deallocate ( id_spc      ,STAT=ierror)
            if( ierror /= 0 ) then
              write(idbg,*) 'stop due to deallocation error <dynamic_reaction_force_2d, id_spc>'
              write(idbg,*) '  rank = ', hecMESH%my_rank,'  ierror = ',ierror
              call flush(idbg)
              call hecmw_abort( hecmw_comm_get_comm())
            end if

      end subroutine dynamic_reaction_force_2d
!C***
!C*** OUTPUT for FSTR solver
!C***


      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE dynamic_output_6d_Shell                              &
                 (hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC) 
!####################################################################
      use m_static_output  
!--------------------------------------------------------------------
      
      TYPE(hecmwST_matrix)     :: hecMAT
      TYPE(hecmwST_local_mesh) :: hecMESH
      TYPE(fstr_solid)         :: fstrSOLID
      TYPE(fstr_param)         :: fstrPARAM
      TYPE(fstr_dynamic)       :: fstrDYNAMIC
      
!--------------------------------------------------------------------
      
      INTEGER(KIND = kint) :: i, j, k
      INTEGER(KIND = kint) :: ii, jj
      INTEGER(KIND = kint) :: ndof
      INTEGER(KIND = kint) :: idx
      
!--------------------------------------------------------------------
      
      ! ndof = hecMESH%n_dof
      ndof = hecMESH%n_dof
      
!--------------------------------------------------------------------
      
      ! Displacement, Velocity and Accelaration
      IF( ( fstrDYNAMIC%idx_eqa .EQ. 1 ) .AND. &
          ( fstrDYNAMIC%i_step .GT. 0 ) ) THEN 
       
       idx = 2
       
      ELSE
       
       idx = 1
       
      END IF
      
!--------------------------------------------------------------------
      
      ! Displacement
      DO i = 1, hecMESH%n_node
       
       DO j = 1, ndof
        
        fstrSOLID%unode( ndof*(i-1)+j )= fstrDYNAMIC%DISP( ndof*(i-1)+j, idx )
        
       END DO
       
      END DO
      
      IF( fstrDYNAMIC%iout_list(1) .EQ. 1 ) THEN
       
       WRITE(ILOG, *) '#### DISPLACEMENT SHELL'
       WRITE(ILOG, '(A, A)')                                         &
                   '     NODE    X-DISP      Y-DISP      Z-DISP   ', &
                   '   X-ROT       Y-ROT       Z-ROT'                
       WRITE(ILOG, '(A, A)')                                         &
                   ' --------+-----------+-----------+-----------+', &
                   '-----------+-----------+-------------'           
       
       DO i= 1, hecMESH%nn_internal
        
        jj = fstrPARAM%global_local_id(1, i)
        ii = fstrPARAM%global_local_id(2, i)
        
        WRITE( ILOG, '(I10, 6E12.4)' )                                     &
               jj, ( fstrDYNAMIC%DISP( ndof*(ii-1)+k, idx ), k = 1, ndof ) 
        
       END DO
       
      END IF

!--------------------------------------------------------------------
      
      ! Velocity
      IF( fstrDYNAMIC%iout_list(2) .EQ. 1 ) THEN
       
       IF( fstrDYNAMIC%idx_eqa .EQ. 11 ) THEN
        
        IF( fstrDYNAMIC%i_step .GT. 0 ) THEN
         
         WRITE( ILOG, '(A, I10, 5X, A, 1PE12.4)' )                            &
                      'STEP-1 =', fstrDYNAMIC%i_step-1,                       &
                      'TIME-dt =', fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step-1) 
         
        END IF
        
       END IF
       
       WRITE(ILOG, *) '#### VELOCITY SHELL'
       WRITE(ILOG, '(A,A)')                                          &
                   '     NODE    X-VELO      Y-VELO      Z-VELO   ', &
                   '   X-ROT-VEL   Y-ROT-VEL   Z-ROT-VEL'            
       WRITE(ILOG, '(A,A)')                                          &
                   ' --------+-----------+-----------+-----------+', &
                   '-----------+-----------+-------------'           
       
       DO i= 1, hecMESH%nn_internal
        
        jj = fstrPARAM%global_local_id(1, i)
        ii = fstrPARAM%global_local_id(2, i)
        
        WRITE( ILOG, '(I10, 6E12.4)' )                                    &
               jj, ( fstrDYNAMIC%VEL( ndof*(ii-1)+k, idx ), k = 1, ndof ) 
        
       END DO
       
      END IF
      
!--------------------------------------------------------------------
      
      ! Acceleration
      IF( fstrDYNAMIC%iout_list(3) .EQ. 1 ) THEN
       
       IF( fstrDYNAMIC%idx_eqa .EQ. 11 ) THEN
        
        IF( fstrDYNAMIC%i_step .GT. 0 ) THEN
         
         WRITE( ILOG, '(A, I10, 5X, A, 1PE12.4)')                       &
                'STEP-1 =' , fstrDYNAMIC%i_step-1,                      &
                'TIME-dt =', fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step-1) 
         
        END IF
        
       END IF
       
       WRITE(ILOG, *) '#### ACCELERATION SHELL'
       WRITE(ILOG, '(A, A)')                                         &
                   '     NODE    X-ACCE      Y-ACCE      Z-ACCE   ', &
                   '   X-ROT-ACC   Y-ROT-ACC   Z-ROT-ACC'            
       WRITE(ILOG, '(A, A)')                                         &
                   ' --------+-----------+-----------+-----------+', &
                   '-----------+-----------+-------------'           
       
       DO i= 1, hecMESH%nn_internal
        
        jj = fstrPARAM%global_local_id(1, i)
        ii = fstrPARAM%global_local_id(2, i)
        
        WRITE( ILOG, '(I10, 6E12.4)' )                                    &
               jj, ( fstrDYNAMIC%ACC( ndof*(ii-1)+k, idx ), k = 1, ndof ) 
        
       END DO
       
      END IF
      
!--------------------------------------------------------------------
      
      ! Reaction force
      WRITE(ILOG, *) '#### REACTION FORCE SHELL'
      
      CALL dynamic_reaction_force_6d_Shell   &
           (hecMESH, fstrSOLID, fstrDYNAMIC) 
      
      CALL flush(ILOG)
      
!--------------------------------------------------------------------
      
      ! Nodal stress
      WRITE(ILOG, *) '#### NODAL STRAIN/STRESS SHELL'
      
      CALL fstr_nodal_stress_6d_Shell                &
           (hecMESH, fstrSOLID, fstrPARAM, fstrDYNAMIC%iout_list) 
      
      CALL flush(ILOG)
      
!--------------------------------------------------------------------
      
      RETURN
      
!####################################################################
      END SUBROUTINE dynamic_output_6d_Shell
!####################################################################
      ! > (Gaku Hashimoto, The University of Tokyo, 2012/11/15)
      
      
      ! (Gaku Hashimoto, The University of Tokyo, 2012/11/15) <
!####################################################################
      SUBROUTINE dynamic_reaction_force_6d_Shell   &
                 (hecMESH, fstrSOLID, fstrDYNAMIC) 
!####################################################################
      
      IMPLICIT NONE
      
!--------------------------------------------------------------------
      
      TYPE(hecmwST_local_mesh) :: hecMESH
      TYPE(fstr_solid)         :: fstrSOLID
      TYPE(fstr_dynamic)       :: fstrDYNAMIC
      
!--------------------------------------------------------------------
      
      REAL(KIND = kreal) :: edisp(6*9), force(6*9)
      REAL(KIND = kreal), ALLOCATABLE   :: spcForce(:)
      REAL(KIND = kreal) :: ecoord(3, 9)
      REAL(KIND = kreal) :: stiff(6*9, 6*9)
      REAL(KIND = kreal) :: thick
      
      INTEGER(KIND = kint) :: i, j, k
      INTEGER(KIND = kint) :: ndof
      INTEGER(KIND = kint) :: nodLOCAL(9)
      INTEGER(KIND = kint) :: icount
      INTEGER(KIND = kint) :: ig, ig0
      INTEGER(KIND = kint) :: itype, ic_type, iS, iE, iS0, iE0
      INTEGER(KIND = kint) :: ik, kk
      INTEGER(KIND = kint) :: in, isid
      INTEGER(KIND = kint) :: nspc
      INTEGER(KIND = kint) :: nn
      INTEGER(KIND = kint) :: icel, jS
      INTEGER(KIND = kint) :: ifspc
      INTEGER(KIND = kint) :: isect, ihead
      INTEGER(KIND = kint), ALLOCATABLE :: id_spc(:)
      
!--------------------------------------------------------------------
      
      ALLOCATE( id_spc( hecMESH%n_node ) )
      
      id_spc = 0
      
!--------------------------------------------------------------------
      
      ndof = hecMESH%n_dof
      
!--------------------------------------------------------------------
      
      ! SPC Nodal Check Flag
      icount = 0
      
      DO ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
       
       ig = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
       iS = hecMESH%node_group%grp_index(ig-1) + 1
       iE = hecMESH%node_group%grp_index(ig  )
       
       DO ik = iS, iE
        
        in = hecMESH%node_group%grp_item(ik)
        
        IF( id_spc(in) .EQ. 0 ) THEN
         
         icount = icount+1
         
         id_spc(in) = icount
         
        END IF
        
       END DO
       
      END DO
      
      icount = 0
      
      DO i = 1, hecMESH%n_node
       
       IF( id_spc(i) .GE. 1 ) icount = icount+1
       
      END DO
      
      nspc = icount
      
      ! Allocate array
      IF( nspc .GT. 0 ) THEN
       
       ALLOCATE( spcForce(nspc*ndof) )
       
       spcForce = 0.0D0
       
      ELSE
       
       DEALLOCATE ( id_spc )
       
       RETURN
       
      END IF
      
!--------------------------------------------------------------------
      
      ! +-------------------------------+
      ! | according to ELEMENT TYPE     |
      ! +-------------------------------+
      DO itype = 1, hecMESH%n_elem_type
       
       !--------------------------------------------------------
       
       iS = hecMESH%elem_type_index(itype-1)+1
       iE = hecMESH%elem_type_index(itype  )
       
       ic_type = hecMESH%elem_type_item(itype)
       
       !--------------------------------------------------------
       
       IF( .NOT. hecmw_is_etype_shell(ic_type) ) CYCLE
       
       !--------------------------------------------------------
       
       !** Set number of nodes
       nn = hecmw_get_max_node(ic_type)
       
       !--------------------------------------------------------
       
       DO icel = iS, iE
        
        jS = hecMESH%elem_node_index(icel-1)
        
        !--------------------------------------------------
        
        DO j = 1, nn
         
         nodLOCAL(j) = hecMESH%elem_node_item(jS+j)
         
         ecoord(1, j) = hecMESH%node( 3*( nodLOCAL(j)-1 )+1 )
         ecoord(2, j) = hecMESH%node( 3*( nodLOCAL(j)-1 )+2 )
         ecoord(3, j) = hecMESH%node( 3*( nodLOCAL(j)-1 )+3 )
         
         DO k = 1, ndof
          
          edisp( ndof*(j-1)+k ) = fstrSOLID%unode( ndof*( nodLOCAL(j)-1 )+k )
          
         END DO
         
        END DO
        
        !--------------------------------------------------
        
        ifspc = 0
        
        DO j = 1, nn
         
         IF( id_spc( nodLOCAL(j) ) .GE. 1 ) ifspc = 1
         
        END DO
        
        !--------------------------------------------------
        
        IF( ifspc .EQ. 1 ) THEN
         
         !--------------------------------------------
         
         isect = hecMESH%section_ID(icel)
         ihead = hecMESH%section%sect_R_index(isect-1)
         thick = hecMESH%section%sect_R_item(ihead+1)
         
         !--------------------------------------------
         
         IF( ( ic_type .EQ. 741 ) .OR. ( ic_type .EQ. 743 ) .OR. &
             ( ic_type .EQ. 731 ) ) THEN                         
          
          CALL STF_Shell_MITC                                                     &
               ( ic_type, nn, ndof, ecoord(:, 1:nn),                              &
                 fstrSOLID%elements(icel)%gausses, stiff(1:6*nn, 1:6*nn), thick ) 
          
         END IF
         
         !--------------------------------------------
         
         force(1:ndof*nn)= MATMUL( stiff(1:ndof*nn,1:ndof*nn), edisp(1:ndof*nn) )
         
         !--------------------------------------------
         
         !*** Add SPC FORCE ****
         DO j = 1, nn
          
          isid = id_spc( nodLOCAL(j) )
          
          IF( isid .GT. 0 ) THEN
           
           DO k = 1, ndof
            
            kk = ndof*(isid-1)+k
            
            spcForce(kk) = spcForce(kk) + force( ndof*(j-1)+k )
            
           END DO
           
          END IF
          
         END DO
         
         !--------------------------------------------
         
        END IF
        
        !--------------------------------------------------
        
       END DO
       
       !--------------------------------------------------------
       
      END DO
      
!--------------------------------------------------------------------
      
      IF( fstrDYNAMIC%iout_list(4) .EQ. 1 ) THEN
       
       WRITE(ILOG, '(A, A)')                                         &
                   '   NODE      X-REAC      Y-REAC      Z-REAC   ', &
                   '  RX-REAC     RY-REAC     RZ-REAC'               
       WRITE(ILOG, '(A, A)')                                         &
                   ' --------+-----------+-----------+-----------+', &
                   '-----------+-----------+-------------'           
       
       DO i = 1, hecMESH%nn_internal
        
        in = hecMESH%global_node_ID(i)
        
        isid = id_spc(i)
        
        IF( isid .GT. 0 ) THEN
         
         WRITE( ILOG, '(I10, 1P6E12.4)' )                        &
                in, ( spcForce( ndof*(isid-1)+k ), k = 1, ndof ) 
         
        END IF
        
       END DO
       
      END IF
      
!--------------------------------------------------------------------
      
      DEALLOCATE( id_spc )
      DEALLOCATE( spcForce )
      
!--------------------------------------------------------------------
      
      RETURN
	  
!####################################################################
      END SUBROUTINE dynamic_reaction_force_6d_Shell
!####################################################################
      

!C 
!C OUTPUT FOR ELEMENTS
!C
      subroutine dynamic_output_elem(hecMESH,fstrSOLID,fstrPARAM, fstrDYNAMIC)
      type (hecmwST_local_mesh) :: hecMESH
      type (fstr_solid)      :: fstrSOLID
      type (fstr_param)      :: fstrPARAM
      type (fstr_dynamic      ) :: fstrDYNAMIC
      integer :: ig, ig0, iS0, iE0, ik, iS, iE, ic_type, nn, jS, icel, isid, isect, num
      integer :: i,j,k,kk, in, itype, ielem, ID_area, ii,jj, NDOF, icount, nspc, ifspc

      if( hecMESH%n_dof .eq. 3 ) then

        if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then

        write(ILOG,*) '#### 3D STRAIN @Element'
        write(ILOG,'(a,a)')                                           &
                    '  ELEMENT     E11         E22         E33     '  &
                             ,'    E12         E23         E13     '
        write(ILOG,'(a,a)')                                            &
                    '  -------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p6e12.4)') jj,(fstrSOLID%ESTRAIN(6*(ii-1)+k),k=1,6)
          end if
        enddo
        end if

        if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then

        write(ILOG,*) '#### 3D STRESS @Element'
        write(ILOG,'(a,a)')                                           &
                    '  ELEMENT     S11         S22         S33     '  &
                 ,'    S12         S23         S13        MISES'
        write(ILOG,'(a,a)')                                           &
                    '  -------+-----------+-----------+-----------+'  &
                 ,'-----------+-----------+-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p7e12.4)') jj,(fstrSOLID%ESTRESS(7*(ii-1)+k),k=1,7)
          end if
        enddo
        end if

      elseif( hecMESH%n_dof .eq. 2 ) then

        if( fstrDYNAMIC%iout_list(5) .eq. 1 ) then
        write(ILOG,*) '#### 2D STRAIN @Element'
        write(ILOG,'(a,a)')                                      &
               'ELEMENT       E11         E22         E33     '  &
                             ,'    E12' 
        write(ILOG,'(a,a)')                                           &
               '  ------------+-----------+-----------+-----------+'  &
                             ,'-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p4e12.4)') jj,(fstrSOLID%ESTRAIN(6*(ii-1)+k),k=1,4)
          end if
        enddo
        end if

        if( fstrDYNAMIC%iout_list(6) .eq. 1 ) then
        write(ILOG,*) '#### 2D STRESS @Element'
        write(ILOG,'(a,a)')                                     &
               'ELEMENT       S11         S22         S33     ' &
                             ,'    S12        MISES' 
        write(ILOG,'(a,a)')                                            &
               '  ------------+-----------+-----------+-----------+'   &
                             ,'-----------+-------------'
        do i=1,hecMESH%n_elem
!!        ii = hecMESH%elem_ID(i*2-1)
          ii = i
          jj = hecMESH%global_elem_ID(i)
          ID_area = hecMESH%elem_ID(i*2)
          if( ID_area.eq.hecMESH%my_rank ) then
            write(ILOG,'(i10,1p5e12.4)') jj,(fstrSOLID%ESTRESS(7*(ii-1)+k),k=1,4)  &
                                   ,fstrSOLID%ESTRESS(7*(ii-1)+7)
          end if
        enddo
        end if

      elseif( hecMESH%n_dof .eq. 6 ) then

        IF( fstrDYNAMIC%iout_list(5) .EQ. 1 ) THEN
         
         WRITE( ILOG, '(a,a)' )                                              &
                '  ------------+-----------+-----------+-----------+',       &
                              '-----------+-----------+-------------'        
         
         WRITE(ILOG, *) '#### SHELL STRAIN @Element'
         WRITE( ILOG, '(A, A)' )                                             &
                '  ELEMENT          E11         E22         E33     ',       &
                               '    E12         E23         E31       EQUIV' 
         WRITE( ILOG, '(A, A)' )                                             &
                '  ------------+-----------+-----------+-----------+',       &
                              '-----------+-----------+-------------'        
         
         DO i = 1, hecMESH%n_elem
          
          ii = i
          jj = hecMESH%global_elem_ID(i)
          
          ID_area = hecMESH%elem_ID(i*2)
          
          IF( ID_area .EQ. hecMESH%my_rank ) then
           
           WRITE( ILOG, '(I10, A5, 1P7E12.4)' )                                &
                  jj, ' (+) ', ( fstrSOLID%ESTRAIN( 14*(ii-1)+k ), k = 1, 7  ) 
           WRITE( ILOG, '(10X, A5, 1P7E12.4)' )                                &
                      ' (-) ', ( fstrSOLID%ESTRAIN( 14*(ii-1)+k ), k = 8, 14 ) 
           
          END IF
          
         END DO
         
        END IF

        IF( fstrDYNAMIC%iout_list(6) .EQ. 1 ) THEN
         
         WRITE(ILOG, *) '#### SHELL STRESS @Element'
         WRITE( ILOG, '(A, A)' )                                             &
                '  ELEMENT          S11         S22         S33     ',       &
                               '    S12         S23         S31       MISES' 
         WRITE( ILOG, '(A, A)' )                                             &
                '  ------------+-----------+-----------+-----------+',       &
                              '-----------+-----------+-------------'        
       
         
         DO i = 1, hecMESH%n_elem
          
          ii = i
          jj = hecMESH%global_elem_ID(i)
          
          ID_area = hecMESH%elem_ID(i*2)
          
          IF( ID_area .EQ. hecMESH%my_rank ) THEN
           
           WRITE( ILOG, '(I10, A5, 1P7E12.4)' )                              &
                  jj, ' (+) ', ( fstrSOLID%ESTRESS(14*(ii-1)+k), k = 1, 7  ) 
           WRITE( ILOG, '(10X, A5, 1P7E12.4)' )                              &
                      ' (-) ', ( fstrSOLID%ESTRESS(14*(ii-1)+k), k = 8, 14 ) 
           
          END IF
          
         END DO
         
        end if

      endif

      end subroutine dynamic_output_elem
	  
!C================================================================C
!C-- subroutine  matvec
!C================================================================C
  subroutine matvec(y,x,hecMAT,ndof,D,AU,AL)
    use hecmw
    type ( hecmwST_matrix      ) :: hecMAT

    integer(kind=kint) :: ndof,i,is,ie,j,icol,ki,kj,ix,iy,ip,nn
    real(kind=kreal) :: D(ndof*ndof*hecMAT%NP)
    real(kind=kreal) :: AU(ndof*ndof*hecMAT%NPU)
    real(kind=kreal) :: AL(ndof*ndof*hecMAT%NPL)
    real(kind=kreal) :: x(ndof*hecMAT%NP)
    real(kind=kreal) :: y(ndof*hecMAT%NP)

    nn=ndof*ndof

    y=0.

    do i=1,hecMAT%NP
      is=hecMAT%indexU(i-1)+1
      ie=hecMAT%indexU(i)
      do j=is,ie
        icol=hecMAT%itemU(j)
        do ki=1,ndof
          iy=ndof*(i-1)+ki
          do kj=1,ndof
            ix=ndof*(icol-1)+kj
            ip=nn*(j-1)+ndof*(ki-1)+kj
            y(iy)=y(iy)+AU(ip)*x(ix)
          enddo
        enddo
      enddo
    enddo

    do i=1,hecMAT%NP
      is=hecMAT%indexL(i-1)+1
      ie=hecMAT%indexL(i)
      do j=is,ie
        icol=hecMAT%itemL(j)
        do ki=1,ndof
          iy=ndof*(i-1)+ki
          do kj=1,ndof
            ix=ndof*(icol-1)+kj
            ip=nn*(j-1)+ndof*(ki-1)+kj
            y(iy)=y(iy)+AL(ip)*x(ix)
          enddo
        enddo
      enddo
    enddo

    do i=1,hecMAT%NP
      do ki=1,ndof
        iy=ndof*(i-1)+ki
        do kj=1,ndof
          ix=ndof*(i-1)+kj
          ip=nn*(i-1)+ndof*(ki-1)+kj
          y(iy)=y(iy)+D(ip)*x(ix)
        enddo
      enddo
    enddo

  end subroutine matvec

  subroutine dynamic_output_and_post(istep,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      use hecmw

      integer (kind = kint) :: istep, interval
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM
      type ( fstr_dynamic        ) :: fstrDYNAMIC

!
!C-- output new displacement, velocity and accelaration
      if( IRESULT.eq.1 .and. &
          (mod(istep,fstrSOLID%output_ctrl(3)%freqency).eq.0 .or. istep.eq.fstrDYNAMIC%n_step) ) then
!
      write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP =',istep,'TIME =',fstrDYNAMIC%t_delta*istep
!C-- output

      call flush(ILOG)
      IF(myrank .EQ. 0) THEN
        IF(hecMAT%Iarray(99) .EQ. 1) THEN
          write(IMSG,*) '*----------------------------*'
          write(IMSG,*) '## No.of ITER:',hecMAT%Iarray(1)
          write(IMSG,*) '*----------------------------*'
        ENDIF
        call flush(IDBG)
        call flush(IMSG)
      ENDIF

!C-- output
      if( hecMESH%n_dof .eq. 3 ) then
        call dynamic_output_3d(hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrDYNAMIC)
      else if( hecMESH%n_dof .eq. 2 ) then
        call dynamic_output_2d(hecMESH,hecMAT,fstrSOLID,fstrPARAM,fstrDYNAMIC)
      else if( hecMESH%n_dof.EQ.6) THEN
        call dynamic_output_6d_Shell(hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC)
      endif

!C-- OUTPUT FOR ELEMENTS
      call dynamic_output_elem(hecMESH,fstrSOLID,fstrPARAM,fstrDYNAMIC)
!
      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
      call flush(ILOG)
!C
!C-- POST PROCESSING
!C
      call dynamic_post(hecMESH,hecMAT,fstrSOLID,istep,fstrDYNAMIC)

      endif

!C
!C-- POST PROCESSING VIA MEMORY
!C
      if( IVISUAL.eq.1 .and. &
          (mod(istep,fstrSOLID%output_ctrl(4)%freqency).eq.0 .or. istep.eq.fstrDYNAMIC%n_step) ) then
        interval = fstrSOLID%output_ctrl(4)%freqency
        call fstr_make_result(hecMESH,fstrSOLID,fstrRESULT)
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize(hecMESH,fstrRESULT,istep,fstrDYNAMIC%n_step,interval)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif

  end subroutine dynamic_output_and_post


!C================================================================C
!C-- subroutine dynamic_output_monit
!C================================================================C
  subroutine dynamic_output_monit(hecMESH, fstrPARAM, fstrDYNAMIC, myEIG, my_rank_monit_1)
      use hecmw
!C
!C-- global variable
!C
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( fstr_dynamic        ) :: fstrDYNAMIC
      type ( fstr_param          ) :: fstrPARAM
      type ( lczparam            ) :: myEIG         
!C
!C-- local variable
!C
    integer(kind=kint) :: ii,jj,ierr                
    integer(kind=kint) :: my_rank_monit_1
    logical :: yes                                      

    if( mod(fstrDYNAMIC%i_step,fstrDYNAMIC%nout_monit) /= 0 ) return
!!
    if( hecMESH%my_rank .eq. my_rank_monit_1) then
!
      jj=fstrPARAM%global_local_id(1,fstrDYNAMIC%node_monit_1)
      ii=fstrPARAM%global_local_id(2,fstrDYNAMIC%node_monit_1)
!!
!!    step = 0
!!
      if(fstrDYNAMIC%i_step .eq. 0) then
!C-- displacement
        if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
!!           write(dynamic_IW4,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW4,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &  
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )      
        end if
!C-- velocity
        if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
!!           write(dynamic_IW5,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW5,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &   
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%VEL ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )                 
        end if
!C-- acceleration
        if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
!!           write(dynamic_IW6,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW6,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &   
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%ACC ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!!
!!    step = 1,2,.......
!!
      else if(fstrDYNAMIC%i_step .gt. 0) then
!
!C-- implicit dynamic analysis
        if(fstrDYNAMIC%idx_eqa .eq. 1) then
!C-- displacement
        if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
!!           write(dynamic_IW4,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW4,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &   
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )         
        end if
!C-- velocity
        if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
!!           write(dynamic_IW5,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW5,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &   
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%VEL ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )    
        end if
!C-- acceleration
        if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
!!           write(dynamic_IW6,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW6,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &   
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%ACC ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )
        end if
!
!C-- explicit dynamic analysis
        else if(fstrDYNAMIC%idx_eqa .eq. 11) then
!C-- displacement
        if( fstrDYNAMIC%iout_list(1) .eq. 1 ) then
!!           write(dynamic_IW4,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step &
           write(fstrDYNAMIC%dynamic_IW4,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step &  
                 ,fstrDYNAMIC%t_curr &
                 ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
!                ,jj ,fstrDYNAMIC%DISP( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 2 )

        end if
!C-- velocity
        if( fstrDYNAMIC%iout_list(2) .eq. 1 ) then
!!           write(dynamic_IW5,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step -1 &
           write(fstrDYNAMIC%dynamic_IW5,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step -1 &  
                 ,fstrDYNAMIC%t_curr-fstrDYNAMIC%t_delta &
                 ,jj ,fstrDYNAMIC%VEL ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!C-- acceleration
        if( fstrDYNAMIC%iout_list(3) .eq. 1 ) then
!!           write(dynamic_IW6,'(i10,1pe12.4,i10,1p6e12.4)') fstrDYNAMIC%i_step -1 &
           write(fstrDYNAMIC%dynamic_IW6,'(i10,1pe13.4e3,i10,1p6e13.4e3)') fstrDYNAMIC%i_step -1 &  
                 ,fstrDYNAMIC%t_curr-fstrDYNAMIC%t_delta &
                 ,jj ,fstrDYNAMIC%ACC ( hecMESH%n_dof*(ii-1)+1 : hecMESH%n_dof*ii , 1 )
        end if
!
        end if
!
      end if
    end if
    if( hecMESH%my_rank == 0 )then     
      if(any(fstrDYNAMIC%iout_list(1:3) == 1)) then 
        inquire(file='dyna_energy.txt',opened=yes)
        if( .not. yes )then
          open(fstrDYNAMIC%dynamic_IW7, file = 'dyna_energy.txt', status = 'replace', iostat=ierr )
          if( ierr /= 0 ) then
            write(*,*) 'stop due to file opening error <dyna_enrgy.txt>'
            call hecmw_abort( hecmw_comm_get_comm()) 
          endif 
          write(fstrDYNAMIC%dynamic_IW7,*) ' time step', '     time    ',   &
                                       '  kinetic energy', '   strain energy', '   total energy'     
        endif
        if(fstrDYNAMIC%i_step == 0) then
           fstrDYNAMIC%kineticEnergy=0.0d0
           do ii=1,hecMESH%n_node*hecMESH%n_dof
             fstrDYNAMIC%kineticEnergy=fstrDYNAMIC%kineticEnergy+&
                                       0.5d0*myEIG%mass(ii)*fstrDYNAMIC%VEL(ii,1)*fstrDYNAMIC%VEL(ii,1)  
           enddo      
        endif
        fstrDYNAMIC%totalEnergy=fstrDYNAMIC%kineticEnergy+fstrDYNAMIC%strainEnergy    
        write(fstrDYNAMIC%dynamic_IW7,'(i10,1pe13.4e3,1p3e16.4e3)') fstrDYNAMIC%i_step &       
              ,fstrDYNAMIC%t_curr &
              ,fstrDYNAMIC%kineticEnergy,fstrDYNAMIC%strainEnergy,fstrDYNAMIC%totalEnergy     
      endif 
      if(fstrDYNAMIC%i_step==fstrDYNAMIC%n_step)close(fstrDYNAMIC%dynamic_IW7)
    endif                              
	
	call flush(fstrDYNAMIC%dynamic_IW4)

  end subroutine dynamic_output_monit
  
  subroutine dynamic_nloutput(istep,hecMESH,hecMAT,fstrSOLID,fstrRESULT,fstrPARAM,fstrDYNAMIC)
      use hecmw

      integer (kind = kint) :: istep, interval
      type ( hecmwST_local_mesh  ) :: hecMESH
      type ( hecmwST_matrix      ) :: hecMAT
      type ( fstr_solid          ) :: fstrSOLID
      type ( hecmwST_result_data ) :: fstrRESULT
      type ( fstr_param          ) :: fstrPARAM
      type ( fstr_dynamic        ) :: fstrDYNAMIC

      integer :: idx
      integer :: ndof, i, ii,jj, k, ID_area, itype, iS, iE, ic_type, icel
      real(kind=kreal) :: s11,s22,s33,s12,s23,s13,ps,smises,ss(6)
!
!C-- output new displacement, velocity and accelaration
      if( IRESULT.eq.1 .and. &
          (mod(istep,fstrSOLID%output_ctrl(3)%freqency).eq.0 .or. istep.eq.fstrDYNAMIC%n_step) ) then
!
      if((fstrDYNAMIC%idx_eqa .eq. 1) .and. (fstrDYNAMIC%i_step .gt. 0)) then
        idx = 2
      else
        idx = 1
      endif
!
      write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP =',istep,'TIME =',fstrDYNAMIC%t_curr   
!C-- output

      call flush(ILOG)
      IF(myrank == 0) THEN
        IF(hecMAT%Iarray(99) == 1) THEN
          write(IMSG,*) '*----------------------------*'
          write(IMSG,*) '## No.of ITER:',hecMAT%Iarray(1)
          write(IMSG,*) '*----------------------------*'
        ENDIF
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
	  
      ndof = hecMESH%n_dof
      
!C
!C-- DISPLACEMENT
!C
      if( fstrDYNAMIC%iout_list(1) == 1 ) then
      if( ndof == 3 ) then
         write(ILOG,*) '#### DISPLACEMENT 3D'
         write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP      Z-DISP'
         write(ILOG,'(a)')'---------+-----------+-----------+------------'
      else if( ndof == 2 ) then
         write(ILOG,*) '#### DISPLACEMENT 2D'
         write(ILOG,'(a)')'    NODE     X-DISP      Y-DISP    '
         write(ILOG,'(a)')'---------+-----------+-----------'
      else if( ndof == 6 ) then
         write(ILOG,*) '#### DISPLACEMENT SHELL'
         write(ILOG,'(a,a)')                                              &
                    '     NODE    X-DISP      Y-DISP      Z-DISP   '   &
                             ,'   X-ROT       Y-ROT       Z-ROT'
         write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      endif

      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

!        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%DISP(ndof*(ii-1)+k, 1),k=1,ndof)   
        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%DISP(ndof*(ii-1)+k, idx),k=1,ndof)           
      enddo
      end if
!C
!C-- VELOCITY   
!C
      if( fstrDYNAMIC%iout_list(2) ==1 ) then
        if(fstrDYNAMIC%i_step > 0) then
 !         write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &    
 !               ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if

      if( ndof == 3 ) then
        write(ILOG,*) '#### VELOCITY 3D'
        write(ILOG,'(a)')'    NODE     X-VELO      Y-VELO      Z-VELO'
        write(ILOG,'(a)')'---------+-----------+-----------+------------'
      else if( ndof == 6 ) then
        write(ILOG,*) '#### VELOCITY SHELL'
        write(ILOG,'(a,a)')                                              &
                    '     NODE    X-VELO      Y-VELO      Z-VELO   '   &
                             ,'   X-ROT-VEL   Y-ROT-VEL   Z-ROT-VEL'
        write(ILOG,'(a,a)')                                              &
                    ' --------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
      endif
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

!        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%VEL(ndof*(ii-1)+k, 1),k=1,ndof)       
        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%VEL(ndof*(ii-1)+k, idx),k=1,ndof)         
      enddo

      end if
!C
!C-- ACCELERATION
!C
      if( fstrDYNAMIC%iout_list(3) == 1 ) then
        if(fstrDYNAMIC%i_step .gt. 0) then
!          write(ILOG,'(a,i10,5x,a,1pe12.4)') 'STEP-1 =' ,fstrDYNAMIC%i_step-1   &    
!                ,'TIME-dt =',fstrDYNAMIC%t_delta*(fstrDYNAMIC%i_step -1)
        end if

      write(ILOG,*) '#### ACCELERATION 3D'
      write(ILOG,'(a)')'    NODE     X-ACCE      Y-ACCE      Z-ACCE'
      write(ILOG,'(a)')'---------+-----------+-----------+------------'
      do i= 1, hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)

!        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%ACC(ndof*(ii-1)+k, 1),k=1,ndof)       
        write(ILOG,'(i10,3e12.4)') jj,(fstrDYNAMIC%ACC(ndof*(ii-1)+k, idx),k=1,ndof)           
      enddo

      end if

!---  Reaction force	  
      if( fstrDYNAMIC%iout_list(4) == 1 ) then
	  
      if( ndof == 3 ) then
               write(ILOG,'(a)') '#### REACTION FORCE 3D'
               write(ILOG,'(a)') '    NODE      X-DIREC       Y-DIREC         Z-DIREC    '
               write(ILOG,'(a)') '---------+--------------+--------------+---------------'
      else if( ndof == 2 ) then
               write(ILOG,'(a)') '#### REACTION FORCE 2D'
               write(ILOG,'(a)') '    NODE     X-DIREC        Y-DIREC     '
               write(ILOG,'(a)') '---------+--------------+---------------'
      else if( ndof == 6 ) then
      endif

      do i=1,hecMESH%nn_internal
        jj=fstrPARAM%global_local_id(1,i)
        ii=fstrPARAM%global_local_id(2,i)                                               
        write(ILOG,'(i10,1p6e15.7)') jj, fstrSOLID%QFORCE( ndof*(ii-1)+1 : ndof*(ii-1)+ndof)
      enddo
      endif

!C-- OUTPUT FOR ELEMENTS
      if( fstrDYNAMIC%iout_list(5) == 1 ) then
        write(ILOG,*) '#### 3D STRAIN @Element'
        write(ILOG,'(a,a)')                                           &
                    '  ELEMENT     E11         E22         E33     '  &
                             ,'    E12         E23         E13     '
        write(ILOG,'(a,a)')                                            &
                    '  -------+-----------+-----------+-----------+'   &
                             ,'-----------+-----------+-------------'
        do itype= 1, hecMESH%n_elem_type
           iS= hecMESH%elem_type_index(itype-1) + 1
           iE= hecMESH%elem_type_index(itype  )
           ic_type= hecMESH%elem_type_item(itype)

           if (hecmw_is_etype_link(ic_type)) cycle

           ii = NumOfQuadPoints( ic_type )
           do icel= iS, iE
	          ss(:) = 0.d0
              do jj=1,ii
                 ss(:) = ss(:) + fstrSOLID%elements(icel)%gausses(jj)%strain(:)
              enddo
              ss(:) = ss(:)/ii
              jj = hecMESH%global_elem_ID(icel)
              write(ILOG,'(i10,1p6e12.4)') jj,(ss(k),k=1,6)
           enddo
        enddo
      end if

      if( fstrDYNAMIC%iout_list(6) == 1 ) then

        write(ILOG,*) '#### 3D STRESS @Element'
        write(ILOG,'(a,a)')                                           &
                    '  ELEMENT     S11         S22         S33     '  &
                 ,'    S12         S23         S13        MISES'
        write(ILOG,'(a,a)')                                           &
                    '  -------+-----------+-----------+-----------+'  &
                 ,'-----------+-----------+-----------+-------------'
        do itype= 1, hecMESH%n_elem_type
           iS= hecMESH%elem_type_index(itype-1) + 1
           iE= hecMESH%elem_type_index(itype  )
           ic_type= hecMESH%elem_type_item(itype)

           if (hecmw_is_etype_link(ic_type)) cycle

           ii = NumOfQuadPoints( ic_type )
           do icel= iS, iE
	          ss(:) = 0.d0
              do jj=1,ii
                 ss(:) = ss(:) + fstrSOLID%elements(icel)%gausses(jj)%stress(:)
              enddo
              ss(:) = ss(:)/ii
              s11=ss(1)
              s22=ss(2)
              s33=ss(3)
              s12=ss(4)
              s23=ss(5)
              s13=ss(6)
              ps=(s11+s22+s33)/3.d0
              smises=0.5d0*( (s11-ps)**2+(s22-ps)**2+(s33-ps)**2 )    &
                            +s12**2+s23**2+s13**2 
              smises=dsqrt(3.d0*smises)
              jj = hecMESH%global_elem_ID(icel)
              write(ILOG,'(i10,1p7e12.4)') jj,(ss(k),k=1,6),smises
           enddo
        enddo
      end if

      IF(myrank == 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
        call flush(IDBG)
        call flush(IMSG)
      ENDIF
      call flush(ILOG)
!C
!C-- POST PROCESSING
!C
      call dynamic_post(hecMESH,hecMAT,fstrSOLID,istep,fstrDYNAMIC)

      endif

!C
!C-- POST PROCESSING VIA MEMORY
!C
      if( IVISUAL.eq.1 .and. &
          (mod(istep,fstrSOLID%output_ctrl(4)%freqency).eq.0 .or. istep.eq.fstrDYNAMIC%n_step) ) then
        interval = fstrSOLID%output_ctrl(4)%freqency
        call fstr_make_result(hecMESH,fstrSOLID,fstrRESULT)
        call fstr2hecmw_mesh_conv(hecMESH)
        call hecmw_visualize_init
        call hecmw_visualize(hecMESH,fstrRESULT,istep,fstrDYNAMIC%n_step,interval)
        call hecmw_visualize_finalize
        call hecmw2fstr_mesh_conv(hecMESH)
        call hecmw_result_free(fstrRESULT)
      endif

  end subroutine dynamic_nloutput

end module m_dynamic_output
