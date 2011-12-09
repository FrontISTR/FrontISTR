!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.0                                   !
!                                                                      !
!     Last Update : 2007/03/02                                         !
!        Category : Dynamic Transit Analysis                           !
!                                                                      !
!                    Written by Noboru Imai (Univ. of Tokyo)           !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

module m_fstr_rcap_io
use m_fstr
use rcapf

      public :: fstr_rcap_initialize ! call after fstr_setup
      public :: fstr_rcap_finalize   ! call before hecmw_finalize
      public :: fstr_rcap_send
      public :: fstr_rcap_get

contains

!------------------------------------------------------------------------------
subroutine fstr_rcap_initialize( hecMESH, fstrPARAM, fstrCPL )
      implicit none
      type( hecmwST_local_mesh ) :: hecMESH
      type( fstr_param  ) :: fstrPARAM
      type( fstr_couple ) :: fstrCPL
      character( len=16)  :: portfile
      integer(kind=kint)  :: myrank
      integer(kind=kint)  :: i,err,nid,nid_old

      write(IDBG,*) "fstr_rcap_initialize: start"

      if( fstrPARAM%fg_couple /= 1 ) return

      write(IDBG,*) "fstr_rcap_initialize: calling rcapf_init_solid_solver"

      portfile='port'//CHAR(0)

      call rcapf_init_solid_solver( hecMESH%my_rank, portfile )

      write(IDBG,*) "fstr_rcap_initialize: calling rcapf_get_num_of_matching_node"

      fstrCPL%dof = 3

      call rcapf_get_num_of_matching_node( fstrCPL%coupled_node_n )

      fstrCPL%ndof = fstrCPL%coupled_node_n * fstrCPL%dof
      allocate( fstrCPL%coupled_node( fstrCPL%coupled_node_n ), stat=err)
      if( err /= 0 ) goto 1000
      allocate( fstrCPL%trac( fstrCPL%ndof ), stat=err)
      if( err /= 0 ) goto 1000
      allocate( fstrCPL%disp( fstrCPL%ndof ), stat=err)
      if( err /= 0 ) goto 1000
      allocate( fstrCPL%velo( fstrCPL%ndof ), stat=err)
      if( err /= 0 ) goto 1000
      allocate( fstrCPL%accel( fstrCPL%ndof ), stat=err)
      if( err /= 0 ) goto 1000
      allocate( fstrCPL%index( hecMESH%n_node_gross ), stat=err)
      if( err /= 0 ) goto 1000
      fstrCPL%index = -1

      write(IDBG,*) "fstr_rcap_initialize: calling rcapf_get_mathcing_node_id"

      call rcapf_get_matching_node_id( fstrCPL%coupled_node, fstrCPL%coupled_node_n )

      write(IDBG,*) "fstr_rcap_initialize: converting to local id: ", fstrCPL%coupled_node_n

      do i=1, fstrCPL%coupled_node_n
            nid = fstrCPL%coupled_node(i)
            if( nid <= 0 .or. nid>hecMESH%n_node_gross ) then
                  write(*,*) '##Fatal error in fstr_rcap_initialize ', i, nid
                  call hecmw_abort( hecmw_comm_get_comm());
            endif
!            if( hecMESH%n_refine > 0 ) then
!                  nid_old = nid
!                  if( associated( hecMESH%node_old2new ) ) then
!                    nid = hecMESH%node_old2new( nid ) + 1
!                  endif
!                  write(IDBG,*) nid_old, nid
!            endif
            fstrCPL%index( nid ) = i
      end do

      write(IDBG,*) "fstr_rcap_initialize: end"

      return
1000  write(*,*) "##Error : not enough memory"
      call hecmw_abort( hecmw_comm_get_comm() )
end subroutine fstr_rcap_initialize

!------------------------------------------------------------------------------
subroutine fstr_rcap_finalize( fstrPARAM, fstrCPL )
      implicit none
      type( fstr_param  ) :: fstrPARAM
      type( fstr_couple ) :: fstrCPL

      write(IDBG,*) "fstr_rcap_finalize: start"

      if( fstrPARAM%fg_couple /= 1 ) return
      deallocate( fstrCPL%coupled_node )
      deallocate( fstrCPL%trac )
      deallocate( fstrCPL%disp )
      deallocate( fstrCPL%velo )
      deallocate( fstrCPL%accel )
      deallocate( fstrCPL%index )

      call rcapf_finalize()

      write(IDBG,*) "fstr_rcap_finalize: end"

end subroutine fstr_rcap_finalize
!------------------------------------------------------------------------------
subroutine fstr_rcap_send( fstrCPL )
      implicit none
      type( fstr_couple ) :: fstrCPL

      write(IDBG,*) "fstr_rcap_send: start"

      call rcapf_set_disp( fstrCPL%coupled_node,       &
                           fstrCPL%coupled_node_n,     &
                           fstrCPL%disp, fstrCPL%ndof )


      call rcapf_set_velo( fstrCPL%coupled_node,       &
                           fstrCPL%coupled_node_n,     &
                           fstrCPL%velo, fstrCPL%ndof )

      call rcapf_set_accel( fstrCPL%coupled_node,       &
                            fstrCPL%coupled_node_n,     &
                            fstrCPL%accel, fstrCPL%ndof )

      write(IDBG,*) "fstr_rcap_send: end"

end subroutine fstr_rcap_send
!------------------------------------------------------------------------------
subroutine fstr_rcap_get( fstrCPL )
      implicit none
      type( fstr_couple ) :: fstrCPL

      write(IDBG,*) "fstr_rcap_get: start"

      call rcapf_get_trac( fstrCPL%coupled_node,       &
                           fstrCPL%coupled_node_n,     &
                           fstrCPL%trac, fstrCPL%ndof )

      write(IDBG,*) "fstr_rcap_get: end"

end subroutine fstr_rcap_get
!------------------------------------------------------------------------------

end module m_fstr_rcap_io
