!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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

      integer(kind=kint), private :: n_iter, n_rcap
      real(kind=kreal), private :: t_start, t_end, t_rcap
      real(kind=kreal), private :: t_start_all, t_end_all, t_rcap_all

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
      real(kind=kreal) :: t_s, t_e

      write(IDBG,*) "fstr_rcap_initialize: start"
      t_start_all = hecmw_Wtime()

      if( fstrPARAM%fg_couple /= 1 ) return

      write(IDBG,*) "fstr_rcap_initialize: calling rcapf_init_solid_solver"

      portfile='port'//CHAR(0)

      t_s = hecmw_Wtime()
      call rcapf_init_solid_solver( hecMESH%my_rank, portfile )
      t_e = hecmw_Wtime()
      t_rcap_all = t_e - t_s

      write(IDBG,*) "fstr_rcap_initialize: calling rcapf_get_num_of_matching_node"

      fstrCPL%dof = 3

      t_s = hecmw_Wtime()
      call rcapf_get_num_of_matching_node( fstrCPL%coupled_node_n )
      t_e = hecmw_Wtime()
      t_rcap_all = t_rcap_all + (t_e - t_s)

      fstrCPL%ndof =       fstrCPL%coupled_node_n * fstrCPL%dof
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
      fstrCPL%trac = 0.0d0
      fstrCPL%index = -1

      write(IDBG,*) "fstr_rcap_initialize: calling rcapf_get_mathcing_node_id"

      t_s = hecmw_Wtime()
        call rcapf_get_matching_node_id( fstrCPL%coupled_node, fstrCPL%coupled_node_n )
      t_e = hecmw_Wtime()
      t_rcap_all = t_rcap_all + (t_e - t_s)

      write(IDBG,*) "fstr_rcap_initialize: converting to local id: ", fstrCPL%coupled_node_n

      do i=1, fstrCPL%coupled_node_n
            nid = fstrCPL%coupled_node(i)
            if( nid <= 0 .or. nid>hecMESH%n_node_gross ) then
                  write(*,*) '##Fatal error in fstr_rcap_initialize ', i, nid
                  call hecmw_abort( hecmw_comm_get_comm());
            endif
            if( hecMESH%n_refine > 0 ) then
                  nid_old = nid
                  if( associated( hecMESH%node_old2new ) ) then
                    nid = hecMESH%node_old2new( nid ) + 1
                  endif
                  write(IDBG,*) nid_old, nid
            endif
            fstrCPL%index( nid ) = i
      end do

      n_iter = 0
      n_rcap = 0

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
      real(kind=kreal) :: t_tot, t_tot_avg, t_rcap_avg, tr
      real(kind=kreal) :: t_tot_all, tr_all
      real(kind=kreal) :: t_s, t_e

      write(IDBG,*) "fstr_rcap_finalize: start"

      if( fstrPARAM%fg_couple /= 1 ) return
      deallocate( fstrCPL%coupled_node )
      deallocate( fstrCPL%trac )
      deallocate( fstrCPL%disp )
      deallocate( fstrCPL%velo )
      deallocate( fstrCPL%accel )
      deallocate( fstrCPL%index )

      write(IDBG,*) "fstr_rcap_finalize: calling rcapf_finalize"

      t_s = hecmw_Wtime()
      call rcapf_finalize()
      t_e = hecmw_Wtime()
      t_rcap_all = t_rcap_all + (t_e - t_s)

      t_tot = t_end - t_start
      t_tot_avg = t_tot / n_iter
      t_rcap_avg = t_rcap / n_rcap
      tr = t_rcap_avg / t_tot_avg * 100.d0

      write(IDBG,'(a,f11.3,a,i0,a,f11.3,a)') &
	& 'fstr + rcap:', t_tot,'s for ',n_iter,' iters / avg. ', t_tot_avg,'s/iter'
      write(IDBG,'(a,f11.3,a,i0,a,f11.3,a)') &
	& '       rcap:',t_rcap,'s for ',n_rcap,' calls / avg. ',t_rcap_avg,'s/call'
      write(IDBG,'(a,f11.3,a)') &
	& ' rcap ratio:',tr,'%/iter'

      t_end_all = hecmw_Wtime()

      t_tot_all = t_end_all - t_start_all
      tr_all = t_rcap_all / t_tot_all * 100.d0

      write(IDBG,'(a,f11.3,a)') 'fstr total:',t_tot_all,'s'
      write(IDBG,'(a,f11.3,a)') 'rcap total:',t_rcap_all,'s'
      write(IDBG,'(a,f11.3,a)') 'rcap ratio:',tr_all,'%'

      write(IDBG,*) "fstr_rcap_finalize: end"

end subroutine fstr_rcap_finalize
!------------------------------------------------------------------------------
subroutine fstr_rcap_send( fstrCPL )
      implicit none
      type( fstr_couple ) :: fstrCPL

      write(IDBG,*) "fstr_rcap_send: start"

!      call rcapf_set_disp( fstrCPL%coupled_node,       &
!                           fstrCPL%coupled_node_n,     &
!                           fstrCPL%disp, fstrCPL%ndof )

      call rcapf_set_velo( fstrCPL%coupled_node,       &
                           fstrCPL%coupled_node_n,     &
                           fstrCPL%velo, fstrCPL%ndof )

!      call rcapf_set_accel( fstrCPL%coupled_node,       &
!                            fstrCPL%coupled_node_n,     &
!                            fstrCPL%accel, fstrCPL%ndof )

      write(IDBG,*) "fstr_rcap_send: end"

end subroutine fstr_rcap_send
!------------------------------------------------------------------------------
subroutine fstr_rcap_get( fstrCPL )
      implicit none
      type( fstr_couple ) :: fstrCPL
      real(kind=kreal) :: t_s, t_e

      write(IDBG,*) "fstr_rcap_get: start"

      if (n_rcap == 0) then
         t_start = hecmw_Wtime()
         t_rcap = 0
      else
         t_end = hecmw_Wtime()
         n_iter = n_iter + 1
      end if

      write(IDBG,*) "fstr_rcap_get: calling rcapf_get_trac"

      t_s = hecmw_Wtime()
      call rcapf_get_trac( fstrCPL%coupled_node,       &
                           fstrCPL%coupled_node_n,     &
                           fstrCPL%trac, fstrCPL%ndof )
      t_e = hecmw_Wtime()

      t_rcap = t_rcap + (t_e - t_s)
      n_rcap = n_rcap + 1

      t_rcap_all = t_rcap_all + (t_e - t_s)

      write(IDBG,*) "fstr_rcap_get: end"

end subroutine fstr_rcap_get
!------------------------------------------------------------------------------
subroutine fstr_get_convergence( revocap_flag )
      implicit none
      integer(kind=kint)  :: revocap_flag

      write(IDBG,*) "fstr_get_convergence: start"

      call rcapf_get_convergence( revocap_flag )

      write(IDBG,*) "fstr_get_convergence: end"

end subroutine fstr_get_convergence
!------------------------------------------------------------------------------

end module m_fstr_rcap_io
