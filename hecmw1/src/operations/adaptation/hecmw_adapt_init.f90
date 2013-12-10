!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Adaptive Mesh Refinement                           !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

      subroutine hecmw_adapt_init (hecMESH)

      use  hecmw_util
      type      (hecmwST_local_mesh) :: hecMESH

      if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#allocate (n_node, nn_int, n_elem, ne_int)'
      call hecmw_adapt_allocate   (hecMESH)

      if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#ACTIVE'
      call hecmw_adapt_ACTIVE (hecMESH)

      if (hecMESH%my_rank.eq.0) write (*,'(/,a)') '#E-COMM.TAB.'
      call hecmw_adapt_EDGE_COMM_TABLE (hecMESH)
      if (hecMESH%my_rank.eq.0) write (*,'(  a)') '#C-COMM.TAB.'
      call hecmw_adapt_CELL_COMM_TABLE (hecMESH)

      end subroutine hecmw_adapt_init
