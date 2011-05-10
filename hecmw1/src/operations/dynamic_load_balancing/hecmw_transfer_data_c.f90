!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Dynamic Load Balancing                             !
!                                                                      !
!            Written by Li Chen (Univ. of Tokyo)                       !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

      subroutine hecmw_transfer_data_f2c (hecMESHnew, adapRES)

      use  hecmw_util
      use  hecmw_io
      use  hecmw_dist_copy_f2c_f
      use  hecmw_result
      type (hecmwST_local_mesh) :: hecMESHnew
      type (hecmwST_result_data):: adapRES

      call hecmw_dlb_f2c_init
      call hecmw_dist_copy_f2c(hecMESHnew, ierr)
      call hecmw_dlb_f2c_finalize
      call test_mesh
      call hecmw_dist_result_copy_f2c(adapRES)
      end subroutine hecmw_transfer_data_f2c

