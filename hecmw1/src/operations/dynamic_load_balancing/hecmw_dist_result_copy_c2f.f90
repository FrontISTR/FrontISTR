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

subroutine hecmw_dist_result_copy_c2f(hecMESHnew, adapRES)
      use   hecmw_util
      use   hecmw_result 
      use   hecmw_io
      use   hecmw_dist_copy_f2c_f
      type (hecmwST_local_mesh) :: hecMESHnew
      type (hecmwST_result_data):: adapRES
   
	  

  if(adapRES%nn_component .gt.0) then
      do i=1, adapRES%nn_component
         in=in+adapRES%nn_dof(i)
      enddo
      allocate (adapRES%node_val_item(in*hecMESHnew%n_node))
      call hecmw_dlb_get_result_node(adapRES%node_val_item)
  endif
  if(adapRES%ne_component .gt.0) then
      do i=1, adapRES%ne_component
         in=in+adapRES%ne_dof(i)
      enddo
      allocate (adapRES%elem_val_item(in*hecMESHnew%n_elem))
      call hecmw_dlb_get_result_elem(adapRES%elem_val_item)
  endif

end subroutine  hecmw_dist_result_copy_c2f
  
