!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.1                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Kazuaki Sakane (RIST)                          !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_io
    use hecmw_util
    use hecmw_dist_copy_f2c_f
    use hecmw_dist_copy_c2f_f
    use hecmw_dist_free_f
    use hecmw_dist_print_f
    use hecmw_result
    use hecmw_restart
    implicit none

    public :: hecmw_get_mesh
    public :: hecmw_put_mesh
    
    contains

!C====================================================================
!C Get HEC-MW dist mesh from file
!C====================================================================

    subroutine hecmw_get_mesh(name_ID, mesh) 
        integer(kind=kint) :: ierr  
        character(len=HECMW_NAME_LEN) :: name_ID
        type(hecmwST_local_mesh) :: mesh

        call hecmw_nullify_mesh(mesh)

        call hecmw_get_mesh_init_if(name_ID,ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_dist_copy_c2f(mesh, ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_get_mesh_finalize_if(ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

    end subroutine hecmw_get_mesh


!C====================================================================
!C Put HEC-MW dist mesh to file
!C====================================================================

    subroutine hecmw_put_mesh(name, mesh) 
        integer(kind=kint) :: ierr  
        character(len=HECMW_NAME_LEN) :: name
        type(hecmwST_local_mesh) :: mesh

        call hecmw_put_mesh_init_if(ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_dist_copy_f2c(mesh, ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_put_mesh_if(name, ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())

        call hecmw_put_mesh_finalize_if(ierr) 
        if(ierr /=0) call hecmw_abort(hecmw_comm_get_comm())
    end subroutine hecmw_put_mesh

end module hecmw_io

