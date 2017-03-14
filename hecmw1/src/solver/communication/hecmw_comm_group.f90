!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.8                                                !
!                                                                      !
!     Last Update : 2016/09/16                                         !
!        Category : I/O and Utility                                    !
!                                                                      !
!            Written by Kazuhisa Inagaki (FATEC)                       !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!

module hecmw_comm_group
    use hecmw_util
    use m_hecmw_comm_f
    implicit none
    
    contains

  function hecmw_ngrp_get_number(hecMESH, ngrpid)
    type(hecmwST_local_mesh) :: hecMESH
    integer(kind=kint) :: ngrpid
    integer(kind=kint) :: hecmw_ngrp_get_number

    integer(kind=kint) :: iS, iE, ik, in 
    integer(kind=kint) :: n_ngrp_internal, nn_internal
    
    nn_internal = hecMESH%nn_internal

    iS = hecMESH%node_group%grp_index(ngrpid-1)+1
    iE = hecMESH%node_group%grp_index(ngrpid)
    
    n_ngrp_internal = 0
    do ik=iS,iE
      in = hecMESH%node_group%grp_item(ik)
      if( in > nn_internal ) cycle
      n_ngrp_internal = n_ngrp_internal + 1
    enddo
    
    call hecmw_allreduce_I1 (hecMESH, n_ngrp_internal, hecmw_sum)
    hecmw_ngrp_get_number = n_ngrp_internal

  end function hecmw_ngrp_get_number    

  function hecmw_ngrp_get_totalvalue(hecMESH, ngrpid, ndof, idof, vector)
    type(hecmwST_local_mesh) :: hecMESH    ! local mesh
    integer(kind=kint) :: ngrpid           ! node group id
    integer(kind=kint) :: idof             ! dof to be returned
    integer(kind=kint) :: ndof             ! dof of noval value
    real(kind=kreal), pointer :: vector(:) ! nodal value vector (size=nodf*hecMESH&n_node)  
    real(kind=kreal) :: hecmw_ngrp_get_totalvalue

    integer(kind=kint) :: iS, iE, ik, in 
    integer(kind=kint) :: nn_internal
    real(kind=kreal)   :: value

    nn_internal = hecMESH%nn_internal
    
    iS = hecMESH%node_group%grp_index(ngrpid-1)+1
    iE = hecMESH%node_group%grp_index(ngrpid)
    
    value = 0.d0
    do ik=iS,iE
      in = hecMESH%node_group%grp_item(ik)
      if( in > nn_internal ) cycle
      value = value + vector(ndof*(in-1)+idof)
    enddo
    
    call hecmw_allreduce_R1 (hecMESH, value, hecmw_sum)
    hecmw_ngrp_get_totalvalue = value

  end function hecmw_ngrp_get_totalvalue    
   
end module hecmw_comm_group
