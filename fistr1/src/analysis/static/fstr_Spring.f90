!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief  This module provides functions to deal with spring force

module m_fstr_spring
  implicit none
contains

  subroutine fstr_AddSPRING(cstep, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrDYNAMIC)
    use m_fstr
    use m_static_lib
    use m_table_dyn
    integer, intent(in)                  :: cstep       !< current step
    type (hecmwST_matrix),intent(inout)  :: hecMAT      !< hecmw matrix
    type (hecmwST_local_mesh),intent(in) :: hecMESH     !< hecmw mesh
    type (fstr_solid),intent(inout)      :: fstrSOLID   !< fstr_solid
    type (fstr_param),intent(inout)      :: fstrPARAM   !< analysis control parameters
    type(fstr_dynamic), optional         :: fstrDYNAMIC

    integer(kind=kint) :: grpid, ndof, ig0, ig, ityp, iS0, iE0, ik, in, idx, num, flag_u
    real(kind=kreal) :: fval, factor

    factor = fstrSOLID%factor(2)

    ndof = hecMAT%NDOF
    do ig0= 1, fstrSOLID%SPRING_ngrp_tot
      grpid = fstrSOLID%SPRING_ngrp_GRPID(ig0)
      if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
      ig= fstrSOLID%SPRING_ngrp_ID(ig0)
      ityp= fstrSOLID%SPRING_ngrp_DOF(ig0)
      fval= fstrSOLID%SPRING_ngrp_val(ig0)
      if( present(fstrDYNAMIC) ) then
        flag_u = 11
        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, factor, flag_u)
        fval = fval*factor
        ! write(*,*)'(@v@)', factor
      else
        if( fval < 0.d0 ) fval = -fval*(1.d0-factor)
      end if

      iS0= hecMESH%node_group%grp_index(ig-1) + 1
      iE0= hecMESH%node_group%grp_index(ig  )
      do ik= iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        idx = ndof**2 * (in - 1) + ndof * (ityp - 1) + ityp
        hecMAT%D(idx) = hecMAT%D(idx) + fval
      enddo
    enddo

  end subroutine fstr_AddSPRING

  subroutine fstr_Update_NDForce_spring( cstep, hecMESH, fstrSOLID, B, fstrDYNAMIC )
    use m_fstr
    use m_table_dyn
    integer(kind=kint), intent(in)       :: cstep      !< current step
    type (hecmwST_local_mesh),intent(in) :: hecMESH    !< mesh information
    type (fstr_solid), intent(in)        :: fstrSOLID  !< we need boundary conditions of curr step
    real(kind=kreal), intent(inout)      :: B(:)       !< right hand side
    type(fstr_dynamic), optional         :: fstrDYNAMIC
    !    Local variables
    integer(kind=kint) ndof,ig0,ig,ityp,iS0,iE0,ik,in,idx,num
    integer(kind=kint) :: grpid, flag_u
    real(kind=kreal) :: fval, factor, damp

    factor = fstrSOLID%factor(2)

    ndof = hecMESH%n_dof
    do ig0= 1, fstrSOLID%SPRING_ngrp_tot
      grpid = fstrSOLID%SPRING_ngrp_GRPID(ig0)
      if( .not. fstr_isLoadActive( fstrSOLID, grpid, cstep ) ) cycle
      damp = fstrSOLID%SPRING_damper(ig0)
      ig= fstrSOLID%SPRING_ngrp_ID(ig0)
      ityp= fstrSOLID%SPRING_ngrp_DOF(ig0)
      fval= fstrSOLID%SPRING_ngrp_val(ig0)
      if( present(fstrDYNAMIC) ) then
        flag_u = 11
        call table_dyn(hecMESH, fstrSOLID, fstrDYNAMIC, ig0, factor, flag_u)
        fval = fval*factor
        ! write(*,*)'(@_@)', factor
      else
        if( fval < 0.d0 ) fval = -fval*(1.d0-factor)
      end if
      iS0= hecMESH%node_group%grp_index(ig-1) + 1
      iE0= hecMESH%node_group%grp_index(ig  )
      do ik= iS0, iE0
        in = hecMESH%node_group%grp_item(ik)
        idx = ndof * (in - 1) + ityp
        B(idx) = B(idx) - fval * ( fstrSOLID%dunode( idx ) + fstrSOLID%unode( idx ) * (1-damp) )
      enddo
    enddo
  end subroutine fstr_Update_NDForce_spring

end module m_fstr_spring
