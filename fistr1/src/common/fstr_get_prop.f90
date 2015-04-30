!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.6                                   !
!                                                                      !
!      Module Name : Static Analysis                                   !
!                                                                      !
!            Written by Toshio Nagashima (Sophia University)           !
!                       Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!

!> This module provide a function to fetch material properties from hecmw
module m_static_get_prop
contains
  !C
  !C***
  !C*** GET_PROP for FSTR solver
  !C***
  !C
  SUBROUTINE fstr_get_prop(hecMESH,fstrSOLID,cid,isect,ee,pp,rho,alpha,thick,alpha_over_mu,  &
       n_totlyr,shell_ortho,  &
       beam_radius,beam_angle1,beam_angle2,beam_angle3,   &
       beam_angle4,beam_angle5,beam_angle6)
    use m_fstr

    implicit none
    TYPE (hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)   :: fstrSOLID
    integer(kind=kint) :: cid
    integer(kind=kint) :: shell_ortho, n_ls, n_mod, n_totlyr, section_type

    n_totlyr = 1
    shell_ortho = -1

    !*EHM June 24 04
    !Get thickness
    n_item = hecMESH%section%sect_R_index(isect)-hecMESH%section%sect_R_index(isect-1)
    ihead = hecMESH%section%sect_R_index(isect-1)

    !do i = 1,n_item
    thick = hecMESH%section%sect_R_item(ihead+1)
    section_type = hecMESH%section%sect_type(isect)

    !IF(thick.LE.0.0) STOP "Zero thickness <= 0 is illegal"
    !Print *,'cval:',cval
    !PAUSE
    !end do
    !C** material ID
    mid=hecMESH%section%sect_mat_ID_item(isect)

    !C** Number of Item
    n_item=hecMESH%material%mat_ITEM_index(mid)-hecMESH%material%mat_ITEM_index(mid-1)

    !C** Head possition
    ihead=hecMESH%material%mat_ITEM_index(mid-1)
    !C Get ITEM of Meterial (Young's Modulus & Possion's Ratio
    pp=0.0
    alpha_over_mu = 0.0
    beam_radius = 0.0
    beam_angle1 = 0.0
    beam_angle2 = 0.0
    beam_angle3 = 0.0
    beam_angle4 = 0.0
    beam_angle5 = 0.0
    beam_angle6 = 0.0
    if( n_item .lt. 1 ) then
      write(IMSG,*) 'n_item=',n_item
      write(IMSG,*) '###Error 1'
      stop
    endif
    if ( n_item.ge.1 ) then
      !C** Number of Sub Item
      n_subitem = hecMESH%material%mat_subITEM_index(ihead+1)                          &
           - hecMESH%material%mat_subITEM_index(ihead)
      !C** Head possition
      mpos=hecMESH%material%mat_subITEM_index(ihead)
      !C Get SUBITEM of Meterial
      alpha_over_mu=1.0D-3
      n_subitem_ls = n_subitem
      !<<********************   shell analysis   ********************
      if (section_type == 2)then
        CALL fstr_get_prop_shell(hecMESH,fstrSOLID,cid,n_subitem,ee,pp,rho,alpha,thick,alpha_over_mu, &
             n_totlyr,shell_ortho,mpos)
      else
        if( n_subitem .lt. 1 ) then
          write(IMSG,*) '###Error 2'
          stop
        endif
        if ( n_subitem.ge.1 ) then
          ee=hecMESH%material%mat_val(mpos+1)
        endif
        if ( n_subitem.ge.2 ) then
          pp=hecMESH%material%mat_val(mpos+2)
        endif
        if ( n_subitem.ge.9 ) then
          beam_radius = hecMESH%material%mat_val(mpos+3)
          beam_angle1 = hecMESH%material%mat_val(mpos+4)
          beam_angle2 = hecMESH%material%mat_val(mpos+5)
          beam_angle3 = hecMESH%material%mat_val(mpos+6)
          beam_angle4 = hecMESH%material%mat_val(mpos+7)
          beam_angle5 = hecMESH%material%mat_val(mpos+8)
          beam_angle6 = hecMESH%material%mat_val(mpos+9)
        endif
      endif
    endif
    !C Get ITEM of Meterial (Density)
    rho=0.0
    if ( n_item.ge.2 ) then
      !C** Number of Sub Item
      n_subitem=hecMESH%material%mat_subITEM_index(ihead+2)                            &
           -hecMESH%material%mat_subITEM_index(ihead+1)
      !C** Head possition
      mpos=hecMESH%material%mat_subITEM_index(ihead+1)
      !C Get SUBITEM of Meterial
      if( n_subitem .lt. 1 ) then
        write(IMSG,*) '###Error 3'
        stop
      endif
      if ( n_subitem.ge.1 ) then
        rho=hecMESH%material%mat_val(mpos+1)
      endif
      if ( n_subitem.ge.2 ) then
        alpha_over_mu=hecMESH%material%mat_val(mpos+2)
      endif
    endif
    !C Get ITEM of Meterial (Thermal Expansion)
    alpha=0.0
    if ( n_item.ge.3 ) then
      !C** Number of Sub Item
      n_subitem=hecMESH%material%mat_subITEM_index(ihead+3)                            &
           -hecMESH%material%mat_subITEM_index(ihead+2)
      !C** Head possition
      mpos=hecMESH%material%mat_subITEM_index(ihead+2)
      !C Get SUBITEM of Meterial
      if( n_subitem .lt. 1 ) then
        write(IMSG,*) '###Error 4'
        stop
      endif
      if ( n_subitem.ge.1 ) then
        alpha=hecMESH%material%mat_val(mpos+1)
      endif
    endif

    fstrSOLID%materials(cid)%name = hecMESH%material%mat_name(cid)
    fstrSOLID%materials(cid)%variables(M_YOUNGS)=ee
    fstrSOLID%materials(cid)%variables(M_POISSON)=pp
    fstrSOLID%materials(cid)%variables(M_DENSITY)=rho
    fstrSOLID%materials(cid)%variables(M_EXAPNSION)=alpha
    fstrSOLID%materials(cid)%variables(M_THICK)=thick
    fstrSOLID%materials(cid)%variables(M_ALPHA_OVER_MU)= alpha_over_mu
    fstrSOLID%materials(cid)%variables(M_BEAM_RADIUS)=beam_radius
    fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE1)=beam_angle1
    fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE2)=beam_angle2
    fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE3)=beam_angle3
    fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE4)=beam_angle4
    fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE5)=beam_angle5
    fstrSOLID%materials(cid)%variables(M_BEAM_ANGLE6)=beam_angle6
    fstrSOLID%materials(cid)%mtype = ELASTIC

  end subroutine fstr_get_prop

  !-----Modifed by N.Morita GSFS, Univ. of Tokyo ------------------------------------------------------
  Subroutine fstr_get_prop_shell(hecMESH,fstrSOLID,cid,n_subitem,ee,pp,rho,alpha,thick,alpha_over_mu, &
       n_totlyr,shell_ortho,mpos)
    use m_fstr
    implicit none
    TYPE (hecmwST_local_mesh) :: hecMESH
    type(fstr_solid)   :: fstrSOLID
    real(kind=kreal) :: ee, pp, rho, alpha, thick, alpha_over_mu
    integer(kind=kint) :: cid, count, i, flag
    integer(kind=kint) :: shell_ortho, n_ls, n_mod, n_totlyr, section_type
    integer(kind=kint) :: n_subitem, mpos

    if( n_subitem .lt. 1 ) then
      write(ISTA,*) '###Error 2'
      stop
    elseif( n_subitem == 2) then
      n_totlyr = 1
      shell_ortho = 0
      ee = hecMESH%material%mat_val(mpos+1)
      pp = hecMESH%material%mat_val(mpos+2)

      allocate(fstrSOLID%materials(cid)%shell_var(1))
      fstrSOLID%materials(cid)%shell_var(1)%ortho = 0
      fstrSOLID%materials(cid)%shell_var(1)%ee = ee
      fstrSOLID%materials(cid)%shell_var(1)%pp = pp
      fstrSOLID%materials(cid)%shell_var(1)%thick = thick

    elseif( n_subitem == 3) then
      n_totlyr = 1
      shell_ortho = 0
      ee = hecMESH%material%mat_val(mpos+1)
      pp = hecMESH%material%mat_val(mpos+2)
      thick=hecMESH%material%mat_val(mpos+3)

      allocate(fstrSOLID%materials(cid)%shell_var(1))
      fstrSOLID%materials(cid)%shell_var(1)%ortho = 0
      fstrSOLID%materials(cid)%shell_var(1)%ee = ee
      fstrSOLID%materials(cid)%shell_var(1)%pp = pp
      fstrSOLID%materials(cid)%shell_var(1)%thick = thick

      write(ISTA,*) '###NOTICE : shell thickness is updated'

    elseif( n_subitem >= 4) then
      n_totlyr=0
      do i=1,n_subitem
        flag=hecMESH%material%mat_val(mpos+i)
        if(flag==0)then
          n_totlyr=n_totlyr+1
          i=i+3
        elseif(flag=1)then
          n_totlyr=n_totlyr+1
          i=i+8
        endif
      enddo
      allocate(fstrSOLID%materials(cid)%shell_var(n_totlyr))
      count=0
      do i=1,n_subitem
        !search material
        flag=hecMESH%material%mat_val(mpos+i)
        if(flag==0)then
          fstrSOLID%materials(cid)%shell_var(count)%ortho = hecMESH%material%mat_val(mpos+i  )
          fstrSOLID%materials(cid)%shell_var(count)%ee    = hecMESH%material%mat_val(mpos+i+1)
          fstrSOLID%materials(cid)%shell_var(count)%pp    = hecMESH%material%mat_val(mpos+i+2)
          fstrSOLID%materials(cid)%shell_var(count)%thick = hecMESH%material%mat_val(mpos+i+3)
          i=i+3
        elseif(flag=1)then
          fstrSOLID%materials(cid)%shell_var(count)%ortho = hecMESH%material%mat_val(mpos+i  )
          fstrSOLID%materials(cid)%shell_var(count)%ee    = hecMESH%material%mat_val(mpos+i+1)
          fstrSOLID%materials(cid)%shell_var(count)%pp    = hecMESH%material%mat_val(mpos+i+2)
          fstrSOLID%materials(cid)%shell_var(count)%ee2   = hecMESH%material%mat_val(mpos+i+3)
          fstrSOLID%materials(cid)%shell_var(count)%g12   = hecMESH%material%mat_val(mpos+i+4)
          fstrSOLID%materials(cid)%shell_var(count)%g23   = hecMESH%material%mat_val(mpos+i+5)
          fstrSOLID%materials(cid)%shell_var(count)%g31   = hecMESH%material%mat_val(mpos+i+6)
          fstrSOLID%materials(cid)%shell_var(count)%angle = hecMESH%material%mat_val(mpos+i+7)
          fstrSOLID%materials(cid)%shell_var(count)%thick = hecMESH%material%mat_val(mpos+i+8)
          i=i+8
        endif
        count=count+1
      enddo
    endif

    fstrSOLID%materials(cid)%totallyr=n_totlyr
    
  end subroutine fstr_get_prop_shell
end module m_static_get_prop
