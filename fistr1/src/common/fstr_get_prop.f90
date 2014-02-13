!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
  SUBROUTINE fstr_get_prop(hecMESH,isect,ee,pp,rho,alpha,thick,alpha_over_mu,  &
       n_totallayer_ls,shell_matltype,shell_variables,  &
       beam_radius,beam_angle1,beam_angle2,beam_angle3,   &
       beam_angle4,beam_angle5,beam_angle6)
    use m_fstr

    IMPLICIT REAL(kind=kreal) (A-H,O-Z)
    real(kind=kreal) :: shell_variables(200)
    integer(kind=kint) :: shell_matltype, n_ls, n_mod, n_totallayer_ls, section_type
    TYPE (hecmwST_local_mesh) :: hecMESH

    n_totallayer_ls = 1
    shell_matltype = -1

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
        CALL fstr_get_prop_shell(hecMESH,n_subitem,ee,pp,rho,alpha,thick,alpha_over_mu, &
             n_totallayer_ls,shell_matltype,shell_variables,mpos)
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
        alpha_over_mu=1.0D-3
        if ( n_subitem.ge.3 ) then
          alpha_over_mu=hecMESH%material%mat_val(mpos+3)
        endif
        if ( n_subitem.ge.9 ) then
          alpha_over_mu = 0.0D0
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
      alpha_over_mu=1.0D-3
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

  end subroutine FSTR_GET_PROP

  !-----Modifed by N.Morita GSFS, Univ. of Tokyo ------------------------------------------------------
  Subroutine FSTR_GET_PROP_SHELL(hecMESH,n_subitem,ee,pp,rho,alpha,thick,alpha_over_mu, &
       n_totallayer_ls,shell_matltype,shell_variables,mpos)
    use m_fstr
    implicit none
    real(kind=kreal) :: shell_variables(200), ee, pp, rho, alpha, thick, alpha_over_mu
    integer(kind=kint) :: shell_matltype, n_ls, n_mod, n_totallayer_ls, section_type
    integer(kind=kint) :: n_subitem, mpos
    TYPE (hecmwST_local_mesh) :: hecMESH

    if( n_subitem .lt. 1 ) then
      write(IMSG,*) '###Error 2'
      stop
    elseif( n_subitem == 2) then
      n_totallayer_ls = 1
      shell_matltype = 0
      ee = hecMESH%material%mat_val(mpos+1)
      pp = hecMESH%material%mat_val(mpos+2)
      shell_variables(1) = hecMESH%material%mat_val(mpos+1)
      shell_variables(2) = hecMESH%material%mat_val(mpos+2)
      shell_variables(3) = thick
      write(IMSG,*) '###Warning : Old format of shell properties'
    elseif( n_subitem >= 4) then
      shell_matltype = int(hecMESH%material%mat_val(mpos+1))
      if (shell_matltype == 0) then
        n_mod = mod(n_subitem-1,3)
        if (n_mod.eq.0) then
          n_totallayer_ls = (n_subitem-1)/3
          thick = 0.0D0
          DO n_ls=1,n_totallayer_ls
            shell_variables(3*n_ls-2) = hecMESH%material%mat_val(mpos+3*(n_ls-1)+2)
            shell_variables(3*n_ls-1) = hecMESH%material%mat_val(mpos+3*(n_ls-1)+3)
            shell_variables(3*n_ls  ) = hecMESH%material%mat_val(mpos+3*(n_ls-1)+4)
            thick = thick + hecMESH%material%mat_val(mpos+3*(n_ls-1)+4)
          END DO
          ee = shell_variables(1)
          pp = shell_variables(2)
        else
          write(IMSG,*) '###Error : shell properties not correct (isotropic)'
          stop
        endif
      elseif(shell_matltype == 1)then
        n_mod = mod(n_subitem-1,8)
        if (n_mod.eq.0) then
          n_totallayer_ls = (n_subitem-1)/8
          thick = 0.0D0
          DO n_ls=1,n_totallayer_ls
            shell_variables(8*n_ls-7) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+2)
            shell_variables(8*n_ls-6) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+3)
            shell_variables(8*n_ls-5) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+4)
            shell_variables(8*n_ls-4) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+5)
            shell_variables(8*n_ls-3) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+6)
            shell_variables(8*n_ls-2) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+7)
            shell_variables(8*n_ls-1) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+8)
            shell_variables(8*n_ls  ) = hecMESH%material%mat_val(mpos+8*(n_ls-1)+9)
            thick = thick + shell_variables(8*n_ls-5)
          END DO
          ee = shell_variables(1)
          pp = shell_variables(2)
        else
          write(IMSG,*) '###Error : shell properties not correct(orthotropic)'
          stop
        endif
      else
        write(IMSG,*) '###Error : shell_matltype not correct'
        stop
      endif
    endif
  end subroutine FSTR_GET_PROP_SHELL
end module m_static_get_prop
