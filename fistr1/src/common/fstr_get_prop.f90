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
  SUBROUTINE fstr_get_prop(hecMESH,shell_var,isect,ee,pp,rho,alpha,thick,n_totlyr,alpha_over_mu,  &
     beam_radius,beam_angle1,beam_angle2,beam_angle3,   &
     beam_angle4,beam_angle5,beam_angle6)
  use m_fstr
  use mMaterial

  implicit none
  type (hecmwST_local_mesh) :: hecMESH
  type (tshellmat), pointer :: shell_var(:)
  integer(kind=kint) :: n_item, n_subitem, ihead, isect, mid, mpos
  integer(kind=kint) :: shell_ortho, n_ls, n_mod, n_totlyr, section_type
  real(kind=kreal) :: ee, pp, rho, alpha, thick, alpha_over_mu
  real(kind=kreal) :: beam_radius,beam_angle1,beam_angle2,beam_angle3,beam_angle4,beam_angle5,beam_angle6

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
    !<<********************   shell analysis   ********************
    if (section_type == 2)then
    CALL fstr_get_prop_shell(hecMESH,shell_var,mid,n_subitem,ee,pp,rho,alpha,thick,alpha_over_mu, &
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

  end subroutine fstr_get_prop

  !-----Modifed by N.Morita GSFS, Univ. of Tokyo ------------------------------------------------------
  Subroutine fstr_get_prop_shell(hecMESH,shell_var,mid,n_subitem,ee,pp,rho,alpha,thick,alpha_over_mu, &
     n_totlyr,shell_ortho,mpos)
  use m_fstr
  implicit none
  type (hecmwST_local_mesh) :: hecMESH
  type (tshellmat),pointer :: shell_var(:)
  real(kind=kreal) :: ee, pp, rho, alpha, thick, alpha_over_mu, tmp
  integer(kind=kint) :: mid, count, i, flag
  integer(kind=kint) :: shell_ortho, n_ls, n_mod, n_totlyr, section_type
  integer(kind=kint) :: n_subitem, mpos

  if( n_subitem .lt. 1 ) then
    write(IMSG,*) '###Error 2'
    stop
  elseif( n_subitem == 2) then
    n_totlyr = 1
    shell_ortho = 0
    ee = hecMESH%material%mat_val(mpos+1)
    pp = hecMESH%material%mat_val(mpos+2)

    allocate(shell_var(1))
    shell_var(1)%ortho = 0
    shell_var(1)%ee = ee
    shell_var(1)%pp = pp
    shell_var(1)%thick = thick
    shell_var(1)%weight= 1.0d0

  elseif( n_subitem == 3) then
    n_totlyr = 1
    shell_ortho = 0
    ee = hecMESH%material%mat_val(mpos+1)
    pp = hecMESH%material%mat_val(mpos+2)
    thick=hecMESH%material%mat_val(mpos+3)

    allocate(shell_var(1))
    shell_var(1)%ortho = 0
    shell_var(1)%ee = ee
    shell_var(1)%pp = pp
    shell_var(1)%thick = thick
    shell_var(1)%weight= 1.0d0

    write(IMSG,*) '###NOTICE : shell thickness is updated'

  elseif( n_subitem >= 4) then
    n_totlyr=0
      
    i=1
    do
      if(i >= n_subitem) exit
      flag=hecMESH%material%mat_val(mpos+i)
      if(flag == 0)then
        n_totlyr=n_totlyr+1
        i=i+4
      elseif(flag == 1)then
        n_totlyr=n_totlyr+1
        i=i+9
      else
        write(IMSG,*) '###Error: Shell property invalid'
        STOP
      endif
    enddo
    allocate(shell_var(n_totlyr))
    count=1
    i=1
    do
    if(i >= n_subitem)exit
    !search material
    flag=hecMESH%material%mat_val(mpos+i)
    if(flag==0)then
      shell_var(count)%ortho = hecMESH%material%mat_val(mpos+i  )
      shell_var(count)%ee    = hecMESH%material%mat_val(mpos+i+1)
      shell_var(count)%pp    = hecMESH%material%mat_val(mpos+i+2)
      shell_var(count)%thick = hecMESH%material%mat_val(mpos+i+3)
      shell_var(count)%weight= hecMESH%material%mat_val(mpos+i+3)/thick
      i=i+4
    elseif(flag==1)then
      shell_var(count)%ortho = hecMESH%material%mat_val(mpos+i  )
      shell_var(count)%ee    = hecMESH%material%mat_val(mpos+i+1)
      shell_var(count)%pp    = hecMESH%material%mat_val(mpos+i+2)
      shell_var(count)%ee2   = hecMESH%material%mat_val(mpos+i+3)
      shell_var(count)%g12   = hecMESH%material%mat_val(mpos+i+4)
      shell_var(count)%g23   = hecMESH%material%mat_val(mpos+i+5)
      shell_var(count)%g31   = hecMESH%material%mat_val(mpos+i+6)
      shell_var(count)%angle = hecMESH%material%mat_val(mpos+i+7)
      shell_var(count)%thick = hecMESH%material%mat_val(mpos+i+8)
      shell_var(count)%weight= hecMESH%material%mat_val(mpos+i+8)/thick
      i=i+9
    else
      write(IMSG,*) '###Error: Shell property invalid'
      STOP
    endif
    count=count+1
    enddo

    !check weight
    tmp = 0.0d0
    do i=1,n_totlyr
      tmp = tmp + shell_var(i)%weight
    enddo
    if(tmp == 1.0d0)then
      write(IMSG,"(a)")"### NOTICCE: Total thickness is not equal to the sum of laminated layers' thickness"
    endif
  endif
  
  end subroutine fstr_get_prop_shell
end module m_static_get_prop
