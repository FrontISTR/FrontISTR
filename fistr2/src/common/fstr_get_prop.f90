!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Static Analysis                                    !
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
   SUBROUTINE fstr_get_prop(hecMESH,isect,ee,pp,rho,alpha,thick)

      use m_fstr
      
      IMPLICIT REAL(kind=kreal) (A-H,O-Z)
      TYPE (hecmwST_local_mesh) :: hecMESH

!*EHM June 24 04
!Get thickness
      n_item = hecMESH%section%sect_R_index(isect)-hecMESH%section%sect_R_index(isect-1)
      !Print *,'n_item for thickness: ',n_item
      !PAUSE
      ihead = hecMESH%section%sect_R_index(isect-1)
      !Print *,'ihead for thickness: ',ihead
      !PAUSE
      !do i = 1,n_item
       thick = hecMESH%section%sect_R_item(ihead+1)
       IF(thick.LE.0.0) STOP "Zero thickness <= 0 is illegal"
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
end module m_static_get_prop
