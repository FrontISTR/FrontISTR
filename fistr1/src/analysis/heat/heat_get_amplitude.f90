!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Heat Analysis                                     !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> This moudle provide a function to get amplitude definition
module m_heat_get_amplitude
   contains
!C***
!C*** GET_AMPLITUDE
!C***
   subroutine heat_get_amplitude ( fstrHEAT,ID,TT,QQ,OutOfRange )

      use m_fstr

      implicit none
      integer(kind=kint) ID,nn,ii,ikk
      real(kind=kreal)    TT,QQ
      type (fstr_heat   ) :: fstrHEAT
      logical, optional :: OutOfRange

      QQ = 1.0
      if (present(OutOfRange)) OutOfRange = .false.
      if( ID>0 ) then
        nn = fstrHEAT%AMPLtab(ID)
        if    ( TT < fstrHEAT%AMPLtime(ID,1) ) then
          ii = 1
          if (present(OutOfRange)) OutOfRange = .true.
        elseif( TT >= fstrHEAT%AMPLtime(ID,nn) ) then
          ii = nn + 1
          if (present(OutOfRange)) OutOfRange = .true.
        else
          ii = 2
          if (present(OutOfRange)) OutOfRange = .false.
          do ikk= 1, nn - 1
            if(       TT .ge. fstrHEAT%AMPLtime(ID,ikk)                        &
               .and.  TT .lt. fstrHEAT%AMPLtime(ID,ikk+1) ) then
              ii = ikk + 1
              exit
            endif
          enddo
        endif
        QQ = fstrHEAT%AMPLfuncA(ID,ii) * TT + fstrHEAT%AMPLfuncB(ID,ii)
      endif

   end subroutine heat_get_amplitude
end module m_heat_get_amplitude
