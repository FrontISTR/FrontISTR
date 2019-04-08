!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> This moudle provide a function to get amplitude definition
module m_heat_get_amplitude
contains
  !C***
  !C*** GET_AMPLITUDE
  !C***
  subroutine heat_get_amplitude ( fstrHEAT,id,TT,QQ,OutOfRange )

    use m_fstr

    implicit none
    integer(kind=kint) :: id, nn, ii, ikk
    real(kind=kreal)   :: TT, QQ
    type(fstr_heat)    :: fstrHEAT
    logical, optional  :: OutOfRange

    QQ = 1.0
    if (present(OutOfRange)) OutOfRange = .false.
    if( id>0 ) then
      nn = fstrHEAT%AMPLtab(id)
      if    ( TT < fstrHEAT%AMPLtime(id,1) ) then
        ii = 1
        if (present(OutOfRange)) OutOfRange = .true.
      elseif( TT >= fstrHEAT%AMPLtime(id,nn) ) then
        ii = nn + 1
        if (present(OutOfRange)) OutOfRange = .true.
      else
        ii = 2
        if (present(OutOfRange)) OutOfRange = .false.
        do ikk= 1, nn - 1
          if(       TT .ge. fstrHEAT%AMPLtime(id,ikk)                        &
              .and.  TT .lt. fstrHEAT%AMPLtime(id,ikk+1) ) then
            ii = ikk + 1
            exit
          endif
        enddo
      endif
      QQ = fstrHEAT%AMPLfuncA(id,ii) * TT + fstrHEAT%AMPLfuncB(id,ii)
    endif

  end subroutine heat_get_amplitude
end module m_heat_get_amplitude
