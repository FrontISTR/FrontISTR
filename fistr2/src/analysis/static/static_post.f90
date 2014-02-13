!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
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
!C***
!>  Prepare and output max and min value of results
!C***
module m_static_post
contains

      subroutine FSTR_POST(fstrSOLID,istep)
      use m_fstr
          use elementInfo
      implicit none
      type (fstr_solid)         :: fstrSOLID
      integer                   :: istep
	  
      include "HEC_MW3_For.h"

      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN) :: label
      character(len=HECMW_NAME_LEN) :: nameID
!C Max Min
      real(kind=kreal)    Umax(6), Umin(6)
      integer(kind=kint) IUmax(6),IUmin(6)

      real(kind=kreal)    EEmax(6), EEmin(6) ,ESmax(6), ESmin(6)
      integer(kind=kint) IEEmax(6),IEEmin(6),IESmax(62),IESmin(6)
	  
      integer :: i, j, k, id, ndof, MDOF, nd


!C Clear
      Umax=0.0d0
      Umin=0.0d0
      EEmax=0.0d0
      EEmin=0.0d0
      ESmax=0.0d0
      ESmin=0.0d0
      IUmax=0
      IUmin=0
      IEEmax=0
      IEEmin=0
      IESmax=0
      IESmin=0
!C*** Update Strain and Stress
      ndof= assDOF(1)
!C*** Show Displacement 
      write(IDBG,*) '#### DISPLACEMENT@POST'
      do i= 1, total_node
        if( i==1) then
          do k=1,NDOF 
            Umax(k)=fstrSOLID%unode(NDOF*(i-1)+k)
            Umin(k)=fstrSOLID%unode(NDOF*(i-1)+k)
            IUmax(k)=i
            IUmin(k)=i
          enddo
        else
          do k=1,NDOF 
            if( fstrSOLID%unode(NDOF*(i-1)+k) .gt. Umax(k) ) then
              Umax(k)=fstrSOLID%unode(NDOF*(i-1)+k)
              IUmax(k)=i
            endif
            if( fstrSOLID%unode(NDOF*(i-1)+k) .lt. Umin(k) ) then
              Umin(k)=fstrSOLID%unode(NDOF*(i-1)+k)
              IUmin(k)=i
            endif
          enddo
        endif
      enddo
!C*** Show Strain
      if( NDOF == 2 ) NDOF = 3
      if( NDOF == 3 ) MDOF = 6
      if( NDOF == 6 ) MDOF = 5
      write(IDBG,*) '#### STRAIN@POST'

!C @element
      do i= 1, total_elem
        if( i==1) then
          do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
            do k=1,MDOF 
              EEmax(k)=fstrSOLID%elements(i)%gausses(j)%strain(k)
              EEmin(k)=fstrSOLID%elements(i)%gausses(j)%strain(k)
              IEEmax(k)=i
              IEEmin(k)=i
            enddo
          enddo
        else
          do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
            do k=1,MDOF
              if( fstrSOLID%elements(i)%gausses(j)%strain(k) > EEmax(k) ) then
                EEmax(k)=fstrSOLID%elements(i)%gausses(j)%strain(k)
                IEEmax(k)=i
              endif
              if( fstrSOLID%elements(i)%gausses(j)%strain(k) < EEmin(k) ) then
                EEmin(k)=fstrSOLID%elements(i)%gausses(j)%strain(k)
                IEEmin(k)=i
              endif
            enddo
          enddo
        endif
      enddo
!C*** Show Stress
      if( NDOF .eq. 3 ) MDOF = 6
      if( NDOF .eq. 6 ) MDOF = 6
      write(IDBG,*) '#### STRESS@POST'
!C @element
      do i= 1, total_elem
        if( i==1) then
          do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
            do k=1,MDOF
              ESmax(k)=fstrSOLID%elements(i)%gausses(j)%stress(k)
              ESmin(k)=fstrSOLID%elements(i)%gausses(j)%stress(k)
              IESmax(k)=i
              IESmin(k)=i
            enddo
          enddo
        else
          do j=1,NumOfQuadPoints( fstrSOLID%elements(i)%etype )
            do k=1,MDOF
              if( fstrSOLID%elements(i)%gausses(j)%stress(k) > ESmax(k) ) then
                ESmax(k)=fstrSOLID%elements(i)%gausses(j)%stress(k)
                IESmax(k)=i
              endif
              if( fstrSOLID%elements(i)%gausses(j)%stress(k) < ESmin(k) ) then
                ESmin(k)=fstrSOLID%elements(i)%gausses(j)%stress(k)
                IESmin(k)=i
              endif
            enddo
          enddo
        endif
      enddo
!C
!
      if( NDOF==2 .or. NDOF==3 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//U1 ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009) '//U2 ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009) '//U3 ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009) '//E11',EEmax(1),IEEmax(1),EEmin(1),IEEmin(1)
        write(ILOG,1009) '//E22',EEmax(2),IEEmax(2),EEmin(2),IEEmin(2)
        write(ILOG,1009) '//E33',EEmax(3),IEEmax(3),EEmin(3),IEEmin(3)
        write(ILOG,1009) '//E12',EEmax(4),IEEmax(4),EEmin(4),IEEmin(4)
        write(ILOG,1009) '//E23',EEmax(5),IEEmax(5),EEmin(5),IEEmin(5)
        write(ILOG,1009) '//E13',EEmax(6),IEEmax(6),EEmin(6),IEEmin(6)
        write(ILOG,1009) '//S11',ESmax(1),IESmax(1),ESmin(1),IESmin(1)
        write(ILOG,1009) '//S22',ESmax(2),IESmax(2),ESmin(2),IESmin(2)
        write(ILOG,1009) '//S33',ESmax(3),IESmax(3),ESmin(3),IESmin(3)
        write(ILOG,1009) '//S12',ESmax(4),IESmax(4),ESmin(4),IESmin(4)
        write(ILOG,1009) '//S23',ESmax(5),IESmax(5),ESmin(5),IESmin(5)
        write(ILOG,1009) '//S13',ESmax(6),IESmax(6),ESmin(6),IESmin(6)


!C*** Show Summary
        call mw_allreduce_r( Umax, 3, mw_mpi_max() )
        call mw_allreduce_r( Umin, 3, mw_mpi_min() )
        call mw_allreduce_r( EEmax, 6, mw_mpi_max() )
        call mw_allreduce_r( EEmin, 6, mw_mpi_min() )
        call mw_allreduce_r( ESmax, 6, mw_mpi_max() )
        call mw_allreduce_r( ESmin, 6, mw_mpi_min() )

        if( myrank == 0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1 ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2 ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3 ',Umax(3),Umin(3)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11',EEmax(1),EEmin(1)
          write(ILOG,1019) '//E22',EEmax(2),EEmin(2)
          write(ILOG,1019) '//E33',EEmax(3),EEmin(3)
          write(ILOG,1019) '//E12',EEmax(4),EEmin(4)
          write(ILOG,1019) '//E23',EEmax(5),EEmin(5)
          write(ILOG,1019) '//E13',EEmax(6),EEmin(6)
          write(ILOG,1019) '//S11',ESmax(1),ESmin(1)
          write(ILOG,1019) '//S22',ESmax(2),ESmin(2)
          write(ILOG,1019) '//S33',ESmax(3),ESmin(3)
          write(ILOG,1019) '//S12',ESmax(4),ESmin(4)
          write(ILOG,1019) '//S23',ESmax(5),ESmin(5)
          write(ILOG,1019) '//S13',ESmax(6),ESmin(6)
        endif
        
!C
!C*for 6D
      elseif( NDOF == 6 ) then
        write(ILOG,*) '##### Local Summary :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//U1    ',Umax(1),IUmax(1),Umin(1),IUmin(1)
        write(ILOG,1009)'//U2    ',Umax(2),IUmax(2),Umin(2),IUmin(2)
        write(ILOG,1009)'//U3    ',Umax(3),IUmax(3),Umin(3),IUmin(3)
        write(ILOG,1009)'//R1    ',Umax(4),IUmax(4),Umin(4),IUmin(4)
        write(ILOG,1009)'//R2    ',Umax(5),IUmax(5),Umin(5),IUmin(5)
        write(ILOG,1009)'//R3    ',Umax(6),IUmax(6),Umin(6),IUmin(6)
        write(ILOG,*) '##### @Element :Max/IdMax/Min/IdMin####'
        write(ILOG,1009)'//E11'                                        &
                              ,EEmax( 1),IEEmax( 1),EEmin( 1),IEEmin( 1)
        write(ILOG,1009)'//E22'                                        &
                              ,EEmax( 2),IEEmax( 2),EEmin( 2),IEEmin( 2)
        write(ILOG,1009)'//E12'                                        &
                              ,EEmax( 3),IEEmax( 3),EEmin( 3),IEEmin( 3)
        write(ILOG,1009)'//E23'                                        &
                              ,EEmax( 4),IEEmax( 4),EEmin( 4),IEEmin( 4)
        write(ILOG,1009)'//E13'                                        &
                              ,EEmax( 5),IEEmax( 5),EEmin( 5),IEEmin( 5)
        write(ILOG,1009)'//S11'                                        &
                              ,ESmax( 1),IESmax( 1),ESmin( 1),IESmin( 1)
        write(ILOG,1009)'//S22'                                        &
                              ,ESmax( 2),IESmax( 2),ESmin( 2),IESmin( 2)
        write(ILOG,1009)'//S12'                                        &
                              ,ESmax( 3),IESmax( 3),ESmin( 3),IESmin( 3)
        write(ILOG,1009)'//S23'                                        &
                              ,ESmax( 4),IESmax( 4),ESmin( 4),IESmin( 4)
        write(ILOG,1009)'//S13'                                        &
                              ,ESmax( 5),IESmax( 5),ESmin( 5),IESmin( 5)
        write(ILOG,1009)'//SMS'                                        &
                              ,ESmax( 6),IESmax( 6),ESmin( 6),IESmin( 6)

!C*** Show Summary
        call mw_allreduce_r( Umax, 6, mw_mpi_max() )
        call mw_allreduce_r( Umin, 6, mw_mpi_min() )
        call mw_allreduce_r( EEmax, 5, mw_mpi_max() )
        call mw_allreduce_r( EEmin, 5, mw_mpi_min() )
        call mw_allreduce_r( ESmax, 6, mw_mpi_max() )
        call mw_allreduce_r( ESmin, 6, mw_mpi_min() )
		

        if( myrank == 0 ) then
          write(ILOG,*) '##### Global Summary :Max/Min####'
          write(ILOG,1019) '//U1    ',Umax(1),Umin(1)
          write(ILOG,1019) '//U2    ',Umax(2),Umin(2)
          write(ILOG,1019) '//U3    ',Umax(3),Umin(3)
          write(ILOG,1019) '//R1    ',Umax(4),Umin(4)
          write(ILOG,1019) '//R2    ',Umax(5),Umin(5)
          write(ILOG,1019) '//R3    ',Umax(6),Umin(6)
          write(ILOG,*) '##### @Element :Max/Min####'
          write(ILOG,1019) '//E11',EEmax( 1),EEmin( 1)
          write(ILOG,1019) '//E22',EEmax( 2),EEmin( 2)
          write(ILOG,1019) '//E12',EEmax( 3),EEmin( 3)
          write(ILOG,1019) '//E23',EEmax( 4),EEmin( 4)
          write(ILOG,1019) '//E13',EEmax( 5),EEmin( 5)

          write(ILOG,1019) '//S11',ESmax( 1),ESmin( 1)
          write(ILOG,1019) '//S22',ESmax( 2),ESmin( 2)
          write(ILOG,1019) '//S12',ESmax( 3),ESmin( 3)
          write(ILOG,1019) '//S23',ESmax( 4),ESmin( 4)
          write(ILOG,1019) '//S13',ESmax( 5),ESmin( 5)
          write(ILOG,1019) '//SMS',ESmax( 6),ESmin( 6)

        endif
      endif
        
 1009 format(a8,1pe12.4,i10,1pe12.4,i10)
 1019 format(a8,1pe12.4,1pe12.4)
      end subroutine FSTR_POST
end module m_static_post
