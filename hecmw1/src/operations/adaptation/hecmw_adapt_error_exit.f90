!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 1.00                                               !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Adaptive Mesh Refinement                           !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using Hight End Computing Middleware (HEC-MW)"      !
!                                                                      !
!======================================================================!

      subroutine hecmw_adapt_ERROR_EXIT (hecMESH, IFLAG)
      use  hecmw_util
      type (hecmwST_local_mesh) :: hecMESH

      write (*,'(/,a, i8)')                                             &
     &        "********** MESSAGE from HEC-MW ADAPTATION **********",   &
     &         hecMESH%my_rank

      if (IFLAG.ge. 1001) then
        write (*,'(/,a)')                                               &
     &        " ### ABORT : unexpected ZERO/minus in the orginal file"
        if (IFLAG.eq.1001) write (*,'(  a,/)')                          &
     &        "       TOTAL NODE and/or ELEMENT NUMBER"
        if (IFLAG.eq.1002) write (*,'(  a,/)')                          &
     &        "       BOUNDARY GROUP NUMBER"
        if (IFLAG.eq.1003) write (*,'(  a,/)')                          &
     &        "       BOUNDARY info ITEMs"
        if (IFLAG.eq.1004) write (*,'(  a,/)')                          &
     &        "       ELEMENT type"
        if (IFLAG.eq.1005) write (*,'(  a,/)')                          &
     &        "       ELEMENT connectivity"
      endif

      if (IFLAG.eq. 1) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE : Parallel Info"
      endif

      if (IFLAG.eq. 2) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE : Parallel Info"
      endif

      if (IFLAG.eq.11) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE"
      endif

      if (IFLAG.eq.12) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : UNEXPECTED EOF in ORIGINAL GRID FILE"
      endif

      if (IFLAG.eq.21) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in MeTiS PARTIRION FILE"
      endif

      if (IFLAG.eq.99) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR allocating ieaddrs in load_mesh.f90"
      endif

      if (IFLAG.eq. 101) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : file not found : GRID FILE"
      endif

      if (IFLAG.eq. 102) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : file not found : RESTART FILE"
      endif

      if (IFLAG.eq. 20) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : invalid file format : CNTL FILE"
      endif

      if (IFLAG.eq. 22) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : invalid file format : RESTART FILE"
      endif

      if (IFLAG.eq. 94) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : invalid value : CNTL. FILE"
      endif

      if (IFLAG.eq.61) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : array overflow (CELL)"
      endif

      if (IFLAG.eq.62) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : array overflow (NODES)"
      endif

      if (IFLAG.eq.63) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : array overflow (LAYER)"
      endif

      if (IFLAG.eq.64) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : array overflow (EDGES)"
      endif

      if (IFLAG.eq.81) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : arithmetic overflow : Zero Volume"
      endif

      if (IFLAG.eq.82) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : arithmetic overflow : Zero Density"
      endif

      if (IFLAG.eq.9) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : adaptation level exceeded"
      endif

      if (IFLAG.eq.95) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : inconsistent grid-result files"
      endif

      if (IFLAG.eq.7) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : iteration failed"
      endif

      if (IFLAG.eq.201) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : initial conditions"
      endif

      if (IFLAG.eq.202) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : element-based partition"
      endif

      if (IFLAG.eq.203) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : multiple DOF group"
      endif

      if (IFLAG.eq.204) then
        write (*,'(/,a,/)')                                             &
     &       " ### ABORT : composite material"
      endif


      call MPI_FINALIZE (                  ier1      )

      stop
      end
  

