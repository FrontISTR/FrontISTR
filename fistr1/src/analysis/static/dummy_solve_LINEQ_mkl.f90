!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module provides functions to stop contact calculations
!> \with Lagrage multipliers.

   subroutine pardiso(pt, maxfct, mnum, mtype, phase, ntdf, values, pointers, indices, &
                       idum, nrhs, iparm, msglvl, ddum1, ddum2, ierr)

   use m_fstr

!< ------------------------------------------------ input parameters for MKL solver(see Intel(R) MKL Reference Manual)
   integer (kind=8)              :: pt(64)
   integer (kind=kint)           :: maxfct
   integer (kind=kint)           :: mnum
   integer (kind=kint)           :: mtype
   integer (kind=kint)           :: phase
   integer (kind=kint)           :: ntdf
   integer (kind=kint)           :: pointers(:)        !< ia
   integer (kind=kint)           :: indices(:)         !< ja
   real     (kind=kreal)          :: values(:)          !< a
   integer (kind=kint)           :: idum
   integer (kind=kint)           :: nrhs
   integer (kind=kint)           :: iparm(64)
   integer (kind=kint)           :: msglvl
   real(kind=kreal)               :: ddum1, ddum2
   integer (kind=kint)           :: ierr

   write(*,*) "paradiso: ERROR was detected. Please install MKL library and try again."
   stop

   end subroutine
