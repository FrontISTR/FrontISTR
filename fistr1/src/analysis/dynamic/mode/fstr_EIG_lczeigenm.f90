!-------------------------------------------------------------------------------
! Copyright (c) 2016 The University of Tokyo
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!> \brief This module contains data definitions for Lanczos eigenvalue solver

      module lczeigen
      use hecmw
      use m_fstr
      IMPLICIT NONE
      PUBLIC

      REAL (KIND=KREAL), DIMENSION (:), pointer :: WORK
      REAL (KIND=KREAL), DIMENSION (:,:), pointer :: EWK

      REAL(KIND=KREAL), pointer ::  EM(:)

      END MODULE lczeigen
