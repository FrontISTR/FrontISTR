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

      REAL (KIND=KREAL), DIMENSION (:), pointer :: WORK, EGV0
      REAL (KIND=KREAL), DIMENSION (:,:), pointer :: EWK
      REAL (KIND=KREAL), DIMENSION (:,:), POINTER :: EVEC
      INTEGER (KIND=KINT), DIMENSION (:), pointer :: NEW
      REAL (KIND=KREAL), DIMENSION (:,:), pointer :: XMODE

      REAL(KIND=KREAL), pointer ::  LVECP(:)   !< Array of P vectors
      REAL(KIND=KREAL), pointer ::  LVECPP(:)  !< Array of modified P vectors
      REAL(KIND=KREAL), pointer ::  LLWRK(:),LLLWRK(:),LWRK(:)
      REAL(KIND=KREAL), pointer ::  EVECQ(:,:),EM(:)


      END MODULE lczeigen
