!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
!                                                                      !
!      Module Name : Steady-state harmonic response analysis           !
!                                                                      !
!            Written by Kuniaki Koike (ASTOM)                          !
!                                                                      !
!                                                                      !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!
!> \brief This module contains data structure for frequency analysis

module m_fstr_freqdata
use hecmw
  implicit none

!C
!C-- Frequency analysis parameter data structure
!C
  type fstr_freqanalysis
    integer(kind=kint)                :: FLOAD_ngrp_tot
    integer(kind=kint), pointer      :: FLOAD_ngrp_GRPID(:) => NULL()
    integer(kind=kint), pointer      :: FLOAD_ngrp_ID(:)    => NULL()
    integer(kind=kint), pointer      :: FLOAD_ngrp_TYPE(:)  => NULL()
    integer(kind=kint), pointer      :: FLOAD_ngrp_DOF(:)   => NULL()
    real(kind=kreal), pointer         :: FLOAD_ngrp_valre(:) => NULL()
    real(kind=kreal), pointer         :: FLOAD_ngrp_valim(:) => NULL()
    character(len=HECMW_FILENAME_LEN) :: eigenlog_filename
    integer(kind=kint)                :: start_mode
    integer(kind=kint)                :: end_mode
  end type
  
  type fstr_freqanalysis_data
    integer(kind=kint)        :: numMode
    integer(kind=kint)        :: numNodeDOF
    real(kind=kreal), pointer :: eigOmega(:)    => NULL()
    real(kind=kreal), pointer :: eigVector(:,:) => NULL()
    real(kind=kreal)           :: rayAlpha, rayBeta
  end type
  
  integer, parameter :: kFLOADTYPE_NODE = 1
  integer, parameter :: kFLOADTYPE_SURF = 2
  
  integer, parameter :: kFLOADCASE_RE = 1
  integer, parameter :: kFLOADCASE_IM = 2
contains

!C
!C-- initialize fstr_freqanalysis structure  
!C
  subroutine fstr_nullify_fstr_freqanalysis( f )
  !---- args
    type( fstr_freqanalysis ), intent(inout) :: f
  !---- body
    f%FLOAD_ngrp_tot = 0
    nullify( f%FLOAD_ngrp_GRPID )
    nullify( f%FLOAD_ngrp_ID )
    nullify( f%FLOAD_ngrp_TYPE )
    nullify( f%FLOAD_ngrp_DOF )
    nullify( f%FLOAD_ngrp_valre )
    nullify( f%FLOAD_ngrp_valim )
    
  end subroutine
end module
