!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 4.3                                   !
!                                                                      !
!      Module Name : Eigen Analysis                                    !
!                                                                      !
!            Written by Yasuji Fukahori (Univ. of Tokyo)               !
!                       Giri Prabhakar (RIST)                          !
!                                                                      !
!      Contact address :  IIS,The University of Tokyo, CISS            !
!                                                                      !
!      "Structural Analysis for Large Scale Assembly"                  !
!                                                                      !
!======================================================================!


!======================================================================!
!                       Description                                    !
!======================================================================!
!> \brief This module contains data definitions for Lanczos eigenvalue solver
!=======================================================================
!----------------------------------------------------------------------*
!     Definitions of module
!----------------------------------------------------------------------*
      module lczeigen
      use hecmw
      USE lczparm
      IMPLICIT NONE
      PUBLIC

!Work arrays (EWK,EVEC carry eigenvectors, EVAL carries eigenvalues)

      INTEGER (KIND=KINT) :: IERR
      REAL (KIND=KREAL), DIMENSION (:), pointer :: MASS, EVAL
      REAL (KIND=KREAL), DIMENSION (:), pointer :: WORK, EGV0
      REAL (KIND=KREAL), DIMENSION (:,:), pointer :: EWK
      REAL (KIND=KREAL), DIMENSION (:,:), POINTER :: EVEC
      INTEGER (KIND=KINT), DIMENSION (:), pointer :: NEW
      REAL (KIND=KREAL), DIMENSION (:,:), pointer :: XMODE

!Lanczos parameters
      TYPE(lczvec) :: lvecq(0:lvecq_size)      !< Array of Q vectors
      REAL(KIND=KREAL), pointer ::  LVECP(:)   !< Array of P vectors
      REAL(KIND=KREAL), pointer ::  LVECPP(:)  !< Array of modified P vectors
      REAL(KIND=KREAL), pointer ::  ALF(:),LWRK(:)
      REAL(KIND=KREAL), pointer ::  BTA(:)     !< Beta parameter
      REAL(KIND=KREAL), pointer ::  LLWRK(:),LLLWRK(:)
      REAL(KIND=KREAL), ALLOCATABLE ::  LLDIAG(:), LNDIAG(:), LSUB(:)
      REAL(KIND=KREAL), ALLOCATABLE ::  LZMAT(:,:), LNZMAT(:,:)
      REAL(KIND=KREAL), pointer ::  EFILT(:)
      INTEGER(kind=kint), pointer ::  modal(:)
      REAL(KIND=KREAL), pointer ::  EVECQ(:,:),EM(:)
      REAL(KIND=KREAL) SGM                     !< Shift
      REAL(KIND=KREAL) CTOL                    !< Tolerance
      REAL(KIND=KREAL) prechk
      LOGICAL CONV                             !< Convergence flag

      INTEGER(KIND=KINT) :: ITLIMIT,NEIG,NGET,ITER,IERROR,LTRIAL
      REAL(KIND=KREAL) :: PRECHK1,CCHK0,CCHK1,CCHK,CERR

!*------------------- Extra variables for FSTR ---------------------------*
      INTEGER (KIND=kint) numnp,NDOF,ntotal,iv

!*------------------- For more than 1 CPU --------------------------------*
      INTEGER(KIND=kint), POINTER :: my_ntotal(:)
      INTEGER(KIND=kint) numn,Gntotal,Gtotal,novl,iovl

!----------------------- For Multiple Eigenvalues ------------------------!
      real(kind=kreal), pointer :: emwk(:,:,:)
      real(kind=kreal), pointer :: ewt(:)
      real(kind=kreal), pointer :: ev(:)
      real(kind=kreal)  res
      integer(kind=kint) isc, nume, nmode, it, kcount
      logical lczmult
      character(len=32) :: fname

      END MODULE lczeigen
