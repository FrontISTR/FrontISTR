!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.5                                   !
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
!> This program provide a function to report timing measurements for the Lanczos Eigensolver routines
module m_time_log
contains


       SUBROUTINE time_log(hecMESH,hecMAT,myEIG)

       USE m_fstr
       USE lczeigen
       USE lczparm

       IMPLICIT NONE
       TYPE (hecmwST_local_mesh) :: hecMESH
       TYPE (hecmwST_matrix)     :: hecMAT
       TYPE (lczparam)           :: myEIG
       


!*Initialize
      IF(tinit) THEN
       tasbts = 0.
       tasbte = 0.
       tasbus = 0.
       tasbue = 0.
       tasbss = 0.
       tasbse = 0.
       tmainte = 0.
       tmaints = 0.
       tmainse = 0.
       tmainss = 0.
       tmainue = 0.
       tmainus = 0.
       tsolavt = 0.
       tsolavu = 0.
       tsolavs = 0.
       tlczavt = 0.
       tlczavu = 0.
       tlczavs = 0.
       treavt  = 0.
       treavu  = 0.
       treavs  = 0.
       tsolts  = 0.
       tsolte  = 0.
       tsolus  = 0.
       tsolue  = 0.
       tsolss  = 0.
       tsolse  = 0.
       tlczts  = 0.
       tlczte  = 0.
       tlczus  = 0.
       tlczue  = 0.
       tlczss  = 0.
       tlczse  = 0.
       treorts = 0.
       treorte = 0.
       treorus = 0.
       treorue = 0.
       treorss = 0.
       treorse = 0.
       tinit = .FALSE.
       RETURN
      ENDIF
      
      IF(myrank .EQ. 0) THEN
       IF(teachiter) THEN
!*------------------------------- Solver timing Results --------------------------------*
       IF(myEIG%eqset.EQ.0) THEN
       WRITE(IMSG,*) '*-----------------------------------------------*'
       WRITE(IMSG,*) '*                TIMING SUMMARY                 *'
       WRITE(IMSG,*) '*-----------------------------------------------*'
       ELSE IF(myEIG%eqset.EQ.1) THEN
       WRITE(IMSG,*) '*-----------------------------------------------*'
       WRITE(IMSG,*) '* TIMING FOR LANCZOS ITERATION NUMBER: ',ITER ,'*'
       WRITE(IMSG,*) '*-----------------------------------------------*'
       ENDIF
       WRITE(IMSG,*) '+=-=-   LINEAR SOLVER TIME INFO.    -=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tsolte - tsolts
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tsolue - tsolus
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tsolse - tsolss
       WRITE(IMSG,*) '+==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-+'
!*---------------------------- Lanczos loop timing Results -----------------------------*
       IF(myEIG%eqset.EQ.1) THEN
       WRITE(IMSG,*) '+=-=- LANCZOS LOOP TIME INFO.  -=-=-=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tlczte - tlczts
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tlczue - tlczus
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tlczse - tlczss
       WRITE(IMSG,*) '+==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-+'
!*------------------------ Reorthogonalization timing Results --------------------------*
       WRITE(IMSG,*) '+=-=-  REORTHOGONALIZATION TIME INFO. -=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)')treorte-treorts
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)')treorue-treorus
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)')treorse-treorss
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-+'
       WRITE(IMSG,*) '*----------------------------------------------*'
       ENDIF
!*--------------------------------------------------------------------------------------*
      ENDIF 

      IF(tenditer) THEN
       WRITE(IMSG,*) 
       WRITE(IMSG,*) 
       WRITE(IMSG,*) '*-----------------------------------------------*'
       WRITE(IMSG,*) '*               TIMING SUMMARY                  *'
       WRITE(IMSG,*) '*-----------------------------------------------*'
!*------------------------------- Stiffness assembly time  -----------------------------*
       WRITE(IMSG,*) '+=-=-STIFFNESS ASSEMBLY TIME INFO. -=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tasbte-tasbts
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tasbue-tasbus
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tasbse-tasbss
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-+'
!*-------------------------- Lanczos Routine Time --------------------------------------*
       IF(myEIG%eqset.EQ.1) THEN
       WRITE(IMSG,*) '+=-=-MAIN LANCZOS ROUTINE TIME INFO. -=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tmainte-tmaints
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tmainue-tmainus
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tmainse-tmainss
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-+'
!*------------------------------- Solver timing Results --------------------------------*
       IF(hecMAT%Iarray(99) .EQ. 1) THEN
       WRITE(IMSG,*) '+=-=- ILU SOLVER AVERAGE TIME INFO. -=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tsolavt/LTRIAL
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tsolavu/LTRIAL
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tsolavs/LTRIAL
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-+'
       ELSE IF(hecMAT%Iarray(99) .EQ. 2) THEN
       WRITE(IMSG,*) '+=-=- DIRECT SOLVER FACTOR+SOLVE TIME INFO.-=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tsolavt
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tsolavu
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tsolavs
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-+'
       ENDIF
!*---------------------------- Lanczos loop timing Results -----------------------------*
       WRITE(IMSG,*) '+=-=- LANCZOS LOOP AVERAGE TIME INFO. -=-+'
       IF (hecMAT%Iarray(99).EQ.2) THEN
         WRITE(IMSG,*) '== For direct solver factorization time is excluded' &
     & //' from Lanczos loop time ==' 
       ENDIF
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') tlczavt/LTRIAL
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') tlczavu/LTRIAL
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') tlczavs/LTRIAL
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-+'
!*------------------------ Reorthogonalization timing Results --------------------------*
       WRITE(IMSG,*) '+=-  REORTHO. AVERAGE TIME INFO.  -=-=-+'
       WRITE(IMSG,'(''CPU TIME TOTAL (SEC)='',F10.2)') treavt/LTRIAL
       WRITE(IMSG,'(''CPU TIME USER   (SEC)='',F10.2)') treavu/LTRIAL
       WRITE(IMSG,'(''CPU TIME SYSTEM  (SEC)='',F10.2)') treavs/LTRIAL
       WRITE(IMSG,*) '+=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-+'
       ENDIF
!*--------------------------------------------------------------------------------------*
      ENDIF
      ENDIF

      RETURN
      END subroutine

end module m_time_log
 
