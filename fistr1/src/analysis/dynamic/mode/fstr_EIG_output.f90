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

!> Display results in .out files
module m_fstr_EIG_output
contains

!C======================================================================
!C----------------------------------------------------------------------
      subroutine fstr_eigen_output(hecMESH,hecMAT,IOUT,myEIG)
!C----------------------------------------------------------------------
!C*
!C* SHOW RESULTS
!C*
      use m_fstr
      use lczparm
      use lczeigen
      use m_fstr_EIG_lanczos
      use m_fstr_EIG_matmult
      implicit none
      integer(kind=kint) :: I, IOUT,IREOR,JJITER,JITER,IITER
      real(kind=kreal) :: prechk0(1)
      type (hecmwST_local_mesh) :: hecMESH
      type (hecmwST_matrix    ) :: hecMAT
      type (fstr_solid       )  :: fstrSOLID
      type (hecmwST_result_data):: fstrRESULT
      type (lczparam) :: myEIG

!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if

      CALL EGLIST(hecMESH,hecMAT,IOUT)
      IF(myrank.EQ.0) THEN
        WRITE(IMSG,*) ''
        WRITE(IMSG,*) '*----------------------------------------------*'
        WRITE(IMSG,*) '*--E I G E N V A L U E  C O N V E R G E N C E--*'
        WRITE(IMSG,*) '*----------------------------------------------*'
        WRITE(IMSG,*) 'Absolute residual =  |(||Kx - lambda*Mx||)|'
      ENDIF
!C
!Normalize eigenvectors
      IF(myrank==0) THEN
        WRITE(IMSG,*) '*------------- Residual Check ----------------*'
        WRITE(IMSG,*) '   Iter.#   Eigenvalue    Abs. Residual  '
        WRITE(IMSG,*) '   *-----*  *---------*  *--------------*'
      ENDIF
!C
      LLLWRK = 0.0D0
      DO JJITER = 1,myEIG%nget
        JITER = NEW(JJITER)

        if( modal(JITER).eq.1 ) then
          LWRK = 0.0D0
!!!          CALL VECPRO( prechk1,ewk(1,JJITER),ewk(1,JJITER),Gntotal,1 )
          prechk1=0.0
          do i = 1, Gntotal
            prechk1=prechk1+ewk(i,JJITER)**2
          enddo
          if (.not. ds) then !In case of Direct Solver prevent MPI
            CALL hecmw_allreduce_R1( hecMESH,prechk1,hecmw_sum )
          end if
          prechk1 = sqrt(prechk1)
          if( prechk1.NE.0.0D0 ) ewk(:,JJITER) = ewk(:,JJITER)/prechk1
!C
          DO IITER = 1, ntotal
            LWRK(IITER) = ewk(IITER,JJITER)
          ENDDO
!C
          IF(NDOF.EQ.3) THEN
            CALL MATMULT3( hecMESH,hecMAT,LLLWRK,1.0D0,LWRK,numnp,NDOF )
          ELSE IF(NDOF.EQ.2) THEN
            CALL MATMULT2( hecMESH,hecMAT,LLLWRK,1.0D0,LWRK,numnp,NDOF )
          ELSE IF(NDOF.EQ.6) THEN
            CALL MATMULT6( hecMESH,hecMAT,LLLWRK,1.0D0,LWRK,numnp,NDOF )
          ENDIF
!C
          LLWRK = 0.0D0
!!!          CALL MATPRO(LLWRK,mass,ewk(1,JJITER),ntotal,1)
          do i = 1, ntotal
            LLWRK(i) = mass(i)*ewk(i,JJITER)
          enddo
!C
          CCHK1 = 0.0D0
          DO IITER = 1,Gntotal
            CCHK1 = CCHK1 + ( LWRK(IITER) - (eval(JITER) &
     &                      - myEIG%lczsgm)*LLWRK(IITER) )**2
          END DO
!C
          if (.not. ds) then !In case of Direct Solver prevent MPI
            CALL hecmw_allreduce_R1(hecMESH,CCHK1,hecmw_sum)
          end if
          CCHK1 = SQRT(CCHK1)
          IF(myrank.EQ.0) THEN
             WRITE(IMSG,'(2X,I5,2X,5E15.6)') JITER,eval(JITER),CCHK1
          ENDIF
        endif
      ENDDO
!C
      IF(myrank==0) THEN
        WRITE(IMSG,*) ''
        WRITE(IMSG,*)'*    ---END Eigenvalue listing---     *'
        WRITE(IMSG,*) ''
      ENDIF

      end subroutine fstr_eigen_output

end module m_fstr_EIG_output


