!======================================================================!
!                                                                      !
! Software Name : FrontISTR Ver. 3.7                                   !
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
!   Record of revision:
!      Date            Programmer           Description of change
!    =========         ==========           =====================
!    Jan 6,2012        YUAN Xi       Eigen analysis considering effects
!                                          of former static step
!======================================================================!
!> This module provides a function to control eigen analysis
module m_fstr_solve_eigen
contains

!C***
!C*** SOLVE EIGENVALUE PROBLEM
!C***
!C
      subroutine fstr_solve_eigen(hecMESH,hecMAT,myEIG,fstrSOLID, &
     &                            fstrRESULT,fstrPARAM, fstrMAT)
      use m_fstr
      use m_fstr_StiffMatrix
      use m_fstr_AddBC
      use fstr_matrix_con_contact
      use hecmw_util
!C*-------- Modules for Lanczos method --------------*
      use lczparm
      use lczeigen
      use m_eigen_lib
      use m_fstr_EIG_getamat
      use m_fstr_EIG_lanczos
      use m_fstr_EIG_matmult
      use m_fstr_EIG_mgs1
      use m_fstr_EIG_output
      use m_fstr_EIG_setMASS
      use m_fstr_EIG_tridiag
      use m_static_lib
      use m_static_mat_ass
      use m_static_make_result
      use m_hecmw2fstr_mesh_conv

      implicit none

      type (hecmwST_local_mesh ) :: hecMESH
      type (hecmwST_matrix     ) :: hecMAT
      type (fstr_solid         ) :: fstrSOLID
      type (hecmwST_result_data) :: fstrRESULT
      type (fstr_param         ) :: fstrPARAM
      type (fstrST_matrix_contact_lagrange)  :: fstrMAT   !< type fstrST_matrix_contact_lagrange

      character(len=HECMW_HEADER_LEN) :: header
      character(len=HECMW_NAME_LEN)   :: label
      character(len=HECMW_NAME_LEN)   :: nameID

!C*-------- Parameters for Lanczos method -----------*
      type (lczparam) :: myEIG

      integer(kind=kint) :: i, j, k , ii, iii, ik, in, in1, in2, in3, nstep, istep
      integer(kind=kint) :: ig, ig0, is0, ie0, its0, ite0, jiter, iiter, kiter
      integer(kind=kint) :: kk, jjiter, ppc
      real(kind=kreal)   :: t1, t2, aalf, tmp, tmp2, gm, gm2, r1, r2, r3, r4, r5, r6

      integer(kind=kint) IOUT,IREOR,eITMAX,itype,iS,iE,ic_type,icel,jS,nn
      integer(kind=kint), allocatable :: isnode33(:)
      real(kind=kreal), allocatable   :: gmass(:)

!C*-------- solver control -----------*
      logical :: ds = .false. !using Direct Solver or not

! in case of direct solver
      if (hecMAT%Iarray(99) .eq. 2) then
        ds = .true.
      end if

!C -------- 22:11 2006/05/27 by n.imai
      do i=0,lvecq_size
        nullify( lvecq(i)%q )
      enddo
!C ------------------------------------

!C*------------ Number of DOF ----------------------*
      ALLOCATE(my_ntotal(0:nprocs), STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"

      IF(myrank.GT.nprocs) THEN
        WRITE(*,*) '+==================================+'
        WRITE(*,*) ' Array my_ntotal needs larger size!'
        WRITE(*,*) '+=================================+'
        WRITE(*,*) ' ...Stopping in fstr_solve_eigen...'
        IF( myrank == 0 ) then
          WRITE(IMSG,*) '+==================================+'
          WRITE(IMSG,*) ' Array my_ntotal needs larger size!'
          WRITE(IMSG,*) '+=================================+'
          WRITE(IMSG,*) ' ...Stopping in fstr_solve_eigen...'
        ENDIF
        CALL hecmw_abort(hecmw_comm_get_comm())
        STOP
      ENDIF

      t1 = hecmw_Wtime() ! DEBUG elap

      numnp  = hecMAT%NP
      numn   = hecMAT%N
      NDOF   = hecMESH%n_dof
      ntotal = numnp*NDOF

!C*------------ Assemble stiffness matrix ----------*
      if( fstrPARAM%solution_type==kstSTATICEIGEN ) then
        fstrSOLID%dunode = 0.d0
        call fstr_StiffMatrix( hecMESH, hecMAT, fstrSOLID, 0.d0 )
        call fstr_AddBC(1, 1, hecMESH, hecMAT, fstrSOLID, fstrPARAM, fstrMAT, 2)
      else
        call fstr_mat_ass(hecMESH, hecMAT, fstrSOLID)
        call fstr_mat_ass_load(hecMESH, hecMAT, myEIG, fstrSOLID)
        call fstr_mat_ass_bc(hecMESH, hecMAT, fstrSOLID)
        call fstr_mat_ass_check_rhs(hecMESH, hecMAT, myEIG)
      endif

      if( myrank == 0 ) then
         write(IMSG,*) 'fstr_mat_ass: OK'
      endif

!C*----------- Check node33 -----*
      allocate(isnode33(numnp))
      isnode33 = 0

      do itype = 1, hecMESH%n_elem_type
        iS = hecMESH%elem_type_index(itype-1) + 1
        iE = hecMESH%elem_type_index(itype  )
        ic_type = hecMESH%elem_type_item(itype)
        if(hecmw_is_etype_33struct(ic_type))then
          nn = HECMW_get_max_node(ic_type)/2
          do icel = iS, iE
            jS = hecMESH%elem_node_index(icel-1)
            do j = 1, nn
              ii = hecMESH%elem_node_item(jS+j+nn)
              isnode33(ii)=1
            enddo
          enddo
        endif
      enddo

!C*----------- Alloc eigenvalue related arrays -----*
      allocate( EFILT( ntotal) )

!C*----------- EHM 7Apr04: Memory monitor ----------*
      CALL memget(mlczi,ntotal,4)

!C* SPC VECTOR FILTER: ignore rotational inertia at present.

      efilt  = 1.0
      kcount = 0

      do ig0 = 1, fstrSOLID%BOUNDARY_ngrp_tot
        ig   = fstrSOLID%BOUNDARY_ngrp_ID(ig0)
        iS0  = hecMESH%node_group%grp_index(ig-1) + 1
        iE0  = hecMESH%node_group%grp_index(ig  )
        it   = fstrSOLID%BOUNDARY_ngrp_type(ig0)
        itS0 = (it - mod(it,10))/10
        itE0 = mod(it,10)

        do ik = iS0, iE0
          in = hecMESH%node_group%grp_item(ik)
          do i = itS0,itE0
            EFILT((in-1)*NDOF+i)=0.0
          enddo
        enddo

      enddo

      if (.not. ds) then !In case of Direct Solver prevent MPI
        IF(NDOF.eq.3) THEN
          CALL hecmw_update_3_R(hecMESH,EFILT,numnp)
        ELSE IF(NDOF.eq.2) THEN
          CALL hecmw_update_2_R(hecMESH,EFILT,numnp)
        ELSE IF(NDOF.eq.6) THEN
          CALL hecmw_update_m_R(hecMESH,EFILT,numnp,NDOF)
        ENDIF
      end if

      Gtotal  = numnp*NDOF
      Gntotal = numn*NDOF
      novl    = ( numnp - numn )*NDOF

      kcount = 0
      DO I = 1, Gntotal
        IF(EFILT(I).EQ.0) kcount = kcount + 1
      END DO

      if (.not. ds) then !In case of Direct Solver prevent MPI
        call hecmw_allreduce_I1(hecMESH,Gtotal,hecmw_sum)
        call hecmw_allreduce_I1(hecMESH,kcount,hecmw_sum)
      end if

      eITMAX  = myEIG%lczmax
      ITLIMIT = Gtotal - kcount
      IF(eITMAX.GT.ITLIMIT) THEN
        IF(myrank .EQ. 0) THEN
          WRITE(IMSG,*) '*-------------------------------------------*'
          WRITE(IMSG,*) '  WARNING: LCZMAX exceeds system matrix size.'
          WRITE(IMSG,*) '  Resetting LCZMAX to system matrix size.'
          WRITE(IMSG,*) '*-------------------------------------------*'
        ENDIF
        eITMAX = ITLIMIT
      ENDIF

      CTOL = myEIG%lcztol
      NGET = myEIG%nget
      neig = eITMAX + NGET

!C* Allocate mass matrix
      ALLOCATE(myEIG%mass(ntotal), STAT=ierror)
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      write(IDBG,*) '*Allocated mass matrix of size: ',ntotal

!C*----------- EHM 7Apr04: Memory monitor ----------*
      CALL memget(mlczr,ntotal,8)

!C*-------------------- Memory allocations -----------------------*
      allocate( mass( ntotal ),            STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( EM( ntotal ),              STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( ewk( ntotal, neig),        STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( eval( neig ),              STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( modal( neig ),             STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( new( neig ),               STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( work( neig*(3*neig + 5) ), STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( LVECP(NTOTAL),             STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( LVECPP(NTOTAL),            STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( lvecq(0)%q(ntotal),        STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( lvecq(1)%q(ntotal),        STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( LWRK(NTOTAL),              STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( LLWRK(NTOTAL),             STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( LLLWRK(NTOTAL),            STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( ALF(NEIG+2),               STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"
      allocate( BTA(NEIG+2),               STAT=ierror )
      IF(ierror.NE.0) STOP "Allocation error, fstr_solve_eigen"

!C*----------- EHM 7Apr04: Memory monitor ----------*
      I = 2*neig
      CALL memget(mlczi,I,4)
      I = (nget+9)*ntotal
      CALL memget(mlczr,I,8)
      I = 10*neig+3*neig**2+4
      CALL memget(mlczr,I,8)
      I = (2+(neig+1)*myEIG%lczrod)*ntotal
      CALL memget(mreorr,I,8)

!C*------------ Zero clear --------------------------------------------*
      mass       = 0.0
      ewk        = 0.0
      eval       = 0.0
      new        = 0.0
      work       = 0.0
      lwrk       = 0.0
      LVECP      = 0.0
      LVECPP     = 0.0
      lvecq(0)%q = 0.0
      lvecq(1)%q = 0.0
      LWRK       = 0.0
      LLWRK      = 0.0
      ALF        = 0.0
      BTA        = 0.0

!C*-------------- Create Mass Matrix -----------------------------------*
      call setMASS(IDBG,fstrSOLID,hecMESH,hecMAT,myEIG)
      do i=1,ntotal
        mass(i) = myEIG%mass(i)
      enddo
!C
!C*------------- Set initial eigenvector guesses -----------------------*
      my_ntotal(myrank) = Gntotal
      if (.not. ds) then !In case of Direct Solver prevent MPI
        DO I = 0,nprocs-1
          CALL hecmw_bcast_I1(hecMESH,my_ntotal(I),I)
        END DO
      end if

      CALL SETIVL(mass,EM,EFILT,ewk,LVECP,lvecq,BTA,ntotal,&
     &                         neig,my_ntotal,hecMESH,hecMAT,NDOF,Gtotal)

!C*------------ Main iteration loop for Lanczos method -----------------*
      CONV = .FALSE.
!C
      do i=1,ntotal
        hecMAT%B(i) = EM(i)
      enddo
      hecMAT%X = 0.
      hecMAT%Iarray(98) = 1   !Assmebly complete
      hecMAT%Iarray(97) = 1   !Need numerical factorization

      hecMAT%Rarray(2)  = 1.0
      hecMAT%Rarray(1)  = myEIG%iluetol

      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *****   STAGE Begin Lanczos loop     **'
      ENDIF
!C
!C*---------- Begin Main loop ----------------
      DO ITER=1,eITMAX

!C*---------- Copy p vector into r
        CALL DUPL(EM,LVECP,ntotal)

!C*---------- Set up RHS for solve
        call hecmw_mat_clear_b(hecMAT)
        do i=1,ntotal
          hecMAT%B(i) = EM(i)
        enddo
        hecMAT%X = 0.0d0

!C*---------- Call solver (direct/iterative)
        CALL solve_LINEQ( hecMESH,hecMAT,IMSG )

!C*---------- Get back RHS VECTOR
        do i=1,ntotal
          EM(i) = hecMAT%X(i)
        enddo

!C*---------- SPC vector filter
        i = 0
        do ii = 1,ntotal
          IF(efilt(ii) .EQ. 0) i = i + 1
        end do

        do ii = 1,ntotal
          EM(ii) = EM(ii)*EFILT(ii)
        end do

!C*---------- Allocate fresh lvecq
        ALLOCATE( lvecq(iter+1)%q(ntotal), STAT=ierror )
        IF(ierror.NE.0) THEN
          IF( myrank==0 ) THEN
             WRITE(IMSG,*)"*------- Lanczos Solver Message --------------*"
             WRITE(IMSG,*)"Allocation failed for additional lanczos vector"
             WRITE(IMSG,*)"Please increase memory, or decrease NGET."
             WRITE(IMSG,*)"*---------------------------------------------*"
          ENDIF
          WRITE(*,*)"*------- Lanczos Solver Message --------------*"
          WRITE(*,*)"Allocation failed for additional lanczos vector"
          WRITE(*,*)"Please increase memory, or decrease NGET."
          WRITE(*,*)"*---------------------------------------------*"
          CALL hecmw_abort(hecmw_comm_get_comm())
          STOP
        ENDIF
        CALL memget(mlczr,ntotal,8)

!C* ///  ^rj = -rj - qj-1*btaj  ///
        CALL SCSHFT(EM,lvecq(iter-1)%q,BTA(ITER),ntotal)

!C* ///  alfj = qj^T[M][^rj] = pj^T[^rj]  ///
        CALL VECPRO1(AALF,LVECP,EM,Gntotal)
        ALF(ITER)=AALF
        if (.not. ds) then !In case of Direct Solver prevent MPI
          CALL hecmw_allreduce_R1(hecMESH,ALF(ITER),hecmw_sum)
        end if

        WRITE(IDBG,*) '+-------------------------------------+'
        WRITE(IDBG,*) ' ITER=',ITER,' ALPHA=',ALF(ITER)
        WRITE(IDBG,*) '+-------------------------------------+'

!C* ///  rj = ^rj - qj alfj  ///
        CALL SCSHFT(EM,lvecq(iter)%q,ALF(ITER),ntotal)

!C*---------- Reorthogonalization
        LWRK = 0.
        CALL MATPRO(LWRK,mass,EM,ntotal,1)
        CALL VECPRO1(prechk,LWRK,LWRK,Gntotal)
        if (.not. ds) then !In case of Direct Solver prevent MPI
          CALL hecmw_allreduce_R1(hecMESH,prechk,hecmw_sum)
        end if
        prechk = sqrt(prechk)
        IF(prechk.NE.0.0D0) LWRK = LWRK/prechk
        IREOR = ITER*(1.0 - myEIG%lczrod)

        DO KK = IREOR,ITER
          prechk = 0.0D0
          CALL VECPRO1(prechk,lvecq(kk)%q,LWRK,Gntotal)
          if (.not. ds) then !In case of Direct Solver prevent MPI
            CALL hecmw_allreduce_R1(hecMESH,prechk,hecmw_sum)
          end if
          prechk1 = 0.0D0
          CALL VECPRO1(prechk1,lvecq(kk)%q,lvecq(kk)%q,Gntotal)
          if (.not. ds) then !In case of Direct Solver prevent MPI
            CALL hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)
          end if
          prechk1 = sqrt(prechk1)
          if(prechk1.ne.0.0D0) prechk = prechk/prechk1
          IF(abs(prechk).GT.myEIG%lczrot) THEN
            CALL MGS1(lvecq(kk)%q,EM,mass,ntotal,my_ntotal,myrank,&
     &                hecMESH,Gntotal)
          ENDIF
        ENDDO

!C* ///  -pj = [M][rj]  ///
        CALL MATPRO(LVECPP,mass,EM,ntotal,1)

!C* ///  btaj+1 = (rj^T[M]rj)^(1/2) = (-pj^Trj)^(1/2)  ///
        CALL VECPRO1(AALF,LVECPP,EM,Gntotal)
        BTA(ITER+1) = AALF
        if (.not. ds) then !In case of Direct Solver prevent MPI
          CALL hecmw_allreduce_R1(hecMESH,BTA(ITER+1),hecmw_sum)
        end if
        BTA(ITER+1) = SQRT(BTA(ITER+1))

        WRITE(IDBG,*) '+-------------------------------------+'
        WRITE(IDBG,*) ' ITER=',ITER,' BTA=',BTA(ITER+1)
        WRITE(IDBG,*) '+-------------------------------------+'

!C*---------- Update Lanczos vectors
        CALL hecmw_barrier(hecMESH)
        CALL UPLCZ(LVECP,lvecq(iter+1)%q,LVECPP,EM,BTA(ITER+1),ntotal)

!C*---------- Set up tridiagonal Lanczos matrix
        LTRIAL = ITER
        ALLOCATE( LLDIAG(LTRIAL)        )
        ALLOCATE( LNDIAG(LTRIAL)        ) !Unordered
        ALLOCATE( LSUB(LTRIAL)          )
        ALLOCATE( LZMAT(LTRIAL,LTRIAL)  )
        ALLOCATE( LNZMAT(LTRIAL,LTRIAL) ) !Unordered

        DO JITER = 1,LTRIAL
          DO IITER = 1,LTRIAL
            LZMAT(JITER,IITER)  = 0.0D0
            LNZMAT(JITER,IITER) = 0.0D0 !Unordered
          END DO
        END DO

        DO IITER = 1,LTRIAL
          LLDIAG(IITER) = ALF(IITER)
          LZMAT(IITER,IITER) = 1.0D0
        END DO

        LSUB(1) = 0.0
        IITER   = 0
        DO IITER = 2,LTRIAL
          LSUB(IITER)  = BTA(IITER)
        END DO

        CALL TRIDIAG(LTRIAL,LTRIAL,LLDIAG,LNDIAG,&
     &                      LSUB,LZMAT,LNZMAT,IERROR) !Unordered

!C*---------- Convergence check
        DO IITER = 1, LTRIAL
          IF( LLDIAG(IITER).NE.0.0D0 ) THEN
            EVAL(IITER) = 1.0D0/LLDIAG(IITER) + myEIG%lczsgm
          ENDIF
        ENDDO
        CALL EVSORT(EVAL,NEW,LTRIAL)
!C
        CCHK  = 0.0
        KITER = NGET+2            !Extra 2 values as a safety feature
        IF(LTRIAL .LT. KITER) KITER = LTRIAL

        DO JJITER = 1,KITER
          JITER = NEW(JJITER)
          CALL VECPRO1(prechk1,LZMAT(1,JITER),LZMAT(1,JITER),LTRIAL)

          IF(prechk1 .GT. 0.) THEN
            prechk1 = SQRT(prechk1)
          ELSE
            prechk1 = 1.
          ENDIF

          CCHK1 = (BTA(JITER))*(LZMAT(LTRIAL,JITER)/prechk)
          IF(cchk .LT. ABS(cchk1)) THEN
            cchk = ABS(cchk1)
            iiter = jiter
            ppc = prechk1
          ENDIF
        END DO
!C
        IF(CCHK.LE.CTOL.AND.ITER.GE.NGET) THEN
          WRITE(IDBG,*) '*=====Desired convergence was obtained =====*'
          CONV = .TRUE.
          GO TO 25
        ENDIF

        IF(ITER.LT.eITMAX) THEN
          if( allocated(LLDIAG) ) DEALLOCATE(LLDIAG)
          if( allocated(LNDIAG) ) DEALLOCATE(LNDIAG)
          if( allocated(LSUB) )   DEALLOCATE(LSUB)
          if( allocated(LZMAT) )  DEALLOCATE(LZMAT)
          if( allocated(LNZMAT) ) DEALLOCATE(LNZMAT)
        ENDIF

!* End of Main Lanczos Loop
      END DO

  25  CONTINUE
      CERR = CCHK

      DO IITER = 1, LTRIAL
        IF(LLDIAG(IITER).NE.0.0D0) THEN
          EVAL(IITER) = 1.0D0/LLDIAG(IITER) + myEIG%lczsgm
        ENDIF
      END DO
      CALL evsort(eval,new,ltrial)

!C*---------- Estimate eigenvectors
      ewk = 0.0
      k = nget
      IF(k .GT. ltrial) k = ltrial
      DO kk = 1,k
        kiter = NEW(kk)
        DO jiter = 1,ltrial
          DO iiter =1,ntotal
            ewk(iiter,kk) = ewk(iiter,kk) &
     &                    + lvecq(jiter)%q(iiter)*LZMAT(jiter,kiter)
          ENDDO
        ENDDO
      ENDDO

!C*---------- Deallocate Lanczos vector array
      DO iiter=0,ltrial
        if( associated(lvecq(iiter)%q) ) DEALLOCATE(lvecq(iiter)%q)
      END DO

!C*---------- Temporarily set modal value to 1
      modal = 0
      do i = 1,LTRIAL
        if(EVAL(i).NE.0) modal(i) = 1
      end do

      lczmult = .FALSE.
      CALL memget(mlczr,LTRIAL*(LTRIAL+3),8)

!C********************* End Lanczos ************************!

      t2 = hecmw_Wtime() !DEBUG elap
      write(idbg,'(a,f10.2)') 'Lanczos loop (sec) :', T2 - T1 ! elap

      IF(myrank .EQ. 0) THEN
        WRITE(IMSG,*)
        WRITE(IMSG,*) ' *     STAGE Output and postprocessing    **'
      ENDIF

      DO JITER=1,NGET
        !prechk1 = 0.0
!C        CALL VECPRO1(prechk1,ewk(1:,JITER:),ewk(1:,JITER:),Gntotal)
        !do iii = 1, Gntotal
        !    prechk1 = prechk1 + mass(iii)*ewk(iii,JITER)*ewk(iii,JITER)
        !enddo
        prechk1 = maxval(ewk(:,JITER))
        if (.not. ds) then !In case of Direct Solver prevent MPI
          CALL hecmw_allreduce_R1(hecMESH,prechk1,hecmw_sum)
        end if
        !prechk1 = sqrt(prechk1)
        if(prechk1.NE.0.0D0)then
          do i = 1, Gntotal
            ewk(i,JITER) = ewk(i,JITER)/prechk1
          enddo
        endif
      END DO

!C***** compute effective mass and participation factor
      allocate(myEIG%effmass(3*NGET))
      allocate(myEIG%partfactor(3*NGET))
      myEIG%effmass    = 0.0d0
      myEIG%partfactor = 0.0d0

      DO i=1,NGET
        r1 = 0.0d0
        r2 = 0.0d0
        r3 = 0.0d0
        gm = 0.0d0
        do j = 1, numn
          !if(isnode33(j)==0)then
            in1 = 3*j-2
            in2 = 3*j-1
            in3 = 3*j
            r1 = r1 + mass(in1)*ewk(in1,i)
            r2 = r2 + mass(in2)*ewk(in2,i)
            r3 = r3 + mass(in3)*ewk(in3,i)
            gm = gm + mass(in1)*ewk(in1,i)*ewk(in1,i) &
            & + mass(in2)*ewk(in2,i)*ewk(in2,i) &
            & + mass(in3)*ewk(in3,i)*ewk(in3,i)
          !endif
        enddo
        CALL hecmw_allreduce_R1(hecMESH,r1,hecmw_sum)
        CALL hecmw_allreduce_R1(hecMESH,r2,hecmw_sum)
        CALL hecmw_allreduce_R1(hecMESH,r3,hecmw_sum)
        CALL hecmw_allreduce_R1(hecMESH,gm,hecmw_sum)
        myEIG%partfactor(3*i-2) = r1/gm
        myEIG%partfactor(3*i-1) = r2/gm
        myEIG%partfactor(3*i  ) = r3/gm
        myEIG%effmass(3*i-2) = r1*r1/gm
        myEIG%effmass(3*i-1) = r2*r2/gm
        myEIG%effmass(3*i  ) = r3*r3/gm
      END DO


!*------------------
!*Eigensolver output
!*------------------
      CALL fstr_eigen_output(hecMESH,hecMAT,ILOG,myEIG)

      nstep = myEIG%nget

      do istep = 1, nstep
        do i=1,ntotal
          hecMAT%X(i) = ewk(i,istep)
        enddo

        IF(NDOF.EQ.3) THEN
          call hecmw_update_3_R(hecMESH,hecMAT%X,hecMAT%NP)
        ELSE IF( NDOF .EQ. 2 ) THEN
          call hecmw_update_2_R(hecMESH,hecMAT%X,hecMAT%NP)
        ELSE IF(ndof.EQ.6) THEN
          call hecmw_update_m_R(hecMESH,hecMAT%X,hecMAT%NP,ndof)
        ENDIF

        if( IRESULT.eq.1 ) then
          header = "*fstrresult"
          call hecmw_result_init(hecMESH,nstep,istep,header)
          label = "DISPLACEMENT"
          call hecmw_result_add(1,NDOF,label,hecMAT%X)
          nameID = "fstrRES"
          call hecmw_result_write_by_name(nameID)
          call hecmw_result_finalize
        endif

        if( IVISUAL.eq.1 ) then
          call hecmw_nullify_result_data(fstrRESULT)
          fstrRESULT%nn_component = 1
          fstrRESULT%ne_component = 0
          allocate(fstrRESULT%nn_dof(1))
          allocate(fstrRESULT%node_label(1))
          allocate(fstrRESULT%node_val_item(NDOF*hecMAT%NP))
          fstrRESULT%nn_dof(1) = NDOF
          fstrRESULT%node_label(1) = 'DISPLACEMENT'
          fstrRESULT%node_val_item = hecMAT%X
          call fstr2hecmw_mesh_conv(hecMESH)
          call hecmw_visualize_init
          call hecmw_visualize(hecMESH,fstrRESULT,istep,nstep,1)
          call hecmw_visualize_finalize
          call hecmw2fstr_mesh_conv(hecMESH)
          call hecmw_result_free(fstrRESULT)
        endif
      enddo

      IF(myrank == 0)THEN
        WRITE(IMSG,'("### FSTR_SOLVE_EIGEN FINISHED!")')
        WRITE(*,'("### FSTR_SOLVE_EIGEN FINISHED!")')
      ENDIF

      end subroutine fstr_solve_eigen

end module m_fstr_solve_eigen

