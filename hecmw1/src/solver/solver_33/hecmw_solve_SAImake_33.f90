!======================================================================!
!                                                                      !
!   Software Name : HEC-MW Library for PC-cluster                      !
!         Version : 2.6                                                !
!                                                                      !
!     Last Update : 2006/06/01                                         !
!        Category : Linear Solver                                      !
!                                                                      !
!            Written by Kengo Nakajima (Univ. of Tokyo)                !
!                                                                      !
!     Contact address :  IIS,The University of Tokyo RSS21 project     !
!                                                                      !
!     "Structural Analysis System for General-purpose Coupling         !
!      Simulations Using High End Computing Middleware (HEC-MW)"       !
!                                                                      !
!======================================================================!
module hecmw_solve_SAI_make_33
contains
!C
!C*** hecmw_solve_SAImake_33
!C
!C    construct sparse approximate inverse (SAI) preconditioner
!C    using LAPACK
!C

      subroutine hecmw_solve_SAImake_33                                 &
     &   (N, NP, NPL, AMAT, AMATsai, INDEX, ITEM, INDEXsai, ITEMsai,    &
     &    B, X,  SYM, THRESH, F_THRESH)

      use hecmw_util
      use m_hecmw_solve_error
      use m_hecmw_comm_f
      use hecmw_solver_SR_33

      implicit none

      integer(kind=kint ), intent(in):: N, NP, NPL, SYM

      real(kind=kreal), dimension(3*NP) , intent(inout) :: B, X
      real(kind=kreal), dimension(9*NPL), intent(inout) :: AMAT
      real(kind=kreal), dimension(9*NPL), intent(inout) :: AMATsai

      integer(kind=kint ), dimension(0:NP) ,intent(in) :: INDEX
      integer(kind=kint ), dimension(  NPL),intent(in) :: ITEM
      integer(kind=kint ), dimension(0:NP) ,intent(inout) :: INDEXsai
      integer(kind=kint ), dimension(  NPL),intent(inout) :: ITEMsai

      integer (kind=kint), dimension(:), allocatable :: ITEMsai_add
      integer (kind=kint), dimension(:), allocatable :: DIAGpos
      integer (kind=kint), dimension(:), allocatable :: WK1, WK2, WK3

      real (kind=kreal), dimension(:),  allocatable :: AHAT , BHAT, WORK
      real (kind=kreal) :: THRESH, F_THRESH

      integer (kind=kint) :: i,j,jj,k,kk,kk1,kk2,ik,jk,iS,iE,iiS
      integer (kind=kint) :: NWORK, IFLAG, info
      integer (kind=kint) :: kdiag_a, kdiag_b,in0
      integer (kind=kint) :: icou, icou0,iadd, length, length3, nn1,npat,npat3
      real (kind=kreal)   :: Da1,Da2,Da3,Db1,Db2,Db3
      real (kind=kreal)   :: rDD11,rDD12,rDD13,rDD21,rDD22,rDD23,rDD31,rDD32,rDD33
      real (kind=kreal)   :: eps11,eps12,eps13,eps21,eps22,eps23,eps31,eps32,eps33
      real (kind=kreal)   :: a11,a12,a13,a21,a22,a23,a31,a32,a33
      real (kind=kreal)   :: b11,b12,b13,b21,b22,b23,b31,b32,b33
      real (kind=kreal)   :: S_time
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      NWORK= 2000*64
      allocate (WORK(NWORK))
      WORK= 0.d0

      allocate (DIAGpos(NP))
      allocate (ITEMsai_add(NPL))
      DIAGpos   = 0
      ITEMsai_add= 0

      do i= 1, NP
        iS= INDEX(i-1) + 1
        iE= INDEX(i)
        do k= iS, iE
          if (ITEM(k).eq.i) DIAGpos(i)= k
        enddo
      enddo
!C===

!C
!C +----------------------------------------------+
!C | drop small components of the original matrix |
!C +----------------------------------------------+
!C===
      do i= 1, NP
        iS= INDEX(i-1) + 1
        iE= INDEX(i)
        kdiag_a= DIAGpos(i)

        Da1= dsqrt(dabs(AMAT(9*kdiag_a-8)))
        Da2= dsqrt(dabs(AMAT(9*kdiag_a-4)))
        Da3= dsqrt(dabs(AMAT(9*kdiag_a  )))

        icou= 0
        do k= iS, iE
          kk= ITEM(k)
          kdiag_b= DIAGpos(kk)
          Db1= dsqrt(dabs(AMAT(9*kdiag_b-8)))
          Db2= dsqrt(dabs(AMAT(9*kdiag_b-4)))
          Db3= dsqrt(dabs(AMAT(9*kdiag_b  )))

          rDD11= 1.d0 / (Da1*Db1)
          rDD12= 1.d0 / (Da1*Db2)
          rDD13= 1.d0 / (Da1*Db3)
          rDD21= 1.d0 / (Da2*Db1)
          rDD22= 1.d0 / (Da2*Db2)
          rDD23= 1.d0 / (Da2*Db3)
          rDD31= 1.d0 / (Da3*Db1)
          rDD32= 1.d0 / (Da3*Db2)
          rDD33= 1.d0 / (Da3*Db3)

          eps11= dabs(AMAT(9*k-8)) * rDD11
          eps12= dabs(AMAT(9*k-7)) * rDD12
          eps13= dabs(AMAT(9*k-6)) * rDD13
          eps21= dabs(AMAT(9*k-5)) * rDD21
          eps22= dabs(AMAT(9*k-4)) * rDD22
          eps23= dabs(AMAT(9*k-3)) * rDD23
          eps31= dabs(AMAT(9*k-2)) * rDD31
          eps32= dabs(AMAT(9*k-1)) * rDD32
          eps33= dabs(AMAT(9*k  )) * rDD33

          IFLAG= 0
          if (eps11.ge.THRESH .or. eps12.ge.THRESH .or. eps13.ge.THRESH &
     &   .or. eps21.ge.THRESH .or. eps22.ge.THRESH .or. eps23.ge.THRESH &
     &   .or. eps31.ge.THRESH .or. eps32.ge.THRESH .or. eps33.ge.THRESH)&
     &        IFLAG= 1

          if (k.eq.kdiag_a) IFLAG= 1
          if (IFLAG.eq.1) then
            icou= icou + 1
            iadd= icou + INDEXsai(i-1)
            INDEXsai(i)       = iadd
            ITEMsai    (iadd)= ITEM(k)
            ITEMsai_add(iadd)= k
          endif
        enddo
      enddo
!C===

!C
!C +----------------------------------+
!C | minimize |MA-E| by least squares |
!C +----------------------------------+
!C===
      allocate (WK1(NP), WK2(NP), WK3(NP))

      do i= 1, N
        WK1= 0
        WK2= 0
        WK3= 0
!C
!C-- matrix size
        length= INDEXsai(i) - INDEXsai(i-1)
        do k= INDEXsai(i-1)+1, INDEXsai(i)
          WK1(ITEMsai(k))= 1
        enddo

        do j= 1, NP
          if (WK1(j).eq.1) then
            do jk= INDEXsai(j-1)+1, INDEXsai(j)
              nn1= ITEMsai(jk)
              if (WK1(nn1).eq.0) WK1(nn1)= 2
            enddo
          endif
        enddo

        icou= 0
        do j= 1, NP
          if (WK1(j).eq.1) then
            icou  = icou + 1
            WK2(j   )= icou
          endif
        enddo

        do j= 1, NP
          if (WK1(j).eq.2) then
            icou  = icou + 1
            WK2(j   )= icou
          endif
        enddo

        npat= icou

!C
!C-- CONSTRUCT MATRIX
        npat3  = npat   * 3
        length3= length * 3
        allocate (AHAT (npat3*length3), BHAT (npat3*3))
        AHAT = 0.d0
        BHAT = 0.d0

        in0= WK2(i)
        DIAGpos(i)= INDEXsai(i-1) + in0

        BHAT(          3*in0-2)= 1.d0
        BHAT(  npat3 + 3*in0-1)= 1.d0
        BHAT(2*npat3 + 3*in0  )= 1.d0

        do k= INDEXsai(i-1)+1, INDEXsai(i)
          j= ITEMsai(k)
          do jk= INDEXsai(j-1)+1, INDEXsai(j)
            kk1= ITEMsai_add(jk )
            kk2= ITEMsai    (jk)
            in0= WK2(kk2)
            a11= AMAT(9*kk1-8)
            a12= AMAT(9*kk1-7)
            a13= AMAT(9*kk1-6)
            a21= AMAT(9*kk1-5)
            a22= AMAT(9*kk1-4)
            a23= AMAT(9*kk1-3)
            a31= AMAT(9*kk1-2)
            a32= AMAT(9*kk1-1)
            a33= AMAT(9*kk1  )

            iiS= (k-INDEXsai(i-1)-1) * npat3 * 3
            AHAT(        3*in0-2+iiS)= a11
            AHAT(        3*in0-1+iiS)= a12
            AHAT(        3*in0  +iiS)= a13
            AHAT(  npat3+3*in0-2+iiS)= a21
            AHAT(  npat3+3*in0-1+iiS)= a22
            AHAT(  npat3+3*in0  +iiS)= a23
            AHAT(2*npat3+3*in0-2+iiS)= a31
            AHAT(2*npat3+3*in0-1+iiS)= a32
            AHAT(2*npat3+3*in0  +iiS)= a33
          enddo
        enddo

        call DGELS ('N', npat3, length3, 3, AHAT, npat3, BHAT, npat3,   &
     &               WORK, NWORK, info)

        do k= INDEXsai(i-1)+1, INDEXsai(i)
          j= ITEMsai(k)
          icou0= WK2(j)
          AMATsai(9*k-8)= BHAT(3*icou0-2          )
          AMATsai(9*k-7)= BHAT(3*icou0-1          )
          AMATsai(9*k-6)= BHAT(3*icou0            )
          AMATsai(9*k-5)= BHAT(3*icou0-2 +   npat3)
          AMATsai(9*k-4)= BHAT(3*icou0-1 +   npat3)
          AMATsai(9*k-3)= BHAT(3*icou0   +   npat3)
          AMATsai(9*k-2)= BHAT(3*icou0-2 + 2*npat3)
          AMATsai(9*k-1)= BHAT(3*icou0-1 + 2*npat3)
          AMATsai(9*k  )= BHAT(3*icou0   + 2*npat3)
        enddo

        deallocate (AHAT , BHAT)
      enddo
!C===

        S_time= HECMW_WTIME()
!C
!C +-----------------------+
!C | DROP small components |
!C +-----------------------+
!C===
      do i= 1, N
        iS= INDEXsai(i-1) + 1
        iE= INDEXsai(i)
        kdiag_a= DIAGpos(i)

        Da1= dsqrt(dabs(AMATsai(9*kdiag_a-8)))
        Da2= dsqrt(dabs(AMATsai(9*kdiag_a-4)))
        Da3= dsqrt(dabs(AMATsai(9*kdiag_a  )))

        icou= 0
        do k= iS, iE
          kk= ITEMsai(k)
          kdiag_b= DIAGpos(kk)
          Db1= dsqrt(dabs(AMATsai(9*kdiag_b-8)))
          Db2= dsqrt(dabs(AMATsai(9*kdiag_b-4)))
          Db3= dsqrt(dabs(AMATsai(9*kdiag_b  )))

          rDD11= 1.d0 / (Da1*Db1)
          rDD12= 1.d0 / (Da1*Db2)
          rDD13= 1.d0 / (Da1*Db3)
          rDD21= 1.d0 / (Da2*Db1)
          rDD22= 1.d0 / (Da2*Db2)
          rDD23= 1.d0 / (Da2*Db3)
          rDD31= 1.d0 / (Da3*Db1)
          rDD32= 1.d0 / (Da3*Db2)
          rDD33= 1.d0 / (Da3*Db3)

          eps11= dabs(AMATsai(9*k-8)) * rDD11
          eps12= dabs(AMATsai(9*k-7)) * rDD12
          eps13= dabs(AMATsai(9*k-6)) * rDD13
          eps21= dabs(AMATsai(9*k-5)) * rDD21
          eps22= dabs(AMATsai(9*k-4)) * rDD22
          eps23= dabs(AMATsai(9*k-3)) * rDD23
          eps31= dabs(AMATsai(9*k-2)) * rDD31
          eps32= dabs(AMATsai(9*k-1)) * rDD32
          eps33= dabs(AMATsai(9*k  )) * rDD33

          IFLAG= 0
          if (eps11.ge.F_THRESH .or. eps12.ge.F_THRESH .or.             &
     &        eps13.ge.F_THRESH .or. eps21.ge.F_THRESH .or.             &
     &        eps22.ge.F_THRESH .or. eps23.ge.F_THRESH .or.             &
     &        eps31.ge.F_THRESH .or. eps32.ge.F_THRESH .or.             &
     &        eps33.ge.F_THRESH)  IFLAG= 1

          if (k.eq.kdiag_a) IFLAG= 1

          if (IFLAG.eq.0) then
            AMATsai(9*k-8)= 0.d0
            AMATsai(9*k-7)= 0.d0
            AMATsai(9*k-6)= 0.d0
            AMATsai(9*k-5)= 0.d0
            AMATsai(9*k-4)= 0.d0
            AMATsai(9*k-3)= 0.d0
            AMATsai(9*k-2)= 0.d0
            AMATsai(9*k-1)= 0.d0
            AMATsai(9*k  )= 0.d0
          endif
        enddo
      enddo
!C===
      deallocate (WK1, WK2, WK3, ITEMsai_add, DIAGpos)

      if (SYM.eq.1) then
!C
!C +-------------+
!C | SYMMETRIC ? |
!C +-------------+
!C===
      do i= 1, N
        do ik= INDEXsai(i-1)+1, INDEXsai(i)
          j= ITEMsai(ik)
          if (j.gt.i) then
            do jk= INDEXsai(j-1)+1, INDEXsai(j)
              jj= ITEMsai(jk)
              if (jj.eq.i) then
                a11= AMATsai(9*ik-8)
                a12= AMATsai(9*ik-7)
                a13= AMATsai(9*ik-6)
                a21= AMATsai(9*ik-5)
                a22= AMATsai(9*ik-4)
                a23= AMATsai(9*ik-3)
                a31= AMATsai(9*ik-2)
                a32= AMATsai(9*ik-1)
                a33= AMATsai(9*ik  )
                b11= AMATsai(9*jk-8)
                b12= AMATsai(9*jk-5)
                b13= AMATsai(9*jk-2)
                b21= AMATsai(9*jk-7)
                b22= AMATsai(9*jk-4)
                b23= AMATsai(9*jk-1)
                b31= AMATsai(9*jk-6)
                b32= AMATsai(9*jk-3)
                b33= AMATsai(9*jk  )
                AMATsai(9*ik-8)= 0.50d0 * ( a11 + b11 )
                AMATsai(9*ik-7)= 0.50d0 * ( a12 + b12 )
                AMATsai(9*ik-6)= 0.50d0 * ( a13 + b13 )
                AMATsai(9*ik-5)= 0.50d0 * ( a21 + b21 )
                AMATsai(9*ik-4)= 0.50d0 * ( a22 + b22 )
                AMATsai(9*ik-3)= 0.50d0 * ( a23 + b23 )
                AMATsai(9*ik-2)= 0.50d0 * ( a31 + b31 )
                AMATsai(9*ik-1)= 0.50d0 * ( a32 + b32 )
                AMATsai(9*ik  )= 0.50d0 * ( a33 + b33 )
                AMATsai(9*jk-8)= 0.50d0 * ( a11 + b11 )
                AMATsai(9*jk-5)= 0.50d0 * ( a12 + b12 )
                AMATsai(9*jk-2)= 0.50d0 * ( a13 + b13 )
                AMATsai(9*jk-7)= 0.50d0 * ( a21 + b21 )
                AMATsai(9*jk-4)= 0.50d0 * ( a22 + b22 )
                AMATsai(9*jk-1)= 0.50d0 * ( a23 + b23 )
                AMATsai(9*jk-6)= 0.50d0 * ( a31 + b31 )
                AMATsai(9*jk-3)= 0.50d0 * ( a32 + b32 )
                AMATsai(9*jk  )= 0.50d0 * ( a33 + b33 )
              endif
            enddo
          endif
        enddo
      enddo
      endif
!C===

      end subroutine hecmw_solve_SAImake_33
end module hecmw_solve_SAI_make_33
